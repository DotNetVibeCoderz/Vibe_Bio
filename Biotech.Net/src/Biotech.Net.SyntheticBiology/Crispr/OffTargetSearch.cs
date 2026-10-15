using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// The entry point for off-target discovery: runs a spacer against a <see cref="CrisprSiteIndex"/>, scores every
/// hit with CFD and MIT, and aggregates the two published specificity scores.
/// </summary>
/// <remarks>
/// <para>
/// The mismatch pass is the seeded index search of <see cref="CrisprSiteIndex.Query(string, int)"/>: exact, and
/// fast enough for a genome. The optional bulge pass is a separate linear scan over the indexed sequences,
/// because a looped-out base shifts every downstream position and so defeats the fixed-width seed buckets. It
/// costs roughly one pass per indexed base and is off by default.
/// </para>
/// <para>
/// <b>Scoring of bulge hits.</b> The published CFD and MIT tables cover substitutions only; neither paper
/// provides a penalty for an RNA or DNA bulge. A bulge hit is therefore scored by collapsing the bulge (the
/// looped-out base is removed and the remaining positions are scored as substitutions) and multiplying by
/// <see cref="BulgePenalty"/>, a single blunt constant chosen here, not a published parameter. For the same
/// reason bulge hits are excluded from the aggregate <see cref="OffTargetProfile.MitSpecificity"/> and
/// <see cref="OffTargetProfile.CfdSpecificity"/>, which are defined in the literature over substitution hits
/// only — mixing the two would produce a number that looks like CRISPOR's but is not comparable to it.
/// </para>
/// </remarks>
public static class OffTargetSearch
{
    /// <summary>
    /// The multiplier applied to a bulge hit's collapsed CFD/MIT score. A blunt in-house constant reflecting
    /// that bulged binding is markedly less efficient than a substitution-only mismatch; it is not a published
    /// value and should not be read as one.
    /// </summary>
    public const double BulgePenalty = 0.2;

    /// <summary>
    /// Searches <paramref name="index"/> for every site <paramref name="spacer"/> can bind, and returns the
    /// scored profile.
    /// </summary>
    /// <param name="spacer">The protospacer, 5′→3′ on the targeted strand, of the index nuclease's length.</param>
    /// <param name="index">The site index to search.</param>
    /// <param name="options">Search limits; defaults to <see cref="OffTargetSearchOptions.Default"/>.</param>
    /// <returns>Every hit found, with the two aggregate specificity scores.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <see cref="OffTargetSearchOptions.MaxMismatches"/> exceeds what the seed index can answer exactly, or
    /// <see cref="OffTargetSearchOptions.MaxBulges"/> is not 0 or 1.
    /// </exception>
    public static OffTargetProfile Search(string spacer, CrisprSiteIndex index, OffTargetSearchOptions? options = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(spacer);
        ArgumentNullException.ThrowIfNull(index);
        options ??= OffTargetSearchOptions.Default;
        if (options.MaxBulges is < 0 or > 1)
            throw new ArgumentOutOfRangeException(nameof(options), "At most one bulge is modelled.");

        string normalised = DnaUtility.NormalizeDna(spacer);
        var hits = new List<OffTarget>(index.Query(normalised, options.MaxMismatches));

        if (options.MaxBulges == 1)
        {
            var seen = new HashSet<(string, int, Strand)>();
            foreach (var hit in hits) seen.Add((hit.SequenceId, hit.Start, hit.Strand));
            foreach (var bulgeHit in BulgeHits(normalised, index, options))
                if (seen.Add((bulgeHit.SequenceId, bulgeHit.Start, bulgeHit.Strand)))
                    hits.Add(bulgeHit);
        }

        MarkOnTarget(hits);

        if (options.SortByCfd)
            hits.Sort(static (a, b) =>
            {
                if (a.IsOnTarget != b.IsOnTarget) return a.IsOnTarget ? -1 : 1;
                double x = double.IsNaN(a.CfdScore) ? -1 : a.CfdScore;
                double y = double.IsNaN(b.CfdScore) ? -1 : b.CfdScore;
                int byScore = y.CompareTo(x);
                return byScore != 0 ? byScore : a.MismatchCount.CompareTo(b.MismatchCount);
            });

        if (options.MaxResults > 0 && hits.Count > options.MaxResults)
            hits.RemoveRange(options.MaxResults, hits.Count - options.MaxResults);

        var substitutionHits = hits.Where(static h => h.Bulge == BulgeKind.None).ToList();
        return new OffTargetProfile(
            normalised,
            hits,
            MitSpecificityScorer.GuideScore(substitutionHits),
            CfdScorer.GuideSpecificity(substitutionHits));
    }

    /// <summary>
    /// Searches an index and marks the hit at a known designed locus as the on-target one, rather than
    /// inferring it from the alignment.
    /// </summary>
    /// <param name="spacer">The protospacer, 5′→3′ on the targeted strand.</param>
    /// <param name="index">The site index to search.</param>
    /// <param name="sequenceId">Identifier of the sequence the guide was designed against.</param>
    /// <param name="start">0-based plus-strand start of the intended protospacer.</param>
    /// <param name="options">Search limits; defaults to <see cref="OffTargetSearchOptions.Default"/>.</param>
    /// <returns>Every hit found, with the designed locus flagged <see cref="OffTarget.IsOnTarget"/>.</returns>
    public static OffTargetProfile SearchAt(string spacer, CrisprSiteIndex index, string sequenceId, int start,
        OffTargetSearchOptions? options = null)
    {
        var profile = Search(spacer, index, options);
        var hits = new List<OffTarget>(profile.Hits.Count);
        foreach (var hit in profile.Hits)
            hits.Add(hit with { IsOnTarget = hit.Start == start && hit.SequenceId == sequenceId });

        var substitutionHits = hits.Where(static h => h.Bulge == BulgeKind.None).ToList();
        return new OffTargetProfile(
            profile.Spacer,
            hits,
            MitSpecificityScorer.GuideScore(substitutionHits),
            CfdScorer.GuideSpecificity(substitutionHits));
    }

    /// <summary>Convenience overload that searches a guide found by <see cref="GuideFinder"/> at its own locus.</summary>
    /// <param name="guide">The guide to profile.</param>
    /// <param name="index">The site index to search.</param>
    /// <param name="options">Search limits; defaults to <see cref="OffTargetSearchOptions.Default"/>.</param>
    /// <returns>The guide's off-target profile.</returns>
    public static OffTargetProfile Search(Guide guide, CrisprSiteIndex index, OffTargetSearchOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(guide);
        return SearchAt(guide.Spacer, index, guide.SequenceId, guide.Start, options);
    }

    /// <summary>
    /// The on-target hit is the perfect match carrying a canonical PAM. When the caller has not said which
    /// locus was designed against, the first such hit is taken; if there is none, nothing is marked.
    /// </summary>
    private static void MarkOnTarget(List<OffTarget> hits)
    {
        for (int i = 0; i < hits.Count; i++)
        {
            if (hits[i].MismatchCount != 0 || hits[i].Bulge != BulgeKind.None) continue;
            hits[i] = hits[i] with { IsOnTarget = true };
            return;
        }
    }

    /// <summary>
    /// Linear bulge pass: for every PAM-adjacent window, try a single looped-out base in the DNA (the genomic
    /// window is one base longer than the spacer) or in the guide RNA (one base shorter).
    /// </summary>
    private static List<OffTarget> BulgeHits(string spacer, CrisprSiteIndex index, OffTargetSearchOptions options)
    {
        var results = new List<OffTarget>();
        var cas = index.Nuclease;
        if (cas.PamSide != PamSide.ThreePrime) return results; // documented limitation

        int pamLength = cas.PamLength;
        var pams = new List<byte[]> { DnaUtility.Ascii(cas.Pam) };
        if (options.IncludeNagPam && cas.Pam == "NGG") pams.Add(DnaUtility.Ascii("NAG"));

        Span<byte> pamRc = stackalloc byte[pamLength];
        for (int s = 0; s < index.SequenceCount; s++)
        {
            var sequence = index.SequenceAt(s);
            string id = index.SequenceIds[s];

            foreach (var pam in pams)
            {
                for (int i = 0; i < pamLength; i++) pamRc[i] = DnaUtility.Complement(pam[pamLength - 1 - i]);

                for (int p = 0; p + pamLength <= sequence.Length; p++)
                {
                    var window = sequence.Slice(p, pamLength);

                    if (DnaUtility.MatchesIupac(window, pam))
                    {
                        AddBulgeCandidates(results, sequence, id, spacer, p, Strand.Forward, cas, options);
                    }
                    if (DnaUtility.MatchesIupac(window, pamRc))
                    {
                        AddBulgeCandidates(results, sequence, id, spacer, p, Strand.Reverse, cas, options);
                    }
                }
            }
        }

        return results;
    }

    private static void AddBulgeCandidates(List<OffTarget> results, ReadOnlySpan<byte> sequence, string sequenceId,
        string spacer, int pamStart, Strand strand, CasNuclease cas, OffTargetSearchOptions options)
    {
        int guideLength = cas.GuideLength;
        int pamLength = cas.PamLength;

        foreach (var kind in (ReadOnlySpan<BulgeKind>)[BulgeKind.Dna, BulgeKind.Rna])
        {
            int windowLength = kind == BulgeKind.Dna ? guideLength + 1 : guideLength - 1;
            int start = strand == Strand.Reverse ? pamStart + pamLength : pamStart - windowLength;
            if (start < 0 || start + windowLength > sequence.Length) continue;

            var raw = sequence.Slice(start, windowLength);
            if (!DnaUtility.IsUnambiguous(raw)) continue;

            string site = strand == Strand.Reverse
                ? DnaUtility.ReverseComplement(raw)
                : Encoding.ASCII.GetString(raw).ToUpperInvariant();

            if (!TryBestBulgeAlignment(spacer, site, kind, options.MaxMismatches, out int mismatches, out string collapsed))
                continue;

            string pam = PamText(sequence, pamStart, pamLength, strand);

            double cfd = BulgePenalty * CfdScorer.Score(spacer, collapsed, pam);
            double mit = BulgePenalty * MitSpecificityScorer.HitScore(spacer, collapsed);

            results.Add(new OffTarget(sequenceId, start, strand, site, pam, mismatches, kind, cfd, mit, IsOnTarget: false));
        }
    }

    private static string PamText(ReadOnlySpan<byte> sequence, int pamStart, int pamLength, Strand strand)
    {
        var window = sequence.Slice(pamStart, pamLength);
        return strand == Strand.Reverse
            ? DnaUtility.ReverseComplement(window)
            : Encoding.ASCII.GetString(window).ToUpperInvariant();
    }

    /// <summary>
    /// Best single-bulge alignment of <paramref name="spacer"/> against <paramref name="site"/>: every possible
    /// looped-out position is tried and the one with the fewest substitutions wins.
    /// </summary>
    /// <param name="spacer">The guide spacer.</param>
    /// <param name="site">The genomic window: one base longer than the spacer for a DNA bulge, one shorter for an RNA bulge.</param>
    /// <param name="kind">Which strand carries the extra base.</param>
    /// <param name="maxMismatches">Substitutions allowed on top of the bulge.</param>
    /// <param name="mismatches">The number of substitutions in the best alignment.</param>
    /// <param name="collapsed">
    /// The site with the bulge removed, aligned one-to-one with the spacer, for substitution-based scoring.
    /// </param>
    /// <returns>True when an alignment within <paramref name="maxMismatches"/> exists.</returns>
    internal static bool TryBestBulgeAlignment(string spacer, string site, BulgeKind kind, int maxMismatches,
        out int mismatches, out string collapsed)
    {
        mismatches = int.MaxValue;
        collapsed = string.Empty;

        if (kind == BulgeKind.Dna)
        {
            // The site has one base too many; drop one of its positions.
            for (int skip = 0; skip < site.Length; skip++)
            {
                int count = 0;
                var buffer = new char[spacer.Length];
                for (int i = 0, j = 0; i < site.Length; i++)
                {
                    if (i == skip) continue;
                    buffer[j] = site[i];
                    if (char.ToUpperInvariant(site[i]) != char.ToUpperInvariant(spacer[j])) count++;
                    j++;
                }
                if (count < mismatches) { mismatches = count; collapsed = new string(buffer); }
            }
        }
        else
        {
            // The guide has one base too many; drop one of its positions and pad the site to match.
            for (int skip = 0; skip < spacer.Length; skip++)
            {
                int count = 0;
                var buffer = new char[spacer.Length];
                for (int i = 0, j = 0; i < spacer.Length; i++)
                {
                    if (i == skip)
                    {
                        // The looped-out guide base pairs with nothing; score it as a mismatch-neutral filler
                        // that never matches, which is what the collapsed scoring needs.
                        buffer[i] = char.ToUpperInvariant(spacer[i]);
                        continue;
                    }
                    buffer[i] = site[j];
                    if (char.ToUpperInvariant(site[j]) != char.ToUpperInvariant(spacer[i])) count++;
                    j++;
                }
                if (count < mismatches) { mismatches = count; collapsed = new string(buffer); }
            }
        }

        return mismatches <= maxMismatches;
    }
}
