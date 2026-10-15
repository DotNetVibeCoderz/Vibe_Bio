using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>Whether a bulge is an extra base in the DNA or in the guide RNA.</summary>
public enum BulgeKind
{
    /// <summary>No bulge; the alignment is mismatches only.</summary>
    None,

    /// <summary>An extra base in the genomic DNA, looped out against the guide.</summary>
    Dna,

    /// <summary>An extra base in the guide RNA, looped out against the DNA.</summary>
    Rna,
}

/// <summary>
/// A genomic site a guide can bind besides its intended target.
/// </summary>
/// <param name="SequenceId">Chromosome or contig the site is on.</param>
/// <param name="Start">0-based plus-strand start of the protospacer.</param>
/// <param name="Strand">Strand the guide would anneal to.</param>
/// <param name="Protospacer">The off-target protospacer, 5′→3′ on the targeted strand.</param>
/// <param name="Pam">The off-target PAM, 5′→3′ on the targeted strand.</param>
/// <param name="MismatchCount">Number of mismatched positions.</param>
/// <param name="Bulge">Bulge type, if the alignment needed one.</param>
/// <param name="CfdScore">Doench 2016 CFD score (0–1).</param>
/// <param name="MitScore">Hsu 2013 single-hit score (0–100).</param>
/// <param name="IsOnTarget">True for the intended site (zero mismatches at the designed locus).</param>
public sealed record OffTarget(
    string SequenceId,
    int Start,
    Strand Strand,
    string Protospacer,
    string Pam,
    int MismatchCount,
    BulgeKind Bulge,
    double CfdScore,
    double MitScore,
    bool IsOnTarget)
{
    /// <summary>The site as a half-open interval.</summary>
    public GenomicInterval Interval => new(SequenceId, Start, Start + Protospacer.Length, Strand);

    /// <summary>A mismatch mask: <c>.</c> where the guide matches, the off-target base where it does not.</summary>
    public string Alignment(string spacer)
    {
        ArgumentNullException.ThrowIfNull(spacer);
        int n = Math.Min(spacer.Length, Protospacer.Length);
        return string.Create(n, (spacer, Protospacer), static (span, state) =>
        {
            var (guide, site) = state;
            for (int i = 0; i < span.Length; i++)
                span[i] = char.ToUpperInvariant(guide[i]) == char.ToUpperInvariant(site[i]) ? '.' : char.ToUpperInvariant(site[i]);
        });
    }

    /// <inheritdoc/>
    public override string ToString()
        => $"{Protospacer} {Pam} {SequenceId}:{Start + 1}({(Strand == Strand.Reverse ? '-' : '+')}) "
           + $"{MismatchCount}mm{(Bulge == BulgeKind.None ? "" : $"+{Bulge}bulge")} CFD {CfdScore:F3}";
}

/// <summary>Controls an off-target search.</summary>
public sealed class OffTargetSearchOptions
{
    /// <summary>Maximum number of mismatches to report. The seed index supports up to 5.</summary>
    public int MaxMismatches { get; init; } = 4;

    /// <summary>
    /// Maximum number of single-base bulges (0 or 1). Bulge search is a separate, slower pass over every PAM
    /// site, so leave it off for genome-wide scans unless you need it.
    /// </summary>
    public int MaxBulges { get; init; }

    /// <summary>
    /// Also search the relaxed NAG PAM, which SpCas9 cleaves at a few percent of NGG efficiency and which the
    /// CFD PAM table scores explicitly.
    /// </summary>
    public bool IncludeNagPam { get; init; } = true;

    /// <summary>Stop after this many hits (0 = unlimited). Hits are not ordered while the cap applies.</summary>
    public int MaxResults { get; init; }

    /// <summary>Sort the results by descending CFD score.</summary>
    public bool SortByCfd { get; init; } = true;

    /// <summary>The default options: up to 4 mismatches, no bulges, NAG included.</summary>
    public static OffTargetSearchOptions Default { get; } = new();
}

/// <summary>A guide's off-target profile and the two aggregate specificity scores computed from it.</summary>
/// <param name="Spacer">The guide that was searched.</param>
/// <param name="Hits">Every site found, including the intended one.</param>
/// <param name="MitSpecificity">Hsu 2013 aggregate specificity (0–100).</param>
/// <param name="CfdSpecificity">CFD aggregate specificity (0–100).</param>
public sealed record OffTargetProfile(string Spacer, IReadOnlyList<OffTarget> Hits, double MitSpecificity, double CfdSpecificity)
{
    /// <summary>Off-target sites only (the intended site removed).</summary>
    public IEnumerable<OffTarget> OffTargets => Hits.Where(h => !h.IsOnTarget);

    /// <summary>Number of off-targets with exactly <paramref name="mismatches"/> mismatches.</summary>
    public int CountWith(int mismatches) => OffTargets.Count(h => h.MismatchCount == mismatches);

    /// <summary>Counts of off-targets by mismatch number, index 0 … <c>MaxMismatches</c>.</summary>
    public int[] MismatchHistogram(int maxMismatches)
    {
        var histogram = new int[maxMismatches + 1];
        foreach (var hit in OffTargets)
            if (hit.MismatchCount <= maxMismatches) histogram[hit.MismatchCount]++;
        return histogram;
    }

    /// <inheritdoc/>
    public override string ToString()
        => $"{Spacer}: {Hits.Count - 1} off-targets, MIT {MitSpecificity:F1}, CFD {CfdSpecificity:F1}";
}
