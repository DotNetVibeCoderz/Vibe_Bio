using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>Filters applied while scanning a target for candidate guides.</summary>
public sealed class GuideSearchOptions
{
    /// <summary>The nuclease to design for. Defaults to <see cref="CasNuclease.SpCas9"/>.</summary>
    public CasNuclease Nuclease { get; init; } = CasNuclease.SpCas9;

    /// <summary>Search the plus strand.</summary>
    public bool SearchForward { get; init; } = true;

    /// <summary>Search the minus strand.</summary>
    public bool SearchReverse { get; init; } = true;

    /// <summary>Minimum spacer GC fraction (0–1).</summary>
    public double MinGcContent { get; init; } = 0.20;

    /// <summary>Maximum spacer GC fraction (0–1).</summary>
    public double MaxGcContent { get; init; } = 0.85;

    /// <summary>Reject spacers containing <c>TTTT</c>, which terminates U6 transcription.</summary>
    public bool RejectPolyT { get; init; } = true;

    /// <summary>Reject spacers with a single-base run longer than this (0 disables the check).</summary>
    public int MaxHomopolymer { get; init; } = 6;

    /// <summary>Reject spacers that contain an ambiguity code such as N.</summary>
    public bool RequireUnambiguous { get; init; } = true;

    /// <summary>Require a 5′ G, as U6-driven transcription prefers (the designer can add one instead).</summary>
    public bool RequireLeadingG { get; init; }

    /// <summary>Compute the Doench 2014 on-target score for every candidate.</summary>
    public bool ScoreOnTarget { get; init; } = true;

    /// <summary>Discard guides whose on-target score is below this (ignored when the score is unavailable).</summary>
    public double MinOnTargetScore { get; init; }

    /// <summary>The default options: SpCas9, both strands, the usual composition filters.</summary>
    public static GuideSearchOptions Default { get; } = new();
}

/// <summary>
/// Scans a target sequence for PAM sites and returns the spacers that pass the composition filters, on both
/// strands.
/// </summary>
/// <remarks>
/// The scan is a single left-to-right pass over the ASCII buffer with an IUPAC PAM match at each offset; it
/// allocates only the candidate list and the spacer strings. A reverse-strand hit is found by matching the
/// reverse complement of the PAM pattern against the plus strand, so the sequence is never copied or
/// reverse-complemented as a whole.
/// </remarks>
public static class GuideFinder
{
    /// <summary>Finds guides on a <see cref="Sequence"/>.</summary>
    public static IReadOnlyList<Guide> Find(Sequence target, GuideSearchOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        return Find(target.Span, target.Id, 0, options);
    }

    /// <summary>Finds guides inside a region of a <see cref="Sequence"/>, with flanking context for scoring.</summary>
    /// <remarks>
    /// Only guides whose <em>cut site</em> falls inside <c>[start, start + length)</c> are returned, but the
    /// scoring context is taken from the full sequence, so a guide at the edge of the region still scores.
    /// </remarks>
    public static IReadOnlyList<Guide> FindInRegion(Sequence target, int start, int length, GuideSearchOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (start + length > target.Length) throw new ArgumentOutOfRangeException(nameof(length));

        var all = Find(target.Span, target.Id, 0, options);
        var filtered = new List<Guide>();
        foreach (var guide in all)
        {
            int cut = guide.CutSite;
            if (cut >= start && cut <= start + length) filtered.Add(guide);
        }
        return filtered;
    }

    /// <summary>Finds guides on a plain DNA string.</summary>
    public static IReadOnlyList<Guide> Find(string target, GuideSearchOptions? options = null, string sequenceId = "target")
    {
        ArgumentNullException.ThrowIfNull(target);
        return Find(DnaUtility.Ascii(target), sequenceId, 0, options);
    }

    /// <summary>
    /// Finds guides on an ASCII span. <paramref name="offset"/> is added to every reported coordinate, so a
    /// caller can scan a slice of a chromosome and still report genomic positions.
    /// </summary>
    public static IReadOnlyList<Guide> Find(ReadOnlySpan<byte> target, string sequenceId, int offset = 0, GuideSearchOptions? options = null)
    {
        options ??= GuideSearchOptions.Default;
        var cas = options.Nuclease;
        var results = new List<Guide>();

        int pamLength = cas.PamLength;
        int guideLength = cas.GuideLength;
        var pam = cas.PamSpan;

        // Reverse complement of the PAM pattern, so a minus-strand PAM can be matched against the plus strand.
        Span<byte> pamRc = stackalloc byte[pamLength];
        for (int i = 0; i < pamLength; i++) pamRc[i] = DnaUtility.Complement(pam[pamLength - 1 - i]);

        for (int p = 0; p + pamLength <= target.Length; p++)
        {
            var window = target.Slice(p, pamLength);

            if (options.SearchForward && DnaUtility.MatchesIupac(window, pam))
            {
                int protoStart = cas.PamSide == PamSide.ThreePrime ? p - guideLength : p + pamLength;
                if (protoStart >= 0 && protoStart + guideLength <= target.Length)
                    TryAdd(results, target, sequenceId, offset, protoStart, Strand.Forward, options);
            }

            if (options.SearchReverse && DnaUtility.MatchesIupac(window, pamRc))
            {
                int protoStart = cas.PamSide == PamSide.ThreePrime ? p + pamLength : p - guideLength;
                if (protoStart >= 0 && protoStart + guideLength <= target.Length)
                    TryAdd(results, target, sequenceId, offset, protoStart, Strand.Reverse, options);
            }
        }

        return results;
    }

    private static void TryAdd(List<Guide> results, ReadOnlySpan<byte> target, string sequenceId, int offset,
        int protoStart, Strand strand, GuideSearchOptions options)
    {
        var cas = options.Nuclease;
        var protospacer = target.Slice(protoStart, cas.GuideLength);

        if (options.RequireUnambiguous && !DnaUtility.IsUnambiguous(protospacer)) return;

        double gc = DnaUtility.GcFraction(protospacer);
        if (gc < options.MinGcContent || gc > options.MaxGcContent) return;

        string spacer = strand == Strand.Reverse
            ? DnaUtility.ReverseComplement(protospacer)
            : Encoding.ASCII.GetString(protospacer).ToUpperInvariant();

        if (options.RejectPolyT && spacer.Contains("TTTT", StringComparison.Ordinal)) return;
        if (options.MaxHomopolymer > 0 && DnaUtility.LongestHomopolymer(DnaUtility.Ascii(spacer)) > options.MaxHomopolymer) return;
        if (options.RequireLeadingG && !spacer.StartsWith('G')) return;

        // The PAM as read on the targeted strand.
        int pamStartPlus = cas.PamSide == PamSide.ThreePrime
            ? (strand == Strand.Reverse ? protoStart - cas.PamLength : protoStart + cas.GuideLength)
            : (strand == Strand.Reverse ? protoStart + cas.GuideLength : protoStart - cas.PamLength);
        if (pamStartPlus < 0 || pamStartPlus + cas.PamLength > target.Length) return;
        var pamPlus = target.Slice(pamStartPlus, cas.PamLength);
        string pamText = strand == Strand.Reverse
            ? DnaUtility.ReverseComplement(pamPlus)
            : Encoding.ASCII.GetString(pamPlus).ToUpperInvariant();

        string? context = ExtractContext(target, protoStart, strand, cas);

        var guide = new Guide(cas, spacer, pamText, sequenceId, protoStart + offset, strand, context);
        if (options.ScoreOnTarget)
        {
            guide.OnTargetScore = OnTargetScorer.Score(guide);
            if (!double.IsNaN(guide.OnTargetScore) && guide.OnTargetScore < options.MinOnTargetScore) return;
        }
        results.Add(guide);
    }

    /// <summary>
    /// The 30-mer Doench context on the targeted strand: 4 nt upstream of the spacer, the spacer, the 3 nt PAM
    /// and 3 nt downstream. Null when the model does not apply or the flanks run off the sequence.
    /// </summary>
    private static string? ExtractContext(ReadOnlySpan<byte> target, int protoStart, Strand strand, CasNuclease cas)
    {
        if (cas.GuideLength != 20 || cas.PamSide != PamSide.ThreePrime || cas.PamLength != 3) return null;

        // On the plus strand the 30-mer runs [protoStart - 4, protoStart + 26); on the minus strand it is the
        // reverse complement of [protoStart - 6, protoStart + 24).
        int windowStart = strand == Strand.Reverse ? protoStart - 6 : protoStart - 4;
        if (windowStart < 0 || windowStart + OnTargetScorer.ContextLength > target.Length) return null;

        var window = target.Slice(windowStart, OnTargetScorer.ContextLength);
        if (!DnaUtility.IsUnambiguous(window)) return null;

        return strand == Strand.Reverse
            ? DnaUtility.ReverseComplement(window)
            : Encoding.ASCII.GetString(window).ToUpperInvariant();
    }
}
