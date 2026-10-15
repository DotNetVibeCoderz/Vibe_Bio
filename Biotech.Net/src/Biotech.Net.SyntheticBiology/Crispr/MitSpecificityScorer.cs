namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// The MIT (crispr.mit.edu) off-target scoring scheme of Hsu et al., <i>Nature Biotechnology</i> 31:827 (2013).
/// </summary>
/// <remarks>
/// <para>
/// A single off-target hit scores
/// <c>100 · Π(1 − M[i]) · 1/(((19 − d̄)/19)·4 + 1) · 1/n²</c>, where <c>M</c> is the published vector of
/// position-specific mismatch penalties (index 0 is PAM-distal), <c>n</c> is the number of mismatches and
/// <c>d̄</c> the mean distance between consecutive mismatches. The two correction terms are 1 when <c>n &lt; 2</c>.
/// </para>
/// <para>
/// A guide's aggregate score is <c>100 · 100/(100 + Σ hit scores)</c> over all off-targets. Higher is more
/// specific; the usual thresholds are &gt; 50 for a usable guide and &gt; 80 for a good one.
/// </para>
/// </remarks>
public static class MitSpecificityScorer
{
    /// <summary>
    /// The published position weights, index 0 being the PAM-distal end of a 20 nt spacer. A mismatch at a
    /// weight of 0 is free; the PAM-proximal "seed" positions carry the largest penalties.
    /// </summary>
    public static ReadOnlySpan<double> PositionWeights =>
    [
        0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
        0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583,
    ];

    /// <summary>
    /// Score (0–100) of one off-target site, given the 0-based positions of the mismatches within a 20 nt
    /// spacer. An empty list (a perfect match) scores 100.
    /// </summary>
    public static double HitScore(ReadOnlySpan<int> mismatchPositions)
    {
        var weights = PositionWeights;
        double score = 1.0;
        for (int i = 0; i < mismatchPositions.Length; i++)
        {
            int p = mismatchPositions[i];
            if ((uint)p >= (uint)weights.Length) throw new ArgumentOutOfRangeException(nameof(mismatchPositions));
            score *= 1.0 - weights[p];
        }

        int n = mismatchPositions.Length;
        if (n > 1)
        {
            // Mean distance between the outermost mismatches, as in the published implementation.
            double meanDistance = (double)(mismatchPositions[n - 1] - mismatchPositions[0]) / (n - 1);
            score *= 1.0 / (((19.0 - meanDistance) / 19.0) * 4.0 + 1.0);
            score *= 1.0 / (n * (double)n);
        }

        return score * 100.0;
    }

    /// <summary>
    /// Score (0–100) of one off-target site from the two 20 nt sequences: the spacer and the off-target
    /// protospacer, both 5′→3′ on the targeted strand.
    /// </summary>
    public static double HitScore(ReadOnlySpan<char> spacer, ReadOnlySpan<char> offTarget)
    {
        if (spacer.Length != offTarget.Length)
            throw new ArgumentException("The spacer and the off-target must be the same length.", nameof(offTarget));

        int skip = Math.Max(0, spacer.Length - 20);
        Span<int> positions = stackalloc int[20];
        int n = 0;
        for (int i = skip; i < spacer.Length && n < 20; i++)
        {
            char a = char.ToUpperInvariant(spacer[i]);
            char b = char.ToUpperInvariant(offTarget[i]);
            if (a == 'U') a = 'T';
            if (b == 'U') b = 'T';
            if (a != b) positions[n++] = i - skip;
        }
        return HitScore(positions[..n]);
    }

    /// <summary>
    /// Aggregate guide specificity (0–100) from a set of off-target hits. The intended site is excluded.
    /// </summary>
    public static double GuideScore(IEnumerable<OffTarget> offTargets)
    {
        ArgumentNullException.ThrowIfNull(offTargets);
        double sum = 0;
        foreach (var hit in offTargets)
        {
            if (hit.IsOnTarget) continue;
            sum += hit.MitScore;
        }
        return 100.0 * 100.0 / (100.0 + sum);
    }
}
