namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// Cutting Frequency Determination (CFD) off-target scoring, Doench et al., <i>Nature Biotechnology</i> 34:184
/// (2016), Supplementary Table 19.
/// </summary>
/// <remarks>
/// <para>
/// CFD is the product of one published penalty per mismatched position and one penalty for the PAM: a perfect
/// match with an NGG PAM scores 1.0, and every mismatch multiplies the score by a factor measured for that
/// specific RNA:DNA pair at that specific position. Positions are numbered 1–20 from the PAM-distal end.
/// </para>
/// <para>
/// The mismatch key is <c>r&lt;guide base, as RNA&gt;:d&lt;complement of the off-target base&gt;</c> — the "d" base is the
/// base on the strand the guide pairs with, which is why it is complemented. This library stores the 240
/// published values in a dense array indexed arithmetically, so scoring a site costs 20 array reads and no
/// allocation.
/// </para>
/// </remarks>
public static class CfdScorer
{
    /// <summary>Spacer length the published table covers.</summary>
    public const int TableLength = 20;

    /// <summary>
    /// The published penalty for a guide base paired against an off-target DNA base at <paramref name="position"/>
    /// (1–20, 1 being PAM-distal). Returns 1.0 for a Watson–Crick match.
    /// </summary>
    public static double MismatchWeight(char guideBase, char offTargetBase, int position)
    {
        if (position is < 1 or > TableLength) throw new ArgumentOutOfRangeException(nameof(position));
        int g = DnaUtility.Code((byte)guideBase);
        int o = DnaUtility.Code((byte)offTargetBase);
        if (g < 0 || o < 0) return 1.0;
        if (g == o) return 1.0;
        int d = DnaUtility.Code(DnaUtility.Complement((byte)char.ToUpperInvariant(offTargetBase)));
        return CfdTables.Mismatch[((g * 4 + d) * TableLength) + position - 1];
    }

    /// <summary>
    /// The published PAM penalty for the two nucleotides that follow the first N of an NGG PAM (so <c>GG</c>
    /// scores 1.0 and <c>AG</c> about 0.26). Unknown input scores 0.
    /// </summary>
    public static double PamWeight(ReadOnlySpan<char> lastTwoPamBases)
    {
        if (lastTwoPamBases.Length != 2) throw new ArgumentException("The CFD PAM table is keyed by two nucleotides.", nameof(lastTwoPamBases));
        int a = DnaUtility.Code((byte)lastTwoPamBases[0]);
        int b = DnaUtility.Code((byte)lastTwoPamBases[1]);
        if (a < 0 || b < 0) return 0;
        return CfdTables.Pam[a * 4 + b];
    }

    /// <summary>
    /// CFD score (0–1) of an off-target site: the 20 nt spacer against the 20 nt off-target protospacer, both
    /// 5′→3′ on the targeted strand, with the off-target's 3 nt PAM.
    /// </summary>
    /// <remarks>
    /// Spacers longer than 20 nt (SaCas9, Cas12a) are scored over their 20 PAM-proximal nucleotides, because the
    /// published table only covers SpCas9 geometry. That is an extrapolation and is documented as such — treat
    /// non-SpCas9 CFD values as a ranking heuristic rather than a calibrated cutting frequency.
    /// </remarks>
    public static double Score(ReadOnlySpan<char> spacer, ReadOnlySpan<char> offTargetProtospacer, ReadOnlySpan<char> pam)
    {
        if (spacer.Length != offTargetProtospacer.Length)
            throw new ArgumentException("The spacer and the off-target protospacer must be the same length.", nameof(offTargetProtospacer));
        if (spacer.Length < TableLength)
            throw new ArgumentException($"CFD needs at least {TableLength} nucleotides.", nameof(spacer));

        // Keep the 20 PAM-proximal nucleotides (the 3′ end for an SpCas9-style guide).
        int skip = spacer.Length - TableLength;
        spacer = spacer[skip..];
        offTargetProtospacer = offTargetProtospacer[skip..];

        double score = 1.0;
        for (int i = 0; i < TableLength; i++)
        {
            char g = char.ToUpperInvariant(spacer[i]);
            char o = char.ToUpperInvariant(offTargetProtospacer[i]);
            if (g == 'U') g = 'T';
            if (o == 'U') o = 'T';
            if (g == o) continue;
            score *= MismatchWeight(g, o, i + 1);
            if (score == 0) return 0;
        }

        if (pam.Length >= 2) score *= PamWeight(pam[^2..]);
        return score;
    }

    /// <summary>CFD score (0–1) for string inputs.</summary>
    public static double Score(string spacer, string offTargetProtospacer, string pam)
    {
        ArgumentNullException.ThrowIfNull(spacer);
        ArgumentNullException.ThrowIfNull(offTargetProtospacer);
        ArgumentNullException.ThrowIfNull(pam);
        return Score(spacer.AsSpan(), offTargetProtospacer.AsSpan(), pam.AsSpan());
    }

    /// <summary>
    /// Aggregate CFD specificity of a guide: <c>100 / (100 + Σ off-target CFD × 100)</c> over every off-target
    /// other than the intended site, expressed as 0–100. This is the "CFD specificity score" reported by CRISPOR.
    /// </summary>
    public static double GuideSpecificity(IEnumerable<OffTarget> offTargets)
    {
        ArgumentNullException.ThrowIfNull(offTargets);
        double sum = 0;
        foreach (var hit in offTargets)
        {
            if (hit.IsOnTarget) continue;
            sum += hit.CfdScore * 100.0;
        }
        return 100.0 * 100.0 / (100.0 + sum);
    }
}
