namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// A breakdown of the heuristic quality checks applied to a spacer, alongside the published on-target model.
/// </summary>
/// <param name="GcPenalty">Penalty (0–1, higher is worse) for GC content outside the 40–70% sweet spot.</param>
/// <param name="PolyTPenalty">1 when the spacer contains a <c>TTTT</c> Pol III terminator, else 0.</param>
/// <param name="HomopolymerPenalty">Penalty (0–1) for single-base runs of 5 or more.</param>
/// <param name="SelfComplementarityPenalty">Penalty (0–1) for a spacer that can fold back on itself.</param>
public readonly record struct GuidePenalties(
    double GcPenalty,
    double PolyTPenalty,
    double HomopolymerPenalty,
    double SelfComplementarityPenalty)
{
    /// <summary>Total penalty, clamped to 0–1.</summary>
    public double Total => Math.Clamp(GcPenalty + PolyTPenalty + HomopolymerPenalty + SelfComplementarityPenalty, 0, 1);

    /// <summary>A 0–1 multiplier (<c>1 − Total</c>) to apply to a model score.</summary>
    public double Multiplier => 1 - Total;
}

/// <summary>
/// On-target (cutting efficiency) scoring for SpCas9 guides.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="RuleSet1(string)"/> is the logistic model of Doench et al., <i>Nature Biotechnology</i> 32:1262 (2014),
/// using the published Supplementary Table 19 weights: an intercept, a GC-count term and 70 positional
/// single-nucleotide and dinucleotide terms over a 30-mer of 4 nt upstream context, the 20 nt spacer, the 3 nt
/// PAM and 3 nt of downstream context. The result is the logistic of the linear predictor, in 0–1.
/// </para>
/// <para>
/// The model was trained on SpCas9 with an NGG PAM and 20 nt spacers. Applying it to another nuclease is not
/// meaningful, so <see cref="Score(Guide)"/> returns <see cref="double.NaN"/> in that case rather than a number
/// that looks authoritative. Rule Set 2 (Azimuth) is a boosted-tree model whose parameters are not published
/// as a table and so is not reproduced here.
/// </para>
/// </remarks>
public static class OnTargetScorer
{
    /// <summary>The 30-mer context length the model needs.</summary>
    public const int ContextLength = 30;

    /// <summary>Offset of the spacer inside the 30-mer context.</summary>
    public const int SpacerOffsetInContext = 4;

    /// <summary>
    /// Doench 2014 Rule Set 1 score (0–1) for a 30-mer context: 4 nt upstream, 20 nt spacer, 3 nt PAM,
    /// 3 nt downstream, all on the targeted strand.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="context30"/> is not 30 characters long.</exception>
    public static double RuleSet1(ReadOnlySpan<char> context30)
    {
        if (context30.Length != ContextLength)
            throw new ArgumentException($"Rule Set 1 needs a {ContextLength}-mer context, got {context30.Length}.", nameof(context30));

        Span<char> seq = stackalloc char[ContextLength];
        for (int i = 0; i < ContextLength; i++)
        {
            char c = char.ToUpperInvariant(context30[i]);
            seq[i] = c == 'U' ? 'T' : c;
        }

        int gc = 0;
        for (int i = SpacerOffsetInContext; i < SpacerOffsetInContext + 20; i++)
            if (seq[i] == 'G' || seq[i] == 'C') gc++;

        double score = Doench2014Weights.Intercept
                       + Math.Abs(10 - gc) * (gc <= 10 ? Doench2014Weights.GcLow : Doench2014Weights.GcHigh);

        var terms = Doench2014Weights.Terms;
        for (int i = 0; i < terms.Length; i++)
        {
            ref readonly var term = ref terms[i];
            string motif = term.Motif;
            int pos = term.Position;
            if (pos + motif.Length > ContextLength) continue;
            bool hit = true;
            for (int j = 0; j < motif.Length; j++)
            {
                if (seq[pos + j] != motif[j]) { hit = false; break; }
            }
            if (hit) score += term.Weight;
        }

        return 1.0 / (1.0 + Math.Exp(-score));
    }

    /// <summary>Doench 2014 Rule Set 1 score (0–1) for a 30-mer context string.</summary>
    public static double RuleSet1(string context30)
    {
        ArgumentNullException.ThrowIfNull(context30);
        return RuleSet1(context30.AsSpan());
    }

    /// <summary>
    /// Rule Set 1 for a guide, or <see cref="double.NaN"/> when the guide is not a 20 nt SpCas9-family NGG guide
    /// or its 30-mer context ran off the end of the target sequence.
    /// </summary>
    public static double Score(Guide guide)
    {
        ArgumentNullException.ThrowIfNull(guide);
        if (guide.Context is not { Length: ContextLength }) return double.NaN;
        if (guide.Nuclease.GuideLength != 20 || guide.Nuclease.PamSide != PamSide.ThreePrime) return double.NaN;
        if (guide.Nuclease.PamLength != 3) return double.NaN;
        return RuleSet1(guide.Context);
    }

    /// <summary>
    /// The sequence-composition penalties that every guide-design tool applies on top of a model score.
    /// </summary>
    /// <remarks>
    /// These follow the rules of thumb published with the major designers (Benchling, CHOPCHOP, CRISPOR):
    /// 40–70% GC is optimal, a <c>TTTT</c> run terminates U6 transcription, long homopolymers synthesise and
    /// sequence poorly, and a spacer with a long self-complementary stem folds instead of loading into Cas9.
    /// </remarks>
    public static GuidePenalties Penalties(string spacer)
    {
        ArgumentException.ThrowIfNullOrEmpty(spacer);
        var bytes = DnaUtility.Ascii(spacer);

        double gc = DnaUtility.GcFraction(bytes);
        double gcPenalty = gc switch
        {
            >= 0.40 and <= 0.70 => 0,
            < 0.40 => Math.Min(0.5, (0.40 - gc) * 2.0),
            _ => Math.Min(0.5, (gc - 0.70) * 2.0),
        };

        double polyT = spacer.Contains("TTTT", StringComparison.OrdinalIgnoreCase) ? 1.0 : 0.0;

        int run = DnaUtility.LongestHomopolymer(bytes);
        double homopolymer = run < 5 ? 0 : Math.Min(0.5, (run - 4) * 0.15);

        double self = Math.Min(0.3, LongestSelfComplementaryStem(spacer) switch
        {
            < 5 => 0,
            int n => (n - 4) * 0.1,
        });

        return new GuidePenalties(gcPenalty, polyT, homopolymer, self);
    }

    /// <summary>
    /// A combined 0–1 efficiency estimate: the Rule Set 1 score reduced by the composition penalties. When the
    /// model does not apply, the penalties alone give a score relative to a neutral 0.5 baseline.
    /// </summary>
    public static double CombinedScore(Guide guide)
    {
        ArgumentNullException.ThrowIfNull(guide);
        double model = Score(guide);
        double baseline = double.IsNaN(model) ? 0.5 : model;
        return Math.Clamp(baseline * Penalties(guide.Spacer).Multiplier, 0, 1);
    }

    /// <summary>
    /// Length of the longest stem the spacer can form with itself (a sub-sequence whose reverse complement also
    /// occurs in the spacer, separated by at least 3 nt of loop).
    /// </summary>
    internal static int LongestSelfComplementaryStem(string spacer)
    {
        var fwd = DnaUtility.Ascii(DnaUtility.NormalizeDna(spacer));
        int n = fwd.Length;
        int best = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = n - 1; j > i + 3; j--)
            {
                int len = 0;
                while (i + len < j - len && DnaUtility.Complement(fwd[i + len]) == fwd[j - len]
                       && (j - len) - (i + len) > 3)
                {
                    len++;
                }
                if (len > best) best = len;
            }
        }
        return best;
    }
}
