namespace Biotech.Net.SyntheticBiology.Circuits;

/// <summary>
/// The decomposed output of the approximate ribosome binding site strength calculation.
/// All energies are in kcal/mol; negative values are stabilising.
/// </summary>
/// <param name="TranslationInitiationRate">
/// Predicted translation initiation rate in arbitrary units, proportional to
/// exp(-beta * <paramref name="TotalFreeEnergy"/>). Only ratios between predictions are meaningful.
/// </param>
/// <param name="TotalFreeEnergy">
/// The summed free energy change of translation initiation, dG_total (kcal/mol). More negative means
/// a stronger ribosome binding site.
/// </param>
/// <param name="MrnaRnaseDuplexEnergy">
/// The mRNA:16S rRNA duplex term, dG_mRNA:rRNA (kcal/mol): the hybridisation energy of the best
/// Shine-Dalgarno-like window in the 5' UTR against the 3' anti-Shine-Dalgarno tail of the 16S rRNA.
/// </param>
/// <param name="SpacingPenalty">
/// The spacing term, dG_spacing (kcal/mol, non-negative): a penalty for an aligned spacing that
/// departs from the optimum of five nucleotides.
/// </param>
/// <param name="SpacingNucleotides">
/// Aligned spacing: the number of nucleotides between the 3' end of the Shine-Dalgarno duplex and the
/// first base of the start codon.
/// </param>
/// <param name="ShineDalgarnoMatch">
/// The 5' UTR subsequence (DNA alphabet, upper case) found to pair with the anti-Shine-Dalgarno tail,
/// or an empty string when no window of at least three paired bases was found.
/// </param>
public sealed record RbsPrediction(
    double TranslationInitiationRate,
    double TotalFreeEnergy,
    double MrnaRnaseDuplexEnergy,
    double SpacingPenalty,
    int SpacingNucleotides,
    string ShineDalgarnoMatch)
{
    /// <summary>
    /// The mRNA self-folding term, dG_mRNA (kcal/mol, non-positive): the stability of the most stable
    /// hairpin detected in the 5' UTR plus the first 35 nucleotides of the coding sequence. Structure
    /// that must be unfolded to admit the ribosome opposes initiation, so this term enters dG_total
    /// with a minus sign.
    /// </summary>
    public double MrnaFoldingEnergy { get; init; }

    /// <summary>
    /// The start-codon term, dG_start (kcal/mol): the initiator tRNA(fMet) anticodon pairs best with
    /// AUG, less well with GUG and UUG, and poorly with anything else.
    /// </summary>
    public double StartCodonEnergy { get; init; }

    /// <summary>The start codon read from the first three bases of the coding sequence, upper case.</summary>
    public string StartCodon { get; init; } = "";
}

/// <summary>
/// An approximate, deterministic heuristic for ribosome binding site strength in <i>E. coli</i>,
/// inspired by the thermodynamic RBS Calculator of Salis, Mirsky and Voigt
/// (<i>Automated design of synthetic ribosome binding sites to control protein expression</i>,
/// Nature Biotechnology 27:946-950, 2009).
/// <para>
/// <b>Read this before using the numbers.</b> This class reproduces the <i>functional form</i> of the
/// published model, namely
/// dG_total = dG_mRNA:rRNA - dG_mRNA + dG_spacing + dG_standby + dG_start and
/// rate proportional to exp(-beta * dG_total) with beta ~ 0.45 mol/kcal, but it does <b>not</b> use the
/// published parameter set. The duplex term here is a crude per-base-pair sum rather than a nearest
/// neighbour RNA thermodynamic calculation, the folding term is a stem-detection scan rather than a
/// partition function over secondary structures, and the spacing and start-codon coefficients are
/// hand-chosen to have the published <i>shape</i> only. Consequently the predictions are useful for
/// <b>relative ranking</b> of candidate ribosome binding sites and for teaching the structure of the
/// model; they are <b>not</b> calibrated expression predictions. Quantitative forward engineering
/// should use the real RBS Calculator (salislab.net), UTR Designer, or RBS Designer, all of which are
/// validated against measured fluorescence data.
/// </para>
/// </summary>
public static class RbsStrengthModel
{
    /// <summary>
    /// The 3' terminal anti-Shine-Dalgarno tail of the <i>E. coli</i> 16S ribosomal RNA, written
    /// 5' to 3' as <c>ACCUCCUUA</c> (Shine &amp; Dalgarno, PNAS 71:1342, 1974). An mRNA Shine-Dalgarno
    /// element pairs with it antiparallel, so the ideal mRNA element is its reverse complement,
    /// <c>UAAGGAGGU</c>.
    /// </summary>
    public const string AntiShineDalgarno = "ACCUCCUUA";

    /// <summary>
    /// Boltzmann-like proportionality constant beta in rate ~ exp(-beta * dG_total), in mol/kcal.
    /// Salis et al. fitted beta = 0.45 mol/kcal across a library of 132 measured ribosome binding sites.
    /// </summary>
    public const double Beta = 0.45;

    /// <summary>Arbitrary proportionality constant K in rate = K * exp(-beta * dG_total).</summary>
    public const double RateScale = 2500.0;

    /// <summary>
    /// Optimal aligned spacing between the 3' end of the Shine-Dalgarno duplex and the start codon,
    /// in nucleotides. Chen et al. (Nucleic Acids Research 22:4953, 1994) measured an optimum of
    /// about five nucleotides, with expression falling off sharply when the spacing is shorter and
    /// more gently when it is longer.
    /// </summary>
    public const int OptimalSpacing = 5;

    /// <summary>
    /// A reference 5' UTR carrying a perfect Shine-Dalgarno element at the optimal spacing; the
    /// normalisation point for <see cref="RelativeStrength(string, string)"/>.
    /// </summary>
    public const string ConsensusUtr = "TAAGGAGGTAAAAA";

    /// <summary>
    /// A reference coding sequence with an AUG start and no strong internal structure, used together
    /// with <see cref="ConsensusUtr"/> as the normalisation point for relative strengths.
    /// </summary>
    public const string ReferenceCodingSequence = "ATGAAAGAAGCAGAAGCAAAAGAAGCATAA";

    private const double EnergyGc = -3.0;
    private const double EnergyAu = -2.0;
    private const double EnergyWobble = -1.0;
    private const double StackingBonus = -0.5;
    private const int MinimumDuplexLength = 3;
    private const int FoldingWindowIntoCds = 35;

    private static readonly double ReferenceRate =
        Predict(ConsensusUtr, ReferenceCodingSequence).TranslationInitiationRate;

    /// <summary>
    /// Predicts the translation initiation strength of a ribosome binding site from its 5' untranslated
    /// region and the coding sequence that follows it.
    /// </summary>
    /// <param name="fivePrimeUtr">
    /// The untranslated region immediately upstream of the start codon (DNA alphabet; T is read as U).
    /// It should contain the Shine-Dalgarno element and the spacer, and must not include the start codon.
    /// </param>
    /// <param name="codingSequence">
    /// The coding sequence beginning with the start codon. Only the first
    /// <see cref="FoldingWindowIntoCds"/> nucleotides influence the folding term.
    /// </param>
    /// <returns>The decomposed prediction; deterministic for a given input.</returns>
    public static RbsPrediction Predict(string fivePrimeUtr, string codingSequence)
    {
        ArgumentNullException.ThrowIfNull(fivePrimeUtr);
        ArgumentNullException.ThrowIfNull(codingSequence);

        string utr = ToRna(fivePrimeUtr);
        string cds = ToRna(codingSequence);

        FindBestDuplex(utr, out double duplexEnergy, out int duplexEnd, out string match);

        int spacing = match.Length == 0 ? utr.Length : utr.Length - duplexEnd;
        double spacingPenalty = SpacingEnergy(spacing);

        double folding = FoldingEnergy(utr, cds);
        string startCodon = cds.Length >= 3 ? cds[..3] : cds;
        double start = StartCodonEnergy(startCodon);

        double total = duplexEnergy - folding + spacingPenalty + start;
        double rate = RateScale * Math.Exp(-Beta * total);

        return new RbsPrediction(rate, total, duplexEnergy, spacingPenalty, spacing, FromRna(match))
        {
            MrnaFoldingEnergy = folding,
            StartCodonEnergy = start,
            StartCodon = FromRna(startCodon),
        };
    }

    /// <summary>
    /// Translation initiation strength normalised so that a perfect Shine-Dalgarno element at the
    /// optimal spacing in front of an unstructured AUG-initiated coding sequence
    /// (<see cref="ConsensusUtr"/> plus <see cref="ReferenceCodingSequence"/>) scores exactly 1.
    /// Values above 1 are therefore rare; typical community ribosome binding sites score well below it.
    /// </summary>
    /// <param name="utr">5' untranslated region upstream of the start codon.</param>
    /// <param name="cds">Coding sequence beginning with the start codon.</param>
    public static double RelativeStrength(string utr, string cds) =>
        Predict(utr, cds).TranslationInitiationRate / ReferenceRate;

    /// <summary>
    /// The spacing term dG_spacing in kcal/mol: zero at the optimum of
    /// <see cref="OptimalSpacing"/> nucleotides, rising quadratically when the Shine-Dalgarno element
    /// sits too close to the start codon (steric clash with the 30S subunit) and linearly when it sits
    /// too far away. This has the shape of the Salis 2009 spacing penalty but not its coefficients.
    /// </summary>
    /// <param name="spacing">Aligned spacing in nucleotides.</param>
    public static double SpacingEnergy(int spacing)
    {
        int delta = spacing - OptimalSpacing;
        if (delta == 0) return 0.0;
        return delta < 0 ? 0.5 * delta * delta : 0.3 * delta;
    }

    /// <summary>
    /// The start-codon term dG_start in kcal/mol. AUG is the canonical initiator read by
    /// tRNA-fMet; GUG and UUG initiate in roughly 10 percent of <i>E. coli</i> genes and are weaker;
    /// anything else is strongly disfavoured.
    /// </summary>
    /// <param name="startCodon">The first three bases of the coding sequence (DNA or RNA alphabet).</param>
    public static double StartCodonEnergy(string startCodon)
    {
        string codon = ToRna(startCodon);
        return codon switch
        {
            "AUG" => -1.19,
            "GUG" => -0.50,
            "UUG" => 0.00,
            "AUU" => 0.40,
            "CUG" => 0.50,
            _ => 2.00,
        };
    }

    // Scans every contiguous ungapped alignment of the UTR against the anti-Shine-Dalgarno tail and
    // keeps the window that minimises duplex energy plus the spacing penalty its position implies.
    private static void FindBestDuplex(string utr, out double energy, out int endExclusive, out string match)
    {
        energy = 0.0;
        endExclusive = utr.Length;
        match = "";

        // The anti-SD read 3' to 5' so that an antiparallel duplex becomes a left-to-right comparison.
        Span<char> reversed = stackalloc char[AntiShineDalgarno.Length];
        for (int i = 0; i < AntiShineDalgarno.Length; i++)
            reversed[i] = AntiShineDalgarno[AntiShineDalgarno.Length - 1 - i];

        double bestScore = double.PositiveInfinity;
        for (int i = 0; i < utr.Length; i++)
        {
            for (int j = 0; j < reversed.Length; j++)
            {
                double running = 0.0;
                int maxLength = Math.Min(utr.Length - i, reversed.Length - j);
                for (int k = 0; k < maxLength; k++)
                {
                    double pair = PairEnergy(utr[i + k], reversed[j + k]);
                    if (pair == 0.0) break;
                    running += pair;
                    if (k > 0) running += StackingBonus;
                    int length = k + 1;
                    if (length < MinimumDuplexLength) continue;

                    int end = i + length;
                    double score = running + SpacingEnergy(utr.Length - end);
                    if (score < bestScore - 1e-12)
                    {
                        bestScore = score;
                        energy = running;
                        endExclusive = end;
                        match = utr.Substring(i, length);
                    }
                }
            }
        }
    }

    // Crude folding term: the most stable perfectly paired stem of 3 to 8 bp with a loop of at least
    // 3 nt, found anywhere in the 5' UTR plus the first 35 nt of the coding sequence.
    private static double FoldingEnergy(string utr, string cds)
    {
        int take = Math.Min(FoldingWindowIntoCds, cds.Length);
        string window = utr + cds[..take];
        double best = 0.0;
        for (int i = 0; i < window.Length; i++)
        {
            for (int j = window.Length - 1; j > i; j--)
            {
                double running = 0.0;
                for (int k = 0; k < 8; k++)
                {
                    int left = i + k;
                    int right = j - k;
                    if (right - left < 4) break;
                    double pair = PairEnergy(window[left], window[right]);
                    if (pair == 0.0) break;
                    running += pair;
                    if (k > 0) running += StackingBonus;
                    if (k + 1 >= 3 && running < best) best = running;
                }
            }
        }
        return best;
    }

    // Per-base-pair stabilisation in kcal/mol: Watson-Crick G:C and A:U plus the G:U wobble.
    private static double PairEnergy(char a, char b) => (a, b) switch
    {
        ('G', 'C') or ('C', 'G') => EnergyGc,
        ('A', 'U') or ('U', 'A') => EnergyAu,
        ('G', 'U') or ('U', 'G') => EnergyWobble,
        _ => 0.0,
    };

    private static string ToRna(string sequence)
    {
        var buffer = new char[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            buffer[i] = c == 'T' ? 'U' : c;
        }
        return new string(buffer);
    }

    private static string FromRna(string sequence)
    {
        var buffer = new char[sequence.Length];
        for (int i = 0; i < sequence.Length; i++) buffer[i] = sequence[i] == 'U' ? 'T' : sequence[i];
        return new string(buffer);
    }
}
