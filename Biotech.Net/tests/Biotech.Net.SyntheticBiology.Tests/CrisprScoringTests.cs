using Biotech.Net.SyntheticBiology.Crispr;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Doench 2014 Rule Set 1, Doench 2016 CFD and Hsu 2013 MIT scoring.
/// </summary>
/// <remarks>
/// The Rule Set 1 cases are hand-derived: for a homopolymeric 30-mer the set of weights that fire is small
/// enough to enumerate by hand from the published Supplementary Table 19, so the expected linear predictor is
/// written out term by term in the test. That checks the parts of the implementation that can actually be
/// wrong — the 0-based positions, the spacer window the GC term is taken over, the low/high GC branch and the
/// dinucleotide offsets — against something other than the table it reads.
/// </remarks>
public class CrisprScoringTests
{
    private static double Logistic(double x) => 1.0 / (1.0 + Math.Exp(-x));

    [Fact]
    public void RuleSet1_PolyA_MatchesHandSummedWeights()
    {
        // Spacer = context[4..24] contains no G or C, so the GC term is |10 - 0| * gcLow.
        const double intercept = 0.59763615, gcLow = -0.2026259;
        double gcTerm = 10 * gcLow;

        // Single-nucleotide 'A' terms at positions 2, 5, 11, 14, 15, 16, 18, 24.
        double singles = -0.3238875 + 0.03644004 - 0.466099 + 0.08537695
                         + 0.27262051 + 0.09745459 + 0.22508903 - 0.0773007;
        // Dinucleotide 'AA' terms at positions 5 and 23.
        double dinucleotides = -0.8348362 - 0.5794924;

        double expected = Logistic(intercept + gcTerm + singles + dinucleotides);

        double actual = OnTargetScorer.RuleSet1(new string('A', 30));

        Assert.Equal(expected, actual, 12);
        // Sanity: the hand value is ~0.0477, a poor guide.
        Assert.InRange(actual, 0.047, 0.048);
    }

    [Fact]
    public void RuleSet1_PolyG_UsesHighGcBranch()
    {
        const double intercept = 0.59763615, gcHigh = -0.1665878;
        double gcTerm = 10 * gcHigh; // 20 G in the spacer window: |10 - 20| * gcHigh

        // Single-nucleotide 'G' terms at positions 1, 4, 6, 16, 17, 19, 20, 23, 27, 29.
        double singles = -0.2753771 + 0.24595663 - 0.3932644 - 0.1755462 - 0.6780964
                         - 0.4173736 + 0.37989937 + 0.27891626 - 0.6890167 + 0.38634258;
        // Dinucleotide 'GG' terms at positions 6, 11, 12, 18, 20, 28.
        double dinucleotides = -0.4908167 - 1.5169074 - 0.3345637 - 0.6668087 - 0.7822076 - 0.69774;

        double expected = Logistic(intercept + gcTerm + singles + dinucleotides);

        Assert.Equal(expected, OnTargetScorer.RuleSet1(new string('G', 30)), 12);
    }

    [Fact]
    public void RuleSet1_IsCaseAndUraciInsensitive()
    {
        string context = "TTTTACGTACGTACGTACGTACGTTGGAAA";
        Assert.Equal(30, context.Length);
        Assert.Equal(OnTargetScorer.RuleSet1(context), OnTargetScorer.RuleSet1(context.ToLowerInvariant()), 12);
        Assert.Equal(OnTargetScorer.RuleSet1(context), OnTargetScorer.RuleSet1(context.Replace('T', 'U')), 12);
    }

    [Fact]
    public void RuleSet1_RejectsWrongLength() =>
        Assert.Throws<ArgumentException>(() => OnTargetScorer.RuleSet1("ACGT"));

    [Fact]
    public void RuleSet1_StaysInUnitInterval()
    {
        var random = new Random(7);
        for (int i = 0; i < 200; i++)
        {
            string context = string.Concat(Enumerable.Range(0, 30).Select(_ => "ACGT"[random.Next(4)]));
            double score = OnTargetScorer.RuleSet1(context);
            Assert.InRange(score, 0.0, 1.0);
        }
    }

    // ---- CFD (Doench et al. 2016, Supplementary Table 19) ----

    [Fact]
    public void Cfd_PerfectMatchWithCanonicalPamScoresOne()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        Assert.Equal(1.0, CfdScorer.Score(spacer, spacer, "TGG"), 12);
        Assert.Equal(1.0, CfdScorer.Score(spacer, spacer, "AGG"), 12);
    }

    [Theory]
    // The published PAM table, keyed by the two nucleotides following the N of NGG.
    [InlineData("GG", 1.0)]
    [InlineData("AG", 0.259259259)]
    [InlineData("CG", 0.107142857)]
    [InlineData("GA", 0.069444444)]
    [InlineData("GC", 0.022222222)]
    [InlineData("GT", 0.016129032)]
    [InlineData("TG", 0.038961039)]
    [InlineData("AA", 0.0)]
    [InlineData("TT", 0.0)]
    public void Cfd_PamWeightsMatchPublishedTable(string lastTwo, double expected) =>
        Assert.Equal(expected, CfdScorer.PamWeight(lastTwo), 9);

    [Fact]
    public void Cfd_NonCanonicalPamScalesThePerfectMatch()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        Assert.Equal(0.259259259, CfdScorer.Score(spacer, spacer, "AAG"), 9);
        Assert.Equal(0.0, CfdScorer.Score(spacer, spacer, "AAA"), 12);
    }

    [Theory]
    // The mismatch key is r<guide base>:d<complement of the off-target base>, so an off-target G opposite a
    // guide C is the rC:dC entry. Values are from the published table, position 1 is PAM-distal.
    [InlineData('C', 'G', 1, 0.913043478)]  // rC:dC, position 1
    [InlineData('G', 'A', 1, 0.9)]          // rG:dT, position 1
    [InlineData('A', 'C', 1, 0.857142857)]  // rA:dG, position 1
    [InlineData('A', 'G', 16, 0.192307692)] // rA:dC, position 16
    [InlineData('C', 'T', 15, 0.066666667)] // rC:dA, position 15
    public void Cfd_MismatchWeightsMatchPublishedTable(char guide, char offTarget, int position, double expected) =>
        Assert.Equal(expected, CfdScorer.MismatchWeight(guide, offTarget, position), 9);

    [Fact]
    public void Cfd_WatsonCrickPairsAreFree()
    {
        foreach (char b in "ACGT") Assert.Equal(1.0, CfdScorer.MismatchWeight(b, b, 7), 12);
    }

    [Fact]
    public void Cfd_IsTheProductOfItsPerPositionWeights()
    {
        // Two mismatches: guide C against off-target G at position 1, guide A against off-target C at 16.
        const string spacer = "CAGTCCGAGCAGAAGAAGAA";
        const string offTarget = "GAGTCCGAGCAGAAGCAGAA";
        double expected = CfdScorer.MismatchWeight('C', 'G', 1) * CfdScorer.MismatchWeight('A', 'C', 16);
        Assert.Equal(expected, CfdScorer.Score(spacer, offTarget, "TGG"), 10);
    }

    [Fact]
    public void Cfd_SeedMismatchesHurtMoreThanDistalOnes()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        // Position 2 (PAM-distal) versus position 18 (seed), same substitution.
        string distal = string.Concat("G", "C", spacer.AsSpan(2));
        string seed = string.Concat(spacer.AsSpan(0, 17), "C", spacer.AsSpan(18));
        Assert.True(CfdScorer.Score(spacer, distal, "TGG") > CfdScorer.Score(spacer, seed, "TGG"));
    }

    [Fact]
    public void Cfd_GuideSpecificityFallsAsOffTargetsAccumulate()
    {
        var one = new[] { OffTargetWithCfd(0.5) };
        var many = new[] { OffTargetWithCfd(0.5), OffTargetWithCfd(0.5), OffTargetWithCfd(0.5) };

        // 100 * 100 / (100 + 50) for a single half-scoring off-target.
        Assert.Equal(100.0 * 100.0 / 150.0, CfdScorer.GuideSpecificity(one), 9);
        Assert.True(CfdScorer.GuideSpecificity(many) < CfdScorer.GuideSpecificity(one));
        Assert.Equal(100.0, CfdScorer.GuideSpecificity([]), 9);
    }

    private static OffTarget OffTargetWithCfd(double cfd) =>
        new("chr", 0, Core.Strand.Forward, new string('A', 20), "TGG", 2, BulgeKind.None, cfd, 0, IsOnTarget: false);

    // ---- MIT / Hsu 2013 ----

    [Fact]
    public void Mit_PerfectMatchScoresOneHundred()
    {
        Assert.Equal(100.0, MitSpecificityScorer.HitScore(ReadOnlySpan<int>.Empty), 12);
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        Assert.Equal(100.0, MitSpecificityScorer.HitScore(spacer, spacer), 12);
    }

    [Fact]
    public void Mit_SingleMismatchIsJustThePositionWeight()
    {
        // One mismatch: no distance or count correction, so the score is 100 * (1 - w[p]).
        var weights = MitSpecificityScorer.PositionWeights;
        for (int p = 0; p < 20; p++)
        {
            double expected = 100.0 * (1.0 - weights[p]);
            Assert.Equal(expected, MitSpecificityScorer.HitScore([p]), 10);
        }
    }

    [Fact]
    public void Mit_PositionWeightsArePublishedValues()
    {
        var weights = MitSpecificityScorer.PositionWeights;
        Assert.Equal(20, weights.Length);
        Assert.Equal(0.000, weights[0], 6);   // PAM-distal, free
        Assert.Equal(0.014, weights[2], 6);
        Assert.Equal(0.851, weights[13], 6);  // the largest penalty
        Assert.Equal(0.583, weights[19], 6);  // PAM-proximal
    }

    [Fact]
    public void Mit_TwoMismatchesApplyTheDistanceAndCountCorrections()
    {
        var weights = MitSpecificityScorer.PositionWeights;
        int[] positions = [5, 15];
        double meanDistance = 15 - 5; // (last - first) / (n - 1)
        double expected = 100.0
                          * (1 - weights[5]) * (1 - weights[15])
                          * (1.0 / (((19.0 - meanDistance) / 19.0) * 4.0 + 1.0))
                          * (1.0 / 4.0);
        Assert.Equal(expected, MitSpecificityScorer.HitScore(positions), 10);
    }

    [Fact]
    public void Mit_SeedMismatchesHurtMoreThanDistalOnes()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        string distal = string.Concat("C", spacer.AsSpan(1));
        string seed = string.Concat(spacer.AsSpan(0, 19), "C");
        Assert.True(MitSpecificityScorer.HitScore(spacer, distal) > MitSpecificityScorer.HitScore(spacer, seed));
    }

    [Fact]
    public void Mit_GuideScoreAggregatesHits()
    {
        var hits = new[]
        {
            new OffTarget("chr", 0, Core.Strand.Forward, new string('A', 20), "TGG", 0, BulgeKind.None, 1, 100, IsOnTarget: true),
            new OffTarget("chr", 50, Core.Strand.Forward, new string('A', 20), "TGG", 2, BulgeKind.None, 0.2, 25, IsOnTarget: false),
        };
        // The on-target hit is excluded, so the denominator is 100 + 25.
        Assert.Equal(100.0 * 100.0 / 125.0, MitSpecificityScorer.GuideScore(hits), 9);
    }

    // ---- composition penalties ----

    [Fact]
    public void Penalties_FlagPolyTAndExtremeGc()
    {
        var polyT = OnTargetScorer.Penalties("AATTTTAAAAAAAAAAAAGG");
        Assert.Equal(1.0, polyT.PolyTPenalty, 12);

        var balanced = OnTargetScorer.Penalties("ACGTACGTACGTACGTACGT");
        Assert.Equal(0.0, balanced.GcPenalty, 12);
        Assert.Equal(0.0, balanced.PolyTPenalty, 12);

        var gcRich = OnTargetScorer.Penalties("GCGCGCGCGCGCGCGCGCGC"); // 100% GC
        Assert.Equal(Math.Min(0.5, (1.0 - 0.70) * 2.0), gcRich.GcPenalty, 12);
        Assert.InRange(gcRich.Total, 0, 1);
        Assert.Equal(1 - gcRich.Total, gcRich.Multiplier, 12);
    }

    [Fact]
    public void Penalties_HomopolymerRunsAreScaled()
    {
        // A run of five is one over the threshold of four: 0.15 per extra base.
        Assert.Equal(0.15, OnTargetScorer.Penalties("ACGTACGTAAAAACGTACGT").HomopolymerPenalty, 12);
        Assert.Equal(0.0, OnTargetScorer.Penalties("ACGTACGTACGTACGTACGT").HomopolymerPenalty, 12);
    }

    [Fact]
    public void SelfComplementaryStem_IsFoundInAHairpinSpacer()
    {
        // A 6 nt stem with a 4 nt loop.
        Assert.True(OnTargetScorer.LongestSelfComplementaryStem("GGGCCCAAAAGGGCCCACGT") >= 5);
        Assert.True(OnTargetScorer.LongestSelfComplementaryStem("AAAAAAAAAAAAAAAAAAAA") < 5);
    }
}
