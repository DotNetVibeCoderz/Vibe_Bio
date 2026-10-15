using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology.Codons;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Codon usage tables, the codon adaptation index (Sharp and Li 1987) and the codon optimiser.
/// </summary>
public class CodonTests
{
    /// <summary>
    /// A deliberately tiny table: two amino-acid families with hand-chosen usage, so every relative
    /// adaptiveness and every CAI below can be worked out on paper.
    /// </summary>
    private static CodonUsageTable ToyTable() => CodonUsageTable.FromCounts(
        "toy", "toy organism",
        new Dictionary<string, double>
        {
            // Leucine: CTG is the family optimum, TTA is half as frequent, the rest are unused.
            ["CTG"] = 100,
            ["TTA"] = 50,
            // Phenylalanine: TTC optimal, TTT a quarter as frequent.
            ["TTC"] = 80,
            ["TTT"] = 20,
        });

    [Fact]
    public void RelativeAdaptiveness_IsFrequencyOverTheFamilyOptimum()
    {
        var table = ToyTable();

        Assert.Equal(1.0, table.RelativeAdaptiveness("CTG"), 12);
        Assert.Equal(0.5, table.RelativeAdaptiveness("TTA"), 12);
        Assert.Equal(0.0, table.RelativeAdaptiveness("CTC"), 12); // never observed
        Assert.Equal(1.0, table.RelativeAdaptiveness("TTC"), 12);
        Assert.Equal(0.25, table.RelativeAdaptiveness("TTT"), 12);

        Assert.Equal(0.1, table.Frequency("CTG"), 12); // 100 per thousand
        Assert.Equal(100.0, table.PerThousand("CTG"), 12);
    }

    [Fact]
    public void Cai_IsTheGeometricMeanOfTheWeights()
    {
        var table = ToyTable();

        // w = 1, 0.5, 0.25 -> (1 * 0.5 * 0.25)^(1/3) = 0.125^(1/3) = 0.5 exactly.
        Assert.Equal(0.5, table.Cai("CTGTTATTT"), 12);
        Assert.Equal(1.0, table.Cai("CTGCTGTTC"), 12);
        Assert.Equal(Math.Pow(0.5 * 0.25, 0.5), table.Cai("TTATTT"), 12);
    }

    [Fact]
    public void Cai_SkipsStopsMethionineAndTryptophan()
    {
        var table = ToyTable();
        // ATG (Met), TGG (Trp) and TAA (stop) carry no synonymous choice, so they must not move the index.
        Assert.Equal(table.Cai("CTGTTATTT"), table.Cai("ATGCTGTTATGGTTTTAA"), 12);
        // Nothing scoreable at all.
        Assert.Equal(0.0, table.Cai("ATGTGGTAA"), 12);
    }

    [Fact]
    public void Cai_UsesTheZeroWeightFloor()
    {
        var table = ToyTable();
        // CTC was never observed, so it takes the 0.001 floor: (1 * 0.001)^(1/2).
        Assert.Equal(Math.Sqrt(CodonUsageTable.ZeroWeightFloor), table.Cai("CTGCTC"), 12);
    }

    [Fact]
    public void Cai_RejectsAPartialCodon() =>
        Assert.Throws<ArgumentException>(() => ToyTable().Cai("CTGTT"));

    [Fact]
    public void SynonymousCodons_AreOrderedByDecreasingUsage()
    {
        var table = ToyTable();
        var leucine = table.SynonymousCodons('L');

        Assert.Equal("CTG", leucine[0]);
        Assert.Equal("TTA", leucine[1]);
        Assert.Equal(6, leucine.Count); // the standard code has six leucine codons
        Assert.Equal("CTG", table.MostFrequentCodon('L'));
        Assert.Equal("TTT", table.RarestCodon('F'));
        Assert.Throws<ArgumentException>(() => table.SynonymousCodons('X'));
    }

    [Fact]
    public void Gc3Content_CountsThirdPositionsOnly()
    {
        var table = ToyTable();
        // CTG TTA TTT -> third bases G, A, T.
        Assert.Equal(1.0 / 3.0, table.Gc3Content("CTGTTATTT"), 12);
        Assert.Equal(1.0, table.Gc3Content("CTGCTCTTC"), 12);
        Assert.Equal(0.0, table.Gc3Content("CTATTATTT"), 12);
    }

    [Fact]
    public void FromCounts_ValidatesItsInput()
    {
        Assert.Throws<ArgumentException>(() => CodonUsageTable.FromCounts("x", "x", new Dictionary<string, double> { ["CT"] = 1 }));
        Assert.Throws<ArgumentException>(() => CodonUsageTable.FromCounts("x", "x", new Dictionary<string, double> { ["CTN"] = 1 }));
        Assert.Throws<ArgumentException>(() => CodonUsageTable.FromCounts("x", "x", new Dictionary<string, double> { ["CTG"] = -1 }));
    }

    [Fact]
    public void FromCodingSequences_CountsWhatItIsGiven()
    {
        // Six codons: CTG three times, TTA once, TTC twice.
        var table = CodonUsageTable.FromCodingSequences("counted", ["CTGCTGCTG", "TTATTCTTC"]);

        Assert.Equal(1.0, table.RelativeAdaptiveness("CTG"), 12);
        Assert.Equal(1.0 / 3.0, table.RelativeAdaptiveness("TTA"), 12);
        Assert.Equal(1000.0 * 3 / 6, table.PerThousand("CTG"), 9);
        Assert.Equal(1000.0 * 2 / 6, table.PerThousand("TTC"), 9);
    }

    [Fact]
    public void BuiltInTables_AreWellFormed()
    {
        foreach (var table in CodonUsageTable.All)
        {
            Assert.Same(table, CodonUsageTable.FromName(table.Name));
            Assert.InRange(table.GcContent, 0.2, 0.8);
            // Every family has an optimum of weight exactly 1.
            foreach (char aa in "ACDEFGHIKLMNPQRSTVWY")
                Assert.Equal(1.0, table.RelativeAdaptiveness(table.MostFrequentCodon(aa)), 12);
        }
        // E. coli and human disagree about synonymous choice, and the distance is symmetric and positive.
        double d = CodonUsageTable.EColiK12.CodonAdaptationDistance(CodonUsageTable.HomoSapiens);
        Assert.True(d > 0);
        Assert.Equal(d, CodonUsageTable.HomoSapiens.CodonAdaptationDistance(CodonUsageTable.EColiK12), 12);
        Assert.Equal(0.0, CodonUsageTable.EColiK12.CodonAdaptationDistance(CodonUsageTable.EColiK12), 12);
    }

    // ---- optimiser ----

    private const string Protein = "MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSR";

    [Fact]
    public void Optimize_MostFrequentGivesCaiOfOneAndPreservesTheProtein()
    {
        var options = new CodonOptimizationOptions
        {
            Target = CodonUsageTable.EColiK12,
            Method = CodonOptimizationMethod.MostFrequent,
        };

        var result = CodonOptimizer.Optimize(Protein, options);

        Assert.Equal(Protein.Length * 3, result.Dna.Length);
        Assert.Equal(1.0, result.Cai, 9);
        Assert.Equal(Protein, Translate(result.Dna));
        Assert.Empty(result.RemainingViolations);

        // Every codon really is the host optimum.
        for (int i = 0; i < Protein.Length; i++)
            Assert.Equal(CodonUsageTable.EColiK12.MostFrequentCodon(Protein[i]), result.Dna.Substring(i * 3, 3));
    }

    [Fact]
    public void Optimize_WeightedRandomIsReproducibleAndStillSynonymous()
    {
        var options = new CodonOptimizationOptions
        {
            Target = CodonUsageTable.HomoSapiens,
            Method = CodonOptimizationMethod.WeightedRandom,
            Seed = 1234,
        };

        var a = CodonOptimizer.Optimize(Protein, options);
        var b = CodonOptimizer.Optimize(Protein, options);
        var c = CodonOptimizer.Optimize(Protein, new CodonOptimizationOptions
        {
            Target = CodonUsageTable.HomoSapiens,
            Method = CodonOptimizationMethod.WeightedRandom,
            Seed = 99,
        });

        Assert.Equal(a.Dna, b.Dna);
        Assert.NotEqual(a.Dna, c.Dna);
        Assert.Equal(Protein, Translate(a.Dna));
        Assert.Equal(Protein, Translate(c.Dna));
        // Sampling the distribution instead of collapsing onto the optimum must cost CAI.
        Assert.True(a.Cai < 1.0);
    }

    [Fact]
    public void OptimizeCoding_KeepsTheProteinAndImprovesAdaptation()
    {
        // A gene written entirely in E. coli's rarest codons, re-coded for E. coli.
        var host = CodonUsageTable.EColiK12;
        string poor = string.Concat(Protein.Select(aa => host.RarestCodon(aa)));
        Assert.Equal(Protein, Translate(poor));

        var result = CodonOptimizer.OptimizeCoding(poor, new CodonOptimizationOptions
        {
            Target = host,
            Method = CodonOptimizationMethod.MostFrequent,
        });

        Assert.Equal(Protein, Translate(result.Dna));
        Assert.True(result.Cai > host.Cai(poor));
        Assert.True(result.CodonsChanged > 0);
    }

    [Fact]
    public void Optimize_RepairsConstraintViolations()
    {
        var constraints = new SequenceConstraints
        {
            ForbiddenSites = ["GAATTC", "GGATCC"],   // EcoRI, BamHI
            MaxHomopolymer = 5,
            AvoidRepeats = false,
            AvoidHairpins = false,
        };

        var result = CodonOptimizer.Optimize(Protein, new CodonOptimizationOptions
        {
            Target = CodonUsageTable.EColiK12,
            Method = CodonOptimizationMethod.MostFrequent,
            Constraints = constraints,
        });

        Assert.Equal(Protein, Translate(result.Dna));
        Assert.DoesNotContain("GAATTC", result.Dna, StringComparison.Ordinal);
        Assert.DoesNotContain("GGATCC", result.Dna, StringComparison.Ordinal);
        Assert.Empty(result.RemainingViolations);
        // Repair costs adaptation, so the result can no longer be a perfect CAI of 1.
        Assert.InRange(result.Cai, 0.5, 1.0);
    }

    [Fact]
    public void Optimize_RejectsUnencodableProteinsAndHarmonisationWithoutASource()
    {
        var options = new CodonOptimizationOptions { Target = CodonUsageTable.EColiK12 };
        Assert.Throws<ArgumentException>(() => CodonOptimizer.Optimize("MKX", options));
        Assert.Throws<ArgumentException>(() => CodonOptimizer.Optimize("", options));
        Assert.Throws<InvalidOperationException>(() => CodonOptimizer.Optimize("MKA",
            new CodonOptimizationOptions
            {
                Target = CodonUsageTable.EColiK12,
                Method = CodonOptimizationMethod.Harmonised,
            }));
    }

    [Fact]
    public void OptimizeCoding_HarmonisationKeepsRareCodonsRare()
    {
        var source = CodonUsageTable.EColiK12;
        var target = CodonUsageTable.HomoSapiens;
        string gene = string.Concat(Protein.Select(aa => source.RarestCodon(aa)));

        var result = CodonOptimizer.OptimizeCoding(gene, new CodonOptimizationOptions
        {
            Target = target,
            Source = source,
            Method = CodonOptimizationMethod.Harmonised,
        });

        Assert.Equal(Protein, Translate(result.Dna));
        // Harmonisation reproduces the source's low adaptation rather than maximising it.
        Assert.True(result.Cai < 0.6, $"harmonised CAI was {result.Cai}");
    }

    [Fact]
    public void CodonFrequencyProfile_SumsToOne()
    {
        var profile = CodonOptimizer.CodonFrequencyProfile("CTGCTGTTATTC");
        Assert.Equal(1.0, profile.Values.Sum(), 12);
        Assert.Equal(0.5, profile["CTG"], 12);
        Assert.Equal(0.25, profile["TTA"], 12);
    }

    [Fact]
    public void BackTranslate_IsSynonymous()
    {
        string dna = CodonOptimizer.BackTranslate("MKAIL*", CodonUsageTable.EColiK12,
            CodonOptimizationMethod.MostFrequent, seed: 1);
        Assert.Equal("MKAIL*", Translate(dna));
    }

    // ---- constraints ----

    [Fact]
    public void Constraints_FindForbiddenSitesWithTheirPositions()
    {
        var constraints = new SequenceConstraints
        {
            ForbiddenSites = ["GAATTC"],
            AvoidHairpins = false,
            AvoidRepeats = false,
            MaxHomopolymer = 0,
            MinGc = 0,
            MaxGc = 1,
        };

        var violations = constraints.Check("AAAGAATTCAAA");
        var violation = Assert.Single(violations);
        Assert.Equal(SequenceConstraints.KindForbiddenSite, violation.Kind);
        Assert.Equal(3, violation.Start);
        Assert.Equal(6, violation.Length);

        Assert.True(constraints.IsSatisfiedBy("AAAGAATTAAAA"));
    }

    [Fact]
    public void Constraints_FindHomopolymersAndGcWindows()
    {
        var homopolymer = new SequenceConstraints
        {
            MaxHomopolymer = 5,
            AvoidHairpins = false,
            AvoidRepeats = false,
            MinGc = 0,
            MaxGc = 1,
        };
        Assert.Contains(homopolymer.Check("ACGTAAAAAAACGT"), v => v.Kind == SequenceConstraints.KindHomopolymer);
        Assert.True(homopolymer.IsSatisfiedBy("ACGTAAAAACGT"));

        var gc = new SequenceConstraints
        {
            MinGc = 0.3,
            MaxGc = 0.7,
            GcWindow = 20,
            MaxHomopolymer = 0,
            AvoidHairpins = false,
            AvoidRepeats = false,
        };
        Assert.Contains(gc.Check(new string('G', 40)), v => v.Kind == SequenceConstraints.KindGcWindow);
        Assert.True(gc.IsSatisfiedBy(string.Concat(Enumerable.Repeat("ACGT", 10))));
    }

    [Fact]
    public void Constraints_FindHairpinsAndDirectRepeats()
    {
        var hairpin = new SequenceConstraints
        {
            AvoidHairpins = true,
            MinHairpinStem = 8,
            MinHairpinLoop = 3,
            AvoidRepeats = false,
            MaxHomopolymer = 0,
            MinGc = 0,
            MaxGc = 1,
        };
        // A 10 bp stem separated by a 4 nt loop.
        const string stem = "ACGTACGTAC";
        string folded = stem + "TTTT" + Biotech.Net.SyntheticBiology.DnaUtility.ReverseComplement(stem);
        Assert.Contains(hairpin.Check(folded), v => v.Kind == SequenceConstraints.KindHairpin);

        var repeat = new SequenceConstraints
        {
            AvoidRepeats = true,
            MaxDirectRepeat = 12,
            AvoidHairpins = false,
            MaxHomopolymer = 0,
            MinGc = 0,
            MaxGc = 1,
        };
        const string unit = "ACGTACGTACGTAC";
        Assert.Contains(repeat.Check(unit + "TTTT" + unit), v => v.Kind == SequenceConstraints.KindDirectRepeat);
    }

    [Fact]
    public void GeneSynthesisPreset_IsStricterThanTheDefault()
    {
        Assert.True(SequenceConstraints.ForGeneSynthesis.ForbiddenSites.Count >= SequenceConstraints.Default.ForbiddenSites.Count);
        Assert.NotEmpty(SequenceConstraints.CommonRestrictionSites);
    }

    private static string Translate(string dna)
    {
        var code = GeneticCode.Standard;
        return string.Concat(Enumerable.Range(0, dna.Length / 3).Select(i => code.Translate(dna.AsSpan(i * 3, 3))));
    }
}
