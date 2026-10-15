using Biotech.Net.Core;
using Biotech.Net.Genomics.Population;

namespace Biotech.Net.Genomics.Tests;

public class PopulationGeneticsTests
{
    [Theory]
    // Reference p-values from direct enumeration of the exact conditional distribution of heterozygotes
    // (Levene 1949 / Wigginton et al. 2005, eq. 2), computed independently in Python with log-gamma.
    [InlineData(57, 14, 50, 0.8422797565707778)]
    [InlineData(10, 40, 50, 3.273808610043725e-17)]
    [InlineData(0, 50, 50, 1.1142241805814346e-30)]
    [InlineData(30, 35, 35, 5.882470534402013e-05)]
    [InlineData(1, 99, 0, 1.0)]
    [InlineData(500, 250, 250, 1.0)]
    [InlineData(3, 10, 2, 0.12337164750957792)]
    public void HardyWeinbergExactTestMatchesEnumeration(int het, int homRef, int homAlt, double expected)
    {
        double p = HardyWeinberg.ExactTest(het, homRef, homAlt);
        Assert.Equal(expected, p, expected * 1e-8 + 1e-15);
    }

    [Fact]
    public void HardyWeinbergExactTestIsSymmetricAndHandlesEdges()
    {
        Assert.Equal(HardyWeinberg.ExactTest(57, 14, 50), HardyWeinberg.ExactTest(57, 50, 14), 12);
        Assert.Equal(1.0, HardyWeinberg.ExactTest(0, 0, 0));
        Assert.Equal(1.0, HardyWeinberg.ExactTest(0, 100, 0));
        var (chi, p) = HardyWeinberg.ChiSquareTest(30, 35, 35);
        // Expected 50 hets at p = 0.5; χ² = (35−25)²/25·2 + (30−50)²/50 = 16; P(χ²₁ > 16) = 6.334e-5.
        Assert.Equal(16.0, chi, 9);
        Assert.Equal(6.334e-5, p, 3e-7);
    }

    [Fact]
    public void GenotypeCountsAndFrequenciesFromVariant()
    {
        var v = new Variant("1", 10, "A", "G");
        foreach (var gt in new[] { "0/0", "0/1", "1/1", "0|1", "./.", "0/0" }) v.Genotypes.Add(Genotype.Parse(gt));
        var c = GenotypeCounts.From(v);
        Assert.Equal(new GenotypeCounts(2, 2, 1, 1), c);
        Assert.Equal(0.4, c.AlternateFrequency, 12);
        Assert.Equal(0.4, c.ObservedHeterozygosity, 12);
        Assert.Equal(0.48, c.ExpectedHeterozygosity, 12);
        Assert.Equal([0.6, 0.4], AlleleFrequencies.Of(v));
        Assert.Equal(0.4, AlleleFrequencies.ObservedHeterozygosity(v), 12);
        // Nei's unbiased gene diversity, 2N = 10 alleles: 10/9 · 0.48.
        Assert.Equal(0.48 * 10 / 9, AlleleFrequencies.ExpectedHeterozygosity([0.6, 0.4], 10), 12);
        Assert.Equal(GenotypeCounts.From(v, [0, 1]), new GenotypeCounts(1, 1, 0, 0));
    }

    private static readonly string[] Alignment =
    [
        "ATGCTAGCTAGGCTAACGTA", "ATGCTAGCTTGGCTAACGTA", "ATGCAAGCTAGGCTTACGTA", "ATGCTAGCTAGGCTAACGGA",
        "TTGCTAGCTAGGCTAACGTA", "ATGCTAGCTAGCCTAACGTA", "ATGCTAGCTTGGCTAACGTA",
    ];

    [Fact]
    public void NucleotideDiversityWattersonAndTajimasD()
    {
        var seqs = Alignment.Select((s, i) => Sequence.Dna(s, $"s{i}")).ToList();
        var d = Diversity.Summarize(seqs);
        // Reference values computed independently (pairwise comparisons, Tajima 1989 eqs. 38).
        Assert.Equal(7, d.SampleSize);
        Assert.Equal(20, d.Sites);
        Assert.Equal(6, d.SegregatingSites);
        Assert.Equal(1.9047619047619047, d.MeanPairwiseDifferences, 12);
        Assert.Equal(0.09523809523809523, d.Pi, 12);
        Assert.Equal(2.4489795918367343, d.ThetaW, 12);
        Assert.Equal(-1.1289806675140472, d.TajimasD, 12);
        // Harmonic numbers for n = 10 (Tajima 1989): a1 = 2.828968, a2 = 1.539768.
        Assert.Equal(2.828968, Diversity.A1(10), 6);
        Assert.Equal(1.539768, Diversity.A2(10), 6);
        Assert.True(double.IsNaN(Diversity.TajimasD(10, 0, 0)));
    }

    [Fact]
    public void SfsStatisticsAgreeWithAlignmentStatistics()
    {
        // Derived-allele counts per segregating site of the alignment above (minor allele = derived here).
        var seqs = Alignment;
        var counts = new List<int>();
        for (int j = 0; j < seqs[0].Length; j++)
        {
            var col = seqs.Select(s => s[j]).ToList();
            var major = col.GroupBy(c => c).MaxBy(g => g.Count())!.Key;
            int derived = col.Count(c => c != major);
            if (derived > 0) counts.Add(derived);
        }
        var sfs = SiteFrequencySpectrum.Unfolded(counts, 7);
        Assert.Equal(6, sfs.Skip(1).Take(6).Sum());
        Assert.Equal(1.9047619047619047, SiteFrequencySpectrum.ThetaPi(sfs), 12);
        Assert.Equal(2.4489795918367343, SiteFrequencySpectrum.ThetaW(sfs), 12);
        Assert.Equal(-1.1289806675140472, SiteFrequencySpectrum.TajimasD(sfs), 12);
        Assert.Equal(SiteFrequencySpectrum.Folded(counts, 7), SiteFrequencySpectrum.Fold(sfs));
        Assert.Equal([0L, 2, 1, 1], SiteFrequencySpectrum.Fold([0, 2, 1, 1, 0, 0, 0, 0]));
    }

    [Fact]
    public void DiversityFromVariantsMatchesHaplotypeStatistics()
    {
        // Encode the alignment's 6 segregating sites as phased haploid-pairs VCF genotypes (7 haplotypes → use 7 haploid samples).
        var variants = new List<Variant>();
        for (int j = 0; j < Alignment[0].Length; j++)
        {
            var col = Alignment.Select(s => s[j]).ToList();
            var major = col.GroupBy(c => c).MaxBy(g => g.Count())!.Key;
            if (col.All(c => c == major)) continue;
            char minor = col.First(c => c != major);
            var v = new Variant("1", j + 1, major.ToString(), minor.ToString());
            foreach (char c in col) v.Genotypes.Add(new Genotype([c == major ? 0 : 1]));
            variants.Add(v);
        }
        Assert.Equal(0.09523809523809523, Diversity.NucleotideDiversity(variants, 20), 12);
        Assert.Equal(-1.1289806675140472, Diversity.TajimasD(variants), 12);
        var sfs = SiteFrequencySpectrum.FromVariants(variants);
        Assert.Equal(8, sfs.Length);
        Assert.Equal(6, sfs[1..7].Sum());
    }

    [Fact]
    public void HudsonFst()
    {
        string[] p1 = ["ACGTACGTAC", "ACGTACGTAA", "ACGAACGTAC", "ACGTACGTAC"];
        string[] p2 = ["TCGTACCTAC", "TCGTACCTGC", "TCGTACCTAC", "ACGTACCTAC"];
        double f = Fst.Hudson(p1.Select(s => Sequence.Dna(s)).ToList(), p2.Select(s => Sequence.Dna(s)).ToList());
        Assert.Equal(0.6, f, 12); // 1 − Hw/Hb, computed independently

        // Bhatia et al. (2013) ratio-of-averages estimator from allele counts.
        (int, int, int, int)[] sites = [(3, 20, 15, 20), (10, 40, 12, 40), (0, 30, 5, 30), (8, 16, 8, 16)];
        Assert.Equal(0.18682330182951057, Fst.Hudson(sites), 12);
        // Fixed differences: Fst = 1.
        Assert.Equal(1.0, Fst.Hudson([(0, 50, 50, 50)]), 12);
    }

    [Fact]
    public void WeirCockerhamFst()
    {
        GenotypeCounts G(int hr, int het, int ha) => new(hr, het, ha);
        IReadOnlyList<GenotypeCounts>[] loci =
        [
            [G(30, 15, 5), G(10, 20, 20)],
            [G(40, 8, 2), G(35, 12, 3)],
            [G(10, 25, 15), G(20, 20, 10), G(5, 10, 35)],
        ];
        // Reference values from an independent implementation of Weir & Cockerham (1984) eqs. 2–4.
        Assert.Equal(0.21335807050092762, Fst.WeirCockerhamComponents(loci[0]).Theta, 12);
        Assert.Equal(0.0018437281454637953, Fst.WeirCockerhamComponents(loci[1]).Theta, 12);
        Assert.Equal(0.14921812880996557, Fst.WeirCockerhamComponents(loci[2]).Theta, 12);
        Assert.Equal(0.14727228129914727, Fst.WeirCockerham(loci), 12);

        // From VCF genotypes: two fully differentiated populations give θ = 1.
        var v = new Variant("1", 1, "A", "T");
        for (int i = 0; i < 10; i++) v.Genotypes.Add(new Genotype([0, 0]));
        for (int i = 0; i < 10; i++) v.Genotypes.Add(new Genotype([1, 1]));
        Assert.Equal(1.0, Fst.WeirCockerham([v], Enumerable.Range(0, 10).ToArray(), Enumerable.Range(10, 10).ToArray()), 12);
        Assert.Equal(1.0, Fst.Hudson([v], Enumerable.Range(0, 10).ToArray(), Enumerable.Range(10, 10).ToArray()), 12);
    }

    [Fact]
    public void LinkageDisequilibrium()
    {
        int[] a = [1, 1, 1, 0, 0, 0, 0, 0, 1, 0];
        int[] b = [1, 1, 0, 0, 0, 0, 0, 0, 1, 0];
        // pA = 0.4, pB = 0.3, pAB = 0.3 → D = 0.18, Dmax = min(0.28, 0.18) → D′ = 1, r² = 0.0324 / 0.0504.
        var ld = Population.LinkageDisequilibrium.FromHaplotypes(a, b);
        Assert.Equal(0.18, ld.D, 12);
        Assert.Equal(1.0, ld.DPrime, 12);
        Assert.Equal(0.0324 / 0.0504, ld.RSquared, 12);

        // Perfect LD and linkage equilibrium.
        Assert.Equal(1.0, Population.LinkageDisequilibrium.FromHaplotypes(a, a).RSquared, 12);
        var eq = Population.LinkageDisequilibrium.FromCounts(100, 50, 50, 25);
        Assert.Equal(0.0, eq.D, 12);
        Assert.Equal(0.0, eq.RSquared, 12);

        // Same data as phased diploid genotypes.
        var va = new Variant("1", 1, "A", "G");
        var vb = new Variant("1", 2, "C", "T");
        for (int i = 0; i < 10; i += 2)
        {
            va.Genotypes.Add(new Genotype([a[i], a[i + 1]], true));
            vb.Genotypes.Add(new Genotype([b[i], b[i + 1]], true));
        }
        var ld2 = Population.LinkageDisequilibrium.FromPhased(va, vb);
        Assert.Equal(ld.RSquared, ld2.RSquared, 12);
        Assert.Equal(10, ld2.Haplotypes);
        Assert.InRange(Population.LinkageDisequilibrium.GenotypeRSquared(va, vb), 0, 1);
        Assert.Equal(1.0, Population.LinkageDisequilibrium.GenotypeRSquared(va, va), 12);
    }
}
