using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics.Tests;

public class NormalizationTests
{
    private static Dataset Make(double[,] values) => Dataset.FromArray("t", values);

    [Fact]
    public void MedianOfRatios_MatchesHandComputation()
    {
        // Geometric means: g1 = 200, g2 = 30, g3 = 20; g4 has a zero and is ignored.
        // Ratios s1: 0.5, 1, 0.5 → 0.5; s2: 1, 1, 2 → 1; s3: 2, 1, 1 → 1.
        var counts = Make(new double[,] { { 100, 200, 400 }, { 30, 30, 30 }, { 10, 40, 20 }, { 0, 5, 7 } });
        var sf = Normalization.MedianOfRatiosSizeFactors(counts);
        Assert.Equal([0.5, 1, 1], sf.Select(v => Math.Round(v, 12)).ToArray());
    }

    [Fact]
    public void MedianOfRatios_IsEquivariantToColumnScaling()
    {
        var (counts, _, _) = TestData.TwoGroupCounts(300, 3, 0, 0, seed: 1);
        var sf = Normalization.MedianOfRatiosSizeFactors(counts);
        double[] k = [1, 2, 0.5, 3, 1, 4];
        var scaled = counts.Map(v => v); // copy
        for (int i = 0; i < scaled.Rows; i++)
            for (int j = 0; j < scaled.Columns; j++) scaled[i, j] *= k[j];
        var sf2 = Normalization.MedianOfRatiosSizeFactors(scaled);
        // Every ratio of sample j is multiplied by k_j / geomean(k).
        double g = Math.Exp(k.Average(Math.Log));
        for (int j = 0; j < sf.Length; j++) Assert.Equal(sf[j] * k[j] / g, sf2[j], 1e-9);
        // The product of DESeq2 size factors is ~1 by construction for balanced data.
        Assert.InRange(Math.Exp(sf.Average(Math.Log)), 0.8, 1.25);
    }

    [Fact]
    public void MedianOfRatios_ThrowsWhenEveryGeneHasAZero()
    {
        var counts = Make(new double[,] { { 0, 1 }, { 2, 0 } });
        Assert.Throws<InvalidOperationException>(() => Normalization.MedianOfRatiosSizeFactors(counts));
    }

    [Fact]
    public void Tmm_IdenticalOrProportionalLibrariesGiveOne()
    {
        var r = new Random(5);
        var values = new double[200, 3];
        for (int i = 0; i < 200; i++)
        {
            double v = 1 + r.Next(1000);
            values[i, 0] = v; values[i, 1] = v; values[i, 2] = 3 * v; // third library is a deeper copy
        }
        var f = Normalization.TmmFactors(Make(values));
        foreach (double x in f) Assert.Equal(1, x, 1e-12);
    }

    [Fact]
    public void Tmm_CorrectsCompositionBias()
    {
        // 100 genes with 100 reads in both samples, but one gene in sample 2 soaks up 10 000 extra reads.
        // Library sizes 10 000 and 20 000; the 99 ordinary genes have M = log2(0.005 / 0.01) = −1, so the raw
        // factor of sample 2 relative to sample 1 is 2^−1, and after scaling to geometric mean 1: √2 and 1/√2.
        var values = new double[100, 2];
        for (int i = 0; i < 100; i++) { values[i, 0] = 100; values[i, 1] = 100; }
        values[0, 1] = 10100;
        var counts = Make(values);
        var f = Normalization.TmmFactors(counts);
        Assert.Equal(Math.Sqrt(2), f[0], 1e-12);
        Assert.Equal(1 / Math.Sqrt(2), f[1], 1e-12);
        Assert.Equal(1, f[0] * f[1], 1e-12);
        var eff = Normalization.EffectiveLibrarySizes(counts, f);
        Assert.Equal(eff[0], eff[1], 1e-6);
        // Plain library-size normalisation (upper quartile is fooled less, but None is fooled completely).
        Assert.Equal([1.0, 1.0], Normalization.Factors(counts, NormalizationMethod.None));
    }

    [Fact]
    public void Tmm_ProductOfFactorsIsOneOnRealisticData()
    {
        var (counts, _, _) = TestData.TwoGroupCounts(500, 4, 50, 2, seed: 2);
        foreach (var method in new[] { NormalizationMethod.Tmm, NormalizationMethod.UpperQuartile, NormalizationMethod.MedianOfRatios })
        {
            var f = Normalization.Factors(counts, method);
            Assert.Equal(0, f.Sum(Math.Log), 1e-10);
            Assert.All(f, x => Assert.InRange(x, 0.5, 2));
        }
    }

    [Fact]
    public void UpperQuartile_MatchesHandComputation()
    {
        // Proportions per sample: s1 = (1,2,3,4)/10, s2 = (2,2,2,2)/8. 75% quantile (type 7): s1 0.325, s2 0.25.
        var counts = Make(new double[,] { { 1, 2 }, { 2, 2 }, { 3, 2 }, { 4, 2 } });
        var f = Normalization.UpperQuartileFactors(counts);
        double g = Math.Sqrt(0.325 * 0.25);
        Assert.Equal(0.325 / g, f[0], 1e-12);
        Assert.Equal(0.25 / g, f[1], 1e-12);
    }

    [Fact]
    public void CpmTpmFpkm_MatchHandComputation()
    {
        // g1 = (10, 0), g2 = (30, 50); lengths 1 kb and 2 kb; libraries 40 and 50.
        var counts = Make(new double[,] { { 10, 0 }, { 30, 50 } });
        double[] len = [1000, 2000];

        var cpm = Normalization.Cpm(counts);
        Assert.Equal(250_000, cpm[0, 0], 1e-6);
        Assert.Equal(750_000, cpm[1, 0], 1e-6);
        Assert.Equal(0, cpm[0, 1]);
        Assert.Equal(1_000_000, cpm[1, 1], 1e-6);

        // RPK s1 = 10, 15 (sum 25); s2 = 0, 25.
        var tpm = Normalization.Tpm(counts, len);
        Assert.Equal(400_000, tpm[0, 0], 1e-6);
        Assert.Equal(600_000, tpm[1, 0], 1e-6);
        Assert.Equal(0, tpm[0, 1]);
        Assert.Equal(1_000_000, tpm[1, 1], 1e-6);

        var fpkm = Normalization.Fpkm(counts, len);
        Assert.Equal(10 / 1.0 / (40 / 1e6), fpkm[0, 0], 1e-6);   // 250 000
        Assert.Equal(30 / 2.0 / (40 / 1e6), fpkm[1, 0], 1e-6);   // 375 000
        Assert.Equal(50 / 2.0 / (50 / 1e6), fpkm[1, 1], 1e-6);   // 500 000
        Assert.Equal(fpkm.Values.ToArray(), Normalization.Rpkm(counts, len).Values.ToArray());

        // With normalisation factors the effective library shrinks/grows accordingly.
        var cpmNf = Normalization.Cpm(counts, [2, 0.5]);
        Assert.Equal(125_000, cpmNf[0, 0], 1e-6);

        // TPM via the length annotation.
        counts.WithGeneLengths(len);
        Assert.Equal(tpm.Values.ToArray(), Normalization.Tpm(counts).Values.ToArray());
    }

    [Fact]
    public void LogCpm_UsesEdgeRPriorCount()
    {
        var counts = Make(new double[,] { { 10, 0 }, { 30, 50 } });
        var lcpm = Normalization.LogCpm(counts, priorCount: 2);
        // lib = (40, 50), mean 45: prior_1 = 2·40/45 = 16/9, adjusted lib_1 = 40 + 32/9.
        double p1 = 16.0 / 9, l1 = 40 + 32.0 / 9, p2 = 20.0 / 9, l2 = 50 + 40.0 / 9;
        Assert.Equal(Math.Log2((10 + p1) / l1 * 1e6), lcpm[0, 0], 1e-12);
        Assert.Equal(Math.Log2((0 + p2) / l2 * 1e6), lcpm[0, 1], 1e-12);
        Assert.Equal(Math.Log2((50 + p2) / l2 * 1e6), lcpm[1, 1], 1e-12);
    }

    [Fact]
    public void NormalizedCounts_DivideBySizeFactor()
    {
        var counts = Make(new double[,] { { 10, 20 }, { 4, 8 } });
        var norm = Normalization.NormalizedCounts(counts, [0.5, 2]);
        Assert.Equal([20.0, 10, 8, 4], norm.Values.ToArray());
    }

    [Fact]
    public void FilterByExpr_FollowsEdgeRRules()
    {
        // Libraries are all 1e6, so min.count = 10 → CPM cutoff 10. Smallest group has 2 samples.
        var values = new double[,]
        {
            { 20, 20, 0, 0, 999_960 / 4.0 },
            { 20, 0, 0, 0, 0 },
            { 5, 5, 5, 5, 0 },
            { 0, 0, 15, 15, 0 },
        };
        // Make all columns sum to 1e6 with a filler gene.
        var v = new double[5, 4];
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) v[i, j] = values[i, j];
        for (int j = 0; j < 4; j++)
        {
            double s = 0;
            for (int i = 0; i < 4; i++) s += v[i, j];
            v[4, j] = 1e6 - s;
        }
        var counts = Make(v);
        var keep = ExpressionFilter.FilterByExpr(counts, ["a", "a", "b", "b"]);
        Assert.Equal([true, false, false, true, true], keep);
        var keepAll = ExpressionFilter.FilterByExpr(counts); // one group of 4
        Assert.Equal([false, false, false, false, true], keepAll);
        var design = DesignMatrix.TwoGroup(4, [0, 1], [2, 3]);
        Assert.Equal(keep, ExpressionFilter.FilterByExpr(counts, design)); // 1 / max leverage = 2
        Assert.Equal(3, ExpressionFilter.Apply(counts, ExpressionFilter.MinimumCount(counts, 15, 2)).Rows);
    }
}
