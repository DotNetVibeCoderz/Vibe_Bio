using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics.Tests;

public class LimmaTests
{
    [Fact]
    public void LmFit_TwoGroupGivesMeanDifferenceAndPooledSigma()
    {
        var y = Dataset.FromArray("y", new double[,] { { 1, 2, 3, 6, 7, 8 }, { 5, 5, 5, 5, 5, 6 } });
        var design = DesignMatrix.TwoGroup(6, [0, 1, 2], [3, 4, 5]);
        var fit = Limma.LmFit(y, design);
        Assert.Equal(2, fit.Coefficients[0], 1e-12);   // intercept = mean of group A
        Assert.Equal(5, fit.Coefficients[1], 1e-12);   // 7 − 2
        Assert.Equal(1, fit.Sigma[0], 1e-12);          // pooled SD: (2 + 2) / 4 = 1
        Assert.Equal(4, fit.DfResidual[0]);
        Assert.Equal(Math.Sqrt(2.0 / 3), fit.StdevUnscaled[1], 1e-12);
        Assert.Equal(4.5, fit.Amean[0], 1e-12);
        // Second gene: group B = (5, 5, 6), RSS = 2/3.
        Assert.Equal(1.0 / 3, fit.Coefficients[3], 1e-12);
        Assert.Equal(Math.Sqrt(2.0 / 3 / 4), fit.Sigma[1], 1e-12);
    }

    [Fact]
    public void FitFDist_RecoversScaledFPrior()
    {
        // σ² ~ s0² d0 / χ²(d0), s² = σ² χ²(d) / d.
        var r = new Random(21);
        const double d0 = 4, s20 = 0.5, d = 6;
        int n = 20000;
        var s2 = new double[n]; var df = new double[n];
        for (int i = 0; i < n; i++)
        {
            double sigma2 = s20 * d0 / TestData.ChiSquare(r, d0);
            s2[i] = sigma2 * TestData.ChiSquare(r, d) / d;
            df[i] = d;
        }
        var (estD0, estS20) = Limma.FitFDist(s2, df);
        Assert.InRange(estD0, 3.2, 5.0);
        Assert.InRange(estS20, 0.45, 0.55);
    }

    [Fact]
    public void FitFDist_NoExtraVariationGivesInfinitePrior()
    {
        var r = new Random(22);
        int n = 5000;
        var s2 = new double[n]; var df = new double[n];
        for (int i = 0; i < n; i++) { s2[i] = 2.0 * TestData.ChiSquare(r, 5) / 5; df[i] = 5; }
        var (d0, s20) = Limma.FitFDist(s2, df);
        Assert.True(d0 > 50, $"d0 = {d0}");
        Assert.InRange(s20, 1.9, 2.1);
    }

    [Fact]
    public void ModeratedT_ShrinksVariancesAndUsesTotalDf()
    {
        var r = new Random(23);
        int n = 500, m = 6;
        var values = new double[n * m];
        for (int i = 0; i < values.Length; i++) values[i] = TestData.Normal(r) * Math.Sqrt(0.3 + 0.2 * (i / m % 3));
        var y = new Dataset("y", Enumerable.Range(0, n).Select(i => $"g{i}").ToArray(), Enumerable.Range(0, m).Select(j => $"s{j}").ToArray(), values);
        var design = DesignMatrix.TwoGroup(m, [0, 1, 2], [3, 4, 5]);
        var fit = Limma.EBayes(Limma.LmFit(y, design));
        Assert.True(fit.DfPrior > 0);
        for (int i = 0; i < n; i++)
        {
            double s2 = fit.Sigma[i] * fit.Sigma[i];
            double post = fit.S2Post[i];
            Assert.InRange(post, Math.Min(s2, fit.S2Prior) - 1e-12, Math.Max(s2, fit.S2Prior) + 1e-12);
        }
        var tt = Limma.TopTable(fit, 1);
        double dfTotal = Math.Min(4 + fit.DfPrior, 4.0 * n);
        Assert.Equal(Biotech.Net.Native.SpecialFunctions.StudentTTwoSided(tt.Stat[0], dfTotal), tt.PValue[0], 1e-12);
        // Under the null p-values are roughly uniform.
        double frac = tt.PValue.Count(p => p < 0.05) / (double)n;
        Assert.InRange(frac, 0.01, 0.1);
    }

    [Fact]
    public void Voom_RecoversPlantedDifferentialGenes()
    {
        var (counts, trueLfc, _) = TestData.TwoGroupCounts(2000, 4, 100, 1.5, seed: 31);
        var design = DesignMatrix.FromAnnotations(counts, "condition", "A");
        var nf = Normalization.TmmFactors(counts);
        var res = Limma.VoomTest(counts, design, normFactors: nf);
        Assert.Equal("limma-voom", res.Method);

        var called = Enumerable.Range(0, 2000).Where(i => res.PAdj[i] < 0.05).ToArray();
        int tp = called.Count(i => i < 100);
        Assert.True(tp >= 80, $"recall {tp}/100");
        Assert.True(called.Length - tp <= 0.1 * called.Length + 2, $"{called.Length - tp} false of {called.Length}");
        foreach (int i in called.Where(i => i < 100)) Assert.Equal(Math.Sign(trueLfc[i]), Math.Sign(res.Log2FoldChange[i]));

        var top100 = Enumerable.Range(0, 2000).OrderBy(i => res.PValue[i]).Take(100).Count(i => i < 100);
        Assert.True(top100 >= 80, $"{top100} of the top 100 are planted");
    }

    [Fact]
    public void Voom_WeightsFollowTheMeanVarianceTrend()
    {
        var (counts, _, _) = TestData.TwoGroupCounts(1000, 3, 0, 0, seed: 32);
        var design = DesignMatrix.FromAnnotations(counts, "condition", "A");
        var v = Limma.Voom(counts, design);
        Assert.Equal(1000 * 6, v.Weights.Length);
        Assert.All(v.Weights, w => Assert.True(w > 0 && double.IsFinite(w)));
        // Low-count observations are noisier on the log scale, so they get smaller weights.
        int lo = Enumerable.Range(0, 1000).OrderBy(i => v.TrendX[i]).First();
        int hi = Enumerable.Range(0, 1000).OrderBy(i => v.TrendX[i]).Last();
        Assert.True(v.Weights[lo * 6] < v.Weights[hi * 6]);
        // log-CPM uses limma's offsets: log2((count + 0.5) / (lib + 1) · 1e6).
        double lib0 = 0;
        for (int i = 0; i < 1000; i++) lib0 += counts[i, 0];
        Assert.Equal(Math.Log2((counts[0, 0] + 0.5) / (lib0 + 1) * 1e6), v.LogCpm[0, 0], 1e-10);
    }

    [Fact]
    public void TTest_WelchAndStudentAgreeForEqualVariances()
    {
        var y = Dataset.FromArray("y", new double[,] { { 1, 2, 3, 6, 7, 8 } });
        var welch = DifferentialExpression.TTest(y, [0, 1, 2], [3, 4, 5]);
        var student = DifferentialExpression.TTest(y, [0, 1, 2], [3, 4, 5], welch: false);
        // t = 5 / sqrt(1/3 + 1/3) with df 4; R: t.test(c(6,7,8), c(1,2,3), var.equal=TRUE)$p.value = 0.002535...
        double t = 5 / Math.Sqrt(2.0 / 3);
        Assert.Equal(t, student.Stat[0], 1e-12);
        Assert.Equal(t, welch.Stat[0], 1e-9);
        Assert.Equal(Biotech.Net.Native.SpecialFunctions.StudentTTwoSided(t, 4), student.PValue[0], 1e-12);
        Assert.Equal(student.PValue[0], welch.PValue[0], 1e-9);
        Assert.Equal(5, welch.Log2FoldChange[0], 1e-12);
    }
}
