using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics.Tests;

public class DESeqTests
{
    internal static readonly Lazy<(Dataset Counts, double[] TrueLfc, DESeqFit Fit, DifferentialExpressionResult Res)> Sim = new(() =>
    {
        var (counts, lfc, _) = TestData.TwoGroupCounts(1000, 3, 100, 2.0, seed: 11);
        var design = DesignMatrix.FromAnnotations(counts, "condition", "A");
        var fit = DESeq.Fit(counts, design);
        return (counts, lfc, fit, fit.Results());
    });

    [Fact]
    public void Fit_UsesMedianOfRatiosSizeFactors()
    {
        var (counts, _, fit, _) = Sim.Value;
        Assert.Equal(Normalization.MedianOfRatiosSizeFactors(counts), fit.SizeFactors);
        for (int i = 0; i < 5; i++)
        {
            double mean = 0;
            for (int j = 0; j < counts.Columns; j++) mean += counts[i, j] / fit.SizeFactors[j];
            Assert.Equal(mean / counts.Columns, fit.BaseMean[i], 1e-9);
        }
    }

    [Fact]
    public void PlantedDifferentialGenes_AreRecovered()
    {
        var (_, trueLfc, _, res) = Sim.Value;
        Assert.Equal("condition_B_vs_A", res.Contrast);
        int de = 0, found = 0, signOk = 0;
        double absErr = 0;
        for (int i = 0; i < 100; i++)
        {
            de++;
            if (res.PAdj[i] < 0.1) found++;
            if (Math.Sign(res.Log2FoldChange[i]) == Math.Sign(trueLfc[i])) signOk++;
            absErr += Math.Abs(res.Log2FoldChange[i] - trueLfc[i]);
        }
        Assert.True(found >= 90, $"Only {found}/{de} planted genes were called.");
        Assert.Equal(de, signOk);
        Assert.True(absErr / de < 0.35, $"Mean |LFC error| {absErr / de:F3}");
        // A clearly changed gene has a tiny p-value.
        Assert.True(res.Top(1)[0].PValue < 1e-8);
    }

    [Fact]
    public void NullGenes_HaveCalibratedPValues()
    {
        var (_, _, _, res) = Sim.Value;
        var nullP = Enumerable.Range(100, 900).Select(i => res.PValue[i]).Where(p => !double.IsNaN(p)).ToArray();
        double below05 = nullP.Count(p => p < 0.05) / (double)nullP.Length;
        double below50 = nullP.Count(p => p < 0.5) / (double)nullP.Length;
        Assert.InRange(below05, 0.01, 0.10);
        Assert.InRange(below50, 0.40, 0.60);
        int falseCalls = Enumerable.Range(100, 900).Count(i => res.PAdj[i] < 0.1);
        Assert.True(falseCalls <= 15, $"{falseCalls} null genes called at FDR 0.1");
    }

    [Fact]
    public void DispersionTrend_RecoversSimulatedParameters()
    {
        var (_, _, fit, _) = Sim.Value;
        Assert.Equal(DispersionFitType.Parametric, fit.TrendType);
        // Simulated α(μ) = 0.05 + 1/μ on the base-mean scale.
        Assert.InRange(fit.Trend.AsymptoticDispersion, 0.02, 0.1);
        Assert.InRange(fit.Trend.ExtraPoisson, 0.4, 2.5);
        Assert.True(fit.DispersionPriorVariance >= 0.25);
        // MAP dispersions lie between gene-wise estimates and the trend (in log space) for non-outliers.
        int checkedGenes = 0;
        for (int i = 0; i < fit.Dispersions.Length; i++)
        {
            if (fit.DispersionOutliers[i] || fit.AllZero[i]) continue;
            double g = Math.Log(fit.GeneDispersions[i]), t = Math.Log(fit.TrendDispersions[i]), d = Math.Log(fit.Dispersions[i]);
            if (fit.GeneDispersions[i] < 1e-6) continue;
            Assert.InRange(d, Math.Min(g, t) - 1e-3, Math.Max(g, t) + 1e-3);
            checkedGenes++;
        }
        Assert.True(checkedGenes > 500);
    }

    [Fact]
    public void GroupDesign_CoefficientsSolveTheNegativeBinomialScoreEquations()
    {
        // The fitted coefficients are a maximum-likelihood solution, so the negative-binomial score
        //     S_k = sum_j x_jk (y_j - mu_j) / (1 + alpha mu_j),   mu_j = s_j exp(x_j . beta)
        // must vanish for every coefficient (up to the 1e-6 ridge DESeq2 also carries). Note the weight
        // 1/(1 + alpha mu_j): it is why the negative-binomial fit is NOT simply the ratio of the normalised
        // group means — that identity is the Poisson one, and it only reappears when the size factors are equal
        // or the dispersion is zero.
        var (counts, _, fit, _) = Sim.Value;
        var design = fit.Design;
        int m = design.Samples, p = design.Coefficients;
        int checkedGenes = 0;

        for (int i = 0; i < 200; i++)
        {
            if (fit.AllZero[i] || !fit.Converged[i]) continue;

            var mu = new double[m];
            for (int j = 0; j < m; j++)
            {
                double eta = 0;
                for (int k = 0; k < p; k++) eta += design[j, k] * fit.Coefficients[i * p + k];
                mu[j] = fit.SizeFactors[j] * Math.Pow(2, eta);
            }
            if (mu.Any(v => v < 1.0)) continue; // the minimum-mu clamp would hold the fit off the score root

            double alpha = fit.Dispersions[i];
            for (int k = 0; k < p; k++)
            {
                double score = 0, scale = 0;
                for (int j = 0; j < m; j++)
                {
                    double w = 1 + alpha * mu[j];
                    score += design[j, k] * (counts[i, j] - mu[j]) / w;
                    scale += Math.Abs(design[j, k]) * mu[j] / w;
                }
                Assert.True(Math.Abs(score) < 1e-5 * Math.Max(scale, 1),
                    $"gene {i} coefficient {k}: score {score} against scale {scale}");
            }
            checkedGenes++;
        }

        Assert.True(checkedGenes > 50, $"only {checkedGenes} genes could be checked");
    }

    [Fact]
    public void GroupDesign_FoldChangeTracksTheRatioOfNormalizedMeans()
    {
        // Not an identity (see the score-equation test), but the two must stay close for well-expressed genes.
        var (counts, _, fit, res) = Sim.Value;
        int checkedGenes = 0;
        for (int i = 0; i < 200; i++)
        {
            double a = 0, b = 0;
            for (int j = 0; j < 3; j++) a += counts[i, j] / fit.SizeFactors[j];
            for (int j = 3; j < 6; j++) b += counts[i, j] / fit.SizeFactors[j];
            if (a / 3 < 20 || b / 3 < 20) continue;
            Assert.Equal(Math.Log2(b / a), res.Log2FoldChange[i], 0.1);
            checkedGenes++;
        }
        Assert.True(checkedGenes > 20, $"only {checkedGenes} genes were well enough expressed");
        Assert.All(fit.Converged, Assert.True);
    }

    [Fact]
    public void Results_WithoutIndependentFilteringIsPlainBh()
    {
        var (_, _, fit, _) = Sim.Value;
        var res = fit.Results(options: new DESeqOptions { IndependentFiltering = false });
        int n = res.Count;
        // Recompute BH by hand.
        var order = Enumerable.Range(0, n).Where(i => !double.IsNaN(res.PValue[i])).OrderBy(i => res.PValue[i]).ToArray();
        int m = order.Length;
        double running = 1;
        var expected = new double[n];
        for (int r = m - 1; r >= 0; r--)
        {
            running = Math.Min(running, res.PValue[order[r]] * m / (r + 1));
            expected[order[r]] = running;
        }
        foreach (int i in order) Assert.Equal(expected[i], res.PAdj[i], 1e-12);
    }

    [Fact]
    public void IndependentFiltering_DropsOnlyTheLowestMeans()
    {
        var (_, _, _, res) = Sim.Value;
        var kept = Enumerable.Range(0, res.Count).Where(i => !double.IsNaN(res.PAdj[i])).ToArray();
        var dropped = Enumerable.Range(0, res.Count).Where(i => double.IsNaN(res.PAdj[i])).ToArray();
        Assert.NotEmpty(kept);
        if (dropped.Length > 0)
            Assert.True(dropped.Max(i => res.BaseMean[i]) <= kept.Min(i => res.BaseMean[i]));
        // Among the kept genes padj is exactly BH over the kept set.
        var p = kept.Select(i => res.PValue[i]).ToArray();
        var order = Enumerable.Range(0, p.Length).OrderBy(i => p[i]).ToArray();
        double running = 1;
        for (int r = p.Length - 1; r >= 0; r--)
        {
            running = Math.Min(running, p[order[r]] * p.Length / (r + 1));
            Assert.Equal(running, res.PAdj[kept[order[r]]], 1e-12);
        }
    }

    [Fact]
    public void AllZeroGenes_AreExcluded()
    {
        var (counts, _, _) = TestData.TwoGroupCounts(200, 3, 10, 2, seed: 4);
        var row = counts.Row(5);
        row.Clear();
        var fit = DESeq.Fit(counts, DesignMatrix.FromAnnotations(counts, "condition", "A"));
        Assert.True(fit.AllZero[5]);
        var res = fit.Results();
        Assert.True(double.IsNaN(res.PValue[5]) && double.IsNaN(res.PAdj[5]));
        Assert.Equal(0, res.BaseMean[5]);
    }

    [Fact]
    public void Run_IsEquivalentToFitAndResults()
    {
        var (counts, _, fit, res) = Sim.Value;
        var run = DESeq.Run(counts, "condition", "A", "B");
        Assert.Equal(res.Log2FoldChange, run.Log2FoldChange);
        Assert.Equal(res.PValue, run.PValue);
        _ = fit;
    }

    [Fact]
    public void Vst_IsLog2LikeAndMatchesClosedForms()
    {
        var counts = Dataset.FromArray("c", new double[,] { { 0, 10, 1e6 } });
        double[] sf = [1, 1, 1];
        var mean = Transforms.Vst(counts, sf, 0.1, 0);
        var para = Transforms.Vst(counts, sf, 0.1, 1e-12);
        for (int j = 0; j < 3; j++) Assert.Equal(mean[0, j], para[0, j], 1e-6);
        Assert.Equal(Math.Log2(1e6), mean[0, 2], 1e-3);
        Assert.True(mean[0, 0] < mean[0, 1]);

        var (sim, _, _) = TestData.TwoGroupCounts(300, 3, 0, 0, seed: 8);
        var vst = Transforms.Vst(sim);
        Assert.Equal(sim.Rows, vst.Rows);
        Assert.All(vst.Values.ToArray(), v => Assert.True(double.IsFinite(v)));
        var rlog = Transforms.RlogApproximation(sim);
        Assert.All(rlog.Values.ToArray(), v => Assert.True(double.IsFinite(v)));
    }
}
