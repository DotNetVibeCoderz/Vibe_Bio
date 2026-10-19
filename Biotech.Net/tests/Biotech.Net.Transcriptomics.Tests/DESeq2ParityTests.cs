using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>
/// Parity of <see cref="DESeq"/> and <see cref="LfcShrink"/> with R (R 4.6.1, DESeq2 1.52.0, apeglm 1.34.0) on the
/// Bioconductor <c>pasilla</c> dataset (every 7th gene, 2,086 genes × 7 samples). Reference values come from
/// Fixtures/DESeq2/make_reference.R. Tolerances are relative unless stated; see each test for the reasoning.
/// </summary>
public class DESeq2ParityTests(ITestOutputHelper output)
{
    private static Dataset Counts => DESeq2Fixtures.Pasilla;

    private static DesignMatrix ConditionDesign => DesignMatrix.FromAnnotations(Counts, "condition", "untreated");
    private static DesignMatrix TypeConditionDesign => DesignMatrix.FromAnnotations(Counts, "condition", "untreated", "type");

    private static readonly Lazy<DESeqFit> FitA = new(() => DESeq.Fit(Counts, ConditionDesign));
    private static readonly Lazy<DESeqFit> FitB = new(() => DESeq.Fit(Counts, TypeConditionDesign));

    private static double[] Col(DESeqFit fit, int k)
    {
        int p = fit.Design.Coefficients, n = fit.Counts.Rows;
        var r = new double[n];
        for (int i = 0; i < n; i++) r[i] = fit.Coefficients[i * p + k];
        return r;
    }

    private static double[] SeCol(DESeqFit fit, int k)
    {
        int p = fit.Design.Coefficients, n = fit.Counts.Rows;
        var r = new double[n];
        for (int i = 0; i < n; i++) r[i] = fit.StandardErrors[i * p + k];
        return r;
    }

    private static double[] D(bool[] b) => b.Select(v => v ? 1.0 : 0.0).ToArray();
    private static double[] D(int[] b) => b.Select(v => (double)v).ToArray();

    private static double[] Mask(double[] v, bool[] allZero) => v.Select((x, i) => allZero[i] ? double.NaN : x).ToArray();

    // Tolerance policy.
    //  * Closed-form stages (size factors, base means/variances) are compared at 1e-12..1e-11 relative: they differ
    //    from R only by floating-point summation order.
    //  * Stages that run DESeq2's iterative optimisers (Armijo line search on log alpha, IRLS) are compared at 1e-6
    //    relative. The line search stops when the log posterior improves by < 1e-6. For high-count genes the log
    //    posterior is of order 1e5 and almost flat around its maximum (it changes by ~1e-10 across a 1e-4 relative change
    //    of alpha), and its accept/reject decisions compare values that differ from R's in the last bit, so R and this
    //    port can stop at different points of that flat top. Up to 1% of genes may therefore exceed 1e-6 if they stay
    //    within 1e-3 (iteration counts are reported and allowed to differ for the same reason).
    //  * Genes whose IRLS did not converge (R: betaIter = maxit) are refitted by R with optim(L-BFGS-B) and here with a
    //    projected Newton method. The optimum is the same but L-BFGS-B stops early on these flat likelihoods, so their
    //    coefficients are compared separately with an absolute tolerance of 1e-3 (log2 scale).
    private const double Tight = 1e-6, Loose = 1e-3;

    private static int Allow(int n) => Math.Max(2, n / 100);

    private static void Opt(ParityReport rep, string col, double[] expected, double[] actual, double abs = 0)
        => rep.Check(col, expected, actual, Tight, abs, Allow(expected.Length), Loose);

    private static double[] Where(double[] v, bool[] keep) => v.Select((x, i) => keep[i] ? x : double.NaN).ToArray();

    /// <summary>Strict check on IRLS-converged genes, loose absolute check on the optim-refitted ones.</summary>
    private static void OptSplit(ParityReport rep, string col, double[] expected, double[] actual, bool[] optimRows, double abs = 0,
        double? outlierAbs = null)
    {
        var irls = optimRows.Select(b => !b).ToArray();
        rep.Check(col, Where(expected, irls), Where(actual, irls), Tight, abs, Allow(expected.Length), Loose, outlierAbs);
        rep.Check(col + " (optim rows)", Where(expected, optimRows), Where(actual, optimRows), 1e-2, 1e-3);
    }

    /// <summary>Genes R refitted with optim because the IRLS reached maxit (any column of the scenario).</summary>
    private static bool[] OptimRows(Dictionary<string, double[]> t) => t["betaIter"].Select(v => v >= 100).ToArray();

    /// <summary>Checks every per-gene column DESeq2 stores in mcols(dds) plus the results() table.</summary>
    private static void CheckFit(ParityReport rep, string scenario, DESeqFit fit, DifferentialExpressionResult res, bool lrt = false,
        bool localTrend = false)
    {
        var t = DESeq2Fixtures.Table(scenario);
        var az = fit.AllZero;
        int n = az.Length;
        rep.Check("sizeFactors", Enumerable.Range(1, fit.SizeFactors.Length).Select(j => DESeq2Fixtures.Scalar(scenario, $"sizeFactor_{j}")).ToArray(),
            fit.SizeFactors, 1e-12);
        rep.Check("baseMean", t["baseMean"], fit.BaseMean, 1e-12);
        rep.Check("baseVar", t["baseVar"], fit.BaseVar, 1e-11);
        rep.CheckExact("allZero", t["allZero"], D(fit.AllZero));
        Opt(rep, "dispGeneEst", t["dispGeneEst"], fit.GeneDispersions);
        rep.CheckExact("dispGeneIter", t["dispGeneIter"], Mask(D(fit.GeneDispersionIterations), az), maxMismatch: n / 10);
        if (localTrend) return;
        Opt(rep, "dispFit", t["dispFit"], fit.TrendDispersions);
        Opt(rep, "dispMAP", t["dispMAP"], fit.MapDispersions);
        rep.CheckExact("dispIter", t["dispIter"], Mask(D(fit.DispersionIterations), az), maxMismatch: Allow(n));
        rep.CheckExact("dispOutlier", t["dispOutlier"], Mask(D(fit.DispersionOutliers), az));
        Opt(rep, "dispersion", t["dispersion"], fit.Dispersions);
        int p = fit.Design.Coefficients;
        var optim = OptimRows(t);
        string[] skip = ["baseVar", "allZero", "betaConv", "fullBetaConv", "reducedBetaConv", "betaIter", "deviance", "maxCooks",
            "replace", "LRTStatistic", "LRTPvalue"];
        var coefNames = t.Keys.Where(k => k != "baseMean" && !k.StartsWith("SE_") && !k.StartsWith("res_") && !k.StartsWith("disp")
            && !k.StartsWith("ape_") && !k.StartsWith("normal_") && !k.StartsWith("flip_") && !k.StartsWith("num_") && !skip.Contains(k))
            .ToArray();
        Assert.Equal(p, coefNames.Length);
        for (int k = 0; k < p; k++)
        {
            OptSplit(rep, $"coef[{coefNames[k]}]", t[coefNames[k]], Col(fit, k), optim, 1e-9);
            OptSplit(rep, $"SE[{coefNames[k]}]", t["SE_" + coefNames[k]], SeCol(fit, k), optim);
        }
        if (t.ContainsKey("betaConv")) rep.CheckExact("betaConv", t["betaConv"], Mask(D(fit.Converged), az));
        rep.CheckExact("betaIter", t["betaIter"], Mask(D(fit.BetaIterations), az), maxMismatch: Allow(n));
        OptSplit(rep, "deviance", t["deviance"], fit.Deviance, optim);
        Opt(rep, "maxCooks", t["maxCooks"], fit.MaxCooks);
        if (lrt)
        {
            OptSplit(rep, "LRTStatistic", t["LRTStatistic"], fit.LrtStatistic!, optim, 1e-8);
            Opt(rep, "LRTPvalue", t["LRTPvalue"], fit.LrtPValue!, 1e-12);
        }
        rep.Check("res.baseMean", t["res_baseMean"], res.BaseMean, 1e-12);
        OptSplit(rep, "res.log2FoldChange", t["res_log2FoldChange"], res.Log2FoldChange, optim, 1e-9);
        OptSplit(rep, "res.lfcSE", t["res_lfcSE"], res.LfcSE, optim);
        OptSplit(rep, "res.stat", t["res_stat"], res.Stat, optim, 1e-8);
        OptSplit(rep, "res.pvalue", t["res_pvalue"], res.PValue, optim, 1e-12);
        OptSplit(rep, "res.padj", t["res_padj"], res.PAdj, optim, 1e-12);
    }

    [Fact]
    public void A_Condition_Wald()
    {
        var fit = FitA.Value;
        var rep = new ParityReport("A_condition");
        CheckFit(rep, "A_condition", fit, fit.Results());
        rep.Check("varLogDispEsts", [DESeq2Fixtures.Scalar("A_condition", "varLogDispEsts")], [fit.VarLogDispersionEstimates], 1e-8);
        rep.Check("dispPriorVar", [DESeq2Fixtures.Scalar("A_condition", "dispPriorVar")], [fit.DispersionPriorVariance], 1e-8);
        rep.Check("asymptDisp", [DESeq2Fixtures.Scalar("A_condition", "asymptDisp")], [fit.Trend.AsymptoticDispersion], 1e-8);
        rep.Check("extraPois", [DESeq2Fixtures.Scalar("A_condition", "extraPois")], [fit.Trend.ExtraPoisson], 1e-8);
        var nofilter = DESeq2Fixtures.Table("A_condition_nofilter");
        var r0 = fit.Results(options: new DESeqOptions { IndependentFiltering = false, CooksFiltering = false });
        Opt(rep, "nofilter.pvalue", nofilter["pvalue"], r0.PValue, 1e-12);
        Opt(rep, "nofilter.padj", nofilter["padj"], r0.PAdj, 1e-12);
        rep.AssertOk(output);
    }

    // apeglm: the MAP is found by a port of the same L-BFGS (LBFGSpp, eps_g = eps_f = 1e-8 on an objective scaled by the
    // null log-likelihood), which reproduces apeglm's C++ iterates to ~1e-13 for identical inputs; differences come from
    // the ~1e-6 differences in the dispersions and MLEs it starts from. The shrunken LFCs are often ~0, so they get an
    // absolute tolerance of 1e-5 (log2 scale). The posterior SD comes from optimHess at the MAP and inherits that accuracy.
    // For low-count genes at the maximum dispersion the apeglm posterior is nearly flat and L-BFGS stops wherever its
    // relative-change test first fires: R's own two starts (apeglm runs two and accepts them if they agree within 0.01)
    // typically end ~1e-4 apart there, and last-bit differences move this port's stopping point by the same amount.
    // Up to 1% of genes may therefore differ by up to 1e-3 (log2 scale) in the shrunken fold change.
    private const double ApeOutlierAbs = 1e-3;

    private static void CheckShrink(ParityReport rep, string scenario, DESeqFit fit, int coef, bool svalue)
    {
        var t = DESeq2Fixtures.Table(scenario);
        var ape = LfcShrink.Apeglm(fit, coef);
        rep.Check("apeglm.priorScale", [DESeq2Fixtures.Scalar(scenario, "ape_prior_scale")], [ape.PriorScale], 1e-9);
        OptSplit(rep, "apeglm.log2FoldChange", t["ape_lfc"], ape.Table.Log2FoldChange, OptimRows(t), 1e-5, ApeOutlierAbs);
        OptSplit(rep, "apeglm.lfcSE", t["ape_se"], ape.Table.LfcSE, OptimRows(t));
        if (svalue) Opt(rep, "apeglm.svalue", t["ape_svalue"], ape.SValue!, 1e-6);
        rep.CheckExact("apeglm.conv", t["ape_conv"], ape.Convergence!.Select((c, i) => fit.AllZero[i] ? double.NaN : c).ToArray(), Allow(t["ape_conv"].Length));
        var nor = LfcShrink.Normal(fit, coef);
        OptSplit(rep, "normal.log2FoldChange", t["normal_lfc"], nor.Table.Log2FoldChange, OptimRows(t), 1e-9);
        OptSplit(rep, "normal.lfcSE", t["normal_se"], nor.Table.LfcSE, OptimRows(t));
    }

    [Fact]
    public void A_Condition_LfcShrink()
    {
        var fit = FitA.Value;
        var rep = new ParityReport("A_condition");
        CheckShrink(rep, "A_condition", fit, 1, svalue: true);
        rep.Check("normal.betaPriorVar", [DESeq2Fixtures.Scalar("A_condition", "normal_betaPriorVar")], [LfcShrink.Normal(fit, 1).PriorVariance], 1e-9);
        rep.AssertOk(output);
    }

    [Fact]
    public void B_TypeCondition_LfcShrink()
    {
        var fit = FitB.Value;
        var rep = new ParityReport("B_type_condition");
        CheckShrink(rep, "B_type_condition", fit, 2, svalue: false);
        rep.Check("normal.betaPriorVar", [DESeq2Fixtures.Scalar("B_type_condition", "normal_betaPriorVar_condition")], [LfcShrink.Normal(fit, 2).PriorVariance], 1e-9);
        rep.AssertOk(output);
    }

    [Fact]
    public void I_Airway_SimulatedPriorVariance()
    {
        // airway: ~ cell + dex leaves m - p = 3 residual degrees of freedom, so DESeq2 estimates the dispersion prior
        // variance by its seeded simulation (set.seed(2), rchisq/rnorm, hist, loess). The argmin lies on a 1,000-point
        // grid, so the prior variance must match exactly.
        var airway = DESeq2Fixtures.Airway;
        var fit = DESeq.Fit(airway, DesignMatrix.FromAnnotations(airway, "dex", "untrt", "cell"));
        var rep = new ParityReport("I_airway");
        rep.Check("dispPriorVar", [DESeq2Fixtures.Scalar("I_airway", "dispPriorVar")], [fit.DispersionPriorVariance], 1e-12);
        CheckFit(rep, "I_airway", fit, fit.Results());
        var t = DESeq2Fixtures.Table("I_airway");
        var ape = LfcShrink.Apeglm(fit);
        OptSplit(rep, "apeglm.log2FoldChange", t["ape_lfc"], ape.Table.Log2FoldChange, OptimRows(t), 1e-5, ApeOutlierAbs);
        OptSplit(rep, "apeglm.lfcSE", t["ape_se"], ape.Table.LfcSE, OptimRows(t));
        rep.AssertOk(output);
    }

    [Fact]
    public void RRandom_ReproducesRsStream()
    {
        // set.seed(2); runif(3); rnorm(3); rchisq(3, 1); rchisq(3, 2); rchisq(3, 3); rnorm(2, 0, 0); rgamma(3, 0.5); rnorm(2)
        var expected = File.ReadAllLines(Path.Combine(DESeq2Fixtures.Directory, "rng.txt"))
            .Select(l => double.Parse(l.Trim(), System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        var r = new Internal.RRandom(2);
        var actual = new List<double>();
        for (int i = 0; i < 3; i++) actual.Add(r.Unif());
        for (int i = 0; i < 3; i++) actual.Add(r.Rnorm(0, 1));
        foreach (double df in new[] { 1.0, 2, 3 }) for (int i = 0; i < 3; i++) actual.Add(r.Rchisq(df));
        for (int i = 0; i < 2; i++) actual.Add(r.Rnorm(0, 0));
        for (int i = 0; i < 3; i++) actual.Add(r.Rgamma(0.5, 1));
        for (int i = 0; i < 2; i++) actual.Add(r.Rnorm(0, 1));
        Assert.Equal(expected.Length, actual.Count);
        for (int i = 0; i < expected.Length; i++) Assert.Equal(expected[i], actual[i], 1e-14 * Math.Max(1, Math.Abs(expected[i])));
    }

    [Fact]
    public void RLoess_MatchesRsInterpolatedSurface()
    {
        var expected = File.ReadAllLines(Path.Combine(DESeq2Fixtures.Directory, "loess.txt"))
            .Select(l => double.Parse(l.Trim(), System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        var x = Internal.DESeqEngine.RSeq(0, 8, 200);
        var y = x.Select(v => Math.Sin(v) + 0.1 * v * v + Math.Cos(7 * v) / 5).ToArray();
        var p = Internal.RLoess.PredictInterpolate(x, y, 0.2, Internal.DESeqEngine.RSeq(0, 8, 1000));
        var rep = new ParityReport("loess");
        rep.Check("predict", expected, p, 1e-10, 1e-12);
        rep.AssertOk(output);
    }

    [Fact]
    public void B_TypeCondition_WaldAndContrasts()
    {
        var fit = FitB.Value;
        var rep = new ParityReport("B_type_condition");
        CheckFit(rep, "B_type_condition", fit, fit.Results());
        var t = DESeq2Fixtures.Table("B_type_condition");
        var flip = fit.Results("condition", "untreated", "treated");
        OptSplit(rep, "flip.log2FoldChange", t["flip_lfc"], flip.Log2FoldChange, OptimRows(t), 1e-9);
        OptSplit(rep, "flip.stat", t["flip_stat"], flip.Stat, OptimRows(t), 1e-8);
        OptSplit(rep, "flip.padj", t["flip_padj"], flip.PAdj, OptimRows(t), 1e-12);
        var num = fit.Results([0, -1, 1]);
        OptSplit(rep, "num.log2FoldChange", t["num_lfc"], num.Log2FoldChange, OptimRows(t), 1e-9);
        OptSplit(rep, "num.lfcSE", t["num_se"], num.LfcSE, OptimRows(t));
        OptSplit(rep, "num.pvalue", t["num_pvalue"], num.PValue, OptimRows(t), 1e-12);
        OptSplit(rep, "num.padj", t["num_padj"], num.PAdj, OptimRows(t), 1e-12);
        rep.Check("dispPriorVar", [DESeq2Fixtures.Scalar("B_type_condition", "dispPriorVar")], [fit.DispersionPriorVariance], 1e-8);
        rep.AssertOk(output);
    }

    [Fact]
    public void C_Lrt_AgainstType()
    {
        var reduced = DesignMatrix.FromAnnotations(Counts, "type", "paired-end");
        var fit = DESeq.FitLrt(Counts, TypeConditionDesign, reduced);
        var rep = new ParityReport("C_lrt_type");
        CheckFit(rep, "C_lrt_type", fit, fit.Results(), lrt: true);
        rep.AssertOk(output);
    }

    [Fact]
    public void D_Lrt_AgainstIntercept()
    {
        var fit = DESeq.FitLrt(Counts, ConditionDesign, DesignMatrix.Intercept(Counts.Columns));
        var rep = new ParityReport("D_lrt_intercept");
        CheckFit(rep, "D_lrt_intercept", fit, fit.Results(), lrt: true);
        rep.AssertOk(output);
    }

    [Fact]
    public void E_PosCountsSizeFactors()
    {
        var fit = DESeq.Fit(Counts, ConditionDesign, new DESeqOptions { SizeFactorType = SizeFactorType.PosCounts });
        var rep = new ParityReport("E_poscounts");
        CheckFit(rep, "E_poscounts", fit, fit.Results());
        rep.AssertOk(output);
    }

    [Fact]
    public void F_MeanTrend()
    {
        var fit = DESeq.Fit(Counts, ConditionDesign, new DESeqOptions { FitType = DispersionFitType.Mean });
        var rep = new ParityReport("F_mean");
        CheckFit(rep, "F_mean", fit, fit.Results());
        rep.Check("meanDisp", [DESeq2Fixtures.Scalar("F_mean", "meanDisp")], [fit.Trend.AsymptoticDispersion], 1e-9);
        rep.AssertOk(output);
    }

    [Fact]
    public void G_LocalTrend_IsApproximate()
    {
        // locfit's adaptive-tree interpolation is not reproduced: the local trend is a direct local-quadratic fit, so only
        // the gene-wise stage is exact and the trend is checked loosely (median and 95th percentile of relative error).
        var fit = DESeq.Fit(Counts, ConditionDesign, new DESeqOptions { FitType = DispersionFitType.Local });
        Assert.Equal(DispersionFitType.Local, fit.TrendType);
        var rep = new ParityReport("G_local");
        CheckFit(rep, "G_local", fit, fit.Results(), localTrend: true);
        var t = DESeq2Fixtures.Table("G_local");
        var res = fit.Results();
        var rel = Enumerable.Range(0, fit.TrendDispersions.Length).Where(i => !fit.AllZero[i])
            .Select(i => Math.Abs(fit.TrendDispersions[i] / t["dispFit"][i] - 1)).Order().ToArray();
        double median = rel[rel.Length / 2], q95 = rel[(int)(rel.Length * 0.95)];
        output.WriteLine($"G_local dispFit relative error: median {median:E2}, 95% {q95:E2}, max {rel[^1]:E2}");
        var lfcErr = Enumerable.Range(0, fit.Counts.Rows).Where(i => !fit.AllZero[i])
            .Select(i => Math.Abs(res.Log2FoldChange[i] - t["res_log2FoldChange"][i])).Max();
        output.WriteLine($"G_local max |log2FC difference| {lfcErr:E2}");
        Assert.True(median < 0.05 && q95 < 0.25, $"local trend: median {median}, q95 {q95}");
        rep.AssertOk(output);
    }

    [Fact]
    public void H_CooksReplacement()
    {
        var fit = DESeq.Fit(Counts, ConditionDesign, new DESeqOptions { MinReplicatesForReplace = 3 });
        var rep = new ParityReport("H_replace");
        CheckFit(rep, "H_replace", fit, fit.Results());
        rep.CheckExact("replace", DESeq2Fixtures.Table("H_replace")["replace"], Mask(D(fit.Replaced), fit.AllZero.Select((z, i) => z && !fit.Replaced[i]).ToArray()));
        rep.AssertOk(output);
    }
}
