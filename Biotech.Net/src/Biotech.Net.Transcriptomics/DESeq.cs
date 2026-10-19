using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>How the mean–dispersion trend is fitted (DESeq2 <c>fitType</c>).</summary>
public enum DispersionFitType
{
    /// <summary>DESeq2's parametric trend α(μ) = asymptDisp + extraPois / μ, fitted by a gamma-family GLM.
    /// When the fit fails DESeq2 (and this port) substitute a <see cref="Local"/> fit.</summary>
    Parametric,

    /// <summary>A constant trend: the mean of the gene-wise estimates above 10·minDisp, trimmed by 0.1% (DESeq2 "mean").</summary>
    Mean,

    /// <summary>A local regression of log dispersion on log mean (DESeq2 "local"; approximates locfit).</summary>
    Local,
}

/// <summary>Size-factor estimator (DESeq2 <c>sfType</c>).</summary>
public enum SizeFactorType
{
    /// <summary>Median-of-ratios over genes without zeros (DESeq2 "ratio", the default).</summary>
    Ratio,

    /// <summary>Median-of-ratios with geometric means over positive counts only, scaled to geometric mean 1
    /// (DESeq2 "poscounts"; use when every gene has a zero, e.g. sparse data).</summary>
    PosCounts,
}

/// <summary>Options for <see cref="DESeq"/> (defaults are DESeq2's).</summary>
public sealed class DESeqOptions
{
    /// <summary>Lower bound for dispersions (DESeq2 <c>minDisp</c>, default 1e-8).</summary>
    public double MinDispersion { get; init; } = 1e-8;

    /// <summary>Trend type (default parametric, falling back to local if the parametric fit fails).</summary>
    public DispersionFitType FitType { get; init; } = DispersionFitType.Parametric;

    /// <summary>Size-factor estimator (default <see cref="SizeFactorType.Ratio"/>). Ignored when <see cref="SizeFactors"/> is given.</summary>
    public SizeFactorType SizeFactorType { get; init; } = SizeFactorType.Ratio;

    /// <summary>Shrink gene-wise dispersions towards the trend (maximum a posteriori). Default true; false uses the
    /// gene-wise estimates directly (not a DESeq2 option, used for blind transformations).</summary>
    public bool ShrinkDispersions { get; init; } = true;

    /// <summary>Gene-wise estimates more than this many SDs (of the log residuals) above the trend are kept unshrunk (default 2).</summary>
    public double OutlierSD { get; init; } = 2;

    /// <summary>Maximum iterations for the dispersion line search and the IRLS GLM fit (DESeq2 <c>maxit</c>, default 100).</summary>
    public int MaxIterations { get; init; } = 100;

    /// <summary>Ridge penalty on log2-scale coefficients (DESeq2 uses 1e-6 when no LFC prior is used).</summary>
    public double Lambda { get; init; } = 1e-6;

    /// <summary>IRLS convergence tolerance on the relative deviance change (DESeq2 <c>betaTol</c>, default 1e-8).</summary>
    public double BetaTolerance { get; init; } = 1e-8;

    /// <summary>Lower bound for fitted means inside the IRLS and dispersion fits (DESeq2 <c>minmu</c>, default 0.5).</summary>
    public double MinMu { get; init; } = 0.5;

    /// <summary>Refit genes whose IRLS did not converge with a bounded optimiser (DESeq2 <c>useOptim</c>, default true).</summary>
    public bool UseOptim { get; init; } = true;

    /// <summary>Optimise the BH threshold by filtering on the mean of normalised counts (default true).</summary>
    public bool IndependentFiltering { get; init; } = true;

    /// <summary>Significance level the independent filter is tuned for (DESeq2 <c>alpha</c>, default 0.1).</summary>
    public double Alpha { get; init; } = 0.1;

    /// <summary>Precomputed size factors; estimated from the counts when null.</summary>
    public double[]? SizeFactors { get; init; }

    /// <summary>Set p-values of genes whose maximum Cook's distance exceeds <see cref="CooksCutoff"/> to NaN (default true).</summary>
    public bool CooksFiltering { get; init; } = true;

    /// <summary>Cook's distance cutoff; null uses DESeq2's default, the 0.99 quantile of F(p, m − p).</summary>
    public double? CooksCutoff { get; init; }

    /// <summary>Minimum number of replicates in a design cell for its outlier counts to be replaced by the trimmed mean and
    /// the gene refitted (DESeq2 <c>minReplicatesForReplace</c>, default 7; <see cref="int.MaxValue"/> disables).</summary>
    public int MinReplicatesForReplace { get; init; } = 7;

    /// <summary>Fixed prior variance of log dispersions (DESeq2 <c>dispPriorVar</c>); estimated when null.</summary>
    public double? DispersionPriorVariance { get; init; }

    /// <summary>Reduced design for a likelihood-ratio test (DESeq2 <c>test = "LRT"</c>, <c>reduced</c>); null runs Wald tests.</summary>
    public DesignMatrix? ReducedDesign { get; init; }
}

/// <summary>
/// The fitted negative-binomial model of every gene: size factors, dispersions at each stage, the GLM coefficients
/// with their standard errors (log2 scale), Cook's distances and, for a likelihood-ratio test, the LRT statistics.
/// Produce result tables with the <c>Results</c> overloads and shrunken fold changes with <see cref="LfcShrink"/>.
/// </summary>
public sealed class DESeqFit
{
    internal DESeqFit(Dataset counts, DesignMatrix design) { Counts = counts; Design = design; }

    /// <summary>The count matrix that was fitted (original counts, before any outlier replacement).</summary>
    public Dataset Counts { get; }

    /// <summary>The (full) model matrix.</summary>
    public DesignMatrix Design { get; }

    /// <summary>The reduced model matrix of a likelihood-ratio test (null for Wald tests).</summary>
    public DesignMatrix? ReducedDesign { get; internal set; }

    /// <summary>"Wald" or "LRT".</summary>
    public string Test => ReducedDesign is null ? "Wald" : "LRT";

    /// <summary>Per-sample size factors.</summary>
    public double[] SizeFactors { get; internal set; } = [];

    /// <summary>Mean of normalised counts per gene (after outlier replacement, as DESeq2 reports it).</summary>
    public double[] BaseMean { get; internal set; } = [];

    /// <summary>Variance of normalised counts per gene.</summary>
    public double[] BaseVar { get; internal set; } = [];

    /// <summary>Gene-wise maximum Cox–Reid adjusted profile likelihood dispersions (<c>dispGeneEst</c>).</summary>
    public double[] GeneDispersions { get; internal set; } = [];

    /// <summary>Line-search iterations of the gene-wise estimate (<c>dispGeneIter</c>).</summary>
    public int[] GeneDispersionIterations { get; internal set; } = [];

    /// <summary>Fitted trend at each gene's mean (<c>dispFit</c>).</summary>
    public double[] TrendDispersions { get; internal set; } = [];

    /// <summary>Maximum a posteriori dispersions (<c>dispMAP</c>).</summary>
    public double[] MapDispersions { get; internal set; } = [];

    /// <summary>Line-search iterations of the MAP estimate (<c>dispIter</c>).</summary>
    public int[] DispersionIterations { get; internal set; } = [];

    /// <summary>Final dispersions used for testing (<c>dispersion</c>): MAP estimates, or gene-wise for outliers.</summary>
    public double[] Dispersions { get; internal set; } = [];

    /// <summary>True where the gene-wise estimate was flagged as a dispersion outlier and kept.</summary>
    public bool[] DispersionOutliers { get; internal set; } = [];

    /// <summary>Parametric trend coefficients (asymptDisp, extraPois); for a mean trend extraPois is 0, for a local trend both are NaN.</summary>
    public (double AsymptoticDispersion, double ExtraPoisson) Trend { get; internal set; }

    /// <summary>Which trend was actually used (a failed parametric fit falls back to local).</summary>
    public DispersionFitType TrendType { get; internal set; }

    /// <summary>The fitted mean–dispersion trend as a function of the mean of normalised counts.</summary>
    public Func<double, double> TrendFunction { get; internal set; } = _ => double.NaN;

    /// <summary>Squared MAD of log(gene-wise / trend) dispersions (<c>varLogDispEsts</c>).</summary>
    public double VarLogDispersionEstimates { get; internal set; }

    /// <summary>Prior variance of log dispersions around the trend (<c>dispPriorVar</c>).</summary>
    public double DispersionPriorVariance { get; internal set; }

    /// <summary>Coefficients on the log2 scale, genes × coefficients (row-major).</summary>
    public double[] Coefficients { get; internal set; } = [];

    /// <summary>Standard errors of the log2 coefficients, genes × coefficients (row-major).</summary>
    public double[] StandardErrors { get; internal set; } = [];

    /// <summary>Whether the GLM fit converged for each gene (<c>betaConv</c>).</summary>
    public bool[] Converged { get; internal set; } = [];

    /// <summary>IRLS iterations per gene (<c>betaIter</c>).</summary>
    public int[] BetaIterations { get; internal set; } = [];

    /// <summary>−2 × log-likelihood of the (full) model (<c>deviance</c>).</summary>
    public double[] Deviance { get; internal set; } = [];

    /// <summary>True for genes with only zero counts (excluded from testing).</summary>
    public bool[] AllZero { get; internal set; } = [];

    /// <summary>Fitted means of the initial fit, genes × samples (DESeq2's <c>mu</c> assay).</summary>
    public double[] Mu { get; internal set; } = [];

    /// <summary>Hat-matrix diagonals, genes × samples (DESeq2's <c>H</c> assay).</summary>
    public double[] HatDiagonals { get; internal set; } = [];

    /// <summary>Cook's distances of the initial fit, genes × samples (DESeq2's <c>cooks</c> assay).</summary>
    public double[] CooksDistances { get; internal set; } = [];

    /// <summary>Maximum Cook's distance per gene over samples in cells with ≥ 3 replicates (<c>maxCooks</c>; NaN when undefined).</summary>
    public double[] MaxCooks { get; internal set; } = [];

    /// <summary>True for genes whose outlier counts were replaced and which were refitted (<c>replace</c>).</summary>
    public bool[] Replaced { get; internal set; } = [];

    /// <summary>Counts after outlier replacement (null when nothing was replaced).</summary>
    public Dataset? ReplacementCounts { get; internal set; }

    /// <summary>Likelihood-ratio statistic 2(ℓ_full − ℓ_reduced) per gene (LRT only).</summary>
    public double[]? LrtStatistic { get; internal set; }

    /// <summary>Likelihood-ratio test p-values (LRT only).</summary>
    public double[]? LrtPValue { get; internal set; }

    /// <summary>Number of genes refitted by the bounded optimiser because the IRLS did not converge.</summary>
    public int OptimRefits { get; internal set; }

    internal bool[] ReplaceableSamples { get; set; } = [];

    /// <summary>
    /// Results for the coefficient at <paramref name="coefficient"/> (default: the last one): Wald test (or the LRT
    /// statistics for a likelihood-ratio fit), Cook's distance filtering, and BH adjustment with independent filtering.
    /// </summary>
    public DifferentialExpressionResult Results(int coefficient = -1, DESeqOptions? options = null)
    {
        options ??= new DESeqOptions();
        int p = Design.Coefficients;
        if (coefficient < 0) coefficient = p - 1;
        if (coefficient >= p) throw new ArgumentOutOfRangeException(nameof(coefficient));
        int n = Counts.Rows;
        var lfc = new double[n]; var se = new double[n]; var stat = new double[n]; var pv = new double[n];
        for (int i = 0; i < n; i++)
        {
            lfc[i] = Coefficients[i * p + coefficient];
            se[i] = StandardErrors[i * p + coefficient];
            stat[i] = lfc[i] / se[i];
            pv[i] = RMath.TwoSidedNormal(stat[i]);
        }
        return Finish(Design.CoefficientNames[coefficient], lfc, se, stat, pv, null, options);
    }

    /// <summary>Results for a named coefficient.</summary>
    public DifferentialExpressionResult Results(string coefficient, DESeqOptions? options = null)
        => Results(Design.IndexOf(coefficient), options);

    /// <summary>
    /// Results for the contrast <paramref name="numerator"/> vs <paramref name="denominator"/> of the design factor
    /// <paramref name="factor"/> (DESeq2 <c>contrast = c(factor, numerator, denominator)</c>). Genes with zero counts in
    /// every sample of both levels get log2FC 0, stat 0 and p = 1, as in DESeq2.
    /// </summary>
    public DifferentialExpressionResult Results(string factor, string numerator, string denominator, DESeqOptions? options = null)
    {
        options ??= new DESeqOptions();
        var f = Design.Factors.FirstOrDefault(d => d.Name == factor)
            ?? throw new KeyNotFoundException($"'{factor}' is not a factor of the design (build it with DesignMatrix.FromAnnotations).");
        int li = IndexOfLevel(f, numerator), lj = IndexOfLevel(f, denominator);
        int ci = f.Columns[li], cj = f.Columns[lj];
        string name = $"{factor}_{numerator}_vs_{denominator}";
        var inLevels = new bool[Design.Samples];
        for (int s = 0; s < Design.Samples; s++) { int lv = LevelOfSample(f, s); inLevels[s] = lv == li || lv == lj; }
        var allZeroInContrast = ContrastAllZero(inLevels);
        int n = Counts.Rows, p = Design.Coefficients;
        DifferentialExpressionResult res;
        if (cj < 0 || ci < 0)
        {
            int col = cj < 0 ? ci : cj;
            double sign = cj < 0 ? 1 : -1;
            var lfc = new double[n]; var se = new double[n]; var stat = new double[n]; var pv = new double[n];
            for (int i = 0; i < n; i++)
            {
                lfc[i] = sign * Coefficients[i * p + col];
                se[i] = StandardErrors[i * p + col];
                double z = Coefficients[i * p + col] / se[i];
                stat[i] = ReducedDesign is null ? sign * z : z;
                pv[i] = RMath.TwoSidedNormal(z);
            }
            res = Finish(name, lfc, se, stat, pv, allZeroInContrast, options);
        }
        else
        {
            var c = new double[p];
            c[ci] = 1; c[cj] = -1;
            res = ContrastResults(c, name, allZeroInContrast, options);
        }
        return res;
    }

    /// <summary>
    /// Results for a numeric contrast vector over the coefficients (DESeq2 <c>contrast = c(...)</c>): the log2 fold
    /// change c'β with standard error sqrt(c'Σc) from the GLM covariance at the fitted coefficients.
    /// </summary>
    public DifferentialExpressionResult Results(double[] contrast, DESeqOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(contrast);
        options ??= new DESeqOptions();
        if (contrast.Length != Design.Coefficients) throw new ArgumentException($"The contrast needs {Design.Coefficients} values.");
        bool[]? zero = null;
        bool anyPos = contrast.Any(v => v > 0), anyNeg = contrast.Any(v => v < 0);
        if (anyPos && anyNeg)
        {
            var which = new bool[Design.Samples];
            for (int s = 0; s < Design.Samples; s++)
            {
                double v = 0;
                for (int k = 0; k < contrast.Length; k++) v += Design[s, k] * (contrast[k] == 0 ? 0 : 1);
                which[s] = v != 0;
            }
            zero = ContrastAllZero(which);
        }
        string name = string.Join(",", contrast.Select(v => (v > 0 ? "+" : "") + v.ToString(System.Globalization.CultureInfo.InvariantCulture)));
        return ContrastResults(contrast, name, zero, options);
    }

    private DifferentialExpressionResult ContrastResults(double[] c, string name, bool[]? allZeroInContrast, DESeqOptions options)
    {
        int n = Counts.Rows, m = Design.Samples, p = Design.Coefficients;
        var lfc = new double[n]; var se = new double[n]; var stat = new double[n]; var pv = new double[n];
        double ln2 = Math.Log(2);
        var lambdaNat = new double[p];
        Array.Fill(lambdaNat, 1 / (ln2 * ln2 * 1e6));
        var y = ReplacementCounts ?? Counts;
        Parallel.For(0, n, i =>
        {
            if (AllZero[i] || double.IsNaN(Coefficients[i * p]))
            {
                lfc[i] = se[i] = stat[i] = pv[i] = double.NaN;
                return;
            }
            var beta = new double[p];
            for (int k = 0; k < p; k++) beta[k] = ln2 * Coefficients[i * p + k];
            var fb = DESeqEngine.FitBeta(y.Row(i), SizeFactors, Design, Dispersions[i], beta, lambdaNat, 0, 1e-8, MinMu);
            double num = 0, den = 0;
            for (int a = 0; a < p; a++)
            {
                num += c[a] * beta[a];
                for (int b = 0; b < p; b++) den += c[a] * fb.Sigma[a * p + b] * c[b];
            }
            const double log2e = 1.4426950408889634;
            lfc[i] = log2e * num;
            se[i] = log2e * Math.Sqrt(den);
            stat[i] = lfc[i] / se[i];
            pv[i] = RMath.TwoSidedNormal(stat[i]);
        });
        return Finish(name, lfc, se, stat, pv, allZeroInContrast, options);
    }

    internal double MinMu { get; set; } = 0.5;

    private bool[] ContrastAllZero(bool[] samples)
    {
        int n = Counts.Rows;
        var res = new bool[n];
        for (int i = 0; i < n; i++)
        {
            if (AllZero[i]) continue;
            var row = Counts.Row(i);
            bool zero = true;
            for (int s = 0; s < samples.Length && zero; s++) if (samples[s] && row[s] != 0) zero = false;
            res[i] = zero;
        }
        return res;
    }

    private static int IndexOfLevel(DesignFactor f, string level)
    {
        for (int i = 0; i < f.Levels.Count; i++) if (f.Levels[i] == level) return i;
        throw new KeyNotFoundException($"'{level}' is not a level of '{f.Name}'. Levels: {string.Join(", ", f.Levels)}.");
    }

    private int LevelOfSample(DesignFactor f, int s)
    {
        for (int l = 1; l < f.Levels.Count; l++) if (Design[s, f.Columns[l]] != 0) return l;
        return 0;
    }

    /// <summary>The shared tail of DESeq2's <c>results()</c>: LRT statistics, Cook's filter, replaced-to-zero genes, padj.</summary>
    private DifferentialExpressionResult Finish(string name, double[] lfc, double[] se, double[] stat, double[] pv,
        bool[]? contrastAllZero, DESeqOptions options)
    {
        int n = Counts.Rows;
        if (contrastAllZero is not null)
            for (int i = 0; i < n; i++)
                if (contrastAllZero[i] && !AllZero[i]) { lfc[i] = 0; stat[i] = 0; pv[i] = 1; }
        if (LrtStatistic is not null)
            for (int i = 0; i < n; i++) { stat[i] = LrtStatistic[i]; pv[i] = LrtPValue![i]; }
        int m = Design.Samples, p = Design.Coefficients;
        for (int i = 0; i < n; i++)
            if (AllZero[i] && !Replaced[i]) { lfc[i] = se[i] = stat[i] = pv[i] = double.NaN; }
        if (options.CooksFiltering && m > p)
        {
            double cutoff = options.CooksCutoff ?? RMath.Qf(0.99, p, m - p);
            var outlier = new bool[n];
            for (int i = 0; i < n; i++) outlier[i] = MaxCooks[i] > cutoff;
            // DESeq2 heuristic for a single two-level factor: keep genes where ≥ 3 samples exceed the outlier count.
            bool twoLevel = Design.Factors.Count == 1 && Design.Factors[0].Levels.Count == 2 && p == 2;
            if (twoLevel)
                for (int i = 0; i < n; i++)
                {
                    if (!outlier[i]) continue;
                    var row = Counts.Row(i);
                    int jmax = 0;
                    for (int j = 1; j < m; j++) if (CooksDistances[i * m + j] > CooksDistances[i * m + jmax]) jmax = j;
                    double outCount = row[jmax];
                    int above = 0;
                    for (int j = 0; j < m; j++) if (row[j] > outCount) above++;
                    if (above >= 3) outlier[i] = false;
                }
            for (int i = 0; i < n; i++) if (outlier[i]) pv[i] = double.NaN;
        }
        for (int i = 0; i < n; i++)
            if (Replaced[i] && BaseMean[i] == 0) { lfc[i] = 0; se[i] = 0; stat[i] = 0; pv[i] = 1; }
        var padj = options.IndependentFiltering ? DESeq.IndependentFilter(BaseMean, pv, options.Alpha) : DESeq.Bh(pv);
        return new DifferentialExpressionResult(ReducedDesign is null ? "DESeq2-Wald" : "DESeq2-LRT", name, Counts.RowNames,
            (double[])BaseMean.Clone(), lfc, se, stat, pv, padj);
    }

    /// <summary>
    /// DESeq2's variance-stabilising transformation using the fitted parametric trend (the closed form of
    /// <c>getVarianceStabilizedData</c> with <c>fitType = "parametric"</c>, applied to normalised counts). With a mean
    /// trend the closed form for constant dispersion is used. Output is on a log2-like scale.
    /// </summary>
    public Dataset VarianceStabilizingTransform()
        => Transforms.Vst(Counts, SizeFactors, Trend.AsymptoticDispersion, Trend.ExtraPoisson);
}

/// <summary>
/// DESeq2 differential expression (Love, Huber &amp; Anders 2014), ported from DESeq2 1.52 step for step:
/// median-of-ratios (or poscounts) size factors; gene-wise dispersions maximising the Cox–Reid adjusted profile
/// likelihood by DESeq2's Armijo line search (grid fallback); a parametric, local or mean trend; the prior variance of
/// log dispersions; MAP shrinkage with dispersion-outlier handling; ridge-penalised IRLS negative-binomial GLMs; Wald or
/// likelihood-ratio tests; Cook's distance flagging and replacement of outliers with refitting; and BH adjustment
/// with genefilter-style independent filtering on the mean of normalised counts.
/// </summary>
/// <remarks>
/// <para>Parity with R (DESeq2 1.52.0, R 4.6.1, Bioconductor 3.23) is checked on the public <c>pasilla</c> and
/// <c>airway</c> datasets by the test suite: closed-form stages agree to ~1e-12 relative, iterative stages (line search,
/// IRLS) to 1e-6 for at least 99% of genes and to 1e-3 for the rest, because the line search's stopping point on flat
/// posteriors is sensitive to last-bit arithmetic differences. Known remaining differences:</para>
/// <list type="bullet">
/// <item>The <see cref="DispersionFitType.Local"/> trend is a direct local-quadratic fit, not locfit's tree interpolation.</item>
/// <item>Genes whose IRLS fails to converge are refitted with a projected Newton method instead of R's L-BFGS-B
/// (same optimum, different stopping point).</item>
/// <item>When the residual degrees of freedom m − p are ≤ 3 DESeq2 estimates the dispersion prior variance by a seeded
/// simulation; this port reproduces R's Mersenne-Twister stream, <c>rchisq</c>/<c>rnorm</c>, <c>hist</c> and
/// <c>loess</c> (k-d tree interpolation), see the airway parity test.</item>
/// <item>Observation weights, normalisation-factor matrices, <c>betaPrior = TRUE</c> and glmGamPoi are not supported.</item>
/// </list>
/// </remarks>
public static class DESeq
{
    /// <summary>Fits the full model for <paramref name="design"/> (Wald tests, or an LRT when
    /// <see cref="DESeqOptions.ReducedDesign"/> is set).</summary>
    public static DESeqFit Fit(Dataset counts, DesignMatrix design, DESeqOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentNullException.ThrowIfNull(design);
        options ??= new DESeqOptions();
        int n = counts.Rows, m = counts.Columns, p = design.Coefficients;
        if (design.Samples != m) throw new ArgumentException($"Design has {design.Samples} rows but counts have {m} samples.");
        if (m <= p) throw new ArgumentException("DESeq needs more samples than coefficients (no residual degrees of freedom).");
        var reduced = options.ReducedDesign;
        if (reduced is not null)
        {
            if (reduced.Samples != m) throw new ArgumentException("The reduced design needs one row per sample.");
            if (reduced.Coefficients >= p) throw new ArgumentException("The reduced design must have fewer coefficients than the full design.");
        }

        var fit = new DESeqFit(counts, design) { ReducedDesign = reduced, MinMu = options.MinMu };
        var sf = options.SizeFactors ?? (options.SizeFactorType == SizeFactorType.PosCounts
            ? Normalization.PosCountsSizeFactors(counts)
            : Normalization.MedianOfRatiosSizeFactors(counts));
        if (sf.Length != m) throw new ArgumentException("Need one size factor per sample.");
        fit.SizeFactors = sf;

        var cm = new CountMatrix(counts.Values.ToArray(), n, m);
        var st = new PipelineState(n, m, p);
        BaseMeansAndVariances(cm, sf, st);
        var rows = Enumerable.Range(0, n).Where(i => !st.AllZero[i]).ToArray();

        // 1. Gene-wise dispersions.
        GeneWiseDispersions(cm, rows, sf, design, options, st);

        // 2. Trend, varLogDispEsts and prior variance (from the initial pass only, reused by the outlier refit).
        FitTrend(rows, options, st, fit);
        // 3. MAP dispersions.
        MapDispersions(cm, rows, design, options, st, fit.DispersionPriorVariance, fit.VarLogDispersionEstimates,
            options.DispersionPriorVariance is not null);
        // 4. GLM fit, tests and Cook's distances.
        TestRows(cm, rows, sf, design, reduced, options, st, fit, recordAssays: true);
        RecordMaxCooks(st, design, st.Cooks);

        // 5. Outlier replacement and refit (DESeq2 refitWithoutOutliers).
        var replaceable = DESeqEngine.NOrMoreInCell(design, options.MinReplicatesForReplace);
        fit.ReplaceableSamples = replaceable;
        if (replaceable.Any(b => b)) ReplaceOutliers(cm, sf, design, reduced, options, st, fit, replaceable);

        fit.BaseMean = st.BaseMean; fit.BaseVar = st.BaseVar; fit.AllZero = st.AllZero;
        fit.GeneDispersions = st.GeneDisp; fit.GeneDispersionIterations = st.GeneIter;
        fit.TrendDispersions = st.DispFit; fit.MapDispersions = st.DispMap; fit.DispersionIterations = st.DispIter;
        fit.Dispersions = st.Disp; fit.DispersionOutliers = st.Outlier;
        fit.Coefficients = st.Coef; fit.StandardErrors = st.Se; fit.Converged = st.Conv; fit.BetaIterations = st.BetaIter;
        fit.Deviance = st.Deviance; fit.MaxCooks = st.MaxCooks; fit.Replaced = st.Replace;
        fit.Mu = st.Mu; fit.HatDiagonals = st.Hat; fit.CooksDistances = st.Cooks;
        if (reduced is not null) { fit.LrtStatistic = st.LrtStat; fit.LrtPValue = st.LrtP; }
        return fit;
    }

    /// <summary>Likelihood-ratio test of <paramref name="full"/> against <paramref name="reduced"/> (DESeq2 <c>test = "LRT"</c>).</summary>
    public static DESeqFit FitLrt(Dataset counts, DesignMatrix full, DesignMatrix reduced, DESeqOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(reduced);
        var o = options ?? new DESeqOptions();
        return Fit(counts, full, new DESeqOptions
        {
            MinDispersion = o.MinDispersion, FitType = o.FitType, SizeFactorType = o.SizeFactorType, ShrinkDispersions = o.ShrinkDispersions,
            OutlierSD = o.OutlierSD, MaxIterations = o.MaxIterations, Lambda = o.Lambda, BetaTolerance = o.BetaTolerance, MinMu = o.MinMu,
            UseOptim = o.UseOptim, IndependentFiltering = o.IndependentFiltering, Alpha = o.Alpha, SizeFactors = o.SizeFactors,
            CooksFiltering = o.CooksFiltering, CooksCutoff = o.CooksCutoff, MinReplicatesForReplace = o.MinReplicatesForReplace,
            DispersionPriorVariance = o.DispersionPriorVariance, ReducedDesign = reduced,
        });
    }

    /// <summary>
    /// One-call two-group analysis using a column annotation: fits <c>~ factor</c> (plus optional covariates) with
    /// <paramref name="reference"/> as baseline and returns the Wald test of <paramref name="treatment"/> vs reference.
    /// </summary>
    public static DifferentialExpressionResult Run(Dataset counts, string factor, string reference, string treatment,
        DESeqOptions? options = null, params string[] covariates)
    {
        var design = DesignMatrix.FromAnnotations(counts, factor, reference, covariates);
        var fit = Fit(counts, design, options);
        return fit.Results($"{factor}_{treatment}_vs_{reference}", options);
    }

    // ---- pipeline ---------------------------------------------------------------------------------------------

    private sealed class PipelineState(int n, int m, int p)
    {
        public double[] BaseMean = new double[n], BaseVar = new double[n];
        public bool[] AllZero = new bool[n];
        public double[] GeneDisp = Nan(n), DispFit = Nan(n), DispMap = Nan(n), Disp = Nan(n);
        public int[] GeneIter = new int[n], DispIter = new int[n], BetaIter = new int[n];
        public bool[] Outlier = new bool[n], Conv = new bool[n], Replace = new bool[n];
        public double[] GeneMu = Nan(n * m);
        public double[] Coef = Nan(n * p), Se = Nan(n * p), Deviance = Nan(n), MaxCooks = Nan(n);
        public double[] Mu = Nan(n * m), Hat = Nan(n * m), Cooks = Nan(n * m);
        public double[] LrtStat = Nan(n), LrtP = Nan(n);
        public Func<double, double> TrendFn = _ => double.NaN;
        private static double[] Nan(int k) { var a = new double[k]; Array.Fill(a, double.NaN); return a; }
    }

    private static void BaseMeansAndVariances(CountMatrix cm, double[] sf, PipelineState st)
    {
        int n = cm.Rows, m = cm.Columns;
        var norm = new double[m];
        for (int i = 0; i < n; i++)
        {
            var y = cm.Row(i);
            double any = 0;
            for (int j = 0; j < m; j++) { norm[j] = y[j] / sf[j]; any += y[j]; }
            double mean = RMath.Mean(norm);
            double v = 0;
            for (int j = 0; j < m; j++) { double d = norm[j] - mean; v += d * d; }
            st.BaseMean[i] = mean;
            st.BaseVar[i] = v / (m - 1);
            st.AllZero[i] = any == 0;
        }
    }

    /// <summary>DESeq2 <c>estimateDispersionsGeneEst</c> (niter = 1) for the given rows.</summary>
    private static void GeneWiseDispersions(CountMatrix cm, int[] rows, double[] sf, DesignMatrix x, DESeqOptions o, PipelineState st)
    {
        int m = x.Samples, p = x.Coefficients;
        double minDisp = o.MinDispersion, maxDisp = Math.Max(10, m);
        if (Math.Log(minDisp / 10) <= -30) throw new ArgumentException("log(minDisp/10) must be above -30.");
        var gramInv = DESeqEngine.Inverse(x.Gram(), p) ?? throw new ArgumentException("The design matrix is not of full rank.");
        double xim = 0;
        foreach (double s in sf) xim += 1 / s;
        xim /= m;
        // Initial dispersions: min(rough, moments), bounded.
        var alphaInit = new double[rows.Length];
        for (int r = 0; r < rows.Length; r++)
        {
            int i = rows[r];
            var y = cm.Row(i);
            var norm = new double[m];
            for (int j = 0; j < m; j++) norm[j] = y[j] / sf[j];
            var muLin = DESeqEngine.LinearModelMu(norm, x, gramInv);
            double est = 0;
            for (int j = 0; j < m; j++)
            {
                double mu = Math.Max(1, muLin[j]);
                est += ((norm[j] - mu) * (norm[j] - mu) - mu) / (mu * mu);
            }
            double rough = Math.Max(est / (m - p), 0);
            double moments = (st.BaseVar[i] - xim * st.BaseMean[i]) / (st.BaseMean[i] * st.BaseMean[i]);
            double a = Math.Min(rough, moments);
            if (double.IsNaN(rough) || double.IsNaN(moments)) a = double.NaN;
            alphaInit[r] = Math.Min(Math.Max(minDisp, a), maxDisp);
        }
        DESeqEngine.RowGroups(x, out int groups);
        bool linearMu = groups == p;
        double[] muRows = new double[rows.Length * m];
        if (linearMu)
        {
            for (int r = 0; r < rows.Length; r++)
            {
                var y = cm.Row(rows[r]);
                var norm = new double[m];
                for (int j = 0; j < m; j++) norm[j] = y[j] / sf[j];
                var mu = DESeqEngine.LinearModelMu(norm, x, gramInv);
                for (int j = 0; j < m; j++) muRows[r * m + j] = mu[j] * sf[j];
            }
        }
        else
        {
            var lambda = Enumerable.Repeat(1e-6, p).ToArray();
            var glm = DESeqEngine.FitNbinomGlms(cm, rows, sf, x, alphaInit, lambda, o.MinMu, o.BetaTolerance, o.MaxIterations, o.UseOptim);
            Array.Copy(glm.Mu, muRows, muRows.Length);
        }
        for (int k = 0; k < muRows.Length; k++) if (muRows[k] < o.MinMu) muRows[k] = o.MinMu;
        double minLogAlpha = Math.Log(minDisp / 10);
        Parallel.For(0, rows.Length, () => (new double[p * p], new double[p * p]), (r, _, ws) =>
        {
            int i = rows[r];
            var y = cm.Row(i);
            var mu = muRows.AsSpan(r * m, m);
            var res = DESeqEngine.FitDisp(y, mu, x, Math.Log(alphaInit[r]), Math.Log(alphaInit[r]), 1, false, minLogAlpha,
                1.0, 1e-6, o.MaxIterations, ws.Item1, ws.Item2);
            double est = Math.Min(Math.Exp(res.LogAlpha), maxDisp);
            if (res.LastLp < res.InitialLp + Math.Abs(res.InitialLp) / 1e6) est = alphaInit[r];
            bool conv = res.Iterations < o.MaxIterations && res.Iterations != 1;
            if (!conv && est > minDisp * 10)
                est = DESeqEngine.FitDispGrid(y, mu, x, 0, 1, false, ws.Item1);
            st.GeneDisp[i] = Math.Min(Math.Max(est, minDisp), maxDisp);
            st.GeneIter[i] = res.Iterations;
            mu.CopyTo(st.GeneMu.AsSpan(i * m, m));
            return ws;
        }, _ => { });
    }

    private static void FitTrend(int[] rows, DESeqOptions o, PipelineState st, DESeqFit fit)
    {
        double minDisp = o.MinDispersion;
        var useMeans = new List<double>(); var useDisps = new List<double>();
        foreach (int i in rows)
            if (st.GeneDisp[i] > 100 * minDisp) { useMeans.Add(st.BaseMean[i]); useDisps.Add(st.GeneDisp[i]); }
        if (useMeans.Count == 0)
            throw new InvalidOperationException("All gene-wise dispersion estimates are within 2 orders of magnitude of the minimum; the trend cannot be fitted.");
        var type = o.FitType;
        Func<double, double> fn = _ => double.NaN;
        (double, double) coefs = (double.NaN, double.NaN);
        if (type == DispersionFitType.Parametric)
        {
            if (DESeqEngine.TryParametricFit(useMeans.ToArray(), useDisps.ToArray(), out double a0, out double a1))
            {
                fn = q => a0 + a1 / q;
                coefs = (a0, a1);
            }
            else type = DispersionFitType.Local;
        }
        if (type == DispersionFitType.Local)
            fn = DESeqEngine.LocalDispersionFit(useMeans.ToArray(), useDisps.ToArray(), minDisp);
        if (type == DispersionFitType.Mean)
        {
            var forMean = rows.Where(i => st.GeneDisp[i] > 10 * minDisp).Select(i => st.GeneDisp[i]).ToArray();
            double meanDisp = RMath.TrimmedMean(forMean, 0.001);
            fn = _ => meanDisp;
            coefs = (meanDisp, 0);
        }
        st.TrendFn = fn;
        foreach (int i in rows) st.DispFit[i] = fn(st.BaseMean[i]);
        var resid = new List<double>();
        foreach (int i in rows)
            if (st.GeneDisp[i] >= 1e-8 * 100) resid.Add(Math.Log(st.GeneDisp[i]) - Math.Log(st.DispFit[i]));
        double varLog = resid.Count > 0 ? Math.Pow(NumericUtil.Mad(resid.ToArray()), 2) : double.NaN;
        fit.Trend = coefs; fit.TrendType = type; fit.TrendFunction = fn; fit.VarLogDispersionEstimates = varLog;

        int m = fit.Design.Samples, p = fit.Design.Coefficients;
        double priorVar;
        if (o.DispersionPriorVariance is double given) priorVar = given;
        else if (m - p <= 3)
        {
            var above = rows.Where(i => st.GeneDisp[i] >= minDisp * 100)
                .Select(i => Math.Log(st.GeneDisp[i]) - Math.Log(st.DispFit[i])).ToArray();
            priorVar = DESeqEngine.SimulatedPriorVar(above, m - p);
        }
        else priorVar = Math.Max(varLog - NumericUtil.Trigamma((m - p) / 2.0), 0.25);
        fit.DispersionPriorVariance = priorVar;
    }

    /// <summary>DESeq2 <c>estimateDispersionsMAP</c> for the given rows.</summary>
    private static void MapDispersions(CountMatrix cm, int[] rows, DesignMatrix x, DESeqOptions o, PipelineState st,
        double priorVar, double varLogDispEsts, bool priorGiven)
    {
        int m = x.Samples, p = x.Coefficients;
        double minDisp = o.MinDispersion, maxDisp = Math.Max(10, m);
        if (!o.ShrinkDispersions)
        {
            foreach (int i in rows) { st.Disp[i] = st.GeneDisp[i]; st.DispMap[i] = st.GeneDisp[i]; }
            return;
        }
        // Only when the prior variance has to be estimated (not in the outlier refit, which passes it in).
        bool anyAbove = rows.Any(i => st.GeneDisp[i] >= minDisp * 100);
        if (!anyAbove && !priorGiven)
        {
            foreach (int i in rows) { st.Disp[i] = minDisp * 10; st.DispMap[i] = double.NaN; }
            return;
        }
        double minLogAlpha = Math.Log(minDisp / 10);
        Parallel.For(0, rows.Length, () => (new double[p * p], new double[p * p]), (r, _, ws) =>
        {
            int i = rows[r];
            var y = cm.Row(i);
            var mu = st.GeneMu.AsSpan(i * m, m);
            double fitv = st.DispFit[i];
            double init = st.GeneDisp[i] > 0.1 * fitv ? st.GeneDisp[i] : fitv;
            if (double.IsNaN(init)) init = fitv;
            var res = DESeqEngine.FitDisp(y, mu, x, Math.Log(init), Math.Log(fitv), priorVar, true, minLogAlpha,
                1.0, 1e-6, o.MaxIterations, ws.Item1, ws.Item2);
            double map = Math.Exp(res.LogAlpha);
            if (res.Iterations >= o.MaxIterations)
                map = DESeqEngine.FitDispGrid(y, mu, x, Math.Log(fitv), priorVar, true, ws.Item1);
            map = Math.Min(Math.Max(map, minDisp), maxDisp);
            st.DispMap[i] = map;
            st.DispIter[i] = res.Iterations;
            bool outlier = Math.Log(st.GeneDisp[i]) > Math.Log(fitv) + o.OutlierSD * Math.Sqrt(varLogDispEsts);
            st.Outlier[i] = outlier;
            st.Disp[i] = outlier ? st.GeneDisp[i] : map;
            return ws;
        }, _ => { });
    }

    /// <summary>nbinomWaldTest / nbinomLRT for the given rows, including Cook's distances.</summary>
    private static void TestRows(CountMatrix cm, int[] rows, double[] sf, DesignMatrix x, DesignMatrix? reduced, DESeqOptions o,
        PipelineState st, DESeqFit fit, bool recordAssays)
    {
        int m = x.Samples, p = x.Coefficients;
        var alpha = rows.Select(i => st.Disp[i]).ToArray();
        var lambda = Enumerable.Repeat(o.Lambda, p).ToArray();
        var full = DESeqEngine.FitNbinomGlms(cm, rows, sf, x, alpha, lambda, o.MinMu, o.BetaTolerance, o.MaxIterations, o.UseOptim);
        fit.OptimRefits += full.OptimRows;
        GlmFitResult? red = null;
        if (reduced is not null)
        {
            var lr = Enumerable.Repeat(o.Lambda, reduced.Coefficients).ToArray();
            red = DESeqEngine.FitNbinomGlms(cm, rows, sf, reduced, alpha, lr, o.MinMu, o.BetaTolerance, o.MaxIterations, o.UseOptim);
            fit.OptimRefits += red.OptimRows;
        }
        int df = reduced is null ? 0 : p - reduced.Coefficients;
        var cells = DESeqEngine.CookCells(x);
        var cooks = new double[rows.Length * m];
        Parallel.For(0, rows.Length, r =>
        {
            int i = rows[r];
            var y = cm.Row(i);
            var norm = new double[m];
            for (int j = 0; j < m; j++) norm[j] = y[j] / sf[j];
            double rd = DESeqEngine.RobustMethodOfMomentsDisp(norm, cells.CellOf, cells.HasThree, cells.Cells, cells.AnyThree);
            for (int j = 0; j < m; j++)
            {
                double mu = full.Mu[r * m + j], h = full.Hat[r * m + j];
                double v = mu + rd * mu * mu;
                double pr = (y[j] - mu) * (y[j] - mu) / v;
                cooks[r * m + j] = pr / p * h / ((1 - h) * (1 - h));
            }
        });
        for (int r = 0; r < rows.Length; r++)
        {
            int i = rows[r];
            for (int k = 0; k < p; k++) { st.Coef[i * p + k] = full.Beta[r * p + k]; st.Se[i * p + k] = full.Se[r * p + k]; }
            st.Conv[i] = full.Converged[r];
            st.BetaIter[i] = full.Iterations[r];
            st.Deviance[i] = -2 * full.LogLike[r];
            if (red is not null)
            {
                double stat = 2 * (full.LogLike[r] - red.LogLike[r]);
                st.LrtStat[i] = stat;
                st.LrtP[i] = RMath.ChiSquareUpper(stat, df);
            }
            if (recordAssays)
                for (int j = 0; j < m; j++)
                {
                    st.Mu[i * m + j] = full.Mu[r * m + j];
                    st.Hat[i * m + j] = full.Hat[r * m + j];
                    st.Cooks[i * m + j] = cooks[r * m + j];
                }
            else
            {
                // Refit rows: DESeq2 records maxCooks from the refit before overwriting it (see ReplaceOutliers).
            }
        }
    }

    private static void RecordMaxCooks(PipelineState st, DesignMatrix x, double[] cooks)
    {
        int n = st.AllZero.Length, m = x.Samples, p = x.Coefficients;
        var forCooks = DESeqEngine.NOrMoreInCell(x, 3);
        bool any = forCooks.Any(b => b);
        for (int i = 0; i < n; i++)
        {
            if (!(m > p && any)) { st.MaxCooks[i] = double.NaN; continue; }
            double mx = double.NegativeInfinity;
            for (int j = 0; j < m; j++)
                if (forCooks[j]) { double c = cooks[i * m + j]; if (double.IsNaN(c)) { mx = double.NaN; break; } mx = Math.Max(mx, c); }
            st.MaxCooks[i] = mx;
        }
    }

    /// <summary>DESeq2 <c>replaceOutliers</c> + <c>refitWithoutOutliers</c>.</summary>
    private static void ReplaceOutliers(CountMatrix cm, double[] sf, DesignMatrix x, DesignMatrix? reduced, DESeqOptions o,
        PipelineState st, DESeqFit fit, bool[] replaceable)
    {
        int n = cm.Rows, m = x.Samples, p = x.Coefficients;
        double cutoff = RMath.Qf(0.99, p, m - p);
        var newCounts = cm.Clone();
        int nrefit = 0;
        for (int i = 0; i < n; i++)
        {
            if (st.AllZero[i]) continue;
            bool rep = false;
            for (int j = 0; j < m; j++) if (st.Cooks[i * m + j] > cutoff) rep = true;
            st.Replace[i] = rep;
            if (!rep) continue;
            nrefit++;
            var norm = new double[m];
            for (int j = 0; j < m; j++) norm[j] = cm[i, j] / sf[j];
            double trimMean = RMath.TrimmedMean(norm, 0.2);
            for (int j = 0; j < m; j++)
                if (replaceable[j] && st.Cooks[i * m + j] > cutoff)
                    newCounts[i, j] = Math.Truncate(trimMean * sf[j]);
        }
        if (nrefit == 0) return;
        var origAllZero = (bool[])st.AllZero.Clone();
        BaseMeansAndVariances(newCounts, sf, st);
        var refitRows = Enumerable.Range(0, n).Where(i => st.Replace[i] && !st.AllZero[i]).ToArray();
        var newAllZero = Enumerable.Range(0, n).Where(i => st.Replace[i] && st.AllZero[i]).ToArray();
        fit.ReplacementCounts = new Dataset(fit.Counts.Name, fit.Counts.RowNames, fit.Counts.ColumnNames, newCounts.Values);
        if (refitRows.Length > 0)
        {
            GeneWiseDispersions(newCounts, refitRows, sf, x, o, st);
            foreach (int i in refitRows) st.DispFit[i] = st.TrendFn(st.BaseMean[i]);
            MapDispersions(newCounts, refitRows, x, o, st, fit.DispersionPriorVariance, fit.VarLogDispersionEstimates, priorGiven: true);
            TestRows(newCounts, refitRows, sf, x, reduced, o, st, fit, recordAssays: false);
        }
        foreach (int i in newAllZero)
        {
            for (int k = 0; k < p; k++) { st.Coef[i * p + k] = double.NaN; st.Se[i * p + k] = double.NaN; }
            st.Deviance[i] = double.NaN; st.LrtStat[i] = double.NaN; st.LrtP[i] = double.NaN; st.Conv[i] = false;
        }
        if (replaceable.All(b => b))
            Array.Fill(st.MaxCooks, double.NaN);
        else
        {
            var rc = (double[])st.Cooks.Clone();
            for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) if (replaceable[j]) rc[i * m + j] = 0;
            RecordMaxCooks(st, x, rc);
        }
        _ = origAllZero;
    }

    // ---- multiple testing -------------------------------------------------------------------------------------

    /// <summary>
    /// Benjamini–Hochberg like R's <c>p.adjust(p, "BH")</c>: NaN p-values stay NaN and do not count towards the
    /// number of tests (they are compacted out before calling the kernel).
    /// </summary>
    internal static double[] Bh(double[] p)
    {
        var q = new double[p.Length];
        int k = 0;
        foreach (double v in p) if (!double.IsNaN(v)) k++;
        if (k == p.Length)
        {
            BioKernels.BenjaminiHochberg(p, q);
            return q;
        }
        var finite = new double[k];
        var at = new int[k];
        k = 0;
        for (int i = 0; i < p.Length; i++)
        {
            if (double.IsNaN(p[i])) { q[i] = double.NaN; continue; }
            finite[k] = p[i]; at[k++] = i;
        }
        var fq = new double[k];
        if (k > 0) BioKernels.BenjaminiHochberg(finite, fq);
        for (int t = 0; t < k; t++) q[at[t]] = fq[t];
        return q;
    }

    /// <summary>
    /// DESeq2's independent filtering (<c>pvalueAdjustment</c> with genefilter's <c>filtered_p</c>): BH-adjusts the
    /// p-values of genes whose base mean is at least each of 50 quantile thresholds, fits a lowess curve to the number of
    /// rejections, and uses the first threshold whose rejections exceed the curve's maximum minus its residual RMSE
    /// (falling back to 90% / 80% of the maximum). Filtered genes get NaN.
    /// </summary>
    internal static double[] IndependentFilter(double[] baseMean, double[] p, double alpha)
        => IndependentFilter(baseMean, p, alpha, out _, out _);

    internal static double[] IndependentFilter(double[] baseMean, double[] p, double alpha, out double threshold, out double theta)
    {
        int n = p.Length;
        int zeros = 0;
        foreach (double b in baseMean) if (b == 0) zeros++;
        double lower = (double)zeros / n, upper = lower < 0.95 ? 0.95 : 1;
        var sortedF = (double[])baseMean.Clone();
        Array.Sort(sortedF);
        const int nTheta = 50;
        var thetas = DESeqEngine.RSeq(lower, upper, nTheta);
        var numRej = new double[nTheta];
        var padjs = new double[nTheta][];
        var cutoffs = new double[nTheta];
        for (int t = 0; t < nTheta; t++)
        {
            double cutoff = NumericUtil.QuantileSorted(sortedF, thetas[t]);
            cutoffs[t] = cutoff;
            var pk = new double[n];
            for (int i = 0; i < n; i++) pk[i] = baseMean[i] >= cutoff ? p[i] : double.NaN;
            var q = Bh(pk);
            padjs[t] = q;
            int rej = 0;
            foreach (double v in q) if (v < alpha) rej++;
            numRej[t] = rej;
        }
        int j = 0;
        double maxRej = numRej.Max();
        if (maxRej > 10)
        {
            var lo = NumericUtil.Lowess(thetas, numRej, 1.0 / 5);
            double maxFit = lo.Max(), ss = 0; int cnt = 0;
            for (int t = 0; t < nTheta; t++)
                if (numRej[t] > 0) { double r = numRej[t] - lo[t]; ss += r * r; cnt++; }
            double thresh = maxFit - Math.Sqrt(cnt > 0 ? ss / cnt : 0);
            int Find(double th) { for (int t = 0; t < nTheta; t++) if (numRej[t] > th) return t; return -1; }
            j = Find(thresh);
            if (j < 0) j = Find(0.9 * maxFit);
            if (j < 0) j = Find(0.8 * maxFit);
            if (j < 0) j = 0;
        }
        threshold = cutoffs[j];
        theta = thetas[j];
        return padjs[j];
    }
}
