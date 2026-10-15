using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>How the mean–dispersion trend is fitted.</summary>
public enum DispersionFitType
{
    /// <summary>DESeq2's parametric trend α(μ) = asymptDisp + extraPois / μ, fitted by a gamma-family GLM.</summary>
    Parametric,

    /// <summary>A constant trend: the (trimmed) mean of the gene-wise estimates.</summary>
    Mean,
}

/// <summary>Options for <see cref="DESeq"/>.</summary>
public sealed class DESeqOptions
{
    /// <summary>Lower bound for dispersions (DESeq2 <c>minDisp</c>, default 1e-8).</summary>
    public double MinDispersion { get; init; } = 1e-8;

    /// <summary>Trend type (default parametric, falling back to mean if the parametric fit fails).</summary>
    public DispersionFitType FitType { get; init; } = DispersionFitType.Parametric;

    /// <summary>Shrink gene-wise dispersions towards the trend (maximum a posteriori). Default true.</summary>
    public bool ShrinkDispersions { get; init; } = true;

    /// <summary>Gene-wise estimates more than this many prior SDs above the trend are kept unshrunk (default 2).</summary>
    public double OutlierSD { get; init; } = 2;

    /// <summary>Maximum IRLS iterations for the GLM (default 100).</summary>
    public int MaxIterations { get; init; } = 100;

    /// <summary>Ridge penalty on log2-scale coefficients (DESeq2 uses 1e-6 when no LFC prior is used).</summary>
    public double Lambda { get; init; } = 1e-6;

    /// <summary>Optimise the BH threshold by filtering on the mean of normalised counts (default true).</summary>
    public bool IndependentFiltering { get; init; } = true;

    /// <summary>Significance level the independent filter is tuned for (DESeq2 <c>alpha</c>, default 0.1).</summary>
    public double Alpha { get; init; } = 0.1;

    /// <summary>Precomputed size factors; estimated by median-of-ratios when null.</summary>
    public double[]? SizeFactors { get; init; }
}

/// <summary>
/// The fitted negative-binomial model of every gene: size factors, dispersions at each stage, and the GLM
/// coefficients with their standard errors (log2 scale). Produce result tables with <see cref="Results(int, DESeqOptions?)"/>.
/// </summary>
public sealed class DESeqFit
{
    internal DESeqFit(Dataset counts, DesignMatrix design) { Counts = counts; Design = design; }

    /// <summary>The count matrix that was fitted.</summary>
    public Dataset Counts { get; }

    /// <summary>The model matrix.</summary>
    public DesignMatrix Design { get; }

    /// <summary>Per-sample size factors.</summary>
    public double[] SizeFactors { get; internal set; } = [];

    /// <summary>Mean of normalised counts per gene.</summary>
    public double[] BaseMean { get; internal set; } = [];

    /// <summary>Variance of normalised counts per gene.</summary>
    public double[] BaseVar { get; internal set; } = [];

    /// <summary>Gene-wise maximum Cox–Reid adjusted profile likelihood dispersions (<c>dispGeneEst</c>).</summary>
    public double[] GeneDispersions { get; internal set; } = [];

    /// <summary>Fitted trend at each gene's mean (<c>dispFit</c>).</summary>
    public double[] TrendDispersions { get; internal set; } = [];

    /// <summary>Final dispersions used for testing (<c>dispersion</c>): MAP estimates, or gene-wise for outliers.</summary>
    public double[] Dispersions { get; internal set; } = [];

    /// <summary>True where the gene-wise estimate was flagged as a dispersion outlier and kept.</summary>
    public bool[] DispersionOutliers { get; internal set; } = [];

    /// <summary>Parametric trend coefficients (asymptDisp, extraPois); for a mean trend extraPois is 0.</summary>
    public (double AsymptoticDispersion, double ExtraPoisson) Trend { get; internal set; }

    /// <summary>Which trend was actually used (a failed parametric fit falls back to mean).</summary>
    public DispersionFitType TrendType { get; internal set; }

    /// <summary>Prior variance of log dispersions around the trend (<c>dispPriorVar</c>).</summary>
    public double DispersionPriorVariance { get; internal set; }

    /// <summary>Coefficients on the log2 scale, genes × coefficients (row-major).</summary>
    public double[] Coefficients { get; internal set; } = [];

    /// <summary>Standard errors of the log2 coefficients, genes × coefficients (row-major).</summary>
    public double[] StandardErrors { get; internal set; } = [];

    /// <summary>Whether IRLS converged for each gene.</summary>
    public bool[] Converged { get; internal set; } = [];

    /// <summary>True for genes with only zero counts (excluded from testing).</summary>
    public bool[] AllZero { get; internal set; } = [];

    /// <summary>Wald test of the coefficient at <paramref name="coefficient"/> (default: the last one).</summary>
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
            if (AllZero[i]) { lfc[i] = se[i] = stat[i] = pv[i] = double.NaN; continue; }
            lfc[i] = Coefficients[i * p + coefficient];
            se[i] = StandardErrors[i * p + coefficient];
            stat[i] = lfc[i] / se[i];
            pv[i] = double.IsNaN(stat[i]) ? double.NaN : 2 * NumericUtil.NormalUpper(Math.Abs(stat[i]));
        }
        var padj = options.IndependentFiltering ? DESeq.IndependentFilter(BaseMean, pv, options.Alpha) : DESeq.Bh(pv);
        return new DifferentialExpressionResult("DESeq2-Wald", Design.CoefficientNames[coefficient], Counts.RowNames,
            (double[])BaseMean.Clone(), lfc, se, stat, pv, padj);
    }

    /// <summary>Wald test of a named coefficient.</summary>
    public DifferentialExpressionResult Results(string coefficient, DESeqOptions? options = null)
        => Results(Design.IndexOf(coefficient), options);

    /// <summary>
    /// DESeq2's variance-stabilising transformation using the fitted parametric trend (the closed form of
    /// <c>getVarianceStabilizedData</c> with <c>fitType = "parametric"</c>, applied to normalised counts). With a mean
    /// trend the closed form for constant dispersion is used. Output is on a log2-like scale.
    /// </summary>
    public Dataset VarianceStabilizingTransform()
        => Transforms.Vst(Counts, SizeFactors, Trend.AsymptoticDispersion, Trend.ExtraPoisson);
}

/// <summary>
/// DESeq2-style differential expression (Love, Huber &amp; Anders 2014): median-of-ratios size factors, gene-wise
/// negative-binomial dispersions by Cox–Reid adjusted profile likelihood, a parametric mean–dispersion trend,
/// empirical-Bayes (MAP) shrinkage of dispersions towards the trend, an IRLS negative-binomial GLM per gene, Wald
/// tests and Benjamini–Hochberg FDR with independent filtering.
/// </summary>
/// <remarks>
/// <para>Differences from DESeq2 (documented simplifications):</para>
/// <list type="bullet">
/// <item>Dispersion likelihoods are maximised by a 30-point grid on log α refined by golden-section search instead
/// of DESeq2's Armijo line search (DESeq2 itself falls back to a grid); a single estimation pass is made (DESeq2's
/// default <c>niter = 1</c>), with μ fitted by the GLM at the rough moment estimate.</item>
/// <item>When residual degrees of freedom are ≤ 3, the expected sampling variance of log dispersions is taken as
/// trigamma((m−p)/2) instead of DESeq2's simulation.</item>
/// <item>No Cook's distance outlier flagging/replacement and no LFC shrinkage (apeglm/ashr); the "mean" trend is an
/// untrimmed mean; local (locfit) trends are not offered.</item>
/// </list>
/// With these exceptions log2 fold changes, standard errors and p-values agree with DESeq2 to a few percent; MLE
/// fold changes for group designs equal the ratio of mean normalised counts exactly.
/// </remarks>
public static class DESeq
{
    /// <summary>Fits the full model for <paramref name="design"/>.</summary>
    public static DESeqFit Fit(Dataset counts, DesignMatrix design, DESeqOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentNullException.ThrowIfNull(design);
        options ??= new DESeqOptions();
        int n = counts.Rows, m = counts.Columns, p = design.Coefficients;
        if (design.Samples != m) throw new ArgumentException($"Design has {design.Samples} rows but counts have {m} samples.");
        if (m <= p) throw new ArgumentException("DESeq needs more samples than coefficients (no residual degrees of freedom).");

        var fit = new DESeqFit(counts, design);
        var sf = options.SizeFactors ?? Normalization.MedianOfRatiosSizeFactors(counts);
        if (sf.Length != m) throw new ArgumentException("Need one size factor per sample.");
        fit.SizeFactors = sf;

        double minDisp = options.MinDispersion, maxDisp = Math.Max(10, m);
        double lambdaNat = options.Lambda / (Math.Log(2) * Math.Log(2));
        var baseMean = new double[n]; var baseVar = new double[n]; var allZero = new bool[n];
        var geneDisp = new double[n];
        var muStore = new double[n * m];
        double xim = 0;
        foreach (double s in sf) xim += 1 / s;
        xim /= m;

        // ---- 1. Gene-wise dispersions (Cox–Reid adjusted profile likelihood, μ held fixed) ----
        Parallel.For(0, n, () => new GlmWorkspace(m, p), (i, _, ws) =>
        {
            var y = counts.Row(i);
            double sum = 0, any = 0;
            for (int j = 0; j < m; j++) { double q = y[j] / sf[j]; ws.Norm[j] = q; sum += q; any += y[j]; }
            double mean = sum / m;
            double v = 0;
            for (int j = 0; j < m; j++) { double d = ws.Norm[j] - mean; v += d * d; }
            v /= m - 1;
            baseMean[i] = mean; baseVar[i] = v;
            if (any == 0) { allZero[i] = true; geneDisp[i] = double.NaN; return ws; }

            double moments = (v - xim * mean) / (mean * mean);
            double rough = RoughDispersion(ws.Norm, design);
            double a0 = Math.Min(rough, moments);
            if (double.IsNaN(a0)) a0 = minDisp;
            a0 = Math.Clamp(a0, minDisp, maxDisp);

            FitGlm(y, sf, design, a0, lambdaNat, options.MaxIterations, ws, out double[] _);
            ws.Mu.AsSpan(0, m).CopyTo(muStore.AsSpan(i * m, m));
            double logA = MaximizeLogPosterior(y, ws.Mu, design, double.NaN, 0, Math.Log(minDisp / 10), Math.Log(maxDisp), ws);
            geneDisp[i] = Math.Clamp(Math.Exp(logA), minDisp, maxDisp);
            return ws;
        }, _ => { });
        fit.BaseMean = baseMean; fit.BaseVar = baseVar; fit.AllZero = allZero; fit.GeneDispersions = geneDisp;

        // ---- 2. Mean–dispersion trend ----
        var trend = new double[n];
        double asym = 0, extra = 0;
        var type = options.FitType;
        if (type == DispersionFitType.Parametric && !TryParametricFit(baseMean, geneDisp, minDisp, out asym, out extra))
            type = DispersionFitType.Mean;
        if (type == DispersionFitType.Mean)
        {
            double s = 0; int c = 0;
            for (int i = 0; i < n; i++) if (!allZero[i] && geneDisp[i] > 10 * minDisp) { s += geneDisp[i]; c++; }
            asym = c > 0 ? s / c : minDisp; extra = 0;
        }
        for (int i = 0; i < n; i++) trend[i] = allZero[i] ? double.NaN : asym + extra / baseMean[i];
        fit.Trend = (asym, extra); fit.TrendType = type; fit.TrendDispersions = trend;

        // ---- 3. Prior variance and MAP shrinkage ----
        var resid = new List<double>(n);
        for (int i = 0; i < n; i++)
            if (!allZero[i] && geneDisp[i] >= 100 * minDisp) resid.Add(Math.Log(geneDisp[i]) - Math.Log(trend[i]));
        double varLogDisp = resid.Count > 0 ? Math.Pow(NumericUtil.Mad(resid.ToArray()), 2) : 0.25;
        double expVar = NumericUtil.Trigamma((m - p) / 2.0);
        double priorVar = Math.Max(varLogDisp - expVar, 0.25);
        fit.DispersionPriorVariance = priorVar;

        var finalDisp = new double[n];
        var outlier = new bool[n];
        Parallel.For(0, n, () => new GlmWorkspace(m, p), (i, _, ws) =>
        {
            if (allZero[i]) { finalDisp[i] = double.NaN; return ws; }
            if (!options.ShrinkDispersions) { finalDisp[i] = geneDisp[i]; return ws; }
            var y = counts.Row(i);
            muStore.AsSpan(i * m, m).CopyTo(ws.Mu);
            double logA = MaximizeLogPosterior(y, ws.Mu, design, Math.Log(trend[i]), priorVar, Math.Log(minDisp / 10), Math.Log(maxDisp), ws);
            double map = Math.Clamp(Math.Exp(logA), minDisp, maxDisp);
            if (Math.Log(geneDisp[i]) > Math.Log(trend[i]) + options.OutlierSD * Math.Sqrt(varLogDisp))
            {
                outlier[i] = true;
                map = geneDisp[i];
            }
            finalDisp[i] = map;
            return ws;
        }, _ => { });
        fit.Dispersions = finalDisp; fit.DispersionOutliers = outlier;

        // ---- 4. Final GLM fit and standard errors ----
        var coef = new double[n * p];
        var ses = new double[n * p];
        var conv = new bool[n];
        double log2e = 1 / Math.Log(2);
        Parallel.For(0, n, () => new GlmWorkspace(m, p), (i, _, ws) =>
        {
            if (allZero[i])
            {
                for (int k = 0; k < p; k++) { coef[i * p + k] = double.NaN; ses[i * p + k] = double.NaN; }
                return ws;
            }
            conv[i] = FitGlm(counts.Row(i), sf, design, finalDisp[i], lambdaNat, options.MaxIterations, ws, out var seNat);
            for (int k = 0; k < p; k++)
            {
                coef[i * p + k] = ws.Beta[k] * log2e;
                ses[i * p + k] = seNat[k] * log2e;
            }
            return ws;
        }, _ => { });
        fit.Coefficients = coef; fit.StandardErrors = ses; fit.Converged = conv;
        return fit;
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

    // ---- GLM machinery ----------------------------------------------------------------------------

    internal sealed class GlmWorkspace(int m, int p)
    {
        public readonly double[] Norm = new double[m];
        public readonly double[] Mu = new double[m];
        public readonly double[] W = new double[m];
        public readonly double[] Z = new double[m];
        public readonly double[] Beta = new double[p];
        public readonly double[] Xtwx = new double[p * p];
        public readonly double[] Xtwz = new double[p];
        public readonly double[] Se = new double[p];
    }

    private static double RoughDispersion(double[] norm, DesignMatrix x)
    {
        // DESeq2 roughDispEstimate: least-squares fit of normalised counts, μ floored at 1.
        int m = x.Samples, p = x.Coefficients;
        var g = x.Gram();
        var xty = new double[p];
        for (int j = 0; j < m; j++) for (int k = 0; k < p; k++) xty[k] += x[j, k] * norm[j];
        var b = NumericUtil.SolveSmall(g, xty, p);
        if (b is null) return double.NaN;
        double s = 0;
        for (int j = 0; j < m; j++)
        {
            double mu = 0;
            for (int k = 0; k < p; k++) mu += x[j, k] * b[k];
            mu = Math.Max(mu, 1);
            double r = norm[j] - mu;
            s += (r * r - mu) / (mu * mu);
        }
        return Math.Max(s / (m - p), 0);
    }

    /// <summary>
    /// IRLS for a negative-binomial GLM with log link and log(sizeFactor) offset (DESeq2's fitBeta). Returns
    /// convergence; on exit ws.Beta holds natural-log coefficients, ws.Mu the fitted means, and
    /// <paramref name="seNat"/> the sandwich standard errors.
    /// </summary>
    internal static bool FitGlm(ReadOnlySpan<double> y, ReadOnlySpan<double> sf, DesignMatrix x, double alpha,
        double lambdaNat, int maxIter, GlmWorkspace ws, out double[] seNat)
    {
        int m = x.Samples, p = x.Coefficients;
        const double minMu = 0.5;
        // Start: least squares of log(normalised + 0.1).
        Array.Clear(ws.Xtwz);
        var g = x.Gram();
        for (int j = 0; j < m; j++)
        {
            double ly = Math.Log(y[j] / sf[j] + 0.1);
            for (int k = 0; k < p; k++) ws.Xtwz[k] += x[j, k] * ly;
        }
        var b0 = NumericUtil.SolveSmall(g, (double[])ws.Xtwz.Clone(), p) ?? new double[p];
        Array.Copy(b0, ws.Beta, p);
        UpdateMu(ws.Beta, x, sf, ws.Mu, minMu);
        double dev = Deviance(y, ws.Mu, alpha), devOld;
        bool converged = false;
        for (int it = 0; it < maxIter; it++)
        {
            Array.Clear(ws.Xtwx); Array.Clear(ws.Xtwz);
            for (int j = 0; j < m; j++)
            {
                double mu = ws.Mu[j];
                double w = mu / (1 + alpha * mu);
                double z = Math.Log(mu / sf[j]) + y[j] / mu - 1;
                for (int a = 0; a < p; a++)
                {
                    double xa = x[j, a] * w;
                    ws.Xtwz[a] += xa * z;
                    for (int c = 0; c < p; c++) ws.Xtwx[a * p + c] += xa * x[j, c];
                }
            }
            var lhs = (double[])ws.Xtwx.Clone();
            for (int a = 0; a < p; a++) lhs[a * p + a] += lambdaNat;
            var nb = NumericUtil.SolveSmall(lhs, (double[])ws.Xtwz.Clone(), p);
            if (nb is null) break;
            bool blewUp = false;
            foreach (double v in nb) if (!double.IsFinite(v) || Math.Abs(v) > 30) blewUp = true;
            Array.Copy(nb, ws.Beta, p);
            UpdateMu(ws.Beta, x, sf, ws.Mu, minMu);
            if (blewUp) break;
            devOld = dev;
            dev = Deviance(y, ws.Mu, alpha);
            if (Math.Abs(dev - devOld) / (Math.Abs(dev) + 0.1) < 1e-8) { converged = true; break; }
        }
        // Sandwich covariance (X'WX + Λ)⁻¹ X'WX (X'WX + Λ)⁻¹ at the final μ.
        Array.Clear(ws.Xtwx);
        for (int j = 0; j < m; j++)
        {
            double w = ws.Mu[j] / (1 + alpha * ws.Mu[j]);
            for (int a = 0; a < p; a++)
                for (int c = 0; c < p; c++) ws.Xtwx[a * p + c] += x[j, a] * w * x[j, c];
        }
        var r = (double[])ws.Xtwx.Clone();
        for (int a = 0; a < p; a++) r[a * p + a] += lambdaNat;
        var rinv = NumericUtil.InvertSpd(r, p);
        seNat = ws.Se;
        if (rinv is null) { Array.Fill(seNat, double.NaN); return false; }
        for (int a = 0; a < p; a++)
        {
            double s = 0;
            for (int u = 0; u < p; u++)
                for (int v = 0; v < p; v++) s += rinv[a * p + u] * ws.Xtwx[u * p + v] * rinv[v * p + a];
            seNat[a] = Math.Sqrt(Math.Max(s, 0));
        }
        return converged;
    }

    private static void UpdateMu(double[] beta, DesignMatrix x, ReadOnlySpan<double> sf, double[] mu, double minMu)
    {
        int p = x.Coefficients;
        for (int j = 0; j < x.Samples; j++)
        {
            double eta = 0;
            for (int k = 0; k < p; k++) eta += x[j, k] * beta[k];
            mu[j] = Math.Max(sf[j] * Math.Exp(eta), minMu);
        }
    }

    private static double Deviance(ReadOnlySpan<double> y, double[] mu, double alpha)
    {
        double ll = 0;
        double r = 1 / alpha;
        for (int j = 0; j < y.Length; j++)
            ll += SpecialFunctions.LogGamma(y[j] + r) - SpecialFunctions.LogGamma(r) - SpecialFunctions.LogGamma(y[j] + 1)
                + r * Math.Log(r / (r + mu[j])) + (y[j] > 0 ? y[j] * Math.Log(mu[j] / (r + mu[j])) : 0);
        return -2 * ll;
    }

    /// <summary>
    /// DESeq2's log posterior of log α with μ fixed: NB log-likelihood terms that depend on α, the Cox–Reid
    /// adjustment −½ log det(XᵀWX) with W = 1/(1/μ + α), and (when <paramref name="priorVar"/> &gt; 0) a normal prior
    /// on log α centred at <paramref name="priorMean"/>.
    /// </summary>
    internal static double LogPosterior(double logAlpha, ReadOnlySpan<double> y, double[] mu, DesignMatrix x,
        double priorMean, double priorVar, GlmWorkspace ws)
    {
        double alpha = Math.Exp(logAlpha);
        double r = 1 / alpha;
        int m = x.Samples, p = x.Coefficients;
        double lgr = SpecialFunctions.LogGamma(r);
        double ll = 0;
        for (int j = 0; j < m; j++)
            ll += SpecialFunctions.LogGamma(y[j] + r) - lgr - y[j] * Math.Log(mu[j] + r) - r * Math.Log(1 + mu[j] * alpha);
        Array.Clear(ws.Xtwx);
        for (int j = 0; j < m; j++)
        {
            double w = 1 / (1 / mu[j] + alpha);
            for (int a = 0; a < p; a++)
                for (int c = 0; c < p; c++) ws.Xtwx[a * p + c] += x[j, a] * w * x[j, c];
        }
        double cr = -0.5 * NumericUtil.LogDetSpd(ws.Xtwx, p);
        if (double.IsNaN(cr)) cr = 0;
        double prior = priorVar > 0 ? -0.5 * (logAlpha - priorMean) * (logAlpha - priorMean) / priorVar : 0;
        return ll + cr + prior;
    }

    private static double MaximizeLogPosterior(ReadOnlySpan<double> y, double[] mu, DesignMatrix x, double priorMean,
        double priorVar, double lo, double hi, GlmWorkspace ws)
    {
        const int grid = 30;
        double best = double.NegativeInfinity; int bestK = 0;
        double step = (hi - lo) / (grid - 1);
        for (int k = 0; k < grid; k++)
        {
            double v = LogPosterior(lo + k * step, y, mu, x, priorMean, priorVar, ws);
            if (v > best) { best = v; bestK = k; }
        }
        double a = lo + Math.Max(bestK - 1, 0) * step, b = lo + Math.Min(bestK + 1, grid - 1) * step;
        // Golden-section refinement.
        double gr = (Math.Sqrt(5) - 1) / 2;
        double c = b - gr * (b - a), d = a + gr * (b - a);
        double fc = LogPosterior(c, y, mu, x, priorMean, priorVar, ws), fd = LogPosterior(d, y, mu, x, priorMean, priorVar, ws);
        for (int it = 0; it < 60 && b - a > 1e-6; it++)
        {
            if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = LogPosterior(c, y, mu, x, priorMean, priorVar, ws); }
            else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = LogPosterior(d, y, mu, x, priorMean, priorVar, ws); }
        }
        double xm = 0.5 * (a + b);
        return LogPosterior(xm, y, mu, x, priorMean, priorVar, ws) >= best ? xm : lo + bestK * step;
    }

    /// <summary>
    /// DESeq2's <c>parametricDispersionFit</c>: iteratively fits α = a₀ + a₁/μ by a gamma GLM with identity link,
    /// dropping genes whose ratio to the current fit lies outside (1e-4, 15).
    /// </summary>
    internal static bool TryParametricFit(double[] means, double[] disps, double minDisp, out double a0, out double a1)
    {
        var xs = new List<double>(); var ys = new List<double>();
        for (int i = 0; i < means.Length; i++)
            if (means[i] > 0 && double.IsFinite(disps[i]) && disps[i] >= 100 * minDisp) { xs.Add(1 / means[i]); ys.Add(disps[i]); }
        a0 = 0.1; a1 = 1;
        if (xs.Count < 3) return false;
        for (int iter = 0; iter <= 10; iter++)
        {
            var gx = new List<double>(); var gy = new List<double>();
            for (int i = 0; i < xs.Count; i++)
            {
                double res = ys[i] / (a0 + a1 * xs[i]);
                if (res > 1e-4 && res < 15) { gx.Add(xs[i]); gy.Add(ys[i]); }
            }
            if (gx.Count < 3) return false;
            if (!GammaIdentityGlm(gx, gy, a0, a1, out double n0, out double n1)) return false;
            if (n0 <= 0 || n1 <= 0) return false;
            double change = Math.Pow(Math.Log(n0 / a0), 2) + Math.Pow(Math.Log(n1 / a1), 2);
            a0 = n0; a1 = n1;
            if (change < 1e-6) return true;
        }
        return false;
    }

    private static bool GammaIdentityGlm(List<double> x, List<double> y, double b0, double b1, out double c0, out double c1)
    {
        c0 = b0; c1 = b1;
        double dev = double.MaxValue;
        for (int it = 0; it < 50; it++)
        {
            double s00 = 0, s01 = 0, s11 = 0, t0 = 0, t1 = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double mu = c0 + c1 * x[i];
                if (mu <= 0) return false;
                double w = 1 / (mu * mu);
                s00 += w; s01 += w * x[i]; s11 += w * x[i] * x[i];
                t0 += w * y[i]; t1 += w * x[i] * y[i];
            }
            double det = s00 * s11 - s01 * s01;
            if (Math.Abs(det) < 1e-300) return false;
            c0 = (s11 * t0 - s01 * t1) / det;
            c1 = (s00 * t1 - s01 * t0) / det;
            double nd = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double mu = c0 + c1 * x[i];
                if (mu <= 0) return false;
                nd += 2 * (-Math.Log(y[i] / mu) + (y[i] - mu) / mu);
            }
            if (Math.Abs(nd - dev) / (Math.Abs(nd) + 0.1) < 1e-8) return true;
            dev = nd;
        }
        return true;
    }

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
    /// DESeq2's independent filtering (<c>pvalueAdjustment</c>): BH-adjusts p-values of genes whose base mean is at least
    /// each of 50 quantile thresholds, fits a lowess curve to the number of rejections, and uses the first threshold
    /// whose rejections come within one residual SD of the curve's maximum. Filtered genes get NaN.
    /// </summary>
    internal static double[] IndependentFilter(double[] baseMean, double[] p, double alpha)
    {
        int n = p.Length;
        bool anyTested = false;
        foreach (double v in p) if (!double.IsNaN(v)) { anyTested = true; break; }
        if (!anyTested) return Bh(p);
        // As DESeq2: the filter statistic covers every gene (all-zero genes included), and genefilter's
        // filtered_p keeps genes with filter >= quantile(filter, theta).
        int zeros = 0;
        foreach (double b in baseMean) if (b == 0) zeros++;
        double lower = (double)zeros / n, upper = lower < 0.95 ? 0.95 : 1;
        var sortedF = (double[])baseMean.Clone();
        Array.Sort(sortedF);
        const int nTheta = 50;
        var theta = new double[nTheta];
        var numRej = new double[nTheta];
        var padjs = new double[nTheta][];
        for (int t = 0; t < nTheta; t++)
        {
            theta[t] = lower + (upper - lower) * t / (nTheta - 1);
            double cutoff = NumericUtil.QuantileSorted(sortedF, theta[t]);
            var pk = new double[n];
            for (int i = 0; i < n; i++) pk[i] = double.NaN;
            for (int i = 0; i < n; i++) if (baseMean[i] >= cutoff) pk[i] = p[i];
            var q = Bh(pk);
            padjs[t] = q;
            int rej = 0;
            foreach (double v in q) if (v < alpha) rej++;
            numRej[t] = rej;
        }
        int j = 0;
        double maxRej = 0;
        foreach (double r in numRej) maxRej = Math.Max(maxRej, r);
        if (maxRej > 10)
        {
            var lo = NumericUtil.Lowess(theta, numRej, 1.0 / 5);
            double maxFit = double.NegativeInfinity, ss = 0; int cnt = 0;
            for (int t = 0; t < nTheta; t++)
            {
                maxFit = Math.Max(maxFit, lo[t]);
                if (numRej[t] > 0) { double r = numRej[t] - lo[t]; ss += r * r; cnt++; }
            }
            double thresh = maxFit - Math.Sqrt(cnt > 0 ? ss / cnt : 0);
            j = 0;
            for (int t = 0; t < nTheta; t++) if (numRej[t] > thresh) { j = t; break; }
        }
        return padjs[j];
    }
}
