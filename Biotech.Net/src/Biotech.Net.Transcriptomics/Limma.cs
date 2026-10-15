using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>Per-gene (weighted) least-squares fits, as returned by limma's <c>lmFit</c>.</summary>
public sealed class LinearModelFit
{
    internal LinearModelFit(IReadOnlyList<string> genes, DesignMatrix design, int genesCount)
    {
        Genes = genes;
        Design = design;
        int p = design.Coefficients;
        Coefficients = new double[genesCount * p];
        StdevUnscaled = new double[genesCount * p];
        Sigma = new double[genesCount];
        DfResidual = new double[genesCount];
        Amean = new double[genesCount];
    }

    /// <summary>Gene identifiers.</summary>
    public IReadOnlyList<string> Genes { get; }

    /// <summary>The design matrix.</summary>
    public DesignMatrix Design { get; }

    /// <summary>Coefficients, genes × coefficients (row-major).</summary>
    public double[] Coefficients { get; }

    /// <summary>Unscaled standard deviations √diag((XᵀWX)⁻¹), genes × coefficients.</summary>
    public double[] StdevUnscaled { get; }

    /// <summary>Residual standard deviation per gene.</summary>
    public double[] Sigma { get; }

    /// <summary>Residual degrees of freedom per gene.</summary>
    public double[] DfResidual { get; }

    /// <summary>Average expression per gene (AveExpr).</summary>
    public double[] Amean { get; }

    /// <summary>Prior degrees of freedom d₀ (set by <see cref="Limma.EBayes"/>).</summary>
    public double DfPrior { get; internal set; } = double.NaN;

    /// <summary>Prior variance s₀² (set by <see cref="Limma.EBayes"/>).</summary>
    public double S2Prior { get; internal set; } = double.NaN;

    /// <summary>Posterior (moderated) variances s̃² per gene (set by <see cref="Limma.EBayes"/>).</summary>
    public double[] S2Post { get; internal set; } = [];
}

/// <summary>The voom transformation: log-CPM values plus observation-level precision weights.</summary>
/// <param name="LogCpm">log2((count + 0.5) / (libSize + 1) × 10⁶).</param>
/// <param name="Weights">Precision weights, genes × samples (row-major).</param>
/// <param name="TrendX">Mean log2 count per gene used for the mean–variance trend (x of the voom plot).</param>
/// <param name="TrendY">√ residual SD per gene (y of the voom plot).</param>
public sealed record VoomResult(Dataset LogCpm, double[] Weights, double[] TrendX, double[] TrendY);

/// <summary>
/// limma-style linear modelling with empirical-Bayes variance moderation (Smyth 2004) and the voom
/// mean–variance weighting for RNA-Seq counts (Law et al. 2014).
/// </summary>
/// <remarks>
/// <c>lmFit</c>, <c>eBayes</c> (non-robust <c>fitFDist</c> with no trend) and <c>voom</c> follow limma's formulas,
/// including R's <c>lowess</c> (ported from its C source) for the voom trend. Not implemented: robust/trended eBayes,
/// array weights, duplicate correlation and voom's quality weights.
/// </remarks>
public static class Limma
{
    /// <summary>
    /// Fits <c>y = Xβ + ε</c> by (weighted) least squares for every row of <paramref name="expression"/>
    /// (log-scale values, genes × samples). <paramref name="weights"/> is genes × samples row-major or null.
    /// </summary>
    public static LinearModelFit LmFit(Dataset expression, DesignMatrix design, double[]? weights = null)
    {
        ArgumentNullException.ThrowIfNull(expression);
        ArgumentNullException.ThrowIfNull(design);
        int n = expression.Rows, m = expression.Columns, p = design.Coefficients;
        if (design.Samples != m) throw new ArgumentException("Design rows must match samples.");
        if (weights is not null && weights.Length != n * m) throw new ArgumentException("Weights must be genes × samples.");
        var fit = new LinearModelFit(expression.RowNames, design, n);
        var values = expression.Values.ToArray();
        Parallel.For(0, n, i =>
        {
            var y = values.AsSpan(i * m, m);
            var xtwx = new double[p * p];
            var xtwy = new double[p];
            double am = 0;
            for (int j = 0; j < m; j++)
            {
                am += y[j];
                double w = weights is null ? 1 : weights[i * m + j];
                for (int a = 0; a < p; a++)
                {
                    double xa = design[j, a] * w;
                    xtwy[a] += xa * y[j];
                    for (int c = 0; c < p; c++) xtwx[a * p + c] += xa * design[j, c];
                }
            }
            fit.Amean[i] = am / m;
            var inv = NumericUtil.InvertSpd(xtwx, p);
            if (inv is null)
            {
                for (int a = 0; a < p; a++) { fit.Coefficients[i * p + a] = double.NaN; fit.StdevUnscaled[i * p + a] = double.NaN; }
                fit.Sigma[i] = double.NaN;
                return;
            }
            for (int a = 0; a < p; a++)
            {
                double b = 0;
                for (int c = 0; c < p; c++) b += inv[a * p + c] * xtwy[c];
                fit.Coefficients[i * p + a] = b;
                fit.StdevUnscaled[i * p + a] = Math.Sqrt(inv[a * p + a]);
            }
            double rss = 0;
            for (int j = 0; j < m; j++)
            {
                double f = 0;
                for (int a = 0; a < p; a++) f += design[j, a] * fit.Coefficients[i * p + a];
                double w = weights is null ? 1 : weights[i * m + j];
                rss += w * (y[j] - f) * (y[j] - f);
            }
            int df = m - p;
            fit.DfResidual[i] = df;
            fit.Sigma[i] = df > 0 ? Math.Sqrt(rss / df) : double.NaN;
        });
        return fit;
    }

    /// <summary>
    /// Empirical-Bayes moderation (limma <c>eBayes</c>): estimates the prior d₀ and s₀² from the sample variances by
    /// the method of moments on log s² (<c>fitFDist</c>) and computes posterior variances
    /// s̃² = (d₀s₀² + d·s²)/(d₀ + d). Returns <paramref name="fit"/> with those properties set.
    /// </summary>
    public static LinearModelFit EBayes(LinearModelFit fit)
    {
        ArgumentNullException.ThrowIfNull(fit);
        int n = fit.Sigma.Length;
        var s2 = new double[n];
        for (int i = 0; i < n; i++) s2[i] = fit.Sigma[i] * fit.Sigma[i];
        var (d0, s20) = FitFDist(s2, fit.DfResidual);
        fit.DfPrior = d0; fit.S2Prior = s20;
        var post = new double[n];
        for (int i = 0; i < n; i++)
        {
            double d = fit.DfResidual[i];
            post[i] = double.IsPositiveInfinity(d0) ? s20 : (d * s2[i] + d0 * s20) / (d + d0);
        }
        fit.S2Post = post;
        return fit;
    }

    /// <summary>
    /// limma's <c>fitFDist</c> for a scaled F (s² ~ s₀² F(d, d₀)) without covariate: method of moments on
    /// e = log s² − ψ(d/2) + log(d/2). Returns (d₀, s₀²); d₀ is +∞ when the variances show no extra spread.
    /// </summary>
    public static (double DfPrior, double S2Prior) FitFDist(ReadOnlySpan<double> variances, ReadOnlySpan<double> df)
    {
        int n = variances.Length;
        if (df.Length != n) throw new ArgumentException("Need one df per variance.", nameof(df));
        // limma: the median (used to floor tiny variances) is taken over the usable values only, after pmax(x, 0).
        var usable = new List<double>(n);
        for (int i = 0; i < n; i++)
        {
            double x = variances[i], d = df[i];
            if (double.IsFinite(x) && d > 1e-15 && x > -1e-15) usable.Add(Math.Max(x, 0));
        }
        var med = NumericUtil.Median(usable.ToArray());
        if (!(med > 0)) med = 1;
        var e = new List<double>(n);
        double triMean = 0;
        for (int i = 0; i < n; i++)
        {
            double x = variances[i], d = df[i];
            if (!double.IsFinite(x) || !(d > 1e-15) || x < -1e-15) continue;
            x = Math.Max(x, 1e-5 * med);
            e.Add(Math.Log(x) - NumericUtil.Digamma(d / 2) + Math.Log(d / 2));
            triMean += NumericUtil.Trigamma(d / 2);
        }
        int k = e.Count;
        if (k == 0) return (double.NaN, double.NaN);
        triMean /= k;
        double emean = 0;
        foreach (double v in e) emean += v;
        emean /= k;
        if (k == 1) return (0, usable[0]); // limma: a single variance is its own scale, with df2 = 0
        double evar = 0;
        foreach (double v in e) evar += (v - emean) * (v - emean);
        evar = evar / (k - 1) - triMean;
        if (evar > 0)
        {
            double d0 = 2 * NumericUtil.TrigammaInverse(evar);
            return (d0, Math.Exp(emean + NumericUtil.Digamma(d0 / 2) - Math.Log(d0 / 2)));
        }
        return (double.PositiveInfinity, Math.Exp(emean));
    }

    /// <summary>Moderated t-test of one coefficient after <see cref="EBayes"/> (limma <c>topTable</c>, unsorted).</summary>
    public static DifferentialExpressionResult TopTable(LinearModelFit fit, int coefficient, string method = "limma")
    {
        ArgumentNullException.ThrowIfNull(fit);
        if (fit.S2Post.Length == 0) EBayes(fit);
        int n = fit.Sigma.Length, p = fit.Design.Coefficients;
        double dfPooled = 0;
        for (int i = 0; i < n; i++) if (double.IsFinite(fit.DfResidual[i])) dfPooled += fit.DfResidual[i];
        var lfc = new double[n]; var se = new double[n]; var t = new double[n]; var pv = new double[n];
        for (int i = 0; i < n; i++)
        {
            lfc[i] = fit.Coefficients[i * p + coefficient];
            se[i] = fit.StdevUnscaled[i * p + coefficient] * Math.Sqrt(fit.S2Post[i]);
            t[i] = lfc[i] / se[i];
            double dfTotal = Math.Min(fit.DfResidual[i] + fit.DfPrior, dfPooled);
            pv[i] = SpecialFunctions.StudentTTwoSided(t[i], dfTotal);
        }
        return new DifferentialExpressionResult(method, fit.Design.CoefficientNames[coefficient], fit.Genes,
            (double[])fit.Amean.Clone(), lfc, se, t, pv, DESeq.Bh(pv));
    }

    /// <summary>
    /// limma <c>voom</c>: log-CPM with offset 0.5, a lowess fit (span 0.5) of √σ against mean log count, and weights
    /// 1 / trend(fitted log count)⁴. Library sizes are multiplied by <paramref name="normFactors"/> (e.g. TMM).
    /// </summary>
    public static VoomResult Voom(Dataset counts, DesignMatrix design, ReadOnlySpan<double> normFactors = default, double span = 0.5)
    {
        ArgumentNullException.ThrowIfNull(counts);
        int n = counts.Rows, m = counts.Columns, p = design.Coefficients;
        var lib = Normalization.EffectiveLibrarySizes(counts, normFactors);
        var y = Normalization.CopyShape(counts, "voom");
        for (int i = 0; i < n; i++)
        {
            var src = counts.Row(i); var dst = y.Row(i);
            for (int j = 0; j < m; j++) dst[j] = Math.Log2((src[j] + 0.5) / (lib[j] + 1) * 1e6);
        }
        var fit = LmFit(y, design);
        double meanLogLib = 0;
        foreach (double l in lib) meanLogLib += Math.Log2(l + 1);
        meanLogLib /= m;
        var sx = new double[n]; var sy = new double[n];
        var xs = new List<double>(n); var ys = new List<double>(n);
        for (int i = 0; i < n; i++)
        {
            sx[i] = fit.Amean[i] + meanLogLib - Math.Log2(1e6);
            sy[i] = Math.Sqrt(fit.Sigma[i]);
            double total = 0;
            foreach (double c in counts.Row(i)) total += c;
            if (total > 0 && double.IsFinite(sy[i])) { xs.Add(sx[i]); ys.Add(sy[i]); }
        }
        // lowess on x-sorted data, then approxfun(rule = 2) with tied x averaged.
        // Stable sort like R's order(), so tied x keep their original order.
        var order = Enumerable.Range(0, xs.Count).ToArray();
        Array.Sort(order, (a, b) => { int c = xs[a].CompareTo(xs[b]); return c != 0 ? c : a.CompareTo(b); });
        var xk = new double[order.Length];
        for (int k = 0; k < order.Length; k++) xk[k] = xs[order[k]];
        var yk = new double[order.Length];
        for (int k = 0; k < order.Length; k++) yk[k] = ys[order[k]];
        var lo = NumericUtil.Lowess(xk, yk, span, 3);
        var (ux, uy) = CollapseTies(xk, lo);
        var w = new double[n * m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                double fv = 0;
                for (int a = 0; a < p; a++) fv += fit.Coefficients[i * p + a] * design[j, a];
                double fittedLogCount = Math.Log2(Math.Pow(2, fv) * 1e-6 * (lib[j] + 1));
                double s = NumericUtil.Interpolate(ux, uy, fittedLogCount);
                w[i * m + j] = 1 / Math.Pow(s, 4);
            }
        return new VoomResult(y, w, sx, sy);
    }

    /// <summary>
    /// limma-voom differential expression in one call: voom, weighted <c>lmFit</c>, <c>eBayes</c>, moderated t-test
    /// of <paramref name="coefficient"/> (default: last).
    /// </summary>
    public static DifferentialExpressionResult VoomTest(Dataset counts, DesignMatrix design, int coefficient = -1,
        ReadOnlySpan<double> normFactors = default)
    {
        if (coefficient < 0) coefficient = design.Coefficients - 1;
        var v = Voom(counts, design, normFactors);
        var fit = EBayes(LmFit(v.LogCpm, design, v.Weights));
        return TopTable(fit, coefficient, "limma-voom");
    }

    /// <summary>limma (trend-free) moderated t-test on log-scale expression (e.g. log-CPM or microarray data).</summary>
    public static DifferentialExpressionResult ModeratedTTest(Dataset logExpression, DesignMatrix design, int coefficient = -1)
    {
        if (coefficient < 0) coefficient = design.Coefficients - 1;
        return TopTable(EBayes(LmFit(logExpression, design)), coefficient);
    }

    private static (double[] X, double[] Y) CollapseTies(double[] x, double[] y)
    {
        var ox = new List<double>(); var oy = new List<double>();
        int s = 0;
        while (s < x.Length)
        {
            int e = s; double sum = y[s];
            while (e + 1 < x.Length && x[e + 1] == x[s]) { e++; sum += y[e]; }
            ox.Add(x[s]); oy.Add(sum / (e - s + 1));
            s = e + 1;
        }
        return (ox.ToArray(), oy.ToArray());
    }
}

/// <summary>Classical two-sample tests on log-scale expression.</summary>
public static class DifferentialExpression
{
    /// <summary>
    /// Welch (default) or Student two-sample t-test per gene on log-scale values (e.g. <see cref="Normalization.LogCpm"/>),
    /// comparing <paramref name="treatment"/> columns with <paramref name="reference"/> columns. Welch runs on the
    /// <see cref="BioKernels.WelchTTest"/> kernel. BaseMean holds the average log expression.
    /// </summary>
    public static DifferentialExpressionResult TTest(Dataset logExpression, int[] reference, int[] treatment, bool welch = true)
    {
        ArgumentNullException.ThrowIfNull(logExpression);
        int n = logExpression.Rows, m = logExpression.Columns;
        var t = new double[n]; var p = new double[n]; var diff = new double[n];
        var values = logExpression.Values.ToArray();
        if (welch) BioKernels.WelchTTest(values, n, m, reference, treatment, t, p, diff);
        else
        {
            if (reference.Length < 2 || treatment.Length < 2) throw new ArgumentException("Each group needs at least two samples.");
            Parallel.For(0, n, i =>
            {
                var row = values.AsSpan(i * m, m);
                MeanVar(row, reference, out double ma, out double va);
                MeanVar(row, treatment, out double mb, out double vb);
                int na = reference.Length, nb = treatment.Length;
                double sp = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2);
                double se = Math.Sqrt(sp * (1.0 / na + 1.0 / nb));
                diff[i] = mb - ma;
                t[i] = se > 0 ? diff[i] / se : double.NaN;
                p[i] = SpecialFunctions.StudentTTwoSided(t[i], na + nb - 2);
            });
        }
        var mean = new double[n]; var se2 = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = 0;
            foreach (int c in reference) s += values[i * m + c];
            foreach (int c in treatment) s += values[i * m + c];
            mean[i] = s / (reference.Length + treatment.Length);
            se2[i] = t[i] != 0 ? diff[i] / t[i] : double.NaN;
        }
        return new DifferentialExpressionResult(welch ? "Welch-t" : "Student-t", "treatment_vs_reference", logExpression.RowNames,
            mean, diff, se2, t, p, DESeq.Bh(p));
    }

    /// <summary>t-test using a column annotation to define groups.</summary>
    public static DifferentialExpressionResult TTest(Dataset logExpression, string factor, string reference, string treatment, bool welch = true)
        => TTest(logExpression, logExpression.ColumnsWhere(factor, reference), logExpression.ColumnsWhere(factor, treatment), welch);

    private static void MeanVar(ReadOnlySpan<double> row, int[] idx, out double mean, out double variance)
    {
        double s = 0;
        foreach (int i in idx) s += row[i];
        mean = s / idx.Length;
        double ss = 0;
        foreach (int i in idx) { double d = row[i] - mean; ss += d * d; }
        variance = ss / (idx.Length - 1);
    }
}
