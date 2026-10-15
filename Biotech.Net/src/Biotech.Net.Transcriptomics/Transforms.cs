using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>
/// Variance-stabilising transformations of counts for visualisation, clustering and PCA (not for testing).
/// </summary>
public static class Transforms
{
    /// <summary>
    /// DESeq2's closed-form variance-stabilising transformation for a dispersion trend α(μ) = a₀ + a₁/μ applied to
    /// normalised counts q = count / sizeFactor:
    /// <c>log2((1 + a₁ + 2a₀q + 2√(a₀q(1 + a₁ + a₀q))) / (4a₀))</c>. When <paramref name="extraPoisson"/> is 0 (mean
    /// trend) DESeq2's <c>(2 asinh(√(a₀q)) − log a₀ − log 4) / log 2</c> is used, which is the same function.
    /// </summary>
    public static Dataset Vst(Dataset counts, ReadOnlySpan<double> sizeFactors, double asymptoticDispersion, double extraPoisson)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (sizeFactors.Length != counts.Columns) throw new ArgumentException("Need one size factor per sample.");
        if (asymptoticDispersion <= 0) throw new ArgumentOutOfRangeException(nameof(asymptoticDispersion));
        double a0 = asymptoticDispersion, a1 = extraPoisson;
        var outp = Normalization.CopyShape(counts, "vst");
        int m = counts.Columns;
        double ln2 = Math.Log(2);
        for (int i = 0; i < counts.Rows; i++)
        {
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++)
            {
                double q = src[j] / sizeFactors[j];
                dst[j] = a1 == 0
                    ? (2 * Math.Asinh(Math.Sqrt(a0 * q)) - Math.Log(a0) - Math.Log(4)) / ln2
                    : Math.Log((1 + a1 + 2 * a0 * q + 2 * Math.Sqrt(a0 * q * (1 + a1 + a0 * q))) / (4 * a0)) / ln2;
            }
        }
        return outp;
    }

    /// <summary>
    /// Fits the DESeq2 dispersion trend (intercept-only design, so it is "blind" to conditions like DESeq2's
    /// <c>vst(blind = TRUE)</c>) and applies <see cref="Vst(Dataset, ReadOnlySpan{double}, double, double)"/>.
    /// </summary>
    public static Dataset Vst(Dataset counts)
    {
        var fit = DESeq.Fit(counts, DesignMatrix.Intercept(counts.Columns), new DESeqOptions { ShrinkDispersions = false });
        return Vst(counts, fit.SizeFactors, fit.Trend.AsymptoticDispersion, fit.Trend.ExtraPoisson);
    }

    /// <summary>
    /// An approximation of DESeq2's regularised log (rlog). <b>Not the exact rlog</b>: DESeq2 fits a ridge-penalised
    /// GLM per gene; here each sample's deviation of log2(q + 0.5) from the gene mean is shrunk by the
    /// normal–normal posterior factor <c>σ²_β / (σ²_β + v)</c>, where v = (1/μ + α)/ln²2 is the sampling variance of
    /// the log count under the fitted dispersion trend and σ²_β is matched to the 95% quantile of observed deviations
    /// (DESeq2's <c>matchUpperQuantileForVariance</c>). Low counts are therefore pulled to the gene mean, high counts
    /// keep their deviations, which reproduces rlog's qualitative behaviour.
    /// </summary>
    public static Dataset RlogApproximation(Dataset counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        var fit = DESeq.Fit(counts, DesignMatrix.Intercept(counts.Columns), new DESeqOptions { ShrinkDispersions = false });
        var sf = fit.SizeFactors;
        int n = counts.Rows, m = counts.Columns;
        var logq = new double[n * m];
        var means = new double[n];
        var devs = new List<double>(n * m);
        for (int i = 0; i < n; i++)
        {
            var row = counts.Row(i);
            double s = 0;
            for (int j = 0; j < m; j++) { double v = Math.Log2(row[j] / sf[j] + 0.5); logq[i * m + j] = v; s += v; }
            means[i] = s / m;
            if (fit.BaseMean[i] > 0)
                for (int j = 0; j < m; j++) devs.Add(Math.Abs(logq[i * m + j] - means[i]));
        }
        double q95 = devs.Count > 0 ? NumericUtil.Quantile(devs.ToArray(), 0.95) : 1;
        double priorVar = Math.Pow(q95 / SpecialFunctions.NormalQuantile(1 - 0.05 / 2), 2);
        priorVar = Math.Max(priorVar, 1e-8);
        double ln2sq = Math.Log(2) * Math.Log(2);
        var outp = Normalization.CopyShape(counts, "rlog");
        for (int i = 0; i < n; i++)
        {
            double alpha = double.IsFinite(fit.TrendDispersions[i]) ? fit.TrendDispersions[i] : fit.Trend.AsymptoticDispersion;
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++)
            {
                double mu = Math.Max(fit.BaseMean[i] * sf[j], 0.5);
                double v = (1 / mu + alpha) / ln2sq;
                double shrink = priorVar / (priorVar + v);
                dst[j] = means[i] + shrink * (logq[i * m + j] - means[i]);
            }
        }
        return outp;
    }
}

/// <summary>Filters of lowly expressed genes.</summary>
public static class ExpressionFilter
{
    /// <summary>
    /// edgeR's <c>filterByExpr</c>: keeps genes with CPM ≥ <c>minCount / median(libSize) × 10⁶</c> in at least
    /// <c>n</c> samples (n = smallest group size, relaxed to <c>largeN + (n − largeN) × minProp</c> for large groups)
    /// and a total count of at least <paramref name="minTotalCount"/>. Returns a keep mask.
    /// </summary>
    /// <param name="counts">Genes × samples counts.</param>
    /// <param name="groups">Group label per sample; null means all samples form one group.</param>
    /// <param name="minCount">Minimum count (converted to a CPM cutoff at the median library size).</param>
    /// <param name="minTotalCount">Minimum total count of a gene across all samples.</param>
    /// <param name="largeN">Group size above which the required number of samples is relaxed.</param>
    /// <param name="minProp">Proportion of samples beyond <paramref name="largeN"/> that must pass in large groups.</param>
    /// <param name="normFactors">Optional normalisation factors applied to the library sizes.</param>
    public static bool[] FilterByExpr(Dataset counts, IReadOnlyList<string>? groups = null, double minCount = 10,
        double minTotalCount = 15, int largeN = 10, double minProp = 0.7, ReadOnlySpan<double> normFactors = default)
    {
        ArgumentNullException.ThrowIfNull(counts);
        int m = counts.Columns;
        double minSample;
        if (groups is null) minSample = m;
        else
        {
            if (groups.Count != m) throw new ArgumentException("Need one group label per sample.");
            var sizes = new Dictionary<string, int>();
            foreach (var g in groups) sizes[g] = sizes.GetValueOrDefault(g) + 1;
            minSample = sizes.Values.Min();
        }
        return FilterCore(counts, minSample, minCount, minTotalCount, largeN, minProp, normFactors);
    }

    /// <summary>
    /// <c>filterByExpr</c> with a design matrix: the minimum sample size is 1 / max leverage of the design.
    /// </summary>
    public static bool[] FilterByExpr(Dataset counts, DesignMatrix design, double minCount = 10, double minTotalCount = 15,
        int largeN = 10, double minProp = 0.7)
    {
        ArgumentNullException.ThrowIfNull(design);
        double maxH = design.Leverages().Max();
        return FilterCore(counts, 1 / maxH, minCount, minTotalCount, largeN, minProp, default);
    }

    /// <summary>Keeps genes with at least <paramref name="minCount"/> counts in at least <paramref name="minSamples"/> samples.</summary>
    public static bool[] MinimumCount(Dataset counts, double minCount = 10, int minSamples = 1)
    {
        ArgumentNullException.ThrowIfNull(counts);
        var keep = new bool[counts.Rows];
        for (int i = 0; i < counts.Rows; i++)
        {
            int c = 0;
            foreach (double x in counts.Row(i)) if (x >= minCount) c++;
            keep[i] = c >= minSamples;
        }
        return keep;
    }

    /// <summary>A dataset with only the rows where <paramref name="keep"/> is true.</summary>
    public static Dataset Apply(Dataset data, ReadOnlySpan<bool> keep)
    {
        ArgumentNullException.ThrowIfNull(data);
        var rows = new List<int>();
        for (int i = 0; i < keep.Length; i++) if (keep[i]) rows.Add(i);
        return data.SelectRows(rows);
    }

    private static bool[] FilterCore(Dataset counts, double minSample, double minCount, double minTotalCount, int largeN,
        double minProp, ReadOnlySpan<double> normFactors)
    {
        if (minSample > largeN) minSample = largeN + (minSample - largeN) * minProp;
        var lib = Normalization.EffectiveLibrarySizes(counts, normFactors);
        double cutoff = minCount / NumericUtil.Median(lib) * 1e6;
        const double tol = 1e-14;
        int m = counts.Columns;
        var keep = new bool[counts.Rows];
        for (int i = 0; i < counts.Rows; i++)
        {
            var row = counts.Row(i);
            int above = 0; double total = 0;
            for (int j = 0; j < m; j++)
            {
                total += row[j];
                if (row[j] / lib[j] * 1e6 >= cutoff) above++;
            }
            keep[i] = above >= minSample - tol && total >= minTotalCount - tol;
        }
        return keep;
    }
}
