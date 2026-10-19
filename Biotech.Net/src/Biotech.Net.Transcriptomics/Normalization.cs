using Biotech.Net.Core;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>Between-sample normalisation methods for count data.</summary>
public enum NormalizationMethod
{
    /// <summary>Raw library size (total counts) only.</summary>
    None,

    /// <summary>edgeR's trimmed mean of M-values (Robinson &amp; Oshlack 2010).</summary>
    Tmm,

    /// <summary>Upper-quartile normalisation (Bullard et al. 2010).</summary>
    UpperQuartile,

    /// <summary>DESeq2's median-of-ratios size factors (Anders &amp; Huber 2010).</summary>
    MedianOfRatios,
}

/// <summary>
/// Library-size and composition normalisation of RNA-Seq counts (genes × samples <see cref="Dataset"/>),
/// plus the gene-length-aware units CPM, TPM and FPKM/RPKM. Formulas follow edgeR and DESeq2 exactly so
/// results can be compared with R.
/// </summary>
public static class Normalization
{
    /// <summary>Total counts per sample.</summary>
    public static double[] LibrarySizes(Dataset counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        return counts.ColumnSums();
    }

    /// <summary>
    /// DESeq2 size factors (<c>estimateSizeFactors</c>, median-of-ratios): for each sample, the median over genes
    /// with no zero count of <c>count / geometric mean of the gene across samples</c>.
    /// </summary>
    /// <exception cref="InvalidOperationException">Every gene contains at least one zero.</exception>
    public static double[] MedianOfRatiosSizeFactors(Dataset counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        int n = counts.Rows, m = counts.Columns;
        var logGeo = new double[n];
        int usable = 0;
        for (int i = 0; i < n; i++)
        {
            var row = counts.Row(i);
            double s = 0;
            bool ok = true;
            foreach (double x in row) { if (x <= 0) { ok = false; break; } s += Math.Log(x); }
            logGeo[i] = ok ? s / m : double.NegativeInfinity;
            if (ok) usable++;
        }
        if (usable == 0) throw new InvalidOperationException("Every gene contains at least one zero; median-of-ratios size factors are undefined.");
        var sf = new double[m];
        var ratios = new double[usable];
        for (int j = 0; j < m; j++)
        {
            int k = 0;
            for (int i = 0; i < n; i++)
                if (!double.IsNegativeInfinity(logGeo[i])) ratios[k++] = Math.Log(counts[i, j]) - logGeo[i];
            sf[j] = Math.Exp(NumericUtil.MedianInPlace(ratios.AsSpan(0, k)));
        }
        return sf;
    }

    /// <summary>
    /// DESeq2 <c>estimateSizeFactors(type = "poscounts")</c>: median-of-ratios against the geometric mean of each gene's
    /// <em>positive</em> counts (divided by the full number of samples), using only positive counts in each sample, then
    /// scaled so the size factors have geometric mean 1. Works when every gene has at least one zero.
    /// </summary>
    public static double[] PosCountsSizeFactors(Dataset counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        int n = counts.Rows, m = counts.Columns;
        var logGeo = new double[n];
        bool any = false;
        for (int i = 0; i < n; i++)
        {
            double s = 0; bool pos = false;
            foreach (double x in counts.Row(i)) if (x > 0) { s += Math.Log(x); pos = true; }
            // geoMeans = exp(sum(log(x[x > 0])) / length(x)); log() of it, -Inf for all-zero genes.
            logGeo[i] = pos ? Math.Log(Math.Exp(s / m)) : double.NegativeInfinity;
            any |= pos;
        }
        if (!any) throw new InvalidOperationException("Every gene is all zero; size factors are undefined.");
        var sf = new double[m];
        var ratios = new double[n];
        for (int j = 0; j < m; j++)
        {
            int k = 0;
            for (int i = 0; i < n; i++)
                if (double.IsFinite(logGeo[i]) && counts[i, j] > 0) ratios[k++] = Math.Log(counts[i, j]) - logGeo[i];
            sf[j] = Math.Exp(NumericUtil.MedianInPlace(ratios.AsSpan(0, k)));
        }
        double ml = 0;
        foreach (double v in sf) ml += Math.Log(v);
        double g = Math.Exp(ml / m);
        for (int j = 0; j < m; j++) sf[j] /= g;
        return sf;
    }

    /// <summary>
    /// edgeR TMM normalisation factors (<c>calcNormFactors(method = "TMM")</c>): weighted trimmed mean of the log
    /// expression ratios against a reference sample, trimming 30% of M-values and 5% of A-values. The reference is
    /// the sample whose upper-quartile factor is closest to the mean, as in edgeR. Factors are scaled to have a
    /// geometric mean of 1; effective library sizes are <c>libSize × factor</c>.
    /// </summary>
    public static double[] TmmFactors(Dataset counts, int? referenceColumn = null, double logRatioTrim = 0.3,
        double sumTrim = 0.05, bool weighting = true, double aCutoff = -1e10)
    {
        ArgumentNullException.ThrowIfNull(counts);
        var x = NonZeroRows(counts);
        int m = counts.Columns;
        var lib = LibraryFromRows(x, m);
        int refCol;
        if (referenceColumn is int rc) refCol = rc;
        else
        {
            var f75 = QuantileFactorsRaw(x, m, lib, 0.75);
            double mean = 0;
            foreach (double v in f75) mean += v;
            mean /= m;
            if (NumericUtil.Median(f75) < 1e-20)
            {
                // edgeR: fall back to the column with the largest sum of square-root counts.
                refCol = 0; double best = double.NegativeInfinity;
                for (int j = 0; j < m; j++)
                {
                    double s = 0;
                    foreach (var r in x) s += Math.Sqrt(r[j]);
                    if (s > best) { best = s; refCol = j; }
                }
            }
            else
            {
                refCol = 0; double best = double.PositiveInfinity;
                for (int j = 0; j < m; j++) { double d = Math.Abs(f75[j] - mean); if (d < best) { best = d; refCol = j; } }
            }
        }
        var f = new double[m];
        for (int j = 0; j < m; j++)
            f[j] = TmmFactor(x, j, refCol, lib[j], lib[refCol], logRatioTrim, sumTrim, weighting, aCutoff);
        return ScaleToUnitGeometricMean(f);
    }

    /// <summary>
    /// Upper-quartile factors (<c>calcNormFactors(method = "upperquartile")</c>): the <paramref name="p"/> quantile
    /// of each sample's count proportions (all-zero genes removed), scaled to geometric mean 1.
    /// </summary>
    public static double[] UpperQuartileFactors(Dataset counts, double p = 0.75)
    {
        ArgumentNullException.ThrowIfNull(counts);
        var x = NonZeroRows(counts);
        var lib = LibraryFromRows(x, counts.Columns);
        return ScaleToUnitGeometricMean(QuantileFactorsRaw(x, counts.Columns, lib, p));
    }

    /// <summary>
    /// Normalisation factors for <paramref name="method"/>, expressed edgeR-style (multipliers of library size).
    /// For <see cref="NormalizationMethod.MedianOfRatios"/> the size factors are converted as <c>sf / libSize</c>,
    /// rescaled to geometric mean 1, so <c>libSize × factor ∝ sizeFactor</c>.
    /// </summary>
    public static double[] Factors(Dataset counts, NormalizationMethod method)
    {
        ArgumentNullException.ThrowIfNull(counts);
        switch (method)
        {
            case NormalizationMethod.Tmm: return TmmFactors(counts);
            case NormalizationMethod.UpperQuartile: return UpperQuartileFactors(counts);
            case NormalizationMethod.MedianOfRatios:
            {
                var sf = MedianOfRatiosSizeFactors(counts);
                var lib = LibrarySizes(counts);
                for (int j = 0; j < sf.Length; j++) sf[j] /= lib[j];
                return ScaleToUnitGeometricMean(sf);
            }
            default:
            {
                var f = new double[counts.Columns];
                Array.Fill(f, 1.0);
                return f;
            }
        }
    }

    /// <summary>Effective library sizes: total counts × normalisation factors (factors default to 1).</summary>
    public static double[] EffectiveLibrarySizes(Dataset counts, ReadOnlySpan<double> normFactors = default)
    {
        var lib = LibrarySizes(counts);
        if (!normFactors.IsEmpty)
        {
            RequireColumns(counts, normFactors.Length);
            for (int j = 0; j < lib.Length; j++) lib[j] *= normFactors[j];
        }
        return lib;
    }

    /// <summary>Counts per million: <c>count / (libSize × factor) × 10⁶</c> (edgeR <c>cpm</c>).</summary>
    public static Dataset Cpm(Dataset counts, ReadOnlySpan<double> normFactors = default)
    {
        var lib = EffectiveLibrarySizes(counts, normFactors);
        var outp = CopyShape(counts, "cpm");
        int m = counts.Columns;
        for (int i = 0; i < counts.Rows; i++)
        {
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) dst[j] = src[j] / lib[j] * 1e6;
        }
        return outp;
    }

    /// <summary>
    /// log2 counts per million with edgeR's library-size-scaled prior count (<c>cpm(y, log = TRUE, prior.count)</c>):
    /// the prior for sample j is <c>prior × lib_j / mean(lib)</c> and the library is enlarged by twice that prior.
    /// </summary>
    public static Dataset LogCpm(Dataset counts, double priorCount = 2, ReadOnlySpan<double> normFactors = default)
    {
        var lib = EffectiveLibrarySizes(counts, normFactors);
        int m = counts.Columns;
        double meanLib = 0;
        foreach (double l in lib) meanLib += l;
        meanLib /= m;
        var prior = new double[m];
        var adj = new double[m];
        for (int j = 0; j < m; j++)
        {
            prior[j] = priorCount * lib[j] / meanLib;
            adj[j] = lib[j] + 2 * prior[j];
        }
        var outp = CopyShape(counts, "logcpm");
        for (int i = 0; i < counts.Rows; i++)
        {
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) dst[j] = Math.Log2((src[j] + prior[j]) / adj[j] * 1e6);
        }
        return outp;
    }

    /// <summary>
    /// Transcripts per million: reads per kilobase, rescaled so each sample sums to 10⁶. Gene lengths are in base
    /// pairs (use featureCounts' <c>Length</c> or the effective length from a quantifier).
    /// </summary>
    public static Dataset Tpm(Dataset counts, ReadOnlySpan<double> geneLengths)
    {
        ArgumentNullException.ThrowIfNull(counts);
        RequireRows(counts, geneLengths.Length);
        int m = counts.Columns;
        var outp = CopyShape(counts, "tpm");
        var totals = new double[m];
        for (int i = 0; i < counts.Rows; i++)
        {
            double kb = geneLengths[i] / 1000.0;
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) { double r = kb > 0 ? src[j] / kb : 0; dst[j] = r; totals[j] += r; }
        }
        for (int i = 0; i < counts.Rows; i++)
        {
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) dst[j] = totals[j] > 0 ? dst[j] / totals[j] * 1e6 : 0;
        }
        return outp;
    }

    /// <summary>TPM using the <c>length</c> row annotation (see <see cref="CountMatrixReader.GeneLengths"/>).</summary>
    public static Dataset Tpm(Dataset counts) => Tpm(counts, CountMatrixReader.GeneLengths(counts));

    /// <summary>
    /// Fragments (or reads) per kilobase per million mapped: <c>count / (length/10³) / (libSize × factor / 10⁶)</c>,
    /// i.e. edgeR's <c>rpkm</c>. FPKM and RPKM are the same formula applied to fragments or reads.
    /// </summary>
    public static Dataset Fpkm(Dataset counts, ReadOnlySpan<double> geneLengths, ReadOnlySpan<double> normFactors = default)
    {
        RequireRows(counts, geneLengths.Length);
        var lib = EffectiveLibrarySizes(counts, normFactors);
        int m = counts.Columns;
        var outp = CopyShape(counts, "fpkm");
        for (int i = 0; i < counts.Rows; i++)
        {
            double kb = geneLengths[i] / 1000.0;
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) dst[j] = kb > 0 ? src[j] / kb / (lib[j] / 1e6) : 0;
        }
        return outp;
    }

    /// <summary>Alias of <see cref="Fpkm(Dataset, ReadOnlySpan{double}, ReadOnlySpan{double})"/> for single-end reads.</summary>
    public static Dataset Rpkm(Dataset counts, ReadOnlySpan<double> geneLengths, ReadOnlySpan<double> normFactors = default)
        => Fpkm(counts, geneLengths, normFactors);

    /// <summary>DESeq2 normalised counts: each column divided by its size factor (<c>counts(dds, normalized = TRUE)</c>).</summary>
    public static Dataset NormalizedCounts(Dataset counts, ReadOnlySpan<double> sizeFactors)
    {
        RequireColumns(counts, sizeFactors.Length);
        var outp = CopyShape(counts, "normalized");
        int m = counts.Columns;
        for (int i = 0; i < counts.Rows; i++)
        {
            var src = counts.Row(i);
            var dst = outp.Row(i);
            for (int j = 0; j < m; j++) dst[j] = src[j] / sizeFactors[j];
        }
        return outp;
    }

    internal static Dataset CopyShape(Dataset src, string suffix)
    {
        var d = new Dataset($"{src.Name}.{suffix}", src.RowNames, src.ColumnNames, new double[src.Rows * src.Columns]);
        foreach (var kv in src.ColumnAnnotations) d.ColumnAnnotations[kv.Key] = kv.Value;
        foreach (var kv in src.RowAnnotations) d.RowAnnotations[kv.Key] = kv.Value;
        return d;
    }

    private static double TmmFactor(List<double[]> x, int obsCol, int refCol, double nO, double nR,
        double logRatioTrim, double sumTrim, bool weighting, double aCutoff)
    {
        var logR = new List<double>(x.Count);
        var absE = new List<double>(x.Count);
        var v = new List<double>(x.Count);
        foreach (var row in x)
        {
            double obs = row[obsCol], rf = row[refCol];
            double lr = Math.Log2(obs / nO / (rf / nR));
            double ae = (Math.Log2(obs / nO) + Math.Log2(rf / nR)) / 2;
            if (!double.IsFinite(lr) || !double.IsFinite(ae) || ae <= aCutoff) continue;
            logR.Add(lr); absE.Add(ae);
            v.Add((nO - obs) / nO / obs + (nR - rf) / nR / rf);
        }
        int n = logR.Count;
        if (n == 0) return 1;
        double maxAbs = 0;
        foreach (double l in logR) maxAbs = Math.Max(maxAbs, Math.Abs(l));
        if (maxAbs < 1e-6) return 1;
        int loL = (int)Math.Floor(n * logRatioTrim) + 1, hiL = n + 1 - loL;
        int loS = (int)Math.Floor(n * sumTrim) + 1, hiS = n + 1 - loS;
        var rankR = NumericUtil.Rank(logR.ToArray());
        var rankE = NumericUtil.Rank(absE.ToArray());
        double num = 0, den = 0, sum = 0;
        int cnt = 0;
        for (int i = 0; i < n; i++)
        {
            if (rankR[i] < loL || rankR[i] > hiL || rankE[i] < loS || rankE[i] > hiS) continue;
            num += logR[i] / v[i]; den += 1 / v[i]; sum += logR[i]; cnt++;
        }
        double f = weighting ? num / den : sum / cnt;
        if (double.IsNaN(f)) f = 0;
        return Math.Pow(2, f);
    }

    private static List<double[]> NonZeroRows(Dataset counts)
    {
        var rows = new List<double[]>(counts.Rows);
        for (int i = 0; i < counts.Rows; i++)
        {
            var r = counts.Row(i);
            bool any = false;
            foreach (double v in r) if (v > 0) { any = true; break; }
            if (any) rows.Add(r.ToArray());
        }
        return rows;
    }

    private static double[] LibraryFromRows(List<double[]> x, int m)
    {
        var lib = new double[m];
        foreach (var r in x) for (int j = 0; j < m; j++) lib[j] += r[j];
        return lib;
    }

    private static double[] QuantileFactorsRaw(List<double[]> x, int m, double[] lib, double p)
    {
        var f = new double[m];
        var col = new double[x.Count];
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < x.Count; i++) col[i] = x[i][j] / lib[j];
            Array.Sort(col);
            f[j] = NumericUtil.QuantileSorted(col, p);
        }
        return f;
    }

    private static double[] ScaleToUnitGeometricMean(double[] f)
    {
        double s = 0;
        foreach (double v in f) s += Math.Log(v);
        double g = Math.Exp(s / f.Length);
        for (int j = 0; j < f.Length; j++) f[j] /= g;
        return f;
    }

    private static void RequireRows(Dataset d, int n)
    {
        ArgumentNullException.ThrowIfNull(d);
        if (n != d.Rows) throw new ArgumentException($"Expected {d.Rows} per-gene values, got {n}.");
    }

    private static void RequireColumns(Dataset d, int n)
    {
        ArgumentNullException.ThrowIfNull(d);
        if (n != d.Columns) throw new ArgumentException($"Expected {d.Columns} per-sample values, got {n}.");
    }
}
