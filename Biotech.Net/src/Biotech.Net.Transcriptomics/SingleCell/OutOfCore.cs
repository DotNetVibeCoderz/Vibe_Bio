using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// Per-cell and per-gene quality-control metrics, named and computed as scanpy's
/// <c>pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=(50, 100, 200, 500), log1p=True)</c>.
/// </summary>
public sealed class QcMetrics
{
    internal QcMetrics(int cells, int genes, int[] percentTop)
    {
        PercentTop = percentTop;
        NGenesByCounts = new double[cells]; TotalCounts = new double[cells];
        TotalCountsMito = new double[cells]; PctCountsMito = new double[cells];
        PctCountsInTop = percentTop.Select(_ => new double[cells]).ToArray();
        NCellsByCounts = new double[genes]; GeneTotalCounts = new double[genes];
        MeanCounts = new double[genes]; PctDropoutByCounts = new double[genes]; Mito = new bool[genes];
    }

    /// <summary>Genes with a non-zero count, per cell (<c>n_genes_by_counts</c>).</summary>
    public double[] NGenesByCounts { get; }

    /// <summary>Total counts per cell (<c>total_counts</c>).</summary>
    public double[] TotalCounts { get; }

    /// <summary>Counts in mitochondrial genes per cell (<c>total_counts_mt</c>).</summary>
    public double[] TotalCountsMito { get; }

    /// <summary>Percentage of counts in mitochondrial genes (<c>pct_counts_mt</c>, 0–100).</summary>
    public double[] PctCountsMito { get; }

    /// <summary>The <c>percent_top</c> values N.</summary>
    public IReadOnlyList<int> PercentTop { get; }

    /// <summary>Percentage of counts in each cell's N most expressed genes, one array per <see cref="PercentTop"/> entry.</summary>
    public double[][] PctCountsInTop { get; }

    /// <summary>Cells in which each gene is detected (<c>n_cells_by_counts</c>).</summary>
    public double[] NCellsByCounts { get; }

    /// <summary>Total counts per gene (<c>total_counts</c> in var).</summary>
    public double[] GeneTotalCounts { get; }

    /// <summary>Mean count per gene over all cells (<c>mean_counts</c>).</summary>
    public double[] MeanCounts { get; }

    /// <summary>Percentage of cells where the gene is not detected (<c>pct_dropout_by_counts</c>).</summary>
    public double[] PctDropoutByCounts { get; }

    /// <summary>Mitochondrial gene flag (<c>mt</c>).</summary>
    public bool[] Mito { get; }

    /// <summary>Adds the metrics as scanpy-named columns (including the <c>log1p_</c> columns).</summary>
    public void WriteTo(AnnDataFrame obs, AnnDataFrame var)
    {
        ArgumentNullException.ThrowIfNull(obs);
        ArgumentNullException.ThrowIfNull(var);
        static double[] L(double[] a) => a.Select(v => Math.Log(1 + v)).ToArray();
        obs.Set("n_genes_by_counts", NGenesByCounts); obs.Set("log1p_n_genes_by_counts", L(NGenesByCounts));
        obs.Set("total_counts", TotalCounts); obs.Set("log1p_total_counts", L(TotalCounts));
        for (int t = 0; t < PercentTop.Count; t++) obs.Set($"pct_counts_in_top_{PercentTop[t]}_genes", PctCountsInTop[t]);
        obs.Set("total_counts_mt", TotalCountsMito); obs.Set("log1p_total_counts_mt", L(TotalCountsMito)); obs.Set("pct_counts_mt", PctCountsMito);
        var.Set("mt", Mito);
        var.Set("n_cells_by_counts", NCellsByCounts); var.Set("mean_counts", MeanCounts); var.Set("log1p_mean_counts", L(MeanCounts));
        var.Set("pct_dropout_by_counts", PctDropoutByCounts); var.Set("total_counts", GeneTotalCounts); var.Set("log1p_total_counts", L(GeneTotalCounts));
    }

    /// <summary>Stores the metrics in a <see cref="SingleCellExperiment"/> under the keys <see cref="Preprocessing.CalculateQcMetrics"/> uses.</summary>
    public void WriteTo(SingleCellExperiment sc)
    {
        ArgumentNullException.ThrowIfNull(sc);
        sc.CellMetrics[Preprocessing.NGenes] = NGenesByCounts; sc.CellMetrics[Preprocessing.TotalCounts] = TotalCounts;
        sc.CellMetrics[Preprocessing.TotalMito] = TotalCountsMito; sc.CellMetrics[Preprocessing.PctMito] = PctCountsMito;
        sc.GeneFlags["mt"] = Mito;
        sc.GeneMetrics["n_cells_by_counts"] = NCellsByCounts; sc.GeneMetrics["total_counts"] = GeneTotalCounts; sc.GeneMetrics["mean_counts"] = MeanCounts;
    }
}

/// <summary>Highly-variable-gene statistics (scanpy's <c>highly_variable_genes(..., inplace=False)</c> frame).</summary>
public sealed class HvgResult
{
    internal HvgResult(HvgFlavor flavor, double[] means, double[] spread, double[] normalized, bool[] hv, double[]? rank)
    {
        Flavor = flavor; Means = means; Spread = spread; Normalized = normalized; HighlyVariable = hv; Rank = rank;
        Selected = Enumerable.Range(0, hv.Length).Where(j => hv[j]).ToArray();
    }

    /// <summary>The flavour used.</summary>
    public HvgFlavor Flavor { get; }

    /// <summary><c>means</c>: log1p of the mean of expm1(X) for Seurat, the mean count for Seurat v3.</summary>
    public double[] Means { get; }

    /// <summary><c>dispersions</c> (log variance/mean, Seurat) or <c>variances</c> (Seurat v3).</summary>
    public double[] Spread { get; }

    /// <summary><c>dispersions_norm</c> (Seurat) or <c>variances_norm</c> (Seurat v3).</summary>
    public double[] Normalized { get; }

    /// <summary><c>highly_variable</c>.</summary>
    public bool[] HighlyVariable { get; }

    /// <summary><c>highly_variable_rank</c> (Seurat v3 only; NaN outside the selection).</summary>
    public double[]? Rank { get; }

    /// <summary>Indices of the selected genes, ascending.</summary>
    public int[] Selected { get; }

    /// <summary>Adds the scanpy-named columns to a var frame.</summary>
    public void WriteTo(AnnDataFrame var)
    {
        ArgumentNullException.ThrowIfNull(var);
        var.Set("highly_variable", HighlyVariable);
        var.Set("means", Means);
        if (Flavor == HvgFlavor.Seurat) { var.Set("dispersions", Spread); var.Set("dispersions_norm", Normalized); }
        else { var.Set("variances", Spread); var.Set("variances_norm", Normalized); var.Set("highly_variable_rank", Rank!); }
    }

    /// <summary>Stores the result in a <see cref="SingleCellExperiment"/> (the keys <see cref="Preprocessing.HighlyVariableGenes"/> uses).</summary>
    public void WriteTo(SingleCellExperiment sc)
    {
        ArgumentNullException.ThrowIfNull(sc);
        sc.GeneFlags["highly_variable"] = HighlyVariable;
        sc.GeneMetrics["means"] = Means;
        if (Flavor == HvgFlavor.Seurat) { sc.GeneMetrics["dispersions"] = Spread; sc.GeneMetrics["dispersions_norm"] = Normalized; }
        else { sc.GeneMetrics["variances"] = Spread; sc.GeneMetrics["variances_norm"] = Normalized; }
    }
}

/// <summary>
/// Out-of-core single-cell preprocessing: QC metrics, library-size normalisation + log1p, highly variable genes,
/// each streaming once (or twice) over row blocks of an <see cref="ISparseMatrixSource"/>, so memory is bounded by
/// the block size plus O(cells + genes) result arrays whatever the matrix size. Results match scanpy
/// (<c>calculate_qc_metrics</c>, <c>normalize_total</c>, <c>log1p</c>, <c>highly_variable_genes</c> with
/// <c>flavor="seurat"</c> and <c>"seurat_v3"</c>); see <see cref="OutOfCorePca"/> for PCA.
/// </summary>
public static class OutOfCore
{
    private static readonly int[] DefaultPercentTop = [50, 100, 200, 500];

    /// <summary>
    /// QC metrics in one pass. Mitochondrial genes are those whose name starts with <paramref name="mitoPrefix"/>
    /// (case-insensitive). <paramref name="percentTop"/> values larger than the gene count are dropped.
    /// </summary>
    public static QcMetrics CalculateQcMetrics(ISparseMatrixSource counts, IReadOnlyList<string> geneNames, string mitoPrefix = "MT-",
        IReadOnlyList<int>? percentTop = null, int blockRows = SparseMatrixSource.DefaultBlockRows, int? maxParallelism = null)
    {
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentNullException.ThrowIfNull(geneNames);
        if (geneNames.Count != counts.Columns) throw new ArgumentException($"Need {counts.Columns} gene names, got {geneNames.Count}.");
        int n = counts.Rows, g = counts.Columns;
        var tops = (percentTop ?? DefaultPercentTop).Where(t => t >= 1 && t <= g).Order().ToArray();
        var q = new QcMetrics(n, g, tops);
        for (int j = 0; j < g; j++) q.Mito[j] = geneNames[j].StartsWith(mitoPrefix, StringComparison.OrdinalIgnoreCase);
        var nCells = new long[g]; var gSum = new double[g];
        var sync = new object();
        SparseMatrixSource.ForEachBlock(counts, blockRows, (start, b) =>
        {
            var lc = new long[g]; var ls = new double[g];
            var buf = new List<float>();
            for (int i = 0; i < b.Rows; i++)
            {
                var idx = b.RowIndices(i); var v = b.RowValues(i);
                double t = 0, m = 0; int det = 0;
                buf.Clear();
                for (int k = 0; k < idx.Length; k++)
                {
                    float x = v[k];
                    if (x == 0) continue;
                    det++; t += x; buf.Add(x);
                    if (q.Mito[idx[k]]) m += x;
                    lc[idx[k]]++; ls[idx[k]] += x;
                }
                int row = start + i;
                q.NGenesByCounts[row] = det; q.TotalCounts[row] = t; q.TotalCountsMito[row] = m;
                q.PctCountsMito[row] = t > 0 ? 100 * m / t : double.NaN;
                if (tops.Length > 0)
                {
                    var span = System.Runtime.InteropServices.CollectionsMarshal.AsSpan(buf);
                    span.Sort((a, c) => c.CompareTo(a));
                    double acc = 0; int at = 0;
                    for (int ti = 0; ti < tops.Length; ti++)
                    {
                        int lim = Math.Min(tops[ti], span.Length);
                        for (; at < lim; at++) acc += span[at];
                        q.PctCountsInTop[ti][row] = t > 0 ? 100 * acc / t : double.NaN;
                    }
                }
            }
            lock (sync) for (int j = 0; j < g; j++) { nCells[j] += lc[j]; gSum[j] += ls[j]; }
        }, maxParallelism);
        for (int j = 0; j < g; j++)
        {
            q.NCellsByCounts[j] = nCells[j];
            q.GeneTotalCounts[j] = gSum[j];
            q.MeanCounts[j] = n > 0 ? gSum[j] / n : 0;
            q.PctDropoutByCounts[j] = n > 0 ? (1 - (double)nCells[j] / n) * 100 : 0;
        }
        return q;
    }

    /// <summary>
    /// Per-gene mean and variance (ddof = 1, implicit zeros included) in one pass, optionally of expm1(X).
    /// </summary>
    public static (double[] Mean, double[] Variance) GeneMeanVariance(ISparseMatrixSource x, bool expm1 = false,
        int blockRows = SparseMatrixSource.DefaultBlockRows, int? maxParallelism = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        var (sum, sq) = ColumnMoments(x, expm1, null, blockRows, maxParallelism);
        return MeanVar(sum, sq, x.Rows);
    }

    private static (double[] Mean, double[] Variance) MeanVar(double[] sum, double[] sq, int n)
    {
        int g = sum.Length;
        var mean = new double[g]; var var = new double[g];
        for (int j = 0; j < g; j++)
        {
            mean[j] = sum[j] / n;
            var[j] = n > 1 ? (sq[j] / n - mean[j] * mean[j]) * n / (n - 1) : 0;
            if (var[j] < 0) var[j] = 0;
        }
        return (mean, var);
    }

    /// <summary>Column sums of x (or expm1(x), or min(x, clip_j)) and of squares, in double precision.</summary>
    private static (double[] Sum, double[] SumSq) ColumnMoments(ISparseMatrixSource x, bool expm1, double[]? clip, int blockRows, int? maxParallelism)
    {
        int g = x.Columns;
        var sum = new double[g]; var sq = new double[g];
        var sync = new object();
        SparseMatrixSource.ForEachBlock(x, blockRows, (_, b) =>
        {
            var ls = new double[g]; var lq = new double[g];
            int nnz = b.NonZeroCount;
            var ci = b.ColumnIndices; var vs = b.Values;
            for (int k = 0; k < nnz; k++)
            {
                int j = ci[k];
                double v = expm1 ? (float)Expm1(vs[k]) : vs[k];
                if (clip is not null && v > clip[j]) v = clip[j];
                ls[j] += v; lq[j] += v * v;
            }
            lock (sync) for (int j = 0; j < g; j++) { sum[j] += ls[j]; sq[j] += lq[j]; }
        }, maxParallelism);
        return (sum, sq);
    }

    internal static double Expm1(double x) => Math.Abs(x) < 1e-5 ? x + x * x / 2 + x * x * x / 6 : Math.Exp(x) - 1;

    /// <summary>
    /// <c>normalize_total</c> (+ <c>log1p</c>) into a new <see cref="MappedSparseMatrix"/> store at
    /// <paramref name="outputDirectory"/>. Each cell is divided by total/<paramref name="targetSum"/> (median of the
    /// non-zero totals when null), in float32 like scanpy. Pass <paramref name="cellTotals"/> (e.g.
    /// <see cref="QcMetrics.TotalCounts"/>) to skip the extra pass that computes them.
    /// </summary>
    public static MappedSparseMatrix NormalizeTotalLog1p(ISparseMatrixSource counts, string outputDirectory, double? targetSum = 1e4, bool log1p = true,
        double[]? cellTotals = null, int blockRows = SparseMatrixSource.DefaultBlockRows)
    {
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentException.ThrowIfNullOrEmpty(outputDirectory);
        var factors = NormalizationFactors(counts, targetSum, cellTotals, blockRows);
        using var w = MappedSparseMatrixWriter.Create(outputDirectory, counts.Columns);
        foreach (var (start, b) in counts.RowBlocks(blockRows))
        {
            Normalize(start, b, factors, log1p);
            w.AppendRows(b);
        }
        return w.Complete();
    }

    /// <summary>
    /// <c>normalize_total</c> (+ <c>log1p</c>) written as a new layer of the same store (4 bytes per stored entry;
    /// the sparsity structure is shared). Returns the new layer.
    /// </summary>
    public static MappedSparseMatrix NormalizeTotalLog1pToLayer(MappedSparseMatrix counts, string layerName = "log1p_norm", double? targetSum = 1e4,
        bool log1p = true, double[]? cellTotals = null, int blockRows = SparseMatrixSource.DefaultBlockRows)
    {
        ArgumentNullException.ThrowIfNull(counts);
        var factors = NormalizationFactors(counts, targetSum, cellTotals, blockRows);
        return counts.WriteLayer(layerName, (start, b) => Normalize(start, b, factors, log1p), blockRows);
    }

    /// <summary>scanpy's divisors: float32(total) / target, as float32 (cells with zero total are left as they are).</summary>
    internal static float[] NormalizationFactors(ISparseMatrixSource counts, double? targetSum, double[]? cellTotals, int blockRows)
    {
        int n = counts.Rows;
        var totals = new float[n];
        if (cellTotals is not null)
        {
            if (cellTotals.Length != n) throw new ArgumentException($"Need {n} cell totals.", nameof(cellTotals));
            for (int i = 0; i < n; i++) totals[i] = (float)cellTotals[i];
        }
        else
        {
            SparseMatrixSource.ForEachBlock(counts, blockRows, (start, b) =>
            {
                for (int i = 0; i < b.Rows; i++) { double t = 0; foreach (float v in b.RowValues(i)) t += v; totals[start + i] = (float)t; }
            });
        }
        float target;
        if (targetSum is double ts) target = (float)ts;
        else
        {
            var pos = totals.Where(t => t > 0).ToArray();
            if (pos.Length == 0) target = 1;
            else
            {
                Array.Sort(pos);
                int h = pos.Length / 2;
                target = pos.Length % 2 == 1 ? pos[h] : (float)((pos[h - 1] + (double)pos[h]) / 2);
            }
        }
        var f = new float[n];
        for (int i = 0; i < n; i++) f[i] = totals[i] > 0 ? (float)(totals[i] / (double)target) : 1f;
        return f;
    }

    private static void Normalize(int start, SparseMatrix b, float[] factors, bool log1p)
    {
        for (int i = 0; i < b.Rows; i++)
        {
            double f = factors[start + i];
            var v = b.RowValues(i);
            for (int k = 0; k < v.Length; k++)
            {
                float y = (float)(v[k] / f);
                v[k] = log1p ? (float)Math.Log(1 + (double)y) : y;
            }
        }
    }

    /// <summary>
    /// Highly variable genes. <see cref="HvgFlavor.Seurat"/> expects log-normalised data and takes one pass;
    /// <see cref="HvgFlavor.SeuratV3"/> expects raw counts and takes two (moments, then clipped moments). Both follow
    /// scanpy 1.12 exactly (pandas-style equal-width mean bins; netlib loess with <paramref name="span"/> = 0.3,
    /// degree 2 and interpolated surface for v3).
    /// </summary>
    public static HvgResult HighlyVariableGenes(ISparseMatrixSource x, int nTopGenes = 2000, HvgFlavor flavor = HvgFlavor.SeuratV3, int bins = 20,
        double span = 0.3, int blockRows = SparseMatrixSource.DefaultBlockRows, int? maxParallelism = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentOutOfRangeException.ThrowIfLessThan(nTopGenes, 1);
        if (flavor == HvgFlavor.Seurat)
        {
            var (mean, var) = GeneMeanVariance(x, expm1: true, blockRows, maxParallelism);
            return SeuratFromMoments(mean, var, nTopGenes, bins);
        }
        var (m, v) = GeneMeanVariance(x, expm1: false, blockRows, maxParallelism);
        var (regStd, clip) = SeuratV3Fit(m, v, x.Rows, span);
        var (cs, cq) = ColumnMoments(x, false, clip, blockRows, maxParallelism);
        return SeuratV3FromMoments(m, v, regStd, cs, cq, x.Rows, nTopGenes);
    }

    internal static HvgResult SeuratFromMoments(double[] meanIn, double[] var, int nTop, int nBins)
    {
        int g = meanIn.Length;
        var mean = new double[g]; var disp = new double[g];
        for (int j = 0; j < g; j++)
        {
            double mu = meanIn[j] == 0 ? 1e-12 : meanIn[j];
            double d = var[j] / mu;
            disp[j] = d == 0 ? double.NaN : Math.Log(d);
            mean[j] = Math.Log(1 + mu);
        }
        var bin = PandasCut(mean, nBins);
        // groupby(mean_bin).agg(mean, std) skipping NaN.
        var cnt = new int[nBins]; var sum = new double[nBins];
        for (int j = 0; j < g; j++) if (!double.IsNaN(disp[j])) { cnt[bin[j]]++; sum[bin[j]] += disp[j]; }
        var avg = new double[nBins]; var dev = new double[nBins];
        for (int b = 0; b < nBins; b++) avg[b] = cnt[b] > 0 ? sum[b] / cnt[b] : double.NaN;
        var ss = new double[nBins];
        for (int j = 0; j < g; j++) if (!double.IsNaN(disp[j])) { double t = disp[j] - avg[bin[j]]; ss[bin[j]] += t * t; }
        for (int b = 0; b < nBins; b++)
        {
            dev[b] = cnt[b] > 1 ? Math.Sqrt(ss[b] / (cnt[b] - 1)) : double.NaN;
            if (double.IsNaN(dev[b])) { dev[b] = avg[b]; avg[b] = 0; } // single-gene bins: z = disp / disp = 1
        }
        var norm = new double[g];
        for (int j = 0; j < g; j++) norm[j] = (disp[j] - avg[bin[j]]) / dev[bin[j]];
        var finite = norm.Where(z => !double.IsNaN(z)).OrderDescending().ToArray();
        int take = Math.Min(nTop, finite.Length);
        var hv = new bool[g];
        if (take > 0)
        {
            double cut = finite[take - 1];
            for (int j = 0; j < g; j++) hv[j] = !double.IsNaN(norm[j]) && norm[j] >= cut;
        }
        return new HvgResult(HvgFlavor.Seurat, mean, disp, norm, hv, null);
    }

    /// <summary>pandas.cut(values, bins=n): equal-width right-closed bins over [min, max], the first edge lowered by 0.1 % of the range.</summary>
    internal static int[] PandasCut(double[] values, int nBins)
    {
        double mn = values.Min(), mx = values.Max();
        var edges = new double[nBins + 1];
        if (mn == mx)
        {
            mn -= mn != 0 ? 0.001 * Math.Abs(mn) : 0.001;
            mx += mx != 0 ? 0.001 * Math.Abs(mx) : 0.001;
            double step0 = (mx - mn) / nBins;
            for (int i = 0; i < nBins; i++) edges[i] = i * step0 + mn;
            edges[nBins] = mx;
        }
        else
        {
            double step = (mx - mn) / nBins;
            for (int i = 0; i < nBins; i++) edges[i] = i * step + mn;
            edges[nBins] = mx;
            edges[0] -= (mx - mn) * 0.001;
        }
        var bin = new int[values.Length];
        for (int j = 0; j < values.Length; j++)
        {
            // searchsorted(edges, x, side="left") - 1  =>  edges[b] < x <= edges[b + 1]
            int lo = 0, hi = edges.Length;
            while (lo < hi) { int mid = (lo + hi) >> 1; if (edges[mid] < values[j]) lo = mid + 1; else hi = mid; }
            bin[j] = Math.Clamp(lo - 1, 0, nBins - 1);
        }
        return bin;
    }

    /// <summary>Loess fit of log10 variance on log10 mean; returns the regularised standard deviation and the clip values.</summary>
    internal static (double[] RegStd, double[] Clip) SeuratV3Fit(double[] mean, double[] var, int n, double span)
    {
        int g = mean.Length;
        var idx = Enumerable.Range(0, g).Where(j => var[j] > 0).ToArray();
        var estimate = new double[g];
        if (idx.Length > 0)
        {
            var lx = idx.Select(j => Math.Log10(mean[j])).ToArray();
            var ly = idx.Select(j => Math.Log10(var[j])).ToArray();
            var fit = RLoess.PredictInterpolate(lx, ly, span, lx);
            for (int k = 0; k < idx.Length; k++) estimate[idx[k]] = fit[k];
        }
        var regStd = new double[g]; var clip = new double[g];
        double sq = Math.Sqrt(n);
        for (int j = 0; j < g; j++)
        {
            regStd[j] = Math.Sqrt(Math.Pow(10, estimate[j]));
            clip[j] = regStd[j] * sq + mean[j];
        }
        return (regStd, clip);
    }

    internal static HvgResult SeuratV3FromMoments(double[] mean, double[] var, double[] regStd, double[] clipSum, double[] clipSq, int n, int nTop)
    {
        int g = mean.Length;
        var norm = new double[g];
        for (int j = 0; j < g; j++)
            norm[j] = 1 / ((n - 1) * regStd[j] * regStd[j]) * (n * mean[j] * mean[j] + clipSq[j] - 2 * clipSum[j] * mean[j]);
        // argsort(argsort(-norm)) with ties broken by index; ranks >= nTop become NaN.
        var order = Enumerable.Range(0, g).OrderByDescending(j => norm[j]).ThenBy(j => j).ToArray();
        var rank = new double[g];
        Array.Fill(rank, double.NaN);
        var hv = new bool[g];
        int take = Math.Min(nTop, g);
        for (int r = 0; r < take; r++) { rank[order[r]] = r; hv[order[r]] = true; }
        return new HvgResult(HvgFlavor.SeuratV3, mean, var, norm, hv, rank);
    }
}
