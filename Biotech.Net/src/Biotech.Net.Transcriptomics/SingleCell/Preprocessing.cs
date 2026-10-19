using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Highly-variable-gene selection flavours (scanpy's <c>flavor</c> argument).</summary>
public enum HvgFlavor
{
    /// <summary>
    /// Seurat (v1–v2) / scanpy <c>"seurat"</c>: dispersion = variance/mean of expm1(log data), normalised as a z-score
    /// within 20 bins of log mean. Expects log-normalised <see cref="SingleCellExperiment.X"/>.
    /// </summary>
    Seurat,

    /// <summary>
    /// Seurat v3 / scanpy <c>"seurat_v3"</c>: variance of raw counts standardised by a loess fit of log variance on
    /// log mean, with values clipped at μ + √N·σ. Uses the <c>counts</c> layer (or X if absent).
    /// </summary>
    SeuratV3,
}

/// <summary>
/// Quality control, filtering, normalisation, highly variable genes and scaling for single-cell data — the
/// scanpy <c>pp</c> pipeline (<c>calculate_qc_metrics</c>, <c>filter_cells</c>, <c>filter_genes</c>,
/// <c>normalize_total</c>, <c>log1p</c>, <c>highly_variable_genes</c>, <c>scale</c>).
/// </summary>
public static class Preprocessing
{
    /// <summary>Cell metric keys written by <see cref="CalculateQcMetrics"/>.</summary>
    public const string NGenes = "n_genes_by_counts", TotalCounts = "total_counts", PctMito = "pct_counts_mt", TotalMito = "total_counts_mt";

    /// <summary>
    /// Per-cell n_genes_by_counts, total_counts, total_counts_mt and pct_counts_mt (percent, 0–100) and per-gene
    /// n_cells_by_counts, mean_counts, total_counts and the <c>mt</c> flag (symbols starting with
    /// <paramref name="mitoPrefix"/>, case-insensitive).
    /// </summary>
    public static void CalculateQcMetrics(SingleCellExperiment sc, string mitoPrefix = "MT-")
    {
        ArgumentNullException.ThrowIfNull(sc);
        var x = sc.X;
        var mt = new bool[sc.GeneCount];
        for (int j = 0; j < mt.Length; j++) mt[j] = sc.GeneNames[j].StartsWith(mitoPrefix, StringComparison.OrdinalIgnoreCase);
        var nGenes = new double[x.Rows]; var total = new double[x.Rows]; var mito = new double[x.Rows]; var pct = new double[x.Rows];
        for (int i = 0; i < x.Rows; i++)
        {
            var idx = x.RowIndices(i); var v = x.RowValues(i);
            double t = 0, m = 0; int g = 0;
            for (int k = 0; k < idx.Length; k++)
            {
                if (v[k] == 0) continue;
                g++; t += v[k];
                if (mt[idx[k]]) m += v[k];
            }
            nGenes[i] = g; total[i] = t; mito[i] = m; pct[i] = t > 0 ? 100 * m / t : 0;
        }
        sc.CellMetrics[NGenes] = nGenes; sc.CellMetrics[TotalCounts] = total;
        sc.CellMetrics[TotalMito] = mito; sc.CellMetrics[PctMito] = pct;
        sc.GeneFlags["mt"] = mt;
        var nCells = x.ColumnNonZeroCounts();
        var gsum = x.ColumnSums();
        sc.GeneMetrics["n_cells_by_counts"] = nCells.Select(c => (double)c).ToArray();
        sc.GeneMetrics["total_counts"] = gsum;
        sc.GeneMetrics["mean_counts"] = gsum.Select(s => s / Math.Max(1, x.Rows)).ToArray();
    }

    /// <summary>
    /// Removes cells outside the given limits (any limit may be null). Recomputes QC metrics first if missing.
    /// Returns the number of cells removed.
    /// </summary>
    public static int FilterCells(SingleCellExperiment sc, int? minGenes = null, int? maxGenes = null, double? minCounts = null,
        double? maxCounts = null, double? maxPctMito = null, string mitoPrefix = "MT-")
    {
        ArgumentNullException.ThrowIfNull(sc);
        if (!sc.CellMetrics.ContainsKey(NGenes)) CalculateQcMetrics(sc, mitoPrefix);
        var ng = sc.CellMetrics[NGenes]; var tc = sc.CellMetrics[TotalCounts]; var pm = sc.CellMetrics[PctMito];
        var keep = new List<int>(sc.CellCount);
        for (int i = 0; i < sc.CellCount; i++)
        {
            if (minGenes is int a && ng[i] < a) continue;
            if (maxGenes is int b && ng[i] > b) continue;
            if (minCounts is double c && tc[i] < c) continue;
            if (maxCounts is double d && tc[i] > d) continue;
            if (maxPctMito is double e && pm[i] > e) continue;
            keep.Add(i);
        }
        int removed = sc.CellCount - keep.Count;
        if (removed > 0) sc.SubsetCells(keep.ToArray());
        return removed;
    }

    /// <summary>Removes genes detected in fewer than <paramref name="minCells"/> cells. Returns the number removed.</summary>
    public static int FilterGenes(SingleCellExperiment sc, int minCells = 3, double? minCounts = null)
    {
        ArgumentNullException.ThrowIfNull(sc);
        var nz = sc.X.ColumnNonZeroCounts();
        var sums = sc.X.ColumnSums();
        var keep = new List<int>(sc.GeneCount);
        for (int j = 0; j < sc.GeneCount; j++)
            if (nz[j] >= minCells && (minCounts is not double mc || sums[j] >= mc)) keep.Add(j);
        int removed = sc.GeneCount - keep.Count;
        if (removed > 0) sc.SubsetGenes(keep.ToArray());
        return removed;
    }

    /// <summary>
    /// Scales every cell to the same total (<c>normalize_total</c>). When <paramref name="targetSum"/> is null the
    /// median total of cells with non-zero counts is used, as in scanpy. The raw counts are kept in the
    /// <c>counts</c> layer (if not already there). Works in place on <see cref="SingleCellExperiment.X"/>.
    /// </summary>
    public static void NormalizeTotal(SingleCellExperiment sc, double? targetSum = 1e4)
    {
        ArgumentNullException.ThrowIfNull(sc);
        if (!sc.Layers.ContainsKey("counts")) sc.Layers["counts"] = sc.X.Clone();
        var x = sc.X.Clone();
        var totals = x.RowSums();
        double target = targetSum ?? NumericUtil.Median(totals.Where(t => t > 0).ToArray());
        for (int i = 0; i < x.Rows; i++)
        {
            if (totals[i] <= 0) continue;
            float f = (float)(target / totals[i]);
            var v = x.RowValues(i);
            for (int k = 0; k < v.Length; k++) v[k] *= f;
        }
        sc.X = x;
    }

    /// <summary>Natural log(1 + x) of <see cref="SingleCellExperiment.X"/> (sparse-preserving), like scanpy's <c>log1p</c>.</summary>
    public static void Log1p(SingleCellExperiment sc)
    {
        ArgumentNullException.ThrowIfNull(sc);
        var x = sc.X.Clone();
        var v = x.Values;
        for (int k = 0; k < v.Length; k++) v[k] = MathF.Log(1 + v[k]);
        sc.X = x;
    }

    /// <summary>
    /// Flags the <paramref name="nTopGenes"/> most variable genes in <c>GeneFlags["highly_variable"]</c> and stores
    /// the per-gene statistics (<c>means</c>, <c>dispersions</c>/<c>variances</c>, <c>dispersions_norm</c>/<c>variances_norm</c>).
    /// Returns the selected gene indices in ascending order.
    /// </summary>
    public static int[] HighlyVariableGenes(SingleCellExperiment sc, int nTopGenes = 2000, HvgFlavor flavor = HvgFlavor.Seurat, int bins = 20)
    {
        ArgumentNullException.ThrowIfNull(sc);
        // Same implementation as the out-of-core path (scanpy 1.12 semantics): Seurat on X, Seurat v3 on the counts layer.
        var source = flavor == HvgFlavor.Seurat ? sc.X : sc.Layers.TryGetValue("counts", out var c) ? c : sc.X;
        var result = OutOfCore.HighlyVariableGenes(source, nTopGenes, flavor, bins);
        result.WriteTo(sc);
        return result.Selected;
    }

    /// <summary>
    /// Dense cells × genes matrix of the selected genes (all when null), each gene centred to mean 0 and scaled to
    /// unit variance (ddof = 1), values clipped to ±<paramref name="maxValue"/> (scanpy <c>scale(max_value=10)</c>).
    /// Genes with zero variance are left centred.
    /// </summary>
    public static DenseMatrix Scale(SingleCellExperiment sc, int[]? genes = null, double? maxValue = 10, bool zeroCenter = true)
    {
        ArgumentNullException.ThrowIfNull(sc);
        genes ??= Enumerable.Range(0, sc.GeneCount).ToArray();
        var d = sc.X.ToDense(genes);
        int n = d.Rows, p = d.Columns;
        var mean = new double[p]; var sd = new double[p];
        for (int i = 0; i < n; i++) { var r = d.Row(i); for (int j = 0; j < p; j++) mean[j] += r[j]; }
        for (int j = 0; j < p; j++) mean[j] /= n;
        for (int i = 0; i < n; i++) { var r = d.Row(i); for (int j = 0; j < p; j++) { double t = r[j] - mean[j]; sd[j] += t * t; } }
        for (int j = 0; j < p; j++) { sd[j] = n > 1 ? Math.Sqrt(sd[j] / (n - 1)) : 0; if (sd[j] == 0) sd[j] = 1; }
        for (int i = 0; i < n; i++)
        {
            var r = d.Row(i);
            for (int j = 0; j < p; j++)
            {
                double v = (zeroCenter ? r[j] - mean[j] : r[j]) / sd[j];
                if (maxValue is double mv) v = Math.Clamp(v, -mv, mv);
                r[j] = v;
            }
        }
        return d;
    }
}
