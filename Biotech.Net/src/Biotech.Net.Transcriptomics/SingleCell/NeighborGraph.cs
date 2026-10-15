namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Distance used to find nearest neighbours.</summary>
public enum NeighborMetric
{
    /// <summary>Euclidean distance.</summary>
    Euclidean,

    /// <summary>Cosine distance 1 − cos θ.</summary>
    Cosine,
}

/// <summary>
/// The k-nearest-neighbour graph of cells and its UMAP fuzzy-simplicial-set connectivities, as scanpy's
/// <c>pp.neighbors</c> computes them (McInnes, Healy &amp; Melville 2018).
/// </summary>
/// <remarks>
/// Neighbours are found exactly by brute force (O(n²·d), parallel over cells), which is fine up to a few
/// tens of thousands of cells in PCA space. As in scanpy and umap-learn, <c>nNeighbors</c> counts the cell itself,
/// so each cell stores <c>nNeighbors − 1</c> true neighbours.
/// </remarks>
public sealed class NeighborGraph
{
    private NeighborGraph(int n, int k, int[] indices, double[] distances, SparseMatrix connectivities)
    {
        CellCount = n; K = k; Indices = indices; Distances = distances; Connectivities = connectivities;
    }

    /// <summary>Number of cells.</summary>
    public int CellCount { get; }

    /// <summary>Stored neighbours per cell (excluding the cell itself).</summary>
    public int K { get; }

    /// <summary>Neighbour indices, cells × K, nearest first.</summary>
    public int[] Indices { get; }

    /// <summary>Neighbour distances, cells × K.</summary>
    public double[] Distances { get; }

    /// <summary>Symmetric fuzzy-union membership strengths (UMAP graph), used for clustering and UMAP.</summary>
    public SparseMatrix Connectivities { get; }

    /// <summary>Builds the graph from a cells × dimensions embedding (typically the first 30–50 PCs).</summary>
    public static NeighborGraph Build(DenseMatrix data, int nNeighbors = 15, NeighborMetric metric = NeighborMetric.Euclidean)
    {
        ArgumentNullException.ThrowIfNull(data);
        int n = data.Rows, d = data.Columns;
        int k = Math.Min(nNeighbors - 1, n - 1);
        if (k < 1) throw new ArgumentException("Need at least two cells and nNeighbors ≥ 2.");
        var x = data.Data;
        double[]? norms = null;
        if (metric == NeighborMetric.Cosine)
        {
            norms = new double[n];
            for (int i = 0; i < n; i++)
            {
                double s = 0;
                for (int t = 0; t < d; t++) s += x[i * d + t] * x[i * d + t];
                norms[i] = Math.Sqrt(s);
            }
        }
        var idx = new int[n * k];
        var dist = new double[n * k];
        Parallel.For(0, n, () => (new double[k], new int[k]), (i, _, buf) =>
        {
            var (bd, bi) = buf;
            int filled = 0;
            var xi = x.AsSpan(i * d, d);
            for (int j = 0; j < n; j++)
            {
                if (j == i) continue;
                var xj = x.AsSpan(j * d, d);
                double s = 0;
                if (metric == NeighborMetric.Euclidean)
                {
                    for (int t = 0; t < d; t++) { double df = xi[t] - xj[t]; s += df * df; }
                }
                else
                {
                    for (int t = 0; t < d; t++) s += xi[t] * xj[t];
                    double den = norms![i] * norms[j];
                    s = den > 0 ? 1 - s / den : 1;
                }
                // Insertion into the sorted top-k buffer (ties broken by index for determinism).
                if (filled < k) { Insert(bd, bi, ref filled, s, j); }
                else if (s < bd[k - 1]) { filled--; Insert(bd, bi, ref filled, s, j); }
            }
            for (int t = 0; t < k; t++)
            {
                idx[i * k + t] = bi[t];
                dist[i * k + t] = metric == NeighborMetric.Euclidean ? Math.Sqrt(bd[t]) : bd[t];
            }
            return buf;
        }, _ => { });
        return new NeighborGraph(n, k, idx, dist, FuzzySimplicialSet(n, k, idx, dist));
    }

    /// <summary>Builds the graph on an embedding of <paramref name="sc"/> (default <c>X_pca</c>) and stores it.</summary>
    public static NeighborGraph Build(SingleCellExperiment sc, int nNeighbors = 15, string embedding = "X_pca", int? nPcs = null)
    {
        ArgumentNullException.ThrowIfNull(sc);
        if (!sc.Embeddings.TryGetValue(embedding, out var e)) throw new InvalidOperationException($"Embedding '{embedding}' not found; run PCA first.");
        if (nPcs is int np && np < e.Columns)
        {
            var sub = new DenseMatrix(e.Rows, np);
            for (int i = 0; i < e.Rows; i++) e.Row(i)[..np].CopyTo(sub.Row(i));
            e = sub;
        }
        var g = Build(e, nNeighbors);
        sc.Neighbors = g;
        return g;
    }

    private static void Insert(double[] bd, int[] bi, ref int filled, double s, int j)
    {
        int pos = filled;
        while (pos > 0 && bd[pos - 1] > s) { bd[pos] = bd[pos - 1]; bi[pos] = bi[pos - 1]; pos--; }
        bd[pos] = s; bi[pos] = j;
        filled++;
    }

    /// <summary>
    /// umap-learn's <c>fuzzy_simplicial_set</c>: per-cell ρ (distance to the nearest neighbour) and σ found by binary
    /// search so that Σ exp(−(d − ρ)/σ) = log₂(k + 1); memberships are symmetrised by fuzzy union a + b − ab.
    /// </summary>
    internal static SparseMatrix FuzzySimplicialSet(int n, int k, int[] idx, double[] dist)
    {
        const double minKDistScale = 1e-3;
        double target = Math.Log2(k + 1);
        double meanAll = 0;
        foreach (double v in dist) meanAll += v;
        meanAll /= dist.Length;
        var rows = new int[n * k]; var cols = new int[n * k]; var vals = new float[n * k];
        Parallel.For(0, n, i =>
        {
            var di = dist.AsSpan(i * k, k);
            double rho = 0;
            foreach (double v in di) if (v > 0) { rho = v; break; }
            double lo = 0, hi = double.PositiveInfinity, mid = 1;
            for (int it = 0; it < 64; it++)
            {
                double psum = 0;
                foreach (double v in di) { double dd = v - rho; psum += dd > 0 ? Math.Exp(-dd / mid) : 1; }
                if (Math.Abs(psum - target) < 1e-5) break;
                if (psum > target) { hi = mid; mid = (lo + hi) / 2; }
                else { lo = mid; mid = double.IsPositiveInfinity(hi) ? mid * 2 : (lo + hi) / 2; }
            }
            double meanI = 0;
            foreach (double v in di) meanI += v;
            meanI /= k;
            double sigma = rho > 0 ? Math.Max(mid, minKDistScale * meanI) : Math.Max(mid, minKDistScale * meanAll);
            for (int t = 0; t < k; t++)
            {
                double dd = di[t] - rho;
                rows[i * k + t] = i; cols[i * k + t] = idx[i * k + t];
                vals[i * k + t] = (float)(dd <= 0 || sigma == 0 ? 1.0 : Math.Exp(-dd / sigma));
            }
        });
        var a = SparseMatrix.FromTriplets(n, n, rows, cols, vals);
        var at = a.Transpose();
        // Fuzzy union: A + Aᵀ − A∘Aᵀ.
        var r2 = new List<int>(2 * n * k); var c2 = new List<int>(2 * n * k); var v2 = new List<float>(2 * n * k);
        for (int i = 0; i < n; i++)
        {
            var ia = a.RowIndices(i); var va = a.RowValues(i);
            var ib = at.RowIndices(i); var vb = at.RowValues(i);
            int p = 0, q = 0;
            while (p < ia.Length || q < ib.Length)
            {
                int ca = p < ia.Length ? ia[p] : int.MaxValue, cb = q < ib.Length ? ib[q] : int.MaxValue;
                float x = 0, y = 0; int c;
                if (ca == cb) { x = va[p++]; y = vb[q++]; c = ca; }
                else if (ca < cb) { x = va[p++]; c = ca; }
                else { y = vb[q++]; c = cb; }
                r2.Add(i); c2.Add(c); v2.Add(x + y - x * y);
            }
        }
        return SparseMatrix.FromTriplets(n, n, r2.ToArray(), c2.ToArray(), v2.ToArray());
    }
}
