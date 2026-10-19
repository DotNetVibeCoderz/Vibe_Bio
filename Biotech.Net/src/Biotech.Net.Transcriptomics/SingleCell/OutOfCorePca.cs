using Gravicode.Science.GraviNum;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Solver for <see cref="OutOfCorePca"/>.</summary>
public enum OutOfCorePcaSolver
{
    /// <summary><see cref="Covariance"/> for up to 4 096 genes, <see cref="Randomized"/> beyond.</summary>
    Auto,

    /// <summary>
    /// Exact: accumulates the genes × genes covariance in one pass (sparse outer products, implicit centring),
    /// eigendecomposes it, then projects in a second pass. Memory O(genes²) per worker thread. Same maths as
    /// scanpy's <c>svd_solver="covariance_eigh"</c>.
    /// </summary>
    Covariance,

    /// <summary>
    /// Randomized SVD (Halko, Martinsson &amp; Tropp 2011) with implicit centring: 3 + 2·powerIterations passes,
    /// memory O(cells × (components + oversampling)). Approximate; accuracy improves with power iterations.
    /// </summary>
    Randomized,
}

/// <summary>
/// Principal components of a large sparse cells × genes matrix without densifying it: the matrix is streamed in
/// row blocks from any <see cref="ISparseMatrixSource"/> (h5ad-backed, memory-mapped or in memory) and centred
/// implicitly, as scanpy does for sparse input. Returns the same <see cref="PcaResult"/> as <see cref="Pca.Fit"/>
/// (scores feed straight into <see cref="NeighborGraph.Build(DenseMatrix, int, NeighborMetric)"/>), with the same
/// sign convention (largest-magnitude loading of each component positive).
/// </summary>
public static class OutOfCorePca
{
    /// <summary>
    /// Fits <paramref name="components"/> principal components on the selected <paramref name="genes"/> (all when
    /// null; typically the highly variable genes of log-normalised data). With <paramref name="scale"/> each gene is
    /// also divided by its standard deviation (ddof = 1), implicitly and without scanpy's clipping at ±10.
    /// </summary>
    public static PcaResult Fit(ISparseMatrixSource x, int[]? genes = null, int components = 50, OutOfCorePcaSolver solver = OutOfCorePcaSolver.Auto,
        bool scale = false, int seed = 0, int powerIterations = 4, int oversamples = 10, int blockRows = SparseMatrixSource.DefaultBlockRows,
        int? maxParallelism = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        int n = x.Rows;
        if (n < 2) throw new ArgumentException("PCA needs at least two cells.");
        var sel = (genes ?? Enumerable.Range(0, x.Columns).ToArray()).Distinct().Order().ToArray();
        if (sel.Length == 0 || sel[0] < 0 || sel[^1] >= x.Columns) throw new ArgumentOutOfRangeException(nameof(genes));
        int p = sel.Length;
        int k = Math.Min(components, Math.Min(n, p));
        var map = new int[x.Columns];
        Array.Fill(map, -1);
        for (int j = 0; j < p; j++) map[sel[j]] = j;
        if (solver == OutOfCorePcaSolver.Auto) solver = p <= 4096 ? OutOfCorePcaSolver.Covariance : OutOfCorePcaSolver.Randomized;
        int workers = Math.Max(1, maxParallelism ?? Environment.ProcessorCount);
        return solver == OutOfCorePcaSolver.Covariance
            ? FitCovariance(x, map, p, k, scale, blockRows, workers)
            : FitRandomized(x, map, p, k, scale, seed, powerIterations, oversamples, blockRows, workers);
    }

    // ------------------------------------------------------------------ covariance

    private static PcaResult FitCovariance(ISparseMatrixSource x, int[] map, int p, int k, bool scale, int blockRows, int workers)
    {
        int n = x.Rows;
        // One p×p accumulator (upper triangle). Blocks are read sequentially; within a block the T threads split the
        // rows of C (row a goes to thread a mod T), so each thread updates a disjoint 1/T of the matrix and no
        // per-thread copies are needed.
        var cov = new double[(long)p * p];
        var sum = new double[p];
        int threads = Math.Max(1, workers);
        foreach (var (_, b) in x.RowBlocks(blockRows))
        {
            // Restrict the block to the selected genes (CSR over 0..p-1).
            var ptr = new int[b.Rows + 1];
            var cols = new int[b.NonZeroCount]; var vals = new double[b.NonZeroCount];
            int m = 0;
            for (int i = 0; i < b.Rows; i++)
            {
                var idx = b.RowIndices(i); var v = b.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    int j = map[idx[t]];
                    if (j >= 0 && v[t] != 0) { cols[m] = j; vals[m] = v[t]; sum[j] += v[t]; m++; }
                }
                ptr[i + 1] = m;
            }
            Parallel.For(0, threads, new ParallelOptions { MaxDegreeOfParallelism = threads }, t =>
            {
                for (int i = 0; i < b.Rows; i++)
                {
                    int s0 = ptr[i], e0 = ptr[i + 1];
                    for (int a = s0; a < e0; a++)
                    {
                        int ca = cols[a];
                        if (ca % threads != t) continue;
                        double va = vals[a];
                        long row = (long)ca * p;
                        for (int bb = a; bb < e0; bb++) cov[row + cols[bb]] += va * vals[bb];
                    }
                }
            });
        }
        var mean = new double[p];
        for (int j = 0; j < p; j++) mean[j] = sum[j] / n;
        // Symmetrise and centre: C = (XᵀX − n μμᵀ) / (n − 1).
        for (int i = 0; i < p; i++)
            for (int j = i; j < p; j++)
            {
                double v = (cov[(long)i * p + j] - n * mean[i] * mean[j]) / (n - 1);
                cov[(long)i * p + j] = v; cov[(long)j * p + i] = v;
            }
        var sd = new double[p];
        for (int j = 0; j < p; j++) { sd[j] = scale ? Math.Sqrt(Math.Max(cov[(long)j * p + j], 0)) : 1; if (sd[j] == 0) sd[j] = 1; }
        if (scale)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) cov[(long)i * p + j] /= sd[i] * sd[j];
        double total = 0;
        for (int j = 0; j < p; j++) total += cov[(long)j * p + j];
        var (vals2, vecs) = TopEigen(cov, p, k, threads);
        var loadings = vecs;
        FixSigns(loadings, p, k);
        var scores = Project(x, map, p, k, loadings, mean, sd, blockRows, workers);
        var variance = vals2.Select(v => Math.Max(v, 0)).ToArray();
        var ratio = variance.Select(v => total > 0 ? v / total : 0).ToArray();
        return new PcaResult(scores, new DenseMatrix(p, k, loadings), variance, ratio, mean);
    }

    /// <summary>Matrices up to this size use the dense Householder/QL eigensolver; larger ones subspace iteration.</summary>
    internal const int DenseEigenLimit = 600;

    /// <summary>
    /// Top-k eigenpairs of a symmetric p×p matrix (descending), as eigenvalues and p×k row-major vectors. Small
    /// matrices use GraviNum's full eigendecomposition; larger ones Chebyshev-filtered block subspace iteration with
    /// Rayleigh–Ritz (block size max(2k, k + 20)), iterated until every wanted Ritz pair has residual
    /// ‖Av − θv‖ ≤ 1e-11·‖A‖ (exact to working precision, unlike the randomized solver's fixed number of passes).
    /// Assumes A is positive semi-definite (a covariance matrix).
    /// </summary>
    internal static (double[] Values, double[] Vectors) TopEigen(double[] a, int p, int k, int threads = 1)
    {
        if (p <= DenseEigenLimit) return DenseTopEigen(a, p, k);
        int l = Math.Min(p, Math.Max(2 * k, k + 20));
        var rng = new Random(12345);
        var q = new double[p * l];
        for (int t = 0; t < q.Length; t++) q[t] = rng.NextDouble() - 0.5;
        Orthonormalize(q, p, l);
        double norm = 0;
        for (int j = 0; j < p; j++) norm = Math.Max(norm, Math.Abs(a[(long)j * p + j]));
        var z = new double[p * l];
        double[] theta = [], ritz = [];
        for (int iter = 0; iter < 3000; iter++)
        {
            MultiplySymmetric(a, p, q, l, z, threads);           // Z = A Q
            // Rayleigh–Ritz: T = Qᵀ A Q, eig T = W Θ Wᵀ, Ritz vectors Q W, A Q W = Z W.
            var tm = new double[l * l];
            for (int r = 0; r < p; r++)
                for (int i = 0; i < l; i++)
                {
                    double qi = q[r * l + i];
                    for (int j = 0; j < l; j++) tm[i * l + j] += qi * z[r * l + j];
                }
            for (int i = 0; i < l; i++) for (int j = i + 1; j < l; j++) { double v = (tm[i * l + j] + tm[j * l + i]) / 2; tm[i * l + j] = v; tm[j * l + i] = v; }
            var (tv, w) = DenseTopEigen(tm, l, l);
            theta = tv;
            ritz = new double[p * l];
            var az = new double[p * l];
            for (int r = 0; r < p; r++)
                for (int c = 0; c < l; c++)
                {
                    double s1 = 0, s2 = 0;
                    for (int i = 0; i < l; i++) { s1 += q[r * l + i] * w[i * l + c]; s2 += z[r * l + i] * w[i * l + c]; }
                    ritz[r * l + c] = s1; az[r * l + c] = s2;
                }
            norm = Math.Max(norm, Math.Abs(theta[0]));
            double worst = 0;
            for (int c = 0; c < k; c++)
            {
                double res = 0;
                for (int r = 0; r < p; r++) { double d = az[r * l + c] - theta[c] * ritz[r * l + c]; res += d * d; }
                worst = Math.Max(worst, Math.Sqrt(res));
            }
            if (worst <= 1e-11 * norm) break;
            // Next iterate: a degree-8 Chebyshev filter that damps [0, θ_l] (the unwanted part of a PSD spectrum)
            // applied to the Ritz vectors (Zhou & Saad 2007); plain A·(Ritz vectors) when θ_l is ~0.
            double upper = theta[l - 1];
            if (upper > 1e-9 * norm) Chebyshev(a, p, ritz, l, 8, -1e-12 * norm, upper, threads, q);
            else Array.Copy(az, q, az.Length);
            Orthonormalize(q, p, l);
        }
        var vecs = new double[p * k];
        for (int r = 0; r < p; r++) for (int c = 0; c < k; c++) vecs[r * k + c] = ritz[r * l + c];
        return (theta[..k], vecs);
    }

    private static (double[] Values, double[] Vectors) DenseTopEigen(double[] a, int p, int k)
    {
        var e = Decomposition.SymmetricEigen(new NdArray(a, [p, p]));
        var vals = e.Values.ToArray(); var vecs = e.Vectors.ToArray();
        var order = Enumerable.Range(0, p).OrderByDescending(i => vals[i]).Take(k).ToArray();
        var outV = new double[p * k];
        for (int c = 0; c < k; c++) for (int j = 0; j < p; j++) outV[j * k + c] = vecs[j * p + order[c]];
        return (order.Select(i => vals[i]).ToArray(), outV);
    }

    /// <summary>
    /// Y = T_m(L(A)) X with L mapping [lower, upper] onto [−1, 1] (three-term recurrence), which leaves eigenvalues in
    /// [lower, upper] bounded by 1 and amplifies those above <paramref name="upper"/> like cosh(m·acosh(L(λ))).
    /// </summary>
    private static void Chebyshev(double[] a, int p, double[] x, int l, int degree, double lower, double upper, int threads, double[] y)
    {
        double e = (upper - lower) / 2, c = (upper + lower) / 2;
        var y0 = (double[])x.Clone();
        var y1 = new double[p * l];
        var ay = new double[p * l];
        MultiplySymmetric(a, p, y0, l, ay, threads);
        for (int t = 0; t < y1.Length; t++) y1[t] = (ay[t] - c * y0[t]) / e;
        for (int d = 2; d <= degree; d++)
        {
            MultiplySymmetric(a, p, y1, l, ay, threads);
            for (int t = 0; t < y1.Length; t++)
            {
                double next = 2 * (ay[t] - c * y1[t]) / e - y0[t];
                y0[t] = y1[t]; y1[t] = next;
            }
            // Rescale to keep magnitudes bounded (the subspace, not the scale, matters).
            double mx = 0;
            for (int t = 0; t < y1.Length; t++) mx = Math.Max(mx, Math.Abs(y1[t]));
            if (mx > 1e100) for (int t = 0; t < y1.Length; t++) { y1[t] /= mx; y0[t] /= mx; }
        }
        Array.Copy(y1, y, y1.Length);
    }

    /// <summary>Z = A Q for symmetric A (p×p, row-major) and Q (p×l, row-major), parallel over rows of A.</summary>
    private static void MultiplySymmetric(double[] a, int p, double[] q, int l, double[] z, int threads)
    {
        Parallel.For(0, p, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, threads) }, r =>
        {
            var zr = z.AsSpan(r * l, l);
            zr.Clear();
            var ar = a.AsSpan((int)((long)r * p), p);
            for (int j = 0; j < p; j++)
            {
                double v = ar[j];
                if (v == 0) continue;
                var qj = q.AsSpan(j * l, l);
                for (int c = 0; c < l; c++) zr[c] += v * qj[c];
            }
        });
    }

    private static void FixSigns(double[] loadings, int p, int k)
    {
        for (int c = 0; c < k; c++)
        {
            double best = 0; int bi = 0;
            for (int j = 0; j < p; j++) if (Math.Abs(loadings[j * k + c]) > best) { best = Math.Abs(loadings[j * k + c]); bi = j; }
            if (loadings[bi * k + c] < 0) for (int j = 0; j < p; j++) loadings[j * k + c] = -loadings[j * k + c];
        }
    }

    /// <summary>Scores = ((X − μ) / sd) V, one streaming pass.</summary>
    private static DenseMatrix Project(ISparseMatrixSource x, int[] map, int p, int k, double[] v, double[] mean, double[] sd, int blockRows, int workers)
    {
        var w = new double[p * k];
        var offset = new double[k];
        for (int j = 0; j < p; j++)
            for (int c = 0; c < k; c++) { w[j * k + c] = v[j * k + c] / sd[j]; offset[c] += mean[j] / sd[j] * v[j * k + c]; }
        var scores = new DenseMatrix(x.Rows, k);
        SparseMatrixSource.ForEachBlock(x, blockRows, (start, b) =>
        {
            for (int i = 0; i < b.Rows; i++)
            {
                var row = scores.Row(start + i);
                for (int c = 0; c < k; c++) row[c] = -offset[c];
                var idx = b.RowIndices(i); var vals = b.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    int j = map[idx[t]];
                    if (j < 0) continue;
                    double val = vals[t];
                    var wj = w.AsSpan(j * k, k);
                    for (int c = 0; c < k; c++) row[c] += val * wj[c];
                }
            }
        }, workers);
        return scores;
    }

    // ------------------------------------------------------------------ randomized

    private static PcaResult FitRandomized(ISparseMatrixSource x, int[] map, int p, int k, bool scale, int seed, int q, int oversamples, int blockRows, int workers)
    {
        int n = x.Rows;
        int l = Math.Min(k + Math.Max(0, oversamples), Math.Min(n, p));
        // Pass 0: moments.
        var sum = new double[p]; var sq = new double[p];
        var sync = new object();
        SparseMatrixSource.ForEachBlock(x, blockRows, (_, b) =>
        {
            var ls = new double[p]; var lq = new double[p];
            for (int t = 0; t < b.NonZeroCount; t++)
            {
                int j = map[b.ColumnIndices[t]];
                if (j < 0) continue;
                double v = b.Values[t]; ls[j] += v; lq[j] += v * v;
            }
            lock (sync) for (int j = 0; j < p; j++) { sum[j] += ls[j]; sq[j] += lq[j]; }
        }, workers);
        var mean = new double[p]; var sd = new double[p];
        double total = 0;
        for (int j = 0; j < p; j++)
        {
            mean[j] = sum[j] / n;
            double var = Math.Max(0, (sq[j] - n * mean[j] * mean[j]) / (n - 1));
            sd[j] = scale && var > 0 ? Math.Sqrt(var) : 1;
            total += var / (sd[j] * sd[j]);
        }
        var rng = new Random(seed);
        var omega = new double[p * l];
        for (int t = 0; t < omega.Length; t++) omega[t] = Gaussian(rng);
        var y = Times(x, map, p, l, omega, mean, sd, blockRows, workers);   // n×l
        Orthonormalize(y, n, l);
        for (int it = 0; it < q; it++)
        {
            var z = TransposeTimes(x, map, p, l, y, mean, sd, blockRows, workers); // p×l
            Orthonormalize(z, p, l);
            y = Times(x, map, p, l, z, mean, sd, blockRows, workers);
            Orthonormalize(y, n, l);
        }
        var bt = TransposeTimes(x, map, p, l, y, mean, sd, blockRows, workers); // p×l = (Qᵀ A)ᵀ
        // Small SVD of B (l×p) through B Bᵀ = U S² Uᵀ.
        var bbt = new double[l * l];
        for (int j = 0; j < p; j++)
            for (int a = 0; a < l; a++)
            {
                double ba = bt[j * l + a];
                for (int c = 0; c < l; c++) bbt[a * l + c] += ba * bt[j * l + c];
            }
        var (s2, u) = TopEigen(bbt, l, k, workers);
        var loadings = new double[p * k];
        for (int c = 0; c < k; c++)
        {
            double s = Math.Sqrt(Math.Max(s2[c], 1e-300));
            for (int j = 0; j < p; j++)
            {
                double acc = 0;
                for (int a = 0; a < l; a++) acc += bt[j * l + a] * u[a * k + c];
                loadings[j * k + c] = acc / s;
            }
        }
        FixSigns(loadings, p, k);
        var scores = Project(x, map, p, k, loadings, mean, sd, blockRows, workers);
        var variance = s2.Select(v => Math.Max(v, 0) / (n - 1)).ToArray();
        var ratio = variance.Select(v => total > 0 ? v / total : 0).ToArray();
        return new PcaResult(scores, new DenseMatrix(p, k, loadings), variance, ratio, mean);
    }

    private static double Gaussian(Random r)
    {
        double u1 = 1 - r.NextDouble(), u2 = r.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }

    /// <summary>Ã M for Ã = (X − 1μᵀ) D⁻¹, M p×l; returns n×l.</summary>
    private static double[] Times(ISparseMatrixSource x, int[] map, int p, int l, double[] m, double[] mean, double[] sd, int blockRows, int workers)
    {
        var w = new double[p * l];
        var offset = new double[l];
        for (int j = 0; j < p; j++)
            for (int c = 0; c < l; c++) { w[j * l + c] = m[j * l + c] / sd[j]; offset[c] += mean[j] / sd[j] * m[j * l + c]; }
        var y = new double[(long)x.Rows * l];
        SparseMatrixSource.ForEachBlock(x, blockRows, (start, b) =>
        {
            for (int i = 0; i < b.Rows; i++)
            {
                var row = y.AsSpan((start + i) * l, l);
                for (int c = 0; c < l; c++) row[c] = -offset[c];
                var idx = b.RowIndices(i); var vals = b.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    int j = map[idx[t]];
                    if (j < 0) continue;
                    double val = vals[t];
                    var wj = w.AsSpan(j * l, l);
                    for (int c = 0; c < l; c++) row[c] += val * wj[c];
                }
            }
        }, workers);
        return y;
    }

    /// <summary>Ãᵀ Y for Y n×l; returns p×l.</summary>
    private static double[] TransposeTimes(ISparseMatrixSource x, int[] map, int p, int l, double[] y, double[] mean, double[] sd, int blockRows, int workers)
    {
        var z = new double[p * l];
        var colSum = new double[l];
        var sync = new object();
        SparseMatrixSource.ForEachBlock(x, blockRows, (start, b) =>
        {
            var lz = new double[p * l]; var lc = new double[l];
            for (int i = 0; i < b.Rows; i++)
            {
                var yr = y.AsSpan((start + i) * l, l);
                for (int c = 0; c < l; c++) lc[c] += yr[c];
                var idx = b.RowIndices(i); var vals = b.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    int j = map[idx[t]];
                    if (j < 0) continue;
                    double val = vals[t];
                    var zj = lz.AsSpan(j * l, l);
                    for (int c = 0; c < l; c++) zj[c] += val * yr[c];
                }
            }
            lock (sync)
            {
                for (int t = 0; t < z.Length; t++) z[t] += lz[t];
                for (int c = 0; c < l; c++) colSum[c] += lc[c];
            }
        }, workers);
        for (int j = 0; j < p; j++)
            for (int c = 0; c < l; c++) z[j * l + c] = (z[j * l + c] - mean[j] * colSum[c]) / sd[j];
        return z;
    }

    /// <summary>Modified Gram–Schmidt on the columns of a rows×l row-major matrix, twice for stability.</summary>
    private static void Orthonormalize(double[] a, int rows, int l)
    {
        for (int pass = 0; pass < 2; pass++)
        {
            for (int c = 0; c < l; c++)
            {
                for (int prev = 0; prev < c; prev++)
                {
                    double dot = 0;
                    for (long r = 0; r < rows; r++) dot += a[r * l + c] * a[r * l + prev];
                    for (long r = 0; r < rows; r++) a[r * l + c] -= dot * a[r * l + prev];
                }
                double norm = 0;
                for (long r = 0; r < rows; r++) norm += a[r * l + c] * a[r * l + c];
                norm = Math.Sqrt(norm);
                if (norm > 0) for (long r = 0; r < rows; r++) a[r * l + c] /= norm;
            }
        }
    }
}
