using Gravicode.Science.GraviNum;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>How UMAP initialises the layout.</summary>
public enum UmapInit
{
    /// <summary>Eigenvectors of the normalised graph Laplacian (umap-learn's default); falls back to random above 3 000 cells.</summary>
    Spectral,

    /// <summary>Uniform random coordinates in [−10, 10].</summary>
    Random,
}

/// <summary>Options for <see cref="Umap"/>.</summary>
public sealed class UmapOptions
{
    /// <summary>Output dimensions (default 2).</summary>
    public int Components { get; init; } = 2;

    /// <summary>Minimum distance between embedded points (default 0.5, scanpy's default; umap-learn uses 0.1).</summary>
    public double MinDist { get; init; } = 0.5;

    /// <summary>Effective scale of embedded points (default 1).</summary>
    public double Spread { get; init; } = 1.0;

    /// <summary>Optimisation epochs; null = 500 for ≤ 10 000 cells, else 200.</summary>
    public int? Epochs { get; init; }

    /// <summary>Initial SGD learning rate (default 1).</summary>
    public double LearningRate { get; init; } = 1.0;

    /// <summary>Negative samples per positive edge sample (default 5).</summary>
    public int NegativeSampleRate { get; init; } = 5;

    /// <summary>Repulsion strength γ (default 1).</summary>
    public double RepulsionStrength { get; init; } = 1.0;

    /// <summary>Initialisation (default spectral).</summary>
    public UmapInit Init { get; init; } = UmapInit.Spectral;

    /// <summary>RNG seed.</summary>
    public int Seed { get; init; }
}

/// <summary>
/// Uniform Manifold Approximation and Projection (McInnes, Healy &amp; Melville 2018), following umap-learn: the
/// fuzzy simplicial set from <see cref="NeighborGraph"/>, curve parameters (a, b) fitted from min_dist/spread,
/// spectral initialisation and the edge-sampled SGD with negative sampling of <c>optimize_layout_euclidean</c>.
/// </summary>
/// <remarks>
/// The optimisation runs single-threaded so a given seed always gives the same layout (umap-learn's parallel mode
/// is not reproducible either). Spectral initialisation uses a dense eigendecomposition and is used up to 3 000
/// cells; larger inputs are initialised randomly.
/// </remarks>
public static class Umap
{
    /// <summary>
    /// Fits a and b of the low-dimensional similarity 1 / (1 + a·d^(2b)) to the target curve that is 1 below
    /// <paramref name="minDist"/> and exp(−(d − minDist)/spread) above (umap-learn's <c>find_ab_params</c>,
    /// nonlinear least squares on 300 points in [0, 3·spread]).
    /// </summary>
    public static (double A, double B) FindAbParams(double spread = 1.0, double minDist = 0.1)
    {
        const int m = 300;
        var x = new double[m]; var y = new double[m];
        for (int i = 0; i < m; i++)
        {
            x[i] = 3 * spread * i / (m - 1);
            y[i] = x[i] < minDist ? 1 : Math.Exp(-(x[i] - minDist) / spread);
        }
        double a = 1, b = 1, lambda = 1e-3;
        double Sse(double aa, double bb)
        {
            double s = 0;
            for (int i = 0; i < m; i++) { double r = 1 / (1 + aa * Math.Pow(x[i], 2 * bb)) - y[i]; s += r * r; }
            return s;
        }
        double cur = Sse(a, b);
        for (int it = 0; it < 500; it++)
        {
            double jaa = 0, jab = 0, jbb = 0, ga = 0, gb = 0;
            for (int i = 0; i < m; i++)
            {
                if (x[i] == 0) continue;
                double p = Math.Pow(x[i], 2 * b);
                double f = 1 / (1 + a * p);
                double r = f - y[i];
                double da = -p * f * f;
                double db = -a * p * 2 * Math.Log(x[i]) * f * f;
                jaa += da * da; jab += da * db; jbb += db * db; ga += da * r; gb += db * r;
            }
            // Levenberg–Marquardt step.
            double maa = jaa * (1 + lambda), mbb = jbb * (1 + lambda);
            double det = maa * mbb - jab * jab;
            if (Math.Abs(det) < 1e-300) break;
            double sa = -(mbb * ga - jab * gb) / det, sb = -(maa * gb - jab * ga) / det;
            double next = Sse(a + sa, b + sb);
            if (next < cur)
            {
                a += sa; b += sb;
                bool done = Math.Abs(cur - next) < 1e-15 * Math.Max(1, cur);
                cur = next; lambda = Math.Max(lambda / 10, 1e-12);
                if (done) break;
            }
            else lambda *= 10;
            if (lambda > 1e12) break;
        }
        return (a, b);
    }

    /// <summary>Embeds the cells of a neighbour graph.</summary>
    public static DenseMatrix Fit(NeighborGraph graph, UmapOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        options ??= new UmapOptions();
        var w = graph.Connectivities;
        int n = w.Rows, dim = options.Components;
        int epochs = options.Epochs ?? (n <= 10000 ? 500 : 200);
        var (a, b) = FindAbParams(options.Spread, options.MinDist);
        var rng = new Random(options.Seed);

        // Edges (both directions, as umap-learn uses the full symmetric matrix), pruned below max/epochs.
        float maxW = 0;
        foreach (float v in w.Values) maxW = Math.Max(maxW, v);
        var head = new List<int>(); var tail = new List<int>(); var ew = new List<double>();
        for (int i = 0; i < n; i++)
        {
            var idx = w.RowIndices(i); var v = w.RowValues(i);
            for (int t = 0; t < idx.Length; t++)
                if (v[t] >= maxW / epochs && idx[t] != i) { head.Add(i); tail.Add(idx[t]); ew.Add(v[t]); }
        }
        int ne = head.Count;
        var eps = new double[ne];
        for (int e = 0; e < ne; e++) eps[e] = maxW / ew[e];

        var y = Initialise(w, n, dim, options, rng);

        var epsNeg = new double[ne];
        var nextSample = new double[ne];
        var nextNeg = new double[ne];
        for (int e = 0; e < ne; e++)
        {
            epsNeg[e] = eps[e] / options.NegativeSampleRate;
            nextSample[e] = eps[e];
            nextNeg[e] = epsNeg[e];
        }
        double gamma = options.RepulsionStrength;
        for (int epoch = 0; epoch < epochs; epoch++)
        {
            double alpha = options.LearningRate * (1 - (double)epoch / epochs);
            for (int e = 0; e < ne; e++)
            {
                if (nextSample[e] > epoch) continue;
                int j = head[e], k = tail[e];
                double d2 = Dist2(y, j, k, dim);
                double coeff = d2 > 0 ? -2 * a * b * Math.Pow(d2, b - 1) / (a * Math.Pow(d2, b) + 1) : 0;
                for (int d = 0; d < dim; d++)
                {
                    double g = Clip(coeff * (y[j * dim + d] - y[k * dim + d]));
                    y[j * dim + d] += g * alpha;
                    y[k * dim + d] -= g * alpha;
                }
                nextSample[e] += eps[e];
                int nNeg = (int)((epoch - nextNeg[e]) / epsNeg[e]);
                for (int s = 0; s < nNeg; s++)
                {
                    int r = rng.Next(n);
                    if (r == j) continue;
                    double dn = Dist2(y, j, r, dim);
                    double rc = dn > 0 ? 2 * gamma * b / ((0.001 + dn) * (a * Math.Pow(dn, b) + 1)) : 0;
                    for (int d = 0; d < dim; d++)
                    {
                        double g = rc > 0 ? Clip(rc * (y[j * dim + d] - y[r * dim + d])) : 4;
                        y[j * dim + d] += g * alpha;
                    }
                }
                nextNeg[e] += nNeg * epsNeg[e];
            }
        }
        return new DenseMatrix(n, dim, y);
    }

    /// <summary>Runs UMAP on <see cref="SingleCellExperiment.Neighbors"/> and stores <c>X_umap</c>.</summary>
    public static DenseMatrix Fit(SingleCellExperiment sc, UmapOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sc);
        var e = Fit(sc.Neighbors ?? throw new InvalidOperationException("Compute the neighbour graph first."), options);
        sc.Embeddings["X_umap"] = e;
        return e;
    }

    private static double[] Initialise(SparseMatrix w, int n, int dim, UmapOptions o, Random rng)
    {
        var y = new double[n * dim];
        if (o.Init == UmapInit.Spectral && n <= 3000 && n > dim + 1)
        {
            // Normalised Laplacian L = I − D^−½ W D^−½; take eigenvectors 2..dim+1 of the smallest eigenvalues,
            // i.e. of the largest eigenvalues of D^−½ W D^−½ after the trivial one.
            var deg = w.RowSums();
            var m = new double[n * n];
            for (int i = 0; i < n; i++)
            {
                var idx = w.RowIndices(i); var v = w.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    double den = Math.Sqrt(deg[i] * deg[idx[t]]);
                    if (den > 0) m[i * n + idx[t]] = v[t] / den;
                }
            }
            var eig = Decomposition.SymmetricEigen(new NdArray(m, [n, n]));
            var vecs = eig.Vectors.ToArray();
            double maxAbs = 0;
            for (int i = 0; i < n; i++)
                for (int d = 0; d < dim; d++)
                {
                    y[i * dim + d] = vecs[i * n + d + 1];
                    maxAbs = Math.Max(maxAbs, Math.Abs(y[i * dim + d]));
                }
            if (maxAbs > 0)
            {
                double expansion = 10 / maxAbs;
                for (int i = 0; i < y.Length; i++) y[i] = y[i] * expansion + 1e-4 * Pca.Gaussian(rng);
                // umap-learn rescales spectral init to [0, 10] per dimension.
                for (int d = 0; d < dim; d++)
                {
                    double lo = double.MaxValue, hi = double.MinValue;
                    for (int i = 0; i < n; i++) { lo = Math.Min(lo, y[i * dim + d]); hi = Math.Max(hi, y[i * dim + d]); }
                    double span = hi - lo > 0 ? hi - lo : 1;
                    for (int i = 0; i < n; i++) y[i * dim + d] = 10 * (y[i * dim + d] - lo) / span;
                }
                return y;
            }
        }
        for (int i = 0; i < y.Length; i++) y[i] = rng.NextDouble() * 20 - 10;
        return y;
    }

    private static double Dist2(double[] y, int i, int j, int dim)
    {
        double s = 0;
        for (int d = 0; d < dim; d++) { double t = y[i * dim + d] - y[j * dim + d]; s += t * t; }
        return s;
    }

    private static double Clip(double v) => v > 4 ? 4 : v < -4 ? -4 : v;
}

/// <summary>Options for <see cref="Tsne"/>.</summary>
public sealed class TsneOptions
{
    /// <summary>Output dimensions (default 2).</summary>
    public int Components { get; init; } = 2;

    /// <summary>Perplexity (default 30; reduced automatically when n is small).</summary>
    public double Perplexity { get; init; } = 30;

    /// <summary>Gradient-descent iterations (default 1000).</summary>
    public int Iterations { get; init; } = 1000;

    /// <summary>Learning rate; null = scikit-learn's "auto": max(n / exaggeration / 4, 50).</summary>
    public double? LearningRate { get; init; }

    /// <summary>Early exaggeration factor (default 12) applied for the first 250 iterations.</summary>
    public double EarlyExaggeration { get; init; } = 12;

    /// <summary>RNG seed for the initial layout.</summary>
    public int Seed { get; init; }
}

/// <summary>
/// Exact t-distributed stochastic neighbour embedding (van der Maaten &amp; Hinton 2008) with perplexity calibration
/// by binary search, early exaggeration, momentum and adaptive gains, following scikit-learn's <c>method="exact"</c>.
/// Cost is O(n²) per iteration, so it is intended for up to a few thousand cells (Barnes–Hut is not implemented).
/// </summary>
public static class Tsne
{
    /// <summary>Embeds the rows of <paramref name="data"/> (e.g. PCA scores). The final KL divergence is returned too.</summary>
    public static (DenseMatrix Embedding, double KlDivergence) Fit(DenseMatrix data, TsneOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(data);
        options ??= new TsneOptions();
        int n = data.Rows, d = data.Columns, dim = options.Components;
        if (n < 4) throw new ArgumentException("t-SNE needs at least four points.");
        double perplexity = Math.Min(options.Perplexity, (n - 1) / 3.0);
        var p = JointProbabilities(data.Data, n, d, perplexity);

        var rng = new Random(options.Seed);
        var y = new double[n * dim];
        for (int i = 0; i < y.Length; i++) y[i] = 1e-4 * Pca.Gaussian(rng);
        var update = new double[n * dim];
        var gains = new double[n * dim];
        Array.Fill(gains, 1.0);
        double lr = options.LearningRate ?? Math.Max(n / options.EarlyExaggeration / 4, 50);
        var num = new double[n * n];
        var grad = new double[n * dim];
        double kl = 0;
        const int exaggerationIters = 250;
        for (int it = 0; it < options.Iterations; it++)
        {
            double exag = it < exaggerationIters ? options.EarlyExaggeration : 1;
            double momentum = it < exaggerationIters ? 0.5 : 0.8;
            // Student-t affinities.
            double z = 0;
            for (int i = 0; i < n; i++)
            {
                num[i * n + i] = 0;
                for (int j = i + 1; j < n; j++)
                {
                    double s = 0;
                    for (int c = 0; c < dim; c++) { double t = y[i * dim + c] - y[j * dim + c]; s += t * t; }
                    double q = 1 / (1 + s);
                    num[i * n + j] = q; num[j * n + i] = q;
                    z += 2 * q;
                }
            }
            Parallel.For(0, n, i =>
            {
                for (int c = 0; c < dim; c++) grad[i * dim + c] = 0;
                for (int j = 0; j < n; j++)
                {
                    if (i == j) continue;
                    double q = num[i * n + j];
                    double mult = 4 * (exag * p[i * n + j] - q / z) * q;
                    for (int c = 0; c < dim; c++) grad[i * dim + c] += mult * (y[i * dim + c] - y[j * dim + c]);
                }
            });
            for (int k = 0; k < y.Length; k++)
            {
                // scikit-learn: gains grow by 0.2 where update·grad < 0 and shrink by 0.8 otherwise (including zeros).
                bool opposite = update[k] * grad[k] < 0;
                gains[k] = opposite ? gains[k] + 0.2 : Math.Max(gains[k] * 0.8, 0.01);
                update[k] = momentum * update[k] - lr * gains[k] * grad[k];
                y[k] += update[k];
            }
            if (it == options.Iterations - 1)
            {
                kl = 0;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        if (i != j && p[i * n + j] > 0) kl += p[i * n + j] * Math.Log(p[i * n + j] / Math.Max(num[i * n + j] / z, 1e-300));
            }
        }
        return (new DenseMatrix(n, dim, y), kl);
    }

    /// <summary>Runs t-SNE on an embedding of <paramref name="sc"/> (default <c>X_pca</c>) and stores <c>X_tsne</c>.</summary>
    public static DenseMatrix Fit(SingleCellExperiment sc, TsneOptions? options = null, string embedding = "X_pca")
    {
        ArgumentNullException.ThrowIfNull(sc);
        if (!sc.Embeddings.TryGetValue(embedding, out var e)) throw new InvalidOperationException($"Embedding '{embedding}' not found.");
        var (y, _) = Fit(e, options);
        sc.Embeddings["X_tsne"] = y;
        return y;
    }

    /// <summary>Symmetrised joint probabilities P with each conditional row calibrated to the perplexity.</summary>
    internal static double[] JointProbabilities(double[] x, int n, int d, double perplexity)
    {
        var d2 = new double[n * n];
        Parallel.For(0, n, i =>
        {
            for (int j = 0; j < n; j++)
            {
                double s = 0;
                for (int c = 0; c < d; c++) { double t = x[i * d + c] - x[j * d + c]; s += t * t; }
                d2[i * n + j] = s;
            }
        });
        var cond = new double[n * n];
        double logU = Math.Log(perplexity);
        Parallel.For(0, n, i =>
        {
            double beta = 1, lo = double.NegativeInfinity, hi = double.PositiveInfinity;
            var row = cond.AsSpan(i * n, n);
            for (int it = 0; it < 100; it++)
            {
                double sum = 0, hsum = 0;
                for (int j = 0; j < n; j++)
                {
                    row[j] = j == i ? 0 : Math.Exp(-d2[i * n + j] * beta);
                    sum += row[j];
                }
                if (sum == 0) sum = 1e-300;
                for (int j = 0; j < n; j++) { row[j] /= sum; if (row[j] > 0) hsum -= row[j] * Math.Log(row[j]); }
                double diff = hsum - logU;
                if (Math.Abs(diff) < 1e-5) break;
                if (diff > 0) { lo = beta; beta = double.IsPositiveInfinity(hi) ? beta * 2 : (beta + hi) / 2; }
                else { hi = beta; beta = double.IsNegativeInfinity(lo) ? beta / 2 : (beta + lo) / 2; }
            }
        });
        var p = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                p[i * n + j] = Math.Max((cond[i * n + j] + cond[j * n + i]) / (2.0 * n), 1e-12);
        return p;
    }
}
