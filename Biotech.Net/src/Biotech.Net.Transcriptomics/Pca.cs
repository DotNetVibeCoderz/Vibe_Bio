using Gravicode.Science.GraviNum;

namespace Biotech.Net.Transcriptomics;

/// <summary>PCA solver.</summary>
public enum PcaSolver
{
    /// <summary>Full decomposition for small problems, randomized otherwise.</summary>
    Auto,

    /// <summary>Exact eigendecomposition of the covariance (or Gram) matrix.</summary>
    Full,

    /// <summary>Randomized range finder with power iterations (Halko, Martinsson &amp; Tropp 2011).</summary>
    Randomized,
}

/// <summary>Principal components of a samples × features matrix.</summary>
public sealed class PcaResult
{
    internal PcaResult(DenseMatrix scores, DenseMatrix loadings, double[] variance, double[] ratio, double[] means)
    {
        Scores = scores; Loadings = loadings; ExplainedVariance = variance; ExplainedVarianceRatio = ratio; Means = means;
    }

    /// <summary>Sample coordinates (samples × components).</summary>
    public DenseMatrix Scores { get; }

    /// <summary>Feature loadings (features × components), unit-norm columns.</summary>
    public DenseMatrix Loadings { get; }

    /// <summary>Variance captured by each component (eigenvalue / (n − 1)).</summary>
    public double[] ExplainedVariance { get; }

    /// <summary>Fraction of total variance captured by each component.</summary>
    public double[] ExplainedVarianceRatio { get; }

    /// <summary>Column means that were subtracted.</summary>
    public double[] Means { get; }

    /// <summary>Number of components.</summary>
    public int Components => ExplainedVariance.Length;
}

/// <summary>
/// Principal component analysis of a dense samples × features matrix (cells × genes, or samples × genes for
/// bulk data). The exact solver and the randomized solver's small eigenproblems use GraviNum's
/// Householder/QL <c>Decomposition.SymmetricEigen</c>. Component signs are fixed so that the largest-magnitude
/// loading of each component is positive, making results deterministic.
/// </summary>
public static class Pca
{
    /// <summary>Computes the first <paramref name="components"/> principal components.</summary>
    public static PcaResult Fit(DenseMatrix data, int components = 50, PcaSolver solver = PcaSolver.Auto, int seed = 0,
        int powerIterations = 4, bool center = true)
    {
        ArgumentNullException.ThrowIfNull(data);
        int n = data.Rows, p = data.Columns;
        if (n < 2) throw new ArgumentException("PCA needs at least two samples.");
        int k = Math.Min(components, Math.Min(n, p));
        var a = (double[])data.Data.Clone();
        var means = new double[p];
        if (center)
        {
            for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) means[j] += a[i * p + j];
            for (int j = 0; j < p; j++) means[j] /= n;
            for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) a[i * p + j] -= means[j];
        }
        double totalVar = 0;
        foreach (double v in a) totalVar += v * v;
        totalVar /= n - 1;

        if (solver == PcaSolver.Auto) solver = Math.Min(n, p) <= 1000 || k > Math.Min(n, p) / 2 ? PcaSolver.Full : PcaSolver.Randomized;
        double[] loadings; double[] eig;
        if (solver == PcaSolver.Full) (loadings, eig) = FullSolve(a, n, p, k);
        else (loadings, eig) = RandomizedSolve(a, n, p, k, seed, powerIterations);

        // Sign convention: largest |loading| positive.
        for (int c = 0; c < k; c++)
        {
            double best = 0; int bi = 0;
            for (int j = 0; j < p; j++) if (Math.Abs(loadings[j * k + c]) > best) { best = Math.Abs(loadings[j * k + c]); bi = j; }
            if (loadings[bi * k + c] < 0) for (int j = 0; j < p; j++) loadings[j * k + c] = -loadings[j * k + c];
        }
        var scores = MatMul(a, n, p, loadings, k);
        var variance = new double[k]; var ratio = new double[k];
        for (int c = 0; c < k; c++)
        {
            variance[c] = Math.Max(eig[c], 0) / (n - 1);
            ratio[c] = totalVar > 0 ? variance[c] / totalVar : 0;
        }
        return new PcaResult(new DenseMatrix(n, k, scores), new DenseMatrix(p, k, loadings), variance, ratio, means);
    }

    private static (double[] Loadings, double[] Eigen) FullSolve(double[] a, int n, int p, int k)
    {
        var loadings = new double[p * k];
        var eigOut = new double[k];
        if (p <= n)
        {
            var cov = AtA(a, n, p);
            var e = Decomposition.SymmetricEigen(new NdArray(cov, [p, p]));
            var vals = e.Values.ToArray(); var vecs = e.Vectors.ToArray();
            for (int c = 0; c < k; c++)
            {
                eigOut[c] = vals[c];
                for (int j = 0; j < p; j++) loadings[j * k + c] = vecs[j * p + c];
            }
        }
        else
        {
            // Gram trick: A Aᵀ = U S² Uᵀ, loadings V = Aᵀ U / s.
            var gram = AAt(a, n, p);
            var e = Decomposition.SymmetricEigen(new NdArray(gram, [n, n]));
            var vals = e.Values.ToArray(); var vecs = e.Vectors.ToArray();
            for (int c = 0; c < k; c++)
            {
                eigOut[c] = vals[c];
                double s = Math.Sqrt(Math.Max(vals[c], 1e-300));
                for (int j = 0; j < p; j++)
                {
                    double v = 0;
                    for (int i = 0; i < n; i++) v += a[i * p + j] * vecs[i * n + c];
                    loadings[j * k + c] = v / s;
                }
            }
        }
        return (loadings, eigOut);
    }

    private static (double[] Loadings, double[] Eigen) RandomizedSolve(double[] a, int n, int p, int k, int seed, int powerIterations)
    {
        int l = Math.Min(k + 10, Math.Min(n, p));
        var rng = new Random(seed);
        var omega = new double[p * l];
        for (int i = 0; i < omega.Length; i++) omega[i] = Gaussian(rng);
        var q = MatMul(a, n, p, omega, l);           // n × l
        Orthonormalize(q, n, l);
        for (int it = 0; it < powerIterations; it++)
        {
            var z = MatTMul(a, n, p, q, l);           // p × l = Aᵀ Q
            Orthonormalize(z, p, l);
            q = MatMul(a, n, p, z, l);
            Orthonormalize(q, n, l);
        }
        var b = MatTMul(q, n, l, a, p);               // l × p = Qᵀ A
        var bbt = AAt(b, l, p);                        // l × l
        var e = Decomposition.SymmetricEigen(new NdArray(bbt, [l, l]));
        var vals = e.Values.ToArray(); var vecs = e.Vectors.ToArray();
        var loadings = new double[p * k];
        var eig = new double[k];
        for (int c = 0; c < k; c++)
        {
            eig[c] = vals[c];
            double s = Math.Sqrt(Math.Max(vals[c], 1e-300));
            for (int j = 0; j < p; j++)
            {
                double v = 0;
                for (int r = 0; r < l; r++) v += b[r * p + j] * vecs[r * l + c];
                loadings[j * k + c] = v / s;
            }
        }
        return (loadings, eig);
    }

    internal static double Gaussian(Random rng)
    {
        double u1 = 1 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }

    /// <summary>C = A (n×p) · B (p×m).</summary>
    internal static double[] MatMul(double[] a, int n, int p, double[] b, int m)
    {
        var c = new double[n * m];
        Parallel.For(0, n, i =>
        {
            var ci = c.AsSpan(i * m, m);
            for (int t = 0; t < p; t++)
            {
                double av = a[i * p + t];
                if (av == 0) continue;
                var bt = b.AsSpan(t * m, m);
                for (int j = 0; j < m; j++) ci[j] += av * bt[j];
            }
        });
        return c;
    }

    /// <summary>C = Aᵀ (A is n×p) · B (n×m) → p×m.</summary>
    internal static double[] MatTMul(double[] a, int n, int p, double[] b, int m)
    {
        var c = new double[p * m];
        int chunks = Math.Max(1, Math.Min(Environment.ProcessorCount, n / 64));
        var partial = new double[chunks][];
        Parallel.For(0, chunks, ch =>
        {
            var acc = new double[p * m];
            int lo = (int)((long)n * ch / chunks), hi = (int)((long)n * (ch + 1) / chunks);
            for (int i = lo; i < hi; i++)
                for (int t = 0; t < p; t++)
                {
                    double av = a[i * p + t];
                    if (av == 0) continue;
                    for (int j = 0; j < m; j++) acc[t * m + j] += av * b[i * m + j];
                }
            partial[ch] = acc;
        });
        foreach (var acc in partial) for (int i = 0; i < c.Length; i++) c[i] += acc[i];
        return c;
    }

    private static double[] AtA(double[] a, int n, int p) => MatTMul(a, n, p, a, p);

    private static double[] AAt(double[] a, int n, int p)
    {
        var g = new double[n * n];
        Parallel.For(0, n, i =>
        {
            var ai = a.AsSpan(i * p, p);
            for (int j = 0; j <= i; j++)
            {
                double s = 0;
                var aj = a.AsSpan(j * p, p);
                for (int t = 0; t < p; t++) s += ai[t] * aj[t];
                g[i * n + j] = s; g[j * n + i] = s;
            }
        });
        return g;
    }

    /// <summary>Modified Gram–Schmidt on the columns of an n×m row-major matrix (twice for stability).</summary>
    private static void Orthonormalize(double[] q, int n, int m)
    {
        for (int pass = 0; pass < 2; pass++)
            for (int c = 0; c < m; c++)
            {
                for (int d = 0; d < c; d++)
                {
                    double dot = 0;
                    for (int i = 0; i < n; i++) dot += q[i * m + c] * q[i * m + d];
                    for (int i = 0; i < n; i++) q[i * m + c] -= dot * q[i * m + d];
                }
                double norm = 0;
                for (int i = 0; i < n; i++) norm += q[i * m + c] * q[i * m + c];
                norm = Math.Sqrt(norm);
                if (norm < 1e-300) continue;
                for (int i = 0; i < n; i++) q[i * m + c] /= norm;
            }
    }
}
