namespace Biotech.Net.Chemistry;

/// <summary>A stereo/tetrahedral constraint: signed volume of (N1−N4)·((N2−N4)×(N3−N4)) kept within [Lower, Upper].</summary>
internal readonly record struct ChiralSet(int Center, int N1, int N2, int N3, int N4, double Lower, double Upper, bool FusedSmallRings);

/// <summary>
/// Distance-geometry error function (RDKit <c>DistGeom::constructForceField</c>): flat-bottom distance violations
/// relative to the bounds, chiral-volume violations and a penalty on the fourth coordinate.
/// </summary>
internal sealed class DgForceField(int dim, List<(int I, int J, double L2, double U2, double W)> dist, List<ChiralSet> chiral, double chiralWeight, double fourthDimWeight)
{
    public int Dim { get; } = dim;

    public double Evaluate(ReadOnlySpan<double> x, Span<double> g)
    {
        bool grad = !g.IsEmpty;
        if (grad) g.Clear();
        double e = 0;
        foreach (var (i, j, l2, u2, w) in dist)
        {
            double d2 = FFMath.Distance2(x, i, j, Dim);
            if (d2 > u2)
            {
                double val = d2 / u2 - 1.0;
                e += w * val * val;
                if (grad)
                {
                    double d = Math.Sqrt(d2);
                    FFMath.AddDistanceGradient(x, g, i, j, Dim, d, 4.0 * w * val * d / u2);
                }
            }
            else if (d2 < l2)
            {
                double s = l2 + d2;
                double val = 2.0 * l2 / s - 1.0;
                e += w * val * val;
                if (grad)
                {
                    double d = Math.Sqrt(d2);
                    if (d < 1e-8) { d = 1e-8; }
                    FFMath.AddDistanceGradient(x, g, i, j, Dim, d, -8.0 * w * val * l2 * d / (s * s));
                }
            }
        }
        if (chiralWeight > 1e-8)
        {
            foreach (var c in chiral)
            {
                double v = FFMath.ChiralVolume(x, c.N1, c.N2, c.N3, c.N4, Dim);
                double dv = v < c.Lower ? v - c.Lower : v > c.Upper ? v - c.Upper : 0;
                if (dv == 0) continue;
                e += chiralWeight * dv * dv;
                if (grad) FFMath.AddChiralVolumeGradient(x, g, c.N1, c.N2, c.N3, c.N4, Dim, 2 * chiralWeight * dv);
            }
        }
        if (Dim == 4 && fourthDimWeight > 1e-8)
        {
            int n = x.Length / 4;
            for (int i = 0; i < n; i++)
            {
                double w4 = x[4 * i + 3];
                e += fourthDimWeight * w4 * w4;
                if (grad) g[4 * i + 3] += 2 * fourthDimWeight * w4;
            }
        }
        return e;
    }
}

/// <summary>
/// The 3-D "ETK" force field of ETKDG (RDKit <c>DistGeom::construct3DForceField</c>): experimental torsion terms,
/// sp² improper (inversion) terms, flat-bottom restraints that hold 1-2 and 1-3 distances, linear-angle restraints and
/// bounds-matrix restraints for all other pairs.
/// </summary>
internal sealed class EtkForceField
{
    public readonly List<(int I, int J, int K, int L, double[] V, int[] S)> Torsions = [];
    public readonly List<(int I, int J, int K, int L, double Kf)> Impropers = [];
    public readonly List<(int I, int J, double Min, double Max, double K)> Distances = [];
    public readonly List<(int I, int J, int K, double MinDeg, double MaxDeg, double Kf)> Angles = [];

    public double Evaluate(ReadOnlySpan<double> x, Span<double> g)
    {
        bool grad = !g.IsEmpty;
        if (grad) g.Clear();
        double e = 0;
        foreach (var (i, j, k, l, v, s) in Torsions)
        {
            double phi = FFMath.Dihedral(x, i, j, k, l, 3);
            double dE = 0;
            for (int m = 0; m < 6; m++)
            {
                if (v[m] == 0) continue;
                e += v[m] * (1 + s[m] * Math.Cos((m + 1) * phi));
                dE -= v[m] * s[m] * (m + 1) * Math.Sin((m + 1) * phi);
            }
            if (grad && dE != 0) FFMath.AddDihedralGradient(x, g, i, j, k, l, 3, dE);
        }
        foreach (var (i, j, k, l, kf) in Impropers)
        {
            double cosY = FFMath.InversionCosY(x, i, j, k, l, 3);
            double sinY2 = 1 - cosY * cosY;
            double sinY = sinY2 > 0 ? Math.Sqrt(sinY2) : 0;
            e += kf * (1 - sinY);
            if (grad) FFMath.AddInversionGradient(x, g, i, j, k, l, 3, kf, -1, 0);
        }
        foreach (var (i, j, min, max, kc) in Distances)
        {
            double d2 = FFMath.Distance2(x, i, j, 3);
            double diff;
            if (d2 < min * min) diff = Math.Sqrt(d2) - min;
            else if (d2 > max * max) diff = Math.Sqrt(d2) - max;
            else continue;
            e += 0.5 * kc * diff * diff;
            if (grad) FFMath.AddDistanceGradient(x, g, i, j, 3, Math.Sqrt(d2), kc * diff);
        }
        foreach (var (i, j, k, minDeg, maxDeg, kf) in Angles)
        {
            double c = FFMath.CosAngle(x, i, j, k, 3);
            double ang = Math.Acos(c) / FFMath.Deg2Rad;
            double t = ang < minDeg ? ang - minDeg : ang > maxDeg ? ang - maxDeg : 0;
            if (t == 0) continue;
            e += kf * t * t;
            if (grad)
            {
                double sinT = Math.Max(Math.Sqrt(1 - c * c), 1e-8);
                double dEdang = 2 * kf * t; // per degree
                double dEdc = dEdang / FFMath.Deg2Rad * (-1.0 / sinT);
                FFMath.AddCosAngleGradient(x, g, i, j, k, 3, dEdc);
            }
        }
        return e;
    }
}

/// <summary>Metric-matrix embedding helpers.</summary>
internal static class MetricMatrixEmbedding
{
    /// <summary>
    /// Embeds a distance matrix (row-major n×n) into <paramref name="dim"/> dimensions via the metric matrix and
    /// power-iteration eigenvectors (RDKit <c>DistGeom::computeInitialCoords</c>). Returns false for degenerate input.
    /// </summary>
    public static bool Embed(double[] dist, int n, int dim, Random rng, double[] coords)
    {
        const double eigTol = 0.001;
        var sq = new double[n * n];
        double sumSq = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double d2 = dist[i * n + j] * dist[i * n + j];
                sq[i * n + j] = sq[j * n + i] = d2;
                sumSq += d2;
            }
        sumSq /= (double)n * n;
        var d0 = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = 0;
            for (int j = 0; j < n; j++) s += sq[i * n + j];
            d0[i] = s / n - sumSq;
            if (d0[i] < eigTol && n > 3) return false;
        }
        var t = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
                t[i * n + j] = t[j * n + i] = 0.5 * (d0[i] + d0[j] - sq[i * n + j]);
        int nEig = Math.Min(dim, n);
        var values = new double[nEig];
        var vectors = new double[nEig][];
        PowerEigen(t, n, nEig, rng, values, vectors);
        int zero = 0;
        bool[] negative = new bool[nEig];
        for (int k = 0; k < nEig; k++)
        {
            if (values[k] > eigTol) values[k] = Math.Sqrt(values[k]);
            else if (Math.Abs(values[k]) < eigTol) { values[k] = 0; zero++; }
            else negative[k] = true;
        }
        if (zero >= 1 && n > 3) return false;
        for (int i = 0; i < n; i++)
            for (int k = 0; k < dim; k++)
                coords[i * dim + k] = k < nEig && !negative[k] ? values[k] * vectors[k][i] : 1.0 - 2.0 * rng.NextDouble();
        return true;
    }

    private static void PowerEigen(double[] a, int n, int count, Random rng, double[] values, double[][] vectors)
    {
        var m = (double[])a.Clone();
        var v = new double[n];
        var w = new double[n];
        for (int k = 0; k < count; k++)
        {
            for (int i = 0; i < n; i++) v[i] = rng.NextDouble() - 0.5;
            Normalize(v);
            double lambda = 0;
            for (int it = 0; it < 1000; it++)
            {
                for (int i = 0; i < n; i++)
                {
                    double s = 0;
                    int row = i * n;
                    for (int j = 0; j < n; j++) s += m[row + j] * v[j];
                    w[i] = s;
                }
                double newLambda = 0;
                for (int i = 0; i < n; i++) newLambda += v[i] * w[i];
                double norm = Normalize(w);
                if (norm < 1e-12) { newLambda = 0; break; }
                double change = 0;
                for (int i = 0; i < n; i++) change = Math.Max(change, Math.Abs(Math.Abs(w[i]) - Math.Abs(v[i])));
                Array.Copy(w, v, n);
                bool done = it > 3 && Math.Abs(newLambda - lambda) <= 1e-6 * Math.Max(1, Math.Abs(newLambda)) && change < 1e-5;
                lambda = newLambda;
                if (done) break;
            }
            values[k] = lambda;
            vectors[k] = (double[])v.Clone();
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m[i * n + j] -= lambda * v[i] * v[j];
        }
    }

    private static double Normalize(double[] v)
    {
        double s = 0;
        foreach (double x in v) s += x * x;
        s = Math.Sqrt(s);
        if (s > 0) for (int i = 0; i < v.Length; i++) v[i] /= s;
        return s;
    }
}
