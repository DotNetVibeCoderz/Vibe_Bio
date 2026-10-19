namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>
/// One-dimensional port of R's <c>loess(y ~ x, span, degree = 2)</c> with the default <c>surface = "interpolate"</c>
/// (Cleveland &amp; Grosse's netlib dloess, R's loessf.f): a k-d tree whose cells hold at most
/// floor(n·span·cell) points, local quadratic fits (tricube weights over the floor(n·span) nearest points) at the
/// cell vertices giving a value and a slope, and cubic Hermite interpolation between vertices. Used by DESeq2's
/// simulation-based dispersion prior variance.
/// </summary>
internal static class RLoess
{
    public static double[] PredictInterpolate(double[] x, double[] y, double span, double[] z, double cell = 0.2)
    {
        int n = x.Length;
        int nf = Math.Min(n, (int)Math.Floor(n * span + 1e-5));
        int fc = (int)Math.Floor(n * (span * cell));
        int nvmax = Math.Max(200, n), ncmax = nvmax;
        var order = Enumerable.Range(0, n).OrderBy(i => x[i]).ToArray();
        var xs = order.Select(i => x[i]).ToArray();

        double alpha = xs.Min(), beta = xs.Max();
        double mu = 0.005 * Math.Max(beta - alpha, 1e-10 * Math.Max(Math.Abs(alpha), Math.Abs(beta)) + 1e-30);
        var vert = new List<double> { alpha - mu, beta + mu };
        // Cells (0-based): split value, children, vertex indices, point range [l, u] in sorted order (1-based like Fortran).
        var split = new List<bool>(); var xi = new List<double>(); var left = new List<int>(); var right = new List<int>();
        var cl = new List<int>(); var cu = new List<int>(); var lo = new List<int>(); var hi = new List<int>();
        void AddCell(int l, int u, int v0, int v1)
        {
            split.Add(false); xi.Add(0); left.Add(-1); right.Add(-1); cl.Add(v0); cu.Add(v1); lo.Add(l); hi.Add(u);
        }
        AddCell(1, n, 0, 1);
        for (int p = 0; p < split.Count; p++)
        {
            int l = lo[p], u = hi[p];
            bool leaf = u - l + 1 <= fc;
            if (!leaf) leaf = ncmax < split.Count + 2 || nvmax < vert.Count + 1.0;
            int m = 0;
            if (!leaf)
            {
                m = (l + u) / 2;
                // All ties go with the hi son (btyner's 2006 fix): search outwards for a non-tied split position.
                int offset = 0;
                while (!(m + offset >= u || m + offset < l))
                {
                    if (xs[m + offset - 1] == xs[m + offset])
                    {
                        offset = -offset;
                        if (offset >= 0) offset++;
                    }
                    else { m += offset; break; }
                }
                double t = xs[m - 1];
                leaf = vert[cl[p]] == t || vert[cu[p]] == t;
            }
            if (leaf) continue;
            double tv = xs[m - 1];
            split[p] = true;
            xi[p] = tv;
            int vnew = vert.IndexOf(tv);
            if (vnew < 0) { vert.Add(tv); vnew = vert.Count - 1; }
            left[p] = split.Count; AddCell(l, m, cl[p], vnew);
            right[p] = split.Count; AddCell(m + 1, u, vnew, cu[p]);
        }

        // Local quadratic fit (value, slope) at every vertex.
        var val = new double[vert.Count]; var slope = new double[vert.Count];
        var dist = new double[n];
        for (int v = 0; v < vert.Count; v++)
        {
            double q = vert[v];
            for (int i = 0; i < n; i++) dist[i] = (x[i] - q) * (x[i] - q);
            var sorted = (double[])dist.Clone();
            Array.Sort(sorted);
            double rho = sorted[nf - 1] * Math.Max(1.0, span);
            double s00 = 0, s01 = 0, s02 = 0, s03 = 0, s04 = 0, t0 = 0, t1 = 0, t2 = 0;
            for (int i = 0; i < n; i++)
            {
                if (dist[i] > sorted[nf - 1]) continue;
                double r = Math.Sqrt(dist[i] / rho);
                double w = 1 - r * r * r;
                w = w <= 0 ? 0 : w * w * w;
                if (w == 0) continue;
                double d = x[i] - q, d2 = d * d;
                s00 += w; s01 += w * d; s02 += w * d2; s03 += w * d2 * d; s04 += w * d2 * d2;
                t0 += w * y[i]; t1 += w * d * y[i]; t2 += w * d2 * y[i];
            }
            var sol = DESeqEngine.Solve([s00, s01, s02, s01, s02, s03, s02, s03, s04], [t0, t1, t2], 3);
            val[v] = sol?[0] ?? double.NaN;
            slope[v] = sol?[1] ?? double.NaN;
        }

        var res = new double[z.Length];
        for (int k = 0; k < z.Length; k++)
        {
            int c = 0;
            while (split[c]) c = z[k] <= xi[c] ? left[c] : right[c];
            double v0 = vert[cl[c]], v1 = vert[cu[c]];
            double h = (z[k] - v0) / (v1 - v0);
            double phi0 = (1 - h) * (1 - h) * (1 + 2 * h), phi1 = h * h * (3 - 2 * h);
            double psi0 = h * (1 - h) * (1 - h), psi1 = h * h * (h - 1);
            res[k] = phi0 * val[cl[c]] + phi1 * val[cu[c]] + (psi0 * slope[cl[c]] + psi1 * slope[cu[c]]) * (v1 - v0);
        }
        return res;
    }
}
