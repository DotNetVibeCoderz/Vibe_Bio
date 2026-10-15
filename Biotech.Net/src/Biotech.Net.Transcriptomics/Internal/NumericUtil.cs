namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>Small numeric helpers shared by the transcriptomics algorithms (R-compatible where it matters).</summary>
internal static class NumericUtil
{
    /// <summary>Median of the finite values (NaN when none). Does not modify <paramref name="values"/>.</summary>
    public static double Median(ReadOnlySpan<double> values)
    {
        var buf = new double[values.Length];
        int n = 0;
        foreach (double v in values) if (!double.IsNaN(v)) buf[n++] = v;
        return MedianInPlace(buf.AsSpan(0, n));
    }

    /// <summary>Median of <paramref name="values"/>, sorting it in place.</summary>
    public static double MedianInPlace(Span<double> values)
    {
        int n = values.Length;
        if (n == 0) return double.NaN;
        values.Sort();
        return n % 2 == 1 ? values[n / 2] : 0.5 * (values[n / 2 - 1] + values[n / 2]);
    }

    /// <summary>R's default (type 7) sample quantile.</summary>
    public static double Quantile(ReadOnlySpan<double> values, double p)
    {
        var buf = values.ToArray();
        Array.Sort(buf);
        return QuantileSorted(buf, p);
    }

    /// <summary>R's type-7 quantile of already sorted data.</summary>
    public static double QuantileSorted(ReadOnlySpan<double> sorted, double p)
    {
        int n = sorted.Length;
        if (n == 0) return double.NaN;
        double h = (n - 1) * p;
        int lo = (int)Math.Floor(h);
        int hi = Math.Min(lo + 1, n - 1);
        return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
    }

    /// <summary>Median absolute deviation scaled by 1.4826 (R's <c>mad</c>).</summary>
    public static double Mad(ReadOnlySpan<double> values)
    {
        double med = Median(values);
        var dev = new double[values.Length];
        for (int i = 0; i < values.Length; i++) dev[i] = Math.Abs(values[i] - med);
        return 1.4826 * Median(dev);
    }

    /// <summary>Ranks with ties averaged (1-based, like R's <c>rank</c>).</summary>
    public static double[] Rank(ReadOnlySpan<double> values)
    {
        int n = values.Length;
        var idx = new int[n];
        var keys = values.ToArray();
        for (int i = 0; i < n; i++) idx[i] = i;
        Array.Sort(keys, idx);
        var ranks = new double[n];
        int s = 0;
        while (s < n)
        {
            int e = s;
            while (e + 1 < n && keys[e + 1] == keys[s]) e++;
            double r = 0.5 * (s + e) + 1;
            for (int k = s; k <= e; k++) ranks[idx[k]] = r;
            s = e + 1;
        }
        return ranks;
    }

    /// <summary>Digamma ψ(x) for x &gt; 0 (recurrence plus asymptotic series).</summary>
    public static double Digamma(double x)
    {
        if (double.IsNaN(x) || x <= 0) return double.NaN;
        double r = 0;
        while (x < 10) { r -= 1 / x; x += 1; }
        double f = 1 / (x * x);
        return r + Math.Log(x) - 0.5 / x
            - f * (1.0 / 12 - f * (1.0 / 120 - f * (1.0 / 252 - f * (1.0 / 240 - f * (1.0 / 132)))));
    }

    /// <summary>Trigamma ψ′(x) for x &gt; 0.</summary>
    public static double Trigamma(double x)
    {
        if (double.IsNaN(x) || x <= 0) return double.NaN;
        double r = 0;
        while (x < 10) { r += 1 / (x * x); x += 1; }
        double f = 1 / (x * x);
        return r + 1 / x + f / 2 + f / x * (1.0 / 6 - f * (1.0 / 30 - f * (1.0 / 42 - f * (1.0 / 30))));
    }

    /// <summary>Tetragamma ψ″(x) for x &gt; 0.</summary>
    public static double Tetragamma(double x)
    {
        if (double.IsNaN(x) || x <= 0) return double.NaN;
        double r = 0;
        while (x < 10) { r -= 2 / (x * x * x); x += 1; }
        double x2 = x * x;
        double f = 1 / x2;
        return r - 1 / x2 - 1 / (x2 * x) - 0.5 * f * f
            + f * f * f * (1.0 / 6 - f * (1.0 / 6 - f * (3.0 / 10 - f * (5.0 / 6))));
    }

    /// <summary>Inverse of the trigamma function (limma's <c>trigammaInverse</c>, Newton iteration).</summary>
    public static double TrigammaInverse(double x)
    {
        if (double.IsNaN(x)) return double.NaN;
        if (x > 1e7) return 1 / Math.Sqrt(x);
        if (x < 1e-6) return 1 / x;
        double y = 0.5 + 1 / x;
        for (int iter = 0; iter < 50; iter++)
        {
            double tri = Trigamma(y);
            double dif = tri * (1 - tri / x) / Tetragamma(y);
            y += dif;
            if (-dif / y < 1e-8) break;
        }
        return y;
    }

    /// <summary>
    /// Port of R's <c>lowess</c> (Cleveland 1979, the <c>clowess</c> C routine) with <paramref name="iterations"/>
    /// robustness steps. <paramref name="x"/> must be sorted ascending; returns the fitted values.
    /// </summary>
    public static double[] Lowess(ReadOnlySpan<double> x, ReadOnlySpan<double> y, double f = 2.0 / 3, int iterations = 3, double? delta = null)
    {
        int n = x.Length;
        var ys = new double[n];
        if (n == 0) return ys;
        if (n < 2) { ys[0] = y[0]; return ys; }
        double d = delta ?? 0.01 * (x[n - 1] - x[0]);
        var rw = new double[n];
        var res = new double[n];
        int ns = Math.Max(2, Math.Min(n, (int)(f * n + 1e-7)));
        for (int iter = 1; iter <= iterations + 1; iter++)
        {
            int nleft = 0, nright = ns - 1, last = -1, i = 0;
            while (true)
            {
                if (nright < n - 1)
                {
                    double d1 = x[i] - x[nleft], d2 = x[nright + 1] - x[i];
                    if (d1 > d2) { nleft++; nright++; continue; }
                }
                bool ok = LowessFit(x, y, x[i], out ys[i], nleft, nright, res, iter > 1, rw);
                if (!ok) ys[i] = y[i];
                if (last < i - 1)
                {
                    double denom = x[i] - x[last];
                    for (int j = last + 1; j < i; j++)
                    {
                        double alpha = (x[j] - x[last]) / denom;
                        ys[j] = alpha * ys[i] + (1 - alpha) * ys[last];
                    }
                }
                last = i;
                double cut = x[last] + d;
                for (i = last + 1; i < n; i++)
                {
                    if (x[i] > cut) break;
                    if (x[i] == x[last]) { ys[i] = ys[last]; last = i; }
                }
                i = Math.Max(last + 1, i - 1);
                if (last >= n - 1) break;
            }
            for (int k = 0; k < n; k++) res[k] = y[k] - ys[k];
            double sc = 0;
            for (int k = 0; k < n; k++) sc += Math.Abs(res[k]);
            sc /= n;
            if (iter > iterations) break;
            for (int k = 0; k < n; k++) rw[k] = Math.Abs(res[k]);
            var sorted = (double[])rw.Clone();
            Array.Sort(sorted);
            int m1 = n / 2;
            double cmad = n % 2 == 0 ? 3 * (sorted[m1] + sorted[n - m1 - 1]) : 6 * sorted[m1];
            if (cmad < 1e-7 * sc) break;
            double c9 = 0.999 * cmad, c1 = 0.001 * cmad;
            for (int k = 0; k < n; k++)
            {
                double r = Math.Abs(res[k]);
                if (r <= c1) rw[k] = 1;
                else if (r <= c9) { double u = r / cmad; u = 1 - u * u; rw[k] = u * u; }
                else rw[k] = 0;
            }
        }
        return ys;
    }

    private static bool LowessFit(ReadOnlySpan<double> x, ReadOnlySpan<double> y, double xs, out double ys,
        int nleft, int nright, double[] w, bool userw, double[] rw)
    {
        int n = x.Length;
        double range = x[n - 1] - x[0];
        double h = Math.Max(xs - x[nleft], x[nright] - xs);
        double h9 = 0.999 * h, h1 = 0.001 * h;
        double a = 0;
        int j = nleft;
        while (j < n)
        {
            w[j] = 0;
            double r = Math.Abs(x[j] - xs);
            if (r <= h9)
            {
                if (r <= h1) w[j] = 1;
                else { double q = r / h; q = 1 - q * q * q; w[j] = q * q * q; }
                if (userw) w[j] *= rw[j];
                a += w[j];
            }
            else if (x[j] > xs) break;
            j++;
        }
        int nrt = j - 1;
        ys = 0;
        if (a <= 0) return false;
        for (j = nleft; j <= nrt; j++) w[j] /= a;
        if (h > 0)
        {
            a = 0;
            for (j = nleft; j <= nrt; j++) a += w[j] * x[j];
            double b = xs - a, c = 0;
            for (j = nleft; j <= nrt; j++) c += w[j] * (x[j] - a) * (x[j] - a);
            if (Math.Sqrt(c) > 0.001 * range)
            {
                b /= c;
                for (j = nleft; j <= nrt; j++) w[j] *= b * (x[j] - a) + 1;
            }
        }
        for (j = nleft; j <= nrt; j++) ys += w[j] * y[j];
        return true;
    }

    /// <summary>
    /// Local polynomial regression (loess, tricube weights, no robustness iterations) evaluated at every x.
    /// Uses exact local fits rather than the kd-tree interpolation of R's default <c>surface = "interpolate"</c>,
    /// so values can differ slightly from R/skmisc in the last digits.
    /// </summary>
    public static double[] Loess(ReadOnlySpan<double> x, ReadOnlySpan<double> y, double span, int degree)
    {
        int n = x.Length;
        var order = new int[n];
        var xs = x.ToArray();
        for (int i = 0; i < n; i++) order[i] = i;
        Array.Sort(xs, order);
        var ysorted = new double[n];
        for (int i = 0; i < n; i++) ysorted[i] = y[order[i]];
        int q = Math.Clamp((int)Math.Floor(span * n), degree + 1, n);
        var fittedSorted = new double[n];
        var xsArr = xs;
        Parallel.For(0, n, () => new double[(degree + 1) * (degree + 2)], (i, _, buf) =>
        {
            double x0 = xsArr[i];
            // Find the q nearest neighbours: a window [lo, lo+q) in sorted order.
            int lo = Math.Clamp(i - q / 2, 0, n - q);
            while (lo > 0 && x0 - xsArr[lo - 1] < xsArr[lo + q - 1] - x0) lo--;
            while (lo + q < n && xsArr[lo + q] - x0 < x0 - xsArr[lo]) lo++;
            double h = Math.Max(x0 - xsArr[lo], xsArr[lo + q - 1] - x0);
            h = Math.Max(h, 1e-12);
            int p = degree + 1;
            Array.Clear(buf);
            Span<double> pw = stackalloc double[3];
            // Normal equations of weighted polynomial regression centred at x0.
            for (int j = lo; j < lo + q; j++)
            {
                double u = Math.Abs(xsArr[j] - x0) / h;
                if (u >= 1) continue;
                double w = 1 - u * u * u; w = w * w * w;
                double dx = xsArr[j] - x0;
                pw[0] = 1; pw[1] = dx; pw[2] = dx * dx;
                for (int r = 0; r < p; r++)
                {
                    for (int c = 0; c < p; c++) buf[r * p + c] += w * pw[r] * pw[c];
                    buf[p * p + r] += w * pw[r] * ysorted[j];
                }
            }
            var sol = SolveSmall(buf.AsSpan(0, p * p).ToArray(), buf.AsSpan(p * p, p).ToArray(), p);
            fittedSorted[i] = sol is null ? ysorted[i] : sol[0];
            return buf;
        }, _ => { });
        var fitted = new double[n];
        for (int i = 0; i < n; i++) fitted[order[i]] = fittedSorted[i];
        return fitted;
    }

    /// <summary>Solves a small dense system A x = b (Gaussian elimination with partial pivoting); null if singular.</summary>
    public static double[]? SolveSmall(double[] a, double[] b, int n)
    {
        for (int col = 0; col < n; col++)
        {
            int piv = col;
            double best = Math.Abs(a[col * n + col]);
            for (int r = col + 1; r < n; r++)
            {
                double v = Math.Abs(a[r * n + col]);
                if (v > best) { best = v; piv = r; }
            }
            if (best < 1e-300) return null;
            if (piv != col)
            {
                for (int c = 0; c < n; c++) (a[col * n + c], a[piv * n + c]) = (a[piv * n + c], a[col * n + c]);
                (b[col], b[piv]) = (b[piv], b[col]);
            }
            for (int r = col + 1; r < n; r++)
            {
                double fct = a[r * n + col] / a[col * n + col];
                if (fct == 0) continue;
                for (int c = col; c < n; c++) a[r * n + c] -= fct * a[col * n + c];
                b[r] -= fct * b[col];
            }
        }
        var x = new double[n];
        for (int r = n - 1; r >= 0; r--)
        {
            double s = b[r];
            for (int c = r + 1; c < n; c++) s -= a[r * n + c] * x[c];
            x[r] = s / a[r * n + r];
        }
        return x;
    }

    /// <summary>Inverts a small symmetric positive-definite matrix (Cholesky); null if not positive definite.</summary>
    public static double[]? InvertSpd(ReadOnlySpan<double> a, int n)
    {
        var l = new double[n * n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s = a[i * n + j];
                for (int k = 0; k < j; k++) s -= l[i * n + k] * l[j * n + k];
                if (i == j)
                {
                    if (s <= 0) return null;
                    l[i * n + i] = Math.Sqrt(s);
                }
                else l[i * n + j] = s / l[j * n + j];
            }
        }
        // inv(L)
        var li = new double[n * n];
        for (int i = 0; i < n; i++)
        {
            li[i * n + i] = 1 / l[i * n + i];
            for (int j = 0; j < i; j++)
            {
                double s = 0;
                for (int k = j; k < i; k++) s -= l[i * n + k] * li[k * n + j];
                li[i * n + j] = s / l[i * n + i];
            }
        }
        var inv = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = 0;
                for (int k = i; k < n; k++) s += li[k * n + i] * li[k * n + j];
                inv[i * n + j] = s; inv[j * n + i] = s;
            }
        return inv;
    }

    /// <summary>log det of a small SPD matrix via Cholesky (NaN when not positive definite).</summary>
    public static double LogDetSpd(ReadOnlySpan<double> a, int n)
    {
        Span<double> l = n <= 16 ? stackalloc double[n * n] : new double[n * n];
        double logDet = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i * n + j];
                for (int k = 0; k < j; k++) s -= l[i * n + k] * l[j * n + k];
                if (i == j)
                {
                    if (s <= 0) return double.NaN;
                    l[i * n + i] = Math.Sqrt(s);
                    logDet += Math.Log(s);
                }
                else l[i * n + j] = s / l[j * n + j];
            }
        return logDet;
    }

    /// <summary>Linear interpolation of sorted knots (<c>approxfun(rule = 2)</c>: constant beyond the ends).</summary>
    public static double Interpolate(ReadOnlySpan<double> xKnots, ReadOnlySpan<double> yKnots, double x)
    {
        int n = xKnots.Length;
        if (x <= xKnots[0]) return yKnots[0];
        if (x >= xKnots[n - 1]) return yKnots[n - 1];
        int lo = 0, hi = n - 1;
        while (hi - lo > 1)
        {
            int mid = (lo + hi) >> 1;
            if (xKnots[mid] <= x) lo = mid; else hi = mid;
        }
        double dx = xKnots[hi] - xKnots[lo];
        return dx <= 0 ? yKnots[lo] : yKnots[lo] + (x - xKnots[lo]) / dx * (yKnots[hi] - yKnots[lo]);
    }

    /// <summary>Upper tail of the standard normal, accurate far into the tail.</summary>
    public static double NormalUpper(double z) => 0.5 * Biotech.Net.Native.SpecialFunctions.Erfc(z / Math.Sqrt(2));
}
