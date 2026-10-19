namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>
/// Port of the L-BFGS solver of LBFGSpp (Yixuan Qiu) as configured by RcppNumerical's <c>optim_lbfgs</c>, which is
/// what apeglm's C++ routine <c>nbinomGLM</c> uses: m = 6 corrections, Nocedal–Wright strong-Wolfe line search
/// (ftol 1e-4, wolfe 0.9, 100 trials), gradient tolerance <c>epsilon = epsilon_rel = eps_g</c> and a one-step
/// relative function-change test <c>delta = eps_f</c>.
/// </summary>
internal sealed class Lbfgs
{
    /// <summary>Objective: returns f(x) and writes the gradient.</summary>
    public delegate double Objective(double[] x, double[] grad);

    private readonly int _m = 6;
    private readonly double _epsilon, _epsilonRel, _delta, _ftol = 1e-4, _wolfe = 0.9;
    private readonly int _maxIterations, _maxLinesearch = 100, _past = 1;

    public Lbfgs(int maxIterations, double epsF, double epsG)
    {
        _maxIterations = maxIterations;
        _epsilon = epsG; _epsilonRel = epsG; _delta = epsF;
    }

    /// <summary>RcppNumerical's <c>optim_lbfgs</c>: minimises in place; returns 0, or -1 when the solver threw.</summary>
    public int Minimize(Objective f, double[] x, out double fx)
    {
        try
        {
            fx = MinimizeCore(f, x);
            return 0;
        }
        catch (LbfgsException e)
        {
            // As in C++, x (updated in place) and fx hold the solver state at the moment it threw.
            fx = e.Fx;
            return -1;
        }
    }

    private sealed class LbfgsException(double fx) : Exception
    {
        public double Fx { get; } = fx;
    }

    private static double Dot(double[] a, double[] b) { double s = 0; for (int i = 0; i < a.Length; i++) s += a[i] * b[i]; return s; }
    private static double Norm(double[] a) => Math.Sqrt(Dot(a, a));

    private double MinimizeCore(Objective f, double[] x)
    {
        int n = x.Length;
        var s = new double[_m][]; var y = new double[_m][];
        for (int i = 0; i < _m; i++) { s[i] = new double[n]; y[i] = new double[n]; }
        var ys = new double[_m]; var alpha = new double[_m];
        double theta = 1;
        int ncorr = 0, ptr = _m;
        var fxHist = new double[_past];
        var grad = new double[n];
        var xp = new double[n]; var gradp = new double[n]; var drt = new double[n];
        var vecs = new double[n]; var vecy = new double[n];

        double fx = f(x, grad);
        double gnorm = Norm(grad);
        fxHist[0] = fx;
        if (gnorm <= _epsilon || gnorm <= _epsilonRel * Norm(x)) return fx;
        for (int i = 0; i < n; i++) drt[i] = -grad[i];
        double step = 1.0 / Norm(drt);
        const double eps = 2.220446049250313e-16;
        int k = 1;
        for (; ; )
        {
            Array.Copy(x, xp, n);
            Array.Copy(grad, gradp, n);
            double dg = Dot(grad, drt);
            LineSearch(f, xp, drt, ref step, ref fx, grad, ref dg, x);
            gnorm = Norm(grad);
            if (gnorm <= _epsilon || gnorm <= _epsilonRel * Norm(x)) return fx;
            double fxd = fxHist[k % _past];
            if (k >= _past && Math.Abs(fxd - fx) <= _delta * Math.Max(Math.Max(Math.Abs(fx), Math.Abs(fxd)), 1.0)) return fx;
            fxHist[k % _past] = fx;
            if (_maxIterations != 0 && k >= _maxIterations) return fx;
            for (int i = 0; i < n; i++) { vecs[i] = x[i] - xp[i]; vecy[i] = grad[i] - gradp[i]; }
            if (Dot(vecs, vecy) > eps * Dot(vecy, vecy))
            {
                int loc = ptr % _m;
                Array.Copy(vecs, s[loc], n); Array.Copy(vecy, y[loc], n);
                double ysv = Dot(s[loc], y[loc]);
                ys[loc] = ysv;
                theta = Dot(y[loc], y[loc]) / ysv;
                if (ncorr < _m) ncorr++;
                ptr = loc + 1;
            }
            // drt = -H * grad (two-loop recursion)
            for (int i = 0; i < n; i++) drt[i] = -grad[i];
            int j = ptr % _m;
            for (int i = 0; i < ncorr; i++)
            {
                j = (j + _m - 1) % _m;
                alpha[j] = Dot(s[j], drt) / ys[j];
                for (int t = 0; t < n; t++) drt[t] -= alpha[j] * y[j][t];
            }
            for (int t = 0; t < n; t++) drt[t] /= theta;
            for (int i = 0; i < ncorr; i++)
            {
                double beta = Dot(y[j], drt) / ys[j];
                for (int t = 0; t < n; t++) drt[t] += (alpha[j] - beta) * s[j][t];
                j = (j + 1) % _m;
            }
            step = 1.0;
            k++;
        }
    }

    private static double QuadInterp(double stepLo, double stepHi, double fxLo, double fxHi, double dgLo)
    {
        double fdiff = fxHi - fxLo, sdiff = stepHi - stepLo, smid = (stepHi + stepLo) / 2;
        double cand = fdiff * stepLo - smid * sdiff * dgLo;
        cand /= fdiff - sdiff * dgLo;
        bool nan = !double.IsFinite(cand);
        double endDist = Math.Min(Math.Abs(cand - stepLo), Math.Abs(cand - stepHi));
        bool nearEnd = endDist < 0.01 * Math.Abs(sdiff);
        bool bisect = nan || cand <= Math.Min(stepLo, stepHi) || cand >= Math.Max(stepLo, stepHi) || nearEnd;
        return bisect ? smid : cand;
    }

    private void LineSearch(Objective f, double[] xp, double[] drt, ref double step, ref double fx, double[] grad, ref double dg, double[] x)
    {
        int n = x.Length;
        if (step <= 0) throw new LbfgsException(fx);
        double fxInit = fx, dgInit = dg;
        if (dgInit > 0) throw new LbfgsException(fx);
        double testDecr = _ftol * dgInit, testCurv = -_wolfe * dgInit;
        double stepHi, fxHi;
        double stepLo = 0, fxLo = fxInit, dgLo = dgInit;
        var xLo = (double[])xp.Clone(); var gradLo = (double[])grad.Clone();
        int iter = 0;
        for (; ; )
        {
            for (int i = 0; i < n; i++) x[i] = xp[i] + step * drt[i];
            fx = f(x, grad);
            dg = Dot(grad, drt);
            if (fx - fxInit > step * testDecr || (0 < stepLo && fx >= fxLo))
            {
                stepHi = step; fxHi = fx;
                break;
            }
            if (Math.Abs(dg) <= testCurv) return;
            stepHi = stepLo; fxHi = fxLo;
            stepLo = step; fxLo = fx; dgLo = dg;
            Swap(xLo, x); Swap(gradLo, grad);
            if (dg >= 0) break;
            iter++;
            if (iter >= _maxLinesearch) { Swap(x, xLo); Swap(grad, gradLo); return; }
            step *= 2;
        }
        for (; ; )
        {
            step = QuadInterp(stepLo, stepHi, fxLo, fxHi, dgLo);
            for (int i = 0; i < n; i++) x[i] = xp[i] + step * drt[i];
            fx = f(x, grad);
            dg = Dot(grad, drt);
            if (fx - fxInit > step * testDecr || fx >= fxLo)
            {
                if (step == stepHi) throw new LbfgsException(fx);
                stepHi = step; fxHi = fx;
            }
            else
            {
                if (Math.Abs(dg) <= testCurv) return;
                if (dg * (stepHi - stepLo) >= 0) { stepHi = stepLo; fxHi = fxLo; }
                if (step == stepLo) throw new LbfgsException(fx);
                stepLo = step; fxLo = fx; dgLo = dg;
                Swap(xLo, x); Swap(gradLo, grad);
            }
            iter++;
            if (iter >= _maxLinesearch)
            {
                if (stepLo <= 0) throw new LbfgsException(fx);
                step = stepLo; fx = fxLo; dg = dgLo;
                Swap(x, xLo); Swap(grad, gradLo);
                return;
            }
        }
    }

    // Eigen's swap exchanges contents; emulate by copying (arrays are owned by the caller).
    private static void Swap(double[] a, double[] b)
    {
        for (int i = 0; i < a.Length; i++) (a[i], b[i]) = (b[i], a[i]);
    }
}
