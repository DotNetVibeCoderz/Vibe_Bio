namespace Biotech.Net.Chemistry;

/// <summary>Energy function: returns E(x) and writes the full gradient into <c>grad</c>.</summary>
internal delegate double EnergyFunction(ReadOnlySpan<double> x, Span<double> grad);

/// <summary>Outcome of a local minimisation.</summary>
/// <param name="Energy">Final energy.</param>
/// <param name="Iterations">Iterations used.</param>
/// <param name="Converged">True when the gradient or energy-change criterion was met before the iteration cap.</param>
public readonly record struct MinimizationResult(double Energy, int Iterations, bool Converged);

/// <summary>
/// Limited-memory BFGS (Nocedal &amp; Wright, Numerical Optimization, alg. 7.4) with a backtracking Armijo line search
/// and a cap on the largest single-coordinate move per step.
/// </summary>
internal static class Lbfgs
{
    public static MinimizationResult Minimize(double[] x, EnergyFunction f, int maxIterations, double gradientTolerance = 1e-4,
        double energyTolerance = 1e-9, int memory = 8, double maxStep = 0.3)
    {
        int n = x.Length;
        var g = new double[n];
        var xNew = new double[n];
        var gNew = new double[n];
        var d = new double[n];
        var s = new double[memory][];
        var y = new double[memory][];
        var rho = new double[memory];
        var alpha = new double[memory];
        for (int i = 0; i < memory; i++) { s[i] = new double[n]; y[i] = new double[n]; }
        int stored = 0, head = 0;
        double e = f(x, g);
        if (!double.IsFinite(e)) return new MinimizationResult(e, 0, false);
        int smallSteps = 0;
        for (int iter = 1; iter <= maxIterations; iter++)
        {
            if (MaxAbs(g) < gradientTolerance) return new MinimizationResult(e, iter - 1, true);
            // two-loop recursion: d = −H g
            for (int k = 0; k < n; k++) d[k] = -g[k];
            for (int t = 0; t < stored; t++)
            {
                int idx = (head - 1 - t + memory) % memory;
                alpha[idx] = rho[idx] * Dot(s[idx], d);
                Axpy(-alpha[idx], y[idx], d);
            }
            if (stored > 0)
            {
                int last = (head - 1 + memory) % memory;
                double gamma = Dot(s[last], y[last]) / Dot(y[last], y[last]);
                for (int k = 0; k < n; k++) d[k] *= gamma;
            }
            for (int t = stored - 1; t >= 0; t--)
            {
                int idx = (head - 1 - t + memory) % memory;
                double beta = rho[idx] * Dot(y[idx], d);
                Axpy(alpha[idx] - beta, s[idx], d);
            }
            double slope = Dot(d, g);
            if (!(slope < 0))
            {
                // not a descent direction: reset to steepest descent
                stored = 0;
                for (int k = 0; k < n; k++) d[k] = -g[k];
                slope = Dot(d, g);
            }
            double big = MaxAbs(d);
            double step = big > maxStep ? maxStep / big : 1.0;
            double eNew = double.NaN;
            bool accepted = false;
            for (int ls = 0; ls < 30; ls++)
            {
                for (int k = 0; k < n; k++) xNew[k] = x[k] + step * d[k];
                eNew = f(xNew, gNew);
                if (double.IsFinite(eNew) && eNew <= e + 1e-4 * step * slope) { accepted = true; break; }
                step *= 0.5;
            }
            if (!accepted)
            {
                if (stored == 0) return new MinimizationResult(e, iter, MaxAbs(g) < gradientTolerance * 10);
                stored = 0; // retry with steepest descent
                continue;
            }
            var sv = s[head];
            var yv = y[head];
            for (int k = 0; k < n; k++) { sv[k] = xNew[k] - x[k]; yv[k] = gNew[k] - g[k]; }
            double sy = Dot(sv, yv);
            if (sy > 1e-12)
            {
                rho[head] = 1 / sy;
                head = (head + 1) % memory;
                if (stored < memory) stored++;
            }
            double de = e - eNew;
            Array.Copy(xNew, x, n);
            Array.Copy(gNew, g, n);
            e = eNew;
            if (de < energyTolerance * Math.Max(1, Math.Abs(e)))
            {
                if (++smallSteps >= 3) return new MinimizationResult(e, iter, true);
            }
            else smallSteps = 0;
        }
        return new MinimizationResult(e, maxIterations, MaxAbs(g) < gradientTolerance);
    }

    private static double Dot(double[] a, double[] b)
    {
        double s = 0;
        for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
        return s;
    }

    private static void Axpy(double a, double[] x, double[] y)
    {
        for (int i = 0; i < x.Length; i++) y[i] += a * x[i];
    }

    private static double MaxAbs(double[] a)
    {
        double m = 0;
        foreach (double v in a) { double t = Math.Abs(v); if (t > m) m = t; }
        return m;
    }
}
