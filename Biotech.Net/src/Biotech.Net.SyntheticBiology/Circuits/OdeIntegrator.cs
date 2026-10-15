namespace Biotech.Net.SyntheticBiology.Circuits;

/// <summary>
/// The right-hand side of an explicit first-order ordinary differential equation system
/// dy/dt = f(t, y), as used for deterministic models of gene regulatory circuits.
/// </summary>
/// <param name="t">Independent variable, normally time in the model's own units.</param>
/// <param name="y">Current state vector; must not be modified by the implementation.</param>
/// <param name="dydt">Output buffer for the derivatives; it has the same length as <paramref name="y"/> and must be fully written.</param>
public delegate void OdeSystem(double t, ReadOnlySpan<double> y, Span<double> dydt);

/// <summary>
/// The result of integrating an <see cref="OdeSystem"/>: the recorded time grid and the state vector
/// at each recorded time.
/// </summary>
/// <param name="Times">Recorded times, strictly increasing, starting at the initial time.</param>
/// <param name="States">Recorded states; <c>States[k]</c> is the state vector at <c>Times[k]</c>.</param>
public sealed record OdeSolution(double[] Times, double[][] States)
{
    /// <summary>Number of recorded output points (not necessarily the number of internal steps taken).</summary>
    public int StepCount => Times.Length;

    /// <summary>The state vector at the last recorded time.</summary>
    public double[] Final => States[^1];

    /// <summary>
    /// Number of adaptive steps whose local error estimate exceeded the tolerance and were rejected
    /// and retried with a smaller step. Always zero for fixed-step integration.
    /// </summary>
    public int RejectedSteps { get; init; }

    /// <summary>
    /// The time series of a single state variable, i.e. <c>States[k][index]</c> for every recorded k.
    /// </summary>
    /// <param name="index">Zero-based index of the state variable.</param>
    public double[] SeriesOf(int index)
    {
        if (index < 0 || index >= States[0].Length)
            throw new ArgumentOutOfRangeException(nameof(index), index, $"The system has {States[0].Length} state variables.");
        var series = new double[States.Length];
        for (int k = 0; k < States.Length; k++) series[k] = States[k][index];
        return series;
    }

    /// <summary>
    /// The state at an arbitrary time, by linear interpolation between the two bracketing recorded
    /// points. Times outside the recorded interval are clamped to the first or last recorded state.
    /// Linear interpolation is first-order accurate only, so output the solution densely if an
    /// accurate value between steps matters.
    /// </summary>
    /// <param name="time">The time at which to evaluate the solution.</param>
    public double[] InterpolateAt(double time)
    {
        if (time <= Times[0]) return (double[])States[0].Clone();
        if (time >= Times[^1]) return (double[])States[^1].Clone();

        int low = 0, high = Times.Length - 1;
        while (high - low > 1)
        {
            int mid = (low + high) / 2;
            if (Times[mid] <= time) low = mid; else high = mid;
        }
        double span = Times[high] - Times[low];
        double w = span <= 0.0 ? 0.0 : (time - Times[low]) / span;
        var result = new double[States[low].Length];
        for (int i = 0; i < result.Length; i++)
            result[i] = States[low][i] + w * (States[high][i] - States[low][i]);
        return result;
    }
}

/// <summary>
/// Explicit Runge-Kutta integrators for the small, smooth, non-stiff ODE systems that arise in
/// deterministic models of genetic circuits: a classical fixed-step fourth-order method and an
/// adaptive Dormand-Prince 5(4) method.
/// <para>
/// Both are explicit, so a stiff system (very fast degradation or binding relative to the simulated
/// horizon) will force the adaptive method to take extremely small steps; in that case reduce the
/// horizon or use a stiff solver.
/// </para>
/// </summary>
public static class OdeIntegrator
{
    /// <summary>
    /// Classical fourth-order Runge-Kutta with a fixed step. Local truncation error is O(h^5) and the
    /// global error O(h^4), so halving the step reduces the error about sixteen-fold. The method is
    /// cheap and perfectly reproducible but has no error control.
    /// </summary>
    /// <param name="f">The system dy/dt = f(t, y).</param>
    /// <param name="y0">Initial state at <paramref name="t0"/>; copied, not modified.</param>
    /// <param name="t0">Initial time.</param>
    /// <param name="t1">Final time; must be greater than <paramref name="t0"/>.</param>
    /// <param name="step">Step size h in time units; must be positive. The last step is shortened to land exactly on <paramref name="t1"/>.</param>
    /// <param name="outputEvery">Record the state every n-th step; the initial and final states are always recorded.</param>
    public static OdeSolution Rk4(OdeSystem f, double[] y0, double t0, double t1, double step, int outputEvery = 1)
    {
        ArgumentNullException.ThrowIfNull(f);
        ArgumentNullException.ThrowIfNull(y0);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(step);
        ArgumentOutOfRangeException.ThrowIfLessThan(outputEvery, 1);
        if (t1 <= t0) throw new ArgumentOutOfRangeException(nameof(t1), t1, "The final time must be greater than the initial time.");

        int n = y0.Length;
        int steps = (int)Math.Ceiling((t1 - t0) / step);
        if (steps < 1) steps = 1;

        var times = new List<double>(steps / outputEvery + 2);
        var states = new List<double[]>(steps / outputEvery + 2);

        var y = (double[])y0.Clone();
        var k1 = new double[n];
        var k2 = new double[n];
        var k3 = new double[n];
        var k4 = new double[n];
        var temp = new double[n];

        double t = t0;
        times.Add(t);
        states.Add((double[])y.Clone());

        for (int s = 0; s < steps; s++)
        {
            double h = Math.Min(step, t1 - t);
            if (h <= 0.0) break;

            f(t, y, k1);
            for (int i = 0; i < n; i++) temp[i] = y[i] + 0.5 * h * k1[i];
            f(t + 0.5 * h, temp, k2);
            for (int i = 0; i < n; i++) temp[i] = y[i] + 0.5 * h * k2[i];
            f(t + 0.5 * h, temp, k3);
            for (int i = 0; i < n; i++) temp[i] = y[i] + h * k3[i];
            f(t + h, temp, k4);
            for (int i = 0; i < n; i++)
                y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

            t += h;
            bool last = s == steps - 1 || t >= t1 - 1e-15;
            if ((s + 1) % outputEvery == 0 || last)
            {
                times.Add(t);
                states.Add((double[])y.Clone());
            }
            if (last) break;
        }

        return new OdeSolution([.. times], [.. states]);
    }

    // Dormand-Prince 5(4) coefficients (Dormand & Prince, J. Comput. Appl. Math. 6:19-26, 1980).
    private static readonly double[] C = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];

    private const double A21 = 1.0 / 5.0;
    private const double A31 = 3.0 / 40.0, A32 = 9.0 / 40.0;
    private const double A41 = 44.0 / 45.0, A42 = -56.0 / 15.0, A43 = 32.0 / 9.0;
    private const double A51 = 19372.0 / 6561.0, A52 = -25360.0 / 2187.0, A53 = 64448.0 / 6561.0, A54 = -212.0 / 729.0;
    private const double A61 = 9017.0 / 3168.0, A62 = -355.0 / 33.0, A63 = 46732.0 / 5247.0, A64 = 49.0 / 176.0, A65 = -5103.0 / 18656.0;
    private const double B1 = 35.0 / 384.0, B3 = 500.0 / 1113.0, B4 = 125.0 / 192.0, B5 = -2187.0 / 6784.0, B6 = 11.0 / 84.0;
    private const double E1 = 5179.0 / 57600.0, E3 = 7571.0 / 16695.0, E4 = 393.0 / 640.0, E5 = -92097.0 / 339200.0, E6 = 187.0 / 2100.0, E7 = 1.0 / 40.0;

    /// <summary>
    /// Adaptive Dormand-Prince 5(4) integration (the DOPRI5 tableau, as used by MATLAB's <c>ode45</c>
    /// and SciPy's <c>RK45</c>). Each step evaluates seven stages, advances with the fifth-order
    /// solution, and estimates the local error from the difference against the embedded fourth-order
    /// solution. Steps whose scaled error exceeds one are rejected and retried; the step size is then
    /// updated by h_new = h * clamp(0.9 * (1/err)^(1/5), 0.2, 5).
    /// </summary>
    /// <param name="f">The system dy/dt = f(t, y).</param>
    /// <param name="y0">Initial state; copied, not modified.</param>
    /// <param name="t0">Initial time.</param>
    /// <param name="t1">Final time; must be greater than <paramref name="t0"/>.</param>
    /// <param name="relativeTolerance">Relative error tolerance per component, as a fraction of the component's magnitude.</param>
    /// <param name="absoluteTolerance">Absolute error tolerance per component, in the units of that state variable; it governs components near zero.</param>
    /// <param name="initialStep">First trial step size in time units; the controller adjusts it immediately.</param>
    /// <param name="maxStep">Upper bound on the step size, useful when a forcing term must not be stepped over.</param>
    /// <param name="maxSteps">Safety limit on the number of attempted steps.</param>
    /// <exception cref="InvalidOperationException">Thrown when <paramref name="maxSteps"/> attempts are exhausted before reaching <paramref name="t1"/>.</exception>
    public static OdeSolution Rk45(
        OdeSystem f,
        double[] y0,
        double t0,
        double t1,
        double relativeTolerance = 1e-6,
        double absoluteTolerance = 1e-9,
        double initialStep = 1e-3,
        double maxStep = double.PositiveInfinity,
        int maxSteps = 1_000_000)
    {
        ArgumentNullException.ThrowIfNull(f);
        ArgumentNullException.ThrowIfNull(y0);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(relativeTolerance);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(absoluteTolerance);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(initialStep);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxSteps);
        if (t1 <= t0) throw new ArgumentOutOfRangeException(nameof(t1), t1, "The final time must be greater than the initial time.");

        int n = y0.Length;
        var y = (double[])y0.Clone();
        var yNext = new double[n];
        var temp = new double[n];
        var k1 = new double[n];
        var k2 = new double[n];
        var k3 = new double[n];
        var k4 = new double[n];
        var k5 = new double[n];
        var k6 = new double[n];
        var k7 = new double[n];

        var times = new List<double>(1024) { t0 };
        var states = new List<double[]>(1024) { (double[])y.Clone() };

        double t = t0;
        double h = Math.Min(initialStep, Math.Min(maxStep, t1 - t0));
        int rejected = 0;
        int attempts = 0;
        bool firstStage = false;

        while (t < t1 - 1e-15)
        {
            if (++attempts > maxSteps)
            {
                throw new InvalidOperationException(
                    $"Rk45 exceeded maxSteps = {maxSteps} at t = {t:R} (step size {h:R}) before reaching t1 = {t1:R}. " +
                    "The system is probably stiff or the tolerances are too tight; relax the tolerances, shorten the horizon, or raise maxSteps.");
            }

            if (h > t1 - t) h = t1 - t;

            if (!firstStage) f(t, y, k1); // FSAL: k1 of this step is k7 of the accepted previous step.

            for (int i = 0; i < n; i++) temp[i] = y[i] + h * A21 * k1[i];
            f(t + C[1] * h, temp, k2);
            for (int i = 0; i < n; i++) temp[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
            f(t + C[2] * h, temp, k3);
            for (int i = 0; i < n; i++) temp[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
            f(t + C[3] * h, temp, k4);
            for (int i = 0; i < n; i++) temp[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
            f(t + C[4] * h, temp, k5);
            for (int i = 0; i < n; i++) temp[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
            f(t + C[5] * h, temp, k6);
            for (int i = 0; i < n; i++)
                yNext[i] = y[i] + h * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i]);
            f(t + h, yNext, k7);

            double error = 0.0;
            for (int i = 0; i < n; i++)
            {
                double fourth = y[i] + h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
                double scale = absoluteTolerance + relativeTolerance * Math.Max(Math.Abs(y[i]), Math.Abs(yNext[i]));
                double e = (yNext[i] - fourth) / scale;
                error += e * e;
            }
            error = n == 0 ? 0.0 : Math.Sqrt(error / n);

            double factor = error <= 0.0 ? 5.0 : Math.Clamp(0.9 * Math.Pow(1.0 / error, 0.2), 0.2, 5.0);

            if (error > 1.0)
            {
                rejected++;
                h *= factor;
                firstStage = false;
                continue;
            }

            t += h;
            Array.Copy(yNext, y, n);
            Array.Copy(k7, k1, n);
            firstStage = true;

            times.Add(t);
            states.Add((double[])y.Clone());

            h = Math.Min(h * factor, maxStep);
            if (h <= 0.0 || double.IsNaN(h))
                throw new InvalidOperationException($"Rk45 step size collapsed to {h:R} at t = {t:R}; the system may be stiff or discontinuous.");
        }

        return new OdeSolution([.. times], [.. states]) { RejectedSteps = rejected };
    }
}
