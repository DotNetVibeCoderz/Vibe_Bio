namespace Biotech.Net.SyntheticBiology.Circuits;

/// <summary>Which integrator a simulation helper should use.</summary>
public enum IntegratorKind
{
    /// <summary>Classical fixed-step fourth-order Runge-Kutta; uniform output grid, no error control.</summary>
    Rk4,

    /// <summary>Adaptive Dormand-Prince 5(4); non-uniform output grid with local error control.</summary>
    Rk45,
}

/// <summary>
/// Parameters of a single Hill-type regulated gene, the standard phenomenological description of
/// transcription factor control (Hill, J. Physiol. 40:iv-vii, 1910; Alon,
/// <i>An Introduction to Systems Biology</i>, 2nd ed., 2019, chapter 2).
/// </summary>
/// <param name="MaxRate">
/// Maximum transcription rate at saturating input, beta_max, in concentration per unit time
/// (for example nM/min or, in dimensionless models, units of the degradation rate).
/// </param>
/// <param name="LeakRate">
/// Basal ("leaky") transcription rate with no input present, in the same units as
/// <paramref name="MaxRate"/>. Real promoters are never fully off; a leak of 0.1-1 percent of the
/// maximum is typical for well-repressed synthetic promoters.
/// </param>
/// <param name="DissociationConstant">
/// The input concentration K giving half-maximal response, in concentration units. For a repressor
/// this is the half-repression threshold; for an inducer the half-induction threshold.
/// </param>
/// <param name="Cooperativity">
/// The Hill coefficient n (dimensionless), reflecting effective cooperative binding: n = 1 for simple
/// binding, n = 2 for a dimeric repressor such as TetR, n up to 4 for a tetramer such as LacI.
/// Larger n gives a sharper, more switch-like transition.
/// </param>
/// <param name="DegradationRate">
/// First-order removal rate constant of the species produced, in inverse time. In growing bacteria
/// this is dominated by dilution, ln(2)/doubling time, unless a degradation tag is present.
/// </param>
public sealed record HillParameters(
    double MaxRate,
    double LeakRate,
    double DissociationConstant,
    double Cooperativity,
    double DegradationRate)
{
    /// <summary>
    /// A conventional parameter set for a moderately strong, slightly leaky inducible promoter:
    /// maximum transcription 10 concentration units per minute, 1 percent leak, half-induction at 1
    /// concentration unit, Hill coefficient 2 (dimeric regulator) and a removal rate of 0.1 per
    /// minute (about a seven-minute half-life, typical of a tagged protein or fast dilution).
    /// </summary>
    public static HillParameters Default { get; } = new(10.0, 0.1, 1.0, 2.0, 0.1);

    /// <summary>
    /// The activating Hill response at an input concentration, LeakRate + MaxRate * x^n / (K^n + x^n),
    /// in units of transcription rate.
    /// </summary>
    /// <param name="input">Inducer or activator concentration; negative values are treated as zero.</param>
    public double Activation(double input)
    {
        double x = input <= 0.0 ? 0.0 : Math.Pow(input, Cooperativity);
        double k = Math.Pow(DissociationConstant, Cooperativity);
        return LeakRate + MaxRate * x / (k + x);
    }

    /// <summary>
    /// The repressing Hill response at a repressor concentration,
    /// LeakRate + MaxRate * K^n / (K^n + x^n), in units of transcription rate.
    /// </summary>
    /// <param name="repressor">Repressor concentration; negative values are treated as zero.</param>
    public double Repression(double repressor)
    {
        double x = repressor <= 0.0 ? 0.0 : Math.Pow(repressor, Cooperativity);
        double k = Math.Pow(DissociationConstant, Cooperativity);
        return LeakRate + MaxRate * k / (k + x);
    }
}

/// <summary>
/// Deterministic ordinary-differential-equation models of the canonical synthetic gene circuits:
/// the genetic toggle switch, the repressilator and a single inducible gene. Each model is exposed
/// both as a bare <see cref="OdeSystem"/>, so it can be composed or integrated by other code, and as
/// a ready-made simulation helper.
/// </summary>
public static class CircuitModels
{
    /// <summary>
    /// The dimensionless genetic toggle switch of Gardner, Cantor and Collins
    /// (<i>Construction of a genetic toggle switch in Escherichia coli</i>, Nature 403:339-342, 2000):
    /// two repressors each transcribed from a promoter repressed by the other,
    /// du/dt = alpha1 / (1 + v^beta) - u and dv/dt = alpha2 / (1 + u^gamma) - v.
    /// Time is scaled by the protein removal rate and concentrations by the repression thresholds, so
    /// all four parameters are dimensionless. The circuit is bistable when both alphas are large
    /// enough and at least one exponent exceeds 1; without cooperativity (beta = gamma = 1) it has a
    /// single stable state.
    /// </summary>
    /// <param name="alpha1">Dimensionless maximum synthesis rate of repressor u.</param>
    /// <param name="alpha2">Dimensionless maximum synthesis rate of repressor v.</param>
    /// <param name="beta">Cooperativity of repression of the u promoter by v.</param>
    /// <param name="gamma">Cooperativity of repression of the v promoter by u.</param>
    /// <returns>A two-state system with state vector [u, v].</returns>
    public static OdeSystem ToggleSwitch(double alpha1, double alpha2, double beta, double gamma)
    {
        return (double t, ReadOnlySpan<double> y, Span<double> dydt) =>
        {
            double u = y[0] <= 0.0 ? 0.0 : y[0];
            double v = y[1] <= 0.0 ? 0.0 : y[1];
            dydt[0] = alpha1 / (1.0 + Math.Pow(v, beta)) - y[0];
            dydt[1] = alpha2 / (1.0 + Math.Pow(u, gamma)) - y[1];
        };
    }

    /// <summary>
    /// Integrates the dimensionless toggle switch from a given initial condition.
    /// </summary>
    /// <param name="alpha1">Dimensionless maximum synthesis rate of repressor u.</param>
    /// <param name="alpha2">Dimensionless maximum synthesis rate of repressor v.</param>
    /// <param name="beta">Cooperativity of repression of the u promoter by v.</param>
    /// <param name="gamma">Cooperativity of repression of the v promoter by u.</param>
    /// <param name="initial">Initial state [u0, v0] in dimensionless concentration units.</param>
    /// <param name="tEnd">End of the simulated interval in dimensionless time (one unit = one protein lifetime).</param>
    /// <param name="integrator">Which integrator to use.</param>
    /// <param name="step">Fixed step for <see cref="IntegratorKind.Rk4"/>, ignored by the adaptive method.</param>
    public static OdeSolution SimulateToggleSwitch(
        double alpha1,
        double alpha2,
        double beta,
        double gamma,
        double[] initial,
        double tEnd,
        IntegratorKind integrator = IntegratorKind.Rk45,
        double step = 0.01)
    {
        ArgumentNullException.ThrowIfNull(initial);
        if (initial.Length != 2) throw new ArgumentException("The toggle switch has exactly two states, [u, v].", nameof(initial));
        var system = ToggleSwitch(alpha1, alpha2, beta, gamma);
        return integrator == IntegratorKind.Rk4
            ? OdeIntegrator.Rk4(system, initial, 0.0, tEnd, step)
            : OdeIntegrator.Rk45(system, initial, 0.0, tEnd, 1e-8, 1e-10, Math.Min(step, 1e-3));
    }

    /// <summary>
    /// Published default maximum transcription rate alpha of the repressilator, in units of the mRNA
    /// decay rate (Elowitz &amp; Leibler, Nature 403:335-338, 2000).
    /// </summary>
    public const double RepressilatorAlpha = 216.0;

    /// <summary>Published default leak alpha0 = alpha / 1000, i.e. 0.1 percent basal transcription.</summary>
    public const double RepressilatorLeak = 0.216;

    /// <summary>
    /// Published default ratio beta of the protein decay rate to the mRNA decay rate; beta = 5 means
    /// proteins are five times more stable than their messages.
    /// </summary>
    public const double RepressilatorBeta = 5.0;

    /// <summary>Published default Hill coefficient n = 2 (dimeric repressors).</summary>
    public const double RepressilatorHill = 2.0;

    /// <summary>
    /// The repressilator of Elowitz and Leibler (<i>A synthetic oscillatory network of
    /// transcriptional regulators</i>, Nature 403:335-338, 2000): three repressor genes in a cyclic
    /// negative feedback loop (lacI represses tetR represses cI represses lacI). The dimensionless
    /// six-state model is
    /// dm_i/dt = -m_i + alpha / (1 + p_j^n) + alpha0 and dp_i/dt = -beta (p_i - m_i),
    /// with j = i - 1 (mod 3), time scaled by the mRNA lifetime and protein concentration scaled by
    /// the repression threshold K. With the published defaults the symmetric fixed point is unstable
    /// and the system settles onto a limit cycle.
    /// </summary>
    /// <param name="alpha">Maximum transcription rate in units of the mRNA decay rate.</param>
    /// <param name="alpha0">Leak transcription rate in the same units.</param>
    /// <param name="beta">Ratio of protein to mRNA decay rate.</param>
    /// <param name="n">Hill coefficient of repression.</param>
    /// <returns>A six-state system with state vector [m1, m2, m3, p1, p2, p3].</returns>
    public static OdeSystem Repressilator(
        double alpha = RepressilatorAlpha,
        double alpha0 = RepressilatorLeak,
        double beta = RepressilatorBeta,
        double n = RepressilatorHill)
    {
        return (double t, ReadOnlySpan<double> y, Span<double> dydt) =>
        {
            for (int i = 0; i < 3; i++)
            {
                int j = (i + 2) % 3;               // the repressor acting on gene i is protein i-1 (mod 3)
                double p = y[3 + j] <= 0.0 ? 0.0 : y[3 + j];
                dydt[i] = -y[i] + alpha / (1.0 + Math.Pow(p, n)) + alpha0;
                dydt[3 + i] = -beta * (y[3 + i] - y[i]);
            }
        };
    }

    /// <summary>
    /// Integrates the repressilator from a given initial condition. The symmetric state is a fixed
    /// point, so a non-symmetric initial condition is required to see oscillations.
    /// </summary>
    /// <param name="initial">Initial state [m1, m2, m3, p1, p2, p3]; defaults to a mildly asymmetric state.</param>
    /// <param name="tEnd">End of the simulated interval in units of the mRNA lifetime.</param>
    /// <param name="alpha">Maximum transcription rate in units of the mRNA decay rate.</param>
    /// <param name="alpha0">Leak transcription rate in the same units.</param>
    /// <param name="beta">Ratio of protein to mRNA decay rate.</param>
    /// <param name="n">Hill coefficient of repression.</param>
    /// <param name="integrator">Which integrator to use.</param>
    /// <param name="step">Fixed step for <see cref="IntegratorKind.Rk4"/>, ignored by the adaptive method.</param>
    /// <param name="outputEvery">For <see cref="IntegratorKind.Rk4"/>, record every n-th step.</param>
    public static OdeSolution SimulateRepressilator(
        double[]? initial = null,
        double tEnd = 300.0,
        double alpha = RepressilatorAlpha,
        double alpha0 = RepressilatorLeak,
        double beta = RepressilatorBeta,
        double n = RepressilatorHill,
        IntegratorKind integrator = IntegratorKind.Rk4,
        double step = 0.005,
        int outputEvery = 10)
    {
        double[] y0 = initial ?? [0.2, 0.1, 0.3, 1.0, 2.0, 3.0];
        if (y0.Length != 6) throw new ArgumentException("The repressilator has exactly six states, [m1, m2, m3, p1, p2, p3].", nameof(initial));
        var system = Repressilator(alpha, alpha0, beta, n);
        return integrator == IntegratorKind.Rk4
            ? OdeIntegrator.Rk4(system, y0, 0.0, tEnd, step, outputEvery)
            : OdeIntegrator.Rk45(system, y0, 0.0, tEnd, 1e-8, 1e-10, 1e-4);
    }

    /// <summary>
    /// Estimates the oscillation period of one state variable by timing successive upward crossings of
    /// its own mean: the series is mean-subtracted, crossings from negative to positive are located by
    /// linear interpolation between samples, and the period is the mean interval between them.
    /// Returns 0 when fewer than two upward crossings occur, which is the expected answer for a
    /// non-oscillating (monotone or damped) trajectory.
    /// <para>
    /// For an adaptive solution the samples are unevenly spaced; the crossing times are taken from the
    /// recorded time grid, so the estimate is accurate as long as the oscillation is well resolved.
    /// Transients bias the estimate, so discard the beginning of the trajectory if a limit-cycle
    /// period is wanted.
    /// </para>
    /// </summary>
    /// <param name="solution">A solution produced by one of the simulate helpers.</param>
    /// <param name="stateIndex">Index of the state variable to analyse.</param>
    /// <returns>The mean period in the time units of the solution, or 0 if no period could be measured.</returns>
    public static double EstimateOscillationPeriod(OdeSolution solution, int stateIndex)
    {
        ArgumentNullException.ThrowIfNull(solution);
        double[] series = solution.SeriesOf(stateIndex);
        if (series.Length < 3) return 0.0;

        double mean = 0.0;
        for (int i = 0; i < series.Length; i++) mean += series[i];
        mean /= series.Length;

        double first = 0.0, last = 0.0;
        int crossings = 0;
        for (int i = 1; i < series.Length; i++)
        {
            double a = series[i - 1] - mean;
            double b = series[i] - mean;
            if (a < 0.0 && b >= 0.0)
            {
                double w = b == a ? 0.0 : -a / (b - a);
                double tc = solution.Times[i - 1] + w * (solution.Times[i] - solution.Times[i - 1]);
                if (crossings == 0) first = tc;
                last = tc;
                crossings++;
            }
        }
        return crossings < 2 ? 0.0 : (last - first) / (crossings - 1);
    }

    /// <summary>
    /// A single gene under an inducible promoter, with explicit mRNA and protein states:
    /// dm/dt = (leak + beta_max * I^n / (K^n + I^n)) - delta_m * m and
    /// dp/dt = k_translate * m - delta_p * p,
    /// where I is the (constant) inducer concentration. Transcription responds to the inducer through
    /// an activating Hill function; translation is first order in the message. This is the standard
    /// two-stage description of expression from, for example, pLac with IPTG or pTet with
    /// anhydrotetracycline.
    /// </summary>
    /// <param name="inducer">Inducer concentration, held constant, in the same units as the dissociation constant.</param>
    /// <param name="hill">Transcription parameters; <see cref="HillParameters.DegradationRate"/> is the mRNA decay rate.</param>
    /// <param name="translationRate">Protein synthesis rate per unit mRNA per unit time.</param>
    /// <param name="proteinDegradation">First-order protein removal rate (degradation plus dilution), in inverse time.</param>
    /// <returns>A two-state system with state vector [mRNA, protein].</returns>
    public static OdeSystem InducibleExpression(
        double inducer,
        HillParameters hill,
        double translationRate,
        double proteinDegradation)
    {
        ArgumentNullException.ThrowIfNull(hill);
        double transcription = hill.Activation(inducer);
        double mrnaDecay = hill.DegradationRate;
        return (double t, ReadOnlySpan<double> y, Span<double> dydt) =>
        {
            dydt[0] = transcription - mrnaDecay * y[0];
            dydt[1] = translationRate * y[0] - proteinDegradation * y[1];
        };
    }

    /// <summary>
    /// Integrates the inducible expression model from a given initial condition at a fixed inducer
    /// concentration.
    /// </summary>
    /// <param name="inducer">Inducer concentration, held constant.</param>
    /// <param name="hill">Transcription parameters.</param>
    /// <param name="translationRate">Protein synthesis rate per unit mRNA per unit time.</param>
    /// <param name="proteinDegradation">First-order protein removal rate, in inverse time.</param>
    /// <param name="initial">Initial state [mRNA, protein]; defaults to an empty cell.</param>
    /// <param name="tEnd">End of the simulated interval, in the time units of the rate constants (minutes by convention).</param>
    /// <param name="integrator">Which integrator to use.</param>
    /// <param name="step">Fixed step for <see cref="IntegratorKind.Rk4"/>, ignored by the adaptive method.</param>
    public static OdeSolution SimulateInduction(
        double inducer,
        HillParameters hill,
        double translationRate = 1.0,
        double proteinDegradation = 0.02,
        double[]? initial = null,
        double tEnd = 1000.0,
        IntegratorKind integrator = IntegratorKind.Rk45,
        double step = 0.05)
    {
        ArgumentNullException.ThrowIfNull(hill);
        double[] y0 = initial ?? [0.0, 0.0];
        if (y0.Length != 2) throw new ArgumentException("The inducible expression model has exactly two states, [mRNA, protein].", nameof(initial));
        var system = InducibleExpression(inducer, hill, translationRate, proteinDegradation);
        return integrator == IntegratorKind.Rk4
            ? OdeIntegrator.Rk4(system, y0, 0.0, tEnd, step)
            : OdeIntegrator.Rk45(system, y0, 0.0, tEnd, 1e-9, 1e-12, Math.Min(step, 1e-3));
    }

    /// <summary>
    /// The analytic steady-state protein concentration of the inducible expression model,
    /// p* = k_translate * (leak + beta_max * I^n / (K^n + I^n)) / (delta_m * delta_p),
    /// obtained by setting both derivatives to zero. Useful as a validation target for the integrators
    /// and as a fast dose-response curve.
    /// </summary>
    /// <param name="inducer">Inducer concentration.</param>
    /// <param name="hill">Transcription parameters.</param>
    /// <param name="translationRate">Protein synthesis rate per unit mRNA per unit time.</param>
    /// <param name="proteinDegradation">First-order protein removal rate, in inverse time; must be positive.</param>
    public static double SteadyStateProtein(
        double inducer,
        HillParameters hill,
        double translationRate = 1.0,
        double proteinDegradation = 0.02)
    {
        ArgumentNullException.ThrowIfNull(hill);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(proteinDegradation);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(hill.DegradationRate);
        double mrna = hill.Activation(inducer) / hill.DegradationRate;
        return translationRate * mrna / proteinDegradation;
    }

    /// <summary>
    /// The analytic steady-state mRNA concentration of the inducible expression model,
    /// m* = (leak + beta_max * I^n / (K^n + I^n)) / delta_m.
    /// </summary>
    /// <param name="inducer">Inducer concentration.</param>
    /// <param name="hill">Transcription parameters.</param>
    public static double SteadyStateMrna(double inducer, HillParameters hill)
    {
        ArgumentNullException.ThrowIfNull(hill);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(hill.DegradationRate);
        return hill.Activation(inducer) / hill.DegradationRate;
    }
}
