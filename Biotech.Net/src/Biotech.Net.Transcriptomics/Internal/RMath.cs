using Biotech.Net.Native;

namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>
/// Ports of the R (nmath / stats) routines that DESeq2 and apeglm depend on, so that the DESeq2 pipeline reproduces
/// R's arithmetic step for step: <c>dnbinom(log = TRUE)</c> (Loader's saddle-point form), <c>mean(x, trim)</c>,
/// <c>qf</c>, <c>pchisq</c>, <c>uniroot</c> (Brent's <c>zeroin</c>), <c>optim(method = "BFGS")</c> (<c>vmmin</c>) and
/// <c>optimHess</c>.
/// </summary>
internal static class RMath
{
    private const double LnSqrt2Pi = 0.918938533204672741780329736406; // log(sqrt(2*pi))
    private const double Ln2Pi = 1.837877066409345483560659472811;     // log(2*pi)

    private static readonly double[] Gamcs =
    [
        +.8571195590989331421920062399942e-2, +.4415381324841006757191315771652e-2, +.5685043681599363378632664588789e-1,
        -.4219835396418560501012500186624e-2, +.1326808181212460220584006796352e-2, -.1893024529798880432523947023886e-3,
        +.3606925327441245256578082217225e-4, -.6056761904460864218485548290365e-5, +.1055829546302283344731823509093e-5,
        -.1811967365542384048291855891166e-6, +.3117724964715322277790254593169e-7, -.5354219639019687140874081024347e-8,
        +.9193275519859588946887786825940e-9, -.1577941280288339761767423273953e-9, +.2707980622934954543266540433089e-10,
        -.4646818653825730144081661058933e-11, +.7973350192007419656460767175359e-12, -.1368078209830916025799499172309e-12,
        +.2347319486563800657233471771688e-13, -.4027432614949066932766570534699e-14, +.6910051747372100912138336975257e-15,
        -.1185584500221992907052387126192e-15,
    ];

    private static readonly double[] Algmcs =
    [
        +.1666389480451863247205729650822e+0, -.1384948176067563840732986059135e-4, +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10, +.6221098041892605227126015543416e-13,
    ];

    private static double ChebyshevEval(double x, double[] a, int n)
    {
        double twox = x * 2, b2 = 0, b1 = 0, b0 = 0;
        for (int i = 1; i <= n; i++) { b2 = b1; b1 = b0; b0 = twox * b1 - b2 + a[n - i]; }
        return (b0 - b2) * 0.5;
    }

    private static double Lgammacor(double x)
    {
        const double xbig = 94906265.62425156;
        if (x < 10) return double.NaN;
        if (x < xbig) { double tmp = 10 / x; return ChebyshevEval(tmp * tmp * 2 - 1, Algmcs, 5) / x; }
        return 1 / (x * 12);
    }

    /// <summary>R's <c>gammafn</c> for 0 &lt; x ≤ 10 (Chebyshev series on [1, 2) and recursion).</summary>
    private static double GammaSmall(double x)
    {
        int n = (int)x;
        if (x < 0) --n;
        double y = x - n;
        --n;
        double value = ChebyshevEval(y * 2 - 1, Gamcs, 22) + .9375;
        if (n == 0) return value;
        if (n < 0)
        {
            n = -n;
            for (int i = 0; i < n; i++) value /= x + i;
            return value;
        }
        for (int i = 1; i <= n; i++) value *= y + i;
        return value;
    }

    /// <summary>R's <c>lgammafn</c> for x &gt; 0 (what DESeq2's C++ <c>lgamma</c> / <c>Rf_lgammafn</c> calls evaluate).</summary>
    public static double LGamma(double x)
    {
        if (double.IsNaN(x)) return x;
        if (x <= 0) return x == Math.Truncate(x) ? double.PositiveInfinity : SpecialFunctions.LogGamma(x);
        if (x < 1e-306) return -Math.Log(x);
        if (x <= 10) return Math.Log(Math.Abs(GammaSmall(x)));
        if (x > 2.5327372760800758e+305) return double.PositiveInfinity;
        if (x > 1e17) return x * (Math.Log(x) - 1.0);
        if (x > 4934720.0) return LnSqrt2Pi + (x - 0.5) * Math.Log(x) - x;
        return LnSqrt2Pi + (x - 0.5) * Math.Log(x) - x + Lgammacor(x);
    }

    private static readonly double[] Bvalues =
    [
        1.00000000000000000e+00, -5.00000000000000000e-01, 1.66666666666666667e-01, -3.33333333333333333e-02,
        2.38095238095238095e-02, -3.33333333333333333e-02, 7.57575757575757576e-02, -2.53113553113553114e-01,
        1.16666666666666667e+00, -7.09215686274509804e+00, 5.49711779448621554e+01, -5.29124242424242424e+02,
        6.19212318840579710e+03, -8.65802531135531136e+04, 1.42551716666666667e+06, -2.72982310678160920e+07,
        6.01580873900642368e+08, -1.51163157670921569e+10, 4.29614643061166667e+11, -1.37116552050883328e+13,
        4.88332318973593167e+14, -1.92965793419400681e+16,
    ];

    /// <summary>R's <c>digamma(x)</c> for x &gt; 0 (Amos' <c>dpsifn</c> with n = 0, kode = 1, m = 1).</summary>
    public static double Digamma(double x)
    {
        if (double.IsNaN(x)) return x;
        if (x <= 0) return NumericUtil.Digamma(x);
        const double eps = 2.220446049250313e-16;
        double xln = Math.Log(x);
        if (x * xln > 1 / (2.0 * eps)) return xln;
        const double wdtol = 1.1102230246251565e-16;
        double elim = 2.302 * (1021 * 0.30102999566398120 - 3.0);
        double rln = Math.Min(0.30102999566398120 * 53, 18.06);
        double fln = Math.Max(rln, 3.0) - 3.0;
        double yint = 3.50 + 0.40 * fln;
        double t = xln;
        if (Math.Abs(t) > elim) return double.NaN;
        if (x < wdtol) return -1 / x;
        double xmin = (int)yint + 1;
        double xdmy = x, xdmln = xln, xinc = 0.0;
        if (x < xmin)
        {
            int nx0 = (int)x;
            xinc = xmin - nx0;
            xdmy = x + xinc;
            xdmln = Math.Log(xdmy);
        }
        // L10, fn = 0
        double tt = 0.5 / xdmy, t1 = tt, tst = wdtol * tt;
        double rxsq = 1.0 / (xdmy * xdmy), ta = 0.5 * rxsq;
        t = ta;
        double s = t * Bvalues[2];
        if (Math.Abs(s) >= tst)
        {
            double tk = 2.0;
            for (int k = 4; k <= 22; k++)
            {
                t = t * ((tk + 1) / (tk + 1.0)) * (tk / (tk + 2.0)) * rxsq;
                double trm = t * Bvalues[k - 1];
                if (Math.Abs(trm) < tst) break;
                s += trm;
                tk += 2.0;
            }
        }
        s += t1;
        if (xinc != 0.0)
        {
            int nx = (int)xinc;
            for (int i = 1; i <= nx; i++) s += 1.0 / (x + (nx - i));
        }
        return -(s - xdmln);
    }

    /// <summary>R's <c>dnbinom(x, size, mu = mu, log = TRUE)</c> (nmath <c>dnbinom_mu</c>).</summary>
    public static double DnbinomMuLog(double x, double size, double mu)
    {
        if (double.IsNaN(x) || double.IsNaN(size) || double.IsNaN(mu)) return x + size + mu;
        if (mu < 0 || size < 0) return double.NaN;
        if (x < 0 || !double.IsFinite(x)) return double.NegativeInfinity;
        if (x == 0 && size == 0) return 0;
        x = Math.Round(x);
        if (!double.IsFinite(size)) return DpoisRawLog(x, mu);
        if (x == 0)
            return size * (size < mu ? Math.Log(size / (size + mu)) : Log1p(-mu / (size + mu)));
        if (x < 1e-10 * size)
        {
            double p = size < mu ? Math.Log(size / (1 + size / mu)) : Math.Log(mu / (1 + mu / size));
            double xx2s = x * (x - 1) / 2 / size;
            return x * p - mu - LGamma(x + 1) + Log1p(xx2s);
        }
        double lp = x < size ? Log1p(-x / (size + x)) : Math.Log(size / (size + x));
        return lp + DbinomRawLog(size, x + size, size / (size + mu), mu / (size + mu));
    }

    private static double DpoisRawLog(double x, double lambda)
    {
        if (lambda == 0) return x == 0 ? 0 : double.NegativeInfinity;
        if (x == 0) return -lambda;
        return -Stirlerr(x) - Bd0(x, lambda) - 0.5 * Math.Log(2 * Math.PI * x);
    }

    private static double DbinomRawLog(double x, double n, double p, double q)
    {
        if (p == 0) return x == 0 ? 0 : double.NegativeInfinity;
        if (q == 0) return x == n ? 0 : double.NegativeInfinity;
        if (x == 0)
        {
            if (n == 0) return 0;
            return p > q ? n * Math.Log(q) : n * Log1p(-p);
        }
        if (x == n) return p > q ? n * Log1p(-q) : n * Math.Log(p);
        if (x < 0 || x > n) return double.NegativeInfinity;
        double lc = Stirlerr(n) - Stirlerr(x) - Stirlerr(n - x) - Bd0(x, n * p) - Bd0(n - x, n * q);
        double lf = Ln2Pi + Math.Log(x) + Log1p(-x / n);
        return lc - 0.5 * lf;
    }

    private static readonly double[] SferrHalves =
    [
        0.0, 0.1534264097200273452913848, 0.0810614667953272582196702, 0.0548141210519176538961390,
        0.0413406959554092940938221, 0.03316287351993628748511048, 0.02767792568499833914878929,
        0.02374616365629749597132920, 0.02079067210376509311152277, 0.01848845053267318523077934,
        0.01664469118982119216319487, 0.01513497322191737887351255, 0.01387612882307074799874573,
        0.01281046524292022692424986, 0.01189670994589177009505572, 0.01110455975820691732662991,
        0.010411265261972096497478567, 0.009799416126158803298389475, 0.009255462182712732917728637,
        0.008768700134139385462952823, 0.008330563433362871256469318, 0.007934114564314020547248100,
        0.007573675487951840794972024, 0.007244554301320383179543912, 0.006942840107209529865664152,
        0.006665247032707682442354394, 0.006408994188004207068439631, 0.006171712263039457647532867,
        0.005951370112758847735624416, 0.005746216513010115682023589, 0.005554733551962801371038690,
    ];

    /// <summary>nmath <c>stirlerr(n) = log(n!) - log(sqrt(2 pi n) (n/e)^n)</c> (R 4.4+ version).</summary>
    private static double Stirlerr(double n)
    {
        const double S0 = 0.083333333333333333333, S1 = 0.00277777777777777777778, S2 = 0.00079365079365079365079365,
            S3 = 0.000595238095238095238095238, S4 = 0.0008417508417508417508417508, S5 = 0.0019175269175269175269175262,
            S6 = 0.0064102564102564102564102561, S7 = 0.029550653594771241830065352, S8 = 0.17964437236883057316493850,
            S9 = 1.3924322169059011164274315, S10 = 13.402864044168391994478957, S11 = 156.84828462600201730636509,
            S12 = 2193.1033333333333333333333, S13 = 36108.771253724989357173269, S14 = 691472.26885131306710839498,
            S15 = 15238221.539407416192283370, S16 = 382900751.39141414141414141;
        double nn;
        if (n <= 23.5)
        {
            nn = n + n;
            if (n <= 15.0 && nn == (int)nn) return SferrHalves[(int)nn];
            if (n <= 5.25)
            {
                if (n >= 1.0)
                {
                    double ln = Math.Log(n);
                    return LGamma(n) + n * (1 - ln) + (ln - Ln2Pi) / 2;
                }
                return LGamma(n + 1) - (n + 0.5) * Math.Log(n) + n - LnSqrt2Pi;
            }
            nn = n * n;
            if (n > 12.8) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - S6 / nn) / nn) / nn) / nn) / nn) / nn) / n;
            if (n > 12.3) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - S7 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
            if (n > 8.9) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - S8 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
            if (n > 7.3) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - S10 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
            if (n > 6.6) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - (S10 - (S11 - S12 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
            if (n > 6.1) return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - (S10 - (S11 - (S12 - (S13 - S14 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
            return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - (S8 - (S9 - (S10 - (S11 - (S12 - (S13 - (S14 - (S15 - S16 / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / nn) / n;
        }
        nn = n * n;
        if (n > 15.7e6) return S0 / n;
        if (n > 6180) return (S0 - S1 / nn) / n;
        if (n > 205) return (S0 - (S1 - S2 / nn) / nn) / n;
        if (n > 86) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
        if (n > 27) return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
        return (S0 - (S1 - (S2 - (S3 - (S4 - S5 / nn) / nn) / nn) / nn) / nn) / n;
    }

    /// <summary>nmath <c>bd0(x, np) = x log(x/np) + np - x</c>, evaluated without cancellation.</summary>
    private static double Bd0(double x, double np)
    {
        if (!double.IsFinite(x) || !double.IsFinite(np) || np == 0.0) return double.NaN;
        if (Math.Abs(x - np) < 0.1 * (x + np))
        {
            double d = x - np, v = d / (x + np);
            if (d != 0.0 && v == 0.0) { double x_ = x / 4, n_ = np / 4; v = (x_ - n_) / (x_ + n_); }
            double s = d / 2 * v;
            if (Math.Abs(2 * s) < 2.2250738585072014e-308) return 2 * s; // |2s| < DBL_MIN
            double ej = x * v;
            v *= v;
            for (int j = 1; j < 1000; j++)
            {
                ej *= v;
                double s0 = s;
                s += ej / ((j << 1) + 1);
                if (s == s0) return 2 * s;
            }
        }
        double lgxn = double.IsFinite(x / np) ? Math.Log(x / np) : Math.Log(x) - Math.Log(np);
        return x > np ? x * (lgxn - 1.0) + np : x * lgxn + np - x;
    }

    /// <summary>log(1 + x) accurate near 0.</summary>
    public static double Log1p(double x)
    {
        double u = 1 + x;
        return u == 1 ? x : Math.Log(u) * x / (u - 1);
    }

    /// <summary>R's <c>mean(x)</c>: long-double style two-pass mean (sum, then a correction pass).</summary>
    public static double Mean(ReadOnlySpan<double> x)
    {
        int n = x.Length;
        if (n == 0) return double.NaN;
        double s = 0;
        foreach (double v in x) s += v;
        s /= n;
        if (!double.IsFinite(s)) return s;
        double t = 0;
        foreach (double v in x) t += v - s;
        return s + t / n;
    }

    /// <summary>R's <c>mean(x, trim = trim)</c>.</summary>
    public static double TrimmedMean(ReadOnlySpan<double> x, double trim)
    {
        int n = x.Length;
        if (n == 0) return double.NaN;
        if (trim <= 0) return Mean(x);
        foreach (double v in x) if (double.IsNaN(v)) return double.NaN;
        if (trim >= 0.5) return NumericUtil.Median(x);
        int lo = (int)Math.Floor(n * trim) + 1;
        int hi = n + 1 - lo;
        var s = x.ToArray();
        Array.Sort(s);
        return Mean(s.AsSpan(lo - 1, hi - lo + 1));
    }

    /// <summary>Upper normal tail <c>pnorm(-|z|)</c> doubled: the two-sided Wald p-value.</summary>
    public static double TwoSidedNormal(double z) => double.IsNaN(z) ? double.NaN : 2 * NumericUtil.NormalUpper(Math.Abs(z));

    /// <summary><c>pchisq(x, df, lower.tail = FALSE)</c>.</summary>
    public static double ChiSquareUpper(double x, double df)
    {
        if (double.IsNaN(x)) return double.NaN;
        return x <= 0 ? 1 : SpecialFunctions.RegularizedGammaQ(df / 2, x / 2);
    }

    /// <summary><c>pnorm(q, mean, sd)</c> (lower tail).</summary>
    public static double Pnorm(double q, double mean, double sd)
    {
        if (double.IsNaN(q) || double.IsNaN(sd)) return double.NaN;
        double z = (q - mean) / sd;
        return z < 0 ? NumericUtil.NormalUpper(-z) : 1 - NumericUtil.NormalUpper(z);
    }

    /// <summary><c>qnorm(p)</c> refined by Newton steps on the accurate normal tail.</summary>
    public static double Qnorm(double p)
    {
        double x = SpecialFunctions.NormalQuantile(p);
        for (int i = 0; i < 3; i++)
        {
            double cdf = Pnorm(x, 0, 1);
            double pdf = Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * Math.PI);
            if (pdf == 0) break;
            x -= (cdf - p) / pdf;
        }
        return x;
    }

    /// <summary><c>qf(p, df1, df2)</c> via inversion of the regularised incomplete beta function.</summary>
    public static double Qf(double p, double df1, double df2)
    {
        double a = df1 / 2, b = df2 / 2;
        // Solve I_x(a, b) = p by bisection refined with Newton steps.
        double lo = 0, hi = 1, x = 0.5;
        for (int it = 0; it < 200; it++)
        {
            x = 0.5 * (lo + hi);
            double f = SpecialFunctions.RegularizedIncompleteBeta(a, b, x) - p;
            if (f > 0) hi = x; else lo = x;
            if (hi - lo < 1e-15) break;
        }
        double lbeta = SpecialFunctions.LogGamma(a) + SpecialFunctions.LogGamma(b) - SpecialFunctions.LogGamma(a + b);
        for (int it = 0; it < 5; it++)
        {
            double f = SpecialFunctions.RegularizedIncompleteBeta(a, b, x) - p;
            double dens = Math.Exp((a - 1) * Math.Log(x) + (b - 1) * Math.Log(1 - x) - lbeta);
            if (dens <= 0 || !double.IsFinite(dens)) break;
            double nx = x - f / dens;
            if (nx <= lo || nx >= hi) break;
            x = nx;
        }
        return df2 / df1 * x / (1 - x);
    }

    /// <summary>R's <c>uniroot</c> (Brent's <c>R_zeroin2</c>) with its default tolerance <c>.Machine$double.eps^0.25</c>.</summary>
    public static double Uniroot(Func<double, double> f, double ax, double bx, double tol = 1.220703125e-4, int maxIter = 1000)
    {
        double fa = f(ax), fb = f(bx);
        if (fa * fb > 0) throw new ArgumentException("f() values at end points not of opposite sign");
        double a = ax, b = bx, c = a, fc = fa;
        int maxit = maxIter + 1;
        const double eps = 2.220446049250313e-16;
        if (fa == 0.0) return a;
        if (fb == 0.0) return b;
        while (maxit-- > 0)
        {
            double prevStep = b - a;
            if (Math.Abs(fc) < Math.Abs(fb))
            {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }
            double tolAct = 2 * eps * Math.Abs(b) + tol / 2;
            double newStep = (c - b) / 2;
            if (Math.Abs(newStep) <= tolAct || fb == 0) return b;
            if (Math.Abs(prevStep) >= tolAct && Math.Abs(fa) > Math.Abs(fb))
            {
                double t1, cb = c - b, t2, p, q;
                if (a == c) { t1 = fb / fa; p = cb * t1; q = 1.0 - t1; }
                else
                {
                    q = fa / fc; t1 = fb / fc; t2 = fb / fa;
                    p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                    q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
                }
                if (p > 0) q = -q; else p = -p;
                if (p < (0.75 * cb * q - Math.Abs(tolAct * q) / 2) && p < Math.Abs(prevStep * q / 2)) newStep = p / q;
            }
            if (Math.Abs(newStep) < tolAct) newStep = newStep > 0 ? tolAct : -tolAct;
            a = b; fa = fb;
            b += newStep; fb = f(b);
            if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; }
        }
        return b;
    }

    /// <summary>
    /// R's <c>optim(par, fn, gr, method = "BFGS")</c> with default control (maxit 100, reltol sqrt(eps)): the
    /// <c>vmmin</c> variable-metric method of Nash. Returns the convergence code (0 = converged, 1 = maxit reached).
    /// </summary>
    public static int Vmmin(double[] b, Func<double[], double> fn, Action<double[], double[]> gr, out double fmin,
        int maxit = 100, double reltol = 1.4901161193847656e-08)
    {
        const double stepredn = 0.2, acctol = 0.0001, reltest = 10.0;
        double abstol = double.NegativeInfinity;
        int n = b.Length;
        var g = new double[n]; var t = new double[n]; var X = new double[n]; var c = new double[n];
        var B = new double[n, n];
        double f = fn(b);
        if (!double.IsFinite(f)) throw new InvalidOperationException("initial value in 'vmmin' is not finite");
        fmin = f;
        int gradcount = 1, iter = 0, count;
        gr(b, g);
        iter++;
        int ilast = gradcount;
        do
        {
            if (ilast == gradcount)
                for (int i = 0; i < n; i++) { for (int j = 0; j < i; j++) B[i, j] = 0.0; B[i, i] = 1.0; }
            for (int i = 0; i < n; i++) { X[i] = b[i]; c[i] = g[i]; }
            double gradproj = 0.0;
            for (int i = 0; i < n; i++)
            {
                double s = 0.0;
                for (int j = 0; j <= i; j++) s -= B[i, j] * g[j];
                for (int j = i + 1; j < n; j++) s -= B[j, i] * g[j];
                t[i] = s;
                gradproj += s * g[i];
            }
            if (gradproj < 0.0)
            {
                double steplength = 1.0;
                bool accpoint = false;
                do
                {
                    count = 0;
                    for (int i = 0; i < n; i++)
                    {
                        b[i] = X[i] + steplength * t[i];
                        if (reltest + X[i] == reltest + b[i]) count++;
                    }
                    if (count < n)
                    {
                        f = fn(b);
                        accpoint = double.IsFinite(f) && f <= fmin + gradproj * steplength * acctol;
                        if (!accpoint) steplength *= stepredn;
                    }
                } while (!(count == n || accpoint));
                bool enough = f > abstol && Math.Abs(f - fmin) > reltol * (Math.Abs(fmin) + reltol);
                if (!enough) { count = n; fmin = f; }
                if (count < n)
                {
                    fmin = f;
                    gr(b, g);
                    gradcount++;
                    iter++;
                    double D1 = 0.0;
                    for (int i = 0; i < n; i++) { t[i] = steplength * t[i]; c[i] = g[i] - c[i]; D1 += t[i] * c[i]; }
                    if (D1 > 0)
                    {
                        double D2 = 0.0;
                        for (int i = 0; i < n; i++)
                        {
                            double s = 0.0;
                            for (int j = 0; j <= i; j++) s += B[i, j] * c[j];
                            for (int j = i + 1; j < n; j++) s += B[j, i] * c[j];
                            X[i] = s;
                            D2 += s * c[i];
                        }
                        D2 = 1.0 + D2 / D1;
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j <= i; j++)
                                B[i, j] += (D2 * t[i] * t[j] - X[i] * t[j] - t[i] * X[j]) / D1;
                    }
                    else ilast = gradcount;
                }
                else if (ilast < gradcount) { count = 0; ilast = gradcount; }
            }
            else
            {
                count = 0;
                if (ilast == gradcount) count = n;
                else ilast = gradcount;
            }
            if (iter >= maxit) break;
            if (gradcount - ilast > 2 * n) ilast = gradcount;
        } while (count != n || ilast != gradcount);
        return iter < maxit ? 0 : 1;
    }

    /// <summary>R's <c>optimHess(par, fn, gr)</c>: central differences of the gradient (ndeps = 1e-3), symmetrised.</summary>
    public static double[,] OptimHess(double[] par, Action<double[], double[]> gr, double ndeps = 1e-3)
    {
        int n = par.Length;
        var h = new double[n, n];
        var dpar = (double[])par.Clone();
        var df1 = new double[n]; var df2 = new double[n];
        for (int i = 0; i < n; i++)
        {
            double eps = ndeps;
            dpar[i] += eps;
            gr(dpar, df1);
            dpar[i] -= 2 * eps;
            gr(dpar, df2);
            for (int j = 0; j < n; j++) h[i, j] = (df1[j] - df2[j]) / (2 * eps);
            dpar[i] += eps;
        }
        var res = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) res[i, j] = 0.5 * (h[i, j] + h[j, i]);
        return res;
    }
}
