namespace Biotech.Net.Native;

/// <summary>
/// Special functions and distribution tails used by the statistics kernels (t-tests, chi-square,
/// Fisher's exact test). Implementations follow Numerical Recipes / Cephes and are accurate to ~1e-12.
/// </summary>
public static class SpecialFunctions
{
    private static readonly double[] LanczosCoefficients =
    [
        0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    /// <summary>Natural log of the gamma function for x &gt; 0 (Lanczos, g = 7).</summary>
    public static double LogGamma(double x)
    {
        if (x <= 0) throw new ArgumentOutOfRangeException(nameof(x), "LogGamma is defined for x > 0.");
        if (x < 0.5)
            return Math.Log(Math.PI / Math.Sin(Math.PI * x)) - LogGamma(1 - x); // reflection
        x -= 1;
        double a = LanczosCoefficients[0];
        double t = x + 7.5;
        for (int i = 1; i < 9; i++) a += LanczosCoefficients[i] / (x + i);
        return 0.5 * Math.Log(2 * Math.PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(a);
    }

    /// <summary>log(n!) for n ≥ 0.</summary>
    public static double LogFactorial(int n) => n < 2 ? 0 : LogGamma(n + 1.0);

    /// <summary>Regularised incomplete beta I_x(a, b).</summary>
    public static double RegularizedIncompleteBeta(double a, double b, double x)
    {
        if (x <= 0) return 0;
        if (x >= 1) return 1;
        double lnFront = LogGamma(a + b) - LogGamma(a) - LogGamma(b) + a * Math.Log(x) + b * Math.Log(1 - x);
        // The continued fraction converges fast for x < (a+1)/(a+b+2); use the symmetry otherwise.
        if (x < (a + 1) / (a + b + 2))
            return Math.Exp(lnFront) * BetaContinuedFraction(a, b, x) / a;
        return 1 - Math.Exp(lnFront) * BetaContinuedFraction(b, a, 1 - x) / b;
    }

    private static double BetaContinuedFraction(double a, double b, double x)
    {
        const double eps = 1e-15, tiny = 1e-300;
        double qab = a + b, qap = a + 1, qam = a - 1;
        double c = 1, d = 1 - qab * x / qap;
        if (Math.Abs(d) < tiny) d = tiny;
        d = 1 / d;
        double h = d;
        for (int m = 1; m <= 300; m++)
        {
            int m2 = 2 * m;
            double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1 + aa * d; if (Math.Abs(d) < tiny) d = tiny;
            c = 1 + aa / c; if (Math.Abs(c) < tiny) c = tiny;
            d = 1 / d; h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1 + aa * d; if (Math.Abs(d) < tiny) d = tiny;
            c = 1 + aa / c; if (Math.Abs(c) < tiny) c = tiny;
            d = 1 / d;
            double del = d * c;
            h *= del;
            if (Math.Abs(del - 1) < eps) break;
        }
        return h;
    }

    /// <summary>Regularised lower incomplete gamma P(a, x).</summary>
    public static double RegularizedGammaP(double a, double x)
    {
        if (x <= 0) return 0;
        if (x < a + 1)
        {
            double sum = 1 / a, term = sum, ap = a;
            for (int n = 0; n < 1000; n++)
            {
                ap += 1; term *= x / ap; sum += term;
                if (Math.Abs(term) < Math.Abs(sum) * 1e-16) break;
            }
            return sum * Math.Exp(-x + a * Math.Log(x) - LogGamma(a));
        }
        return 1 - RegularizedGammaQContinued(a, x);
    }

    /// <summary>Regularised upper incomplete gamma Q(a, x) = 1 − P(a, x).</summary>
    public static double RegularizedGammaQ(double a, double x)
    {
        if (x <= 0) return 1;
        return x < a + 1 ? 1 - RegularizedGammaP(a, x) : RegularizedGammaQContinued(a, x);
    }

    private static double RegularizedGammaQContinued(double a, double x)
    {
        const double tiny = 1e-300;
        double b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
        for (int i = 1; i < 1000; i++)
        {
            double an = -i * (i - a);
            b += 2;
            d = an * d + b; if (Math.Abs(d) < tiny) d = tiny;
            c = b + an / c; if (Math.Abs(c) < tiny) c = tiny;
            d = 1 / d;
            double del = d * c;
            h *= del;
            if (Math.Abs(del - 1) < 1e-16) break;
        }
        return Math.Exp(-x + a * Math.Log(x) - LogGamma(a)) * h;
    }

    /// <summary>Two-sided p-value of a Student t statistic with <paramref name="df"/> degrees of freedom.</summary>
    public static double StudentTTwoSided(double t, double df)
    {
        if (double.IsNaN(t) || double.IsNaN(df) || df <= 0) return double.NaN;
        if (double.IsInfinity(t)) return 0;
        return RegularizedIncompleteBeta(df / 2, 0.5, df / (df + t * t));
    }

    /// <summary>Upper tail P(X ≥ x) of a chi-square distribution with <paramref name="df"/> degrees of freedom.</summary>
    public static double ChiSquareUpper(double x, double df) => x <= 0 ? 1 : RegularizedGammaQ(df / 2, x / 2);

    /// <summary>Standard normal cumulative distribution Φ(z).</summary>
    public static double NormalCdf(double z) => 0.5 * Erfc(-z / Math.Sqrt(2));

    /// <summary>Complementary error function, via the incomplete gamma function (full double precision).</summary>
    public static double Erfc(double x)
    {
        // erfc(x) = Q(1/2, x^2) for x ≥ 0; the gamma route keeps full double precision.
        if (x >= 0) return RegularizedGammaQ(0.5, x * x);
        return 2 - RegularizedGammaQ(0.5, x * x);
    }

    /// <summary>Inverse standard normal CDF (Acklam's algorithm with one Halley refinement step).</summary>
    public static double NormalQuantile(double p)
    {
        if (p <= 0) return double.NegativeInfinity;
        if (p >= 1) return double.PositiveInfinity;
        double[] a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
        double[] b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
        double[] c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
        double[] d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
        const double plow = 0.02425, phigh = 1 - plow;
        double x;
        if (p < plow)
        {
            double q = Math.Sqrt(-2 * Math.Log(p));
            x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        }
        else if (p <= phigh)
        {
            double q = p - 0.5, r = q * q;
            x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
        }
        else
        {
            double q = Math.Sqrt(-2 * Math.Log(1 - p));
            x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        }
        double e = NormalCdf(x) - p;
        double u = e * Math.Sqrt(2 * Math.PI) * Math.Exp(x * x / 2);
        return x - u / (1 + x * u / 2);
    }
}
