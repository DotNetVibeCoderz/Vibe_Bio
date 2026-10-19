namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>
/// R's default random number generation, reproduced bit for bit: <c>set.seed(seed)</c> with the Mersenne-Twister
/// generator (R's seeding scramble), <c>unif_rand</c>, inversion <c>norm_rand</c> (two uniforms, AS 241 <c>qnorm</c>),
/// <c>exp_rand</c>, <c>rgamma</c> (Ahrens–Dieter GD / GS) and <c>rchisq</c>. Used by DESeq2's simulation-based estimate of
/// the dispersion prior variance, which calls <c>set.seed(2)</c>.
/// </summary>
internal sealed class RRandom
{
    private const int N = 624, M = 397;
    private const uint MatrixA = 0x9908b0df, UpperMask = 0x80000000, LowerMask = 0x7fffffff;
    private readonly uint[] _mt = new uint[N];
    private int _mti;
    private double _aa, _aaa, _s, _s2, _d, _q0, _b, _si, _c;

    /// <summary><c>set.seed(seed)</c>.</summary>
    public RRandom(int seed)
    {
        uint s = unchecked((uint)seed);
        for (int j = 0; j < 50; j++) s = unchecked(69069 * s + 1);
        // i_seed[0] is mti (overwritten by FixupSeeds), i_seed[1..624] the state.
        s = unchecked(69069 * s + 1);
        for (int j = 0; j < N; j++) { s = unchecked(69069 * s + 1); _mt[j] = s; }
        _mti = N;
    }

    private double MtGenrand()
    {
        uint y;
        if (_mti >= N)
        {
            int kk;
            for (kk = 0; kk < N - M; kk++)
            {
                y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                _mt[kk] = _mt[kk + M] ^ (y >> 1) ^ ((y & 1) != 0 ? MatrixA : 0);
            }
            for (; kk < N - 1; kk++)
            {
                y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                _mt[kk] = _mt[kk + (M - N)] ^ (y >> 1) ^ ((y & 1) != 0 ? MatrixA : 0);
            }
            y = (_mt[N - 1] & UpperMask) | (_mt[0] & LowerMask);
            _mt[N - 1] = _mt[M - 1] ^ (y >> 1) ^ ((y & 1) != 0 ? MatrixA : 0);
            _mti = 0;
        }
        y = _mt[_mti++];
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        return y * 2.3283064365386963e-10;
    }

    /// <summary><c>unif_rand()</c>.</summary>
    public double Unif()
    {
        const double i2_32m1 = 2.328306437080797e-10;
        double x = MtGenrand();
        if (x <= 0.0) return 0.5 * i2_32m1;
        if (1.0 - x <= 0.0) return 1.0 - 0.5 * i2_32m1;
        return x;
    }

    /// <summary><c>norm_rand()</c> with R's default "Inversion".</summary>
    public double Norm()
    {
        const double big = 134217728; // 2^27
        double u = Unif();
        u = (int)(big * u) + Unif();
        return Qnorm5(u / big);
    }

    /// <summary><c>rnorm(1, mu, sigma)</c> (no draw when sigma is 0).</summary>
    public double Rnorm(double mu, double sigma)
    {
        if (sigma == 0 || !double.IsFinite(mu)) return mu;
        return mu + sigma * Norm();
    }

    private static readonly double[] ExpQ =
    [
        0.6931471805599453, 0.9333736875190459, 0.9888777961838675, 0.9984959252914960, 0.9998292811061389,
        0.9999833164100727, 0.9999985691438767, 0.9999998906925558, 0.9999999924734159, 0.9999999995283275,
        0.9999999999728814, 0.9999999999985598, 0.9999999999999289, 0.9999999999999968, 0.9999999999999999,
        1.0000000000000000,
    ];

    /// <summary><c>exp_rand()</c> (Ahrens &amp; Dieter 1972).</summary>
    public double Exp()
    {
        double a = 0.0;
        double u = Unif();
        while (u <= 0.0 || u >= 1.0) u = Unif();
        for (; ; )
        {
            u += u;
            if (u > 1.0) break;
            a += ExpQ[0];
        }
        u -= 1.0;
        if (u <= ExpQ[0]) return a + u;
        int i = 0;
        double ustar = Unif(), umin = ustar;
        do
        {
            ustar = Unif();
            if (umin > ustar) umin = ustar;
            i++;
        } while (u > ExpQ[i]);
        return a + umin * ExpQ[0];
    }

    /// <summary><c>rgamma(1, shape = a, scale = scale)</c>.</summary>
    public double Rgamma(double a, double scale)
    {
        const double sqrt32 = 5.656854, expM1 = 0.36787944117144232159;
        const double q1 = 0.04166669, q2 = 0.02083148, q3 = 0.00801191, q4 = 0.00144121, q5 = -7.388e-5, q6 = 2.4511e-4, q7 = 2.424e-4;
        const double a1 = 0.3333333, a2 = -0.250003, a3 = 0.2000062, a4 = -0.1662921, a5 = 0.1423657, a6 = -0.1367177, a7 = 0.1233795;
        double e, p, q, r, t, u, v, w, x, retVal;
        if (double.IsNaN(a) || double.IsNaN(scale)) return double.NaN;
        if (a <= 0.0 || scale <= 0.0) return scale == 0 || a == 0 ? 0 : double.NaN;
        if (!double.IsFinite(a) || !double.IsFinite(scale)) return double.PositiveInfinity;
        if (a < 1.0)
        {
            e = 1.0 + expM1 * a;
            for (; ; )
            {
                p = e * Unif();
                if (p >= 1.0)
                {
                    x = -Math.Log((e - p) / a);
                    if (Exp() >= (1.0 - a) * Math.Log(x)) break;
                }
                else
                {
                    x = Math.Exp(Math.Log(p) / a);
                    if (Exp() >= x) break;
                }
            }
            return scale * x;
        }
        if (a != _aa)
        {
            _aa = a;
            _s2 = a - 0.5;
            _s = Math.Sqrt(_s2);
            _d = sqrt32 - _s * 12.0;
        }
        t = Norm();
        x = _s + 0.5 * t;
        retVal = x * x;
        if (t >= 0.0) return scale * retVal;
        u = Unif();
        if (_d * u <= t * t * t) return scale * retVal;
        if (a != _aaa)
        {
            _aaa = a;
            r = 1.0 / a;
            _q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r + q1) * r;
            if (a <= 3.686) { _b = 0.463 + _s + 0.178 * _s2; _si = 1.235; _c = 0.195 / _s - 0.079 + 0.16 * _s; }
            else if (a <= 13.022) { _b = 1.654 + 0.0076 * _s2; _si = 1.68 / _s + 0.275; _c = 0.062 / _s + 0.024; }
            else { _b = 1.77; _si = 0.75; _c = 0.1515 / _s; }
        }
        if (x > 0.0)
        {
            v = t / (_s + _s);
            if (Math.Abs(v) <= 0.25)
                q = _q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
            else
                q = _q0 - _s * t + 0.25 * t * t + (_s2 + _s2) * Math.Log(1.0 + v);
            if (Math.Log(1.0 - u) <= q) return scale * retVal;
        }
        for (; ; )
        {
            e = Exp();
            u = Unif();
            u = u + u - 1.0;
            t = u < 0.0 ? _b - _si * e : _b + _si * e;
            if (t >= -0.71874483771719)
            {
                v = t / (_s + _s);
                if (Math.Abs(v) <= 0.25)
                    q = _q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
                else
                    q = _q0 - _s * t + 0.25 * t * t + (_s2 + _s2) * Math.Log(1.0 + v);
                if (q > 0.0)
                {
                    double th = Math.Tanh(q / 2);
                    w = q < 1 ? 2 * th / (1 - th) : Math.Exp(q) - 1; // expm1(q) without cancellation
                    if (_c * Math.Abs(u) <= w * Math.Exp(e - 0.5 * t * t)) break;
                }
            }
        }
        x = _s + 0.5 * t;
        return scale * x * x;
    }

    /// <summary><c>rchisq(1, df)</c> = <c>rgamma(df / 2, 2)</c>.</summary>
    public double Rchisq(double df) => !double.IsFinite(df) || df < 0.0 ? double.NaN : Rgamma(df / 2.0, 2.0);

    /// <summary>nmath <c>qnorm5(p, 0, 1, lower = TRUE, log = FALSE)</c> (Wichura's AS 241).</summary>
    public static double Qnorm5(double p)
    {
        if (double.IsNaN(p)) return p;
        if (p <= 0) return p == 0 ? double.NegativeInfinity : double.NaN;
        if (p >= 1) return p == 1 ? double.PositiveInfinity : double.NaN;
        double q = p - 0.5, r, val;
        if (Math.Abs(q) <= .425)
        {
            r = .180625 - q * q;
            return q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r +
                             45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
                         133.14166789178437745) * r + 3.387132872796366608)
                   / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r + 39307.89580009271061) * r +
                          21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
                       42.313330701600911252) * r + 1.0);
        }
        double lp = Math.Log(q > 0 ? 1 - p : p);
        r = Math.Sqrt(-lp);
        if (r <= 5.0)
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 + .0227238449892691845833) * r + .24178072517745061177) * r +
                        1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
                     4.6303378461565452959) * r + 1.42343711074968357734)
                  / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + .0151986665636164571966) * r +
                         .14810397642748007459) * r + .68976733498510000455) * r + 1.6763848301838038494) * r +
                      2.05319162663775882187) * r + 1.0);
        }
        else
        {
            r += -5.0;
            val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + .0012426609473880784386) * r +
                        .026532189526576123093) * r + .29656057182850489123) * r + 1.7848265399172913358) * r +
                     5.4637849111641143699) * r + 6.6579046435011037772)
                  / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
                         7.868691311456132591e-4) * r + .0148753612908506148525) * r + .13692988092273580531) * r +
                      .59983220655588793769) * r + 1.0);
        }
        return q < 0.0 ? -val : val;
    }
}
