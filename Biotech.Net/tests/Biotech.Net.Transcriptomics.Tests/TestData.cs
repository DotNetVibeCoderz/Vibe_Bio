using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>Seeded simulators shared by the tests (gamma, Poisson and negative-binomial draws).</summary>
internal static class TestData
{
    public static double Normal(Random r)
    {
        double u1 = 1 - r.NextDouble(), u2 = r.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }

    /// <summary>Gamma(shape, scale) by Marsaglia–Tsang.</summary>
    public static double Gamma(Random r, double shape, double scale = 1)
    {
        if (shape < 1) return Gamma(r, shape + 1, scale) * Math.Pow(r.NextDouble(), 1 / shape);
        double d = shape - 1.0 / 3, c = 1 / Math.Sqrt(9 * d);
        while (true)
        {
            double x, v;
            do { x = Normal(r); v = 1 + c * x; } while (v <= 0);
            v = v * v * v;
            double u = r.NextDouble();
            if (u < 1 - 0.0331 * x * x * x * x) return d * v * scale;
            if (Math.Log(u) < 0.5 * x * x + d * (1 - v + Math.Log(v))) return d * v * scale;
        }
    }

    public static double ChiSquare(Random r, double df) => Gamma(r, df / 2, 2);

    /// <summary>Poisson draw: Knuth for small means, numpy's PTRS (Hörmann 1993) for larger ones.</summary>
    public static int Poisson(Random r, double lambda)
    {
        if (lambda <= 0) return 0;
        if (lambda < 10)
        {
            double l = Math.Exp(-lambda), p = 1;
            int k = 0;
            do { k++; p *= r.NextDouble(); } while (p > l);
            return k - 1;
        }
        double slam = Math.Sqrt(lambda), loglam = Math.Log(lambda);
        double b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
        double invAlpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2);
        while (true)
        {
            double u = r.NextDouble() - 0.5, v = r.NextDouble();
            double us = 0.5 - Math.Abs(u);
            int k = (int)Math.Floor((2 * a / us + b) * u + lambda + 0.43);
            if (us >= 0.07 && v <= vr) return k;
            if (k < 0 || (us < 0.013 && v > us)) continue;
            if (Math.Log(v) + Math.Log(invAlpha) - Math.Log(a / (us * us) + b) <= -lambda + k * loglam - SpecialFunctions.LogGamma(k + 1.0))
                return k;
        }
    }

    /// <summary>Negative binomial with mean <paramref name="mu"/> and dispersion α (variance μ + αμ²) as gamma–Poisson.</summary>
    public static int NegativeBinomial(Random r, double mu, double alpha)
        => alpha <= 0 ? Poisson(r, mu) : Poisson(r, Gamma(r, 1 / alpha, mu * alpha));

    /// <summary>
    /// A two-group RNA-Seq experiment: <paramref name="genes"/> genes, <paramref name="perGroup"/> samples per group,
    /// NB counts with dispersion 0.05 + 1/μ. The first <paramref name="deGenes"/> genes change by ±<paramref name="lfc"/>
    /// log2 units (alternating sign). Column annotation "condition" = A (reference) / B.
    /// </summary>
    public static (Dataset Counts, double[] TrueLfc, double[] SizeFactors) TwoGroupCounts(int genes, int perGroup, int deGenes,
        double lfc, int seed, double asymptDisp = 0.05, double extraPois = 1.0)
    {
        var r = new Random(seed);
        int m = 2 * perGroup;
        var sf = new double[m];
        for (int j = 0; j < m; j++) sf[j] = Math.Exp(0.4 * Normal(r));
        var trueLfc = new double[genes];
        var values = new double[genes * m];
        for (int i = 0; i < genes; i++)
        {
            double baseMu = Math.Exp(Math.Log(20) + 2.0 * r.NextDouble() * Math.Log(10) * 0.8); // ~20 .. ~1600
            if (i < deGenes) trueLfc[i] = i % 2 == 0 ? lfc : -lfc;
            for (int j = 0; j < m; j++)
            {
                double mu = baseMu * sf[j] * (j >= perGroup ? Math.Pow(2, trueLfc[i] / 2) : Math.Pow(2, -trueLfc[i] / 2));
                double alpha = asymptDisp + extraPois / (baseMu);
                values[i * m + j] = NegativeBinomial(r, mu, alpha);
            }
        }
        var geneNames = Enumerable.Range(0, genes).Select(i => $"g{i}").ToArray();
        var samples = Enumerable.Range(0, m).Select(j => $"s{j}").ToArray();
        var d = new Dataset("sim", geneNames, samples, values);
        d.ColumnAnnotations["condition"] = Enumerable.Range(0, m).Select(j => j < perGroup ? "A" : "B").ToArray();
        return (d, trueLfc, sf);
    }

    /// <summary>Adjusted Rand index of two labelings.</summary>
    public static double AdjustedRandIndex(IReadOnlyList<int> a, IReadOnlyList<int> b)
    {
        int n = a.Count;
        var table = new Dictionary<(int, int), long>();
        var ra = new Dictionary<int, long>(); var rb = new Dictionary<int, long>();
        for (int i = 0; i < n; i++)
        {
            table[(a[i], b[i])] = table.GetValueOrDefault((a[i], b[i])) + 1;
            ra[a[i]] = ra.GetValueOrDefault(a[i]) + 1;
            rb[b[i]] = rb.GetValueOrDefault(b[i]) + 1;
        }
        static double C2(long x) => x * (x - 1) / 2.0;
        double index = table.Values.Sum(C2), sa = ra.Values.Sum(C2), sb = rb.Values.Sum(C2);
        double expected = sa * sb / C2(n), max = 0.5 * (sa + sb);
        return (index - expected) / (max - expected);
    }
}
