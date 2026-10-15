using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Population;

/// <summary>Diploid genotype counts at a bi-allelic site.</summary>
/// <param name="HomRef">Homozygous reference (0/0).</param>
/// <param name="Het">Heterozygous (0/1).</param>
/// <param name="HomAlt">Homozygous alternate (1/1).</param>
/// <param name="Missing">Calls with a missing allele.</param>
public readonly record struct GenotypeCounts(int HomRef, int Het, int HomAlt, int Missing = 0)
{
    /// <summary>Called individuals.</summary>
    public int Called => HomRef + Het + HomAlt;

    /// <summary>Called alleles (2 per individual).</summary>
    public int AlleleNumber => 2 * Called;

    /// <summary>Alternate allele count.</summary>
    public int AlternateCount => Het + 2 * HomAlt;

    /// <summary>Alternate allele frequency (NaN when nothing is called).</summary>
    public double AlternateFrequency => Called == 0 ? double.NaN : (double)AlternateCount / AlleleNumber;

    /// <summary>Reference allele frequency.</summary>
    public double ReferenceFrequency => 1 - AlternateFrequency;

    /// <summary>Observed heterozygosity, Het / Called.</summary>
    public double ObservedHeterozygosity => Called == 0 ? double.NaN : (double)Het / Called;

    /// <summary>Expected heterozygosity under HWE, 2pq.</summary>
    public double ExpectedHeterozygosity => 2 * AlternateFrequency * ReferenceFrequency;

    /// <summary>Inbreeding coefficient F = 1 − Hobs/Hexp.</summary>
    public double InbreedingCoefficient => 1 - ObservedHeterozygosity / ExpectedHeterozygosity;

    /// <summary>Hardy–Weinberg exact test p-value (Wigginton et al. 2005).</summary>
    public double HardyWeinbergP() => HardyWeinberg.ExactTest(Het, HomRef, HomAlt);

    /// <summary>
    /// Counts genotypes of <paramref name="variant"/> (diploid; any non-reference allele counts as alternate).
    /// <paramref name="samples"/> restricts the count to a sub-population (sample indices).
    /// </summary>
    public static GenotypeCounts From(Variant variant, IEnumerable<int>? samples = null)
    {
        int hr = 0, het = 0, ha = 0, miss = 0;
        var gts = variant.Genotypes;
        foreach (int i in samples ?? Enumerable.Range(0, gts.Count))
        {
            var a = gts[i].Alleles;
            if (a.Length < 2 || a[0] < 0 || a[1] < 0) { miss++; continue; }
            int alt = (a[0] > 0 ? 1 : 0) + (a[1] > 0 ? 1 : 0);
            if (alt == 0) hr++; else if (alt == 1) het++; else ha++;
        }
        return new(hr, het, ha, miss);
    }
}

/// <summary>Hardy–Weinberg equilibrium tests.</summary>
public static class HardyWeinberg
{
    /// <summary>
    /// Exact test of Hardy–Weinberg equilibrium for a bi-allelic marker (Wigginton, Cutler &amp; Abecasis 2005,
    /// <i>Am J Hum Genet</i> 76:887). Returns the two-sided p-value: the summed probability of all heterozygote
    /// counts no more likely than the observed one, given the allele counts.
    /// </summary>
    public static double ExactTest(int hets, int homRef, int homAlt)
    {
        if (hets < 0 || homRef < 0 || homAlt < 0) throw new ArgumentOutOfRangeException(nameof(hets), "Counts must be non-negative.");
        int obsHomr = Math.Min(homRef, homAlt), obsHomc = Math.Max(homRef, homAlt);
        int rare = 2 * obsHomr + hets;
        int genotypes = hets + obsHomc + obsHomr;
        if (genotypes == 0) return 1.0;
        var probs = new double[rare + 1];
        int mid = (int)((long)rare * (2L * genotypes - rare) / (2L * genotypes));
        if ((rare & 1) != (mid & 1)) mid++;

        int currHomr = (rare - mid) / 2, currHomc = genotypes - mid - currHomr;
        probs[mid] = 1.0;
        double sum = 1.0;
        for (int h = mid; h > 1; h -= 2)
        {
            probs[h - 2] = probs[h] * h * (h - 1.0) / (4.0 * (currHomr + 1.0) * (currHomc + 1.0));
            sum += probs[h - 2];
            currHomr++; currHomc++;
        }
        currHomr = (rare - mid) / 2; currHomc = genotypes - mid - currHomr;
        for (int h = mid; h <= rare - 2; h += 2)
        {
            probs[h + 2] = probs[h] * 4.0 * currHomr * currHomc / ((h + 2.0) * (h + 1.0));
            sum += probs[h + 2];
            currHomr--; currHomc--;
        }
        double target = probs[hets], p = 0;
        // Relative tolerance guards against floating-point ties.
        for (int i = 0; i <= rare; i++)
            if (probs[i] <= target * (1 + 1e-7)) p += probs[i];
        return Math.Min(1.0, p / sum);
    }

    /// <summary>Exact test from <see cref="GenotypeCounts"/>.</summary>
    public static double ExactTest(GenotypeCounts counts) => ExactTest(counts.Het, counts.HomRef, counts.HomAlt);

    /// <summary>Pearson χ² goodness-of-fit test (1 degree of freedom); returns the statistic and p-value.</summary>
    public static (double ChiSquare, double P) ChiSquareTest(int hets, int homRef, int homAlt)
    {
        double n = hets + homRef + homAlt;
        if (n == 0) return (0, 1);
        double p = (2.0 * homRef + hets) / (2 * n), q = 1 - p;
        double eRR = p * p * n, eRA = 2 * p * q * n, eAA = q * q * n;
        double chi = 0;
        if (eRR > 0) chi += Sq(homRef - eRR) / eRR;
        if (eRA > 0) chi += Sq(hets - eRA) / eRA;
        if (eAA > 0) chi += Sq(homAlt - eAA) / eAA;
        return (chi, SpecialFunctions.ChiSquareSurvival1(chi));
    }

    private static double Sq(double x) => x * x;
}

/// <summary>Allele frequencies and heterozygosity.</summary>
public static class AlleleFrequencies
{
    /// <summary>Frequency of every allele (index 0 = reference) over called alleles of <paramref name="variant"/>.</summary>
    public static double[] Of(Variant variant, IEnumerable<int>? samples = null)
    {
        var (counts, total) = Counts(variant, samples);
        var f = new double[counts.Length];
        for (int i = 0; i < f.Length; i++) f[i] = total == 0 ? double.NaN : (double)counts[i] / total;
        return f;
    }

    /// <summary>Allele counts (index 0 = reference) and the number of called alleles.</summary>
    public static (int[] Counts, int Called) Counts(Variant variant, IEnumerable<int>? samples = null)
    {
        var counts = new int[variant.Alternates.Count + 1];
        int total = 0;
        var gts = variant.Genotypes;
        foreach (int i in samples ?? Enumerable.Range(0, gts.Count))
            foreach (int a in gts[i].Alleles)
            {
                if (a < 0 || a >= counts.Length) continue;
                counts[a]++;
                total++;
            }
        return (counts, total);
    }

    /// <summary>
    /// Expected heterozygosity (gene diversity) 1 − Σp². With <paramref name="alleleNumber"/> (2N sampled alleles)
    /// Nei's unbiased estimator 2N/(2N−1)·(1 − Σp²) is returned.
    /// </summary>
    public static double ExpectedHeterozygosity(ReadOnlySpan<double> frequencies, int alleleNumber = 0)
    {
        double s = 0;
        foreach (double p in frequencies) s += p * p;
        double h = 1 - s;
        return alleleNumber > 1 ? h * alleleNumber / (alleleNumber - 1.0) : h;
    }

    /// <summary>Observed heterozygosity: fraction of called diploid genotypes that are heterozygous.</summary>
    public static double ObservedHeterozygosity(Variant variant, IEnumerable<int>? samples = null)
    {
        int het = 0, called = 0;
        foreach (int i in samples ?? Enumerable.Range(0, variant.Genotypes.Count))
        {
            var g = variant.Genotypes[i];
            if (g.IsMissing || g.Alleles.Any(a => a < 0)) continue;
            called++;
            if (g.IsHeterozygous) het++;
        }
        return called == 0 ? double.NaN : (double)het / called;
    }
}

/// <summary>Numerical helpers.</summary>
internal static class SpecialFunctions
{
    /// <summary>Complementary error function (Numerical Recipes Chebyshev approximation, |error| &lt; 1.2e-7).</summary>
    public static double Erfc(double x)
    {
        double z = Math.Abs(x);
        double t = 1.0 / (1.0 + 0.5 * z);
        double ans = t * Math.Exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
            t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
            t * (-0.82215223 + t * 0.17087277)))))))));
        return x >= 0 ? ans : 2.0 - ans;
    }

    /// <summary>Upper tail of χ² with one degree of freedom.</summary>
    public static double ChiSquareSurvival1(double x) => x <= 0 ? 1.0 : Erfc(Math.Sqrt(x / 2));

    /// <summary>log10 of the sum of 10^values (stable).</summary>
    public static double Log10SumExp(ReadOnlySpan<double> values)
    {
        double max = double.NegativeInfinity;
        foreach (double v in values) if (v > max) max = v;
        if (double.IsNegativeInfinity(max)) return max;
        double s = 0;
        foreach (double v in values) s += Math.Pow(10, v - max);
        return max + Math.Log10(s);
    }
}
