using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Population;

/// <summary>Pairwise linkage disequilibrium between two bi-allelic loci.</summary>
/// <param name="D">Coefficient D = p<sub>AB</sub> − p<sub>A</sub>p<sub>B</sub> (alternate alleles as A and B).</param>
/// <param name="DPrime">Lewontin's |D′| = |D| / D<sub>max</sub> (0–1).</param>
/// <param name="RSquared">r² = D² / (p<sub>A</sub>(1−p<sub>A</sub>)p<sub>B</sub>(1−p<sub>B</sub>)).</param>
/// <param name="Haplotypes">Haplotypes used.</param>
public readonly record struct LdResult(double D, double DPrime, double RSquared, int Haplotypes);

/// <summary>Linkage disequilibrium statistics.</summary>
public static class LinkageDisequilibrium
{
    /// <summary>
    /// LD from haplotypes: <paramref name="locusA"/>[i] and <paramref name="locusB"/>[i] are the alleles of haplotype i
    /// (0 = reference, &gt;0 = alternate, &lt;0 = missing and skipped).
    /// </summary>
    public static LdResult FromHaplotypes(ReadOnlySpan<int> locusA, ReadOnlySpan<int> locusB)
    {
        if (locusA.Length != locusB.Length) throw new ArgumentException("Both loci need the same number of haplotypes.");
        int n = 0, a = 0, b = 0, ab = 0;
        for (int i = 0; i < locusA.Length; i++)
        {
            if (locusA[i] < 0 || locusB[i] < 0) continue;
            n++;
            bool x = locusA[i] > 0, y = locusB[i] > 0;
            if (x) a++;
            if (y) b++;
            if (x && y) ab++;
        }
        return FromCounts(n, a, b, ab);
    }

    /// <summary>LD from haplotype counts: total haplotypes, carriers of A, carriers of B, and carriers of both.</summary>
    public static LdResult FromCounts(int haplotypes, int countA, int countB, int countAB)
    {
        if (haplotypes == 0) return new(double.NaN, double.NaN, double.NaN, 0);
        double pA = (double)countA / haplotypes, pB = (double)countB / haplotypes, pAB = (double)countAB / haplotypes;
        double d = pAB - pA * pB;
        double denom = pA * (1 - pA) * pB * (1 - pB);
        double r2 = denom == 0 ? double.NaN : d * d / denom;
        double dmax = d >= 0 ? Math.Min(pA * (1 - pB), (1 - pA) * pB) : Math.Min(pA * pB, (1 - pA) * (1 - pB));
        double dp = dmax == 0 ? double.NaN : Math.Abs(d) / dmax;
        return new(d, dp, r2, haplotypes);
    }

    /// <summary>
    /// LD between two variants from their phased genotypes (each allele slot of each sample is one haplotype).
    /// Unphased heterozygotes are skipped.
    /// </summary>
    public static LdResult FromPhased(Variant a, Variant b, IEnumerable<int>? samples = null)
    {
        var ha = new List<int>();
        var hb = new List<int>();
        foreach (int i in samples ?? Enumerable.Range(0, Math.Min(a.Genotypes.Count, b.Genotypes.Count)))
        {
            var ga = a.Genotypes[i];
            var gb = b.Genotypes[i];
            if (ga.Alleles.Length != gb.Alleles.Length) continue;
            bool phased = (ga.Phased || !ga.IsHeterozygous) && (gb.Phased || !gb.IsHeterozygous);
            if (!phased) continue;
            for (int k = 0; k < ga.Alleles.Length; k++) { ha.Add(ga.Alleles[k]); hb.Add(gb.Alleles[k]); }
        }
        return FromHaplotypes(ha.ToArray(), hb.ToArray());
    }

    /// <summary>
    /// Genotype (composite) r²: the squared Pearson correlation of alternate-allele dosages across samples. Works on
    /// unphased data; equals haplotype r² under HWE in expectation.
    /// </summary>
    public static double GenotypeRSquared(Variant a, Variant b, IEnumerable<int>? samples = null)
    {
        double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        int n = 0;
        foreach (int i in samples ?? Enumerable.Range(0, Math.Min(a.Genotypes.Count, b.Genotypes.Count)))
        {
            var ga = a.Genotypes[i];
            var gb = b.Genotypes[i];
            if (ga.Alleles.Any(x => x < 0) || gb.Alleles.Any(x => x < 0)) continue;
            double x = ga.AlternateCount, y = gb.AlternateCount;
            n++; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
        }
        if (n < 2) return double.NaN;
        double cov = sxy - sx * sy / n, vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
        return vx <= 0 || vy <= 0 ? double.NaN : cov * cov / (vx * vy);
    }
}

/// <summary>Site frequency spectrum (SFS) utilities.</summary>
public static class SiteFrequencySpectrum
{
    /// <summary>
    /// Unfolded SFS: <c>sfs[j]</c> = number of sites where the derived allele is carried by j of the
    /// <paramref name="sampleSize"/> chromosomes (j = 0..n).
    /// </summary>
    public static long[] Unfolded(IEnumerable<int> derivedCounts, int sampleSize)
    {
        var sfs = new long[sampleSize + 1];
        foreach (int j in derivedCounts)
        {
            if (j < 0 || j > sampleSize) throw new ArgumentOutOfRangeException(nameof(derivedCounts), $"Count {j} outside 0..{sampleSize}.");
            sfs[j]++;
        }
        return sfs;
    }

    /// <summary>Folded (minor-allele) SFS: <c>sfs[j]</c> for j = 0..⌊n/2⌋.</summary>
    public static long[] Folded(IEnumerable<int> alleleCounts, int sampleSize)
    {
        var sfs = new long[sampleSize / 2 + 1];
        foreach (int j in alleleCounts)
        {
            if (j < 0 || j > sampleSize) throw new ArgumentOutOfRangeException(nameof(alleleCounts));
            sfs[Math.Min(j, sampleSize - j)]++;
        }
        return sfs;
    }

    /// <summary>Folds an unfolded spectrum.</summary>
    public static long[] Fold(ReadOnlySpan<long> unfolded)
    {
        int n = unfolded.Length - 1;
        var f = new long[n / 2 + 1];
        for (int j = 0; j <= n; j++) f[Math.Min(j, n - j)] += unfolded[j];
        return f;
    }

    /// <summary>
    /// SFS from VCF variants, treating the alternate allele as derived (bi-allelic sites with every chromosome called;
    /// others are skipped). <paramref name="folded"/> returns the minor-allele spectrum.
    /// </summary>
    public static long[] FromVariants(IEnumerable<Variant> variants, bool folded = false, IEnumerable<int>? samples = null)
    {
        var list = samples?.ToArray();
        var counts = new List<int>();
        int n = -1;
        foreach (var v in variants)
        {
            if (v.Alternates.Count != 1) continue;
            var (c, called) = AlleleFrequencies.Counts(v, list);
            int expected = (list ?? Enumerable.Range(0, v.Genotypes.Count)).Sum(i => v.Genotypes[i].Alleles.Length);
            if (called != expected) continue;
            if (n < 0) n = called;
            if (called != n) continue;
            counts.Add(c[1]);
        }
        if (n < 0) return [];
        return folded ? Folded(counts, n) : Unfolded(counts, n);
    }

    /// <summary>θπ (mean pairwise differences) from an unfolded SFS: Σ j(n−j)·ξⱼ / C(n,2).</summary>
    public static double ThetaPi(ReadOnlySpan<long> unfolded)
    {
        int n = unfolded.Length - 1;
        double s = 0;
        for (int j = 1; j < n; j++) s += (double)j * (n - j) * unfolded[j];
        return s / (n * (n - 1) / 2.0);
    }

    /// <summary>Watterson's θ from an SFS (segregating sites / a₁).</summary>
    public static double ThetaW(ReadOnlySpan<long> unfolded)
    {
        int n = unfolded.Length - 1;
        long seg = 0;
        for (int j = 1; j < n; j++) seg += unfolded[j];
        return seg / Diversity.A1(n);
    }

    /// <summary>Tajima's D from an unfolded SFS.</summary>
    public static double TajimasD(ReadOnlySpan<long> unfolded)
    {
        int n = unfolded.Length - 1;
        long seg = 0;
        for (int j = 1; j < n; j++) seg += unfolded[j];
        return Diversity.TajimasD(n, (int)seg, ThetaPi(unfolded));
    }
}
