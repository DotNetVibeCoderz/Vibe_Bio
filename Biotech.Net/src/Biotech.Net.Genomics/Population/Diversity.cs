using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Population;

/// <summary>Summary of nucleotide variation in a sample of aligned sequences.</summary>
/// <param name="SampleSize">Number of sequences n.</param>
/// <param name="Sites">Alignment columns analysed (columns with a gap or ambiguity code are excluded).</param>
/// <param name="SegregatingSites">Polymorphic columns S.</param>
/// <param name="MeanPairwiseDifferences">Average number of differences between two sequences, k (Tajima's θπ).</param>
/// <param name="ThetaW">Watterson's θ = S / a₁ (per sequence).</param>
/// <param name="TajimasD">Tajima's D (NaN when S = 0).</param>
public sealed record DiversitySummary(int SampleSize, int Sites, int SegregatingSites, double MeanPairwiseDifferences, double ThetaW, double TajimasD)
{
    /// <summary>Nucleotide diversity π per site (Nei &amp; Li 1979).</summary>
    public double Pi => Sites == 0 ? double.NaN : MeanPairwiseDifferences / Sites;

    /// <summary>Watterson's θ per site.</summary>
    public double ThetaWPerSite => Sites == 0 ? double.NaN : ThetaW / Sites;
}

/// <summary>Nucleotide diversity, Watterson's θ and Tajima's D.</summary>
public static class Diversity
{
    /// <summary>a₁ = Σ_{i=1}^{n−1} 1/i.</summary>
    public static double A1(int n)
    {
        double s = 0;
        for (int i = 1; i < n; i++) s += 1.0 / i;
        return s;
    }

    /// <summary>a₂ = Σ_{i=1}^{n−1} 1/i².</summary>
    public static double A2(int n)
    {
        double s = 0;
        for (int i = 1; i < n; i++) s += 1.0 / ((double)i * i);
        return s;
    }

    /// <summary>Watterson's estimator θ_W = S / a₁ (per sequence; divide by length for per site).</summary>
    public static double WattersonTheta(int segregatingSites, int sampleSize)
    {
        if (sampleSize < 2) throw new ArgumentOutOfRangeException(nameof(sampleSize), "Need at least two sequences.");
        return segregatingSites / A1(sampleSize);
    }

    /// <summary>
    /// Tajima's D (Tajima 1989, <i>Genetics</i> 123:585) from the sample size, the number of segregating sites and
    /// the mean number of pairwise differences. NaN when S = 0.
    /// </summary>
    public static double TajimasD(int sampleSize, int segregatingSites, double meanPairwiseDifferences)
    {
        int n = sampleSize, s = segregatingSites;
        if (n < 4) throw new ArgumentOutOfRangeException(nameof(sampleSize), "Tajima's D needs at least four sequences.");
        if (s == 0) return double.NaN;
        double a1 = A1(n), a2 = A2(n);
        double b1 = (n + 1.0) / (3.0 * (n - 1));
        double b2 = 2.0 * (n * (double)n + n + 3) / (9.0 * n * (n - 1));
        double c1 = b1 - 1 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1, e2 = c2 / (a1 * a1 + a2);
        return (meanPairwiseDifferences - s / a1) / Math.Sqrt(e1 * s + e2 * s * (s - 1.0));
    }

    /// <summary>
    /// Diversity summary of aligned sequences of equal length. Columns containing anything other than A/C/G/T in any
    /// sequence are excluded ("complete deletion").
    /// </summary>
    public static DiversitySummary Summarize(IReadOnlyList<Sequence> alignment)
    {
        int n = alignment.Count;
        if (n < 2) throw new ArgumentException("Need at least two sequences.", nameof(alignment));
        int len = alignment[0].Length;
        foreach (var s in alignment)
            if (s.Length != len) throw new ArgumentException("Sequences must be aligned (equal length).", nameof(alignment));
        int sites = 0, seg = 0;
        double diffs = 0;
        Span<int> c = stackalloc int[4];
        double pairs = n * (n - 1) / 2.0;
        for (int j = 0; j < len; j++)
        {
            c.Clear();
            bool ok = true;
            for (int i = 0; i < n; i++)
            {
                int code = BaseCode(alignment[i][j]);
                if (code < 0) { ok = false; break; }
                c[code]++;
            }
            if (!ok) continue;
            sites++;
            int distinct = 0;
            double same = 0;
            for (int b = 0; b < 4; b++)
            {
                if (c[b] > 0) distinct++;
                same += c[b] * (c[b] - 1) / 2.0;
            }
            if (distinct > 1) seg++;
            diffs += pairs - same;
        }
        double k = diffs / pairs;
        double tw = WattersonTheta(seg, n);
        double d = n >= 4 ? TajimasD(n, seg, k) : double.NaN;
        return new DiversitySummary(n, sites, seg, k, tw, d);
    }

    /// <summary>Nucleotide diversity π per site of an alignment.</summary>
    public static double NucleotideDiversity(IReadOnlyList<Sequence> alignment) => Summarize(alignment).Pi;

    /// <summary>
    /// Mean pairwise differences summed over bi-allelic sites from allele counts: Σ 2·j·(n−j)/(n(n−1)) where j is the
    /// alternate count and n the number of called alleles at the site. Divide by the sequence length for π per site.
    /// </summary>
    public static double MeanPairwiseDifferences(IEnumerable<(int AlternateCount, int AlleleNumber)> sites)
    {
        double k = 0;
        foreach (var (j, n) in sites)
            if (n > 1) k += 2.0 * j * (n - j) / (n * (n - 1.0));
        return k;
    }

    /// <summary>
    /// π per site from VCF variants (bi-allelic sites, called alleles only) over a region of
    /// <paramref name="sequenceLength"/> bases.
    /// </summary>
    public static double NucleotideDiversity(IEnumerable<Variant> variants, long sequenceLength, IEnumerable<int>? samples = null)
    {
        var list = samples?.ToArray();
        return MeanPairwiseDifferences(variants.Select(v =>
        {
            var (counts, called) = AlleleFrequencies.Counts(v, list);
            return (called - counts[0], called);
        })) / sequenceLength;
    }

    /// <summary>Tajima's D from VCF variants (every sample must be called at every site; incomplete sites are skipped).</summary>
    public static double TajimasD(IEnumerable<Variant> variants, IEnumerable<int>? samples = null)
    {
        var list = samples?.ToArray();
        int n = -1, s = 0;
        double k = 0;
        foreach (var v in variants)
        {
            var (counts, called) = AlleleFrequencies.Counts(v, list);
            int expected = (list ?? Enumerable.Range(0, v.Genotypes.Count)).Sum(i => v.Genotypes[i].Alleles.Length);
            if (called != expected || called < 2) continue;
            if (n < 0) n = called;
            if (called != n) continue;
            int j = called - counts[0];
            if (j == 0 || j == called) continue;
            s++;
            k += 2.0 * j * (called - j) / (called * (called - 1.0));
        }
        return n < 4 ? double.NaN : TajimasD(n, s, k);
    }

    /// <summary>Number of differences between two aligned sequences over columns where both have A/C/G/T.</summary>
    public static int PairwiseDifferences(Sequence a, Sequence b, out int comparedSites)
    {
        int len = Math.Min(a.Length, b.Length), d = 0;
        comparedSites = 0;
        for (int i = 0; i < len; i++)
        {
            int x = BaseCode(a[i]), y = BaseCode(b[i]);
            if (x < 0 || y < 0) continue;
            comparedSites++;
            if (x != y) d++;
        }
        return d;
    }

    internal static int BaseCode(byte b) => (char)(b & 0xDF) switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' or 'U' => 3, _ => -1 };
}

/// <summary>Population differentiation (F<sub>ST</sub>) estimators.</summary>
public static class Fst
{
    /// <summary>
    /// Hudson, Slatkin &amp; Maddison (1992) F<sub>ST</sub> = 1 − H<sub>w</sub>/H<sub>b</sub> from aligned sequences of two
    /// populations: H<sub>w</sub> is the mean of the within-population mean pairwise differences and H<sub>b</sub> the mean
    /// pairwise difference between sequences drawn from different populations (per site, gap/N columns excluded per pair).
    /// </summary>
    public static double Hudson(IReadOnlyList<Sequence> population1, IReadOnlyList<Sequence> population2)
    {
        double w1 = MeanDiff(population1, population1, within: true);
        double w2 = MeanDiff(population2, population2, within: true);
        double b = MeanDiff(population1, population2, within: false);
        return b == 0 ? 0 : 1 - (w1 + w2) / 2 / b;
    }

    private static double MeanDiff(IReadOnlyList<Sequence> x, IReadOnlyList<Sequence> y, bool within)
    {
        double sum = 0;
        int pairs = 0;
        for (int i = 0; i < x.Count; i++)
            for (int j = within ? i + 1 : 0; j < y.Count; j++)
            {
                int d = Diversity.PairwiseDifferences(x[i], y[j], out int sites);
                if (sites == 0) continue;
                sum += (double)d / sites;
                pairs++;
            }
        return pairs == 0 ? 0 : sum / pairs;
    }

    /// <summary>
    /// Hudson F<sub>ST</sub> from bi-allelic allele counts, as formulated by Bhatia et al. (2013, <i>Genome Res</i> 23:1514),
    /// combined over sites as a ratio of averages: N = (p₁−p₂)² − p₁(1−p₁)/(n₁−1) − p₂(1−p₂)/(n₂−1),
    /// D = p₁(1−p₂) + p₂(1−p₁), F<sub>ST</sub> = ΣN / ΣD.
    /// </summary>
    public static double Hudson(IEnumerable<(int Ac1, int An1, int Ac2, int An2)> sites)
    {
        double num = 0, den = 0;
        foreach (var (ac1, an1, ac2, an2) in sites)
        {
            if (an1 < 2 || an2 < 2) continue;
            double p1 = (double)ac1 / an1, p2 = (double)ac2 / an2;
            num += (p1 - p2) * (p1 - p2) - p1 * (1 - p1) / (an1 - 1) - p2 * (1 - p2) / (an2 - 1);
            den += p1 * (1 - p2) + p2 * (1 - p1);
        }
        return den == 0 ? double.NaN : num / den;
    }

    /// <summary>Hudson F<sub>ST</sub> between two sample groups of VCF variants (bi-allelic, any alternate counted).</summary>
    public static double Hudson(IEnumerable<Variant> variants, IReadOnlyList<int> population1, IReadOnlyList<int> population2)
        => Hudson(variants.Select(v =>
        {
            var (c1, n1) = AlleleFrequencies.Counts(v, population1);
            var (c2, n2) = AlleleFrequencies.Counts(v, population2);
            return (n1 - c1[0], n1, n2 - c2[0], n2);
        }));

    /// <summary>Per-locus Weir–Cockerham variance components.</summary>
    /// <param name="A">Between-population component.</param>
    /// <param name="B">Between individuals within populations.</param>
    /// <param name="C">Between gametes within individuals.</param>
    public readonly record struct VarianceComponents(double A, double B, double C)
    {
        /// <summary>θ for this locus alone.</summary>
        public double Theta => A / (A + B + C);
    }

    /// <summary>
    /// Weir &amp; Cockerham (1984, <i>Evolution</i> 38:1358) variance components for one bi-allelic locus from the diploid
    /// genotype counts of each population.
    /// </summary>
    public static VarianceComponents WeirCockerhamComponents(IReadOnlyList<GenotypeCounts> populations)
    {
        int r = populations.Count;
        if (r < 2) throw new ArgumentException("Need at least two populations.", nameof(populations));
        double sumN = 0, sumN2 = 0;
        foreach (var p in populations) { sumN += p.Called; sumN2 += (double)p.Called * p.Called; }
        double nBar = sumN / r;
        double nC = (r * nBar - sumN2 / (r * nBar)) / (r - 1);
        double pBar = 0, hBar = 0;
        foreach (var p in populations)
        {
            if (p.Called == 0) continue;
            pBar += p.Called * p.AlternateFrequency;
            hBar += p.Called * p.ObservedHeterozygosity;
        }
        pBar /= r * nBar;
        hBar /= r * nBar;
        double s2 = 0;
        foreach (var pop in populations)
        {
            if (pop.Called == 0) continue;
            double d = pop.AlternateFrequency - pBar;
            s2 += pop.Called * d * d;
        }
        s2 /= (r - 1) * nBar;
        double pq = pBar * (1 - pBar);
        double a = nBar / nC * (s2 - 1 / (nBar - 1) * (pq - (r - 1.0) / r * s2 - hBar / 4));
        double b = nBar / (nBar - 1) * (pq - (r - 1.0) / r * s2 - (2 * nBar - 1) / (4 * nBar) * hBar);
        double c = hBar / 2;
        return new VarianceComponents(a, b, c);
    }

    /// <summary>
    /// Weir &amp; Cockerham θ over loci: Σa / Σ(a + b + c). Each element of <paramref name="loci"/> holds the genotype
    /// counts of every population at one locus.
    /// </summary>
    public static double WeirCockerham(IEnumerable<IReadOnlyList<GenotypeCounts>> loci)
    {
        double num = 0, den = 0;
        foreach (var locus in loci)
        {
            var vc = WeirCockerhamComponents(locus);
            if (double.IsNaN(vc.A) || double.IsNaN(vc.B)) continue;
            num += vc.A;
            den += vc.A + vc.B + vc.C;
        }
        return den == 0 ? double.NaN : num / den;
    }

    /// <summary>Weir &amp; Cockerham θ between sample groups of VCF variants.</summary>
    public static double WeirCockerham(IEnumerable<Variant> variants, params IReadOnlyList<int>[] populations)
        => WeirCockerham(variants.Select(v => (IReadOnlyList<GenotypeCounts>)populations.Select(p => GenotypeCounts.From(v, p)).ToArray()));
}
