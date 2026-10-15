using System.Numerics;
using Biotech.Net.Transcriptomics.Enrichment;

namespace Biotech.Net.Transcriptomics.Tests;

public class EnrichmentTests
{
    private static BigInteger Choose(int n, int k)
    {
        if (k < 0 || k > n) return BigInteger.Zero;
        BigInteger r = BigInteger.One;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }

    /// <summary>Exact P(X ≥ k) by summing binomial coefficients with arbitrary-precision integers.</summary>
    private static double ExactUpperTail(int k, int bigK, int bigN, int n)
    {
        BigInteger num = BigInteger.Zero;
        for (int x = Math.Max(k, 0); x <= Math.Min(n, bigK); x++) num += Choose(bigK, x) * Choose(bigN - bigK, n - x);
        if (num.IsZero) return 0;
        return Math.Exp(BigInteger.Log(num) - BigInteger.Log(Choose(bigN, n)));
    }

    [Theory]
    [InlineData(5, 50, 1000, 40)]
    [InlineData(10, 200, 20000, 300)]
    [InlineData(3, 4, 10, 5)]
    [InlineData(1, 30, 100, 10)]
    [InlineData(25, 40, 15000, 200)]   // extreme tail
    [InlineData(0, 30, 100, 10)]
    public void Hypergeometric_MatchesExactSums(int k, int bigK, int bigN, int n)
    {
        double exact = ExactUpperTail(k, bigK, bigN, n);
        double p = OverRepresentation.HypergeometricUpperTail(k, bigK, bigN, n);
        Assert.True(Math.Abs(p - exact) <= 1e-8 * Math.Max(exact, 1e-300) + 1e-300, $"{p} vs {exact}");
    }

    [Fact]
    public void Hypergeometric_MatchesPhyper()
    {
        // phyper(4, 50, 950, 40, lower.tail = FALSE)
        Assert.Equal(0.04446989715343539, OverRepresentation.HypergeometricUpperTail(5, 50, 1000, 40), 1e-13);
        // phyper(2, 4, 6, 5, lower.tail = FALSE) = 66/252
        Assert.Equal(66.0 / 252, OverRepresentation.HypergeometricUpperTail(3, 4, 10, 5), 1e-14);
        Assert.Equal(0, OverRepresentation.HypergeometricUpperTail(5, 4, 10, 5));
    }

    [Fact]
    public void Ora_ComputesOverlapExpectedAndBh()
    {
        var universe = Enumerable.Range(0, 100).Select(i => $"G{i}").ToArray();
        var hits = Enumerable.Range(0, 10).Select(i => $"G{i}").Append("NOT_IN_UNIVERSE").ToArray();
        var sets = new[]
        {
            new GeneSet("enriched", Enumerable.Range(0, 8).Select(i => $"G{i}").Concat(["G50", "G51", "X1"])),   // K = 10, k = 8
            new GeneSet("neutral", Enumerable.Range(40, 20).Select(i => $"G{i}")),                                // K = 20, k = 0
            new GeneSet("partial", ["G9", "G60", "G61", "G62", "G63"]),                                           // K = 5,  k = 1
        };
        var res = OverRepresentation.Test(hits, universe, sets);
        Assert.Equal("enriched", res[0].Set);
        var e = res[0];
        Assert.Equal(10, e.SetSize);
        Assert.Equal(8, e.Overlap);
        Assert.Equal(1.0, e.Expected, 1e-12);          // n·K/N = 10·10/100
        Assert.Equal(8.0, e.FoldEnrichment, 1e-12);
        Assert.Equal(8.0 * (100 - 10 - 10 + 8) / (2.0 * 2), e.OddsRatio, 1e-12);
        Assert.Equal(ExactUpperTail(8, 10, 100, 10), e.PValue, 1e-15);
        Assert.Equal(["G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7"], e.Genes);

        var neutral = res.Single(r => r.Set == "neutral");
        Assert.Equal(1.0, neutral.PValue);
        var partial = res.Single(r => r.Set == "partial");
        Assert.Equal(ExactUpperTail(1, 5, 100, 10), partial.PValue, 1e-12);

        // BH over three tests.
        var ps = res.Select(r => r.PValue).ToArray();
        Assert.Equal(Math.Min(ps[0] * 3, Math.Min(ps[1] * 3 / 2, 1)), res[0].PAdj, 1e-14);
        Assert.Equal(Math.Min(ps[1] * 3 / 2, 1), res[1].PAdj, 1e-14);
        Assert.Equal(1.0, res[2].PAdj, 1e-14);
    }

    [Fact]
    public void Gmt_RoundTrips()
    {
        const string gmt = "SET_A\thttp://a\tTP53\tMDM2\t\tCDKN1A \n# comment\n\nSET_B\tdesc\tEGFR\n";
        var sets = Gmt.Read(new StringReader(gmt));
        Assert.Equal(2, sets.Count);
        Assert.Equal(3, sets[0].Count);
        Assert.Contains("CDKN1A", sets[0].Genes);
        Assert.Equal("http://a", sets[0].Description);
        var sw = new StringWriter();
        Gmt.Write(sw, sets);
        var again = Gmt.Read(new StringReader(sw.ToString()));
        Assert.Equal(sets[0].Genes.Order(), again[0].Genes.Order());
        Assert.Equal("SET_B", again[1].Name);
    }

    [Fact]
    public void EnrichmentScore_MatchesHandComputedRunningSum()
    {
        // Ranked scores 3, 2, 1, −1, −2, −3 for genes A..F; |r|^1 weights.
        double[] w = [3, 2, 1, 1, 2, 3];
        // {A, C}: N_R = 4. Running sum: +0.75, 0.5, 0.75, 0.5, 0.25, 0 → ES = 0.75 at A.
        Assert.Equal(0.75, Gsea.EnrichmentScore([0, 2], w, 6, out int peak), 1e-15);
        Assert.Equal(0, peak);
        // {A, E}: N_R = 5, miss step 1/4: 0.6, 0.35, 0.1, −0.15, 0.25, 0 → ES = 0.6.
        Assert.Equal(0.6, Gsea.EnrichmentScore([0, 4], w, 6, out _), 1e-15);
        // {D, F}: N_R = 4: −0.25, −0.5, −0.75, −0.5, −0.75, 0 → ES = −0.75, first reached before D.
        Assert.Equal(-0.75, Gsea.EnrichmentScore([3, 5], w, 6, out peak), 1e-15);
        Assert.Equal(3, peak);
        // Unweighted (p = 0) {A, E}: 0.5, 0.25, 0, −0.25, 0.25, 0 → 0.5.
        Assert.Equal(0.5, Gsea.EnrichmentScore([0, 4], [1, 1, 1, 1, 1, 1], 6, out _), 1e-15);
    }

    [Fact]
    public void Gsea_DetectsPlantedSetsAtBothEnds()
    {
        var r = new Random(41);
        var ranking = new Dictionary<string, double>();
        for (int i = 0; i < 1000; i++) ranking[$"G{i}"] = TestData.Normal(r);
        var up = Enumerable.Range(0, 30).Select(i => $"G{i}").ToArray();
        var down = Enumerable.Range(30, 30).Select(i => $"G{i}").ToArray();
        foreach (var g in up) ranking[g] = 3 + r.NextDouble();
        foreach (var g in down) ranking[g] = -3 - r.NextDouble();
        var randomSet = Enumerable.Range(0, 30).Select(i => $"G{100 + 7 * i}").ToArray();
        var sets = new[] { new GeneSet("UP", up), new GeneSet("DOWN", down), new GeneSet("RANDOM", randomSet) };
        var res = Gsea.Run(ranking, sets, permutations: 1000, seed: 1);
        var u = res.Single(x => x.Set == "UP");
        var d = res.Single(x => x.Set == "DOWN");
        var rnd = res.Single(x => x.Set == "RANDOM");
        Assert.True(u.EnrichmentScore > 0.9 && u.NormalizedEnrichmentScore > 1.5);
        Assert.True(d.EnrichmentScore < -0.9 && d.NormalizedEnrichmentScore < -1.5);
        Assert.True(u.PValue < 0.005 && d.PValue < 0.005);
        Assert.True(rnd.PValue > 0.01);
        Assert.Equal(30, u.LeadingEdge.Count);
        Assert.Equal(30, u.Size);
        Assert.InRange(u.FdrQ, 0, 0.05);
        // Deterministic for a seed.
        var again = Gsea.Run(ranking, sets, permutations: 1000, seed: 1);
        Assert.Equal(res.Select(x => x.PValue), again.Select(x => x.PValue));
    }

    [Fact]
    public void BenjaminiHochberg_MatchesRPAdjust()
    {
        // p.adjust(c(0.01, 0.02, 0.03, 0.04, 0.05), "BH") = rep(0.05, 5)
        Assert.All(DESeq.Bh([0.01, 0.02, 0.03, 0.04, 0.05]), q => Assert.Equal(0.05, q, 1e-15));

        double[] p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.34, 0.341,
            0.384, 0.569, 0.594, 0.696, 0.762, 0.94, 0.942, 0.975, 0.986];
        double[] expected = [0.025, 0.1, 0.21, 0.21, 0.21, 0.25, 0.26428571428571429, 0.49107142857142857, 0.49107142857142857,
            0.49107142857142857, 0.49107142857142857, 0.49107142857142857, 0.49107142857142857, 0.49107142857142857,
            0.5328125, 0.5328125, 0.56470588235294118, 0.78157894736842105, 0.78157894736842105, 0.87, 0.90714285714285714,
            0.986, 0.986, 0.986, 0.986];
        var q = DESeq.Bh(p);
        for (int i = 0; i < p.Length; i++) Assert.Equal(expected[i], q[i], 1e-12);

        // NaNs are ignored when counting tests, like p.adjust's n = sum(!is.na(p)).
        var withNa = DESeq.Bh([0.01, double.NaN, 0.02, double.NaN]);
        Assert.Equal(0.02, withNa[0], 1e-15);
        Assert.True(double.IsNaN(withNa[1]));
        Assert.Equal(0.02, withNa[2], 1e-15);
    }
}
