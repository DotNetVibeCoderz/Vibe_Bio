using Biotech.Net.Core;
using Biotech.Net.Genomics.Assembly;

namespace Biotech.Net.Genomics.Tests;

public class AssemblyTests
{
    private static bool SameOrReverse(string a, string b) => a == b || a == Simulation.ReverseComplement(b);

    private static bool Within(string contig, string genome) => genome.Contains(contig) || genome.Contains(Simulation.ReverseComplement(contig));

    [Fact]
    public void AssemblesErrorFreeReadsIntoTheGenome()
    {
        var rnd = new Random(11);
        string genome = Simulation.RandomDna(rnd, 3000);
        var reads = Simulation.Reads(rnd, genome, 600, 100, 0); // 20× coverage
        var result = DeBruijnAssembler.Assemble(reads.Select(r => r.Read), new AssemblyOptions { K = 31, MinKmerCount = 2 });
        var contig = Assert.Single(result.Contigs);
        // Terminal k-mers covered by fewer than MinKmerCount reads are dropped, so a few end bases may be missing.
        Assert.True(Within(contig.ToString(), genome) && contig.Length >= 2990, $"contig of {contig.Length} bp");
        Assert.Equal(contig.Length, result.Stats.N50);
        Assert.Equal(1, result.Stats.L50);
        Assert.Equal(contig.Length - 31 + 1, result.SolidKmers);
    }

    [Fact]
    public void SolidKmerThresholdRemovesSequencingErrors()
    {
        var rnd = new Random(12);
        string genome = Simulation.RandomDna(rnd, 4000);
        var reads = Simulation.Reads(rnd, genome, 1200, 100, 0.002); // ~0.2 errors per read, 30×
        var noisy = DeBruijnGraph.Build(reads.Select(r => r.Read.Sequence), 31, minCount: 1);
        var solid = DeBruijnGraph.Build(reads.Select(r => r.Read.Sequence), 31, minCount: 3);
        Assert.True(noisy.KmerCount > solid.KmerCount);
        var result = DeBruijnAssembler.Assemble(reads.Select(r => r.Read.Sequence), new AssemblyOptions { MinKmerCount = 3 });
        Assert.True(result.Stats.Longest >= 3900, result.Stats.ToString());
        Assert.Contains(result.Contigs, c => genome.Contains(c.ToString()) || genome.Contains(Simulation.ReverseComplement(c.ToString())));
    }

    [Fact]
    public void TipClippingAndBubblePoppingCleanTheGraph()
    {
        var rnd = new Random(13);
        string genome = Simulation.RandomDna(rnd, 2000);
        var reads = Simulation.Reads(rnd, genome, 400, 100, 0).Select(r => r.Read.Sequence).ToList();
        // An error in the middle of one read makes a bubble; an error at a read end makes a tip.
        string mid = genome.Substring(500, 100);
        mid = mid[..50] + (mid[50] == 'A' ? 'C' : 'A') + mid[51..];
        string end = genome.Substring(1200, 100);
        end = end[..97] + (end[97] == 'A' ? 'C' : 'A') + end[98..];
        reads.Add(Sequence.Dna(mid, "bubble"));
        reads.Add(Sequence.Dna(end, "tip"));

        var raw = DeBruijnAssembler.Assemble(reads, new AssemblyOptions { MinKmerCount = 1, RemoveTips = false, PopBubbles = false, MinContigLength = 1 });
        Assert.True(raw.Contigs.Count > 1);

        var clean = DeBruijnAssembler.Assemble(reads, new AssemblyOptions { MinKmerCount = 1 });
        Assert.True(clean.TipsRemoved >= 1, $"tips {clean.TipsRemoved}");
        Assert.True(clean.BubblesPopped >= 1, $"bubbles {clean.BubblesPopped}");
        var contig = Assert.Single(clean.Contigs);
        Assert.True(Within(contig.ToString(), genome) && contig.Length >= 1990, $"contig of {contig.Length} bp, tips {clean.TipsRemoved}, bubbles {clean.BubblesPopped}");
    }

    [Fact]
    public void UnitigsOfALinearGraph()
    {
        var g = DeBruijnGraph.Build([Sequence.Dna("CTGGGTAGGATACGGCGGA")], 5, minCount: 1);
        var u = Assert.Single(g.Unitigs());
        Assert.True(SameOrReverse(u.Sequence.ToString(), "CTGGGTAGGATACGGCGGA"));
        Assert.Equal(15, u.KmerCount);
        Assert.Equal(1, g.Count("CTGGG"));
        Assert.Equal(1, g.Count("CCCAG")); // reverse complement of CTGGG
        Assert.Equal(0, g.Count("GGGGG"));
    }

    [Fact]
    public void AssemblyStatistics()
    {
        var s = AssemblyStats.Compute([100L, 200, 300, 400, 500]);
        Assert.Equal(5, s.Count);
        Assert.Equal(1500, s.TotalLength);
        Assert.Equal(400, s.N50);
        Assert.Equal(2, s.L50);
        Assert.Equal(200, s.N90);
        Assert.Equal(4, s.L90);
        Assert.Equal(300, s.MeanLength);
        var withSeq = AssemblyStats.Compute([Sequence.Dna("GGCC"), Sequence.Dna("ATAT")]);
        Assert.Equal(0.5, withSeq.GcContent, 12);
    }
}
