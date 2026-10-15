using Biotech.Net.Alignments.Phylogenetics;

namespace Biotech.Net.Alignments.Tests;

public class PhylogeneticsTests
{
    private static DistanceMatrix Matrix(string names, double[,] d) => new(names.Select(c => c.ToString()), d);

    [Fact]
    public void NeighborJoining_WikipediaWorkedExample()
    {
        // https://en.wikipedia.org/wiki/Neighbor_joining#Example: a–u 2, b–u 3, u–v 3, c–v 4, v–w 2, d–w 2, e–w 1.
        var dm = Matrix("abcde", new double[,]
        {
            { 0, 5, 9, 9, 8 },
            { 5, 0, 10, 10, 9 },
            { 9, 10, 0, 8, 7 },
            { 9, 10, 8, 0, 3 },
            { 8, 9, 7, 3, 0 },
        });
        var tree = TreeBuilder.NeighborJoining(dm);
        Assert.Equal(2, tree.FindLeaf("a")!.BranchLength!.Value, 10);
        Assert.Equal(3, tree.FindLeaf("b")!.BranchLength!.Value, 10);
        Assert.Equal(4, tree.FindLeaf("c")!.BranchLength!.Value, 10);
        Assert.Equal(2, tree.FindLeaf("d")!.BranchLength!.Value, 10);
        Assert.Equal(1, tree.FindLeaf("e")!.BranchLength!.Value, 10);
        Assert.Equal(17, tree.TotalBranchLength, 10); // 2+3+3+4+2+2+1
        // The matrix is additive, so NJ reproduces every pairwise distance exactly.
        foreach (string x in dm.Names)
            foreach (string y in dm.Names)
                Assert.Equal(dm[x, y], tree.Distance(x, y), 10);
        Assert.Equal(0, tree.RobinsonFoulds(PhylogeneticTree.ParseNewick("((a,b),c,(d,e));")));
    }

    [Fact]
    public void Upgma_WikipediaWorkedExample()
    {
        // https://en.wikipedia.org/wiki/UPGMA#Working_example (5S rRNA of Bacillus subtilis, B. stearothermophilus,
        // Lactobacillus viridescens, Acholeplasma modicum, Micrococcus luteus): (((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);
        var dm = Matrix("abcde", new double[,]
        {
            { 0, 17, 21, 31, 23 },
            { 17, 0, 30, 34, 21 },
            { 21, 30, 0, 28, 39 },
            { 31, 34, 28, 0, 43 },
            { 23, 21, 39, 43, 0 },
        });
        var tree = TreeBuilder.Upgma(dm);
        var expected = PhylogeneticTree.ParseNewick("(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);");
        Assert.Equal(0, tree.RobinsonFoulds(expected));
        Assert.Equal(8.5, tree.FindLeaf("a")!.BranchLength!.Value, 10);
        Assert.Equal(11, tree.FindLeaf("e")!.BranchLength!.Value, 10);
        Assert.Equal(14, tree.FindLeaf("c")!.BranchLength!.Value, 10);
        Assert.Equal(2.5, tree.FindLeaf("a")!.Parent!.BranchLength!.Value, 10);
        Assert.Equal(5.5, tree.FindLeaf("e")!.Parent!.BranchLength!.Value, 10);
        Assert.Equal(2.5, tree.FindLeaf("c")!.Parent!.BranchLength!.Value, 10);
        // Ultrametric: every leaf at height 16.5.
        foreach (var leaf in tree.Leaves) Assert.Equal(16.5, PhylogeneticTree.Distance(leaf, tree.Root), 10);
    }

    [Theory]
    [InlineData("((A:0.1,B:0.2):0.05,C:0.3);")]
    [InlineData("(A,B,(C,D)E)F;")]
    [InlineData("((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);")]
    [InlineData("('quoted name':1,'it''s':2.5e-3);")]
    public void Newick_RoundTrips(string newick)
    {
        var tree = PhylogeneticTree.ParseNewick(newick);
        string written = tree.ToNewick();
        var again = PhylogeneticTree.ParseNewick(written);
        Assert.Equal(written, again.ToNewick());
        Assert.Equal(tree.LeafNames, again.LeafNames);
    }

    [Fact]
    public void Newick_ParsesLabelsLengthsAndComments()
    {
        var tree = PhylogeneticTree.ParseNewick("((A:1,B:2)95:0.5,[a comment]C_x:6);");
        Assert.Equal(["A", "B", "C x"], tree.LeafNames);
        Assert.Equal("95", tree.FindLeaf("A")!.Parent!.Name);
        Assert.Equal(9.5, tree.TotalBranchLength, 12);
        Assert.Equal("((A:1,B:2)95:0.5,'C x':6);", tree.ToNewick());
        Assert.Equal("((A,B)95,'C x');", tree.ToNewick(branchLengths: false));
        Assert.Throws<FormatException>(() => PhylogeneticTree.ParseNewick("((A,B);"));
    }

    [Fact]
    public void RootAtMidpoint_SplitsLongestPath()
    {
        var tree = PhylogeneticTree.ParseNewick("((A:1,B:2):0.5,C:6);");
        // Longest path B–C = 2 + 0.5 + 6 = 8.5; midpoint lies 4.25 from C.
        tree.RootAtMidpoint();
        Assert.Equal(2, tree.Root.Children.Count);
        Assert.Equal(4.25, PhylogeneticTree.Distance(tree.Root, tree.FindLeaf("C")!), 12);
        Assert.Equal(4.25, PhylogeneticTree.Distance(tree.Root, tree.FindLeaf("B")!), 12);
        Assert.Equal(8.5, tree.Distance("B", "C"), 12);
        Assert.Equal(3, tree.Distance("A", "B"), 12);
        Assert.Equal(9.5, tree.TotalBranchLength, 12);
        Assert.Equal(3, tree.Leaves.Count);
    }

    [Fact]
    public void RootAtMidpoint_OnNjTree_PreservesDistances()
    {
        var dm = Matrix("abcde", new double[,]
        {
            { 0, 5, 9, 9, 8 }, { 5, 0, 10, 10, 9 }, { 9, 10, 0, 8, 7 }, { 9, 10, 8, 0, 3 }, { 8, 9, 7, 3, 0 },
        });
        var tree = TreeBuilder.NeighborJoining(dm).RootAtMidpoint();
        foreach (string x in dm.Names)
            foreach (string y in dm.Names)
                Assert.Equal(dm[x, y], tree.Distance(x, y), 10);
        // Longest path is b–c (10): the root sits 5 from both.
        Assert.Equal(5, PhylogeneticTree.Distance(tree.Root, tree.FindLeaf("b")!), 10);
        Assert.Equal(5, PhylogeneticTree.Distance(tree.Root, tree.FindLeaf("c")!), 10);
    }

    [Fact]
    public void RobinsonFoulds_KnownValues()
    {
        var t1 = PhylogeneticTree.ParseNewick("((A,B),(C,D));");
        var t2 = PhylogeneticTree.ParseNewick("((A,C),(B,D));");
        var t1Rerooted = PhylogeneticTree.ParseNewick("(A,(B,(C,D)));");
        Assert.Equal(2, t1.RobinsonFoulds(t2));
        Assert.Equal(1.0, t1.RobinsonFoulds(t2, normalized: true));
        Assert.Equal(0, t1.RobinsonFoulds(t1Rerooted));
        var t5a = PhylogeneticTree.ParseNewick("(((A,B),C),(D,E));");
        var t5b = PhylogeneticTree.ParseNewick("(((A,C),B),(D,E));");
        Assert.Equal(2, t5a.RobinsonFoulds(t5b)); // {A,B}|{C,D,E} vs {A,C}|{B,D,E}
        Assert.Throws<ArgumentException>(() => t1.RobinsonFoulds(t5a));
    }

    [Fact]
    public void Distances_PDistance_JukesCantor_Kimura()
    {
        // 8 sites: 1 transition (A↔G) and 1 transversion (C↔A).
        string a = "ACGTACGT", b = "GCGTAAGT";
        Assert.Equal(0.25, DistanceMatrix.PairDistance(a, b, DistanceModel.PDistance, true), 12);
        Assert.Equal(-0.75 * Math.Log(1 - 4.0 / 3 * 0.25), DistanceMatrix.PairDistance(a, b, DistanceModel.JukesCantor, true), 12);
        double P = 0.125, Q = 0.125;
        Assert.Equal(-0.5 * Math.Log(1 - 2 * P - Q) - 0.25 * Math.Log(1 - 2 * Q), DistanceMatrix.PairDistance(a, b, DistanceModel.Kimura2P, true), 12);
        // Gap columns are skipped pairwise.
        Assert.Equal(0.0, DistanceMatrix.PairDistance("AC-T", "ACGT", DistanceModel.PDistance, true));
        // Saturation.
        Assert.True(double.IsPositiveInfinity(DistanceMatrix.PairDistance("AAAA", "CGTC", DistanceModel.JukesCantor, true)));
        // Protein JC (Poisson-corrected, b = 19/20).
        Assert.Equal(-0.95 * Math.Log(1 - 0.5 / 0.95), DistanceMatrix.PairDistance("MKVL", "MRVA", DistanceModel.JukesCantor, false), 12);
    }

    [Fact]
    public void DistanceMatrix_FromAlignment_AndPhylip()
    {
        var msa = new MultipleSequenceAlignment(["s1", "s2", "s3"], ["ACGTACGT", "ACGTACGA", "TCGTACGA"]);
        var dm = DistanceMatrix.FromAlignment(msa);
        Assert.Equal(0.125, dm["s1", "s2"], 12);
        Assert.Equal(0.25, dm["s1", "s3"], 12);
        Assert.Equal(dm[0, 2], dm[2, 0]);
        Assert.StartsWith("3", dm.ToPhylip());
        Assert.Throws<ArgumentException>(() => DistanceMatrix.FromAlignment(new MultipleSequenceAlignment(["p", "q"], ["MKLV", "MRLV"]), DistanceModel.Kimura2P));
    }
}
