using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Core;
using Microsoft.DotNet.Interactive.Formatting;

namespace Biotech.Net.Notebooks.Tests;

public class AlignmentFormatterTests
{
    [Fact]
    public void PairwiseShowsMatchLineAndStats()
    {
        var aligner = new PairwiseAligner { MatchScore = 2, MismatchScore = -1, GapScore = -2 };
        var a = aligner.Align(Sequence.Dna("GATTACAGATTACA", "t"), Sequence.Dna("GATTCAGATTTACA", "q"));
        var html = AlignmentHtml.Pairwise(a);
        Assert.Contains("t × q", html);
        Assert.Contains($"score {a.Score}", html);
        Assert.Contains("|||", html);
        Assert.Contains("CIGAR", html);
        Assert.Contains("btn-A", html);
    }

    [Fact]
    public void MultipleAlignmentHasConsensus()
    {
        var msa = new MultipleSequenceAlignment(["s1", "s2", "s3"], ["ACGT-A", "ACGTTA", "ACCT-A"]);
        var html = AlignmentHtml.Multiple(msa);
        Assert.Contains("3 sequences", html);
        Assert.Contains("consensus", html);
        Assert.Contains("s2", html);
    }

    [Fact]
    public void TreeRendersSvgAndNewick()
    {
        var tree = PhylogeneticTree.ParseNewick("((A:0.1,B:0.2):0.3,(C:0.25,D:0.05):0.1);");
        var html = AlignmentHtml.Tree(tree);
        Assert.Contains("<svg", html);
        Assert.Contains(">A</text>", html);
        Assert.Contains("4 leaves", html);
        Assert.Contains("Newick", html);
    }

    [Fact]
    public void FormatterRegistryCoversAlignments()
    {
        BiotechFormatters.Register();
        var tree = PhylogeneticTree.ParseNewick("(A,B,(C,D));");
        Assert.Contains("phylogenetic tree", tree.ToDisplayString(HtmlFormatter.MimeType));
    }
}
