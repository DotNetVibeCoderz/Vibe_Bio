using Biotech.Net.Chemistry;
using Biotech.Net.Core;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Transcriptomics;
using Biotech.Net.Visualization;
using Microsoft.DotNet.Interactive.Formatting;

namespace Biotech.Net.Notebooks.Tests;

public class BiotechHtmlTests
{
    [Fact]
    public void SequenceIsColouredWrappedWithGc()
    {
        var s = Sequence.Dna(string.Concat(Enumerable.Repeat("ACGTGC", 30)), "seq1", "test sequence");
        var html = BiotechHtml.Sequence(s);
        Assert.Contains("seq1", html);
        Assert.Contains("180 bp", html);
        Assert.Contains("GC 66.7%", html);
        Assert.Contains("btn-A", html);
        Assert.Contains("btn-G", html);
        Assert.Contains("test sequence", html);
        Assert.Contains("\n 61  ", html.Replace("<span class=\"n\">", "\n")); // line 2 starts at residue 61
        Assert.Contains(BiotechInfo.CreditsId, html);
    }

    [Fact]
    public void LongSequenceIsTruncated()
    {
        var s = Sequence.Dna(new string('A', 5000), "long");
        var html = BiotechHtml.Sequence(s);
        Assert.Contains("2,000 more residues", html);
    }

    [Fact]
    public void ProteinUsesAminoAcidColours()
    {
        var p = new Protein("P1", Sequence.Protein("MKDEWLLR", "P1"), "Test protein") { Organism = "E. coli", GeneName = "tst" };
        var html = BiotechHtml.Protein(p);
        Assert.Contains("Test protein (P1)", html);
        Assert.Contains("8 aa", html);
        Assert.Contains("kDa", html);
        Assert.Contains("E. coli", html);
        Assert.Contains("btn-pos", html);
        Assert.Contains("btn-neg", html);
    }

    [Fact]
    public void ReadShowsQualityBars()
    {
        var r = SequenceRead.FromAscii("read1", "ACGTACGT", "IIII####");
        var html = BiotechHtml.Read(r);
        Assert.Contains("read1", html);
        Assert.Contains("<svg", html);
        Assert.Contains("#1E9E5A", html); // high quality
        Assert.Contains("#D1414A", html); // low quality
    }

    [Fact]
    public void VariantTableIsOneBased()
    {
        var v = new Variant("chr7", 117559590, "A", "G", "rs113993960") { Quality = 50 };
        v.Info["DP"] = "42";
        var html = BiotechHtml.Variant(v);
        Assert.Contains("117,559,590", html);
        Assert.Contains("rs113993960", html);
        Assert.Contains("transition", html);
        Assert.Contains("DP=42", html);
        Assert.Contains("PASS", html);
    }

    [Fact]
    public void IntervalShowsBothConventions()
    {
        var html = BiotechHtml.Interval(new GenomicInterval("chr1", 100, 200, Strand.Reverse));
        Assert.Contains("chr1:101-200(-)", html);
        Assert.Contains("[100, 200)", html);
        Assert.Contains("− strand", html);
    }

    [Fact]
    public void DatasetTableTruncates()
    {
        var ds = Dataset.Zeros("counts", Enumerable.Range(0, 50).Select(i => $"g{i}").ToArray(), Enumerable.Range(0, 20).Select(i => $"s{i}").ToArray());
        ds[0, 0] = 1234.5;
        ds.Annotate("condition", Enumerable.Range(0, 20).Select(i => i < 10 ? "ctl" : "trt").ToArray());
        var html = BiotechHtml.Dataset(ds);
        Assert.Contains("50 × 20", html);
        Assert.Contains("<table", html);
        Assert.Contains("1234.5", html);
        Assert.Contains("30 more rows", html);
        Assert.Contains("condition", html);
        Assert.DoesNotContain(">s15<", html);
    }

    [Fact]
    public void DifferentialExpressionTopTable()
    {
        var res = new DifferentialExpressionResult("DESeq2-Wald", "trt vs ctl", ["a", "b", "c"],
            [10, 20, 30], [2, -3, 0.1], [0.2, 0.3, 0.4], [10, -10, 0.25], [1e-8, 1e-6, 0.8], [3e-8, 1.5e-6, 0.8]);
        var html = BiotechHtml.DifferentialExpression(res);
        Assert.Contains("1 up", html);
        Assert.Contains("1 down", html);
        Assert.Contains("class=\"up\"", html);
        Assert.True(html.IndexOf(">a<", StringComparison.Ordinal) < html.IndexOf(">c<", StringComparison.Ordinal));
    }

    [Fact]
    public void MoleculeHasSvgDepictionAndDescriptors()
    {
        var html = BiotechHtml.Molecule(Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin"));
        Assert.Contains("aspirin", html);
        Assert.Contains("<svg", html);
        Assert.DoesNotContain("<?xml", html);
        Assert.Contains("C9H8O4", html);
        Assert.Contains("SMILES", html);
        Assert.Contains("TPSA", html);
    }

    [Fact]
    public void StructureEmbeds3DViewer()
    {
        var html = BiotechHtml.Structure(SampleStructures.Crambin());
        Assert.Contains("1CRN", html);
        Assert.Contains("3dmol@" + ProteinViewer.ThreeDmolVersion, html);
        Assert.Contains("TTCCPSIVARSNFNVCRLPGTPEA", html);
    }

    [Fact]
    public void PlotsEmbedInEveryMode()
    {
        var logo = new SequenceLogo(["ACGT", "ACGA", "ACTT"]);
        Assert.Contains("<svg", BiotechHtml.Plot(logo, PlotEmbedding.Svg));
        Assert.Contains("data:image/png;base64,", BiotechHtml.Plot(logo, PlotEmbedding.Png));
        Assert.Contains("<svg", BiotechHtml.Plot(logo)); // no interactive form: falls back to SVG

        var volcano = new VolcanoPlot(["g1", "g2", "g3"], [2, -2, 0], [1e-6, 1e-5, 0.5]);
        var html = BiotechHtml.Plot(volcano);
        Assert.Contains("plotly", html);
        Assert.Contains(BiotechInfo.CreditsId, html);
    }

    [Fact]
    public void GeneListsTranscripts()
    {
        var g = new Gene("ENSG1", new GenomicInterval("chr1", 0, 1000, Strand.Forward), "TP");
        g.Transcripts.Add(new Transcript("ENST1", "ENSG1", [new("chr1", 0, 200, Strand.Forward), new("chr1", 600, 1000, Strand.Forward)]));
        var html = BiotechHtml.Gene(g);
        Assert.Contains("ENST1", html);
        Assert.Contains("600", html);
    }

    [Fact]
    public void HtmlIsEscaped()
    {
        var html = BiotechHtml.Sequence(Sequence.Dna("ACGT", "<script>x</script>"));
        Assert.DoesNotContain("<script>x", html);
        Assert.Contains("&lt;script&gt;", html);
    }
}

[Collection("formatter registry")]
public class FormatterRegistrationTests
{
    [Fact]
    public void RegisteredFormattersProduceHtml()
    {
        BiotechFormatters.Register();
        BiotechFormatters.Register(); // idempotent
        Assert.True(BiotechFormatters.IsRegistered);

        string seq = Sequence.Dna("ACGTACGT", "fmt").ToDisplayString(HtmlFormatter.MimeType);
        Assert.Contains("btn-card", seq);
        Assert.Contains("fmt", seq);

        // A derived plot type uses the BioPlot formatter.
        string plot = new SequenceLogo(["ACGT", "ACGA"]).ToDisplayString(HtmlFormatter.MimeType);
        Assert.Contains("biotech-plot-output", plot);

        string iv = new GenomicInterval("chr2", 0, 10).ToDisplayString(HtmlFormatter.MimeType);
        Assert.Contains("chr2:1-10", iv);

        string mol = Molecule.FromSmiles("CCO", "ethanol").ToDisplayString(HtmlFormatter.MimeType);
        Assert.Contains("ethanol", mol);

        Assert.Contains(HtmlFormatter.MimeType, Formatter.GetPreferredMimeTypesFor(typeof(Sequence)));
        Assert.Contains("notebook formatters loaded", BiotechFormatters.BannerHtml());
    }
}
