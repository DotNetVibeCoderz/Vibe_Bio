using Biotech.Net.Chemistry;
using Biotech.Net.Core;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Transcriptomics;

namespace Biotech.Net.Visualization.Tests;

public class GenomeBrowserTests
{
    private static GenomeBrowser Build()
    {
        var seq = Sequence.Dna(string.Concat(Enumerable.Repeat("ACGTTGCA", 1000)), "chr1");
        var feats = new List<SequenceFeature>
        {
            new("CDS", FeatureLocation.Span(1300, 2200, Strand.Forward)),
            new("promoter", FeatureLocation.Span(1050, 1200, Strand.Forward)),
            new("misc_feature", new FeatureLocation([(4100, 4300), (4500, 4700)], Strand.Reverse)),
        };
        feats[0].With("gene", "alphaCDS");
        return new GenomeBrowser(TestData.Region) { Title = "Test region" }
            .AddRuler()
            .AddGenes(TestData.Genes())
            .AddFeatures(feats)
            .AddCoverage(TestData.Coverage(), 1000)
            .AddVariants(TestData.Variants())
            .AddSequence(seq);
    }

    [Fact]
    public void RendersPngAndSvg()
    {
        var b = Build();
        var (w, h) = ImageAssert.Png(b.ToPng());
        Assert.Equal(b.Width, w);
        Assert.Equal(b.Height, h);
        var svg = b.ToSvg();
        ImageAssert.Svg(svg);
        Assert.Contains("ALPHA", svg);
        Assert.Contains("rs1", svg);
        ImageAssert.Dump(b, "genome-browser");
    }

    [Fact]
    public void LayoutIsDeterministicAndStacked()
    {
        var b = Build();
        var l1 = b.ComputeLayout();
        var l2 = b.ComputeLayout();
        Assert.Equal(6, l1.Count);
        Assert.Equal(l1.Select(p => p.Row), l2.Select(p => p.Row));
        for (int i = 1; i < l1.Count; i++) Assert.True(l1[i].Row.Y >= l1[i - 1].Row.Bottom);
        Assert.All(l1, p => Assert.Equal(b.Padding + b.NameGutter, p.DataArea.X));
        Assert.True(b.Height >= l1[^1].Row.Bottom);
    }

    [Fact]
    public void GeneRowsPackWithoutOverlap()
    {
        var b = new GenomeBrowser(TestData.Region);
        var track = new GeneTrack(TestData.Genes());
        var rows = track.PackRows(b, 800);
        var models = track.Models();
        Assert.Equal(4, models.Count); // two isoforms of ALPHA + BETA + GAMMA
        Assert.NotEqual(rows[0], rows[1]); // overlapping isoforms go to different rows
        for (int i = 0; i < models.Count; i++)
            for (int j = i + 1; j < models.Count; j++)
                if (rows[i] == rows[j]) Assert.False(models[i].Span.Overlaps(models[j].Span));
    }

    [Fact]
    public void ToXMapsRegionEdges()
    {
        var b = new GenomeBrowser(TestData.Region);
        var area = new Matplotlib.Net.Rendering.Rect(100, 0, 500, 10);
        Assert.Equal(100, b.ToX(1000, area), 6);
        Assert.Equal(600, b.ToX(6000, area), 6);
        Assert.Equal(350, b.ToX(3500, area), 6);
    }

    [Fact]
    public void CoverageBinsTakeMaximum()
    {
        var t = new CoverageTrack([1, 5, 2, 8, 3, 1], 10);
        var bins = t.Bin(new GenomicInterval("c", 10, 16), 3);
        Assert.Equal([5, 8, 3], bins);
    }

    [Fact]
    public void DarkThemeRendersDarkBackground()
    {
        var b = Build();
        b.Theme = BioTheme.Dark;
        using var bmp = SkiaSharp.SKBitmap.Decode(b.ToPng());
        var px = bmp.GetPixel(2, 2);
        Assert.True(px.Red < 40 && px.Green < 40 && px.Blue < 50);
    }

    [Fact]
    public void ScaledPngDoublesSize()
    {
        var b = Build();
        var (w, h) = ImageAssert.Png(b.ToPng(2));
        Assert.Equal(b.Width * 2, w);
        Assert.Equal(b.Height * 2, h);
    }
}

public class ClusteringTests
{
    [Fact]
    public void ClustersTwoObviousGroups()
    {
        double[][] obs = [[0, 0], [0.1, 0], [0, 0.2], [10, 10], [10.2, 10], [10, 10.1]];
        foreach (var link in Enum.GetValues<LinkageMethod>())
        {
            var d = HierarchicalClustering.Cluster(obs, DistanceMetric.Euclidean, link);
            Assert.Equal(5, d.Merges.Count);
            Assert.True(d.Merges.Zip(d.Merges.Skip(1)).All(p => p.First.Height <= p.Second.Height));
            Assert.Equal(6, d.Merges[^1].Size);
            var labels = d.Cut(2);
            Assert.Equal(labels[0], labels[1]);
            Assert.Equal(labels[0], labels[2]);
            Assert.Equal(labels[3], labels[5]);
            Assert.NotEqual(labels[0], labels[3]);
            Assert.Equal(Enumerable.Range(0, 6), d.Order.Order());
        }
    }

    [Fact]
    public void MatchesKnownSingleLinkageHeights()
    {
        // 1-D points 0, 1, 3, 7: single linkage merges at 1, 2, 4.
        double[][] obs = [[0], [1], [3], [7]];
        var d = HierarchicalClustering.Cluster(obs, DistanceMetric.Euclidean, LinkageMethod.Single);
        Assert.Equal([1.0, 2.0, 4.0], d.Merges.Select(m => m.Height));
        var avg = HierarchicalClustering.Cluster(obs, DistanceMetric.Euclidean, LinkageMethod.Average);
        Assert.Equal(1.0, avg.Merges[0].Height, 9);
        Assert.Equal(2.5, avg.Merges[1].Height, 9);   // (3 + 2) / 2
        Assert.Equal(17.0 / 3, avg.Merges[2].Height, 9); // (7 + 6 + 4) / 3
    }

    [Fact]
    public void SegmentsCoverEveryMerge()
    {
        double[][] obs = [[0], [1], [3], [7], [8]];
        var d = HierarchicalClustering.Cluster(obs);
        Assert.Equal(d.Merges.Count * 3, d.Segments().Count);
    }
}

public class HeatmapTests
{
    [Fact]
    public void RowZScoresHaveZeroMeanUnitSd()
    {
        var z = Heatmap.Standardise(TestData.Expression(), HeatmapScale.Rows);
        for (int i = 0; i < z.GetLength(0); i++)
        {
            var row = Enumerable.Range(0, z.GetLength(1)).Select(j => z[i, j]).ToArray();
            Assert.Equal(0, row.Average(), 9);
            Assert.Equal(1, Math.Sqrt(row.Select(v => v * v).Average()), 9);
        }
    }

    [Fact]
    public void ClusteredLayoutGroupsBlocks()
    {
        var h = new Heatmap(TestData.Expression(), columnLabels: Enumerable.Range(1, 8).Select(i => $"S{i}").ToArray())
        {
            ColumnGroups = ["ctl", "ctl", "ctl", "ctl", "trt", "trt", "trt", "trt"],
            Title = "Expression",
        };
        var layout = h.ComputeLayout();
        Assert.NotNull(layout.RowTree);
        Assert.NotNull(layout.ColumnDendrogram);
        Assert.NotNull(layout.GroupStrip);
        // The first four columns stay together.
        var firstHalf = layout.ColumnOrder.Take(4).Order().ToArray();
        Assert.True(firstHalf.SequenceEqual([0, 1, 2, 3]) || firstHalf.SequenceEqual([4, 5, 6, 7]));
        Assert.Equal(layout, layout with { }); // record value semantics
        Assert.Equal(layout.RowOrder, h.ComputeLayout().RowOrder);

        ImageAssert.Png(h.ToPng());
        var svg = h.ToSvg();
        ImageAssert.Svg(svg);
        Assert.Contains("z-score", svg);
        ImageAssert.Dump(h, "heatmap");
    }

    [Fact]
    public void LargeMatrixUsesImagePath()
    {
        var h = new Heatmap(TestData.Expression(400, 120)) { ClusterColumns = false };
        var (w, hh) = ImageAssert.Png(h.ToPng());
        Assert.Equal(h.Width, w);
        Assert.Equal(h.Height, hh);
    }

    [Fact]
    public void PlotlyHeatmapContainsClusteredLabels()
    {
        var ds = new Dataset("expr", Enumerable.Range(0, 30).Select(i => $"g{i}").ToArray(), Enumerable.Range(0, 8).Select(i => $"s{i}").ToArray(),
            Flatten(TestData.Expression()));
        ds.Annotate("condition", "a", "a", "a", "a", "b", "b", "b", "b");
        var h = Heatmap.FromDataset(ds);
        Assert.NotNull(h.ColumnGroups);
        var html = h.ToHtml();
        Assert.Contains("heatmap", html);
        Assert.Contains("g0", html);
    }

    private static double[] Flatten(double[,] m)
    {
        var v = new double[m.Length];
        int k = 0;
        foreach (double x in m) v[k++] = x;
        return v;
    }
}

public class ExpressionPlotTests
{
    [Fact]
    public void VolcanoClassifiesAndRenders()
    {
        var res = TestData.DeResult();
        var v = new VolcanoPlot(res) { LabelTop = 5, Highlight = ["gene100"] };
        var (up, down, ns) = v.Counts();
        Assert.Equal(res.Count, up + down + ns);
        Assert.True(up >= 15 && down >= 10);
        var calls = v.Classify();
        Assert.All(Enumerable.Range(0, res.Count).Where(i => calls[i] == Regulation.Up), i => Assert.True(res.Log2FoldChange[i] >= 1 && res.PAdj[i] < 0.05));
        Assert.Contains(100, v.LabelledGenes());

        var (w, h) = ImageAssert.Png(v.ToPng());
        Assert.Equal((v.Width, v.Height), (w, h));
        var svg = v.ToSvg();
        ImageAssert.Svg(svg);
        Assert.Contains("log2 fold change", svg);
        var html = v.ToHtml();
        Assert.Contains("scatter", html);
        Assert.Contains("gene0", html);
        ImageAssert.Dump(v, "volcano");
    }

    [Fact]
    public void MaPlotRenders()
    {
        var ma = new MaPlot(TestData.DeResult());
        Assert.False(ma.MeanIsLog);
        ImageAssert.Png(ma.ToPng());
        ImageAssert.Svg(ma.ToSvg());
        ImageAssert.Dump(ma, "ma");
    }

    [Fact]
    public void BenjaminiHochbergMatchesReference()
    {
        // R: p.adjust(c(0.01, 0.04, 0.03, 0.2), "BH") = 0.04 0.0533 0.0533 0.2
        var q = DifferentialPlot.BenjaminiHochberg([0.01, 0.04, 0.03, 0.2]);
        Assert.Equal(0.04, q[0], 9);
        Assert.Equal(0.04 * 4 / 3, q[1], 9);
        Assert.Equal(0.04 * 4 / 3, q[2], 9);
        Assert.Equal(0.2, q[3], 9);
    }

    [Fact]
    public void PcaPlotLabelsVariance()
    {
        var rng = new Random(1);
        var m = new DenseMatrix(12, 50);
        for (int i = 0; i < 12; i++) for (int j = 0; j < 50; j++) m[i, j] = rng.NextDouble() + (i < 6 && j < 10 ? 4 : 0);
        var pca = Pca.Fit(m, 3);
        var groups = Enumerable.Range(0, 12).Select(i => i < 6 ? "treated" : "control").ToArray();
        var p = EmbeddingPlot.FromPca(pca, groups, Enumerable.Range(1, 12).Select(i => $"S{i}").ToArray());
        Assert.StartsWith("PC1 (", p.XLabel);
        Assert.Equal(["treated", "control"], p.GroupLevels());
        ImageAssert.Png(p.ToPng());
        var svg = p.ToSvg();
        ImageAssert.Svg(svg);
        Assert.Contains("treated", svg);
        Assert.Contains("PC2", p.ToHtml());
        ImageAssert.Dump(p, "pca");
    }

    [Fact]
    public void UmapEmbeddingWithContinuousColour()
    {
        var rng = new Random(2);
        var e = new DenseMatrix(300, 2);
        for (int i = 0; i < 300; i++) { e[i, 0] = rng.NextDouble() * 10 + (i % 3) * 12; e[i, 1] = rng.NextDouble() * 10; }
        var p = EmbeddingPlot.FromEmbedding(e, Enumerable.Range(0, 300).Select(i => $"cluster {i % 3}").ToArray());
        Assert.Equal("UMAP 1", p.XLabel);
        ImageAssert.Png(p.ToPng());
        p.Values = Enumerable.Range(0, 300).Select(i => (double)i).ToArray();
        ImageAssert.Png(p.ToPng());
        ImageAssert.Dump(p, "umap");
    }
}

public class SequencePlotTests
{
    private static readonly string[] Sites = ["TATAAT", "TATAAT", "TATGAT", "TACAAT", "TATAAT", "GATAAT", "TATACT", "TATAAT"];

    [Fact]
    public void LogoColumnsMatchInformationContent()
    {
        var logo = new SequenceLogo(Sites);
        Assert.True(logo.IsNucleotide);
        var cols = logo.ComputeColumns();
        Assert.Equal(6, cols.Count);
        Assert.Equal(2.0, cols[1].Information, 9); // A in every sequence
        Assert.Equal('A', cols[1].Letters.Single().Symbol);
        Assert.True(cols[0].Information < 2 && cols[0].Information > 1);
        Assert.Equal('T', cols[0].Letters[^1].Symbol); // the most frequent letter is on top
        Assert.All(cols, c => Assert.Equal(c.Information, c.Letters.Sum(l => l.Height), 9));

        logo.Mode = LogoMode.Probability;
        Assert.All(logo.ComputeColumns(), c => Assert.Equal(1, c.Letters.Sum(l => l.Height), 9));
        logo.Mode = LogoMode.Bits;
        logo.SmallSampleCorrection = true;
        Assert.True(logo.ComputeColumns()[1].Information < 2);
    }

    [Fact]
    public void LogoRendersGlyphs()
    {
        var logo = new SequenceLogo(Sites) { Title = "-10 box" };
        Assert.NotEmpty(LetterGlyphs.Outline('A'));
        ImageAssert.Png(logo.ToPng());
        var svg = logo.ToSvg();
        ImageAssert.Svg(svg);
        Assert.Contains("<path", svg);
        ImageAssert.Dump(logo, "logo");
    }

    [Fact]
    public void EmbeddedLetterOutlinesCoverTheUnitSquare()
    {
        for (char c = 'A'; c <= 'Z'; c++)
        {
            var outline = LetterGlyphs.Outline(c);
            Assert.NotEmpty(outline);
            var pts = outline.SelectMany(s => s.Points).ToArray();
            Assert.All(pts, p => Assert.InRange(p.X, -1e-9, 1 + 1e-9));
            Assert.All(pts, p => Assert.InRange(p.Y, -1e-9, 1 + 1e-9));
            Assert.Equal(1, pts.Max(p => p.X) - pts.Min(p => p.X), 2);
            Assert.Equal(1, pts.Max(p => p.Y) - pts.Min(p => p.Y), 2);
        }
        Assert.Equal(LetterGlyphs.Outline('A').Length, LetterGlyphs.Outline('a').Length);
    }

    [Fact]
    public void ProteinLogoFromSequences()
    {
        var seqs = new[] { "MKVLA", "MKILA", "MRVLG", "MKVLA" }.Select((s, i) => Sequence.Protein(s, $"p{i}"));
        var logo = new SequenceLogo(seqs);
        Assert.False(logo.IsNucleotide);
        Assert.Equal(Math.Log2(20), logo.MaxBits, 9);
        ImageAssert.Png(logo.ToPng());
    }

    [Fact]
    public void DotPlotFindsDiagonalAndInversion()
    {
        var a = Sequence.Dna("ACGTTGCATGCCATAGGCTAGCTAGGATCCGATCGATTACGGA", "a");
        var self = new DotPlot(a, a, window: 8);
        var m = self.FindMatches();
        Assert.Contains(m, d => d.X == 0 && d.Y == 0 && !d.ReverseComplement);
        Assert.True(m.Count(d => !d.ReverseComplement && d.X == d.Y) == a.Length - 8 + 1);

        var rc = a.ReverseComplement().WithId("rc");
        var inv = new DotPlot(a, rc, window: 8);
        // The inversion lies on the anti-diagonal: window i of a matches window (n - i - w) of its reverse complement.
        Assert.Contains(inv.FindMatches(), d => d.ReverseComplement && d.X == 0 && d.Y == a.Length - 8);

        ImageAssert.Png(self.ToPng());
        ImageAssert.Svg(self.ToSvg());
        ImageAssert.Dump(inv, "dotplot");
    }
}

public class StructurePlotTests
{
    [Fact]
    public void CaTraceRendersCrambin()
    {
        var s = SampleStructures.Crambin();
        var p = ProteinViewer.CaTrace(s);
        Assert.Equal(46, p.ResidueCount);
        var pts = p.Project();
        Assert.Equal(46, pts.Length);
        Assert.Equal(pts, p.Project()); // deterministic
        ImageAssert.Png(p.ToPng());
        ImageAssert.Svg(p.ToSvg());
        foreach (var c in Enum.GetValues<ResidueColoring>())
        {
            p.Coloring = c;
            ImageAssert.Png(p.ToPng());
        }
        p.Coloring = ResidueColoring.SecondaryStructure;
        ImageAssert.Dump(p, "catrace");
    }

    [Fact]
    public void ContactMapIsSymmetric()
    {
        var p = ProteinViewer.ContactMap(SampleStructures.Crambin());
        var m = p.Matrix();
        int n = m.GetLength(0);
        Assert.Equal(46, n);
        for (int i = 0; i < n; i++) { Assert.Equal(1, m[i, i]); for (int j = 0; j < n; j++) Assert.Equal(m[i, j], m[j, i]); }
        ImageAssert.Png(p.ToPng());
        ImageAssert.Svg(p.ToSvg());
        p.ShowDistances = true;
        ImageAssert.Png(p.ToPng());
        ImageAssert.Dump(p, "contactmap");
    }

    [Fact]
    public void ThreeDViewerHtmlPinsVersion()
    {
        var html = ProteinViewer.To3DHtml(SampleStructures.Crambin(), new Viewer3DOptions { Coloring = ResidueColoring.SecondaryStructure });
        Assert.Contains("3dmol@" + ProteinViewer.ThreeDmolVersion, html);
        Assert.Contains("createViewer", html);
        Assert.Contains("ATOM", html);
        Assert.Contains("ssJmol", html);
    }
}

public class MoleculeViewerTests
{
    [Theory]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "aspirin")]
    [InlineData("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "caffeine")]
    [InlineData("C#CC[NH3+]", "propargylammonium")]
    public void DepictsMolecule(string smiles, string name)
    {
        var mol = Molecule.FromSmiles(smiles, name);
        var v = new MoleculeViewer(mol);
        var c = v.Coordinates2D();
        Assert.Equal(mol.AtomCount, c.Length);
        Assert.Null(mol.Coordinates); // the input molecule is not modified
        ImageAssert.Png(v.ToPng());
        var svg = v.ToSvg();
        ImageAssert.Svg(svg, 500);
        Assert.Contains("<text", svg);
        if (smiles.Contains('O')) Assert.Contains(">O<", svg.Replace(" ", ""), StringComparison.Ordinal);
        Assert.Contains("3dmol", v.To3DHtml());
        ImageAssert.Dump(v, "molecule-" + name);
    }

    [Fact]
    public void LabelsHeteroatomHydrogens()
    {
        var v = new MoleculeViewer(Molecule.FromSmiles("OCCN", "ethanolamine"));
        var svg = v.ToSvg();
        Assert.Contains("OH", svg);
        Assert.Contains("NH₂", svg);
    }
}

public class ThemeTests
{
    [Fact]
    public void ThemesAreConsistent()
    {
        foreach (var t in new[] { BioTheme.Light, BioTheme.Dark })
        {
            Assert.True(t.Palette.Count >= 8);
            Assert.NotEqual(t.Up, t.Down);
            Assert.Equal(t.BaseT, t.NucleotideColor('U'));
            Assert.True(Visual.Net.Core.Colors.Color.ContrastRatio(t.Foreground, t.Background) > 7);
            var style = t.ToPlotStyle();
            Assert.Equal(t.Background, style.FigureFaceColor);
        }
        BiotechColormaps.EnsureRegistered();
        Assert.NotNull(Visual.Net.Core.Colors.Colormaps.Get(BiotechColormaps.Diverging));
    }

    [Fact]
    public void SaveWritesEveryFormat()
    {
        var dir = Directory.CreateTempSubdirectory("biotechviz");
        try
        {
            var logo = new SequenceLogo(["ACGT", "ACGA"]);
            foreach (var ext in new[] { "png", "svg", "html", "pdf" })
            {
                string path = Path.Combine(dir.FullName, "logo." + ext);
                logo.Save(path);
                Assert.True(new FileInfo(path).Length > 200);
            }
            Assert.StartsWith("%PDF", System.Text.Encoding.ASCII.GetString(File.ReadAllBytes(Path.Combine(dir.FullName, "logo.pdf"))[..4]));
            Assert.Throws<NotSupportedException>(() => logo.Save(Path.Combine(dir.FullName, "x.bmp")));
        }
        finally { dir.Delete(true); }
    }
}
