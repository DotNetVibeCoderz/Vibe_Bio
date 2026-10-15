using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

namespace Biotech.Net.Visualization.Tests;

internal static class TestData
{
    public static readonly GenomicInterval Region = new("chr1", 1000, 6000);

    public static List<Gene> Genes()
    {
        var g1 = new Gene("G1", new GenomicInterval("chr1", 1200, 3400, Strand.Forward), "ALPHA");
        g1.Transcripts.Add(new Transcript("T1.1", "G1",
        [
            new("chr1", 1200, 1500, Strand.Forward), new("chr1", 2000, 2300, Strand.Forward), new("chr1", 3000, 3400, Strand.Forward),
        ]) { CodingRegion = new GenomicInterval("chr1", 1350, 3200, Strand.Forward) });
        g1.Transcripts.Add(new Transcript("T1.2", "G1",
        [
            new("chr1", 1200, 1500, Strand.Forward), new("chr1", 3000, 3400, Strand.Forward),
        ]) { CodingRegion = new GenomicInterval("chr1", 1350, 3200, Strand.Forward) });
        var g2 = new Gene("G2", new GenomicInterval("chr1", 3800, 5600, Strand.Reverse), "BETA", "lncRNA");
        g2.Transcripts.Add(new Transcript("T2.1", "G2", [new("chr1", 3800, 4200, Strand.Reverse), new("chr1", 5000, 5600, Strand.Reverse)]));
        var g3 = new Gene("G3", new GenomicInterval("chr1", 2500, 2900, Strand.Reverse), "GAMMA");
        return [g1, g2, g3];
    }

    public static List<Variant> Variants() =>
    [
        new("chr1", 1400, "A", "G", "rs1"),
        new("chr1", 2100, "C", "T"),
        new("chr1", 2101, "C", "A"),
        new("chr1", 4000, "AT", "A"),
        new("chr1", 5200, "G", "GTT") { Filter = "LowQual" },
    ];

    public static double[] Coverage()
    {
        var v = new double[5000];
        for (int i = 0; i < v.Length; i++) v[i] = 20 + 15 * Math.Sin(i / 300.0) + (i is > 1200 and < 1500 ? 30 : 0);
        return v;
    }

    public static DifferentialExpressionResult DeResult(int n = 400)
    {
        var rng = new Random(7);
        var genes = Enumerable.Range(0, n).Select(i => $"gene{i}").ToArray();
        var lfc = new double[n]; var p = new double[n]; var bm = new double[n]; var se = new double[n]; var stat = new double[n];
        for (int i = 0; i < n; i++)
        {
            lfc[i] = (rng.NextDouble() - 0.5) * 2 + (i < 20 ? 3 : i < 35 ? -3 : 0);
            p[i] = i < 35 ? Math.Pow(10, -3 - rng.NextDouble() * 10) : rng.NextDouble();
            bm[i] = Math.Exp(rng.NextDouble() * 8);
            se[i] = 0.3; stat[i] = lfc[i] / se[i];
        }
        var padj = VisualizationHelpers.Bh(p);
        return new DifferentialExpressionResult("DESeq2-Wald", "treated vs control", genes, bm, lfc, se, stat, p, padj);
    }

    public static double[,] Expression(int rows = 30, int cols = 8)
    {
        var rng = new Random(3);
        var m = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i, j] = rng.NextDouble() + (i < rows / 2 == j < cols / 2 ? 3 : 0);
        return m;
    }
}

internal static class VisualizationHelpers
{
    public static double[] Bh(double[] p) => DifferentialPlot.BenjaminiHochberg(p);
}

internal static class ImageAssert
{
    public static (int Width, int Height) Png(byte[] png, int minBytes = 2000)
    {
        Assert.True(png.Length >= minBytes, $"PNG too small: {png.Length} bytes");
        Assert.Equal(new byte[] { 0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A }, png[..8]);
        int w = (png[16] << 24) | (png[17] << 16) | (png[18] << 8) | png[19];
        int h = (png[20] << 24) | (png[21] << 16) | (png[22] << 8) | png[23];
        using var bmp = SkiaSharp.SKBitmap.Decode(png);
        Assert.NotNull(bmp);
        // Non-trivial content: more than a handful of distinct colours.
        var colours = new HashSet<uint>();
        for (int y = 0; y < bmp.Height; y += 3)
            for (int x = 0; x < bmp.Width; x += 3)
                colours.Add((uint)bmp.GetPixel(x, y));
        Assert.True(colours.Count > 4, $"Image looks blank ({colours.Count} colours)");
        return (w, h);
    }

    public static void Svg(string svg, int minLength = 1000)
    {
        Assert.StartsWith("<?xml", svg);
        Assert.Contains("<svg", svg);
        Assert.EndsWith("</svg>", svg.TrimEnd());
        Assert.True(svg.Length >= minLength, $"SVG too short: {svg.Length}");
        _ = System.Xml.Linq.XDocument.Parse(svg);
    }

    /// <summary>Writes the plot to $BIOTECH_VIZ_OUT for visual inspection (no-op otherwise).</summary>
    public static void Dump(BioPlot plot, string name)
    {
        var dir = Environment.GetEnvironmentVariable("BIOTECH_VIZ_OUT");
        if (string.IsNullOrEmpty(dir)) return;
        Directory.CreateDirectory(dir);
        plot.Save(Path.Combine(dir, name + ".png"));
        plot.Save(Path.Combine(dir, name + ".svg"));
        plot.Save(Path.Combine(dir, name + ".html"));
    }
}
