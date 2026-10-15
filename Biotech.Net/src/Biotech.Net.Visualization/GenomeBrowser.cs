using Biotech.Net.Core;
using Matplotlib.Net.Rendering;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>Where a track is drawn: its full row and the data area to the right of the name gutter.</summary>
public readonly record struct TrackPlacement(GenomeTrack Track, Rect Row, Rect DataArea);

/// <summary>
/// A static, multi-track genome browser view of one region (UCSC / IGV style): ruler, gene and transcript
/// models, sequence features, coverage, variants and bases. Coordinates follow Core: 0-based half-open
/// <see cref="GenomicInterval"/>s, 1-based <see cref="Variant.Position"/>.
/// </summary>
/// <example>
/// <code>
/// var view = new GenomeBrowser(GenomicInterval.Parse("chr1:1,000-6,000"))
///     .AddRuler().AddGenes(genome.GenesIn(region)).AddCoverage(depth, start: 999).AddVariants(variants);
/// view.Save("region.svg");
/// </code>
/// </example>
public sealed class GenomeBrowser : BioPlot
{
    private readonly List<GenomeTrack> _tracks = [];

    public GenomeBrowser(GenomicInterval region, int width = 1000) : base(width, 200)
    {
        if (region.Length <= 0) throw new ArgumentException("The region must cover at least one base.", nameof(region));
        Region = region;
    }

    /// <summary>The region shown.</summary>
    public GenomicInterval Region { get; }

    /// <summary>Tracks from top to bottom.</summary>
    public IReadOnlyList<GenomeTrack> Tracks => _tracks;

    /// <summary>Width of the left gutter holding track names.</summary>
    public double NameGutter { get; set; } = 120;

    /// <summary>Outer padding.</summary>
    public double Padding { get; set; } = 12;

    /// <summary>Vertical gap between tracks.</summary>
    public double TrackGap { get; set; } = 6;

    /// <summary>Height follows the tracks; setting it has no effect.</summary>
    public override int Height
    {
        get => (int)Math.Ceiling(ComputeLayout().Sum(p => p.Row.Height) + _tracks.Count * TrackGap + HeaderHeight + 2 * Padding);
        set { }
    }

    private double HeaderHeight => string.IsNullOrEmpty(Title) ? 20 : 44;

    /// <summary>Adds any track.</summary>
    public GenomeBrowser AddTrack(GenomeTrack track)
    {
        _tracks.Add(track ?? throw new ArgumentNullException(nameof(track)));
        return this;
    }

    /// <summary>Adds a coordinate ruler.</summary>
    public GenomeBrowser AddRuler() => AddTrack(new RulerTrack());

    /// <summary>Adds gene models: each transcript is drawn (or the gene span when it has none).</summary>
    public GenomeBrowser AddGenes(IEnumerable<Gene> genes, string name = "Genes", bool showTranscripts = true)
        => AddTrack(new GeneTrack(genes, name) { ShowTranscripts = showTranscripts });

    /// <summary>Adds transcript models.</summary>
    public GenomeBrowser AddTranscripts(IEnumerable<Transcript> transcripts, string name = "Transcripts")
        => AddTrack(GeneTrack.FromTranscripts(transcripts, name));

    /// <summary>Adds annotated features; their <see cref="FeatureLocation"/> coordinates are on the region's sequence.</summary>
    public GenomeBrowser AddFeatures(IEnumerable<SequenceFeature> features, string name = "Features")
        => AddTrack(new FeatureTrack(features, name));

    /// <summary>Adds labelled intervals (peaks, BED records, ...).</summary>
    public GenomeBrowser AddIntervals(IEnumerable<(GenomicInterval Interval, string Label)> intervals, string name = "Intervals")
        => AddTrack(FeatureTrack.FromIntervals(intervals, name));

    /// <summary>Adds a per-base signal (read depth, conservation). <paramref name="start"/> is the 0-based position of <c>values[0]</c>.</summary>
    public GenomeBrowser AddCoverage(double[] values, long start, string name = "Coverage", Color? color = null)
        => AddTrack(new CoverageTrack(values, start, name) { Color = color });

    /// <summary>Adds variants as a lollipop track.</summary>
    public GenomeBrowser AddVariants(IEnumerable<Variant> variants, string name = "Variants")
        => AddTrack(new VariantTrack(variants, name));

    /// <summary>Adds the bases of <paramref name="sequence"/>; <paramref name="offset"/> is the 0-based position of its first base.</summary>
    public GenomeBrowser AddSequence(Sequence sequence, long offset = 0, string name = "Sequence")
        => AddTrack(new SequenceTrack(sequence, offset, name));

    /// <summary>Data-area x of a 0-based position (left edge of the base).</summary>
    public double ToX(double position, Rect dataArea)
        => dataArea.X + (position - Region.Start) / Region.Length * dataArea.Width;

    /// <summary>Bases per pixel at the current width.</summary>
    public double BasesPerPixel => Region.Length / Math.Max(1, Width - NameGutter - 2 * Padding);

    /// <summary>Deterministic placement of each track.</summary>
    public IReadOnlyList<TrackPlacement> ComputeLayout()
    {
        var result = new List<TrackPlacement>(_tracks.Count);
        double y = Padding + HeaderHeight;
        double dataX = Padding + NameGutter, dataW = Math.Max(10, Width - dataX - Padding);
        foreach (var t in _tracks)
        {
            double h = t.MeasureHeight(this, dataW);
            var row = new Rect(Padding, y, Width - 2 * Padding, h);
            result.Add(new TrackPlacement(t, row, new Rect(dataX, y, dataW, h)));
            y += h + TrackGap;
        }
        return result;
    }

    /// <inheritdoc />
    public override void Render(IRenderer r)
    {
        var th = Theme;
        r.Clear(th.Background);
        var layout = ComputeLayout();
        double top = Padding;
        if (!string.IsNullOrEmpty(Title))
        {
            r.DrawText(Title, new Pt(Padding, top), th.TitleFont, th.Foreground, TextAlign.Left, TextBaseline.Top);
            top += 24;
        }
        string where = $"{Region.Chromosome}:{Region.Start + 1:N0}-{Region.End:N0}  ·  {Region.Length:N0} bp";
        r.DrawText(where, new Pt(Padding, top), th.SmallFont, th.Muted, TextAlign.Left, TextBaseline.Top);

        for (int i = 0; i < layout.Count; i++)
        {
            var p = layout[i];
            if (i % 2 == 1) r.DrawRect(p.Row.Inflate(0, TrackGap / 2), new Brush(th.Grid.WithAlpha(th.IsDark ? 0.25 : 0.35)), Pen.None);
            string label = Paint.Fit(r, p.Track.Name, th.Font(th.FontSize * 0.9, bold: true), NameGutter - 10);
            r.DrawText(label, new Pt(p.Row.X, p.Row.Y + Math.Min(p.Row.Height / 2, 14)), th.Font(th.FontSize * 0.9, bold: true), th.Foreground, TextAlign.Left, TextBaseline.Middle);
            r.PushClip(p.DataArea.Inflate(0, 2));
            p.Track.Draw(new TrackContext(r, this, p.DataArea, th));
            r.PopClip();
        }
    }
}

/// <summary>Everything a track needs to draw itself.</summary>
public readonly record struct TrackContext(IRenderer Renderer, GenomeBrowser Browser, Rect Area, BioTheme Theme)
{
    /// <summary>x of a 0-based position.</summary>
    public double X(double position) => Browser.ToX(position, Area);

    /// <summary>The region shown.</summary>
    public GenomicInterval Region => Browser.Region;
}

/// <summary>A horizontal track of a <see cref="GenomeBrowser"/>. Derive to add custom tracks.</summary>
public abstract class GenomeTrack(string name)
{
    /// <summary>Name shown in the gutter.</summary>
    public string Name { get; set; } = name;

    /// <summary>Height needed for the current view.</summary>
    public abstract double MeasureHeight(GenomeBrowser browser, double dataWidth);

    /// <summary>Draws the track inside <see cref="TrackContext.Area"/>.</summary>
    public abstract void Draw(TrackContext context);
}

/// <summary>Coordinate ruler with 1-2-5 ticks.</summary>
public sealed class RulerTrack() : GenomeTrack("")
{
    public override double MeasureHeight(GenomeBrowser browser, double dataWidth) => 30;

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme; var a = c.Area;
        double y = a.Y + 6;
        Paint.Line(r, a.X, y, a.Right, y, th.Border, 1);
        double span = c.Region.Length;
        int maxTicks = Math.Max(2, (int)(a.Width / 90));
        foreach (double t in Paint.NiceTicks(c.Region.Start, c.Region.End, maxTicks))
        {
            if (t < c.Region.Start || t > c.Region.End) continue;
            double x = c.X(t);
            Paint.Line(r, x, y, x, y + 6, th.Border, 1);
            string label = Paint.Coordinate(t + 1, span);
            double half = r.MeasureText(label, th.SmallFont).Width / 2;
            var align = x - half < a.X ? TextAlign.Left : x + half > a.Right ? TextAlign.Right : TextAlign.Center;
            r.DrawText(label, new Pt(x, y + 9), th.SmallFont, th.Muted, align, TextBaseline.Top);
        }
    }
}

/// <summary>Gene / transcript models: thick CDS exons, thin UTRs, intron lines with strand chevrons, packed in rows.</summary>
public sealed class GeneTrack : GenomeTrack
{
    private readonly List<GeneModel> _models;

    /// <summary>One drawable model: a transcript (or a gene span) with its label.</summary>
    public sealed record GeneModel(string Label, GenomicInterval Span, IReadOnlyList<GenomicInterval> Exons, GenomicInterval? CodingRegion, Strand Strand, string Biotype);

    public GeneTrack(IEnumerable<Gene> genes, string name = "Genes") : base(name)
    {
        Genes = genes.ToList();
        _models = [];
    }

    private GeneTrack(List<GeneModel> models, string name) : base(name)
    {
        Genes = [];
        _models = models;
    }

    /// <summary>A track of transcripts without gene objects.</summary>
    public static GeneTrack FromTranscripts(IEnumerable<Transcript> transcripts, string name = "Transcripts")
        => new(transcripts.Select(t => new GeneModel(t.Id, t.Span, t.Exons, t.CodingRegion, t.Exons[0].Strand, "transcript")).ToList(), name);

    /// <summary>Genes of the track.</summary>
    public IReadOnlyList<Gene> Genes { get; }

    /// <summary>Draw every transcript (true) or one model per gene (false: the first transcript or the gene span).</summary>
    public bool ShowTranscripts { get; set; } = true;

    /// <summary>Height of one row.</summary>
    public double RowHeight { get; set; } = 30;

    /// <summary>The models to draw, in input order.</summary>
    public IReadOnlyList<GeneModel> Models()
    {
        if (Genes.Count == 0) return _models;
        var list = new List<GeneModel>();
        foreach (var g in Genes)
        {
            if (g.Transcripts.Count == 0)
            {
                list.Add(new GeneModel(g.Name, g.Location, [g.Location], null, g.Strand, g.Biotype));
                continue;
            }
            foreach (var t in ShowTranscripts ? g.Transcripts : g.Transcripts.Take(1))
            {
                string label = ShowTranscripts && g.Transcripts.Count > 1 ? $"{g.Name} · {t.Id}" : g.Name;
                var strand = t.Exons[0].Strand == Strand.Unknown ? g.Strand : t.Exons[0].Strand;
                list.Add(new GeneModel(label, t.Span, t.Exons, t.CodingRegion, strand, g.Biotype));
            }
        }
        return list;
    }

    /// <summary>
    /// Greedy row packing in pixel space: models sorted by start go into the first row whose last model (plus its
    /// label) ends before they begin. Returns the row of each model in <see cref="Models"/> order.
    /// </summary>
    public int[] PackRows(GenomeBrowser browser, double dataWidth)
    {
        var models = Models();
        var rows = new int[models.Count];
        var rowEnds = new List<double>();
        double pxPerBase = dataWidth / browser.Region.Length;
        foreach (int i in Enumerable.Range(0, models.Count).OrderBy(i => models[i].Span.Start))
        {
            var m = models[i];
            if (m.Span.Chromosome != browser.Region.Chromosome || !m.Span.Overlaps(browser.Region)) { rows[i] = -1; continue; }
            double x0 = (Math.Max(m.Span.Start, browser.Region.Start) - browser.Region.Start) * pxPerBase;
            double x1 = (Math.Min(m.Span.End, browser.Region.End) - browser.Region.Start) * pxPerBase;
            double labelEnd = x0 + m.Label.Length * 7 + 4;
            double end = Math.Max(x1, labelEnd) + 8;
            int row = rowEnds.FindIndex(e => e <= x0);
            if (row < 0) { row = rowEnds.Count; rowEnds.Add(end); } else rowEnds[row] = end;
            rows[i] = row;
        }
        return rows;
    }

    public override double MeasureHeight(GenomeBrowser browser, double dataWidth)
    {
        var rows = PackRows(browser, dataWidth);
        int n = rows.Length == 0 ? 1 : Math.Max(1, rows.Max() + 1);
        return n * RowHeight + 4;
    }

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme;
        var models = Models();
        var rows = PackRows(c.Browser, c.Area.Width);
        for (int i = 0; i < models.Count; i++)
        {
            if (rows[i] < 0) continue;
            var m = models[i];
            double top = c.Area.Y + rows[i] * RowHeight;
            double mid = top + 19;
            var color = m.Biotype is "protein_coding" or "transcript" or "mRNA" ? th.Accent : th.CategoryColor(4);
            double x0 = Math.Max(c.Area.X, c.X(m.Span.Start)), x1 = Math.Min(c.Area.Right, c.X(m.Span.End));
            r.DrawText(Paint.Fit(r, m.Label, th.SmallFont, Math.Max(40, c.Area.Right - x0)), new Pt(x0, top + 1), th.SmallFont, th.Foreground, TextAlign.Left, TextBaseline.Top);

            // Intron line and strand chevrons.
            Paint.Line(r, x0, mid, x1, mid, color, 1.2);
            if (m.Strand is Strand.Forward or Strand.Reverse)
            {
                int dir = m.Strand == Strand.Forward ? 1 : -1;
                for (double x = x0 + 10; x < x1 - 6; x += 22)
                {
                    r.DrawLine(new Pt(x - 2.5 * dir, mid - 3.5), new Pt(x + 2.5 * dir, mid), Paint.Pen(color, 1.2));
                    r.DrawLine(new Pt(x + 2.5 * dir, mid), new Pt(x - 2.5 * dir, mid + 3.5), Paint.Pen(color, 1.2));
                }
            }

            // Exons: UTR thin, CDS thick.
            foreach (var e in m.Exons)
            {
                DrawBlock(r, c, e.Start, e.End, mid, m.CodingRegion is null ? 12 : 7, color);
                if (m.CodingRegion is { } cds)
                {
                    long s = Math.Max(e.Start, cds.Start), en = Math.Min(e.End, cds.End);
                    if (en > s) DrawBlock(r, c, s, en, mid, 12, color);
                }
            }
        }
    }

    private static void DrawBlock(IRenderer r, TrackContext c, long start, long end, double mid, double height, Color color)
    {
        double x0 = c.X(start), x1 = c.X(end);
        if (x1 < c.Area.X || x0 > c.Area.Right) return;
        double w = Math.Max(1, x1 - x0);
        r.DrawRect(new Rect(x0, mid - height / 2, w, height), new Brush(color), Pen.None);
    }
}

/// <summary>Annotated features as boxes (joined segments connected by a line) with a strand arrow tip.</summary>
public sealed class FeatureTrack : GenomeTrack
{
    /// <summary>A drawable feature.</summary>
    public sealed record FeatureItem(string Label, string Type, IReadOnlyList<(long Start, long End)> Segments, Strand Strand);

    private readonly List<FeatureItem> _items;

    public FeatureTrack(IEnumerable<SequenceFeature> features, string name = "Features") : base(name)
        => _items = features.Select(f => new FeatureItem(f.Name, f.Type,
            f.Location.Segments.Select(s => ((long)s.Start, (long)s.End)).OrderBy(s => s.Item1).ToList(), f.Location.Strand)).ToList();

    private FeatureTrack(List<FeatureItem> items, string name) : base(name) => _items = items;

    /// <summary>A track of labelled intervals.</summary>
    public static FeatureTrack FromIntervals(IEnumerable<(GenomicInterval Interval, string Label)> intervals, string name = "Intervals")
        => new(intervals.Select(i => new FeatureItem(i.Label, "region", [(i.Interval.Start, i.Interval.End)], i.Interval.Strand)).ToList(), name);

    /// <summary>The features.</summary>
    public IReadOnlyList<FeatureItem> Items => _items;

    /// <summary>Height of one row.</summary>
    public double RowHeight { get; set; } = 28;

    /// <summary>Row of each item (−1 when outside the view).</summary>
    public int[] PackRows(GenomeBrowser browser, double dataWidth)
    {
        var rows = new int[_items.Count];
        var rowEnds = new List<double>();
        double pxPerBase = dataWidth / browser.Region.Length;
        foreach (int i in Enumerable.Range(0, _items.Count).OrderBy(i => _items[i].Segments[0].Start))
        {
            var it = _items[i];
            long s = it.Segments[0].Start, e = it.Segments[^1].End;
            if (e <= browser.Region.Start || s >= browser.Region.End) { rows[i] = -1; continue; }
            double x0 = (Math.Max(s, browser.Region.Start) - browser.Region.Start) * pxPerBase;
            double x1 = (Math.Min(e, browser.Region.End) - browser.Region.Start) * pxPerBase;
            double end = Math.Max(x1, x0 + it.Label.Length * 7) + 8;
            int row = rowEnds.FindIndex(v => v <= x0);
            if (row < 0) { row = rowEnds.Count; rowEnds.Add(end); } else rowEnds[row] = end;
            rows[i] = row;
        }
        return rows;
    }

    public override double MeasureHeight(GenomeBrowser browser, double dataWidth)
    {
        var rows = PackRows(browser, dataWidth);
        return (rows.Length == 0 ? 1 : Math.Max(1, rows.Max() + 1)) * RowHeight + 4;
    }

    /// <summary>Colour of a feature type.</summary>
    public static Color TypeColor(BioTheme th, string type) => type.ToLowerInvariant() switch
    {
        "cds" => th.CategoryColor(0),
        "gene" => th.CategoryColor(9),
        "mrna" or "exon" or "transcript" => th.CategoryColor(2),
        "promoter" or "regulatory" or "tata_signal" => th.CategoryColor(1),
        "rrna" or "trna" or "ncrna" => th.CategoryColor(4),
        "repeat_region" or "misc_feature" => th.CategoryColor(7),
        "primer_bind" or "protein_bind" => th.CategoryColor(3),
        _ => th.CategoryColor(6),
    };

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme;
        var rows = PackRows(c.Browser, c.Area.Width);
        for (int i = 0; i < _items.Count; i++)
        {
            if (rows[i] < 0) continue;
            var it = _items[i];
            var color = TypeColor(th, it.Type);
            double top = c.Area.Y + rows[i] * RowHeight, mid = top + 18, h = 10;
            double sx = Math.Max(c.Area.X, c.X(it.Segments[0].Start)), ex = Math.Min(c.Area.Right, c.X(it.Segments[^1].End));
            r.DrawText(Paint.Fit(r, it.Label, th.SmallFont, Math.Max(40, c.Area.Right - sx)), new Pt(sx, top + 1), th.SmallFont, th.Foreground, TextAlign.Left, TextBaseline.Top);
            if (it.Segments.Count > 1) Paint.Line(r, sx, mid, ex, mid, color, 1);
            for (int k = 0; k < it.Segments.Count; k++)
            {
                var (s, e) = it.Segments[k];
                double x0 = c.X(s), x1 = c.X(e);
                double w = Math.Max(1.5, x1 - x0);
                bool tipRight = it.Strand == Strand.Forward && k == it.Segments.Count - 1;
                bool tipLeft = it.Strand == Strand.Reverse && k == 0;
                double tip = Math.Min(8, w / 2);
                if (tipRight)
                    Paint.Polygon(r, color, new Pt(x0, mid - h / 2), new Pt(x0 + w - tip, mid - h / 2), new Pt(x0 + w, mid), new Pt(x0 + w - tip, mid + h / 2), new Pt(x0, mid + h / 2));
                else if (tipLeft)
                    Paint.Polygon(r, color, new Pt(x0 + tip, mid - h / 2), new Pt(x0 + w, mid - h / 2), new Pt(x0 + w, mid + h / 2), new Pt(x0 + tip, mid + h / 2), new Pt(x0, mid));
                else
                    r.DrawRect(new Rect(x0, mid - h / 2, w, h), new Brush(color), Pen.None);
            }
        }
    }
}

/// <summary>A filled signal track (read depth, conservation), reduced to the maximum per pixel column.</summary>
public sealed class CoverageTrack : GenomeTrack
{
    public CoverageTrack(double[] values, long start, string name = "Coverage") : base(name)
    {
        Values = values ?? throw new ArgumentNullException(nameof(values));
        Start = start;
    }

    /// <summary>Per-base values.</summary>
    public double[] Values { get; }

    /// <summary>0-based position of <c>Values[0]</c>.</summary>
    public long Start { get; }

    /// <summary>Track height.</summary>
    public double TrackHeight { get; set; } = 70;

    /// <summary>Fill colour (theme accent palette colour by default).</summary>
    public Color? Color { get; set; }

    /// <summary>Fixed y maximum (null = data maximum in view).</summary>
    public double? Max { get; set; }

    public override double MeasureHeight(GenomeBrowser browser, double dataWidth) => TrackHeight;

    /// <summary>Maximum value in each of <paramref name="columns"/> pixel columns across the region (NaN where no data).</summary>
    public double[] Bin(GenomicInterval region, int columns)
    {
        var bins = new double[columns];
        Array.Fill(bins, double.NaN);
        double basesPerCol = (double)region.Length / columns;
        for (int col = 0; col < columns; col++)
        {
            long a = region.Start + (long)Math.Floor(col * basesPerCol);
            long b = Math.Max(a + 1, region.Start + (long)Math.Floor((col + 1) * basesPerCol));
            for (long p = a; p < b; p++)
            {
                long k = p - Start;
                if (k < 0 || k >= Values.Length) continue;
                double v = Values[k];
                if (double.IsNaN(v)) continue;
                bins[col] = double.IsNaN(bins[col]) ? v : Math.Max(bins[col], v);
            }
        }
        return bins;
    }

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme; var a = c.Area;
        int cols = Math.Max(1, (int)a.Width);
        var bins = Bin(c.Region, cols);
        double max = Max ?? bins.Where(v => !double.IsNaN(v)).DefaultIfEmpty(0).Max();
        if (max <= 0) max = 1;
        var color = Color ?? th.CategoryColor(2);
        double baseY = a.Bottom - 2, h = a.Height - 14;
        var path = new GraphicsPath();
        path.MoveTo(a.X, baseY);
        for (int i = 0; i < cols; i++)
        {
            double v = double.IsNaN(bins[i]) ? 0 : Math.Clamp(bins[i] / max, 0, 1);
            double x = a.X + i * a.Width / cols;
            path.LineTo(x, baseY - v * h);
            path.LineTo(x + a.Width / cols, baseY - v * h);
        }
        path.LineTo(a.Right, baseY);
        path.Close();
        r.FillAndDrawPath(path, new Brush(color.WithAlpha(0.55)), Paint.Pen(color, 1));
        Paint.Line(r, a.X, baseY, a.Right, baseY, th.Border, 1);
        r.DrawText($"0–{max:G4}", new Pt(a.Right - 2, a.Y + 1), th.SmallFont, th.Muted, TextAlign.Right, TextBaseline.Top);
    }
}

/// <summary>Variants as lollipops: stem height by alternate-allele frequency (or stacked count), head colour by type.</summary>
public sealed class VariantTrack : GenomeTrack
{
    public VariantTrack(IEnumerable<Variant> variants, string name = "Variants") : base(name) => Variants = variants.ToList();

    /// <summary>The variants.</summary>
    public IReadOnlyList<Variant> Variants { get; }

    /// <summary>Track height.</summary>
    public double TrackHeight { get; set; } = 80;

    /// <summary>Show "Ref&gt;Alt" or ID labels when there is room.</summary>
    public bool ShowLabels { get; set; } = true;

    public override double MeasureHeight(GenomeBrowser browser, double dataWidth) => TrackHeight;

    /// <summary>Head colour of a variant.</summary>
    public static Color TypeColor(BioTheme th, Variant v) => v.Type switch
    {
        VariantType.Snv => v.IsTransition ? th.CategoryColor(0) : th.CategoryColor(1),
        VariantType.Insertion => th.CategoryColor(2),
        VariantType.Deletion => th.CategoryColor(3),
        _ => th.CategoryColor(4),
    };

    /// <summary>Relative stem height (0..1) of each variant, in input order.</summary>
    public double[] StemHeights()
    {
        var h = new double[Variants.Count];
        bool anyGenotypes = Variants.Any(v => v.Genotypes.Count > 0);
        for (int i = 0; i < h.Length; i++)
        {
            double af = anyGenotypes ? Variants[i].AlternateAlleleFrequency() : double.NaN;
            h[i] = double.IsNaN(af) ? 0.6 : 0.2 + 0.8 * Math.Clamp(af, 0, 1);
        }
        return h;
    }

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme; var a = c.Area;
        double baseY = a.Bottom - 6, maxStem = a.Height - 26;
        Paint.Line(r, a.X, baseY, a.Right, baseY, th.Border, 1.5);
        var heights = StemHeights();
        var stackAt = new Dictionary<long, int>();
        double lastLabelX = double.NegativeInfinity;
        foreach (int i in Enumerable.Range(0, Variants.Count).OrderBy(i => Variants[i].Position))
        {
            var v = Variants[i];
            if (v.Chromosome != c.Region.Chromosome) continue;
            long p0 = v.Position - 1;
            if (p0 < c.Region.Start || p0 >= c.Region.End) continue;
            int stack = stackAt.TryGetValue(p0, out int s) ? s : 0;
            stackAt[p0] = stack + 1;
            double x = c.X(p0 + 0.5);
            double top = baseY - heights[i] * maxStem - stack * 10;
            var col = TypeColor(th, v);
            Paint.Line(r, x, baseY, x, top, th.Muted, 1);
            double rad = v.IsPass ? 4.5 : 3.5;
            r.DrawEllipse(new Rect(x - rad, top - rad, 2 * rad, 2 * rad), new Brush(v.IsPass ? col : th.Panel), Paint.Pen(col, 1.4));
            if (ShowLabels && x - lastLabelX > 46)
            {
                string label = v.Id ?? $"{Short(v.Reference)}>{Short(v.Alternates.Count > 0 ? v.Alternates[0] : ".")}";
                r.DrawText(label, new Pt(x, top - rad - 2), th.Font(th.FontSize * 0.75), th.Muted, TextAlign.Center, TextBaseline.Bottom);
                lastLabelX = x;
            }
        }
    }

    private static string Short(string allele) => allele.Length <= 4 ? allele : allele[..3] + "…";
}

/// <summary>Bases of a sequence: letters when zoomed in, coloured bars when moderately zoomed, a note otherwise.</summary>
public sealed class SequenceTrack : GenomeTrack
{
    public SequenceTrack(Sequence sequence, long offset = 0, string name = "Sequence") : base(name)
    {
        Sequence = sequence ?? throw new ArgumentNullException(nameof(sequence));
        Offset = offset;
    }

    /// <summary>The sequence.</summary>
    public Sequence Sequence { get; }

    /// <summary>0-based position of the first base.</summary>
    public long Offset { get; }

    public override double MeasureHeight(GenomeBrowser browser, double dataWidth) => 20;

    public override void Draw(TrackContext c)
    {
        var r = c.Renderer; var th = c.Theme; var a = c.Area;
        double pxPerBase = a.Width / c.Region.Length;
        if (pxPerBase < 1)
        {
            r.DrawText("zoom in to see bases", new Pt(a.X + a.Width / 2, a.Y + a.Height / 2), th.SmallFont, th.Muted, TextAlign.Center, TextBaseline.Middle);
            return;
        }
        var font = th.Font(Math.Min(13, pxPerBase * 0.95), bold: true);
        for (long p = c.Region.Start; p < c.Region.End; p++)
        {
            long k = p - Offset;
            if (k < 0 || k >= Sequence.Length) continue;
            char b = (char)Sequence[(int)k];
            double x = c.X(p);
            var col = th.NucleotideColor(b);
            if (pxPerBase >= 7)
                r.DrawText(b.ToString(), new Pt(x + pxPerBase / 2, a.Y + a.Height / 2), font, col, TextAlign.Center, TextBaseline.Middle);
            else
                r.DrawRect(new Rect(x, a.Y + 4, Math.Max(1, pxPerBase), a.Height - 8), new Brush(col), Pen.None);
        }
    }
}
