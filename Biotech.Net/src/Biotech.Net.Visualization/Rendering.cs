using System.Collections.Concurrent;
using System.Globalization;
using Matplotlib.Net.Rendering;
using SkiaSharp;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>An <see cref="IRenderer"/> decorator that scales all coordinates and font sizes (high-DPI PNG output).</summary>
internal sealed class ScaledRenderer(IRenderer inner, double scale) : IRenderer
{
    public SizeD Size => new(inner.Size.Width / scale, inner.Size.Height / scale);
    public double Dpi => inner.Dpi;

    private Pt S(Pt p) => new(p.X * scale, p.Y * scale);
    private Rect S(Rect r) => new(r.X * scale, r.Y * scale, r.Width * scale, r.Height * scale);
    private Pen S(Pen p) => p.IsInvisible ? p : p with { Width = p.Width * scale };
    private Font S(Font f) => f.Scaled(scale);

    private GraphicsPath S(GraphicsPath path)
    {
        var result = new GraphicsPath();
        int k = 0;
        var pts = path.Points;
        foreach (var verb in path.Verbs)
        {
            switch (verb)
            {
                case PathVerb.MoveTo: result.MoveTo(S(pts[k++])); break;
                case PathVerb.LineTo: result.LineTo(S(pts[k++])); break;
                case PathVerb.QuadTo: result.QuadTo(S(pts[k]), S(pts[k + 1])); k += 2; break;
                case PathVerb.CubicTo: result.CubicTo(S(pts[k]), S(pts[k + 1]), S(pts[k + 2])); k += 3; break;
                case PathVerb.Close: result.Close(); break;
            }
        }
        return result;
    }

    public void Clear(Color color) => inner.Clear(color);
    public void DrawPath(GraphicsPath path, Pen pen) => inner.DrawPath(S(path), S(pen));
    public void FillPath(GraphicsPath path, Brush brush) => inner.FillPath(S(path), brush);
    public void FillAndDrawPath(GraphicsPath path, Brush brush, Pen pen) => inner.FillAndDrawPath(S(path), brush, S(pen));
    public void DrawLine(Pt from, Pt to, Pen pen) => inner.DrawLine(S(from), S(to), S(pen));

    public void DrawPolyline(ReadOnlySpan<Pt> points, Pen pen)
    {
        var scaled = new Pt[points.Length];
        for (int i = 0; i < points.Length; i++) scaled[i] = S(points[i]);
        inner.DrawPolyline(scaled, S(pen));
    }

    public void DrawRect(Rect rect, Brush brush, Pen pen) => inner.DrawRect(S(rect), brush, S(pen));
    public void DrawEllipse(Rect bounds, Brush brush, Pen pen) => inner.DrawEllipse(S(bounds), brush, S(pen));

    public void DrawText(string text, Pt anchor, Font font, Color color, TextAlign align = TextAlign.Left,
        TextBaseline baseline = TextBaseline.Baseline, double rotationDegrees = 0)
        => inner.DrawText(text, S(anchor), S(font), color, align, baseline, rotationDegrees);

    public TextMetrics MeasureText(string text, Font font)
    {
        var m = inner.MeasureText(text, S(font));
        return new TextMetrics(m.Width / scale, m.Height / scale, m.Ascent / scale, m.Descent / scale);
    }

    public void DrawImage(ReadOnlySpan<uint> pixels, int width, int height, Rect destination, bool smooth = false)
        => inner.DrawImage(pixels, width, height, S(destination), smooth);

    public void PushClip(Rect rect) => inner.PushClip(S(rect));
    public void PopClip() => inner.PopClip();
    public void Flush() => inner.Flush();
    public void Dispose() { }
}

/// <summary>Drawing helpers shared by the custom-rendered plots.</summary>
internal static class Paint
{
    public static Pen Pen(Color c, double width = 1, LineStyle style = LineStyle.Solid) => new(c, width, style, LineCap.Round, LineJoin.Round);

    public static Brush Fill(Color c) => new(c);

    public static void Rect(IRenderer r, double x, double y, double w, double h, Color fill)
        => r.DrawRect(new Rect(x, y, w, h), new Brush(fill), Matplotlib.Net.Rendering.Pen.None);

    public static void Line(IRenderer r, double x1, double y1, double x2, double y2, Color c, double width = 1, LineStyle style = LineStyle.Solid)
        => r.DrawLine(new Pt(x1, y1), new Pt(x2, y2), Pen(c, width, style));

    public static void Polygon(IRenderer r, Color fill, params Pt[] points)
    {
        var path = new GraphicsPath().AddPolygon(points);
        r.FillPath(path, new Brush(fill));
    }

    /// <summary>Text truncated with an ellipsis to fit <paramref name="maxWidth"/>.</summary>
    public static string Fit(IRenderer r, string text, Font font, double maxWidth)
    {
        if (maxWidth <= 0) return string.Empty;
        if (r.MeasureText(text, font).Width <= maxWidth) return text;
        for (int n = text.Length - 1; n > 0; n--)
        {
            string t = text[..n] + "…";
            if (r.MeasureText(t, font).Width <= maxWidth) return t;
        }
        return string.Empty;
    }

    /// <summary>Formats a genomic coordinate compactly (bp, kb, Mb) for a view of <paramref name="span"/> bases.</summary>
    public static string Coordinate(double position, double span)
    {
        if (span >= 5_000_000) return $"{position / 1e6:0.##} Mb";
        if (span >= 50_000) return $"{position / 1e3:0.#} kb";
        return ((long)Math.Round(position)).ToString("N0", System.Globalization.CultureInfo.InvariantCulture);
    }

    /// <summary>"Nice" tick positions (1-2-5 steps) covering [min, max].</summary>
    public static double[] NiceTicks(double min, double max, int maxTicks = 8)
    {
        if (!(max > min)) return [min];
        double raw = (max - min) / Math.Max(1, maxTicks);
        double mag = Math.Pow(10, Math.Floor(Math.Log10(raw)));
        double norm = raw / mag;
        double step = (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 2.5 ? 2.5 : norm <= 5 ? 5 : 10) * mag;
        var ticks = new List<double>();
        for (double t = Math.Ceiling(min / step) * step; t <= max + step * 1e-9; t += step) ticks.Add(Math.Round(t / step) * step);
        return ticks.ToArray();
    }
}

/// <summary>
/// Letter outlines for sequence logos: the glyph is taken from a bold system font through SkiaSharp
/// (<see cref="SKFont.GetTextPath(string, SKPoint)"/>) and normalised to the unit square, so it can be stretched to any
/// height and drawn by both the PNG and the SVG renderer as a filled path.
/// </summary>
internal static class LetterGlyphs
{
    private static readonly ConcurrentDictionary<char, (PathVerb Verb, Pt[] Points)[]> _cache = new();

    /// <summary>Fills <paramref name="letter"/> stretched to exactly cover <paramref name="box"/>.</summary>
    public static void Fill(IRenderer r, char letter, Rect box, Color color)
    {
        if (box.Width <= 0 || box.Height <= 0) return;
        var glyph = _cache.GetOrAdd(letter, Outline);
        if (glyph.Length == 0)
        {
            r.DrawRect(box, new Brush(color), Pen.None);
            return;
        }
        var path = new GraphicsPath();
        Pt Map(Pt p) => new(box.X + p.X * box.Width, box.Y + p.Y * box.Height);
        foreach (var (verb, pts) in glyph)
        {
            switch (verb)
            {
                case PathVerb.MoveTo: path.MoveTo(Map(pts[0])); break;
                case PathVerb.LineTo: path.LineTo(Map(pts[0])); break;
                case PathVerb.QuadTo: path.QuadTo(Map(pts[0]), Map(pts[1])); break;
                case PathVerb.CubicTo: path.CubicTo(Map(pts[0]), Map(pts[1]), Map(pts[2])); break;
                case PathVerb.Close: path.Close(); break;
            }
        }
        r.FillPath(path, new Brush(color));
    }

    /// <summary>
    /// The unit-square outline of a letter: the embedded DejaVu Sans Bold outline for A–Z (so logos look the same
    /// on every machine, fonts installed or not), otherwise the system font (empty when no font is available).
    /// </summary>
    internal static (PathVerb Verb, Pt[] Points)[] Outline(char letter)
    {
        if (LetterOutlines.Paths.TryGetValue(char.ToUpperInvariant(letter), out var embedded)) return Parse(embedded);
        try
        {
            using var typeface = SKTypeface.FromFamilyName("Arial", SKFontStyle.Bold) ?? SKTypeface.Default;
            using var font = new SKFont(typeface, 100);
            using var skPath = font.GetTextPath(letter.ToString(), new SKPoint(0, 0));
            if (skPath is null || skPath.IsEmpty) return [];
            var b = skPath.TightBounds;
            if (b.Width <= 0 || b.Height <= 0) return [];
            Pt N(SKPoint p) => new((p.X - b.Left) / b.Width, (p.Y - b.Top) / b.Height);
            var result = new List<(PathVerb, Pt[])>();
            using var it = skPath.CreateRawIterator();
            var pts = new SKPoint[4];
            SKPathVerb verb;
            while ((verb = it.Next(pts)) != SKPathVerb.Done)
            {
                switch (verb)
                {
                    case SKPathVerb.Move: result.Add((PathVerb.MoveTo, [N(pts[0])])); break;
                    case SKPathVerb.Line: result.Add((PathVerb.LineTo, [N(pts[1])])); break;
                    case SKPathVerb.Quad: result.Add((PathVerb.QuadTo, [N(pts[1]), N(pts[2])])); break;
                    case SKPathVerb.Conic: result.Add((PathVerb.QuadTo, [N(pts[1]), N(pts[2])])); break;
                    case SKPathVerb.Cubic: result.Add((PathVerb.CubicTo, [N(pts[1]), N(pts[2]), N(pts[3])])); break;
                    case SKPathVerb.Close: result.Add((PathVerb.Close, [])); break;
                }
            }
            return result.ToArray();
        }
        catch (Exception)
        {
            return [];
        }
    }

    private static (PathVerb Verb, Pt[] Points)[] Parse(string commands)
    {
        var tokens = commands.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        var result = new List<(PathVerb, Pt[])>();
        int i = 0;
        Pt[] Read(int n)
        {
            var pts = new Pt[n];
            for (int k = 0; k < n; k++, i += 2)
                pts[k] = new Pt(double.Parse(tokens[i], CultureInfo.InvariantCulture), double.Parse(tokens[i + 1], CultureInfo.InvariantCulture));
            return pts;
        }
        while (i < tokens.Length)
        {
            switch (tokens[i++])
            {
                case "M": result.Add((PathVerb.MoveTo, Read(1))); break;
                case "L": result.Add((PathVerb.LineTo, Read(1))); break;
                case "Q": result.Add((PathVerb.QuadTo, Read(2))); break;
                case "C": result.Add((PathVerb.CubicTo, Read(3))); break;
                case "Z": result.Add((PathVerb.Close, [])); break;
                default: throw new FormatException($"Bad outline command '{tokens[i - 1]}'.");
            }
        }
        return result.ToArray();
    }
}
