using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace Biotech.Net.Gallery.Controls;

internal static class Draw
{
    public static readonly Typeface MonoFace = new(Palette.Mono);
    public static readonly Typeface DisplayFace = new(Palette.Display, FontStyle.Normal, FontWeight.SemiBold);

    public static FormattedText Text(string s, double size, Color color, Typeface? face = null)
        => new(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face ?? MonoFace, size, new SolidColorBrush(color));

    public static void TextAt(DrawingContext ctx, string s, double x, double y, double size, Color color, double align = 0, Typeface? face = null)
    {
        var ft = Text(s, size, color, face);
        ctx.DrawText(ft, new Point(x - ft.Width * align, y));
    }

    public static Pen Pen(Color c, double w = 1, bool dashed = false)
        => new(new SolidColorBrush(c), w) { DashStyle = dashed ? new DashStyle([3, 3], 0) : null };
}

/// <summary>A labelled vertical bar chart (base or residue composition).</summary>
public sealed class BarChart : Control
{
    public sealed record Bar(string Label, double Value, Color Color, string? ValueText = null);

    private IReadOnlyList<Bar> _bars = [];

    public IReadOnlyList<Bar> Bars
    {
        get => _bars;
        set { _bars = value; InvalidateVisual(); }
    }

    /// <summary>Optional fixed maximum (else the data maximum).</summary>
    public double? Max { get; set; }

    public override void Render(DrawingContext ctx)
    {
        double w = Bounds.Width, h = Bounds.Height;
        if (_bars.Count == 0 || w < 20 || h < 30) return;
        double labelH = 18, valueH = 16;
        double plotTop = valueH, plotBottom = h - labelH;
        double max = Max ?? Math.Max(1e-9, _bars.Max(b => b.Value));
        double slot = w / _bars.Count;
        double bw = Math.Min(46, slot * 0.62);
        for (int k = 1; k <= 3; k++)
        {
            double y = plotBottom - (plotBottom - plotTop) * k / 4.0;
            ctx.DrawLine(Draw.Pen(Palette.Grid), new Point(0, y), new Point(w, y));
        }
        ctx.DrawLine(Draw.Pen(Palette.Rule), new Point(0, plotBottom + 0.5), new Point(w, plotBottom + 0.5));
        for (int i = 0; i < _bars.Count; i++)
        {
            var bar = _bars[i];
            double cx = slot * i + slot / 2;
            double bh = (plotBottom - plotTop) * Math.Clamp(bar.Value / max, 0, 1);
            var rect = new Rect(cx - bw / 2, plotBottom - bh, bw, bh);
            ctx.FillRectangle(new SolidColorBrush(bar.Color, 0.22), rect);
            ctx.FillRectangle(new SolidColorBrush(bar.Color), new Rect(rect.X, rect.Y, rect.Width, Math.Min(2, bh)));
            Draw.TextAt(ctx, bar.ValueText ?? bar.Value.ToString("0.#", CultureInfo.InvariantCulture), cx, plotBottom - bh - 15, 10.5, Palette.Ink2, 0.5);
            Draw.TextAt(ctx, bar.Label, cx, plotBottom + 3, 11, bar.Color, 0.5, Draw.DisplayFace);
        }
    }
}

/// <summary>A per-position profile (hydropathy, GC content) with a zero/threshold line and highlighted segments.</summary>
public sealed class ProfileChart : Control
{
    private double[] _values = [];

    public double[] Values
    {
        get => _values;
        set { _values = value; InvalidateVisual(); }
    }

    /// <summary>Horizontal reference lines (value, label).</summary>
    public IReadOnlyList<(double Value, string Label)> References { get; set; } = [];

    /// <summary>Highlighted x ranges (index start, end exclusive).</summary>
    public IReadOnlyList<(int Start, int End)> Segments { get; set; } = [];

    public Color Positive { get; set; } = Palette.DyeT;
    public Color Negative { get; set; } = Palette.DyeC;
    public double? Baseline { get; set; } = 0;
    public double? YMin { get; set; }
    public double? YMax { get; set; }
    public string XLabel { get; set; } = "position";
    public int XOffset { get; set; } = 1;
    public string Format { get; set; } = "0.0";

    public override void Render(DrawingContext ctx)
    {
        double w = Bounds.Width, h = Bounds.Height;
        if (_values.Length < 2 || w < 60 || h < 40) return;
        double left = 40, right = 8, top = 8, bottom = 22;
        double pw = w - left - right, ph = h - top - bottom;
        var finite = _values.Where(double.IsFinite).ToArray();
        if (finite.Length == 0) return;
        double min = YMin ?? Math.Min(finite.Min(), Baseline ?? finite.Min());
        double max = YMax ?? Math.Max(finite.Max(), Baseline ?? finite.Max());
        if (max - min < 1e-9) { max += 1; min -= 1; }
        double pad = (max - min) * 0.06;
        if (YMin is null) min -= pad;
        if (YMax is null) max += pad;
        double X(double i) => left + pw * i / (_values.Length - 1);
        double Y(double v) => top + ph * (1 - (v - min) / (max - min));

        foreach (var (s, e) in Segments)
            ctx.FillRectangle(new SolidColorBrush(Palette.DyeG, 0.10), new Rect(X(s), top, Math.Max(2, X(Math.Max(s, e - 1)) - X(s)), ph));

        // y ticks
        for (int k = 0; k <= 4; k++)
        {
            double v = min + (max - min) * k / 4;
            double y = Y(v);
            ctx.DrawLine(Draw.Pen(Palette.Grid), new Point(left, y), new Point(w - right, y));
            Draw.TextAt(ctx, v.ToString(Format, CultureInfo.InvariantCulture), left - 6, y - 7, 10, Palette.Ink3, 1);
        }
        // x ticks
        int step = NiceStep(_values.Length);
        for (int i = 0; i < _values.Length; i += step)
            Draw.TextAt(ctx, (i + XOffset).ToString(CultureInfo.InvariantCulture), X(i), top + ph + 5, 10, Palette.Ink3, 0.5);

        double baseV = Baseline ?? min;
        double by = Y(Math.Clamp(baseV, min, max));
        // filled area split by sign
        for (int sign = 0; sign < 2; sign++)
        {
            var geo = new StreamGeometry();
            using (var g = geo.Open())
            {
                g.BeginFigure(new Point(X(0), by), true);
                for (int i = 0; i < _values.Length; i++)
                {
                    double v = double.IsFinite(_values[i]) ? _values[i] : baseV;
                    v = sign == 0 ? Math.Max(v, baseV) : Math.Min(v, baseV);
                    g.LineTo(new Point(X(i), Y(v)));
                }
                g.LineTo(new Point(X(_values.Length - 1), by));
                g.EndFigure(true);
            }
            ctx.DrawGeometry(new SolidColorBrush(sign == 0 ? Positive : Negative, 0.20), null, geo);
        }
        var line = new StreamGeometry();
        using (var g = line.Open())
        {
            g.BeginFigure(new Point(X(0), Y(double.IsFinite(_values[0]) ? _values[0] : baseV)), false);
            for (int i = 1; i < _values.Length; i++) g.LineTo(new Point(X(i), Y(double.IsFinite(_values[i]) ? _values[i] : baseV)));
            g.EndFigure(false);
        }
        ctx.DrawGeometry(null, Draw.Pen(Palette.Ink, 1.2), line);
        if (Baseline is not null) ctx.DrawLine(Draw.Pen(Palette.Ink3), new Point(left, by), new Point(w - right, by));
        foreach (var (v, label) in References)
        {
            double y = Y(v);
            ctx.DrawLine(Draw.Pen(Palette.DyeG, 1, dashed: true), new Point(left, y), new Point(w - right, y));
            Draw.TextAt(ctx, label, left + 6, top, 10, Palette.DyeG, 0);
        }
    }

    private static int NiceStep(int n)
    {
        foreach (int s in new[] { 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000 })
            if (n / s <= 10) return s;
        return 10000;
    }
}

/// <summary>DSSP secondary-structure strip: one cell per residue with its one-letter code underneath.</summary>
public sealed class SecondaryStructureStrip : Control
{
    private string _ss = "", _residues = "";

    public void Set(string dssp, string residues)
    {
        _ss = dssp;
        _residues = residues;
        InvalidateVisual();
    }

    public static Color ColorOf(char c) => c switch
    {
        'H' or 'G' or 'I' => Palette.DyeT,
        'E' or 'B' => Palette.DyeG,
        'T' or 'S' => Palette.Ink3,
        _ => Palette.Rule,
    };

    public override void Render(DrawingContext ctx)
    {
        double w = Bounds.Width, h = Bounds.Height;
        int n = _ss.Length;
        if (n == 0 || w < 20) return;
        double cw = w / n;
        double barY = 16, barH = 14;
        for (int i = 0; i < n; i++)
        {
            char c = _ss[i];
            var col = ColorOf(c);
            bool regular = c is 'H' or 'G' or 'I' or 'E' or 'B';
            double hh = regular ? barH : c is 'T' or 'S' ? 6 : 2;
            ctx.FillRectangle(new SolidColorBrush(col, regular ? 0.85 : 1), new Rect(i * cw + 0.5, barY + (barH - hh) / 2, Math.Max(1, cw - 1), hh));
            if (cw >= 9 && i < _residues.Length)
                Draw.TextAt(ctx, _residues[i].ToString(), i * cw + cw / 2, barY + barH + 4, 10, Palette.Ink2, 0.5);
            if ((i + 1) % 10 == 0)
                Draw.TextAt(ctx, (i + 1).ToString(CultureInfo.InvariantCulture), i * cw + cw / 2, 0, 9.5, Palette.Ink3, 0.5);
        }
    }
}

/// <summary>A linear sequence map: ruler plus rows of features (ORFs, cut sites, guides).</summary>
public sealed class SequenceMap : Control
{
    public sealed record Mark(int Start, int End, string Label, Color Color, int Row, bool Tick = false, bool Reverse = false);

    private int _length;
    private IReadOnlyList<Mark> _marks = [];

    public void Set(int length, IReadOnlyList<Mark> marks)
    {
        _length = length;
        _marks = marks;
        InvalidateVisual();
    }

    public IReadOnlyList<string> RowLabels { get; set; } = [];

    public override void Render(DrawingContext ctx)
    {
        double w = Bounds.Width, h = Bounds.Height;
        if (_length <= 0 || w < 80) return;
        double left = RowLabels.Count > 0 ? 70 : 6, right = 10;
        double pw = w - left - right;
        double X(double p) => left + pw * p / _length;
        // ruler
        double ry = 18;
        ctx.DrawLine(Draw.Pen(Palette.Ink3), new Point(left, ry), new Point(left + pw, ry));
        int step = Step(_length);
        for (int p = 0; p <= _length; p += step)
        {
            ctx.DrawLine(Draw.Pen(Palette.Ink3), new Point(X(p), ry - 4), new Point(X(p), ry));
            Draw.TextAt(ctx, p.ToString("N0", CultureInfo.InvariantCulture), X(p), 0, 9.5, Palette.Ink3, p == 0 ? 0 : 0.5);
        }
        double rowH = 22, rowTop = ry + 12;
        for (int r = 0; r < RowLabels.Count; r++)
            Draw.TextAt(ctx, RowLabels[r], 0, rowTop + r * rowH + 3, 10, Palette.Ink3, 0, Draw.DisplayFace);
        // ticks with labels are packed greedily into rows starting at their requested row
        var rowEnds = new Dictionary<int, double>();
        foreach (var m in _marks.OrderBy(m => m.Tick).ThenBy(m => m.Start))
        {
            double y = rowTop + m.Row * rowH;
            if (m.Tick)
            {
                double x = X(m.Start);
                double labelW = Draw.Text(m.Label, 9.5, m.Color).Width + 8;
                int row = m.Row;
                while (rowEnds.TryGetValue(row, out double end) && end > x && row < m.Row + 4) row++;
                rowEnds[row] = x + labelW;
                y = rowTop + row * rowH;
                ctx.DrawLine(Draw.Pen(m.Color, 1.4), new Point(x, y), new Point(x, y + 12));
                Draw.TextAt(ctx, m.Label, x + 3, y - 1, 9.5, m.Color);
            }
            else
            {
                double x0 = X(m.Start), x1 = Math.Max(X(m.End), x0 + 2);
                ctx.FillRectangle(new SolidColorBrush(m.Color, 0.25), new Rect(x0, y + 2, x1 - x0, 10));
                ctx.FillRectangle(new SolidColorBrush(m.Color), new Rect(m.Reverse ? x0 : x1 - 2, y + 2, 2, 10));
                if (x1 - x0 > 40) Draw.TextAt(ctx, m.Label, x0 + 4, y + 13, 9.5, Palette.Ink2);
            }
        }
    }

    private static int Step(int n)
    {
        foreach (int s in new[] { 10, 20, 50, 100, 200, 250, 500, 1000, 2000, 5000, 10000 })
            if (n / s <= 8) return s;
        return 50000;
    }
}
