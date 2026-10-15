using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace Biotech.Net.Gallery.Controls;

/// <summary>
/// The Gallery's signature: a Sanger-style four-channel trace of a DNA sequence, one Gaussian peak per base in its
/// dye colour, with the base calls printed above the peaks. Peak heights are deterministic per position so the
/// same sequence always draws the same trace.
/// </summary>
public sealed class Chromatogram : Control
{
    public static readonly StyledProperty<string?> SequenceProperty =
        AvaloniaProperty.Register<Chromatogram, string?>(nameof(Sequence));

    public static readonly StyledProperty<double> BaseSpacingProperty =
        AvaloniaProperty.Register<Chromatogram, double>(nameof(BaseSpacing), 12.5);

    static Chromatogram() => AffectsRender<Chromatogram>(SequenceProperty, BaseSpacingProperty);

    public string? Sequence { get => GetValue(SequenceProperty); set => SetValue(SequenceProperty, value); }

    public double BaseSpacing { get => GetValue(BaseSpacingProperty); set => SetValue(BaseSpacingProperty, value); }

    private static readonly Typeface CallFace = new(Palette.Mono);

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        double w = b.Width, h = b.Height;
        if (w < 10 || h < 10) return;
        ctx.FillRectangle(Brushes.Transparent, new Rect(0, 0, w, h));

        string seq = (Sequence ?? "").ToUpperInvariant();
        double sp = BaseSpacing;
        int n = Math.Min(seq.Length, (int)(w / sp));
        double callBand = 14, baseY = h - 2, top = callBand + 2;
        double amp = baseY - top;

        // baseline
        ctx.DrawLine(new Pen(Palette.RuleBrush, 1), new Point(0, baseY + 0.5), new Point(w, baseY + 0.5));
        if (n == 0) return;

        // peak heights
        var heights = new double[n];
        for (int i = 0; i < n; i++)
        {
            uint hsh = (uint)(i * 2654435761u) ^ seq[i] * 40503u;
            heights[i] = 0.55 + 0.45 * ((hsh >> 8) % 1000) / 1000.0;
        }

        const string channels = "ACGT";
        double sigma = sp * 0.40;
        int samples = (int)(n * sp);
        foreach (char ch in channels)
        {
            var color = Palette.Channel(ch);
            var geo = new StreamGeometry();
            using (var g = geo.Open())
            {
                bool started = false;
                for (int px = 0; px <= samples; px += 1)
                {
                    double x = px;
                    double y = 0;
                    int center = (int)(x / sp);
                    for (int k = Math.Max(0, center - 2); k <= Math.Min(n - 1, center + 2); k++)
                    {
                        double cx = k * sp + sp / 2;
                        double d = (x - cx) / sigma;
                        double peak = seq[k] == ch ? heights[k] : 0.06 * heights[k]; // small cross-talk
                        y += peak * Math.Exp(-0.5 * d * d);
                    }
                    var p = new Point(x, baseY - Math.Min(1.05, y) * amp);
                    if (!started) { g.BeginFigure(p, false); started = true; }
                    else g.LineTo(p);
                }
                g.EndFigure(false);
            }
            ctx.DrawGeometry(null, new Pen(new SolidColorBrush(color, 0.95), 1.3), geo);
        }

        // base calls
        for (int i = 0; i < n; i++)
        {
            var ft = new FormattedText(seq[i].ToString(), CultureInfo.InvariantCulture, FlowDirection.LeftToRight, CallFace, 9.5,
                new SolidColorBrush(Palette.Channel(seq[i])));
            ctx.DrawText(ft, new Point(i * sp + sp / 2 - ft.Width / 2, 0));
        }
        // fade the tail so the ribbon ends softly
        var fade = new LinearGradientBrush
        {
            StartPoint = new RelativePoint(0, 0, RelativeUnit.Relative),
            EndPoint = new RelativePoint(1, 0, RelativeUnit.Relative),
            GradientStops = { new GradientStop(Color.FromArgb(0, 11, 16, 21), 0), new GradientStop(Palette.Bench, 1) },
        };
        double fw = Math.Min(140, w * 0.15);
        ctx.FillRectangle(fade, new Rect(w - fw, 0, fw, h));
    }
}

/// <summary>Small square dye mark for a channel (navigation rail).</summary>
public sealed class ChannelMark : Control
{
    public static readonly StyledProperty<char> ChannelProperty = AvaloniaProperty.Register<ChannelMark, char>(nameof(Channel), 'A');

    static ChannelMark() => AffectsRender<ChannelMark>(ChannelProperty);

    public char Channel { get => GetValue(ChannelProperty); set => SetValue(ChannelProperty, value); }

    public override void Render(DrawingContext ctx)
        => ctx.FillRectangle(new SolidColorBrush(Palette.Channel(Channel)), new Rect(Bounds.Size), 1);
}
