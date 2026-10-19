using System.Globalization;
using System.Net;
using System.Text;

namespace Biotech.Net.WebGallery.Services;

/// <summary>
/// Small inline-SVG instruments drawn in the chassis style (bar meter, profile trace, linear map). Library figures
/// (genome browser, volcano, molecule, Cα trace) come from Biotech.Net.Visualization instead.
/// </summary>
public static class Charts
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public sealed record Bar(string Label, double Value, string Color, string Caption);

    /// <summary>Vertical bars with a caption above each; the bars are tinted wells with a bright cap, like the desktop meter.</summary>
    public static string Bars(IReadOnlyList<Bar> bars, double width = 460, double height = 180)
    {
        var sb = Begin(width, height, "Bar chart");
        double top = 18, bottom = height - 18, max = Math.Max(1e-9, bars.Count == 0 ? 1 : bars.Max(b => b.Value));
        double slot = width / Math.Max(1, bars.Count), bw = Math.Min(46, slot * 0.62);
        for (int i = 1; i <= 3; i++)
        {
            double y = bottom - (bottom - top) * i / 3.0;
            sb.Append(Inv, $"<line x1=\"0\" x2=\"{width}\" y1=\"{y:0.#}\" y2=\"{y:0.#}\" stroke=\"#1A2530\"/>");
        }
        for (int i = 0; i < bars.Count; i++)
        {
            var b = bars[i];
            double x = slot * i + (slot - bw) / 2, h = (bottom - top) * b.Value / max, y = bottom - h;
            sb.Append(Inv, $"<rect x=\"{x:0.#}\" y=\"{y:0.#}\" width=\"{bw:0.#}\" height=\"{Math.Max(0, h):0.#}\" fill=\"{b.Color}\" fill-opacity=\"0.26\"/>");
            sb.Append(Inv, $"<rect x=\"{x:0.#}\" y=\"{y:0.#}\" width=\"{bw:0.#}\" height=\"2\" fill=\"{b.Color}\"/>");
            sb.Append(Inv, $"<text x=\"{x + bw / 2:0.#}\" y=\"{y - 5:0.#}\" text-anchor=\"middle\" class=\"t-cap\">{E(b.Caption)}</text>");
            sb.Append(Inv, $"<text x=\"{x + bw / 2:0.#}\" y=\"{height - 3:0.#}\" text-anchor=\"middle\" fill=\"{b.Color}\" class=\"t-lab\">{E(b.Label)}</text>");
        }
        sb.Append(Inv, $"<line x1=\"0\" x2=\"{width}\" y1=\"{bottom + 0.5:0.#}\" y2=\"{bottom + 0.5:0.#}\" stroke=\"#233140\"/>");
        return End(sb);
    }

    /// <summary>A profile trace around a baseline: above in <paramref name="positive"/>, below in <paramref name="negative"/>.</summary>
    public static string Profile(IReadOnlyList<double> values, double baseline, string positive, string negative,
        string format = "0", double width = 460, double height = 180, (double Value, string Label)? reference = null,
        IReadOnlyList<(int Start, int End)>? segments = null)
    {
        var sb = Begin(width, height, "Profile");
        if (values.Count < 2) return End(sb);
        double left = 34, right = width - 6, top = 10, bottom = height - 20;
        double min = Math.Min(values.Min(), baseline), max = Math.Max(values.Max(), baseline);
        if (reference is { } r0) { min = Math.Min(min, r0.Value); max = Math.Max(max, r0.Value); }
        double pad = (max - min) * 0.08 + 1e-9;
        min -= pad; max += pad;
        double X(double i) => left + (right - left) * i / (values.Count - 1);
        double Y(double v) => bottom - (bottom - top) * (v - min) / (max - min);
        for (int k = 0; k <= 4; k++)
        {
            double v = min + (max - min) * k / 4, y = Y(v);
            sb.Append(Inv, $"<line x1=\"{left}\" x2=\"{right}\" y1=\"{y:0.#}\" y2=\"{y:0.#}\" stroke=\"#1A2530\"/>");
            sb.Append(Inv, $"<text x=\"{left - 6}\" y=\"{y + 3:0.#}\" text-anchor=\"end\" class=\"t-cap\">{v.ToString(format, Inv)}</text>");
        }
        if (segments is not null)
            foreach (var (s, e) in segments)
            {
                double x0 = X(Math.Clamp(s, 0, values.Count - 1)), x1 = X(Math.Clamp(e, 0, values.Count - 1));
                sb.Append(Inv, $"<rect x=\"{x0:0.#}\" y=\"{top}\" width=\"{Math.Max(1, x1 - x0):0.#}\" height=\"{bottom - top}\" fill=\"#F0B03F\" fill-opacity=\"0.10\"/>");
            }
        double yb = Y(baseline);
        // Filled areas above / below the baseline (clip by drawing the same polygon twice with clip rects).
        var poly = new StringBuilder();
        poly.Append(Inv, $"M{X(0):0.#} {yb:0.#}");
        for (int i = 0; i < values.Count; i++) poly.Append(Inv, $"L{X(i):0.#} {Y(values[i]):0.#}");
        poly.Append(Inv, $"L{X(values.Count - 1):0.#} {yb:0.#}Z");
        string id = "p" + Guid.NewGuid().ToString("N")[..8];
        sb.Append(Inv, $"<defs><clipPath id=\"{id}a\"><rect x=\"0\" y=\"0\" width=\"{width}\" height=\"{yb:0.#}\"/></clipPath>");
        sb.Append(Inv, $"<clipPath id=\"{id}b\"><rect x=\"0\" y=\"{yb:0.#}\" width=\"{width}\" height=\"{height - yb:0.#}\"/></clipPath></defs>");
        sb.Append(Inv, $"<path d=\"{poly}\" fill=\"{positive}\" fill-opacity=\"0.22\" clip-path=\"url(#{id}a)\"/>");
        sb.Append(Inv, $"<path d=\"{poly}\" fill=\"{negative}\" fill-opacity=\"0.22\" clip-path=\"url(#{id}b)\"/>");
        var line = new StringBuilder();
        for (int i = 0; i < values.Count; i++) line.Append(i == 0 ? 'M' : 'L').Append(Inv, $"{X(i):0.#} {Y(values[i]):0.#}");
        sb.Append(Inv, $"<path d=\"{line}\" fill=\"none\" stroke=\"#E2E8EE\" stroke-width=\"1.1\" stroke-opacity=\"0.85\"/>");
        sb.Append(Inv, $"<line x1=\"{left}\" x2=\"{right}\" y1=\"{yb:0.#}\" y2=\"{yb:0.#}\" stroke=\"#93A3B2\" stroke-dasharray=\"3 3\"/>");
        if (reference is { } r)
        {
            double yr = Y(r.Value);
            sb.Append(Inv, $"<line x1=\"{left}\" x2=\"{right}\" y1=\"{yr:0.#}\" y2=\"{yr:0.#}\" stroke=\"#F0B03F\" stroke-dasharray=\"4 3\"/>");
            sb.Append(Inv, $"<text x=\"{left + 6}\" y=\"{yr - 4:0.#}\" fill=\"#F0B03F\" class=\"t-cap\">{E(r.Label)}</text>");
        }
        for (int k = 0; k <= 5; k++)
        {
            int i = (int)Math.Round((values.Count - 1) * k / 5.0);
            sb.Append(Inv, $"<text x=\"{X(i):0.#}\" y=\"{height - 4}\" text-anchor=\"middle\" class=\"t-cap\">{i + 1}</text>");
        }
        return End(sb);
    }

    public sealed record Mark(int Start, int End, string Label, string Color, int Row, bool Tick = false, bool Reverse = false);

    /// <summary>A linear map: coordinate ruler, rows of features (bars with direction) and cut-site ticks.</summary>
    public static string LinearMap(int length, IReadOnlyList<Mark> marks, IReadOnlyList<string> rowLabels, double width = 980)
    {
        double gutter = 64, right = width - 12, rowH = 30, top = 30;
        double height = top + rowLabels.Count * rowH + 40;
        var sb = Begin(width, height, "Linear map");
        double X(double p) => gutter + (right - gutter) * p / Math.Max(1, length);
        // ruler
        sb.Append(Inv, $"<line x1=\"{gutter}\" x2=\"{right}\" y1=\"18.5\" y2=\"18.5\" stroke=\"#233140\"/>");
        foreach (double t in NiceTicks(0, length, 10))
        {
            double x = X(t);
            sb.Append(Inv, $"<line x1=\"{x:0.#}\" x2=\"{x:0.#}\" y1=\"14\" y2=\"19\" stroke=\"#5E6F7F\"/>");
            sb.Append(Inv, $"<text x=\"{x:0.#}\" y=\"10\" text-anchor=\"middle\" class=\"t-cap\">{t.ToString("N0", Inv)}</text>");
        }
        for (int r = 0; r < rowLabels.Count; r++)
            sb.Append(Inv, $"<text x=\"0\" y=\"{top + r * rowH + 12:0.#}\" class=\"t-cap\">{E(rowLabels[r])}</text>");
        foreach (var m in marks.Where(m => !m.Tick))
        {
            double y = top + m.Row * rowH + 2, x0 = X(m.Start), x1 = Math.Max(x0 + 2, X(m.End));
            sb.Append(Inv, $"<rect x=\"{x0:0.#}\" y=\"{y:0.#}\" width=\"{x1 - x0:0.#}\" height=\"10\" fill=\"{m.Color}\" fill-opacity=\"0.32\"/>");
            sb.Append(Inv, $"<rect x=\"{(m.Reverse ? x0 : x1 - 3):0.#}\" y=\"{y:0.#}\" width=\"3\" height=\"10\" fill=\"{m.Color}\"/>");
            sb.Append(Inv, $"<text x=\"{x0 + 3:0.#}\" y=\"{y + 22:0.#}\" fill=\"{m.Color}\" class=\"t-cap\">{E(m.Label)}</text>");
        }
        // cut ticks: stagger labels into lanes so neighbours do not collide
        var lanes = new List<double>();
        foreach (var m in marks.Where(m => m.Tick).OrderBy(m => m.Start))
        {
            double x = X(m.Start), labelW = 7 * m.Label.Length + 8;
            int lane = 0;
            while (lane < lanes.Count && lanes[lane] > x) lane++;
            if (lane == lanes.Count) lanes.Add(0);
            lanes[lane] = x + labelW;
            double y = top + m.Row * rowH + lane * 13;
            if (y > height - 14) continue;
            sb.Append(Inv, $"<line x1=\"{x:0.#}\" x2=\"{x:0.#}\" y1=\"{y:0.#}\" y2=\"{y + 11:0.#}\" stroke=\"{m.Color}\" stroke-width=\"1.4\"/>");
            sb.Append(Inv, $"<text x=\"{x + 3:0.#}\" y=\"{y + 9:0.#}\" fill=\"{m.Color}\" class=\"t-cap\">{E(m.Label)}</text>");
        }
        return End(sb);
    }

    public static double[] NiceTicks(double min, double max, int maxTicks)
    {
        double span = Math.Max(1e-9, max - min), raw = span / Math.Max(1, maxTicks);
        double mag = Math.Pow(10, Math.Floor(Math.Log10(raw))), norm = raw / mag;
        double step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
        var ticks = new List<double>();
        for (double t = Math.Ceiling(min / step) * step; t <= max + 1e-9; t += step) ticks.Add(t);
        return ticks.ToArray();
    }

    private static StringBuilder Begin(double w, double h, string label)
    {
        var sb = new StringBuilder();
        sb.Append(Inv, $"<svg class=\"chart\" viewBox=\"0 0 {w:0.#} {h:0.#}\" role=\"img\" aria-label=\"{label}\" preserveAspectRatio=\"xMinYMin meet\">");
        return sb;
    }

    private static string End(StringBuilder sb) => sb.Append("</svg>").ToString();

    private static string E(string s) => WebUtility.HtmlEncode(s);
}
