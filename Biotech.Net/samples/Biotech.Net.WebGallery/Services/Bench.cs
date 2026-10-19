using System.Globalization;
using System.Text;
using Biotech.Net.Visualization;
using Microsoft.AspNetCore.Components;

namespace Biotech.Net.WebGallery.Services;

/// <summary>
/// The "bench instrument" design tokens shared with the desktop Gallery (samples/Biotech.Net.Gallery/Controls/Palette.cs),
/// the matching Biotech.Net.Visualization theme, and safe SVG embedding of library figures.
/// </summary>
public static class Bench
{
    public const string DyeA = "#3FD69B", DyeC = "#4EA3FF", DyeG = "#F0B03F", DyeT = "#FF6363";
    public const string Ink = "#E2E8EE", Ink2 = "#93A3B2", Ink3 = "#5E6F7F", Rule = "#233140";

    public static string Channel(char c) => char.ToUpperInvariant(c) switch
    {
        'A' => DyeA,
        'C' => DyeC,
        'G' => DyeG,
        'T' or 'U' => DyeT,
        _ => Ink3,
    };

    /// <summary>CSS class of a dye channel (a, c, g, t).</summary>
    public static string ChannelClass(char c) => char.ToUpperInvariant(c) switch
    {
        'A' => "a",
        'C' => "c",
        'G' => "g",
        'T' or 'U' => "t",
        _ => "n",
    };

    private const string BenchColormap = "gallery_bench";
    private static BioTheme? _theme;

    /// <summary>Plot theme matching the chassis; used for every library-rendered figure.</summary>
    public static BioTheme PlotTheme => _theme ??= CreateTheme();

    private static BioTheme CreateTheme()
    {
        BiotechColormaps.EnsureRegistered();
        Visual.Net.Core.Colors.Colormaps.Register(Visual.Net.Core.Colors.LinearSegmentedColormap.FromHex(BenchColormap,
            ["#FFE2A6", "#F0B03F", "#C46F3A", "#3C6FB0", "#1D3550", "#111921"]));
        static Visual.Net.Core.Colors.Color C(string hex) => Visual.Net.Core.Colors.Color.Parse(hex);
        return BioTheme.Dark with
        {
            Name = "webgallery-bench",
            Background = C("#111921"),
            Panel = C("#111921"),
            Foreground = C(Ink),
            Muted = C(Ink2),
            Grid = C("#1C2833"),
            Border = C("#2A3A4A"),
            Accent = C(DyeC),
            BaseA = C(DyeA),
            BaseC = C(DyeC),
            BaseG = C(DyeG),
            BaseT = C(DyeT),
            Up = C(DyeT),
            Down = C(DyeC),
            FontFamily = "Inter",
            FontSize = 12,
            SequentialColormap = BenchColormap,
        };
    }

    /// <summary>
    /// Renders a library figure to inline SVG that scales to its container. A failure (for example a build without
    /// the native SkiaSharp text metrics) becomes a readable notice instead of an exception.
    /// </summary>
    public static MarkupString Svg(BioPlot plot)
    {
        try
        {
            string svg = plot.ToSvg();
            int start = svg.IndexOf("<svg", StringComparison.Ordinal);
            if (start > 0) svg = svg[start..];
            return new MarkupString(Responsive(svg, plot.Width, plot.Height));
        }
        catch (Exception ex)
        {
            return Notice("Figure unavailable in this build", Explain(ex));
        }
    }

    /// <summary>Makes an SVG document scale: adds a viewBox when missing and lets CSS size it.</summary>
    public static string Responsive(string svg, double width, double height)
    {
        int end = svg.IndexOf('>');
        if (end < 0) return svg;
        string head = svg[..end];
        if (!head.Contains("viewBox", StringComparison.Ordinal))
            head += string.Create(CultureInfo.InvariantCulture, $" viewBox=\"0 0 {width:0.##} {height:0.##}\"");
        head += " class=\"figure\" preserveAspectRatio=\"xMinYMin meet\"";
        return head + svg[end..];
    }

    public static string Explain(Exception ex)
    {
        var root = ex;
        while (root.InnerException is not null) root = root.InnerException;
        bool skia = root is DllNotFoundException or EntryPointNotFoundException || root.Message.Contains("SkiaSharp", StringComparison.OrdinalIgnoreCase);
        return skia
            ? "Figures need SkiaSharp's native library linked into the WebAssembly runtime. Publish with the wasm-tools workload installed (dotnet workload install wasm-tools)."
            : root.GetType().Name + ": " + root.Message;
    }

    public static MarkupString Notice(string title, string body)
        => new($"<div class=\"notice\"><b>{Enc(title)}</b><span>{Enc(body)}</span></div>");

    public static string Enc(string s) => System.Net.WebUtility.HtmlEncode(s);

    public static string F(double v, string format) => double.IsNaN(v) ? "—" : v.ToString(format, CultureInfo.InvariantCulture);

    public static string N(long v) => v.ToString("N0", CultureInfo.InvariantCulture);

    /// <summary>
    /// The signature: a Sanger-style four-channel trace of a DNA sequence, one Gaussian peak per base in its dye colour,
    /// base calls above the peaks. Same algorithm as the desktop Gallery's Chromatogram control.
    /// </summary>
    public static string Chromatogram(string? sequence, int maxBases = 150)
    {
        string seq = new((sequence ?? "").ToUpperInvariant().Where(char.IsLetter).Take(maxBases).ToArray());
        const double sp = 12.5, callBand = 14, h = 62, baseY = h - 2, top = callBand + 2;
        double amp = baseY - top;
        int n = seq.Length;
        double w = Math.Max(1, n) * sp;
        var sb = new StringBuilder();
        sb.Append(CultureInfo.InvariantCulture, $"<svg class=\"trace\" viewBox=\"0 0 {w:0.#} {h}\" preserveAspectRatio=\"xMinYMid slice\" aria-hidden=\"true\">");
        sb.Append(CultureInfo.InvariantCulture, $"<line x1=\"0\" y1=\"{baseY + 0.5}\" x2=\"{w:0.#}\" y2=\"{baseY + 0.5}\" stroke=\"{Rule}\" stroke-width=\"1\"/>");
        if (n > 0)
        {
            var heights = new double[n];
            for (int i = 0; i < n; i++)
            {
                uint hsh = (uint)(i * 2654435761u) ^ seq[i] * 40503u;
                heights[i] = 0.55 + 0.45 * ((hsh >> 8) % 1000) / 1000.0;
            }
            double sigma = sp * 0.40;
            int samples = (int)(n * sp);
            sb.Append("<g class=\"peaks\" fill=\"none\" stroke-width=\"1.3\" stroke-linejoin=\"round\">");
            foreach (char ch in "ACGT")
            {
                sb.Append(CultureInfo.InvariantCulture, $"<path stroke=\"{Channel(ch)}\" stroke-opacity=\"0.95\" d=\"");
                for (int px = 0; px <= samples; px += 2)
                {
                    double y = 0;
                    int center = (int)(px / sp);
                    for (int k = Math.Max(0, center - 2); k <= Math.Min(n - 1, center + 2); k++)
                    {
                        double d = (px - (k * sp + sp / 2)) / sigma;
                        y += (seq[k] == ch ? heights[k] : 0.06 * heights[k]) * Math.Exp(-0.5 * d * d);
                    }
                    sb.Append(px == 0 ? 'M' : 'L').Append(CultureInfo.InvariantCulture, $"{px} {baseY - Math.Min(1.05, y) * amp:0.#}");
                }
                sb.Append("\"/>");
            }
            sb.Append("</g><g class=\"calls\" font-size=\"9.5\" text-anchor=\"middle\">");
            for (int i = 0; i < n; i++)
                sb.Append(CultureInfo.InvariantCulture, $"<text x=\"{i * sp + sp / 2:0.##}\" y=\"10\" fill=\"{Channel(seq[i])}\">{seq[i]}</text>");
            sb.Append("</g>");
        }
        sb.Append("</svg>");
        return sb.ToString();
    }
}
