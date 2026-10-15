using System.Net;
using System.Text;
using Matplotlib.Net;
using Matplotlib.Net.Rendering;
using Plotly.Net;

namespace Biotech.Net.Visualization;

/// <summary>Output format for <see cref="BioPlot.Save(string, PlotFormat)"/> and notebook embedding.</summary>
public enum PlotFormat
{
    /// <summary>Raster PNG.</summary>
    Png,

    /// <summary>Vector SVG (text stays text).</summary>
    Svg,

    /// <summary>Interactive HTML (Plotly.Net) when the plot supports it, otherwise inline SVG.</summary>
    Html,

    /// <summary>Single-page PDF.</summary>
    Pdf,
}

/// <summary>
/// Base class of every Biotech.Net plot. A plot draws itself once onto a Matplotlib.Net
/// <see cref="IRenderer"/>, so the same geometry is written to PNG (SkiaSharp), SVG and PDF. Plots that have a
/// natural interactive form also return a Plotly.Net figure from <see cref="ToPlotly"/>.
/// </summary>
public abstract class BioPlot
{
    private int _width;
    private int _height;

    protected BioPlot(int width, int height)
    {
        _width = width;
        _height = height;
        BiotechColormaps.EnsureRegistered();
    }

    /// <summary>Width in pixels.</summary>
    public virtual int Width
    {
        get => _width;
        set => _width = value > 0 ? value : throw new ArgumentOutOfRangeException(nameof(value));
    }

    /// <summary>Height in pixels.</summary>
    public virtual int Height
    {
        get => _height;
        set => _height = value > 0 ? value : throw new ArgumentOutOfRangeException(nameof(value));
    }

    /// <summary>Plot title (null for none).</summary>
    public string? Title { get; set; }

    /// <summary>Colours and fonts; <see cref="BioTheme.Default"/> unless set.</summary>
    public BioTheme Theme { get; set; } = BioTheme.Default;

    /// <summary>Draws the plot in pixel coordinates (origin top-left, <see cref="Width"/> × <see cref="Height"/>).</summary>
    public abstract void Render(IRenderer renderer);

    /// <summary>Renders a PNG. <paramref name="scale"/> 2 gives a sharp image for high-DPI screens.</summary>
    public virtual byte[] ToPng(double scale = 1)
    {
        int w = Math.Max(1, (int)Math.Round(Width * scale)), h = Math.Max(1, (int)Math.Round(Height * scale));
        using var skia = SkiaRenderer.CreateRaster(w, h);
        Render(scale == 1 ? skia : new ScaledRenderer(skia, scale));
        skia.Flush();
        return skia.Encode(ImageFormat.Png);
    }

    /// <summary>Renders an SVG document.</summary>
    public virtual string ToSvg()
    {
        using var svg = new SvgRenderer(Width, Height);
        Render(svg);
        svg.Flush();
        return svg.ToSvg();
    }

    /// <summary>Renders a single-page PDF.</summary>
    public virtual byte[] ToPdf()
    {
        using var stream = new MemoryStream();
        using (var pdf = SkiaRenderer.CreatePdf(stream, Width, Height))
        {
            Render(pdf);
            pdf.Flush();
            pdf.ClosePdf();
        }
        return stream.ToArray();
    }

    /// <summary>PNG as a <c>data:image/png;base64,...</c> URI.</summary>
    public string ToPngDataUri(double scale = 1) => "data:image/png;base64," + Convert.ToBase64String(ToPng(scale));

    /// <summary>The interactive Plotly.Net version of this plot, or null when it has none.</summary>
    public virtual PlotlyFigure? ToPlotly() => null;

    /// <summary>Whether <see cref="ToPlotly"/> returns a figure.</summary>
    public virtual bool SupportsInteractive => false;

    /// <summary>
    /// An HTML fragment: the Plotly chart when <paramref name="interactive"/> and supported, otherwise inline SVG.
    /// </summary>
    public virtual string ToHtml(bool interactive = true)
    {
        if (interactive && ToPlotly() is { } figure) return figure.ToNotebookHtml();
        var svg = ToSvg();
        int start = svg.IndexOf("<svg", StringComparison.Ordinal);
        return $"<div class=\"biotech-plot\" style=\"max-width:100%;overflow:auto\">{(start >= 0 ? svg[start..] : svg)}</div>";
    }

    /// <summary>A complete standalone HTML page containing the plot.</summary>
    public string ToHtmlDocument(bool interactive = true)
    {
        if (interactive && ToPlotly() is { } figure) return figure.ToHtml();
        var sb = new StringBuilder();
        sb.Append("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">");
        sb.Append("<title>").Append(WebUtility.HtmlEncode(Title ?? GetType().Name)).Append("</title>");
        sb.Append("<style>body{margin:0;padding:16px;background:").Append(Theme.Background.ToHex()).Append(";color:")
          .Append(Theme.Foreground.ToHex()).Append(";font-family:'Segoe UI',system-ui,sans-serif}</style></head><body>");
        sb.Append(ToHtml(interactive: false));
        sb.Append("</body></html>");
        return sb.ToString();
    }

    /// <summary>Saves to a file; the format follows the extension (.png, .svg, .html/.htm, .pdf).</summary>
    public void Save(string path)
    {
        string ext = Path.GetExtension(path).ToLowerInvariant();
        switch (ext)
        {
            case ".png": File.WriteAllBytes(path, ToPng()); break;
            case ".svg": File.WriteAllText(path, ToSvg()); break;
            case ".html" or ".htm": File.WriteAllText(path, ToHtmlDocument()); break;
            case ".pdf": File.WriteAllBytes(path, ToPdf()); break;
            default: throw new NotSupportedException($"Unsupported plot format '{ext}'. Use .png, .svg, .html or .pdf.");
        }
    }

    /// <summary>Saves in an explicit format.</summary>
    public void Save(string path, PlotFormat format)
    {
        switch (format)
        {
            case PlotFormat.Png: File.WriteAllBytes(path, ToPng()); break;
            case PlotFormat.Svg: File.WriteAllText(path, ToSvg()); break;
            case PlotFormat.Html: File.WriteAllText(path, ToHtmlDocument()); break;
            case PlotFormat.Pdf: File.WriteAllBytes(path, ToPdf()); break;
        }
    }

    /// <summary>Draws the title centred at the top and returns the height it used.</summary>
    protected double DrawTitle(IRenderer r, double top = 10)
    {
        if (string.IsNullOrEmpty(Title)) return 0;
        var font = Theme.TitleFont;
        string title = Paint.Fit(r, Title, font, Width - 24);
        r.DrawText(title, new Pt(Width / 2.0, top), font, Theme.Foreground, TextAlign.Center, TextBaseline.Top);
        return r.MeasureText(Title, font).Height + top + 6;
    }

    /// <summary>A Plotly figure with this theme's template, size and title.</summary>
    protected PlotlyFigure NewPlotly()
    {
        var f = new PlotlyFigure();
        f.WithTemplate(Theme.PlotlyTemplate).WithSize(Width, Height);
        if (Title is not null) f.WithTitle(Title);
        f.Layout.PaperBackground = Theme.Background;
        f.Layout.PlotBackground = Theme.Panel;
        return f;
    }

    /// <inheritdoc />
    public override string ToString() => $"{GetType().Name} {Width}×{Height}{(Title is null ? "" : $" \"{Title}\"")}";
}

/// <summary>
/// A plot built as a Matplotlib.Net <see cref="Figure"/> (axes, ticks, legend, colour bar). Call
/// <see cref="ToFigure"/> to get the figure and adjust it with the Matplotlib.Net API before saving.
/// </summary>
public abstract class FigurePlot : BioPlot
{
    protected FigurePlot(int width, int height) : base(width, height) { }

    /// <summary>Axis label for x (null keeps the plot's default).</summary>
    public string? XLabel { get; set; }

    /// <summary>Axis label for y (null keeps the plot's default).</summary>
    public string? YLabel { get; set; }

    /// <summary>Builds a fresh figure for the current settings.</summary>
    public Figure ToFigure()
    {
        var fig = new Figure(Width / 100.0, Height / 100.0, 100, Theme.ToPlotStyle());
        Build(fig);
        return fig;
    }

    /// <summary>Adds the plot's axes and artists to <paramref name="figure"/>.</summary>
    protected abstract void Build(Figure figure);

    /// <summary>A constant marker-size array for Matplotlib.Net scatter calls.</summary>
    protected static double[] Sizes(int count, double size)
    {
        var s = new double[count];
        Array.Fill(s, size);
        return s;
    }

    /// <inheritdoc />
    public override void Render(IRenderer renderer) => ToFigure().Draw(renderer);

    /// <inheritdoc />
    public override byte[] ToPng(double scale = 1) => ToFigure().ToPng(100 * scale);

    /// <inheritdoc />
    public override string ToSvg() => ToFigure().ToSvg();

    /// <inheritdoc />
    public override byte[] ToPdf() => ToFigure().ToPdf();

}
