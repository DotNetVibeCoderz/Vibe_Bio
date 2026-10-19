using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Visualization;
using Biotech.Net.WebGallery.Services;

namespace Biotech.Net.WebGallery.Playground;

/// <summary>One rich output of a snippet: inline SVG, a sandboxed HTML document (3-D viewers) or plain text.</summary>
public sealed record SnippetOutput(string Kind, string Content, string Caption);

/// <summary>
/// <c>Show(...)</c> for Playground snippets (imported with <c>using static</c>): figures render as SVG under the console,
/// structures and molecules get a 3Dmol.js viewer, anything else is printed.
/// </summary>
public static class Display
{
    internal static List<SnippetOutput>? Sink;

    /// <summary>Renders <paramref name="value"/> below the console output.</summary>
    public static void Show(object? value, string? caption = null)
    {
        var sink = Sink;
        if (sink is null) { Console.WriteLine(value); return; }
        switch (value)
        {
            case null:
                sink.Add(new("text", "null", caption ?? ""));
                break;
            case BioPlot plot:
                string svg = plot.ToSvg();
                int start = svg.IndexOf("<svg", StringComparison.Ordinal);
                sink.Add(new("svg", Bench.Responsive(start > 0 ? svg[start..] : svg, plot.Width, plot.Height), caption ?? plot.ToString()));
                break;
            case Molecule mol:
                Show(new MoleculeViewer(mol) { Theme = Bench.PlotTheme }, caption ?? mol.Name);
                break;
            case Structure structure:
                sink.Add(new("html", HtmlDocument(ProteinViewer.To3DHtml(structure, new Viewer3DOptions { Width = 720, Height = 420, Theme = Bench.PlotTheme })),
                    caption ?? structure.Id));
                break;
            case string s when s.TrimStart().StartsWith("<svg", StringComparison.OrdinalIgnoreCase):
                sink.Add(new("svg", s, caption ?? "SVG"));
                break;
            default:
                sink.Add(new("text", value.ToString() ?? "", caption ?? value.GetType().Name));
                break;
        }
    }

    /// <summary>A full HTML page for an iframe's srcdoc (the library's 3Dmol.js fragments include their own script).</summary>
    public static string HtmlDocument(string fragment)
        => "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><style>html,body{margin:0;background:#111921;color:#E2E8EE;font-family:Inter,system-ui,sans-serif}" +
           ".biotech-3d>div{border:0!important;border-radius:0!important;width:100%!important}</style></head><body>" + fragment + "</body></html>";
}
