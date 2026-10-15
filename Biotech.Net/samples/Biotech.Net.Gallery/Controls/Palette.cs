using Avalonia.Media;
using Biotech.Net.Visualization;

namespace Biotech.Net.Gallery.Controls;

/// <summary>The Gallery's design tokens in code (mirrors App.axaml) and the matching plot theme.</summary>
public static class Palette
{
    public static readonly Color Bench = Color.Parse("#0B1015");
    public static readonly Color Deck = Color.Parse("#111921");
    public static readonly Color Raised = Color.Parse("#17212B");
    public static readonly Color Rule = Color.Parse("#233140");
    public static readonly Color Grid = Color.Parse("#1A2530");
    public static readonly Color Ink = Color.Parse("#E2E8EE");
    public static readonly Color Ink2 = Color.Parse("#93A3B2");
    public static readonly Color Ink3 = Color.Parse("#5E6F7F");
    public static readonly Color DyeA = Color.Parse("#3FD69B");
    public static readonly Color DyeC = Color.Parse("#4EA3FF");
    public static readonly Color DyeG = Color.Parse("#F0B03F");
    public static readonly Color DyeT = Color.Parse("#FF6363");

    public static IBrush Brush(Color c) => new SolidColorBrush(c);

    public static readonly IBrush InkBrush = Brush(Ink);
    public static readonly IBrush Ink2Brush = Brush(Ink2);
    public static readonly IBrush Ink3Brush = Brush(Ink3);
    public static readonly IBrush RuleBrush = Brush(Rule);
    public static readonly IBrush LampBrush = Brush(DyeG);

    public static Color Channel(char c) => char.ToUpperInvariant(c) switch
    {
        'A' => DyeA,
        'C' => DyeC,
        'G' => DyeG,
        'T' or 'U' => DyeT,
        _ => Ink3,
    };

    public static readonly FontFamily Display = new("Bahnschrift, Inter, Segoe UI");
    public static readonly FontFamily Mono = new("Cascadia Mono, Consolas, Menlo, monospace");
    public static readonly FontFamily Body = new("Inter, Segoe UI");

    /// <summary>Biotech.Net.Visualization theme matching the chassis, used for every library-rendered figure.</summary>
    private static BioTheme? _plotTheme;
    public static BioTheme PlotTheme => _plotTheme ??= CreatePlotTheme();

    /// <summary>Sequential map from lamp amber (near) through dye blue to the chassis (far).</summary>
    public const string BenchColormap = "gallery_bench";

    private static BioTheme CreatePlotTheme()
    {
        BiotechColormaps.EnsureRegistered();
        Visual.Net.Core.Colors.Colormaps.Register(Visual.Net.Core.Colors.LinearSegmentedColormap.FromHex(BenchColormap,
            ["#FFE2A6", "#F0B03F", "#C46F3A", "#3C6FB0", "#1D3550", "#111921"]));
        return BaseTheme;
    }

    private static readonly BioTheme BaseTheme = BioTheme.Dark with
    {
        Name = "gallery-bench",
        Background = Visual.Net.Core.Colors.Color.Parse("#111921"),
        Panel = Visual.Net.Core.Colors.Color.Parse("#111921"),
        Foreground = Visual.Net.Core.Colors.Color.Parse("#E2E8EE"),
        Muted = Visual.Net.Core.Colors.Color.Parse("#93A3B2"),
        Grid = Visual.Net.Core.Colors.Color.Parse("#1C2833"),
        Border = Visual.Net.Core.Colors.Color.Parse("#2A3A4A"),
        Accent = Visual.Net.Core.Colors.Color.Parse("#4EA3FF"),
        BaseA = Visual.Net.Core.Colors.Color.Parse("#3FD69B"),
        BaseC = Visual.Net.Core.Colors.Color.Parse("#4EA3FF"),
        BaseG = Visual.Net.Core.Colors.Color.Parse("#F0B03F"),
        BaseT = Visual.Net.Core.Colors.Color.Parse("#FF6363"),
        Up = Visual.Net.Core.Colors.Color.Parse("#FF6363"),
        Down = Visual.Net.Core.Colors.Color.Parse("#4EA3FF"),
        FontFamily = "Segoe UI",
        FontSize = 12,
        SequentialColormap = BenchColormap,
    };
}
