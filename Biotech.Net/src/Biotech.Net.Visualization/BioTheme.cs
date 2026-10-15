using Matplotlib.Net.Rendering;
using Matplotlib.Net.Styling;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>
/// The shared look of every Biotech.Net plot: background and ink colours, a colour-blind-aware categorical
/// palette, regulation colours (up / down / not significant), nucleotide and amino-acid colours, fonts, and the
/// matching Matplotlib.Net <see cref="PlotStyle"/> and Plotly template. Use <see cref="Light"/> or
/// <see cref="Dark"/>, or <see cref="With"/> to derive a variant.
/// </summary>
public sealed record BioTheme
{
    /// <summary>Theme name.</summary>
    public required string Name { get; init; }

    /// <summary>True for dark backgrounds.</summary>
    public bool IsDark { get; init; }

    /// <summary>Page background.</summary>
    public Color Background { get; init; }

    /// <summary>Panel / plot-area background.</summary>
    public Color Panel { get; init; }

    /// <summary>Primary ink (titles, labels).</summary>
    public Color Foreground { get; init; }

    /// <summary>Secondary ink (tick labels, captions).</summary>
    public Color Muted { get; init; }

    /// <summary>Grid lines.</summary>
    public Color Grid { get; init; }

    /// <summary>Axis lines and panel borders.</summary>
    public Color Border { get; init; }

    /// <summary>Accent colour (gene models, single-series plots).</summary>
    public Color Accent { get; init; }

    /// <summary>Categorical palette for groups and series.</summary>
    public required IReadOnlyList<Color> Palette { get; init; }

    /// <summary>Up-regulated / positive.</summary>
    public Color Up { get; init; }

    /// <summary>Down-regulated / negative.</summary>
    public Color Down { get; init; }

    /// <summary>Not significant / background points.</summary>
    public Color Neutral { get; init; }

    /// <summary>Colours of A, C, G, T/U and other symbols.</summary>
    public Color BaseA { get; init; }
    public Color BaseC { get; init; }
    public Color BaseG { get; init; }
    public Color BaseT { get; init; }
    public Color BaseOther { get; init; }

    /// <summary>Font family for all text.</summary>
    public string FontFamily { get; init; } = "Segoe UI";

    /// <summary>Base font size in pixels.</summary>
    public double FontSize { get; init; } = 12;

    /// <summary>Name of the diverging colormap (z-scores, fold changes).</summary>
    public string DivergingColormap { get; init; } = BiotechColormaps.Diverging;

    /// <summary>Name of the sequential colormap (expression, distance, density).</summary>
    public string SequentialColormap { get; init; } = "viridis";

    /// <summary>Plotly layout template.</summary>
    public string PlotlyTemplate { get; init; } = "plotly_white";

    /// <summary>The light theme (default).</summary>
    public static BioTheme Light { get; } = new()
    {
        Name = "light",
        Background = Hex("#FFFFFF"),
        Panel = Hex("#FFFFFF"),
        Foreground = Hex("#1C2430"),
        Muted = Hex("#5B6573"),
        Grid = Hex("#E6E9EE"),
        Border = Hex("#C3C9D2"),
        Accent = Hex("#1F5FAD"),
        Palette = HexList("#2F6FDB", "#E8833A", "#2BA58A", "#D64F6E", "#8A63D2", "#A07C2F", "#3FA7D6", "#6B7280", "#B5A00E", "#1E4E79"),
        Up = Hex("#D1414A"),
        Down = Hex("#2F6FDB"),
        Neutral = Hex("#BAC1CB"),
        BaseA = Hex("#1E9E5A"),
        BaseC = Hex("#2F6FDB"),
        BaseG = Hex("#EBA300"),
        BaseT = Hex("#D1414A"),
        BaseOther = Hex("#8A93A0"),
        PlotlyTemplate = "plotly_white",
    };

    /// <summary>The dark theme.</summary>
    public static BioTheme Dark { get; } = new()
    {
        Name = "dark",
        IsDark = true,
        Background = Hex("#11161D"),
        Panel = Hex("#161C25"),
        Foreground = Hex("#E6EAF0"),
        Muted = Hex("#9AA4B2"),
        Grid = Hex("#26303C"),
        Border = Hex("#3A4553"),
        Accent = Hex("#6EA8FF"),
        Palette = HexList("#6EA8FF", "#FFA45C", "#4FD1A8", "#FF7A98", "#B596FF", "#D9B25F", "#6CCDF2", "#A3ACB9", "#E3D04A", "#5C8FD6"),
        Up = Hex("#FF6B6B"),
        Down = Hex("#5B9BFF"),
        Neutral = Hex("#4A5462"),
        BaseA = Hex("#3FCB7E"),
        BaseC = Hex("#5B9BFF"),
        BaseG = Hex("#FFC233"),
        BaseT = Hex("#FF6B6B"),
        BaseOther = Hex("#8A93A0"),
        PlotlyTemplate = "plotly_dark",
    };

    /// <summary>Theme used when a plot does not set one. Defaults to <see cref="Light"/>.</summary>
    public static BioTheme Default { get; set; } = Light;

    /// <summary>A copy with changes applied (convenience for <c>with</c> from other .NET languages).</summary>
    public BioTheme With(Func<BioTheme, BioTheme> change) => change(this);

    /// <summary>Palette colour <paramref name="index"/> (cycled).</summary>
    public Color CategoryColor(int index) => Palette[((index % Palette.Count) + Palette.Count) % Palette.Count];

    /// <summary>Colour of a nucleotide symbol.</summary>
    public Color NucleotideColor(char symbol) => char.ToUpperInvariant(symbol) switch
    {
        'A' => BaseA,
        'C' => BaseC,
        'G' => BaseG,
        'T' or 'U' => BaseT,
        _ => BaseOther,
    };

    /// <summary>
    /// Colour of an amino acid by physico-chemical class (hydrophobic, positive, negative, polar, glycine, proline,
    /// aromatic-polar), close to the Clustal X scheme.
    /// </summary>
    public Color AminoAcidColor(char residue) => char.ToUpperInvariant(residue) switch
    {
        'A' or 'V' or 'L' or 'I' or 'M' or 'F' or 'W' or 'C' => IsDark ? Hex("#6EA8FF") : Hex("#2F6FDB"),
        'K' or 'R' => IsDark ? Hex("#FF6B6B") : Hex("#D1414A"),
        'D' or 'E' => IsDark ? Hex("#E58BFF") : Hex("#B23CC8"),
        'S' or 'T' or 'N' or 'Q' => IsDark ? Hex("#4FD1A8") : Hex("#1E9E5A"),
        'G' => IsDark ? Hex("#FFA45C") : Hex("#E07A1F"),
        'P' => IsDark ? Hex("#E3D04A") : Hex("#B59A00"),
        'H' or 'Y' => IsDark ? Hex("#6CCDF2") : Hex("#1597A5"),
        _ => BaseOther,
    };

    /// <summary>Colour of a residue: nucleotide scheme for nucleotide alphabets, amino-acid scheme otherwise.</summary>
    public Color ResidueColor(char symbol, bool nucleotide) => nucleotide ? NucleotideColor(symbol) : AminoAcidColor(symbol);

    /// <summary>CPK-style colour of a chemical element (carbon follows the foreground ink).</summary>
    public Color ElementColor(string symbol) => symbol switch
    {
        "C" => Foreground,
        "H" => Muted,
        "N" => IsDark ? Hex("#7FA6FF") : Hex("#2B50C8"),
        "O" => IsDark ? Hex("#FF6B6B") : Hex("#D02A2A"),
        "S" => IsDark ? Hex("#F2D04A") : Hex("#B8960A"),
        "P" => Hex("#E8833A"),
        "F" or "Cl" => IsDark ? Hex("#4FD18A") : Hex("#1E9E3A"),
        "Br" => Hex("#A0522D"),
        "I" => Hex("#8A2BE2"),
        "B" => Hex("#E0868A"),
        _ => IsDark ? Hex("#E58BFF") : Hex("#9B3FB0"),
    };

    /// <summary>Regular text font.</summary>
    public Font Font(double? size = null, bool bold = false)
        => new(FontFamily, size ?? FontSize, bold ? FontWeight.Bold : FontWeight.Normal);

    /// <summary>Title font.</summary>
    public Font TitleFont => Font(FontSize * 1.3, bold: true);

    /// <summary>Small font for tick labels and annotations.</summary>
    public Font SmallFont => Font(FontSize * 0.85);

    /// <summary>A Matplotlib.Net style matching this theme (white-grid base, no top/right spines).</summary>
    public PlotStyle ToPlotStyle()
    {
        var s = (IsDark ? PlotStyle.DarkBackground() : PlotStyle.SeabornWhiteGrid()).Clone();
        s.FigureFaceColor = Background;
        s.AxesFaceColor = Panel;
        s.AxesEdgeColor = Border;
        s.AxesLineWidth = 1;
        // Spines are (Left, Bottom, Top, Right).
        s.Spines = (true, true, false, false);
        s.GridVisible = true;
        s.GridColor = Grid;
        s.GridLineWidth = 0.8;
        s.GridLineStyle = LineStyle.Solid;
        s.GridAlpha = 1;
        s.TickColor = Border;
        s.TickLabelColor = Muted;
        s.TextColor = Foreground;
        s.FontFamily = FontFamily;
        s.FontSize = FontSize;
        s.TitleWeight = FontWeight.Bold;
        s.ColorCycle = Palette;
        s.Colormap = SequentialColormap;
        s.LegendFaceColor = Panel;
        s.LegendEdgeColor = Grid;
        s.LegendFrame = true;
        s.PatchEdgeColor = Panel;
        return s;
    }

    internal static Color Hex(string hex) => Color.Parse(hex);

    private static IReadOnlyList<Color> HexList(params string[] hex) => hex.Select(Color.Parse).ToArray();
}

/// <summary>Colormaps registered by Biotech.Net with Visual.Net's <see cref="Colormaps"/> registry.</summary>
public static class BiotechColormaps
{
    /// <summary>Blue–white–red diverging map for z-scores and fold changes.</summary>
    public const string Diverging = "biotech_diverging";

    /// <summary>White-to-accent-blue sequential map for densities and contact maps.</summary>
    public const string Blues = "biotech_blues";

    /// <summary>N-to-C rainbow used for protein chains (blue → green → yellow → red).</summary>
    public const string Spectrum = "biotech_spectrum";

    private static int _registered;

    /// <summary>Registers the maps (idempotent). Called automatically by every plot.</summary>
    public static void EnsureRegistered()
    {
        if (Interlocked.Exchange(ref _registered, 1) == 1) return;
        Colormaps.Register(LinearSegmentedColormap.FromHex(Diverging, ["#1E4E9C", "#5B8FD6", "#BFD3EE", "#F7F7F7", "#F2C4C0", "#D9665F", "#9E1B22"]));
        Colormaps.Register(LinearSegmentedColormap.FromHex(Blues, ["#F7FAFD", "#C9DCF2", "#7FAADF", "#2F6FDB", "#123E85"]));
        Colormaps.Register(LinearSegmentedColormap.FromHex(Spectrum, ["#2F3DDB", "#2F9FDB", "#2BC48A", "#B8D432", "#F2B51C", "#E05A2A", "#C8202F"]));
    }
}
