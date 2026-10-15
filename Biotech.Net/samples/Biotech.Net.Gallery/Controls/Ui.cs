using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Primitives;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using Biotech.Net.Visualization;

namespace Biotech.Net.Gallery.Controls;

/// <summary>Small builders so code-built pages share the type scale and surfaces defined in App.axaml.</summary>
public static class Ui
{
    public static TextBlock Text(string text, string cls = "body")
    {
        var t = new TextBlock { Text = text };
        foreach (var c in cls.Split(' ', StringSplitOptions.RemoveEmptyEntries)) t.Classes.Add(c);
        return t;
    }

    public static TextBlock Eyebrow(string text) => Text(text.ToUpperInvariant(), "eyebrow");

    public static SelectableTextBlock Mono(string text, double size = 12, IBrush? fg = null)
    {
        var t = new SelectableTextBlock { Text = text, FontSize = size, TextWrapping = TextWrapping.Wrap, Foreground = fg ?? Palette.InkBrush };
        t.Classes.Add("mono");
        return t;
    }

    /// <summary>A titled panel (Border.panel) with an optional right-aligned note.</summary>
    public static Border Panel(string heading, Control content, string? note = null, Control? headerExtra = null)
    {
        var head = new Grid { ColumnDefinitions = new ColumnDefinitions("*,Auto"), Margin = new Thickness(0, 0, 0, 10) };
        head.Children.Add(Text(heading.ToUpperInvariant(), "h2"));
        if (headerExtra is not null) { Grid.SetColumn(headerExtra, 1); head.Children.Add(headerExtra); }
        else if (note is not null)
        {
            var n = Text(note, "unit");
            n.VerticalAlignment = VerticalAlignment.Center;
            Grid.SetColumn(n, 1);
            head.Children.Add(n);
        }
        var dock = new DockPanel();
        DockPanel.SetDock(head, Dock.Top);
        dock.Children.Add(head);
        dock.Children.Add(content);
        var b = new Border { Child = dock };
        b.Classes.Add("panel");
        return b;
    }

    /// <summary>An instrument readout: small label, big monospaced value, unit.</summary>
    public static Control Readout(string label, string value, string unit = "", Color? accent = null)
    {
        var sp = new StackPanel { Spacing = 2 };
        sp.Children.Add(Text(label.ToUpperInvariant(), "eyebrow"));
        var row = new StackPanel { Orientation = Orientation.Horizontal, Spacing = 5 };
        var v = Text(value, "readout");
        if (accent is { } a) v.Foreground = new SolidColorBrush(a);
        row.Children.Add(v);
        if (unit.Length > 0)
        {
            var u = Text(unit, "unit");
            u.VerticalAlignment = VerticalAlignment.Bottom;
            u.Margin = new Thickness(0, 0, 0, 4);
            row.Children.Add(u);
        }
        sp.Children.Add(row);
        return sp;
    }

    /// <summary>A row of readouts separated by hairlines.</summary>
    public static Control ReadoutRow(params Control[] readouts)
    {
        var g = new Grid();
        for (int i = 0; i < readouts.Length; i++)
        {
            g.ColumnDefinitions.Add(new ColumnDefinition(GridLength.Star));
            var cell = new Border
            {
                Child = readouts[i],
                Padding = new Thickness(i == 0 ? 0 : 16, 2, 8, 2),
                BorderBrush = Palette.RuleBrush,
                BorderThickness = new Thickness(i == 0 ? 0 : 1, 0, 0, 0),
            };
            Grid.SetColumn(cell, i);
            g.Children.Add(cell);
        }
        return g;
    }

    /// <summary>A key/value list in two columns (label left, mono value right).</summary>
    public static Control KeyValues(IEnumerable<(string Key, string Value, Color? Color)> rows)
    {
        var g = new Grid { ColumnDefinitions = new ColumnDefinitions("Auto,*") };
        int r = 0;
        foreach (var (k, v, c) in rows)
        {
            g.RowDefinitions.Add(new RowDefinition(GridLength.Auto));
            var kt = Text(k, "caption");
            kt.FontSize = 12;
            kt.Foreground = Palette.Ink2Brush;
            kt.Margin = new Thickness(0, 4);
            var vt = Text(v, "td");
            vt.HorizontalAlignment = HorizontalAlignment.Right;
            vt.TextAlignment = TextAlignment.Right;
            vt.TextTrimming = TextTrimming.CharacterEllipsis;
            vt.Margin = new Thickness(16, 0, 0, 0);
            if (c is { } col) vt.Foreground = new SolidColorBrush(col);
            Grid.SetRow(kt, r); Grid.SetRow(vt, r); Grid.SetColumn(vt, 1);
            var line = new Border { Height = 1, Background = new SolidColorBrush(Palette.Grid), VerticalAlignment = VerticalAlignment.Bottom };
            Grid.SetRow(line, r); Grid.SetColumnSpan(line, 2);
            g.Children.Add(line); g.Children.Add(kt); g.Children.Add(vt);
            r++;
        }
        return g;
    }

    /// <summary>Cell content for <see cref="Table"/>: plain text, or text in a colour.</summary>
    public readonly record struct Cell(string Text, Color? Color = null, bool Bar = false, double BarValue = 0);

    /// <summary>A hairline table. Column widths use Grid syntax ("Auto", "*", "2*", "80").</summary>
    public static Control Table(string[] headers, string[] widths, IEnumerable<Cell[]> rows, int maxRows = 200)
    {
        var g = new Grid();
        foreach (var w in widths) g.ColumnDefinitions.Add(new ColumnDefinition(GridLength.Parse(w)));
        g.RowDefinitions.Add(new RowDefinition(GridLength.Auto));
        for (int c = 0; c < headers.Length; c++)
        {
            var t = Text(headers[c].ToUpperInvariant(), "th");
            t.Margin = new Thickness(c == 0 ? 0 : 10, 0, 6, 6);
            Grid.SetColumn(t, c);
            g.Children.Add(t);
        }
        var hl = new Border { Height = 1, Background = Palette.RuleBrush, VerticalAlignment = VerticalAlignment.Bottom };
        Grid.SetColumnSpan(hl, headers.Length);
        g.Children.Add(hl);
        int r = 1;
        foreach (var row in rows.Take(maxRows))
        {
            g.RowDefinitions.Add(new RowDefinition(new GridLength(26)));
            for (int c = 0; c < row.Length; c++)
            {
                Control cell;
                if (row[c].Bar)
                {
                    var stack = new Grid { ColumnDefinitions = new ColumnDefinitions("38,*"), VerticalAlignment = VerticalAlignment.Center };
                    var label = Text(row[c].Text, "td");
                    var track = new Border { Height = 4, Background = new SolidColorBrush(Palette.Grid), VerticalAlignment = VerticalAlignment.Center, Margin = new Thickness(4, 0, 0, 0) };
                    var fill = new Border
                    {
                        Height = 4, HorizontalAlignment = HorizontalAlignment.Left, VerticalAlignment = VerticalAlignment.Center,
                        Background = new SolidColorBrush(row[c].Color ?? Palette.DyeC), Margin = new Thickness(4, 0, 0, 0),
                    };
                    fill.Width = 0;
                    double frac = Math.Clamp(row[c].BarValue, 0, 1);
                    track.SizeChanged += (_, e) => fill.Width = Math.Max(1, e.NewSize.Width * frac);
                    Grid.SetColumn(track, 1); Grid.SetColumn(fill, 1);
                    stack.Children.Add(label); stack.Children.Add(track); stack.Children.Add(fill);
                    cell = stack;
                }
                else
                {
                    var t = Text(row[c].Text, "td");
                    if (row[c].Color is { } col) t.Foreground = new SolidColorBrush(col);
                    t.TextTrimming = TextTrimming.CharacterEllipsis;
                    cell = t;
                }
                cell.Margin = new Thickness(c == 0 ? 0 : 10, 0, 6, 0);
                Grid.SetRow(cell, r); Grid.SetColumn(cell, c);
                g.Children.Add(cell);
            }
            var line = new Border { Height = 1, Background = new SolidColorBrush(Palette.Grid), VerticalAlignment = VerticalAlignment.Bottom };
            Grid.SetRow(line, r); Grid.SetColumnSpan(line, headers.Length);
            g.Children.Add(line);
            r++;
        }
        return g;
    }

    public static Button Button(string text, Action onClick, string cls = "quiet")
    {
        var b = new Button { Content = text };
        foreach (var c in cls.Split(' ', StringSplitOptions.RemoveEmptyEntries)) b.Classes.Add(c);
        b.Click += (_, _) => onClick();
        return b;
    }

    public static Control Labelled(string label, Control control)
    {
        var sp = new StackPanel { Spacing = 5 };
        sp.Children.Add(Text(label.ToUpperInvariant(), "eyebrow"));
        sp.Children.Add(control);
        return sp;
    }

    public static StackPanel Row(double spacing, params Control[] children)
    {
        var sp = new StackPanel { Orientation = Orientation.Horizontal, Spacing = spacing };
        foreach (var c in children) sp.Children.Add(c);
        return sp;
    }

    public static StackPanel Column(double spacing, params Control[] children)
    {
        var sp = new StackPanel { Spacing = spacing };
        foreach (var c in children) sp.Children.Add(c);
        return sp;
    }

    /// <summary>Grid with star/auto columns from a spec, children placed left to right.</summary>
    public static Grid Columns(string spec, double gap, params Control[] children)
    {
        var g = new Grid { ColumnDefinitions = new ColumnDefinitions(spec) };
        for (int i = 0; i < children.Length; i++)
        {
            if (i > 0) children[i].Margin = new Thickness(gap, children[i].Margin.Top, children[i].Margin.Right, children[i].Margin.Bottom);
            Grid.SetColumn(children[i], i);
            g.Children.Add(children[i]);
        }
        return g;
    }

    /// <summary>Renders a Biotech.Net.Visualization plot to a crisp bitmap (2× for high DPI).</summary>
    public static Bitmap PlotBitmap(BioPlot plot, double scale = 2)
    {
        using var ms = new MemoryStream(plot.ToPng(scale));
        return new Bitmap(ms);
    }

    public static Image PlotImage(BioPlot plot, double scale = 2) => new()
    {
        Source = PlotBitmap(plot, scale),
        Stretch = Stretch.Uniform,
        HorizontalAlignment = HorizontalAlignment.Left,
        VerticalAlignment = VerticalAlignment.Top,
    };

    /// <summary>A small coloured legend item.</summary>
    public static Control Legend(params (string Label, Color Color)[] items)
    {
        var sp = new StackPanel { Orientation = Orientation.Horizontal, Spacing = 12 };
        foreach (var (label, color) in items)
        {
            var row = new StackPanel { Orientation = Orientation.Horizontal, Spacing = 5 };
            row.Children.Add(new Border { Width = 8, Height = 8, Background = new SolidColorBrush(color), VerticalAlignment = VerticalAlignment.Center });
            row.Children.Add(Text(label, "unit"));
            sp.Children.Add(row);
        }
        return sp;
    }

    public static ScrollViewer Scroll(Control content, bool horizontal = false) => new()
    {
        Content = content,
        HorizontalScrollBarVisibility = horizontal ? ScrollBarVisibility.Auto : ScrollBarVisibility.Disabled,
        VerticalScrollBarVisibility = ScrollBarVisibility.Auto,
    };

    public static TextBox Input(string text, string cls = "", double height = double.NaN)
    {
        var t = new TextBox { Text = text, Height = height };
        foreach (var c in cls.Split(' ', StringSplitOptions.RemoveEmptyEntries)) t.Classes.Add(c);
        return t;
    }

    /// <summary>Notice block for empty or offline states: what is missing and how to enable it.</summary>
    public static Control Notice(string title, string body, Color accent)
    {
        var sp = new StackPanel { Spacing = 4 };
        var t = Text(title, "body");
        t.FontWeight = FontWeight.SemiBold;
        sp.Children.Add(t);
        sp.Children.Add(Text(body, "lede"));
        return new Border
        {
            Child = sp,
            BorderBrush = new SolidColorBrush(accent),
            BorderThickness = new Thickness(2, 0, 0, 0),
            Background = new SolidColorBrush(accent, 0.06),
            Padding = new Thickness(12, 8),
        };
    }
}
