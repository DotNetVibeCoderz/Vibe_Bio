using Biotech.Net.Core;
using Matplotlib.Net.Rendering;
using Plotly.Net;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>Standardisation applied before colouring a heatmap.</summary>
public enum HeatmapScale
{
    /// <summary>Raw values.</summary>
    None,

    /// <summary>Z-score each row (gene-wise, the usual expression heatmap).</summary>
    Rows,

    /// <summary>Z-score each column.</summary>
    Columns,
}

/// <summary>Computed geometry of a <see cref="Heatmap"/>.</summary>
public sealed record HeatmapLayout(
    int[] RowOrder, int[] ColumnOrder, Dendrogram? RowTree, Dendrogram? ColumnTree,
    Rect Cells, Rect? RowDendrogram, Rect? ColumnDendrogram, Rect? GroupStrip, Rect Colorbar, double VMin, double VMax);

/// <summary>
/// A clustered heatmap (pheatmap / seaborn clustermap style): optional row or column z-scoring, hierarchical
/// clustering with dendrograms, a column-group annotation strip and a colour bar. Rendered directly on the
/// Matplotlib.Net renderer with Visual.Net colormaps; <see cref="ToPlotly"/> gives an interactive version in the
/// same clustered order.
/// </summary>
public sealed class Heatmap : BioPlot
{
    private readonly double[,] _values;

    public Heatmap(double[,] values, IReadOnlyList<string>? rowLabels = null, IReadOnlyList<string>? columnLabels = null)
        : base(820, 640)
    {
        _values = values ?? throw new ArgumentNullException(nameof(values));
        int r = values.GetLength(0), c = values.GetLength(1);
        if (r == 0 || c == 0) throw new ArgumentException("The heatmap needs at least one row and one column.");
        RowLabels = rowLabels ?? Enumerable.Range(1, r).Select(i => $"row {i}").ToArray();
        ColumnLabels = columnLabels ?? Enumerable.Range(1, c).Select(i => $"col {i}").ToArray();
        if (RowLabels.Count != r || ColumnLabels.Count != c) throw new ArgumentException("Label counts must match the matrix shape.");
    }

    /// <summary>A heatmap of a Core <see cref="Dataset"/> (rows = features, columns = samples). A "group" or "condition" column annotation becomes the group strip.</summary>
    public static Heatmap FromDataset(Dataset dataset, string? groupAnnotation = null)
    {
        var v = new double[dataset.Rows, dataset.Columns];
        for (int i = 0; i < dataset.Rows; i++)
            for (int j = 0; j < dataset.Columns; j++) v[i, j] = dataset[i, j];
        var h = new Heatmap(v, dataset.RowNames, dataset.ColumnNames) { Title = dataset.Name };
        string? key = groupAnnotation ?? new[] { "group", "condition", "treatment" }.FirstOrDefault(k => dataset.ColumnAnnotations.ContainsKey(k));
        if (key is not null && dataset.ColumnAnnotations.TryGetValue(key, out var groups) && groups.Length == dataset.Columns)
            h.ColumnGroups = groups;
        return h;
    }

    public IReadOnlyList<string> RowLabels { get; }
    public IReadOnlyList<string> ColumnLabels { get; }

    /// <summary>Number of rows / columns.</summary>
    public int Rows => _values.GetLength(0);
    public int Columns => _values.GetLength(1);

    /// <summary>Standardisation (default: z-score rows).</summary>
    public HeatmapScale Scale { get; set; } = HeatmapScale.Rows;

    /// <summary>Cluster and reorder rows / columns.</summary>
    public bool ClusterRows { get; set; } = true;
    public bool ClusterColumns { get; set; } = true;

    /// <summary>Show dendrograms for clustered dimensions.</summary>
    public bool ShowDendrograms { get; set; } = true;

    public DistanceMetric Metric { get; set; } = DistanceMetric.Euclidean;
    public LinkageMethod Linkage { get; set; } = LinkageMethod.Average;

    /// <summary>Colormap name (null: diverging when scaled, sequential otherwise).</summary>
    public string? Colormap { get; set; }

    /// <summary>Colour limits (null: symmetric ±<see cref="ZClip"/> when scaled, data range otherwise).</summary>
    public double? Min { get; set; }
    public double? Max { get; set; }

    /// <summary>Clip for z-scores.</summary>
    public double ZClip { get; set; } = 3;

    /// <summary>Optional group of each column (e.g. condition), drawn as a coloured strip.</summary>
    public IReadOnlyList<string>? ColumnGroups { get; set; }

    /// <summary>Row labels are shown when there are at most this many rows.</summary>
    public int MaxRowLabels { get; set; } = 60;

    /// <summary>Column labels are shown when there are at most this many columns.</summary>
    public int MaxColumnLabels { get; set; } = 80;

    /// <summary>Colour-bar caption (default "z-score" or "value").</summary>
    public string? ColorbarLabel { get; set; }

    /// <summary>The matrix after <see cref="Scale"/> (input order).</summary>
    public double[,] ScaledValues() => Standardise(_values, Scale);

    /// <summary>Z-scores rows or columns (population SD; constant vectors become 0).</summary>
    public static double[,] Standardise(double[,] values, HeatmapScale scale)
    {
        int r = values.GetLength(0), c = values.GetLength(1);
        var z = (double[,])values.Clone();
        if (scale == HeatmapScale.None) return z;
        bool rows = scale == HeatmapScale.Rows;
        int outer = rows ? r : c, inner = rows ? c : r;
        for (int o = 0; o < outer; o++)
        {
            double sum = 0; int n = 0;
            for (int k = 0; k < inner; k++) { double v = rows ? values[o, k] : values[k, o]; if (!double.IsNaN(v)) { sum += v; n++; } }
            double mean = n > 0 ? sum / n : 0, ss = 0;
            for (int k = 0; k < inner; k++) { double v = rows ? values[o, k] : values[k, o]; if (!double.IsNaN(v)) ss += (v - mean) * (v - mean); }
            double sd = n > 0 ? Math.Sqrt(ss / n) : 0;
            for (int k = 0; k < inner; k++)
            {
                double v = rows ? values[o, k] : values[k, o];
                double s = double.IsNaN(v) ? double.NaN : sd > 0 ? (v - mean) / sd : 0;
                if (rows) z[o, k] = s; else z[k, o] = s;
            }
        }
        return z;
    }

    private string EffectiveColormap => Colormap ?? (Scale != HeatmapScale.None ? Theme.DivergingColormap : Theme.SequentialColormap);

    private (double Min, double Max) Limits(double[,] z)
    {
        if (Min is { } a && Max is { } b) return (a, b);
        if (Scale != HeatmapScale.None)
        {
            double m = 0;
            foreach (double v in z) if (!double.IsNaN(v)) m = Math.Max(m, Math.Abs(v));
            m = Math.Min(Math.Max(m, 1e-9), ZClip);
            return (Min ?? -m, Max ?? m);
        }
        double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
        foreach (double v in z) if (!double.IsNaN(v)) { lo = Math.Min(lo, v); hi = Math.Max(hi, v); }
        if (double.IsInfinity(lo)) (lo, hi) = (0, 1);
        if (hi <= lo) hi = lo + 1;
        return (Min ?? lo, Max ?? hi);
    }

    /// <summary>Clusters (when enabled) and lays out the plot deterministically.</summary>
    public HeatmapLayout ComputeLayout()
    {
        var z = ScaledValues();
        Dendrogram? rowTree = null, colTree = null;
        if (ClusterRows && Rows > 1)
            rowTree = HierarchicalClustering.Cluster(Enumerable.Range(0, Rows).Select(i => RowVector(z, i)).ToList(), Metric, Linkage);
        if (ClusterColumns && Columns > 1)
            colTree = HierarchicalClustering.Cluster(Enumerable.Range(0, Columns).Select(j => ColumnVector(z, j)).ToList(), Metric, Linkage);
        int[] rowOrder = rowTree?.Order ?? Enumerable.Range(0, Rows).ToArray();
        int[] colOrder = colTree?.Order ?? Enumerable.Range(0, Columns).ToArray();

        double pad = 14, top = pad + (string.IsNullOrEmpty(Title) ? 0 : 30);
        bool rowLabels = Rows <= MaxRowLabels, colLabels = Columns <= MaxColumnLabels;
        double rowLabelW = rowLabels ? Math.Min(170, 12 + 6.6 * RowLabels.Max(l => l.Length)) : 0;
        double colLabelH = colLabels ? Math.Min(130, 14 + 6.0 * ColumnLabels.Max(l => l.Length)) : 0;
        double dendroLeft = rowTree is not null && ShowDendrograms ? 90 : 0;
        double dendroTop = colTree is not null && ShowDendrograms ? 70 : 0;
        double stripH = ColumnGroups is not null ? 16 : 0;
        double legendH = ColumnGroups is not null ? 22 : 0;
        double colorbarW = 70;

        double x0 = pad + dendroLeft;
        double y0 = top + dendroTop + stripH + (stripH > 0 ? 4 : 0);
        double w = Math.Max(20, Width - x0 - rowLabelW - colorbarW - pad);
        double h = Math.Max(20, Height - y0 - colLabelH - legendH - pad);
        var cells = new Rect(x0, y0, w, h);
        Rect? rd = dendroLeft > 0 ? new Rect(pad, y0, dendroLeft - 6, h) : null;
        Rect? cd = dendroTop > 0 ? new Rect(x0, top, w, dendroTop - 6) : null;
        Rect? strip = stripH > 0 ? new Rect(x0, y0 - stripH - 4, w, stripH) : null;
        var cb = new Rect(cells.Right + rowLabelW + 18, y0, 14, Math.Min(h, 220));
        var (vmin, vmax) = Limits(z);
        return new HeatmapLayout(rowOrder, colOrder, rowTree, colTree, cells, rd, cd, strip, cb, vmin, vmax);
    }

    /// <inheritdoc />
    public override void Render(IRenderer r)
    {
        var th = Theme;
        r.Clear(th.Background);
        DrawTitle(r, 12);
        var L = ComputeLayout();
        var z = ScaledValues();
        var cmap = Visual.Net.Core.Colors.Colormaps.Get(EffectiveColormap);
        var norm = new Normalize(L.VMin, L.VMax, clip: true);
        double cw = L.Cells.Width / Columns, ch = L.Cells.Height / Rows;

        Color CellColor(double v) => double.IsNaN(v) ? th.Grid : cmap.At(norm.Map(v));

        if ((long)Rows * Columns <= 40_000)
        {
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                {
                    var rect = new Rect(L.Cells.X + j * cw, L.Cells.Y + i * ch, cw + 0.35, ch + 0.35);
                    r.DrawRect(rect, new Brush(CellColor(z[L.RowOrder[i], L.ColumnOrder[j]])), Pen.None);
                }
        }
        else
        {
            var pixels = new uint[Rows * Columns];
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    pixels[i * Columns + j] = CellColor(z[L.RowOrder[i], L.ColumnOrder[j]]).ToArgb();
            r.DrawImage(pixels, Columns, Rows, L.Cells, smooth: false);
        }
        r.DrawRect(L.Cells, Brush.None, Paint.Pen(th.Border, 1));

        // Dendrograms.
        var dpen = Paint.Pen(th.Muted, 1);
        if (L.RowDendrogram is { } rd && L.RowTree is { } rt && rt.RootHeight > 0)
            foreach (var (x1, y1, x2, y2) in rt.Segments())
                r.DrawLine(new Pt(rd.Right - y1 / rt.RootHeight * rd.Width, rd.Y + x1 * ch), new Pt(rd.Right - y2 / rt.RootHeight * rd.Width, rd.Y + x2 * ch), dpen);
        if (L.ColumnDendrogram is { } cd && L.ColumnTree is { } ct && ct.RootHeight > 0)
            foreach (var (x1, y1, x2, y2) in ct.Segments())
                r.DrawLine(new Pt(cd.X + x1 * cw, cd.Bottom - y1 / ct.RootHeight * cd.Height), new Pt(cd.X + x2 * cw, cd.Bottom - y2 / ct.RootHeight * cd.Height), dpen);

        // Group strip and legend.
        if (L.GroupStrip is { } gs && ColumnGroups is { } groups)
        {
            var distinct = groups.Distinct().ToList();
            for (int j = 0; j < Columns; j++)
                r.DrawRect(new Rect(gs.X + j * cw, gs.Y, cw + 0.35, gs.Height), new Brush(th.CategoryColor(distinct.IndexOf(groups[L.ColumnOrder[j]]))), Pen.None);
            double lx = L.Cells.X, ly = Height - 26;
            for (int g = 0; g < distinct.Count; g++)
            {
                r.DrawRect(new Rect(lx, ly, 12, 12), new Brush(th.CategoryColor(g)), Pen.None);
                r.DrawText(distinct[g], new Pt(lx + 16, ly + 6), th.SmallFont, th.Foreground, TextAlign.Left, TextBaseline.Middle);
                lx += 32 + r.MeasureText(distinct[g], th.SmallFont).Width;
            }
        }

        // Labels.
        var small = th.Font(Math.Clamp(Math.Min(ch, 12) * 0.9, 7, 11));
        if (Rows <= MaxRowLabels)
            for (int i = 0; i < Rows; i++)
                r.DrawText(RowLabels[L.RowOrder[i]], new Pt(L.Cells.Right + 6, L.Cells.Y + (i + 0.5) * ch), small, th.Foreground, TextAlign.Left, TextBaseline.Middle);
        var csmall = th.Font(Math.Clamp(Math.Min(cw, 12) * 0.9, 7, 11));
        if (Columns <= MaxColumnLabels)
            for (int j = 0; j < Columns; j++)
                r.DrawText(ColumnLabels[L.ColumnOrder[j]], new Pt(L.Cells.X + (j + 0.5) * cw, L.Cells.Bottom + 6), csmall, th.Foreground, TextAlign.Right, TextBaseline.Middle, -90);

        // Colour bar.
        var cb = L.Colorbar;
        int steps = 64;
        for (int s = 0; s < steps; s++)
        {
            double t = 1 - (s + 0.5) / steps;
            r.DrawRect(new Rect(cb.X, cb.Y + s * cb.Height / steps, cb.Width, cb.Height / steps + 0.5), new Brush(cmap.At(t)), Pen.None);
        }
        r.DrawRect(cb, Brush.None, Paint.Pen(th.Border, 1));
        foreach (double t in Paint.NiceTicks(L.VMin, L.VMax, 4))
        {
            double y = cb.Bottom - (t - L.VMin) / (L.VMax - L.VMin) * cb.Height;
            Paint.Line(r, cb.Right, y, cb.Right + 4, y, th.Border);
            r.DrawText(t.ToString("G3", System.Globalization.CultureInfo.InvariantCulture), new Pt(cb.Right + 6, y), th.SmallFont, th.Muted, TextAlign.Left, TextBaseline.Middle);
        }
        r.DrawText(ColorbarLabel ?? (Scale == HeatmapScale.None ? "value" : "z-score"), new Pt(cb.X, cb.Y - 6), th.SmallFont, th.Muted, TextAlign.Left, TextBaseline.Bottom);
    }

    /// <inheritdoc />
    public override bool SupportsInteractive => true;

    /// <inheritdoc />
    public override PlotlyFigure ToPlotly()
    {
        var L = ComputeLayout();
        var z = ScaledValues();
        var flat = new double[Rows, Columns];
        // Plotly draws row 0 at the bottom; reverse so the first clustered row is on top.
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Columns; j++)
                flat[Rows - 1 - i, j] = z[L.RowOrder[i], L.ColumnOrder[j]];
        var f = NewPlotly();
        var trace = f.Add(new Plotly.Net.Traces.Heatmap
        {
            Z = flat,
            X = L.ColumnOrder.Select(j => ColumnLabels[j]).ToArray(),
            Y = L.RowOrder.Reverse().Select(i => RowLabels[i]).ToArray(),
            ColorScale = Scale != HeatmapScale.None ? "RdBu" : "Viridis",
            Reversescale = Scale != HeatmapScale.None,
            ZMin = L.VMin,
            ZMax = L.VMax,
        });
        trace.SetColorbarTitle(ColorbarLabel ?? (Scale == HeatmapScale.None ? "value" : "z-score"));
        return f;
    }

    private static double[] RowVector(double[,] z, int i)
    {
        var v = new double[z.GetLength(1)];
        for (int j = 0; j < v.Length; j++) v[j] = double.IsNaN(z[i, j]) ? 0 : z[i, j];
        return v;
    }

    private static double[] ColumnVector(double[,] z, int j)
    {
        var v = new double[z.GetLength(0)];
        for (int i = 0; i < v.Length; i++) v[i] = double.IsNaN(z[i, j]) ? 0 : z[i, j];
        return v;
    }
}
