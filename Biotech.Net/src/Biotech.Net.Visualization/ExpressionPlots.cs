using Biotech.Net.Transcriptomics;
using Biotech.Net.Transcriptomics.SingleCell;
using Matplotlib.Net;
using Matplotlib.Net.Artists;
using Matplotlib.Net.Rendering;
using Plotly.Net;
using Plotly.Net.Traces;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>Regulation call of one gene.</summary>
public enum Regulation
{
    NotSignificant,
    Up,
    Down,
}

/// <summary>
/// Shared logic of volcano and MA plots: significance calls from adjusted p-values and a fold-change threshold,
/// colouring, and labelling of the top genes.
/// </summary>
public abstract class DifferentialPlot : FigurePlot
{
    protected DifferentialPlot(IReadOnlyList<string> genes, double[] log2FoldChange, double[] pValue, double[]? pAdj, int width, int height)
        : base(width, height)
    {
        int n = genes.Count;
        if (log2FoldChange.Length != n || pValue.Length != n || (pAdj is not null && pAdj.Length != n))
            throw new ArgumentException("All columns need one value per gene.");
        Genes = genes;
        Log2FoldChange = log2FoldChange;
        PValue = pValue;
        PAdj = pAdj ?? BenjaminiHochberg(pValue);
    }

    public IReadOnlyList<string> Genes { get; }
    public double[] Log2FoldChange { get; }
    public double[] PValue { get; }

    /// <summary>Adjusted p-values (computed with Benjamini–Hochberg when not supplied).</summary>
    public double[] PAdj { get; }

    /// <summary>Significance level on the adjusted p-value.</summary>
    public double Alpha { get; set; } = 0.05;

    /// <summary>Minimum |log2 fold change| for a gene to be called up or down.</summary>
    public double FoldChangeThreshold { get; set; } = 1;

    /// <summary>Number of most significant genes to label.</summary>
    public int LabelTop { get; set; } = 10;

    /// <summary>Genes to always label and outline.</summary>
    public IReadOnlyCollection<string>? Highlight { get; set; }

    /// <summary>Point size (pt²-like Matplotlib.Net scatter size).</summary>
    public double PointSize { get; set; } = 14;

    /// <summary>Regulation call of each gene.</summary>
    public Regulation[] Classify()
    {
        var c = new Regulation[Genes.Count];
        for (int i = 0; i < c.Length; i++)
        {
            double q = PAdj[i], lfc = Log2FoldChange[i];
            if (double.IsNaN(q) || double.IsNaN(lfc) || q >= Alpha || Math.Abs(lfc) < FoldChangeThreshold) continue;
            c[i] = lfc > 0 ? Regulation.Up : Regulation.Down;
        }
        return c;
    }

    /// <summary>Counts per call.</summary>
    public (int Up, int Down, int NotSignificant) Counts()
    {
        var c = Classify();
        return (c.Count(x => x == Regulation.Up), c.Count(x => x == Regulation.Down), c.Count(x => x == Regulation.NotSignificant));
    }

    /// <summary>Indices of genes to label: highlighted ones plus the <see cref="LabelTop"/> most significant called genes.</summary>
    public int[] LabelledGenes()
    {
        var calls = Classify();
        var set = new HashSet<string>(Highlight ?? [], StringComparer.OrdinalIgnoreCase);
        var top = Enumerable.Range(0, Genes.Count)
            .Where(i => calls[i] != Regulation.NotSignificant)
            .OrderBy(i => double.IsNaN(PValue[i]) ? double.MaxValue : PValue[i])
            .Take(Math.Max(0, LabelTop));
        return top.Concat(Enumerable.Range(0, Genes.Count).Where(i => set.Contains(Genes[i]))).Distinct().ToArray();
    }

    /// <summary>Plot coordinates of gene <paramref name="i"/> (NaN when not drawable).</summary>
    protected abstract (double X, double Y) Coordinates(int i);

    protected abstract string DefaultXLabel { get; }
    protected abstract string DefaultYLabel { get; }

    protected Color ColorOf(Regulation r) => r switch { Regulation.Up => Theme.Up, Regulation.Down => Theme.Down, _ => Theme.Neutral };

    protected override void Build(Figure figure)
    {
        var ax = figure.AddSubplot();
        var calls = Classify();
        var (up, down, ns) = (calls.Count(c => c == Regulation.Up), calls.Count(c => c == Regulation.Down), calls.Count(c => c == Regulation.NotSignificant));
        foreach (var (call, label) in new[] { (Regulation.NotSignificant, $"NS ({ns:N0})"), (Regulation.Down, $"Down ({down:N0})"), (Regulation.Up, $"Up ({up:N0})") })
        {
            var idx = Enumerable.Range(0, Genes.Count).Where(i => calls[i] == call).Select(i => (i, p: Coordinates(i)))
                .Where(t => double.IsFinite(t.p.X) && double.IsFinite(t.p.Y)).ToList();
            if (idx.Count == 0) continue;
            ax.Scatter(idx.Select(t => t.p.X).ToArray(), idx.Select(t => t.p.Y).ToArray(), label, ColorOf(call),
                sizes: Sizes(idx.Count, PointSize), alpha: call == Regulation.NotSignificant ? 0.55 : 0.85, edgeWidth: 0);
        }
        AddGuides(ax);
        var highlight = new HashSet<string>(Highlight ?? [], StringComparer.OrdinalIgnoreCase);
        foreach (int i in LabelledGenes())
        {
            var (x, y) = Coordinates(i);
            if (!double.IsFinite(x) || !double.IsFinite(y)) continue;
            if (highlight.Contains(Genes[i]))
                ax.Scatter([x], [y], color: Color.Transparent, sizes: [PointSize * 3], alpha: 1, edgeColor: Theme.Foreground, edgeWidth: 1.2,
                    marker: MarkerStyle.Circle);
            ax.Text(x, y, " " + Genes[i], Theme.Foreground, Theme.FontSize * 0.8, x >= 0 ? TextAlign.Left : TextAlign.Right, TextBaseline.Bottom);
        }
        ax.SetXLabel(XLabel ?? DefaultXLabel);
        ax.SetYLabel(YLabel ?? DefaultYLabel);
        if (Title is not null) ax.SetTitle(Title);
        ax.Legend(LegendLocation.Best);
    }

    protected abstract void AddGuides(Axes ax);

    /// <inheritdoc />
    public override bool SupportsInteractive => true;

    /// <inheritdoc />
    public override PlotlyFigure ToPlotly()
    {
        var f = NewPlotly();
        var calls = Classify();
        foreach (var call in new[] { Regulation.NotSignificant, Regulation.Down, Regulation.Up })
        {
            var idx = Enumerable.Range(0, Genes.Count).Where(i => calls[i] == call && double.IsFinite(Coordinates(i).X) && double.IsFinite(Coordinates(i).Y)).ToArray();
            if (idx.Length == 0) continue;
            f.Add(new Scatter(webgl: idx.Length > 5000)
            {
                Name = call switch { Regulation.Up => "Up", Regulation.Down => "Down", _ => "NS" },
                X = idx.Select(i => Coordinates(i).X).ToArray(),
                Y = idx.Select(i => Coordinates(i).Y).ToArray(),
                Mode = ScatterMode.Markers,
                Text = idx.Select(i => Genes[i]).ToArray(),
                HoverTemplate = "<b>%{text}</b><br>x=%{x:.3f}<br>y=%{y:.3f}<extra></extra>",
                Marker = new Marker { Color = ColorOf(call), Size = 6, Opacity = call == Regulation.NotSignificant ? 0.5 : 0.85 },
            });
        }
        f.Layout.XAxis.Title = XLabel ?? DefaultXLabel;
        f.Layout.YAxis.Title = YLabel ?? DefaultYLabel;
        return f;
    }

    /// <summary>Benjamini–Hochberg adjusted p-values (NaN stays NaN).</summary>
    public static double[] BenjaminiHochberg(double[] p)
    {
        var idx = Enumerable.Range(0, p.Length).Where(i => !double.IsNaN(p[i])).OrderBy(i => p[i]).ToArray();
        var q = new double[p.Length];
        Array.Fill(q, double.NaN);
        double min = 1;
        for (int k = idx.Length - 1; k >= 0; k--)
        {
            min = Math.Min(min, p[idx[k]] * idx.Length / (k + 1));
            q[idx[k]] = Math.Min(1, min);
        }
        return q;
    }
}

/// <summary>Volcano plot: log2 fold change against −log10 p, with threshold guides and labelled top genes.</summary>
public sealed class VolcanoPlot : DifferentialPlot
{
    public VolcanoPlot(IReadOnlyList<string> genes, double[] log2FoldChange, double[] pValue, double[]? pAdj = null)
        : base(genes, log2FoldChange, pValue, pAdj, 760, 560) { }

    /// <summary>From a Transcriptomics result (DESeq2, limma, t-test).</summary>
    public VolcanoPlot(DifferentialExpressionResult result)
        : this(result.Genes, result.Log2FoldChange, result.PValue, result.PAdj) => Title = $"{result.Method}: {result.Contrast}";

    /// <summary>Plot −log10 of the adjusted p-value instead of the raw p-value.</summary>
    public bool UseAdjustedPValue { get; set; }

    protected override (double X, double Y) Coordinates(int i)
    {
        double p = UseAdjustedPValue ? PAdj[i] : PValue[i];
        if (double.IsNaN(p) || double.IsNaN(Log2FoldChange[i])) return (double.NaN, double.NaN);
        return (Log2FoldChange[i], -Math.Log10(Math.Max(p, 1e-300)));
    }

    protected override string DefaultXLabel => "log2 fold change";
    protected override string DefaultYLabel => UseAdjustedPValue ? "−log10 adjusted p" : "−log10 p";

    protected override void AddGuides(Axes ax)
    {
        var g = Theme.Muted;
        if (FoldChangeThreshold > 0)
        {
            ax.AxVLine(-FoldChangeThreshold, g, LineStyle.Dashed, 1);
            ax.AxVLine(FoldChangeThreshold, g, LineStyle.Dashed, 1);
        }
        // The p-value at which the adjusted p crosses alpha (largest raw p among significant genes).
        double cut = double.NaN;
        for (int i = 0; i < PValue.Length; i++)
            if (PAdj[i] < Alpha && !double.IsNaN(PValue[i])) cut = double.IsNaN(cut) ? PValue[i] : Math.Max(cut, PValue[i]);
        double y = UseAdjustedPValue ? -Math.Log10(Alpha) : double.IsNaN(cut) ? double.NaN : -Math.Log10(Math.Max(cut, 1e-300));
        if (double.IsFinite(y)) ax.AxHLine(y, g, LineStyle.Dashed, 1);
    }
}

/// <summary>MA plot: mean expression (A) against log2 fold change (M).</summary>
public sealed class MaPlot : DifferentialPlot
{
    private readonly double[] _a;

    public MaPlot(IReadOnlyList<string> genes, double[] meanExpression, double[] log2FoldChange, double[] pValue, double[]? pAdj = null, bool meanIsLog = false)
        : base(genes, log2FoldChange, pValue, pAdj, 760, 520)
    {
        if (meanExpression.Length != genes.Count) throw new ArgumentException("meanExpression needs one value per gene.");
        _a = meanIsLog ? meanExpression : meanExpression.Select(v => v > 0 ? Math.Log2(v) : double.NaN).ToArray();
        MeanIsLog = meanIsLog;
    }

    /// <summary>From a Transcriptomics result; <paramref name="baseMeanIsLog"/> for limma / t-test tables whose base mean is already log-scale.</summary>
    public MaPlot(DifferentialExpressionResult result, bool? baseMeanIsLog = null)
        : this(result.Genes, result.BaseMean, result.Log2FoldChange, result.PValue, result.PAdj,
            baseMeanIsLog ?? !result.Method.StartsWith("DESeq", StringComparison.OrdinalIgnoreCase))
        => Title = $"{result.Method}: {result.Contrast}";

    /// <summary>Whether the mean was supplied on the log scale.</summary>
    public bool MeanIsLog { get; }

    protected override (double X, double Y) Coordinates(int i) => (_a[i], Log2FoldChange[i]);

    protected override string DefaultXLabel => MeanIsLog ? "average log2 expression (A)" : "log2 mean of normalised counts (A)";
    protected override string DefaultYLabel => "log2 fold change (M)";

    protected override void AddGuides(Axes ax)
    {
        ax.AxHLine(0, Theme.Muted, LineStyle.Solid, 1);
        if (FoldChangeThreshold > 0)
        {
            ax.AxHLine(FoldChangeThreshold, Theme.Muted, LineStyle.Dashed, 1);
            ax.AxHLine(-FoldChangeThreshold, Theme.Muted, LineStyle.Dashed, 1);
        }
    }
}

/// <summary>
/// A 2-D embedding scatter (PCA, UMAP, t-SNE) coloured by group or by a continuous value, with axis labels that
/// carry the explained variance for PCA.
/// </summary>
public sealed class EmbeddingPlot : FigurePlot
{
    public EmbeddingPlot(double[] x, double[] y, IReadOnlyList<string>? groups = null, IReadOnlyList<string>? labels = null)
        : base(720, 560)
    {
        if (x.Length != y.Length) throw new ArgumentException("x and y need the same length.");
        if (groups is not null && groups.Count != x.Length) throw new ArgumentException("groups needs one entry per point.");
        if (labels is not null && labels.Count != x.Length) throw new ArgumentException("labels needs one entry per point.");
        X = x; Y = y; Groups = groups; Labels = labels;
    }

    /// <summary>PCA scores of two components; axis labels include the variance explained.</summary>
    public static EmbeddingPlot FromPca(PcaResult pca, IReadOnlyList<string>? groups = null, IReadOnlyList<string>? sampleNames = null,
        int xComponent = 0, int yComponent = 1)
    {
        if (xComponent >= pca.Components || yComponent >= pca.Components) throw new ArgumentOutOfRangeException(nameof(yComponent));
        return new EmbeddingPlot(pca.Scores.Column(xComponent), pca.Scores.Column(yComponent), groups, sampleNames)
        {
            Title = "PCA",
            XLabel = $"PC{xComponent + 1} ({pca.ExplainedVarianceRatio[xComponent] * 100:F1}%)",
            YLabel = $"PC{yComponent + 1} ({pca.ExplainedVarianceRatio[yComponent] * 100:F1}%)",
            ShowLabels = sampleNames is not null && pca.Scores.Rows <= 40,
        };
    }

    /// <summary>A UMAP / t-SNE embedding (first two columns).</summary>
    public static EmbeddingPlot FromEmbedding(DenseMatrix embedding, IReadOnlyList<string>? groups = null, string method = "UMAP")
    {
        if (embedding.Columns < 2) throw new ArgumentException("The embedding needs at least two columns.");
        return new EmbeddingPlot(embedding.Column(0), embedding.Column(1), groups)
        {
            Title = method, XLabel = $"{method} 1", YLabel = $"{method} 2", PointSize = embedding.Rows > 2000 ? 4 : 10,
        };
    }

    /// <summary>
    /// An embedding stored on a single-cell experiment (<c>X_umap</c>, <c>X_tsne</c>, <c>X_pca</c>), coloured by a cell
    /// label (e.g. <c>leiden</c>) or a cell metric.
    /// </summary>
    public static EmbeddingPlot FromSingleCell(SingleCellExperiment sc, string embedding = "X_umap", string? colorBy = null)
    {
        if (!sc.Embeddings.TryGetValue(embedding, out var e)) throw new KeyNotFoundException($"No embedding '{embedding}'. Available: {string.Join(", ", sc.Embeddings.Keys)}.");
        string method = embedding.Replace("X_", "").ToUpperInvariant();
        IReadOnlyList<string>? groups = null;
        double[]? values = null;
        if (colorBy is not null)
        {
            if (sc.CellLabels.TryGetValue(colorBy, out var labels)) groups = labels;
            else if (sc.CellMetrics.TryGetValue(colorBy, out var metric)) values = metric;
            else throw new KeyNotFoundException($"No cell label or metric '{colorBy}'.");
        }
        var plot = FromEmbedding(e, groups, method);
        plot.Values = values;
        plot.ValueLabel = colorBy;
        return plot;
    }

    public double[] X { get; }
    public double[] Y { get; }
    public IReadOnlyList<string>? Groups { get; }
    public IReadOnlyList<string>? Labels { get; }

    /// <summary>Continuous colouring (overrides <see cref="Groups"/>).</summary>
    public double[]? Values { get; set; }

    /// <summary>Caption of the colour bar for <see cref="Values"/>.</summary>
    public string? ValueLabel { get; set; }

    /// <summary>Draw point labels.</summary>
    public bool ShowLabels { get; set; }

    /// <summary>Marker size.</summary>
    public double PointSize { get; set; } = 36;

    /// <summary>Distinct groups in first-seen order.</summary>
    public IReadOnlyList<string> GroupLevels() => Groups?.Distinct().ToList() ?? [];

    protected override void Build(Figure figure)
    {
        var ax = figure.AddSubplot();
        if (Values is { } values)
        {
            ax.Scatter(X, Y, sizes: Sizes(X.Length, PointSize), colorValues: values, colormap: Theme.SequentialColormap, alpha: 0.9, edgeWidth: 0);
            ax.AddColorbar(ValueLabel);
        }
        else if (Groups is not null)
        {
            var levels = GroupLevels();
            for (int g = 0; g < levels.Count; g++)
            {
                var idx = Enumerable.Range(0, X.Length).Where(i => Groups[i] == levels[g]).ToArray();
                ax.Scatter(idx.Select(i => X[i]).ToArray(), idx.Select(i => Y[i]).ToArray(), levels[g], Theme.CategoryColor(g),
                    sizes: Sizes(idx.Length, PointSize), alpha: 0.88, edgeColor: Theme.Panel, edgeWidth: PointSize >= 20 ? 0.8 : 0);
            }
            ax.Legend(levels.Count > 8 ? LegendLocation.OutsideRight : LegendLocation.Best, columns: levels.Count > 16 ? 2 : 1);
        }
        else
        {
            ax.Scatter(X, Y, color: Theme.Accent, sizes: Sizes(X.Length, PointSize), alpha: 0.85, edgeWidth: 0);
        }
        if (ShowLabels && Labels is not null)
            for (int i = 0; i < X.Length; i++)
                ax.Text(X[i], Y[i], " " + Labels[i], Theme.Muted, Theme.FontSize * 0.75, TextAlign.Left, TextBaseline.Bottom);
        ax.SetXLabel(XLabel ?? "dimension 1");
        ax.SetYLabel(YLabel ?? "dimension 2");
        if (Title is not null) ax.SetTitle(Title);
    }

    /// <inheritdoc />
    public override bool SupportsInteractive => true;

    /// <inheritdoc />
    public override PlotlyFigure ToPlotly()
    {
        var f = NewPlotly();
        bool webgl = X.Length > 5000;
        if (Values is { } values)
        {
            f.Add(new Scatter(webgl)
            {
                Name = ValueLabel ?? "value", X = X, Y = Y, Mode = ScatterMode.Markers, Text = Labels?.ToArray(),
                Marker = new Marker { Colors = values, ColorScale = "Viridis", ShowScale = true, Size = 6 },
            });
        }
        else if (Groups is not null)
        {
            var levels = GroupLevels();
            for (int g = 0; g < levels.Count; g++)
            {
                var idx = Enumerable.Range(0, X.Length).Where(i => Groups[i] == levels[g]).ToArray();
                f.Add(new Scatter(webgl)
                {
                    Name = levels[g], X = idx.Select(i => X[i]).ToArray(), Y = idx.Select(i => Y[i]).ToArray(), Mode = ScatterMode.Markers,
                    Text = Labels is null ? null : idx.Select(i => Labels[i]).ToArray(),
                    Marker = new Marker { Color = Theme.CategoryColor(g), Size = 7, Opacity = 0.88 },
                });
            }
        }
        else
        {
            f.Add(new Scatter(webgl) { X = X, Y = Y, Mode = ScatterMode.Markers, Text = Labels?.ToArray(), Marker = new Marker { Color = Theme.Accent, Size = 7 } });
        }
        f.Layout.XAxis.Title = XLabel ?? "dimension 1";
        f.Layout.YAxis.Title = YLabel ?? "dimension 2";
        return f;
    }
}
