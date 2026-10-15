using Biotech.Net.Core;
using Matplotlib.Net;
using Matplotlib.Net.Artists;
using Matplotlib.Net.Rendering;
using Plotly.Net;
using Plotly.Net.Traces;

namespace Biotech.Net.Visualization;

/// <summary>Letter-height scale of a <see cref="SequenceLogo"/>.</summary>
public enum LogoMode
{
    /// <summary>Information content in bits (letter height = frequency × column information).</summary>
    Bits,

    /// <summary>Plain frequency (columns sum to 1).</summary>
    Probability,
}

/// <summary>One stacked letter of a logo column.</summary>
public readonly record struct LogoLetter(char Symbol, double Frequency, double Height);

/// <summary>One position of a logo: letters sorted from smallest (bottom) to largest (top).</summary>
public sealed record LogoColumn(int Position, double Information, IReadOnlyList<LogoLetter> Letters);

/// <summary>
/// A sequence logo (Schneider &amp; Stephens 1990) of aligned sequences or a position frequency matrix, in bits or
/// probability, with optional small-sample correction. Letter outlines are real font glyphs stretched to their
/// height, so the logo is crisp in both PNG and SVG.
/// </summary>
public sealed class SequenceLogo : BioPlot
{
    private readonly double[,] _counts;

    /// <summary>A logo from aligned sequences of equal length (gaps '-' and '.' are not counted).</summary>
    public SequenceLogo(IEnumerable<Sequence> alignedSequences) : this(alignedSequences.Select(s => s.ToString()).ToList(), DetectNucleotide(alignedSequences)) { }

    /// <summary>A logo from aligned strings; the alphabet (nucleotide or protein) is detected when not given.</summary>
    public SequenceLogo(IReadOnlyList<string> alignedSequences, bool? nucleotide = null) : base(700, 240)
    {
        if (alignedSequences.Count == 0) throw new ArgumentException("No sequences.", nameof(alignedSequences));
        int len = alignedSequences[0].Length;
        if (alignedSequences.Any(s => s.Length != len)) throw new ArgumentException("Aligned sequences must have equal length.");
        IsNucleotide = nucleotide ?? alignedSequences.All(s => s.All(ch => "ACGTUNacgtun-.".Contains(ch)));
        Symbols = IsNucleotide
            ? (alignedSequences.Any(s => s.Contains('U') || s.Contains('u')) ? "ACGU" : "ACGT")
            : "ACDEFGHIKLMNPQRSTVWY";
        _counts = new double[len, Symbols.Length];
        foreach (var s in alignedSequences)
            for (int i = 0; i < len; i++)
            {
                int k = Symbols.IndexOf(char.ToUpperInvariant(s[i]));
                if (k >= 0) _counts[i, k]++;
            }
        SequenceCount = alignedSequences.Count;
        Width = Math.Clamp(80 + len * 34, 240, 1600);
    }

    /// <summary>
    /// A logo from a position × symbol count or frequency matrix whose columns follow <paramref name="symbols"/>
    /// (e.g. <c>"ACGT"</c>). <paramref name="sampleSize"/> is used for the small-sample correction.
    /// </summary>
    public static SequenceLogo FromMatrix(double[,] matrix, string symbols, int sampleSize = 0)
        => new(matrix, symbols, sampleSize);

    private SequenceLogo(double[,] matrix, string symbols, int sampleSize) : base(700, 240)
    {
        if (matrix.GetLength(1) != symbols.Length) throw new ArgumentException("Matrix columns must match the symbols.");
        _counts = (double[,])matrix.Clone();
        Symbols = symbols.ToUpperInvariant();
        IsNucleotide = Symbols.Length <= 4 || Symbols.All(c => "ACGTUN".Contains(c));
        SequenceCount = sampleSize;
        Width = Math.Clamp(80 + matrix.GetLength(0) * 34, 240, 1600);
    }

    /// <summary>Alphabet symbols, in matrix column order.</summary>
    public string Symbols { get; }

    /// <summary>Nucleotide (2-bit maximum) or protein (log2 20 ≈ 4.32 bits) logo.</summary>
    public bool IsNucleotide { get; }

    /// <summary>Number of sequences (0 when built from frequencies without a sample size).</summary>
    public int SequenceCount { get; }

    /// <summary>Number of positions.</summary>
    public int Length => _counts.GetLength(0);

    /// <summary>Height scale.</summary>
    public LogoMode Mode { get; set; } = LogoMode.Bits;

    /// <summary>Subtract the small-sample error e(n) = (s − 1) / (2 ln 2 · n) from each column's information.</summary>
    public bool SmallSampleCorrection { get; set; }

    /// <summary>Number printed under the first column.</summary>
    public int FirstPosition { get; set; } = 1;

    /// <summary>Maximum information per column (log2 of the alphabet size).</summary>
    public double MaxBits => Math.Log2(Symbols.Length);

    /// <summary>Letter stacks of every column.</summary>
    public IReadOnlyList<LogoColumn> ComputeColumns()
    {
        var cols = new List<LogoColumn>(Length);
        int s = Symbols.Length;
        for (int i = 0; i < Length; i++)
        {
            double total = 0;
            for (int k = 0; k < s; k++) total += Math.Max(0, _counts[i, k]);
            var freq = new double[s];
            if (total > 0) for (int k = 0; k < s; k++) freq[k] = Math.Max(0, _counts[i, k]) / total;
            double h = 0;
            foreach (double f in freq) if (f > 0) h -= f * Math.Log2(f);
            double info = total > 0 ? MaxBits - h : 0;
            int n = SequenceCount > 0 ? SequenceCount : (int)Math.Round(total);
            if (SmallSampleCorrection && n > 0) info -= (s - 1) / (2 * Math.Log(2) * n);
            info = Math.Max(0, info);
            var letters = Enumerable.Range(0, s)
                .Where(k => freq[k] > 0)
                .Select(k => new LogoLetter(Symbols[k], freq[k], Mode == LogoMode.Bits ? freq[k] * info : freq[k]))
                .OrderBy(l => l.Height).ThenByDescending(l => l.Symbol)
                .ToList();
            cols.Add(new LogoColumn(FirstPosition + i, info, letters));
        }
        return cols;
    }

    /// <inheritdoc />
    public override void Render(IRenderer r)
    {
        var th = Theme;
        r.Clear(th.Background);
        double top = 14 + DrawTitle(r, 10);
        double left = 56, right = 14, bottom = 34;
        var area = new Rect(left, top, Width - left - right, Height - top - bottom);
        double yMax = Mode == LogoMode.Bits ? MaxBits : 1;

        // Axis.
        Paint.Line(r, area.X - 6, area.Y, area.X - 6, area.Bottom, th.Border);
        foreach (double t in Paint.NiceTicks(0, yMax, 4))
        {
            double y = area.Bottom - t / yMax * area.Height;
            Paint.Line(r, area.X - 10, y, area.X - 6, y, th.Border);
            r.DrawText(t.ToString("0.#", System.Globalization.CultureInfo.InvariantCulture), new Pt(area.X - 12, y), th.SmallFont, th.Muted, TextAlign.Right, TextBaseline.Middle);
        }
        r.DrawText(Mode == LogoMode.Bits ? "bits" : "probability", new Pt(16, area.Y + area.Height / 2), th.SmallFont, th.Muted, TextAlign.Center, TextBaseline.Middle, -90);
        Paint.Line(r, area.X, area.Bottom + 0.5, area.Right, area.Bottom + 0.5, th.Border);

        var cols = ComputeColumns();
        double cw = area.Width / Math.Max(1, cols.Count);
        double gap = Math.Min(3, cw * 0.08);
        for (int i = 0; i < cols.Count; i++)
        {
            double x = area.X + i * cw + gap / 2, y = area.Bottom;
            foreach (var l in cols[i].Letters)
            {
                double h = l.Height / yMax * area.Height;
                if (h < 0.25) continue;
                y -= h;
                LetterGlyphs.Fill(r, l.Symbol, new Rect(x, y, cw - gap, h), th.ResidueColor(l.Symbol, IsNucleotide));
            }
            if (cols.Count <= 60 || i % 5 == 0)
                r.DrawText(cols[i].Position.ToString(System.Globalization.CultureInfo.InvariantCulture), new Pt(area.X + (i + 0.5) * cw, area.Bottom + 6), th.SmallFont, th.Muted, TextAlign.Center, TextBaseline.Top);
        }
    }

    private static bool? DetectNucleotide(IEnumerable<Sequence> sequences)
    {
        var first = sequences.FirstOrDefault();
        return first is null ? null : first.Alphabet.IsNucleotide;
    }
}

/// <summary>A dot match of a <see cref="DotPlot"/>: window starting at <c>X</c> in the first sequence matches the window at <c>Y</c> in the second.</summary>
public readonly record struct DotMatch(int X, int Y, bool ReverseComplement);

/// <summary>
/// A dot plot of two sequences: a dot where a window of <see cref="Window"/> residues has at least
/// <see cref="MinMatches"/> identities; nucleotide plots also show reverse-complement matches (inversions).
/// Computed with running diagonal sums in O(n·m).
/// </summary>
public sealed class DotPlot : FigurePlot
{
    public DotPlot(Sequence first, Sequence second, int window = 10, int? minMatches = null) : base(640, 640)
    {
        First = first ?? throw new ArgumentNullException(nameof(first));
        Second = second ?? throw new ArgumentNullException(nameof(second));
        Window = Math.Max(1, window);
        MinMatches = Math.Clamp(minMatches ?? Window, 1, Window);
        IncludeReverseComplement = first.Alphabet.IsNucleotide && second.Alphabet.IsNucleotide;
        Title = $"{first.Id} vs {second.Id}";
    }

    public Sequence First { get; }
    public Sequence Second { get; }

    /// <summary>Window length.</summary>
    public int Window { get; set; }

    /// <summary>Identities required in a window.</summary>
    public int MinMatches { get; set; }

    /// <summary>Also plot matches against the reverse complement of the second sequence.</summary>
    public bool IncludeReverseComplement { get; set; }

    /// <summary>Cap on the number of dots (the densest plots are thinned evenly).</summary>
    public int MaxDots { get; set; } = 250_000;

    /// <summary>All window matches.</summary>
    public IReadOnlyList<DotMatch> FindMatches()
    {
        var result = new List<DotMatch>();
        Scan(First.Span, Second.Span, false, result);
        if (IncludeReverseComplement) Scan(First.Span, Second.ReverseComplement().Span, true, result);
        return result;
    }

    private void Scan(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, bool reverse, List<DotMatch> sink)
    {
        int n = a.Length, m = b.Length, w = Window;
        if (n < w || m < w) return;
        // Each diagonal d = j - i, walk with a running window sum.
        for (int d = -(n - w); d <= m - w; d++)
        {
            int i0 = Math.Max(0, -d), j0 = i0 + d;
            int len = Math.Min(n - i0, m - j0);
            if (len < w) continue;
            int sum = 0;
            for (int k = 0; k < len; k++)
            {
                if (Same(a[i0 + k], b[j0 + k])) sum++;
                if (k >= w && Same(a[i0 + k - w], b[j0 + k - w])) sum--;
                if (k >= w - 1 && sum >= MinMatches)
                {
                    int i = i0 + k - w + 1, j = j0 + k - w + 1;
                    sink.Add(new DotMatch(i, reverse ? m - j - w : j, reverse));
                    if (sink.Count >= MaxDots * 4) return;
                }
            }
        }
    }

    private static bool Same(byte x, byte y) => (x | 0x20) == (y | 0x20) && x != (byte)'N' && x != (byte)'n' && x != (byte)'-';

    protected override void Build(Figure figure)
    {
        var ax = figure.AddSubplot();
        var matches = FindMatches();
        int stride = Math.Max(1, matches.Count / MaxDots);
        foreach (bool rev in new[] { false, true })
        {
            var pts = matches.Where((mt, k) => mt.ReverseComplement == rev && k % stride == 0).ToList();
            if (pts.Count == 0) continue;
            ax.Scatter(pts.Select(p => (double)p.X + 1).ToArray(), pts.Select(p => (double)p.Y + 1).ToArray(),
                rev ? "reverse complement" : "forward", rev ? Theme.Up : Theme.Accent,
                sizes: Sizes(pts.Count, Math.Clamp(6000.0 / Math.Max(First.Length, Second.Length), 2, 30)), marker: MarkerStyle.Square, alpha: 0.9, edgeWidth: 0);
        }
        ax.SetXLim(0, First.Length + 1);
        ax.SetYLim(0, Second.Length + 1);
        ax.InvertYAxis();
        ax.SetXLabel(XLabel ?? $"{First.Id} (position)");
        ax.SetYLabel(YLabel ?? $"{Second.Id} (position)");
        if (Title is not null) ax.SetTitle(Title);
        if (IncludeReverseComplement && matches.Any(m => m.ReverseComplement)) ax.Legend(LegendLocation.Best);
    }

    /// <inheritdoc />
    public override bool SupportsInteractive => true;

    /// <inheritdoc />
    public override PlotlyFigure ToPlotly()
    {
        var f = NewPlotly();
        var matches = FindMatches();
        int stride = Math.Max(1, matches.Count / MaxDots);
        foreach (bool rev in new[] { false, true })
        {
            var pts = matches.Where((mt, k) => mt.ReverseComplement == rev && k % stride == 0).ToList();
            if (pts.Count == 0) continue;
            f.Add(new Scatter(pts.Count > 5000)
            {
                Name = rev ? "reverse complement" : "forward",
                X = pts.Select(p => (double)p.X + 1).ToArray(), Y = pts.Select(p => (double)p.Y + 1).ToArray(),
                Mode = ScatterMode.Markers,
                Marker = new Marker { Color = rev ? Theme.Up : Theme.Accent, Size = 3, Symbol = "square" },
            });
        }
        f.Layout.XAxis.Title = $"{First.Id} (position)";
        f.Layout.YAxis.Title = $"{Second.Id} (position)";
        f.Layout.YAxis.AutoRangeReversed = true;
        return f;
    }
}
