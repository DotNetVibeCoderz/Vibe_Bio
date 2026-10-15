using System.Globalization;
using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics;

/// <summary>One gene of a differential-expression table.</summary>
/// <param name="Gene">Gene identifier.</param>
/// <param name="BaseMean">DESeq2: mean of normalised counts. t-test / limma: average log2 expression (AveExpr).</param>
/// <param name="Log2FoldChange">Estimated log2 fold change (treatment vs reference).</param>
/// <param name="LfcSE">Standard error of the log2 fold change.</param>
/// <param name="Stat">Test statistic (Wald z, Welch/Student t, or moderated t).</param>
/// <param name="PValue">Two-sided p-value.</param>
/// <param name="PAdj">Benjamini–Hochberg adjusted p-value (NaN when filtered or untestable).</param>
public readonly record struct DifferentialExpressionRow(
    string Gene, double BaseMean, double Log2FoldChange, double LfcSE, double Stat, double PValue, double PAdj);

/// <summary>
/// A differential-expression result table in DESeq2's column layout (baseMean, log2FoldChange, lfcSE, stat,
/// pvalue, padj), stored column-wise. Provides top-table selection and MA/volcano plot coordinates.
/// </summary>
public sealed class DifferentialExpressionResult
{
    /// <summary>Creates a result table from parallel columns (not copied).</summary>
    public DifferentialExpressionResult(string method, string contrast, IReadOnlyList<string> genes, double[] baseMean,
        double[] log2FoldChange, double[] lfcSE, double[] stat, double[] pValue, double[] pAdj)
    {
        int n = genes.Count;
        if (baseMean.Length != n || log2FoldChange.Length != n || lfcSE.Length != n || stat.Length != n || pValue.Length != n || pAdj.Length != n)
            throw new ArgumentException("All result columns need one value per gene.");
        Method = method;
        Contrast = contrast;
        Genes = genes;
        BaseMean = baseMean;
        Log2FoldChange = log2FoldChange;
        LfcSE = lfcSE;
        Stat = stat;
        PValue = pValue;
        PAdj = pAdj;
    }

    /// <summary>The test that produced the table (e.g. <c>DESeq2-Wald</c>, <c>limma-voom</c>, <c>Welch</c>).</summary>
    public string Method { get; }

    /// <summary>Description of the tested contrast.</summary>
    public string Contrast { get; }

    /// <summary>Gene identifiers.</summary>
    public IReadOnlyList<string> Genes { get; }

    /// <summary>DESeq2: mean of normalised counts; t-test / limma: average log2 expression.</summary>
    public double[] BaseMean { get; }

    /// <summary>Estimated log2 fold changes.</summary>
    public double[] Log2FoldChange { get; }

    /// <summary>Standard errors of the log2 fold changes.</summary>
    public double[] LfcSE { get; }

    /// <summary>Test statistics.</summary>
    public double[] Stat { get; }

    /// <summary>Two-sided p-values.</summary>
    public double[] PValue { get; }

    /// <summary>Benjamini–Hochberg adjusted p-values (NaN when filtered or untestable).</summary>
    public double[] PAdj { get; }

    /// <summary>Number of genes.</summary>
    public int Count => Genes.Count;

    /// <summary>One gene's row.</summary>
    public DifferentialExpressionRow this[int index]
        => new(Genes[index], BaseMean[index], Log2FoldChange[index], LfcSE[index], Stat[index], PValue[index], PAdj[index]);

    /// <summary>All rows in gene order.</summary>
    public IEnumerable<DifferentialExpressionRow> Rows()
    {
        for (int i = 0; i < Count; i++) yield return this[i];
    }

    /// <summary>The <paramref name="n"/> rows with the smallest p-values (NaN last), like limma's <c>topTable</c>.</summary>
    public IReadOnlyList<DifferentialExpressionRow> Top(int n = 10)
    {
        var order = OrderByPValue();
        var rows = new List<DifferentialExpressionRow>(Math.Min(n, order.Length));
        for (int k = 0; k < order.Length && rows.Count < n; k++) rows.Add(this[order[k]]);
        return rows;
    }

    /// <summary>
    /// Rows with <c>padj &lt; alpha</c> and <c>|log2FC| ≥ minAbsLog2FoldChange</c>, sorted by p-value.
    /// </summary>
    public IReadOnlyList<DifferentialExpressionRow> Significant(double alpha = 0.05, double minAbsLog2FoldChange = 0)
    {
        var rows = new List<DifferentialExpressionRow>();
        foreach (int i in OrderByPValue())
            if (PAdj[i] < alpha && Math.Abs(Log2FoldChange[i]) >= minAbsLog2FoldChange) rows.Add(this[i]);
        return rows;
    }

    /// <summary>Number of genes with <c>padj &lt; alpha</c>, split into up- and down-regulated.</summary>
    public (int Up, int Down) CountSignificant(double alpha = 0.05)
    {
        int up = 0, down = 0;
        for (int i = 0; i < Count; i++)
            if (PAdj[i] < alpha) { if (Log2FoldChange[i] > 0) up++; else if (Log2FoldChange[i] < 0) down++; }
        return (up, down);
    }

    /// <summary>Volcano-plot coordinates: x = log2 fold change, y = −log10 p (genes with NaN p are skipped).</summary>
    public (string[] Genes, double[] X, double[] Y) Volcano()
    {
        var g = new List<string>(); var x = new List<double>(); var y = new List<double>();
        for (int i = 0; i < Count; i++)
        {
            if (double.IsNaN(PValue[i]) || double.IsNaN(Log2FoldChange[i])) continue;
            g.Add(Genes[i]); x.Add(Log2FoldChange[i]); y.Add(-Math.Log10(Math.Max(PValue[i], 1e-300)));
        }
        return (g.ToArray(), x.ToArray(), y.ToArray());
    }

    /// <summary>
    /// MA-plot coordinates: A = log2 base mean for count-based tests (log2(baseMean) of DESeq2), or the average log
    /// expression itself when <see cref="BaseMean"/> already is on the log scale; M = log2 fold change.
    /// </summary>
    public (string[] Genes, double[] A, double[] M) MaPlot(bool baseMeanIsLog = false)
    {
        var g = new List<string>(); var a = new List<double>(); var m = new List<double>();
        for (int i = 0; i < Count; i++)
        {
            if (double.IsNaN(Log2FoldChange[i]) || (!baseMeanIsLog && BaseMean[i] <= 0)) continue;
            g.Add(Genes[i]); a.Add(baseMeanIsLog ? BaseMean[i] : Math.Log2(BaseMean[i])); m.Add(Log2FoldChange[i]);
        }
        return (g.ToArray(), a.ToArray(), m.ToArray());
    }

    /// <summary>The table as a genes × 6 <see cref="Dataset"/>.</summary>
    public Dataset ToDataset()
    {
        var v = new double[Count * 6];
        for (int i = 0; i < Count; i++)
        {
            v[i * 6] = BaseMean[i]; v[i * 6 + 1] = Log2FoldChange[i]; v[i * 6 + 2] = LfcSE[i];
            v[i * 6 + 3] = Stat[i]; v[i * 6 + 4] = PValue[i]; v[i * 6 + 5] = PAdj[i];
        }
        return new Dataset($"{Method}: {Contrast}", Genes, ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], v);
    }

    /// <summary>Writes a TSV with DESeq2-style column names.</summary>
    public void WriteTsv(TextWriter writer)
    {
        ArgumentNullException.ThrowIfNull(writer);
        writer.WriteLine("gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj");
        var ci = CultureInfo.InvariantCulture;
        for (int i = 0; i < Count; i++)
            writer.WriteLine(string.Join('\t', Genes[i], BaseMean[i].ToString("G10", ci), Log2FoldChange[i].ToString("G10", ci),
                LfcSE[i].ToString("G10", ci), Stat[i].ToString("G10", ci), PValue[i].ToString("G10", ci), PAdj[i].ToString("G10", ci)));
    }

    /// <inheritdoc />
    public override string ToString()
    {
        var (up, down) = CountSignificant(0.1);
        return $"{Method} [{Contrast}]: {Count} genes, padj<0.1: {up} up, {down} down";
    }

    private int[] OrderByPValue()
    {
        var idx = new int[Count];
        var keys = new double[Count];
        for (int i = 0; i < Count; i++) { idx[i] = i; keys[i] = double.IsNaN(PValue[i]) ? double.PositiveInfinity : PValue[i]; }
        Array.Sort(keys, idx);
        return idx;
    }
}
