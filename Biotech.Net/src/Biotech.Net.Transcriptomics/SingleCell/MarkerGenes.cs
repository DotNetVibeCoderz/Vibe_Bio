using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Statistical test used by <see cref="MarkerGenes"/>.</summary>
public enum MarkerTest
{
    /// <summary>Welch's t-test, one group against the rest (scanpy <c>method="t-test"</c>).</summary>
    TTest,

    /// <summary>Wilcoxon rank-sum test with normal approximation, no tie correction (scanpy <c>method="wilcoxon"</c>).</summary>
    Wilcoxon,
}

/// <summary>One gene's statistics for one group in a marker-gene ranking.</summary>
/// <param name="Gene">Gene symbol.</param>
/// <param name="GeneIndex">Column index of the gene in <see cref="SingleCellExperiment.X"/>.</param>
/// <param name="Score">Test statistic (Welch t or Wilcoxon z); positive means higher in the group.</param>
/// <param name="Log2FoldChange">scanpy's <c>log2((expm1(mean_group) + 1e-9) / (expm1(mean_rest) + 1e-9))</c> on log1p data.</param>
/// <param name="PValue">Two-sided p-value.</param>
/// <param name="PAdj">Benjamini–Hochberg adjusted p-value over all genes of this group.</param>
/// <param name="PctInGroup">Fraction of cells of the group expressing the gene.</param>
/// <param name="PctRest">Fraction of the other cells expressing the gene.</param>
public readonly record struct MarkerGene(string Gene, int GeneIndex, double Score, double Log2FoldChange, double PValue,
    double PAdj, double PctInGroup, double PctRest);

/// <summary>
/// Marker genes of cell clusters — scanpy's <c>tl.rank_genes_groups</c> with <c>reference="rest"</c>: each group
/// is compared with all other cells on the (log-normalised) <see cref="SingleCellExperiment.X"/>.
/// </summary>
public static class MarkerGenes
{
    /// <summary>
    /// Ranks genes for every group of <paramref name="labels"/> (one label per cell). Returns, per group, all genes
    /// sorted by decreasing score.
    /// </summary>
    public static Dictionary<string, List<MarkerGene>> Rank(SingleCellExperiment sc, IReadOnlyList<string> labels,
        MarkerTest test = MarkerTest.Wilcoxon)
    {
        ArgumentNullException.ThrowIfNull(sc);
        ArgumentNullException.ThrowIfNull(labels);
        int n = sc.CellCount, g = sc.GeneCount;
        if (labels.Count != n) throw new ArgumentException($"Need {n} labels, got {labels.Count}.");
        var groups = labels.Distinct().Order(StringComparer.Ordinal).ToArray();
        if (groups.Length < 2) throw new ArgumentException("Need at least two groups.");
        var groupIndex = groups.Select((name, i) => (name, i)).ToDictionary(t => t.name, t => t.i, StringComparer.Ordinal);
        var cellGroup = labels.Select(l => groupIndex[l]).ToArray();
        var groupSize = new int[groups.Length];
        foreach (int c in cellGroup) groupSize[c]++;

        // Per gene (CSC), per group: mean, variance, detection fraction and (Wilcoxon) rank sum.
        var csc = sc.X.Transpose();
        int ng = groups.Length;
        var score = new double[ng * g]; var lfc = new double[ng * g]; var pv = new double[ng * g];
        var pctIn = new double[ng * g]; var pctOut = new double[ng * g];
        Parallel.For(0, g, j =>
        {
            var column = new double[n];
            var idx = csc.RowIndices(j); var v = csc.RowValues(j);
            for (int k = 0; k < idx.Length; k++) column[idx[k]] = v[k];
            var sum = new double[ng]; var sq = new double[ng]; var nz = new int[ng];
            for (int i = 0; i < n; i++)
            {
                int c = cellGroup[i];
                sum[c] += column[i]; sq[c] += column[i] * column[i];
                if (column[i] != 0) nz[c]++;
            }
            double totSum = sum.Sum(), totSq = sq.Sum();
            int totNz = nz.Sum();
            double[]? rankSum = null;
            if (test == MarkerTest.Wilcoxon)
            {
                var ranks = NumericUtil.Rank(column);
                rankSum = new double[ng];
                for (int i = 0; i < n; i++) rankSum[cellGroup[i]] += ranks[i];
            }
            for (int c = 0; c < ng; c++)
            {
                int n1 = groupSize[c], n2 = n - n1;
                double m1 = sum[c] / n1, m2 = (totSum - sum[c]) / n2;
                int at = c * g + j;
                pctIn[at] = (double)nz[c] / n1;
                pctOut[at] = (double)(totNz - nz[c]) / n2;
                lfc[at] = Math.Log2((Math.Exp(m1) - 1 + 1e-9) / (Math.Exp(m2) - 1 + 1e-9));
                if (test == MarkerTest.TTest)
                {
                    double v1 = n1 > 1 ? Math.Max(0, (sq[c] - n1 * m1 * m1) / (n1 - 1)) : 0;
                    double v2 = n2 > 1 ? Math.Max(0, (totSq - sq[c] - n2 * m2 * m2) / (n2 - 1)) : 0;
                    double a = v1 / n1, b = v2 / n2;
                    double se = Math.Sqrt(a + b);
                    double t = se > 0 ? (m1 - m2) / se : double.NaN;
                    double df = (a + b) * (a + b) / (a * a / Math.Max(n1 - 1, 1) + b * b / Math.Max(n2 - 1, 1));
                    // scanpy sets NaN scores (constant genes) to 0 and their p-values to 1.
                    score[at] = double.IsNaN(t) ? 0 : t;
                    pv[at] = double.IsNaN(t) ? 1 : SpecialFunctions.StudentTTwoSided(t, df);
                }
                else
                {
                    double expected = n1 * (n + 1.0) / 2;
                    double sd = Math.Sqrt(n1 * (double)n2 * (n + 1.0) / 12);
                    double z = (rankSum![c] - expected) / sd;
                    score[at] = z;
                    pv[at] = 2 * NumericUtil.NormalUpper(Math.Abs(z));
                }
            }
        });

        var result = new Dictionary<string, List<MarkerGene>>(StringComparer.Ordinal);
        for (int c = 0; c < ng; c++)
        {
            var p = pv.AsSpan(c * g, g).ToArray();
            var q = DESeq.Bh(p);
            var list = new List<MarkerGene>(g);
            for (int j = 0; j < g; j++)
            {
                int at = c * g + j;
                list.Add(new MarkerGene(sc.GeneNames[j], j, score[at], lfc[at], pv[at], q[j], pctIn[at], pctOut[at]));
            }
            list.Sort((x, y) => y.Score.CompareTo(x.Score));
            result[groups[c]] = list;
        }
        return result;
    }

    /// <summary>Ranks genes for the groups stored in <c>CellLabels[key]</c> (e.g. <c>leiden</c>).</summary>
    public static Dictionary<string, List<MarkerGene>> Rank(SingleCellExperiment sc, string key = "leiden",
        MarkerTest test = MarkerTest.Wilcoxon)
    {
        ArgumentNullException.ThrowIfNull(sc);
        if (!sc.CellLabels.TryGetValue(key, out var labels)) throw new KeyNotFoundException($"No cell labels '{key}'.");
        return Rank(sc, labels, test);
    }
}
