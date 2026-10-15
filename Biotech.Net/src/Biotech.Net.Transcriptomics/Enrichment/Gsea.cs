using Biotech.Net.Native;

namespace Biotech.Net.Transcriptomics.Enrichment;

/// <summary>Result of preranked GSEA for one gene set.</summary>
/// <param name="Set">Gene-set name.</param>
/// <param name="Size">Members present in the ranked list.</param>
/// <param name="EnrichmentScore">Maximum signed deviation of the running sum (ES).</param>
/// <param name="NormalizedEnrichmentScore">ES divided by the mean of same-signed null ES (NES).</param>
/// <param name="PValue">Permutation p-value: (#null with same sign and |ES| ≥ |ES_obs| + 1) / (#same-sign null + 1).</param>
/// <param name="PAdj">Benjamini–Hochberg adjusted p-value across sets (as fgsea reports).</param>
/// <param name="FdrQ">GSEA-style FDR q-value computed from the NES null distribution (Subramanian et al. 2005).</param>
/// <param name="LeadingEdge">Genes contributing to the ES (before the peak for positive ES, after it for negative).</param>
public sealed record GseaResult(string Set, int Size, double EnrichmentScore, double NormalizedEnrichmentScore,
    double PValue, double PAdj, double FdrQ, IReadOnlyList<string> LeadingEdge);

/// <summary>
/// Preranked gene-set enrichment analysis (Subramanian et al. 2005): a weighted Kolmogorov–Smirnov-like
/// running sum over a ranked gene list, with significance from gene-set permutations (random sets of the same
/// size, as in GSEA-preranked and fgsea's simple method). Permutations use a seeded RNG, so results are reproducible.
/// </summary>
public static class Gsea
{
    /// <summary>
    /// Runs preranked GSEA. <paramref name="ranking"/> maps genes to a score (e.g. Wald stat or signed −log10 p);
    /// genes are ranked by decreasing score. <paramref name="weight"/> is the exponent p of the hit increments |r|^p
    /// (1 = classic weighted GSEA, 0 = unweighted KS).
    /// </summary>
    public static List<GseaResult> Run(IReadOnlyDictionary<string, double> ranking, IEnumerable<GeneSet> sets,
        int permutations = 1000, double weight = 1, int minSize = 15, int maxSize = 500, int seed = 42)
    {
        ArgumentNullException.ThrowIfNull(ranking);
        ArgumentNullException.ThrowIfNull(sets);
        var genes = ranking.Keys.ToArray();
        var scores = genes.Select(g => ranking[g]).ToArray();
        // Sort by decreasing score; break ties by gene name for determinism.
        var order = Enumerable.Range(0, genes.Length).ToArray();
        Array.Sort(order, (a, b) =>
        {
            int c = scores[b].CompareTo(scores[a]);
            return c != 0 ? c : string.CompareOrdinal(genes[a], genes[b]);
        });
        int n = genes.Length;
        var rankedGenes = new string[n];
        var absW = new double[n];
        var position = new Dictionary<string, int>(n, StringComparer.Ordinal);
        for (int r = 0; r < n; r++)
        {
            rankedGenes[r] = genes[order[r]];
            absW[r] = Math.Pow(Math.Abs(scores[order[r]]), weight);
            position[rankedGenes[r]] = r;
        }

        var tested = new List<(GeneSet Set, int[] Pos)>();
        foreach (var s in sets)
        {
            var pos = new List<int>();
            foreach (var g in s.Genes) if (position.TryGetValue(g, out int p)) pos.Add(p);
            if (pos.Count < minSize || pos.Count > maxSize || pos.Count >= n) continue;
            pos.Sort();
            tested.Add((s, pos.ToArray()));
        }

        int t = tested.Count;
        var es = new double[t]; var peak = new int[t];
        var nullEs = new double[t][];
        Parallel.For(0, t, k =>
        {
            var pos = tested[k].Pos;
            es[k] = EnrichmentScore(pos, absW, n, out peak[k]);
            var rng = new Random(unchecked(seed * 7919 + k * 104729));
            var nul = new double[permutations];
            var sample = new int[pos.Length];
            var chosen = new HashSet<int>();
            for (int b = 0; b < permutations; b++)
            {
                chosen.Clear();
                int c = 0;
                while (c < pos.Length) { int r = rng.Next(n); if (chosen.Add(r)) sample[c++] = r; }
                Array.Sort(sample);
                nul[b] = EnrichmentScore(sample, absW, n, out _);
            }
            nullEs[k] = nul;
        });

        var nes = new double[t]; var pv = new double[t];
        var nullNes = new List<double>(t * permutations);
        for (int k = 0; k < t; k++)
        {
            double posMean = 0, negMean = 0; int posN = 0, negN = 0, posGe = 0, negLe = 0;
            foreach (double v in nullEs[k])
            {
                if (v >= 0) { posMean += v; posN++; if (v >= es[k]) posGe++; }
                else { negMean += v; negN++; if (v <= es[k]) negLe++; }
            }
            posMean = posN > 0 ? posMean / posN : double.NaN;
            negMean = negN > 0 ? -negMean / negN : double.NaN;
            if (es[k] >= 0) { nes[k] = es[k] / posMean; pv[k] = (posGe + 1.0) / (posN + 1.0); }
            else { nes[k] = es[k] / negMean; pv[k] = (negLe + 1.0) / (negN + 1.0); }
            foreach (double v in nullEs[k]) nullNes.Add(v >= 0 ? v / posMean : v / negMean);
        }
        var padj = new double[t];
        BioKernels.BenjaminiHochberg(pv, padj);

        // GSEA FDR: (fraction of null NES at least as extreme) / (fraction of observed NES at least as extreme), per sign.
        var nullPos = nullNes.Where(v => v >= 0 && double.IsFinite(v)).Order().ToArray();
        var nullNeg = nullNes.Where(v => v < 0 && double.IsFinite(v)).Order().ToArray();
        var obsPos = nes.Where(v => v >= 0).Order().ToArray();
        var obsNeg = nes.Where(v => v < 0).Order().ToArray();
        var fdr = new double[t];
        for (int k = 0; k < t; k++)
        {
            if (nes[k] >= 0)
            {
                double fNull = nullPos.Length > 0 ? (double)(nullPos.Length - LowerBound(nullPos, nes[k])) / nullPos.Length : 0;
                double fObs = (double)(obsPos.Length - LowerBound(obsPos, nes[k])) / obsPos.Length;
                fdr[k] = Math.Min(1, fNull / fObs);
            }
            else
            {
                double fNull = nullNeg.Length > 0 ? (double)UpperBound(nullNeg, nes[k]) / nullNeg.Length : 0;
                double fObs = (double)UpperBound(obsNeg, nes[k]) / obsNeg.Length;
                fdr[k] = Math.Min(1, fNull / fObs);
            }
        }

        var results = new List<GseaResult>(t);
        for (int k = 0; k < t; k++)
        {
            var pos = tested[k].Pos;
            var le = new List<string>();
            foreach (int p in pos)
                if (es[k] >= 0 ? p <= peak[k] : p >= peak[k]) le.Add(rankedGenes[p]);
            results.Add(new GseaResult(tested[k].Set.Name, pos.Length, es[k], nes[k], pv[k], padj[k], fdr[k], le));
        }
        results.Sort((a, b) => a.PValue != b.PValue ? a.PValue.CompareTo(b.PValue) : Math.Abs(b.NormalizedEnrichmentScore).CompareTo(Math.Abs(a.NormalizedEnrichmentScore)));
        return results;
    }

    /// <summary>
    /// The enrichment score of a set whose members sit at the sorted 0-based ranks <paramref name="hitPositions"/>
    /// of a list of length <paramref name="n"/> with hit weights <paramref name="absWeights"/> (|r|^p in rank order).
    /// Evaluates the running sum only at hits, so the cost is O(set size).
    /// </summary>
    public static double EnrichmentScore(ReadOnlySpan<int> hitPositions, ReadOnlySpan<double> absWeights, int n, out int peakPosition)
    {
        int k = hitPositions.Length;
        double nr = 0;
        foreach (int p in hitPositions) nr += absWeights[p];
        bool unweighted = nr == 0;
        if (unweighted) nr = k;
        double missStep = 1.0 / (n - k);
        double hit = 0, max = 0, min = 0;
        int maxPos = -1, minPos = -1;
        for (int i = 0; i < k; i++)
        {
            int p = hitPositions[i];
            double misses = (p - i) * missStep;
            double before = hit - misses; // running sum at rank p − 1 (lowest point before this hit)
            if (before < min) { min = before; minPos = p; }
            hit += (unweighted ? 1 : absWeights[p]) / nr;
            double after = hit - misses;
            if (after > max) { max = after; maxPos = p; }
        }
        if (max >= -min) { peakPosition = maxPos; return max; }
        peakPosition = minPos;
        return min;
    }

    private static int LowerBound(double[] sorted, double v)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi) { int mid = (lo + hi) >> 1; if (sorted[mid] < v) lo = mid + 1; else hi = mid; }
        return lo;
    }

    private static int UpperBound(double[] sorted, double v)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi) { int mid = (lo + hi) >> 1; if (sorted[mid] <= v) lo = mid + 1; else hi = mid; }
        return lo;
    }
}
