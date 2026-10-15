using Biotech.Net.Native;

namespace Biotech.Net.Transcriptomics.Enrichment;

/// <summary>A named set of genes (pathway, GO term, MSigDB signature).</summary>
public sealed class GeneSet
{
    /// <summary>Creates a gene set; duplicate genes are removed.</summary>
    public GeneSet(string name, IEnumerable<string> genes, string? description = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentNullException.ThrowIfNull(genes);
        Name = name;
        Description = description ?? string.Empty;
        Genes = new HashSet<string>(genes, StringComparer.Ordinal);
    }

    /// <summary>Set name.</summary>
    public string Name { get; }

    /// <summary>Description (second GMT column; often a URL).</summary>
    public string Description { get; }

    /// <summary>Member genes.</summary>
    public IReadOnlySet<string> Genes { get; }

    /// <summary>Number of member genes.</summary>
    public int Count => Genes.Count;

    /// <inheritdoc />
    public override string ToString() => $"{Name} ({Count} genes)";
}

/// <summary>Reader and writer for the Gene Matrix Transposed (GMT) format used by MSigDB and GSEA.</summary>
public static class Gmt
{
    /// <summary>Reads GMT: one set per line, <c>name&lt;TAB&gt;description&lt;TAB&gt;gene1&lt;TAB&gt;gene2...</c>.</summary>
    public static List<GeneSet> Read(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var sets = new List<GeneSet>();
        string? line;
        int lineNo = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (string.IsNullOrWhiteSpace(line) || line[0] == '#') continue;
            var parts = line.Split('\t');
            if (parts.Length < 2) throw new InvalidDataException($"GMT line {lineNo}: expected name, description and genes.");
            var genes = new List<string>(parts.Length - 2);
            for (int i = 2; i < parts.Length; i++)
            {
                var gene = parts[i].Trim();
                if (gene.Length > 0) genes.Add(gene);
            }
            sets.Add(new GeneSet(parts[0].Trim(), genes, parts[1]));
        }
        return sets;
    }

    /// <summary>Reads a GMT file.</summary>
    public static List<GeneSet> ReadFile(string path)
    {
        using var r = CountMatrixReader.OpenText(path);
        return Read(r);
    }

    /// <summary>Writes gene sets in GMT format.</summary>
    public static void Write(TextWriter writer, IEnumerable<GeneSet> sets)
    {
        ArgumentNullException.ThrowIfNull(writer);
        foreach (var s in sets) writer.WriteLine(string.Join('\t', new[] { s.Name, s.Description }.Concat(s.Genes.Order(StringComparer.Ordinal))));
    }
}

/// <summary>One gene set tested for over-representation.</summary>
/// <param name="Set">Gene-set name.</param>
/// <param name="SetSize">Set members present in the universe (K).</param>
/// <param name="Overlap">Hits in the set (k).</param>
/// <param name="Expected">Expected overlap n·K/N.</param>
/// <param name="FoldEnrichment">k / expected.</param>
/// <param name="OddsRatio">Sample odds ratio (k(N−K−n+k)) / ((K−k)(n−k)).</param>
/// <param name="PValue">Hypergeometric upper tail P(X ≥ k).</param>
/// <param name="PAdj">Benjamini–Hochberg adjusted p-value over the tested sets.</param>
/// <param name="Genes">The overlapping genes.</param>
public sealed record OverRepresentationResult(string Set, int SetSize, int Overlap, double Expected, double FoldEnrichment,
    double OddsRatio, double PValue, double PAdj, IReadOnlyList<string> Genes);

/// <summary>
/// Over-representation analysis (ORA): the one-sided hypergeometric test (equivalently, Fisher's exact test,
/// greater) of a hit list against each gene set within a universe of measured genes, as in clusterProfiler's
/// <c>enricher</c> or R's <c>phyper(k − 1, K, N − K, n, lower.tail = FALSE)</c>.
/// </summary>
public static class OverRepresentation
{
    /// <summary>
    /// Tests every set with at least <paramref name="minSetSize"/> and at most <paramref name="maxSetSize"/> members in
    /// the universe. Hits outside the universe are ignored; results are sorted by p-value.
    /// </summary>
    public static List<OverRepresentationResult> Test(IEnumerable<string> hits, IEnumerable<string> universe,
        IEnumerable<GeneSet> sets, int minSetSize = 1, int maxSetSize = int.MaxValue)
    {
        ArgumentNullException.ThrowIfNull(hits);
        ArgumentNullException.ThrowIfNull(universe);
        ArgumentNullException.ThrowIfNull(sets);
        var u = new HashSet<string>(universe, StringComparer.Ordinal);
        var h = new HashSet<string>(StringComparer.Ordinal);
        foreach (var g in hits) if (u.Contains(g)) h.Add(g);
        int bigN = u.Count, n = h.Count;
        var raw = new List<(GeneSet Set, int K, List<string> Overlap, double P)>();
        foreach (var s in sets)
        {
            int k = 0, bigK = 0;
            var ov = new List<string>();
            foreach (var g in s.Genes)
            {
                if (!u.Contains(g)) continue;
                bigK++;
                if (h.Contains(g)) { k++; ov.Add(g); }
            }
            if (bigK < minSetSize || bigK > maxSetSize) continue;
            ov.Sort(StringComparer.Ordinal);
            raw.Add((s, bigK, ov, HypergeometricUpperTail(k, bigK, bigN, n)));
        }
        var pv = raw.Select(r => r.P).ToArray();
        var q = new double[pv.Length];
        BioKernels.BenjaminiHochberg(pv, q);
        var results = new List<OverRepresentationResult>(raw.Count);
        for (int i = 0; i < raw.Count; i++)
        {
            var (s, bigK, ov, p) = raw[i];
            int k = ov.Count;
            double expected = (double)n * bigK / bigN;
            double or = (double)k * (bigN - bigK - n + k) / ((double)(bigK - k) * (n - k));
            results.Add(new OverRepresentationResult(s.Name, bigK, k, expected, expected > 0 ? k / expected : double.NaN, or, p, q[i], ov));
        }
        results.Sort((a, b) => a.PValue.CompareTo(b.PValue));
        return results;
    }

    /// <summary>
    /// P(X ≥ k) for X ~ Hypergeometric(population N, successes K, draws n), summed in log space for stability.
    /// </summary>
    public static double HypergeometricUpperTail(int k, int successes, int population, int draws)
    {
        int bigK = successes, bigN = population, n = draws;
        int lo = Math.Max(0, n - (bigN - bigK)), hi = Math.Min(n, bigK);
        if (k <= lo) return 1;
        if (k > hi) return 0;
        double logDenom = LogChoose(bigN, n);
        double max = double.NegativeInfinity;
        var terms = new double[hi - k + 1];
        for (int x = k; x <= hi; x++)
        {
            double t = LogChoose(bigK, x) + LogChoose(bigN - bigK, n - x) - logDenom;
            terms[x - k] = t;
            max = Math.Max(max, t);
        }
        double s = 0;
        foreach (double t in terms) s += Math.Exp(t - max);
        return Math.Min(1, Math.Exp(max) * s);
    }

    private static double LogChoose(int n, int k)
        => SpecialFunctions.LogFactorial(n) - SpecialFunctions.LogFactorial(k) - SpecialFunctions.LogFactorial(n - k);
}
