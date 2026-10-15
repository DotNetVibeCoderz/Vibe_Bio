using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>A named library compound.</summary>
/// <param name="Id">Identifier.</param>
/// <param name="Molecule">Structure.</param>
public sealed record LibraryCompound(string Id, Molecule Molecule)
{
    /// <summary>Parses a SMILES into a compound.</summary>
    public static LibraryCompound FromSmiles(string id, string smiles) => new(id, Molecule.FromSmiles(smiles, id));
}

/// <summary>How scores from several fingerprints and/or several query molecules are combined.</summary>
public enum FusionRule
{
    /// <summary>Highest similarity (MAX fusion — the usual best choice for group fusion).</summary>
    Max,
    /// <summary>Arithmetic mean of similarities.</summary>
    Mean,
    /// <summary>Lowest similarity (conservative consensus).</summary>
    Min,
    /// <summary>Sum of ranks (lower is better; reported as negative rank sum so higher is still better).</summary>
    RankSum,
    /// <summary>Reciprocal-rank fusion: Σ 1 / (60 + rank).</summary>
    ReciprocalRank,
}

/// <summary>Options for <see cref="SimilarityScreen"/>.</summary>
public sealed record SimilarityScreenOptions
{
    /// <summary>Fingerprints to use (default ECFP4, FCFP4, MACCS and atom pairs).</summary>
    public IReadOnlyList<FingerprintType> Fingerprints { get; init; } =
        [FingerprintType.Ecfp4, FingerprintType.Fcfp4, FingerprintType.Maccs, FingerprintType.AtomPair];

    /// <summary>Similarity coefficient (default Tanimoto).</summary>
    public SimilarityMetric Metric { get; init; } = SimilarityMetric.Tanimoto;

    /// <summary>Fusion across fingerprints (default Mean).</summary>
    public FusionRule FingerprintFusion { get; init; } = FusionRule.Mean;

    /// <summary>Fusion across query molecules (default Max).</summary>
    public FusionRule QueryFusion { get; init; } = FusionRule.Max;

    /// <summary>Keep at most this many hits (default all).</summary>
    public int TopK { get; init; } = int.MaxValue;
}

/// <summary>One ranked screening result.</summary>
/// <param name="Id">Compound id.</param>
/// <param name="Rank">1-based rank.</param>
/// <param name="Score">Fused score (higher is better).</param>
/// <param name="PerFingerprint">Best similarity to any query, per fingerprint type.</param>
public sealed record ScreeningHit(string Id, int Rank, double Score, IReadOnlyDictionary<FingerprintType, double> PerFingerprint);

/// <summary>
/// Ligand-based virtual screening by fingerprint similarity with data fusion (Willett, "Enhancing the effectiveness
/// of ligand-based virtual screening using data fusion", QSAR Comb. Sci. 2006): several fingerprints per pair are
/// fused into one score, and several query actives are fused by group fusion (MAX by default).
/// </summary>
public static class SimilarityScreen
{
    /// <summary>Ranks <paramref name="library"/> by fused similarity to <paramref name="queries"/>.</summary>
    public static IReadOnlyList<ScreeningHit> Run(IReadOnlyList<Molecule> queries, IReadOnlyList<LibraryCompound> library, SimilarityScreenOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(queries);
        if (queries.Count == 0) throw new ArgumentException("At least one query molecule is required.", nameof(queries));
        options ??= new SimilarityScreenOptions();
        var fpTypes = options.Fingerprints;
        int nQ = queries.Count, nL = library.Count, nF = fpTypes.Count;
        var qfp = new BitVector[nF][];
        for (int f = 0; f < nF; f++) qfp[f] = queries.Select(q => Fingerprint.Compute(q, fpTypes[f])).ToArray();
        // sim[f][q][l]
        var sim = new double[nF][][];
        for (int f = 0; f < nF; f++)
        {
            var lfp = new BitVector[nL];
            int ff = f;
            Parallel.For(0, nL, l => lfp[l] = Fingerprint.Compute(library[l].Molecule, fpTypes[ff]));
            sim[f] = new double[nQ][];
            for (int q = 0; q < nQ; q++) sim[f][q] = Similarity.Bulk(qfp[f][q], lfp, options.Metric);
        }
        // fuse fingerprints for each (q, l)
        var perQuery = new double[nQ][];
        for (int q = 0; q < nQ; q++)
        {
            var lists = Enumerable.Range(0, nF).Select(f => sim[f][q]).ToArray();
            perQuery[q] = Fuse(lists, options.FingerprintFusion);
        }
        var fused = Fuse(perQuery, options.QueryFusion);
        var order = Enumerable.Range(0, nL).OrderByDescending(l => fused[l]).ThenBy(l => l).Take(options.TopK).ToList();
        var hits = new List<ScreeningHit>(order.Count);
        int rank = 0;
        foreach (int l in order)
        {
            var per = new Dictionary<FingerprintType, double>();
            for (int f = 0; f < nF; f++) per[fpTypes[f]] = Enumerable.Range(0, nQ).Max(q => sim[f][q][l]);
            hits.Add(new ScreeningHit(library[l].Id, ++rank, fused[l], per));
        }
        return hits;
    }

    /// <summary>Single-query convenience overload.</summary>
    public static IReadOnlyList<ScreeningHit> Run(Molecule query, IReadOnlyList<LibraryCompound> library, SimilarityScreenOptions? options = null)
        => Run([query], library, options);

    /// <summary>Fuses several score lists (each over the same items) with <paramref name="rule"/>.</summary>
    public static double[] Fuse(IReadOnlyList<double[]> lists, FusionRule rule)
    {
        if (lists.Count == 0) return [];
        int n = lists[0].Length;
        var r = new double[n];
        switch (rule)
        {
            case FusionRule.Max:
                for (int i = 0; i < n; i++) r[i] = lists.Max(l => l[i]);
                break;
            case FusionRule.Min:
                for (int i = 0; i < n; i++) r[i] = lists.Min(l => l[i]);
                break;
            case FusionRule.Mean:
                for (int i = 0; i < n; i++) r[i] = lists.Average(l => l[i]);
                break;
            case FusionRule.RankSum:
            case FusionRule.ReciprocalRank:
                foreach (var l in lists)
                {
                    var ranks = Ranks(l);
                    for (int i = 0; i < n; i++) r[i] += rule == FusionRule.RankSum ? -ranks[i] : 1.0 / (60 + ranks[i]);
                }
                break;
        }
        return r;
    }

    /// <summary>1-based ranks (highest score = 1; ties share the average rank).</summary>
    internal static double[] Ranks(double[] scores)
    {
        int n = scores.Length;
        var idx = Enumerable.Range(0, n).OrderByDescending(i => scores[i]).ToArray();
        var ranks = new double[n];
        int k = 0;
        while (k < n)
        {
            int j = k;
            while (j + 1 < n && scores[idx[j + 1]] == scores[idx[k]]) j++;
            double avg = (k + j) / 2.0 + 1;
            for (int t = k; t <= j; t++) ranks[idx[t]] = avg;
            k = j + 1;
        }
        return ranks;
    }

    /// <summary>
    /// Enrichment factor at fraction <paramref name="fraction"/> of a ranked list:
    /// (actives found in the top x%) / (actives expected at random).
    /// </summary>
    public static double EnrichmentFactor(IReadOnlyList<ScreeningHit> ranked, ISet<string> actives, double fraction = 0.01)
    {
        if (ranked.Count == 0 || actives.Count == 0) return 0;
        int top = Math.Max(1, (int)Math.Ceiling(ranked.Count * fraction));
        int found = ranked.Take(top).Count(h => actives.Contains(h.Id));
        int totalActives = ranked.Count(h => actives.Contains(h.Id));
        if (totalActives == 0) return 0;
        return (double)found / top / ((double)totalActives / ranked.Count);
    }

    /// <summary>Area under the ROC curve of a ranked list (probability an active ranks above a decoy).</summary>
    public static double RocAuc(IReadOnlyList<ScreeningHit> ranked, ISet<string> actives)
    {
        long pos = ranked.Count(h => actives.Contains(h.Id)), neg = ranked.Count - pos;
        if (pos == 0 || neg == 0) return double.NaN;
        long decoysSeen = 0;
        double sum = 0;
        foreach (var h in ranked)
        {
            if (actives.Contains(h.Id)) sum += neg - decoysSeen;
            else decoysSeen++;
        }
        return sum / (pos * neg);
    }
}

