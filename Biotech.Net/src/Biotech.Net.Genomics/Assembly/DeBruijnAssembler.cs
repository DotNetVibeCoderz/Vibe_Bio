using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Assembly;

/// <summary>Options for <see cref="DeBruijnAssembler"/>.</summary>
public sealed record AssemblyOptions
{
    /// <summary>k-mer length (odd, 3–31; default 31).</summary>
    public int K { get; init; } = 31;

    /// <summary>Minimum occurrences for a k-mer to be kept (solid k-mers; default 2).</summary>
    public int MinKmerCount { get; init; } = 2;

    /// <summary>Clip tips (default true).</summary>
    public bool RemoveTips { get; init; } = true;

    /// <summary>Maximum tip length in k-mers (default 2k).</summary>
    public int? MaxTipKmers { get; init; }

    /// <summary>Pop bubbles (default true).</summary>
    public bool PopBubbles { get; init; } = true;

    /// <summary>Maximum bubble-arm length in k-mers (default 3k).</summary>
    public int? MaxBubbleKmers { get; init; }

    /// <summary>Maximum clean-up rounds (tips then bubbles; default 4).</summary>
    public int CleaningRounds { get; init; } = 4;

    /// <summary>Contigs shorter than this are dropped (default 2k).</summary>
    public int? MinContigLength { get; init; }
}

/// <summary>Contiguity statistics of an assembly.</summary>
public sealed record AssemblyStats(int Count, long TotalLength, long Longest, long Shortest, double MeanLength,
    long N50, int L50, long N90, int L90, double GcContent)
{
    /// <summary>Statistics of a set of contigs.</summary>
    public static AssemblyStats Compute(IEnumerable<Sequence> contigs)
    {
        var list = contigs.ToList();
        long gc = 0, acgt = 0;
        foreach (var c in list)
        {
            var n = c.NucleotideCounts();
            gc += n.G + n.C;
            acgt += n.A + n.C + n.G + n.T;
        }
        var s = Compute(list.Select(c => (long)c.Length));
        return s with { GcContent = acgt == 0 ? double.NaN : (double)gc / acgt };
    }

    /// <summary>Statistics from contig lengths (GC content NaN).</summary>
    public static AssemblyStats Compute(IEnumerable<long> lengths)
    {
        var l = lengths.OrderDescending().ToArray();
        if (l.Length == 0) return new(0, 0, 0, 0, 0, 0, 0, 0, 0, double.NaN);
        long total = l.Sum();
        var (n50, l50) = Nx(l, total, 0.5);
        var (n90, l90) = Nx(l, total, 0.9);
        return new(l.Length, total, l[0], l[^1], (double)total / l.Length, n50, l50, n90, l90, double.NaN);
    }

    private static (long N, int L) Nx(long[] sortedDesc, long total, double fraction)
    {
        long acc = 0;
        for (int i = 0; i < sortedDesc.Length; i++)
        {
            acc += sortedDesc[i];
            if (acc >= total * fraction) return (sortedDesc[i], i + 1);
        }
        return (sortedDesc[^1], sortedDesc.Length);
    }

    public override string ToString()
        => $"{Count} contigs, {TotalLength:N0} bp, longest {Longest:N0}, N50 {N50:N0} (L50 {L50}), N90 {N90:N0}, GC {GcContent:P1}";
}

/// <summary>Result of a de Bruijn assembly.</summary>
public sealed record AssemblyResult(IReadOnlyList<Sequence> Contigs, AssemblyStats Stats, int SolidKmers, int TipsRemoved, int BubblesPopped);

/// <summary>
/// Short-read de Bruijn graph assembler: solid k-mer counting → graph → iterative tip clipping and bubble popping →
/// unitigs as contigs.
/// </summary>
public static class DeBruijnAssembler
{
    /// <summary>Assembles <paramref name="reads"/>.</summary>
    public static AssemblyResult Assemble(IEnumerable<Sequence> reads, AssemblyOptions? options = null)
    {
        options ??= new AssemblyOptions();
        var graph = DeBruijnGraph.Build(reads, options.K, options.MinKmerCount);
        int solid = graph.KmerCount, tips = 0, bubbles = 0;
        int maxTip = options.MaxTipKmers ?? 2 * options.K;
        int maxBubble = options.MaxBubbleKmers ?? 3 * options.K;
        for (int round = 0; round < options.CleaningRounds; round++)
        {
            int t = options.RemoveTips ? graph.RemoveTips(maxTip) : 0;
            int b = options.PopBubbles ? graph.PopBubbles(maxBubble) : 0;
            tips += t; bubbles += b;
            if (t == 0 && b == 0) break;
        }
        int minLen = options.MinContigLength ?? 2 * options.K;
        var contigs = new List<Sequence>();
        foreach (var u in graph.Unitigs())
        {
            if (u.Length < minLen) continue;
            contigs.Add(u.Sequence.WithId($"contig_{contigs.Count + 1}", $"length={u.Length} coverage={u.MeanCoverage:F1}"));
        }
        return new AssemblyResult(contigs, AssemblyStats.Compute(contigs), solid, tips, bubbles);
    }

    /// <summary>Assembles reads (qualities ignored).</summary>
    public static AssemblyResult Assemble(IEnumerable<SequenceRead> reads, AssemblyOptions? options = null)
        => Assemble(reads.Select(r => r.Sequence), options);
}
