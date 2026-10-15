using System.Buffers.Binary;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Biotech.Net.Core;
using Biotech.Net.Sequences;

namespace Biotech.Net.Genomics.Assembly;

/// <summary>A maximal non-branching path of the de Bruijn graph.</summary>
/// <param name="Sequence">The unitig bases.</param>
/// <param name="KmerCount">Number of k-mers on the path.</param>
/// <param name="MeanCoverage">Mean k-mer count.</param>
public sealed record Unitig(Sequence Sequence, int KmerCount, double MeanCoverage)
{
    /// <summary>Length in bases.</summary>
    public int Length => Sequence.Length;
}

/// <summary>
/// A bidirected de Bruijn graph over canonical k-mers (odd k ≤ 31, 2-bit packed). Edges are implicit: two k-mers are
/// adjacent when they overlap by k−1 bases in some orientation. Supports solid-k-mer filtering, compaction into
/// unitigs, tip clipping and bubble popping.
/// </summary>
public sealed class DeBruijnGraph
{
    private readonly Dictionary<ulong, int> _counts;
    private readonly ulong _mask;

    private DeBruijnGraph(int k, Dictionary<ulong, int> counts)
    {
        K = k;
        _counts = counts;
        _mask = (1UL << (2 * k)) - 1;
    }

    /// <summary>k-mer length.</summary>
    public int K { get; }

    /// <summary>Number of (canonical) k-mers in the graph.</summary>
    public int KmerCount => _counts.Count;

    /// <summary>Count of a k-mer (either orientation), 0 when absent.</summary>
    public int Count(string kmer)
        => kmer.Length == K && Kmer.TryEncode(kmer, out ulong code) ? _counts.GetValueOrDefault(Canonical(code)) : 0;

    /// <summary>
    /// Counts the canonical k-mers of <paramref name="reads"/> and keeps those seen at least
    /// <paramref name="minCount"/> times ("solid" k-mers; sequencing errors mostly create singletons).
    /// </summary>
    public static DeBruijnGraph Build(IEnumerable<Sequence> reads, int k, int minCount = 2)
    {
        if (k is < 3 or > 31 || k % 2 == 0) throw new ArgumentOutOfRangeException(nameof(k), "k must be odd and between 3 and 31.");
        var counts = new Dictionary<ulong, int>();
        ulong mask = (1UL << (2 * k)) - 1;
        int shift = 2 * (k - 1);
        foreach (var read in reads)
        {
            var s = read.Span;
            ulong fwd = 0, rev = 0;
            int valid = 0;
            for (int i = 0; i < s.Length; i++)
            {
                int c = Kmer.Code(s[i]);
                if (c < 0) { valid = 0; fwd = rev = 0; continue; }
                fwd = ((fwd << 2) | (uint)c) & mask;
                rev = (rev >> 2) | ((ulong)(3 - c) << shift);
                if (++valid < k) continue;
                CollectionsMarshal.GetValueRefOrAddDefault(counts, Math.Min(fwd, rev), out _)++;
            }
        }
        if (minCount > 1)
        {
            var solid = new Dictionary<ulong, int>(counts.Count / 2);
            foreach (var kv in counts) if (kv.Value >= minCount) solid[kv.Key] = kv.Value;
            counts = solid;
        }
        return new DeBruijnGraph(k, counts);
    }

    /// <summary>k-mer count histogram of the graph (<c>h[m]</c> = k-mers seen m times).</summary>
    public long[] Spectrum()
    {
        int max = 0;
        foreach (int c in _counts.Values) max = Math.Max(max, c);
        var h = new long[max + 1];
        foreach (int c in _counts.Values) h[c]++;
        return h;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong ReverseComplement(ulong x)
    {
        x = ~x;
        x = ((x >> 2) & 0x3333333333333333UL) | ((x & 0x3333333333333333UL) << 2);
        x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FUL) | ((x & 0x0F0F0F0F0F0F0F0FUL) << 4);
        x = BinaryPrimitives.ReverseEndianness(x);
        return x >> (64 - 2 * K);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ulong Canonical(ulong x) => Math.Min(x, ReverseComplement(x));

    private bool Has(ulong oriented) => _counts.ContainsKey(Canonical(oriented));

    private int CountOf(ulong oriented) => _counts.GetValueOrDefault(Canonical(oriented));

    /// <summary>Successors of an oriented k-mer; returns how many, the last one in <paramref name="next"/>.</summary>
    private int OutDegree(ulong x, out ulong next)
    {
        int d = 0;
        next = 0;
        ulong baseCode = (x << 2) & _mask;
        for (uint b = 0; b < 4; b++)
        {
            ulong y = baseCode | b;
            if (Has(y)) { d++; next = y; }
        }
        return d;
    }

    private int InDegree(ulong x, out ulong prev)
    {
        int d = OutDegree(ReverseComplement(x), out ulong p);
        prev = ReverseComplement(p);
        return d;
    }

    private int OutDegree(ulong x) => OutDegree(x, out _);
    private int InDegree(ulong x) => InDegree(x, out _);

    private sealed class Path
    {
        public required List<ulong> Kmers;   // oriented, in path order
        public double Coverage;
    }

    private List<Path> BuildPaths()
    {
        var visited = new HashSet<ulong>(_counts.Count);
        var paths = new List<Path>();
        var keys = _counts.Keys.ToArray();
        Array.Sort(keys);
        foreach (ulong c in keys)
        {
            if (!visited.Add(c)) continue;
            var forward = Extend(c, visited);
            var backward = Extend(ReverseComplement(c), visited);
            var kmers = new List<ulong>(forward.Count + backward.Count + 1);
            for (int i = backward.Count - 1; i >= 0; i--) kmers.Add(ReverseComplement(backward[i]));
            kmers.Add(c);
            kmers.AddRange(forward);
            double cov = 0;
            foreach (var x in kmers) cov += CountOf(x);
            paths.Add(new Path { Kmers = kmers, Coverage = cov / kmers.Count });
        }
        return paths;
    }

    private List<ulong> Extend(ulong start, HashSet<ulong> visited)
    {
        var list = new List<ulong>();
        ulong cur = start;
        while (true)
        {
            if (OutDegree(cur, out ulong next) != 1) break;
            if (InDegree(next) != 1) break;
            if (!visited.Add(Canonical(next))) break; // cycle or already used
            list.Add(next);
            cur = next;
        }
        return list;
    }

    private Sequence ToSequence(List<ulong> kmers, string id)
    {
        var buf = new byte[K + kmers.Count - 1];
        string first = Kmer.Decode(kmers[0], K);
        for (int i = 0; i < K; i++) buf[i] = (byte)first[i];
        for (int i = 1; i < kmers.Count; i++) buf[K + i - 1] = (byte)"ACGT"[(int)(kmers[i] & 3)];
        return Sequence.FromBuffer(id, buf, Alphabet.Dna, validate: false);
    }

    /// <summary>Compacts the graph into unitigs (maximal non-branching paths).</summary>
    public IReadOnlyList<Unitig> Unitigs()
    {
        var paths = BuildPaths();
        var result = new List<Unitig>(paths.Count);
        int n = 0;
        foreach (var p in paths.OrderByDescending(p => p.Kmers.Count))
            result.Add(new Unitig(ToSequence(p.Kmers, $"unitig_{++n}"), p.Kmers.Count, p.Coverage));
        return result;
    }

    private void Remove(Path p)
    {
        foreach (var x in p.Kmers) _counts.Remove(Canonical(x));
    }

    /// <summary>
    /// Clips tips: unitigs of at most <paramref name="maxKmers"/> k-mers that dead-end on one side and hang off a
    /// branching node on the other (typical of sequencing errors near read ends). Lower-coverage tips go first and a
    /// tip is only removed while its junction still has another branch, so genuine sequence ends survive.
    /// </summary>
    /// <returns>Number of tips removed.</returns>
    public int RemoveTips(int maxKmers)
    {
        int removed = 0;
        foreach (var p in BuildPaths().Where(p => p.Kmers.Count <= maxKmers).OrderBy(p => p.Coverage))
        {
            ulong s = p.Kmers[0], e = p.Kmers[^1];
            int inS = InDegree(s, out _), outE = OutDegree(e, out ulong after);
            InDegree(s, out ulong before);
            bool tip;
            if (inS == 0 && outE == 1) tip = InDegree(after) > 1;
            else if (outE == 0 && inS == 1) tip = OutDegree(before) > 1;
            else if (inS == 0 && outE > 1 || outE == 0 && inS > 1) tip = true;
            else tip = false;
            if (!tip) continue;
            Remove(p);
            removed++;
        }
        return removed;
    }

    /// <summary>
    /// Pops simple bubbles: two or more unitigs of at most <paramref name="maxKmers"/> k-mers that leave the same node
    /// and rejoin at the same node (SNP-like errors or heterozygous sites). The highest-coverage path is kept.
    /// </summary>
    /// <returns>Number of paths removed.</returns>
    public int PopBubbles(int maxKmers)
    {
        var groups = new Dictionary<(ulong, ulong), List<Path>>();
        foreach (var p in BuildPaths())
        {
            if (p.Kmers.Count > maxKmers) continue;
            ulong s = p.Kmers[0], e = p.Kmers[^1];
            if (InDegree(s, out ulong pred) != 1 || OutDegree(e, out ulong succ) != 1) continue;
            var key1 = (pred, succ);
            var key2 = (ReverseComplement(succ), ReverseComplement(pred));
            var key = key1.CompareTo(key2) <= 0 ? key1 : key2;
            if (!groups.TryGetValue(key, out var list)) groups[key] = list = [];
            list.Add(p);
        }
        int removed = 0;
        foreach (var list in groups.Values)
        {
            if (list.Count < 2) continue;
            var keep = list.MaxBy(p => p.Coverage)!;
            foreach (var p in list)
            {
                if (ReferenceEquals(p, keep)) continue;
                if (Math.Abs(p.Kmers.Count - keep.Kmers.Count) > Math.Max(2, keep.Kmers.Count / 5)) continue;
                Remove(p);
                removed++;
            }
        }
        return removed;
    }
}
