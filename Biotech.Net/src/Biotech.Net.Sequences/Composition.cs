using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>A value computed over the window <c>[Start, End)</c>.</summary>
public readonly record struct WindowValue(int Start, int End, double Value)
{
    /// <summary>Window midpoint (for plotting).</summary>
    public double Center => (Start + End) / 2.0;
}

/// <summary>
/// Sliding-window composition profiles. Counts are updated incrementally as the window moves (one base in, one
/// out), so a profile costs O(N) regardless of window size.
/// </summary>
public static class SlidingWindow
{
    /// <summary>GC fraction (over unambiguous bases) per window.</summary>
    public static IReadOnlyList<WindowValue> GcContent(Sequence sequence, int window, int step = 1)
        => Profile(sequence, window, step, (a, c, g, t) => a + c + g + t == 0 ? 0 : (double)(g + c) / (a + c + g + t));

    /// <summary>GC skew (G − C)/(G + C) per window; its sign flips at replication origins and termini.</summary>
    public static IReadOnlyList<WindowValue> GcSkew(Sequence sequence, int window, int step = 1)
        => Profile(sequence, window, step, (a, c, g, t) => g + c == 0 ? 0 : (double)(g - c) / (g + c));

    /// <summary>AT skew (A − T)/(A + T) per window.</summary>
    public static IReadOnlyList<WindowValue> AtSkew(Sequence sequence, int window, int step = 1)
        => Profile(sequence, window, step, (a, c, g, t) => a + t == 0 ? 0 : (double)(a - t) / (a + t));

    /// <summary>
    /// Cumulative GC skew over non-overlapping windows. In bacteria the minimum marks the replication origin and the
    /// maximum the terminus.
    /// </summary>
    public static IReadOnlyList<WindowValue> CumulativeGcSkew(Sequence sequence, int window)
    {
        var skew = GcSkew(sequence, window, window);
        var result = new List<WindowValue>(skew.Count);
        double acc = 0;
        foreach (var w in skew) { acc += w.Value; result.Add(w with { Value = acc }); }
        return result;
    }

    /// <summary>Shannon entropy (bits, 0–2) of the base composition per window.</summary>
    public static IReadOnlyList<WindowValue> Entropy(Sequence sequence, int window, int step = 1)
        => Profile(sequence, window, step, (a, c, g, t) =>
        {
            double n = a + c + g + t, h = 0;
            if (n == 0) return 0;
            foreach (int x in new[] { a, c, g, t }) if (x > 0) { double p = x / n; h -= p * Math.Log2(p); }
            return h;
        });

    /// <summary>A custom per-window statistic computed from the A, C, G, T counts of each window.</summary>
    public static IReadOnlyList<WindowValue> Profile(Sequence sequence, int window, int step, Func<int, int, int, int, double> statistic)
    {
        if (window < 1) throw new ArgumentOutOfRangeException(nameof(window));
        if (step < 1) throw new ArgumentOutOfRangeException(nameof(step));
        var s = sequence.Span;
        var result = new List<WindowValue>(Math.Max(0, (s.Length - window) / step + 1));
        if (s.Length < window) return result;
        Span<int> counts = stackalloc int[5];
        counts.Clear();
        for (int i = 0; i < window; i++) counts[Slot(s[i])]++;
        int start = 0;
        while (true)
        {
            result.Add(new WindowValue(start, start + window, statistic(counts[0], counts[1], counts[2], counts[3])));
            int next = start + step;
            if (next + window > s.Length) break;
            if (step < window)
            {
                for (int i = start; i < next; i++) counts[Slot(s[i])]--;
                for (int i = start + window; i < next + window; i++) counts[Slot(s[i])]++;
            }
            else
            {
                counts.Clear();
                for (int i = next; i < next + window; i++) counts[Slot(s[i])]++;
            }
            start = next;
        }
        return result;
    }

    private static int Slot(byte b) => (b & 0xDF) switch { (byte)'A' => 0, (byte)'C' => 1, (byte)'G' => 2, (byte)'T' or (byte)'U' => 3, _ => 4 };
}

/// <summary>How low-complexity regions are masked.</summary>
public enum MaskStyle
{
    /// <summary>Lower-case the residues (soft-masking, reversible).</summary>
    Soft,

    /// <summary>Replace residues with <c>N</c> (hard-masking).</summary>
    Hard,
}

/// <summary>
/// Low-complexity detection with the DUST score (Tatusov &amp; Lipman; Morgulis et al. 2006): within a window of
/// w bases the triplet counts c_t give score = Σ c_t(c_t − 1)/2 / (l − 1), with l = w − 2 triplets. A random
/// sequence scores about 0.5–1.5; tandem repeats and homopolymers score far higher. The window score is maintained
/// incrementally in O(N).
/// </summary>
public static class Complexity
{
    /// <summary>DUST score of the whole sequence (one window).</summary>
    public static double DustScore(Sequence sequence) => DustScore(sequence.Span);

    /// <summary>DUST score of a residue span.</summary>
    public static double DustScore(ReadOnlySpan<byte> s)
    {
        Span<int> counts = stackalloc int[64];
        counts.Clear();
        long sum = 0; int l = 0;
        for (int i = 0; i + 3 <= s.Length; i++)
        {
            int t = Triplet(s, i);
            if (t < 0) continue;
            sum += counts[t]++;
            l++;
        }
        return l <= 1 ? 0 : (double)sum / (l - 1);
    }

    /// <summary>DUST score of every window of <paramref name="window"/> bases (step 1).</summary>
    public static IReadOnlyList<WindowValue> DustProfile(Sequence sequence, int window = 64)
    {
        var s = sequence.Span;
        var result = new List<WindowValue>();
        if (window < 4 || s.Length < window) return result;
        Span<int> counts = stackalloc int[64];
        counts.Clear();
        long sum = 0;
        int l = 0;
        int triplets = window - 2;
        // Prime with the first window's triplets.
        for (int i = 0; i < triplets; i++)
        {
            int t = Triplet(s, i);
            if (t >= 0) { sum += counts[t]++; l++; }
        }
        for (int start = 0; ; start++)
        {
            result.Add(new WindowValue(start, start + window, l <= 1 ? 0 : (double)sum / (l - 1)));
            int add = start + triplets; // next triplet entering
            if (add + 3 > s.Length) break;
            int outT = Triplet(s, start);
            if (outT >= 0) { sum -= --counts[outT]; l--; }
            int inT = Triplet(s, add);
            if (inT >= 0) { sum += counts[inT]++; l++; }
        }
        return result;
    }

    /// <summary>Intervals (0-based half-open, merged) whose DUST window score exceeds <paramref name="threshold"/> (dustmasker's default level is 20).</summary>
    public static IReadOnlyList<(int Start, int End)> LowComplexityRegions(Sequence sequence, int window = 64, double threshold = 20)
    {
        var regions = new List<(int Start, int End)>();
        IReadOnlyList<WindowValue> profile = sequence.Length >= window ? DustProfile(sequence, window)
            : [new WindowValue(0, sequence.Length, DustScore(sequence))];
        foreach (var w in profile)
        {
            if (w.Value <= threshold) continue;
            if (regions.Count > 0 && w.Start <= regions[^1].End) regions[^1] = (regions[^1].Start, Math.Max(regions[^1].End, w.End));
            else regions.Add((w.Start, w.End));
        }
        return regions;
    }

    /// <summary>A copy with low-complexity regions masked.</summary>
    public static Sequence Mask(Sequence sequence, int window = 64, double threshold = 20, MaskStyle style = MaskStyle.Soft)
    {
        var dst = sequence.Span.ToArray();
        foreach (var (start, end) in LowComplexityRegions(sequence, window, threshold))
            for (int i = start; i < end; i++)
                dst[i] = style == MaskStyle.Soft ? (byte)(dst[i] | 0x20) : (byte)'N';
        return Sequence.FromBuffer(sequence.Id, dst, sequence.Alphabet, sequence.Description, validate: false);
    }

    /// <summary>
    /// Linguistic complexity: the product over word sizes 1..<paramref name="maxWordLength"/> of observed distinct
    /// words / maximum possible distinct words (Trifonov 1990). 1 for maximally complex sequences.
    /// </summary>
    public static double LinguisticComplexity(Sequence sequence, int maxWordLength = 6)
    {
        double product = 1;
        int n = sequence.Length;
        for (int k = 1; k <= maxWordLength && k <= n; k++)
        {
            var counter = new KmerCounter(k).Add(sequence);
            double possible = Math.Min(Math.Pow(4, k), n - k + 1);
            product *= counter.Distinct / possible;
        }
        return product;
    }

    private static int Triplet(ReadOnlySpan<byte> s, int i)
    {
        int a = Kmer.Code(s[i]), b = Kmer.Code(s[i + 1]), c = Kmer.Code(s[i + 2]);
        return (a | b | c) < 0 ? -1 : (a << 4) | (b << 2) | c;
    }
}

/// <summary>Summary statistics for a set of sequences (reads, contigs, scaffolds).</summary>
public sealed record SequenceSetSummary(
    int Count, long TotalLength, int MinLength, int MaxLength, double MeanLength, double MedianLength,
    long N50, int L50, long N90, int L90, double GcContent, long AmbiguousBases)
{
    public override string ToString()
        => $"{Count:N0} sequences, {TotalLength:N0} residues (min {MinLength:N0}, max {MaxLength:N0}, mean {MeanLength:F1}), N50 {N50:N0}, L50 {L50}, GC {GcContent:P2}";
}

/// <summary>Length statistics: Nx/Lx, distributions and set summaries.</summary>
public static class SequenceStatistics
{
    /// <summary>
    /// Nx: the length L such that sequences of length ≥ L contain at least x% of all residues (N50 for x = 50).
    /// Returns (Nx, Lx) where Lx is the number of sequences needed.
    /// </summary>
    public static (long Nx, int Lx) Nx(IEnumerable<long> lengths, double x = 50)
    {
        if (x is <= 0 or > 100) throw new ArgumentOutOfRangeException(nameof(x));
        var sorted = lengths.Where(l => l > 0).ToArray();
        Array.Sort(sorted);
        Array.Reverse(sorted);
        long total = 0;
        foreach (var l in sorted) total += l;
        if (total == 0) return (0, 0);
        double target = total * x / 100.0;
        long acc = 0;
        for (int i = 0; i < sorted.Length; i++)
        {
            acc += sorted[i];
            if (acc >= target) return (sorted[i], i + 1);
        }
        return (sorted[^1], sorted.Length);
    }

    /// <summary>N50 of a set of lengths.</summary>
    public static long N50(IEnumerable<long> lengths) => Nx(lengths, 50).Nx;

    /// <summary>Summarises a set of sequences in one pass over their compositions.</summary>
    public static SequenceSetSummary Summarize(IEnumerable<Sequence> sequences)
    {
        var lengths = new List<long>();
        long gc = 0, acgt = 0, ambiguous = 0;
        foreach (var s in sequences)
        {
            lengths.Add(s.Length);
            if (s.Alphabet.IsNucleotide)
            {
                var c = s.NucleotideCounts();
                gc += c.G + c.C; acgt += c.Acgt; ambiguous += c.N + c.Other;
            }
        }
        if (lengths.Count == 0) return new SequenceSetSummary(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        var sorted = lengths.ToArray();
        Array.Sort(sorted);
        long total = 0;
        foreach (var l in sorted) total += l;
        double median = sorted.Length % 2 == 1 ? sorted[sorted.Length / 2] : (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]) / 2.0;
        var (n50, l50) = Nx(sorted, 50);
        var (n90, l90) = Nx(sorted, 90);
        return new SequenceSetSummary(sorted.Length, total, (int)sorted[0], (int)sorted[^1], (double)total / sorted.Length, median,
            n50, l50, n90, l90, acgt == 0 ? 0 : (double)gc / acgt, ambiguous);
    }

    /// <summary>Length histogram with bins of <paramref name="binWidth"/>: (bin start, count), ascending.</summary>
    public static IReadOnlyList<(long BinStart, int Count)> LengthHistogram(IEnumerable<long> lengths, long binWidth)
    {
        if (binWidth < 1) throw new ArgumentOutOfRangeException(nameof(binWidth));
        var bins = new SortedDictionary<long, int>();
        foreach (var l in lengths)
        {
            long b = l / binWidth * binWidth;
            bins[b] = bins.GetValueOrDefault(b) + 1;
        }
        return bins.Select(kv => (kv.Key, kv.Value)).ToList();
    }
}

/// <summary>Seeded random sequences and composition-preserving shuffles (null models for motif and ORF statistics).</summary>
public static class RandomSequence
{
    /// <summary>A random DNA sequence with the given expected GC fraction.</summary>
    public static Sequence Dna(int length, int seed, double gcContent = 0.5, string id = "random")
    {
        if (gcContent is < 0 or > 1) throw new ArgumentOutOfRangeException(nameof(gcContent));
        var rng = new Random(seed);
        var dst = new byte[length];
        for (int i = 0; i < length; i++)
        {
            bool gc = rng.NextDouble() < gcContent;
            bool second = rng.Next(2) == 1;
            dst[i] = gc ? (second ? (byte)'G' : (byte)'C') : (second ? (byte)'T' : (byte)'A');
        }
        return Sequence.FromBuffer(id, dst, Alphabet.Dna, validate: false);
    }

    /// <summary>A random sequence drawing each residue from <paramref name="weights"/> (symbol → relative weight).</summary>
    public static Sequence FromComposition(int length, int seed, IReadOnlyDictionary<char, double> weights, string id = "random")
    {
        var symbols = weights.Keys.ToArray();
        var cumulative = new double[symbols.Length];
        double acc = 0;
        for (int i = 0; i < symbols.Length; i++) { acc += weights[symbols[i]]; cumulative[i] = acc; }
        if (acc <= 0) throw new ArgumentException("Weights must sum to a positive value.", nameof(weights));
        var rng = new Random(seed);
        var dst = new byte[length];
        for (int i = 0; i < length; i++)
        {
            double u = rng.NextDouble() * acc;
            int j = Array.BinarySearch(cumulative, u);
            if (j < 0) j = ~j;
            dst[i] = (byte)symbols[Math.Min(j, symbols.Length - 1)];
        }
        return Sequence.FromBuffer(id, dst, Alphabet.Detect(dst), validate: false);
    }

    /// <summary>A random protein with uniform amino-acid usage over the 20 standard residues.</summary>
    public static Sequence Protein(int length, int seed, string id = "random")
    {
        const string aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
        var rng = new Random(seed);
        var dst = new byte[length];
        for (int i = 0; i < length; i++) dst[i] = (byte)aminoAcids[rng.Next(aminoAcids.Length)];
        return Sequence.FromBuffer(id, dst, Alphabet.ProteinExtended, validate: false);
    }

    /// <summary>A Fisher–Yates shuffle: exactly the same residue composition in random order.</summary>
    public static Sequence Shuffle(Sequence sequence, int seed)
    {
        var dst = sequence.Span.ToArray();
        new Random(seed).Shuffle(dst);
        return Sequence.FromBuffer(sequence.Id, dst, sequence.Alphabet, sequence.Description, validate: false);
    }

    /// <summary>
    /// Shuffles non-overlapping blocks of <paramref name="blockLength"/> residues (e.g. 3 to keep codons intact),
    /// preserving the composition of the blocks; a trailing partial block stays at the end.
    /// </summary>
    public static Sequence ShuffleBlocks(Sequence sequence, int blockLength, int seed)
    {
        if (blockLength < 1) throw new ArgumentOutOfRangeException(nameof(blockLength));
        int blocks = sequence.Length / blockLength;
        var order = Enumerable.Range(0, blocks).ToArray();
        new Random(seed).Shuffle(order);
        var src = sequence.Span;
        var dst = new byte[sequence.Length];
        for (int i = 0; i < blocks; i++) src.Slice(order[i] * blockLength, blockLength).CopyTo(dst.AsSpan(i * blockLength));
        src[(blocks * blockLength)..].CopyTo(dst.AsSpan(blocks * blockLength));
        return Sequence.FromBuffer(sequence.Id, dst, sequence.Alphabet, sequence.Description, validate: false);
    }
}
