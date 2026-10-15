using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Sequences;

/// <summary>
/// Dense k-mer counts (1 ≤ k ≤ 12) indexed by 2-bit code (A=0, C=1, G=2, T=3, first base most significant), filled
/// by <see cref="BioKernels.CountKmers"/>. Windows containing a non-ACGT symbol are skipped.
/// </summary>
public sealed class KmerProfile
{
    private readonly int[] _counts;

    private KmerProfile(int k, int[] counts, long total)
    {
        K = k;
        _counts = counts;
        Total = total;
    }

    /// <summary>Counts the k-mers of one sequence.</summary>
    public static KmerProfile Count(Sequence sequence, int k) => Count([sequence], k);

    /// <summary>Counts the k-mers of many sequences (reads, contigs) into one profile.</summary>
    /// <param name="canonical">Count each k-mer together with its reverse complement (strand-independent spectra).</param>
    public static KmerProfile Count(IEnumerable<Sequence> sequences, int k, bool canonical = false)
    {
        if (k is < 1 or > 12) throw new ArgumentOutOfRangeException(nameof(k), "Dense k-mer profiles support 1 ≤ k ≤ 12; use KmerCounter for larger k.");
        var counts = new int[1 << (2 * k)];
        long total = 0;
        foreach (var s in sequences) total += BioKernels.CountKmers(s.Span, k, counts);
        if (canonical)
        {
            // A k-mer on one strand is its reverse complement on the other: fold each pair onto the smaller code.
            for (int code = 0; code < counts.Length; code++)
            {
                int rc = (int)Kmer.ReverseComplement((ulong)code, k);
                if (rc > code) { counts[code] += counts[rc]; counts[rc] = 0; }
            }
        }
        return new KmerProfile(k, counts, total);
    }

    /// <summary>k.</summary>
    public int K { get; }

    /// <summary>Number of k-mers counted.</summary>
    public long Total { get; }

    /// <summary>Counts indexed by k-mer code (4^k entries).</summary>
    public ReadOnlySpan<int> Counts => _counts;

    /// <summary>Count of one k-mer (any case, DNA or RNA).</summary>
    public int this[string kmer] => kmer.Length == K && Kmer.TryEncode(kmer, out ulong code) ? _counts[code] : 0;

    /// <summary>Number of distinct k-mers seen.</summary>
    public int Distinct
    {
        get { int d = 0; foreach (int c in _counts) if (c > 0) d++; return d; }
    }

    /// <summary>Relative frequencies (sum 1).</summary>
    public double[] Frequencies()
    {
        var f = new double[_counts.Length];
        if (Total == 0) return f;
        long sum = 0;
        foreach (int c in _counts) sum += c;
        for (int i = 0; i < f.Length; i++) f[i] = (double)_counts[i] / sum;
        return f;
    }

    /// <summary>
    /// The k-mer spectrum: <c>spectrum[m]</c> is the number of distinct k-mers seen exactly m times. Used to estimate
    /// coverage and genome size and to pick an error threshold for assembly.
    /// </summary>
    public long[] Spectrum()
    {
        int max = 0;
        foreach (int c in _counts) if (c > max) max = c;
        var h = new long[max + 1];
        foreach (int c in _counts) if (c > 0) h[c]++;
        return h;
    }

    /// <summary>The <paramref name="count"/> most frequent k-mers.</summary>
    public IReadOnlyList<(string Kmer, int Count)> Top(int count)
    {
        var pq = new PriorityQueue<int, int>();
        for (int i = 0; i < _counts.Length; i++)
        {
            if (_counts[i] == 0) continue;
            if (pq.Count < count) pq.Enqueue(i, _counts[i]);
            else if (pq.TryPeek(out _, out int min) && _counts[i] > min) { pq.Dequeue(); pq.Enqueue(i, _counts[i]); }
        }
        var list = new List<(string, int)>(pq.Count);
        while (pq.TryDequeue(out int code, out int c)) list.Add((Kmer.Decode((ulong)code, K), c));
        list.Reverse();
        return list;
    }

    /// <summary>Euclidean distance between the frequency vectors of two profiles with the same k.</summary>
    public double EuclideanDistance(KmerProfile other)
    {
        RequireSameK(other);
        var a = Frequencies(); var b = other.Frequencies();
        double s = 0;
        for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
        return Math.Sqrt(s);
    }

    /// <summary>Cosine similarity of the count vectors (1 = identical composition).</summary>
    public double CosineSimilarity(KmerProfile other)
    {
        RequireSameK(other);
        double dot = 0, na = 0, nb = 0;
        for (int i = 0; i < _counts.Length; i++)
        {
            double x = _counts[i], y = other._counts[i];
            dot += x * y; na += x * x; nb += y * y;
        }
        return na == 0 || nb == 0 ? 0 : dot / Math.Sqrt(na * nb);
    }

    private void RequireSameK(KmerProfile other)
    {
        if (other.K != K) throw new ArgumentException("Profiles must use the same k.");
    }
}

/// <summary>2-bit k-mer encoding helpers (k ≤ 32 fits in a <see cref="ulong"/>).</summary>
public static class Kmer
{
    private static readonly sbyte[] Codes = BuildCodes();

    private static sbyte[] BuildCodes()
    {
        var c = new sbyte[256];
        Array.Fill(c, (sbyte)-1);
        c['A'] = c['a'] = 0; c['C'] = c['c'] = 1; c['G'] = c['g'] = 2; c['T'] = c['t'] = c['U'] = c['u'] = 3;
        return c;
    }

    /// <summary>2-bit code of a base, or −1 for anything but A/C/G/T/U.</summary>
    public static int Code(byte b) => Codes[b];

    /// <summary>Encodes a k-mer; false when it contains a non-ACGT symbol or is longer than 32.</summary>
    public static bool TryEncode(ReadOnlySpan<char> kmer, out ulong code)
    {
        code = 0;
        if (kmer.Length > 32) return false;
        foreach (char ch in kmer)
        {
            int c = ch < 256 ? Codes[ch] : -1;
            if (c < 0) return false;
            code = (code << 2) | (uint)c;
        }
        return true;
    }

    /// <summary>Encodes a k-mer given as bytes.</summary>
    public static bool TryEncode(ReadOnlySpan<byte> kmer, out ulong code)
    {
        code = 0;
        if (kmer.Length > 32) return false;
        foreach (byte b in kmer)
        {
            int c = Codes[b];
            if (c < 0) return false;
            code = (code << 2) | (uint)c;
        }
        return true;
    }

    /// <summary>Decodes a k-mer code.</summary>
    public static string Decode(ulong code, int k)
        => string.Create(k, code, (span, c) => { for (int i = k - 1; i >= 0; i--) { span[i] = "ACGT"[(int)(c & 3)]; c >>= 2; } });

    /// <summary>Reverse complement of a k-mer code.</summary>
    public static ulong ReverseComplement(ulong code, int k)
    {
        ulong rc = 0;
        for (int i = 0; i < k; i++) { rc = (rc << 2) | (3 - (code & 3)); code >>= 2; }
        return rc;
    }

    /// <summary>The smaller of a code and its reverse complement.</summary>
    public static ulong Canonical(ulong code, int k) => Math.Min(code, ReverseComplement(code, k));

    /// <summary>
    /// Enumerates (position, code) for every valid k-mer of <paramref name="sequence"/> with a rolling 2-bit encoding;
    /// windows containing a non-ACGT symbol are skipped.
    /// </summary>
    public static void ForEach(ReadOnlySpan<byte> sequence, int k, Action<int, ulong> action)
    {
        if (k is < 1 or > 32) throw new ArgumentOutOfRangeException(nameof(k));
        ulong mask = k == 32 ? ulong.MaxValue : (1UL << (2 * k)) - 1;
        ulong code = 0;
        int valid = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            int c = Codes[sequence[i]];
            if (c < 0) { valid = 0; code = 0; continue; }
            code = ((code << 2) | (uint)c) & mask;
            if (++valid >= k) action(i - k + 1, code);
        }
    }
}

/// <summary>Sparse k-mer counter for any k up to 32 (hash table keyed by the 2-bit code).</summary>
public sealed class KmerCounter
{
    private readonly Dictionary<ulong, int> _counts = [];

    /// <summary>Creates a counter for k-mers of length <paramref name="k"/>.</summary>
    public KmerCounter(int k, bool canonical = false)
    {
        if (k is < 1 or > 32) throw new ArgumentOutOfRangeException(nameof(k));
        K = k;
        Canonical = canonical;
    }

    public int K { get; }

    /// <summary>Whether k-mers are merged with their reverse complements.</summary>
    public bool Canonical { get; }

    /// <summary>Number of distinct k-mers.</summary>
    public int Distinct => _counts.Count;

    /// <summary>Adds every k-mer of <paramref name="sequence"/>.</summary>
    public KmerCounter Add(Sequence sequence) => Add(sequence.Span);

    /// <summary>Adds every k-mer of a residue span.</summary>
    public KmerCounter Add(ReadOnlySpan<byte> sequence)
    {
        int k = K;
        ulong mask = k == 32 ? ulong.MaxValue : (1UL << (2 * k)) - 1;
        ulong code = 0;
        int valid = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            int c = Kmer.Code(sequence[i]);
            if (c < 0) { valid = 0; code = 0; continue; }
            code = ((code << 2) | (uint)c) & mask;
            if (++valid < k) continue;
            ulong key = Canonical ? Kmer.Canonical(code, k) : code;
            System.Runtime.InteropServices.CollectionsMarshal.GetValueRefOrAddDefault(_counts, key, out _)++;
        }
        return this;
    }

    /// <summary>Count of one k-mer.</summary>
    public int this[string kmer]
    {
        get
        {
            if (kmer.Length != K || !Kmer.TryEncode(kmer, out ulong code)) return 0;
            if (Canonical) code = Kmer.Canonical(code, K);
            return _counts.GetValueOrDefault(code);
        }
    }

    /// <summary>All (k-mer, count) pairs.</summary>
    public IEnumerable<(string Kmer, int Count)> Entries() => _counts.Select(kv => (Kmer.Decode(kv.Key, K), kv.Value));

    /// <summary>Raw (code, count) pairs.</summary>
    public IReadOnlyDictionary<ulong, int> Codes => _counts;

    /// <summary><c>spectrum[m]</c> = number of distinct k-mers seen exactly m times.</summary>
    public long[] Spectrum()
    {
        int max = 0;
        foreach (int c in _counts.Values) if (c > max) max = c;
        var h = new long[max + 1];
        foreach (int c in _counts.Values) h[c]++;
        return h;
    }
}
