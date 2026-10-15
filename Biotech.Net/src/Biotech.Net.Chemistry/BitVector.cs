using System.Numerics;
using System.Text;

namespace Biotech.Net.Chemistry;

/// <summary>
/// A fixed-length bit vector backed by 64-bit words, used for molecular fingerprints. Population counts and the
/// similarity coefficients in <see cref="Similarity"/> use hardware popcount.
/// </summary>
public sealed class BitVector : IEquatable<BitVector>
{
    private readonly ulong[] _words;

    /// <summary>Creates an all-zero vector of <paramref name="length"/> bits.</summary>
    public BitVector(int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        Length = length;
        _words = new ulong[(length + 63) >> 6];
    }

    /// <summary>Creates a vector of <paramref name="length"/> bits with the given bits set.</summary>
    public BitVector(int length, IEnumerable<int> onBits) : this(length)
    {
        foreach (int b in onBits) Set(b);
    }

    private BitVector(int length, ulong[] words)
    {
        Length = length;
        _words = words;
    }

    /// <summary>Number of bits.</summary>
    public int Length { get; }

    /// <summary>The backing words (bit <c>i</c> is bit <c>i % 64</c> of word <c>i / 64</c>).</summary>
    public ReadOnlySpan<ulong> Words => _words;

    /// <summary>Gets or sets bit <paramref name="index"/>.</summary>
    public bool this[int index]
    {
        get
        {
            if ((uint)index >= (uint)Length) throw new ArgumentOutOfRangeException(nameof(index));
            return (_words[index >> 6] & (1UL << (index & 63))) != 0;
        }
        set
        {
            if (value) Set(index); else Clear(index);
        }
    }

    /// <summary>Sets bit <paramref name="index"/>.</summary>
    public void Set(int index)
    {
        if ((uint)index >= (uint)Length) throw new ArgumentOutOfRangeException(nameof(index));
        _words[index >> 6] |= 1UL << (index & 63);
    }

    /// <summary>Clears bit <paramref name="index"/>.</summary>
    public void Clear(int index)
    {
        if ((uint)index >= (uint)Length) throw new ArgumentOutOfRangeException(nameof(index));
        _words[index >> 6] &= ~(1UL << (index & 63));
    }

    /// <summary>Number of set bits.</summary>
    public int Count
    {
        get
        {
            int n = 0;
            foreach (ulong w in _words) n += BitOperations.PopCount(w);
            return n;
        }
    }

    /// <summary>Indices of the set bits in ascending order.</summary>
    public IEnumerable<int> OnBits()
    {
        for (int w = 0; w < _words.Length; w++)
        {
            ulong word = _words[w];
            while (word != 0)
            {
                int bit = BitOperations.TrailingZeroCount(word);
                yield return (w << 6) + bit;
                word &= word - 1;
            }
        }
    }

    /// <summary>A copy of this vector.</summary>
    public BitVector Clone() => new(Length, (ulong[])_words.Clone());

    /// <summary>Bitwise AND (vectors must have equal length).</summary>
    public static BitVector operator &(BitVector a, BitVector b) => Combine(a, b, static (x, y) => x & y);

    /// <summary>Bitwise OR (vectors must have equal length).</summary>
    public static BitVector operator |(BitVector a, BitVector b) => Combine(a, b, static (x, y) => x | y);

    /// <summary>Bitwise XOR (vectors must have equal length).</summary>
    public static BitVector operator ^(BitVector a, BitVector b) => Combine(a, b, static (x, y) => x ^ y);

    private static BitVector Combine(BitVector a, BitVector b, Func<ulong, ulong, ulong> op)
    {
        CheckLength(a, b);
        var w = new ulong[a._words.Length];
        for (int i = 0; i < w.Length; i++) w[i] = op(a._words[i], b._words[i]);
        return new BitVector(a.Length, w);
    }

    internal static void CheckLength(BitVector a, BitVector b)
    {
        if (a.Length != b.Length) throw new ArgumentException($"Bit vectors differ in length ({a.Length} vs {b.Length}).");
    }

    /// <summary>Number of bits set in both vectors.</summary>
    public static int IntersectionCount(BitVector a, BitVector b)
    {
        CheckLength(a, b);
        int n = 0;
        var x = a._words;
        var y = b._words;
        for (int i = 0; i < x.Length; i++) n += BitOperations.PopCount(x[i] & y[i]);
        return n;
    }

    /// <summary>
    /// Folds the vector to <paramref name="newLength"/> bits by OR-ing bit <c>i</c> into bit <c>i % newLength</c>
    /// (RDKit <c>FoldFingerprint</c> semantics when <paramref name="newLength"/> divides <see cref="Length"/>).
    /// </summary>
    public BitVector Fold(int newLength)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(newLength);
        var r = new BitVector(newLength);
        foreach (int b in OnBits()) r.Set(b % newLength);
        return r;
    }

    /// <summary>The vector as a string of '0'/'1' characters (bit 0 first).</summary>
    public string ToBitString()
    {
        var sb = new StringBuilder(Length);
        for (int i = 0; i < Length; i++) sb.Append((_words[i >> 6] & (1UL << (i & 63))) != 0 ? '1' : '0');
        return sb.ToString();
    }

    /// <summary>Parses a string of '0'/'1' characters (bit 0 first).</summary>
    public static BitVector FromBitString(string bits)
    {
        var v = new BitVector(bits.Length);
        for (int i = 0; i < bits.Length; i++)
        {
            if (bits[i] == '1') v.Set(i);
            else if (bits[i] != '0') throw new FormatException($"Invalid bit character '{bits[i]}' at {i}.");
        }
        return v;
    }

    /// <summary>Compact Base64 encoding of the words (little-endian); see <see cref="FromBase64"/>.</summary>
    public string ToBase64()
    {
        var bytes = new byte[_words.Length * 8];
        Buffer.BlockCopy(_words, 0, bytes, 0, bytes.Length);
        return $"{Length}:{Convert.ToBase64String(bytes)}";
    }

    /// <summary>Decodes a vector written by <see cref="ToBase64"/>.</summary>
    public static BitVector FromBase64(string text)
    {
        int colon = text.IndexOf(':');
        if (colon < 0) throw new FormatException("Expected '<length>:<base64>'.");
        int length = int.Parse(text.AsSpan(0, colon), System.Globalization.CultureInfo.InvariantCulture);
        var bytes = Convert.FromBase64String(text[(colon + 1)..]);
        var v = new BitVector(length);
        if (bytes.Length != v._words.Length * 8) throw new FormatException("Length does not match the encoded data.");
        Buffer.BlockCopy(bytes, 0, v._words, 0, bytes.Length);
        return v;
    }

    /// <inheritdoc/>
    public bool Equals(BitVector? other) => other is not null && other.Length == Length && _words.AsSpan().SequenceEqual(other._words);

    /// <inheritdoc/>
    public override bool Equals(object? obj) => obj is BitVector b && Equals(b);

    /// <inheritdoc/>
    public override int GetHashCode()
    {
        var h = new HashCode();
        h.Add(Length);
        foreach (ulong w in _words) h.Add(w);
        return h.ToHashCode();
    }

    /// <inheritdoc/>
    public override string ToString() => $"BitVector({Length} bits, {Count} set)";
}

/// <summary>
/// A sparse count vector (feature id → count), e.g. unfolded Morgan or atom-pair fingerprints with counts.
/// </summary>
public sealed class CountVector
{
    private readonly Dictionary<long, int> _counts = [];

    /// <summary>Creates an empty count vector over a feature space of <paramref name="length"/> ids (informational).</summary>
    public CountVector(long length = uint.MaxValue + 1L) => Length = length;

    /// <summary>Size of the feature space.</summary>
    public long Length { get; }

    /// <summary>Feature counts.</summary>
    public IReadOnlyDictionary<long, int> Counts => _counts;

    /// <summary>Number of distinct non-zero features.</summary>
    public int NonZeroCount => _counts.Count;

    /// <summary>Sum of all counts.</summary>
    public int TotalCount => _counts.Values.Sum();

    /// <summary>Count of feature <paramref name="id"/> (0 when absent).</summary>
    public int this[long id] => _counts.GetValueOrDefault(id);

    /// <summary>Adds <paramref name="count"/> to feature <paramref name="id"/>.</summary>
    public void Add(long id, int count = 1)
    {
        int c = _counts.GetValueOrDefault(id) + count;
        if (c == 0) _counts.Remove(id); else _counts[id] = c;
    }

    /// <summary>Folds into a bit vector: feature <c>id</c> sets bit <c>id mod nBits</c>.</summary>
    public BitVector ToBitVector(int nBits)
    {
        var v = new BitVector(nBits);
        foreach (var id in _counts.Keys) v.Set((int)(((id % nBits) + nBits) % nBits));
        return v;
    }

    /// <inheritdoc/>
    public override string ToString() => $"CountVector({NonZeroCount} features, total {TotalCount})";
}

/// <summary>Similarity coefficients for binary fingerprints.</summary>
public enum SimilarityMetric
{
    /// <summary>Tanimoto / Jaccard: c / (a + b − c).</summary>
    Tanimoto,
    /// <summary>Dice / Sørensen: 2c / (a + b).</summary>
    Dice,
    /// <summary>Cosine / Ochiai: c / √(a·b).</summary>
    Cosine,
    /// <summary>Sokal–Sneath: c / (2a + 2b − 3c).</summary>
    Sokal,
    /// <summary>Kulczynski: c(a + b) / (2ab).</summary>
    Kulczynski,
    /// <summary>Russell–Rao: c / n.</summary>
    Russel,
}

/// <summary>
/// Fingerprint similarity coefficients (definitions as in RDKit's <c>DataStructs</c>). For two vectors with
/// <c>a</c> and <c>b</c> bits set and <c>c</c> bits in common. Two all-zero vectors have similarity 0 (RDKit returns 0 for Tanimoto too).
/// </summary>
public static class Similarity
{
    /// <summary>Tanimoto coefficient c / (a + b − c).</summary>
    public static double Tanimoto(BitVector x, BitVector y)
    {
        var (a, b, c) = Counts(x, y);
        int d = a + b - c;
        return d == 0 ? 0 : (double)c / d;
    }

    /// <summary>Dice coefficient 2c / (a + b).</summary>
    public static double Dice(BitVector x, BitVector y)
    {
        var (a, b, c) = Counts(x, y);
        return a + b == 0 ? 0 : 2.0 * c / (a + b);
    }

    /// <summary>Cosine (Ochiai) coefficient c / √(ab).</summary>
    public static double Cosine(BitVector x, BitVector y)
    {
        var (a, b, c) = Counts(x, y);
        return a == 0 || b == 0 ? 0 : c / Math.Sqrt((double)a * b);
    }

    /// <summary>
    /// Tversky index c / (α(a − c) + β(b − c) + c). α = β = 1 gives Tanimoto, α = β = 0.5 gives Dice; α = 1, β = 0
    /// measures how much of the query <paramref name="x"/> is contained in <paramref name="y"/>.
    /// </summary>
    public static double Tversky(BitVector x, BitVector y, double alpha, double beta)
    {
        var (a, b, c) = Counts(x, y);
        double d = alpha * (a - c) + beta * (b - c) + c;
        return d <= 0 ? 0 : c / d;
    }

    /// <summary>Sokal–Sneath coefficient c / (2a + 2b − 3c).</summary>
    public static double Sokal(BitVector x, BitVector y)
    {
        var (a, b, c) = Counts(x, y);
        int d = 2 * a + 2 * b - 3 * c;
        return d == 0 ? 0 : (double)c / d;
    }

    /// <summary>Kulczynski coefficient c(a + b) / (2ab).</summary>
    public static double Kulczynski(BitVector x, BitVector y)
    {
        var (a, b, c) = Counts(x, y);
        return a == 0 || b == 0 ? 0 : c * (double)(a + b) / (2.0 * a * b);
    }

    /// <summary>Russell–Rao coefficient c / n.</summary>
    public static double Russel(BitVector x, BitVector y)
    {
        var (_, _, c) = Counts(x, y);
        return x.Length == 0 ? 0 : (double)c / x.Length;
    }

    /// <summary>Computes the coefficient selected by <paramref name="metric"/>.</summary>
    public static double Compute(SimilarityMetric metric, BitVector x, BitVector y) => metric switch
    {
        SimilarityMetric.Tanimoto => Tanimoto(x, y),
        SimilarityMetric.Dice => Dice(x, y),
        SimilarityMetric.Cosine => Cosine(x, y),
        SimilarityMetric.Sokal => Sokal(x, y),
        SimilarityMetric.Kulczynski => Kulczynski(x, y),
        SimilarityMetric.Russel => Russel(x, y),
        _ => throw new ArgumentOutOfRangeException(nameof(metric)),
    };

    /// <summary>Tanimoto similarity of count vectors: Σmin / Σmax (RDKit <c>TanimotoSimilarity</c> on sparse int vectors).</summary>
    public static double Tanimoto(CountVector x, CountVector y)
    {
        long mins = 0, maxs = 0;
        foreach (var (id, cx) in x.Counts)
        {
            int cy = y[id];
            mins += Math.Min(cx, cy);
            maxs += Math.Max(cx, cy);
        }
        foreach (var (id, cy) in y.Counts) if (x[id] == 0) maxs += cy;
        return maxs == 0 ? 0 : (double)mins / maxs;
    }

    /// <summary>Dice similarity of count vectors: 2Σmin / (Σx + Σy) (RDKit <c>DiceSimilarity</c>).</summary>
    public static double Dice(CountVector x, CountVector y)
    {
        long mins = 0;
        foreach (var (id, cx) in x.Counts) mins += Math.Min(cx, y[id]);
        long total = (long)x.TotalCount + y.TotalCount;
        return total == 0 ? 0 : 2.0 * mins / total;
    }

    /// <summary>Similarity of <paramref name="query"/> to every target (parallel for large lists).</summary>
    public static double[] Bulk(BitVector query, IReadOnlyList<BitVector> targets, SimilarityMetric metric = SimilarityMetric.Tanimoto)
    {
        var r = new double[targets.Count];
        if (targets.Count < 2048)
        {
            for (int i = 0; i < r.Length; i++) r[i] = Compute(metric, query, targets[i]);
        }
        else
        {
            Parallel.For(0, r.Length, i => r[i] = Compute(metric, query, targets[i]));
        }
        return r;
    }

    /// <summary>Bulk Tanimoto similarity (RDKit <c>BulkTanimotoSimilarity</c>).</summary>
    public static double[] BulkTanimoto(BitVector query, IReadOnlyList<BitVector> targets) => Bulk(query, targets, SimilarityMetric.Tanimoto);

    /// <summary>Bulk Dice similarity.</summary>
    public static double[] BulkDice(BitVector query, IReadOnlyList<BitVector> targets) => Bulk(query, targets, SimilarityMetric.Dice);

    /// <summary>The full symmetric similarity matrix of a set of fingerprints.</summary>
    public static double[,] Matrix(IReadOnlyList<BitVector> fps, SimilarityMetric metric = SimilarityMetric.Tanimoto)
    {
        int n = fps.Count;
        var m = new double[n, n];
        Parallel.For(0, n, i =>
        {
            m[i, i] = fps[i].Count == 0 ? 0 : 1;
            for (int j = i + 1; j < n; j++)
            {
                double s = Compute(metric, fps[i], fps[j]);
                m[i, j] = s;
                m[j, i] = s;
            }
        });
        return m;
    }

    private static (int A, int B, int C) Counts(BitVector x, BitVector y)
    {
        BitVector.CheckLength(x, y);
        var wx = x.Words;
        var wy = y.Words;
        int a = 0, b = 0, c = 0;
        for (int i = 0; i < wx.Length; i++)
        {
            a += BitOperations.PopCount(wx[i]);
            b += BitOperations.PopCount(wy[i]);
            c += BitOperations.PopCount(wx[i] & wy[i]);
        }
        return (a, b, c);
    }
}

/// <summary>One hit of a similarity search.</summary>
/// <param name="Index">Position of the hit in the searched collection.</param>
/// <param name="Id">Identifier of the hit (molecule name or supplied id).</param>
/// <param name="Score">Similarity to the query.</param>
public sealed record SimilarityHit(int Index, string Id, double Score);

/// <summary>
/// An in-memory fingerprint database for fast similarity searching (threshold and top-k), with optional
/// popcount-bound pruning for Tanimoto (Swamidass &amp; Baldi bound: T ≤ min(a, b) / max(a, b)).
/// </summary>
public sealed class FingerprintIndex
{
    private readonly List<BitVector> _fps = [];
    private readonly List<string> _ids = [];
    private readonly List<int> _counts = [];

    /// <summary>Creates an empty index; <paramref name="fingerprinter"/> is used by the molecule overloads.</summary>
    public FingerprintIndex(Func<Molecule, BitVector>? fingerprinter = null)
        => Fingerprinter = fingerprinter ?? (m => MorganFingerprint.Compute(m));

    /// <summary>The fingerprint function used for molecules (default ECFP4, 2048 bits).</summary>
    public Func<Molecule, BitVector> Fingerprinter { get; }

    /// <summary>Number of entries.</summary>
    public int Count => _fps.Count;

    /// <summary>Fingerprint of entry <paramref name="index"/>.</summary>
    public BitVector this[int index] => _fps[index];

    /// <summary>Identifier of entry <paramref name="index"/>.</summary>
    public string IdAt(int index) => _ids[index];

    /// <summary>Adds a precomputed fingerprint.</summary>
    public void Add(string id, BitVector fp)
    {
        _fps.Add(fp);
        _ids.Add(id);
        _counts.Add(fp.Count);
    }

    /// <summary>Fingerprints and adds a molecule (its <see cref="Molecule.Name"/> is the id unless given).</summary>
    public void Add(Molecule mol, string? id = null) => Add(id ?? mol.Name, Fingerprinter(mol));

    /// <summary>Adds many molecules, fingerprinting in parallel.</summary>
    public void AddRange(IEnumerable<Molecule> molecules)
    {
        var list = molecules.ToList();
        var fps = new BitVector[list.Count];
        Parallel.For(0, list.Count, i => fps[i] = Fingerprinter(list[i]));
        for (int i = 0; i < list.Count; i++) Add(list[i].Name, fps[i]);
    }

    /// <summary>All entries with similarity ≥ <paramref name="threshold"/>, best first.</summary>
    public IReadOnlyList<SimilarityHit> Search(BitVector query, double threshold = 0.0, SimilarityMetric metric = SimilarityMetric.Tanimoto)
    {
        int qa = query.Count;
        var hits = new System.Collections.Concurrent.ConcurrentBag<SimilarityHit>();
        Parallel.For(0, _fps.Count, i =>
        {
            if (metric == SimilarityMetric.Tanimoto && threshold > 0)
            {
                int cb = _counts[i];
                double bound = Math.Max(qa, cb) == 0 ? 0 : (double)Math.Min(qa, cb) / Math.Max(qa, cb);
                if (bound < threshold) return;
            }
            double s = Similarity.Compute(metric, query, _fps[i]);
            if (s >= threshold) hits.Add(new SimilarityHit(i, _ids[i], s));
        });
        return hits.OrderByDescending(h => h.Score).ThenBy(h => h.Index).ToList();
    }

    /// <summary>All entries with similarity ≥ <paramref name="threshold"/> to a query molecule, best first.</summary>
    public IReadOnlyList<SimilarityHit> Search(Molecule query, double threshold = 0.0, SimilarityMetric metric = SimilarityMetric.Tanimoto)
        => Search(Fingerprinter(query), threshold, metric);

    /// <summary>The <paramref name="k"/> most similar entries, best first.</summary>
    public IReadOnlyList<SimilarityHit> TopK(BitVector query, int k, SimilarityMetric metric = SimilarityMetric.Tanimoto)
    {
        var scores = Similarity.Bulk(query, _fps, metric);
        return Enumerable.Range(0, scores.Length)
            .OrderByDescending(i => scores[i]).ThenBy(i => i)
            .Take(k)
            .Select(i => new SimilarityHit(i, _ids[i], scores[i]))
            .ToList();
    }

    /// <summary>The <paramref name="k"/> entries most similar to a query molecule, best first.</summary>
    public IReadOnlyList<SimilarityHit> TopK(Molecule query, int k, SimilarityMetric metric = SimilarityMetric.Tanimoto)
        => TopK(Fingerprinter(query), k, metric);
}
