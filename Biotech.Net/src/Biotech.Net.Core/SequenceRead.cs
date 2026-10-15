namespace Biotech.Net.Core;

/// <summary>
/// A sequencing read: a nucleotide <see cref="Sequence"/> with one Phred quality score per base (FASTQ).
/// Qualities are stored as raw Phred values (0–93), not ASCII.
/// </summary>
public sealed class SequenceRead
{
    private readonly byte[] _quality;

    /// <summary>Creates a read; <paramref name="phredQualities"/> must have one value per base.</summary>
    public SequenceRead(Sequence sequence, ReadOnlySpan<byte> phredQualities)
        : this(sequence, phredQualities.ToArray()) { }

    private SequenceRead(Sequence sequence, byte[] owned)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (owned.Length != sequence.Length)
            throw new ArgumentException($"Read '{sequence.Id}' has {sequence.Length} bases but {owned.Length} qualities.");
        Sequence = sequence;
        _quality = owned;
    }

    /// <summary>Wraps an owned quality buffer without copying (for parsers).</summary>
    public static SequenceRead FromBuffer(Sequence sequence, byte[] phredQualities) => new(sequence, phredQualities);

    /// <summary>Builds a read from residues and an ASCII quality string with the given Phred offset (33 for Sanger/Illumina 1.8+).</summary>
    public static SequenceRead FromAscii(string id, string residues, string quality, int offset = 33)
    {
        var q = new byte[quality.Length];
        for (int i = 0; i < q.Length; i++) q[i] = (byte)(quality[i] - offset);
        return new(Sequence.Dna(residues, id), q);
    }

    /// <summary>The bases.</summary>
    public Sequence Sequence { get; }

    /// <summary>Read identifier.</summary>
    public string Id => Sequence.Id;

    /// <summary>Read length.</summary>
    public int Length => Sequence.Length;

    /// <summary>Phred quality per base.</summary>
    public ReadOnlySpan<byte> Quality => _quality;

    /// <summary>Mean Phred quality (0 for an empty read).</summary>
    public double MeanQuality()
    {
        if (_quality.Length == 0) return 0;
        long s = 0;
        foreach (byte q in _quality) s += q;
        return (double)s / _quality.Length;
    }

    /// <summary>Expected number of base-call errors, Σ 10^(−Q/10).</summary>
    public double ExpectedErrors()
    {
        double e = 0;
        foreach (byte q in _quality) e += Math.Pow(10, -q / 10.0);
        return e;
    }

    /// <summary>The quality string encoded with <paramref name="offset"/>.</summary>
    public string QualityString(int offset = 33) => string.Create(_quality.Length, (_quality, offset),
        static (s, st) => { for (int i = 0; i < s.Length; i++) s[i] = (char)(st._quality[i] + st.offset); });

    /// <summary>A copy of the read between <paramref name="start"/> and <paramref name="start"/> + <paramref name="length"/>.</summary>
    public SequenceRead Slice(int start, int length)
        => new(Sequence.Subsequence(start, length), _quality.AsSpan(start, length).ToArray());

    /// <summary>
    /// Trims both ends with a sliding window: drops leading/trailing bases until a window of
    /// <paramref name="window"/> bases reaches a mean quality of <paramref name="minQuality"/>, then drops any
    /// remaining edge bases that are individually below <paramref name="minQuality"/>.
    /// </summary>
    public SequenceRead TrimByQuality(int minQuality = 20, int window = 4)
    {
        int n = _quality.Length;
        if (n == 0) return this;
        window = Math.Clamp(window, 1, n);
        int start = 0, end = n;
        while (start + window <= end && WindowMean(start, window) < minQuality) start++;
        while (end - window >= start && WindowMean(end - window, window) < minQuality) end--;
        while (start < end && _quality[start] < minQuality) start++;
        while (end > start && _quality[end - 1] < minQuality) end--;
        if (end <= start) return Slice(0, 0);
        return start == 0 && end == n ? this : Slice(start, end - start);
    }

    private double WindowMean(int start, int window)
    {
        int s = 0;
        for (int i = start; i < start + window; i++) s += _quality[i];
        return (double)s / window;
    }

    public override string ToString() => $"@{Id} {Length} bp, mean Q{MeanQuality():F1}";
}
