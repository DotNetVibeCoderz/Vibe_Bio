using System.Text;
using Biotech.Net.Native;

namespace Biotech.Net.Core;

/// <summary>
/// An immutable biological sequence (DNA, RNA or protein) stored as ASCII bytes. Residues are kept in the
/// case they were given (lower case is commonly soft-masking); every operation treats case insensitively.
/// </summary>
/// <remarks>
/// Byte storage halves memory against <see cref="string"/> and lets the Rust kernels read the buffer in place.
/// Operations that return a new sequence allocate exactly one array of the result length.
/// </remarks>
public sealed class Sequence : IEquatable<Sequence>
{
    private readonly byte[] _data;
    private List<SequenceFeature>? _features;
    private Dictionary<string, string>? _annotations;

    /// <summary>Creates a sequence by copying <paramref name="residues"/>.</summary>
    /// <exception cref="ArgumentException">A residue is not in <paramref name="alphabet"/> and <paramref name="validate"/> is true.</exception>
    public Sequence(string id, ReadOnlySpan<byte> residues, Alphabet alphabet, string? description = null, bool validate = true)
        : this(id, residues.ToArray(), alphabet, description, validate)
    {
    }

    private Sequence(string id, byte[] owned, Alphabet alphabet, string? description, bool validate)
    {
        ArgumentNullException.ThrowIfNull(id);
        ArgumentNullException.ThrowIfNull(alphabet);
        if (validate)
        {
            int bad = alphabet.FirstInvalid(owned);
            if (bad >= 0)
                throw new ArgumentException($"Symbol '{(char)owned[bad]}' at position {bad} is not valid for {alphabet.Name}.", nameof(alphabet));
        }
        _data = owned;
        Id = id;
        Alphabet = alphabet;
        Description = description ?? string.Empty;
    }

    /// <summary>
    /// Wraps <paramref name="buffer"/> without copying. The caller hands over ownership and must not modify the
    /// array afterwards. Intended for parsers that build the residue buffer themselves.
    /// </summary>
    public static Sequence FromBuffer(string id, byte[] buffer, Alphabet alphabet, string? description = null, bool validate = true)
        => new(id, buffer, alphabet, description, validate);

    /// <summary>A DNA sequence (IUPAC codes allowed).</summary>
    public static Sequence Dna(string residues, string id = "dna", string? description = null)
        => new(id, Encoding.ASCII.GetBytes(residues), Alphabet.DnaIupac, description, validate: true);

    /// <summary>An RNA sequence (IUPAC codes allowed).</summary>
    public static Sequence Rna(string residues, string id = "rna", string? description = null)
        => new(id, Encoding.ASCII.GetBytes(residues), Alphabet.RnaIupac, description, validate: true);

    /// <summary>A protein sequence (one-letter codes, extended set allowed).</summary>
    public static Sequence Protein(string residues, string id = "protein", string? description = null)
        => new(id, Encoding.ASCII.GetBytes(residues), Alphabet.ProteinExtended, description, validate: true);

    /// <summary>A sequence whose type is detected from its residues (see <see cref="Alphabet.Detect"/>).</summary>
    public static Sequence Parse(string residues, string id = "seq", string? description = null)
    {
        var bytes = Encoding.ASCII.GetBytes(residues);
        return new(id, bytes, Alphabet.Detect(bytes), description, validate: false);
    }

    /// <summary>Identifier (FASTA header up to the first whitespace, GenBank accession, ...).</summary>
    public string Id { get; }

    /// <summary>Free-text description.</summary>
    public string Description { get; }

    /// <summary>The alphabet the residues belong to.</summary>
    public Alphabet Alphabet { get; }

    /// <summary>Polymer type.</summary>
    public SequenceType Type => Alphabet.Type;

    /// <summary>Number of residues.</summary>
    public int Length => _data.Length;

    /// <summary>Residue at <paramref name="index"/> as an ASCII byte.</summary>
    public byte this[int index] => _data[index];

    /// <summary>The residues, without copying.</summary>
    public ReadOnlySpan<byte> Span => _data;

    /// <summary>The residues as memory, for async and pipeline APIs.</summary>
    public ReadOnlyMemory<byte> Memory => _data;

    /// <summary>Annotated features (genes, CDS, ...), typically from GenBank/EMBL.</summary>
    public IList<SequenceFeature> Features => _features ??= [];

    /// <summary>Free-form key/value annotations (organism, molecule type, source file, ...).</summary>
    public IDictionary<string, string> Annotations => _annotations ??= new(StringComparer.OrdinalIgnoreCase);

    /// <summary>True when the sequence carries features.</summary>
    public bool HasFeatures => _features is { Count: > 0 };

    /// <summary>A copy of <paramref name="length"/> residues starting at <paramref name="start"/> (0-based).</summary>
    public Sequence Subsequence(int start, int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (start + length > _data.Length) throw new ArgumentOutOfRangeException(nameof(length));
        return new(Id, _data.AsSpan(start, length).ToArray(), Alphabet, Description, validate: false);
    }

    /// <summary>Same residues and metadata under a different identifier.</summary>
    public Sequence WithId(string id, string? description = null)
    {
        var s = new Sequence(id, _data, Alphabet, description ?? Description, validate: false);
        CopyMetadataTo(s);
        return s;
    }

    /// <summary>The IUPAC reverse complement (nucleotides only).</summary>
    public Sequence ReverseComplement()
    {
        RequireNucleotide();
        var dst = new byte[_data.Length];
        BioKernels.ReverseComplement(_data, dst);
        if (Type == SequenceType.Rna) TToU(dst);
        return new(Id, dst, Alphabet, Description, validate: false);
    }

    /// <summary>The IUPAC complement, not reversed (nucleotides only).</summary>
    public Sequence Complement()
    {
        RequireNucleotide();
        var dst = new byte[_data.Length];
        BioKernels.Complement(_data, dst);
        if (Type == SequenceType.Rna) TToU(dst);
        return new(Id, dst, Alphabet, Description, validate: false);
    }

    /// <summary>DNA → RNA (T → U).</summary>
    public Sequence Transcribe()
    {
        if (Type != SequenceType.Dna) throw new InvalidOperationException("Only DNA can be transcribed.");
        var dst = (byte[])_data.Clone();
        TToU(dst);
        return new(Id, dst, Alphabet.RnaIupac, Description, validate: false);
    }

    /// <summary>RNA → DNA (U → T).</summary>
    public Sequence BackTranscribe()
    {
        if (Type != SequenceType.Rna) throw new InvalidOperationException("Only RNA can be back-transcribed.");
        var dst = (byte[])_data.Clone();
        for (int i = 0; i < dst.Length; i++)
        {
            if (dst[i] == (byte)'U') dst[i] = (byte)'T';
            else if (dst[i] == (byte)'u') dst[i] = (byte)'t';
        }
        return new(Id, dst, Alphabet.DnaIupac, Description, validate: false);
    }

    /// <summary>
    /// Translates the sequence in reading frame <paramref name="frame"/> (0, 1 or 2) with <paramref name="code"/>
    /// (standard by default). Trailing partial codons are ignored. With <paramref name="toStop"/> translation ends
    /// before the first stop codon; otherwise stops appear as <c>*</c>.
    /// </summary>
    public Sequence Translate(GeneticCode? code = null, int frame = 0, bool toStop = false)
    {
        RequireNucleotide();
        if (frame is < 0 or > 2) throw new ArgumentOutOfRangeException(nameof(frame));
        code ??= GeneticCode.Standard;
        var src = _data.AsSpan(Math.Min(frame, _data.Length));
        var protein = new byte[src.Length / 3];
        BioKernels.Translate(src, code.CodonTable, protein);
        if (toStop)
        {
            int stop = Array.IndexOf(protein, (byte)'*');
            if (stop >= 0) Array.Resize(ref protein, stop);
        }
        return new(Id, protein, Alphabet.ProteinExtended, Description, validate: false);
    }

    /// <summary>Base composition (nucleotides).</summary>
    public NucleotideCounts NucleotideCounts() => BioKernels.CountNucleotides(_data);

    /// <summary>GC fraction over unambiguous bases (0–1).</summary>
    public double GcContent() => NucleotideCounts().GcFraction;

    /// <summary>Occurrences of each residue symbol (upper-cased).</summary>
    public IReadOnlyDictionary<char, long> Composition()
    {
        Span<long> h = stackalloc long[256];
        h.Clear();
        BioKernels.ByteHistogram(_data, h);
        var result = new SortedDictionary<char, long>();
        for (int b = 0; b < 256; b++)
        {
            if (h[b] == 0) continue;
            char c = char.ToUpperInvariant((char)b);
            result[c] = result.GetValueOrDefault(c) + h[b];
        }
        return result;
    }

    /// <summary>First 0-based index of <paramref name="motif"/> (exact, case-sensitive bytes) at or after <paramref name="start"/>, or −1.</summary>
    public int IndexOf(ReadOnlySpan<byte> motif, int start = 0)
    {
        int i = _data.AsSpan(start).IndexOf(motif);
        return i < 0 ? -1 : i + start;
    }

    /// <summary>First 0-based index of <paramref name="motif"/> (case-insensitive) at or after <paramref name="start"/>, or −1.</summary>
    public int IndexOf(string motif, int start = 0)
    {
        var haystack = _data.AsSpan(start);
        Span<byte> pattern = motif.Length <= 256 ? stackalloc byte[motif.Length] : new byte[motif.Length];
        Encoding.ASCII.GetBytes(motif, pattern);
        if (!HasLowerCase(haystack) && !HasLowerCase(pattern))
        {
            int j = haystack.IndexOf(pattern);
            return j < 0 ? -1 : j + start;
        }
        for (int i = 0; i + pattern.Length <= haystack.Length; i++)
            if (EqualsIgnoreCase(haystack.Slice(i, pattern.Length), pattern)) return i + start;
        return -1;
    }

    /// <summary>Whether the residues contain <paramref name="motif"/> (case-insensitive).</summary>
    public bool Contains(string motif) => IndexOf(motif) >= 0;

    /// <summary>Concatenates two sequences of the same type.</summary>
    public Sequence Concat(Sequence other, string? id = null)
    {
        if (other.Type != Type) throw new ArgumentException("Cannot join sequences of different types.", nameof(other));
        var dst = new byte[_data.Length + other._data.Length];
        _data.CopyTo(dst, 0);
        other._data.CopyTo(dst, _data.Length);
        return new(id ?? Id, dst, Alphabet, Description, validate: false);
    }

    /// <summary>Upper-cased copy (removes soft-masking).</summary>
    public Sequence ToUpper()
    {
        var dst = (byte[])_data.Clone();
        Ascii.ToUpperInPlace(dst, out _);
        return new(Id, dst, Alphabet, Description, validate: false);
    }

    /// <summary>The residues as a string.</summary>
    public override string ToString() => Encoding.ASCII.GetString(_data);

    /// <summary>A short preview, e.g. <c>ACGTACGT…TTGA (1,204 bp)</c>.</summary>
    public string Preview(int head = 30)
    {
        string unit = Type == SequenceType.Protein ? "aa" : "bp";
        if (_data.Length <= head * 2) return $"{this} ({_data.Length:N0} {unit})";
        return $"{Encoding.ASCII.GetString(_data, 0, head)}…{Encoding.ASCII.GetString(_data, _data.Length - 10, 10)} ({_data.Length:N0} {unit})";
    }

    /// <summary>Residue equality, case-insensitive; identifiers are ignored.</summary>
    public bool Equals(Sequence? other) => other is not null && other.Type == Type && EqualsIgnoreCase(_data, other._data);

    public override bool Equals(object? obj) => obj is Sequence s && Equals(s);

    public override int GetHashCode()
    {
        var h = new HashCode();
        h.Add(Type);
        foreach (byte b in _data) h.Add(b & 0xDF);
        return h.ToHashCode();
    }

    private void CopyMetadataTo(Sequence target)
    {
        if (_features is not null) foreach (var f in _features) target.Features.Add(f);
        if (_annotations is not null) foreach (var kv in _annotations) target.Annotations[kv.Key] = kv.Value;
    }

    private void RequireNucleotide()
    {
        if (!Alphabet.IsNucleotide)
            throw new InvalidOperationException($"This operation needs a nucleotide sequence; '{Id}' is {Type}.");
    }

    private static void TToU(Span<byte> s)
    {
        for (int i = 0; i < s.Length; i++)
        {
            if (s[i] == (byte)'T') s[i] = (byte)'U';
            else if (s[i] == (byte)'t') s[i] = (byte)'u';
        }
    }

    private static bool HasLowerCase(ReadOnlySpan<byte> s) => s.IndexOfAnyInRange((byte)'a', (byte)'z') >= 0;

    private static bool EqualsIgnoreCase(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b)
        => a.Length == b.Length && Ascii.EqualsIgnoreCase(a, b);
}
