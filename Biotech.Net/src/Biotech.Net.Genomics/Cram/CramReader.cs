using System.Buffers.Binary;
using System.Globalization;
using System.IO.Compression;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Cram;

/// <summary>Options for <see cref="CramReader"/>.</summary>
public sealed record CramReaderOptions
{
    /// <summary>Recreate MD and NM tags that the writer dropped because they matched the reference (default true, as samtools).</summary>
    public bool RegenerateMdNm { get; init; } = true;

    /// <summary>Check each slice's reference MD5 against the reference bases used (default true).</summary>
    public bool ValidateReferenceMd5 { get; init; } = true;

    /// <summary>Prefix of generated read names when names were not stored (default: the file name, as htslib).</summary>
    public string? ReadNamePrefix { get; init; }
}

/// <summary>One line of a CRAM index (<c>.crai</c>): a slice and the reference span it covers.</summary>
/// <param name="ReferenceId">Reference id (-1 for unplaced reads).</param>
/// <param name="AlignmentStart">1-based start of the alignments in the slice.</param>
/// <param name="AlignmentSpan">Reference span.</param>
/// <param name="ContainerOffset">File offset of the container.</param>
/// <param name="SliceOffset">Offset of the slice after the container header (the landmark).</param>
/// <param name="SliceSize">Bytes of the slice.</param>
public sealed record CramIndexEntry(int ReferenceId, long AlignmentStart, long AlignmentSpan, long ContainerOffset, int SliceOffset, int SliceSize);

/// <summary>A CRAM index (<c>.crai</c>, gzipped text).</summary>
public sealed class CramIndex
{
    /// <summary>Creates an index from entries.</summary>
    public CramIndex(IEnumerable<CramIndexEntry> entries) => Entries = entries.ToList();

    /// <summary>The slices, in file order.</summary>
    public IReadOnlyList<CramIndexEntry> Entries { get; }

    /// <summary>Loads a <c>.crai</c> file.</summary>
    public static CramIndex Load(string path)
    {
        using var fs = File.OpenRead(path);
        return Read(fs);
    }

    /// <summary>Reads CRAI data (gzip or plain text).</summary>
    public static CramIndex Read(Stream stream)
    {
        var entries = new List<CramIndexEntry>();
        using var r = new StreamReader(BioFile.WrapDecompression(stream), Encoding.ASCII);
        string? line;
        var inv = CultureInfo.InvariantCulture;
        while ((line = r.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            if (f.Length < 6) throw new InvalidDataException($"Malformed CRAI line '{line}'.");
            entries.Add(new CramIndexEntry(int.Parse(f[0], inv), long.Parse(f[1], inv), long.Parse(f[2], inv),
                long.Parse(f[3], inv), int.Parse(f[4], inv), int.Parse(f[5], inv)));
        }
        return new CramIndex(entries);
    }

    /// <summary>Slices that may hold reads of <paramref name="refId"/> overlapping 1-based closed <c>[start, end]</c>.</summary>
    public IEnumerable<CramIndexEntry> Overlapping(int refId, long start, long end)
        => Entries.Where(e => e.ReferenceId == refId && e.AlignmentStart <= end && e.AlignmentStart + e.AlignmentSpan > start);
}

/// <summary>
/// CRAM 2.1/3.0/3.1 reader: file definition, containers, compression headers, slices and records, with reference
/// based sequence reconstruction from a FASTA (via its <c>.fai</c>) or from embedded references, mate resolution and
/// MD/NM regeneration as htslib does. Region queries use a <c>.crai</c> index when present.
/// </summary>
/// <remarks>
/// Block codecs: raw, gzip, bzip2, rANS 4x8, rANS Nx16, the adaptive arithmetic coder, fqzcomp and the name
/// tokeniser. LZMA blocks raise <see cref="NotSupportedException"/>. One enumeration at a time.
/// </remarks>
public sealed class CramReader : IAlignmentReader
{
    private readonly Stream _stream;
    private readonly bool _leaveOpen;
    private readonly CramReaderOptions _options;
    private readonly Func<string, byte[]?>? _referenceSource;
    private readonly Dictionary<int, RefWindow?> _refCache = [];
    private readonly string _prefix;
    private readonly long _firstContainer;

    /// <summary>Reads the file definition and SAM header from <paramref name="stream"/> (must be seekable).</summary>
    /// <param name="stream">CRAM data.</param>
    /// <param name="referenceSource">Returns the bases of a reference by name (null when unknown); needed for reference-based CRAM.</param>
    /// <param name="options">Decoding options.</param>
    /// <param name="leaveOpen">Leave <paramref name="stream"/> open on dispose.</param>
    public CramReader(Stream stream, Func<string, byte[]?>? referenceSource = null, CramReaderOptions? options = null, bool leaveOpen = false)
    {
        _stream = stream ?? throw new ArgumentNullException(nameof(stream));
        if (!stream.CanSeek) throw new ArgumentException("CRAM decoding needs a seekable stream.", nameof(stream));
        _leaveOpen = leaveOpen;
        _options = options ?? new CramReaderOptions();
        _referenceSource = referenceSource;
        _prefix = _options.ReadNamePrefix ?? "";
        Span<byte> def = stackalloc byte[26];
        if (stream.ReadAtLeast(def, 26, throwOnEndOfStream: false) != 26 || !def[..4].SequenceEqual("CRAM"u8))
            throw new InvalidDataException("Not a CRAM file (magic 'CRAM' missing).");
        MajorVersion = def[4];
        MinorVersion = def[5];
        if (MajorVersion is < 2 or > 3) throw new NotSupportedException($"CRAM version {MajorVersion}.{MinorVersion} is not supported (2.1, 3.0 and 3.1 are).");
        var h = ContainerHeader.Read(stream, MajorVersion) ?? throw new InvalidDataException("CRAM file has no header container.");
        var body = new byte[h.Length];
        stream.ReadExactly(body);
        var c = new ByteCursor(body);
        var b = CramBlock.Read(ref c, MajorVersion);
        if (b.Data.Length < 4) throw new InvalidDataException("CRAM header block is too short.");
        int textLen = BinaryPrimitives.ReadInt32LittleEndian(b.Data);
        Header = new SamHeader(Encoding.UTF8.GetString(b.Data, 4, Math.Min(textLen, b.Data.Length - 4)));
        _firstContainer = stream.Position;
    }

    /// <summary>
    /// Opens a CRAM file. <paramref name="referencePath"/> is a FASTA (indexed with <c>.fai</c>, built if missing)
    /// holding the references named in the header. The index is <c>file.cram.crai</c> when it exists.
    /// </summary>
    public static CramReader Open(string path, string? referencePath = null, CramReaderOptions? options = null)
    {
        Func<string, byte[]?>? source = null;
        IndexedFastaReader? fasta = null;
        if (referencePath is not null)
        {
            fasta = IndexedFastaReader.Open(referencePath, saveIndex: false);
            source = name => fasta.Index.TryGet(name, out _) ? fasta.FetchAll(name).Span.ToArray() : null;
        }
        options ??= new CramReaderOptions();
        if (options.ReadNamePrefix is null) options = options with { ReadNamePrefix = Path.GetFileName(path) };
        var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16, FileOptions.RandomAccess);
        var reader = new CramReader(fs, source, options) { _fasta = fasta };
        if (File.Exists(path + ".crai")) reader.Index = CramIndex.Load(path + ".crai");
        return reader;
    }

    /// <summary>Opens a CRAM stream whose references are given as sequences.</summary>
    public static CramReader Open(Stream stream, IEnumerable<Sequence> references, CramReaderOptions? options = null)
    {
        var map = references.ToDictionary(r => r.Id, r => r.Span.ToArray(), StringComparer.Ordinal);
        return new CramReader(stream, n => map.TryGetValue(n, out var b) ? b : null, options);
    }

    private IndexedFastaReader? _fasta;

    /// <summary>CRAM major version.</summary>
    public int MajorVersion { get; }

    /// <summary>CRAM minor version.</summary>
    public int MinorVersion { get; }

    /// <inheritdoc/>
    public SamHeader Header { get; }

    /// <summary>The <c>.crai</c> index used for queries, or null.</summary>
    public CramIndex? Index { get; set; }

    /// <inheritdoc/>
    public bool IsIndexed => Index is not null;

    private RefWindow? GetReference(int refId)
    {
        if (_refCache.TryGetValue(refId, out var w)) return w;
        string name = Header.GetReferenceName(refId);
        var bases = _referenceSource?.Invoke(name);
        if (bases is not null)
        {
            var up = new byte[bases.Length];
            for (int i = 0; i < up.Length; i++) up[i] = bases[i] is >= (byte)'a' and <= (byte)'z' ? (byte)(bases[i] - 32) : bases[i];
            long len = Header.References[refId].Length;
            w = new RefWindow(up, 1, Math.Min(len > 0 ? len : up.Length, up.Length));
        }
        // keep only the most recent reference in memory (sorted files move through them in order)
        if (_refCache.Count > 0 && w is not null) _refCache.Clear();
        _refCache[refId] = w;
        return w;
    }

    private (ContainerHeader Header, byte[] Body, CompressionHeader? Comp)? ReadContainer(long offset)
    {
        _stream.Position = offset;
        var h = ContainerHeader.Read(_stream, MajorVersion);
        if (h is null || h.IsEof) return null;
        var body = new byte[h.Length];
        _stream.ReadExactly(body);
        if (h.NRecords == 0 && h.Landmarks.Length == 0) return (h, body, null);
        var c = new ByteCursor(body);
        var cb = CramBlock.Read(ref c, MajorVersion);
        if (cb.ContentType != CramBlock.CompressionHeader) throw new InvalidDataException("CRAM: container does not start with a compression header.");
        return (h, body, CompressionHeader.Parse(cb.Data));
    }

    private SliceDecoder Decoder(CompressionHeader ch) => new(Header, ch, GetReference, _options, _prefix);

    /// <inheritdoc/>
    public IEnumerable<AlignmentRecord> Records()
    {
        long offset = _firstContainer;
        while (true)
        {
            var ct = ReadContainer(offset);
            if (ct is null) yield break;
            var (h, body, ch) = ct.Value;
            offset = h.Offset + h.HeaderSize + h.Length;
            if (ch is null) continue;
            var dec = Decoder(ch);
            foreach (int landmark in h.Landmarks)
                foreach (var rec in dec.Decode(body, landmark, MajorVersion))
                    yield return rec;
        }
    }

    /// <inheritdoc/>
    /// <remarks>With a <c>.crai</c> index only the overlapping slices are decoded; otherwise containers whose span cannot overlap are skipped.</remarks>
    public IEnumerable<AlignmentRecord> Query(string referenceName, long start, long end)
    {
        int tid = Header.GetReferenceId(referenceName);
        if (tid < 0 || end <= start) yield break;
        if (start < 0) start = 0;
        if (Index is not null)
        {
            var slices = Index.Overlapping(tid, start + 1, end == long.MaxValue ? long.MaxValue : end)
                .Select(e => (e.ContainerOffset, e.SliceOffset)).Distinct().ToList();
            long current = -1;
            SliceDecoder? dec = null;
            byte[] body = [];
            foreach (var (co, so) in slices)
            {
                if (co != current)
                {
                    var ct = ReadContainer(co) ?? throw new InvalidDataException($"CRAM index points at a missing container ({co}).");
                    body = ct.Body;
                    dec = Decoder(ct.Comp ?? throw new InvalidDataException("CRAM index points at an empty container."));
                    current = co;
                }
                foreach (var rec in dec!.Decode(body, so, MajorVersion))
                {
                    if (rec.ReferenceId != tid || rec.Position >= end) continue;
                    if (rec.Overlaps(start, end)) yield return rec;
                }
            }
            yield break;
        }
        long offset = _firstContainer;
        while (true)
        {
            var ct = ReadContainer(offset);
            if (ct is null) yield break;
            var (h, body, ch) = ct.Value;
            offset = h.Offset + h.HeaderSize + h.Length;
            if (ch is null) continue;
            if (h.RefId >= 0 && (h.RefId != tid || h.Start > end || h.Start + h.Span < start + 1)) continue;
            if (h.RefId == -1) continue;
            var dec = Decoder(ch);
            foreach (int landmark in h.Landmarks)
                foreach (var rec in dec.Decode(body, landmark, MajorVersion))
                    if (rec.ReferenceId == tid && rec.Overlaps(start, end)) yield return rec;
        }
    }

    /// <summary>Builds a <c>.crai</c>-equivalent index by scanning container and slice headers (multi-reference slices are expanded per reference).</summary>
    public CramIndex BuildIndex()
    {
        var entries = new List<CramIndexEntry>();
        long offset = _firstContainer;
        while (true)
        {
            var ct = ReadContainer(offset);
            if (ct is null) break;
            var (h, body, ch) = ct.Value;
            offset = h.Offset + h.HeaderSize + h.Length;
            if (ch is null) continue;
            for (int i = 0; i < h.Landmarks.Length; i++)
            {
                int so = h.Landmarks[i];
                int size = (i + 1 < h.Landmarks.Length ? h.Landmarks[i + 1] : h.Length) - so;
                var c = new ByteCursor(body) { Pos = so };
                var sb = CramBlock.Read(ref c, MajorVersion);
                var sh = new ByteCursor(sb.Data);
                int refId = sh.Itf8();
                long st = sh.Itf8(), sp = sh.Itf8();
                if (refId != -2)
                {
                    // htslib writes unplaced slices as "-1 0 1"
                    entries.Add(refId < 0 ? new CramIndexEntry(-1, 0, 1, h.Offset, so, size) : new CramIndexEntry(refId, st, sp, h.Offset, so, size));
                    continue;
                }
                // multi-reference slice: one entry per reference, spans from the decoded records
                var spans = new SortedDictionary<int, (long Start, long End)>();
                foreach (var rec in Decoder(ch).Decode(body, so, MajorVersion))
                {
                    int id = rec.ReferenceId;
                    long a = rec.Position + 1, e = rec.IsUnmapped ? a : rec.Position + Math.Max(1, rec.ReferenceLength);
                    if (id < 0) { a = 0; e = 0; }
                    spans[id] = spans.TryGetValue(id, out var x) ? (Math.Min(x.Start, a), Math.Max(x.End, e)) : (a, e);
                }
                foreach (var (id, (a, e)) in spans.Where(kv => kv.Key >= 0).Concat(spans.Where(kv => kv.Key < 0)))
                    entries.Add(id < 0 ? new CramIndexEntry(-1, 0, 1, h.Offset, so, size) : new CramIndexEntry(id, a, e - a + 1, h.Offset, so, size));
            }
        }
        return new CramIndex(entries);
    }

    /// <summary>Closes the file (and the reference FASTA).</summary>
    public void Dispose()
    {
        if (!_leaveOpen) _stream.Dispose();
        _fasta?.Dispose();
    }
}
