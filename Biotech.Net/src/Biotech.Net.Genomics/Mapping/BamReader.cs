using System.Buffers.Binary;
using Biotech.Net.Core;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>A source of alignments (SAM, BAM or CRAM) with a header and, when indexed, region queries.</summary>
public interface IAlignmentReader : IDisposable
{
    /// <summary>The file header.</summary>
    SamHeader Header { get; }

    /// <summary>Whether <see cref="Query(string, long, long)"/> can use an index.</summary>
    bool IsIndexed { get; }

    /// <summary>Streams the records from the current position to the end.</summary>
    IEnumerable<AlignmentRecord> Records();

    /// <summary>Records on reference <paramref name="referenceName"/> overlapping 0-based <c>[start, end)</c>.</summary>
    IEnumerable<AlignmentRecord> Query(string referenceName, long start, long end);
}

/// <summary>Region-query conveniences shared by every <see cref="IAlignmentReader"/>.</summary>
public static class AlignmentReaderExtensions
{
    /// <summary>Records overlapping a samtools-style region (<c>chr1</c>, <c>chr1:1000-2000</c>, 1-based inclusive).</summary>
    public static IEnumerable<AlignmentRecord> Query(this IAlignmentReader reader, string region)
    {
        var r = HtsRegion.Parse(region, n => reader.Header.GetReferenceId(n) >= 0);
        return reader.Query(r.Name, r.Start, r.End);
    }

    /// <summary>Records overlapping <paramref name="interval"/> (0-based half-open; strand ignored).</summary>
    public static IEnumerable<AlignmentRecord> Query(this IAlignmentReader reader, GenomicInterval interval)
        => reader.Query(interval.Chromosome, interval.Start, interval.End);
}

/// <summary>
/// Streaming BAM reader: decodes the header and records from BGZF, and answers region queries through a BAI or
/// CSI index (bins plus linear index, as htslib). One enumeration at a time: <see cref="Records"/> and the query
/// methods share the underlying stream.
/// </summary>
public sealed class BamReader : IAlignmentReader
{
    private readonly BgzfReader _bgzf;
    private byte[] _buffer = new byte[4096];

    /// <summary>Reads the header of the BAM data in <paramref name="stream"/>.</summary>
    public BamReader(Stream stream, bool leaveOpen = false)
    {
        _bgzf = new BgzfReader(stream, leaveOpen);
        Header = BamCodec.ReadHeader(_bgzf);
        FirstRecordOffset = _bgzf.VirtualPosition;
    }

    /// <summary>
    /// Opens a BAM file. The index is <paramref name="indexPath"/> when given, otherwise the first of
    /// <c>file.bam.bai</c>, <c>file.bai</c> and <c>file.bam.csi</c> that exists (none is fine for streaming).
    /// </summary>
    public static BamReader Open(string path, string? indexPath = null)
    {
        var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16, FileOptions.RandomAccess);
        var reader = new BamReader(fs);
        indexPath ??= FindIndex(path);
        if (indexPath is not null) reader.Index = BinningIndex.Load(indexPath);
        return reader;
    }

    internal static string? FindIndex(string path)
    {
        foreach (var candidate in new[] { path + ".bai", Path.ChangeExtension(path, ".bai"), path + ".csi", Path.ChangeExtension(path, ".csi") })
            if (File.Exists(candidate)) return candidate;
        return null;
    }

    /// <summary>The file header.</summary>
    public SamHeader Header { get; }

    /// <summary>The BAI or CSI index used for queries, or null.</summary>
    public BinningIndex? Index { get; set; }

    /// <inheritdoc/>
    public bool IsIndexed => Index is not null;

    /// <summary>Virtual offset of the first record (just after the header).</summary>
    public VirtualOffset FirstRecordOffset { get; }

    /// <summary>Virtual offset of the next record.</summary>
    public VirtualOffset Position => _bgzf.VirtualPosition;

    /// <summary>Moves to a record boundary (e.g. one saved from <see cref="Position"/>).</summary>
    public void Seek(VirtualOffset offset) => _bgzf.Seek(offset);

    /// <summary>Reads the next record, or null at the end of the file.</summary>
    public AlignmentRecord? Read()
    {
        var raw = ReadRaw();
        return raw.IsEmpty ? null : BamCodec.Decode(raw, Header);
    }

    /// <summary>The next record's raw bytes (after <c>block_size</c>), or empty at the end; valid until the next read.</summary>
    internal ReadOnlySpan<byte> ReadRaw()
    {
        Span<byte> four = stackalloc byte[4];
        int n = _bgzf.Read(four);
        if (n == 0) return default;
        if (n < 4 && _bgzf.Read(four[n..]) + n < 4) throw new EndOfStreamException("Truncated BAM record length.");
        int size = BinaryPrimitives.ReadInt32LittleEndian(four);
        if (size < 32) throw new InvalidDataException($"Invalid BAM record size {size}.");
        if (_buffer.Length < size) _buffer = new byte[Math.Max(size, _buffer.Length * 2)];
        if (_bgzf.ReadAtLeast(_buffer.AsSpan(0, size), size, throwOnEndOfStream: false) != size)
            throw new EndOfStreamException("Truncated BAM record.");
        return _buffer.AsSpan(0, size);
    }

    /// <inheritdoc/>
    public IEnumerable<AlignmentRecord> Records()
    {
        AlignmentRecord? r;
        while ((r = Read()) is not null) yield return r;
    }

    /// <summary>Streams every record from the start of the file.</summary>
    public IEnumerable<AlignmentRecord> AllRecords()
    {
        Seek(FirstRecordOffset);
        return Records();
    }

    /// <inheritdoc/>
    /// <exception cref="InvalidOperationException">No index is loaded.</exception>
    public IEnumerable<AlignmentRecord> Query(string referenceName, long start, long end)
    {
        var index = Index ?? throw new InvalidOperationException("Region queries need a BAI or CSI index (BamReader.Open finds file.bam.bai; BamIndexer.Build creates one).");
        int tid = Header.GetReferenceId(referenceName);
        if (tid < 0) return [];
        return QueryCore(index, tid, Math.Max(0, start), end);
    }

    private IEnumerable<AlignmentRecord> QueryCore(BinningIndex index, int tid, long start, long end)
    {
        foreach (var chunk in index.Query(tid, start, end))
        {
            if (_bgzf.VirtualPosition != chunk.Start) _bgzf.Seek(chunk.Start);
            while (_bgzf.VirtualPosition < chunk.End)
            {
                var raw = ReadRaw();
                if (raw.IsEmpty) yield break;
                int refId = BinaryPrimitives.ReadInt32LittleEndian(raw);
                int pos = BinaryPrimitives.ReadInt32LittleEndian(raw[4..]);
                if (refId != tid || pos >= end) yield break;
                var rec = BamCodec.Decode(raw, Header);
                if (rec.Overlaps(start, end)) yield return rec;
            }
        }
    }

    /// <summary>The unplaced unmapped records stored after all placed ones (needs an index).</summary>
    public IEnumerable<AlignmentRecord> QueryUnmapped()
    {
        var index = Index ?? throw new InvalidOperationException("QueryUnmapped needs a BAI or CSI index.");
        var start = index.UnplacedStart();
        if (start is null) yield break;
        _bgzf.Seek(start.Value);
        AlignmentRecord? r;
        while ((r = Read()) is not null)
            if (r.ReferenceId < 0) yield return r;
    }

    /// <summary>Disposes the underlying stream.</summary>
    public void Dispose() => _bgzf.Dispose();
}
