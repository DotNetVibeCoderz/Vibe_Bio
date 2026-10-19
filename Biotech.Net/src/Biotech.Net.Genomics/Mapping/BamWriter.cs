using System.Buffers.Binary;
using System.IO.Compression;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>Options for <see cref="BamWriter"/>.</summary>
public sealed record BamWriterOptions
{
    /// <summary>Deflate level of the BGZF blocks (default <see cref="CompressionLevel.Optimal"/>).</summary>
    public CompressionLevel CompressionLevel { get; init; } = CompressionLevel.Optimal;

    /// <summary>
    /// Build an index while writing (records must then be coordinate-sorted). <see cref="BamWriter.Create"/> saves
    /// it next to the file (<c>.bai</c>, or <c>.csi</c> for <see cref="BinningIndexFormat.Csi"/>); null for none.
    /// </summary>
    public BinningIndexFormat? Index { get; init; }

    /// <summary>CSI bottom-level bin shift (default 14).</summary>
    public int CsiMinShift { get; init; } = BinningIndex.DefaultMinShift;

    /// <summary>CSI depth (default 5).</summary>
    public int CsiDepth { get; init; } = BinningIndex.DefaultDepth;
}

/// <summary>
/// BAM writer: the header is written on construction, records are encoded to the BAM binary layout (4-bit bases,
/// typed optional fields, bins computed as htslib) and compressed into BGZF blocks that never split a record
/// unnecessarily. Optionally builds a BAI/CSI index on the fly.
/// </summary>
public sealed class BamWriter : IDisposable
{
    private readonly BgzfWriter _bgzf;
    private readonly ByteBuffer _record = new(1024);
    private readonly string? _indexPath;
    private BinningIndexBuilder? _indexer;
    private bool _disposed;

    /// <summary>Creates a writer over <paramref name="stream"/> and writes the header.</summary>
    public BamWriter(Stream stream, SamHeader header, BamWriterOptions? options = null, bool leaveOpen = false)
    {
        Header = header ?? throw new ArgumentNullException(nameof(header));
        options ??= new BamWriterOptions();
        _bgzf = new BgzfWriter(stream, options.CompressionLevel, leaveOpen);
        BamCodec.WriteHeader(_bgzf, header);
        _bgzf.FlushBlock(); // records start in a fresh block, as htslib does
        if (options.Index is { } fmt)
        {
            if (fmt == BinningIndexFormat.Tbi) throw new ArgumentException("BAM files use BAI or CSI indexes.", nameof(options));
            _indexer = fmt == BinningIndexFormat.Csi
                ? new BinningIndexBuilder(fmt, _bgzf.VirtualPosition, header.References.Count, options.CsiMinShift, options.CsiDepth)
                : new BinningIndexBuilder(fmt, _bgzf.VirtualPosition, header.References.Count);
        }
    }

    private BamWriter(Stream stream, SamHeader header, BamWriterOptions options, string? indexPath) : this(stream, header, options)
        => _indexPath = indexPath;

    /// <summary>Creates a BAM file; with <see cref="BamWriterOptions.Index"/> the index is saved beside it on dispose.</summary>
    public static BamWriter Create(string path, SamHeader header, BamWriterOptions? options = null)
    {
        options ??= new BamWriterOptions();
        string? idx = options.Index switch
        {
            BinningIndexFormat.Bai => path + ".bai",
            BinningIndexFormat.Csi => path + ".csi",
            _ => null,
        };
        return new BamWriter(new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 16), header, options, idx);
    }

    /// <summary>The header written.</summary>
    public SamHeader Header { get; }

    /// <summary>The index built while writing (available after <see cref="Dispose"/>), or null.</summary>
    public BinningIndex? Index { get; private set; }

    /// <summary>Writes one record.</summary>
    public void Write(AlignmentRecord record)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        _record.Clear();
        _record.Append(stackalloc byte[4]);
        BamCodec.Encode(record, Header, _record, out int refId, out int end);
        _record.SetInt32(0, _record.Length - 4);
        var span = _record.Span;
        _bgzf.EnsureSpace(span.Length);
        _bgzf.Write(span);
        _indexer?.Push(refId, record.Position, end, _bgzf.VirtualPosition, !record.IsUnmapped);
    }

    /// <summary>Writes records.</summary>
    public void Write(IEnumerable<AlignmentRecord> records)
    {
        foreach (var r in records) Write(r);
    }

    /// <summary>Writes a <see cref="ReadAlignment"/> from <see cref="ReadMapper"/>.</summary>
    public void Write(ReadAlignment alignment) => Write(alignment.ToAlignmentRecord(Header));

    /// <summary>Writes mapper output.</summary>
    public void Write(IEnumerable<ReadAlignment> alignments)
    {
        foreach (var a in alignments) Write(a);
    }

    /// <summary>Finishes the BGZF stream (EOF marker), completes and saves the index if requested.</summary>
    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _bgzf.FlushBlock();
        var final = _bgzf.VirtualPosition;
        _bgzf.Dispose();
        if (_indexer is not null)
        {
            Index = _indexer.Finish(final);
            _indexer = null;
            if (_indexPath is not null) Index.Save(_indexPath);
        }
    }
}

/// <summary>Builds BAI/CSI indexes for coordinate-sorted BAM files (like <c>samtools index</c>).</summary>
public static class BamIndexer
{
    /// <summary>Indexes the BAM data in <paramref name="bam"/> (read from its start).</summary>
    public static BinningIndex Build(Stream bam, BinningIndexFormat format = BinningIndexFormat.Bai,
        int minShift = BinningIndex.DefaultMinShift, int depth = BinningIndex.DefaultDepth)
    {
        if (format == BinningIndexFormat.Tbi) throw new ArgumentException("BAM files use BAI or CSI indexes.", nameof(format));
        using var reader = new BamReader(bam, leaveOpen: true);
        var builder = format == BinningIndexFormat.Csi
            ? new BinningIndexBuilder(format, reader.Position, reader.Header.References.Count, minShift, depth)
            : new BinningIndexBuilder(format, reader.Position, reader.Header.References.Count);
        while (true)
        {
            var raw = reader.ReadRaw();
            if (raw.IsEmpty) break;
            int refId = BinaryPrimitives.ReadInt32LittleEndian(raw);
            int pos = BinaryPrimitives.ReadInt32LittleEndian(raw[4..]);
            int flag = BinaryPrimitives.ReadUInt16LittleEndian(raw[14..]);
            int end = pos + Math.Max(1, (flag & 4) != 0 ? 0 : RawReferenceLength(raw));
            builder.Push(refId, pos, end, reader.Position, (flag & 4) == 0);
        }
        return builder.Finish(reader.Position);
    }

    /// <summary>Indexes <paramref name="bamPath"/> and saves <c>.bai</c> (or <c>.csi</c>) next to it; returns the index path.</summary>
    public static string Index(string bamPath, BinningIndexFormat format = BinningIndexFormat.Bai, int minShift = BinningIndex.DefaultMinShift, int depth = BinningIndex.DefaultDepth)
    {
        BinningIndex idx;
        using (var fs = File.OpenRead(bamPath)) idx = Build(fs, format, minShift, depth);
        string path = bamPath + (format == BinningIndexFormat.Csi ? ".csi" : ".bai");
        idx.Save(path);
        return path;
    }

    private static int RawReferenceLength(ReadOnlySpan<byte> raw)
    {
        int lName = raw[8];
        int nCigar = BinaryPrimitives.ReadUInt16LittleEndian(raw[12..]);
        int rlen = 0;
        for (int i = 0, p = 32 + lName; i < nCigar; i++, p += 4)
        {
            uint c = BinaryPrimitives.ReadUInt32LittleEndian(raw[p..]);
            // M D N = X consume the reference (codes 0 2 3 7 8)
            if (((0x18D >> (int)(c & 0xF)) & 1) != 0) rlen += (int)(c >> 4);
        }
        return rlen;
    }
}
