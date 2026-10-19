using System.Buffers.Binary;
using System.IO.Compression;

namespace Biotech.Net.Genomics.IO;

/// <summary>
/// Writer for BGZF: data is cut into blocks of at most 65 280 uncompressed bytes (as htslib), each compressed as an
/// independent gzip member carrying its size in a <c>BC</c> extra field, and the stream ends with the standard
/// 28-byte empty EOF block. The output is valid gzip, so any gzip reader can decompress it.
/// </summary>
public sealed class BgzfWriter : Stream
{
    /// <summary>Uncompressed bytes per block (htslib's <c>BGZF_BLOCK_SIZE</c>).</summary>
    public const int BlockSize = 0xff00;

    /// <summary>The 28-byte BGZF end-of-file marker.</summary>
    public static ReadOnlySpan<byte> EofMarker =>
    [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    private readonly Stream _stream;
    private readonly bool _leaveOpen;
    private readonly CompressionLevel _level;
    private readonly byte[] _buffer = new byte[BlockSize];
    private readonly MemoryStream _compressed = new(BgzfReader.MaxBlockSize);
    private int _length;
    private long _address;
    private bool _disposed;

    /// <summary>Creates a writer over <paramref name="stream"/>.</summary>
    /// <param name="stream">Destination.</param>
    /// <param name="level">Deflate level (default <see cref="CompressionLevel.Optimal"/>).</param>
    /// <param name="leaveOpen">Leave <paramref name="stream"/> open when disposed.</param>
    public BgzfWriter(Stream stream, CompressionLevel level = CompressionLevel.Optimal, bool leaveOpen = false)
    {
        _stream = stream ?? throw new ArgumentNullException(nameof(stream));
        _level = level;
        _leaveOpen = leaveOpen;
        _address = stream.CanSeek ? stream.Position : 0;
    }

    /// <summary>Creates a BGZF file.</summary>
    public static BgzfWriter Create(string path, CompressionLevel level = CompressionLevel.Optimal)
        => new(new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 16), level);

    /// <summary>Virtual offset of the next byte to be written.</summary>
    public VirtualOffset VirtualPosition => new(_address, _length);

    /// <summary>Starts a new block when fewer than <paramref name="bytes"/> bytes fit in the current one, so a
    /// record of that size is not split across blocks (as htslib does for BAM records).</summary>
    public void EnsureSpace(int bytes)
    {
        if (_length > 0 && _length + bytes > BlockSize) FlushBlock();
    }

    /// <inheritdoc/>
    public override void Write(byte[] buffer, int offset, int count) => Write(buffer.AsSpan(offset, count));

    /// <inheritdoc/>
    public override void Write(ReadOnlySpan<byte> data)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        while (!data.IsEmpty)
        {
            int n = Math.Min(BlockSize - _length, data.Length);
            data[..n].CopyTo(_buffer.AsSpan(_length));
            _length += n;
            data = data[n..];
            if (_length == BlockSize) FlushBlock();
        }
    }

    /// <inheritdoc/>
    public override void WriteByte(byte value)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        _buffer[_length++] = value;
        if (_length == BlockSize) FlushBlock();
    }

    /// <summary>Ends the current block (if not empty) and flushes the underlying stream.</summary>
    public override void Flush()
    {
        FlushBlock();
        _stream.Flush();
    }

    /// <summary>Compresses and writes the buffered bytes as one block.</summary>
    public void FlushBlock()
    {
        if (_length == 0) return;
        var data = _buffer.AsSpan(0, _length);
        Compress(data, _level);
        if (_compressed.Length > BgzfReader.MaxBlockSize - 26) Compress(data, CompressionLevel.NoCompression);
        int cdata = (int)_compressed.Length;
        int total = cdata + 26;
        Span<byte> header = stackalloc byte[18];
        header[0] = 0x1f; header[1] = 0x8b; header[2] = 8; header[3] = 4;
        header[4] = header[5] = header[6] = header[7] = 0;
        header[8] = 0; header[9] = 0xff;
        BinaryPrimitives.WriteUInt16LittleEndian(header[10..], 6);
        header[12] = (byte)'B'; header[13] = (byte)'C';
        BinaryPrimitives.WriteUInt16LittleEndian(header[14..], 2);
        BinaryPrimitives.WriteUInt16LittleEndian(header[16..], (ushort)(total - 1));
        _stream.Write(header);
        _stream.Write(_compressed.GetBuffer().AsSpan(0, cdata));
        Span<byte> trailer = stackalloc byte[8];
        BinaryPrimitives.WriteUInt32LittleEndian(trailer, Crc32.Compute(data));
        BinaryPrimitives.WriteUInt32LittleEndian(trailer[4..], (uint)_length);
        _stream.Write(trailer);
        _address += total;
        _length = 0;
    }

    private void Compress(ReadOnlySpan<byte> data, CompressionLevel level)
    {
        _compressed.SetLength(0);
        using var ds = new DeflateStream(_compressed, level, leaveOpen: true);
        ds.Write(data);
    }

    public override bool CanRead => false;
    public override bool CanSeek => false;
    public override bool CanWrite => !_disposed;
    public override long Length => throw new NotSupportedException();
    public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
    public override int Read(byte[] buffer, int offset, int count) => throw new NotSupportedException();
    public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
    public override void SetLength(long value) => throw new NotSupportedException();

    /// <summary>Writes the last block and the EOF marker.</summary>
    protected override void Dispose(bool disposing)
    {
        if (!_disposed && disposing)
        {
            FlushBlock();
            _stream.Write(EofMarker);
            _address += EofMarker.Length;
            _stream.Flush();
            if (!_leaveOpen) _stream.Dispose();
        }
        _disposed = true;
        base.Dispose(disposing);
    }
}
