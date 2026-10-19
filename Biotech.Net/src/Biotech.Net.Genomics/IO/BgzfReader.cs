using System.Buffers.Binary;
using System.IO.Compression;

namespace Biotech.Net.Genomics.IO;

/// <summary>
/// Reader for BGZF (blocked gzip, the container of BAM, tabix-indexed VCF/BED/GFF and bgzipped FASTA).
/// Decompresses one block at a time, checks every block's CRC32 and ISIZE, and supports random access through
/// <see cref="VirtualOffset"/>s with a small LRU cache of decompressed blocks.
/// </summary>
/// <remarks>
/// <see cref="VirtualPosition"/> follows htslib's <c>bgzf_tell</c>: once a block is fully consumed the position is
/// reported as the start of the next block, so offsets recorded while reading match those in indexes built by htslib.
/// </remarks>
public sealed class BgzfReader : Stream
{
    /// <summary>Largest BGZF block, compressed or not.</summary>
    public const int MaxBlockSize = 65536;

    private readonly Stream _stream;
    private readonly bool _leaveOpen;
    private readonly int _cacheCapacity;
    private readonly Dictionary<long, LinkedListNode<Block>> _cache = [];
    private readonly LinkedList<Block> _lru = new();
    private readonly byte[] _header = new byte[18];
    private byte[] _compressed = new byte[MaxBlockSize];
    private Block? _block;
    private int _pos;
    private long _nextAddress;
    private long _streamPos;
    private bool _disposed;

    private sealed class Block
    {
        public long Address;
        public long Next;
        public int Length;
        public byte[] Data = [];
    }

    /// <summary>Wraps <paramref name="stream"/>, which must be positioned at the start of a BGZF block.</summary>
    /// <param name="stream">BGZF data; must be seekable for <see cref="Seek(VirtualOffset)"/>.</param>
    /// <param name="leaveOpen">Leave <paramref name="stream"/> open when this reader is disposed.</param>
    /// <param name="cacheBlocks">Decompressed blocks kept for random access (default 32, about 2 MB).</param>
    public BgzfReader(Stream stream, bool leaveOpen = false, int cacheBlocks = 32)
    {
        _stream = stream ?? throw new ArgumentNullException(nameof(stream));
        _leaveOpen = leaveOpen;
        _cacheCapacity = Math.Max(2, cacheBlocks);
        _streamPos = stream.CanSeek ? stream.Position : 0;
        _nextAddress = _streamPos;
    }

    /// <summary>Opens a BGZF file for reading.</summary>
    public static BgzfReader Open(string path, int cacheBlocks = 32)
        => new(new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16, FileOptions.RandomAccess), cacheBlocks: cacheBlocks);

    /// <summary>Whether the stream starts with a BGZF block header (gzip with a <c>BC</c> extra subfield). The position is restored.</summary>
    public static bool IsBgzf(Stream stream)
    {
        if (!stream.CanSeek) throw new ArgumentException("The stream must be seekable.", nameof(stream));
        long pos = stream.Position;
        try
        {
            Span<byte> h = stackalloc byte[18];
            return stream.ReadAtLeast(h, 18, throwOnEndOfStream: false) == 18 && IsBgzfHeader(h);
        }
        finally { stream.Position = pos; }
    }

    /// <summary>Whether <paramref name="path"/> is a BGZF file.</summary>
    public static bool IsBgzf(string path)
    {
        using var fs = File.OpenRead(path);
        return IsBgzf(fs);
    }

    private static bool IsBgzfHeader(ReadOnlySpan<byte> h)
        => h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && (h[3] & 4) != 0 && h[12] == (byte)'B' && h[13] == (byte)'C' && h[14] == 2 && h[15] == 0;

    /// <summary>Check the CRC32 of every block (default true).</summary>
    public bool VerifyChecksums { get; init; } = true;

    /// <summary>The current virtual offset (htslib <c>bgzf_tell</c> semantics).</summary>
    public VirtualOffset VirtualPosition
    {
        get
        {
            if (_block is null) return new VirtualOffset(_nextAddress, 0);
            if (_pos >= _block.Length) return new VirtualOffset(_block.Next, 0);
            return new VirtualOffset(_block.Address, _pos);
        }
    }

    /// <summary>Moves to <paramref name="offset"/> (requires a seekable stream).</summary>
    public void Seek(VirtualOffset offset)
    {
        var b = LoadBlock(offset.BlockAddress) ?? throw new EndOfStreamException($"No BGZF block at file offset {offset.BlockAddress}.");
        if (offset.BlockOffset > b.Length)
            throw new InvalidDataException($"Virtual offset {offset} lies beyond the {b.Length}-byte block.");
        _block = b;
        _pos = offset.BlockOffset;
        _nextAddress = b.Next;
    }

    /// <summary>Whether the end of the data has been reached (no further non-empty block).</summary>
    public bool AtEnd => (_block is null || _pos >= _block.Length) && !Advance(peekOnly: true);

    /// <inheritdoc/>
    public override int Read(byte[] buffer, int offset, int count) => Read(buffer.AsSpan(offset, count));

    /// <inheritdoc/>
    public override int Read(Span<byte> buffer)
    {
        int total = 0;
        while (!buffer.IsEmpty)
        {
            if (_block is null || _pos >= _block.Length)
                if (!Advance(peekOnly: false)) break;
            int n = Math.Min(buffer.Length, _block!.Length - _pos);
            _block.Data.AsSpan(_pos, n).CopyTo(buffer);
            _pos += n;
            total += n;
            buffer = buffer[n..];
        }
        return total;
    }

    /// <inheritdoc/>
    public override int ReadByte()
    {
        if (_block is null || _pos >= _block.Length)
            if (!Advance(peekOnly: false)) return -1;
        return _block!.Data[_pos++];
    }

    /// <summary>
    /// Reads one line (without the <c>\n</c>) into <paramref name="line"/>; returns false at end of data.
    /// A final line without a terminator is returned too.
    /// </summary>
    public bool ReadLine(ByteBuffer line)
    {
        line.Clear();
        bool any = false;
        while (true)
        {
            if (_block is null || _pos >= _block.Length)
            {
                if (!Advance(peekOnly: false)) return any;
            }
            any = true;
            var span = _block!.Data.AsSpan(_pos, _block.Length - _pos);
            int nl = span.IndexOf((byte)'\n');
            if (nl >= 0)
            {
                line.Append(span[..nl]);
                _pos += nl + 1;
                return true;
            }
            line.Append(span);
            _pos = _block.Length;
        }
    }

    /// <summary>Moves to the next non-empty block; false at end of data.</summary>
    private bool Advance(bool peekOnly)
    {
        long addr = _block is null ? _nextAddress : _block.Next;
        while (true)
        {
            var b = LoadBlock(addr);
            if (b is null) return false;
            if (b.Length > 0)
            {
                if (!peekOnly)
                {
                    _block = b;
                    _pos = 0;
                    _nextAddress = b.Next;
                }
                return true;
            }
            addr = b.Next;
        }
    }

    private Block? LoadBlock(long address)
    {
        if (_cache.TryGetValue(address, out var node))
        {
            _lru.Remove(node);
            _lru.AddFirst(node);
            return node.Value;
        }
        if (_streamPos != address)
        {
            if (!_stream.CanSeek) throw new NotSupportedException("Random access needs a seekable BGZF stream.");
            _stream.Position = address;
            _streamPos = address;
        }
        int got = _stream.ReadAtLeast(_header, 18, throwOnEndOfStream: false);
        _streamPos += got;
        if (got == 0) return null;
        if (got < 18) throw new InvalidDataException($"Truncated BGZF block header at file offset {address}.");
        var h = _header.AsSpan();
        if (h[0] != 0x1f || h[1] != 0x8b || h[2] != 8 || (h[3] & 4) == 0)
            throw new InvalidDataException($"Not a BGZF block at file offset {address} (gzip header with extra field expected).");
        int xlen = BinaryPrimitives.ReadUInt16LittleEndian(h[10..]);
        // The fixed header read above holds the first 6 bytes of the extra field; read the rest.
        int extraRest = xlen - 6;
        if (extraRest < 0) throw new InvalidDataException($"BGZF block at file offset {address} has a too-short extra field.");
        Span<byte> extra = xlen <= 256 ? stackalloc byte[xlen] : new byte[xlen];
        h[12..18].CopyTo(extra);
        if (extraRest > 0)
        {
            ReadFully(extra[6..]);
        }
        int bsize = -1;
        for (int i = 0; i + 4 <= xlen;)
        {
            int slen = BinaryPrimitives.ReadUInt16LittleEndian(extra[(i + 2)..]);
            if (extra[i] == (byte)'B' && extra[i + 1] == (byte)'C' && slen == 2 && i + 6 <= xlen)
                bsize = BinaryPrimitives.ReadUInt16LittleEndian(extra[(i + 4)..]) + 1;
            i += 4 + slen;
        }
        if (bsize < 0) throw new InvalidDataException($"Gzip member at file offset {address} has no BGZF 'BC' block-size field (plain gzip, not BGZF).");
        int cdata = bsize - xlen - 20;
        if (cdata < 0) throw new InvalidDataException($"Invalid BGZF block size {bsize} at file offset {address}.");
        if (_compressed.Length < cdata + 8) _compressed = new byte[cdata + 8];
        ReadFully(_compressed.AsSpan(0, cdata + 8));
        uint crc = BinaryPrimitives.ReadUInt32LittleEndian(_compressed.AsSpan(cdata));
        int isize = (int)BinaryPrimitives.ReadUInt32LittleEndian(_compressed.AsSpan(cdata + 4));
        if (isize is < 0 or > MaxBlockSize) throw new InvalidDataException($"Invalid BGZF ISIZE {isize} at file offset {address}.");

        var block = NewBlock();
        block.Address = address;
        block.Next = address + bsize;
        block.Length = isize;
        if (isize > 0)
        {
            using (var ms = new MemoryStream(_compressed, 0, cdata, writable: false))
            using (var ds = new DeflateStream(ms, CompressionMode.Decompress))
            {
                int n = ds.ReadAtLeast(block.Data.AsSpan(0, isize), isize, throwOnEndOfStream: false);
                if (n != isize || ds.ReadByte() >= 0)
                    throw new InvalidDataException($"BGZF block at file offset {address}: inflated size does not match ISIZE {isize}.");
            }
            if (VerifyChecksums && Crc32.Compute(block.Data.AsSpan(0, isize)) != crc)
                throw new InvalidDataException($"BGZF block at file offset {address}: CRC32 mismatch.");
        }
        else if (crc != 0)
            throw new InvalidDataException($"BGZF block at file offset {address}: CRC32 mismatch for an empty block.");
        var entry = new LinkedListNode<Block>(block);
        _cache[address] = entry;
        _lru.AddFirst(entry);
        return block;
    }

    /// <summary>A block object to fill: a fresh one while the cache is filling, else the least recently used one.</summary>
    private Block NewBlock()
    {
        if (_cache.Count < _cacheCapacity) return new Block { Data = new byte[MaxBlockSize] };
        var last = _lru.Last!;
        if (ReferenceEquals(last.Value, _block)) { _lru.RemoveLast(); _lru.AddFirst(last); last = _lru.Last!; }
        _lru.Remove(last);
        _cache.Remove(last.Value.Address);
        return last.Value;
    }

    private void ReadFully(Span<byte> dest)
    {
        int n = _stream.ReadAtLeast(dest, dest.Length, throwOnEndOfStream: false);
        _streamPos += n;
        if (n < dest.Length) throw new EndOfStreamException("Truncated BGZF block.");
    }

    public override bool CanRead => !_disposed;
    public override bool CanSeek => false;
    public override bool CanWrite => false;
    public override long Length => throw new NotSupportedException();
    public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
    public override void Flush() { }
    public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException("Use Seek(VirtualOffset).");
    public override void SetLength(long value) => throw new NotSupportedException();
    public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();

    protected override void Dispose(bool disposing)
    {
        if (!_disposed && disposing && !_leaveOpen) _stream.Dispose();
        _disposed = true;
        base.Dispose(disposing);
    }
}

/// <summary>A growable byte buffer used for line and record scratch space.</summary>
public sealed class ByteBuffer
{
    private byte[] _data;

    /// <summary>Creates a buffer with the given initial capacity.</summary>
    public ByteBuffer(int capacity = 256) => _data = new byte[Math.Max(16, capacity)];

    /// <summary>Number of bytes held.</summary>
    public int Length { get; private set; }

    /// <summary>The bytes held.</summary>
    public ReadOnlySpan<byte> Span => _data.AsSpan(0, Length);

    /// <summary>The bytes held, as memory.</summary>
    public ReadOnlyMemory<byte> Memory => _data.AsMemory(0, Length);

    /// <summary>Empties the buffer (capacity is kept).</summary>
    public void Clear() => Length = 0;

    /// <summary>Appends bytes.</summary>
    public void Append(ReadOnlySpan<byte> bytes)
    {
        EnsureCapacity(Length + bytes.Length);
        bytes.CopyTo(_data.AsSpan(Length));
        Length += bytes.Length;
    }

    /// <summary>Appends one byte.</summary>
    public void Append(byte b)
    {
        if (Length == _data.Length) EnsureCapacity(Length + 1);
        _data[Length++] = b;
    }

    /// <summary>Overwrites four bytes at <paramref name="offset"/> with a little-endian integer.</summary>
    public void SetInt32(int offset, int value)
    {
        if ((uint)offset > (uint)(Length - 4)) throw new ArgumentOutOfRangeException(nameof(offset));
        BinaryPrimitives.WriteInt32LittleEndian(_data.AsSpan(offset), value);
    }

    private void EnsureCapacity(int n)
    {
        if (n <= _data.Length) return;
        Array.Resize(ref _data, Math.Max(n, _data.Length * 2));
    }

    /// <summary>A copy of the bytes held.</summary>
    public byte[] ToArray() => Span.ToArray();

    /// <summary>The bytes held as ASCII text.</summary>
    public override string ToString() => System.Text.Encoding.ASCII.GetString(Span);
}
