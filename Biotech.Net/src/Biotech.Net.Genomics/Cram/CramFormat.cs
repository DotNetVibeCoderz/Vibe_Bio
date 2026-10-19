using System.Buffers.Binary;
using System.IO.Compression;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Cram;

/// <summary>CRAM block compression methods (block header byte).</summary>
public enum CramBlockMethod
{
    /// <summary>Uncompressed.</summary>
    Raw = 0,
    /// <summary>gzip / deflate.</summary>
    Gzip = 1,
    /// <summary>bzip2.</summary>
    Bzip2 = 2,
    /// <summary>LZMA (xz container).</summary>
    Lzma = 3,
    /// <summary>rANS 4x8 (CRAM 3.0).</summary>
    Rans4x8 = 4,
    /// <summary>rANS Nx16 (CRAM 3.1).</summary>
    RansNx16 = 5,
    /// <summary>Adaptive arithmetic coder (CRAM 3.1).</summary>
    Arith = 6,
    /// <summary>fqzcomp quality codec (CRAM 3.1).</summary>
    FqzComp = 7,
    /// <summary>Name tokeniser (CRAM 3.1).</summary>
    NameTokenizer = 8,
}

/// <summary>A decompressed CRAM block.</summary>
internal sealed class CramBlock
{
    public CramBlockMethod Method;
    public int ContentType;
    public int ContentId;
    public byte[] Data = [];

    public const int FileHeader = 0, CompressionHeader = 1, SliceHeader = 2, ExternalData = 4, CoreData = 5;

    /// <summary>Reads (and decompresses) the block at the cursor; CRC32 is checked for CRAM 3+.</summary>
    public static CramBlock Read(ref ByteCursor c, int major)
    {
        int start = c.Pos;
        var b = new CramBlock { Method = (CramBlockMethod)c.U8(), ContentType = c.U8(), ContentId = c.Itf8() };
        int compSize = c.Itf8();
        int rawSize = c.Itf8();
        var payload = c.Take(compSize);
        if (major >= 3)
        {
            uint crc = c.U32();
            if (Crc32.Compute(c.Data[start..(c.Pos - 4)]) != crc)
                throw new InvalidDataException($"CRAM block (content type {b.ContentType}, id {b.ContentId}): CRC32 mismatch.");
        }
        b.Data = rawSize == 0 ? [] : Decompress(b.Method, payload, rawSize);
        if (b.Data.Length != rawSize)
            throw new InvalidDataException($"CRAM block (content id {b.ContentId}): decompressed {b.Data.Length} bytes, header says {rawSize}.");
        return b;
    }

    private static byte[] Decompress(CramBlockMethod method, ReadOnlySpan<byte> data, int rawSize)
    {
        switch (method)
        {
            case CramBlockMethod.Raw:
                return data.ToArray();
            case CramBlockMethod.Gzip:
            {
                var o = new byte[rawSize];
                using var gz = new GZipStream(new MemoryStream(data.ToArray()), CompressionMode.Decompress);
                int n = gz.ReadAtLeast(o, rawSize, throwOnEndOfStream: false);
                return n == rawSize ? o : o[..n];
            }
            case CramBlockMethod.Bzip2: return Bzip2.Decompress(data, rawSize);
            case CramBlockMethod.Lzma: return Lzma.DecompressXz(data, rawSize);
            case CramBlockMethod.Rans4x8: return Rans4x8.Decode(data);
            case CramBlockMethod.RansNx16: return RansNx16.Decode(data, rawSize);
            case CramBlockMethod.Arith: return ArithDynamic.Decode(data, rawSize);
            case CramBlockMethod.FqzComp: return FqzComp.Decode(data);
            case CramBlockMethod.NameTokenizer: return NameTokenizer3.Decode(data);
            default: throw new NotSupportedException($"CRAM block compression method {(int)method} is not supported.");
        }
    }
}

/// <summary>A CRAM container header.</summary>
internal sealed class ContainerHeader
{
    public long Offset;       // file offset of the container
    public int HeaderSize;    // bytes of this header
    public int Length;        // bytes of blocks after the header
    public int RefId, Start, Span, NRecords, NBlocks;
    public long RecordCounter, Bases;
    public int[] Landmarks = [];

    public bool IsEof => NRecords == 0 && RefId == -1 && Start == 4542278;

    /// <summary>Reads a container header at the current stream position; null at end of file.</summary>
    public static ContainerHeader? Read(Stream s, int major)
    {
        long offset = s.Position;
        Span<byte> buf = stackalloc byte[4];
        int n = s.ReadAtLeast(buf, 4, throwOnEndOfStream: false);
        if (n == 0) return null;
        if (n < 4) throw new EndOfStreamException("Truncated CRAM container header.");
        var raw = new List<byte>(64);
        foreach (var b in buf) raw.Add(b);
        var h = new ContainerHeader { Offset = offset, Length = BinaryPrimitives.ReadInt32LittleEndian(buf) };
        h.RefId = Itf8(s, raw);
        h.Start = Itf8(s, raw);
        h.Span = Itf8(s, raw);
        h.NRecords = Itf8(s, raw);
        h.RecordCounter = major >= 3 ? Ltf8(s, raw) : Itf8(s, raw);
        h.Bases = major >= 3 ? Ltf8(s, raw) : Itf8(s, raw);
        h.NBlocks = Itf8(s, raw);
        int nl = Itf8(s, raw);
        h.Landmarks = new int[nl];
        for (int i = 0; i < nl; i++) h.Landmarks[i] = Itf8(s, raw);
        if (major >= 3)
        {
            if (s.ReadAtLeast(buf, 4, throwOnEndOfStream: false) != 4) throw new EndOfStreamException("Truncated CRAM container header.");
            uint crc = BinaryPrimitives.ReadUInt32LittleEndian(buf);
            if (Crc32.Compute(raw.ToArray()) != crc) throw new InvalidDataException($"CRAM container at offset {offset}: header CRC32 mismatch.");
            raw.AddRange(buf.ToArray());
        }
        h.HeaderSize = raw.Count;
        return h;
    }

    private static int Next(Stream s, List<byte> raw)
    {
        int b = s.ReadByte();
        if (b < 0) throw new EndOfStreamException("Truncated CRAM container header.");
        raw.Add((byte)b);
        return b;
    }

    private static int Itf8(Stream s, List<byte> raw)
    {
        int b = Next(s, raw);
        if (b < 0x80) return b;
        if (b < 0xc0) return ((b & 0x3f) << 8) | Next(s, raw);
        if (b < 0xe0) return ((b & 0x1f) << 16) | (Next(s, raw) << 8) | Next(s, raw);
        if (b < 0xf0) return ((b & 0x0f) << 24) | (Next(s, raw) << 16) | (Next(s, raw) << 8) | Next(s, raw);
        return (int)(((uint)(b & 0x0f) << 28) | ((uint)Next(s, raw) << 20) | ((uint)Next(s, raw) << 12) | ((uint)Next(s, raw) << 4) | ((uint)Next(s, raw) & 0x0f));
    }

    private static long Ltf8(Stream s, List<byte> raw)
    {
        int b = Next(s, raw);
        int extra;
        long v;
        if (b < 0x80) return b;
        if (b < 0xc0) { extra = 1; v = b & 0x3f; }
        else if (b < 0xe0) { extra = 2; v = b & 0x1f; }
        else if (b < 0xf0) { extra = 3; v = b & 0x0f; }
        else if (b < 0xf8) { extra = 4; v = b & 0x07; }
        else if (b < 0xfc) { extra = 5; v = b & 0x03; }
        else if (b < 0xfe) { extra = 6; v = b & 0x01; }
        else if (b < 0xff) { extra = 7; v = 0; }
        else { extra = 8; v = 0; }
        for (int i = 0; i < extra; i++) v = (v << 8) | (uint)Next(s, raw);
        return v;
    }
}

/// <summary>LTF8 from a byte cursor.</summary>
internal static class CramVarint
{
    public static long Ltf8(ref ByteCursor c)
    {
        int b = c.U8();
        int extra;
        long v;
        if (b < 0x80) return b;
        if (b < 0xc0) { extra = 1; v = b & 0x3f; }
        else if (b < 0xe0) { extra = 2; v = b & 0x1f; }
        else if (b < 0xf0) { extra = 3; v = b & 0x0f; }
        else if (b < 0xf8) { extra = 4; v = b & 0x07; }
        else if (b < 0xfc) { extra = 5; v = b & 0x03; }
        else if (b < 0xfe) { extra = 6; v = b & 0x01; }
        else if (b < 0xff) { extra = 7; v = 0; }
        else { extra = 8; v = 0; }
        for (int i = 0; i < extra; i++) v = (v << 8) | c.U8();
        return v;
    }
}

// ---------------------------------------------------------------------------------------------- slice data streams

/// <summary>MSB-first bit reader over the core data block.</summary>
internal sealed class CoreBits(byte[] data)
{
    private readonly byte[] _data = data;
    private long _bit;

    public int ReadBit()
    {
        long byteIdx = _bit >> 3;
        if (byteIdx >= _data.Length) throw new InvalidDataException("CRAM core block exhausted.");
        int b = (_data[byteIdx] >> (7 - (int)(_bit & 7))) & 1;
        _bit++;
        return b;
    }

    public int ReadBits(int n)
    {
        int v = 0;
        for (int i = 0; i < n; i++) v = (v << 1) | ReadBit();
        return v;
    }
}

/// <summary>An external data block being consumed.</summary>
internal sealed class ExternalStream(byte[] data)
{
    public readonly byte[] Data = data;
    public int Pos;

    public byte ReadByte() => Pos < Data.Length ? Data[Pos++] : throw new InvalidDataException("CRAM external block exhausted.");

    public int ReadItf8()
    {
        var c = new ByteCursor(Data) { Pos = Pos };
        int v = c.Itf8();
        Pos = c.Pos;
        return v;
    }

    public void ReadBytes(Span<byte> dst)
    {
        if (Pos + dst.Length > Data.Length) throw new InvalidDataException("CRAM external block exhausted.");
        Data.AsSpan(Pos, dst.Length).CopyTo(dst);
        Pos += dst.Length;
    }
}

/// <summary>The core and external streams of one slice.</summary>
internal sealed class SliceStreams
{
    public CoreBits Core = new([]);
    public readonly Dictionary<int, ExternalStream> External = [];
    public readonly HashSet<int> TokenizedBlocks = [];

    public ExternalStream Ext(int id) => External.TryGetValue(id, out var e) ? e
        : External[id] = new ExternalStream([]); // a series may point at a block that was omitted because it is empty
}

// ---------------------------------------------------------------------------------------------- encodings

/// <summary>A CRAM data-series encoding (codec + parameters).</summary>
internal abstract class CramEncoding
{
    public virtual int ReadInt(SliceStreams s) => throw new InvalidDataException($"{GetType().Name} cannot decode integers.");
    public virtual byte ReadByte(SliceStreams s) => (byte)ReadInt(s);

    /// <summary>Fills <paramref name="dst"/> with successive byte values.</summary>
    public virtual void ReadBytes(SliceStreams s, Span<byte> dst)
    {
        for (int i = 0; i < dst.Length; i++) dst[i] = ReadByte(s);
    }

    /// <summary>A variable-length byte array (BYTE_ARRAY_LEN / BYTE_ARRAY_STOP).</summary>
    public virtual byte[] ReadArray(SliceStreams s) => throw new InvalidDataException($"{GetType().Name} cannot decode byte arrays.");

    public static CramEncoding Parse(ref ByteCursor c)
    {
        int id = c.Itf8();
        int n = c.Itf8();
        var p = new ByteCursor(c.Take(n));
        switch (id)
        {
            case 0: return new NullEncoding();
            case 1: return new ExternalEncoding(p.Itf8());
            case 3:
            {
                int na = p.Itf8();
                var alphabet = new int[na];
                for (int i = 0; i < na; i++) alphabet[i] = p.Itf8();
                int nl = p.Itf8();
                var lengths = new int[nl];
                for (int i = 0; i < nl; i++) lengths[i] = p.Itf8();
                return new HuffmanEncoding(alphabet, lengths);
            }
            case 4:
            {
                var len = Parse(ref p);
                var val = Parse(ref p);
                return new ByteArrayLenEncoding(len, val);
            }
            case 5:
            {
                byte stop = p.U8();
                return new ByteArrayStopEncoding(stop, p.Itf8());
            }
            case 6: return new BetaEncoding(p.Itf8(), p.Itf8());
            case 7: return new SubexpEncoding(p.Itf8(), p.Itf8());
            case 9: return new GammaEncoding(p.Itf8());
            case 2: return new GolombEncoding(p.Itf8(), p.Itf8(), rice: false);
            case 8: return new GolombEncoding(p.Itf8(), p.Itf8(), rice: true);
            default: throw new NotSupportedException($"CRAM encoding (codec id {id}) is not supported.");
        }
    }
}

internal sealed class NullEncoding : CramEncoding
{
    public override int ReadInt(SliceStreams s) => 0;
    public override byte[] ReadArray(SliceStreams s) => [];
}

internal sealed class ExternalEncoding(int id) : CramEncoding
{
    public readonly int Id = id;
    public override int ReadInt(SliceStreams s) => s.Ext(Id).ReadItf8();
    public override byte ReadByte(SliceStreams s) => s.Ext(Id).ReadByte();
    public override void ReadBytes(SliceStreams s, Span<byte> dst) => s.Ext(Id).ReadBytes(dst);
}

internal sealed class HuffmanEncoding : CramEncoding
{
    private readonly int _constant;
    private readonly bool _isConstant;
    private readonly int[] _count = new int[33];
    private readonly int[] _first = new int[33];
    private readonly int[] _offset = new int[33];
    private readonly int[] _symbols;
    private readonly int _maxLen;

    public HuffmanEncoding(int[] alphabet, int[] lengths)
    {
        if (alphabet.Length != lengths.Length) throw new InvalidDataException("CRAM HUFFMAN: alphabet and lengths differ in size.");
        if (alphabet.Length == 0) throw new InvalidDataException("CRAM HUFFMAN: empty alphabet.");
        if (alphabet.Length == 1 && lengths[0] == 0) { _isConstant = true; _constant = alphabet[0]; _symbols = []; return; }
        // canonical codes: ordered by (length, symbol value)
        var order = Enumerable.Range(0, alphabet.Length).OrderBy(i => lengths[i]).ThenBy(i => alphabet[i]).ToArray();
        _symbols = order.Select(i => alphabet[i]).ToArray();
        foreach (int l in lengths)
        {
            if (l is < 0 or > 32) throw new InvalidDataException("CRAM HUFFMAN: invalid code length.");
            _count[l]++;
            _maxLen = Math.Max(_maxLen, l);
        }
        int code = 0, off = 0;
        for (int len = 1; len <= _maxLen; len++)
        {
            _first[len] = code;
            _offset[len] = off;
            off += _count[len];
            code = (code + _count[len]) << 1;
        }
        _offset[0] = 0;
    }

    public override int ReadInt(SliceStreams s)
    {
        if (_isConstant) return _constant;
        int code = 0;
        for (int len = 1; len <= _maxLen; len++)
        {
            code = (code << 1) | s.Core.ReadBit();
            int d = code - _first[len];
            if (d >= 0 && d < _count[len]) return _symbols[_offset[len] + d];
        }
        throw new InvalidDataException("CRAM HUFFMAN: invalid code.");
    }

    public override void ReadBytes(SliceStreams s, Span<byte> dst)
    {
        if (_isConstant) { dst.Fill((byte)_constant); return; }
        base.ReadBytes(s, dst);
    }
}

internal sealed class ByteArrayLenEncoding(CramEncoding length, CramEncoding values) : CramEncoding
{
    public override byte[] ReadArray(SliceStreams s)
    {
        int n = length.ReadInt(s);
        if (n < 0) throw new InvalidDataException("CRAM BYTE_ARRAY_LEN: negative length.");
        var a = new byte[n];
        values.ReadBytes(s, a);
        return a;
    }

    public override void ReadBytes(SliceStreams s, Span<byte> dst)
    {
        var a = ReadArray(s);
        a.AsSpan(0, Math.Min(a.Length, dst.Length)).CopyTo(dst);
    }
}

internal sealed class ByteArrayStopEncoding(byte stop, int id) : CramEncoding
{
    public override byte[] ReadArray(SliceStreams s)
    {
        var e = s.Ext(id);
        // the name tokeniser always terminates with NUL, whatever stop byte was declared (as htslib)
        byte st = s.TokenizedBlocks.Contains(id) ? (byte)0 : stop;
        int start = e.Pos;
        int end = Array.IndexOf(e.Data, st, start);
        if (end < 0) end = e.Data.Length;
        var a = e.Data.AsSpan(start, end - start).ToArray();
        e.Pos = Math.Min(e.Data.Length, end + 1);
        return a;
    }

    public override void ReadBytes(SliceStreams s, Span<byte> dst)
    {
        var a = ReadArray(s);
        a.AsSpan(0, Math.Min(a.Length, dst.Length)).CopyTo(dst);
    }
}

internal sealed class BetaEncoding(int offset, int bits) : CramEncoding
{
    public override int ReadInt(SliceStreams s) => (bits == 0 ? 0 : s.Core.ReadBits(bits)) - offset;
}

internal sealed class SubexpEncoding(int offset, int k) : CramEncoding
{
    public override int ReadInt(SliceStreams s)
    {
        int u = 0;
        while (s.Core.ReadBit() == 1) u++;
        int n;
        if (u == 0) n = s.Core.ReadBits(k);
        else
        {
            int b = u + k - 1;
            n = (1 << b) | s.Core.ReadBits(b);
        }
        return n - offset;
    }
}

internal sealed class GammaEncoding(int offset) : CramEncoding
{
    public override int ReadInt(SliceStreams s)
    {
        int nz = 0;
        while (s.Core.ReadBit() == 0) nz++;
        int v = 1;
        for (int i = 0; i < nz; i++) v = (v << 1) | s.Core.ReadBit();
        return v - offset;
    }
}

/// <summary>The deprecated GOLOMB (2) and GOLOMB_RICE (8) codecs, from the specification.</summary>
internal sealed class GolombEncoding(int offset, int m, bool rice) : CramEncoding
{
    public override int ReadInt(SliceStreams s)
    {
        int q = 0;
        while (s.Core.ReadBit() == 1) q++;
        if (rice)
        {
            // m is log2 of the divisor
            return (q << m | s.Core.ReadBits(m)) - offset;
        }
        int b = 0;
        while (1 << b < m) b++;
        int r = b > 0 ? s.Core.ReadBits(b - 1) : 0;
        int threshold = (1 << b) - m;
        if (r >= threshold) r = r * 2 + s.Core.ReadBit() - threshold;
        return q * m + r - offset;
    }
}
