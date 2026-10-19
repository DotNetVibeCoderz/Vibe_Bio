using System.Buffers.Binary;

namespace Biotech.Net.Genomics.Cram;

/// <summary>Sequential little-endian reader over a byte span, with the variable-length integers used by CRAM codecs.</summary>
internal ref struct ByteCursor(ReadOnlySpan<byte> data)
{
    public readonly ReadOnlySpan<byte> Data = data;
    public int Pos;

    public readonly int Remaining => Data.Length - Pos;
    public readonly bool End => Pos >= Data.Length;

    public byte U8()
    {
        if (Pos >= Data.Length) throw new InvalidDataException("Unexpected end of CRAM codec data.");
        return Data[Pos++];
    }

    /// <summary>A byte, or 0 past the end (entropy decoders may peek beyond the last renormalisation byte).</summary>
    public byte U8OrZero() => Pos < Data.Length ? Data[Pos++] : (byte)0;

    public ushort U16()
    {
        var v = BinaryPrimitives.ReadUInt16LittleEndian(Take(2));
        return v;
    }

    public uint U32() => BinaryPrimitives.ReadUInt32LittleEndian(Take(4));

    public ReadOnlySpan<byte> Take(int n)
    {
        if (n < 0 || Pos + n > Data.Length) throw new InvalidDataException("Unexpected end of CRAM codec data.");
        var s = Data.Slice(Pos, n);
        Pos += n;
        return s;
    }

    /// <summary>7 bits per byte, most significant group first, top bit = continuation (htscodecs <c>var_get_u32</c>).</summary>
    public uint Uint7()
    {
        uint v = 0;
        byte c;
        int n = 0;
        do
        {
            c = U8();
            v = (v << 7) | (uint)(c & 0x7f);
        } while ((c & 0x80) != 0 && ++n < 5);
        return v;
    }

    /// <summary>CRAM ITF8 integer.</summary>
    public int Itf8()
    {
        int b = U8();
        if (b < 0x80) return b;
        if (b < 0xc0) return ((b & 0x3f) << 8) | U8();
        if (b < 0xe0) return ((b & 0x1f) << 16) | (U8() << 8) | U8();
        if (b < 0xf0) return ((b & 0x0f) << 24) | (U8() << 16) | (U8() << 8) | U8();
        return (int)(((uint)(b & 0x0f) << 28) | ((uint)U8() << 20) | ((uint)U8() << 12) | ((uint)U8() << 4) | ((uint)U8() & 0x0f));
    }
}

/// <summary>rANS 4x8 (CRAM 3.0 block method 4): 4 interleaved states, 8-bit renormalisation, order 0 or 1.</summary>
internal static class Rans4x8
{
    private const int TotFreqBits = 12;
    private const uint Lower = 1u << 23;

    public static byte[] Decode(ReadOnlySpan<byte> input)
    {
        var c = new ByteCursor(input);
        int order = c.U8();
        c.U32(); // compressed size
        int nOut = checked((int)c.U32());
        var output = new byte[nOut];
        if (nOut == 0) return output;
        if (order == 0) DecodeO0(ref c, output);
        else if (order == 1) DecodeO1(ref c, output);
        else throw new InvalidDataException($"rANS 4x8: invalid order {order}.");
        return output;
    }

    private static void ReadFreqs0(ref ByteCursor c, int[] f, int[] cum, byte[] lut)
    {
        Array.Clear(f);
        int s = c.U8(), last = s, rle = 0;
        do
        {
            f[s] = c.Itf8();
            if (rle > 0) { rle--; s++; }
            else
            {
                s = c.U8();
                if (s == last + 1) rle = c.U8();
            }
            last = s;
        } while (s != 0);
        int x = 0;
        for (int i = 0; i < 256; i++)
        {
            cum[i] = x;
            if (f[i] < 0 || x + f[i] > (1 << TotFreqBits)) throw new InvalidDataException("rANS 4x8: frequencies exceed 4096.");
            for (int k = 0; k < f[i]; k++) lut[x + k] = (byte)i;
            x += f[i];
        }
        cum[256] = x;
    }

    private static uint Renorm(uint r, ref ByteCursor c)
    {
        while (r < Lower) r = (r << 8) | c.U8OrZero();
        return r;
    }

    private static void DecodeO0(ref ByteCursor c, byte[] output)
    {
        var f = new int[256];
        var cum = new int[257];
        var lut = new byte[1 << TotFreqBits];
        ReadFreqs0(ref c, f, cum, lut);
        Span<uint> r = stackalloc uint[4];
        for (int j = 0; j < 4; j++) r[j] = c.U32();
        const uint mask = (1u << TotFreqBits) - 1;
        for (int i = 0; i < output.Length; i++)
        {
            int j = i & 3;
            uint x = r[j];
            byte s = lut[x & mask];
            output[i] = s;
            x = (uint)f[s] * (x >> TotFreqBits) + (x & mask) - (uint)cum[s];
            r[j] = Renorm(x, ref c);
        }
    }

    private static void DecodeO1(ref ByteCursor c, byte[] output)
    {
        var f = new int[256][];
        var cum = new int[256][];
        var lut = new byte[256][];
        int ctx = c.U8(), last = ctx, rle = 0;
        do
        {
            f[ctx] = new int[256];
            cum[ctx] = new int[257];
            lut[ctx] = new byte[1 << TotFreqBits];
            ReadFreqs0(ref c, f[ctx], cum[ctx], lut[ctx]);
            if (rle > 0) { rle--; ctx++; }
            else
            {
                ctx = c.U8();
                if (ctx == last + 1) rle = c.U8();
            }
            last = ctx;
        } while (ctx != 0);

        Span<uint> r = stackalloc uint[4];
        Span<int> l = stackalloc int[4];
        for (int j = 0; j < 4; j++) r[j] = c.U32();
        const uint mask = (1u << TotFreqBits) - 1;
        int n = output.Length, isz = n >> 2;
        for (int i = 0; i < isz; i++)
            for (int j = 0; j < 4; j++)
            {
                output[i + j * isz] = Step(ref r[j], ref l[j], ref c);
            }
        for (int i = 4 * isz; i < n; i++) output[i] = Step(ref r[3], ref l[3], ref c);

        byte Step(ref uint state, ref int prev, ref ByteCursor cur)
        {
            var lt = lut[prev] ?? throw new InvalidDataException("rANS 4x8: missing order-1 context.");
            uint x = state;
            byte s = lt[x & mask];
            x = (uint)f[prev][s] * (x >> TotFreqBits) + (x & mask) - (uint)cum[prev][s];
            state = Renorm(x, ref cur);
            prev = s;
            return s;
        }
    }
}

/// <summary>Byte-stream transforms shared by rANS Nx16 and the adaptive arithmetic coder.</summary>
internal static class CodecTransforms
{
    public const int Order = 1, N32 = 4, Stripe = 8, NoSize = 16, Cat = 32, Rle = 64, Pack = 128;
    public const int ArithExt = 4;

    public static byte[] Unpack(ReadOnlySpan<byte> data, byte[] map, int nsym, int len)
    {
        var o = new byte[len];
        if (nsym <= 1) { Array.Fill(o, map[0]); return o; }
        if (nsym > 16)
        {
            // htscodecs: more than 16 symbols means the data was not packed at all
            data[..Math.Min(len, data.Length)].CopyTo(o);
            return o;
        }
        int bits = nsym <= 2 ? 1 : nsym <= 4 ? 2 : nsym <= 16 ? 4 : throw new InvalidDataException("Pack: more than 16 symbols.");
        int per = 8 / bits, mask = (1 << bits) - 1, j = 0, v = 0;
        for (int i = 0; i < len; i++)
        {
            if (i % per == 0) v = j < data.Length ? data[j++] : throw new InvalidDataException("Pack: data too short.");
            o[i] = map[v & mask];
            v >>= bits;
        }
        return o;
    }

    /// <summary>Reads pack meta-data: returns the symbol map and count; <paramref name="packedLength"/> is the packed stream size.</summary>
    public static (byte[] Map, int NSym) ReadPackMeta(ref ByteCursor c, out int packedLength)
    {
        int n = c.U8();
        if (n == 0) n = 256;
        var map = new byte[Math.Max(n, 16)];
        if (n <= 16) for (int i = 0; i < n; i++) map[i] = c.U8();
        packedLength = (int)c.Uint7();
        return (map, n);
    }

    public static byte[] Unstripe(byte[][] parts, int len)
    {
        int n = parts.Length;
        var o = new byte[len];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < parts[j].Length; i++) o[i * n + j] = parts[j][i];
        return o;
    }

    public static int StripeLength(int len, int n, int j) => len / n + (len % n > j ? 1 : 0);
}

/// <summary>rANS Nx16 (CRAM 3.1 block method 5) with its stripe, RLE, pack and cat transforms.</summary>
internal static class RansNx16
{
    public static byte[] Decode(ReadOnlySpan<byte> input, int expectedLength = -1)
    {
        var c = new ByteCursor(input);
        return Decode(ref c, expectedLength);
    }

    public static byte[] Decode(ref ByteCursor c, int len)
    {
        int flags = c.U8();
        if ((flags & CodecTransforms.NoSize) == 0) len = (int)c.Uint7();
        if (len < 0) throw new InvalidDataException("rANS Nx16: unknown output size.");
        if ((flags & CodecTransforms.Stripe) != 0)
        {
            int n = c.U8();
            var clen = new int[n];
            for (int j = 0; j < n; j++) clen[j] = (int)c.Uint7();
            var parts = new byte[n][];
            for (int j = 0; j < n; j++)
            {
                var sub = new ByteCursor(c.Take(clen[j]));
                parts[j] = Decode(ref sub, CodecTransforms.StripeLength(len, n, j));
            }
            return CodecTransforms.Unstripe(parts, len);
        }
        int nStates = (flags & CodecTransforms.N32) != 0 ? 32 : 4;
        int packLen = len, rleLen = len;
        byte[]? map = null;
        int nsym = 0;
        if ((flags & CodecTransforms.Pack) != 0)
        {
            (map, nsym) = CodecTransforms.ReadPackMeta(ref c, out len);
        }
        bool[]? rleSyms = null;
        byte[]? runs = null;
        int runPos = 0;
        if ((flags & CodecTransforms.Rle) != 0)
        {
            rleLen = len;
            uint metaLen = c.Uint7();
            len = (int)c.Uint7();
            byte[] meta;
            if ((metaLen & 1) != 0) meta = c.Take((int)(metaLen / 2)).ToArray();
            else
            {
                int compLen = (int)c.Uint7();
                var mc = new ByteCursor(c.Take(compLen));
                meta = DecodeO0(ref mc, (int)(metaLen / 2), nStates);
            }
            int ns = meta[0] == 0 ? 256 : meta[0];
            rleSyms = new bool[256];
            for (int i = 0; i < ns; i++) rleSyms[meta[1 + i]] = true;
            runs = meta;
            runPos = 1 + ns;
        }
        byte[] data;
        if (c.End && len > 0 && (flags & CodecTransforms.Cat) == 0) data = new byte[0];
        else if ((flags & CodecTransforms.Cat) != 0) data = c.Take(len).ToArray();
        else if ((flags & CodecTransforms.Order) != 0) data = DecodeO1(ref c, len, nStates);
        else data = DecodeO0(ref c, len, nStates);
        if (rleSyms is not null)
        {
            var o = new byte[rleLen];
            var rc = new ByteCursor(runs!) { Pos = runPos };
            int j = 0;
            foreach (byte b in data)
            {
                if (rleSyms[b])
                {
                    int run = (int)rc.Uint7();
                    if (j + run + 1 > o.Length) throw new InvalidDataException("rANS Nx16: RLE overflows the output.");
                    o.AsSpan(j, run + 1).Fill(b);
                    j += run + 1;
                }
                else
                {
                    if (j >= o.Length) throw new InvalidDataException("rANS Nx16: RLE overflows the output.");
                    o[j++] = b;
                }
            }
            data = j == o.Length ? o : o[..j];
        }
        if (map is not null) data = CodecTransforms.Unpack(data, map, nsym, packLen);
        return data;
    }

    private static int[] ReadAlphabet(ref ByteCursor c)
    {
        var a = new List<int>();
        int s = c.U8(), last = s, rle = 0;
        do
        {
            a.Add(s);
            if (rle > 0) { rle--; s++; }
            else
            {
                s = c.U8();
                if (s == last + 1) rle = c.U8();
            }
            last = s;
        } while (s != 0);
        return a.ToArray();
    }

    private static void Normalise(int[] f, int bits)
    {
        long tot = 0;
        foreach (int v in f) tot += v;
        if (tot == 0 || tot == 1L << bits) return;
        int shift = 0;
        while (tot < 1L << bits) { tot *= 2; shift++; }
        for (int i = 0; i < 256; i++) f[i] <<= shift;
    }

    private static void BuildTables(int[] f, int bits, int[] cum, byte[] lut)
    {
        int x = 0;
        for (int i = 0; i < 256; i++)
        {
            cum[i] = x;
            if (x + f[i] > 1 << bits) throw new InvalidDataException("rANS Nx16: frequencies exceed the table size.");
            for (int k = 0; k < f[i]; k++) lut[x + k] = (byte)i;
            x += f[i];
        }
    }

    private static uint Renorm(uint r, ref ByteCursor c)
    {
        if (r < 1u << 15)
        {
            uint lo = c.U8OrZero();
            uint hi = c.U8OrZero();
            r = (r << 16) | (hi << 8) | lo;
        }
        return r;
    }

    internal static byte[] DecodeO0(ref ByteCursor c, int len, int n)
    {
        var f = new int[256];
        foreach (int s in ReadAlphabet(ref c)) f[s] = (int)c.Uint7();
        Normalise(f, 12);
        var cum = new int[256];
        var lut = new byte[4096];
        BuildTables(f, 12, cum, lut);
        var r = new uint[n];
        for (int j = 0; j < n; j++) r[j] = c.U32();
        var o = new byte[len];
        for (int i = 0; i < len; i++)
        {
            int j = i % n;
            uint x = r[j];
            byte s = lut[x & 0xfff];
            o[i] = s;
            x = (uint)f[s] * (x >> 12) + (x & 0xfff) - (uint)cum[s];
            r[j] = Renorm(x, ref c);
        }
        return o;
    }

    private static byte[] DecodeO1(ref ByteCursor c, int len, int n)
    {
        int comp = c.U8();
        int bits = comp >> 4;
        byte[]? freqData = null;
        if ((comp & 1) != 0)
        {
            int uSize = (int)c.Uint7();
            int cSize = (int)c.Uint7();
            var fc = new ByteCursor(c.Take(cSize));
            freqData = DecodeO0(ref fc, uSize, 4);
        }
        var f = new int[256][];
        var cum = new int[256][];
        var lut = new byte[256][];
        if (freqData is not null)
        {
            var src = new ByteCursor(freqData);
            ReadO1Freqs(ref src, bits, f, cum, lut);
        }
        else ReadO1Freqs(ref c, bits, f, cum, lut);
        var r = new uint[n];
        var l = new int[n];
        for (int j = 0; j < n; j++) r[j] = c.U32();
        var o = new byte[len];
        int isz = len / n;
        uint mask = (1u << bits) - 1;
        for (int i = 0; i < isz; i++)
            for (int j = 0; j < n; j++)
                o[i + j * isz] = Step(j, ref c);
        for (int i = n * isz; i < len; i++) o[i] = Step(n - 1, ref c);
        return o;

        byte Step(int j, ref ByteCursor cur)
        {
            int prev = l[j];
            var lt = lut[prev] ?? throw new InvalidDataException("rANS Nx16: missing order-1 context.");
            uint x = r[j];
            byte s = lt[x & mask];
            x = (uint)f[prev][s] * (x >> bits) + (x & mask) - (uint)cum[prev][s];
            r[j] = Renorm(x, ref cur);
            l[j] = s;
            return s;
        }
    }

    private static void ReadO1Freqs(ref ByteCursor src, int bits, int[][] f, int[][] cum, byte[][] lut)
    {
        if (bits is not (10 or 12)) throw new InvalidDataException($"rANS Nx16: invalid order-1 frequency bits {bits}.");
        var a = ReadAlphabet(ref src);
        foreach (int i in a)
        {
            var fi = new int[256];
            int run = 0;
            foreach (int j in a)
            {
                if (run > 0) { run--; continue; }
                fi[j] = (int)src.Uint7();
                if (fi[j] == 0) run = src.U8();
            }
            Normalise(fi, bits);
            f[i] = fi;
            cum[i] = new int[256];
            lut[i] = new byte[1 << bits];
            BuildTables(fi, bits, cum[i], lut[i]);
        }
    }
}

/// <summary>The byte-wise range coder of htscodecs (Eugene Shelwien's carry-less coder, decode side).</summary>
internal ref struct RangeDecoder
{
    private ByteCursor _in;
    private uint _range, _code;

    public RangeDecoder(ReadOnlySpan<byte> data)
    {
        _in = new ByteCursor(data);
        _range = 0xFFFFFFFF;
        _code = 0;
        for (int i = 0; i < 5; i++) _code = (_code << 8) | _in.U8OrZero();
    }

    public uint GetFreq(uint totFreq) => totFreq != 0 && _range >= totFreq ? _code / (_range /= totFreq) : 0;

    public void Decode(uint cumFreq, uint freq)
    {
        _code -= cumFreq * _range;
        _range *= freq;
        while (_range < 1u << 24)
        {
            _code = (_code << 8) | _in.U8OrZero();
            _range <<= 8;
        }
    }
}

/// <summary>Adaptive frequency model (htscodecs <c>c_simple_model.h</c>): +16 per hit, halving above 65519, bubble towards the front.</summary>
internal sealed class SimpleModel
{
    private const uint MaxFreq = (1 << 16) - 17;
    private const int Step = 16;
    private readonly ushort[] _freq;
    private readonly ushort[] _sym;
    private uint _total;

    public SimpleModel(int maxSym, int capacity = -1)
    {
        int n = Math.Max(maxSym, capacity);
        _freq = new ushort[n + 1];
        _sym = new ushort[n + 1];
        for (int i = 0; i < n; i++) { _sym[i] = (ushort)i; _freq[i] = (ushort)(i < maxSym ? 1 : 0); }
        _total = (uint)maxSym;
    }

    public int Decode(ref RangeDecoder rc)
    {
        uint f = rc.GetFreq(_total);
        if (f > MaxFreq) return 0;
        uint acc = 0;
        int x = 0;
        while (x < _freq.Length - 1 && acc + _freq[x] <= f) acc += _freq[x++];
        if (_freq[x] == 0) throw new InvalidDataException("Arithmetic decoder: corrupt data (symbol out of range).");
        rc.Decode(acc, _freq[x]);
        _freq[x] += Step;
        _total += Step;
        if (_total > MaxFreq)
        {
            _total = 0;
            for (int i = 0; i < _freq.Length && _freq[i] != 0; i++)
            {
                _freq[i] -= (ushort)(_freq[i] >> 1);
                _total += _freq[i];
            }
        }
        int sym = _sym[x];
        if (x > 0 && _freq[x] > _freq[x - 1])
        {
            (_freq[x], _freq[x - 1]) = (_freq[x - 1], _freq[x]);
            (_sym[x], _sym[x - 1]) = (_sym[x - 1], _sym[x]);
        }
        return sym;
    }
}

/// <summary>The CRAM 3.1 adaptive arithmetic coder (block method 6), order 0/1, with RLE, pack, stripe and cat.</summary>
internal static class ArithDynamic
{
    public static byte[] Decode(ReadOnlySpan<byte> input, int expectedLength = -1)
    {
        var c = new ByteCursor(input);
        return Decode(ref c, expectedLength);
    }

    public static byte[] Decode(ref ByteCursor c, int len)
    {
        int flags = c.U8();
        if ((flags & CodecTransforms.NoSize) == 0) len = (int)c.Uint7();
        if ((flags & CodecTransforms.Stripe) != 0)
        {
            int n = c.U8();
            var clen = new int[n];
            for (int j = 0; j < n; j++) clen[j] = (int)c.Uint7();
            var parts = new byte[n][];
            for (int j = 0; j < n; j++)
            {
                var sub = new ByteCursor(c.Take(clen[j]));
                parts[j] = Decode(ref sub, CodecTransforms.StripeLength(len, n, j));
            }
            return CodecTransforms.Unstripe(parts, len);
        }
        int packLen = len;
        byte[]? map = null;
        int nsym = 0;
        if ((flags & CodecTransforms.Pack) != 0) (map, nsym) = CodecTransforms.ReadPackMeta(ref c, out len);
        byte[] data;
        var rest = c.Data[c.Pos..];
        if ((flags & CodecTransforms.Cat) != 0) data = c.Take(len).ToArray();
        else if ((flags & CodecTransforms.ArithExt) != 0) data = Bzip2.Decompress(rest, len);
        else if ((flags & CodecTransforms.Rle) != 0) data = DecodeRle(rest, len, (flags & CodecTransforms.Order) != 0);
        else data = DecodeOrder(rest, len, (flags & CodecTransforms.Order) != 0);
        if (map is not null) data = CodecTransforms.Unpack(data, map, nsym, packLen);
        return data;
    }

    private static byte[] DecodeOrder(ReadOnlySpan<byte> data, int len, bool order1)
    {
        var o = new byte[len];
        if (len == 0) return o;
        int maxSym = data[0] == 0 ? 256 : data[0];
        var rc = new RangeDecoder(data[1..]);
        if (!order1)
        {
            var m = new SimpleModel(maxSym, 256);
            for (int i = 0; i < len; i++) o[i] = (byte)m.Decode(ref rc);
        }
        else
        {
            var models = new SimpleModel?[256];
            int last = 0;
            for (int i = 0; i < len; i++)
            {
                o[i] = (byte)(models[last] ??= new SimpleModel(maxSym, 256)).Decode(ref rc);
                last = o[i];
            }
        }
        return o;
    }

    private static byte[] DecodeRle(ReadOnlySpan<byte> data, int len, bool order1)
    {
        var o = new byte[len];
        if (len == 0) return o;
        int maxSym = data[0] == 0 ? 256 : data[0];
        var rc = new RangeDecoder(data[1..]);
        var lit = new SimpleModel?[order1 ? 256 : 1];
        var runModels = new SimpleModel[258];
        for (int i = 0; i < 258; i++) runModels[i] = new SimpleModel(4);
        int last = 0, k = 0;
        while (k < len)
        {
            var m = lit[order1 ? last : 0] ??= new SimpleModel(maxSym, 256);
            int sym = m.Decode(ref rc);
            last = sym;
            int part = runModels[sym].Decode(ref rc);
            int run = part, rctx = 256;
            while (part == 3 && run < len)
            {
                part = runModels[rctx].Decode(ref rc);
                rctx = 257;
                run += part;
            }
            int n = Math.Min(run + 1, len - k); // htscodecs clips runs at the output size
            o.AsSpan(k, n).Fill((byte)sym);
            k += n;
        }
        return o;
    }
}
