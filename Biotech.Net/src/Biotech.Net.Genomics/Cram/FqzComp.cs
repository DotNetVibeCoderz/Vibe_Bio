namespace Biotech.Net.Genomics.Cram;

/// <summary>
/// The fqzcomp quality codec (CRAM 3.1 block method 7), decode side. Follows the htscodecs reference decoder
/// (<c>fqzcomp_qual.c</c>), which differs from the published pseudocode in a few details: later contexts do not
/// re-add the parameter's starting context, the <c>DO_LEN</c> flag means "fixed length", and the selector,
/// context update and quality map use the first parameter block.
/// </summary>
internal static class FqzComp
{
    private const int GMultiParam = 1, GHaveStab = 2, GDoRev = 4;
    private const int PDoDedup = 2, PDoLen = 4, PDoSel = 8, PHaveQmap = 16, PHavePtab = 32, PHaveDtab = 64, PHaveQtab = 128;

    private sealed class Param
    {
        public int Context, Flags, MaxSym, QBits, QShift, QLoc, SLoc, PLoc, DLoc;
        public uint QMask;
        public readonly int[] QMap = new int[256];
        public readonly int[] QTab = new int[256];
        public readonly int[] PTab = new int[1024];
        public readonly int[] DTab = new int[256];
        public bool FixedLen => (Flags & PDoLen) != 0;
        public bool DoDedup => (Flags & PDoDedup) != 0;
        public bool DoSel => (Flags & PDoSel) != 0;
    }

    public static byte[] Decode(ReadOnlySpan<byte> input)
    {
        var c = new ByteCursor(input);
        int len = (int)c.Uint7();
        int vers = c.U8();
        if (vers != 5) throw new InvalidDataException($"fqzcomp: unsupported version {vers}.");
        int gflags = c.U8();
        int nparam = (gflags & GMultiParam) != 0 ? c.U8() : 1;
        if (nparam <= 0) throw new InvalidDataException("fqzcomp: no parameter blocks.");
        int maxSel = nparam > 1 ? nparam : 0;
        var stab = new int[256];
        if ((gflags & GHaveStab) != 0)
        {
            maxSel = c.U8();
            ReadArray(ref c, stab, 256);
        }
        else
        {
            for (int i = 0; i < 256; i++) stab[i] = i < nparam ? i : nparam - 1;
        }
        var ps = new Param[nparam];
        int maxSym = 0;
        for (int i = 0; i < nparam; i++)
        {
            ps[i] = ReadParam(ref c);
            if (ps[i].DoSel && maxSel == 0) throw new InvalidDataException("fqzcomp: selector used without selector values.");
            maxSym = Math.Max(maxSym, ps[i].MaxSym);
        }
        foreach (var p in ps)
        {
            for (int j = 0; j < 1024; j++) p.PTab[j] <<= p.PLoc;
            for (int j = 0; j < 256; j++) p.DTab[j] <<= p.DLoc;
        }

        var qual = new SimpleModel?[1 << 16];
        var lenModels = new SimpleModel[4];
        for (int i = 0; i < 4; i++) lenModels[i] = new SimpleModel(256);
        var revModel = new SimpleModel(2);
        var dupModel = new SimpleModel(2);
        var selModel = maxSel > 0 ? new SimpleModel(maxSel + 1, 256) : null;
        var rc = new RangeDecoder(c.Data[c.Pos..]);

        var output = new byte[len];
        var revs = new List<(bool Rev, int Len)>();
        var p0 = ps[0];
        int pos = 0, delta = 0, prevq = 0, sel = 0, lastLen = 0;
        uint qctx = 0, last = 0;
        bool firstLen = true;
        int i0 = 0;
        while (i0 < len)
        {
            if (pos == 0)
            {
                // new record
                sel = p0.DoSel ? selModel!.Decode(ref rc) : 0;
                int x = (gflags & GHaveStab) != 0 ? stab[Math.Min(255, sel)] : sel;
                if (x >= nparam) throw new InvalidDataException("fqzcomp: parameter selector out of range.");
                var pm = ps[x];
                int rlen = lastLen;
                if (!pm.FixedLen || firstLen)
                {
                    rlen = lenModels[0].Decode(ref rc);
                    rlen |= lenModels[1].Decode(ref rc) << 8;
                    rlen |= lenModels[2].Decode(ref rc) << 16;
                    rlen |= lenModels[3].Decode(ref rc) << 24;
                    firstLen = false;
                    lastLen = rlen;
                }
                if (rlen <= 0 || rlen > len - i0) throw new InvalidDataException("fqzcomp: invalid record length.");
                if ((gflags & GDoRev) != 0) revs.Add((revModel.Decode(ref rc) != 0, rlen));
                if (pm.DoDedup && dupModel.Decode(ref rc) != 0)
                {
                    if (rlen > i0) throw new InvalidDataException("fqzcomp: duplicate of a missing record.");
                    output.AsSpan(i0 - rlen, rlen).CopyTo(output.AsSpan(i0));
                    i0 += rlen;
                    continue;
                }
                pos = rlen;
                delta = 0;
                prevq = 0;
                qctx = 0;
                last = (uint)pm.Context;
            }
            do
            {
                var m = qual[last] ??= new SimpleModel(maxSym + 1, 256);
                int q = m.Decode(ref rc);
                // context update (fqz_update_ctx) with the first parameter block, as htscodecs does
                qctx = (qctx << p0.QShift) + (uint)p0.QTab[q];
                uint ctx = (qctx & p0.QMask) << p0.QLoc;
                ctx += (uint)p0.PTab[Math.Min(1023, pos)];
                ctx += (uint)p0.DTab[Math.Min(255, delta)];
                ctx += (uint)(sel << p0.SLoc);
                delta += prevq != q ? 1 : 0;
                prevq = q;
                pos--;
                last = ctx & 0xFFFF;
                output[i0++] = (byte)p0.QMap[q];
            } while (pos != 0 && i0 < len);
        }
        if ((gflags & GDoRev) != 0)
        {
            int at = 0;
            foreach (var (rev, rl) in revs)
            {
                if (rev) Array.Reverse(output, at, rl);
                at += rl;
            }
        }
        return output;
    }

    private static Param ReadParam(ref ByteCursor c)
    {
        var p = new Param
        {
            Context = c.U16(),
            Flags = c.U8(),
            MaxSym = c.U8(),
        };
        int x = c.U8();
        p.QBits = x >> 4; p.QShift = x & 15;
        x = c.U8();
        p.QLoc = x >> 4; p.SLoc = x & 15;
        x = c.U8();
        p.PLoc = x >> 4; p.DLoc = x & 15;
        p.QMask = (1u << p.QBits) - 1;
        if ((p.Flags & PHaveQmap) != 0)
            for (int i = 0; i < p.MaxSym; i++) p.QMap[i] = c.U8();
        else
            for (int i = 0; i < 256; i++) p.QMap[i] = i;
        if (p.QBits > 0)
        {
            if ((p.Flags & PHaveQtab) != 0) ReadArray(ref c, p.QTab, 256);
            else for (int i = 0; i < 256; i++) p.QTab[i] = i;
        }
        if ((p.Flags & PHavePtab) != 0) ReadArray(ref c, p.PTab, 1024);
        if ((p.Flags & PHaveDtab) != 0) ReadArray(ref c, p.DTab, 256);
        return p;
    }

    /// <summary>Two-level run-length decoded monotonic table (htscodecs <c>read_array</c>).</summary>
    private static void ReadArray(ref ByteCursor c, int[] array, int size)
    {
        var r = new List<int>(size);
        int z = 0, last = -1;
        while (z < size)
        {
            int run = c.U8();
            r.Add(run);
            z += run;
            if (run == last)
            {
                int copy = c.U8();
                z += run * copy;
                for (int k = 0; k < copy && z <= size; k++) r.Add(run);
            }
            last = run;
        }
        int j = 0, zi = 0;
        for (int i = 0; j < size; i++)
        {
            int runLen = 0, part;
            do
            {
                if (zi >= r.Count) throw new InvalidDataException("fqzcomp: truncated lookup table.");
                part = r[zi++];
                runLen += part;
            } while (part == 255);
            while (runLen-- > 0 && j < size) array[j++] = i;
        }
    }
}

/// <summary>The CRAM 3.1 read-name tokeniser (block method 8, "tok3"), decode side.</summary>
internal static class NameTokenizer3
{
    private const int TType = 0, TAlpha = 1, TChar = 2, TDigits0 = 3, TDzLen = 4, TDup = 5, TDiff = 6, TDigits = 7,
        TDelta = 8, TDelta0 = 9, TMatch = 10, TNop = 11, TEnd = 12;

    private sealed class Stream(byte[] data)
    {
        public readonly byte[] Data = data;
        public int Pos;
    }

    private struct Token
    {
        public int Type;
        public uint Int;
        public int Str; // offset of an ALPHA token in its name, or DIGITS0 width
        public int Len; // ALPHA length
    }

    /// <summary>Decodes the names; each is followed by a NUL byte, as htscodecs returns them.</summary>
    public static byte[] Decode(ReadOnlySpan<byte> input)
    {
        if (input.Length < 9) throw new InvalidDataException("tok3: data too short.");
        var c = new ByteCursor(input);
        int ulen = (int)c.U32();
        int nreads = (int)c.U32();
        bool arith = c.U8() != 0;
        var streams = new Dictionary<int, Stream>();
        int tnum = -1;
        while (!c.End)
        {
            int ttype = c.U8();
            if ((ttype & 128) != 0)
            {
                tnum++;
                if (tnum >= 128) throw new InvalidDataException("tok3: too many token positions.");
                for (int k = 0; k < 16; k++) streams.Remove(tnum << 4 | k);
                if ((ttype & 15) != 0)
                {
                    // a type stream of "this type, then MATCH" was elided
                    var buf = new byte[Math.Max(1, nreads)];
                    buf[0] = (byte)(ttype & 15);
                    buf.AsSpan(1).Fill(TMatch);
                    streams[tnum << 4] = new Stream(buf);
                }
            }
            if (tnum < 0) throw new InvalidDataException("tok3: stream before the first token position.");
            int id = tnum << 4 | (ttype & 15);
            if ((ttype & 64) != 0)
            {
                int j = c.U8() << 4;
                j += c.U8();
                if (j >= id || !streams.TryGetValue(j, out var src)) throw new InvalidDataException("tok3: duplicate of a missing stream.");
                streams[id] = new Stream(src.Data);
                continue;
            }
            int clen = (int)c.Uint7();
            var cdata = c.Take(clen);
            streams[id] = new Stream(arith ? ArithDynamic.Decode(cdata) : RansNx16.Decode(cdata));
        }

        var output = new byte[ulen];
        int outLen = 0;
        var names = new (int Start, int Len, Token[] Tokens, int NTok)[nreads];
        for (int n = 0; n < nreads; n++)
        {
            int t0 = ReadType(streams, 0);
            if (t0 < 0) break;
            uint dist = ReadU32(streams, 0, t0);
            if (dist > n) throw new InvalidDataException("tok3: name refers beyond the first name.");
            int pnum = n - (int)dist;
            int start = outLen;
            if (t0 == TDup)
            {
                if (pnum == n) throw new InvalidDataException("tok3: name duplicates itself.");
                var pn = names[pnum];
                Append(ref output, ref outLen, output.AsSpan(pn.Start, pn.Len).ToArray());
                names[n] = (start, pn.Len, pn.Tokens, pn.NTok);
                continue;
            }
            var toks = new Token[128];
            int ntok;
            bool ended = false;
            for (ntok = 1; ntok < 128; ntok++)
            {
                int tok = ReadType(streams, ntok);
                var prev = names[pnum];
                switch (tok)
                {
                    case TChar:
                    {
                        byte ch = ReadU8(streams, ntok, TChar);
                        AppendByte(ref output, ref outLen, ch);
                        toks[ntok] = new Token { Type = TChar, Int = ch };
                        break;
                    }
                    case TAlpha:
                    {
                        var s = Get(streams, ntok, TAlpha);
                        int st = s.Pos;
                        while (s.Pos < s.Data.Length && s.Data[s.Pos] != 0) s.Pos++;
                        int l = s.Pos - st;
                        if (s.Pos < s.Data.Length) s.Pos++;
                        int off = outLen - start;
                        Append(ref output, ref outLen, s.Data.AsSpan(st, l).ToArray());
                        toks[ntok] = new Token { Type = TAlpha, Str = off, Len = l };
                        break;
                    }
                    case TDigits0:
                    {
                        uint w = ReadU8(streams, ntok, TDzLen);
                        uint v = ReadU32(streams, ntok, TDigits0);
                        AppendNumber(ref output, ref outLen, v, (int)w);
                        toks[ntok] = new Token { Type = TDigits0, Int = v, Str = (int)w };
                        break;
                    }
                    case TDelta0:
                    {
                        if (ntok >= prev.NTok) throw new InvalidDataException("tok3: DELTA0 without a previous token.");
                        uint v = ReadU8(streams, ntok, TDelta0) + prev.Tokens[ntok].Int;
                        int w = prev.Tokens[ntok].Str;
                        AppendNumber(ref output, ref outLen, v, w);
                        toks[ntok] = new Token { Type = TDigits0, Int = v, Str = w };
                        break;
                    }
                    case TDigits:
                    {
                        uint v = ReadU32(streams, ntok, TDigits);
                        AppendNumber(ref output, ref outLen, v, 0);
                        toks[ntok] = new Token { Type = TDigits, Int = v };
                        break;
                    }
                    case TDelta:
                    {
                        if (ntok >= prev.NTok) throw new InvalidDataException("tok3: DELTA without a previous token.");
                        uint v = ReadU8(streams, ntok, TDelta) + prev.Tokens[ntok].Int;
                        AppendNumber(ref output, ref outLen, v, 0);
                        toks[ntok] = new Token { Type = TDigits, Int = v };
                        break;
                    }
                    case TNop:
                        toks[ntok] = new Token { Type = TNop };
                        break;
                    case TMatch:
                    {
                        if (ntok >= prev.NTok) throw new InvalidDataException("tok3: MATCH without a previous token.");
                        var pt = prev.Tokens[ntok];
                        switch (pt.Type)
                        {
                            case TChar: AppendByte(ref output, ref outLen, (byte)pt.Int); toks[ntok] = pt; break;
                            case TAlpha:
                            {
                                int off = outLen - start;
                                Append(ref output, ref outLen, output.AsSpan(prev.Start + pt.Str, pt.Len).ToArray());
                                toks[ntok] = new Token { Type = TAlpha, Str = off, Len = pt.Len };
                                break;
                            }
                            case TDigits: AppendNumber(ref output, ref outLen, pt.Int, 0); toks[ntok] = pt; break;
                            case TDigits0: AppendNumber(ref output, ref outLen, pt.Int, pt.Str); toks[ntok] = pt; break;
                            default: throw new InvalidDataException("tok3: MATCH of an unmatchable token.");
                        }
                        break;
                    }
                    default: // END, or an elided END
                        AppendByte(ref output, ref outLen, 0);
                        toks[ntok] = new Token { Type = TEnd };
                        ended = true;
                        break;
                }
                if (ended) break;
            }
            if (!ended) throw new InvalidDataException("tok3: name has too many tokens.");
            names[n] = (start, outLen - start, toks, ntok);
        }
        return outLen == output.Length ? output : output[..outLen];
    }

    private static Stream Get(Dictionary<int, Stream> s, int ntok, int type)
        => s.TryGetValue(ntok << 4 | type, out var st) ? st : throw new InvalidDataException($"tok3: missing stream {ntok}/{type}.");

    private static int ReadType(Dictionary<int, Stream> s, int ntok)
    {
        if (!s.TryGetValue(ntok << 4, out var st) || st.Pos >= st.Data.Length) return -1;
        return st.Data[st.Pos++];
    }

    private static byte ReadU8(Dictionary<int, Stream> s, int ntok, int type)
    {
        var st = Get(s, ntok, type);
        if (st.Pos >= st.Data.Length) throw new InvalidDataException("tok3: stream exhausted.");
        return st.Data[st.Pos++];
    }

    private static uint ReadU32(Dictionary<int, Stream> s, int ntok, int type)
    {
        var st = Get(s, ntok, type);
        if (st.Pos + 4 > st.Data.Length) throw new InvalidDataException("tok3: stream exhausted.");
        uint v = (uint)(st.Data[st.Pos] | st.Data[st.Pos + 1] << 8 | st.Data[st.Pos + 2] << 16 | st.Data[st.Pos + 3] << 24);
        st.Pos += 4;
        return v;
    }

    private static void Grow(ref byte[] o, int need)
    {
        if (need > o.Length) Array.Resize(ref o, Math.Max(need, o.Length * 2 + 64));
    }

    private static void AppendByte(ref byte[] o, ref int len, byte b)
    {
        Grow(ref o, len + 1);
        o[len++] = b;
    }

    private static void Append(ref byte[] o, ref int len, byte[] data)
    {
        Grow(ref o, len + data.Length);
        data.CopyTo(o, len);
        len += data.Length;
    }

    private static void AppendNumber(ref byte[] o, ref int len, uint v, int width)
    {
        var s = v.ToString(System.Globalization.CultureInfo.InvariantCulture);
        if (s.Length < width) s = new string('0', width - s.Length) + s;
        Grow(ref o, len + s.Length);
        foreach (char ch in s) o[len++] = (byte)ch;
    }
}
