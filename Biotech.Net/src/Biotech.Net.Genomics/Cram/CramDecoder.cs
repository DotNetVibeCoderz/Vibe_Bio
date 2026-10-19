using System.Security.Cryptography;
using System.Text;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Genomics.Mapping;

namespace Biotech.Net.Genomics.Cram;

/// <summary>The compression header of a CRAM container: preservation map, data-series and tag encodings.</summary>
internal sealed class CompressionHeader
{
    public bool ReadNamesIncluded = true, ApDelta = true, ReferenceRequired = true;
    public readonly byte[,] Substitution = new byte[5, 4];
    public readonly List<string[]> TagDictionary = [];
    public readonly Dictionary<string, CramEncoding> Series = new(StringComparer.Ordinal);
    public readonly Dictionary<int, CramEncoding> Tags = [];

    private static readonly byte[] Bases = "ACGTN"u8.ToArray();

    public CramEncoding Get(string key) => Series.TryGetValue(key, out var e) ? e
        : throw new InvalidDataException($"CRAM data series {key} is used but has no encoding in the compression header.");

    public static CompressionHeader Parse(byte[] data)
    {
        var h = new CompressionHeader();
        // default substitution matrix: codes 0-3 are the other bases in ACGTN order
        for (int r = 0; r < 5; r++)
        {
            int k = 0;
            for (int b = 0; b < 5; b++) if (b != r) h.Substitution[r, k++] = Bases[b];
        }
        var c = new ByteCursor(data);
        c.Itf8(); // size in bytes
        int n = c.Itf8();
        for (int i = 0; i < n; i++)
        {
            string key = Encoding.ASCII.GetString(c.Take(2));
            switch (key)
            {
                case "RN": h.ReadNamesIncluded = c.U8() != 0; break;
                case "AP": h.ApDelta = c.U8() != 0; break;
                case "RR": h.ReferenceRequired = c.U8() != 0; break;
                case "SM":
                {
                    var sm = c.Take(5);
                    for (int r = 0; r < 5; r++)
                    {
                        int k = 0;
                        for (int b = 0; b < 5; b++)
                        {
                            if (b == r) continue;
                            int code = (sm[r] >> (6 - 2 * k)) & 3;
                            h.Substitution[r, code] = Bases[b];
                            k++;
                        }
                    }
                    break;
                }
                case "TD":
                {
                    int len = c.Itf8();
                    var td = c.Take(len);
                    int start = 0;
                    for (int p = 0; p <= td.Length; p++)
                    {
                        if (p < td.Length && td[p] != 0) continue;
                        if (p == td.Length && start == p && h.TagDictionary.Count > 0) break;
                        var entry = td[start..p];
                        var keys = new string[entry.Length / 3];
                        for (int k = 0; k < keys.Length; k++) keys[k] = Encoding.ASCII.GetString(entry.Slice(3 * k, 3));
                        h.TagDictionary.Add(keys);
                        start = p + 1;
                        if (p == td.Length) break;
                    }
                    break;
                }
                default: throw new InvalidDataException($"Unknown CRAM preservation map key '{key}'.");
            }
        }
        c.Itf8();
        n = c.Itf8();
        for (int i = 0; i < n; i++)
        {
            string key = Encoding.ASCII.GetString(c.Take(2));
            h.Series[key] = CramEncoding.Parse(ref c);
        }
        c.Itf8();
        n = c.Itf8();
        for (int i = 0; i < n; i++)
        {
            int key = c.Itf8();
            h.Tags[key] = CramEncoding.Parse(ref c);
        }
        return h;
    }
}

/// <summary>The reference bases visible to a slice (1-based coordinates as in CRAM).</summary>
internal sealed class RefWindow(byte[] bases, long start, long end)
{
    public readonly byte[] Bases = bases;
    public readonly long Start = start; // 1-based position of Bases[0]
    public readonly long End = end;     // last usable 1-based position

    /// <summary>Base at 0-based position <paramref name="p0"/> (N outside the window).</summary>
    public byte At0(long p0)
    {
        long i = p0 + 1 - Start;
        return i >= 0 && i < Bases.Length && p0 + 1 <= End ? Bases[i] : (byte)'N';
    }
}

/// <summary>Decodes one CRAM slice into alignment records, following htslib's decoder (cram_decode.c).</summary>
internal sealed class SliceDecoder
{
    private readonly SamHeader _header;
    private readonly CompressionHeader _ch;
    private readonly Func<int, RefWindow?> _externalRef;
    private readonly CramReaderOptions _options;
    private readonly string _namePrefix;
    private readonly IReadOnlyList<string> _readGroups;
    private readonly int _unknownRg;

    public SliceDecoder(SamHeader header, CompressionHeader ch, Func<int, RefWindow?> externalRef, CramReaderOptions options, string namePrefix)
    {
        _header = header;
        _ch = ch;
        _externalRef = externalRef;
        _options = options;
        _namePrefix = namePrefix;
        _readGroups = header.ReadGroupIds;
        _unknownRg = _readGroups.Count > 0 && _readGroups[^1] == "UNKNOWN" ? _readGroups.Count - 1 : -2;
    }

    private sealed class Rec
    {
        public int Flags, CramFlags, RefId, Len, Rg = -1, MateFlags, MapQ;
        public long Apos, Aend, MatePos;
        public int MateLine = -1, MateRefId = -1;
        public long Tlen = long.MinValue, ExplicitTlen = long.MinValue;
        public string? Name;
        public byte[] Seq = [];
        public byte[]? Qual;
        public List<CigarElement> Cigar = [];
        public ByteBuffer Aux = new(32);
        public int HasMd, HasNm; // 0 absent, 1 verbatim, <0 placeholder at -offset
    }

    private const int CfPreserveQual = 1, CfDetached = 2, CfMateDownstream = 4, CfNoSeq = 8, CfExplicitTlen = 16;

    public List<AlignmentRecord> Decode(ReadOnlySpan<byte> container, int sliceOffset, int major)
    {
        var c = new ByteCursor(container) { Pos = sliceOffset };
        var hb = CramBlock.Read(ref c, major);
        if (hb.ContentType != CramBlock.SliceHeader) throw new InvalidDataException("CRAM: expected a slice header block.");
        var sh = new ByteCursor(hb.Data);
        int refId = sh.Itf8();
        long sliceStart = sh.Itf8();
        long sliceSpan = sh.Itf8();
        int nRecords = sh.Itf8();
        long recordCounter = major >= 3 ? CramVarint.Ltf8(ref sh) : sh.Itf8();
        int nBlocks = sh.Itf8();
        int nIds = sh.Itf8();
        for (int i = 0; i < nIds; i++) sh.Itf8();
        int embeddedId = sh.Itf8();
        var md5 = sh.Take(16).ToArray();

        var streams = new SliceStreams();
        for (int i = 0; i < nBlocks; i++)
        {
            var b = CramBlock.Read(ref c, major);
            if (b.ContentType == CramBlock.CoreData) streams.Core = new CoreBits(b.Data);
            else
            {
                streams.External[b.ContentId] = new ExternalStream(b.Data);
                if (b.Method == CramBlockMethod.NameTokenizer) streams.TokenizedBlocks.Add(b.ContentId);
            }
        }

        RefWindow? sliceRef = null;
        bool noRef = !_ch.ReferenceRequired;
        if (refId >= 0)
        {
            if (embeddedId >= 0)
            {
                var eb = streams.External.TryGetValue(embeddedId, out var e) ? e.Data : throw new InvalidDataException("CRAM: embedded reference block missing.");
                var up = ToUpper(eb);
                sliceRef = new RefWindow(up, sliceStart, sliceStart + sliceSpan - 1);
            }
            else if (!noRef)
            {
                sliceRef = _externalRef(refId) ?? throw new InvalidOperationException(
                    $"CRAM slice on '{_header.GetReferenceName(refId)}' needs the reference sequence; pass a FASTA (with .fai) to CramReader.Open.");
            }
            if (sliceRef is not null && _options.ValidateReferenceMd5 && md5.Any(x => x != 0))
            {
                long from = sliceStart - sliceRef.Start;
                long len = Math.Min(sliceSpan, sliceRef.Bases.Length - from);
                if (from < 0 || len < 0) throw new InvalidDataException("CRAM: slice lies outside its reference.");
                var digest = MD5.HashData(sliceRef.Bases.AsSpan((int)from, (int)len));
                if (!digest.AsSpan().SequenceEqual(md5))
                    throw new InvalidDataException($"CRAM: reference MD5 mismatch for '{_header.GetReferenceName(refId)}':{sliceStart}-{sliceStart + sliceSpan - 1} (wrong reference FASTA?).");
            }
        }

        var recs = new Rec[nRecords];
        long lastApos = sliceStart;
        int lastRefId = -9;
        RefWindow? recRef = sliceRef;
        for (int r = 0; r < nRecords; r++)
        {
            var cr = recs[r] = new Rec();
            cr.Flags = _ch.Get("BF").ReadInt(streams);
            cr.CramFlags = _ch.Get("CF").ReadInt(streams);
            int cf = cr.CramFlags;
            if (refId == -2)
            {
                cr.RefId = _ch.Get("RI").ReadInt(streams);
                if (cr.RefId >= 0 && cr.RefId != lastRefId)
                {
                    recRef = noRef ? null : _externalRef(cr.RefId) ?? throw new InvalidOperationException(
                        $"CRAM records on '{_header.GetReferenceName(cr.RefId)}' need the reference sequence; pass a FASTA (with .fai) to CramReader.Open.");
                    lastRefId = cr.RefId;
                }
            }
            else cr.RefId = refId;
            if (cr.RefId < -1 || cr.RefId >= _header.References.Count) throw new InvalidDataException($"CRAM record refers to reference id {cr.RefId}.");
            cr.Len = _ch.Get("RL").ReadInt(streams);
            if (cr.Len < 0) throw new InvalidDataException("CRAM: negative read length.");
            int ap = _ch.Get("AP").ReadInt(streams);
            cr.Apos = _ch.ApDelta ? lastApos + ap : ap;
            lastApos = cr.Apos;
            cr.Rg = _ch.Get("RG").ReadInt(streams);
            if (cr.Rg == _unknownRg) cr.Rg = -1;
            if (_ch.ReadNamesIncluded) cr.Name = Encoding.ASCII.GetString(_ch.Get("RN").ReadArray(streams));
            if ((cf & CfDetached) != 0)
            {
                cr.MateFlags = _ch.Series.ContainsKey("MF") ? _ch.Get("MF").ReadInt(streams) : 0;
                if (!_ch.ReadNamesIncluded) cr.Name = Encoding.ASCII.GetString(_ch.Get("RN").ReadArray(streams));
                cr.MateRefId = _ch.Get("NS").ReadInt(streams);
                cr.MatePos = _ch.Get("NP").ReadInt(streams);
                cr.Tlen = _ch.Get("TS").ReadInt(streams);
            }
            else if ((cf & CfMateDownstream) != 0)
            {
                cr.MateLine = _ch.Get("NF").ReadInt(streams) + r + 1;
                if ((cf & CfExplicitTlen) != 0) cr.ExplicitTlen = _ch.Get("TS").ReadInt(streams);
            }
            else if ((cf & CfExplicitTlen) != 0) cr.ExplicitTlen = _ch.Get("TS").ReadInt(streams);

            DecodeTags(streams, cr);

            var rref = refId == -2 ? (cr.RefId >= 0 ? recRef : null) : sliceRef;
            if ((cr.Flags & 4) == 0) DecodeMapped(streams, cr, rref);
            else
            {
                cr.Aend = cr.Apos;
                cr.MapQ = 0;
                cr.Seq = new byte[cr.Len];
                if (cr.Len > 0) _ch.Get("BA").ReadBytes(streams, cr.Seq);
                if ((cf & CfPreserveQual) != 0 && cr.Len > 0)
                {
                    cr.Qual = new byte[cr.Len];
                    _ch.Get("QS").ReadBytes(streams, cr.Qual);
                }
            }
            if ((cf & CfNoSeq) != 0) { cr.Seq = []; cr.Qual = null; }
        }

        ResolveMates(recs);

        var output = new List<AlignmentRecord>(nRecords);
        for (int r = 0; r < nRecords; r++)
        {
            var cr = recs[r];
            string name = cr.Name ?? GenerateName(recs, r, recordCounter);
            if (cr.Rg >= 0)
            {
                if (cr.Rg >= _readGroups.Count) throw new InvalidDataException($"CRAM read group index {cr.Rg} is not in the header.");
                cr.Aux.Append("RGZ"u8);
                cr.Aux.Append(Encoding.ASCII.GetBytes(_readGroups[cr.Rg]));
                cr.Aux.Append(0);
            }
            var qual = cr.Qual;
            if (qual is not null && (qual.Length == 0 || qual[0] == 0xFF)) qual = null;
            output.Add(new AlignmentRecord
            {
                Name = name,
                Flags = (SamFlags)cr.Flags,
                ReferenceId = cr.RefId,
                ReferenceName = _header.GetReferenceName(cr.RefId),
                Position = (int)cr.Apos - 1,
                MapQ = cr.MapQ,
                Cigar = cr.Cigar.ToArray(),
                MateReferenceId = cr.MateRefId,
                MateReferenceName = _header.GetReferenceName(cr.MateRefId),
                MatePosition = (int)cr.MatePos - 1,
                TemplateLength = (int)cr.Tlen,
                Sequence = cr.Seq,
                Qualities = qual,
                AuxData = cr.Aux.ToArray(),
            });
        }
        return output;
    }

    private string GenerateName(Rec[] recs, int r, long counter)
    {
        var cr = recs[r];
        if (cr.MateLine >= 0 && cr.MateLine < recs.Length && recs[cr.MateLine].Name is { Length: > 0 } mateName) return mateName;
        long n = cr.MateLine >= 0 && cr.MateLine < r ? counter + cr.MateLine + 1 : counter + r + 1;
        return _namePrefix + ":" + n.ToString(System.Globalization.CultureInfo.InvariantCulture);
    }

    private static byte[] ToUpper(byte[] b)
    {
        var o = new byte[b.Length];
        for (int i = 0; i < b.Length; i++) o[i] = b[i] is >= (byte)'a' and <= (byte)'z' ? (byte)(b[i] - 32) : b[i];
        return o;
    }

    private void DecodeTags(SliceStreams s, Rec cr)
    {
        int tl = _ch.Get("TL").ReadInt(s);
        if (tl < 0 || tl >= _ch.TagDictionary.Count) throw new InvalidDataException($"CRAM tag line {tl} is not in the tag dictionary.");
        foreach (var key in _ch.TagDictionary[tl])
        {
            int id = key[0] << 16 | key[1] << 8 | key[2];
            if (!_ch.Tags.TryGetValue(id, out var enc)) throw new InvalidDataException($"CRAM tag {key} has no encoding.");
            bool isMd = key[0] == 'M' && key[1] == 'D', isNm = key[0] == 'N' && key[1] == 'M';
            if (isMd) cr.HasMd = key[2] == '*' ? -(cr.Aux.Length + 3) : 1;
            if (isNm) cr.HasNm = key[2] == '*' ? -(cr.Aux.Length + 3) : 1;
            byte[] value = enc is ByteArrayLenEncoding or ByteArrayStopEncoding or NullEncoding ? enc.ReadArray(s) : ReadFixed(s, enc, key[2]);
            if (key == "cFC" && value.Length == 1)
            {
                // htslib private flags: bit 1 = do not regenerate MD, bit 2 = do not regenerate NM
                if ((value[0] & 1) != 0 && cr.HasMd == 0) cr.HasMd = 1;
                if ((value[0] & 2) != 0 && cr.HasNm == 0) cr.HasNm = 1;
                continue;
            }
            cr.Aux.Append((byte)key[0]);
            cr.Aux.Append((byte)key[1]);
            cr.Aux.Append((byte)key[2]);
            cr.Aux.Append(value);
        }
    }

    private static byte[] ReadFixed(SliceStreams s, CramEncoding enc, char type)
    {
        int n = type switch { 'A' or 'c' or 'C' => 1, 's' or 'S' => 2, 'i' or 'I' or 'f' => 4, _ => -1 };
        if (n < 0) throw new NotSupportedException($"CRAM tag of type {type} with a fixed-width encoding is not supported.");
        var v = new byte[n];
        enc.ReadBytes(s, v);
        return v;
    }

    private static int L1(byte b) => b switch
    {
        (byte)'A' or (byte)'a' => 0, (byte)'C' or (byte)'c' => 1, (byte)'G' or (byte)'g' => 2, (byte)'T' or (byte)'t' => 3, _ => 4,
    };

    private void DecodeMapped(SliceStreams s, Rec cr, RefWindow? rf)
    {
        int cf = cr.CramFlags;
        int len = cr.Len;
        var seq = new byte[len];
        var qual = new byte[len];
        if ((cf & CfPreserveQual) == 0) Array.Fill(qual, (byte)0xFF);
        if (rf is null) Array.Fill(seq, (byte)'=');
        bool doMd = _options.RegenerateMdNm;
        bool decodeMd = rf is not null && cr.RefId >= 0 && ((doMd && cr.HasMd == 0) || cr.HasMd < 0);
        bool decodeNm = rf is not null && cr.RefId >= 0 && ((doMd && cr.HasNm == 0) || cr.HasNm < 0);
        if ((cf & CfNoSeq) != 0) decodeMd = decodeNm = false;
        long contigLen = cr.RefId >= 0 ? _header.References[cr.RefId].Length : 0;
        var md = new StringBuilder();
        long mdDist = 0;
        long nm = 0;

        int fn = _ch.Get("FN").ReadInt(s);
        int prevPos = 0, seqPos = 1;
        long refPos = cr.Apos - 1;
        var cigar = cr.Cigar;
        char cigOp = 'M';
        long cigLen = 0;

        void Flush(char op)
        {
            if (cigLen != 0 && cigOp != op) { cigar.Add(new CigarElement((int)cigLen, cigOp)); cigLen = 0; }
        }
        void AddMd(byte refBase)
        {
            if (!decodeMd) return;
            md.Append(mdDist);
            md.Append((char)refBase);
            mdDist = 0;
        }
        void Put(byte[] a, int idx, byte v) { if (idx >= 0 && idx < a.Length) a[idx] = v; }
        void QualDefault() { if ((cf & CfPreserveQual) == 0 && len > 0 && qual[0] == 0xFF) Array.Fill(qual, (byte)30); }

        for (int f = 0; f < fn; f++)
        {
            byte op = _ch.Get("FC").ReadByte(s);
            int pos = _ch.Get("FP").ReadInt(s) + prevPos;
            if (pos > seqPos)
            {
                int n = pos - seqPos;
                if (rf is not null && cr.RefId >= 0)
                {
                    if (refPos + n > contigLen)
                    {
                        long rlen = contigLen - refPos;
                        if (rlen > 0)
                        {
                            for (int i = 0; i < rlen; i++) Put(seq, seqPos - 1 + i, rf.At0(refPos + i));
                            for (long i = rlen; i < n; i++) Put(seq, (int)(seqPos - 1 + i), (byte)'N');
                        }
                        else for (int i = seqPos - 1; i < len; i++) seq[i] = (byte)'N';
                        if (mdDist >= 0) mdDist += n;
                    }
                    else
                    {
                        for (int i = 0; i < n; i++)
                        {
                            byte b = rf.At0(refPos + i);
                            if (decodeMd || decodeNm)
                            {
                                if (b == 'N') { AddMd((byte)'N'); nm++; }
                                else mdDist++;
                            }
                            Put(seq, seqPos - 1 + i, b);
                        }
                    }
                }
                Flush('M');
                cigOp = 'M';
                cigLen += n;
                refPos += n;
                seqPos = pos;
            }
            prevPos = pos;
            switch ((char)op)
            {
                case 'S':
                {
                    var sc = _ch.Get("SC").ReadArray(s);
                    if (cigLen != 0) { cigar.Add(new CigarElement((int)cigLen, cigOp)); cigLen = 0; }
                    for (int i = 0; i < sc.Length; i++) Put(seq, pos - 1 + i, sc[i]);
                    cigar.Add(new CigarElement(sc.Length, 'S'));
                    cigOp = 'S';
                    seqPos += sc.Length;
                    break;
                }
                case 'X':
                {
                    Flush('M');
                    int bs = _ch.Get("BS").ReadByte(s);
                    if (bs is < 0 or > 3) throw new InvalidDataException("CRAM: invalid base substitution code.");
                    if (cr.RefId < 0 || refPos >= contigLen || rf is null)
                    {
                        Put(seq, pos - 1, _ch.Substitution[4, bs]);
                        if (decodeMd || decodeNm)
                        {
                            if (mdDist >= 0 && decodeMd) md.Append(mdDist);
                            mdDist = -1;
                            nm--;
                        }
                    }
                    else
                    {
                        byte refCall = refPos + 1 <= rf.End ? rf.At0(refPos) : (byte)'N';
                        Put(seq, pos - 1, _ch.Substitution[L1(refCall), bs]);
                        AddMd(refCall);
                    }
                    cigOp = 'M';
                    nm++;
                    cigLen++;
                    seqPos++;
                    refPos++;
                    break;
                }
                case 'D':
                {
                    Flush('D');
                    int dl = _ch.Get("DL").ReadInt(s);
                    if (dl < 0) throw new InvalidDataException("CRAM: negative deletion length.");
                    if (decodeMd || decodeNm)
                    {
                        if (mdDist >= 0 && decodeMd) md.Append(mdDist);
                        if (refPos + dl <= contigLen)
                        {
                            if (decodeMd)
                            {
                                md.Append('^');
                                for (int i = 0; i < dl; i++) md.Append((char)rf!.At0(refPos + i));
                                mdDist = 0;
                            }
                            nm += dl;
                        }
                        else
                        {
                            long dlen;
                            if (contigLen >= refPos)
                            {
                                if (decodeMd)
                                {
                                    md.Append('^');
                                    for (long i = refPos; i < contigLen; i++) md.Append((char)rf!.At0(i));
                                    md.Append(0);
                                }
                                dlen = dl - (contigLen - refPos);
                                nm += dl - dlen;
                            }
                            mdDist = -1;
                        }
                    }
                    cigOp = 'D';
                    cigLen += dl;
                    refPos += dl;
                    break;
                }
                case 'I':
                {
                    Flush('I');
                    var ins = _ch.Get("IN").ReadArray(s);
                    for (int i = 0; i < ins.Length; i++) Put(seq, pos - 1 + i, ins[i]);
                    cigOp = 'I';
                    cigLen += ins.Length;
                    seqPos += ins.Length;
                    nm += ins.Length;
                    break;
                }
                case 'i':
                {
                    Flush('I');
                    Put(seq, pos - 1, _ch.Get("BA").ReadByte(s));
                    cigOp = 'I';
                    cigLen++;
                    seqPos++;
                    nm++;
                    break;
                }
                case 'b':
                {
                    Flush('M');
                    var bb = _ch.Get("BB").ReadArray(s);
                    for (int i = 0; i < bb.Length; i++) Put(seq, pos - 1 + i, bb[i]);
                    if (decodeMd || decodeNm)
                    {
                        if (mdDist >= 0 && decodeMd) md.Append(mdDist);
                        int x;
                        for (x = 0; x < bb.Length; x++)
                        {
                            if (x > 0 && decodeMd) md.Append('0');
                            if (refPos + x >= contigLen || rf is null) { mdDist = -1; break; }
                            if (decodeMd) md.Append((char)rf.At0(refPos + x));
                        }
                        nm += x;
                        mdDist = 0;
                    }
                    cigOp = 'M';
                    cigLen += bb.Length;
                    seqPos += bb.Length;
                    refPos += bb.Length;
                    break;
                }
                case 'q':
                {
                    Flush('M');
                    var qq = _ch.Get("QQ").ReadArray(s);
                    QualDefault();
                    for (int i = 0; i < qq.Length; i++) Put(qual, pos - 1 + i, qq[i]);
                    cigOp = 'M';
                    break;
                }
                case 'B':
                {
                    Flush('M');
                    Put(seq, pos - 1, _ch.Get("BA").ReadByte(s));
                    if (decodeMd || decodeNm)
                    {
                        if (mdDist >= 0 && decodeMd) md.Append(mdDist);
                        if (refPos >= contigLen || rf is null) mdDist = -1;
                        else
                        {
                            if (decodeMd) md.Append((char)rf.At0(refPos));
                            nm++;
                            mdDist = 0;
                        }
                    }
                    byte q = _ch.Get("QS").ReadByte(s);
                    QualDefault();
                    Put(qual, pos - 1, q);
                    cigOp = 'M';
                    cigLen++;
                    seqPos++;
                    refPos++;
                    break;
                }
                case 'Q':
                {
                    byte q = _ch.Get("QS").ReadByte(s);
                    QualDefault();
                    Put(qual, pos - 1, q);
                    break;
                }
                case 'H':
                {
                    Flush('H');
                    int hc = _ch.Get("HC").ReadInt(s);
                    if (hc < 0) throw new InvalidDataException("CRAM: negative hard-clip length.");
                    cigOp = 'H';
                    cigLen += hc;
                    break;
                }
                case 'P':
                {
                    Flush('P');
                    int pd = _ch.Get("PD").ReadInt(s);
                    if (pd < 0) throw new InvalidDataException("CRAM: negative padding length.");
                    cigOp = 'P';
                    cigLen += pd;
                    break;
                }
                case 'N':
                {
                    Flush('N');
                    int rs = _ch.Get("RS").ReadInt(s);
                    if (rs < 0) throw new InvalidDataException("CRAM: negative reference-skip length.");
                    cigOp = 'N';
                    cigLen += rs;
                    refPos += rs;
                    break;
                }
                default:
                    throw new InvalidDataException($"CRAM: unknown read feature code '{(char)op}'.");
            }
        }

        // implicit match for the rest of the read
        if (len >= seqPos)
        {
            int rem = len - seqPos + 1;
            if (rf is not null && cr.RefId >= 0)
            {
                if (refPos + rem > contigLen)
                {
                    long rlen = contigLen - refPos;
                    if (rlen > 0)
                    {
                        if (seqPos - 1 + rlen < len)
                            for (int i = 0; i < rlen; i++) seq[seqPos - 1 + i] = rf.At0(refPos + i);
                        for (long i = rlen; i < rem; i++) seq[seqPos - 1 + i] = (byte)'N';
                    }
                    else for (int i = seqPos - 1; i < len; i++) seq[i] = (byte)'N';
                    if (mdDist >= 0) mdDist += rem;
                }
                else
                {
                    for (int i = 0; i < rem; i++)
                    {
                        byte b = rf.At0(refPos + i);
                        if (decodeMd || decodeNm)
                        {
                            if (b == 'N') { AddMd((byte)'N'); nm++; }
                            else mdDist++;
                        }
                        seq[seqPos - 1 + i] = b;
                    }
                    refPos += rem;
                }
            }
            else if (cr.RefId >= 0) refPos += rem;
            Flush('M');
            cigOp = 'M';
            cigLen += rem;
        }
        if (decodeMd && mdDist >= 0) md.Append(mdDist);
        if (cigLen != 0) cigar.Add(new CigarElement((int)cigLen, cigOp));
        cr.Aend = refPos > cr.Apos ? refPos : cr.Apos;
        cr.MapQ = _ch.Get("MQ").ReadInt(s);
        if ((cf & CfPreserveQual) != 0 && len > 0) _ch.Get("QS").ReadBytes(s, qual);
        cr.Seq = seq;
        cr.Qual = qual;

        if (decodeMd)
        {
            var mdBytes = new ByteBuffer(md.Length + 4);
            if (cr.HasMd == 0) mdBytes.Append("MDZ"u8);
            mdBytes.Append(Encoding.ASCII.GetBytes(md.ToString()));
            mdBytes.Append(0);
            if (cr.HasMd < 0) InsertAux(cr, -cr.HasMd, mdBytes.Span, ref cr.HasNm);
            else cr.Aux.Append(mdBytes.Span);
        }
        if (decodeNm)
        {
            uint v = (uint)nm;
            if (cr.HasNm == 0)
            {
                cr.Aux.Append((byte)'N');
                cr.Aux.Append((byte)'M');
                if (v <= byte.MaxValue) { cr.Aux.Append((byte)'C'); cr.Aux.Append((byte)v); }
                else if (v <= ushort.MaxValue) { cr.Aux.Append((byte)'S'); cr.Aux.Append((byte)v); cr.Aux.Append((byte)(v >> 8)); }
                else { cr.Aux.Append((byte)'I'); cr.Aux.Append((byte)v); cr.Aux.Append((byte)(v >> 8)); cr.Aux.Append((byte)(v >> 16)); cr.Aux.Append((byte)(v >> 24)); }
            }
            else
            {
                // placeholder NM*: fill a 32-bit value at its position
                int at = -cr.HasNm;
                var raw = cr.Aux.ToArray();
                raw[at - 1] = (byte)'I';
                var tmp = new ByteBuffer(raw.Length + 4);
                tmp.Append(raw.AsSpan(0, at));
                tmp.Append([(byte)v, (byte)(v >> 8), (byte)(v >> 16), (byte)(v >> 24)]);
                tmp.Append(raw.AsSpan(at));
                cr.Aux = tmp;
            }
        }
    }

    private static void InsertAux(Rec cr, int at, ReadOnlySpan<byte> bytes, ref int hasNm)
    {
        var raw = cr.Aux.ToArray();
        raw[at - 1] = (byte)'Z';
        var tmp = new ByteBuffer(raw.Length + bytes.Length);
        tmp.Append(raw.AsSpan(0, at));
        tmp.Append(bytes);
        tmp.Append(raw.AsSpan(at));
        cr.Aux = tmp;
        if (-hasNm > at) hasNm -= bytes.Length;
    }

    /// <summary>Mate cross-references inside the slice (htslib cram_decode_slice_xref).</summary>
    private static void ResolveMates(Rec[] recs)
    {
        int n = recs.Length;
        for (int rec = 0; rec < n; rec++)
        {
            var cr = recs[rec];
            if (cr.MateLine >= 0)
            {
                if (cr.MateLine >= n) throw new InvalidDataException("CRAM: mate line beyond the slice.");
                if (cr.Tlen == long.MinValue)
                {
                    int id1 = rec, id2 = rec;
                    long aleft = cr.Apos, aright = cr.Aend;
                    int refx = cr.RefId, leftCnt = 0, rightCnt = 0;
                    do
                    {
                        if (aleft > recs[id2].Apos) { aleft = recs[id2].Apos; leftCnt = 1; }
                        else if (aleft == recs[id2].Apos) leftCnt++;
                        if (aright < recs[id2].Aend) { aright = recs[id2].Aend; rightCnt = 1; }
                        else if (aright == recs[id2].Aend) rightCnt++;
                        if (recs[id2].MateLine == -1) { recs[id2].MateLine = rec; break; }
                        if (recs[id2].MateLine <= id2 || recs[id2].MateLine >= n) throw new InvalidDataException("CRAM: invalid mate chain.");
                        id2 = recs[id2].MateLine;
                        if (recs[id2].RefId != refx) refx = -1;
                    } while (id2 != id1);
                    if (refx != -1)
                    {
                        long tlen = aright - aleft + 1;
                        id1 = id2 = rec;
                        if (recs[id2].Apos == aleft && (recs[id2].Aend < aright || leftCnt <= 1))
                        {
                            recs[id2].Tlen = tlen;
                            tlen = -tlen;
                        }
                        else if (recs[id2].Apos == aleft && recs[id2].Aend == aright && leftCnt > 1 && rightCnt > 1)
                        {
                            if ((recs[id2].Flags & 0x40) != 0) { recs[id2].Tlen = tlen; tlen = -tlen; }
                            else recs[id2].Tlen = -tlen;
                        }
                        else recs[id2].Tlen = -tlen;
                        id2 = recs[id2].MateLine;
                        while (id2 != id1)
                        {
                            recs[id2].Tlen = tlen;
                            id2 = recs[id2].MateLine;
                        }
                    }
                    else
                    {
                        id1 = id2 = rec;
                        recs[id2].Tlen = 0;
                        id2 = recs[id2].MateLine;
                        while (id2 != id1)
                        {
                            recs[id2].Tlen = 0;
                            id2 = recs[id2].MateLine;
                        }
                    }
                }
                var mate = recs[cr.MateLine];
                cr.MatePos = mate.Apos;
                cr.MateRefId = mate.RefId;
                cr.Flags |= 0x1;
                if ((mate.Flags & 0x4) != 0) { cr.Flags |= 0x8; cr.Tlen = 0; }
                if ((cr.Flags & 0x4) != 0) cr.Tlen = 0;
                if ((mate.Flags & 0x10) != 0) cr.Flags |= 0x20;
            }
            else
            {
                if ((cr.MateFlags & 1) != 0) cr.Flags |= 0x1 | 0x20;
                if ((cr.MateFlags & 2) != 0) cr.Flags |= 0x8;
                if ((cr.Flags & 0x1) == 0) cr.MateRefId = -1;
            }
            if (cr.Tlen == long.MinValue) cr.Tlen = 0;
        }
        foreach (var cr in recs)
            if (cr.ExplicitTlen != long.MinValue) cr.Tlen = cr.ExplicitTlen;
    }
}
