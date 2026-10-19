using System.Buffers.Binary;
using System.Text;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>Conversion between <see cref="AlignmentRecord"/> and the BAM binary record layout.</summary>
internal static class BamCodec
{
    public const string CigarOps = "MIDNSHP=X";
    private const string Nt16 = "=ACMGRSVTWYHKDBN";

    /// <summary>htslib's <c>seq_nt16_table</c>: ASCII to 4-bit base code (unknown characters become N).</summary>
    private static readonly byte[] NtToCode = BuildNtToCode();

    /// <summary>A byte of packed bases to its two ASCII bases (high nibble first).</summary>
    private static readonly ushort[] PairToAscii = BuildPairs();

    private static readonly int[] OpToCode = BuildOpCodes();

    private static byte[] BuildNtToCode()
    {
        var t = new byte[256];
        Array.Fill(t, (byte)15);
        for (int i = 0; i < Nt16.Length; i++)
        {
            t[Nt16[i]] = (byte)i;
            t[char.ToLowerInvariant(Nt16[i])] = (byte)i;
        }
        t['U'] = t['u'] = 8;
        t['0'] = 1; t['1'] = 2; t['2'] = 4; t['3'] = 8; // colour-space digits, as htslib
        return t;
    }

    private static ushort[] BuildPairs()
    {
        var t = new ushort[256];
        for (int i = 0; i < 256; i++) t[i] = (ushort)(Nt16[i >> 4] | Nt16[i & 15] << 8);
        return t;
    }

    private static int[] BuildOpCodes()
    {
        var t = new int[128];
        Array.Fill(t, -1);
        for (int i = 0; i < CigarOps.Length; i++) t[CigarOps[i]] = i;
        return t;
    }

    /// <summary>Decodes one record (the bytes after <c>block_size</c>).</summary>
    public static AlignmentRecord Decode(ReadOnlySpan<byte> b, SamHeader header)
    {
        if (b.Length < 32) throw new InvalidDataException("Truncated BAM record.");
        int refId = BinaryPrimitives.ReadInt32LittleEndian(b);
        int pos = BinaryPrimitives.ReadInt32LittleEndian(b[4..]);
        int lName = b[8];
        int mapq = b[9];
        int nCigar = BinaryPrimitives.ReadUInt16LittleEndian(b[12..]);
        int flag = BinaryPrimitives.ReadUInt16LittleEndian(b[14..]);
        int lSeq = BinaryPrimitives.ReadInt32LittleEndian(b[16..]);
        int mRef = BinaryPrimitives.ReadInt32LittleEndian(b[20..]);
        int mPos = BinaryPrimitives.ReadInt32LittleEndian(b[24..]);
        int tlen = BinaryPrimitives.ReadInt32LittleEndian(b[28..]);
        int p = 32;
        int seqBytes = (lSeq + 1) >> 1;
        if (lSeq < 0 || p + lName + 4 * nCigar + seqBytes + lSeq > b.Length) throw new InvalidDataException("Corrupt BAM record (field lengths exceed the record).");
        if (refId < -1 || refId >= header.References.Count || mRef < -1 || mRef >= header.References.Count)
            throw new InvalidDataException($"BAM record refers to reference id {refId}/{mRef}, but the header has {header.References.Count}.");
        var r = new AlignmentRecord
        {
            Name = lName > 1 ? Encoding.ASCII.GetString(b.Slice(p, lName - 1)) : "",
            Flags = (SamFlags)flag,
            ReferenceId = refId,
            ReferenceName = header.GetReferenceName(refId),
            Position = pos,
            MapQ = mapq,
            MateReferenceId = mRef,
            MateReferenceName = header.GetReferenceName(mRef),
            MatePosition = mPos,
            TemplateLength = tlen,
        };
        p += lName;
        var cigar = new CigarElement[nCigar];
        for (int i = 0; i < nCigar; i++, p += 4)
        {
            uint c = BinaryPrimitives.ReadUInt32LittleEndian(b[p..]);
            int op = (int)(c & 0xF);
            if (op >= CigarOps.Length) throw new InvalidDataException($"Invalid BAM CIGAR operation code {op}.");
            cigar[i] = new CigarElement((int)(c >> 4), CigarOps[op]);
        }
        r.Cigar = cigar;
        if (lSeq > 0)
        {
            var seq = new byte[lSeq];
            var packed = b.Slice(p, seqBytes);
            int k = 0;
            for (int i = 0; i < lSeq >> 1; i++)
            {
                ushort two = PairToAscii[packed[i]];
                seq[k++] = (byte)two;
                seq[k++] = (byte)(two >> 8);
            }
            if ((lSeq & 1) != 0) seq[k] = (byte)Nt16[packed[seqBytes - 1] >> 4];
            r.Sequence = seq;
            p += seqBytes;
            var qual = b.Slice(p, lSeq);
            if (qual[0] != 0xFF) r.Qualities = qual.ToArray();
            p += lSeq;
        }
        r.AuxData = p < b.Length ? b[p..].ToArray() : [];
        // CIGARs with more than 65535 operations are stored in a CG:B,I tag with a placeholder kSmN CIGAR
        if (nCigar == 2 && cigar[0].Op == 'S' && cigar[0].Length == lSeq && cigar[1].Op == 'N' && r.GetTag("CG") is { } cg && cg.Value is uint[] ops)
        {
            var full = new CigarElement[ops.Length];
            for (int i = 0; i < ops.Length; i++) full[i] = new CigarElement((int)(ops[i] >> 4), CigarOps[(int)(ops[i] & 0xF)]);
            r.Cigar = full;
            r.RemoveTag("CG");
        }
        return r;
    }

    /// <summary>Encodes <paramref name="r"/> (without the leading <c>block_size</c>) into <paramref name="buf"/>.</summary>
    public static void Encode(AlignmentRecord r, SamHeader header, ByteBuffer buf, out int refId, out int end)
    {
        refId = ResolveId(header, r.ReferenceName, r.ReferenceId, "reference");
        int mRef = r.MateReferenceName == r.ReferenceName && r.MateReferenceName != "*" ? refId
            : ResolveId(header, r.MateReferenceName, r.MateReferenceId, "mate reference");
        var nameBytes = Encoding.ASCII.GetBytes(r.Name.Length == 0 ? "*" : r.Name);
        if (nameBytes.Length > 254) throw new InvalidDataException($"Read name '{r.Name}' is longer than 254 characters.");
        int lSeq = r.Sequence.Length;
        if (r.Qualities is not null && r.Qualities.Length != lSeq)
            throw new InvalidDataException($"Read '{r.Name}': {lSeq} bases but {r.Qualities.Length} qualities.");
        var cigar = r.Cigar;
        var aux = r.AuxData;
        if (cigar.Length > 0xFFFF)
        {
            // too many operations for the 16-bit count: move the CIGAR to a CG tag (as htslib)
            var ops = new uint[cigar.Length];
            for (int i = 0; i < ops.Length; i++) ops[i] = EncodeOp(cigar[i]);
            var tmp = r.Clone();
            tmp.SetTag(AuxField.Array("CG", ops));
            aux = tmp.AuxData;
            cigar = [new CigarElement(lSeq, 'S'), new CigarElement(r.ReferenceLength, 'N')];
        }
        int rlen = r.IsUnmapped ? 0 : Cigar.ReferenceLength(r.Cigar);
        end = r.Position + (rlen == 0 ? 1 : rlen);
        int bin = BinningIndex.RegionToBin(r.Position, end, BinningIndex.DefaultMinShift, BinningIndex.DefaultDepth);

        Span<byte> h = stackalloc byte[32];
        BinaryPrimitives.WriteInt32LittleEndian(h, refId);
        BinaryPrimitives.WriteInt32LittleEndian(h[4..], r.Position);
        h[8] = (byte)(nameBytes.Length + 1);
        h[9] = (byte)Math.Clamp(r.MapQ, 0, 255);
        BinaryPrimitives.WriteUInt16LittleEndian(h[10..], (ushort)bin);
        BinaryPrimitives.WriteUInt16LittleEndian(h[12..], (ushort)cigar.Length);
        BinaryPrimitives.WriteUInt16LittleEndian(h[14..], (ushort)r.Flags);
        BinaryPrimitives.WriteInt32LittleEndian(h[16..], lSeq);
        BinaryPrimitives.WriteInt32LittleEndian(h[20..], mRef);
        BinaryPrimitives.WriteInt32LittleEndian(h[24..], r.MatePosition);
        BinaryPrimitives.WriteInt32LittleEndian(h[28..], r.TemplateLength);
        buf.Append(h);
        buf.Append(nameBytes);
        buf.Append(0);
        Span<byte> four = stackalloc byte[4];
        foreach (var e in cigar)
        {
            BinaryPrimitives.WriteUInt32LittleEndian(four, EncodeOp(e));
            buf.Append(four);
        }
        var seq = r.Sequence;
        for (int i = 0; i + 1 < lSeq; i += 2) buf.Append((byte)(NtToCode[seq[i]] << 4 | NtToCode[seq[i + 1]]));
        if ((lSeq & 1) != 0) buf.Append((byte)(NtToCode[seq[lSeq - 1]] << 4));
        if (r.Qualities is null) for (int i = 0; i < lSeq; i++) buf.Append(0xFF);
        else buf.Append(r.Qualities);
        buf.Append(aux);
    }

    private static uint EncodeOp(CigarElement e)
    {
        int code = e.Op < 128 ? OpToCode[e.Op] : -1;
        if (code < 0) throw new InvalidDataException($"Invalid CIGAR operation '{e.Op}'.");
        return (uint)e.Length << 4 | (uint)code;
    }

    private static int ResolveId(SamHeader header, string name, int id, string what)
    {
        if (name is null or "*") return id >= 0 && id < header.References.Count ? id : -1;
        int i = header.GetReferenceId(name);
        if (i < 0) throw new InvalidDataException($"The {what} '{name}' is not in the header.");
        return i;
    }

    /// <summary>Reads the BAM header (magic, text, reference list) from <paramref name="s"/>.</summary>
    public static SamHeader ReadHeader(Stream s)
    {
        Span<byte> four = stackalloc byte[4];
        ReadExactly(s, four);
        if (!four.SequenceEqual("BAM\u0001"u8)) throw new InvalidDataException("Not a BAM file (magic 'BAM\\1' missing).");
        ReadExactly(s, four);
        int lText = BinaryPrimitives.ReadInt32LittleEndian(four);
        if (lText < 0) throw new InvalidDataException("Invalid BAM header text length.");
        var text = new byte[lText];
        ReadExactly(s, text);
        int nul = Array.IndexOf(text, (byte)0);
        string headerText = Encoding.UTF8.GetString(text, 0, nul < 0 ? lText : nul);
        ReadExactly(s, four);
        int nRef = BinaryPrimitives.ReadInt32LittleEndian(four);
        var refs = new List<SamReference>(nRef);
        for (int i = 0; i < nRef; i++)
        {
            ReadExactly(s, four);
            int lName = BinaryPrimitives.ReadInt32LittleEndian(four);
            var name = new byte[lName];
            ReadExactly(s, name);
            ReadExactly(s, four);
            int len = BinaryPrimitives.ReadInt32LittleEndian(four);
            refs.Add(new SamReference(Encoding.ASCII.GetString(name, 0, Math.Max(0, lName - 1)), (uint)len));
        }
        return new SamHeader(headerText, refs);
    }

    /// <summary>Writes the BAM header.</summary>
    public static void WriteHeader(Stream s, SamHeader header)
    {
        var w = new BinaryWriter(s, Encoding.ASCII, leaveOpen: true);
        w.Write("BAM\u0001"u8);
        var text = Encoding.UTF8.GetBytes(header.Text);
        w.Write(text.Length);
        w.Write(text);
        w.Write(header.References.Count);
        foreach (var r in header.References)
        {
            var name = Encoding.ASCII.GetBytes(r.Name);
            w.Write(name.Length + 1);
            w.Write(name);
            w.Write((byte)0);
            w.Write(checked((int)r.Length));
        }
        w.Flush();
    }

    private static void ReadExactly(Stream s, Span<byte> dest)
    {
        if (s.ReadAtLeast(dest, dest.Length, throwOnEndOfStream: false) != dest.Length)
            throw new EndOfStreamException("Truncated BAM header.");
    }
}
