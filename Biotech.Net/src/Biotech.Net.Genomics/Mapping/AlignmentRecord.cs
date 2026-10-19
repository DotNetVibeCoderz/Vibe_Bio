using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>
/// One alignment as read from or written to SAM, BAM or CRAM. Coordinates are 0-based internally
/// (<see cref="Position"/> and <see cref="MatePosition"/> are -1 when absent); the SAM text form prints them 1-based.
/// Bases are ASCII in reference orientation (as stored in the file), qualities are Phred values without offset, and
/// the optional fields are kept in their compact BAM encoding and decoded on demand.
/// </summary>
public sealed class AlignmentRecord
{
    /// <summary>Read name (<c>QNAME</c>; <c>*</c> when unknown).</summary>
    public string Name { get; set; } = "*";

    /// <summary>SAM flags.</summary>
    public SamFlags Flags { get; set; }

    /// <summary>Index of the reference in the header, or -1.</summary>
    public int ReferenceId { get; set; } = -1;

    /// <summary>Reference name (<c>*</c> when unplaced). Writers resolve the reference by this name when it is set.</summary>
    public string ReferenceName { get; set; } = "*";

    /// <summary>0-based leftmost aligned position, or -1.</summary>
    public int Position { get; set; } = -1;

    /// <summary>Mapping quality (255 = unavailable).</summary>
    public int MapQ { get; set; } = 255;

    /// <summary>CIGAR operations (empty for <c>*</c>).</summary>
    public CigarElement[] Cigar { get; set; } = [];

    /// <summary>Reference index of the mate / next segment, or -1.</summary>
    public int MateReferenceId { get; set; } = -1;

    /// <summary>Reference name of the mate / next segment (<c>*</c> when none; never <c>=</c>).</summary>
    public string MateReferenceName { get; set; } = "*";

    /// <summary>0-based position of the mate / next segment, or -1.</summary>
    public int MatePosition { get; set; } = -1;

    /// <summary>Observed template length (<c>TLEN</c>).</summary>
    public int TemplateLength { get; set; }

    /// <summary>Bases (ASCII); empty for <c>*</c>.</summary>
    public byte[] Sequence { get; set; } = [];

    /// <summary>Phred base qualities (no offset), or null for <c>*</c>.</summary>
    public byte[]? Qualities { get; set; }

    /// <summary>Optional fields in BAM binary encoding.</summary>
    public byte[] AuxData { get; set; } = [];

    // -------------------------------------------------------------------------------------- flag helpers

    /// <summary>Template has multiple segments (0x1).</summary>
    public bool IsPaired => (Flags & SamFlags.Paired) != 0;

    /// <summary>Each segment properly aligned (0x2).</summary>
    public bool IsProperPair => (Flags & SamFlags.ProperPair) != 0;

    /// <summary>Segment unmapped (0x4).</summary>
    public bool IsUnmapped => (Flags & SamFlags.Unmapped) != 0;

    /// <summary>Next segment unmapped (0x8).</summary>
    public bool IsMateUnmapped => (Flags & SamFlags.MateUnmapped) != 0;

    /// <summary>Sequence reverse-complemented (0x10).</summary>
    public bool IsReverse => (Flags & SamFlags.Reverse) != 0;

    /// <summary>Next segment reverse-complemented (0x20).</summary>
    public bool IsMateReverse => (Flags & SamFlags.MateReverse) != 0;

    /// <summary>First segment in the template (0x40).</summary>
    public bool IsFirstOfPair => (Flags & SamFlags.First) != 0;

    /// <summary>Last segment in the template (0x80).</summary>
    public bool IsSecondOfPair => (Flags & SamFlags.Last) != 0;

    /// <summary>Secondary alignment (0x100).</summary>
    public bool IsSecondary => (Flags & SamFlags.Secondary) != 0;

    /// <summary>Fails quality checks (0x200).</summary>
    public bool IsQcFail => (Flags & SamFlags.QcFail) != 0;

    /// <summary>PCR or optical duplicate (0x400).</summary>
    public bool IsDuplicate => (Flags & SamFlags.Duplicate) != 0;

    /// <summary>Supplementary alignment (0x800).</summary>
    public bool IsSupplementary => (Flags & SamFlags.Supplementary) != 0;

    /// <summary>Neither secondary nor supplementary.</summary>
    public bool IsPrimary => (Flags & (SamFlags.Secondary | SamFlags.Supplementary)) == 0;

    /// <summary>Whether the record is placed on a reference (it may still be unmapped).</summary>
    public bool IsPlaced => ReferenceId >= 0 || ReferenceName != "*";

    // -------------------------------------------------------------------------------------- coordinates

    /// <summary>Reference bases covered by the CIGAR (M, D, N, =, X).</summary>
    public int ReferenceLength => Mapping.Cigar.ReferenceLength(Cigar);

    /// <summary>Read bases consumed by the CIGAR (M, I, S, =, X).</summary>
    public int QueryLength => Mapping.Cigar.QueryLength(Cigar);

    /// <summary>
    /// Exclusive 0-based end on the reference, as htslib's <c>bam_endpos</c>: unmapped reads and reads whose CIGAR
    /// consumes no reference count as one base.
    /// </summary>
    public int End
    {
        get
        {
            int rlen = IsUnmapped ? 0 : ReferenceLength;
            return Position + (rlen == 0 ? 1 : rlen);
        }
    }

    /// <summary>The 0-based half-open reference interval covered.</summary>
    public GenomicInterval Interval => new(ReferenceName, Position, End, IsReverse ? Strand.Reverse : Strand.Forward);

    /// <summary>Whether the record overlaps <c>[start, end)</c> on its reference (htslib region semantics).</summary>
    public bool Overlaps(long start, long end) => Position < end && End > start;

    /// <summary>The aligned (query position, reference position) pairs; insertions/soft clips have a null reference position, deletions/skips a null query position.</summary>
    public IEnumerable<(int? QueryPosition, int? ReferencePosition)> AlignedPairs(bool matchesOnly = false)
        => Mapping.Cigar.AlignedPairs(Cigar, Position, matchesOnly);

    // -------------------------------------------------------------------------------------- tags

    /// <summary>The decoded optional fields, in file order.</summary>
    public IReadOnlyList<AuxField> Tags => AuxCodec.Decode(AuxData);

    /// <summary>Field <paramref name="tag"/>, or null.</summary>
    public AuxField? GetTag(string tag)
    {
        int p = AuxCodec.Find(AuxData, tag);
        return p < 0 ? null : AuxCodec.DecodeOne(AuxData.AsSpan(p, AuxCodec.FieldLength(AuxData, p)));
    }

    /// <summary>Whether field <paramref name="tag"/> is present.</summary>
    public bool HasTag(string tag) => AuxCodec.Find(AuxData, tag) >= 0;

    /// <summary>Integer field <paramref name="tag"/>, or null.</summary>
    public long? GetInt(string tag) => GetTag(tag) is { IsInteger: true } f ? f.AsInt64() : null;

    /// <summary>String (<c>Z</c>/<c>H</c>/<c>A</c>) field <paramref name="tag"/>, or null.</summary>
    public string? GetString(string tag) => GetTag(tag) is { } f && f.Type is 'Z' or 'H' or 'A' ? f.AsString() : null;

    /// <summary>Adds or replaces a field (a replaced field keeps its position).</summary>
    public void SetTag(AuxField field)
    {
        var buf = new ByteBuffer(AuxData.Length + 16);
        int p = AuxCodec.Find(AuxData, field.Tag);
        if (p < 0)
        {
            buf.Append(AuxData);
            AuxCodec.Encode(buf, field);
        }
        else
        {
            int len = AuxCodec.FieldLength(AuxData, p);
            buf.Append(AuxData.AsSpan(0, p));
            AuxCodec.Encode(buf, field);
            buf.Append(AuxData.AsSpan(p + len));
        }
        AuxData = buf.ToArray();
    }

    /// <summary>Removes field <paramref name="tag"/>; returns whether it was present.</summary>
    public bool RemoveTag(string tag)
    {
        int p = AuxCodec.Find(AuxData, tag);
        if (p < 0) return false;
        int len = AuxCodec.FieldLength(AuxData, p);
        var a = new byte[AuxData.Length - len];
        AuxData.AsSpan(0, p).CopyTo(a);
        AuxData.AsSpan(p + len).CopyTo(a.AsSpan(p));
        AuxData = a;
        return true;
    }

    /// <summary>Replaces all fields.</summary>
    public void SetTags(IEnumerable<AuxField> fields)
    {
        var buf = new ByteBuffer(64);
        foreach (var f in fields) AuxCodec.Encode(buf, f);
        AuxData = buf.ToArray();
    }

    // -------------------------------------------------------------------------------------- conversions

    /// <summary>The bases as text (<c>*</c> when absent).</summary>
    public string SequenceString => Sequence.Length == 0 ? "*" : Encoding.ASCII.GetString(Sequence);

    /// <summary>The qualities as Phred+33 text (<c>*</c> when absent).</summary>
    public string QualityString => Qualities is null ? "*"
        : string.Create(Qualities.Length, Qualities, static (s, q) => { for (int i = 0; i < s.Length; i++) s[i] = (char)(q[i] + 33); });

    /// <summary>
    /// The read as a <see cref="SequenceRead"/>; with <paramref name="originalOrientation"/> a reverse-strand read is
    /// reverse-complemented back to the orientation it was sequenced in.
    /// </summary>
    public SequenceRead ToSequenceRead(bool originalOrientation = true)
    {
        var bases = (byte[])Sequence.Clone();
        var quals = Qualities is null ? new byte[bases.Length] : (byte[])Qualities.Clone();
        if (Qualities is null) Array.Fill(quals, (byte)0);
        if (originalOrientation && IsReverse)
        {
            Array.Reverse(bases);
            for (int i = 0; i < bases.Length; i++) bases[i] = Complement(bases[i]);
            Array.Reverse(quals);
        }
        var seq = Core.Sequence.FromBuffer(Name, bases, Alphabet.Dna, validate: false);
        return SequenceRead.FromBuffer(seq, quals);
    }

    private static byte Complement(byte b) => b switch
    {
        (byte)'A' => (byte)'T', (byte)'C' => (byte)'G', (byte)'G' => (byte)'C', (byte)'T' => (byte)'A',
        (byte)'a' => (byte)'t', (byte)'c' => (byte)'g', (byte)'g' => (byte)'c', (byte)'t' => (byte)'a',
        _ => b,
    };

    /// <summary>The string-based <see cref="SamRecord"/> (1-based positions).</summary>
    public SamRecord ToSamRecord()
    {
        var tags = AuxCodec.Decode(AuxData).Select(t => t.ToString()).ToList();
        return new SamRecord(Name, Flags, ReferenceName, Position + 1, MapQ, Mapping.Cigar.Format(Cigar),
            MateReferenceName == "*" ? "*" : MateReferenceName == ReferenceName ? "=" : MateReferenceName,
            MatePosition + 1, TemplateLength, SequenceString, QualityString) { Tags = tags };
    }

    /// <summary>Converts a string-based <see cref="SamRecord"/>; reference ids are resolved through <paramref name="header"/> when given.</summary>
    public static AlignmentRecord FromSamRecord(SamRecord record, SamHeader? header = null) => Parse(record.ToString(), header);

    /// <summary>A deep copy.</summary>
    public AlignmentRecord Clone()
    {
        var r = (AlignmentRecord)MemberwiseClone();
        r.Cigar = (CigarElement[])Cigar.Clone();
        r.Sequence = (byte[])Sequence.Clone();
        r.Qualities = (byte[]?)Qualities?.Clone();
        r.AuxData = (byte[])AuxData.Clone();
        return r;
    }

    // -------------------------------------------------------------------------------------- SAM text

    /// <summary>The SAM line (no newline).</summary>
    public override string ToString()
    {
        var sb = new StringBuilder(256);
        AppendSam(sb);
        return sb.ToString();
    }

    /// <summary>Appends the SAM line (no newline).</summary>
    public void AppendSam(StringBuilder sb)
    {
        var inv = CultureInfo.InvariantCulture;
        sb.Append(Name).Append('\t').Append(((int)Flags).ToString(inv)).Append('\t').Append(ReferenceName).Append('\t')
          .Append((Position + 1).ToString(inv)).Append('\t').Append(MapQ.ToString(inv)).Append('\t');
        if (Cigar.Length == 0) sb.Append('*');
        else foreach (var e in Cigar) sb.Append(e.Length.ToString(inv)).Append(e.Op);
        sb.Append('\t');
        if (MateReferenceName == "*") sb.Append('*');
        else if (MateReferenceName == ReferenceName) sb.Append('=');
        else sb.Append(MateReferenceName);
        sb.Append('\t').Append((MatePosition + 1).ToString(inv)).Append('\t').Append(TemplateLength.ToString(inv)).Append('\t');
        if (Sequence.Length == 0) sb.Append('*');
        else foreach (byte b in Sequence) sb.Append((char)b);
        sb.Append('\t');
        if (Qualities is null) sb.Append('*');
        else foreach (byte q in Qualities) sb.Append((char)(q + 33));
        AuxCodec.AppendSamText(sb, AuxData);
    }

    /// <summary>Parses a SAM alignment line; reference ids are resolved through <paramref name="header"/> when given.</summary>
    public static AlignmentRecord Parse(string line, SamHeader? header = null) => Parse(line.AsSpan(), header);

    /// <summary>Parses a SAM alignment line.</summary>
    public static AlignmentRecord Parse(ReadOnlySpan<char> line, SamHeader? header = null)
    {
        if (!line.IsEmpty && line[^1] == '\r') line = line[..^1];
        if (line.Count('\t') < 10) throw new InvalidDataException("A SAM line has at least 11 tab-separated columns.");
        var s = line;
        var qname = Next(ref s);
        var flag = Next(ref s);
        var rname = Next(ref s);
        var pos = Next(ref s);
        var mapq = Next(ref s);
        var cigar = Next(ref s);
        var rnext = Next(ref s);
        var pnext = Next(ref s);
        var tlen = Next(ref s);
        var seq = Next(ref s);
        var qual = Next(ref s);

        var r = new AlignmentRecord
        {
            Name = qname.ToString(),
            Flags = (SamFlags)ParseInt(flag, "FLAG"),
            Position = (int)(ParseLong(pos, "POS") - 1),
            MapQ = ParseInt(mapq, "MAPQ"),
            MatePosition = (int)(ParseLong(pnext, "PNEXT") - 1),
            TemplateLength = ParseInt(tlen, "TLEN"),
        };
        r.ReferenceName = rname is "*" ? "*" : header?.References is { Count: > 0 } ? ResolveName(header, rname, out _) : rname.ToString();
        r.ReferenceId = header is null || r.ReferenceName == "*" ? -1 : header.GetReferenceId(r.ReferenceName);
        if (rnext is "*") r.MateReferenceName = "*";
        else if (rnext is "=") { r.MateReferenceName = r.ReferenceName; r.MateReferenceId = r.ReferenceId; }
        else
        {
            r.MateReferenceName = header?.References is { Count: > 0 } ? ResolveName(header, rnext, out _) : rnext.ToString();
            r.MateReferenceId = header?.GetReferenceId(r.MateReferenceName) ?? -1;
        }
        r.Cigar = ParseCigar(cigar);
        if (!(seq is "*"))
        {
            var b = new byte[seq.Length];
            for (int i = 0; i < b.Length; i++) b[i] = (byte)seq[i];
            r.Sequence = b;
        }
        if (!(qual is "*"))
        {
            var q = new byte[qual.Length];
            for (int i = 0; i < q.Length; i++)
            {
                int v = qual[i] - 33;
                if (v is < 0 or > 93) throw new InvalidDataException($"Invalid quality character '{qual[i]}'.");
                q[i] = (byte)v;
            }
            if (r.Sequence.Length > 0 && q.Length != r.Sequence.Length)
                throw new InvalidDataException($"Read '{r.Name}': SEQ and QUAL lengths differ ({r.Sequence.Length} vs {q.Length}).");
            r.Qualities = q;
        }
        if (!s.IsEmpty)
        {
            var buf = new ByteBuffer(64);
            while (!s.IsEmpty)
            {
                var f = Next(ref s);
                if (!f.IsEmpty) AuxCodec.AppendSamField(buf, f);
            }
            r.AuxData = buf.ToArray();
        }
        return r;

        static int ParseInt(ReadOnlySpan<char> v, string what)
            => int.TryParse(v, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int x) ? x
               : throw new InvalidDataException($"Invalid SAM {what} '{v}'.");
        static long ParseLong(ReadOnlySpan<char> v, string what)
            => long.TryParse(v, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long x) ? x
               : throw new InvalidDataException($"Invalid SAM {what} '{v}'.");
    }

    private static string ResolveName(SamHeader header, ReadOnlySpan<char> name, out int id)
    {
        // reuse the header's string instance (no per-record allocation for reference names)
        for (id = 0; id < header.References.Count; id++)
            if (name.SequenceEqual(header.References[id].Name)) return header.References[id].Name;
        id = -1;
        return name.ToString();
    }

    private static ReadOnlySpan<char> Next(ref ReadOnlySpan<char> s)
    {
        int tab = s.IndexOf('\t');
        ReadOnlySpan<char> f;
        if (tab < 0) { f = s; s = default; }
        else { f = s[..tab]; s = s[(tab + 1)..]; }
        return f;
    }

    internal static CigarElement[] ParseCigar(ReadOnlySpan<char> cigar)
    {
        if (cigar.IsEmpty || cigar is "*") return [];
        int n = 0;
        foreach (char c in cigar) if (c is < '0' or > '9') n++;
        var ops = new CigarElement[n];
        int len = 0, k = 0;
        bool digits = false;
        foreach (char c in cigar)
        {
            if (c is >= '0' and <= '9') { len = checked(len * 10 + (c - '0')); digits = true; continue; }
            if ("MIDNSHP=X".IndexOf(c) < 0 || !digits) throw new InvalidDataException($"Invalid CIGAR '{cigar}'.");
            ops[k++] = new CigarElement(len, c);
            len = 0;
            digits = false;
        }
        if (digits) throw new InvalidDataException($"Invalid CIGAR '{cigar}' (trailing length).");
        return ops;
    }
}
