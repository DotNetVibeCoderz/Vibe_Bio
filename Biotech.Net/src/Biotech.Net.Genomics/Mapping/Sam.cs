using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>One CIGAR operation (<c>M</c>, <c>I</c>, <c>D</c>, <c>N</c>, <c>S</c>, <c>H</c>, <c>P</c>, <c>=</c>, <c>X</c>).</summary>
public readonly record struct CigarElement(int Length, char Op)
{
    /// <summary>Whether the operation consumes reference bases.</summary>
    public bool ConsumesReference => Op is 'M' or 'D' or 'N' or '=' or 'X';

    /// <summary>Whether the operation consumes read bases.</summary>
    public bool ConsumesQuery => Op is 'M' or 'I' or 'S' or '=' or 'X';

    public override string ToString() => $"{Length}{Op}";
}

/// <summary>CIGAR helpers.</summary>
public static class Cigar
{
    /// <summary>Parses a CIGAR string (<c>*</c> gives an empty array).</summary>
    public static CigarElement[] Parse(string cigar)
    {
        if (cigar is "*" or "") return [];
        var list = new List<CigarElement>();
        int n = 0;
        foreach (char c in cigar)
        {
            if (c is >= '0' and <= '9') { n = n * 10 + (c - '0'); continue; }
            list.Add(new CigarElement(n, c));
            n = 0;
        }
        return list.ToArray();
    }

    /// <summary>Formats a CIGAR (<c>*</c> when empty).</summary>
    public static string Format(ReadOnlySpan<CigarElement> cigar)
    {
        if (cigar.IsEmpty) return "*";
        var sb = new StringBuilder();
        foreach (var e in cigar) sb.Append(e.Length).Append(e.Op);
        return sb.ToString();
    }

    /// <summary>Reference bases spanned.</summary>
    public static int ReferenceLength(ReadOnlySpan<CigarElement> cigar)
    {
        int n = 0;
        foreach (var e in cigar) if (e.ConsumesReference) n += e.Length;
        return n;
    }

    /// <summary>Read bases consumed.</summary>
    public static int QueryLength(ReadOnlySpan<CigarElement> cigar)
    {
        int n = 0;
        foreach (var e in cigar) if (e.ConsumesQuery) n += e.Length;
        return n;
    }

    /// <summary>Run-length encodes single-character operations.</summary>
    internal static CigarElement[] Compress(IEnumerable<char> ops)
    {
        var list = new List<CigarElement>();
        foreach (char op in ops)
        {
            if (list.Count > 0 && list[^1].Op == op) list[^1] = list[^1] with { Length = list[^1].Length + 1 };
            else list.Add(new CigarElement(1, op));
        }
        return list.ToArray();
    }
}

/// <summary>SAM flag bits.</summary>
[Flags]
public enum SamFlags
{
    None = 0,
    Paired = 0x1,
    ProperPair = 0x2,
    Unmapped = 0x4,
    MateUnmapped = 0x8,
    Reverse = 0x10,
    MateReverse = 0x20,
    First = 0x40,
    Last = 0x80,
    Secondary = 0x100,
    QcFail = 0x200,
    Duplicate = 0x400,
    Supplementary = 0x800,
}

/// <summary>One SAM alignment line. <see cref="Position"/> is 1-based (0 for unmapped), as in the file.</summary>
public sealed record SamRecord(string QueryName, SamFlags Flag, string ReferenceName, long Position, int MapQ, string Cigar,
    string MateReferenceName, long MatePosition, long TemplateLength, string Sequence, string Quality)
{
    /// <summary>Optional fields verbatim (<c>NM:i:1</c>, <c>MD:Z:10A5</c>, ...).</summary>
    public List<string> Tags { get; init; } = [];

    /// <summary>Whether the read is unmapped.</summary>
    public bool IsUnmapped => (Flag & SamFlags.Unmapped) != 0;

    /// <summary>Whether the read maps to the reverse strand.</summary>
    public bool IsReverse => (Flag & SamFlags.Reverse) != 0;

    /// <summary>The value of tag <paramref name="name"/> (e.g. <c>NM</c>), or null.</summary>
    public string? Tag(string name)
    {
        foreach (var t in Tags)
            if (t.Length > 5 && t.StartsWith(name, StringComparison.Ordinal) && t[2] == ':') return t[5..];
        return null;
    }

    /// <summary>The parsed CIGAR.</summary>
    public CigarElement[] CigarElements => Mapping.Cigar.Parse(Cigar);

    /// <summary>The SAM line (no newline).</summary>
    public override string ToString()
    {
        var sb = new StringBuilder(256);
        sb.Append(QueryName).Append('\t').Append((int)Flag).Append('\t').Append(ReferenceName).Append('\t')
          .Append(Position).Append('\t').Append(MapQ).Append('\t').Append(Cigar).Append('\t')
          .Append(MateReferenceName).Append('\t').Append(MatePosition).Append('\t').Append(TemplateLength).Append('\t')
          .Append(Sequence).Append('\t').Append(Quality);
        foreach (var t in Tags) sb.Append('\t').Append(t);
        return sb.ToString();
    }

    /// <summary>Parses a SAM line.</summary>
    public static SamRecord Parse(string line)
    {
        var c = line.Split('\t');
        if (c.Length < 11) throw new InvalidDataException("A SAM line has at least 11 columns.");
        return new SamRecord(c[0], (SamFlags)int.Parse(c[1], CultureInfo.InvariantCulture), c[2], long.Parse(c[3], CultureInfo.InvariantCulture),
            int.Parse(c[4], CultureInfo.InvariantCulture), c[5], c[6], long.Parse(c[7], CultureInfo.InvariantCulture),
            long.Parse(c[8], CultureInfo.InvariantCulture), c[9], c[10]) { Tags = c.Skip(11).ToList() };
    }
}

/// <summary>Writes SAM text.</summary>
public sealed class SamWriter : IDisposable
{
    private readonly TextWriter _w;
    private readonly bool _leaveOpen;

    /// <summary>Creates a writer and writes the header (<c>@HD</c>, one <c>@SQ</c> per reference, <c>@PG</c>).</summary>
    public SamWriter(TextWriter writer, IEnumerable<(string Name, long Length)> references, string program = "Biotech.Net", bool leaveOpen = false)
    {
        _w = writer;
        _leaveOpen = leaveOpen;
        _w.Write("@HD\tVN:1.6\tSO:unsorted\n");
        foreach (var (name, len) in references) _w.Write($"@SQ\tSN:{name}\tLN:{len}\n");
        _w.Write($"@PG\tID:{program}\tPN:{program}\n");
    }

    /// <summary>Creates a writer whose header lists <paramref name="reference"/>.</summary>
    public SamWriter(TextWriter writer, IEnumerable<Sequence> reference, bool leaveOpen = false)
        : this(writer, reference.Select(s => (s.Id, (long)s.Length)), leaveOpen: leaveOpen) { }

    /// <summary>Creates a SAM file (gzip when the path ends in <c>.gz</c>).</summary>
    public static SamWriter Create(string path, IEnumerable<Sequence> reference) => new(BioFile.CreateText(path), reference);

    /// <summary>Writes one record.</summary>
    public void Write(SamRecord record) { _w.Write(record.ToString()); _w.Write('\n'); }

    /// <summary>Writes one alignment.</summary>
    public void Write(ReadAlignment alignment) => Write(alignment.ToSam());

    /// <summary>Writes alignments.</summary>
    public void Write(IEnumerable<ReadAlignment> alignments) { foreach (var a in alignments) Write(a); }

    /// <summary>Flushes and closes (unless leaveOpen).</summary>
    public void Dispose()
    {
        _w.Flush();
        if (!_leaveOpen) _w.Dispose();
    }
}

/// <summary>Reads SAM text.</summary>
public static class SamReader
{
    /// <summary>Streams the alignment lines (header lines are skipped; use <see cref="ReadHeader"/>).</summary>
    public static IEnumerable<SamRecord> Read(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '@') continue;
            yield return SamRecord.Parse(line);
        }
    }

    /// <summary>Streams a SAM file (gzip-aware).</summary>
    public static IEnumerable<SamRecord> Read(string path)
    {
        using var r = BioFile.OpenText(path);
        foreach (var rec in Read(r)) yield return rec;
    }

    /// <summary>Reference names and lengths from the <c>@SQ</c> lines of a SAM text.</summary>
    public static List<(string Name, long Length)> ReadHeader(TextReader reader)
    {
        var refs = new List<(string, long)>();
        string? line;
        while ((line = reader.ReadLine()) is not null && line.StartsWith('@'))
        {
            if (!line.StartsWith("@SQ", StringComparison.Ordinal)) continue;
            string? name = null;
            long len = 0;
            foreach (var f in line.Split('\t'))
            {
                if (f.StartsWith("SN:", StringComparison.Ordinal)) name = f[3..];
                else if (f.StartsWith("LN:", StringComparison.Ordinal)) len = long.Parse(f[3..], CultureInfo.InvariantCulture);
            }
            if (name is not null) refs.Add((name, len));
        }
        return refs;
    }
}
