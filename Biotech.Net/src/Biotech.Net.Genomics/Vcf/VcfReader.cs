using System.Globalization;
using Biotech.Net.Core;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Vcf;

/// <summary>Options for <see cref="VcfReader"/>.</summary>
public sealed record VcfReaderOptions
{
    /// <summary>Parse the GT field of every sample into <see cref="Variant.Genotypes"/> (default true).</summary>
    public bool ParseGenotypes { get; init; } = true;

    /// <summary>Parse the INFO column into <see cref="Variant.Info"/> (default true).</summary>
    public bool ParseInfo { get; init; } = true;
}

/// <summary>
/// Streaming VCF 4.x reader (plain or gzip/BGZF). The header is read on construction; records are produced lazily
/// by <see cref="Records"/>. Chromosome, FILTER and FORMAT strings are interned, genotypes are parsed with a
/// span scanner and the other per-sample fields are left unparsed until first requested.
/// </summary>
public sealed class VcfReader : IDisposable
{
    private readonly TextReader _reader;
    private readonly VcfReaderOptions _options;
    private readonly Dictionary<string, string> _strings = new(StringComparer.Ordinal);
    private readonly Dictionary<string, List<string>> _formats = new(StringComparer.Ordinal);
    private string? _pending;
    private long _lineNumber;

    /// <summary>Reads the header from <paramref name="reader"/>.</summary>
    public VcfReader(TextReader reader, VcfReaderOptions? options = null)
    {
        _reader = reader ?? throw new ArgumentNullException(nameof(reader));
        _options = options ?? new VcfReaderOptions();
        Header = ReadHeader();
    }

    /// <summary>Opens a VCF file; gzip/BGZF compression is detected from the content.</summary>
    /// <remarks>
    /// When the file is BGZF-compressed and a tabix index (<c>path.tbi</c> or <c>path.csi</c>) sits beside it,
    /// <see cref="Query(GenomicInterval)"/> uses the index instead of scanning.
    /// </remarks>
    public static VcfReader Open(string path, VcfReaderOptions? options = null)
    {
        var reader = new VcfReader(BioFile.OpenText(path), options);
        if (TabixReader.FindIndex(path) is { } idx && BgzfReader.IsBgzf(path)) reader._indexed = (path, idx);
        return reader;
    }

    /// <summary>Opens a BGZF-compressed VCF with an explicit tabix index (<c>.tbi</c> or <c>.csi</c>).</summary>
    public static VcfReader Open(string path, string indexPath, VcfReaderOptions? options = null)
    {
        var reader = new VcfReader(BioFile.OpenText(path), options) { _indexed = (path, indexPath) };
        return reader;
    }

    private (string Path, string Index)? _indexed;
    private TabixReader? _tabix;

    /// <summary>Whether region queries use a tabix index.</summary>
    public bool IsIndexed => _indexed is not null;

    /// <summary>Records overlapping a samtools-style region (<c>chr1:1000-2000</c>, 1-based inclusive, or <c>chr1</c>).</summary>
    public IEnumerable<VcfRecord> Query(string region)
    {
        var r = HtsRegion.Parse(region, n => _indexed is not null && (_tabix ??= TabixReader.Open(_indexed.Value.Path, _indexed.Value.Index)).Index.GetId(n) >= 0);
        return Query(new GenomicInterval(r.Name, r.Start, r.End));
    }

    /// <summary>Opens a VCF stream (gzip detected).</summary>
    public static VcfReader Open(Stream stream, VcfReaderOptions? options = null)
        => new(new StreamReader(BioFile.WrapDecompression(stream), System.Text.Encoding.UTF8, true, 1 << 16), options);

    /// <summary>Parses VCF text held in memory.</summary>
    public static (VcfHeader Header, List<VcfRecord> Records) Parse(string text, VcfReaderOptions? options = null)
    {
        using var r = new VcfReader(new StringReader(text), options);
        return (r.Header, r.Records().ToList());
    }

    /// <summary>Reads all records of a file lazily (the file stays open while enumerating).</summary>
    public static IEnumerable<VcfRecord> ReadRecords(string path, VcfReaderOptions? options = null)
    {
        using var r = Open(path, options);
        foreach (var rec in r.Records()) yield return rec;
    }

    /// <summary>The file header.</summary>
    public VcfHeader Header { get; }

    private VcfHeader ReadHeader()
    {
        var header = new VcfHeader();
        string? line;
        while ((line = _reader.ReadLine()) is not null)
        {
            _lineNumber++;
            if (line.Length == 0) continue;
            if (line.StartsWith("##", StringComparison.Ordinal)) { header.ParseMetaLine(line); continue; }
            if (line[0] == '#')
            {
                var cols = line.Split('\t');
                for (int i = 9; i < cols.Length; i++) header.Samples.Add(cols[i]);
                return header;
            }
            _pending = line; // header without #CHROM line
            return header;
        }
        return header;
    }

    /// <summary>Streams the data records.</summary>
    public IEnumerable<VcfRecord> Records()
    {
        string? line;
        if (_pending is not null) { line = _pending; _pending = null; yield return ParseLine(line); }
        while ((line = _reader.ReadLine()) is not null)
        {
            _lineNumber++;
            if (line.Length == 0 || line[0] == '#') continue;
            yield return ParseLine(line);
        }
    }

    /// <summary>Streams the variants (records without their extra FORMAT fields).</summary>
    public IEnumerable<Variant> Variants() => Records().Select(r => r.Variant);

    /// <summary>
    /// Streams the records whose reference span overlaps <paramref name="region"/> (a linear scan; lines on other
    /// chromosomes or outside the window are skipped before full parsing).
    /// </summary>
    public IEnumerable<VcfRecord> Query(GenomicInterval region)
    {
        if (_indexed is { } ix)
        {
            // tabix semantics: REF span, INFO END and SVLEN, exactly as the tabix command; repeatable
            _tabix ??= TabixReader.Open(ix.Path, ix.Index);
            foreach (var l in _tabix.Query(region.Chromosome, region.Start, region.End)) yield return ParseLine(l);
            yield break;
        }
        string? line;
        if (_pending is not null)
        {
            line = _pending; _pending = null;
            if (QuickOverlaps(line, region)) yield return ParseLine(line);
        }
        while ((line = _reader.ReadLine()) is not null)
        {
            _lineNumber++;
            if (line.Length == 0 || line[0] == '#') continue;
            if (QuickOverlaps(line, region)) yield return ParseLine(line);
        }
    }

    private static bool QuickOverlaps(string line, GenomicInterval region)
    {
        var s = line.AsSpan();
        int t1 = s.IndexOf('\t');
        if (t1 < 0 || !s[..t1].SequenceEqual(region.Chromosome)) return false;
        s = s[(t1 + 1)..];
        int t2 = s.IndexOf('\t');
        if (t2 < 0 || !long.TryParse(s[..t2], NumberStyles.Integer, CultureInfo.InvariantCulture, out long pos)) return false;
        s = s[(t2 + 1)..];
        int t3 = s.IndexOf('\t');
        if (t3 < 0) return false;
        s = s[(t3 + 1)..];
        int t4 = s.IndexOf('\t');
        long refLen = t4 < 0 ? 1 : Math.Max(1, t4);
        long start = pos - 1;
        return start < region.End && start + refLen > region.Start;
    }

    private string Intern(ReadOnlySpan<char> s)
    {
        var lookup = _strings.GetAlternateLookup<ReadOnlySpan<char>>();
        if (lookup.TryGetValue(s, out var existing)) return existing;
        string str = s.ToString();
        _strings[str] = str;
        return str;
    }

    private List<string> FormatKeys(ReadOnlySpan<char> s)
    {
        var lookup = _formats.GetAlternateLookup<ReadOnlySpan<char>>();
        if (lookup.TryGetValue(s, out var existing)) return existing;
        string str = s.ToString();
        var keys = str.Split(':').Select(k => Intern(k)).ToList();
        _formats[str] = keys;
        return keys;
    }

    private static ReadOnlySpan<char> NextField(ref ReadOnlySpan<char> s)
    {
        int tab = s.IndexOf('\t');
        ReadOnlySpan<char> field;
        if (tab < 0) { field = s; s = default; }
        else { field = s[..tab]; s = s[(tab + 1)..]; }
        return field;
    }

    private VcfRecord ParseLine(string line)
    {
        var s = line.AsSpan();
        var chrom = NextField(ref s);
        var posField = NextField(ref s);
        var id = NextField(ref s);
        var refField = NextField(ref s);
        var alt = NextField(ref s);
        var qual = NextField(ref s);
        var filter = NextField(ref s);
        var info = NextField(ref s);
        if (refField.IsEmpty || chrom.IsEmpty)
            throw new InvalidDataException($"VCF line {_lineNumber}: expected at least 8 tab-separated columns.");
        if (!long.TryParse(posField, NumberStyles.Integer, CultureInfo.InvariantCulture, out long pos))
            throw new InvalidDataException($"VCF line {_lineNumber}: invalid POS '{posField}'.");

        string[] alts;
        if (alt.IsEmpty || alt is ".") alts = [];
        else
        {
            alts = new string[alt.Count(',') + 1];
            for (int i = 0; i < alts.Length; i++)
            {
                int c = alt.IndexOf(',');
                var a = c < 0 ? alt : alt[..c];
                alts[i] = a.Length == 1 ? Intern(a) : a.ToString();
                if (c >= 0) alt = alt[(c + 1)..];
            }
        }

        double q = double.NaN;
        if (!(qual.IsEmpty || qual is "."))
            q = double.Parse(qual, NumberStyles.Float, CultureInfo.InvariantCulture);

        string refAllele = refField.Length <= 2 ? Intern(refField) : refField.ToString();
        var v = new Variant(Intern(chrom), pos, refAllele, alts, id.IsEmpty || id is "." ? null : id.ToString())
        {
            Quality = q,
            Filter = filter.IsEmpty ? "." : Intern(filter),
        };

        if (_options.ParseInfo && !(info.IsEmpty || info is "."))
        {
            while (!info.IsEmpty)
            {
                int sc = info.IndexOf(';');
                var item = sc < 0 ? info : info[..sc];
                info = sc < 0 ? default : info[(sc + 1)..];
                if (item.IsEmpty) continue;
                int eq = item.IndexOf('=');
                if (eq < 0) v.Info[Intern(item)] = "";
                else v.Info[Intern(item[..eq])] = item[(eq + 1)..].ToString();
            }
        }

        if (s.IsEmpty) return new VcfRecord(v);
        var formatField = NextField(ref s);
        var format = FormatKeys(formatField);
        if (s.IsEmpty) return new VcfRecord(v, format, null, 0);

        int sampleCount = s.Count('\t') + 1;
        int gtIndex = format.Count > 0 && format[0] == "GT" ? 0 : format.IndexOf("GT");
        if (_options.ParseGenotypes && gtIndex >= 0)
        {
            v.Genotypes.Capacity = sampleCount;
            var rest = s;
            while (true)
            {
                var col = NextField(ref rest);
                var gt = col;
                for (int f = 0; f < gtIndex; f++)
                {
                    int c = gt.IndexOf(':');
                    gt = c < 0 ? default : gt[(c + 1)..];
                }
                int colon = gt.IndexOf(':');
                if (colon >= 0) gt = gt[..colon];
                v.Genotypes.Add(ParseGenotype(gt));
                if (rest.IsEmpty) break;
            }
        }
        return new VcfRecord(v, format, s.ToString(), sampleCount);
    }

    /// <summary>Span-based GT parser (<c>0/1</c>, <c>1|0</c>, <c>./.</c>, <c>2</c>, <c>.</c>).</summary>
    internal static Genotype ParseGenotype(ReadOnlySpan<char> gt)
    {
        if (gt.IsEmpty) return new Genotype([-1]);
        int n = 1;
        bool phased = false;
        foreach (char c in gt)
        {
            if (c == '/') n++;
            else if (c == '|') { n++; phased = true; }
        }
        var alleles = new int[n];
        int k = 0, value = 0;
        bool missing = false, any = false;
        foreach (char c in gt)
        {
            if (c is '/' or '|')
            {
                alleles[k++] = missing || !any ? -1 : value;
                value = 0; missing = false; any = false;
            }
            else if (c == '.') missing = true;
            else if (c is >= '0' and <= '9') { value = value * 10 + (c - '0'); any = true; }
            else throw new InvalidDataException($"Invalid genotype '{gt}'.");
        }
        alleles[k] = missing || !any ? -1 : value;
        return new Genotype(alleles, phased);
    }

    /// <summary>Closes the underlying reader.</summary>
    public void Dispose()
    {
        _reader.Dispose();
        _tabix?.Dispose();
    }
}
