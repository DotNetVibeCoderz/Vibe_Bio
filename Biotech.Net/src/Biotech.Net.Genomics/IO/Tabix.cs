using System.Buffers.Binary;
using System.Globalization;
using System.IO.Compression;
using System.Text;

namespace Biotech.Net.Genomics.IO;

/// <summary>
/// How tabix finds the interval of a line: column numbers are 1-based; <see cref="EndColumn"/> 0 means the end is
/// derived from the format (VCF: REF length, INFO END/SVLEN; SAM: CIGAR).
/// </summary>
/// <param name="Format">0 generic, 1 SAM, 2 VCF, plus 0x10000 when coordinates are 0-based half-open (UCSC/BED).</param>
/// <param name="SequenceColumn">Column of the sequence name.</param>
/// <param name="BeginColumn">Column of the start.</param>
/// <param name="EndColumn">Column of the end (0: derived).</param>
/// <param name="MetaChar">Lines starting with this character are headers/comments.</param>
/// <param name="SkipLines">Leading lines to skip.</param>
public sealed record TabixConfig(int Format, int SequenceColumn, int BeginColumn, int EndColumn, char MetaChar = '#', int SkipLines = 0)
{
    /// <summary>Format code for generic tab-separated files.</summary>
    public const int Generic = 0;
    /// <summary>Format code for SAM.</summary>
    public const int Sam = 1;
    /// <summary>Format code for VCF.</summary>
    public const int Vcf = 2;
    /// <summary>Flag for 0-based half-open coordinates.</summary>
    public const int Ucsc = 0x10000;

    /// <summary>The <c>tabix -p vcf</c> preset.</summary>
    public static TabixConfig VcfPreset { get; } = new(Vcf, 1, 2, 0, '#', 0);
    /// <summary>The <c>tabix -p bed</c> preset.</summary>
    public static TabixConfig BedPreset { get; } = new(Ucsc, 1, 2, 3, '#', 0);
    /// <summary>The <c>tabix -p gff</c> preset (also GTF).</summary>
    public static TabixConfig GffPreset { get; } = new(Generic, 1, 4, 5, '#', 0);
    /// <summary>The <c>tabix -p sam</c> preset.</summary>
    public static TabixConfig SamPreset { get; } = new(Sam, 3, 4, 0, '@', 0);

    /// <summary>The format without the UCSC flag.</summary>
    public int Preset => Format & 0xFFFF;

    /// <summary>Whether coordinates are 0-based half-open.</summary>
    public bool ZeroBased => (Format & Ucsc) != 0;

    /// <summary>A preset by name (<c>vcf</c>, <c>bed</c>, <c>gff</c>, <c>gtf</c>, <c>sam</c>).</summary>
    public static TabixConfig FromPreset(string name) => name.ToLowerInvariant() switch
    {
        "vcf" => VcfPreset,
        "bed" => BedPreset,
        "gff" or "gtf" or "gff3" => GffPreset,
        "sam" => SamPreset,
        _ => throw new ArgumentException($"Unknown tabix preset '{name}'.", nameof(name)),
    };

    /// <summary>
    /// Parses the sequence name and 0-based half-open interval of a data line exactly as htslib's <c>tbx_parse1</c>.
    /// Returns false when the line has no valid interval.
    /// </summary>
    public bool TryParse(ReadOnlySpan<byte> line, out ReadOnlySpan<byte> name, out long beg, out long end)
    {
        name = default;
        beg = end = -1;
        bool haveName = false;
        int id = 1, b = 0, alcnt = 0, lenpos = -1;
        bool getlen = false, useSvlen = false;
        long reflen = 0, svlen = 0, fmtlen = 0;
        Span<bool> svAllele = stackalloc bool[64];
        int preset = Preset;
        for (int i = 0; i <= line.Length; i++)
        {
            if (i < line.Length && line[i] != (byte)'\t') continue;
            var field = line[b..i];
            if (id == SequenceColumn)
            {
                name = field;
                haveName = true;
            }
            else if (id == BeginColumn)
            {
                if (!TryLong(field, out beg)) return false;
                if (BeginColumn <= EndColumn) end = beg;
                if (!ZeroBased) --beg;
                else if (BeginColumn <= EndColumn) ++end;
                if (beg < 0) beg = 0;
                if (end < 1) end = 1;
            }
            else if (preset == Generic)
            {
                if (id == EndColumn && !TryLong(field, out end)) return false;
            }
            else if (preset == Sam)
            {
                if (id == 6)
                {
                    long l = 0, x = 0;
                    foreach (byte c in field)
                    {
                        if (c is >= (byte)'0' and <= (byte)'9') { x = x * 10 + (c - '0'); continue; }
                        char op = char.ToUpperInvariant((char)c);
                        if (op is 'M' or 'D' or 'N') l += x;
                        x = 0;
                    }
                    if (l == 0) l = 1;
                    end = beg + l;
                }
            }
            else if (preset == Vcf)
            {
                if (id == 4)
                {
                    if (field.Length > 0) end = beg + field.Length;
                    ++alcnt;
                    reflen = field.Length;
                }
                if (id == 5)
                {
                    var rest = field;
                    while (true)
                    {
                        int comma = rest.IndexOf((byte)',');
                        var alt = comma < 0 ? rest : rest[..comma];
                        ++alcnt;
                        if (SvlenOnRef(alt))
                        {
                            if (alcnt - 1 < svAllele.Length) svAllele[alcnt - 1] = true;
                            useSvlen = true;
                        }
                        else if (alt.SequenceEqual("<*>"u8) || alt.SequenceEqual("<NON_REF>"u8)) getlen = true;
                        if (comma < 0 || alcnt >= 65536) break;
                        rest = rest[(comma + 1)..];
                    }
                }
                else if (id == 8)
                {
                    if (TryFindInfo(field, "END="u8, out var ev) && !(ev.Length > 0 && ev[0] == (byte)'.'))
                    {
                        TryLeadingLong(ev, out long e);
                        if (e > beg) end = e;
                    }
                    if (TryFindInfo(field, "SVLEN="u8, out var rest))
                    {
                        for (int d = 1; d < alcnt; d++)
                        {
                            long tmp = 1;
                            if (useSvlen && d < svAllele.Length && svAllele[d])
                            {
                                TryLeadingLong(rest, out tmp);
                                tmp = Math.Abs(tmp);
                            }
                            svlen = Math.Max(svlen, tmp);
                            int comma = rest.IndexOf((byte)',');
                            if (comma < 0) break;
                            rest = rest[(comma + 1)..];
                        }
                    }
                }
                else if (getlen && id == 9)
                {
                    int pos = 0;
                    var rest = field;
                    while (true)
                    {
                        int colon = rest.IndexOf((byte)':');
                        var key = colon < 0 ? rest : rest[..colon];
                        if (key.SequenceEqual("LEN"u8)) { lenpos = pos; break; }
                        if (colon < 0) break;
                        rest = rest[(colon + 1)..];
                        ++pos;
                    }
                    if (lenpos == -1) break; // htslib stops scanning the line here
                }
                else if (id > 9 && getlen && lenpos != -1)
                {
                    var rest = field;
                    long tmp = 0;
                    for (int k = 0; k <= lenpos; k++)
                    {
                        if (k == lenpos) { TryLeadingLong(rest, out tmp); break; }
                        int colon = rest.IndexOf((byte)':');
                        if (colon < 0) break;
                        rest = rest[(colon + 1)..];
                    }
                    fmtlen = Math.Max(fmtlen, tmp);
                }
            }
            b = i + 1;
            ++id;
        }
        if (preset == Vcf)
        {
            long tmp = Math.Max(reflen, Math.Max(svlen, fmtlen)) + beg;
            end = Math.Max(end, tmp);
        }
        return haveName && beg >= 0 && end >= 0;
    }

    /// <summary>htslib's lookup of an INFO key: at the start of the field or after a <c>;</c>.</summary>
    private static bool TryFindInfo(ReadOnlySpan<byte> info, ReadOnlySpan<byte> key, out ReadOnlySpan<byte> value)
    {
        value = default;
        int p = info.IndexOf(key);
        if (p < 0) return false;
        if (p == 0) { value = info[key.Length..]; return true; }
        Span<byte> semi = stackalloc byte[key.Length + 1];
        semi[0] = (byte)';';
        key.CopyTo(semi[1..]);
        p = info.IndexOf(semi);
        if (p < 0) return false;
        value = info[(p + semi.Length)..];
        return true;
    }

    private static bool SvlenOnRef(ReadOnlySpan<byte> alt)
    {
        if (alt.Length < 5 || alt[0] != (byte)'<') return false;
        if (alt[4] != (byte)'>' && alt[4] != (byte)':') return false;
        var k = alt[..4];
        if (!k.SequenceEqual("<CNV"u8) && !k.SequenceEqual("<DEL"u8) && !k.SequenceEqual("<DUP"u8) && !k.SequenceEqual("<INV"u8)) return false;
        return alt[^1] == (byte)'>';
    }

    /// <summary>strtoll-like: optional sign and leading digits; false when no digit.</summary>
    private static bool TryLeadingLong(ReadOnlySpan<byte> s, out long v)
    {
        v = 0;
        int i = 0;
        bool neg = false;
        while (i < s.Length && (s[i] == (byte)' ')) i++;
        if (i < s.Length && (s[i] == (byte)'-' || s[i] == (byte)'+')) { neg = s[i] == (byte)'-'; i++; }
        int start = i;
        while (i < s.Length && s[i] is >= (byte)'0' and <= (byte)'9') v = v * 10 + (s[i++] - '0');
        if (neg) v = -v;
        return i > start;
    }

    private static bool TryLong(ReadOnlySpan<byte> s, out long v) => TryLeadingLong(s, out v);
}

/// <summary>
/// A tabix index (<c>.tbi</c>, or tabix-flavoured <c>.csi</c>) over a BGZF-compressed, position-sorted text file:
/// the column configuration, the sequence names and the binning index.
/// </summary>
public sealed class TabixIndex
{
    private readonly Dictionary<string, int> _ids;

    internal TabixIndex(TabixConfig config, List<string> names, BinningIndex index)
    {
        Config = config;
        Names = names;
        Index = index;
        _ids = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < names.Count; i++) _ids.TryAdd(names[i], i);
    }

    /// <summary>Column configuration.</summary>
    public TabixConfig Config { get; }

    /// <summary>Sequence names in index order.</summary>
    public IReadOnlyList<string> Names { get; }

    /// <summary>Bins, chunks and linear offsets.</summary>
    public BinningIndex Index { get; }

    /// <summary>The index of sequence <paramref name="name"/>, or -1.</summary>
    public int GetId(string name) => _ids.TryGetValue(name, out int i) ? i : -1;

    /// <summary>Loads a <c>.tbi</c> or <c>.csi</c> tabix index.</summary>
    public static TabixIndex Load(string path)
    {
        using var fs = File.OpenRead(path);
        return Read(fs);
    }

    /// <summary>Reads a TBI or tabix-CSI index (both BGZF-compressed).</summary>
    public static TabixIndex Read(Stream stream)
    {
        var raw = BinningIndex.ReadAll(stream);
        var data = raw.Length >= 2 && raw[0] == 0x1f && raw[1] == 0x8b ? BinningIndex.Decompress(raw) : raw;
        if (data.Length >= 4 && data.AsSpan(0, 4).SequenceEqual("TBI\u0001"u8))
        {
            var r = new SpanReader(data);
            r.Bytes(4);
            int nRef = r.Int32();
            var (config, names) = ReadMeta(ref r);
            var idx = BinningIndex.ReadBaiBody(ref r, BinningIndexFormat.Tbi, nRef);
            return new TabixIndex(config, names, idx);
        }
        if (data.Length >= 4 && data.AsSpan(0, 4).SequenceEqual("CSI\u0001"u8))
        {
            var idx = BinningIndex.Read(new MemoryStream(data));
            if (idx.Aux.Length < 28) throw new InvalidDataException("CSI index has no tabix header (it indexes BAM, not a text file).");
            var r = new SpanReader(idx.Aux);
            var (config, names) = ReadMeta(ref r);
            return new TabixIndex(config, names, idx);
        }
        throw new InvalidDataException("Not a tabix index (TBI or CSI magic expected).");
    }

    private static (TabixConfig, List<string>) ReadMeta(ref SpanReader r)
    {
        int format = r.Int32(), sc = r.Int32(), bc = r.Int32(), ec = r.Int32(), meta = r.Int32(), skip = r.Int32();
        int lNm = r.Int32();
        var nm = r.Bytes(lNm);
        var names = new List<string>();
        int start = 0;
        for (int i = 0; i < nm.Length; i++)
            if (nm[i] == 0) { names.Add(Encoding.UTF8.GetString(nm[start..i])); start = i + 1; }
        return (new TabixConfig(format, sc, bc, ec, (char)meta, skip), names);
    }

    private byte[] MetaBytes()
    {
        using var ms = new MemoryStream();
        var w = new BinaryWriter(ms);
        w.Write(Config.Format); w.Write(Config.SequenceColumn); w.Write(Config.BeginColumn); w.Write(Config.EndColumn);
        w.Write((int)Config.MetaChar); w.Write(Config.SkipLines);
        var nm = new MemoryStream();
        foreach (var n in Names) { nm.Write(Encoding.UTF8.GetBytes(n)); nm.WriteByte(0); }
        w.Write((int)nm.Length);
        w.Write(nm.ToArray());
        w.Flush();
        return ms.ToArray();
    }

    /// <summary>Saves as TBI or CSI according to <see cref="BinningIndex.Format"/>.</summary>
    public void Save(string path)
    {
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None);
        Write(fs);
    }

    /// <summary>Writes as TBI or CSI according to <see cref="BinningIndex.Format"/>.</summary>
    public void Write(Stream stream)
    {
        if (Index.Format == BinningIndexFormat.Csi)
        {
            Index.Aux = MetaBytes();
            Index.Write(stream);
            return;
        }
        using var bg = new BgzfWriter(stream, CompressionLevel.Optimal, leaveOpen: true);
        var w = new BinaryWriter(bg, Encoding.ASCII, leaveOpen: true);
        w.Write("TBI\u0001"u8);
        bool vcf = Config.Preset == TabixConfig.Vcf;
        int nids = vcf ? Index.References.Count(r => r.Bins.Count > 0 || r.HasMetadata) : Index.References.Count;
        w.Write(nids);
        w.Write(MetaBytes());
        Index.WriteBaiBody(w, skipEmptyReferences: vcf);
        w.Flush();
    }

    /// <summary>
    /// Indexes a BGZF-compressed, position-sorted file (like <c>tabix -p</c>). <paramref name="format"/> TBI uses
    /// 14/5 bins; CSI uses <paramref name="minShift"/> with htslib's choice of depth.
    /// </summary>
    public static TabixIndex Build(string bgzfPath, TabixConfig config, BinningIndexFormat format = BinningIndexFormat.Tbi, int minShift = 14)
    {
        using var bg = BgzfReader.Open(bgzfPath);
        return Build(bg, config, format, minShift);
    }

    /// <summary>Indexes BGZF data read from the start of <paramref name="bgzf"/>.</summary>
    public static TabixIndex Build(BgzfReader bgzf, TabixConfig config, BinningIndexFormat format = BinningIndexFormat.Tbi, int minShift = 14)
    {
        if (format == BinningIndexFormat.Bai) throw new ArgumentException("Tabix uses TBI or CSI.", nameof(format));
        var names = new List<string>();
        var ids = new Dictionary<string, int>(StringComparer.Ordinal);
        var line = new ByteBuffer(1024);
        long lineNo = 0, maxRefLen = 0;
        var lastOff = bgzf.VirtualPosition;
        BinningIndexBuilder? builder = null;
        int depth = 5;
        if (format == BinningIndexFormat.Csi) depth = (31 - minShift + 2) / 3; // TBX_MAX_SHIFT = 31
        while (bgzf.ReadLine(line))
        {
            ++lineNo;
            var s = line.Span;
            if (format == BinningIndexFormat.Csi && s.Length > 0 && s[0] == config.MetaChar) maxRefLen = Math.Max(maxRefLen, HeaderLength(s, config.Preset));
            if (lineNo <= config.SkipLines || (s.Length > 0 && s[0] == config.MetaChar))
            {
                lastOff = bgzf.VirtualPosition;
                continue;
            }
            if (builder is null)
            {
                if (format == BinningIndexFormat.Csi)
                {
                    if (maxRefLen > 0) AdjustCsi(maxRefLen, ref minShift, ref depth);
                    else depth = minShift < 10 ? 9 : minShift < 25 ? 9 - (minShift - 10) / 3 : 4;
                    builder = new BinningIndexBuilder(format, lastOff, 0, minShift, depth);
                }
                else builder = new BinningIndexBuilder(format, lastOff);
            }
            if (!config.TryParse(s, out var nameBytes, out long beg, out long end))
                throw new InvalidDataException($"Line {lineNo}: cannot parse a {(config.Preset == TabixConfig.Vcf ? "VCF" : "tab-separated")} interval: \"{line}\"");
            string name = Encoding.UTF8.GetString(nameBytes);
            if (!ids.TryGetValue(name, out int tid)) { tid = names.Count; ids[name] = tid; names.Add(name); }
            builder.Push(tid, beg, end, bgzf.VirtualPosition, true);
        }
        builder ??= format == BinningIndexFormat.Csi ? new BinningIndexBuilder(format, lastOff, 0, minShift, depth) : new BinningIndexBuilder(format, lastOff);
        var idx = builder.Finish(bgzf.VirtualPosition);
        return new TabixIndex(config, names, idx);
    }

    private static void AdjustCsi(long maxRefLen, ref int minShift, ref int depth)
    {
        const int maxDepth = 9;
        long maxLen = maxRefLen + 256;
        if (maxLen <= 1L << (minShift + 3 * maxDepth))
        {
            long maxPos = 1L << (minShift + 3 * depth);
            while (maxLen > maxPos) { ++depth; maxPos *= 8; }
        }
        else
        {
            depth = maxDepth;
            long maxPos = 1L << (minShift + 3 * depth);
            while (maxLen > maxPos) { ++minShift; maxPos *= 2; }
        }
    }

    private static long HeaderLength(ReadOnlySpan<byte> s, int preset)
    {
        // ##contig=<...length=N> (VCF) or @SQ ... LN:N (SAM): htslib uses these to size CSI bins
        if (preset == TabixConfig.Vcf && s.StartsWith("##contig"u8))
        {
            int p = s[8..].IndexOf("length"u8);
            if (p < 0) return 0;
            var t = s[(8 + p + 6)..];
            int i = 0;
            while (i < t.Length && (t[i] == (byte)' ' || t[i] == (byte)'=')) i++;
            long v = 0;
            while (i < t.Length && t[i] is >= (byte)'0' and <= (byte)'9') v = v * 10 + (t[i++] - '0');
            return v;
        }
        if (preset == TabixConfig.Sam && s.StartsWith("@SQ"u8))
        {
            int p = s.IndexOf("\tLN:"u8);
            if (p < 0) return 0;
            long v = 0;
            for (int i = p + 4; i < s.Length && s[i] is >= (byte)'0' and <= (byte)'9'; i++) v = v * 10 + (s[i] - '0');
            return v;
        }
        return 0;
    }
}

/// <summary>
/// Region queries on a BGZF-compressed, tabix-indexed text file (VCF, BED, GFF, ...), returning the matching lines
/// exactly as the <c>tabix</c> command does. One enumeration at a time.
/// </summary>
public sealed class TabixReader : IDisposable
{
    private readonly BgzfReader _bgzf;
    private readonly ByteBuffer _line = new(1024);

    /// <summary>Opens a reader over BGZF data and its index.</summary>
    public TabixReader(BgzfReader bgzf, TabixIndex index)
    {
        _bgzf = bgzf ?? throw new ArgumentNullException(nameof(bgzf));
        Index = index ?? throw new ArgumentNullException(nameof(index));
    }

    /// <summary>Opens <paramref name="path"/> with <paramref name="indexPath"/> (default <c>path.tbi</c>, else <c>path.csi</c>).</summary>
    public static TabixReader Open(string path, string? indexPath = null)
    {
        indexPath ??= FindIndex(path) ?? throw new FileNotFoundException($"No tabix index (.tbi or .csi) found for '{path}'.", path + ".tbi");
        return new TabixReader(BgzfReader.Open(path), TabixIndex.Load(indexPath));
    }

    /// <summary>The index file next to <paramref name="path"/> (<c>.tbi</c> then <c>.csi</c>), or null.</summary>
    public static string? FindIndex(string path)
    {
        if (File.Exists(path + ".tbi")) return path + ".tbi";
        if (File.Exists(path + ".csi")) return path + ".csi";
        return null;
    }

    /// <summary>The index in use.</summary>
    public TabixIndex Index { get; }

    /// <summary>The leading header/comment lines (meta character or skipped lines).</summary>
    public IEnumerable<string> HeaderLines()
    {
        _bgzf.Seek(new VirtualOffset(0));
        long n = 0;
        while (_bgzf.ReadLine(_line))
        {
            ++n;
            var s = _line.Span;
            if (n <= Index.Config.SkipLines || (s.Length > 0 && s[0] == Index.Config.MetaChar)) yield return Encoding.UTF8.GetString(s);
            else yield break;
        }
    }

    /// <summary>Lines whose interval overlaps 0-based <c>[start, end)</c> on <paramref name="name"/>.</summary>
    public IEnumerable<string> Query(string name, long start, long end)
    {
        int tid = Index.GetId(name);
        if (tid < 0) yield break;
        if (start < 0) start = 0;
        foreach (var chunk in Index.Index.Query(tid, start, end))
        {
            if (_bgzf.VirtualPosition != chunk.Start) _bgzf.Seek(chunk.Start);
            while (_bgzf.VirtualPosition < chunk.End)
            {
                if (!_bgzf.ReadLine(_line)) yield break;
                var s = _line.Span;
                if (s.Length > 0 && s[0] == Index.Config.MetaChar) continue;
                if (!Index.Config.TryParse(s, out var nm, out long b, out long e))
                    throw new InvalidDataException($"Cannot parse an interval from line \"{_line}\".");
                if (!NameEquals(nm, Index.Names[tid]) || b >= end) yield break;
                if (e > start && end > b) yield return Encoding.UTF8.GetString(s);
            }
        }
    }

    /// <summary>Lines overlapping a samtools-style region (<c>chr1</c>, <c>chr1:100-200</c>, 1-based inclusive).</summary>
    public IEnumerable<string> Query(string region)
    {
        var r = HtsRegion.Parse(region, n => Index.GetId(n) >= 0);
        return Query(r.Name, r.Start, r.End);
    }

    /// <summary>All data lines in file order.</summary>
    public IEnumerable<string> Lines()
    {
        _bgzf.Seek(new VirtualOffset(0));
        long n = 0;
        while (_bgzf.ReadLine(_line))
        {
            ++n;
            var s = _line.Span;
            if (n <= Index.Config.SkipLines || (s.Length > 0 && s[0] == Index.Config.MetaChar)) continue;
            yield return Encoding.UTF8.GetString(s);
        }
    }

    private static bool NameEquals(ReadOnlySpan<byte> a, string b)
    {
        if (a.Length != b.Length) return Encoding.UTF8.GetString(a) == b;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) return false;
        return true;
    }

    /// <summary>Closes the file.</summary>
    public void Dispose() => _bgzf.Dispose();
}
