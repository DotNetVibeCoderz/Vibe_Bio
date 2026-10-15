using System.Text;

namespace Biotech.Net.Genomics.Vcf;

/// <summary>Definition of an INFO or FORMAT field from a VCF meta-line (<c>##INFO=&lt;ID=DP,Number=1,Type=Integer,Description="..."&gt;</c>).</summary>
/// <param name="Id">Field key.</param>
/// <param name="Number">Cardinality: an integer, <c>A</c> (one per alternate), <c>R</c> (one per allele), <c>G</c> (one per genotype) or <c>.</c>.</param>
/// <param name="Type">Integer, Float, Flag, Character or String.</param>
/// <param name="Description">Free-text description.</param>
public sealed record VcfFieldDefinition(string Id, string Number, string Type, string Description);

/// <summary>A <c>##FILTER</c> definition.</summary>
public sealed record VcfFilterDefinition(string Id, string Description);

/// <summary>A <c>##contig</c> definition; <paramref name="Length"/> is null when the header omits it.</summary>
public sealed record VcfContig(string Id, long? Length = null);

/// <summary>
/// The header of a VCF 4.x file: file format, structured meta-lines (INFO, FORMAT, FILTER, contig), any other
/// meta-lines verbatim, and the sample names from the <c>#CHROM</c> line.
/// </summary>
public sealed class VcfHeader
{
    private readonly List<VcfFieldDefinition> _info = [];
    private readonly List<VcfFieldDefinition> _format = [];
    private readonly List<VcfFilterDefinition> _filters = [];
    private readonly List<VcfContig> _contigs = [];

    /// <summary>The <c>##fileformat</c> value, e.g. <c>VCFv4.2</c>.</summary>
    public string FileFormat { get; set; } = "VCFv4.2";

    /// <summary>INFO definitions in header order.</summary>
    public IReadOnlyList<VcfFieldDefinition> Info => _info;

    /// <summary>FORMAT definitions in header order.</summary>
    public IReadOnlyList<VcfFieldDefinition> Format => _format;

    /// <summary>FILTER definitions in header order.</summary>
    public IReadOnlyList<VcfFilterDefinition> Filters => _filters;

    /// <summary>Contig definitions in header order.</summary>
    public IReadOnlyList<VcfContig> Contigs => _contigs;

    /// <summary>Other meta-lines (<c>##source=...</c>, <c>##reference=...</c>, ...) without the leading <c>##</c>.</summary>
    public List<string> OtherLines { get; } = [];

    /// <summary>Sample names (columns after FORMAT).</summary>
    public List<string> Samples { get; } = [];

    /// <summary>Looks up an INFO definition.</summary>
    public VcfFieldDefinition? GetInfo(string id) => _info.Find(d => d.Id == id);

    /// <summary>Looks up a FORMAT definition.</summary>
    public VcfFieldDefinition? GetFormat(string id) => _format.Find(d => d.Id == id);

    /// <summary>Index of a sample by name, or −1.</summary>
    public int SampleIndex(string name) => Samples.IndexOf(name);

    /// <summary>Adds or replaces an INFO definition.</summary>
    public VcfHeader AddInfo(string id, string number, string type, string description)
    {
        Replace(_info, new VcfFieldDefinition(id, number, type, description), d => d.Id == id);
        return this;
    }

    /// <summary>Adds or replaces a FORMAT definition.</summary>
    public VcfHeader AddFormat(string id, string number, string type, string description)
    {
        Replace(_format, new VcfFieldDefinition(id, number, type, description), d => d.Id == id);
        return this;
    }

    /// <summary>Adds or replaces a FILTER definition.</summary>
    public VcfHeader AddFilter(string id, string description)
    {
        Replace(_filters, new VcfFilterDefinition(id, description), d => d.Id == id);
        return this;
    }

    /// <summary>Adds or replaces a contig definition.</summary>
    public VcfHeader AddContig(string id, long? length = null)
    {
        Replace(_contigs, new VcfContig(id, length), d => d.Id == id);
        return this;
    }

    /// <summary>
    /// The standard FORMAT/INFO definitions used by <c>SnvCaller</c> output: GT, AD, DP, GQ, PL and INFO DP/AF/AC/AN.
    /// </summary>
    public VcfHeader AddStandardGenotypeFields()
    {
        AddInfo("DP", "1", "Integer", "Total read depth");
        AddInfo("AC", "A", "Integer", "Allele count in genotypes");
        AddInfo("AN", "1", "Integer", "Total number of alleles in called genotypes");
        AddInfo("AF", "A", "Float", "Allele frequency");
        AddFormat("GT", "1", "String", "Genotype");
        AddFormat("AD", "R", "Integer", "Allelic depths for the ref and alt alleles");
        AddFormat("DP", "1", "Integer", "Read depth");
        AddFormat("GQ", "1", "Integer", "Genotype quality");
        AddFormat("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
        return this;
    }

    private static void Replace<T>(List<T> list, T item, Predicate<T> match)
    {
        int i = list.FindIndex(match);
        if (i >= 0) list[i] = item; else list.Add(item);
    }

    /// <summary>Parses one meta-line (with or without the leading <c>##</c>).</summary>
    public void ParseMetaLine(string line)
    {
        var s = line.AsSpan();
        if (s.StartsWith("##")) s = s[2..];
        int eq = s.IndexOf('=');
        if (eq < 0) { OtherLines.Add(s.ToString()); return; }
        var key = s[..eq];
        var value = s[(eq + 1)..];
        if (key.SequenceEqual("fileformat")) { FileFormat = value.ToString(); return; }
        if (value.Length < 2 || value[0] != '<' || value[^1] != '>'
            || !(key is "INFO" or "FORMAT" or "FILTER" or "contig"))
        {
            OtherLines.Add(s.ToString());
            return;
        }
        var attrs = ParseStructured(value[1..^1]);
        string id = attrs.GetValueOrDefault("ID") ?? throw new InvalidDataException($"Meta-line without ID: ##{s}");
        string desc = attrs.GetValueOrDefault("Description") ?? "";
        switch (key)
        {
            case "INFO": AddInfo(id, attrs.GetValueOrDefault("Number") ?? ".", attrs.GetValueOrDefault("Type") ?? "String", desc); break;
            case "FORMAT": AddFormat(id, attrs.GetValueOrDefault("Number") ?? ".", attrs.GetValueOrDefault("Type") ?? "String", desc); break;
            case "FILTER": AddFilter(id, desc); break;
            default:
                AddContig(id, long.TryParse(attrs.GetValueOrDefault("length"), out long len) ? len : null);
                break;
        }
    }

    /// <summary>Parses <c>ID=x,Number=1,Description="a, b"</c> respecting quotes.</summary>
    internal static Dictionary<string, string> ParseStructured(ReadOnlySpan<char> body)
    {
        var result = new Dictionary<string, string>(StringComparer.Ordinal);
        int i = 0;
        while (i < body.Length)
        {
            int eq = body[i..].IndexOf('=');
            if (eq < 0) break;
            string k = body.Slice(i, eq).Trim().ToString();
            i += eq + 1;
            string v;
            if (i < body.Length && body[i] == '"')
            {
                var sb = new StringBuilder();
                i++;
                while (i < body.Length && body[i] != '"')
                {
                    if (body[i] == '\\' && i + 1 < body.Length) i++;
                    sb.Append(body[i++]);
                }
                i++; // closing quote
                v = sb.ToString();
                int comma = body[Math.Min(i, body.Length)..].IndexOf(',');
                i = comma < 0 ? body.Length : i + comma + 1;
            }
            else
            {
                int comma = body[i..].IndexOf(',');
                int end = comma < 0 ? body.Length : i + comma;
                v = body[i..end].ToString();
                i = end + 1;
            }
            result[k] = v;
        }
        return result;
    }

    /// <summary>Writes the header (meta-lines and the <c>#CHROM</c> line).</summary>
    public void Write(TextWriter w)
    {
        w.Write("##fileformat="); w.Write(FileFormat); w.Write('\n');
        foreach (var l in OtherLines) { w.Write("##"); w.Write(l); w.Write('\n'); }
        foreach (var f in _filters) w.Write($"##FILTER=<ID={f.Id},Description=\"{Escape(f.Description)}\">\n");
        foreach (var d in _info) w.Write($"##INFO=<ID={d.Id},Number={d.Number},Type={d.Type},Description=\"{Escape(d.Description)}\">\n");
        foreach (var d in _format) w.Write($"##FORMAT=<ID={d.Id},Number={d.Number},Type={d.Type},Description=\"{Escape(d.Description)}\">\n");
        foreach (var c in _contigs) w.Write(c.Length is { } len ? $"##contig=<ID={c.Id},length={len}>\n" : $"##contig=<ID={c.Id}>\n");
        w.Write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if (Samples.Count > 0)
        {
            w.Write("\tFORMAT");
            foreach (var s in Samples) { w.Write('\t'); w.Write(s); }
        }
        w.Write('\n');
    }

    private static string Escape(string s) => s.Replace("\\", "\\\\").Replace("\"", "\\\"");

    /// <summary>The header as text.</summary>
    public override string ToString()
    {
        var sw = new StringWriter { NewLine = "\n" };
        Write(sw);
        return sw.ToString();
    }
}
