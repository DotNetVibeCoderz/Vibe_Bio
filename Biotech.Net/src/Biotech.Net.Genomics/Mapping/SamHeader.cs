using System.Globalization;
using System.Text;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>A reference sequence of a SAM/BAM/CRAM header (<c>@SQ</c> line).</summary>
/// <param name="Name">Reference name (<c>SN</c>).</param>
/// <param name="Length">Reference length (<c>LN</c>).</param>
public sealed record SamReference(string Name, long Length);

/// <summary>One header line split into its record type and <c>TAG:value</c> fields (<c>@CO</c> keeps its text as the single field <c>""</c>).</summary>
public sealed record SamHeaderLine(string Type, IReadOnlyList<KeyValuePair<string, string>> Fields)
{
    /// <summary>The value of field <paramref name="tag"/>, or null.</summary>
    public string? this[string tag]
    {
        get
        {
            foreach (var kv in Fields) if (kv.Key == tag) return kv.Value;
            return null;
        }
    }

    public override string ToString()
    {
        var sb = new StringBuilder("@").Append(Type);
        foreach (var kv in Fields) sb.Append('\t').Append(kv.Key.Length > 0 ? kv.Key + ":" : "").Append(kv.Value);
        return sb.ToString();
    }
}

/// <summary>
/// A SAM/BAM/CRAM header: the header text plus the reference dictionary. For BAM the dictionary comes from the
/// binary reference list (which wins over <c>@SQ</c> lines, as in htslib); for SAM and CRAM it is parsed from the
/// <c>@SQ</c> lines.
/// </summary>
public sealed class SamHeader
{
    private readonly List<SamReference> _references;
    private readonly Dictionary<string, int> _index = new(StringComparer.Ordinal);
    private List<SamHeaderLine>? _lines;

    /// <summary>Creates a header from its text; when <paramref name="references"/> is null they are read from the <c>@SQ</c> lines.</summary>
    public SamHeader(string text, IEnumerable<SamReference>? references = null)
    {
        Text = text ?? "";
        _references = references?.ToList() ?? ParseReferences(Text);
        for (int i = 0; i < _references.Count; i++) _index.TryAdd(_references[i].Name, i);
    }

    /// <summary>An empty header (no references).</summary>
    public static SamHeader Empty { get; } = new("");

    /// <summary>Builds a header with <c>@HD</c> and one <c>@SQ</c> per reference.</summary>
    public static SamHeader Create(IEnumerable<(string Name, long Length)> references, string sortOrder = "unsorted", IEnumerable<string>? extraLines = null)
    {
        var sb = new StringBuilder();
        sb.Append("@HD\tVN:1.6\tSO:").Append(sortOrder).Append('\n');
        var refs = new List<SamReference>();
        foreach (var (name, len) in references)
        {
            sb.Append("@SQ\tSN:").Append(name).Append("\tLN:").Append(len.ToString(CultureInfo.InvariantCulture)).Append('\n');
            refs.Add(new SamReference(name, len));
        }
        if (extraLines is not null) foreach (var l in extraLines) sb.Append(l.TrimEnd('\n')).Append('\n');
        return new SamHeader(sb.ToString(), refs);
    }

    /// <summary>Builds a header listing <paramref name="references"/> (e.g. the reference genome given to a mapper).</summary>
    public static SamHeader Create(IEnumerable<Core.Sequence> references, string sortOrder = "unsorted")
        => Create(references.Select(s => (s.Id, (long)s.Length)), sortOrder);

    /// <summary>The header text (lines separated by <c>\n</c>).</summary>
    public string Text { get; }

    /// <summary>The reference dictionary, in order (record reference ids index into it).</summary>
    public IReadOnlyList<SamReference> References => _references;

    /// <summary>The index of reference <paramref name="name"/>, or -1.</summary>
    public int GetReferenceId(string name) => name is null or "*" ? -1 : _index.TryGetValue(name, out int i) ? i : -1;

    /// <summary>The name of reference <paramref name="id"/>, or <c>*</c> for -1.</summary>
    public string GetReferenceName(int id) => id >= 0 && id < _references.Count ? _references[id].Name : "*";

    /// <summary>The parsed header lines.</summary>
    public IReadOnlyList<SamHeaderLine> Lines => _lines ??= ParseLines(Text);

    /// <summary>Lines of one record type (<c>"RG"</c>, <c>"PG"</c>, ...).</summary>
    public IEnumerable<SamHeaderLine> GetLines(string type) => Lines.Where(l => l.Type == type);

    /// <summary>The <c>@HD SO</c> sort order, or null.</summary>
    public string? SortOrder => GetLines("HD").FirstOrDefault()?["SO"];

    /// <summary>The <c>@HD VN</c> format version, or null.</summary>
    public string? Version => GetLines("HD").FirstOrDefault()?["VN"];

    /// <summary>Read-group identifiers (<c>@RG ID</c>) in header order.</summary>
    public IReadOnlyList<string> ReadGroupIds => GetLines("RG").Select(l => l["ID"]).OfType<string>().ToList();

    /// <summary>A copy with the <c>@HD SO</c> field set to <paramref name="sortOrder"/> (an <c>@HD</c> line is added if missing).</summary>
    public SamHeader WithSortOrder(string sortOrder)
    {
        var lines = Text.Split('\n').Where(l => l.Length > 0).ToList();
        int hd = lines.FindIndex(l => l.StartsWith("@HD", StringComparison.Ordinal));
        if (hd < 0) lines.Insert(0, "@HD\tVN:1.6\tSO:" + sortOrder);
        else
        {
            var f = lines[hd].Split('\t').ToList();
            int so = f.FindIndex(x => x.StartsWith("SO:", StringComparison.Ordinal));
            if (so < 0) f.Add("SO:" + sortOrder); else f[so] = "SO:" + sortOrder;
            lines[hd] = string.Join('\t', f);
        }
        return new SamHeader(string.Join('\n', lines) + "\n", _references);
    }

    private static List<SamReference> ParseReferences(string text)
    {
        var refs = new List<SamReference>();
        foreach (var raw in text.Split('\n'))
        {
            var line = raw.TrimEnd('\r');
            if (!line.StartsWith("@SQ\t", StringComparison.Ordinal)) continue;
            string? name = null;
            long len = 0;
            foreach (var f in line.Split('\t'))
            {
                if (f.StartsWith("SN:", StringComparison.Ordinal)) name = f[3..];
                else if (f.StartsWith("LN:", StringComparison.Ordinal)) len = long.Parse(f.AsSpan(3), NumberStyles.Integer, CultureInfo.InvariantCulture);
            }
            if (name is not null) refs.Add(new SamReference(name, len));
        }
        return refs;
    }

    private static List<SamHeaderLine> ParseLines(string text)
    {
        var lines = new List<SamHeaderLine>();
        foreach (var raw in text.Split('\n'))
        {
            var line = raw.TrimEnd('\r');
            if (line.Length < 3 || line[0] != '@') continue;
            int tab = line.IndexOf('\t');
            string type = tab < 0 ? line[1..] : line[1..tab];
            var fields = new List<KeyValuePair<string, string>>();
            if (tab >= 0)
            {
                if (type == "CO") fields.Add(new("", line[(tab + 1)..]));
                else
                    foreach (var f in line[(tab + 1)..].Split('\t'))
                    {
                        int c = f.IndexOf(':');
                        fields.Add(c < 0 ? new(f, "") : new(f[..c], f[(c + 1)..]));
                    }
            }
            lines.Add(new SamHeaderLine(type, fields));
        }
        return lines;
    }

    /// <summary>The header text.</summary>
    public override string ToString() => Text;
}
