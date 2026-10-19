using System.Collections;
using System.Text;
using System.Text.RegularExpressions;

namespace Biotech.Net.Workflows;

/// <summary>
/// The wildcard values of one job, such as <c>sample=A</c>, in the order the rule's output pattern names them.
/// </summary>
public sealed class Wildcards : IReadOnlyDictionary<string, string>
{
    private readonly List<KeyValuePair<string, string>> _items;

    /// <summary>No wildcards.</summary>
    public static Wildcards Empty { get; } = new([]);

    /// <summary>Creates wildcards from name/value pairs (order is kept).</summary>
    public Wildcards(IEnumerable<KeyValuePair<string, string>> values)
    {
        ArgumentNullException.ThrowIfNull(values);
        _items = [];
        foreach (var kv in values)
        {
            if (_items.Exists(i => i.Key == kv.Key)) throw new ArgumentException($"Wildcard '{kv.Key}' is given twice.", nameof(values));
            _items.Add(kv);
        }
    }

    /// <summary>Creates wildcards from <c>("sample", "A")</c> tuples.</summary>
    public static Wildcards Of(params (string Name, string Value)[] values)
        => new(values.Select(v => new KeyValuePair<string, string>(v.Name, v.Value)));

    /// <summary>The value of wildcard <paramref name="name"/>.</summary>
    /// <exception cref="KeyNotFoundException">The job has no such wildcard.</exception>
    public string this[string name] => TryGetValue(name, out var v)
        ? v
        : throw new KeyNotFoundException($"No wildcard '{name}'; this job has {(Count == 0 ? "none" : string.Join(", ", Keys))}.");

    /// <inheritdoc />
    public int Count => _items.Count;

    /// <inheritdoc />
    public IEnumerable<string> Keys => _items.Select(i => i.Key);

    /// <inheritdoc />
    public IEnumerable<string> Values => _items.Select(i => i.Value);

    /// <inheritdoc />
    public bool ContainsKey(string key) => _items.Exists(i => i.Key == key);

    /// <inheritdoc />
    public bool TryGetValue(string key, out string value)
    {
        foreach (var kv in _items)
        {
            if (kv.Key == key)
            {
                value = kv.Value;
                return true;
            }
        }

        value = "";
        return false;
    }

    /// <inheritdoc />
    public IEnumerator<KeyValuePair<string, string>> GetEnumerator() => _items.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <summary><c>sample=A,chunk=1</c>.</summary>
    public override string ToString() => string.Join(",", _items.Select(i => $"{i.Key}={i.Value}"));

    internal bool SameAs(Wildcards other)
        => Count == other.Count && _items.All(i => other.TryGetValue(i.Key, out var v) && v == i.Value);
}

/// <summary>
/// A path or value-key pattern with <c>{name}</c> wildcards, Snakemake style: <c>work/{sample}.trimmed.fastq</c>.
/// </summary>
/// <remarks>
/// A wildcard matches one or more characters other than <c>/</c>; <c>{name,regex}</c> gives it its own regular
/// expression, for example <c>{chunk,\d+}</c>. A name used twice must match the same text both times.
/// <c>{{</c> and <c>}}</c> are literal braces.
/// </remarks>
public sealed partial class WildcardPattern
{
    private readonly List<(string? Literal, string? Name)> _parts = [];
    private readonly Regex _regex;

    /// <summary>Parses a pattern.</summary>
    /// <exception cref="FormatException">An unbalanced brace or an invalid wildcard name.</exception>
    public WildcardPattern(string pattern)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        Text = pattern;
        var names = new List<string>();
        var regex = new StringBuilder("^");
        var literal = new StringBuilder();
        for (int i = 0; i < pattern.Length; i++)
        {
            char c = pattern[i];
            if (c == '{' && i + 1 < pattern.Length && pattern[i + 1] == '{')
            {
                literal.Append('{');
                i++;
                continue;
            }

            if (c == '}' && i + 1 < pattern.Length && pattern[i + 1] == '}')
            {
                literal.Append('}');
                i++;
                continue;
            }

            if (c == '}') throw new FormatException($"Unbalanced '}}' at {i} in pattern '{pattern}' (write '}}}}' for a literal brace).");
            if (c != '{')
            {
                literal.Append(c);
                continue;
            }

            int close = FindClose(pattern, i);
            if (close < 0) throw new FormatException($"Unclosed '{{' at {i} in pattern '{pattern}'.");
            string body = pattern[(i + 1)..close];
            int comma = body.IndexOf(',');
            string name = (comma < 0 ? body : body[..comma]).Trim();
            string? constraint = comma < 0 ? null : body[(comma + 1)..];
            if (!NameRegex().IsMatch(name)) throw new FormatException($"'{name}' is not a valid wildcard name in pattern '{pattern}'.");

            Flush(literal, regex);
            if (names.Contains(name))
            {
                regex.Append(@"\k<").Append(name).Append('>');
            }
            else
            {
                names.Add(name);
                regex.Append("(?<").Append(name).Append('>').Append(constraint ?? "[^/]+").Append(')');
            }

            _parts.Add((null, name));
            i = close;
        }

        Flush(literal, regex);
        regex.Append('$');
        WildcardNames = names;
        _regex = new Regex(regex.ToString(), RegexOptions.CultureInvariant);
    }

    private void Flush(StringBuilder literal, StringBuilder regex)
    {
        if (literal.Length == 0) return;
        _parts.Add((literal.ToString(), null));
        regex.Append(Regex.Escape(literal.ToString()));
        literal.Clear();
    }

    private static int FindClose(string s, int open)
    {
        int depth = 0;
        for (int i = open; i < s.Length; i++)
        {
            if (s[i] == '{') depth++;
            else if (s[i] == '}' && --depth == 0) return i;
        }

        return -1;
    }

    /// <summary>The pattern text.</summary>
    public string Text { get; }

    /// <summary>Wildcard names in order of first appearance.</summary>
    public IReadOnlyList<string> WildcardNames { get; }

    /// <summary>Whether the pattern has no wildcards.</summary>
    public bool IsConcrete => WildcardNames.Count == 0;

    /// <summary>Matches <paramref name="text"/> and returns the wildcard values, or null.</summary>
    public Wildcards? Match(string text)
    {
        var m = _regex.Match(text);
        if (!m.Success) return null;
        return new Wildcards(WildcardNames.Select(n => new KeyValuePair<string, string>(n, m.Groups[n].Value)));
    }

    /// <summary>Substitutes wildcard values.</summary>
    /// <exception cref="KeyNotFoundException">A wildcard in the pattern has no value.</exception>
    public string Format(IReadOnlyDictionary<string, string> values)
    {
        ArgumentNullException.ThrowIfNull(values);
        var sb = new StringBuilder();
        foreach (var (lit, name) in _parts)
        {
            if (lit is not null)
            {
                sb.Append(lit);
                continue;
            }

            if (!values.TryGetValue(name!, out var v))
                throw new KeyNotFoundException($"Pattern '{Text}' needs wildcard '{name}', which is not defined here ({(values.Count == 0 ? "no wildcards" : string.Join(", ", values.Keys))}).");
            sb.Append(v);
        }

        return sb.ToString();
    }

    /// <inheritdoc />
    public override string ToString() => Text;

    [GeneratedRegex(@"^[A-Za-z_][A-Za-z0-9_]*$")]
    private static partial Regex NameRegex();
}

/// <summary>
/// A table of samples (one row per sample, one column per attribute), read from CSV/TSV or built in code.
/// Rows double as wildcard values for <see cref="Expand.Rows(string, IEnumerable{IReadOnlyDictionary{string, string}})"/>.
/// </summary>
public sealed class SampleSheet : IReadOnlyList<IReadOnlyDictionary<string, string>>
{
    private readonly List<Dictionary<string, string>> _rows;

    /// <summary>Creates a sheet from rows; every row should have the same columns.</summary>
    public SampleSheet(IEnumerable<IReadOnlyDictionary<string, string>> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        _rows = [.. rows.Select(r => new Dictionary<string, string>(r, StringComparer.Ordinal))];
        Columns = _rows.Count == 0 ? [] : [.. _rows[0].Keys];
    }

    /// <summary>A one-column sheet, for example <c>FromColumn("sample", "A", "B")</c>.</summary>
    public static SampleSheet FromColumn(string column, params string[] values)
        => new(values.Select(v => (IReadOnlyDictionary<string, string>)new Dictionary<string, string> { [column] = v }));

    /// <summary>Parses CSV or TSV text with a header row (the separator is a tab when the header has one, else a comma). Lines starting with <c>#</c> are skipped.</summary>
    /// <exception cref="FormatException">A row has a different number of fields than the header.</exception>
    public static SampleSheet Parse(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = text.Split('\n').Select(l => l.TrimEnd('\r')).Where(l => l.Trim().Length > 0 && !l.StartsWith('#')).ToList();
        if (lines.Count == 0) return new([]);
        char sep = lines[0].Contains('\t') ? '\t' : ',';
        var header = lines[0].Split(sep).Select(h => h.Trim()).ToArray();
        var rows = new List<IReadOnlyDictionary<string, string>>();
        for (int i = 1; i < lines.Count; i++)
        {
            var fields = lines[i].Split(sep);
            if (fields.Length != header.Length)
                throw new FormatException($"Sample sheet row {i + 1} has {fields.Length} fields; the header has {header.Length}.");
            var row = new Dictionary<string, string>(StringComparer.Ordinal);
            for (int j = 0; j < header.Length; j++) row[header[j]] = fields[j].Trim();
            rows.Add(row);
        }

        return new(rows);
    }

    /// <summary>Reads a CSV/TSV sample sheet.</summary>
    public static SampleSheet Load(string path) => Parse(File.ReadAllText(path));

    /// <summary>Column names, from the first row.</summary>
    public IReadOnlyList<string> Columns { get; }

    /// <inheritdoc />
    public int Count => _rows.Count;

    /// <inheritdoc />
    public IReadOnlyDictionary<string, string> this[int index] => _rows[index];

    /// <summary>All values of one column, in row order.</summary>
    public IReadOnlyList<string> Column(string name) => [.. _rows.Select(r => r.TryGetValue(name, out var v)
        ? v : throw new KeyNotFoundException($"The sample sheet has no column '{name}' (columns: {string.Join(", ", Columns)})."))];

    /// <summary>The row whose <paramref name="column"/> equals <paramref name="value"/>.</summary>
    /// <exception cref="KeyNotFoundException">No such row.</exception>
    public IReadOnlyDictionary<string, string> Find(string column, string value)
        => _rows.FirstOrDefault(r => r.TryGetValue(column, out var v) && v == value)
           ?? throw new KeyNotFoundException($"No sample sheet row has {column} = '{value}'.");

    /// <inheritdoc />
    public IEnumerator<IReadOnlyDictionary<string, string>> GetEnumerator() => _rows.Cast<IReadOnlyDictionary<string, string>>().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}

/// <summary>Expands patterns into concrete paths: the gather side of scatter/gather, like Snakemake's <c>expand</c>.</summary>
public static class Expand
{
    /// <summary>The cross product of the given wildcard values: <c>Expand.Product("calls/{sample}.{chunk}.vcf", ("sample", s), ("chunk", c))</c>.</summary>
    public static IReadOnlyList<string> Product(string pattern, params (string Name, IEnumerable<string> Values)[] wildcards)
    {
        var p = new WildcardPattern(pattern);
        IEnumerable<Dictionary<string, string>> combos = [new Dictionary<string, string>()];
        foreach (var (name, values) in wildcards)
        {
            var list = values.ToList();
            combos = combos.SelectMany(c => list.Select(v => new Dictionary<string, string>(c) { [name] = v }));
        }

        return [.. combos.Select(c => p.Format(c))];
    }

    /// <summary>One path per value of a single wildcard: <c>Expand.Each("calls/{sample}.vcf", "sample", samples)</c>.</summary>
    public static IReadOnlyList<string> Each(string pattern, string name, IEnumerable<string> values) => Product(pattern, (name, values));

    /// <summary>One path per row (sample-sheet rows or any wildcard dictionaries).</summary>
    public static IReadOnlyList<string> Rows(string pattern, IEnumerable<IReadOnlyDictionary<string, string>> rows)
    {
        var p = new WildcardPattern(pattern);
        return [.. rows.Select(p.Format)];
    }

}
