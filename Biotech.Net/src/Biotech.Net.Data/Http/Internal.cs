using System.Globalization;
using System.Text;
using System.Text.Json;

namespace Biotech.Net.Data.Http;

/// <summary>Builds <c>name=value&amp;...</c> strings with RFC 3986 escaping, skipping null values.</summary>
internal sealed class QueryString
{
    private readonly List<KeyValuePair<string, string>> _pairs = [];

    public QueryString Add(string name, string? value)
    {
        if (value is not null) _pairs.Add(new(name, value));
        return this;
    }

    public QueryString Add(string name, int? value) => Add(name, value?.ToString(CultureInfo.InvariantCulture));

    public QueryString Add(string name, double? value) => Add(name, value?.ToString("R", CultureInfo.InvariantCulture));

    public int Count => _pairs.Count;

    public override string ToString()
    {
        var sb = new StringBuilder();
        foreach (var (k, v) in _pairs)
        {
            if (sb.Length > 0) sb.Append('&');
            sb.Append(Uri.EscapeDataString(k)).Append('=').Append(Uri.EscapeDataString(v));
        }

        return sb.ToString();
    }

    public Uri Apply(string baseUrl) => new(_pairs.Count == 0 ? baseUrl : baseUrl + (baseUrl.Contains('?') ? "&" : "?") + this);
}

/// <summary>Tolerant accessors over <see cref="JsonElement"/>: missing, null or wrongly typed members give null/empty.</summary>
internal static class Json
{
    public static JsonElement? Prop(this JsonElement e, string name) =>
        e.ValueKind == JsonValueKind.Object && e.TryGetProperty(name, out var v) && v.ValueKind != JsonValueKind.Null ? v : null;

    public static JsonElement? Path(this JsonElement e, params string[] names)
    {
        JsonElement? cur = e;
        foreach (string n in names)
        {
            cur = cur?.Prop(n);
            if (cur is null) return null;
        }

        return cur;
    }

    public static string? Str(this JsonElement e, string name) => e.Prop(name) is { } v ? AsString(v) : null;

    public static string? AsString(this JsonElement v) => v.ValueKind switch
    {
        JsonValueKind.String => v.GetString(),
        JsonValueKind.Number or JsonValueKind.True or JsonValueKind.False => v.GetRawText(),
        _ => null,
    };

    public static long? Long(this JsonElement e, string name) => e.Prop(name) is { } v ? AsLong(v) : null;

    public static long? AsLong(this JsonElement v) => v.ValueKind switch
    {
        JsonValueKind.Number when v.TryGetInt64(out long l) => l,
        JsonValueKind.Number => (long)v.GetDouble(),
        JsonValueKind.String when long.TryParse(v.GetString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out long l) => l,
        _ => null,
    };

    public static int? Int(this JsonElement e, string name) => e.Long(name) is { } l ? (int)l : null;

    public static double? Double(this JsonElement e, string name) => e.Prop(name) is { } v ? AsDouble(v) : null;

    public static double? AsDouble(this JsonElement v) => v.ValueKind switch
    {
        JsonValueKind.Number => v.GetDouble(),
        JsonValueKind.String when double.TryParse(v.GetString(), NumberStyles.Float, CultureInfo.InvariantCulture, out double d) => d,
        _ => null,
    };

    public static bool Bool(this JsonElement e, string name) => e.Prop(name) is { } v && v.ValueKind switch
    {
        JsonValueKind.True => true,
        JsonValueKind.Number => v.GetDouble() != 0,
        JsonValueKind.String => v.GetString() is "1" or "true" or "Y" or "y",
        _ => false,
    };

    public static IEnumerable<JsonElement> Arr(this JsonElement e, string name) =>
        e.Prop(name) is { ValueKind: JsonValueKind.Array } a ? a.EnumerateArray() : [];

    public static IEnumerable<JsonElement> Items(this JsonElement? e) =>
        e is { ValueKind: JsonValueKind.Array } a ? a.EnumerateArray() : [];

    public static IReadOnlyList<string> Strings(this JsonElement e, string name) =>
        e.Arr(name).Select(AsString).OfType<string>().ToList();
}
