using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>A parsed INSDC (GenBank/EMBL/DDBJ) location with the details <see cref="FeatureLocation"/> cannot hold.</summary>
/// <param name="Location">Segments and strand (0-based half-open).</param>
/// <param name="PartialStart">The 5' end is incomplete (<c>&lt;</c> on the first base).</param>
/// <param name="PartialEnd">The 3' end is incomplete (<c>&gt;</c> on the last base).</param>
/// <param name="IsOrder">The segments were combined with <c>order(...)</c> rather than <c>join(...)</c>.</param>
/// <param name="HasRemoteParts">Some segments referenced another entry (<c>J00194.1:100..202</c>) and were dropped.</param>
public sealed record InsdcLocation(FeatureLocation Location, bool PartialStart = false, bool PartialEnd = false,
    bool IsOrder = false, bool HasRemoteParts = false);

/// <summary>
/// Parses and formats INSDC feature locations: <c>467</c>, <c>340..565</c>, <c>&lt;345..500</c>, <c>102.110</c>,
/// <c>123^124</c>, <c>join(...)</c>, <c>order(...)</c>, <c>complement(...)</c> and nestings of them.
/// </summary>
/// <remarks>
/// Core's <see cref="FeatureLocation"/> concatenates segments and then reverse-complements for the reverse strand,
/// so reverse-strand segments are stored in ascending genomic order whichever way the file spelled them
/// (<c>complement(join(a,b))</c> and <c>join(complement(b),complement(a))</c> both become <c>[a, b]</c>, reverse).
/// Mixed-strand (trans-spliced) locations keep the strand of the majority of bases.
/// </remarks>
public static class InsdcLocationParser
{
    private static readonly ConditionalWeakTable<SequenceFeature, InsdcLocation> Details = new();

    /// <summary>Parses a location string.</summary>
    /// <exception cref="FormatException">The string is not a valid location.</exception>
    public static InsdcLocation Parse(string text)
    {
        var p = new Parser(text);
        var parts = p.ParseLocation(false);
        p.ExpectEnd();
        bool remote = p.Remote;
        if (parts.Count == 0) throw new FormatException($"Location '{text}' has no local segments.");
        long rev = 0, fwd = 0;
        foreach (var s in parts) { if (s.Reverse) rev += s.End - s.Start; else fwd += s.End - s.Start; }
        var strand = rev > fwd || (rev == fwd && parts.All(s => s.Reverse)) ? Strand.Reverse : Strand.Forward;
        // parts are in biological order; reverse-strand storage is ascending (see remarks).
        var ordered = strand == Strand.Reverse ? Enumerable.Reverse(parts).ToList() : parts;
        var segments = ordered.Select(s => (s.Start, s.End)).ToList();
        // Partial flags refer to the biological 5'/3' ends.
        var first = parts[0];
        var last = parts[^1];
        bool partialStart = first.Reverse ? first.FuzzyEnd : first.FuzzyStart;
        bool partialEnd = last.Reverse ? last.FuzzyStart : last.FuzzyEnd;
        return new InsdcLocation(new FeatureLocation(segments, strand), partialStart, partialEnd, p.SawOrder, remote);
    }

    /// <summary>Formats a location in INSDC syntax, 1-based.</summary>
    public static string Format(InsdcLocation location)
    {
        var loc = location.Location;
        var segs = loc.Segments;
        var sb = new StringBuilder();
        bool rev = loc.Strand == Strand.Reverse;
        for (int i = 0; i < segs.Count; i++)
        {
            if (i > 0) sb.Append(',');
            // Fuzzy markers sit on the genomic left/right ends; for the reverse strand the 5' end is on the right.
            bool fuzzyLeft = rev ? location.PartialEnd && i == 0 : location.PartialStart && i == 0;
            bool fuzzyRight = rev ? location.PartialStart && i == segs.Count - 1 : location.PartialEnd && i == segs.Count - 1;
            var (s, e) = segs[i];
            if (e == s) sb.Append(CultureInfo.InvariantCulture, $"{s}^{s + 1}");
            else if (e - s == 1 && !fuzzyLeft && !fuzzyRight) sb.Append(CultureInfo.InvariantCulture, $"{s + 1}");
            else sb.Append(CultureInfo.InvariantCulture, $"{(fuzzyLeft ? "<" : "")}{s + 1}..{(fuzzyRight ? ">" : "")}{e}");
        }
        string body = sb.ToString();
        if (segs.Count > 1) body = $"{(location.IsOrder ? "order" : "join")}({body})";
        return rev ? $"complement({body})" : body;
    }

    /// <summary>Formats a plain <see cref="FeatureLocation"/>.</summary>
    public static string Format(FeatureLocation location) => Format(new InsdcLocation(location));

    /// <summary>The parsed location details recorded for <paramref name="feature"/> by a reader, if any.</summary>
    public static InsdcLocation? GetDetails(SequenceFeature feature) => Details.TryGetValue(feature, out var d) ? d : null;

    /// <summary>Records partial/order details for <paramref name="feature"/> so writers can reproduce them.</summary>
    public static void SetDetails(SequenceFeature feature, InsdcLocation details) => Details.AddOrUpdate(feature, details);

    /// <summary>Whether the feature's 5' or 3' end is marked incomplete.</summary>
    public static bool IsPartial(SequenceFeature feature) => GetDetails(feature) is { } d && (d.PartialStart || d.PartialEnd);

    private readonly record struct Part(int Start, int End, bool Reverse, bool FuzzyStart, bool FuzzyEnd);

    private ref struct Parser(string text)
    {
        private readonly ReadOnlySpan<char> _s = text.AsSpan();
        private int _i;
        public bool SawOrder;
        public bool Remote;

        public List<Part> ParseLocation(bool reverse)
        {
            SkipWs();
            if (TryKeyword("complement"))
            {
                var inner = ParseLocation(false);
                Expect(')');
                var flipped = new List<Part>(inner.Count);
                for (int k = inner.Count - 1; k >= 0; k--)
                    flipped.Add(inner[k] with { Reverse = !inner[k].Reverse });
                return reverse ? Flip(flipped) : flipped;
            }
            bool isJoin = TryKeyword("join");
            bool isOrder = !isJoin && TryKeyword("order");
            if (isJoin || isOrder)
            {
                if (isOrder) SawOrder = true;
                var list = new List<Part>();
                while (true)
                {
                    list.AddRange(ParseLocation(false));
                    SkipWs();
                    if (Peek() == ',') { _i++; continue; }
                    Expect(')');
                    break;
                }
                return reverse ? Flip(list) : list;
            }
            if (TryKeyword("bond") || TryKeyword("gap"))
            {
                // gap(...) / bond(...) carry no residues we can represent; skip their content.
                int depth = 1;
                while (_i < _s.Length && depth > 0) { if (_s[_i] == '(') depth++; else if (_s[_i] == ')') depth--; _i++; }
                return [];
            }
            return ParseRange();
        }

        private static List<Part> Flip(List<Part> l)
        {
            l.Reverse();
            for (int k = 0; k < l.Count; k++) l[k] = l[k] with { Reverse = !l[k].Reverse };
            return l;
        }

        private List<Part> ParseRange()
        {
            SkipWs();
            // Remote reference: ACCESSION.V:range
            int save = _i;
            while (_i < _s.Length && _s[_i] != ':' && _s[_i] != ',' && _s[_i] != ')') _i++;
            if (_i < _s.Length && _s[_i] == ':')
            {
                _i++;
                ParseSimple();
                Remote = true;
                return [];
            }
            _i = save;
            return [ParseSimple()];
        }

        private Part ParseSimple()
        {
            SkipWs();
            bool fuzzyStart = false, fuzzyEnd = false;
            if (Peek() == '<') { fuzzyStart = true; _i++; }
            else if (Peek() == '>') { fuzzyStart = true; _i++; }
            int a = ParseInt();
            if (Peek() == '.' && PeekAt(1) == '.')
            {
                _i += 2;
                if (Peek() == '>') { fuzzyEnd = true; _i++; }
                else if (Peek() == '<') { fuzzyEnd = true; _i++; }
                int b = ParseInt();
                if (b < a) throw new FormatException($"Range {a}..{b} is reversed.");
                return new Part(a - 1, b, false, fuzzyStart, fuzzyEnd);
            }
            if (Peek() == '^')
            {
                _i++;
                ParseInt();
                return new Part(a, a, false, false, false); // site between two bases: zero-length
            }
            if (Peek() == '.')
            {
                // Old "102.110": a single base somewhere in the range; keep the whole range.
                _i++;
                int b = ParseInt();
                return new Part(a - 1, b, false, true, true);
            }
            return new Part(a - 1, a, false, fuzzyStart, fuzzyStart);
        }

        private int ParseInt()
        {
            SkipWs();
            int start = _i;
            while (_i < _s.Length && char.IsAsciiDigit(_s[_i])) _i++;
            if (_i == start) throw new FormatException($"Expected a number at position {start} of location '{_s.ToString()}'.");
            return int.Parse(_s[start.._i], CultureInfo.InvariantCulture);
        }

        private bool TryKeyword(string kw)
        {
            if (_s[_i..].StartsWith(kw, StringComparison.OrdinalIgnoreCase))
            {
                int j = _i + kw.Length;
                while (j < _s.Length && _s[j] == ' ') j++;
                if (j < _s.Length && _s[j] == '(') { _i = j + 1; return true; }
            }
            return false;
        }

        private readonly char Peek() => _i < _s.Length ? _s[_i] : '\0';
        private readonly char PeekAt(int o) => _i + o < _s.Length ? _s[_i + o] : '\0';
        private void SkipWs() { while (_i < _s.Length && char.IsWhiteSpace(_s[_i])) _i++; }

        private void Expect(char c)
        {
            SkipWs();
            if (Peek() != c) throw new FormatException($"Expected '{c}' at position {_i} of location '{_s.ToString()}'.");
            _i++;
        }

        public void ExpectEnd()
        {
            SkipWs();
            if (_i != _s.Length) throw new FormatException($"Unexpected text at position {_i} of location '{_s.ToString()}'.");
        }
    }
}

/// <summary>
/// The feature-table block shared by GenBank and EMBL. After removing the five-character line prefix
/// (five spaces in GenBank, <c>FT   </c> in EMBL) both use the same layout: key in columns 1–16, location
/// and qualifiers from column 17.
/// </summary>
internal sealed class FeatureTableBuilder
{
    private readonly List<SequenceFeature> _features = [];
    private string? _key;
    private readonly StringBuilder _location = new();
    private readonly List<(string Key, StringBuilder Value, bool Quoted)> _qualifiers = [];
    private bool _openQuote;

    // Values of these qualifiers are wrapped without spaces (sequence-like).
    private static readonly HashSet<string> NoSpaceJoin = new(StringComparer.OrdinalIgnoreCase) { "translation", "peptide", "rpt_unit_seq" };

    public void AddLine(ReadOnlySpan<char> content)
    {
        if (content.IsWhiteSpace()) return;
        if (!_openQuote && content.Length > 0 && content[0] != ' ')
        {
            Flush();
            content = content.TrimEnd();
            int ws = content.IndexOf(' ');
            _key = (ws < 0 ? content : content[..ws]).ToString();
            if (ws >= 0) _location.Append(content[ws..].Trim());
            return;
        }
        var t = content.Trim();
        if (_key is null) return;
        if (_openQuote)
        {
            var (_, v, _) = _qualifiers[^1];
            if (!NoSpaceJoin.Contains(_qualifiers[^1].Key) && v.Length > 0) v.Append(' ');
            v.Append(t);
            _openQuote = !QuoteClosed(v);
            return;
        }
        if (t.Length > 0 && t[0] == '/')
        {
            t = t[1..];
            int eq = t.IndexOf('=');
            if (eq < 0) { _qualifiers.Add((t.ToString(), new StringBuilder(), false)); return; }
            string key = t[..eq].ToString();
            var value = new StringBuilder().Append(t[(eq + 1)..]);
            bool quoted = value.Length > 0 && value[0] == '"';
            _qualifiers.Add((key, value, quoted));
            _openQuote = quoted && !QuoteClosed(value);
            return;
        }
        if (_qualifiers.Count == 0) { _location.Append(t); return; }
        // Continuation of an unquoted value.
        _qualifiers[^1].Value.Append(t);
    }

    // A quoted value is closed when it has an even number of quote characters (doubled quotes are escapes).
    private static bool QuoteClosed(StringBuilder v)
    {
        int q = 0;
        for (int i = 0; i < v.Length; i++) if (v[i] == '"') q++;
        return q >= 2 && q % 2 == 0;
    }

    private void Flush()
    {
        if (_key is null) return;
        string locText = _location.ToString();
        InsdcLocation loc;
        try { loc = InsdcLocationParser.Parse(locText); }
        catch (FormatException ex) { throw new InvalidDataException($"Feature '{_key}': {ex.Message}", ex); }
        var f = new SequenceFeature(_key, loc.Location);
        if (loc.PartialStart || loc.PartialEnd || loc.IsOrder || loc.HasRemoteParts) InsdcLocationParser.SetDetails(f, loc);
        foreach (var (k, v, quoted) in _qualifiers)
        {
            string value = v.ToString();
            if (quoted && value.Length >= 2 && value[0] == '"' && value[^1] == '"') value = value[1..^1];
            if (quoted) value = value.Replace("\"\"", "\"");
            f.With(k, value);
        }
        _features.Add(f);
        _key = null;
        _location.Clear();
        _qualifiers.Clear();
        _openQuote = false;
    }

    public List<SequenceFeature> Finish()
    {
        Flush();
        return _features;
    }
}

/// <summary>Writes the feature table shared by GenBank and EMBL.</summary>
internal static class FeatureTableWriter
{
    private const int Width = 79;
    private const int QualifierIndent = 21;

    // Qualifiers whose values the INSDC specification writes without quotes.
    private static readonly HashSet<string> Unquoted = new(StringComparer.OrdinalIgnoreCase)
    {
        "codon_start", "transl_table", "citation", "number", "estimated_length", "direction", "rpt_type",
        "rpt_unit_range", "tag_peptide", "transl_except", "anticodon", "mod_base", "compare", "codon",
    };

    public static void Write(TextWriter w, string prefix, IEnumerable<SequenceFeature> features)
    {
        string indent = prefix + new string(' ', QualifierIndent - prefix.Length);
        foreach (var f in features)
        {
            var details = InsdcLocationParser.GetDetails(f);
            string loc = details is not null && details.Location == f.Location
                ? InsdcLocationParser.Format(details)
                : InsdcLocationParser.Format(f.Location);
            string head = prefix + f.Type.PadRight(QualifierIndent - prefix.Length);
            WrapAt(w, head, indent, loc, ',');
            foreach (var (key, values) in f.Qualifiers)
                foreach (var value in values)
                {
                    string text = value.Length == 0 && !Unquoted.Contains(key) && IsFlag(key)
                        ? "/" + key
                        : Unquoted.Contains(key) ? $"/{key}={value}" : $"/{key}=\"{value.Replace("\"", "\"\"")}\"";
                    bool hard = key.Equals("translation", StringComparison.OrdinalIgnoreCase);
                    WrapText(w, indent, text, hard);
                }
        }
    }

    private static bool IsFlag(string key) => key is "pseudo" or "partial" or "environmental_sample" or "focus" or "germline"
        or "macronuclear" or "proviral" or "rearranged" or "ribosomal_slippage" or "transgenic" or "trans_splicing" or "circular_RNA";

    // Locations wrap after commas so each line stays within 79 columns.
    private static void WrapAt(TextWriter w, string head, string indent, string text, char breakChar)
    {
        int avail = Width - head.Length;
        string first = head;
        while (text.Length > avail)
        {
            int cut = text.LastIndexOf(breakChar, Math.Min(text.Length - 1, avail - 1));
            if (cut <= 0) cut = avail - 1;
            w.Write(first); w.Write(text.AsSpan(0, cut + 1)); w.Write('\n');
            text = text[(cut + 1)..];
            first = indent;
            avail = Width - indent.Length;
        }
        w.Write(first); w.Write(text); w.Write('\n');
    }

    // Qualifier text wraps at spaces where possible (always hard-wrapped for translations).
    private static void WrapText(TextWriter w, string indent, string text, bool hard)
    {
        int avail = Width - indent.Length;
        while (text.Length > avail)
        {
            int cut = hard ? -1 : text.LastIndexOf(' ', avail);
            if (cut <= 0)
            {
                w.Write(indent); w.Write(text.AsSpan(0, avail)); w.Write('\n');
                text = text[avail..];
            }
            else
            {
                w.Write(indent); w.Write(text.AsSpan(0, cut)); w.Write('\n');
                text = text[(cut + 1)..];
            }
        }
        w.Write(indent); w.Write(text); w.Write('\n');
    }

    /// <summary>Wraps free text (DEFINITION, ORGANISM lineage, ...) to the given width.</summary>
    public static void WrapKeyword(TextWriter w, string firstPrefix, string nextPrefix, string text, int width = 79)
    {
        string prefix = firstPrefix;
        if (text.Length == 0) { w.Write(prefix.TrimEnd()); w.Write('\n'); return; }
        while (text.Length > 0)
        {
            int avail = width - prefix.Length;
            if (text.Length <= avail) { w.Write(prefix); w.Write(text); w.Write('\n'); break; }
            int cut = text.LastIndexOf(' ', avail);
            if (cut <= 0) cut = avail;
            w.Write(prefix); w.Write(text.AsSpan(0, cut).TrimEnd()); w.Write('\n');
            text = text[cut..].TrimStart();
            prefix = nextPrefix;
        }
    }
}
