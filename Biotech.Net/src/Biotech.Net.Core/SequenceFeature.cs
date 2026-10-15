namespace Biotech.Net.Core;

/// <summary>
/// The location of a feature on its parent sequence: one or more 0-based half-open segments
/// (several for <c>join(...)</c>), on one strand. Segments are kept in biological order.
/// </summary>
public sealed record FeatureLocation(IReadOnlyList<(int Start, int End)> Segments, Strand Strand = Strand.Forward)
{
    /// <summary>A single contiguous segment.</summary>
    public static FeatureLocation Span(int start, int end, Strand strand = Strand.Forward) => new([(start, end)], strand);

    /// <summary>Leftmost 0-based position.</summary>
    public int Start => Segments.Min(s => s.Start);

    /// <summary>Rightmost position (exclusive).</summary>
    public int End => Segments.Max(s => s.End);

    /// <summary>Total bases covered by all segments.</summary>
    public int Length => Segments.Sum(s => s.End - s.Start);

    /// <summary>Whether the location is split over several segments.</summary>
    public bool IsJoin => Segments.Count > 1;

    /// <summary>
    /// The residues of this location extracted from <paramref name="parent"/>: segments concatenated, and the
    /// result reverse-complemented for features on the reverse strand.
    /// </summary>
    public Sequence Extract(Sequence parent, string? id = null)
    {
        var buffer = new byte[Length];
        int o = 0;
        foreach (var (s, e) in Segments)
        {
            parent.Span[s..e].CopyTo(buffer.AsSpan(o));
            o += e - s;
        }
        var seq = Sequence.FromBuffer(id ?? parent.Id, buffer, parent.Alphabet, parent.Description, validate: false);
        return Strand == Strand.Reverse ? seq.ReverseComplement() : seq;
    }

    /// <summary>GenBank-style location string, 1-based (e.g. <c>complement(join(1..10,20..30))</c>).</summary>
    public override string ToString()
    {
        string body = string.Join(",", Segments.Select(s => $"{s.Start + 1}..{s.End}"));
        if (IsJoin) body = $"join({body})";
        return Strand == Strand.Reverse ? $"complement({body})" : body;
    }
}

/// <summary>An annotated feature of a sequence (gene, CDS, mRNA, promoter, ...).</summary>
public sealed class SequenceFeature
{
    public SequenceFeature(string type, FeatureLocation location)
    {
        ArgumentException.ThrowIfNullOrEmpty(type);
        Type = type;
        Location = location ?? throw new ArgumentNullException(nameof(location));
    }

    /// <summary>Feature key, e.g. <c>gene</c>, <c>CDS</c>, <c>mRNA</c>.</summary>
    public string Type { get; }

    /// <summary>Where the feature lies on its parent sequence.</summary>
    public FeatureLocation Location { get; }

    /// <summary>Qualifiers (<c>/gene="lacZ"</c>, <c>/product=...</c>); a key may repeat, so values are lists.</summary>
    public Dictionary<string, List<string>> Qualifiers { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>First value of a qualifier, or null.</summary>
    public string? this[string qualifier] => Qualifiers.TryGetValue(qualifier, out var v) && v.Count > 0 ? v[0] : null;

    /// <summary>Adds a qualifier value.</summary>
    public SequenceFeature With(string qualifier, string value)
    {
        if (!Qualifiers.TryGetValue(qualifier, out var list)) Qualifiers[qualifier] = list = [];
        list.Add(value);
        return this;
    }

    /// <summary>The best available name: gene, locus_tag, product or label.</summary>
    public string Name => this["gene"] ?? this["locus_tag"] ?? this["product"] ?? this["label"] ?? Type;

    public override string ToString() => $"{Type} {Location} {Name}";
}
