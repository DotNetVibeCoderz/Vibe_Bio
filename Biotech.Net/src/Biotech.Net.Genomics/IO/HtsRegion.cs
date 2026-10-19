using System.Globalization;
using Biotech.Net.Core;

namespace Biotech.Net.Genomics.IO;

/// <summary>A query region: reference name and 0-based half-open <c>[Start, End)</c>.</summary>
/// <param name="Name">Reference name.</param>
/// <param name="Start">0-based start.</param>
/// <param name="End">Exclusive end (<see cref="long.MaxValue"/> for "to the end of the reference").</param>
public readonly record struct HtsRegion(string Name, long Start, long End)
{
    /// <summary>
    /// Parses a samtools/tabix region: <c>chr1</c> (whole reference), <c>chr1:1000</c> (from base 1000 to the end) or
    /// <c>chr1:1,000-2,000</c> (1-based, inclusive). When <paramref name="isReference"/> recognises the whole text
    /// as a name (names may contain <c>:</c>), it is taken as a whole reference.
    /// </summary>
    public static HtsRegion Parse(string region, Func<string, bool>? isReference = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(region);
        if (isReference?.Invoke(region) == true) return new HtsRegion(region, 0, long.MaxValue);
        int colon = region.LastIndexOf(':');
        if (colon <= 0) return new HtsRegion(region, 0, long.MaxValue);
        string name = region[..colon];
        var range = region.AsSpan(colon + 1).ToString().Replace(",", "");
        if (range.Length == 0) return new HtsRegion(name, 0, long.MaxValue);
        int dash = range.IndexOf('-');
        var inv = CultureInfo.InvariantCulture;
        if (dash < 0)
        {
            if (!long.TryParse(range, NumberStyles.None, inv, out long b)) return new HtsRegion(region, 0, long.MaxValue);
            return new HtsRegion(name, Math.Max(0, b - 1), long.MaxValue);
        }
        long beg = dash == 0 ? 1 : long.Parse(range.AsSpan(0, dash), NumberStyles.None, inv);
        long end = dash == range.Length - 1 ? long.MaxValue : long.Parse(range.AsSpan(dash + 1), NumberStyles.None, inv);
        if (end < beg) throw new FormatException($"Region '{region}' ends before it starts.");
        return new HtsRegion(name, Math.Max(0, beg - 1), end);
    }

    /// <summary>The region of a <see cref="GenomicInterval"/>.</summary>
    public static HtsRegion From(GenomicInterval interval) => new(interval.Chromosome, interval.Start, interval.End);

    public override string ToString() => End == long.MaxValue ? $"{Name}:{Start + 1}-" : $"{Name}:{Start + 1}-{End}";
}
