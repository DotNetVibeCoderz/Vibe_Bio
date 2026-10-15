namespace Biotech.Net.Core;

/// <summary>Strand of a genomic feature.</summary>
public enum Strand : sbyte
{
    Reverse = -1,
    Unknown = 0,
    Forward = 1,
}

/// <summary>
/// A region on a named sequence, 0-based and half-open (<c>[Start, End)</c>) like BED and most
/// libraries; <see cref="ToString"/> prints the 1-based, closed UCSC form (<c>chr1:101-200</c>).
/// </summary>
public readonly record struct GenomicInterval(string Chromosome, long Start, long End, Strand Strand = Strand.Unknown)
    : IComparable<GenomicInterval>
{
    /// <summary>Number of bases covered.</summary>
    public long Length => End - Start;

    /// <summary>Validates and creates an interval.</summary>
    public static GenomicInterval Create(string chromosome, long start, long end, Strand strand = Strand.Unknown)
    {
        ArgumentException.ThrowIfNullOrEmpty(chromosome);
        if (start < 0 || end < start) throw new ArgumentOutOfRangeException(nameof(end), $"Invalid interval [{start}, {end}).");
        return new(chromosome, start, end, strand);
    }

    /// <summary>
    /// Parses <c>chr1:101-200</c> (1-based closed, commas allowed), <c>chr1:101-200:-</c> or a bare <c>chr1</c>
    /// (which needs <paramref name="chromosomeLength"/>).
    /// </summary>
    public static GenomicInterval Parse(string region, long chromosomeLength = long.MaxValue)
    {
        var parts = region.Replace(",", "").Split(':');
        string chrom = parts[0];
        if (parts.Length == 1) return Create(chrom, 0, chromosomeLength);
        var range = parts[1].Split('-');
        long start = long.Parse(range[0]) - 1;
        long end = range.Length > 1 ? long.Parse(range[1]) : start + 1;
        var strand = parts.Length > 2 ? (parts[2] == "-" ? Strand.Reverse : parts[2] == "+" ? Strand.Forward : Strand.Unknown) : Strand.Unknown;
        return Create(chrom, start, end, strand);
    }

    /// <summary>Whether the two intervals share at least one base (strand ignored).</summary>
    public bool Overlaps(GenomicInterval other)
        => Chromosome == other.Chromosome && Start < other.End && other.Start < End;

    /// <summary>Number of shared bases (0 when disjoint).</summary>
    public long OverlapLength(GenomicInterval other)
        => Chromosome != other.Chromosome ? 0 : Math.Max(0, Math.Min(End, other.End) - Math.Max(Start, other.Start));

    /// <summary>Whether a 0-based position lies inside the interval.</summary>
    public bool Contains(long position) => position >= Start && position < End;

    /// <summary>Whether <paramref name="other"/> lies entirely inside this interval.</summary>
    public bool Contains(GenomicInterval other) => Chromosome == other.Chromosome && other.Start >= Start && other.End <= End;

    /// <summary>The smallest interval covering both.</summary>
    public GenomicInterval Union(GenomicInterval other)
    {
        if (Chromosome != other.Chromosome) throw new ArgumentException("Intervals are on different sequences.");
        return this with { Start = Math.Min(Start, other.Start), End = Math.Max(End, other.End) };
    }

    /// <summary>Distance between two non-overlapping intervals on the same sequence (0 when they overlap, −1 across sequences).</summary>
    public long DistanceTo(GenomicInterval other)
    {
        if (Chromosome != other.Chromosome) return -1;
        if (Overlaps(other)) return 0;
        return other.Start >= End ? other.Start - End : Start - other.End;
    }

    /// <summary>Extends both sides by <paramref name="upstream"/>/<paramref name="downstream"/> bases, clamped at 0.</summary>
    public GenomicInterval Expand(long upstream, long downstream)
        => this with { Start = Math.Max(0, Start - upstream), End = End + downstream };

    public int CompareTo(GenomicInterval other)
    {
        int c = string.CompareOrdinal(Chromosome, other.Chromosome);
        if (c != 0) return c;
        c = Start.CompareTo(other.Start);
        return c != 0 ? c : End.CompareTo(other.End);
    }

    /// <summary>Merges overlapping or touching intervals (per chromosome, strand ignored).</summary>
    public static IReadOnlyList<GenomicInterval> Merge(IEnumerable<GenomicInterval> intervals)
    {
        var sorted = intervals.OrderBy(i => i).ToList();
        var result = new List<GenomicInterval>(sorted.Count);
        foreach (var iv in sorted)
        {
            if (result.Count > 0 && result[^1].Chromosome == iv.Chromosome && iv.Start <= result[^1].End)
                result[^1] = result[^1] with { End = Math.Max(result[^1].End, iv.End), Strand = Strand.Unknown };
            else
                result.Add(iv);
        }
        return result;
    }

    public override string ToString()
    {
        string s = $"{Chromosome}:{Start + 1:N0}-{End:N0}";
        return Strand switch { Strand.Forward => s + "(+)", Strand.Reverse => s + "(-)", _ => s };
    }
}
