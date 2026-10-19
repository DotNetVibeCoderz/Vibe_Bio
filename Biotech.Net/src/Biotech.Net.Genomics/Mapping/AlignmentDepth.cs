namespace Biotech.Net.Genomics.Mapping;

/// <summary>Filters for <see cref="AlignmentDepth"/>; the defaults are those of <c>samtools depth</c>.</summary>
public sealed record DepthOptions
{
    /// <summary>Reads with any of these flags are skipped (default UNMAP, SECONDARY, QCFAIL, DUP).</summary>
    public SamFlags ExcludeFlags { get; init; } = SamFlags.Unmapped | SamFlags.Secondary | SamFlags.QcFail | SamFlags.Duplicate;

    /// <summary>Minimum mapping quality (default 0).</summary>
    public int MinMapQ { get; init; }

    /// <summary>Minimum base quality (default 0; reads without qualities count as 255, as in BAM).</summary>
    public int MinBaseQuality { get; init; }

    /// <summary>Count deleted bases (D) as covered, like <c>samtools depth -J</c> (default false).</summary>
    public bool CountDeletions { get; init; }
}

/// <summary>Per-base read depth over a region (like <c>samtools depth -a</c>).</summary>
public static class AlignmentDepth
{
    /// <summary>Depth at every 0-based position of <c>[start, end)</c> on <paramref name="referenceName"/>, using the reader's index when available.</summary>
    public static int[] Compute(IAlignmentReader reader, string referenceName, long start, long end, DepthOptions? options = null)
        => Compute(reader.Query(referenceName, start, end), referenceName, start, end, options);

    /// <summary>Depth at every 0-based position of <c>[start, end)</c> from <paramref name="records"/>.</summary>
    public static int[] Compute(IEnumerable<AlignmentRecord> records, string referenceName, long start, long end, DepthOptions? options = null)
    {
        options ??= new DepthOptions();
        if (end <= start) return [];
        var depth = new int[checked((int)(end - start))];
        foreach (var r in records)
        {
            if ((r.Flags & options.ExcludeFlags) != 0 || r.MapQ < options.MinMapQ || r.ReferenceName != referenceName || r.Position < 0) continue;
            long rp = r.Position;
            int q = 0;
            foreach (var e in r.Cigar)
            {
                switch (e.Op)
                {
                    case 'M' or '=' or 'X':
                        for (int i = 0; i < e.Length; i++, rp++, q++)
                        {
                            if (rp < start || rp >= end) continue;
                            int bq = r.Qualities is null ? 255 : q < r.Qualities.Length ? r.Qualities[q] : 0;
                            if (bq >= options.MinBaseQuality) depth[rp - start]++;
                        }
                        break;
                    case 'D':
                        if (options.CountDeletions)
                        {
                            // samtools depth -J judges a deletion by the quality of the next read base
                            int dq = r.Qualities is null ? 255 : q < r.Qualities.Length ? r.Qualities[q] : 0;
                            if (dq >= options.MinBaseQuality)
                                for (int i = 0; i < e.Length; i++) if (rp + i >= start && rp + i < end) depth[rp + i - start]++;
                        }
                        rp += e.Length;
                        break;
                    case 'N': rp += e.Length; break;
                    case 'I' or 'S': q += e.Length; break;
                }
            }
        }
        return depth;
    }
}

/// <summary>Orders records by reference id then position (unplaced last), as <c>samtools sort</c> does.</summary>
public sealed class CoordinateComparer : IComparer<AlignmentRecord>
{
    /// <summary>The shared instance.</summary>
    public static CoordinateComparer Instance { get; } = new();

    /// <summary>A stable coordinate sort (ties keep their input order, as <c>samtools sort</c>).</summary>
    public static List<AlignmentRecord> Sort(IEnumerable<AlignmentRecord> records) => records.OrderBy(r => r, Instance).ToList();

    /// <inheritdoc/>
    public int Compare(AlignmentRecord? x, AlignmentRecord? y)
    {
        if (ReferenceEquals(x, y)) return 0;
        if (x is null) return -1;
        if (y is null) return 1;
        uint a = (uint)x.ReferenceId, b = (uint)y.ReferenceId; // -1 sorts last
        int c = a.CompareTo(b);
        if (c != 0) return c;
        c = x.Position.CompareTo(y.Position);
        return c != 0 ? c : ((int)x.Flags & 0x10).CompareTo((int)y.Flags & 0x10);
    }
}
