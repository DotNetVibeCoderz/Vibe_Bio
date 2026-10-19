using System.Globalization;
using Biotech.Net.Core;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Vcf;

/// <summary>Writes VCF 4.x text (plain, or gzip when created from a <c>.gz</c> path).</summary>
public sealed class VcfWriter : IDisposable
{
    private readonly TextWriter _writer;
    private readonly bool _leaveOpen;
    private bool _headerWritten;

    /// <summary>Creates a writer; the header is written before the first record (or on <see cref="Flush"/>).</summary>
    public VcfWriter(TextWriter writer, VcfHeader header, bool leaveOpen = false)
    {
        _writer = writer ?? throw new ArgumentNullException(nameof(writer));
        Header = header ?? throw new ArgumentNullException(nameof(header));
        _leaveOpen = leaveOpen;
    }

    /// <summary>Creates a VCF file (gzip when the name ends in <c>.gz</c>).</summary>
    /// <remarks>A <c>.gz</c>/<c>.bgz</c> path is written as BGZF (gzip-compatible), so it can be indexed with <see cref="TabixIndex.Build(string, TabixConfig, BinningIndexFormat, int)"/>.</remarks>
    public static VcfWriter Create(string path, VcfHeader header)
    {
        if (!BioFile.IsGzipPath(path)) return new(BioFile.CreateText(path), header);
        var w = new StreamWriter(BgzfWriter.Create(path), new System.Text.UTF8Encoding(false), 1 << 16) { NewLine = "\n" };
        return new(w, header);
    }

    /// <summary>The header being written.</summary>
    public VcfHeader Header { get; }

    private void EnsureHeader()
    {
        if (_headerWritten) return;
        Header.Write(_writer);
        _headerWritten = true;
    }

    /// <summary>Writes one record.</summary>
    public void Write(VcfRecord record)
    {
        EnsureHeader();
        WriteRecord(_writer, record);
        _writer.Write('\n');
    }

    /// <summary>Writes a variant (genotypes only).</summary>
    public void Write(Variant variant) => Write(new VcfRecord(variant));

    /// <summary>Writes several records.</summary>
    public void Write(IEnumerable<VcfRecord> records)
    {
        foreach (var r in records) Write(r);
    }

    /// <summary>Writes a header and records to a string.</summary>
    public static string ToString(VcfHeader header, IEnumerable<VcfRecord> records)
    {
        var sw = new StringWriter { NewLine = "\n" };
        using (var w = new VcfWriter(sw, header, leaveOpen: true))
        {
            w.Write(records);
            w.Flush();
        }
        return sw.ToString();
    }

    /// <summary>Writes a header and records to a file.</summary>
    public static void Write(string path, VcfHeader header, IEnumerable<VcfRecord> records)
    {
        using var w = Create(path, header);
        w.Write(records);
    }

    /// <summary>Flushes (writing the header if no record was written yet).</summary>
    public void Flush()
    {
        EnsureHeader();
        _writer.Flush();
    }

    internal static void WriteRecord(TextWriter w, VcfRecord record)
    {
        var v = record.Variant;
        w.Write(v.Chromosome); w.Write('\t');
        w.Write(v.Position.ToString(CultureInfo.InvariantCulture)); w.Write('\t');
        w.Write(v.Id ?? "."); w.Write('\t');
        w.Write(v.Reference); w.Write('\t');
        if (v.Alternates.Count == 0) w.Write('.');
        else
            for (int i = 0; i < v.Alternates.Count; i++)
            {
                if (i > 0) w.Write(',');
                w.Write(v.Alternates[i]);
            }
        w.Write('\t');
        w.Write(double.IsNaN(v.Quality) ? "." : FormatDouble(v.Quality)); w.Write('\t');
        w.Write(string.IsNullOrEmpty(v.Filter) ? "." : v.Filter); w.Write('\t');
        if (v.Info.Count == 0) w.Write('.');
        else
        {
            bool first = true;
            foreach (var kv in v.Info)
            {
                if (!first) w.Write(';');
                first = false;
                w.Write(kv.Key);
                if (kv.Value.Length > 0) { w.Write('='); w.Write(kv.Value); }
            }
        }

        int samples = record.SampleCount;
        if (record.Format.Count == 0 || samples == 0) return;
        w.Write('\t');
        w.Write(string.Join(':', record.Format));
        if (v.Genotypes.Count == 0 && record.HasUnparsedSamples)
        {
            w.Write('\t');
            w.Write(record.RawSamples);
            return;
        }
        var values = record.Values();
        for (int s = 0; s < samples; s++)
        {
            w.Write('\t');
            // Trailing missing fields may be dropped (VCF 4.x), but GT is always written.
            int last = record.Format.Count - 1;
            while (last > 0 && Field(record, values, s, last) == ".") last--;
            for (int f = 0; f <= last; f++)
            {
                if (f > 0) w.Write(':');
                w.Write(Field(record, values, s, f));
            }
        }
    }

    private static string Field(VcfRecord record, List<string[]> values, int sample, int f)
    {
        if (record.Format[f] == "GT")
            return sample < record.Variant.Genotypes.Count ? record.Variant.Genotypes[sample].ToString() : ".";
        if (sample >= values.Count) return ".";
        var row = values[sample];
        return f < row.Length ? row[f] : ".";
    }

    internal static string FormatDouble(double d)
    {
        if (d == Math.Floor(d) && Math.Abs(d) < 1e15) return ((long)d).ToString(CultureInfo.InvariantCulture);
        return d.ToString("0.######", CultureInfo.InvariantCulture);
    }

    /// <summary>Flushes and, unless <c>leaveOpen</c>, closes the underlying writer.</summary>
    public void Dispose()
    {
        Flush();
        if (!_leaveOpen) _writer.Dispose();
    }
}

/// <summary>Streaming filters over VCF records and variants.</summary>
public static class VcfFilters
{
    /// <summary>Records with QUAL ≥ <paramref name="minQuality"/> (missing QUAL fails).</summary>
    public static IEnumerable<VcfRecord> MinQuality(this IEnumerable<VcfRecord> records, double minQuality)
        => records.Where(r => !double.IsNaN(r.Quality) && r.Quality >= minQuality);

    /// <summary>Records whose FILTER is <c>PASS</c> (or <c>.</c>, i.e. unfiltered).</summary>
    public static IEnumerable<VcfRecord> PassOnly(this IEnumerable<VcfRecord> records) => records.Where(r => r.Variant.IsPass);

    /// <summary>Records whose reference allele overlaps <paramref name="region"/> (0-based half-open).</summary>
    public static IEnumerable<VcfRecord> InRegion(this IEnumerable<VcfRecord> records, GenomicInterval region)
        => records.Where(r => r.Variant.Interval.Overlaps(region));

    /// <summary>Records overlapping any of <paramref name="regions"/> (merged and binary-searched per chromosome).</summary>
    public static IEnumerable<VcfRecord> InRegions(this IEnumerable<VcfRecord> records, IEnumerable<GenomicInterval> regions)
    {
        var byChrom = GenomicInterval.Merge(regions).GroupBy(r => r.Chromosome)
            .ToDictionary(g => g.Key, g => g.ToArray(), StringComparer.Ordinal);
        foreach (var r in records)
        {
            if (!byChrom.TryGetValue(r.Chromosome, out var list)) continue;
            var iv = r.Variant.Interval;
            // First merged interval whose End > iv.Start.
            int lo = 0, hi = list.Length;
            while (lo < hi)
            {
                int mid = (lo + hi) >>> 1;
                if (list[mid].End <= iv.Start) lo = mid + 1; else hi = mid;
            }
            if (lo < list.Length && list[lo].Start < iv.End) yield return r;
        }
    }

    /// <summary>Bi-allelic SNVs only.</summary>
    public static IEnumerable<VcfRecord> BiallelicSnvs(this IEnumerable<VcfRecord> records)
        => records.Where(r => r.Alternates.Count == 1 && r.Variant.Type == VariantType.Snv);

    /// <summary>Variants with QUAL ≥ <paramref name="minQuality"/>.</summary>
    public static IEnumerable<Variant> MinQuality(this IEnumerable<Variant> variants, double minQuality)
        => variants.Where(v => !double.IsNaN(v.Quality) && v.Quality >= minQuality);

    /// <summary>Variants that pass all filters.</summary>
    public static IEnumerable<Variant> PassOnly(this IEnumerable<Variant> variants) => variants.Where(v => v.IsPass);

    /// <summary>Variants overlapping <paramref name="region"/>.</summary>
    public static IEnumerable<Variant> InRegion(this IEnumerable<Variant> variants, GenomicInterval region)
        => variants.Where(v => v.Interval.Overlaps(region));
}
