using System.Globalization;
using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Vcf;

/// <summary>
/// One VCF data line: a <see cref="Core.Variant"/> (site columns, INFO and per-sample genotypes) plus the
/// remaining per-sample FORMAT fields (AD, DP, GQ, PL, ...).
/// </summary>
/// <remarks>
/// Genotypes (<c>GT</c>) are authoritative in <see cref="Variant"/>.<see cref="Core.Variant.Genotypes"/>. The other
/// sample fields are parsed lazily from the source line the first time one is requested, so streaming through a
/// large multi-sample file costs little when only genotypes are used.
/// </remarks>
public sealed class VcfRecord
{
    private List<string> _format;
    private string? _rawSamples;
    private int _rawSampleCount;
    private List<string[]>? _values;

    /// <summary>Wraps a variant; <paramref name="format"/> lists the FORMAT keys (GT is added first when genotypes exist).</summary>
    public VcfRecord(Variant variant, IEnumerable<string>? format = null)
    {
        Variant = variant ?? throw new ArgumentNullException(nameof(variant));
        _format = format?.ToList() ?? [];
        if (variant.Genotypes.Count > 0 && !_format.Contains("GT")) _format.Insert(0, "GT");
    }

    internal VcfRecord(Variant variant, List<string> format, string? rawSamples, int sampleCount)
    {
        Variant = variant;
        _format = format;
        _rawSamples = rawSamples;
        _rawSampleCount = sampleCount;
    }

    /// <summary>The site, INFO and genotypes.</summary>
    public Variant Variant { get; }

    /// <summary>FORMAT keys in column order (GT first when present).</summary>
    public IReadOnlyList<string> Format => _format;

    /// <summary>Chromosome name.</summary>
    public string Chromosome => Variant.Chromosome;

    /// <summary>1-based position.</summary>
    public long Position => Variant.Position;

    /// <summary>Reference allele.</summary>
    public string Reference => Variant.Reference;

    /// <summary>Alternate alleles.</summary>
    public IReadOnlyList<string> Alternates => Variant.Alternates;

    /// <summary>Phred-scaled site quality (NaN when missing).</summary>
    public double Quality => Variant.Quality;

    /// <summary>FILTER column.</summary>
    public string Filter => Variant.Filter;

    /// <summary>INFO fields.</summary>
    public Dictionary<string, string> Info => Variant.Info;

    /// <summary>Number of samples carried by the record.</summary>
    public int SampleCount => Math.Max(Variant.Genotypes.Count, _values?.Count ?? _rawSampleCount);

    /// <summary>Genotype of a sample (missing when the record has no GT).</summary>
    public Genotype GetGenotype(int sample) => sample < Variant.Genotypes.Count ? Variant.Genotypes[sample] : Genotype.Missing;

    /// <summary>Raw value of FORMAT <paramref name="key"/> for a sample, or null when absent or <c>.</c>.</summary>
    public string? GetField(int sample, string key)
    {
        if (key == "GT") return sample < Variant.Genotypes.Count ? Variant.Genotypes[sample].ToString() : null;
        int f = _format.IndexOf(key);
        if (f < 0) return null;
        var values = Values();
        if (sample < 0 || sample >= values.Count) return null;
        var row = values[sample];
        if (f >= row.Length) return null;
        string v = row[f];
        return v is "." or "" ? null : v;
    }

    /// <summary>Sets FORMAT <paramref name="key"/> for a sample (adding the key to FORMAT when new). Use <see cref="Core.Variant.Genotypes"/> for GT.</summary>
    public VcfRecord SetField(int sample, string key, string value)
    {
        if (key == "GT") throw new ArgumentException("Set genotypes through Variant.Genotypes.", nameof(key));
        var values = Values();
        int f = _format.IndexOf(key);
        if (f < 0) { _format.Add(key); f = _format.Count - 1; }
        while (values.Count <= sample) values.Add([]);
        var row = values[sample];
        if (row.Length < _format.Count)
        {
            int old = row.Length;
            Array.Resize(ref row, _format.Count);
            for (int i = old; i < row.Length; i++) row[i] = ".";
            values[sample] = row;
        }
        row[f] = value;
        return this;
    }

    /// <summary>Sets an integer-list FORMAT field (e.g. AD, PL).</summary>
    public VcfRecord SetField(int sample, string key, ReadOnlySpan<int> values)
    {
        var parts = new string[values.Length];
        for (int i = 0; i < values.Length; i++) parts[i] = values[i].ToString(CultureInfo.InvariantCulture);
        return SetField(sample, key, string.Join(',', parts));
    }

    /// <summary>Allelic depths (<c>AD</c>), one per allele, or null.</summary>
    public int[]? AlleleDepths(int sample) => IntList(GetField(sample, "AD"));

    /// <summary>Read depth (<c>DP</c>), or null.</summary>
    public int? Depth(int sample) => Int(GetField(sample, "DP"));

    /// <summary>Genotype quality (<c>GQ</c>), or null.</summary>
    public int? GenotypeQuality(int sample) => Int(GetField(sample, "GQ"));

    /// <summary>Phred-scaled genotype likelihoods (<c>PL</c>) in VCF genotype order, or null.</summary>
    public int[]? PhredLikelihoods(int sample) => IntList(GetField(sample, "PL"));

    /// <summary>
    /// VCF ordering index of genotype <c>j/k</c> (<c>j ≤ k</c>) in a G-numbered field: <c>k(k+1)/2 + j</c>.
    /// </summary>
    public static int GenotypeIndex(int j, int k) => j > k ? GenotypeIndex(k, j) : k * (k + 1) / 2 + j;

    /// <summary>
    /// Splits a multi-allelic record into one bi-allelic record per alternate (like <c>bcftools norm -m-</c>):
    /// genotype alleles of other alternates become reference, AD keeps ref plus the chosen alternate, PL keeps the
    /// three genotypes over {ref, alt}. Other FORMAT fields are copied unchanged; INFO fields of Number=A/R type
    /// are not split (they are copied as-is).
    /// </summary>
    public IEnumerable<VcfRecord> SplitMultiallelic()
    {
        if (Variant.Alternates.Count <= 1) { yield return this; yield break; }
        for (int a = 1; a <= Variant.Alternates.Count; a++)
        {
            var v = new Variant(Variant.Chromosome, Variant.Position, Variant.Reference, [Variant.Alternates[a - 1]], Variant.Id)
            {
                Quality = Variant.Quality,
                Filter = Variant.Filter,
            };
            foreach (var kv in Variant.Info) v.Info[kv.Key] = kv.Value;
            foreach (var g in Variant.Genotypes)
            {
                var al = new int[g.Alleles.Length];
                for (int i = 0; i < al.Length; i++) al[i] = g.Alleles[i] < 0 ? -1 : g.Alleles[i] == a ? 1 : 0;
                v.Genotypes.Add(new Genotype(al, g.Phased));
            }
            var rec = new VcfRecord(v, _format.Where(f => f != "GT"));
            for (int s = 0; s < SampleCount; s++)
            {
                foreach (var key in _format)
                {
                    if (key == "GT") continue;
                    string? raw = GetField(s, key);
                    if (raw is null) { rec.SetField(s, key, "."); continue; }
                    if (key == "AD" && IntList(raw) is { } ad && ad.Length > a)
                        rec.SetField(s, key, [ad[0], ad[a]]);
                    else if (key == "PL" && IntList(raw) is { } pl && pl.Length > GenotypeIndex(a, a))
                        rec.SetField(s, key, [pl[0], pl[GenotypeIndex(0, a)], pl[GenotypeIndex(a, a)]]);
                    else
                        rec.SetField(s, key, raw);
                }
            }
            yield return rec;
        }
    }

    internal List<string[]> Values()
    {
        if (_values is not null) return _values;
        _values = new List<string[]>(_rawSampleCount);
        if (_rawSamples is not null)
        {
            var s = _rawSamples.AsSpan();
            while (true)
            {
                int tab = s.IndexOf('\t');
                var col = tab < 0 ? s : s[..tab];
                _values.Add(SplitColon(col));
                if (tab < 0) break;
                s = s[(tab + 1)..];
            }
            _rawSamples = null;
        }
        return _values;
    }

    private static string[] SplitColon(ReadOnlySpan<char> col)
    {
        int n = col.Count(':') + 1;
        var parts = new string[n];
        for (int i = 0; i < n; i++)
        {
            int c = col.IndexOf(':');
            var p = c < 0 ? col : col[..c];
            parts[i] = p.Length == 1 && p[0] == '.' ? "." : p.ToString();
            if (c >= 0) col = col[(c + 1)..];
        }
        return parts;
    }

    internal bool HasUnparsedSamples => _rawSamples is not null;

    internal string? RawSamples => _rawSamples;

    private static int? Int(string? s) => s is not null && int.TryParse(s, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) ? v : null;

    private static int[]? IntList(string? s)
    {
        if (s is null) return null;
        var span = s.AsSpan();
        var result = new int[span.Count(',') + 1];
        for (int i = 0; i < result.Length; i++)
        {
            int c = span.IndexOf(',');
            var p = c < 0 ? span : span[..c];
            result[i] = int.TryParse(p, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) ? v : -1;
            if (c >= 0) span = span[(c + 1)..];
        }
        return result;
    }

    /// <summary>The record as a VCF data line (no trailing newline).</summary>
    public override string ToString()
    {
        var sw = new StringWriter();
        VcfWriter.WriteRecord(sw, this);
        return sw.ToString();
    }
}
