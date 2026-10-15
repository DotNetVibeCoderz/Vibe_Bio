using System.Globalization;
using Biotech.Net.Core;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Genomics.Population;
using Biotech.Net.Genomics.Vcf;

namespace Biotech.Net.Genomics.Calling;

/// <summary>Read and base filters for <see cref="Pileup"/>.</summary>
public sealed record PileupOptions
{
    /// <summary>Bases below this Phred quality are ignored (default 13, as samtools mpileup).</summary>
    public int MinBaseQuality { get; init; } = 13;

    /// <summary>Reads below this mapping quality are ignored (default 20).</summary>
    public int MinMapQ { get; init; } = 20;

    /// <summary>Quality assumed for reads without qualities (default 30).</summary>
    public byte DefaultQuality { get; init; } = 30;
}

/// <summary>Base counts at one reference position.</summary>
/// <param name="Reference">Reference name.</param>
/// <param name="Position">0-based position.</param>
/// <param name="ReferenceBase">Reference base (upper case).</param>
/// <param name="A">Reads showing A.</param>
/// <param name="C">Reads showing C.</param>
/// <param name="G">Reads showing G.</param>
/// <param name="T">Reads showing T.</param>
public readonly record struct PileupColumn(string Reference, int Position, char ReferenceBase, int A, int C, int G, int T)
{
    /// <summary>Total depth over A/C/G/T.</summary>
    public int Depth => A + C + G + T;

    /// <summary>Count of base <paramref name="b"/>.</summary>
    public int this[char b] => char.ToUpperInvariant(b) switch { 'A' => A, 'C' => C, 'G' => G, 'T' => T, _ => 0 };
}

/// <summary>
/// Base-quality-aware pileup of aligned reads. For every position it keeps the reference-supporting depth and the
/// three per-read log-likelihood sums needed for diploid genotype likelihoods; positions where any non-reference base
/// is seen get full per-base accumulators (sparse), so memory stays ~28 bytes per reference base.
/// </summary>
public sealed class Pileup
{
    private const string Bases = "ACGT";
    private readonly Sequence[] _reference;
    private readonly Dictionary<string, int> _refIndex;
    private readonly RefTrack?[] _tracks;

    private sealed class RefTrack(int length)
    {
        public readonly int[] Depth = new int[length];   // reads matching the reference base
        public readonly double[] S1 = new double[length]; // Σ log10(1 − e)
        public readonly double[] S2 = new double[length]; // Σ log10(e / 3)
        public readonly double[] S3 = new double[length]; // Σ log10(1/2 − e/3)
        public readonly Dictionary<int, AltSite> Alt = [];
    }

    internal sealed class AltSite
    {
        public readonly int[] Count = new int[4];
        public readonly double[] S1 = new double[4], S2 = new double[4], S3 = new double[4];
    }

    /// <summary>Creates an empty pileup over <paramref name="reference"/>.</summary>
    public Pileup(IReadOnlyList<Sequence> reference, PileupOptions? options = null)
    {
        _reference = reference.ToArray();
        _refIndex = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < _reference.Length; i++) _refIndex[_reference[i].Id] = i;
        _tracks = new RefTrack?[_reference.Length];
        Options = options ?? new PileupOptions();
    }

    /// <summary>The filters in use.</summary>
    public PileupOptions Options { get; }

    /// <summary>The reference sequences.</summary>
    public IReadOnlyList<Sequence> Reference => _reference;

    /// <summary>Adds a mapped read.</summary>
    public void Add(ReadAlignment alignment)
    {
        if (!alignment.IsMapped) return;
        Add(alignment.ReferenceName, alignment.Position, alignment.Cigar, alignment.Bases, alignment.Qualities, alignment.MapQ);
    }

    /// <summary>Adds a SAM record (unmapped, secondary, QC-fail and duplicate reads are skipped).</summary>
    public void Add(SamRecord record)
    {
        const SamFlags skip = SamFlags.Unmapped | SamFlags.Secondary | SamFlags.QcFail | SamFlags.Duplicate | SamFlags.Supplementary;
        if ((record.Flag & skip) != 0 || record.Sequence == "*") return;
        var bases = System.Text.Encoding.ASCII.GetBytes(record.Sequence);
        byte[]? q = null;
        if (record.Quality != "*")
        {
            q = new byte[record.Quality.Length];
            for (int i = 0; i < q.Length; i++) q[i] = (byte)(record.Quality[i] - 33);
        }
        Add(record.ReferenceName, (int)record.Position - 1, record.CigarElements, bases, q, record.MapQ);
    }

    /// <summary>Adds many alignments.</summary>
    public Pileup AddRange(IEnumerable<ReadAlignment> alignments)
    {
        foreach (var a in alignments) Add(a);
        return this;
    }

    /// <summary>Adds an aligned read given its 0-based position, CIGAR, bases and qualities (reference orientation).</summary>
    public void Add(string referenceName, int position, ReadOnlySpan<CigarElement> cigar, ReadOnlySpan<byte> bases, ReadOnlySpan<byte> qualities, int mapQ)
    {
        if (mapQ < Options.MinMapQ || !_refIndex.TryGetValue(referenceName, out int ri)) return;
        var refSeq = _reference[ri].Span;
        var track = _tracks[ri] ??= new RefTrack(refSeq.Length);
        int r = position, q = 0;
        foreach (var e in cigar)
        {
            if (e.Op is 'M' or '=' or 'X')
            {
                for (int k = 0; k < e.Length; k++, r++, q++)
                {
                    if (r < 0 || r >= refSeq.Length || q >= bases.Length) continue;
                    int b = Code(bases[q]);
                    if (b < 0) continue;
                    byte qual = qualities.IsEmpty ? Options.DefaultQuality : qualities[q];
                    if (qual < Options.MinBaseQuality) continue;
                    double err = Math.Pow(10, -Math.Max((int)qual, 1) / 10.0);
                    double s1 = Math.Log10(1 - err), s2 = Math.Log10(err / 3), s3 = Math.Log10(0.5 - err / 3);
                    if (b == Code(refSeq[r]))
                    {
                        track.Depth[r]++;
                        track.S1[r] += s1; track.S2[r] += s2; track.S3[r] += s3;
                    }
                    else
                    {
                        if (!track.Alt.TryGetValue(r, out var site)) track.Alt[r] = site = new AltSite();
                        site.Count[b]++;
                        site.S1[b] += s1; site.S2[b] += s2; site.S3[b] += s3;
                    }
                }
            }
            else if (e.Op is 'I' or 'S') q += e.Length;
            else if (e.Op is 'D' or 'N') r += e.Length;
        }
    }

    private static int Code(byte b) => (char)(b & 0xDF) switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => -1 };

    /// <summary>The base counts at a 0-based position.</summary>
    public PileupColumn Column(string referenceName, int position)
    {
        int ri = _refIndex[referenceName];
        char refBase = char.ToUpperInvariant((char)_reference[ri][position]);
        Span<int> c = stackalloc int[4];
        if (_tracks[ri] is { } t)
        {
            int rb = Code((byte)refBase);
            if (rb >= 0) c[rb] = t.Depth[position];
            if (t.Alt.TryGetValue(position, out var site)) for (int b = 0; b < 4; b++) c[b] += site.Count[b];
        }
        return new PileupColumn(referenceName, position, refBase, c[0], c[1], c[2], c[3]);
    }

    /// <summary>Depth (A/C/G/T after filtering) at a 0-based position.</summary>
    public int Depth(string referenceName, int position) => Column(referenceName, position).Depth;

    /// <summary>Columns where a non-reference base was observed, in reference/position order.</summary>
    public IEnumerable<PileupColumn> VariableColumns()
    {
        for (int ri = 0; ri < _reference.Length; ri++)
        {
            if (_tracks[ri] is not { } t) continue;
            foreach (int pos in t.Alt.Keys.Order()) yield return Column(_reference[ri].Id, pos);
        }
    }

    internal IEnumerable<(int Ref, int Pos, int RefCode, int RefDepth, double R1, double R2, double R3, AltSite Site)> AltSites()
    {
        for (int ri = 0; ri < _reference.Length; ri++)
        {
            if (_tracks[ri] is not { } t) continue;
            foreach (int pos in t.Alt.Keys.Order())
                yield return (ri, pos, Code(_reference[ri][pos]), t.Depth[pos], t.S1[pos], t.S2[pos], t.S3[pos], t.Alt[pos]);
        }
    }
}

/// <summary>Options for <see cref="SnvCaller"/>.</summary>
public sealed record SnvCallerOptions
{
    /// <summary>Prior probability of a heterozygous site, θ (default 0.001; hom-alt prior θ/2).</summary>
    public double Heterozygosity { get; init; } = 0.001;

    /// <summary>Minimum reads supporting the alternate base (default 2).</summary>
    public int MinAltReads { get; init; } = 2;

    /// <summary>Minimum fraction of reads supporting the alternate base (default 0.1).</summary>
    public double MinAltFraction { get; init; } = 0.1;

    /// <summary>Minimum site QUAL (default 20).</summary>
    public double MinQuality { get; init; } = 20;

    /// <summary>Sample name for the VCF column (default "sample").</summary>
    public string SampleName { get; init; } = "sample";
}

/// <summary>
/// Diploid SNV caller over a <see cref="Pileup"/>, using GATK-style genotype likelihoods: for genotype {a₁, a₂} each read
/// base b with error rate e contributes P(b | G) = ½P(b | a₁) + ½P(b | a₂), P(b | a) = 1 − e if b = a, else e/3.
/// With priors P(het) = θ, P(hom-alt) = θ/2, QUAL is the Phred-scaled posterior probability of hom-ref, GQ the
/// Phred-scaled probability that the called genotype is wrong, and PL the normalized Phred likelihoods (0/0, 0/1, 1/1).
/// </summary>
public static class SnvCaller
{
    /// <summary>Calls SNVs from a pileup.</summary>
    public static IEnumerable<VcfRecord> Call(Pileup pileup, SnvCallerOptions? options = null)
    {
        options ??= new SnvCallerOptions();
        const string bases = "ACGT";
        double theta = options.Heterozygosity;
        double lpRR = Math.Log10(1 - 1.5 * theta), lpRA = Math.Log10(theta), lpAA = Math.Log10(theta / 2);
        foreach (var (ri, pos, refCode, refDepth, r1, r2, r3, site) in pileup.AltSites())
        {
            if (refCode < 0) continue;
            int alt = -1;
            for (int b = 0; b < 4; b++)
                if (b != refCode && (alt < 0 || site.Count[b] > site.Count[alt])) alt = b;
            int altCount = site.Count[alt];
            int depth = refDepth;
            for (int b = 0; b < 4; b++) depth += site.Count[b];
            if (altCount < options.MinAltReads || (double)altCount / depth < options.MinAltFraction) continue;

            // Per-base sums: index the reference-base sums from the track and the others from the site.
            double S1(int b) => b == refCode ? r1 + site.S1[b] : site.S1[b];
            double S2(int b) => b == refCode ? r2 + site.S2[b] : site.S2[b];
            double S3(int b) => b == refCode ? r3 + site.S3[b] : site.S3[b];
            double others = 0;
            for (int b = 0; b < 4; b++) if (b != refCode && b != alt) others += S2(b);
            double lRR = S1(refCode) + S2(alt) + others;
            double lAA = S1(alt) + S2(refCode) + others;
            double lRA = S3(refCode) + S3(alt) + others;

            Span<double> post = [lRR + lpRR, lRA + lpRA, lAA + lpAA];
            double norm = SpecialFunctions.Log10SumExp(post);
            double pRR = Math.Pow(10, post[0] - norm);
            double qual = pRR <= 0 ? 999 : Math.Min(999, -10 * Math.Log10(pRR));
            if (qual < options.MinQuality) continue;
            int gt = post[1] >= post[2] ? (post[1] >= post[0] ? 1 : 0) : (post[2] >= post[0] ? 2 : 0);
            if (gt == 0) continue;
            double pBest = Math.Pow(10, post[gt] - norm);
            int gq = pBest >= 1 ? 99 : (int)Math.Min(99, Math.Round(-10 * Math.Log10(1 - pBest)));
            double maxL = Math.Max(lRR, Math.Max(lRA, lAA));
            int[] pl = [Pl(lRR, maxL), Pl(lRA, maxL), Pl(lAA, maxL)];

            var refSeq = pileup.Reference[ri];
            var v = new Variant(refSeq.Id, pos + 1, bases[refCode].ToString(), bases[alt].ToString())
            {
                Quality = Math.Round(qual, 2),
                Filter = "PASS",
            };
            int ac = gt == 1 ? 1 : 2;
            v.Info["DP"] = depth.ToString(CultureInfo.InvariantCulture);
            v.Info["AC"] = ac.ToString(CultureInfo.InvariantCulture);
            v.Info["AN"] = "2";
            v.Info["AF"] = (ac / 2.0).ToString("0.###", CultureInfo.InvariantCulture);
            v.Genotypes.Add(new Genotype(gt == 1 ? [0, 1] : [1, 1]));
            var rec = new VcfRecord(v, ["GT"]);
            rec.SetField(0, "AD", [depth - altCount - Others(site, refCode, alt), altCount]);
            rec.SetField(0, "DP", depth.ToString(CultureInfo.InvariantCulture));
            rec.SetField(0, "GQ", gq.ToString(CultureInfo.InvariantCulture));
            rec.SetField(0, "PL", pl);
            yield return rec;
        }
    }

    private static int Others(Pileup.AltSite site, int refCode, int alt)
    {
        int n = 0;
        for (int b = 0; b < 4; b++) if (b != refCode && b != alt) n += site.Count[b];
        return n;
    }

    private static int Pl(double l, double max) => (int)Math.Min(int.MaxValue, Math.Round(-10 * (l - max)));

    /// <summary>Builds a pileup from alignments and calls SNVs.</summary>
    public static IReadOnlyList<VcfRecord> Call(IReadOnlyList<Sequence> reference, IEnumerable<ReadAlignment> alignments,
        SnvCallerOptions? options = null, PileupOptions? pileupOptions = null)
        => Call(new Pileup(reference, pileupOptions).AddRange(alignments), options).ToList();

    /// <summary>A VCF header for caller output: standard GT/AD/DP/GQ/PL fields, contigs and the sample column.</summary>
    public static VcfHeader CreateHeader(IReadOnlyList<Sequence> reference, string sampleName = "sample")
    {
        var h = new VcfHeader().AddStandardGenotypeFields();
        h.OtherLines.Add("source=Biotech.Net.Genomics.SnvCaller");
        foreach (var s in reference) h.AddContig(s.Id, s.Length);
        h.Samples.Add(sampleName);
        return h;
    }
}
