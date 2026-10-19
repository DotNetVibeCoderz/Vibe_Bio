using BenchmarkDotNet.Attributes;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Genomics.Mapping;

namespace Biotech.Net.Benchmarks;

/// <summary>BAM decode throughput from an in-memory file: raw BGZF inflate, and full record decoding.</summary>
[MemoryDiagnoser]
public class BamBenchmarks
{
    public const int Records = 200_000, ReadLength = 150;

    private byte[] _bam = [];

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(11);
        var header = SamHeader.Create([("chr1", 50_000_000L)], sortOrder: "coordinate");
        var ms = new MemoryStream();
        using (var w = new BamWriter(ms, header, leaveOpen: true))
        {
            int pos = 0;
            var seq = new byte[ReadLength];
            var qual = new byte[ReadLength];
            for (int i = 0; i < Records; i++)
            {
                pos += rng.Next(0, 200);
                for (int k = 0; k < ReadLength; k++) { seq[k] = (byte)"ACGT"[rng.Next(4)]; qual[k] = (byte)rng.Next(2, 41); }
                var r = new AlignmentRecord
                {
                    Name = "read" + i, Flags = SamFlags.Paired | SamFlags.ProperPair, ReferenceName = "chr1", Position = pos, MapQ = 60,
                    Cigar = [new CigarElement(ReadLength, 'M')], MateReferenceName = "chr1", MatePosition = pos + 200, TemplateLength = 350,
                    Sequence = (byte[])seq.Clone(), Qualities = (byte[])qual.Clone(),
                };
                r.SetTags([AuxField.Int("NM", rng.Next(0, 4)), AuxField.String("RG", "grp1"), AuxField.Int("AS", 140)]);
                w.Write(r);
            }
        }
        _bam = ms.ToArray();
        Console.WriteLine($"BAM: {Records:N0} records, {_bam.Length / 1e6:F1} MB compressed");
    }

    [Benchmark(Description = "BGZF inflate (whole BAM)")]
    public long BgzfInflate()
    {
        using var bg = new BgzfReader(new MemoryStream(_bam));
        var buf = new byte[1 << 16];
        long n = 0;
        int k;
        while ((k = bg.Read(buf)) > 0) n += k;
        return n;
    }

    [Benchmark(Description = "BAM decode to AlignmentRecord")]
    public long BamDecode()
    {
        using var reader = new BamReader(new MemoryStream(_bam));
        long bases = 0;
        foreach (var r in reader.Records()) bases += r.Sequence.Length;
        return bases;
    }
}
