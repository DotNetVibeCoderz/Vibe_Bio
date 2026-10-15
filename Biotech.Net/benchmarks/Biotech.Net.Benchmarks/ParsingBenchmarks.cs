using System.Text;
using BenchmarkDotNet.Attributes;
using Biotech.Net.Genomics.Vcf;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Benchmarks;

/// <summary>FASTA, FASTQ and VCF parsing from in-memory text (no disk I/O in the measurement).</summary>
public class ParsingBenchmarks
{
    public const int FastaRecords = 2_000, FastaRecordLength = 5_000;   // 10 Mb of sequence
    public const int FastqReads = 50_000, ReadLength = 150;             // 7.5 Mb of bases
    public const int VcfRecords = 20_000, VcfSamples = 10;

    private string _fasta = "", _fastq = "", _vcf = "";

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(7);
        const string acgt = "ACGT";

        var sb = new StringBuilder();
        for (int r = 0; r < FastaRecords; r++)
        {
            sb.Append(">contig").Append(r).Append(" synthetic record\n");
            for (int i = 0; i < FastaRecordLength; i++)
            {
                sb.Append(acgt[rng.Next(4)]);
                if ((i + 1) % 60 == 0) sb.Append('\n');
            }
            sb.Append('\n');
        }
        _fasta = sb.ToString();

        sb.Clear();
        for (int r = 0; r < FastqReads; r++)
        {
            sb.Append("@read").Append(r).Append('\n');
            for (int i = 0; i < ReadLength; i++) sb.Append(acgt[rng.Next(4)]);
            sb.Append("\n+\n");
            for (int i = 0; i < ReadLength; i++) sb.Append((char)(33 + 20 + rng.Next(21)));
            sb.Append('\n');
        }
        _fastq = sb.ToString();

        sb.Clear();
        sb.Append("##fileformat=VCFv4.2\n");
        sb.Append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
        sb.Append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
        sb.Append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        sb.Append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
        sb.Append("##contig=<ID=chr1,length=248956422>\n");
        sb.Append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int s = 0; s < VcfSamples; s++) sb.Append("\tS").Append(s);
        sb.Append('\n');
        string[] gts = ["0/0", "0/1", "1/1", "0|1", "./."];
        for (int r = 0; r < VcfRecords; r++)
        {
            char refBase = acgt[rng.Next(4)], alt = acgt[(acgt.IndexOf(refBase) + 1 + rng.Next(3)) % 4];
            sb.Append("chr1\t").Append(1000 + r * 50).Append("\trs").Append(r).Append('\t').Append(refBase).Append('\t').Append(alt)
              .Append('\t').Append(30 + rng.Next(70)).Append("\tPASS\tDP=").Append(20 + rng.Next(80)).Append(";AF=0.")
              .Append(rng.Next(1, 99)).Append("\tGT:DP");
            for (int s = 0; s < VcfSamples; s++) sb.Append('\t').Append(gts[rng.Next(gts.Length)]).Append(':').Append(rng.Next(5, 60));
            sb.Append('\n');
        }
        _vcf = sb.ToString();
    }

    [Benchmark(Description = "FASTA 2,000 x 5 kb (10 Mb)")]
    public long Fasta()
    {
        long total = 0;
        foreach (var s in FastaReader.Read(new StringReader(_fasta))) total += s.Length;
        return total;
    }

    [Benchmark(Description = "FASTQ 50,000 x 150 bp")]
    public long Fastq()
    {
        long total = 0;
        foreach (var r in FastqReader.Read(new StringReader(_fastq), FastqEncoding.Phred33)) total += r.Length;
        return total;
    }

    [Benchmark(Description = "VCF 20,000 records x 10 samples")]
    public int Vcf()
    {
        using var reader = new VcfReader(new StringReader(_vcf));
        int n = 0;
        foreach (var rec in reader.Records()) n += rec.Variant.Genotypes.Count;
        return n;
    }
}
