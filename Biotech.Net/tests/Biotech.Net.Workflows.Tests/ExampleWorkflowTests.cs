using Biotech.Net.Core;
using Biotech.Net.Genomics.Calling;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Genomics.Vcf;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Workflows.Tests;

// The variant-calling example from docs/en/tutorials/workflows.md, kept identical to it.
public static class VariantWorkflow
{
    public sealed record QcParams(int MinQuality, int MinLength);

    public sealed record FilterParams(double MinQuality);

    public sealed record QcStats(string Sample, int ReadsIn, int ReadsKept, double MeanQuality);

    public static Workflow Create(string directory, SampleSheet samples)
    {
        var wf = new Workflow(directory, "variant-calling");
        var qc = wf.Value<QcStats>("qc", "qc/{sample}");

        // Map side (one job per sample): quality-trim the reads and record QC statistics.
        wf.Rule("qc").Describe("trim low-quality 3' ends, drop short reads")
            .Input("reads", w => samples.Find("sample", w["sample"])["fastq"])
            .Output("trimmed", "work/{sample}.trimmed.fastq")
            .Output(qc)
            .Params(new QcParams(MinQuality: 20, MinLength: 50))
            .Run(ctx =>
            {
                var p = ctx.Params<QcParams>();
                var reads = FastqReader.Read(ctx.Input("reads")).ToList();
                var kept = reads.Select(r => r.TrimByQuality(p.MinQuality)).Where(r => r.Length >= p.MinLength).ToList();
                FastqWriter.Write(ctx.Output("trimmed"), kept);
                ctx.Set(qc, new QcStats(ctx.Wildcards["sample"], reads.Count, kept.Count,
                    kept.Count == 0 ? 0 : Math.Round(kept.Average(r => r.MeanQuality()), 2)));
            });

        wf.Rule("map").Describe("map reads to the reference (minimizer seeds + banded alignment)")
            .Input("reference", "ref/genome.fa")
            .Input("reads", "work/{sample}.trimmed.fastq")
            .Output("sam", "work/{sample}.sam")
            .Threads(2)
            .Run(ctx =>
            {
                var reference = FastaReader.Read(ctx.Input("reference")).ToList();
                var mapper = new ReadMapper(reference);
                using var sam = SamWriter.Create(ctx.Output("sam"), reference);
                sam.Write(mapper.Map(FastqReader.Read(ctx.Input("reads"))));
            });

        wf.Rule("call").Describe("pile up the alignments and call SNVs")
            .Input("reference", "ref/genome.fa")
            .Input("sam", "work/{sample}.sam")
            .Output("vcf", "calls/{sample}.raw.vcf")
            .Run(ctx =>
            {
                var reference = FastaReader.Read(ctx.Input("reference")).ToList();
                using var alignments = SamAlignmentReader.Open(ctx.Input("sam"));
                var pileup = new Pileup(reference).AddRange(alignments.Records());
                string sample = ctx.Wildcards["sample"];
                VcfWriter.Write(ctx.Output("vcf"), SnvCaller.CreateHeader(reference, sample),
                    SnvCaller.Call(pileup, new SnvCallerOptions { SampleName = sample }));
            });

        wf.Rule("filter").Describe("keep PASS calls with QUAL >= MinQuality")
            .Input("vcf", "calls/{sample}.raw.vcf")
            .Output("vcf", "calls/{sample}.filtered.vcf")
            .Params(new FilterParams(MinQuality: 30))
            .Run(ctx =>
            {
                using var reader = VcfReader.Open(ctx.Input("vcf"));
                var kept = reader.Records().MinQuality(ctx.Params<FilterParams>().MinQuality).PassOnly().ToList();
                VcfWriter.Write(ctx.Output("vcf"), reader.Header, kept);
            });

        // Collect side: one job over every sample's VCF and QC value.
        wf.Rule("summary").Describe("one row per sample")
            .Input("vcfs", Expand.Rows("calls/{sample}.filtered.vcf", samples))
            .Input(qc, samples)
            .Output("table", "results/summary.tsv")
            .Run(ctx =>
            {
                var stats = ctx.GetAll(qc);
                var vcfs = ctx.Inputs("vcfs");
                var lines = new List<string> { "sample\treads_in\treads_kept\tmean_q\tsnvs" };
                for (int i = 0; i < stats.Count; i++)
                {
                    var s = stats[i];
                    lines.Add($"{s.Sample}\t{s.ReadsIn}\t{s.ReadsKept}\t{s.MeanQuality:F2}\t{VcfReader.ReadRecords(vcfs[i]).Count()}");
                }

                File.WriteAllLines(ctx.Output("table"), lines);
            });

        wf.Target("results/summary.tsv");
        return wf;
    }

    /// <summary>Writes a random 3 kb reference, and for each sample 100 bp reads (20x) carrying that sample's SNVs.</summary>
    public static SampleSheet SimulateData(string directory, IReadOnlyDictionary<string, int[]> snvsBySample, int seed = 7)
    {
        var rnd = new Random(seed);
        var genome = new char[3000];
        for (int i = 0; i < genome.Length; i++) genome[i] = "ACGT"[rnd.Next(4)];
        Directory.CreateDirectory(Path.Combine(directory, "ref"));
        Directory.CreateDirectory(Path.Combine(directory, "data"));
        FastaWriter.Write(Path.Combine(directory, "ref", "genome.fa"), [Sequence.Dna(new string(genome), "chr1")]);

        var rows = new List<IReadOnlyDictionary<string, string>>();
        foreach (var (sample, snvs) in snvsBySample)
        {
            var mutated = (char[])genome.Clone();
            foreach (int pos in snvs) mutated[pos] = mutated[pos] == 'A' ? 'G' : 'A';
            string text = new(mutated);
            var reads = new List<SequenceRead>();
            for (int start = 0; start + 100 <= text.Length; start += 5)
            {
                // Q40 bases with a low-quality 3' tail for the QC step to trim.
                string quality = new string('I', 92) + new string('#', 8);
                reads.Add(SequenceRead.FromAscii($"{sample}_{start}", text.Substring(start, 100), quality));
            }

            string fastq = $"data/{sample}.fastq";
            FastqWriter.Write(Path.Combine(directory, fastq), reads);
            rows.Add(new Dictionary<string, string> { ["sample"] = sample, ["fastq"] = fastq });
        }

        return new SampleSheet(rows);
    }
}

public class ExampleWorkflowTests
{
    private static readonly Dictionary<string, int[]> Snvs = new()
    {
        ["A"] = [500, 1500],
        ["B"] = [800],
        ["C"] = [],
    };

    [Fact]
    public async Task VariantCalling_OverThreeSamples_EndToEnd_AndIncremental()
    {
        using var dir = new TempDir();
        var samples = VariantWorkflow.SimulateData(dir.Path, Snvs);

        var plan = VariantWorkflow.Create(dir.Path, samples).Plan();
        Assert.Equal(13, plan.ToRun.Count);                          // 4 rules x 3 samples + summary

        var report = await VariantWorkflow.Create(dir.Path, samples).RunAsync(new RunOptions { MaxParallelism = 3 });
        Assert.True(report.Succeeded, report.ToMarkdown());
        var table = File.ReadAllLines(dir.File("results/summary.tsv"));
        Assert.Equal("sample\treads_in\treads_kept\tmean_q\tsnvs", table[0]);
        Assert.Equal(["A", "B", "C"], table.Skip(1).Select(l => l.Split('\t')[0]));
        Assert.Equal(["2", "1", "0"], table.Skip(1).Select(l => l.Split('\t')[4]));
        Assert.All(table.Skip(1), l => Assert.Equal("581", l.Split('\t')[1]));
        Assert.Equal([501L, 1501L], VcfReader.ReadRecords(dir.File("calls/A.filtered.vcf")).Select(r => r.Position));

        // Nothing changed: every job is skipped.
        var again = await VariantWorkflow.Create(dir.Path, samples).RunAsync();
        Assert.Equal(13, again.Count(JobStatus.Skipped));

        // New reads for sample B only: B's chain and the summary rerun; A and C stay.
        VariantWorkflow.SimulateData(dir.Path, new Dictionary<string, int[]> { ["B"] = [800, 2200] });
        var rerun = await VariantWorkflow.Create(dir.Path, samples).RunAsync();
        Assert.Equal(["qc[sample=B]", "map[sample=B]", "call[sample=B]", "filter[sample=B]", "summary"],
            rerun.Jobs.Where(j => j.Status == JobStatus.Succeeded).Select(j => j.Id));
        Assert.Equal("2", File.ReadAllLines(dir.File("results/summary.tsv"))[2].Split('\t')[4]);
    }
}
