# Tutorial: Workflows

Package: `Gravicode.Biotech.Net.Workflows`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

A workflow engine in the spirit of Snakemake and Nextflow, written as plain C#: no DSL to learn and no parser between you and the code. You declare **rules** with named inputs and outputs (files, or typed in-memory values); the engine connects them into a DAG through the data they exchange, runs ready jobs in parallel within thread, memory and custom-resource limits, retries and times out failing steps, and on the next run skips every job whose inputs, parameters and code version are unchanged. Steps can be C# lambdas calling any Biotech.Net API, shell commands or external programs.

The package depends only on `Biotech.Net.Core`, so it can drive any module.

## 1. Set-up

```text
dotnet add package Gravicode.Biotech.Net.Workflows
```

A `Workflow` is rooted at a working directory. Relative paths are relative to it, shell and process steps run in it, and its state lives in `<working directory>/.biotech-workflow/`: one record per finished job, stored values, a hash cache, per-job logs and one report folder per run.

## 2. A first workflow

Two rules: `upper` reads `data/greeting.txt` and writes `work/greeting.upper.txt`; `count` reads that file. Nothing says `count` comes after `upper` — the engine finds out because one's output is the other's input.

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-hello-" + Guid.NewGuid().ToString("N")[..8]);
Directory.CreateDirectory(Path.Combine(dir, "data"));
File.WriteAllText(Path.Combine(dir, "data", "greeting.txt"), "hello workflows");

var wf = new Workflow(dir, "hello");
wf.Rule("upper")
    .Input("text", "data/greeting.txt")
    .Output("upper", "work/greeting.upper.txt")
    .Run(ctx => File.WriteAllText(ctx.Output("upper"), File.ReadAllText(ctx.Input("text")).ToUpperInvariant()));
wf.Rule("count")
    .Input("text", "work/greeting.upper.txt")
    .Output("count", "results/length.txt")
    .Run(ctx => File.WriteAllText(ctx.Output("count"), File.ReadAllText(ctx.Input("text")).Length.ToString()));
wf.Target("results/length.txt");

var first = await wf.RunAsync();
Console.WriteLine($"{first.Status}: {first.Count(JobStatus.Succeeded)} ran");          // Succeeded: 2 ran
var second = await wf.RunAsync();
Console.WriteLine($"{second.Status}: {second.Count(JobStatus.Skipped)} up to date");   // Succeeded: 2 up to date

File.WriteAllText(Path.Combine(dir, "data", "greeting.txt"), "hello again");
var third = await wf.RunAsync();
foreach (var job in third.Jobs) Console.WriteLine($"{job.Id}: {job.Status} ({job.Reason})");
// upper: Succeeded (input changed: data/greeting.txt)
// count: Succeeded (input changed: work/greeting.upper.txt)
```

`ctx.Input(name)` is the full path of an input. `ctx.Output(name)` is a **temporary** path: the engine moves it to `work/greeting.upper.txt` only after the step returns without error, so a crash never leaves a half-written file that a later run would mistake for a result. `RunAsync` does not throw when a job fails or the run is cancelled; check `report.Status`, or call `report.EnsureSuccess()`.

## 3. Wildcards, sample sheets, scatter and gather

A pattern such as `counts/{sample}.txt` makes one job per value of `{sample}`. Jobs are found backwards from the targets, as in Snakemake: the target `results/words.tsv` needs every `counts/{sample}.txt` (`Expand` writes the list), and each of those matches the `words` rule with `sample` bound. The rule with wildcards is the map (scatter) side; the rule whose input is an `Expand` list is the collect (gather) side.

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-sheet-" + Guid.NewGuid().ToString("N")[..8]);
Directory.CreateDirectory(Path.Combine(dir, "data"));
File.WriteAllText(Path.Combine(dir, "data", "liver.txt"), "albumin transferrin fibrinogen");
File.WriteAllText(Path.Combine(dir, "data", "brain.txt"), "synapsin tubulin");
File.WriteAllText(Path.Combine(dir, "samples.tsv"), "sample\tfile\nliver\tdata/liver.txt\nbrain\tdata/brain.txt\n");

var samples = SampleSheet.Load(Path.Combine(dir, "samples.tsv"));
var wf = new Workflow(dir, "scatter-gather");
wf.Rule("words")                                            // map: one job per sample
    .Input("text", w => samples.Find("sample", w["sample"])["file"])
    .Output("n", "counts/{sample}.txt")
    .Run(ctx => File.WriteAllText(ctx.Output("n"), File.ReadAllText(ctx.Input("text")).Split(' ').Length.ToString()));
wf.Rule("table")                                            // collect: one job over all samples
    .Input("counts", Expand.Rows("counts/{sample}.txt", samples))
    .Output("tsv", "results/words.tsv")
    .Run(ctx => File.WriteAllLines(ctx.Output("tsv"),
        ctx.Inputs("counts").Select(p => $"{Path.GetFileNameWithoutExtension(p)}\t{File.ReadAllText(p)}")));
wf.Target("results/words.tsv");

Console.Write(wf.Plan().ToText());       // a dry run: what would run, and why
Console.Write(wf.Build().ToMermaid());   // the DAG, for a Markdown page
await wf.RunAsync();
Console.Write(File.ReadAllText(Path.Combine(dir, "results", "words.tsv")));
```

```text
run   words[sample=liver]  # missing output counts/liver.txt
run   words[sample=brain]  # missing output counts/brain.txt
run   table  # missing output results/words.tsv

Job counts:
      2 to run,     0 up to date  words
      1 to run,     0 up to date  table
      3 to run,     0 up to date  total
flowchart LR
  n0["words<br/>sample=liver"]
  n1["words<br/>sample=brain"]
  n2["table"]
  n0 --> n2
  n1 --> n2
liver	3
brain	2
```

Wildcards match one or more characters other than `/`. Give one a regular expression with `{name,regex}`, for example `chunks/{sample}.{n,\d+}.fastq`. `Expand.Product` builds cross products (`("sample", samples), ("chunk", chunks)`), which is how a scatter into chunks and a gather back are written. `Input(name, w => ...)` computes an input from the wildcards, which is how a sample sheet's columns become paths.

## 4. Typed in-memory values

Not everything between steps needs to be a file. `wf.Value<T>(name, pattern)` declares a typed value; a rule outputs it with `ctx.Set` and another reads it with `ctx.Get`, or gathers one per sample with `ctx.GetAll`. Values are stored as JSON in the state directory, so a skipped producer's value is still there for its consumers, and a change in a value reruns what reads it.

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-values-" + Guid.NewGuid().ToString("N")[..8]);
Directory.CreateDirectory(Path.Combine(dir, "data"));
File.WriteAllText(Path.Combine(dir, "data", "s1.seq"), "ATGCGCGTA");
File.WriteAllText(Path.Combine(dir, "data", "s2.seq"), "ATATATGC");

var samples = SampleSheet.FromColumn("sample", "s1", "s2");
var wf = new Workflow(dir, "values");
var gc = wf.Value<GcStats>("gc", "gc/{sample}");
wf.Rule("gc").Input("seq", "data/{sample}.seq").Output(gc).Run(ctx =>
{
    string s = File.ReadAllText(ctx.Input("seq"));
    ctx.Set(gc, new GcStats(ctx.Wildcards["sample"], s.Length, (double)s.Count(c => c is 'G' or 'C') / s.Length));
});
wf.Rule("report").Input(gc, samples).Output("tsv", "results/gc.tsv").Run(ctx =>
    File.WriteAllLines(ctx.Output("tsv"), ctx.GetAll(gc).Select(g => $"{g.Sample}\t{g.Length}\t{g.Gc:F3}")));

var report = await wf.RunAsync(new RunOptions { Targets = ["results/gc.tsv"] });
Console.WriteLine(report["gc[sample=s1]"].Outputs[0]);                        // value:gc/s1
Console.Write(File.ReadAllText(Path.Combine(dir, "results", "gc.tsv")));
// s1	9	0.556
// s2	8	0.250

record GcStats(string Sample, int Length, double Gc);
```

## 5. Shell and process steps

`Shell` runs a command line through `cmd /c` on Windows and `/bin/sh -c` elsewhere; `Process` starts a program with an argument list and no shell. Both take placeholders: `{input}`, `{input.NAME}`, `{output}`, `{output.NAME}` (the temporary paths), `{wildcards.NAME}`, `{threads}` and `{log}`; `{{` and `}}` are literal braces. Paths that contain spaces are quoted. Standard output and standard error go to the job's log line by line, a non-zero exit code fails the job, and the rendered command is part of the job's signature, so editing it reruns the rule. Inside a C# step, `ctx.RunProcessAsync` and `ctx.RunShellAsync` do the same.

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-shell-" + Guid.NewGuid().ToString("N")[..8]);
Directory.CreateDirectory(Path.Combine(dir, "data"));
File.WriteAllText(Path.Combine(dir, "data", "reads.txt"), "ACGT\n");

var wf = new Workflow(dir, "shell");
string cat = OperatingSystem.IsWindows() ? "type" : "cat";
wf.Rule("copy").Input("src", "data/{name}.txt").Output("dst", "out/{name}.copy").Shell(cat + " {input.src} > {output.dst}");
wf.Rule("sdk").Output("v", "out/sdk.txt").Shell("dotnet --version > {output.v}");
wf.Rule("fail").Output("x", "out/never.txt").Shell("exit 3");

var report = await wf.RunAsync(new RunOptions { Targets = ["out/reads.copy", "out/sdk.txt", "out/never.txt"], KeepGoing = true });
Console.WriteLine(File.ReadAllText(Path.Combine(dir, "out", "reads.copy")).Trim());   // ACGT
Console.WriteLine(report["fail"].Error);       // Windows: 'cmd.exe /d /s /c "exit 3"' exited with code 3.  Linux: '/bin/sh -c exit 3' ...
Console.WriteLine(report["fail"].ExitCode);                                          // 3
Console.WriteLine(File.Exists(Path.Combine(dir, "out", "never.txt")));               // False
```

## 6. What makes a job out of date

When all of a job's dependencies have finished, the engine compares the job with the record of its last success and runs it if any of these holds (the first that applies is the reported reason):

| Reason | Meaning |
|---|---|
| `forced (...)` | `RunOptions.ForceAll`, or the rule is in `ForceRules` |
| `incomplete: ...` | a previous run started the job and did not finish it (failure, crash, cancellation) |
| `missing output X` | an output file or value is absent |
| `no record of a previous run` | the job has never succeeded here |
| `rule definition changed (...)` | its parameters (`Params`, serialized to JSON), `Version`, shell command or input list changed |
| `input changed: X` | an input's content hash differs from the one recorded |
| `output modified outside the workflow: X` | someone edited a result by hand |

Inputs are fingerprinted by SHA-256 of their content (`HashMode.Content`, cached by size and modification time), so touching a file without changing it reruns nothing, and a job whose upstream reran but produced byte-identical output is skipped. `HashMode.Timestamp` uses size and modification time only. C# bodies cannot be hashed: bump `.Version("2")` when a step's code changes.

`wf.Plan(options)` is the dry run. A job whose upstream will run is listed with `upstream ... will run`; whether it really runs is decided when it becomes ready, from the actual content of its inputs.

## 7. Parallelism, resources, retries, timeouts and cancellation

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-limits-" + Guid.NewGuid().ToString("N")[..8]);
var wf = new Workflow(dir, "limits");
int calls = 0;
wf.Rule("align").Output("bam", "aln/{chunk}.txt")
    .Threads(4).Memory(2_000).Resource("gpu", 1)            // counted against the RunOptions limits
    .Retries(2, delay: TimeSpan.FromMilliseconds(200), backoffFactor: 2)
    .Timeout(TimeSpan.FromMinutes(30))
    .Run(async ctx =>
    {
        if (Interlocked.Increment(ref calls) == 1) throw new IOException("transient network error");
        await File.WriteAllTextAsync(ctx.Output("bam"), $"{ctx.Wildcards["chunk"]} on {ctx.Threads} threads", ctx.CancellationToken);
    });

var options = new RunOptions
{
    Targets = Expand.Each("aln/{chunk}.txt", "chunk", ["1", "2", "3"]),
    MaxParallelism = 8,                                    // jobs at once
    MaxThreads = 8,                                        // two 4-thread jobs at a time
    MaxMemoryMb = 16_000,
    Resources = new Dictionary<string, int> { ["gpu"] = 1 }, // ...but only one GPU, so one at a time
    KeepGoing = true,                                      // keep running independent jobs after a failure
};
using var cts = new CancellationTokenSource(TimeSpan.FromHours(1));
WorkflowRunReport? report = null;
await foreach (var e in wf.RunStreamingAsync(options, cts.Token))   // or RunAsync(options, IProgress<WorkflowEvent>, token)
{
    if (e is JobRetryingEvent r) Console.WriteLine($"retry {r.JobId} in {r.Delay.TotalMilliseconds} ms: {r.Error}");
    if (e is RunCompletedEvent done) report = done.Report;
}

Console.WriteLine($"{report!.Status}, attempts: {string.Join(", ", report.Jobs.Select(j => $"{j.Id}={j.Attempts}"))}");
// retry align[chunk=1] in 200 ms: IOException: transient network error
// Succeeded, attempts: align[chunk=1]=2, align[chunk=2]=1, align[chunk=3]=1
```

A request larger than the limit is capped to it (a 64-thread rule on an 8-thread run gets 8), so the job still runs alone. Retry `n` waits `delay × backoffFactor^(n-1)`. A timed-out attempt counts as a failure and may be retried. Cancelling the token stops running jobs (external processes are killed with their child processes), deletes their temporary outputs, leaves them marked incomplete, and returns a report with `RunStatus.Cancelled`. Two runs of the same working directory at once are refused. Delays and timeouts use `RunOptions.TimeProvider`, which tests replace with a fake clock.

## 8. Events, logs, reports and graphs

- **Events.** `RunAsync(options, IProgress<WorkflowEvent>)` reports `RunStartedEvent`, `JobStartedEvent`, `JobLogEvent`, `JobSkippedEvent`, `JobRetryingEvent`, `JobSucceededEvent`, `JobFailedEvent`, `JobCancelledEvent`, `RunProgressEvent` and finally `RunCompletedEvent` with the report. `await foreach (var e in wf.RunStreamingAsync(options))` yields the same events as an `IAsyncEnumerable`.
- **Logs.** Each job writes `.biotech-workflow/logs/<run id>/<job>.log`: timestamped attempt lines, everything `ctx.Log` writes, and a tool's stdout/stderr prefixed `[out]`/`[err]`.
- **Reports.** Every run writes `report.json`, `report.md`, `report.html` and `dag.dot` (coloured by outcome) to `.biotech-workflow/runs/<run id>/`, plus `last-run.json`; `wf.LastReport()` reads it back. `WorkflowRunReport` has `ToJson`, `FromJson`, `ToMarkdown` and `ToHtml`.
- **Graphs.** `wf.Build().ToDot()` and `ToMermaid()` export the job DAG; `rules: true` collapses it to one node per rule; `wf.Plan().ToDot()` colours the jobs that would run.

## 9. End to end: variant calling over three samples

FASTQ quality control, read mapping, SNV calling and VCF filtering per sample, then a summary table over all of them — using only Biotech.Net APIs. The QC statistics travel as a typed value; everything else is files. The same code is a test in `tests/Biotech.Net.Workflows.Tests/ExampleWorkflowTests.cs`.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Calling;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Genomics.Vcf;
using Biotech.Net.Sequences.IO;
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-variants-" + Guid.NewGuid().ToString("N")[..8]);
var samples = VariantWorkflow.SimulateData(dir, new Dictionary<string, int[]>
{
    ["A"] = [500, 1500],
    ["B"] = [800],
    ["C"] = [],
});

var report = await VariantWorkflow.Create(dir, samples).RunAsync(new RunOptions { MaxParallelism = 3 });
Console.WriteLine($"{report.Status}: {report.Count(JobStatus.Succeeded)} jobs ran");   // Succeeded: 13 jobs ran
Console.Write(File.ReadAllText(Path.Combine(dir, "results", "summary.tsv")));
// sample	reads_in	reads_kept	mean_q	snvs
// A	581	581	40.00	2
// B	581	581	40.00	1
// C	581	581	40.00	0

// New reads for sample B only: B's four jobs and the summary rerun, A and C are skipped.
VariantWorkflow.SimulateData(dir, new Dictionary<string, int[]> { ["B"] = [800, 2200] });
var rerun = await VariantWorkflow.Create(dir, samples).RunAsync();
Console.WriteLine(string.Join(" ", rerun.Jobs.Where(j => j.Status == JobStatus.Succeeded).Select(j => j.Id)));
// qc[sample=B] map[sample=B] call[sample=B] filter[sample=B] summary

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
```

The run's `report.html` lists the 13 jobs with their timings; `dag.dot` beside it renders the DAG with Graphviz (`dot -Tsvg dag.dot -o dag.svg`).

## Limitations

- **One machine.** Jobs run as tasks and local processes; there is no cluster, cloud or container executor (Slurm, Kubernetes, Docker), and no remote file systems or object stores as inputs.
- **No DSL and no rule inheritance.** Workflows are C# objects; there is no Snakefile or Nextflow-script parser, no `ruleorder` (two rules that can produce the same file are an error), no checkpoints that change the DAG during a run, and no streaming channels: data flows by file or by stored value once a job has finished.
- **C# bodies are not hashed.** Changing the code of a `Run` lambda does not rerun anything until you bump `.Version(...)`; parameters, shell commands and input lists are hashed automatically.
- **Values must be JSON-serializable** with `System.Text.Json` and are kept whole in memory and on disk; use files for large data.
- **Atomic outputs are per file.** Each output is written to a temporary folder beside it (same volume) and renamed into place. A job's outputs are moved one after another, so a crash during those moves leaves the job marked incomplete: the next run redoes it rather than trusting a mixed set.
- **Directories are not outputs.** Declare the files inside them.
- **Shell portability is yours.** `Shell` runs `cmd.exe` on Windows and `/bin/sh` elsewhere; write commands for both (as the examples do) or use `Process` with an argument list.

---
Previous: [Cloud data connectors](data.md) · [Tutorial index](index.md)
