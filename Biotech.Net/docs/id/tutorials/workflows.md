# Tutorial: Workflow

Paket: `Gravicode.Biotech.Net.Workflows`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Mesin workflow bergaya Snakemake dan Nextflow yang ditulis sebagai C# biasa: tanpa DSL yang harus dipelajari dan tanpa parser di antara Anda dan kode. Anda mendeklarasikan **rule** dengan input dan output bernama (berkas, atau nilai in-memory bertipe); mesin menghubungkannya menjadi DAG melalui data yang mereka pertukarkan, menjalankan job yang siap secara paralel dalam batas thread, memori, dan sumber daya kustom, mengulang dan memberi batas waktu pada langkah yang gagal, dan pada run berikutnya melewati setiap job yang input, parameter, dan versi kodenya tidak berubah. Langkah bisa berupa lambda C# yang memanggil API Biotech.Net apa pun, perintah shell, atau program eksternal.

Paket ini hanya bergantung pada `Biotech.Net.Core`, sehingga bisa menggerakkan modul mana pun.

## 1. Persiapan

```text
dotnet add package Gravicode.Biotech.Net.Workflows
```

Sebuah `Workflow` berakar pada direktori kerja. Path relatif dihitung dari direktori itu, langkah shell dan proses berjalan di dalamnya, dan state-nya disimpan di `<direktori kerja>/.biotech-workflow/`: satu rekaman per job yang selesai, nilai tersimpan, cache hash, log per job, dan satu folder laporan per run.

## 2. Workflow pertama

Dua rule: `upper` membaca `data/greeting.txt` dan menulis `work/greeting.upper.txt`; `count` membaca berkas itu. Tidak ada yang menyatakan `count` berjalan setelah `upper` — mesin mengetahuinya karena output yang satu adalah input yang lain.

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

`ctx.Input(name)` adalah path lengkap sebuah input. `ctx.Output(name)` adalah path **sementara**: mesin memindahkannya ke `work/greeting.upper.txt` hanya setelah langkah selesai tanpa galat, sehingga crash tidak pernah meninggalkan berkas setengah tertulis yang kelak dikira hasil. `RunAsync` tidak melempar exception saat job gagal atau run dibatalkan; periksa `report.Status`, atau panggil `report.EnsureSuccess()`.

## 3. Wildcard, sample sheet, scatter dan gather

Pola seperti `counts/{sample}.txt` menghasilkan satu job per nilai `{sample}`. Job dicari mundur dari target, seperti di Snakemake: target `results/words.tsv` membutuhkan setiap `counts/{sample}.txt` (daftarnya ditulis oleh `Expand`), dan masing-masing cocok dengan rule `words` dengan `sample` terikat. Rule yang memakai wildcard adalah sisi map (scatter); rule yang inputnya berupa daftar `Expand` adalah sisi collect (gather).

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-sheet-" + Guid.NewGuid().ToString("N")[..8]);
Directory.CreateDirectory(Path.Combine(dir, "data"));
File.WriteAllText(Path.Combine(dir, "data", "liver.txt"), "albumin transferrin fibrinogen");
File.WriteAllText(Path.Combine(dir, "data", "brain.txt"), "synapsin tubulin");
File.WriteAllText(Path.Combine(dir, "samples.tsv"), "sample\tfile\nliver\tdata/liver.txt\nbrain\tdata/brain.txt\n");

var samples = SampleSheet.Load(Path.Combine(dir, "samples.tsv"));
var wf = new Workflow(dir, "scatter-gather");
wf.Rule("words")                                            // map: satu job per sampel
    .Input("text", w => samples.Find("sample", w["sample"])["file"])
    .Output("n", "counts/{sample}.txt")
    .Run(ctx => File.WriteAllText(ctx.Output("n"), File.ReadAllText(ctx.Input("text")).Split(' ').Length.ToString()));
wf.Rule("table")                                            // collect: satu job untuk semua sampel
    .Input("counts", Expand.Rows("counts/{sample}.txt", samples))
    .Output("tsv", "results/words.tsv")
    .Run(ctx => File.WriteAllLines(ctx.Output("tsv"),
        ctx.Inputs("counts").Select(p => $"{Path.GetFileNameWithoutExtension(p)}\t{File.ReadAllText(p)}")));
wf.Target("results/words.tsv");

Console.Write(wf.Plan().ToText());       // dry run: apa yang akan dijalankan, dan alasannya
Console.Write(wf.Build().ToMermaid());   // DAG-nya, untuk halaman Markdown
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

Wildcard cocok dengan satu karakter atau lebih selain `/`. Beri ekspresi reguler sendiri dengan `{name,regex}`, misalnya `chunks/{sample}.{n,\d+}.fastq`. `Expand.Product` membangun hasil kali silang (`("sample", samples), ("chunk", chunks)`), yaitu cara menulis scatter ke beberapa potongan lalu gather kembali. `Input(name, w => ...)` menghitung input dari wildcard, yaitu cara kolom sample sheet menjadi path.

## 4. Nilai in-memory bertipe

Tidak semua data antarlangkah harus berupa berkas. `wf.Value<T>(name, pattern)` mendeklarasikan nilai bertipe; sebuah rule mengeluarkannya dengan `ctx.Set` dan rule lain membacanya dengan `ctx.Get`, atau mengumpulkan satu per sampel dengan `ctx.GetAll`. Nilai disimpan sebagai JSON di direktori state, sehingga nilai dari produsen yang dilewati tetap tersedia bagi konsumennya, dan perubahan sebuah nilai menjalankan ulang apa pun yang membacanya.

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

## 5. Langkah shell dan proses

`Shell` menjalankan baris perintah melalui `cmd /c` di Windows dan `/bin/sh -c` di sistem lain; `Process` memulai program dengan daftar argumen tanpa shell. Keduanya menerima placeholder: `{input}`, `{input.NAME}`, `{output}`, `{output.NAME}` (path sementara), `{wildcards.NAME}`, `{threads}`, dan `{log}`; `{{` dan `}}` adalah kurung kurawal literal. Path yang mengandung spasi diberi tanda kutip. Standard output dan standard error masuk ke log job baris per baris, exit code bukan nol menggagalkan job, dan perintah yang sudah dirender menjadi bagian dari signature job, sehingga mengubahnya menjalankan ulang rule. Di dalam langkah C#, `ctx.RunProcessAsync` dan `ctx.RunShellAsync` melakukan hal yang sama.

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

## 6. Kapan sebuah job dianggap usang

Ketika semua dependensi sebuah job selesai, mesin membandingkan job itu dengan rekaman keberhasilan terakhirnya dan menjalankannya bila salah satu kondisi berikut berlaku (kondisi pertama yang berlaku menjadi alasan yang dilaporkan):

| Alasan | Arti |
|---|---|
| `forced (...)` | `RunOptions.ForceAll`, atau rule-nya ada di `ForceRules` |
| `incomplete: ...` | run sebelumnya memulai job ini dan tidak menyelesaikannya (gagal, crash, dibatalkan) |
| `missing output X` | sebuah berkas atau nilai output tidak ada |
| `no record of a previous run` | job ini belum pernah berhasil di sini |
| `rule definition changed (...)` | parameternya (`Params`, diserialisasi ke JSON), `Version`, perintah shell, atau daftar inputnya berubah |
| `input changed: X` | hash konten sebuah input berbeda dari yang tercatat |
| `output modified outside the workflow: X` | seseorang menyunting hasil secara manual |

Input diberi sidik jari SHA-256 dari kontennya (`HashMode.Content`, di-cache berdasarkan ukuran dan waktu modifikasi), sehingga menyentuh berkas tanpa mengubahnya tidak menjalankan ulang apa pun, dan job yang upstream-nya dijalankan ulang tetapi menghasilkan output yang identik byte demi byte akan dilewati. `HashMode.Timestamp` hanya memakai ukuran dan waktu modifikasi. Badan C# tidak bisa di-hash: naikkan `.Version("2")` ketika kode sebuah langkah berubah.

`wf.Plan(options)` adalah dry run. Job yang upstream-nya akan berjalan dicantumkan dengan `upstream ... will run`; apakah job itu benar-benar berjalan diputuskan saat ia siap, dari konten aktual inputnya.

## 7. Paralelisme, sumber daya, pengulangan, batas waktu, dan pembatalan

```csharp
using Biotech.Net.Workflows;

string dir = Path.Combine(Path.GetTempPath(), "wf-limits-" + Guid.NewGuid().ToString("N")[..8]);
var wf = new Workflow(dir, "limits");
int calls = 0;
wf.Rule("align").Output("bam", "aln/{chunk}.txt")
    .Threads(4).Memory(2_000).Resource("gpu", 1)            // dihitung terhadap batas RunOptions
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
    MaxParallelism = 8,                                    // job sekaligus
    MaxThreads = 8,                                        // dua job 4-thread sekaligus
    MaxMemoryMb = 16_000,
    Resources = new Dictionary<string, int> { ["gpu"] = 1 }, // ...tetapi hanya satu GPU, jadi satu per satu
    KeepGoing = true,                                      // tetap jalankan job independen setelah ada kegagalan
};
using var cts = new CancellationTokenSource(TimeSpan.FromHours(1));
WorkflowRunReport? report = null;
await foreach (var e in wf.RunStreamingAsync(options, cts.Token))   // atau RunAsync(options, IProgress<WorkflowEvent>, token)
{
    if (e is JobRetryingEvent r) Console.WriteLine($"retry {r.JobId} in {r.Delay.TotalMilliseconds} ms: {r.Error}");
    if (e is RunCompletedEvent done) report = done.Report;
}

Console.WriteLine($"{report!.Status}, attempts: {string.Join(", ", report.Jobs.Select(j => $"{j.Id}={j.Attempts}"))}");
// retry align[chunk=1] in 200 ms: IOException: transient network error
// Succeeded, attempts: align[chunk=1]=2, align[chunk=2]=1, align[chunk=3]=1
```

Permintaan yang lebih besar dari batas dipotong ke batas itu (rule 64 thread pada run 8 thread mendapat 8), sehingga job tetap berjalan sendirian. Pengulangan ke-`n` menunggu `delay × backoffFactor^(n-1)`. Percobaan yang melewati batas waktu dihitung sebagai kegagalan dan dapat diulang. Membatalkan token menghentikan job yang sedang berjalan (proses eksternal dimatikan beserta proses anaknya), menghapus output sementaranya, menandainya incomplete, dan mengembalikan laporan dengan `RunStatus.Cancelled`. Dua run pada direktori kerja yang sama secara bersamaan ditolak. Jeda dan batas waktu memakai `RunOptions.TimeProvider`, yang di pengujian diganti jam palsu.

## 8. Event, log, laporan, dan graf

- **Event.** `RunAsync(options, IProgress<WorkflowEvent>)` melaporkan `RunStartedEvent`, `JobStartedEvent`, `JobLogEvent`, `JobSkippedEvent`, `JobRetryingEvent`, `JobSucceededEvent`, `JobFailedEvent`, `JobCancelledEvent`, `RunProgressEvent`, dan terakhir `RunCompletedEvent` beserta laporannya. `await foreach (var e in wf.RunStreamingAsync(options))` menghasilkan event yang sama sebagai `IAsyncEnumerable`.
- **Log.** Setiap job menulis `.biotech-workflow/logs/<run id>/<job>.log`: baris percobaan bertanda waktu, semua yang ditulis `ctx.Log`, serta stdout/stderr tool dengan awalan `[out]`/`[err]`.
- **Laporan.** Setiap run menulis `report.json`, `report.md`, `report.html`, dan `dag.dot` (diwarnai menurut hasil) ke `.biotech-workflow/runs/<run id>/`, ditambah `last-run.json`; `wf.LastReport()` membacanya kembali. `WorkflowRunReport` memiliki `ToJson`, `FromJson`, `ToMarkdown`, dan `ToHtml`.
- **Graf.** `wf.Build().ToDot()` dan `ToMermaid()` mengekspor DAG job; `rules: true` meringkasnya menjadi satu simpul per rule; `wf.Plan().ToDot()` mewarnai job yang akan berjalan.

## 9. Ujung ke ujung: pemanggilan varian pada tiga sampel

Kontrol kualitas FASTQ, pemetaan read, pemanggilan SNV, dan penyaringan VCF per sampel, lalu tabel ringkasan atas semuanya — hanya memakai API Biotech.Net. Statistik QC berpindah sebagai nilai bertipe; selebihnya berupa berkas. Kode yang sama menjadi pengujian di `tests/Biotech.Net.Workflows.Tests/ExampleWorkflowTests.cs`.

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

// Read baru hanya untuk sampel B: empat job B dan summary dijalankan ulang, A dan C dilewati.
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

        // Sisi map (satu job per sampel): pangkas read berdasarkan kualitas dan catat statistik QC.
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

        // Sisi collect: satu job atas VCF dan nilai QC setiap sampel.
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

    /// <summary>Menulis referensi acak 3 kb, dan untuk setiap sampel read 100 bp (20x) yang membawa SNV sampel itu.</summary>
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
                // Basa Q40 dengan ekor 3' berkualitas rendah untuk dipangkas langkah QC.
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

`report.html` dari run ini mencantumkan 13 job beserta waktunya; `dag.dot` di sebelahnya merender DAG dengan Graphviz (`dot -Tsvg dag.dot -o dag.svg`).

## Batasan

- **Satu mesin.** Job berjalan sebagai task dan proses lokal; tidak ada eksekutor klaster, cloud, atau container (Slurm, Kubernetes, Docker), dan tidak ada sistem berkas jarak jauh atau object store sebagai input.
- **Tanpa DSL dan tanpa pewarisan rule.** Workflow adalah objek C#; tidak ada parser Snakefile atau skrip Nextflow, tidak ada `ruleorder` (dua rule yang dapat menghasilkan berkas yang sama dianggap galat), tidak ada checkpoint yang mengubah DAG selama run, dan tidak ada channel streaming: data mengalir lewat berkas atau nilai tersimpan setelah sebuah job selesai.
- **Badan C# tidak di-hash.** Mengubah kode lambda `Run` tidak menjalankan ulang apa pun sampai Anda menaikkan `.Version(...)`; parameter, perintah shell, dan daftar input di-hash secara otomatis.
- **Nilai harus bisa diserialisasi JSON** dengan `System.Text.Json` dan disimpan utuh di memori dan di disk; gunakan berkas untuk data besar.
- **Output atomik berlaku per berkas.** Setiap output ditulis ke folder sementara di sebelahnya (volume yang sama) lalu di-rename ke tempatnya. Output sebuah job dipindahkan satu per satu, jadi crash di tengah pemindahan meninggalkan job bertanda incomplete: run berikutnya mengulangnya alih-alih memercayai kumpulan yang tercampur.
- **Direktori bukan output.** Deklarasikan berkas-berkas di dalamnya.
- **Portabilitas shell adalah tanggung jawab Anda.** `Shell` menjalankan `cmd.exe` di Windows dan `/bin/sh` di sistem lain; tulis perintah untuk keduanya (seperti contoh di atas) atau gunakan `Process` dengan daftar argumen.

---
Sebelumnya: [Konektor data cloud](data.md) · [Indeks tutorial](index.md)
