# Referensi API: Workflows

Paket `Gravicode.Biotech.Net.Workflows` · namespace `Biotech.Net.Workflows`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Hanya bergantung pada Core. Mesin bergaya Snakemake/Nextflow: rule dengan input/output berkas dan nilai yang dideklarasikan, wildcard, DAG yang dibangun mundur dari target, penjadwalan paralel dengan sumber daya, pengulangan, batas waktu, pembatalan, run ulang inkremental dari hash konten, output atomik, event, dan laporan. Tutorial: [workflows.md](../tutorials/workflows.md).

## Mendefinisikan workflow

| Tipe | Anggota |
|---|---|
| `Workflow` | `new Workflow(workingDirectory = null, name = "workflow")`; `Name`, `WorkingDirectory`, `StateDirectory` (`.biotech-workflow`), `HashMode`, `Rules`, `Targets`; `Rule(name)`, `Value<T>(name, pattern)`, `Target(params string[])`, `Target(IEnumerable<string>)`; `Build(targets)` → `WorkflowGraph`, `Plan(options)` → `WorkflowPlan`, `RunAsync(options, progress, cancellationToken)` → `WorkflowRunReport`, `RunStreamingAsync(options, cancellationToken)` → `IAsyncEnumerable<WorkflowEvent>`, `LastReport()`, `ToDot(rules)`, `ToMermaid(rules)`; `StateDirectoryName` |
| `Rule` (fluent; setiap metode mengembalikan rule-nya) | `Input(name, params string[] patterns)`, `Input(name, IEnumerable<string>)`, `Input(name, Func<Wildcards, string>)`, `Input(name, Func<Wildcards, IEnumerable<string>>)`, `Input<T>(ValueRef<T>)`, `Input<T>(ValueRef<T>, rows)` (gather); `Output(name, pattern)`, `Output<T>(ValueRef<T>)`; `Params(object)`, `Params(Func<Wildcards, object>)`, `Version(string)`, `Describe(string)`; `Threads(n)`, `Memory(mb)`, `Resource(name, amount)`; `Retries(n, delay, backoffFactor)`, `Timeout(limit)`; badan `Run(Func<JobContext, Task>)`, `Run(Action<JobContext>)`, `Shell(command)`, `Process(fileName, params string[] args)`; properti `Name`, `Kind` (`RuleKind`: `Aggregate`, `Code`, `Shell`, `Process`), `CodeVersion`, `ThreadCount`, `MemoryMb`, `Resources`, `RetryCount`, `RetryDelay`, `RetryBackoffFactor`, `TimeLimit`, `WildcardNames`, `Description` |
| `ValueRef<T>` | `Name`, `Pattern`, `Key(wildcards)`; nilai diserialisasi ke JSON di direktori state |
| `WildcardPattern` | `new WildcardPattern("work/{sample}.{n,\\d+}.fq")`; `Text`, `WildcardNames`, `IsConcrete`, `Match(text)` → `Wildcards?`, `Format(values)`; `{name}` cocok dengan `[^/]+`, `{name,regex}` dengan ekspresi kustom, `{{`/`}}` kurung kurawal literal |
| `Wildcards` | kamus baca-saja nama → nilai menurut urutan pola; `Of(("sample", "A"))`, `Empty`, indexer dengan `KeyNotFoundException` yang jelas, `ToString()` = `sample=A,chunk=1` |
| `SampleSheet` | `Parse(text)` (TSV atau CSV dengan header, komentar `#`), `Load(path)`, `FromColumn(column, values)`, `new SampleSheet(rows)`; `Columns`, `Column(name)`, `Find(column, value)`, indexer, enumerable baris |
| `Expand` | `Product(pattern, params (name, values)[])` (hasil kali silang), `Each(pattern, name, values)`, `Rows(pattern, rows)` |
| `ArtifactKey` | `(Kind, Id)`; `ArtifactKind.File` (path relatif terhadap direktori kerja, pemisah `/`) atau `ArtifactKind.Value` (`value:key` dalam target) |

Target berupa path berkas, `value:key`, atau nama rule tanpa wildcard. Tanpa target, `Build` memakai `Target(...)`, atau bila tidak ada, setiap rule tanpa wildcard. Rule tanpa badan adalah agregat (seperti `rule all` di Snakemake) dan tidak boleh mendeklarasikan output.

## DAG

| Tipe | Anggota |
|---|---|
| `WorkflowGraph` | `Jobs` (urutan topologis), `Targets`, `SourceArtifacts` (input yang tidak dihasilkan rule mana pun), indexer menurut id job; `ToDot(rules, status)`, `ToMermaid(rules, status)` |
| `WorkflowJob` | `Id` (`rule` atau `rule[sample=A]`), `Rule`, `Wildcards`, `Index`, `Inputs` (nama → artefak), `Outputs` (nama → artefak), `AllInputs`, `Dependencies`, `Dependents` |
| `WorkflowDefinitionException` | siklus (pesannya memuat jalurnya, `a -> b -> a`), dua rule yang menghasilkan satu artefak, target yang tidak dihasilkan apa pun, pola yang membutuhkan wildcard yang tidak didefinisikan, output tanpa badan |

## Menjalankan

| Tipe | Anggota |
|---|---|
| `RunOptions` (record) | `Targets`, `ForceAll`, `ForceRules`, `MaxParallelism` (jumlah prosesor), `MaxThreads` (jumlah prosesor), `MaxMemoryMb` (null = tak terbatas), `Resources` (batas kustom; yang tidak tercantum = tak terbatas), `KeepGoing`, `TimeProvider`, `WriteReport` (true) |
| `JobContext` | `Job`, `JobId`, `Wildcards`, `Attempt`, `Threads`, `LogPath`, `WorkingDirectory`, `CancellationToken`, `ExitCode`; `Input(name)`, `Inputs(name)`, `Output(name)` (path sementara), `FinalOutput(name)`, `Get<T>(ValueRef<T>)`, `GetAll<T>`, `GetAllByKey<T>`, `Set<T>(ValueRef<T>, value)`, `Params<T>()`, `Log(message)`, `RunProcessAsync(fileName, args, throwOnError)`, `RunShellAsync(command, throwOnError)`, `Render(template, quote)` |
| `HashMode` | `Content` (SHA-256, di-cache menurut ukuran dan mtime; bawaan), `Timestamp` (ukuran dan mtime) |
| `ProcessFailedException` | `Command`, `ExitCode`; pesannya memuat baris stderr terakhir |
| `WorkflowStepException` | sebuah langkah tidak menulis output yang dideklarasikan atau tidak mengisi nilai yang dideklarasikan |

Placeholder shell: `{input}`, `{input.NAME}`, `{output}`, `{output.NAME}`, `{wildcards.NAME}`, `{threads}`, `{log}`. `Shell` memakai `cmd.exe /d /s /c` di Windows (path dirender dengan `\`) dan `/bin/sh -c` di sistem lain; path berspasi diberi tanda kutip. Proses berjalan di direktori kerja dengan `BIOTECH_WORKFLOW_THREADS` terisi, dan dimatikan beserta proses anaknya saat pembatalan atau batas waktu.

Aturan up-to-date (kecocokan pertama menjadi alasannya): forced → penanda incomplete → output hilang → tidak ada rekaman → definisi berubah (JSON parameter, versi, perintah, daftar input) → sidik jari input berubah → output disunting di luar workflow → selain itu dilewati. Keputusan diambil saat job siap, sehingga job yang upstream-nya dijalankan ulang dengan output identik akan dilewati.

## Hasil

| Tipe | Anggota |
|---|---|
| `WorkflowRunReport` (record) | `Workflow`, `RunId`, `Status` (`RunStatus`: `Succeeded`, `Failed`, `Cancelled`), `StartedAt`, `FinishedAt`, `Targets`, `Jobs`, `Succeeded`, `Count(status)`, indexer menurut id job; `EnsureSuccess()` (melempar `WorkflowFailedException`), `ToJson()`, `FromJson(json)`, `ToMarkdown()`, `ToHtml()` |
| `JobReport` (record) | `Id`, `Rule`, `Wildcards`, `Status` (`JobStatus`: `NotRun`, `WillRun`, `Skipped`, `Succeeded`, `Failed`, `Cancelled`), `Reason`, `Attempts`, `StartedAt`, `FinishedAt`, `Seconds`, `Outputs`, `Log`, `Error`, `ExitCode` |
| `WorkflowPlan` | `Graph`, `Jobs` (`PlannedJob(Job, WillRun, Reason)`), `ToRun`, indexer; `ToText(includeUpToDate)`, `ToDot()`, `ToMermaid()` |
| Event (`WorkflowEvent(Timestamp)`) | `RunStartedEvent(RunId, JobCount)`, `JobSkippedEvent(JobId, Reason)`, `JobStartedEvent(JobId, Attempt, Reason)`, `JobLogEvent(JobId, Message)`, `JobRetryingEvent(JobId, FailedAttempt, Delay, Error)`, `JobSucceededEvent(JobId, Duration)`, `JobFailedEvent(JobId, Error)`, `JobCancelledEvent(JobId)`, `RunProgressEvent(Finished, Total)`, `RunCompletedEvent(Report)` |

## Direktori state

```text
.biotech-workflow/
  jobs/<hash>.json         run sukses terakhir setiap job: hash definisi, sidik jari input dan output
  incomplete/<hash>        penanda selama job berjalan; tertinggal bila gagal, crash, atau dibatalkan
  values/<hash>.json       nilai ValueRef<T> yang tersimpan
  hashes.json              cache hash konten menurut path, ukuran, dan mtime
  logs/<run id>/<job>.log  log per job
  runs/<run id>/           report.json, report.md, report.html, dag.dot
  last-run.json            laporan terbaru
  lock                     dipegang selama run; run kedua yang bersamaan ditolak
```

Output sementara berada di folder `.biotech-tmp-<hash job>-<percobaan>/` di sebelah setiap output dan dihapus setelah setiap percobaan.
