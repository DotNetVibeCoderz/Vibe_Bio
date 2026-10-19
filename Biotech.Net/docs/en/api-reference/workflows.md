# API reference: Workflows

Package `Gravicode.Biotech.Net.Workflows` · namespace `Biotech.Net.Workflows`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Depends only on Core. A Snakemake/Nextflow-style engine: rules with declared file and value inputs/outputs, wildcards, a DAG built backwards from targets, parallel scheduling with resources, retries, timeouts, cancellation, incremental re-runs from content hashes, atomic outputs, events and reports. Tutorial: [workflows.md](../tutorials/workflows.md).

## Defining a workflow

| Type | Members |
|---|---|
| `Workflow` | `new Workflow(workingDirectory = null, name = "workflow")`; `Name`, `WorkingDirectory`, `StateDirectory` (`.biotech-workflow`), `HashMode`, `Rules`, `Targets`; `Rule(name)`, `Value<T>(name, pattern)`, `Target(params string[])`, `Target(IEnumerable<string>)`; `Build(targets)` → `WorkflowGraph`, `Plan(options)` → `WorkflowPlan`, `RunAsync(options, progress, cancellationToken)` → `WorkflowRunReport`, `RunStreamingAsync(options, cancellationToken)` → `IAsyncEnumerable<WorkflowEvent>`, `LastReport()`, `ToDot(rules)`, `ToMermaid(rules)`; `StateDirectoryName` |
| `Rule` (fluent; every method returns the rule) | `Input(name, params string[] patterns)`, `Input(name, IEnumerable<string>)`, `Input(name, Func<Wildcards, string>)`, `Input(name, Func<Wildcards, IEnumerable<string>>)`, `Input<T>(ValueRef<T>)`, `Input<T>(ValueRef<T>, rows)` (gather); `Output(name, pattern)`, `Output<T>(ValueRef<T>)`; `Params(object)`, `Params(Func<Wildcards, object>)`, `Version(string)`, `Describe(string)`; `Threads(n)`, `Memory(mb)`, `Resource(name, amount)`; `Retries(n, delay, backoffFactor)`, `Timeout(limit)`; bodies `Run(Func<JobContext, Task>)`, `Run(Action<JobContext>)`, `Shell(command)`, `Process(fileName, params string[] args)`; properties `Name`, `Kind` (`RuleKind`: `Aggregate`, `Code`, `Shell`, `Process`), `CodeVersion`, `ThreadCount`, `MemoryMb`, `Resources`, `RetryCount`, `RetryDelay`, `RetryBackoffFactor`, `TimeLimit`, `WildcardNames`, `Description` |
| `ValueRef<T>` | `Name`, `Pattern`, `Key(wildcards)`; values are JSON-serialized into the state directory |
| `WildcardPattern` | `new WildcardPattern("work/{sample}.{n,\\d+}.fq")`; `Text`, `WildcardNames`, `IsConcrete`, `Match(text)` → `Wildcards?`, `Format(values)`; `{name}` matches `[^/]+`, `{name,regex}` a custom expression, `{{`/`}}` literal braces |
| `Wildcards` | read-only dictionary of name → value in pattern order; `Of(("sample", "A"))`, `Empty`, indexer with a helpful `KeyNotFoundException`, `ToString()` = `sample=A,chunk=1` |
| `SampleSheet` | `Parse(text)` (TSV or CSV with header, `#` comments), `Load(path)`, `FromColumn(column, values)`, `new SampleSheet(rows)`; `Columns`, `Column(name)`, `Find(column, value)`, indexer, enumerable of rows |
| `Expand` | `Product(pattern, params (name, values)[])` (cross product), `Each(pattern, name, values)`, `Rows(pattern, rows)` |
| `ArtifactKey` | `(Kind, Id)`; `ArtifactKind.File` (path relative to the working directory, `/` separators) or `ArtifactKind.Value` (`value:key` in targets) |

A target is a file path, `value:key`, or the name of a rule without wildcards. With no targets, `Build` uses `Target(...)`, or failing that every rule without wildcards. A rule without a body is an aggregate (like Snakemake's `rule all`) and may not declare outputs.

## The DAG

| Type | Members |
|---|---|
| `WorkflowGraph` | `Jobs` (topological order), `Targets`, `SourceArtifacts` (inputs no rule produces), indexer by job id; `ToDot(rules, status)`, `ToMermaid(rules, status)` |
| `WorkflowJob` | `Id` (`rule` or `rule[sample=A]`), `Rule`, `Wildcards`, `Index`, `Inputs` (name → artifacts), `Outputs` (name → artifact), `AllInputs`, `Dependencies`, `Dependents` |
| `WorkflowDefinitionException` | cycles (message lists the path, `a -> b -> a`), two rules producing one artifact, a target nothing produces, a pattern needing an undefined wildcard, outputs without a body |

## Running

| Type | Members |
|---|---|
| `RunOptions` (record) | `Targets`, `ForceAll`, `ForceRules`, `MaxParallelism` (processor count), `MaxThreads` (processor count), `MaxMemoryMb` (null = unlimited), `Resources` (custom limits; unlisted = unlimited), `KeepGoing`, `TimeProvider`, `WriteReport` (true) |
| `JobContext` | `Job`, `JobId`, `Wildcards`, `Attempt`, `Threads`, `LogPath`, `WorkingDirectory`, `CancellationToken`, `ExitCode`; `Input(name)`, `Inputs(name)`, `Output(name)` (temporary path), `FinalOutput(name)`, `Get<T>(ValueRef<T>)`, `GetAll<T>`, `GetAllByKey<T>`, `Set<T>(ValueRef<T>, value)`, `Params<T>()`, `Log(message)`, `RunProcessAsync(fileName, args, throwOnError)`, `RunShellAsync(command, throwOnError)`, `Render(template, quote)` |
| `HashMode` | `Content` (SHA-256, cached by size and mtime; default), `Timestamp` (size and mtime) |
| `ProcessFailedException` | `Command`, `ExitCode`; message includes the last stderr lines |
| `WorkflowStepException` | a step did not write a declared output or set a declared value |

Shell placeholders: `{input}`, `{input.NAME}`, `{output}`, `{output.NAME}`, `{wildcards.NAME}`, `{threads}`, `{log}`. `Shell` uses `cmd.exe /d /s /c` on Windows (paths rendered with `\`) and `/bin/sh -c` elsewhere; paths with spaces are quoted. Processes run in the working directory with `BIOTECH_WORKFLOW_THREADS` set, and are killed with their child processes on cancellation or timeout.

Up-to-date rule (first match is the reason): forced → incomplete marker → missing output → no record → definition changed (parameters JSON, version, command, input list) → input fingerprint changed → output edited outside the workflow → otherwise skipped. The decision is made when the job becomes ready, so a job whose upstream reran with identical output is skipped.

## Results

| Type | Members |
|---|---|
| `WorkflowRunReport` (record) | `Workflow`, `RunId`, `Status` (`RunStatus`: `Succeeded`, `Failed`, `Cancelled`), `StartedAt`, `FinishedAt`, `Targets`, `Jobs`, `Succeeded`, `Count(status)`, indexer by job id; `EnsureSuccess()` (throws `WorkflowFailedException`), `ToJson()`, `FromJson(json)`, `ToMarkdown()`, `ToHtml()` |
| `JobReport` (record) | `Id`, `Rule`, `Wildcards`, `Status` (`JobStatus`: `NotRun`, `WillRun`, `Skipped`, `Succeeded`, `Failed`, `Cancelled`), `Reason`, `Attempts`, `StartedAt`, `FinishedAt`, `Seconds`, `Outputs`, `Log`, `Error`, `ExitCode` |
| `WorkflowPlan` | `Graph`, `Jobs` (`PlannedJob(Job, WillRun, Reason)`), `ToRun`, indexer; `ToText(includeUpToDate)`, `ToDot()`, `ToMermaid()` |
| Events (`WorkflowEvent(Timestamp)`) | `RunStartedEvent(RunId, JobCount)`, `JobSkippedEvent(JobId, Reason)`, `JobStartedEvent(JobId, Attempt, Reason)`, `JobLogEvent(JobId, Message)`, `JobRetryingEvent(JobId, FailedAttempt, Delay, Error)`, `JobSucceededEvent(JobId, Duration)`, `JobFailedEvent(JobId, Error)`, `JobCancelledEvent(JobId)`, `RunProgressEvent(Finished, Total)`, `RunCompletedEvent(Report)` |

## State directory

```text
.biotech-workflow/
  jobs/<hash>.json         last successful run of each job: definition hash, input and output fingerprints
  incomplete/<hash>        marker while a job runs; left behind by a failure, crash or cancellation
  values/<hash>.json       stored ValueRef<T> values
  hashes.json              content-hash cache keyed by path, size and mtime
  logs/<run id>/<job>.log  per-job logs
  runs/<run id>/           report.json, report.md, report.html, dag.dot
  last-run.json            the latest report
  lock                     held during a run; a second concurrent run is refused
```

Temporary outputs live in `.biotech-tmp-<job hash>-<attempt>/` folders next to each output and are removed after every attempt.
