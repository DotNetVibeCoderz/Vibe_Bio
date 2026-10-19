using System.Globalization;
using System.Text;
using System.Text.Json;

namespace Biotech.Net.Workflows;

/// <summary>Settings for one run (or dry run) of a workflow.</summary>
public sealed record RunOptions
{
    /// <summary>Files (<c>results/x.vcf</c>), values (<c>value:key</c>) or rule names to build; empty uses <see cref="Workflow.Targets"/>.</summary>
    public IReadOnlyList<string> Targets { get; init; } = [];

    /// <summary>Rerun every job, up to date or not.</summary>
    public bool ForceAll { get; init; }

    /// <summary>Rerun every job of these rules (downstream jobs rerun only if the outputs change).</summary>
    public IReadOnlyCollection<string> ForceRules { get; init; } = [];

    /// <summary>Jobs running at once. Defaults to the processor count.</summary>
    public int MaxParallelism { get; init; } = Environment.ProcessorCount;

    /// <summary>Total threads shared by running jobs (<see cref="Rule.Threads"/>). Defaults to the processor count.</summary>
    public int MaxThreads { get; init; } = Environment.ProcessorCount;

    /// <summary>Total memory in MB shared by running jobs (<see cref="Rule.Memory"/>); null does not limit memory.</summary>
    public int? MaxMemoryMb { get; init; }

    /// <summary>Limits of custom resources (<see cref="Rule.Resource"/>); a resource not listed is unlimited.</summary>
    public IReadOnlyDictionary<string, int> Resources { get; init; } = new Dictionary<string, int>();

    /// <summary>Keep running independent jobs after a failure (Snakemake's <c>--keep-going</c>). Off by default: no new jobs start after a failure.</summary>
    public bool KeepGoing { get; init; }

    /// <summary>The clock for timestamps, retry backoff and timeouts; tests pass a fake one.</summary>
    public TimeProvider TimeProvider { get; init; } = TimeProvider.System;

    /// <summary>Write <c>report.json</c>, <c>report.md</c> and <c>report.html</c> under <c>.biotech-workflow/runs/&lt;run id&gt;/</c>. On by default.</summary>
    public bool WriteReport { get; init; } = true;
}

internal sealed record Decision(bool Run, string Reason, string? Error = null);

internal sealed record Grant(int Threads, int Memory, Dictionary<string, int> Resources);

internal sealed record Outcome(WorkflowJob Job, JobReport Report);

/// <summary>Runs a <see cref="WorkflowGraph"/>: schedules ready jobs within the limits, decides up-to-date-ness from hashes, retries, commits outputs atomically.</summary>
internal sealed class Executor
{
    private static readonly JsonSerializerOptions ParamsJson = new() { WriteIndented = false };
    private readonly Action<WorkflowEvent> _emit;
    private readonly Dictionary<string, StreamWriter> _logs = new(StringComparer.Ordinal);
    private readonly Lock _logGate = new();
    private int _usedThreads, _usedMemory;
    private readonly Dictionary<string, int> _usedResources = new(StringComparer.Ordinal);

    public Executor(Workflow workflow, WorkflowGraph graph, RunOptions options, Action<WorkflowEvent> emit)
    {
        Workflow = workflow;
        Graph = graph;
        Options = options;
        Store = new StateStore(workflow);
        _emit = emit;
        var now = options.TimeProvider.GetUtcNow();
        RunId = now.ToString("yyyyMMdd-HHmmss", CultureInfo.InvariantCulture) + "-" + Guid.NewGuid().ToString("N")[..6];
    }

    public Workflow Workflow { get; }

    public WorkflowGraph Graph { get; }

    public RunOptions Options { get; }

    public StateStore Store { get; }

    public string RunId { get; }

    private DateTimeOffset Now => Options.TimeProvider.GetUtcNow();

    // ------------------------------------------------------------------ up-to-date decisions

    public string Definition(WorkflowJob job)
    {
        var r = job.Rule;
        object? p = r.ParamsFactory?.Invoke(job.Wildcards);
        var def = new
        {
            rule = r.Name,
            version = r.CodeVersion,
            kind = r.Kind.ToString(),
            command = r.CommandTemplate,
            parameters = p is null ? null : JsonSerializer.Serialize(p, p.GetType(), ParamsJson),
            outputs = job.Outputs.Select(o => o.Key + "=" + o.Value).ToArray(),
            inputs = job.Inputs.Select(i => i.Key + "=" + string.Join(";", i.Value)).ToArray(),
        };
        return "sha256:" + StateStore.Sha(JsonSerializer.Serialize(def, ParamsJson));
    }

    public Decision Decide(WorkflowJob job, Func<WorkflowJob, bool>? upstreamWillRun)
    {
        var r = job.Rule;
        if (r.Body is null) return new(false, "aggregate rule (no body)");
        bool planning = upstreamWillRun is not null;

        if (!planning)
        {
            foreach (var input in job.AllInputs)
            {
                if (Store.Fingerprint(input, Workflow.HashMode) is null)
                {
                    string why = Graph.SourceArtifacts.Contains(input) ? "no rule produces it and it does not exist" : "its producer did not create it";
                    return new(true, "missing input", $"missing input {input} ({why})");
                }
            }
        }

        if (Options.ForceAll) return new(true, "forced (all jobs)");
        if (Options.ForceRules.Contains(r.Name)) return new(true, $"forced (rule {r.Name})");
        if (Store.IsIncomplete(job)) return new(true, "incomplete: a previous run started this job and did not finish it");
        if (job.Outputs.Count == 0) return new(true, "rule has no outputs, so it always runs");
        foreach (var o in job.Outputs.Values)
        {
            if (Store.Fingerprint(o, Workflow.HashMode) is null) return new(true, $"missing output {o}");
        }

        var record = Store.Load(job);
        if (record is null) return new(true, "no record of a previous run");
        if (planning && job.Dependencies.FirstOrDefault(upstreamWillRun!) is { } up) return new(true, $"upstream {up.Id} will run");
        if (record.Definition != Definition(job)) return new(true, "rule definition changed (parameters, version, command or inputs list)");

        foreach (var input in job.AllInputs.Distinct())
        {
            string? fp = Store.Fingerprint(input, Workflow.HashMode);
            if (fp is null) return new(true, $"missing input {input}");
            if (!record.Inputs.TryGetValue(input.ToString(), out var old) || old != fp) return new(true, $"input changed: {input}");
        }

        foreach (var (_, o) in job.Outputs)
        {
            if (!record.Outputs.TryGetValue(o.ToString(), out var old) || old != Store.Fingerprint(o, Workflow.HashMode))
                return new(true, $"output modified outside the workflow: {o}");
        }

        return new(false, "up to date");
    }

    public WorkflowPlan Plan()
    {
        var willRun = new HashSet<WorkflowJob>();
        var planned = new List<PlannedJob>();
        foreach (var job in Graph.Jobs)
        {
            var d = Decide(job, j => willRun.Contains(j));
            if (d.Run) willRun.Add(job);
            planned.Add(new PlannedJob(job, d.Run, d.Reason));
        }

        return new WorkflowPlan(Graph, planned);
    }

    // ------------------------------------------------------------------ scheduling

    public async Task<WorkflowRunReport> RunAsync(CancellationToken cancellationToken)
    {
        var started = Now;
        _emit(new RunStartedEvent(started, RunId, Graph.Jobs.Count));
        var reports = new Dictionary<WorkflowJob, JobReport>();
        var remaining = Graph.Jobs.ToDictionary(j => j, j => j.Dependencies.Count);
        var ready = new SortedSet<int>(Graph.Jobs.Where(j => j.Dependencies.Count == 0).Select(j => j.Index));
        var running = new Dictionary<Task<Outcome>, Grant>();
        bool stop = false;
        int finished = 0;

        while (true)
        {
            if (!stop && !cancellationToken.IsCancellationRequested)
            {
                foreach (int index in ready.ToList())
                {
                    if (running.Count >= Math.Max(1, Options.MaxParallelism)) break;
                    var job = Graph.Jobs[index];
                    if (TryAcquire(job) is not { } grant) continue; // backfill: a smaller job may still fit
                    ready.Remove(index);
                    running[Task.Run(() => ExecuteAsync(job, grant, cancellationToken), CancellationToken.None)] = grant;
                }
            }

            if (running.Count == 0) break;
            var done = await Task.WhenAny(running.Keys).ConfigureAwait(false);
            Release(running[done]);
            running.Remove(done);
            var outcome = await done.ConfigureAwait(false);
            reports[outcome.Job] = outcome.Report;
            _emit(new RunProgressEvent(Now, ++finished, Graph.Jobs.Count));

            switch (outcome.Report.Status)
            {
                case JobStatus.Succeeded or JobStatus.Skipped:
                    foreach (var d in outcome.Job.Dependents)
                    {
                        if (--remaining[d] == 0) ready.Add(d.Index);
                    }

                    break;
                case JobStatus.Failed:
                    if (!Options.KeepGoing) stop = true;
                    break;
                default:
                    stop = true;
                    break;
            }
        }

        Store.SaveHashes();
        CloseLogs();
        var jobs = Graph.Jobs.Select(j => reports.TryGetValue(j, out var r) ? r : new JobReport
        {
            Id = j.Id,
            Rule = j.Rule.Name,
            Wildcards = j.Wildcards.ToDictionary(),
            Status = JobStatus.NotRun,
            Reason = cancellationToken.IsCancellationRequested ? "run cancelled" : "an upstream job failed",
            Outputs = [.. j.Outputs.Values.Select(o => o.ToString())],
        }).ToList();

        var status = cancellationToken.IsCancellationRequested || jobs.Any(j => j.Status == JobStatus.Cancelled) ? RunStatus.Cancelled
            : jobs.Any(j => j.Status is JobStatus.Failed or JobStatus.NotRun) ? RunStatus.Failed
            : RunStatus.Succeeded;
        var report = new WorkflowRunReport
        {
            Workflow = Workflow.Name,
            RunId = RunId,
            Status = status,
            StartedAt = started,
            FinishedAt = Now,
            Targets = [.. Graph.Targets.Select(t => t.ToString())],
            Jobs = jobs,
        };
        if (Options.WriteReport) SaveReport(report);
        _emit(new RunCompletedEvent(Now, report));
        return report;
    }

    private void SaveReport(WorkflowRunReport report)
    {
        string dir = Store.RunDirectory(RunId);
        File.WriteAllText(Path.Combine(dir, "report.json"), report.ToJson());
        File.WriteAllText(Path.Combine(dir, "report.md"), report.ToMarkdown());
        File.WriteAllText(Path.Combine(dir, "report.html"), report.ToHtml());
        File.WriteAllText(Path.Combine(dir, "dag.dot"), Graph.ToDot(status: report.Jobs.ToDictionary(j => j.Id, j => j.Status)));
        StateStore.WriteAtomic(Path.Combine(Store.Root, "last-run.json"), Encoding.UTF8.GetBytes(report.ToJson()));
    }

    private Grant? TryAcquire(WorkflowJob job)
    {
        var r = job.Rule;
        int threads = Math.Min(r.ThreadCount, Math.Max(1, Options.MaxThreads));
        int memory = Options.MaxMemoryMb is int maxMem ? Math.Min(r.MemoryMb, maxMem) : r.MemoryMb;
        if (_usedThreads + threads > Math.Max(1, Options.MaxThreads)) return null;
        if (Options.MaxMemoryMb is int limit && _usedMemory + memory > limit) return null;
        var custom = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (var (name, amount) in r.Resources)
        {
            if (!Options.Resources.TryGetValue(name, out int cap)) continue;
            int take = Math.Min(amount, cap);
            if (_usedResources.GetValueOrDefault(name) + take > cap) return null;
            custom[name] = take;
        }

        _usedThreads += threads;
        _usedMemory += memory;
        foreach (var (name, take) in custom) _usedResources[name] = _usedResources.GetValueOrDefault(name) + take;
        return new Grant(threads, memory, custom);
    }

    private void Release(Grant g)
    {
        _usedThreads -= g.Threads;
        _usedMemory -= g.Memory;
        foreach (var (name, take) in g.Resources) _usedResources[name] -= take;
    }

    // ------------------------------------------------------------------ one job

    private JobReport Report(WorkflowJob job, JobStatus status, string? reason, int attempts = 0, DateTimeOffset? start = null,
        string? log = null, string? error = null, int? exit = null) => new()
    {
        Id = job.Id,
        Rule = job.Rule.Name,
        Wildcards = job.Wildcards.ToDictionary(),
        Status = status,
        Reason = reason,
        Attempts = attempts,
        StartedAt = start,
        FinishedAt = start is null ? null : Now,
        Outputs = [.. job.Outputs.Values.Select(o => o.ToString())],
        Log = log is null ? null : Workflow.Relative(log),
        Error = error,
        ExitCode = exit,
    };

    private async Task<Outcome> ExecuteAsync(WorkflowJob job, Grant grant, CancellationToken runToken)
    {
        Decision decision;
        try
        {
            decision = Decide(job, null);
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException or JsonException)
        {
            decision = new(true, "state unreadable", ex.Message);
        }

        if (decision.Error is not null)
        {
            _emit(new JobFailedEvent(Now, job.Id, decision.Error));
            return new(job, Report(job, JobStatus.Failed, decision.Reason, error: decision.Error));
        }

        if (!decision.Run)
        {
            _emit(new JobSkippedEvent(Now, job.Id, decision.Reason));
            return new(job, Report(job, JobStatus.Skipped, decision.Reason));
        }

        var rule = job.Rule;
        var start = Now;
        string logPath = Store.LogPath(RunId, job);
        Store.MarkIncomplete(job);
        object? parameters = rule.ParamsFactory?.Invoke(job.Wildcards);
        string? error = null;
        int? exit = null;
        int attempts = rule.RetryCount + 1;
        for (int attempt = 1; attempt <= attempts; attempt++)
        {
            var temps = TempPaths(job, attempt);
            using var timeout = rule.TimeLimit is { } limit ? new CancellationTokenSource(limit, Options.TimeProvider) : new CancellationTokenSource();
            using var linked = CancellationTokenSource.CreateLinkedTokenSource(runToken, timeout.Token);
            var ctx = new JobContext(this, job, attempt, grant.Threads, logPath, temps, parameters, linked.Token);
            _emit(new JobStartedEvent(Now, job.Id, attempt, decision.Reason));
            Log(ctx, $"attempt {attempt}/{attempts} ({decision.Reason}); threads {grant.Threads}");
            try
            {
                await rule.Body!(ctx).ConfigureAwait(false);
                linked.Token.ThrowIfCancellationRequested();
                Commit(job, ctx, temps);
                Log(ctx, "succeeded");
                _emit(new JobSucceededEvent(Now, job.Id, Now - start));
                return new(job, Report(job, JobStatus.Succeeded, decision.Reason, attempt, start, logPath, exit: ctx.ExitCode));
            }
            catch (OperationCanceledException) when (runToken.IsCancellationRequested)
            {
                Log(ctx, "cancelled");
                _emit(new JobCancelledEvent(Now, job.Id));
                return new(job, Report(job, JobStatus.Cancelled, decision.Reason, attempt, start, logPath, "cancelled", ctx.ExitCode));
            }
            catch (OperationCanceledException) when (timeout.IsCancellationRequested)
            {
                error = $"timed out after {rule.TimeLimit!.Value.TotalSeconds.ToString("0.###", CultureInfo.InvariantCulture)} s";
            }
            catch (Exception ex)
            {
                error = ex is ProcessFailedException or WorkflowStepException ? ex.Message : $"{ex.GetType().Name}: {ex.Message}";
            }
            finally
            {
                DeleteTemps(temps);
            }

            exit = ctx.ExitCode;
            Log(ctx, "failed: " + error);
            if (attempt < attempts)
            {
                var delay = TimeSpan.FromTicks((long)(rule.RetryDelay.Ticks * Math.Pow(rule.RetryBackoffFactor, attempt - 1)));
                _emit(new JobRetryingEvent(Now, job.Id, attempt, delay, error));
                try
                {
                    await Task.Delay(delay, Options.TimeProvider, runToken).ConfigureAwait(false);
                }
                catch (OperationCanceledException)
                {
                    _emit(new JobCancelledEvent(Now, job.Id));
                    return new(job, Report(job, JobStatus.Cancelled, decision.Reason, attempt, start, logPath, "cancelled", exit));
                }
            }
        }

        _emit(new JobFailedEvent(Now, job.Id, error!));
        return new(job, Report(job, JobStatus.Failed, decision.Reason, attempts, start, logPath, error, exit));
    }

    private Dictionary<string, string> TempPaths(WorkflowJob job, int attempt)
    {
        var temps = new Dictionary<string, string>(StringComparer.Ordinal);
        string tag = StateStore.Sha(job.Id)[..10] + "-" + attempt.ToString(CultureInfo.InvariantCulture);
        foreach (var (name, key) in job.Outputs)
        {
            if (key.Kind != ArtifactKind.File) continue;
            string final = Workflow.FullPath(key);
            string dir = Path.Combine(Path.GetDirectoryName(final)!, ".biotech-tmp-" + tag);
            if (Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            Directory.CreateDirectory(dir);
            temps[name] = Path.Combine(dir, Path.GetFileName(final));
        }

        return temps;
    }

    private static void DeleteTemps(Dictionary<string, string> temps)
    {
        foreach (var dir in temps.Values.Select(Path.GetDirectoryName).Distinct())
        {
            try
            {
                if (dir is not null && Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            }
            catch (IOException)
            {
                // A tool may still hold a handle; the next attempt or run removes it.
            }
        }
    }

    private void Commit(WorkflowJob job, JobContext ctx, Dictionary<string, string> temps)
    {
        foreach (var (name, key) in job.Outputs)
        {
            if (key.Kind == ArtifactKind.File && !File.Exists(temps[name]))
                throw new WorkflowStepException($"{job.Id} did not write output '{name}' ({key}); write it to ctx.Output(\"{name}\").");
            if (key.Kind == ArtifactKind.Value && !ctx.PendingValues.ContainsKey(name))
                throw new WorkflowStepException($"{job.Id} did not set value output '{name}' ({key}); call ctx.Set.");
        }

        // Everything exists: move files into place (same directory, so each move is a rename), then store values.
        foreach (var (name, key) in job.Outputs)
        {
            if (key.Kind == ArtifactKind.File) File.Move(temps[name], Workflow.FullPath(key), overwrite: true);
            else Store.WriteValue(key.Id, ctx.PendingValues[name]);
        }

        var record = new JobRecord
        {
            Job = job.Id,
            Definition = Definition(job),
            Signature = "",
            Inputs = job.AllInputs.Distinct().ToDictionary(i => i.ToString(), i => Store.Fingerprint(i, Workflow.HashMode) ?? ""),
            Outputs = job.Outputs.Values.ToDictionary(o => o.ToString(), o => Store.Fingerprint(o, Workflow.HashMode) ?? ""),
            FinishedAt = Now,
        };
        record = record with { Signature = "sha256:" + StateStore.Sha(record.Definition + string.Concat(record.Inputs.Select(kv => kv.Key + kv.Value))) };
        if (job.Outputs.Count > 0) Store.Save(job, record);
        Store.ClearIncomplete(job);
    }

    // ------------------------------------------------------------------ logs

    public void Log(JobContext ctx, string message)
    {
        var now = Now;
        lock (_logGate)
        {
            if (!_logs.TryGetValue(ctx.LogPath, out var w))
            {
                w = new StreamWriter(ctx.LogPath, append: true, new UTF8Encoding(false)) { AutoFlush = true };
                _logs[ctx.LogPath] = w;
            }

            w.WriteLine($"{now.ToString("O", CultureInfo.InvariantCulture)} {message}");
        }

        _emit(new JobLogEvent(now, ctx.JobId, message));
    }

    private void CloseLogs()
    {
        lock (_logGate)
        {
            foreach (var w in _logs.Values) w.Dispose();
            _logs.Clear();
        }
    }
}

/// <summary>A step broke the workflow's contract, for example by not writing a declared output.</summary>
public sealed class WorkflowStepException : Exception
{
    /// <summary>Creates the exception.</summary>
    public WorkflowStepException(string message) : base(message) { }
}
