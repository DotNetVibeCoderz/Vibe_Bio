using System.Globalization;
using System.Net;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;
using Biotech.Net.Core;

namespace Biotech.Net.Workflows;

/// <summary>What happened (or, in a plan, will happen) to a job.</summary>
[JsonConverter(typeof(JsonStringEnumConverter<JobStatus>))]
public enum JobStatus
{
    /// <summary>Not reached: an upstream job failed, or the run stopped first.</summary>
    NotRun,

    /// <summary>Plan only: the job is out of date and will run.</summary>
    WillRun,

    /// <summary>Up to date (or an aggregate): nothing to do.</summary>
    Skipped,

    /// <summary>Ran and succeeded; outputs are in place.</summary>
    Succeeded,

    /// <summary>Failed after all attempts.</summary>
    Failed,

    /// <summary>Cancelled while running.</summary>
    Cancelled,
}

/// <summary>Outcome of a run.</summary>
[JsonConverter(typeof(JsonStringEnumConverter<RunStatus>))]
public enum RunStatus
{
    /// <summary>Every job succeeded or was up to date.</summary>
    Succeeded,

    /// <summary>At least one job failed.</summary>
    Failed,

    /// <summary>The run was cancelled.</summary>
    Cancelled,
}

/// <summary>A run failed; thrown by <see cref="WorkflowRunReport.EnsureSuccess"/>.</summary>
public sealed class WorkflowFailedException : Exception
{
    /// <summary>Creates the exception.</summary>
    public WorkflowFailedException(WorkflowRunReport report)
        : base($"Workflow '{report.Workflow}' {report.Status.ToString().ToLowerInvariant()}: " +
               string.Join("; ", report.Jobs.Where(j => j.Status == JobStatus.Failed).Select(j => $"{j.Id}: {j.Error}")))
        => Report = report;

    /// <summary>The run's report.</summary>
    public WorkflowRunReport Report { get; }
}

/// <summary>One job in a run report.</summary>
public sealed record JobReport
{
    /// <summary>Job id, for example <c>trim[sample=A]</c>.</summary>
    public required string Id { get; init; }

    /// <summary>Rule name.</summary>
    public required string Rule { get; init; }

    /// <summary>Wildcard values.</summary>
    public IReadOnlyDictionary<string, string> Wildcards { get; init; } = new Dictionary<string, string>();

    /// <summary>Outcome.</summary>
    public JobStatus Status { get; init; }

    /// <summary>Why the job ran or was skipped (for example <c>input changed: data/A.fastq</c>).</summary>
    public string? Reason { get; init; }

    /// <summary>Attempts made (0 when skipped).</summary>
    public int Attempts { get; init; }

    /// <summary>Start of the first attempt.</summary>
    public DateTimeOffset? StartedAt { get; init; }

    /// <summary>End of the last attempt.</summary>
    public DateTimeOffset? FinishedAt { get; init; }

    /// <summary>Wall-clock seconds from <see cref="StartedAt"/> to <see cref="FinishedAt"/>.</summary>
    public double? Seconds => StartedAt is { } s && FinishedAt is { } f ? (f - s).TotalSeconds : null;

    /// <summary>Output artifacts.</summary>
    public IReadOnlyList<string> Outputs { get; init; } = [];

    /// <summary>The job's log file, relative to the working directory.</summary>
    public string? Log { get; init; }

    /// <summary>The error of the last failed attempt.</summary>
    public string? Error { get; init; }

    /// <summary>Exit code of a shell or process step's last attempt.</summary>
    public int? ExitCode { get; init; }
}

/// <summary>The record of one run: every job's outcome, timing, reason and log.</summary>
public sealed record WorkflowRunReport
{
    private static readonly JsonSerializerOptions Json = new() { WriteIndented = true, DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull };

    /// <summary>Workflow name.</summary>
    public required string Workflow { get; init; }

    /// <summary>Run id (also the folder under <c>.biotech-workflow/runs</c>).</summary>
    public required string RunId { get; init; }

    /// <summary>Outcome.</summary>
    public RunStatus Status { get; init; }

    /// <summary>Start of the run.</summary>
    public DateTimeOffset StartedAt { get; init; }

    /// <summary>End of the run.</summary>
    public DateTimeOffset FinishedAt { get; init; }

    /// <summary>Requested targets.</summary>
    public IReadOnlyList<string> Targets { get; init; } = [];

    /// <summary>Jobs in topological order.</summary>
    public IReadOnlyList<JobReport> Jobs { get; init; } = [];

    /// <summary>Whether the run succeeded.</summary>
    [JsonIgnore]
    public bool Succeeded => Status == RunStatus.Succeeded;

    /// <summary>Number of jobs with <paramref name="status"/>.</summary>
    public int Count(JobStatus status) => Jobs.Count(j => j.Status == status);

    /// <summary>The report of job <paramref name="id"/>.</summary>
    public JobReport this[string id] => Jobs.FirstOrDefault(j => j.Id == id) ?? throw new KeyNotFoundException($"No job '{id}' in the report.");

    /// <summary>Throws <see cref="WorkflowFailedException"/> unless the run succeeded.</summary>
    public WorkflowRunReport EnsureSuccess() => Succeeded ? this : throw new WorkflowFailedException(this);

    /// <summary>The report as indented JSON.</summary>
    public string ToJson() => JsonSerializer.Serialize(this, Json);

    /// <summary>Reads a report written by <see cref="ToJson"/>.</summary>
    public static WorkflowRunReport FromJson(string json) => JsonSerializer.Deserialize<WorkflowRunReport>(json, Json)
        ?? throw new FormatException("Not a workflow run report.");

    /// <summary>A Markdown summary: totals, then one table row per job.</summary>
    public string ToMarkdown()
    {
        var sb = new StringBuilder();
        sb.Append("# Workflow run: ").Append(Workflow).Append("\n\n");
        sb.Append(Inv($"- Run `{RunId}`: **{Status}** in {(FinishedAt - StartedAt).TotalSeconds:F1} s\n"));
        sb.Append(Inv($"- Jobs: {Jobs.Count} ({Count(JobStatus.Succeeded)} ran, {Count(JobStatus.Skipped)} up to date, {Count(JobStatus.Failed)} failed, {Count(JobStatus.Cancelled)} cancelled, {Count(JobStatus.NotRun)} not run)\n\n"));
        sb.Append("| Job | Status | Attempts | Seconds | Reason / error |\n|---|---|---|---|---|\n");
        foreach (var j in Jobs)
        {
            string note = (j.Error ?? j.Reason ?? "").Replace("|", "\\|").Replace("\n", " ");
            sb.Append(Inv($"| `{j.Id}` | {j.Status} | {j.Attempts} | {(j.Seconds is { } s ? s.ToString("F2", CultureInfo.InvariantCulture) : "")} | {note} |\n"));
        }

        sb.Append("\n_").Append(BiotechInfo.Product).Append(".Workflows · ").Append(BiotechInfo.CreditsId).Append("_\n");
        return sb.ToString();
    }

    /// <summary>A self-contained HTML page with the same content as <see cref="ToMarkdown"/>.</summary>
    public string ToHtml()
    {
        static string E(string? s) => WebUtility.HtmlEncode(s ?? "");
        var sb = new StringBuilder();
        sb.Append("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
        sb.Append("<title>").Append(E(Workflow)).Append(" run</title><style>");
        sb.Append("body{font:14px/1.45 system-ui,sans-serif;margin:24px;color:#1d1b16;background:#fbfaf6}table{border-collapse:collapse;width:100%}");
        sb.Append("th,td{text-align:left;padding:6px 10px;border-bottom:1px solid #e4e0d6;vertical-align:top}code{font:12px ui-monospace,monospace}");
        sb.Append(".Succeeded{color:#2f7a3a}.Skipped{color:#77736a}.Failed{color:#b3261e;font-weight:600}.Cancelled{color:#a26100}.NotRun{color:#77736a}");
        sb.Append("footer{margin-top:24px;color:#77736a;font-size:12px}</style></head><body>\n");
        sb.Append("<h1>").Append(E(Workflow)).Append("</h1>\n");
        sb.Append(Inv($"<p>Run <code>{E(RunId)}</code>: <strong class=\"{Status}\">{Status}</strong> in {(FinishedAt - StartedAt).TotalSeconds:F1} s · {Jobs.Count} jobs: {Count(JobStatus.Succeeded)} ran, {Count(JobStatus.Skipped)} up to date, {Count(JobStatus.Failed)} failed, {Count(JobStatus.Cancelled)} cancelled, {Count(JobStatus.NotRun)} not run.</p>\n"));
        sb.Append("<table><thead><tr><th>Job</th><th>Status</th><th>Attempts</th><th>Seconds</th><th>Reason / error</th><th>Log</th></tr></thead><tbody>\n");
        foreach (var j in Jobs)
        {
            sb.Append(Inv($"<tr><td><code>{E(j.Id)}</code></td><td class=\"{j.Status}\">{j.Status}</td><td>{j.Attempts}</td><td>{(j.Seconds is { } s ? s.ToString("F2", CultureInfo.InvariantCulture) : "")}</td>"));
            sb.Append("<td>").Append(E(j.Error ?? j.Reason)).Append("</td><td>").Append(j.Log is null ? "" : $"<code>{E(j.Log)}</code>").Append("</td></tr>\n");
        }

        sb.Append("</tbody></table>\n<footer>").Append(E(BiotechInfo.Product)).Append(".Workflows · ").Append(E(BiotechInfo.CreditsId)).Append("</footer>\n</body></html>\n");
        return sb.ToString();
    }

    private static string Inv(FormattableString s) => s.ToString(CultureInfo.InvariantCulture);
}

/// <summary>One job of a dry-run plan.</summary>
/// <param name="Job">The job.</param>
/// <param name="WillRun">Whether it is out of date.</param>
/// <param name="Reason">Why it will run, or why it is up to date.</param>
public sealed record PlannedJob(WorkflowJob Job, bool WillRun, string Reason);

/// <summary>A dry run: which jobs would run and why, without running anything.</summary>
public sealed class WorkflowPlan
{
    internal WorkflowPlan(WorkflowGraph graph, IReadOnlyList<PlannedJob> jobs)
    {
        Graph = graph;
        Jobs = jobs;
    }

    /// <summary>The DAG.</summary>
    public WorkflowGraph Graph { get; }

    /// <summary>Every job in topological order.</summary>
    public IReadOnlyList<PlannedJob> Jobs { get; }

    /// <summary>The jobs that would run.</summary>
    public IReadOnlyList<PlannedJob> ToRun => [.. Jobs.Where(j => j.WillRun)];

    /// <summary>The plan of job <paramref name="id"/>.</summary>
    public PlannedJob this[string id] => Jobs.FirstOrDefault(j => j.Job.Id == id) ?? throw new KeyNotFoundException($"No job '{id}' in the plan.");

    /// <summary>A <c>--dry-run</c>-style listing: one line per job to run with its reason, then totals per rule.</summary>
    public string ToText(bool includeUpToDate = false)
    {
        var sb = new StringBuilder();
        foreach (var p in Jobs.Where(j => j.WillRun || includeUpToDate))
            sb.Append(p.WillRun ? "run   " : "skip  ").Append(p.Job.Id).Append("  # ").Append(p.Reason).Append('\n');
        sb.Append("\nJob counts:\n");
        foreach (var g in Jobs.GroupBy(j => j.Job.Rule.Name))
            sb.Append(CultureInfo.InvariantCulture, $"  {g.Count(j => j.WillRun),5} to run, {g.Count(j => !j.WillRun),5} up to date  {g.Key}\n");
        sb.Append(CultureInfo.InvariantCulture, $"  {ToRun.Count,5} to run, {Jobs.Count - ToRun.Count,5} up to date  total\n");
        return sb.ToString();
    }

    /// <summary>The DAG as Graphviz DOT, coloured by plan status.</summary>
    public string ToDot() => Graph.ToDot(status: Status());

    /// <summary>The DAG as a Mermaid flowchart, coloured by plan status.</summary>
    public string ToMermaid() => Graph.ToMermaid(status: Status());

    private Dictionary<string, JobStatus> Status() => Jobs.ToDictionary(j => j.Job.Id, j => j.WillRun ? JobStatus.WillRun : JobStatus.Skipped);
}
