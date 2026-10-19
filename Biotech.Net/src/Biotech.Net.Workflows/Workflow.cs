using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Channels;

namespace Biotech.Net.Workflows;

/// <summary>
/// A workflow: rules whose declared inputs and outputs (files and in-memory values) connect them into a DAG,
/// in the manner of Snakemake and Nextflow but in plain C#.
/// </summary>
/// <remarks>
/// <para>
/// Jobs are found backwards from the targets: each requested file is matched against the rules' output patterns,
/// which binds the wildcards (<c>{sample}</c>) and gives the job's inputs, and so on down to files that exist.
/// A rule with wildcards is the <i>map</i> (scatter) side; a rule whose input is an <see cref="Expand"/> list is the
/// <i>collect</i> (gather) side.
/// </para>
/// <para>
/// A run executes ready jobs in parallel within <see cref="RunOptions"/> limits. Before a job runs, its inputs'
/// content hashes and its definition (parameters, code version, command) are compared with what
/// <c>.biotech-workflow/</c> recorded for its last success; an unchanged job is skipped, so a re-run after an edit
/// redoes only what the edit affects. Outputs are written to temporary paths and moved into place only when the job
/// succeeds.
/// </para>
/// </remarks>
public sealed class Workflow
{
    /// <summary>The state directory's name inside the working directory.</summary>
    public const string StateDirectoryName = ".biotech-workflow";

    private readonly List<Rule> _rules = [];
    private readonly List<string> _targets = [];
    private readonly Dictionary<string, object> _values = new(StringComparer.Ordinal);

    /// <summary>Creates a workflow rooted at <paramref name="workingDirectory"/> (the current directory by default).</summary>
    public Workflow(string? workingDirectory = null, string name = "workflow")
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        Name = name;
        WorkingDirectory = Path.TrimEndingDirectorySeparator(Path.GetFullPath(workingDirectory ?? Directory.GetCurrentDirectory()));
        StateDirectory = Path.Combine(WorkingDirectory, StateDirectoryName);
    }

    /// <summary>The workflow's name, for reports and graphs.</summary>
    public string Name { get; }

    /// <summary>Relative paths are relative to this directory; shell and process steps run in it.</summary>
    public string WorkingDirectory { get; }

    /// <summary><c>&lt;working directory&gt;/.biotech-workflow</c>: job records, values, hash cache, logs and reports.</summary>
    public string StateDirectory { get; }

    /// <summary>How files are fingerprinted. Defaults to <see cref="HashMode.Content"/>.</summary>
    public HashMode HashMode { get; set; } = HashMode.Content;

    /// <summary>The rules, in declaration order.</summary>
    public IReadOnlyList<Rule> Rules => _rules;

    /// <summary>Default targets set with <see cref="Target(string[])"/>.</summary>
    public IReadOnlyList<string> Targets => _targets;

    /// <summary>Adds a rule; configure it with the returned builder.</summary>
    /// <exception cref="ArgumentException">A rule with that name exists.</exception>
    public Rule Rule(string name)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        if (_rules.Exists(r => r.Name == name)) throw new ArgumentException($"The workflow already has a rule '{name}'.", nameof(name));
        var rule = new Rule(this, name);
        _rules.Add(rule);
        return rule;
    }

    /// <summary>A typed in-memory value handle; <paramref name="pattern"/> (default: the name) may contain wildcards, such as <c>qc/{sample}</c>.</summary>
    public ValueRef<T> Value<T>(string name, string? pattern = null)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        if (_values.TryGetValue(name, out var existing))
            return existing as ValueRef<T> ?? throw new ArgumentException($"Value '{name}' is already declared with another type.", nameof(name));
        var v = new ValueRef<T>(name, pattern ?? name);
        _values[name] = v;
        return v;
    }

    /// <summary>Adds default targets: file paths, <c>value:key</c>, or names of rules without wildcards.</summary>
    public Workflow Target(params string[] targets)
    {
        _targets.AddRange(targets);
        return this;
    }

    /// <summary>Adds default targets from a list (for example an <see cref="Expand"/> result).</summary>
    public Workflow Target(IEnumerable<string> targets) => Target([.. targets]);

    /// <summary>Builds the DAG for <paramref name="targets"/> (or the default targets).</summary>
    /// <exception cref="WorkflowDefinitionException">A cycle, an ambiguous producer, a missing input, or a rule that cannot be expanded.</exception>
    public WorkflowGraph Build(params string[] targets)
    {
        foreach (var r in _rules)
        {
            if (r.Body is null && r.OutputSpecs.Count > 0)
                throw new WorkflowDefinitionException($"Rule '{r.Name}' declares outputs but has no body (Run, Shell or Process).");
        }

        IReadOnlyList<string> list = targets.Length > 0 ? targets
            : _targets.Count > 0 ? _targets
            : [.. _rules.Where(r => r.WildcardNames.Count == 0).Select(r => r.Name)];
        if (list.Count == 0)
            throw new WorkflowDefinitionException("No targets: call Target(...), pass RunOptions.Targets, or add a rule without wildcards.");
        return new GraphBuilder(this).Build(list);
    }

    /// <summary>A dry run: which jobs would run, and why, without running anything.</summary>
    public WorkflowPlan Plan(RunOptions? options = null)
    {
        options ??= new RunOptions();
        return new Executor(this, Build([.. options.Targets]), options, _ => { }).Plan();
    }

    /// <summary>
    /// Runs the workflow. Failures and cancellation do not throw: check <see cref="WorkflowRunReport.Status"/>
    /// (or call <see cref="WorkflowRunReport.EnsureSuccess"/>).
    /// </summary>
    /// <exception cref="WorkflowDefinitionException">The DAG cannot be built.</exception>
    /// <exception cref="InvalidOperationException">Another run holds the state directory's lock.</exception>
    public async Task<WorkflowRunReport> RunAsync(RunOptions? options = null, IProgress<WorkflowEvent>? progress = null,
        CancellationToken cancellationToken = default)
    {
        options ??= new RunOptions();
        var graph = Build([.. options.Targets]);
        var gate = new Lock();
        using var runLock = AcquireLock();
        var executor = new Executor(this, graph, options, e =>
        {
            lock (gate) progress?.Report(e);
        });
        return await executor.RunAsync(cancellationToken).ConfigureAwait(false);
    }

    /// <summary>Runs the workflow and streams its events; the last one is a <see cref="RunCompletedEvent"/> with the report.</summary>
    public async IAsyncEnumerable<WorkflowEvent> RunStreamingAsync(RunOptions? options = null,
        [EnumeratorCancellation] CancellationToken cancellationToken = default)
    {
        var channel = Channel.CreateUnbounded<WorkflowEvent>(new UnboundedChannelOptions { SingleReader = true });
        var run = Task.Run(async () =>
        {
            try
            {
                await RunAsync(options, new ChannelProgress(channel.Writer), cancellationToken).ConfigureAwait(false);
                channel.Writer.TryComplete();
            }
            catch (Exception ex)
            {
                channel.Writer.TryComplete(ex);
            }
        }, CancellationToken.None);

        await foreach (var e in channel.Reader.ReadAllAsync(CancellationToken.None).ConfigureAwait(false)) yield return e;
        await run.ConfigureAwait(false);
    }

    private sealed class ChannelProgress(ChannelWriter<WorkflowEvent> writer) : IProgress<WorkflowEvent>
    {
        public void Report(WorkflowEvent value) => writer.TryWrite(value);
    }

    /// <summary>The report of the most recent run, or null.</summary>
    public WorkflowRunReport? LastReport()
    {
        string path = Path.Combine(StateDirectory, "last-run.json");
        return File.Exists(path) ? WorkflowRunReport.FromJson(File.ReadAllText(path)) : null;
    }

    /// <summary>Graphviz DOT of the DAG for the default targets (or of the rule graph).</summary>
    public string ToDot(bool rules = false) => Build().ToDot(rules);

    /// <summary>Mermaid flowchart of the DAG for the default targets (or of the rule graph).</summary>
    public string ToMermaid(bool rules = false) => Build().ToMermaid(rules);

    private FileStream AcquireLock()
    {
        Directory.CreateDirectory(StateDirectory);
        try
        {
            return new FileStream(Path.Combine(StateDirectory, "lock"), FileMode.OpenOrCreate, FileAccess.ReadWrite, FileShare.None, 1, FileOptions.DeleteOnClose);
        }
        catch (IOException ex)
        {
            throw new InvalidOperationException($"Another run of the workflow in '{WorkingDirectory}' is in progress (the state directory is locked).", ex);
        }
    }

    // ------------------------------------------------------------------ paths

    internal ArtifactKey Key(string pathOrKey)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(pathOrKey);
        if (pathOrKey.StartsWith("value:", StringComparison.Ordinal)) return new ArtifactKey(ArtifactKind.Value, pathOrKey[6..]);
        string p = pathOrKey.Replace('\\', '/');
        string full = Path.GetFullPath(Path.IsPathRooted(p) ? p : Path.Combine(WorkingDirectory, p));
        return new ArtifactKey(ArtifactKind.File, Relative(full));
    }

    internal string FullPath(ArtifactKey key) => Path.GetFullPath(Path.Combine(WorkingDirectory, key.Id));

    internal bool Exists(ArtifactKey key) => key.Kind == ArtifactKind.File
        ? File.Exists(FullPath(key))
        : new StateStore(this).ValueExists(key.Id);

    /// <summary>A path relative to the working directory with forward slashes, or the full path when outside it.</summary>
    internal string Relative(string fullPath)
    {
        string rel = Path.GetRelativePath(WorkingDirectory, fullPath);
        bool outside = rel == ".." || rel.StartsWith(".." + Path.DirectorySeparatorChar, StringComparison.Ordinal) || Path.IsPathRooted(rel);
        return (outside ? fullPath : rel).Replace('\\', '/');
    }

    /// <summary>Normalizes a file output pattern the way <see cref="Key"/> normalizes paths, leaving <c>{...}</c> untouched.</summary>
    internal string NormalizePattern(string pattern)
    {
        var sb = new StringBuilder(pattern.Length);
        int depth = 0;
        foreach (char c in pattern)
        {
            if (c == '{') depth++;
            else if (c == '}') depth = Math.Max(0, depth - 1);
            sb.Append(depth == 0 && c == '\\' ? '/' : c);
        }

        string s = sb.ToString();
        while (s.StartsWith("./", StringComparison.Ordinal)) s = s[2..];
        string root = WorkingDirectory.Replace('\\', '/') + "/";
        var comparison = OperatingSystem.IsWindows() ? StringComparison.OrdinalIgnoreCase : StringComparison.Ordinal;
        if (s.StartsWith(root, comparison)) s = s[root.Length..];
        return s;
    }

    /// <inheritdoc />
    public override string ToString() => $"workflow {Name}: {_rules.Count} rules in {WorkingDirectory}";
}
