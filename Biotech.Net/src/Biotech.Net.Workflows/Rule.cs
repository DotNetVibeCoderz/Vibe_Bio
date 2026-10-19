namespace Biotech.Net.Workflows;

/// <summary>Whether an artifact is a file or an in-memory value.</summary>
public enum ArtifactKind
{
    /// <summary>A file, identified by its path (relative to the working directory when inside it).</summary>
    File,

    /// <summary>An in-memory value, identified by its key; persisted as JSON in the state directory.</summary>
    Value,
}

/// <summary>One concrete input or output of a job: a file path or a value key.</summary>
/// <param name="Kind">File or value.</param>
/// <param name="Id">The normalized path (forward slashes) or the value key.</param>
public readonly record struct ArtifactKey(ArtifactKind Kind, string Id)
{
    /// <summary><c>path</c> for files, <c>value:key</c> for values.</summary>
    public override string ToString() => Kind == ArtifactKind.File ? Id : "value:" + Id;
}

/// <summary>
/// A typed handle to in-memory values that rules produce and consume, such as per-sample QC statistics.
/// The key pattern may contain wildcards (<c>qc/{sample}</c>), so one handle names one value per job.
/// </summary>
/// <typeparam name="T">The value type; it must round-trip through <see cref="System.Text.Json.JsonSerializer"/> so that up-to-date jobs can be skipped.</typeparam>
public sealed class ValueRef<T>
{
    internal ValueRef(string name, string pattern)
    {
        Name = name;
        Pattern = new WildcardPattern(pattern);
    }

    /// <summary>The handle's name, also the input/output name it is declared under.</summary>
    public string Name { get; }

    /// <summary>The key pattern.</summary>
    public WildcardPattern Pattern { get; }

    /// <summary>The concrete key for one set of wildcard values.</summary>
    public string Key(IReadOnlyDictionary<string, string> wildcards) => Pattern.Format(wildcards);

    /// <inheritdoc />
    public override string ToString() => $"value {Name} ({Pattern.Text})";
}

/// <summary>How a rule does its work.</summary>
public enum RuleKind
{
    /// <summary>No body: the rule only gathers its inputs (like Snakemake's <c>rule all</c>).</summary>
    Aggregate,

    /// <summary>An in-process C# step.</summary>
    Code,

    /// <summary>A shell command line (<c>cmd /c</c> on Windows, <c>/bin/sh -c</c> elsewhere).</summary>
    Shell,

    /// <summary>An executable started directly with an argument list.</summary>
    Process,
}

internal sealed record InputSpec(string Name, ArtifactKind Kind, Func<Wildcards, IEnumerable<string>> Resolve);

internal sealed record OutputSpec(string Name, ArtifactKind Kind, WildcardPattern Pattern);

/// <summary>
/// A rule: a step with declared inputs and outputs, a body (C# code, a shell command or a process), parameters,
/// resources and failure policy. Configure it fluently; every method returns the rule.
/// </summary>
/// <example>
/// <code>
/// wf.Rule("trim")
///   .Input("reads", "data/{sample}.fastq")
///   .Output("trimmed", "work/{sample}.trimmed.fastq")
///   .Params(new { MinQuality = 20 })
///   .Run(ctx => { /* read ctx.Input("reads"), write ctx.Output("trimmed") */ });
/// </code>
/// </example>
public sealed class Rule
{
    private readonly List<InputSpec> _inputs = [];
    private readonly List<OutputSpec> _outputs = [];
    private readonly Dictionary<string, int> _resources = new(StringComparer.Ordinal);

    internal Rule(Workflow workflow, string name)
    {
        Workflow = workflow;
        Name = name;
    }

    internal Workflow Workflow { get; }

    /// <summary>The rule's name, unique in the workflow.</summary>
    public string Name { get; }

    /// <summary>What the rule does, for reports.</summary>
    public string? Description { get; private set; }

    /// <summary>How the rule does its work.</summary>
    public RuleKind Kind { get; private set; } = RuleKind.Aggregate;

    /// <summary>The code or command version; part of the signature, so bumping it reruns the rule.</summary>
    public string CodeVersion { get; private set; } = "1";

    /// <summary>Threads the rule asks for (default 1).</summary>
    public int ThreadCount { get; private set; } = 1;

    /// <summary>Memory the rule asks for in MB (default 0, i.e. not tracked).</summary>
    public int MemoryMb { get; private set; }

    /// <summary>Custom resources (for example <c>gpu = 1</c>), limited by <see cref="RunOptions.Resources"/>.</summary>
    public IReadOnlyDictionary<string, int> Resources => _resources;

    /// <summary>Extra attempts after a failure (default 0).</summary>
    public int RetryCount { get; private set; }

    /// <summary>Delay before the first retry; each later retry waits <see cref="RetryBackoffFactor"/> times longer.</summary>
    public TimeSpan RetryDelay { get; private set; } = TimeSpan.FromSeconds(1);

    /// <summary>Multiplier between consecutive retry delays (default 2, exponential backoff).</summary>
    public double RetryBackoffFactor { get; private set; } = 2;

    /// <summary>Time limit per attempt, or null for none.</summary>
    public TimeSpan? TimeLimit { get; private set; }

    /// <summary>Wildcard names of the rule (from its outputs).</summary>
    public IReadOnlyList<string> WildcardNames => _outputs.Count == 0 ? [] : _outputs[0].Pattern.WildcardNames;

    internal IReadOnlyList<InputSpec> InputSpecs => _inputs;

    internal IReadOnlyList<OutputSpec> OutputSpecs => _outputs;

    internal Func<Wildcards, object?>? ParamsFactory { get; private set; }

    internal Func<JobContext, Task>? Body { get; private set; }

    internal string? CommandTemplate { get; private set; }

    internal string? ProcessFileName { get; private set; }

    internal IReadOnlyList<string> ProcessArguments { get; private set; } = [];

    /// <summary>Sets the description shown in reports.</summary>
    public Rule Describe(string description)
    {
        Description = description;
        return this;
    }

    /// <summary>
    /// Declares named file inputs. Each pattern may use the rule's wildcards; a pattern list without wildcards
    /// (from <see cref="Expand"/>) gathers many files into one input.
    /// </summary>
    public Rule Input(string name, params string[] patterns) => Input(name, (IEnumerable<string>)patterns);

    /// <summary>Declares named file inputs from a list of patterns (for example the result of <see cref="Expand.Each"/>).</summary>
    public Rule Input(string name, IEnumerable<string> patterns)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        var parsed = patterns.Select(p => new WildcardPattern(p)).ToList();
        return AddInput(new InputSpec(name, ArtifactKind.File, w => parsed.Select(p => p.Format(w))));
    }

    /// <summary>Declares a file input computed from the job's wildcards, for example a sample-sheet lookup.</summary>
    public Rule Input(string name, Func<Wildcards, string> path)
    {
        ArgumentNullException.ThrowIfNull(path);
        return AddInput(new InputSpec(name, ArtifactKind.File, w => [path(w)]));
    }

    /// <summary>Declares file inputs computed from the job's wildcards.</summary>
    public Rule Input(string name, Func<Wildcards, IEnumerable<string>> paths)
    {
        ArgumentNullException.ThrowIfNull(paths);
        return AddInput(new InputSpec(name, ArtifactKind.File, paths));
    }

    /// <summary>Consumes the value <paramref name="value"/> for this job's wildcards; read it with <see cref="JobContext.Get{T}(ValueRef{T})"/>.</summary>
    public Rule Input<T>(ValueRef<T> value)
    {
        ArgumentNullException.ThrowIfNull(value);
        return AddInput(new InputSpec(value.Name, ArtifactKind.Value, w => [value.Key(w)]));
    }

    /// <summary>Gathers one value per row (for example every sample); read them with <see cref="JobContext.GetAll{T}(ValueRef{T})"/>.</summary>
    public Rule Input<T>(ValueRef<T> value, IEnumerable<IReadOnlyDictionary<string, string>> rows)
    {
        ArgumentNullException.ThrowIfNull(value);
        var keys = rows.Select(value.Key).ToList();
        return AddInput(new InputSpec(value.Name, ArtifactKind.Value, _ => keys));
    }

    private Rule AddInput(InputSpec spec)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(spec.Name);
        if (_inputs.Exists(i => i.Name == spec.Name)) throw new ArgumentException($"Rule '{Name}' already has an input named '{spec.Name}'.");
        _inputs.Add(spec);
        return this;
    }

    /// <summary>Declares a named file output. The output patterns' wildcards are the rule's wildcards.</summary>
    public Rule Output(string name, string pattern) => AddOutput(new OutputSpec(name, ArtifactKind.File, new WildcardPattern(Workflow.NormalizePattern(pattern))));

    /// <summary>Declares a value output; set it with <see cref="JobContext.Set{T}(ValueRef{T}, T)"/>.</summary>
    public Rule Output<T>(ValueRef<T> value)
    {
        ArgumentNullException.ThrowIfNull(value);
        return AddOutput(new OutputSpec(value.Name, ArtifactKind.Value, value.Pattern));
    }

    private Rule AddOutput(OutputSpec spec)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(spec.Name);
        if (_outputs.Exists(o => o.Name == spec.Name)) throw new ArgumentException($"Rule '{Name}' already has an output named '{spec.Name}'.");
        if (_outputs.Count > 0)
        {
            var first = _outputs[0].Pattern.WildcardNames.Order(StringComparer.Ordinal);
            if (!first.SequenceEqual(spec.Pattern.WildcardNames.Order(StringComparer.Ordinal)))
                throw new ArgumentException(
                    $"Every output of rule '{Name}' must use the same wildcards: '{_outputs[0].Pattern.Text}' and '{spec.Pattern.Text}' differ.");
        }

        _outputs.Add(spec);
        return this;
    }

    /// <summary>Parameters of the rule; serialized to JSON and hashed, so changing them reruns the rule. Read them with <see cref="JobContext.Params{T}"/>.</summary>
    public Rule Params(object parameters)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        ParamsFactory = _ => parameters;
        return this;
    }

    /// <summary>Parameters that depend on the job's wildcards.</summary>
    public Rule Params(Func<Wildcards, object> parameters)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        ParamsFactory = parameters;
        return this;
    }

    /// <summary>The code version: C# bodies cannot be hashed, so bump this when the code changes.</summary>
    public Rule Version(string version)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(version);
        CodeVersion = version;
        return this;
    }

    /// <summary>Threads the rule needs; counted against <see cref="RunOptions.MaxThreads"/>.</summary>
    public Rule Threads(int threads)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(threads, 1);
        ThreadCount = threads;
        return this;
    }

    /// <summary>Memory the rule needs in MB; counted against <see cref="RunOptions.MaxMemoryMb"/>.</summary>
    public Rule Memory(int megabytes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(megabytes);
        MemoryMb = megabytes;
        return this;
    }

    /// <summary>A custom resource (for example a GPU or a licence seat); counted against <see cref="RunOptions.Resources"/>.</summary>
    public Rule Resource(string name, int amount)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        ArgumentOutOfRangeException.ThrowIfNegative(amount);
        _resources[name] = amount;
        return this;
    }

    /// <summary>Retries a failed attempt up to <paramref name="retries"/> times, waiting <paramref name="delay"/> × <paramref name="backoffFactor"/>^n between attempts.</summary>
    public Rule Retries(int retries, TimeSpan? delay = null, double backoffFactor = 2)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(retries);
        ArgumentOutOfRangeException.ThrowIfLessThan(backoffFactor, 1);
        RetryCount = retries;
        RetryDelay = delay ?? RetryDelay;
        RetryBackoffFactor = backoffFactor;
        return this;
    }

    /// <summary>Fails an attempt that runs longer than <paramref name="limit"/>.</summary>
    public Rule Timeout(TimeSpan limit)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(limit, TimeSpan.Zero);
        TimeLimit = limit;
        return this;
    }

    /// <summary>An in-process C# step.</summary>
    public Rule Run(Func<JobContext, Task> body)
    {
        ArgumentNullException.ThrowIfNull(body);
        return SetBody(RuleKind.Code, body);
    }

    /// <summary>An in-process synchronous C# step.</summary>
    public Rule Run(Action<JobContext> body)
    {
        ArgumentNullException.ThrowIfNull(body);
        return SetBody(RuleKind.Code, ctx =>
        {
            body(ctx);
            return Task.CompletedTask;
        });
    }

    /// <summary>
    /// A shell command. Placeholders: <c>{input}</c>, <c>{input.NAME}</c>, <c>{output}</c>, <c>{output.NAME}</c>
    /// (temporary paths, moved into place on success), <c>{wildcards.NAME}</c>, <c>{threads}</c>, <c>{log}</c>.
    /// Paths are quoted when they contain spaces. The rendered command is part of the signature.
    /// </summary>
    public Rule Shell(string command)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(command);
        CommandTemplate = command;
        return SetBody(RuleKind.Shell, ctx => ProcessRunner.RunShellAsync(ctx, ctx.Render(command, quote: true)));
    }

    /// <summary>
    /// Starts <paramref name="fileName"/> with an argument list (no shell); arguments take the same placeholders as
    /// <see cref="Shell"/> (unquoted, one argument each).
    /// </summary>
    public Rule Process(string fileName, params string[] arguments)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(fileName);
        ProcessFileName = fileName;
        ProcessArguments = [.. arguments];
        CommandTemplate = fileName + " " + string.Join(" ", arguments);
        return SetBody(RuleKind.Process, ctx =>
            ProcessRunner.RunProcessAsync(ctx, ctx.Render(fileName, quote: false), [.. arguments.Select(a => ctx.Render(a, quote: false))]));
    }

    private Rule SetBody(RuleKind kind, Func<JobContext, Task> body)
    {
        if (Body is not null) throw new InvalidOperationException($"Rule '{Name}' already has a body ({Kind}).");
        Kind = kind;
        Body = body;
        return this;
    }

    /// <inheritdoc />
    public override string ToString() => $"rule {Name} ({Kind})";
}
