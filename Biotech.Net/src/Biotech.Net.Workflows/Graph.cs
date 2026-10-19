using System.Text;

namespace Biotech.Net.Workflows;

/// <summary>The workflow definition or its DAG is invalid: a cycle, an ambiguous or missing producer, a bad pattern.</summary>
public sealed class WorkflowDefinitionException : Exception
{
    /// <summary>Creates the exception.</summary>
    public WorkflowDefinitionException(string message, Exception? inner = null) : base(message, inner) { }
}

/// <summary>One job of the DAG: a rule applied to one set of wildcard values.</summary>
public sealed class WorkflowJob
{
    private readonly List<WorkflowJob> _dependencies = [];
    private readonly List<WorkflowJob> _dependents = [];

    internal WorkflowJob(Rule rule, Wildcards wildcards)
    {
        Rule = rule;
        Wildcards = wildcards;
        Id = wildcards.Count == 0 ? rule.Name : $"{rule.Name}[{wildcards}]";
    }

    /// <summary><c>rule</c> or <c>rule[sample=A]</c>; unique in the graph.</summary>
    public string Id { get; }

    /// <summary>The rule.</summary>
    public Rule Rule { get; }

    /// <summary>The wildcard values.</summary>
    public Wildcards Wildcards { get; }

    /// <summary>Position in topological order (dependencies first).</summary>
    public int Index { get; internal set; } = -1;

    /// <summary>Concrete inputs by input name.</summary>
    public IReadOnlyDictionary<string, IReadOnlyList<ArtifactKey>> Inputs { get; internal set; } = new Dictionary<string, IReadOnlyList<ArtifactKey>>();

    /// <summary>Concrete outputs by output name.</summary>
    public IReadOnlyDictionary<string, ArtifactKey> Outputs { get; internal set; } = new Dictionary<string, ArtifactKey>();

    /// <summary>Jobs whose outputs this job reads.</summary>
    public IReadOnlyList<WorkflowJob> Dependencies => _dependencies;

    /// <summary>Jobs that read this job's outputs.</summary>
    public IReadOnlyList<WorkflowJob> Dependents => _dependents;

    /// <summary>Every input artifact, in declaration order.</summary>
    public IEnumerable<ArtifactKey> AllInputs => Inputs.Values.SelectMany(v => v);

    internal void DependOn(WorkflowJob producer)
    {
        if (_dependencies.Contains(producer)) return;
        _dependencies.Add(producer);
        producer._dependents.Add(this);
    }

    /// <inheritdoc />
    public override string ToString() => Id;
}

/// <summary>
/// The jobs needed for a set of targets and the data dependencies between them, in topological order.
/// </summary>
public sealed class WorkflowGraph
{
    internal WorkflowGraph(Workflow workflow, IReadOnlyList<WorkflowJob> jobs, IReadOnlyList<ArtifactKey> targets,
        IReadOnlyList<ArtifactKey> sources)
    {
        Workflow = workflow;
        Jobs = jobs;
        Targets = targets;
        SourceArtifacts = sources;
    }

    /// <summary>The workflow.</summary>
    public Workflow Workflow { get; }

    /// <summary>Jobs, dependencies before dependents.</summary>
    public IReadOnlyList<WorkflowJob> Jobs { get; }

    /// <summary>The requested target artifacts.</summary>
    public IReadOnlyList<ArtifactKey> Targets { get; }

    /// <summary>Inputs that no rule produces (they must exist before the run).</summary>
    public IReadOnlyList<ArtifactKey> SourceArtifacts { get; }

    /// <summary>The job with id <paramref name="id"/>.</summary>
    public WorkflowJob this[string id] => Jobs.FirstOrDefault(j => j.Id == id) ?? throw new KeyNotFoundException($"No job '{id}' in the graph.");

    /// <summary>Graphviz DOT of the job DAG (or of the rule graph when <paramref name="rules"/> is true).</summary>
    /// <param name="rules">Collapse jobs into one node per rule.</param>
    /// <param name="status">Optional label per job id (for example a plan's reasons); colours the node.</param>
    public string ToDot(bool rules = false, IReadOnlyDictionary<string, JobStatus>? status = null)
    {
        var sb = new StringBuilder();
        sb.Append("digraph \"").Append(Escape(Workflow.Name)).Append("\" {\n");
        sb.Append("  rankdir=LR;\n  node [shape=box, style=\"rounded,filled\", fillcolor=\"#f4f1ea\", fontname=\"Helvetica\"];\n");
        foreach (var (id, label, s) in Nodes(rules, status))
        {
            sb.Append("  \"").Append(Escape(id)).Append("\" [label=\"").Append(Escape(label)).Append('"');
            if (s is { } st) sb.Append(", fillcolor=\"").Append(Colour(st)).Append('"');
            sb.Append("];\n");
        }

        foreach (var (from, to) in Edges(rules)) sb.Append("  \"").Append(Escape(from)).Append("\" -> \"").Append(Escape(to)).Append("\";\n");
        sb.Append("}\n");
        return sb.ToString();
    }

    /// <summary>A Mermaid <c>flowchart</c> of the job DAG (or of the rule graph when <paramref name="rules"/> is true).</summary>
    public string ToMermaid(bool rules = false, IReadOnlyDictionary<string, JobStatus>? status = null)
    {
        var nodes = Nodes(rules, status).ToList();
        var ids = nodes.Select((n, i) => (n.Id, Key: "n" + i)).ToDictionary(x => x.Id, x => x.Key);
        var sb = new StringBuilder("flowchart LR\n");
        foreach (var (id, label, _) in nodes)
            sb.Append("  ").Append(ids[id]).Append("[\"").Append(label.Replace("\"", "#quot;").Replace("\n", "<br/>")).Append("\"]\n");
        foreach (var (from, to) in Edges(rules)) sb.Append("  ").Append(ids[from]).Append(" --> ").Append(ids[to]).Append('\n');
        if (status is not null)
        {
            foreach (var st in Enum.GetValues<JobStatus>())
                sb.Append("  classDef ").Append(st.ToString().ToLowerInvariant()).Append(" fill:").Append(Colour(st)).Append('\n');
            foreach (var (id, _, s) in nodes)
            {
                if (s is { } v) sb.Append("  class ").Append(ids[id]).Append(' ').Append(v.ToString().ToLowerInvariant()).Append('\n');
            }
        }

        return sb.ToString();
    }

    private IEnumerable<(string Id, string Label, JobStatus? Status)> Nodes(bool rules, IReadOnlyDictionary<string, JobStatus>? status)
    {
        if (rules)
        {
            foreach (var g in Jobs.GroupBy(j => j.Rule.Name))
                yield return (g.Key, g.Count() == 1 ? g.Key : $"{g.Key}\n×{g.Count()}", null);
            yield break;
        }

        foreach (var j in Jobs)
        {
            JobStatus? s = status is not null && status.TryGetValue(j.Id, out var v) ? v : null;
            yield return (j.Id, j.Wildcards.Count == 0 ? j.Rule.Name : $"{j.Rule.Name}\n{j.Wildcards}", s);
        }
    }

    private IEnumerable<(string From, string To)> Edges(bool rules)
    {
        if (!rules) return Jobs.SelectMany(j => j.Dependencies.Select(d => (d.Id, j.Id)));
        return Jobs.SelectMany(j => j.Dependencies.Select(d => (d.Rule.Name, j.Rule.Name))).Distinct();
    }

    private static string Escape(string s) => s.Replace("\\", "\\\\").Replace("\"", "\\\"").Replace("\n", "\\n");

    private static string Colour(JobStatus s) => s switch
    {
        JobStatus.Succeeded => "#cde8c9",
        JobStatus.Skipped => "#e3e3e3",
        JobStatus.Failed => "#f3b9b0",
        JobStatus.Cancelled => "#f5deb0",
        JobStatus.WillRun => "#c9dcf2",
        _ => "#f4f1ea",
    };
}

/// <summary>Builds the job DAG backwards from target artifacts, Snakemake style.</summary>
internal sealed class GraphBuilder
{
    private const int MaxDepth = 2000;
    private readonly Workflow _wf;
    private readonly Dictionary<string, WorkflowJob> _jobs = new(StringComparer.Ordinal);
    private readonly Dictionary<ArtifactKey, WorkflowJob?> _producers = [];
    private readonly HashSet<WorkflowJob> _done = [];
    private readonly List<WorkflowJob> _stack = [];
    private readonly List<WorkflowJob> _order = [];
    private readonly List<ArtifactKey> _sources = [];

    public GraphBuilder(Workflow wf) => _wf = wf;

    public WorkflowGraph Build(IReadOnlyList<string> targets)
    {
        var keys = new List<ArtifactKey>();
        foreach (var t in targets)
        {
            var rule = _wf.Rules.FirstOrDefault(r => r.Name == t);
            if (rule is not null)
            {
                if (rule.WildcardNames.Count > 0)
                    throw new WorkflowDefinitionException($"Target '{t}' is a rule with wildcards ({string.Join(", ", rule.WildcardNames)}); request its output files instead.");
                var job = Job(rule, Wildcards.Empty);
                Visit(job, 0);
                keys.AddRange(job.Outputs.Values);
                continue;
            }

            var key = _wf.Key(t);
            keys.Add(key);
            if (Resolve(key, 0) is null && !_wf.Exists(key))
                throw new WorkflowDefinitionException($"Target '{t}' is not produced by any rule and does not exist.");
        }

        for (int i = 0; i < _order.Count; i++) _order[i].Index = i;
        return new WorkflowGraph(_wf, [.. _order], keys, [.. _sources]);
    }

    private WorkflowJob? Resolve(ArtifactKey key, int depth)
    {
        if (_producers.TryGetValue(key, out var known)) return known is null ? null : Visit(known, depth);

        string text = key.Id;
        var matches = new List<(Rule Rule, Wildcards Wildcards)>();
        foreach (var rule in _wf.Rules)
        {
            foreach (var o in rule.OutputSpecs)
            {
                if (o.Kind != key.Kind) continue;
                if (o.Pattern.Match(text) is { } w)
                {
                    matches.Add((rule, w));
                    break;
                }
            }
        }

        if (matches.Count > 1)
            throw new WorkflowDefinitionException(
                $"'{key}' can be produced by more than one rule: {string.Join(", ", matches.Select(m => $"{m.Rule.Name}[{m.Wildcards}]"))}. Make the output patterns distinct, for example with wildcard constraints like {{sample,[A-Za-z0-9]+}}.");
        if (matches.Count == 0)
        {
            _producers[key] = null;
            _sources.Add(key);
            return null;
        }

        var job = Job(matches[0].Rule, matches[0].Wildcards);
        if (!job.Outputs.Values.Contains(key))
            throw new WorkflowDefinitionException($"Rule '{job.Rule.Name}' matched '{key}' but formats its outputs differently; check the output patterns.");
        return Visit(job, depth);
    }

    private WorkflowJob Job(Rule rule, Wildcards wildcards)
    {
        string id = wildcards.Count == 0 ? rule.Name : $"{rule.Name}[{wildcards}]";
        if (_jobs.TryGetValue(id, out var existing)) return existing;
        var job = new WorkflowJob(rule, wildcards);
        try
        {
            job.Outputs = rule.OutputSpecs.ToDictionary(o => o.Name,
                o => o.Kind == ArtifactKind.File ? _wf.Key(o.Pattern.Format(wildcards)) : new ArtifactKey(ArtifactKind.Value, o.Pattern.Format(wildcards)));
            job.Inputs = rule.InputSpecs.ToDictionary(i => i.Name,
                i => (IReadOnlyList<ArtifactKey>)[.. i.Resolve(wildcards).Select(p => i.Kind == ArtifactKind.File ? _wf.Key(p) : new ArtifactKey(ArtifactKind.Value, p))]);
        }
        catch (Exception ex) when (ex is KeyNotFoundException or FormatException or ArgumentException)
        {
            throw new WorkflowDefinitionException($"Cannot expand job {id}: {ex.Message}", ex);
        }

        _jobs[id] = job;
        foreach (var o in job.Outputs.Values)
        {
            if (_producers.TryGetValue(o, out var other) && other is not null && other != job)
                throw new WorkflowDefinitionException($"'{o}' is an output of both {other.Id} and {job.Id}.");
            _producers[o] = job;
        }

        return job;
    }

    private WorkflowJob Visit(WorkflowJob job, int depth)
    {
        if (_done.Contains(job)) return job;
        int at = _stack.IndexOf(job);
        if (at >= 0)
        {
            var cycle = _stack.Skip(at).Append(job).Select(j => j.Id);
            throw new WorkflowDefinitionException($"The workflow has a cycle: {string.Join(" -> ", cycle)}.");
        }

        if (depth > MaxDepth)
            throw new WorkflowDefinitionException($"Dependency chain deeper than {MaxDepth} at {job.Id}; an input pattern probably matches its own rule's output (pattern recursion).");

        _stack.Add(job);
        foreach (var input in job.AllInputs)
        {
            if (Resolve(input, depth + 1) is { } producer) job.DependOn(producer);
        }

        _stack.RemoveAt(_stack.Count - 1);
        _done.Add(job);
        _order.Add(job);
        return job;
    }
}
