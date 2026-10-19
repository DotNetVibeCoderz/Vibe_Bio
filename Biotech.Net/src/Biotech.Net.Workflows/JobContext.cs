using System.Diagnostics;
using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;

namespace Biotech.Net.Workflows;

/// <summary>A shell or process step exited with a non-zero code.</summary>
public sealed class ProcessFailedException : Exception
{
    /// <summary>Creates the exception.</summary>
    public ProcessFailedException(string command, int exitCode, string stderrTail)
        : base($"'{command}' exited with code {exitCode}." + (stderrTail.Length > 0 ? " stderr: " + stderrTail : ""))
    {
        Command = command;
        ExitCode = exitCode;
    }

    /// <summary>The command line.</summary>
    public string Command { get; }

    /// <summary>The exit code.</summary>
    public int ExitCode { get; }
}

/// <summary>
/// What a running job sees: its wildcards, input paths and values, temporary output paths, parameters, threads, a log
/// and a cancellation token. Outputs written here are moved into place only when the job succeeds.
/// </summary>
public sealed partial class JobContext
{
    private readonly Executor _executor;
    private readonly Dictionary<string, string> _tempPaths;
    private readonly Dictionary<string, byte[]> _values = new(StringComparer.Ordinal);
    private readonly object? _params;

    internal JobContext(Executor executor, WorkflowJob job, int attempt, int threads, string logPath,
        Dictionary<string, string> tempPaths, object? parameters, CancellationToken cancellationToken)
    {
        _executor = executor;
        Job = job;
        Attempt = attempt;
        Threads = threads;
        LogPath = logPath;
        _tempPaths = tempPaths;
        _params = parameters;
        CancellationToken = cancellationToken;
    }

    /// <summary>The job.</summary>
    public WorkflowJob Job { get; }

    /// <summary>The job id, for example <c>trim[sample=A]</c>.</summary>
    public string JobId => Job.Id;

    /// <summary>The job's wildcard values.</summary>
    public Wildcards Wildcards => Job.Wildcards;

    /// <summary>1 for the first attempt, 2 for the first retry, …</summary>
    public int Attempt { get; }

    /// <summary>Threads granted to the job (its request, capped by <see cref="RunOptions.MaxThreads"/>).</summary>
    public int Threads { get; }

    /// <summary>The job's log file.</summary>
    public string LogPath { get; }

    /// <summary>The workflow's working directory (also the working directory of shell and process steps).</summary>
    public string WorkingDirectory => _executor.Workflow.WorkingDirectory;

    /// <summary>Cancelled when the run is cancelled or the attempt times out.</summary>
    public CancellationToken CancellationToken { get; }

    /// <summary>Exit code of the last process this job ran, if any.</summary>
    public int? ExitCode { get; internal set; }

    internal IReadOnlyDictionary<string, byte[]> PendingValues => _values;

    /// <summary>The full path of the single file of input <paramref name="name"/>.</summary>
    /// <exception cref="KeyNotFoundException">No such input.</exception>
    /// <exception cref="InvalidOperationException">The input holds several files (use <see cref="Inputs"/>) or values.</exception>
    public string Input(string name)
    {
        var list = Inputs(name);
        return list.Count == 1 ? list[0] : throw new InvalidOperationException($"Input '{name}' of {JobId} has {list.Count} files; use Inputs(\"{name}\").");
    }

    /// <summary>Full paths of every file of input <paramref name="name"/>, in declaration order.</summary>
    public IReadOnlyList<string> Inputs(string name)
    {
        var keys = InputKeys(name);
        if (keys.Any(k => k.Kind != ArtifactKind.File)) throw new InvalidOperationException($"Input '{name}' of {JobId} is a value; use Get or GetAll.");
        return [.. keys.Select(_executor.Workflow.FullPath)];
    }

    private IReadOnlyList<ArtifactKey> InputKeys(string name)
        => Job.Inputs.TryGetValue(name, out var keys)
            ? keys
            : throw new KeyNotFoundException($"{JobId} has no input '{name}' (inputs: {string.Join(", ", Job.Inputs.Keys)}).");

    /// <summary>
    /// The temporary full path to write output <paramref name="name"/> to. It is moved to the real path when the job
    /// succeeds and deleted when it fails, so a failure never leaves a half-written output behind.
    /// </summary>
    public string Output(string name)
        => _tempPaths.TryGetValue(name, out var p)
            ? p
            : throw new KeyNotFoundException($"{JobId} has no file output '{name}' (outputs: {string.Join(", ", Job.Outputs.Keys)}).");

    /// <summary>The final (not temporary) path of output <paramref name="name"/>, for messages and reports.</summary>
    public string FinalOutput(string name) => _executor.Workflow.FullPath(Job.Outputs[name]);

    /// <summary>The value of input <paramref name="value"/> for this job.</summary>
    public T Get<T>(ValueRef<T> value)
    {
        var all = GetAll(value);
        return all.Count == 1 ? all[0] : throw new InvalidOperationException($"Input '{value.Name}' of {JobId} has {all.Count} values; use GetAll.");
    }

    /// <summary>Every value of input <paramref name="value"/>, in declaration order (the gather side).</summary>
    public IReadOnlyList<T> GetAll<T>(ValueRef<T> value)
    {
        ArgumentNullException.ThrowIfNull(value);
        var keys = InputKeys(value.Name);
        return [.. keys.Select(k => JsonSerializer.Deserialize<T>(_executor.Store.ReadValue(k.Id))!)];
    }

    /// <summary>Every value of input <paramref name="value"/> by its key.</summary>
    public IReadOnlyDictionary<string, T> GetAllByKey<T>(ValueRef<T> value)
    {
        var keys = InputKeys(value.Name);
        return keys.ToDictionary(k => k.Id, k => JsonSerializer.Deserialize<T>(_executor.Store.ReadValue(k.Id))!);
    }

    /// <summary>Sets value output <paramref name="value"/>; it is stored when the job succeeds.</summary>
    public void Set<T>(ValueRef<T> value, T result)
    {
        ArgumentNullException.ThrowIfNull(value);
        if (!Job.Outputs.TryGetValue(value.Name, out var key) || key.Kind != ArtifactKind.Value)
            throw new KeyNotFoundException($"{JobId} has no value output '{value.Name}'.");
        _values[value.Name] = JsonSerializer.SerializeToUtf8Bytes(result);
    }

    /// <summary>The rule's parameters, converted to <typeparamref name="T"/> through JSON (or cast when already a <typeparamref name="T"/>).</summary>
    public T Params<T>()
    {
        if (_params is T t) return t;
        if (_params is null) throw new InvalidOperationException($"Rule '{Job.Rule.Name}' has no parameters.");
        return JsonSerializer.Deserialize<T>(JsonSerializer.Serialize(_params))!;
    }

    /// <summary>Writes a line to the job's log and raises a <see cref="JobLogEvent"/>.</summary>
    public void Log(string message) => _executor.Log(this, message);

    /// <summary>Runs an executable (no shell) in the working directory, logging its output; throws <see cref="ProcessFailedException"/> on a non-zero exit unless <paramref name="throwOnError"/> is false.</summary>
    public Task<int> RunProcessAsync(string fileName, IEnumerable<string> arguments, bool throwOnError = true)
        => ProcessRunner.RunAsync(this, ProcessRunner.Direct(fileName, arguments), throwOnError);

    /// <summary>Runs a shell command line (<c>cmd /c</c> or <c>/bin/sh -c</c>) in the working directory, logging its output.</summary>
    public Task<int> RunShellAsync(string command, bool throwOnError = true)
        => ProcessRunner.RunAsync(this, ProcessRunner.ShellStart(command), throwOnError);

    /// <summary>Substitutes <c>{input}</c>, <c>{input.NAME}</c>, <c>{output}</c>, <c>{output.NAME}</c>, <c>{wildcards.NAME}</c>, <c>{threads}</c>, <c>{log}</c> and <c>{{ }}</c>.</summary>
    /// <param name="template">The template.</param>
    /// <param name="quote">For shell command lines: quote paths that contain spaces, and use backslashes on Windows.</param>
    /// <exception cref="FormatException">An unknown placeholder.</exception>
    public string Render(string template, bool quote)
    {
        ArgumentNullException.ThrowIfNull(template);
        string Q(string p)
        {
            if (quote && OperatingSystem.IsWindows()) p = p.Replace('/', '\\');   // cmd.exe reads '/' as a switch
            return quote && p.Any(char.IsWhiteSpace) ? "\"" + p + "\"" : p;
        }

        string Rel(string full) => _executor.Workflow.Relative(full);
        string Files(IEnumerable<ArtifactKey> keys) => string.Join(" ", keys.Select(k => k.Kind == ArtifactKind.File
            ? Q(Rel(_executor.Workflow.FullPath(k)))
            : throw new FormatException($"Value input '{k.Id}' cannot be passed to a command; read it in a C# step.")));

        string text = template.Replace("{{", "\u0001").Replace("}}", "\u0002");
        text = Placeholder().Replace(text, m =>
        {
            string head = m.Groups[1].Value, tail = m.Groups[2].Value;
            return (head, tail) switch
            {
                ("input", "") => Files(Job.AllInputs),
                ("input", _) => Files(InputKeys(tail)),
                ("output", "") => string.Join(" ", _tempPaths.Values.Select(p => Q(Rel(p)))),
                ("output", _) => Q(Rel(Output(tail))),
                ("wildcards", _) => Wildcards[tail],
                ("threads", "") => Threads.ToString(CultureInfo.InvariantCulture),
                ("log", "") => Q(LogPath),
                _ => throw new FormatException($"Unknown placeholder '{m.Value}' in '{template}'."),
            };
        });
        return text.Replace('\u0001', '{').Replace('\u0002', '}');
    }

    [GeneratedRegex(@"\{([A-Za-z_]\w*)(?:\.([A-Za-z_]\w*))?\}")]
    private static partial Regex Placeholder();
}

/// <summary>Starts external tools, streams their output to the job log, and kills them on cancellation.</summary>
internal static class ProcessRunner
{
    public static Task RunShellAsync(JobContext ctx, string command) => RunAsync(ctx, ShellStart(command), throwOnError: true);

    public static Task RunProcessAsync(JobContext ctx, string fileName, IReadOnlyList<string> arguments)
        => RunAsync(ctx, Direct(fileName, arguments), throwOnError: true);

    public static ProcessStartInfo ShellStart(string command)
    {
        if (OperatingSystem.IsWindows())
            return new ProcessStartInfo("cmd.exe") { Arguments = "/d /s /c \"" + command + "\"" };
        var psi = new ProcessStartInfo("/bin/sh");
        psi.ArgumentList.Add("-c");
        psi.ArgumentList.Add(command);
        return psi;
    }

    public static ProcessStartInfo Direct(string fileName, IEnumerable<string> arguments)
    {
        var psi = new ProcessStartInfo(fileName);
        foreach (var a in arguments) psi.ArgumentList.Add(a);
        return psi;
    }

    public static async Task<int> RunAsync(JobContext ctx, ProcessStartInfo psi, bool throwOnError)
    {
        psi.WorkingDirectory = ctx.WorkingDirectory;
        psi.UseShellExecute = false;
        psi.RedirectStandardOutput = true;
        psi.RedirectStandardError = true;
        psi.CreateNoWindow = true;
        psi.StandardOutputEncoding = Encoding.UTF8;
        psi.StandardErrorEncoding = Encoding.UTF8;
        psi.Environment["BIOTECH_WORKFLOW_THREADS"] = ctx.Threads.ToString(CultureInfo.InvariantCulture);
        string display = psi.ArgumentList.Count > 0 ? psi.FileName + " " + string.Join(" ", psi.ArgumentList) : psi.FileName + " " + psi.Arguments;
        ctx.Log("$ " + display);

        var tail = new Queue<string>();
        using var process = new Process { StartInfo = psi };
        process.OutputDataReceived += (_, e) =>
        {
            if (e.Data is not null) ctx.Log("[out] " + e.Data);
        };
        process.ErrorDataReceived += (_, e) =>
        {
            if (e.Data is null) return;
            ctx.Log("[err] " + e.Data);
            lock (tail)
            {
                tail.Enqueue(e.Data);
                if (tail.Count > 5) tail.Dequeue();
            }
        };

        try
        {
            process.Start();
        }
        catch (System.ComponentModel.Win32Exception ex)
        {
            throw new ProcessFailedException(display, -1, "could not start: " + ex.Message);
        }

        process.BeginOutputReadLine();
        process.BeginErrorReadLine();
        try
        {
            await process.WaitForExitAsync(ctx.CancellationToken).ConfigureAwait(false);
        }
        catch (OperationCanceledException)
        {
            try
            {
                process.Kill(entireProcessTree: true);
            }
            catch (InvalidOperationException)
            {
                // Already exited.
            }

            throw;
        }

        int code = process.ExitCode;
        ctx.ExitCode = code;
        ctx.Log($"exit code {code.ToString(CultureInfo.InvariantCulture)}");
        if (code != 0 && throwOnError)
        {
            string err;
            lock (tail) err = string.Join(" | ", tail);
            throw new ProcessFailedException(display, code, err);
        }

        return code;
    }
}
