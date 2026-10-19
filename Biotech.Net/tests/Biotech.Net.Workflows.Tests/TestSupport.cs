namespace Biotech.Net.Workflows.Tests;

/// <summary>A fresh working directory per test, deleted afterwards.</summary>
internal sealed class TempDir : IDisposable
{
    public TempDir()
    {
        Path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "biotechnet-wf-" + Guid.NewGuid().ToString("N")[..12]);
        Directory.CreateDirectory(Path);
    }

    public string Path { get; }

    public string File(string relative) => System.IO.Path.Combine(Path, relative);

    public void Write(string relative, string text)
    {
        string full = File(relative);
        Directory.CreateDirectory(System.IO.Path.GetDirectoryName(full)!);
        System.IO.File.WriteAllText(full, text);
    }

    public string Read(string relative) => System.IO.File.ReadAllText(File(relative));

    public bool Exists(string relative) => System.IO.File.Exists(File(relative));

    /// <summary>Temporary output folders the engine creates while a job runs.</summary>
    public IEnumerable<string> TempFolders()
        => Directory.EnumerateDirectories(Path, ".biotech-tmp-*", SearchOption.AllDirectories);

    public void Dispose()
    {
        try
        {
            Directory.Delete(Path, recursive: true);
        }
        catch (IOException)
        {
            // A process may still be releasing a handle; the temp folder is cleaned up eventually.
        }
    }
}

/// <summary>
/// A clock that never waits: every timer fires at once and the clock jumps forward by its due time, so retry
/// backoff is deterministic and recorded without sleeping.
/// </summary>
internal sealed class ManualTimeProvider : TimeProvider
{
    private readonly Lock _gate = new();
    private DateTimeOffset _now = new(2026, 1, 1, 0, 0, 0, TimeSpan.Zero);

    public List<TimeSpan> Delays { get; } = [];

    public override DateTimeOffset GetUtcNow()
    {
        lock (_gate) return _now;
    }

    public override ITimer CreateTimer(TimerCallback callback, object? state, TimeSpan dueTime, TimeSpan period)
    {
        if (dueTime != global::System.Threading.Timeout.InfiniteTimeSpan)
        {
            lock (_gate)
            {
                Delays.Add(dueTime);
                _now += dueTime;
            }

            ThreadPool.QueueUserWorkItem(_ => callback(state));
        }

        return new NoTimer();
    }

    private sealed class NoTimer : ITimer
    {
        public bool Change(TimeSpan dueTime, TimeSpan period) => true;
        public void Dispose() { }
        public ValueTask DisposeAsync() => ValueTask.CompletedTask;
    }
}

/// <summary>Collects events synchronously (unlike <see cref="Progress{T}"/>, which posts to the thread pool).</summary>
internal sealed class EventLog : IProgress<WorkflowEvent>
{
    private readonly Lock _gate = new();
    private readonly List<WorkflowEvent> _events = [];

    public IReadOnlyList<WorkflowEvent> Events
    {
        get
        {
            lock (_gate) return [.. _events];
        }
    }

    public void Report(WorkflowEvent value)
    {
        lock (_gate) _events.Add(value);
    }
}
