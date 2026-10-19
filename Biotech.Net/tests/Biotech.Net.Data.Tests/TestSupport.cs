using System.Collections.Concurrent;
using System.Net;
using System.Net.Http.Headers;
using System.Text;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.Tests;

/// <summary>A request the fake handler saw.</summary>
public sealed record SeenRequest(HttpMethod Method, Uri Uri, string? Body, string? ContentType, string? UserAgent)
{
    public string Path => Uri.AbsolutePath;
    public IReadOnlyDictionary<string, List<string>> Query => FakeHttpHandler.ParseQuery(Uri);
    public string? Param(string name) => Query.TryGetValue(name, out var v) ? v[0] : null;
}

/// <summary>
/// Replays recorded responses. Routes are tried in order; an unmatched request gets a 404 that names the URL,
/// so a missing fixture fails loudly instead of hitting the network.
/// </summary>
public sealed class FakeHttpHandler : HttpMessageHandler
{
    private readonly List<(Func<HttpRequestMessage, bool> Match, Func<HttpRequestMessage, CancellationToken, Task<HttpResponseMessage>> Respond)> _routes = [];

    public ConcurrentQueue<SeenRequest> Requests { get; } = new();

    public int Count => Requests.Count;

    public FakeHttpHandler On(Func<HttpRequestMessage, bool> match, Func<HttpRequestMessage, CancellationToken, Task<HttpResponseMessage>> respond)
    {
        _routes.Add((match, respond));
        return this;
    }

    public FakeHttpHandler On(Func<HttpRequestMessage, bool> match, Func<HttpResponseMessage> respond) =>
        On(match, (_, _) => Task.FromResult(respond()));

    /// <summary>Answers the matching requests with the given responses in order (the last one repeats).</summary>
    public FakeHttpHandler OnSequence(Func<HttpRequestMessage, bool> match, params Func<HttpResponseMessage>[] responses)
    {
        int i = 0;
        return On(match, () => responses[Math.Min(Interlocked.Increment(ref i) - 1, responses.Length - 1)]());
    }

    /// <summary>Answers with a recorded fixture file (plus its <c>.headers</c> file when one exists).</summary>
    public FakeHttpHandler OnFixture(Func<HttpRequestMessage, bool> match, string fixture, string contentType = "application/json") =>
        On(match, () => Fixtures.Response(fixture, contentType));

    public HttpClient Client() => new(this, disposeHandler: false);

    protected override async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
    {
        string? body = request.Content is null ? null : await request.Content.ReadAsStringAsync(cancellationToken);
        Requests.Enqueue(new SeenRequest(request.Method, request.RequestUri!, body, request.Content?.Headers.ContentType?.MediaType,
            request.Headers.UserAgent.ToString()));
        foreach (var (match, respond) in _routes)
            if (match(request)) return await respond(request, cancellationToken);
        return new HttpResponseMessage(HttpStatusCode.NotFound) { Content = new StringContent($"{{\"error\":\"no fixture for {request.RequestUri}\"}}") };
    }

    public static Dictionary<string, List<string>> ParseQuery(Uri uri)
    {
        var d = new Dictionary<string, List<string>>(StringComparer.Ordinal);
        foreach (string pair in uri.Query.TrimStart('?').Split(['&', ';'], StringSplitOptions.RemoveEmptyEntries))
        {
            var kv = pair.Split('=', 2);
            string k = Uri.UnescapeDataString(kv[0]);
            string v = kv.Length > 1 ? Uri.UnescapeDataString(kv[1].Replace('+', ' ')) : "";
            if (!d.TryGetValue(k, out var list)) d[k] = list = [];
            list.Add(v);
        }

        return d;
    }

    public static Func<HttpRequestMessage, bool> Path(string pathSuffix, Func<Dictionary<string, List<string>>, bool>? query = null) =>
        r => r.RequestUri!.AbsolutePath.EndsWith(pathSuffix, StringComparison.Ordinal) && (query is null || query(ParseQuery(r.RequestUri)));
}

public static class Fixtures
{
    public static string Dir => System.IO.Path.Combine(AppContext.BaseDirectory, "Fixtures");

    public static string File(string name) => System.IO.Path.Combine(Dir, name);

    public static byte[] Bytes(string name) => System.IO.File.ReadAllBytes(File(name));

    public static string Text(string name) => System.IO.File.ReadAllText(File(name));

    public static HttpResponseMessage Response(string name, string contentType = "application/json")
    {
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new ByteArrayContent(Bytes(name)) };
        response.Content.Headers.ContentType = MediaTypeHeaderValue.Parse(contentType);
        string headers = File(System.IO.Path.ChangeExtension(name, ".headers"));
        if (System.IO.File.Exists(headers))
        {
            foreach (string line in System.IO.File.ReadAllLines(headers))
            {
                int colon = line.IndexOf(':');
                if (colon <= 0) continue;
                string key = line[..colon].Trim();
                if (key.Equals("Link", StringComparison.OrdinalIgnoreCase) || key.Equals("X-Total-Results", StringComparison.OrdinalIgnoreCase))
                    response.Headers.TryAddWithoutValidation(key, line[(colon + 1)..].Trim());
            }
        }

        return response;
    }

    public static HttpResponseMessage Json(string json, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(json, Encoding.UTF8, "application/json") };

    public static HttpResponseMessage Status(HttpStatusCode status, string body = "") =>
        new(status) { Content = new StringContent(body) };
}

/// <summary>
/// A deterministic clock. Timers whose due time is at most <see cref="AutoFireLimit"/> fire on the thread pool right
/// away and move the clock forward by their due time, so backoff and rate-limit waits take no real time but are
/// recorded in <see cref="Delays"/>. Longer timers (e.g. a one-hour request timeout) never fire on their own.
/// </summary>
public sealed class FakeClock : TimeProvider
{
    private readonly object _gate = new();
    private readonly List<Action> _paused = [];
    private DateTimeOffset _now = new(2026, 1, 1, 0, 0, 0, TimeSpan.Zero);
    private bool _isPaused;

    public TimeSpan AutoFireLimit { get; set; } = TimeSpan.FromMinutes(10);

    /// <summary>Due times of every auto-fired timer, in creation order.</summary>
    public ConcurrentQueue<TimeSpan> Delays { get; } = new();

    public override DateTimeOffset GetUtcNow()
    {
        lock (_gate) return _now;
    }

    public void Advance(TimeSpan by)
    {
        lock (_gate) _now += by;
    }

    /// <summary>Holds auto-fired timers until <see cref="Resume"/>.</summary>
    public void Pause()
    {
        lock (_gate) _isPaused = true;
    }

    public void Resume()
    {
        List<Action> pending;
        lock (_gate)
        {
            _isPaused = false;
            pending = [.. _paused];
            _paused.Clear();
        }

        foreach (var a in pending) ThreadPool.QueueUserWorkItem(_ => a());
    }

    public override ITimer CreateTimer(TimerCallback callback, object? state, TimeSpan dueTime, TimeSpan period)
    {
        var t = new FakeTimer(this, callback, state);
        t.Change(dueTime, period);
        return t;
    }

    private void Schedule(FakeTimer timer, int generation, DateTimeOffset due, TimeSpan delay)
    {
        Delays.Enqueue(delay);
        void Fire()
        {
            lock (_gate)
                if (due > _now) _now = due;
            timer.Fire(generation);
        }

        lock (_gate)
        {
            if (_isPaused)
            {
                _paused.Add(Fire);
                return;
            }
        }

        ThreadPool.QueueUserWorkItem(_ => Fire());
    }

    private sealed class FakeTimer(FakeClock clock, TimerCallback callback, object? state) : ITimer
    {
        private int _generation;
        private bool _disposed;

        public bool Change(TimeSpan dueTime, TimeSpan period)
        {
            int gen = Interlocked.Increment(ref _generation);
            if (_disposed || dueTime == Timeout.InfiniteTimeSpan || dueTime > clock.AutoFireLimit) return true;
            clock.Schedule(this, gen, clock.GetUtcNow() + dueTime, dueTime);
            return true;
        }

        public void Fire(int generation)
        {
            if (_disposed || generation != Volatile.Read(ref _generation)) return;
            callback(state);
        }

        public void Dispose() => _disposed = true;

        public ValueTask DisposeAsync()
        {
            _disposed = true;
            return ValueTask.CompletedTask;
        }
    }
}

/// <summary>Builds options wired to a fake clock, a private limiter and an in-memory cache.</summary>
public static class TestOptions
{
    public static DataClientOptions Create(FakeClock clock, IDataCache? cache = null, DataCacheMode mode = DataCacheMode.Default) => new()
    {
        Cache = cache ?? new MemoryDataCache(),
        CacheMode = mode,
        TimeProvider = clock,
        RateLimiter = new HostRateLimiter(),
        Timeout = TimeSpan.FromHours(1),
        Random = new Random(1234),
    };
}

/// <summary>A temporary directory deleted at dispose.</summary>
public sealed class TempDir : IDisposable
{
    public TempDir() => System.IO.Directory.CreateDirectory(Path);

    public string Path { get; } = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "biotechnet-data-tests", Guid.NewGuid().ToString("N"));

    public void Dispose()
    {
        try
        {
            System.IO.Directory.Delete(Path, recursive: true);
        }
        catch (IOException)
        {
        }
    }
}
