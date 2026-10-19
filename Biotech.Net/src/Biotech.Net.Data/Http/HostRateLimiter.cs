using System.Collections.Concurrent;

namespace Biotech.Net.Data.Http;

/// <summary>
/// Spaces requests to each host evenly (a request may start every <c>1 / rate</c> seconds) and honours
/// server back-off signals (<c>Retry-After</c>, Ensembl's <c>X-RateLimit-Reset</c>) by pausing the whole host.
/// </summary>
/// <remarks>
/// Slots are reserved under a per-host lock, so concurrent callers queue fairly instead of bursting.
/// All waiting goes through the supplied <see cref="TimeProvider"/>, which makes the behaviour testable
/// with a fake clock.
/// </remarks>
public sealed class HostRateLimiter
{
    private readonly ConcurrentDictionary<string, HostState> _hosts = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>The process-wide limiter used when <see cref="DataClientOptions.RateLimiter"/> is not set.</summary>
    public static HostRateLimiter Shared { get; } = new();

    /// <summary>Waits until a request to <paramref name="host"/> may start at <paramref name="requestsPerSecond"/>.</summary>
    /// <returns>The time spent waiting.</returns>
    public async ValueTask<TimeSpan> WaitAsync(string host, double requestsPerSecond, TimeProvider timeProvider, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(timeProvider);
        var state = _hosts.GetOrAdd(host, static _ => new HostState());
        TimeSpan wait;
        lock (state)
        {
            var now = timeProvider.GetUtcNow();
            var slot = Max(now, Max(state.NextSlot, state.BlockedUntil));
            var spacing = requestsPerSecond > 0 && !double.IsInfinity(requestsPerSecond)
                ? TimeSpan.FromTicks((long)Math.Round(TimeSpan.TicksPerSecond / requestsPerSecond))
                : TimeSpan.Zero;
            state.NextSlot = slot + spacing;
            wait = slot - now;
        }

        // Timers have millisecond resolution and Task.Delay truncates; round up so spacing never undershoots.
        if (wait > TimeSpan.Zero) wait = TimeSpan.FromMilliseconds(Math.Ceiling(wait.TotalMilliseconds));

        if (wait > TimeSpan.Zero) await Task.Delay(wait, timeProvider, cancellationToken).ConfigureAwait(false);
        return wait > TimeSpan.Zero ? wait : TimeSpan.Zero;
    }

    /// <summary>Stops every request to <paramref name="host"/> from starting before <paramref name="until"/>.</summary>
    public void BlockUntil(string host, DateTimeOffset until)
    {
        var state = _hosts.GetOrAdd(host, static _ => new HostState());
        lock (state)
        {
            if (until > state.BlockedUntil) state.BlockedUntil = until;
        }
    }

    /// <summary>When <paramref name="host"/> is blocked by a server back-off signal, the time it is released; otherwise null.</summary>
    public DateTimeOffset? BlockedUntil(string host)
        => _hosts.TryGetValue(host, out var s) && s.BlockedUntil != default ? s.BlockedUntil : null;

    private static DateTimeOffset Max(DateTimeOffset a, DateTimeOffset b) => a > b ? a : b;

    private sealed class HostState
    {
        public DateTimeOffset NextSlot;
        public DateTimeOffset BlockedUntil;
    }
}
