using Biotech.Net.Data.Caching;

namespace Biotech.Net.Data.Http;

/// <summary>How a data client uses its <see cref="IDataCache"/>.</summary>
public enum DataCacheMode
{
    /// <summary>Serve a cached response while it is younger than the time-to-live, otherwise fetch and store.</summary>
    Default,

    /// <summary>Always fetch from the network and overwrite the cached copy.</summary>
    Refresh,

    /// <summary>Never read or write the cache.</summary>
    Bypass,

    /// <summary>
    /// Never touch the network. Serve any cached copy, however old; throw <see cref="DataOfflineException"/>
    /// when nothing is cached.
    /// </summary>
    Offline,
}

/// <summary>
/// Settings shared by every Biotech.Net data client: cache, retries, timeouts and rate limiting.
/// A <c>record</c>, so <c>options with { CacheMode = DataCacheMode.Offline }</c> derives a variant.
/// </summary>
public sealed record DataClientOptions
{
    /// <summary>Environment variable that, when set to <c>1</c> or <c>true</c>, makes <see cref="DataCacheMode.Offline"/> the default mode.</summary>
    public const string OfflineEnvironmentVariable = "BIOTECH_NET_OFFLINE";

    /// <summary>The defaults: the shared on-disk cache, four retries, a 100 s timeout and the shared host rate limiter.</summary>
    public static DataClientOptions Default { get; } = new();

    /// <summary>
    /// The response cache. <c>null</c> (the default) means <see cref="DataCache.Default"/>, the on-disk cache
    /// under <see cref="DiskDataCache.DefaultDirectory"/>.
    /// </summary>
    public IDataCache? Cache { get; init; }

    /// <summary>How the cache is used. Defaults to <see cref="DataCacheMode.Offline"/> when <c>BIOTECH_NET_OFFLINE=1</c>.</summary>
    public DataCacheMode CacheMode { get; init; } = DefaultMode();

    /// <summary>
    /// Overrides every per-endpoint time-to-live (sequence records 30 days, searches 1 day, ...). <c>null</c> keeps
    /// the per-endpoint defaults.
    /// </summary>
    public TimeSpan? CacheTtl { get; init; }

    /// <summary>Clock used for TTLs, backoff delays, timeouts and rate limiting. Replace it in tests.</summary>
    public TimeProvider TimeProvider { get; init; } = TimeProvider.System;

    /// <summary>Retries after the first attempt for 408/429/5xx responses, timeouts and network errors.</summary>
    public int MaxRetries { get; init; } = 4;

    /// <summary>Backoff before the first retry; doubles on each further retry, with jitter.</summary>
    public TimeSpan InitialBackoff { get; init; } = TimeSpan.FromSeconds(1);

    /// <summary>Upper bound on the exponential backoff (before jitter).</summary>
    public TimeSpan MaxBackoff { get; init; } = TimeSpan.FromSeconds(30);

    /// <summary>A server <c>Retry-After</c> longer than this fails the request instead of waiting.</summary>
    public TimeSpan MaxRetryAfter { get; init; } = TimeSpan.FromMinutes(2);

    /// <summary>Timeout of one attempt (connect, send and read the whole body).</summary>
    public TimeSpan Timeout { get; init; } = TimeSpan.FromSeconds(100);

    /// <summary>
    /// Per-host request spacing. <c>null</c> means <see cref="HostRateLimiter.Shared"/>, so every client in the
    /// process shares one budget per host (NCBI counts requests per IP address, not per client object).
    /// </summary>
    public HostRateLimiter? RateLimiter { get; init; }

    /// <summary>Overrides the per-service request rate (NCBI 3/s or 10/s with a key, Ensembl 15/s, others 10/s).</summary>
    public double? RequestsPerSecond { get; init; }

    /// <summary>Source of backoff jitter. Defaults to <see cref="Random.Shared"/>.</summary>
    public Random? Random { get; init; }

    /// <summary>The <c>User-Agent</c> sent with every request.</summary>
    public string UserAgent { get; init; } = $"Biotech.Net/{typeof(DataClientOptions).Assembly.GetName().Version?.ToString(3) ?? "0"} (+https://github.com/DotNetVibeCoderz/Vibe_Bio)";

    internal IDataCache ResolveCache() => Cache ?? DataCache.Default;

    internal HostRateLimiter ResolveRateLimiter() => RateLimiter ?? HostRateLimiter.Shared;

    private static DataCacheMode DefaultMode()
    {
        string? v = Environment.GetEnvironmentVariable(OfflineEnvironmentVariable);
        return v is "1" || string.Equals(v, "true", StringComparison.OrdinalIgnoreCase) ? DataCacheMode.Offline : DataCacheMode.Default;
    }
}
