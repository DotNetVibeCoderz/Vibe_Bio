namespace Biotech.Net.Data.Caching;

/// <summary>One cached HTTP response body plus the few headers the clients need (pagination links, content type).</summary>
public sealed class CacheEntry
{
    public CacheEntry(byte[] content, DateTimeOffset storedAt)
    {
        Content = content ?? throw new ArgumentNullException(nameof(content));
        StoredAt = storedAt;
    }

    /// <summary>The response body, exactly as received (still gzip-compressed for <c>.gz</c> downloads).</summary>
    public byte[] Content { get; }

    /// <summary>When the response was fetched. Freshness is judged against the caller's TTL at read time.</summary>
    public DateTimeOffset StoredAt { get; }

    /// <summary>The request URL, with credentials removed (informational).</summary>
    public string? Url { get; init; }

    /// <summary>The response <c>Content-Type</c>.</summary>
    public string? ContentType { get; init; }

    /// <summary>Selected response headers (<c>Link</c>, <c>X-Total-Results</c>, ...).</summary>
    public IReadOnlyDictionary<string, string> Headers { get; init; } = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

    /// <summary>Age of the entry at <paramref name="now"/>.</summary>
    public TimeSpan Age(DateTimeOffset now) => now - StoredAt;
}

/// <summary>
/// A store of HTTP responses keyed by <see cref="CacheKey"/>. Implementations must be safe for concurrent use;
/// a failed write must never fail the request that produced it.
/// </summary>
public interface IDataCache
{
    /// <summary>The entry for <paramref name="key"/>, or null.</summary>
    ValueTask<CacheEntry?> GetAsync(string key, CancellationToken cancellationToken = default);

    /// <summary>Stores (or replaces) the entry for <paramref name="key"/>.</summary>
    ValueTask SetAsync(string key, CacheEntry entry, CancellationToken cancellationToken = default);

    /// <summary>Removes the entry for <paramref name="key"/>; returns whether one existed.</summary>
    ValueTask<bool> RemoveAsync(string key, CancellationToken cancellationToken = default);

    /// <summary>Removes every entry.</summary>
    ValueTask ClearAsync(CancellationToken cancellationToken = default);
}

/// <summary>The process-wide default cache.</summary>
public static class DataCache
{
    private static readonly Lazy<IDataCache> s_default = new(() => new DiskDataCache(), LazyThreadSafetyMode.ExecutionAndPublication);

    /// <summary>
    /// The on-disk cache at <see cref="DiskDataCache.DefaultDirectory"/> (<c>BIOTECH_NET_CACHE</c>, or
    /// <c>LocalApplicationData/Biotech.Net/cache</c>), created on first use.
    /// </summary>
    public static IDataCache Default => s_default.Value;
}
