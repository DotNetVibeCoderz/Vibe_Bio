namespace Biotech.Net.Data.Caching;

/// <summary>
/// An in-process, least-recently-used response cache bounded by entry count and total bytes.
/// Useful for tests, short-lived tools and as a fast front for <see cref="DiskDataCache"/>.
/// </summary>
public sealed class MemoryDataCache : IDataCache
{
    private readonly object _gate = new();
    private readonly Dictionary<string, LinkedListNode<(string Key, CacheEntry Entry)>> _map = new(StringComparer.Ordinal);
    private readonly LinkedList<(string Key, CacheEntry Entry)> _lru = new();
    private long _bytes;

    /// <param name="maxEntries">Most entries kept (default 1,024).</param>
    /// <param name="maxBytes">Most body bytes kept (default 256 MB).</param>
    public MemoryDataCache(int maxEntries = 1024, long maxBytes = 256L * 1024 * 1024)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxEntries);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxBytes);
        MaxEntries = maxEntries;
        MaxBytes = maxBytes;
    }

    public int MaxEntries { get; }
    public long MaxBytes { get; }

    /// <summary>Entries currently held.</summary>
    public int Count
    {
        get { lock (_gate) return _map.Count; }
    }

    /// <summary>Body bytes currently held.</summary>
    public long SizeBytes
    {
        get { lock (_gate) return _bytes; }
    }

    public ValueTask<CacheEntry?> GetAsync(string key, CancellationToken cancellationToken = default)
    {
        lock (_gate)
        {
            if (!_map.TryGetValue(key, out var node)) return ValueTask.FromResult<CacheEntry?>(null);
            _lru.Remove(node);
            _lru.AddFirst(node);
            return ValueTask.FromResult<CacheEntry?>(node.Value.Entry);
        }
    }

    public ValueTask SetAsync(string key, CacheEntry entry, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(entry);
        if (entry.Content.LongLength > MaxBytes) return ValueTask.CompletedTask; // would evict everything and still not fit
        lock (_gate)
        {
            if (_map.Remove(key, out var old))
            {
                _lru.Remove(old);
                _bytes -= old.Value.Entry.Content.LongLength;
            }

            _map[key] = _lru.AddFirst((key, entry));
            _bytes += entry.Content.LongLength;
            while ((_map.Count > MaxEntries || _bytes > MaxBytes) && _lru.Last is { } last)
            {
                _lru.RemoveLast();
                _map.Remove(last.Value.Key);
                _bytes -= last.Value.Entry.Content.LongLength;
            }
        }

        return ValueTask.CompletedTask;
    }

    public ValueTask<bool> RemoveAsync(string key, CancellationToken cancellationToken = default)
    {
        lock (_gate)
        {
            if (!_map.Remove(key, out var node)) return ValueTask.FromResult(false);
            _lru.Remove(node);
            _bytes -= node.Value.Entry.Content.LongLength;
            return ValueTask.FromResult(true);
        }
    }

    public ValueTask ClearAsync(CancellationToken cancellationToken = default)
    {
        lock (_gate)
        {
            _map.Clear();
            _lru.Clear();
            _bytes = 0;
        }

        return ValueTask.CompletedTask;
    }
}
