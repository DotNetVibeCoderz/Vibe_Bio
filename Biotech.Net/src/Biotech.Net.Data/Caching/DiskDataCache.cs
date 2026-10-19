using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace Biotech.Net.Data.Caching;

/// <summary>
/// A response cache on disk, one file per entry, safe for concurrent use by several threads and processes.
/// </summary>
/// <remarks>
/// <para>
/// Layout: <c>&lt;directory&gt;/&lt;first two hex digits&gt;/&lt;key&gt;.entry</c>. Each file holds a 4-byte magic
/// (<c>BND1</c>), a 4-byte little-endian metadata length, the UTF-8 JSON metadata (URL, content type, stored-at
/// time, headers) and the raw body.
/// </para>
/// <para>
/// Writes are atomic: the entry is written to a uniquely named temporary file in the same directory, flushed,
/// then moved over the final name, so a reader sees the old entry or the new one and never a torn file.
/// Readers open with <see cref="FileShare.Delete"/> so a replace or an eviction never fails because a file is
/// being read. A file that is unreadable or corrupt counts as a miss.
/// </para>
/// <para>
/// Size cap: when the total grows past <see cref="MaxBytes"/>, the least recently used entries (by file write time,
/// which a cache hit refreshes) are deleted until the total is under 90 % of the cap.
/// </para>
/// </remarks>
public sealed class DiskDataCache : IDataCache
{
    /// <summary>Environment variable that overrides <see cref="DefaultDirectory"/>.</summary>
    public const string DirectoryEnvironmentVariable = "BIOTECH_NET_CACHE";

    private const string Extension = ".entry";
    private static ReadOnlySpan<byte> Magic => "BND1"u8;
    private readonly SemaphoreSlim _trimGate = new(1, 1);
    private long _approxBytes = -1;

    /// <param name="directory">Cache root; <c>null</c> means <see cref="DefaultDirectory"/>.</param>
    /// <param name="maxBytes">Size cap in bytes (default 2 GB).</param>
    public DiskDataCache(string? directory = null, long maxBytes = 2L * 1024 * 1024 * 1024)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(maxBytes);
        Directory = Path.GetFullPath(directory ?? DefaultDirectory);
        MaxBytes = maxBytes;
    }

    /// <summary>
    /// <c>$BIOTECH_NET_CACHE</c> when set; otherwise <c>Biotech.Net/cache</c> under the user's local application data
    /// folder (<c>%LOCALAPPDATA%</c> on Windows, <c>~/.local/share</c> on Linux, <c>~/Library/Application Support</c> on macOS).
    /// </summary>
    public static string DefaultDirectory
    {
        get
        {
            string? env = Environment.GetEnvironmentVariable(DirectoryEnvironmentVariable);
            if (!string.IsNullOrWhiteSpace(env)) return env;
            string root = Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData, Environment.SpecialFolderOption.DoNotVerify);
            if (string.IsNullOrEmpty(root)) root = Path.GetTempPath();
            return Path.Combine(root, "Biotech.Net", "cache");
        }
    }

    /// <summary>The cache root directory.</summary>
    public string Directory { get; }

    /// <summary>The size cap in bytes.</summary>
    public long MaxBytes { get; }

    /// <summary>The file that holds (or would hold) <paramref name="key"/>.</summary>
    public string PathFor(string key)
    {
        string name = IsSafeKey(key) ? key : Convert.ToHexStringLower(SHA256.HashData(Encoding.UTF8.GetBytes(key)));
        return Path.Combine(Directory, name[..2], name + Extension);
    }

    public async ValueTask<CacheEntry?> GetAsync(string key, CancellationToken cancellationToken = default)
    {
        string path = PathFor(key);
        byte[] bytes;
        try
        {
            await using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete, 81920, useAsync: true);
            bytes = new byte[fs.Length];
            await fs.ReadExactlyAsync(bytes, cancellationToken).ConfigureAwait(false);
        }
        catch (Exception ex) when (ex is FileNotFoundException or DirectoryNotFoundException)
        {
            return null;
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
        {
            return null; // locked or vanished mid-read: a miss, not an error
        }

        var entry = Decode(bytes);
        if (entry is null)
        {
            TryDelete(path); // corrupt or foreign file
            return null;
        }

        try
        {
            File.SetLastWriteTimeUtc(path, DateTime.UtcNow); // LRU: a hit makes the entry young again
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
        {
            // Best effort only.
        }

        return entry;
    }

    public async ValueTask SetAsync(string key, CacheEntry entry, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(entry);
        string path = PathFor(key);
        string dir = Path.GetDirectoryName(path)!;
        string temp = Path.Combine(dir, $"{Path.GetFileName(path)}.{Guid.NewGuid():N}.tmp");
        byte[] encoded = Encode(entry);
        try
        {
            System.IO.Directory.CreateDirectory(dir);
            await using (var fs = new FileStream(temp, FileMode.CreateNew, FileAccess.Write, FileShare.None, 81920, useAsync: true))
            {
                await fs.WriteAsync(encoded, cancellationToken).ConfigureAwait(false);
                await fs.FlushAsync(cancellationToken).ConfigureAwait(false);
            }

            long previous = TryLength(path);
            if (!await MoveWithRetryAsync(temp, path, cancellationToken).ConfigureAwait(false))
            {
                TryDelete(temp);
                return;
            }

            long total = Interlocked.Read(ref _approxBytes);
            if (total < 0) total = Interlocked.Exchange(ref _approxBytes, ScanSize());
            else Interlocked.Add(ref _approxBytes, encoded.LongLength - Math.Max(0, previous));
            if (Interlocked.Read(ref _approxBytes) > MaxBytes) await TrimAsync(cancellationToken).ConfigureAwait(false);
        }
        catch (OperationCanceledException)
        {
            TryDelete(temp);
            throw;
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
        {
            TryDelete(temp); // a cache write failure never fails the request
        }
    }

    public ValueTask<bool> RemoveAsync(string key, CancellationToken cancellationToken = default)
    {
        string path = PathFor(key);
        long len = TryLength(path);
        bool removed = len >= 0 && TryDelete(path);
        if (removed && Interlocked.Read(ref _approxBytes) >= 0) Interlocked.Add(ref _approxBytes, -len);
        return ValueTask.FromResult(removed);
    }

    public ValueTask ClearAsync(CancellationToken cancellationToken = default)
    {
        foreach (var f in EnumerateFiles("*" + Extension)) TryDelete(f.FullName);
        foreach (var f in EnumerateFiles("*.tmp")) TryDelete(f.FullName);
        Interlocked.Exchange(ref _approxBytes, 0);
        return ValueTask.CompletedTask;
    }

    /// <summary>Number of entries on disk.</summary>
    public int Count => EnumerateFiles("*" + Extension).Count();

    /// <summary>Total bytes of the entries on disk.</summary>
    public long SizeBytes => ScanSize();

    /// <summary>Deletes least recently used entries until the cache is under 90 % of <see cref="MaxBytes"/>.</summary>
    public async Task TrimAsync(CancellationToken cancellationToken = default)
    {
        await _trimGate.WaitAsync(cancellationToken).ConfigureAwait(false);
        try
        {
            // Temporary files left behind by a crashed writer.
            foreach (var t in EnumerateFiles("*.tmp"))
                if (DateTime.UtcNow - t.LastWriteTimeUtc > TimeSpan.FromHours(1)) TryDelete(t.FullName);

            var files = EnumerateFiles("*" + Extension).ToList();
            long total = files.Sum(f => f.Length);
            long target = MaxBytes / 10 * 9;
            if (total > MaxBytes)
            {
                foreach (var f in files.OrderBy(f => f.LastWriteTimeUtc))
                {
                    if (total <= target) break;
                    if (TryDelete(f.FullName)) total -= f.Length;
                }
            }

            Interlocked.Exchange(ref _approxBytes, total);
        }
        finally
        {
            _trimGate.Release();
        }
    }

    private static bool IsSafeKey(string key)
    {
        if (key.Length is < 16 or > 128) return false;
        foreach (char c in key)
            if (c is not (>= '0' and <= '9' or >= 'a' and <= 'f')) return false;
        return true;
    }

    private IEnumerable<FileInfo> EnumerateFiles(string pattern)
    {
        var root = new DirectoryInfo(Directory);
        if (!root.Exists) return [];
        try
        {
            return root.EnumerateFiles(pattern, SearchOption.AllDirectories).ToList();
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
        {
            return [];
        }
    }

    private long ScanSize() => EnumerateFiles("*" + Extension).Sum(f =>
    {
        try { return f.Length; }
        catch (IOException) { return 0L; }
    });

    private static async Task<bool> MoveWithRetryAsync(string from, string to, CancellationToken ct)
    {
        for (int attempt = 0; ; attempt++)
        {
            try
            {
                File.Move(from, to, overwrite: true);
                return true;
            }
            catch (Exception ex) when ((ex is IOException or UnauthorizedAccessException) && attempt < 5)
            {
                // Windows: another process may hold the target without FILE_SHARE_DELETE for a moment.
                await Task.Delay(10 * (attempt + 1), ct).ConfigureAwait(false);
            }
            catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
            {
                return false;
            }
        }
    }

    private static long TryLength(string path)
    {
        try { return new FileInfo(path) is { Exists: true } fi ? fi.Length : -1; }
        catch (IOException) { return -1; }
    }

    private static bool TryDelete(string path)
    {
        try
        {
            if (!File.Exists(path)) return false;
            File.Delete(path);
            return true;
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException)
        {
            return false;
        }
    }

    internal static byte[] Encode(CacheEntry entry)
    {
        using var meta = new MemoryStream();
        using (var w = new Utf8JsonWriter(meta))
        {
            w.WriteStartObject();
            w.WriteString("url", entry.Url);
            w.WriteString("contentType", entry.ContentType);
            w.WriteString("storedAt", entry.StoredAt);
            w.WriteStartObject("headers");
            foreach (var (k, v) in entry.Headers) w.WriteString(k, v);
            w.WriteEndObject();
            w.WriteEndObject();
        }

        var result = new byte[8 + meta.Length + entry.Content.Length];
        Magic.CopyTo(result);
        BinaryPrimitives.WriteInt32LittleEndian(result.AsSpan(4), (int)meta.Length);
        meta.GetBuffer().AsSpan(0, (int)meta.Length).CopyTo(result.AsSpan(8));
        entry.Content.CopyTo(result.AsSpan(8 + (int)meta.Length));
        return result;
    }

    internal static CacheEntry? Decode(byte[] bytes)
    {
        if (bytes.Length < 8 || !bytes.AsSpan(0, 4).SequenceEqual(Magic)) return null;
        int metaLength = BinaryPrimitives.ReadInt32LittleEndian(bytes.AsSpan(4));
        if (metaLength <= 0 || metaLength > bytes.Length - 8) return null;
        try
        {
            using var doc = JsonDocument.Parse(bytes.AsMemory(8, metaLength));
            var root = doc.RootElement;
            var headers = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (root.TryGetProperty("headers", out var h) && h.ValueKind == JsonValueKind.Object)
                foreach (var p in h.EnumerateObject()) headers[p.Name] = p.Value.GetString() ?? "";
            return new CacheEntry(bytes.AsSpan(8 + metaLength).ToArray(), root.GetProperty("storedAt").GetDateTimeOffset())
            {
                Url = root.TryGetProperty("url", out var u) ? u.GetString() : null,
                ContentType = root.TryGetProperty("contentType", out var ct) ? ct.GetString() : null,
                Headers = headers,
            };
        }
        catch (Exception ex) when (ex is JsonException or FormatException or InvalidOperationException or KeyNotFoundException)
        {
            return null;
        }
    }
}
