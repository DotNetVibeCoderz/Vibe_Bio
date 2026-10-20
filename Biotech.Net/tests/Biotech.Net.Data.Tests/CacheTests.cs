using System.Text;
using Biotech.Net.Data.Caching;

namespace Biotech.Net.Data.Tests;

public class CacheKeyTests
{
    [Fact]
    public void NormalizeUrl_SortsParameters_LowercasesHost_AndDropsCredentials()
    {
        var uri = new Uri("HTTPS://EUtils.NCBI.nlm.nih.gov:443/entrez/eutils/esearch.fcgi?term=BRCA1&db=gene&api_key=SECRET&email=a%40b.org&tool=x");
        Assert.Equal("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=BRCA1", CacheKey.NormalizeUrl(uri));
    }

    [Fact]
    public void Compute_IgnoresParameterOrderSeparatorsAndCredentials()
    {
        string a = CacheKey.Compute("GET", new Uri("https://rest.ensembl.org/lookup/id/X?expand=1&content-type=application/json"));
        string b = CacheKey.Compute("get", new Uri("https://rest.ensembl.org/lookup/id/X?content-type=application/json;expand=1;api_key=k"));
        Assert.Equal(a, b);
        Assert.Matches("^[0-9a-f]{64}$", a);
    }

    [Fact]
    public void Compute_DependsOnMethodPathAndBody()
    {
        var uri = new Uri("https://search.rcsb.org/rcsbsearch/v2/query");
        string get = CacheKey.Compute("GET", uri);
        string post1 = CacheKey.Compute("POST", uri, "{\"q\":1}"u8);
        string post2 = CacheKey.Compute("POST", uri, "{\"q\":2}"u8);
        Assert.NotEqual(get, post1);
        Assert.NotEqual(post1, post2);
        Assert.NotEqual(CacheKey.Compute("GET", new Uri("https://x.org/a")), CacheKey.Compute("GET", new Uri("https://x.org/b")));
    }
}

public class MemoryDataCacheTests
{
    private static CacheEntry Entry(int size, byte fill = 1) => new(Enumerable.Repeat(fill, size).ToArray(), DateTimeOffset.UnixEpoch);

    [Fact]
    public async Task EvictsLeastRecentlyUsed_ByCount()
    {
        var cache = new MemoryDataCache(maxEntries: 2);
        await cache.SetAsync("a", Entry(1));
        await cache.SetAsync("b", Entry(1));
        Assert.NotNull(await cache.GetAsync("a")); // a is now more recent than b
        await cache.SetAsync("c", Entry(1));
        Assert.NotNull(await cache.GetAsync("a"));
        Assert.Null(await cache.GetAsync("b"));
        Assert.NotNull(await cache.GetAsync("c"));
        Assert.Equal(2, cache.Count);
    }

    [Fact]
    public async Task EvictsByBytes_AndRejectsOversizedEntries()
    {
        var cache = new MemoryDataCache(maxEntries: 100, maxBytes: 100);
        await cache.SetAsync("a", Entry(60));
        await cache.SetAsync("b", Entry(60));
        Assert.Null(await cache.GetAsync("a"));
        Assert.Equal(60, cache.SizeBytes);
        await cache.SetAsync("huge", Entry(101));
        Assert.Null(await cache.GetAsync("huge"));
        Assert.NotNull(await cache.GetAsync("b"));
    }

    [Fact]
    public async Task ReplaceRemoveAndClear_KeepTheByteCountRight()
    {
        var cache = new MemoryDataCache();
        await cache.SetAsync("a", Entry(10));
        await cache.SetAsync("a", Entry(30, 2));
        Assert.Equal(30, cache.SizeBytes);
        Assert.Equal(2, (await cache.GetAsync("a"))!.Content[0]);
        Assert.True(await cache.RemoveAsync("a"));
        Assert.False(await cache.RemoveAsync("a"));
        await cache.SetAsync("b", Entry(5));
        await cache.ClearAsync();
        Assert.Equal(0, cache.Count);
        Assert.Equal(0, cache.SizeBytes);
    }
}

public class DiskDataCacheTests
{
    private static readonly string Key1 = CacheKey.Compute("GET", new Uri("https://example.org/1"));
    private static readonly string Key2 = CacheKey.Compute("GET", new Uri("https://example.org/2"));
    private static readonly string Key3 = CacheKey.Compute("GET", new Uri("https://example.org/3"));

    [Fact]
    public async Task RoundTrip_KeepsBodyMetadataAndHeaders()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path);
        var stored = new DateTimeOffset(2026, 5, 4, 3, 2, 1, TimeSpan.Zero);
        var entry = new CacheEntry([0x1f, 0x8b, 0, 255], stored)
        {
            Url = "https://example.org/1",
            ContentType = "application/octet-stream",
            Headers = new Dictionary<string, string> { ["Link"] = "<https://example.org/2>; rel=\"next\"" },
        };
        await cache.SetAsync(Key1, entry);
        var back = await new DiskDataCache(dir.Path).GetAsync(Key1); // a second instance, as another process would see it
        Assert.NotNull(back);
        Assert.Equal(entry.Content, back.Content);
        Assert.Equal(stored, back.StoredAt);
        Assert.Equal("application/octet-stream", back.ContentType);
        Assert.Equal("https://example.org/1", back.Url);
        Assert.Equal("<https://example.org/2>; rel=\"next\"", back.Headers["link"]);
        Assert.StartsWith(Path.Combine(dir.Path, Key1[..2]), cache.PathFor(Key1));
        Assert.Equal(1, cache.Count);
    }

    [Fact]
    public async Task Overwrite_ReplacesAtomically_AndLeavesNoTemporaryFiles()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path);
        await cache.SetAsync(Key1, new CacheEntry("old"u8.ToArray(), DateTimeOffset.UnixEpoch));
        await cache.SetAsync(Key1, new CacheEntry("new"u8.ToArray(), DateTimeOffset.UnixEpoch));
        Assert.Equal("new", Encoding.ASCII.GetString((await cache.GetAsync(Key1))!.Content));
        Assert.Empty(Directory.GetFiles(dir.Path, "*.tmp", SearchOption.AllDirectories));
        Assert.Single(Directory.GetFiles(dir.Path, "*.entry", SearchOption.AllDirectories));
    }

    [Fact]
    public async Task CorruptFile_IsAMiss_AndIsDeleted()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path);
        string path = cache.PathFor(Key1);
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        await File.WriteAllBytesAsync(path, "BND1\xff\xff\xff\x7fgarbage"u8.ToArray());
        Assert.Null(await cache.GetAsync(Key1));
        Assert.False(File.Exists(path));
        Assert.Null(await cache.GetAsync(Key2)); // never written
    }

    [Fact]
    public async Task SizeCap_EvictsLeastRecentlyUsedFirst()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path, maxBytes: 2500);
        await cache.SetAsync(Key1, new CacheEntry(new byte[1000], DateTimeOffset.UnixEpoch));
        await cache.SetAsync(Key2, new CacheEntry(new byte[1000], DateTimeOffset.UnixEpoch));
        File.SetLastWriteTimeUtc(cache.PathFor(Key1), DateTime.UtcNow.AddHours(-2));
        File.SetLastWriteTimeUtc(cache.PathFor(Key2), DateTime.UtcNow.AddHours(-1));
        await cache.SetAsync(Key3, new CacheEntry(new byte[1000], DateTimeOffset.UnixEpoch)); // > 2500 bytes: trim to <= 2250
        Assert.Null(await cache.GetAsync(Key1));
        Assert.NotNull(await cache.GetAsync(Key2));
        Assert.NotNull(await cache.GetAsync(Key3));
        Assert.True(cache.SizeBytes <= 2500);
    }

    [Fact]
    public async Task ConcurrentWritersAndReaders_NeverSeeATornEntry()
    {
        using var dir = new TempDir();
        var writerA = new DiskDataCache(dir.Path);
        var writerB = new DiskDataCache(dir.Path); // two instances = two "processes" sharing a directory
        byte[] small = Enumerable.Repeat((byte)'a', 1_000).ToArray();
        byte[] large = Enumerable.Repeat((byte)'b', 200_000).ToArray();
        int reads = 0, hits = 0;
        // Seed the entry so readers that start before the first concurrent write still hit; otherwise a fast
        // scheduler can finish every read before any writer runs and the test would prove nothing.
        await writerA.SetAsync(Key1, new CacheEntry(small, DateTimeOffset.UnixEpoch));
        var tasks = new List<Task>();
        for (int w = 0; w < 4; w++)
        {
            var cache = w % 2 == 0 ? writerA : writerB;
            byte[] payload = w % 2 == 0 ? small : large;
            tasks.Add(Task.Run(async () =>
            {
                for (int i = 0; i < 40; i++) await cache.SetAsync(Key1, new CacheEntry(payload, DateTimeOffset.UnixEpoch));
            }));
        }

        for (int r = 0; r < 4; r++)
        {
            var cache = r % 2 == 0 ? writerA : writerB;
            tasks.Add(Task.Run(async () =>
            {
                for (int i = 0; i < 80; i++)
                {
                    var e = await cache.GetAsync(Key1);
                    Interlocked.Increment(ref reads);
                    if (e is null) continue;
                    Interlocked.Increment(ref hits);
                    // Either the whole small payload or the whole large one: never a mix or a prefix.
                    Assert.True(e.Content.AsSpan().SequenceEqual(small) || e.Content.AsSpan().SequenceEqual(large),
                        $"torn read of {e.Content.Length} bytes");
                }
            }));
        }

        await Task.WhenAll(tasks);
        Assert.Equal(320, reads);
        Assert.True(hits > 0);
        Assert.Empty(Directory.GetFiles(dir.Path, "*.tmp", SearchOption.AllDirectories));
    }

    [Fact]
    public async Task RemoveAndClear()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path);
        await cache.SetAsync(Key1, new CacheEntry([1], DateTimeOffset.UnixEpoch));
        await cache.SetAsync(Key2, new CacheEntry([2], DateTimeOffset.UnixEpoch));
        Assert.True(await cache.RemoveAsync(Key1));
        Assert.False(await cache.RemoveAsync(Key1));
        await cache.ClearAsync();
        Assert.Equal(0, cache.Count);
    }

    [Fact]
    public async Task ArbitraryKeys_AreHashedToSafeFileNames()
    {
        using var dir = new TempDir();
        var cache = new DiskDataCache(dir.Path);
        await cache.SetAsync("../../not a path?*", new CacheEntry([7], DateTimeOffset.UnixEpoch));
        Assert.Equal([7], (await cache.GetAsync("../../not a path?*"))!.Content);
        Assert.StartsWith(dir.Path, cache.PathFor("../../not a path?*"));
    }

    [Fact]
    public void DefaultDirectory_HonoursTheEnvironmentVariable()
    {
        string? old = Environment.GetEnvironmentVariable(DiskDataCache.DirectoryEnvironmentVariable);
        try
        {
            Environment.SetEnvironmentVariable(DiskDataCache.DirectoryEnvironmentVariable, null);
            Assert.EndsWith(Path.Combine("Biotech.Net", "cache"), DiskDataCache.DefaultDirectory);
            Environment.SetEnvironmentVariable(DiskDataCache.DirectoryEnvironmentVariable, Path.Combine(Path.GetTempPath(), "bn-cache-x"));
            Assert.Equal(Path.Combine(Path.GetTempPath(), "bn-cache-x"), DiskDataCache.DefaultDirectory);
        }
        finally
        {
            Environment.SetEnvironmentVariable(DiskDataCache.DirectoryEnvironmentVariable, old);
        }
    }
}
