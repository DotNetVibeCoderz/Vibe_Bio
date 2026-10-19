using System.Net;
using System.Net.Http.Headers;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.Tests;

public class HttpEngineTests
{
    private static readonly Uri Url = new("https://data.example.org/item/1?b=2&a=1");

    private static DataRequest Get(TimeSpan? ttl = null) => DataRequest.Get(Url) with { CacheTtl = ttl ?? TimeSpan.FromHours(1), RequestsPerSecond = 1000 };

    private static (FakeHttpHandler Handler, FakeClock Clock, DataHttp Http, MemoryDataCache Cache) Setup(
        DataCacheMode mode = DataCacheMode.Default, Func<DataClientOptions, DataClientOptions>? tweak = null)
    {
        var handler = new FakeHttpHandler();
        var clock = new FakeClock();
        var cache = new MemoryDataCache();
        var options = TestOptions.Create(clock, cache, mode);
        if (tweak is not null) options = tweak(options);
        return (handler, clock, new DataHttp("Test", handler.Client(), options), cache);
    }

    [Fact]
    public async Task CacheMissThenHit_FetchesOnce()
    {
        var (handler, _, http, cache) = Setup();
        handler.On(_ => true, () => Fixtures.Json("{\"v\":1}"));
        var first = await http.SendAsync(Get(), default);
        var second = await http.SendAsync(Get(), default);
        Assert.False(first.FromCache);
        Assert.True(second.FromCache);
        Assert.Equal("{\"v\":1}", second.Text());
        Assert.Equal(1, handler.Count);
        Assert.Equal(1, cache.Count);
    }

    [Fact]
    public async Task ExpiredEntry_IsRefetched()
    {
        var (handler, clock, http, _) = Setup();
        int n = 0;
        handler.On(_ => true, () => Fixtures.Json($"{{\"v\":{++n}}}"));
        await http.SendAsync(Get(TimeSpan.FromHours(1)), default);
        clock.Advance(TimeSpan.FromMinutes(59));
        Assert.Equal("{\"v\":1}", (await http.SendAsync(Get(TimeSpan.FromHours(1)), default)).Text());
        clock.Advance(TimeSpan.FromMinutes(2));
        var refreshed = await http.SendAsync(Get(TimeSpan.FromHours(1)), default);
        Assert.False(refreshed.FromCache);
        Assert.Equal("{\"v\":2}", refreshed.Text());
        Assert.Equal(2, handler.Count);
    }

    [Fact]
    public async Task OptionsTtl_OverridesTheEndpointTtl()
    {
        var (handler, clock, http, _) = Setup(tweak: o => o with { CacheTtl = TimeSpan.FromSeconds(10) });
        handler.On(_ => true, () => Fixtures.Json("{}"));
        await http.SendAsync(Get(TimeSpan.FromDays(30)), default);
        clock.Advance(TimeSpan.FromSeconds(11));
        await http.SendAsync(Get(TimeSpan.FromDays(30)), default);
        Assert.Equal(2, handler.Count);
    }

    [Fact]
    public async Task RefreshMode_AlwaysFetches_AndOverwrites()
    {
        var (handler, clock, _, cache) = Setup();
        int n = 0;
        handler.On(_ => true, () => Fixtures.Json($"{{\"v\":{++n}}}"));
        var refresh = new DataHttp("Test", handler.Client(), TestOptions.Create(clock, cache, DataCacheMode.Refresh));
        var normal = new DataHttp("Test", handler.Client(), TestOptions.Create(clock, cache));
        await refresh.SendAsync(Get(), default);
        await refresh.SendAsync(Get(), default);
        Assert.Equal(2, handler.Count);
        Assert.Equal("{\"v\":2}", (await normal.SendAsync(Get(), default)).Text());
        Assert.Equal(2, handler.Count);
    }

    [Fact]
    public async Task BypassMode_AndUncacheableRequests_NeverTouchTheCache()
    {
        var (handler, _, http, cache) = Setup(DataCacheMode.Bypass);
        handler.On(_ => true, () => Fixtures.Json("{}"));
        await http.SendAsync(Get(), default);
        await http.SendAsync(Get(), default);
        Assert.Equal(2, handler.Count);
        Assert.Equal(0, cache.Count);

        var (handler2, _, http2, cache2) = Setup();
        handler2.On(_ => true, () => Fixtures.Json("{}"));
        await http2.SendAsync(DataRequest.Get(Url) with { CacheTtl = null }, default);
        Assert.Equal(0, cache2.Count);
    }

    [Fact]
    public async Task OfflineMode_ServesStaleEntries_WithoutTheNetwork()
    {
        var (handler, clock, http, cache) = Setup();
        handler.On(_ => true, () => Fixtures.Json("{\"v\":1}"));
        await http.SendAsync(Get(), default);
        clock.Advance(TimeSpan.FromDays(365)); // far past the TTL
        var offline = new DataHttp("Test", handler.Client(), TestOptions.Create(clock, cache, DataCacheMode.Offline));
        var r = await offline.SendAsync(Get(), default);
        Assert.True(r.FromCache);
        Assert.Equal("{\"v\":1}", r.Text());
        Assert.Equal(1, handler.Count);
    }

    [Fact]
    public async Task OfflineMode_MissingEntry_ThrowsAClearException()
    {
        var (handler, _, http, _) = Setup(DataCacheMode.Offline);
        var ex = await Assert.ThrowsAsync<DataOfflineException>(() => http.SendAsync(Get(), default));
        Assert.Contains("no cached response", ex.Message);
        Assert.Contains("data.example.org/item/1", ex.Message);
        Assert.Contains(DataClientOptions.OfflineEnvironmentVariable, ex.Message);
        var ex2 = await Assert.ThrowsAsync<DataOfflineException>(() => http.SendAsync(DataRequest.Get(Url), default));
        Assert.Contains("never cached", ex2.Message);
        Assert.Equal(0, handler.Count);
    }

    [Fact]
    public async Task TransientFailure_IsRetried_WithJitteredExponentialBackoff()
    {
        var (handler, clock, http, _) = Setup();
        handler.OnSequence(_ => true,
            () => Fixtures.Status(HttpStatusCode.ServiceUnavailable),
            () => Fixtures.Status(HttpStatusCode.BadGateway),
            () => Fixtures.Json("{\"ok\":true}"));
        var r = await http.SendAsync(Get(), default);
        Assert.Equal("{\"ok\":true}", r.Text());
        Assert.Equal(3, handler.Count);
        var delays = clock.Delays.ToArray();
        Assert.Equal(2, delays.Length);
        Assert.InRange(delays[0].TotalMilliseconds, 500, 1000);   // initial 1 s, equal jitter
        Assert.InRange(delays[1].TotalMilliseconds, 1000, 2000);  // doubled
    }

    [Fact]
    public async Task TooManyRequests_HonoursRetryAfter_AndPausesTheHost()
    {
        var limiter = new HostRateLimiter();
        var (handler, clock, http, _) = Setup(tweak: o => o with { RateLimiter = limiter });
        handler.OnSequence(_ => true,
            () =>
            {
                var r = Fixtures.Status(HttpStatusCode.TooManyRequests);
                r.Headers.RetryAfter = new RetryConditionHeaderValue(TimeSpan.FromSeconds(7));
                return r;
            },
            () => Fixtures.Json("{}"));
        var start = clock.GetUtcNow();
        await http.SendAsync(Get(), default);
        Assert.Equal(TimeSpan.FromSeconds(7), clock.Delays.First());
        Assert.Equal(start + TimeSpan.FromSeconds(7), limiter.BlockedUntil("data.example.org"));
    }

    [Fact]
    public async Task RetryAfterLongerThanTheLimit_FailsFast()
    {
        var (handler, _, http, _) = Setup(tweak: o => o with { MaxRetryAfter = TimeSpan.FromSeconds(30) });
        handler.On(_ => true, () =>
        {
            var r = Fixtures.Status(HttpStatusCode.ServiceUnavailable);
            r.Headers.RetryAfter = new RetryConditionHeaderValue(TimeSpan.FromHours(1));
            return r;
        });
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => http.SendAsync(Get(), default));
        Assert.Contains("MaxRetryAfter", ex.Message);
        Assert.Equal(1, handler.Count);
    }

    [Fact]
    public async Task ExhaustedRetries_ThrowWithStatusAndAttemptCount()
    {
        var (handler, clock, http, _) = Setup(tweak: o => o with { MaxRetries = 3 });
        handler.On(_ => true, () => Fixtures.Json("{\"error\":\"database is down\"}", HttpStatusCode.InternalServerError));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => http.SendAsync(Get(), default));
        Assert.Equal(HttpStatusCode.InternalServerError, ex.StatusCode);
        Assert.Equal("database is down", ex.ServiceMessage);
        Assert.Contains("after 4 attempts", ex.Message);
        Assert.Equal(4, handler.Count);
        Assert.Equal(3, clock.Delays.Count);
    }

    [Fact]
    public async Task ClientErrors_AreNotRetried_AndCarryTheServiceMessage()
    {
        var (handler, _, http, _) = Setup();
        handler.On(_ => true, () => Fixtures.Json("{\"messages\":[\"The 'accession' value has invalid format.\"]}", HttpStatusCode.BadRequest));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => http.SendAsync(Get(), default));
        Assert.Equal(1, handler.Count);
        Assert.Equal(HttpStatusCode.BadRequest, ex.StatusCode);
        Assert.Contains("invalid format", ex.Message);
        Assert.False(ex.IsNotFound);
    }

    [Fact]
    public async Task HtmlErrorPages_AreSummarisedByTitle()
    {
        var (handler, _, http, _) = Setup(tweak: o => o with { MaxRetries = 0 });
        handler.On(_ => true, () => new HttpResponseMessage(HttpStatusCode.BadGateway)
        {
            Content = new StringContent("<!doctype html><html><head><title>Error: 502 | EBI</title></head><body>...</body></html>", System.Text.Encoding.UTF8, "text/html"),
        });
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => http.SendAsync(Get(), default));
        Assert.Equal("HTML error page: Error: 502 | EBI", ex.ServiceMessage);
    }

    [Fact]
    public async Task Timeout_IsRetried_ThenReportedAsATimeout()
    {
        var (handler, _, http, _) = Setup(tweak: o => o with { Timeout = TimeSpan.FromSeconds(5), MaxRetries = 1 });
        handler.On(_ => true, async (_, ct) =>
        {
            await Task.Delay(Timeout.Infinite, ct); // never answers
            return Fixtures.Json("{}");
        });
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => http.SendAsync(Get(), default));
        Assert.IsType<TimeoutException>(ex.InnerException);
        Assert.Contains("did not answer within 5 s (2 attempts)", ex.Message);
        Assert.Equal(2, handler.Count);
    }

    [Fact]
    public async Task NetworkErrors_AreRetried()
    {
        var (handler, _, http, _) = Setup();
        int n = 0;
        handler.On(_ => true, (_, _) => ++n == 1 ? throw new HttpRequestException("connection reset") : Task.FromResult(Fixtures.Json("{}")));
        await http.SendAsync(Get(), default);
        Assert.Equal(2, handler.Count);
    }

    [Fact]
    public async Task CallerCancellation_IsNotWrappedOrRetried()
    {
        var (handler, _, http, _) = Setup();
        using var cts = new CancellationTokenSource();
        handler.On(_ => true, async (_, ct) =>
        {
            cts.Cancel();
            await Task.Delay(Timeout.Infinite, ct);
            return Fixtures.Json("{}");
        });
        await Assert.ThrowsAnyAsync<OperationCanceledException>(() => http.SendAsync(Get(), cts.Token));
        Assert.Equal(1, handler.Count);
    }

    [Fact]
    public async Task ExhaustedQuotaHeaders_PauseTheHostUntilReset()
    {
        var limiter = new HostRateLimiter();
        var (handler, clock, http, _) = Setup(tweak: o => o with { RateLimiter = limiter });
        handler.On(_ => true, () =>
        {
            var r = Fixtures.Json("{}");
            r.Headers.Add("X-RateLimit-Limit", "55000");
            r.Headers.Add("X-RateLimit-Remaining", "0");
            r.Headers.Add("X-RateLimit-Reset", "3");
            return r;
        });
        var t0 = clock.GetUtcNow();
        await http.SendAsync(Get() with { CacheTtl = null }, default);
        Assert.Equal(t0 + TimeSpan.FromSeconds(3), limiter.BlockedUntil("data.example.org"));
        await http.SendAsync(Get() with { CacheTtl = null }, default);
        Assert.Contains(TimeSpan.FromSeconds(3), clock.Delays);
    }

    [Fact]
    public async Task SendsUserAgent_AndBodies()
    {
        var (handler, _, http, _) = Setup();
        handler.On(_ => true, () => Fixtures.Json("{}"));
        await http.SendAsync(DataRequest.PostJson(new Uri("https://data.example.org/q"), "{\"a\":1}") with { CacheTtl = TimeSpan.FromHours(1) }, default);
        var seen = handler.Requests.Single();
        Assert.StartsWith("Biotech.Net/", seen.UserAgent);
        Assert.Equal("{\"a\":1}", seen.Body);
        Assert.Equal("application/json", seen.ContentType);
    }

    [Fact]
    public async Task Exceptions_NeverCarryCredentials()
    {
        var (handler, _, http, _) = Setup(tweak: o => o with { MaxRetries = 0 });
        handler.On(_ => true, () => Fixtures.Status(HttpStatusCode.Forbidden, "denied"));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() =>
            http.SendAsync(DataRequest.Get(new Uri("https://data.example.org/x?id=1&api_key=TOPSECRET&email=me%40x.org")), default));
        Assert.DoesNotContain("TOPSECRET", ex.Message);
        Assert.DoesNotContain("TOPSECRET", ex.RequestUri!.ToString());
        Assert.DoesNotContain("me%40x.org", ex.Message);
        Assert.Contains("id=1", ex.Message);
    }

    [Fact]
    public void Backoff_GrowsExponentially_UpToTheCap()
    {
        var (_, _, http, _) = Setup(tweak: o => o with { InitialBackoff = TimeSpan.FromSeconds(1), MaxBackoff = TimeSpan.FromSeconds(8) });
        for (int attempt = 0; attempt < 8; attempt++)
        {
            double d = Math.Min(Math.Pow(2, attempt), 8) * 1000;
            Assert.InRange(http.Backoff(attempt).TotalMilliseconds, d / 2, d);
        }
    }

    [Fact]
    public async Task DiskCache_PersistsAcrossClientInstances()
    {
        using var dir = new TempDir();
        var clock = new FakeClock();
        var handler = new FakeHttpHandler();
        handler.On(_ => true, () => Fixtures.Json("{\"persisted\":true}"));
        var options = TestOptions.Create(clock, new DiskDataCache(dir.Path));
        await new DataHttp("Test", handler.Client(), options).SendAsync(Get(), default);
        var offline = TestOptions.Create(clock, new DiskDataCache(dir.Path), DataCacheMode.Offline);
        var r = await new DataHttp("Test", handler.Client(), offline).SendAsync(Get(), default);
        Assert.True(r.FromCache);
        Assert.Equal("{\"persisted\":true}", r.Text());
        Assert.Equal(1, handler.Count);
    }
}

public class HostRateLimiterTests
{
    [Fact]
    public async Task SequentialRequests_AreSpacedByTheRate()
    {
        var clock = new FakeClock();
        var limiter = new HostRateLimiter();
        var waits = new List<TimeSpan>();
        for (int i = 0; i < 4; i++) waits.Add(await limiter.WaitAsync("eutils.ncbi.nlm.nih.gov", 3, clock));
        Assert.Equal(TimeSpan.Zero, waits[0]);
        // Request k starts at k/3 s: never early, at most 1 ms late (timers have millisecond resolution).
        double elapsed = 0;
        for (int k = 1; k < waits.Count; k++)
        {
            elapsed += waits[k].TotalMilliseconds;
            Assert.InRange(elapsed, k * 1000.0 / 3, k * 1000.0 / 3 + 1);
        }

        Assert.InRange((clock.GetUtcNow() - new DateTimeOffset(2026, 1, 1, 0, 0, 0, TimeSpan.Zero)).TotalMilliseconds, 1000, 1001);
    }

    [Fact]
    public async Task ConcurrentRequests_GetDistinctSlots()
    {
        var clock = new FakeClock();
        var limiter = new HostRateLimiter();
        clock.Pause();
        var tasks = Enumerable.Range(0, 10).Select(_ => limiter.WaitAsync("rest.ensembl.org", 10, clock).AsTask()).ToList();
        var reserved = clock.Delays.Select(d => d.TotalMilliseconds).Order().ToArray();
        Assert.Equal(Enumerable.Range(1, 9).Select(i => i * 100.0).ToArray(), reserved, (a, b) => Math.Abs(a - b) < 1e-9);
        clock.Resume();
        await Task.WhenAll(tasks);
    }

    [Fact]
    public async Task Hosts_AreIndependent_AndBlocksDelayOnlyTheirHost()
    {
        var clock = new FakeClock();
        var limiter = new HostRateLimiter();
        await limiter.WaitAsync("a.org", 1, clock);
        Assert.Equal(TimeSpan.Zero, await limiter.WaitAsync("b.org", 1, clock));
        limiter.BlockUntil("b.org", clock.GetUtcNow() + TimeSpan.FromSeconds(30));
        Assert.Equal(TimeSpan.FromSeconds(30), await limiter.WaitAsync("b.org", 1, clock));
        Assert.Null(limiter.BlockedUntil("c.org"));
    }
}
