using System.Globalization;
using System.IO.Compression;
using System.Net;
using System.Net.Http.Headers;
using System.Text;
using System.Text.Json;
using Biotech.Net.Data.Caching;

namespace Biotech.Net.Data.Http;

/// <summary>One logical request to a data service.</summary>
internal sealed record DataRequest(HttpMethod Method, Uri Uri)
{
    /// <summary>Request body (POST), or null.</summary>
    public byte[]? Body { get; init; }

    /// <summary>Content type of <see cref="Body"/>.</summary>
    public string? BodyContentType { get; init; }

    /// <summary>The <c>Accept</c> header, or null.</summary>
    public string? Accept { get; init; }

    /// <summary>Requests per second allowed for the host.</summary>
    public double RequestsPerSecond { get; init; } = 10;

    /// <summary>How long a cached response stays fresh; null means "never cache" (job submission, polling).</summary>
    public TimeSpan? CacheTtl { get; init; }

    public static DataRequest Get(Uri uri) => new(HttpMethod.Get, uri);

    public static DataRequest PostJson(Uri uri, string json) =>
        new(HttpMethod.Post, uri) { Body = Encoding.UTF8.GetBytes(json), BodyContentType = "application/json", Accept = "application/json" };

    public static DataRequest PostForm(Uri uri, string form) =>
        new(HttpMethod.Post, uri) { Body = Encoding.UTF8.GetBytes(form), BodyContentType = "application/x-www-form-urlencoded" };
}

/// <summary>A successful (2xx) response, from the network or the cache.</summary>
internal sealed class DataResponse
{
    public required Uri RequestUri { get; init; }
    public required HttpStatusCode StatusCode { get; init; }
    public required byte[] Content { get; init; }
    public string? ContentType { get; init; }
    public IReadOnlyDictionary<string, string> Headers { get; init; } = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
    public bool FromCache { get; init; }

    public string? Header(string name) => Headers.TryGetValue(name, out var v) ? v : null;

    /// <summary>The body as text, transparently un-gzipping <c>.gz</c> downloads.</summary>
    public string Text()
    {
        var bytes = Content;
        if (bytes.Length > 2 && bytes[0] == 0x1f && bytes[1] == 0x8b)
        {
            using var gz = new GZipStream(new MemoryStream(bytes), CompressionMode.Decompress);
            using var ms = new MemoryStream();
            gz.CopyTo(ms);
            bytes = ms.ToArray();
        }

        return Encoding.UTF8.GetString(bytes);
    }

    public JsonDocument Json() => JsonDocument.Parse(Content);
}

/// <summary>
/// The engine every client shares: cache lookup, per-host rate limiting, per-attempt timeout, retries with
/// exponential backoff and jitter honouring <c>Retry-After</c>, and conversion of failures to
/// <see cref="DataServiceException"/>.
/// </summary>
internal sealed class DataHttp
{
    /// <summary>Headers kept with cached entries (the rest are dropped).</summary>
    internal static readonly string[] KeptHeaders = ["Link", "X-Total-Results", "X-UniProt-Release", "ETag", "Last-Modified"];

    private static readonly Lazy<HttpClient> s_sharedClient = new(() => new HttpClient(new SocketsHttpHandler
    {
        AutomaticDecompression = DecompressionMethods.All,
        PooledConnectionLifetime = TimeSpan.FromMinutes(5),
    })
    {
        Timeout = System.Threading.Timeout.InfiniteTimeSpan, // per-attempt timeouts are ours
    });

    private readonly HttpClient _http;

    public DataHttp(string service, HttpClient? httpClient, DataClientOptions? options)
    {
        Service = service;
        _http = httpClient ?? s_sharedClient.Value;
        Options = options ?? DataClientOptions.Default;
    }

    public string Service { get; }
    public DataClientOptions Options { get; }
    public TimeProvider Time => Options.TimeProvider;

    public async Task<DataResponse> SendAsync(DataRequest request, CancellationToken cancellationToken)
    {
        var mode = Options.CacheMode;
        var cache = mode == DataCacheMode.Bypass || request.CacheTtl is null ? null : Options.ResolveCache();
        string? key = cache is null ? null : CacheKey.Compute(request.Method.Method, request.Uri, request.Body);

        if (mode == DataCacheMode.Offline)
        {
            if (cache is null || key is null)
                throw new DataOfflineException(Service, request.Uri, "this request is never cached (it starts or polls a server-side job)");
            var stale = await cache.GetAsync(key, cancellationToken).ConfigureAwait(false);
            return stale is not null ? FromEntry(request, stale) : throw new DataOfflineException(Service, request.Uri, "no cached response");
        }

        if (cache is not null && key is not null && mode == DataCacheMode.Default)
        {
            var entry = await cache.GetAsync(key, cancellationToken).ConfigureAwait(false);
            var ttl = Options.CacheTtl ?? request.CacheTtl!.Value;
            if (entry is not null && entry.Age(Time.GetUtcNow()) <= ttl) return FromEntry(request, entry);
        }

        var response = await SendWithRetriesAsync(request, cancellationToken).ConfigureAwait(false);
        if (cache is not null && key is not null)
        {
            var entry = new CacheEntry(response.Content, Time.GetUtcNow())
            {
                Url = DataServiceException.Redact(request.Uri).ToString(),
                ContentType = response.ContentType,
                Headers = response.Headers,
            };
            await cache.SetAsync(key, entry, cancellationToken).ConfigureAwait(false);
        }

        return response;
    }

    private static DataResponse FromEntry(DataRequest request, CacheEntry entry) => new()
    {
        RequestUri = request.Uri,
        StatusCode = HttpStatusCode.OK,
        Content = entry.Content,
        ContentType = entry.ContentType,
        Headers = entry.Headers,
        FromCache = true,
    };

    private async Task<DataResponse> SendWithRetriesAsync(DataRequest request, CancellationToken ct)
    {
        var limiter = Options.ResolveRateLimiter();
        string host = request.Uri.Host;
        double rate = Options.RequestsPerSecond ?? request.RequestsPerSecond;

        for (int attempt = 0; ; attempt++)
        {
            bool canRetry = attempt < Options.MaxRetries;
            await limiter.WaitAsync(host, rate, Time, ct).ConfigureAwait(false);

            TimeSpan? serverDelay = null;
            using (var timeout = new CancellationTokenSource(Options.Timeout, Time))
            using (var linked = CancellationTokenSource.CreateLinkedTokenSource(ct, timeout.Token))
            using (var message = BuildMessage(request))
            {
                try
                {
                    using var response = await _http.SendAsync(message, HttpCompletionOption.ResponseContentRead, linked.Token).ConfigureAwait(false);
                    byte[] body = await response.Content.ReadAsByteArrayAsync(linked.Token).ConfigureAwait(false);
                    ApplyRateLimitHeaders(limiter, host, response);
                    if (response.IsSuccessStatusCode) return ToResponse(request, response, body);

                    string? contentType = response.Content.Headers.ContentType?.ToString();
                    if (!canRetry || !IsTransient(response.StatusCode))
                        throw DataServiceException.FromResponse(Service, request.Uri, response.StatusCode, body, contentType, attempt + 1);

                    serverDelay = RetryAfter(response);
                    if (serverDelay > Options.MaxRetryAfter)
                        throw new DataServiceException(Service,
                            $"{Service} asked to retry after {serverDelay.Value.TotalSeconds:F0} s (HTTP {(int)response.StatusCode}) for {DataServiceException.Redact(request.Uri)}, longer than MaxRetryAfter.",
                            DataServiceException.Redact(request.Uri), response.StatusCode, DataServiceException.ExtractMessage(body, contentType));
                    if (response.StatusCode == HttpStatusCode.TooManyRequests)
                        limiter.BlockUntil(host, Time.GetUtcNow() + (serverDelay ?? Backoff(attempt)));
                }
                catch (OperationCanceledException ex) when (!ct.IsCancellationRequested)
                {
                    // Our per-attempt timeout (or the HttpClient's own Timeout) fired.
                    if (!canRetry)
                        throw new DataServiceException(Service,
                            $"{Service} did not answer within {Options.Timeout.TotalSeconds:F0} s ({attempt + 1} attempts) for {DataServiceException.Redact(request.Uri)}.",
                            DataServiceException.Redact(request.Uri), innerException: new TimeoutException(ex.Message, ex));
                }
                catch (HttpRequestException ex)
                {
                    if (!canRetry)
                        throw new DataServiceException(Service,
                            $"Could not reach {Service} ({attempt + 1} attempts) for {DataServiceException.Redact(request.Uri)}: {ex.Message}",
                            DataServiceException.Redact(request.Uri), ex.StatusCode, innerException: ex);
                }
            }

            var delay = serverDelay ?? Backoff(attempt);
            if (delay > TimeSpan.Zero) await Task.Delay(delay, Time, ct).ConfigureAwait(false);
        }
    }

    private HttpRequestMessage BuildMessage(DataRequest request)
    {
        var message = new HttpRequestMessage(request.Method, request.Uri);
        message.Headers.UserAgent.ParseAdd(Options.UserAgent);
        if (request.Accept is not null) message.Headers.Accept.ParseAdd(request.Accept);
        if (request.Body is not null)
        {
            message.Content = new ByteArrayContent(request.Body);
            if (request.BodyContentType is not null) message.Content.Headers.ContentType = MediaTypeHeaderValue.Parse(request.BodyContentType);
        }

        return message;
    }

    private static DataResponse ToResponse(DataRequest request, HttpResponseMessage response, byte[] body)
    {
        var headers = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        foreach (string name in KeptHeaders)
        {
            if (response.Headers.TryGetValues(name, out var values) || response.Content.Headers.TryGetValues(name, out values))
                headers[name] = string.Join(", ", values);
        }

        return new DataResponse
        {
            RequestUri = request.Uri,
            StatusCode = response.StatusCode,
            Content = body,
            ContentType = response.Content.Headers.ContentType?.ToString(),
            Headers = headers,
        };
    }

    internal static bool IsTransient(HttpStatusCode status) => status is HttpStatusCode.RequestTimeout or HttpStatusCode.TooManyRequests
        or HttpStatusCode.InternalServerError or HttpStatusCode.BadGateway or HttpStatusCode.ServiceUnavailable or HttpStatusCode.GatewayTimeout;

    /// <summary>Exponential backoff with "equal jitter": a uniform delay in [d/2, d], d = min(max, initial * 2^attempt).</summary>
    internal TimeSpan Backoff(int attempt)
    {
        double baseMs = Options.InitialBackoff.TotalMilliseconds * Math.Pow(2, Math.Min(attempt, 30));
        baseMs = Math.Min(baseMs, Options.MaxBackoff.TotalMilliseconds);
        double jitter = (Options.Random ?? Random.Shared).NextDouble();
        return TimeSpan.FromMilliseconds(baseMs / 2 + jitter * baseMs / 2);
    }

    private TimeSpan? RetryAfter(HttpResponseMessage response)
    {
        var ra = response.Headers.RetryAfter;
        if (ra?.Delta is { } delta) return delta < TimeSpan.Zero ? TimeSpan.Zero : delta;
        if (ra?.Date is { } date)
        {
            var d = date - Time.GetUtcNow();
            return d < TimeSpan.Zero ? TimeSpan.Zero : d;
        }

        return null;
    }

    /// <summary>Ensembl-style quota headers: when the window is exhausted, pause the host until it resets.</summary>
    private void ApplyRateLimitHeaders(HostRateLimiter limiter, string host, HttpResponseMessage response)
    {
        if (!response.Headers.TryGetValues("X-RateLimit-Remaining", out var remaining)) return;
        if (!long.TryParse(remaining.FirstOrDefault(), NumberStyles.Integer, CultureInfo.InvariantCulture, out long left) || left > 0) return;
        if (response.Headers.TryGetValues("X-RateLimit-Reset", out var reset)
            && double.TryParse(reset.FirstOrDefault(), NumberStyles.Float, CultureInfo.InvariantCulture, out double seconds) && seconds > 0)
            limiter.BlockUntil(host, Time.GetUtcNow() + TimeSpan.FromSeconds(Math.Min(seconds, Options.MaxRetryAfter.TotalSeconds)));
    }
}
