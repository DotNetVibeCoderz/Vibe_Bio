using System.Net;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;

namespace Biotech.Net.Data.Http;

/// <summary>A request to a remote data service failed: an HTTP error, a service-reported error, a timeout or a network fault.</summary>
public class DataServiceException : Exception
{
    public DataServiceException(string service, string message, Uri? requestUri = null, HttpStatusCode? statusCode = null,
        string? serviceMessage = null, Exception? innerException = null)
        : base(message, innerException)
    {
        Service = service;
        RequestUri = requestUri;
        StatusCode = statusCode;
        ServiceMessage = serviceMessage;
    }

    /// <summary>The service name, e.g. <c>NCBI</c>, <c>UniProt</c>, <c>RCSB PDB</c>, <c>Ensembl</c>.</summary>
    public string Service { get; }

    /// <summary>The URL that failed, with credentials (<c>api_key</c>) removed.</summary>
    public Uri? RequestUri { get; }

    /// <summary>The HTTP status, when the server answered.</summary>
    public HttpStatusCode? StatusCode { get; }

    /// <summary>The error text the service returned, when it could be extracted.</summary>
    public string? ServiceMessage { get; }

    /// <summary>True for HTTP 404 (unknown accession, identifier or file).</summary>
    public bool IsNotFound => StatusCode == HttpStatusCode.NotFound;

    internal static DataServiceException FromResponse(string service, Uri uri, HttpStatusCode status, byte[] body, string? contentType, int attempts)
    {
        string? detail = ExtractMessage(body, contentType);
        var sb = new StringBuilder();
        sb.Append(service).Append(" returned HTTP ").Append((int)status).Append(' ').Append(status).Append(" for ").Append(Redact(uri));
        if (attempts > 1) sb.Append(" after ").Append(attempts).Append(" attempts");
        if (!string.IsNullOrWhiteSpace(detail)) sb.Append(": ").Append(detail);
        return new DataServiceException(service, sb.ToString(), Redact(uri), status, detail);
    }

    /// <summary>Pulls a human-readable message from a JSON, HTML or plain-text error body.</summary>
    internal static string? ExtractMessage(byte[] body, string? contentType)
    {
        if (body.Length == 0) return null;
        string text = Encoding.UTF8.GetString(body, 0, Math.Min(body.Length, 64 * 1024)).Trim();
        if (text.StartsWith('{') || text.StartsWith('['))
        {
            try
            {
                using var doc = JsonDocument.Parse(text);
                var root = doc.RootElement;
                if (root.ValueKind == JsonValueKind.Object)
                {
                    foreach (string name in (string[])["error", "message", "ERROR", "messages", "errors"])
                    {
                        if (!root.TryGetProperty(name, out var v)) continue;
                        string? m = v.ValueKind switch
                        {
                            JsonValueKind.String => v.GetString(),
                            JsonValueKind.Array => string.Join("; ", v.EnumerateArray().Select(e =>
                                e.ValueKind == JsonValueKind.String ? e.GetString() : e.TryGetProperty("message", out var mm) ? mm.ToString() : e.ToString())),
                            _ => v.ToString(),
                        };
                        if (!string.IsNullOrWhiteSpace(m)) return Truncate(m);
                    }
                }
            }
            catch (JsonException)
            {
                // Not JSON after all: fall through to the text form.
            }
        }

        if (text.StartsWith('<') || contentType?.Contains("html", StringComparison.OrdinalIgnoreCase) == true)
        {
            var title = Regex.Match(text, "<title>(.*?)</title>", RegexOptions.IgnoreCase | RegexOptions.Singleline);
            return title.Success ? "HTML error page: " + WebUtility.HtmlDecode(title.Groups[1].Value.Trim()) : "HTML error page";
        }

        // NCBI sends some errors URL-encoded ("Cannot+retrieve+history+data%3A+...").
        if (text.Contains('+') && text.Contains('%') && !text.Contains(' '))
            text = WebUtility.UrlDecode(text);
        return Truncate(text);
    }

    /// <summary>Removes credential query parameters so that exception messages and logs never carry an API key.</summary>
    internal static Uri Redact(Uri uri)
    {
        if (string.IsNullOrEmpty(uri.Query)) return uri;
        var parts = uri.Query.TrimStart('?').Split('&', ';')
            .Where(p => !Caching.CacheKey.IsCredential(p.Split('=', 2)[0]));
        var b = new UriBuilder(uri) { Query = string.Join('&', parts) };
        return b.Uri;
    }

    private static string Truncate(string s)
    {
        s = Regex.Replace(s, @"\s+", " ").Trim();
        return s.Length <= 400 ? s : s[..400] + "...";
    }
}

/// <summary>
/// Thrown in <see cref="DataCacheMode.Offline"/> mode when the response is not in the cache.
/// </summary>
public sealed class DataOfflineException : DataServiceException
{
    public DataOfflineException(string service, Uri requestUri, string reason)
        : base(service, $"Offline mode: {reason} for {Redact(requestUri)}. Run the same call once online (or switch " +
                        $"{nameof(DataClientOptions.CacheMode)} away from {nameof(DataCacheMode.Offline)}, or unset {DataClientOptions.OfflineEnvironmentVariable}).",
            Redact(requestUri))
    {
    }
}
