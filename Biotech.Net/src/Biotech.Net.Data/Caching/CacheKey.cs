using System.Security.Cryptography;
using System.Text;

namespace Biotech.Net.Data.Caching;

/// <summary>
/// Cache keys: a SHA-256 over the HTTP method, the normalised URL and a hash of the request body.
/// </summary>
/// <remarks>
/// Normalisation lower-cases the scheme and host, drops the default port, sorts query parameters (so
/// <c>?a=1&amp;b=2</c> and <c>?b=2;a=1</c> share an entry) and removes credentials and identity parameters
/// (<c>api_key</c>, <c>email</c>, <c>tool</c>), so the key never depends on, or reveals, who asked.
/// </remarks>
public static class CacheKey
{
    private static readonly HashSet<string> s_credentials = new(StringComparer.OrdinalIgnoreCase) { "api_key", "email", "tool" };

    /// <summary>True for query parameters that identify the caller rather than the resource.</summary>
    public static bool IsCredential(string parameterName) => s_credentials.Contains(Uri.UnescapeDataString(parameterName));

    /// <summary>The URL in canonical form (see remarks on <see cref="CacheKey"/>).</summary>
    public static string NormalizeUrl(Uri uri)
    {
        ArgumentNullException.ThrowIfNull(uri);
        if (!uri.IsAbsoluteUri) throw new ArgumentException("An absolute URL is required.", nameof(uri));
        var sb = new StringBuilder();
        sb.Append(uri.Scheme.ToLowerInvariant()).Append("://").Append(uri.Host.ToLowerInvariant());
        if (!uri.IsDefaultPort) sb.Append(':').Append(uri.Port);
        sb.Append(uri.AbsolutePath);
        var query = uri.Query.TrimStart('?')
            .Split(['&', ';'], StringSplitOptions.RemoveEmptyEntries)
            .Where(p => !IsCredential(p.Split('=', 2)[0]))
            .Order(StringComparer.Ordinal)
            .ToList();
        if (query.Count > 0) sb.Append('?').Append(string.Join('&', query));
        return sb.ToString();
    }

    /// <summary>The cache key (64 lower-case hex digits) for a request.</summary>
    public static string Compute(string method, Uri uri, ReadOnlySpan<byte> body = default)
    {
        string bodyHash = body.IsEmpty ? "-" : Convert.ToHexStringLower(SHA256.HashData(body));
        string material = $"{method.ToUpperInvariant()} {NormalizeUrl(uri)}\n{bodyHash}";
        return Convert.ToHexStringLower(SHA256.HashData(Encoding.UTF8.GetBytes(material)));
    }
}
