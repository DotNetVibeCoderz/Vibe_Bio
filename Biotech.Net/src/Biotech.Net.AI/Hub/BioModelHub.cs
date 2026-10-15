using System.Net.Http;
using System.Net.Sockets;
using Gravicode.HFNet.GraviHub;

namespace Biotech.Net.AI.Hub;

/// <summary>
/// Thrown when a model or file is not in the local cache and cannot be downloaded — offline mode, no network,
/// a gated repository without a token, or a wrong repository id.
/// </summary>
public sealed class BioModelUnavailableException : Exception
{
    /// <summary>Creates the exception.</summary>
    public BioModelUnavailableException(string repoId, string message, Exception? inner = null)
        : base(message, inner) => RepoId = repoId;

    /// <summary>The repository that could not be obtained.</summary>
    public string RepoId { get; }
}

/// <summary>Settings for <see cref="BioModelHub"/>.</summary>
public sealed record BioHubOptions
{
    /// <summary>Where checkpoints are cached. Null keeps HF.Net's default (<c>HF_HUB_CACHE</c>, <c>HF_HOME</c>, then <c>%LOCALAPPDATA%/huggingface/hub</c>).</summary>
    public string? CacheDirectory { get; init; }

    /// <summary>Hugging Face access token. Null falls back to the <c>HF_TOKEN</c> environment variable.</summary>
    public string? Token { get; init; }

    /// <summary>Never touch the network: use cached files or fail with <see cref="BioModelUnavailableException"/>.</summary>
    public bool Offline { get; init; }

    /// <summary>Per-request timeout; large checkpoints need a generous value.</summary>
    public TimeSpan Timeout { get; init; } = TimeSpan.FromMinutes(30);
}

/// <summary>
/// Downloads and caches checkpoints through HF.Net's GraviHub, with an offline mode and errors that say what to do.
/// </summary>
public static class BioModelHub
{
    private static readonly object Gate = new();
    private static BioHubOptions _options = new() { Offline = IsOfflineByEnvironment() };

    /// <summary>The settings in use.</summary>
    public static BioHubOptions Options
    {
        get { lock (Gate) return _options; }
    }

    /// <summary>The cache root currently in use.</summary>
    public static string CacheDirectory => Options.CacheDirectory ?? HubCache.DefaultRoot();

    /// <summary>
    /// Replaces the hub settings. Call before the first download. <c>HF_HUB_OFFLINE=1</c> in the environment also
    /// turns offline mode on.
    /// </summary>
    public static void Configure(BioHubOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        lock (Gate)
        {
            _options = options;
            var defaults = new HubOptions();
            var hub = new HubOptions
            {
                OfflineMode = options.Offline,
                Timeout = options.Timeout,
                CacheRoot = options.CacheDirectory ?? defaults.CacheRoot,
                Token = options.Token ?? defaults.Token,
            };
            Gravicode.HFNet.GraviHub.Hub.Configure(hub);
        }
    }

    /// <summary>The cache directory a repository snapshot lives in (it may not exist yet).</summary>
    public static string CachedDirectory(string repoId, string revision = "main")
    {
        ArgumentException.ThrowIfNullOrEmpty(repoId);
        return Gravicode.HFNet.GraviHub.Hub.Cache.DirectoryFor(repoId, RepoKind.Model, revision);
    }

    /// <summary>Whether a repository's <c>config.json</c> is already cached.</summary>
    public static bool IsCached(string repoId, string revision = "main")
        => File.Exists(Path.Combine(CachedDirectory(repoId, revision), "config.json"));

    /// <summary>
    /// Returns a local directory holding the model's config, weights and tokenizer files, downloading them when needed.
    /// </summary>
    /// <param name="repoId">Hub repository id, or an existing local directory.</param>
    /// <param name="revision">Branch, tag or commit.</param>
    /// <exception cref="BioModelUnavailableException">Offline and not cached, or the download failed.</exception>
    public static string EnsureModel(string repoId, string revision = "main")
    {
        ArgumentException.ThrowIfNullOrEmpty(repoId);
        if (Directory.Exists(repoId) && File.Exists(Path.Combine(repoId, "config.json"))) return repoId;
        string cached = CachedDirectory(repoId, revision);
        bool haveCache = File.Exists(Path.Combine(cached, "config.json")) && HasWeights(cached);
        if (Options.Offline)
        {
            return haveCache
                ? cached
                : throw new BioModelUnavailableException(repoId,
                    $"Model '{repoId}' is not in the cache ({cached}) and offline mode is on. Download it once while online, " +
                    "or point BioHubOptions.CacheDirectory at a cache that holds it.");
        }

        try
        {
            return Gravicode.HFNet.GraviHub.Hub.DownloadModel(repoId, revision);
        }
        catch (Exception ex) when (IsTransportFailure(ex))
        {
            if (haveCache) return cached;
            throw new BioModelUnavailableException(repoId, Explain(repoId, ex), ex);
        }
    }

    /// <summary>Returns the local path of one file of a repository, downloading it when needed.</summary>
    /// <exception cref="BioModelUnavailableException">Offline and not cached, or the download failed.</exception>
    public static string DownloadFile(string repoId, string fileName, string revision = "main")
    {
        ArgumentException.ThrowIfNullOrEmpty(repoId);
        ArgumentException.ThrowIfNullOrEmpty(fileName);
        string cached = Gravicode.HFNet.GraviHub.Hub.Cache.PathFor(repoId, RepoKind.Model, revision, fileName);
        if (Options.Offline)
        {
            return File.Exists(cached)
                ? cached
                : throw new BioModelUnavailableException(repoId, $"'{fileName}' of '{repoId}' is not cached ({cached}) and offline mode is on.");
        }

        try
        {
            return Gravicode.HFNet.GraviHub.Hub.DownloadFile(repoId, fileName, revision);
        }
        catch (Exception ex) when (IsTransportFailure(ex))
        {
            if (File.Exists(cached)) return cached;
            throw new BioModelUnavailableException(repoId, Explain(repoId, ex), ex);
        }
    }

    private static bool HasWeights(string dir)
        => Directory.EnumerateFiles(dir, "*.safetensors").Any() || File.Exists(Path.Combine(dir, "pytorch_model.bin"))
           || File.Exists(Path.Combine(dir, "model.safetensors.index.json"));

    private static bool IsTransportFailure(Exception ex)
        => ex is HubException or HttpRequestException or SocketException or TaskCanceledException or TimeoutException or IOException
           || ex.InnerException is not null && IsTransportFailure(ex.InnerException);

    private static string Explain(string repoId, Exception ex)
    {
        var hub = ex as HubException ?? ex.InnerException as HubException;
        string hint = hub?.Status switch
        {
            System.Net.HttpStatusCode.Unauthorized or System.Net.HttpStatusCode.Forbidden =>
                "The repository is gated or private: set HF_TOKEN (or BioHubOptions.Token) to a token that has accepted its terms.",
            System.Net.HttpStatusCode.NotFound => "Check the repository id and revision.",
            _ => "Check the network connection, or run offline against a cache that already holds the model.",
        };
        return $"Could not obtain '{repoId}' from the Hugging Face Hub: {ex.Message} {hint}";
    }

    private static bool IsOfflineByEnvironment()
    {
        string? v = Environment.GetEnvironmentVariable("HF_HUB_OFFLINE");
        return v is "1" or "true" or "TRUE" or "True" or "yes";
    }
}
