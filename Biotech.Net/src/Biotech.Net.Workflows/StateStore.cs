using System.Security.Cryptography;
using System.Text;
using System.Text.Json;

namespace Biotech.Net.Workflows;

/// <summary>How file inputs and outputs are fingerprinted for up-to-date checks.</summary>
public enum HashMode
{
    /// <summary>SHA-256 of the content (cached by size and modification time, so unchanged files are not re-read). Touching a file without changing it does not rerun anything.</summary>
    Content,

    /// <summary>Size and modification time only: cheaper for very large files, but a touch reruns downstream jobs.</summary>
    Timestamp,
}

/// <summary>What the state store remembers about a job's last successful run.</summary>
internal sealed record JobRecord
{
    public required string Job { get; init; }
    public required string Signature { get; init; }
    public required string Definition { get; init; }
    public Dictionary<string, string> Inputs { get; init; } = [];
    public Dictionary<string, string> Outputs { get; init; } = [];
    public DateTimeOffset FinishedAt { get; init; }
}

/// <summary>
/// The persistent state in <c>.biotech-workflow/</c>: job records, incomplete markers, stored values, the hash cache,
/// logs and run reports.
/// </summary>
internal sealed class StateStore
{
    private static readonly JsonSerializerOptions Json = new() { WriteIndented = true };
    private readonly Lock _gate = new();
    private readonly Workflow _wf;
    private Dictionary<string, CachedHash> _hashes = [];
    private bool _hashesDirty;

    private sealed record CachedHash(long Length, long Ticks, string Sha);

    public StateStore(Workflow wf)
    {
        _wf = wf;
        Root = wf.StateDirectory;
        string cache = Path.Combine(Root, "hashes.json");
        if (File.Exists(cache))
        {
            try
            {
                _hashes = JsonSerializer.Deserialize<Dictionary<string, CachedHash>>(File.ReadAllText(cache)) ?? [];
            }
            catch (JsonException)
            {
                _hashes = [];
            }
        }
    }

    public string Root { get; }

    private string Dir(string name)
    {
        string d = Path.Combine(Root, name);
        Directory.CreateDirectory(d);
        return d;
    }

    public static string Sha(string text) => Sha(Encoding.UTF8.GetBytes(text));

    public static string Sha(ReadOnlySpan<byte> bytes) => Convert.ToHexStringLower(SHA256.HashData(bytes));

    private static string FileKey(string jobId) => Sha(jobId)[..24];

    public JobRecord? Load(WorkflowJob job)
    {
        string path = Path.Combine(Root, "jobs", FileKey(job.Id) + ".json");
        if (!File.Exists(path)) return null;
        try
        {
            var r = JsonSerializer.Deserialize<JobRecord>(File.ReadAllText(path));
            return r?.Job == job.Id ? r : null;
        }
        catch (JsonException)
        {
            return null;
        }
    }

    public void Save(WorkflowJob job, JobRecord record)
        => WriteAtomic(Path.Combine(Dir("jobs"), FileKey(job.Id) + ".json"), Encoding.UTF8.GetBytes(JsonSerializer.Serialize(record, Json)));

    public void Forget(WorkflowJob job)
    {
        string path = Path.Combine(Root, "jobs", FileKey(job.Id) + ".json");
        if (File.Exists(path)) File.Delete(path);
    }

    private string MarkerPath(WorkflowJob job) => Path.Combine(Root, "incomplete", FileKey(job.Id));

    public bool IsIncomplete(WorkflowJob job) => File.Exists(MarkerPath(job));

    public void MarkIncomplete(WorkflowJob job)
    {
        Dir("incomplete");
        File.WriteAllText(MarkerPath(job), job.Id);
    }

    public void ClearIncomplete(WorkflowJob job)
    {
        if (File.Exists(MarkerPath(job))) File.Delete(MarkerPath(job));
    }

    private string ValuePath(string key) => Path.Combine(Root, "values", Sha(key)[..32] + ".json");

    public bool ValueExists(string key) => File.Exists(ValuePath(key));

    public byte[] ReadValue(string key)
        => File.Exists(ValuePath(key)) ? File.ReadAllBytes(ValuePath(key)) : throw new FileNotFoundException($"Value '{key}' has not been produced.");

    public void WriteValue(string key, byte[] json)
    {
        Dir("values");
        WriteAtomic(ValuePath(key), json);
    }

    public void DeleteValue(string key)
    {
        if (File.Exists(ValuePath(key))) File.Delete(ValuePath(key));
    }

    /// <summary>The fingerprint of an artifact, or null when it does not exist.</summary>
    public string? Fingerprint(ArtifactKey key, HashMode mode)
    {
        if (key.Kind == ArtifactKind.Value)
            return ValueExists(key.Id) ? "sha256:" + Sha(File.ReadAllBytes(ValuePath(key.Id))) : null;

        var info = new FileInfo(_wf.FullPath(key));
        if (!info.Exists) return null;
        long ticks = info.LastWriteTimeUtc.Ticks;
        if (mode == HashMode.Timestamp) return $"size:{info.Length};mtime:{ticks}";
        lock (_gate)
        {
            if (_hashes.TryGetValue(info.FullName, out var c) && c.Length == info.Length && c.Ticks == ticks) return "sha256:" + c.Sha;
        }

        string sha;
        using (var s = info.OpenRead()) sha = Convert.ToHexStringLower(SHA256.HashData(s));
        lock (_gate)
        {
            _hashes[info.FullName] = new CachedHash(info.Length, ticks, sha);
            _hashesDirty = true;
        }

        return "sha256:" + sha;
    }

    public void SaveHashes()
    {
        byte[] bytes;
        lock (_gate)
        {
            if (!_hashesDirty) return;
            bytes = JsonSerializer.SerializeToUtf8Bytes(_hashes);
            _hashesDirty = false;
        }

        Directory.CreateDirectory(Root);
        WriteAtomic(Path.Combine(Root, "hashes.json"), bytes);
    }

    public string LogPath(string runId, WorkflowJob job)
    {
        var safe = new StringBuilder();
        foreach (char c in job.Id) safe.Append(char.IsLetterOrDigit(c) || c is '-' or '_' or '.' ? c : '_');
        return Path.Combine(Dir(Path.Combine("logs", runId)), safe + ".log");
    }

    public string RunDirectory(string runId) => Dir(Path.Combine("runs", runId));

    /// <summary>Writes to a temporary file beside <paramref name="path"/> and renames it into place.</summary>
    public static void WriteAtomic(string path, byte[] bytes)
    {
        string tmp = path + "." + Guid.NewGuid().ToString("N")[..8] + ".tmp";
        File.WriteAllBytes(tmp, bytes);
        File.Move(tmp, path, overwrite: true);
    }
}
