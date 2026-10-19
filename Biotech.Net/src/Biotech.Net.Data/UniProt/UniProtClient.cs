using System.Runtime.CompilerServices;
using System.Text.Json;
using System.Text.RegularExpressions;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Data.UniProt;

/// <summary>
/// UniProt REST API (<c>rest.uniprot.org</c>): UniProtKB entries (JSON and FASTA), searches with the UniProt query
/// syntax and cursor pagination, and the asynchronous ID-mapping service.
/// </summary>
public sealed partial class UniProtClient
{
    /// <summary>The REST base URL.</summary>
    public const string DefaultBaseUrl = "https://rest.uniprot.org/";

    private static readonly TimeSpan EntryTtl = TimeSpan.FromDays(7);
    private static readonly TimeSpan SearchTtl = TimeSpan.FromDays(1);
    private readonly DataHttp _http;

    /// <param name="httpClient">An <see cref="HttpClient"/> to use (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public UniProtClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        _http = new DataHttp("UniProt", httpClient, options);
    }

    /// <summary>The REST base URL (ending in <c>/</c>).</summary>
    public string BaseUrl { get; init; } = DefaultBaseUrl;

    /// <summary>Interval between ID-mapping status polls (grows by half each poll, up to 5 s).</summary>
    public TimeSpan IdMappingPollInterval { get; init; } = TimeSpan.FromSeconds(1);

    /// <summary>Gives up waiting for an ID-mapping job after this long.</summary>
    public TimeSpan IdMappingTimeout { get; init; } = TimeSpan.FromMinutes(5);

    // ------------------------------------------------------------------ entries

    /// <summary>A UniProtKB entry by accession (or entry name), parsed from JSON.</summary>
    public async Task<UniProtEntry> GetEntryAsync(string accession, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(accession);
        var url = new QueryString().Add("format", "json").Apply($"{BaseUrl}uniprotkb/{Uri.EscapeDataString(accession.Trim())}");
        var response = await SendAsync(DataRequest.Get(url) with { Accept = "application/json" }, EntryTtl, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        if (doc.RootElement.Str("entryType") is "Inactive")
            throw new DataServiceException("UniProt", $"UniProtKB entry {accession} is inactive ({doc.RootElement.Path("inactiveReason")?.Str("inactiveReasonType")}).",
                url, System.Net.HttpStatusCode.Gone);
        return UniProtEntry.FromJson(doc.RootElement);
    }

    /// <summary>The canonical sequence of an entry, from the FASTA endpoint.</summary>
    public async Task<Sequence> GetFastaAsync(string accession, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(accession);
        var url = new QueryString().Add("format", "fasta").Apply($"{BaseUrl}uniprotkb/{Uri.EscapeDataString(accession.Trim())}");
        var response = await SendAsync(DataRequest.Get(url), EntryTtl, cancellationToken).ConfigureAwait(false);
        var list = FastaReader.Parse(response.Text(), Alphabet.ProteinExtended);
        return list.Count > 0 ? list[0] : throw new DataServiceException("UniProt", $"No FASTA record for '{accession}'.", url);
    }

    /// <summary>The entry as a Core <see cref="Protein"/>.</summary>
    public async Task<Protein> GetProteinAsync(string accession, CancellationToken cancellationToken = default) =>
        (await GetEntryAsync(accession, cancellationToken).ConfigureAwait(false)).ToProtein();

    // ------------------------------------------------------------------ search

    /// <summary>One page of a UniProtKB search.</summary>
    /// <param name="query">UniProt query syntax, e.g. <c>gene:BRCA1 AND organism_id:9606 AND reviewed:true</c>.</param>
    /// <param name="fields">Return fields (e.g. <c>accession,protein_name,sequence</c>); null returns full entries.</param>
    /// <param name="pageSize">Entries per page (UniProt allows up to 500).</param>
    public async Task<UniProtSearchPage> SearchPageAsync(string query, string? fields = null, int pageSize = 25, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(query);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(pageSize);
        var url = new QueryString().Add("query", query).Add("fields", fields).Add("format", "json").Add("size", Math.Min(pageSize, 500))
            .Apply($"{BaseUrl}uniprotkb/search");
        return await PageAsync(url, cancellationToken).ConfigureAwait(false);
    }

    /// <summary>The page after <paramref name="page"/>, or null when it was the last.</summary>
    public async Task<UniProtSearchPage?> NextPageAsync(UniProtSearchPage page, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(page);
        return page.NextPage is null ? null : await PageAsync(page.NextPage, cancellationToken).ConfigureAwait(false);
    }

    /// <summary>Streams search results across pages, following the <c>Link</c> header.</summary>
    /// <param name="maxResults">Stop after this many entries (null: all of them).</param>
    public async IAsyncEnumerable<UniProtEntry> SearchAsync(string query, string? fields = null, int pageSize = 100, int? maxResults = null,
        [EnumeratorCancellation] CancellationToken cancellationToken = default)
    {
        int yielded = 0;
        if (maxResults is <= 0) yield break;
        UniProtSearchPage? page = await SearchPageAsync(query, fields, pageSize, cancellationToken).ConfigureAwait(false);
        while (page is not null)
        {
            foreach (var e in page.Entries)
            {
                yielded++;
                yield return e;
                if (yielded >= maxResults) yield break; // no request for a page nobody will read
            }

            page = await NextPageAsync(page, cancellationToken).ConfigureAwait(false);
        }
    }

    private async Task<UniProtSearchPage> PageAsync(Uri url, CancellationToken ct)
    {
        var response = await SendAsync(DataRequest.Get(url) with { Accept = "application/json" }, SearchTtl, ct).ConfigureAwait(false);
        using var doc = response.Json();
        var entries = doc.RootElement.Arr("results").Select(UniProtEntry.FromJson).ToList();
        long? total = long.TryParse(response.Header("X-Total-Results"), out long t) ? t : null;
        return new UniProtSearchPage(entries, total, NextLink(response.Header("Link")));
    }

    /// <summary>Extracts the <c>rel="next"</c> target of a <c>Link</c> header.</summary>
    internal static Uri? NextLink(string? linkHeader)
    {
        if (string.IsNullOrEmpty(linkHeader)) return null;
        var m = NextLinkRegex().Match(linkHeader);
        return m.Success ? new Uri(m.Groups[1].Value) : null;
    }

    [GeneratedRegex("<([^>]+)>\\s*;\\s*rel=\"?next\"?", RegexOptions.IgnoreCase)]
    private static partial Regex NextLinkRegex();

    // ------------------------------------------------------------------ ID mapping

    /// <summary>Submits an ID-mapping job and returns its job id.</summary>
    /// <param name="from">Source database, e.g. <c>UniProtKB_AC-ID</c>, <c>Gene_Name</c>, <c>Ensembl</c>, <c>RefSeq_Protein</c>.</param>
    /// <param name="to">Target database, e.g. <c>GeneID</c>, <c>Ensembl</c>, <c>PDB</c>, <c>UniProtKB</c>.</param>
    /// <param name="ids">Identifiers to map.</param>
    /// <param name="taxonId">Restricts gene-name mappings to one organism.</param>
    public async Task<string> SubmitIdMappingAsync(string from, string to, IEnumerable<string> ids, long? taxonId = null,
        CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(from);
        ArgumentException.ThrowIfNullOrWhiteSpace(to);
        var list = ids.Select(i => i.Trim()).Where(i => i.Length > 0).ToList();
        if (list.Count == 0) throw new ArgumentException("At least one id is required.", nameof(ids));
        var form = new QueryString().Add("from", from).Add("to", to).Add("ids", string.Join(',', list))
            .Add("taxId", taxonId?.ToString(System.Globalization.CultureInfo.InvariantCulture));
        var response = await SendAsync(DataRequest.PostForm(new Uri($"{BaseUrl}idmapping/run"), form.ToString()), null, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        return doc.RootElement.Str("jobId") ?? throw new DataServiceException("UniProt", "ID mapping submission returned no jobId.");
    }

    /// <summary>The state of an ID-mapping job.</summary>
    public async Task<IdMappingJobState> GetIdMappingStatusAsync(string jobId, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(jobId);
        var response = await SendAsync(DataRequest.Get(new Uri($"{BaseUrl}idmapping/status/{Uri.EscapeDataString(jobId)}")), null, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        var root = doc.RootElement;
        // A finished job may redirect straight to its results.
        if (root.Prop("results") is not null || root.Prop("failedIds") is not null) return IdMappingJobState.Finished;
        string? status = root.Str("jobStatus");
        return status?.ToUpperInvariant() switch
        {
            "NEW" => IdMappingJobState.New,
            "RUNNING" => IdMappingJobState.Running,
            "FINISHED" => IdMappingJobState.Finished,
            "ERROR" or "FAILURE" => IdMappingJobState.Error,
            _ => root.Prop("errors") is not null || root.Prop("messages") is not null
                ? throw new DataServiceException("UniProt", $"ID mapping job {jobId} failed: {DataServiceException.ExtractMessage(response.Content, response.ContentType)}")
                : IdMappingJobState.Running,
        };
    }

    /// <summary>The mappings of a finished job (all pages).</summary>
    public async Task<IdMappingResult> GetIdMappingResultsAsync(string jobId, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(jobId);
        var mappings = new List<IdMapping>();
        var failed = new List<string>();
        Uri? url = new QueryString().Add("size", 500).Apply($"{BaseUrl}idmapping/results/{Uri.EscapeDataString(jobId)}");
        while (url is not null)
        {
            var response = await SendAsync(DataRequest.Get(url) with { Accept = "application/json" }, null, cancellationToken).ConfigureAwait(false);
            using var doc = response.Json();
            foreach (var r in doc.RootElement.Arr("results"))
            {
                string from = r.Str("from") ?? "";
                var to = r.Prop("to");
                // Mapping to UniProtKB returns whole entries; other targets return plain identifiers.
                string? target = to is { ValueKind: JsonValueKind.Object } o ? o.Str("primaryAccession") : to?.AsString();
                if (target is not null) mappings.Add(new IdMapping(from, target));
            }

            failed.AddRange(doc.RootElement.Strings("failedIds"));
            url = NextLink(response.Header("Link"));
        }

        return new IdMappingResult(jobId, mappings, failed);
    }

    /// <summary>
    /// Maps identifiers end to end: submits the job, polls its status (through the configured
    /// <see cref="TimeProvider"/>) until it finishes, then collects every page of results.
    /// </summary>
    public async Task<IdMappingResult> MapIdsAsync(string from, string to, IEnumerable<string> ids, long? taxonId = null,
        CancellationToken cancellationToken = default)
    {
        string jobId = await SubmitIdMappingAsync(from, to, ids, taxonId, cancellationToken).ConfigureAwait(false);
        var time = _http.Time;
        var deadline = time.GetUtcNow() + IdMappingTimeout;
        var interval = IdMappingPollInterval;
        while (true)
        {
            var state = await GetIdMappingStatusAsync(jobId, cancellationToken).ConfigureAwait(false);
            if (state == IdMappingJobState.Finished) break;
            if (state == IdMappingJobState.Error) throw new DataServiceException("UniProt", $"ID mapping job {jobId} ended in an error.");
            if (time.GetUtcNow() >= deadline)
                throw new DataServiceException("UniProt", $"ID mapping job {jobId} did not finish within {IdMappingTimeout.TotalSeconds:F0} s.",
                    innerException: new TimeoutException());
            await Task.Delay(interval, time, cancellationToken).ConfigureAwait(false);
            interval = TimeSpan.FromTicks(Math.Min(interval.Ticks * 3 / 2, TimeSpan.FromSeconds(5).Ticks));
        }

        return await GetIdMappingResultsAsync(jobId, cancellationToken).ConfigureAwait(false);
    }

    private Task<DataResponse> SendAsync(DataRequest request, TimeSpan? ttl, CancellationToken ct) =>
        _http.SendAsync(request with { RequestsPerSecond = 10, CacheTtl = ttl }, ct);
}
