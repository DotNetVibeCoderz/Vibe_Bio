using System.Runtime.CompilerServices;
using System.Text.Json;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Data.Ncbi;

/// <summary>
/// NCBI Entrez Programming Utilities (E-utilities): EInfo, ESearch, ESummary, EFetch and ELink, with the history
/// server for large result sets.
/// </summary>
/// <remarks>
/// <para>
/// NCBI allows 3 requests per second per IP address, or 10 with an API key. The client spaces requests accordingly
/// through the shared <see cref="HostRateLimiter"/>. <see cref="ApiKey"/> and <see cref="Email"/> default to the
/// <c>NCBI_API_KEY</c> and <c>NCBI_EMAIL</c> environment variables; NCBI asks every tool to identify itself
/// with <c>tool</c> and <c>email</c>. These parameters are never part of a cache key, a cached file or an exception
/// message.
/// </para>
/// <para>Id lists longer than 200 are sent as a POST body, as NCBI recommends.</para>
/// </remarks>
public sealed class EntrezClient
{
    /// <summary>The E-utilities base URL.</summary>
    public const string DefaultBaseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

    private const int PostThreshold = 200;
    private static readonly TimeSpan SearchTtl = TimeSpan.FromDays(1);
    private static readonly TimeSpan HistoryTtl = TimeSpan.FromHours(1); // WebEnv sessions expire after a few hours
    private static readonly TimeSpan SummaryTtl = TimeSpan.FromDays(7);
    private static readonly TimeSpan RecordTtl = TimeSpan.FromDays(30);
    private readonly DataHttp _http;

    /// <param name="httpClient">An <see cref="HttpClient"/> to use (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public EntrezClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        _http = new DataHttp("NCBI", httpClient, options);
    }

    /// <summary>NCBI API key (raises the limit to 10 requests/s). Defaults to <c>$NCBI_API_KEY</c>.</summary>
    public string? ApiKey { get; init; } = NullIfEmpty(Environment.GetEnvironmentVariable("NCBI_API_KEY"));

    /// <summary>Contact e-mail sent as <c>email</c>. Defaults to <c>$NCBI_EMAIL</c>; never hard-code a personal address.</summary>
    public string? Email { get; init; } = NullIfEmpty(Environment.GetEnvironmentVariable("NCBI_EMAIL"));

    /// <summary>Tool name sent as <c>tool</c>.</summary>
    public string Tool { get; init; } = "Biotech.Net";

    /// <summary>The E-utilities base URL (ending in <c>/</c>).</summary>
    public string BaseUrl { get; init; } = DefaultBaseUrl;

    /// <summary>The request rate in effect: 10/s with an API key, 3/s without.</summary>
    public double RequestsPerSecond => _http.Options.RequestsPerSecond ?? (ApiKey is null ? 3 : 10);

    // ------------------------------------------------------------------ EInfo

    /// <summary>The names of all Entrez databases.</summary>
    public async Task<IReadOnlyList<string>> EInfoDatabasesAsync(CancellationToken cancellationToken = default)
    {
        using var doc = await GetJsonAsync("einfo.fcgi", new QueryString().Add("retmode", "json"), SearchTtl, cancellationToken).ConfigureAwait(false);
        return doc.RootElement.Path("einforesult")?.Strings("dblist") ?? [];
    }

    /// <summary>Statistics, searchable fields and links of one database.</summary>
    public async Task<EInfoDatabase> EInfoAsync(string database, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(database);
        using var doc = await GetJsonAsync("einfo.fcgi", new QueryString().Add("db", database).Add("retmode", "json"), SearchTtl, cancellationToken).ConfigureAwait(false);
        var info = doc.RootElement.Path("einforesult")?.Arr("dbinfo").FirstOrDefault()
                   ?? throw new DataServiceException("NCBI", $"EInfo returned no information for database '{database}'.");
        return new EInfoDatabase(
            info.Str("dbname") ?? database,
            info.Str("menuname") ?? database,
            info.Str("description") ?? "",
            info.Long("count") ?? 0,
            info.Str("lastupdate"),
            info.Arr("fieldlist").Select(f => new EInfoField(f.Str("name") ?? "", f.Str("fullname") ?? "", f.Str("description") ?? "")).ToList(),
            info.Arr("linklist").Select(l => new EInfoLink(l.Str("name") ?? "", l.Str("menu") ?? "", l.Str("description") ?? "", l.Str("dbto") ?? "")).ToList());
    }

    // ------------------------------------------------------------------ ESearch

    /// <summary>Searches <paramref name="database"/> with an Entrez query (<c>BRCA1[gene] AND human[orgn]</c>).</summary>
    /// <param name="retMax">Ids to return (Entrez caps this at 10,000; use the history server beyond that).</param>
    /// <param name="retStart">Offset of the first id.</param>
    /// <param name="useHistory">Keep the result set on the history server and return its handle.</param>
    /// <param name="sort">Database-specific sort order, e.g. <c>relevance</c> or <c>pub_date</c> for PubMed.</param>
    public async Task<ESearchResult> ESearchAsync(string database, string term, int retMax = 20, int retStart = 0,
        bool useHistory = false, string? sort = null, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(database);
        ArgumentException.ThrowIfNullOrWhiteSpace(term);
        var q = new QueryString().Add("db", database).Add("term", term).Add("retmax", retMax).Add("retstart", retStart)
            .Add("usehistory", useHistory ? "y" : null).Add("sort", sort).Add("retmode", "json");
        using var doc = await GetJsonAsync("esearch.fcgi", q, useHistory ? HistoryTtl : SearchTtl, cancellationToken).ConfigureAwait(false);
        var r = doc.RootElement.Path("esearchresult")
                ?? throw new DataServiceException("NCBI", "ESearch response has no 'esearchresult'.");
        if (r.Str("ERROR") is { Length: > 0 } err) throw new DataServiceException("NCBI", $"ESearch failed: {err}", serviceMessage: err);
        string? webEnv = r.Str("webenv");
        string? queryKey = r.Str("querykey");
        return new ESearchResult(database, r.Int("count") ?? 0, r.Int("retstart") ?? retStart, r.Int("retmax") ?? retMax,
            r.Strings("idlist"), r.Str("querytranslation"),
            webEnv is not null && queryKey is not null ? new EntrezHistory(webEnv, queryKey) : null);
    }

    // ------------------------------------------------------------------ ESummary

    /// <summary>Document summaries for a list of UIDs.</summary>
    public async Task<IReadOnlyList<ESummaryDocument>> ESummaryAsync(string database, IEnumerable<string> ids, CancellationToken cancellationToken = default)
    {
        var idList = IdList(ids);
        var q = new QueryString().Add("db", database).Add("id", string.Join(',', idList)).Add("retmode", "json");
        using var doc = await SendJsonAsync("esummary.fcgi", q, post: idList.Count > PostThreshold, SummaryTtl, cancellationToken).ConfigureAwait(false);
        return ParseSummaries(doc);
    }

    /// <summary>Document summaries for a slice of a history-server result set.</summary>
    public async Task<IReadOnlyList<ESummaryDocument>> ESummaryAsync(string database, EntrezHistory history, int retStart = 0, int retMax = 500,
        CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(history);
        var q = new QueryString().Add("db", database).Add("query_key", history.QueryKey).Add("WebEnv", history.WebEnv)
            .Add("retstart", retStart).Add("retmax", retMax).Add("retmode", "json");
        using var doc = await GetJsonAsync("esummary.fcgi", q, SummaryTtl, cancellationToken).ConfigureAwait(false);
        return ParseSummaries(doc);
    }

    /// <summary>PubMed citation summaries (title, journal, date, authors, DOI).</summary>
    public async Task<IReadOnlyList<PubMedSummary>> PubMedSummariesAsync(IEnumerable<string> pmids, CancellationToken cancellationToken = default)
    {
        var docs = await ESummaryAsync("pubmed", pmids, cancellationToken).ConfigureAwait(false);
        return docs.Select(d =>
        {
            var f = d.Fields;
            string? doi = f.Arr("articleids").FirstOrDefault(a => a.Str("idtype") == "doi") is { ValueKind: JsonValueKind.Object } a ? a.Str("value") : null;
            return new PubMedSummary(d.Uid, f.Str("title") ?? "", f.Str("source") ?? f.Str("fulljournalname") ?? "", f.Str("pubdate") ?? "",
                f.Arr("authors").Select(x => x.Str("name")).OfType<string>().ToList(),
                string.IsNullOrEmpty(doi) ? null : doi, NullIfEmpty(f.Str("volume")), NullIfEmpty(f.Str("issue")), NullIfEmpty(f.Str("pages")));
        }).ToList();
    }

    // ------------------------------------------------------------------ EFetch

    /// <summary>Raw EFetch output for a list of UIDs or accessions.</summary>
    /// <param name="retType">E.g. <c>fasta</c>, <c>gb</c>, <c>gp</c>, <c>abstract</c>, <c>xml</c>.</param>
    /// <param name="retMode">E.g. <c>text</c> or <c>xml</c>.</param>
    public async Task<string> EFetchAsync(string database, IEnumerable<string> ids, string retType, string retMode = "text",
        CancellationToken cancellationToken = default)
    {
        var idList = IdList(ids);
        var q = new QueryString().Add("db", database).Add("id", string.Join(',', idList)).Add("rettype", retType).Add("retmode", retMode);
        var response = await SendAsync("efetch.fcgi", q, post: idList.Count > PostThreshold, RecordTtl, cancellationToken).ConfigureAwait(false);
        return response.Text();
    }

    /// <summary>Raw EFetch output for a slice of a history-server result set.</summary>
    public async Task<string> EFetchAsync(string database, EntrezHistory history, string retType, string retMode = "text",
        int retStart = 0, int retMax = 200, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(history);
        var q = new QueryString().Add("db", database).Add("query_key", history.QueryKey).Add("WebEnv", history.WebEnv)
            .Add("retstart", retStart).Add("retmax", retMax).Add("rettype", retType).Add("retmode", retMode);
        var response = await SendAsync("efetch.fcgi", q, post: false, RecordTtl, cancellationToken).ConfigureAwait(false);
        return response.Text();
    }

    /// <summary>
    /// Fetches nucleotide (<c>nuccore</c>) or protein (<c>protein</c>) records and parses them into
    /// <see cref="Sequence"/> objects (GenBank records keep their annotations and features).
    /// </summary>
    public async Task<IReadOnlyList<Sequence>> FetchSequencesAsync(string database, IEnumerable<string> ids,
        EntrezSequenceFormat format = EntrezSequenceFormat.Fasta, CancellationToken cancellationToken = default)
    {
        string text = await EFetchAsync(database, ids, RetType(format), "text", cancellationToken).ConfigureAwait(false);
        return ParseSequences(database, text, format);
    }

    /// <summary>Fetches one record by accession or UID.</summary>
    public async Task<Sequence> FetchSequenceAsync(string database, string id, EntrezSequenceFormat format = EntrezSequenceFormat.Fasta,
        CancellationToken cancellationToken = default)
    {
        var list = await FetchSequencesAsync(database, [id], format, cancellationToken).ConfigureAwait(false);
        return list.Count > 0 ? list[0] : throw new DataServiceException("NCBI", $"EFetch returned no {database} record for '{id}'.");
    }

    /// <summary>
    /// Streams every record of a history-server search (see <see cref="ESearchAsync"/> with <c>useHistory: true</c>)
    /// in batches of <paramref name="batchSize"/>.
    /// </summary>
    public async IAsyncEnumerable<Sequence> FetchSequencesAsync(ESearchResult search, EntrezSequenceFormat format = EntrezSequenceFormat.Fasta,
        int batchSize = 200, [EnumeratorCancellation] CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(search);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(batchSize);
        var history = search.History ?? throw new ArgumentException("The search was run without useHistory: true.", nameof(search));
        for (int start = 0; start < search.Count; start += batchSize)
        {
            string text = await EFetchAsync(search.Database, history, RetType(format), "text", start, batchSize, cancellationToken).ConfigureAwait(false);
            foreach (var s in ParseSequences(search.Database, text, format)) yield return s;
        }
    }

    // ------------------------------------------------------------------ ELink

    /// <summary>Records in <paramref name="database"/> linked to <paramref name="ids"/> in <paramref name="dbFrom"/>.</summary>
    /// <param name="linkName">Restrict to one link, e.g. <c>protein_gene</c>; null returns every link to the target database.</param>
    public async Task<IReadOnlyList<ELinkSet>> ELinkAsync(string dbFrom, string database, IEnumerable<string> ids, string? linkName = null,
        CancellationToken cancellationToken = default)
    {
        var idList = IdList(ids);
        var q = new QueryString().Add("dbfrom", dbFrom).Add("db", database).Add("id", string.Join(',', idList)).Add("linkname", linkName).Add("retmode", "json");
        using var doc = await SendJsonAsync("elink.fcgi", q, post: idList.Count > PostThreshold, SummaryTtl, cancellationToken).ConfigureAwait(false);
        return doc.RootElement.Arr("linksets").Select(set => new ELinkSet(
            set.Str("dbfrom") ?? dbFrom,
            set.Strings("ids"),
            set.Arr("linksetdbs").Select(l => new ELinkResult(l.Str("dbto") ?? database, l.Str("linkname") ?? "", l.Strings("links"))).ToList())).ToList();
    }

    // ------------------------------------------------------------------ plumbing

    private static string RetType(EntrezSequenceFormat format) => format == EntrezSequenceFormat.GenBank ? "gb" : "fasta";

    private static IReadOnlyList<Sequence> ParseSequences(string database, string text, EntrezSequenceFormat format)
    {
        if (format == EntrezSequenceFormat.GenBank) return GenBankReader.Parse(text);
        var alphabet = database.Equals("protein", StringComparison.OrdinalIgnoreCase) ? Alphabet.ProteinExtended : null;
        return FastaReader.Parse(text, alphabet);
    }

    private static List<string> IdList(IEnumerable<string> ids)
    {
        ArgumentNullException.ThrowIfNull(ids);
        var list = ids.Select(i => i.Trim()).Where(i => i.Length > 0).ToList();
        if (list.Count == 0) throw new ArgumentException("At least one id is required.", nameof(ids));
        return list;
    }

    private static IReadOnlyList<ESummaryDocument> ParseSummaries(JsonDocument doc)
    {
        var root = doc.RootElement;
        if (root.Str("error") is { Length: > 0 } err) throw new DataServiceException("NCBI", $"ESummary failed: {err}", serviceMessage: err);
        var result = root.Path("result");
        if (result is null) return [];
        var list = new List<ESummaryDocument>();
        foreach (string uid in result.Value.Strings("uids"))
        {
            if (result.Value.Prop(uid) is not { } d) continue;
            if (d.Str("error") is { Length: > 0 } e) throw new DataServiceException("NCBI", $"ESummary failed for UID {uid}: {e}", serviceMessage: e);
            list.Add(new ESummaryDocument(uid, d.Clone()));
        }

        return list;
    }

    private QueryString WithIdentity(QueryString q) =>
        q.Add("tool", Tool).Add("email", Email).Add("api_key", ApiKey);

    private async Task<JsonDocument> GetJsonAsync(string endpoint, QueryString q, TimeSpan ttl, CancellationToken ct) =>
        await SendJsonAsync(endpoint, q, post: false, ttl, ct).ConfigureAwait(false);

    private async Task<JsonDocument> SendJsonAsync(string endpoint, QueryString q, bool post, TimeSpan ttl, CancellationToken ct)
    {
        var response = await SendAsync(endpoint, q, post, ttl, ct).ConfigureAwait(false);
        try
        {
            return response.Json();
        }
        catch (JsonException ex)
        {
            throw new DataServiceException("NCBI", $"{endpoint} did not return JSON: {DataServiceException.ExtractMessage(response.Content, response.ContentType)}",
                DataServiceException.Redact(response.RequestUri), innerException: ex);
        }
    }

    private Task<DataResponse> SendAsync(string endpoint, QueryString q, bool post, TimeSpan ttl, CancellationToken ct)
    {
        DataRequest request;
        if (post)
        {
            // Identity parameters stay in the URL so the POST body (and thus the cache key) does not depend on them.
            var url = WithIdentity(new QueryString()).Apply(BaseUrl + endpoint);
            request = DataRequest.PostForm(url, q.ToString());
        }
        else
        {
            request = DataRequest.Get(WithIdentity(q).Apply(BaseUrl + endpoint));
        }

        return _http.SendAsync(request with { RequestsPerSecond = ApiKey is null ? 3 : 10, CacheTtl = ttl }, ct);
    }

    private static string? NullIfEmpty(string? s) => string.IsNullOrWhiteSpace(s) ? null : s;
}
