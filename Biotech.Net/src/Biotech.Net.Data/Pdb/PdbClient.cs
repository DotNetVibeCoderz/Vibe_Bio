using System.Text.Json;
using System.Text.Json.Nodes;
using Biotech.Net.Data.Http;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.Data.Pdb;

/// <summary>
/// RCSB Protein Data Bank: coordinate downloads (<c>files.rcsb.org</c>, gzip) parsed into
/// <see cref="Structure"/>, the Data API (REST <c>core/entry</c>, <c>core/polymer_entity</c> and GraphQL) and the
/// Search API (full text, sequence, attribute queries).
/// </summary>
public sealed class PdbClient
{
    private static readonly TimeSpan FileTtl = TimeSpan.FromDays(30);
    private static readonly TimeSpan DataTtl = TimeSpan.FromDays(7);
    private static readonly TimeSpan SearchTtl = TimeSpan.FromDays(1);
    private readonly DataHttp _http;

    /// <param name="httpClient">An <see cref="HttpClient"/> to use (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public PdbClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        _http = new DataHttp("RCSB PDB", httpClient, options);
    }

    /// <summary>Coordinate file server.</summary>
    public string FilesUrl { get; init; } = "https://files.rcsb.org/download/";

    /// <summary>Data API server.</summary>
    public string DataUrl { get; init; } = "https://data.rcsb.org/";

    /// <summary>Search API endpoint.</summary>
    public string SearchUrl { get; init; } = "https://search.rcsb.org/rcsbsearch/v2/query";

    // ------------------------------------------------------------------ coordinates

    /// <summary>Downloads a structure (gzip-compressed on the wire) and parses it.</summary>
    /// <remarks>Very large entries (more than 99,999 atoms or 62 chains) exist only as mmCIF; asking for PDB format throws a 404 <see cref="DataServiceException"/>.</remarks>
    public async Task<Structure> GetStructureAsync(string pdbId, StructureFileFormat format = StructureFileFormat.MmCif,
        CancellationToken cancellationToken = default)
    {
        string text = await GetStructureTextAsync(pdbId, format, cancellationToken).ConfigureAwait(false);
        string id = NormalizeId(pdbId);
        return format == StructureFileFormat.Pdb ? PdbReader.Parse(text, id) : MmCifReader.Parse(text, id);
    }

    /// <summary>Downloads a coordinate file and returns its (decompressed) text.</summary>
    public async Task<string> GetStructureTextAsync(string pdbId, StructureFileFormat format = StructureFileFormat.MmCif,
        CancellationToken cancellationToken = default)
    {
        string id = NormalizeId(pdbId);
        string ext = format == StructureFileFormat.Pdb ? "pdb" : "cif";
        var url = new Uri($"{FilesUrl}{id}.{ext}.gz");
        try
        {
            var response = await SendAsync(DataRequest.Get(url), FileTtl, cancellationToken).ConfigureAwait(false);
            return response.Text();
        }
        catch (DataServiceException ex) when (ex.IsNotFound && format == StructureFileFormat.Pdb)
        {
            throw new DataServiceException(ex.Service, $"{ex.Message} (Large entries have no PDB-format file; request StructureFileFormat.MmCif.)",
                ex.RequestUri, ex.StatusCode, ex.ServiceMessage, ex);
        }
    }

    // ------------------------------------------------------------------ Data API

    /// <summary>Entry-level metadata: title, method, resolution, dates, citation.</summary>
    public async Task<PdbEntryInfo> GetEntryAsync(string pdbId, CancellationToken cancellationToken = default)
    {
        var response = await SendAsync(DataRequest.Get(new Uri($"{DataUrl}rest/v1/core/entry/{NormalizeId(pdbId)}")), DataTtl, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        return PdbEntryInfo.FromJson(doc.RootElement);
    }

    /// <summary>One polymer entity: sequence, chains, UniProt mapping, source organism.</summary>
    public async Task<PdbPolymerEntity> GetPolymerEntityAsync(string pdbId, string entityId = "1", CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(entityId);
        var url = new Uri($"{DataUrl}rest/v1/core/polymer_entity/{NormalizeId(pdbId)}/{Uri.EscapeDataString(entityId)}");
        var response = await SendAsync(DataRequest.Get(url), DataTtl, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        return PdbPolymerEntity.FromJson(doc.RootElement);
    }

    /// <summary>
    /// Runs a GraphQL query against <c>data.rcsb.org/graphql</c> and returns its <c>data</c> member (cloned).
    /// GraphQL errors become a <see cref="DataServiceException"/>.
    /// </summary>
    public async Task<JsonElement> GraphQlAsync(string query, JsonObject? variables = null, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(query);
        var body = new JsonObject { ["query"] = query };
        if (variables is not null) body["variables"] = variables.DeepClone();
        var response = await SendAsync(DataRequest.PostJson(new Uri($"{DataUrl}graphql"), body.ToJsonString()), DataTtl, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        if (doc.RootElement.Prop("errors") is { ValueKind: JsonValueKind.Array } errors && errors.GetArrayLength() > 0)
        {
            string msg = string.Join("; ", errors.EnumerateArray().Select(e => e.Str("message") ?? e.ToString()));
            throw new DataServiceException("RCSB PDB", $"GraphQL query failed: {msg}", serviceMessage: msg);
        }

        return doc.RootElement.Prop("data")?.Clone() ?? default;
    }

    // ------------------------------------------------------------------ Search API

    /// <summary>Runs a Search API query.</summary>
    public Task<PdbSearchResult> SearchAsync(PdbSearchQuery query, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(query);
        return SearchAsync(query.ToJson(), cancellationToken);
    }

    /// <summary>Runs a Search API query given as raw JSON.</summary>
    public async Task<PdbSearchResult> SearchAsync(string queryJson, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(queryJson);
        var response = await SendAsync(DataRequest.PostJson(new Uri(SearchUrl), queryJson), SearchTtl, cancellationToken).ConfigureAwait(false);
        if (response.Content.Length == 0) return new PdbSearchResult(0, "", []); // 204 No Content: no hits
        using var doc = response.Json();
        var root = doc.RootElement;
        return new PdbSearchResult(
            root.Long("total_count") ?? 0,
            root.Str("result_type") ?? "",
            root.Arr("result_set").Select(h => new PdbSearchHit(h.Str("identifier") ?? "", h.Double("score") ?? 0)).ToList());
    }

    private Task<DataResponse> SendAsync(DataRequest request, TimeSpan ttl, CancellationToken ct) =>
        _http.SendAsync(request with { RequestsPerSecond = 10, CacheTtl = ttl }, ct);

    internal static string NormalizeId(string pdbId)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(pdbId);
        string id = pdbId.Trim().ToUpperInvariant();
        if (id.Length < 4 || !id.All(char.IsAsciiLetterOrDigit))
            throw new ArgumentException($"'{pdbId}' is not a PDB id.", nameof(pdbId));
        return id;
    }
}

/// <summary>
/// AlphaFold Protein Structure Database (<c>alphafold.ebi.ac.uk</c>): prediction metadata by UniProt accession and
/// model downloads parsed into <see cref="Structure"/> (per-residue pLDDT is in each atom's B-factor).
/// </summary>
public sealed class AlphaFoldClient
{
    private static readonly TimeSpan Ttl = TimeSpan.FromDays(30);
    private readonly DataHttp _http;

    /// <param name="httpClient">An <see cref="HttpClient"/> to use (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public AlphaFoldClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        _http = new DataHttp("AlphaFold DB", httpClient, options);
    }

    /// <summary>The API base URL.</summary>
    public string BaseUrl { get; init; } = "https://alphafold.ebi.ac.uk/api/";

    /// <summary>All predictions for a UniProt accession (usually one; long proteins may have fragments).</summary>
    public async Task<IReadOnlyList<AlphaFoldPrediction>> GetPredictionsAsync(string uniProtAccession, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(uniProtAccession);
        var url = new Uri($"{BaseUrl}prediction/{Uri.EscapeDataString(uniProtAccession.Trim())}");
        var response = await _http.SendAsync(DataRequest.Get(url) with { RequestsPerSecond = 10, CacheTtl = Ttl, Accept = "application/json" }, cancellationToken).ConfigureAwait(false);
        using var doc = response.Json();
        var root = doc.RootElement;
        var items = root.ValueKind == JsonValueKind.Array ? root.EnumerateArray().ToList() : [root];
        return items.Select(AlphaFoldPrediction.FromJson).ToList();
    }

    /// <summary>Downloads the (first) predicted model for a UniProt accession.</summary>
    public async Task<Structure> GetStructureAsync(string uniProtAccession, StructureFileFormat format = StructureFileFormat.Pdb,
        CancellationToken cancellationToken = default)
    {
        var predictions = await GetPredictionsAsync(uniProtAccession, cancellationToken).ConfigureAwait(false);
        var p = predictions.FirstOrDefault() ?? throw new DataServiceException("AlphaFold DB", $"No AlphaFold prediction for {uniProtAccession}.");
        return await GetStructureAsync(p, format, cancellationToken).ConfigureAwait(false);
    }

    /// <summary>Downloads the model of one prediction.</summary>
    public async Task<Structure> GetStructureAsync(AlphaFoldPrediction prediction, StructureFileFormat format = StructureFileFormat.Pdb,
        CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(prediction);
        var url = (format == StructureFileFormat.Pdb ? prediction.PdbUrl : prediction.CifUrl)
                  ?? throw new DataServiceException("AlphaFold DB", $"Prediction {prediction.EntryId} has no {format} file URL.");
        var response = await _http.SendAsync(DataRequest.Get(url) with { RequestsPerSecond = 10, CacheTtl = Ttl }, cancellationToken).ConfigureAwait(false);
        string text = response.Text();
        return format == StructureFileFormat.Pdb ? PdbReader.Parse(text, prediction.EntryId) : MmCifReader.Parse(text, prediction.EntryId);
    }
}
