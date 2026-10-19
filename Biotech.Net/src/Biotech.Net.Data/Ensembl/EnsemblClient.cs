using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.Ensembl;

/// <summary>
/// Ensembl REST API (<c>rest.ensembl.org</c>): lookup by id or symbol, sequences, region overlap, the Variant
/// Effect Predictor and cross-references.
/// </summary>
/// <remarks>
/// Ensembl allows 15 requests per second and reports its quota in <c>X-RateLimit-*</c> headers; the client spaces
/// requests at that rate and pauses the host when the quota is exhausted or a 429 carries <c>Retry-After</c>.
/// Set <see cref="Server"/> to <see cref="Grch37Server"/> for the GRCh37 archive.
/// </remarks>
public sealed class EnsemblClient
{
    /// <summary>The current-release server (GRCh38 for human).</summary>
    public const string DefaultServer = "https://rest.ensembl.org/";

    /// <summary>The GRCh37 archive server.</summary>
    public const string Grch37Server = "https://grch37.rest.ensembl.org/";

    private const int VepBatchLimit = 200;
    private static readonly TimeSpan Ttl = TimeSpan.FromDays(7);
    private readonly DataHttp _http;

    /// <param name="httpClient">An <see cref="HttpClient"/> to use (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public EnsemblClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        _http = new DataHttp("Ensembl", httpClient, options);
    }

    /// <summary>The REST server (ending in <c>/</c>).</summary>
    public string Server { get; init; } = DefaultServer;

    // ------------------------------------------------------------------ lookup

    /// <summary>Looks up a stable id (gene, transcript, exon or translation).</summary>
    /// <param name="expand">Include transcripts, exons and translations.</param>
    public async Task<EnsemblFeature> LookupIdAsync(string id, bool expand = false, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(id);
        var q = Json().Add("expand", expand ? "1" : null);
        using var doc = await GetJsonAsync($"lookup/id/{Escape(id)}", q, cancellationToken).ConfigureAwait(false);
        return EnsemblFeature.FromJson(doc.RootElement);
    }

    /// <summary>Looks up a gene by symbol, e.g. <c>("homo_sapiens", "BRCA2")</c>.</summary>
    public async Task<EnsemblFeature> LookupSymbolAsync(string species, string symbol, bool expand = false, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        ArgumentException.ThrowIfNullOrWhiteSpace(symbol);
        var q = Json().Add("expand", expand ? "1" : null);
        using var doc = await GetJsonAsync($"lookup/symbol/{Escape(species)}/{Escape(symbol)}", q, cancellationToken).ConfigureAwait(false);
        return EnsemblFeature.FromJson(doc.RootElement);
    }

    // ------------------------------------------------------------------ sequence

    /// <summary>The sequence of a feature. A gene asked for cDNA/CDS/protein has several; use <see cref="GetSequencesAsync"/>.</summary>
    public async Task<Sequence> GetSequenceAsync(string id, EnsemblSequenceType type = EnsemblSequenceType.Genomic, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(id);
        using var doc = await GetJsonAsync($"sequence/id/{Escape(id)}", Json().Add("type", TypeName(type)), cancellationToken).ConfigureAwait(false);
        return ToSequence(doc.RootElement, type);
    }

    /// <summary>Every sequence of a feature of the given type (e.g. all transcripts' cDNA of a gene).</summary>
    public async Task<IReadOnlyList<Sequence>> GetSequencesAsync(string id, EnsemblSequenceType type, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(id);
        var q = Json().Add("type", TypeName(type)).Add("multiple_sequences", "1");
        using var doc = await GetJsonAsync($"sequence/id/{Escape(id)}", q, cancellationToken).ConfigureAwait(false);
        var root = doc.RootElement;
        return root.ValueKind == JsonValueKind.Array ? root.EnumerateArray().Select(e => ToSequence(e, type)).ToList() : [ToSequence(root, type)];
    }

    /// <summary>
    /// The reference sequence of a region. <paramref name="region"/> is 0-based half-open (Biotech.Net convention);
    /// a <see cref="Strand.Reverse"/> region returns the reverse complement.
    /// </summary>
    public async Task<Sequence> GetRegionSequenceAsync(string species, GenomicInterval region, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        string r = $"{region.Chromosome}:{region.Start + 1}..{region.End}:{(region.Strand == Strand.Reverse ? -1 : 1)}";
        using var doc = await GetJsonAsync($"sequence/region/{Escape(species)}/{Escape(r)}", Json(), cancellationToken).ConfigureAwait(false);
        return ToSequence(doc.RootElement, EnsemblSequenceType.Genomic);
    }

    // ------------------------------------------------------------------ overlap

    /// <summary>
    /// Features overlapping a region (0-based half-open <paramref name="region"/>; Ensembl caps regions at 5 Mb).
    /// </summary>
    /// <param name="features">Feature types: <c>gene</c> (default), <c>transcript</c>, <c>exon</c>, <c>cds</c>, <c>variation</c>, <c>regulatory</c>, ...</param>
    public async Task<IReadOnlyList<EnsemblOverlapFeature>> OverlapRegionAsync(string species, GenomicInterval region, IEnumerable<string>? features = null,
        CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        var q = new QueryString();
        foreach (string f in features ?? ["gene"]) q.Add("feature", f);
        q.Add("content-type", "application/json");
        string r = $"{region.Chromosome}:{region.Start + 1}-{region.End}";
        using var doc = await GetJsonAsync($"overlap/region/{Escape(species)}/{Escape(r)}", q, cancellationToken).ConfigureAwait(false);
        return doc.RootElement.ValueKind == JsonValueKind.Array ? doc.RootElement.EnumerateArray().Select(EnsemblOverlapFeature.FromJson).ToList() : [];
    }

    // ------------------------------------------------------------------ VEP

    /// <summary>Variant effects for an HGVS notation, e.g. <c>ENST00000335295:c.20A&gt;T</c> or <c>9:g.22125504G&gt;C</c>.</summary>
    public async Task<IReadOnlyList<VepResult>> VepHgvsAsync(string species, string hgvs, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        ArgumentException.ThrowIfNullOrWhiteSpace(hgvs);
        using var doc = await GetJsonAsync($"vep/{Escape(species)}/hgvs/{Escape(hgvs)}", Json(), cancellationToken).ConfigureAwait(false);
        return ParseVep(doc);
    }

    /// <summary>
    /// Variant effects for a VCF-style variant (1-based position, first alternate allele). Anchor bases of indels are
    /// handled as VEP expects.
    /// </summary>
    public async Task<IReadOnlyList<VepResult>> VepRegionAsync(string species, Variant variant, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        ArgumentNullException.ThrowIfNull(variant);
        if (variant.Alternates.Count == 0) throw new ArgumentException("The variant has no alternate allele.", nameof(variant));
        string region = ToVepRegion(variant.Chromosome, variant.Position, variant.Reference, variant.Alternates[0]);
        int slash = region.LastIndexOf('/');
        string path = $"vep/{Escape(species)}/region/{Escape(region[..slash])}/{Escape(region[(slash + 1)..])}";
        using var doc = await GetJsonAsync(path, Json(), cancellationToken).ConfigureAwait(false);
        return ParseVep(doc);
    }

    /// <summary>Variant effects for many variants, POSTed in batches of 200 VCF lines.</summary>
    public async Task<IReadOnlyList<VepResult>> VepRegionsAsync(string species, IEnumerable<Variant> variants, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        ArgumentNullException.ThrowIfNull(variants);
        var lines = variants.Select(v => string.Join(' ', v.Chromosome, v.Position.ToString(CultureInfo.InvariantCulture),
            v.Id ?? ".", v.Reference, v.Alternates.Count == 0 ? "." : string.Join(',', v.Alternates), ".", ".", ".")).ToList();
        var results = new List<VepResult>();
        foreach (var chunk in lines.Chunk(VepBatchLimit))
        {
            var body = new JsonObject { ["variants"] = new JsonArray(chunk.Select(l => (JsonNode)JsonValue.Create(l)).ToArray()) };
            var url = Json().Apply($"{Server}vep/{Escape(species)}/region");
            var response = await _http.SendAsync(DataRequest.PostJson(url, body.ToJsonString()) with { RequestsPerSecond = 15, CacheTtl = Ttl }, cancellationToken).ConfigureAwait(false);
            using var doc = response.Json();
            results.AddRange(ParseVep(doc));
        }

        return results;
    }

    /// <summary>The Ensembl VEP region notation (<c>chr:start-end:strand/allele</c>) for a VCF-style variant.</summary>
    public static string ToVepRegion(string chromosome, long position, string reference, string alternate)
    {
        // VCF writes indels with a shared leading anchor base; VEP region notation omits it.
        if (reference.Length != alternate.Length && reference.Length > 0 && alternate.Length > 0 && reference[0] == alternate[0])
        {
            reference = reference[1..];
            alternate = alternate[1..];
            position++;
        }

        long end = position + reference.Length - 1; // an insertion gives end = start - 1, as VEP expects
        string allele = alternate.Length == 0 ? "-" : alternate;
        return $"{chromosome}:{position}-{end}:1/{allele}";
    }

    // ------------------------------------------------------------------ xrefs

    /// <summary>External references of an Ensembl object.</summary>
    /// <param name="externalDb">Restrict to one database, e.g. <c>HGNC</c>, <c>Uniprot/SWISSPROT</c>, <c>EntrezGene</c>, <c>RefSeq_mRNA</c>.</param>
    public async Task<IReadOnlyList<EnsemblXref>> GetXrefsAsync(string id, string? externalDb = null, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(id);
        using var doc = await GetJsonAsync($"xrefs/id/{Escape(id)}", Json().Add("external_db", externalDb), cancellationToken).ConfigureAwait(false);
        return doc.RootElement.ValueKind == JsonValueKind.Array ? doc.RootElement.EnumerateArray().Select(EnsemblXref.FromJson).ToList() : [];
    }

    /// <summary>Ensembl objects linked to a symbol or external name.</summary>
    public async Task<IReadOnlyList<EnsemblIdentity>> XrefsBySymbolAsync(string species, string symbol, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(species);
        ArgumentException.ThrowIfNullOrWhiteSpace(symbol);
        using var doc = await GetJsonAsync($"xrefs/symbol/{Escape(species)}/{Escape(symbol)}", Json(), cancellationToken).ConfigureAwait(false);
        return doc.RootElement.ValueKind == JsonValueKind.Array
            ? doc.RootElement.EnumerateArray().Select(e => new EnsemblIdentity(e.Str("id") ?? "", e.Str("type") ?? "")).ToList()
            : [];
    }

    // ------------------------------------------------------------------ plumbing

    private static QueryString Json() => new QueryString().Add("content-type", "application/json");

    /// <summary>Escapes one path segment, keeping <c>:</c> readable (Ensembl regions and HGVS use it).</summary>
    private static string Escape(string s) => Uri.EscapeDataString(s.Trim()).Replace("%3A", ":", StringComparison.Ordinal);

    private static string TypeName(EnsemblSequenceType type) => type switch
    {
        EnsemblSequenceType.Cdna => "cdna",
        EnsemblSequenceType.Cds => "cds",
        EnsemblSequenceType.Protein => "protein",
        _ => "genomic",
    };

    private static Sequence ToSequence(JsonElement e, EnsemblSequenceType type)
    {
        string residues = e.Str("seq") ?? throw new DataServiceException("Ensembl", "Sequence response has no 'seq'.");
        string id = e.Str("id") ?? e.Str("query") ?? "ensembl";
        var alphabet = type == EnsemblSequenceType.Protein || e.Str("molecule") == "protein" ? Alphabet.ProteinExtended : Alphabet.DnaIupac;
        var seq = new Sequence(id, Encoding.ASCII.GetBytes(residues), alphabet, e.Str("desc") ?? e.Str("query"), validate: false);
        if (e.Str("query") is { } query) seq.Annotations["query"] = query;
        if (e.Str("version") is { } version) seq.Annotations["version"] = version;
        return seq;
    }

    private static IReadOnlyList<VepResult> ParseVep(JsonDocument doc) =>
        doc.RootElement.ValueKind == JsonValueKind.Array ? doc.RootElement.EnumerateArray().Select(VepResult.FromJson).ToList() : [];

    private async Task<JsonDocument> GetJsonAsync(string path, QueryString q, CancellationToken ct)
    {
        var request = DataRequest.Get(q.Apply(Server + path)) with { RequestsPerSecond = 15, CacheTtl = Ttl, Accept = "application/json" };
        var response = await _http.SendAsync(request, ct).ConfigureAwait(false);
        return response.Json();
    }
}
