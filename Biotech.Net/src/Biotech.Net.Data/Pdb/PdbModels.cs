using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.Pdb;

/// <summary>Coordinate file formats.</summary>
public enum StructureFileFormat
{
    /// <summary>Legacy PDB format (not available for very large entries).</summary>
    Pdb,

    /// <summary>PDBx/mmCIF, the archive's master format.</summary>
    MmCif,
}

/// <summary>Summary of a PDB entry from the RCSB Data API (<c>/rest/v1/core/entry</c>).</summary>
public sealed record PdbEntryInfo
{
    public required string Id { get; init; }
    public string Title { get; init; } = "";
    public IReadOnlyList<string> ExperimentalMethods { get; init; } = [];

    /// <summary>Best resolution in Å, when the method reports one.</summary>
    public double? Resolution { get; init; }

    public DateOnly? DepositionDate { get; init; }
    public DateOnly? ReleaseDate { get; init; }
    public DateOnly? RevisionDate { get; init; }
    public string? Keywords { get; init; }
    public int PolymerEntityCount { get; init; }
    public int DepositedAtomCount { get; init; }
    public int DepositedModelCount { get; init; }

    /// <summary>Molecular weight of the deposited model in kDa.</summary>
    public double? MolecularWeightKDa { get; init; }

    public string? CitationTitle { get; init; }
    public string? CitationJournal { get; init; }
    public int? CitationYear { get; init; }
    public string? CitationDoi { get; init; }
    public long? CitationPubMedId { get; init; }
    public IReadOnlyList<string> CitationAuthors { get; init; } = [];

    internal static PdbEntryInfo FromJson(JsonElement e)
    {
        var info = e.Prop("rcsb_entry_info");
        var acc = e.Prop("rcsb_accession_info");
        var cit = e.Prop("rcsb_primary_citation");
        double? resolution = info?.Arr("resolution_combined").Select(Json.AsDouble).FirstOrDefault(d => d is not null)
                             ?? e.Arr("refine").Select(r => r.Double("ls_d_res_high")).FirstOrDefault(d => d is not null);
        return new PdbEntryInfo
        {
            Id = e.Str("rcsb_id") ?? e.Path("entry")?.Str("id") ?? "",
            Title = e.Path("struct")?.Str("title") ?? "",
            ExperimentalMethods = e.Arr("exptl").Select(x => x.Str("method")).OfType<string>().ToList(),
            Resolution = resolution,
            DepositionDate = Date(acc?.Str("deposit_date")),
            ReleaseDate = Date(acc?.Str("initial_release_date")),
            RevisionDate = Date(acc?.Str("revision_date")),
            Keywords = e.Path("struct_keywords")?.Str("pdbx_keywords"),
            PolymerEntityCount = info?.Int("polymer_entity_count") ?? 0,
            DepositedAtomCount = info?.Int("deposited_atom_count") ?? 0,
            DepositedModelCount = info?.Int("deposited_model_count") ?? 0,
            MolecularWeightKDa = info?.Double("molecular_weight"),
            CitationTitle = cit?.Str("title"),
            CitationJournal = cit?.Str("rcsb_journal_abbrev") ?? cit?.Str("journal_abbrev"),
            CitationYear = cit?.Int("year"),
            CitationDoi = cit?.Str("pdbx_database_id_DOI"),
            CitationPubMedId = cit?.Long("pdbx_database_id_PubMed"),
            CitationAuthors = cit?.Strings("rcsb_authors") ?? [],
        };
    }

    private static DateOnly? Date(string? iso) =>
        iso is { Length: >= 10 } && DateOnly.TryParseExact(iso[..10], "yyyy-MM-dd", CultureInfo.InvariantCulture, DateTimeStyles.None, out var d) ? d : null;
}

/// <summary>A polymer entity (one distinct chain sequence) of a PDB entry.</summary>
public sealed record PdbPolymerEntity
{
    public required string EntryId { get; init; }
    public required string EntityId { get; init; }
    public string Description { get; init; } = "";

    /// <summary>Canonical one-letter sequence.</summary>
    public string SequenceResidues { get; init; } = "";

    /// <summary><c>Protein</c>, <c>DNA</c>, <c>RNA</c>, ...</summary>
    public string? PolymerType { get; init; }

    /// <summary>Author chain ids carrying this entity.</summary>
    public IReadOnlyList<string> ChainIds { get; init; } = [];

    public IReadOnlyList<string> UniProtAccessions { get; init; } = [];
    public IReadOnlyList<string> Organisms { get; init; } = [];
    public IReadOnlyList<long> TaxonomyIds { get; init; } = [];

    /// <summary>Formula weight in kDa.</summary>
    public double? FormulaWeightKDa { get; init; }

    /// <summary>The entity sequence (id <c>1CRN_1</c>).</summary>
    public Sequence ToSequence()
    {
        var alphabet = PolymerType switch
        {
            "DNA" => Alphabet.DnaIupac,
            "RNA" => Alphabet.RnaIupac,
            _ => Alphabet.ProteinExtended,
        };
        return new Sequence($"{EntryId}_{EntityId}", System.Text.Encoding.ASCII.GetBytes(SequenceResidues), alphabet, Description, validate: false);
    }

    internal static PdbPolymerEntity FromJson(JsonElement e)
    {
        var ids = e.Prop("rcsb_polymer_entity_container_identifiers");
        var poly = e.Prop("entity_poly");
        var orgs = e.Arr("rcsb_entity_source_organism").ToList();
        return new PdbPolymerEntity
        {
            EntryId = ids?.Str("entry_id") ?? "",
            EntityId = ids?.Str("entity_id") ?? "",
            Description = e.Path("rcsb_polymer_entity")?.Str("pdbx_description") ?? "",
            SequenceResidues = (poly?.Str("pdbx_seq_one_letter_code_can") ?? "").Replace("\n", "", StringComparison.Ordinal),
            PolymerType = poly?.Str("rcsb_entity_polymer_type"),
            ChainIds = ids?.Strings("auth_asym_ids") ?? [],
            UniProtAccessions = ids?.Strings("uniprot_ids") ?? [],
            Organisms = orgs.Select(o => o.Str("scientific_name") ?? o.Str("ncbi_scientific_name")).OfType<string>().Distinct().ToList(),
            TaxonomyIds = orgs.Select(o => o.Long("ncbi_taxonomy_id")).OfType<long>().Distinct().ToList(),
            FormulaWeightKDa = e.Path("rcsb_polymer_entity")?.Double("formula_weight"),
        };
    }
}

/// <summary>One hit of an RCSB Search API query.</summary>
public sealed record PdbSearchHit(string Identifier, double Score);

/// <summary>The result of an RCSB Search API query.</summary>
public sealed record PdbSearchResult(long TotalCount, string ResultType, IReadOnlyList<PdbSearchHit> Hits)
{
    /// <summary>The hit identifiers (entry ids, or <c>ENTRY_ENTITY</c> ids for polymer-entity searches).</summary>
    public IEnumerable<string> Identifiers => Hits.Select(h => h.Identifier);
}

/// <summary>Sequence types understood by the RCSB sequence search.</summary>
public enum PdbSequenceType
{
    Protein,
    Dna,
    Rna,
}

/// <summary>
/// Builds RCSB Search API v2 queries (<c>search.rcsb.org/rcsbsearch/v2/query</c>). For anything the builders do not
/// cover, write the JSON yourself and call <see cref="PdbClient.SearchAsync(string, CancellationToken)"/>.
/// </summary>
public sealed class PdbSearchQuery
{
    private readonly JsonObject _node;

    private PdbSearchQuery(JsonObject node) => _node = node;

    /// <summary><c>entry</c> (default), <c>polymer_entity</c>, <c>polymer_instance</c>, <c>assembly</c>, <c>non_polymer_entity</c>.</summary>
    public string ReturnType { get; init; } = "entry";

    /// <summary>Offset of the first hit.</summary>
    public int Start { get; init; }

    /// <summary>Hits per request (the API's default is 10).</summary>
    public int Rows { get; init; } = 10;

    /// <summary>Also search computed structure models (AlphaFold DB, ModelArchive).</summary>
    public bool IncludeComputedModels { get; init; }

    /// <summary>Free-text search over the whole archive.</summary>
    public static PdbSearchQuery FullText(string text)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(text);
        return new(Terminal("full_text", new JsonObject { ["value"] = text }));
    }

    /// <summary>Sequence similarity search (MMseqs2). Returns polymer entities (<c>1CRN_1</c>), since sequences belong to entities.</summary>
    /// <param name="residues">Query sequence (one-letter codes).</param>
    /// <param name="type">Protein, DNA or RNA.</param>
    /// <param name="evalueCutoff">Largest E-value kept.</param>
    /// <param name="identityCutoff">Smallest sequence identity kept (0-1).</param>
    public static PdbSearchQuery Sequence(string residues, PdbSequenceType type = PdbSequenceType.Protein,
        double evalueCutoff = 0.1, double identityCutoff = 0)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(residues);
        var p = new JsonObject
        {
            ["value"] = residues,
            ["sequence_type"] = type switch { PdbSequenceType.Dna => "dna", PdbSequenceType.Rna => "rna", _ => "protein" },
            ["evalue_cutoff"] = evalueCutoff,
            ["identity_cutoff"] = identityCutoff,
        };
        return new(Terminal("sequence", p)) { ReturnType = "polymer_entity" };
    }

    /// <summary>Attribute search, e.g. <c>Attribute("rcsb_entry_info.resolution_combined", "less", 1.0)</c>.</summary>
    /// <param name="attribute">A schema attribute path.</param>
    /// <param name="operator">E.g. <c>exact_match</c>, <c>equals</c>, <c>less</c>, <c>greater</c>, <c>contains_words</c>, <c>in</c>.</param>
    /// <param name="value">A string, number, bool, or a collection for <c>in</c>.</param>
    public static PdbSearchQuery Attribute(string attribute, string @operator, object value)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(attribute);
        var p = new JsonObject { ["attribute"] = attribute, ["operator"] = @operator, ["value"] = JsonSerializer.SerializeToNode(value) };
        return new(Terminal("text", p));
    }

    /// <summary>All of the queries must match (their paging options are ignored; set them on the result).</summary>
    public static PdbSearchQuery And(params PdbSearchQuery[] queries) => Group("and", queries);

    /// <summary>Any of the queries may match.</summary>
    public static PdbSearchQuery Or(params PdbSearchQuery[] queries) => Group("or", queries);

    /// <summary>A copy with different paging/return options.</summary>
    public PdbSearchQuery With(string? returnType = null, int? start = null, int? rows = null, bool? includeComputedModels = null) =>
        new((JsonObject)_node.DeepClone())
        {
            ReturnType = returnType ?? ReturnType,
            Start = start ?? Start,
            Rows = rows ?? Rows,
            IncludeComputedModels = includeComputedModels ?? IncludeComputedModels,
        };

    /// <summary>The request JSON.</summary>
    public string ToJson()
    {
        var content = new JsonArray("experimental");
        if (IncludeComputedModels) content.Add("computational");
        var request = new JsonObject
        {
            ["query"] = _node.DeepClone(),
            ["return_type"] = ReturnType,
            ["request_options"] = new JsonObject
            {
                ["paginate"] = new JsonObject { ["start"] = Start, ["rows"] = Rows },
                ["results_content_type"] = content,
            },
        };
        return request.ToJsonString();
    }

    public override string ToString() => ToJson();

    private static JsonObject Terminal(string service, JsonObject parameters) =>
        new() { ["type"] = "terminal", ["service"] = service, ["parameters"] = parameters };

    private static PdbSearchQuery Group(string op, PdbSearchQuery[] queries)
    {
        if (queries.Length == 0) throw new ArgumentException("At least one query is required.", nameof(queries));
        var nodes = new JsonArray();
        foreach (var q in queries) nodes.Add(q._node.DeepClone());
        return new(new JsonObject { ["type"] = "group", ["logical_operator"] = op, ["nodes"] = nodes });
    }
}

/// <summary>An AlphaFold DB prediction for one UniProt accession.</summary>
public sealed record AlphaFoldPrediction
{
    /// <summary>E.g. <c>AF-P69905-F1</c>.</summary>
    public required string EntryId { get; init; }

    public required string UniProtAccession { get; init; }
    public string? UniProtId { get; init; }
    public string? Description { get; init; }
    public string? Gene { get; init; }
    public string? OrganismScientificName { get; init; }
    public long? TaxonId { get; init; }
    public string SequenceResidues { get; init; } = "";
    public int SequenceStart { get; init; }
    public int SequenceEnd { get; init; }
    public int? LatestVersion { get; init; }
    public DateTimeOffset? ModelCreatedDate { get; init; }

    /// <summary>Mean pLDDT over the model (0-100).</summary>
    public double? MeanPlddt { get; init; }

    public Uri? PdbUrl { get; init; }
    public Uri? CifUrl { get; init; }
    public Uri? PaeUrl { get; init; }

    internal static AlphaFoldPrediction FromJson(JsonElement e) => new()
    {
        EntryId = e.Str("entryId") ?? e.Str("modelEntityId") ?? "",
        UniProtAccession = e.Str("uniprotAccession") ?? "",
        UniProtId = e.Str("uniprotId"),
        Description = e.Str("uniprotDescription"),
        Gene = e.Str("gene"),
        OrganismScientificName = e.Str("organismScientificName"),
        TaxonId = e.Long("taxId"),
        SequenceResidues = e.Str("uniprotSequence") ?? e.Str("sequence") ?? "",
        SequenceStart = e.Int("uniprotStart") ?? e.Int("sequenceStart") ?? 1,
        SequenceEnd = e.Int("uniprotEnd") ?? e.Int("sequenceEnd") ?? 0,
        LatestVersion = e.Int("latestVersion"),
        ModelCreatedDate = DateTimeOffset.TryParse(e.Str("modelCreatedDate"), CultureInfo.InvariantCulture, DateTimeStyles.AssumeUniversal, out var d) ? d : null,
        MeanPlddt = e.Double("globalMetricValue"),
        PdbUrl = Url(e.Str("pdbUrl")),
        CifUrl = Url(e.Str("cifUrl")),
        PaeUrl = Url(e.Str("paeDocUrl")),
    };

    private static Uri? Url(string? s) => Uri.TryCreate(s, UriKind.Absolute, out var u) ? u : null;
}
