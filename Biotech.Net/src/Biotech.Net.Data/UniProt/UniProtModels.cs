using System.Text.Json;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.UniProt;

/// <summary>A UniProtKB comment block (FUNCTION, SUBCELLULAR LOCATION, ...), flattened to text.</summary>
public sealed record UniProtComment(string Type, string Text);

/// <summary>A sequence feature (Chain, Domain, Disulfide bond, ...). Positions are 1-based and inclusive, as in UniProt.</summary>
public sealed record UniProtFeature(string Type, string Description, int? Start, int? End, string? FeatureId);

/// <summary>A UniProtKB keyword.</summary>
public sealed record UniProtKeyword(string Id, string Name, string Category);

/// <summary>A cross-reference to another database (PDB, Ensembl, RefSeq, GO, ...).</summary>
public sealed record UniProtCrossReference(string Database, string Id, IReadOnlyDictionary<string, string> Properties);

/// <summary>
/// A UniProtKB entry, parsed from the REST JSON format. Searches restricted with <c>fields</c> fill only the
/// requested members; the rest are empty.
/// </summary>
public sealed class UniProtEntry
{
    /// <summary>Primary accession, e.g. <c>P69905</c>.</summary>
    public required string Accession { get; init; }

    public IReadOnlyList<string> SecondaryAccessions { get; init; } = [];

    /// <summary>Entry name (mnemonic), e.g. <c>HBA_HUMAN</c>.</summary>
    public string? EntryName { get; init; }

    /// <summary>True for Swiss-Prot (reviewed), false for TrEMBL.</summary>
    public bool IsReviewed { get; init; }

    /// <summary>Recommended (or submitted) full protein name.</summary>
    public string? ProteinName { get; init; }

    public IReadOnlyList<string> AlternativeNames { get; init; } = [];

    /// <summary>Gene names: the primary name first, then synonyms and ordered-locus names.</summary>
    public IReadOnlyList<string> GeneNames { get; init; } = [];

    public string? OrganismName { get; init; }
    public string? OrganismCommonName { get; init; }
    public long? TaxonId { get; init; }
    public IReadOnlyList<string> Lineage { get; init; } = [];

    /// <summary>Protein existence level, e.g. <c>1: Evidence at protein level</c>.</summary>
    public string? ProteinExistence { get; init; }

    /// <summary>Annotation score (1-5).</summary>
    public double? AnnotationScore { get; init; }

    public int? EntryVersion { get; init; }
    public int? SequenceVersion { get; init; }

    /// <summary>The canonical sequence (one-letter codes).</summary>
    public string SequenceResidues { get; init; } = "";

    public int Length { get; init; }

    /// <summary>Average molecular weight in daltons, as UniProt computes it.</summary>
    public double? MolecularWeight { get; init; }

    public string? Crc64 { get; init; }

    public IReadOnlyList<UniProtComment> Comments { get; init; } = [];
    public IReadOnlyList<UniProtFeature> Features { get; init; } = [];
    public IReadOnlyList<UniProtKeyword> Keywords { get; init; } = [];
    public IReadOnlyList<UniProtCrossReference> CrossReferences { get; init; } = [];

    /// <summary>The FUNCTION comment, when present.</summary>
    public string? Function => Comments.FirstOrDefault(c => c.Type == "FUNCTION")?.Text;

    /// <summary>Cross-references to one database, e.g. <c>"PDB"</c> or <c>"Ensembl"</c>.</summary>
    public IEnumerable<UniProtCrossReference> CrossReferencesTo(string database) =>
        CrossReferences.Where(x => x.Database.Equals(database, StringComparison.OrdinalIgnoreCase));

    /// <summary>The canonical sequence as a protein <see cref="Sequence"/> (id = accession, description = name).</summary>
    public Sequence ToSequence()
    {
        if (SequenceResidues.Length == 0) throw new InvalidOperationException($"Entry {Accession} carries no sequence (was 'sequence' left out of the fields?).");
        string description = $"{ProteinName}{(OrganismName is null ? "" : $" OS={OrganismName}")}{(GeneNames.Count > 0 ? $" GN={GeneNames[0]}" : "")}".Trim();
        var seq = new Sequence(Accession, System.Text.Encoding.ASCII.GetBytes(SequenceResidues), Alphabet.ProteinExtended, description, validate: false);
        if (EntryName is not null) seq.Annotations["entry_name"] = EntryName;
        if (OrganismName is not null) seq.Annotations["organism"] = OrganismName;
        if (TaxonId is { } tax) seq.Annotations["taxon_id"] = tax.ToString(System.Globalization.CultureInfo.InvariantCulture);
        return seq;
    }

    /// <summary>The entry as a Core <see cref="Protein"/> with name, organism, gene and function filled in.</summary>
    public Protein ToProtein()
    {
        var p = new Protein(Accession, ToSequence(), ProteinName ?? EntryName ?? Accession)
        {
            Organism = OrganismName,
            GeneName = GeneNames.Count > 0 ? GeneNames[0] : null,
            Function = Function,
            StructureId = CrossReferencesTo("PDB").FirstOrDefault()?.Id,
        };
        if (EntryName is not null) p.Annotations["entry_name"] = EntryName;
        p.Annotations["reviewed"] = IsReviewed ? "true" : "false";
        if (Keywords.Count > 0) p.Annotations["keywords"] = string.Join("; ", Keywords.Select(k => k.Name));
        return p;
    }

    public override string ToString() => $"{Accession} {EntryName} {ProteinName} ({OrganismName}, {Length} aa)";

    /// <summary>Parses one entry of the UniProtKB REST JSON format.</summary>
    public static UniProtEntry FromJson(JsonElement e)
    {
        var desc = e.Prop("proteinDescription");
        string? name = desc?.Path("recommendedName", "fullName")?.Str("value")
                       ?? desc?.Arr("submissionNames").FirstOrDefault().Path("fullName")?.Str("value");
        var alt = desc?.Arr("alternativeNames").Select(a => a.Path("fullName")?.Str("value")).OfType<string>().ToList() ?? [];

        var genes = new List<string>();
        foreach (var g in e.Arr("genes"))
        {
            if (g.Path("geneName")?.Str("value") is { } gn) genes.Add(gn);
            foreach (string list in (string[])["synonyms", "orderedLocusNames", "orfNames"])
                genes.AddRange(g.Arr(list).Select(s => s.Str("value")).OfType<string>());
        }

        var org = e.Prop("organism");
        var seq = e.Prop("sequence");
        var audit = e.Prop("entryAudit");
        return new UniProtEntry
        {
            Accession = e.Str("primaryAccession") ?? throw new DataServiceException("UniProt", "Entry JSON has no primaryAccession."),
            SecondaryAccessions = e.Strings("secondaryAccessions"),
            EntryName = e.Str("uniProtkbId"),
            IsReviewed = e.Str("entryType")?.Contains("Swiss-Prot", StringComparison.OrdinalIgnoreCase) == true
                         || (e.Str("entryType")?.StartsWith("UniProtKB reviewed", StringComparison.OrdinalIgnoreCase) ?? false),
            ProteinName = name,
            AlternativeNames = alt,
            GeneNames = genes,
            OrganismName = org?.Str("scientificName"),
            OrganismCommonName = org?.Str("commonName"),
            TaxonId = org?.Long("taxonId"),
            Lineage = org?.Strings("lineage") ?? [],
            ProteinExistence = e.Str("proteinExistence"),
            AnnotationScore = e.Double("annotationScore"),
            EntryVersion = audit?.Int("entryVersion"),
            SequenceVersion = audit?.Int("sequenceVersion"),
            SequenceResidues = seq?.Str("value") ?? "",
            Length = seq?.Int("length") ?? 0,
            MolecularWeight = seq?.Double("molWeight"),
            Crc64 = seq?.Str("crc64"),
            Comments = e.Arr("comments").Select(ParseComment).Where(c => c.Text.Length > 0).ToList(),
            Features = e.Arr("features").Select(f => new UniProtFeature(
                f.Str("type") ?? "",
                f.Str("description") ?? "",
                f.Path("location", "start")?.Int("value"),
                f.Path("location", "end")?.Int("value"),
                f.Str("featureId"))).ToList(),
            Keywords = e.Arr("keywords").Select(k => new UniProtKeyword(k.Str("id") ?? "", k.Str("name") ?? "", k.Str("category") ?? "")).ToList(),
            CrossReferences = e.Arr("uniProtKBCrossReferences").Select(x => new UniProtCrossReference(
                x.Str("database") ?? "",
                x.Str("id") ?? "",
                x.Arr("properties").Where(p => p.Str("key") is not null)
                    .GroupBy(p => p.Str("key")!)
                    .ToDictionary(g => g.Key, g => g.First().Str("value") ?? ""))).ToList(),
        };
    }

    private static UniProtComment ParseComment(JsonElement c)
    {
        string type = c.Str("commentType") ?? "";
        var parts = new List<string>();
        parts.AddRange(c.Arr("texts").Select(t => t.Str("value")).OfType<string>());
        parts.AddRange(c.Arr("subcellularLocations").Select(s => s.Path("location")?.Str("value")).OfType<string>());
        if (c.Path("disease") is { } d && d.Str("diseaseId") is { } did) parts.Add(did + (d.Str("description") is { } dd ? ": " + dd : ""));
        if (c.Path("reaction")?.Str("name") is { } reaction) parts.Add(reaction);
        if (c.Str("note") is { } note) parts.Add(note);
        else if (c.Path("note") is { } n) parts.AddRange(n.Arr("texts").Select(t => t.Str("value")).OfType<string>());
        return new UniProtComment(type, string.Join(" ", parts));
    }
}

/// <summary>One page of a UniProtKB search.</summary>
/// <param name="Entries">Entries on this page.</param>
/// <param name="TotalResults">Total hits (from <c>X-Total-Results</c>), when the server reports it.</param>
/// <param name="NextPage">URL of the next page (from the <c>Link: rel="next"</c> header), or null on the last page.</param>
public sealed record UniProtSearchPage(IReadOnlyList<UniProtEntry> Entries, long? TotalResults, Uri? NextPage);

/// <summary>State of a UniProt ID-mapping job.</summary>
public enum IdMappingJobState
{
    New,
    Running,
    Finished,
    Error,
}

/// <summary>One source-to-target identifier pair.</summary>
public sealed record IdMapping(string From, string To);

/// <summary>The result of an ID-mapping job.</summary>
public sealed record IdMappingResult(string JobId, IReadOnlyList<IdMapping> Mappings, IReadOnlyList<string> FailedIds)
{
    /// <summary>Targets grouped by source id (one source may map to several targets).</summary>
    public ILookup<string, string> BySource() => Mappings.ToLookup(m => m.From, m => m.To);
}
