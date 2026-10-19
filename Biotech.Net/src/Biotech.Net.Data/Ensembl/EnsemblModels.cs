using System.Text.Json;
using Biotech.Net.Core;
using Biotech.Net.Data.Http;
using Biotech.Net.Genomics.Annotation;

namespace Biotech.Net.Data.Ensembl;

/// <summary>Sequence types of the Ensembl <c>/sequence/id</c> endpoint.</summary>
public enum EnsemblSequenceType
{
    /// <summary>The genomic span of the feature (introns included).</summary>
    Genomic,

    /// <summary>Spliced transcript (exons, UTRs included).</summary>
    Cdna,

    /// <summary>Coding sequence only.</summary>
    Cds,

    /// <summary>Translated protein.</summary>
    Protein,
}

/// <summary>
/// An Ensembl feature from <c>/lookup</c>: a gene, transcript, exon or translation. Coordinates are as Ensembl
/// gives them: 1-based and inclusive; <see cref="ToInterval"/> converts to Biotech.Net's 0-based half-open form.
/// </summary>
public sealed record EnsemblFeature
{
    public required string Id { get; init; }

    /// <summary><c>Gene</c>, <c>Transcript</c>, <c>Exon</c> or <c>Translation</c>.</summary>
    public string ObjectType { get; init; } = "";

    public string? DisplayName { get; init; }
    public string? Species { get; init; }
    public string? Biotype { get; init; }
    public string? Description { get; init; }
    public string? AssemblyName { get; init; }
    public string? SeqRegionName { get; init; }

    /// <summary>1-based first base.</summary>
    public long Start { get; init; }

    /// <summary>1-based last base (inclusive).</summary>
    public long End { get; init; }

    /// <summary>1 or -1.</summary>
    public int Strand { get; init; }

    public int? Version { get; init; }
    public string? Parent { get; init; }
    public string? CanonicalTranscript { get; init; }
    public bool IsCanonical { get; init; }
    public int? Length { get; init; }

    /// <summary>Transcripts of a gene (lookup with <c>expand</c>).</summary>
    public IReadOnlyList<EnsemblFeature> Transcripts { get; init; } = [];

    /// <summary>Exons of a transcript (lookup with <c>expand</c>).</summary>
    public IReadOnlyList<EnsemblFeature> Exons { get; init; } = [];

    /// <summary>The translation of a coding transcript (lookup with <c>expand</c>).</summary>
    public EnsemblFeature? Translation { get; init; }

    /// <summary>The stable id with its version, e.g. <c>ENSG00000139618.19</c>.</summary>
    public string VersionedId => Version is { } v ? $"{Id}.{v}" : Id;

    /// <summary>The span as a 0-based, half-open <see cref="GenomicInterval"/>.</summary>
    public GenomicInterval ToInterval() =>
        GenomicInterval.Create(SeqRegionName ?? "?", Start - 1, End, Strand switch { 1 => Core.Strand.Forward, -1 => Core.Strand.Reverse, _ => Core.Strand.Unknown });

    /// <summary>
    /// A Core <see cref="Gene"/> with its transcripts, exons and coding regions (from a gene looked up with <c>expand</c>).
    /// </summary>
    public Gene ToGene()
    {
        if (!ObjectType.Equals("Gene", StringComparison.OrdinalIgnoreCase))
            throw new InvalidOperationException($"{Id} is a {ObjectType}, not a Gene.");
        var gene = new Gene(Id, ToInterval(), DisplayName, Biotype ?? "unknown");
        if (Description is not null) gene.Attributes["description"] = Description;
        if (CanonicalTranscript is not null) gene.Attributes["canonical_transcript"] = CanonicalTranscript;
        if (AssemblyName is not null) gene.Attributes["assembly"] = AssemblyName;
        foreach (var t in Transcripts)
        {
            var exons = t.Exons.Count > 0 ? t.Exons.Select(e => e.ToInterval()) : [t.ToInterval()];
            gene.Transcripts.Add(new Transcript(t.Id, Id, exons)
            {
                CodingRegion = t.Translation?.ToInterval(),
            });
        }

        return gene;
    }

    public override string ToString() => $"{ObjectType} {DisplayName ?? Id} ({Id}) {SeqRegionName}:{Start}-{End}:{Strand}";

    internal static EnsemblFeature FromJson(JsonElement e) => new()
    {
        Id = e.Str("id") ?? "",
        ObjectType = e.Str("object_type") ?? "",
        DisplayName = e.Str("display_name"),
        Species = e.Str("species"),
        Biotype = e.Str("biotype"),
        Description = e.Str("description"),
        AssemblyName = e.Str("assembly_name"),
        SeqRegionName = e.Str("seq_region_name"),
        Start = e.Long("start") ?? 0,
        End = e.Long("end") ?? 0,
        Strand = e.Int("strand") ?? 0,
        Version = e.Int("version"),
        Parent = e.Str("Parent"),
        CanonicalTranscript = e.Str("canonical_transcript"),
        IsCanonical = e.Bool("is_canonical"),
        Length = e.Int("length"),
        Transcripts = e.Arr("Transcript").Select(FromJson).ToList(),
        Exons = e.Arr("Exon").Select(FromJson).ToList(),
        // A translation carries no region or strand of its own: it inherits the transcript's.
        Translation = e.Prop("Translation") is { ValueKind: JsonValueKind.Object } t
            ? FromJson(t) with { SeqRegionName = e.Str("seq_region_name"), Strand = e.Int("strand") ?? 0 }
            : null,
    };
}

/// <summary>A feature returned by <c>/overlap/region</c> (gene, transcript, exon, variation, regulatory, ...).</summary>
public sealed record EnsemblOverlapFeature
{
    /// <summary><c>gene</c>, <c>transcript</c>, <c>exon</c>, <c>cds</c>, <c>variation</c>, <c>regulatory</c>, ...</summary>
    public required string FeatureType { get; init; }

    public required string Id { get; init; }
    public string? SeqRegionName { get; init; }

    /// <summary>1-based first base.</summary>
    public long Start { get; init; }

    /// <summary>1-based last base (inclusive).</summary>
    public long End { get; init; }

    public int Strand { get; init; }
    public string? Biotype { get; init; }
    public string? ExternalName { get; init; }
    public string? Description { get; init; }
    public string? Parent { get; init; }
    public bool IsCanonical { get; init; }

    /// <summary>For variations: the most severe consequence.</summary>
    public string? ConsequenceType { get; init; }

    /// <summary>For variations: the alleles.</summary>
    public IReadOnlyList<string> Alleles { get; init; } = [];

    /// <summary>For variations: clinical significance terms.</summary>
    public IReadOnlyList<string> ClinicalSignificance { get; init; } = [];

    public GenomicInterval ToInterval() =>
        GenomicInterval.Create(SeqRegionName ?? "?", Start - 1, End, Strand switch { 1 => Core.Strand.Forward, -1 => Core.Strand.Reverse, _ => Core.Strand.Unknown });

    public override string ToString() => $"{FeatureType} {ExternalName ?? Id} {SeqRegionName}:{Start}-{End}";

    internal static EnsemblOverlapFeature FromJson(JsonElement e) => new()
    {
        FeatureType = e.Str("feature_type") ?? "",
        Id = e.Str("id") ?? e.Str("transcript_id") ?? e.Str("gene_id") ?? "",
        SeqRegionName = e.Str("seq_region_name"),
        Start = e.Long("start") ?? 0,
        End = e.Long("end") ?? 0,
        Strand = e.Int("strand") ?? 0,
        Biotype = e.Str("biotype"),
        ExternalName = e.Str("external_name"),
        Description = e.Str("description"),
        Parent = e.Str("Parent"),
        IsCanonical = e.Bool("is_canonical"),
        ConsequenceType = e.Str("consequence_type"),
        Alleles = e.Strings("alleles"),
        ClinicalSignificance = e.Strings("clinical_significance"),
    };
}

/// <summary>The predicted effect of a variant on one transcript (VEP).</summary>
public sealed record VepTranscriptConsequence
{
    public required string TranscriptId { get; init; }
    public string? GeneId { get; init; }
    public string? GeneSymbol { get; init; }
    public string? Biotype { get; init; }
    public string? VariantAllele { get; init; }

    /// <summary>Sequence Ontology terms, e.g. <c>missense_variant</c>.</summary>
    public IReadOnlyList<string> ConsequenceTerms { get; init; } = [];

    /// <summary><c>HIGH</c>, <c>MODERATE</c>, <c>LOW</c> or <c>MODIFIER</c>.</summary>
    public string? Impact { get; init; }

    /// <summary>Reference/alternate amino acids, e.g. <c>E/V</c>.</summary>
    public string? AminoAcids { get; init; }

    /// <summary>Reference/alternate codons, e.g. <c>gAg/gTg</c>.</summary>
    public string? Codons { get; init; }

    public int? ProteinStart { get; init; }
    public int? ProteinEnd { get; init; }
    public int? CdsStart { get; init; }
    public int? CdsEnd { get; init; }
    public string? SiftPrediction { get; init; }
    public double? SiftScore { get; init; }
    public string? PolyPhenPrediction { get; init; }
    public double? PolyPhenScore { get; init; }
    public bool IsCanonical { get; init; }

    /// <summary>
    /// The consequence terms that Biotech.Net's own <see cref="Consequence"/> enumeration models (terms it does not
    /// model, such as <c>NMD_transcript_variant</c>, are left out).
    /// </summary>
    public IReadOnlyList<Consequence> Consequences => VepResult.ToConsequences(ConsequenceTerms);

    /// <summary>HGVS-style protein change, e.g. <c>E7V</c>, when VEP reports amino acids and a position.</summary>
    public string? ProteinChange => AminoAcids is { } aa && aa.Split('/') is [var r, var a] && ProteinStart is { } p ? $"{r}{p}{a}" : null;

    internal static VepTranscriptConsequence FromJson(JsonElement e) => new()
    {
        TranscriptId = e.Str("transcript_id") ?? "",
        GeneId = e.Str("gene_id"),
        GeneSymbol = e.Str("gene_symbol"),
        Biotype = e.Str("biotype"),
        VariantAllele = e.Str("variant_allele"),
        ConsequenceTerms = e.Strings("consequence_terms"),
        Impact = e.Str("impact"),
        AminoAcids = e.Str("amino_acids"),
        Codons = e.Str("codons"),
        ProteinStart = e.Int("protein_start"),
        ProteinEnd = e.Int("protein_end"),
        CdsStart = e.Int("cds_start"),
        CdsEnd = e.Int("cds_end"),
        SiftPrediction = e.Str("sift_prediction"),
        SiftScore = e.Double("sift_score"),
        PolyPhenPrediction = e.Str("polyphen_prediction"),
        PolyPhenScore = e.Double("polyphen_score"),
        IsCanonical = e.Bool("canonical"),
    };
}

/// <summary>One Variant Effect Predictor result (one input variant).</summary>
public sealed record VepResult
{
    /// <summary>The input as submitted (HGVS, region string or VCF line).</summary>
    public string? Input { get; init; }

    public string? Id { get; init; }
    public string? AssemblyName { get; init; }
    public string? SeqRegionName { get; init; }

    /// <summary>1-based start.</summary>
    public long Start { get; init; }

    /// <summary>1-based end (inclusive).</summary>
    public long End { get; init; }

    public int Strand { get; init; }

    /// <summary>E.g. <c>A/T</c>.</summary>
    public string? AlleleString { get; init; }

    public string? MostSevereConsequence { get; init; }
    public IReadOnlyList<VepTranscriptConsequence> TranscriptConsequences { get; init; } = [];

    /// <summary>Known variants at the same position (dbSNP rs ids, COSMIC, HGMD, ...).</summary>
    public IReadOnlyList<string> ColocatedVariantIds { get; init; } = [];

    /// <summary>The canonical-transcript consequence if any, else the first one.</summary>
    public VepTranscriptConsequence? Canonical =>
        TranscriptConsequences.FirstOrDefault(c => c.IsCanonical) ?? TranscriptConsequences.FirstOrDefault();

    /// <summary><see cref="MostSevereConsequence"/> as a Biotech.Net <see cref="Consequence"/>, when modelled.</summary>
    public Consequence? MostSevere => MostSevereConsequence is { } t && ToConsequences([t]) is [var c] ? c : null;

    public override string ToString() => $"{Input ?? Id} {SeqRegionName}:{Start} {AlleleString} {MostSevereConsequence}";

    internal static IReadOnlyList<Consequence> ToConsequences(IEnumerable<string> terms)
    {
        var list = new List<Consequence>();
        foreach (string term in terms)
        {
            try
            {
                list.Add(ConsequenceExtensions.ParseSoTerm(term));
            }
            catch (ArgumentException)
            {
                // A term Biotech.Net does not model.
            }
        }

        return list;
    }

    internal static VepResult FromJson(JsonElement e) => new()
    {
        Input = e.Str("input"),
        Id = e.Str("id"),
        AssemblyName = e.Str("assembly_name"),
        SeqRegionName = e.Str("seq_region_name"),
        Start = e.Long("start") ?? 0,
        End = e.Long("end") ?? 0,
        Strand = e.Int("strand") ?? 0,
        AlleleString = e.Str("allele_string"),
        MostSevereConsequence = e.Str("most_severe_consequence"),
        TranscriptConsequences = e.Arr("transcript_consequences").Select(VepTranscriptConsequence.FromJson).ToList(),
        ColocatedVariantIds = e.Arr("colocated_variants").Select(c => c.Str("id")).OfType<string>().ToList(),
    };
}

/// <summary>An external reference (<c>/xrefs/id</c>): HGNC, UniProt, RefSeq, EntrezGene, ...</summary>
public sealed record EnsemblXref(string DbName, string DbDisplayName, string PrimaryId, string DisplayId, string? Description,
    IReadOnlyList<string> Synonyms)
{
    internal static EnsemblXref FromJson(JsonElement e) => new(
        e.Str("dbname") ?? "", e.Str("db_display_name") ?? "", e.Str("primary_id") ?? "", e.Str("display_id") ?? "",
        e.Str("description"), e.Strings("synonyms"));
}

/// <summary>An Ensembl object found by symbol (<c>/xrefs/symbol</c>).</summary>
public sealed record EnsemblIdentity(string Id, string Type);
