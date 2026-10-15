namespace Biotech.Net.Genomics.Annotation;

/// <summary>
/// Sequence Ontology variant consequence terms, declared from most to least severe (Ensembl VEP ranking), so the
/// enum order can be used to pick the most severe consequence.
/// </summary>
public enum Consequence
{
    /// <summary>SO:0001574 — changes the 2 bases at the 3' end of an intron.</summary>
    SpliceAcceptorVariant,
    /// <summary>SO:0001575 — changes the 2 bases at the 5' end of an intron.</summary>
    SpliceDonorVariant,
    /// <summary>SO:0001587 — introduces a premature stop codon.</summary>
    StopGained,
    /// <summary>SO:0001589 — indel whose length is not a multiple of three.</summary>
    FrameshiftVariant,
    /// <summary>SO:0001578 — the stop codon becomes a sense codon.</summary>
    StopLost,
    /// <summary>SO:0002012 — the initiator codon is changed to a non-start codon.</summary>
    StartLost,
    /// <summary>SO:0001821 — in-frame insertion of amino acids.</summary>
    InframeInsertion,
    /// <summary>SO:0001822 — in-frame deletion of amino acids.</summary>
    InframeDeletion,
    /// <summary>SO:0001583 — changes one or more amino acids, length unchanged.</summary>
    MissenseVariant,
    /// <summary>SO:0001818 — changes the protein in a way not covered by a more specific term.</summary>
    ProteinAlteringVariant,
    /// <summary>SO:0001630 — within 1–3 exonic or 3–8 intronic bases of a splice site.</summary>
    SpliceRegionVariant,
    /// <summary>SO:0002019 — the initiator codon changes but remains a start codon.</summary>
    StartRetainedVariant,
    /// <summary>SO:0001567 — the stop codon changes but remains a stop codon.</summary>
    StopRetainedVariant,
    /// <summary>SO:0001819 — coding change that leaves the amino acid unchanged.</summary>
    SynonymousVariant,
    /// <summary>SO:0001580 — in the coding sequence, effect not determined (e.g. spans an exon boundary).</summary>
    CodingSequenceVariant,
    /// <summary>SO:0001623 — in the 5' untranslated region.</summary>
    FivePrimeUtrVariant,
    /// <summary>SO:0001624 — in the 3' untranslated region.</summary>
    ThreePrimeUtrVariant,
    /// <summary>SO:0001792 — in an exon of a non-coding transcript.</summary>
    NonCodingTranscriptExonVariant,
    /// <summary>SO:0001627 — in an intron.</summary>
    IntronVariant,
    /// <summary>SO:0001631 — 5' of a gene (within the upstream window).</summary>
    UpstreamGeneVariant,
    /// <summary>SO:0001632 — 3' of a gene (within the downstream window).</summary>
    DownstreamGeneVariant,
    /// <summary>SO:0001628 — not near any gene.</summary>
    IntergenicVariant,
}

/// <summary>Putative impact class of a consequence (as used by VEP and SnpEff).</summary>
public enum VariantImpact
{
    /// <summary>Probably disruptive (loss of function).</summary>
    High,
    /// <summary>Non-disruptive but may change protein effectiveness.</summary>
    Moderate,
    /// <summary>Mostly harmless.</summary>
    Low,
    /// <summary>Non-coding, effect hard to predict.</summary>
    Modifier,
}

/// <summary>Sequence Ontology names, accessions and impacts for <see cref="Consequence"/>.</summary>
public static class ConsequenceExtensions
{
    private static readonly (string Term, string Accession, VariantImpact Impact)[] Table =
    [
        ("splice_acceptor_variant", "SO:0001574", VariantImpact.High),
        ("splice_donor_variant", "SO:0001575", VariantImpact.High),
        ("stop_gained", "SO:0001587", VariantImpact.High),
        ("frameshift_variant", "SO:0001589", VariantImpact.High),
        ("stop_lost", "SO:0001578", VariantImpact.High),
        ("start_lost", "SO:0002012", VariantImpact.High),
        ("inframe_insertion", "SO:0001821", VariantImpact.Moderate),
        ("inframe_deletion", "SO:0001822", VariantImpact.Moderate),
        ("missense_variant", "SO:0001583", VariantImpact.Moderate),
        ("protein_altering_variant", "SO:0001818", VariantImpact.Moderate),
        ("splice_region_variant", "SO:0001630", VariantImpact.Low),
        ("start_retained_variant", "SO:0002019", VariantImpact.Low),
        ("stop_retained_variant", "SO:0001567", VariantImpact.Low),
        ("synonymous_variant", "SO:0001819", VariantImpact.Low),
        ("coding_sequence_variant", "SO:0001580", VariantImpact.Modifier),
        ("5_prime_UTR_variant", "SO:0001623", VariantImpact.Modifier),
        ("3_prime_UTR_variant", "SO:0001624", VariantImpact.Modifier),
        ("non_coding_transcript_exon_variant", "SO:0001792", VariantImpact.Modifier),
        ("intron_variant", "SO:0001627", VariantImpact.Modifier),
        ("upstream_gene_variant", "SO:0001631", VariantImpact.Modifier),
        ("downstream_gene_variant", "SO:0001632", VariantImpact.Modifier),
        ("intergenic_variant", "SO:0001628", VariantImpact.Modifier),
    ];

    /// <summary>The SO term name, e.g. <c>missense_variant</c>.</summary>
    public static string ToSoTerm(this Consequence c) => Table[(int)c].Term;

    /// <summary>The SO accession, e.g. <c>SO:0001583</c>.</summary>
    public static string SoAccession(this Consequence c) => Table[(int)c].Accession;

    /// <summary>The impact class.</summary>
    public static VariantImpact Impact(this Consequence c) => Table[(int)c].Impact;

    /// <summary>Parses an SO term name.</summary>
    public static Consequence ParseSoTerm(string term)
    {
        for (int i = 0; i < Table.Length; i++)
            if (Table[i].Term.Equals(term, StringComparison.OrdinalIgnoreCase)) return (Consequence)i;
        throw new ArgumentException($"Unknown consequence term '{term}'.", nameof(term));
    }
}

/// <summary>The predicted effect of one alternate allele on one transcript (or intergenic).</summary>
public sealed record VariantEffect
{
    /// <summary>The alternate allele as given in the VCF.</summary>
    public required string Allele { get; init; }

    /// <summary>Consequences, most severe first.</summary>
    public required IReadOnlyList<Consequence> Consequences { get; init; }

    /// <summary>The most severe consequence.</summary>
    public Consequence MostSevere => Consequences[0];

    /// <summary>Impact of the most severe consequence.</summary>
    public VariantImpact Impact => MostSevere.Impact();

    /// <summary>Gene identifier (null for intergenic).</summary>
    public string? GeneId { get; init; }

    /// <summary>Gene symbol.</summary>
    public string? GeneName { get; init; }

    /// <summary>Transcript identifier.</summary>
    public string? TranscriptId { get; init; }

    /// <summary>HGVS coding (<c>c.</c>) or non-coding (<c>n.</c>) notation, without the transcript prefix.</summary>
    public string? HgvsC { get; init; }

    /// <summary>HGVS protein notation (three-letter, e.g. <c>p.Arg12Trp</c>), without the protein prefix.</summary>
    public string? HgvsP { get; init; }

    /// <summary>1-based position of the first affected CDS base.</summary>
    public int? CdsPosition { get; init; }

    /// <summary>1-based position of the first affected amino acid.</summary>
    public int? ProteinPosition { get; init; }

    /// <summary>Reference/alternate amino acids, VEP style (e.g. <c>R/W</c>).</summary>
    public string? AminoAcids { get; init; }

    /// <summary>Reference/alternate codons, changed bases in upper case (e.g. <c>Cgg/Tgg</c>).</summary>
    public string? Codons { get; init; }

    /// <summary>Distance to the transcript for upstream/downstream variants.</summary>
    public long? Distance { get; init; }

    /// <summary>The consequences joined with <c>&amp;</c>, VEP style.</summary>
    public string ConsequenceTerms => string.Join('&', Consequences.Select(c => c.ToSoTerm()));

    /// <summary>VEP-style CSQ entry: <c>Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|Protein_position|Amino_acids|Codons</c>.</summary>
    public string ToCsq() => string.Join('|', Allele, ConsequenceTerms, Impact.ToString().ToUpperInvariant(), GeneName ?? "", GeneId ?? "",
        TranscriptId ?? "", HgvsC ?? "", HgvsP ?? "", ProteinPosition?.ToString() ?? "", AminoAcids ?? "", Codons ?? "");

    public override string ToString() => $"{Allele} {ConsequenceTerms} {GeneName}/{TranscriptId} {HgvsC} {HgvsP}".Trim();
}
