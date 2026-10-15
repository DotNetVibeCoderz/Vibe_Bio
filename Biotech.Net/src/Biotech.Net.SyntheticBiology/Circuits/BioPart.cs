using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Circuits;

/// <summary>
/// The functional class of a standard biological part, following the part categories used by the
/// Registry of Standard Biological Parts (parts.igem.org) and the BioBrick assembly standard
/// (Knight, "Idempotent Vector Design for Standard Assembly of BioBricks", MIT Synthetic Biology
/// Working Group Technical Report, 2003).
/// </summary>
public enum PartType
{
    /// <summary>A DNA element bound by RNA polymerase holoenzyme that initiates transcription.</summary>
    Promoter,

    /// <summary>
    /// A ribosome binding site: the 5' untranslated region carrying the Shine-Dalgarno sequence that
    /// base-pairs with the 3' end of the 16S rRNA and positions the ribosome on the start codon
    /// (Shine &amp; Dalgarno, PNAS 71:1342, 1974).
    /// </summary>
    RibosomeBindingSite,

    /// <summary>A protein coding sequence, read from a start codon to a stop codon in a single frame.</summary>
    CodingSequence,

    /// <summary>
    /// A transcriptional terminator: in bacteria usually an intrinsic (Rho-independent) GC-rich hairpin
    /// followed by a poly-U tract that releases RNA polymerase.
    /// </summary>
    Terminator,

    /// <summary>A protein-binding operator site that modulates transcription of a nearby promoter.</summary>
    Operator,

    /// <summary>A neutral spacer or insulator sequence with no intended regulatory activity.</summary>
    Spacer,

    /// <summary>
    /// An assembly scar: the short sequence left between parts after BioBrick (RFC[10]) ligation,
    /// canonically the 8 bp scar <c>tactagag</c> between two assembled parts.
    /// </summary>
    Scar,
}

/// <summary>
/// A single standard biological part: a named DNA sequence with a declared function, optionally
/// annotated with regulatory behaviour and a Registry of Standard Biological Parts identifier.
/// </summary>
/// <param name="Name">Short human-readable name, e.g. <c>J23119</c> or <c>GFP</c>.</param>
/// <param name="Type">Functional class of the part.</param>
/// <param name="Sequence">The DNA sequence of the part, sense strand, 5' to 3'.</param>
/// <param name="Description">Free text: provenance, caveats, and whether the sequence is authentic or representative.</param>
public sealed record BioPart(string Name, PartType Type, string Sequence, string Description = "")
{
    /// <summary>Length of the part in base pairs.</summary>
    public int Length => Sequence.Length;

    /// <summary>
    /// Fraction of G or C bases (0..1). GC content correlates with duplex stability and is a
    /// practical proxy for synthesis and PCR difficulty.
    /// </summary>
    public double GcContent => Sequence.Length == 0 ? 0.0 : DnaUtility.GcFraction(Sequence);

    /// <summary>
    /// Relative activity of the part in <b>arbitrary units</b>, with no physical dimension.
    /// For a promoter this is transcription initiation activity relative to a reference promoter
    /// (the Anderson family is conventionally scaled to BBa_J23100 = 1); for a ribosome binding site
    /// it is a relative translation initiation rate. These numbers are catalogue-style estimates
    /// intended for <i>ranking</i> parts, not for predicting absolute expression: measured part
    /// activity depends strongly on host strain, growth phase, plasmid copy number and the
    /// neighbouring sequence context (Kelly et al., Journal of Biological Engineering 3:4, 2009).
    /// Defaults to 1.0.
    /// </summary>
    public double RelativeStrength { get; init; } = 1.0;

    /// <summary>
    /// Name of the transcription factor that acts on this part, if any (for example <c>LacI</c>,
    /// <c>TetR</c> or <c>cI</c>). Null for constitutive or unregulated parts.
    /// </summary>
    public string? RegulatedBy { get; init; }

    /// <summary>True when binding of <see cref="RegulatedBy"/> reduces the activity of the part.</summary>
    public bool IsRepressible { get; init; }

    /// <summary>True when a small-molecule inducer increases the activity of the part.</summary>
    public bool IsInducible { get; init; }

    /// <summary>
    /// The small molecule that relieves repression or activates the part, e.g. <c>IPTG</c> for the
    /// <c>lac</c> operator or <c>aTc</c> (anhydrotetracycline) for the <c>tet</c> operator.
    /// </summary>
    public string? Inducer { get; init; }

    /// <summary>
    /// Registry of Standard Biological Parts identifier, e.g. <c>BBa_R0010</c>, when the part
    /// corresponds to a catalogued BioBrick.
    /// </summary>
    public string? RegistryId { get; init; }

    /// <summary>
    /// Produces a <see cref="Core.SequenceFeature"/> for this part placed at <paramref name="start"/>
    /// (0-based, half-open) in a larger construct, with the feature key taken from the part type and
    /// <c>label</c> / <c>standard</c> qualifiers describing the part and its Registry identifier.
    /// </summary>
    /// <param name="start">0-based offset of the first base of the part within the parent construct.</param>
    public SequenceFeature ToFeature(int start)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        var feature = new SequenceFeature(FeatureKey(Type), FeatureLocation.Span(start, start + Length))
            .With("label", Name)
            .With("standard", RegistryId is null ? "BioBrick RFC[10]" : $"BioBrick RFC[10]:{RegistryId}");
        if (Description.Length > 0) feature.With("note", Description);
        if (RegulatedBy is not null) feature.With("regulated_by", RegulatedBy);
        if (Inducer is not null) feature.With("inducer", Inducer);
        return feature;
    }

    /// <summary>Maps a <see cref="PartType"/> to a GenBank-style feature key.</summary>
    /// <param name="type">The part class.</param>
    public static string FeatureKey(PartType type) => type switch
    {
        PartType.Promoter => "promoter",
        PartType.RibosomeBindingSite => "RBS",
        PartType.CodingSequence => "CDS",
        PartType.Terminator => "terminator",
        PartType.Operator => "protein_bind",
        PartType.Spacer => "misc_feature",
        PartType.Scar => "misc_feature",
        _ => "misc_feature",
    };

    /// <summary>A compact description such as <c>J23119 (Promoter, 35 bp)</c>.</summary>
    public override string ToString() => $"{Name} ({Type}, {Length} bp)";
}

/// <summary>
/// A transcriptional unit: an ordered list of parts assembled into one construct, conventionally
/// promoter - RBS - coding sequence - terminator. This is the unit of composition in the BioBrick
/// standard, where parts are concatenated head to tail and each junction carries an assembly scar.
/// </summary>
/// <param name="Name">Name of the assembled device.</param>
/// <param name="Parts">The parts in 5' to 3' order.</param>
public sealed record TranscriptionalUnit(string Name, IReadOnlyList<BioPart> Parts)
{
    /// <summary>
    /// Concatenates the part sequences in order, giving the DNA of the assembled construct.
    /// No scars are inserted automatically; include <see cref="PartType.Scar"/> parts explicitly
    /// if a scarred assembly standard is being modelled.
    /// </summary>
    public string Assemble()
    {
        int total = 0;
        for (int i = 0; i < Parts.Count; i++) total += Parts[i].Sequence.Length;
        var buffer = new char[total];
        int offset = 0;
        for (int i = 0; i < Parts.Count; i++)
        {
            Parts[i].Sequence.AsSpan().CopyTo(buffer.AsSpan(offset));
            offset += Parts[i].Sequence.Length;
        }
        return new string(buffer);
    }

    /// <summary>Total length of the assembled construct in base pairs.</summary>
    public int Length
    {
        get
        {
            int total = 0;
            for (int i = 0; i < Parts.Count; i++) total += Parts[i].Sequence.Length;
            return total;
        }
    }

    /// <summary>
    /// Annotation of the assembled construct: one feature per part, with contiguous 0-based
    /// half-open coordinates that tile the whole construct in order.
    /// </summary>
    public IReadOnlyList<SequenceFeature> Features()
    {
        var features = new List<SequenceFeature>(Parts.Count);
        int offset = 0;
        for (int i = 0; i < Parts.Count; i++)
        {
            features.Add(Parts[i].ToFeature(offset));
            offset += Parts[i].Sequence.Length;
        }
        return features;
    }

    /// <summary>The assembled construct as a Core <see cref="Sequence"/> identified by <see cref="Name"/>.</summary>
    public Sequence ToSequence() => Sequence.Dna(Assemble(), Name);

    /// <summary>A compact description such as <c>gfp_device: 5 parts, 1024 bp</c>.</summary>
    public override string ToString() => $"{Name}: {Parts.Count} parts, {Length} bp";
}
