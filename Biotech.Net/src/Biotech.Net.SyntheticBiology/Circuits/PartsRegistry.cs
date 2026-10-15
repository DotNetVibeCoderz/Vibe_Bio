namespace Biotech.Net.SyntheticBiology.Circuits;

/// <summary>
/// A small built-in catalogue of iGEM-style standard biological parts, modelled on the Registry of
/// Standard Biological Parts (parts.igem.org).
/// <para>
/// <b>Provenance and honesty.</b> Short, widely reproduced sequences (the 35 bp Anderson constitutive
/// promoter family, the 12 bp BBa_B0034 ribosome binding site, the lac and tet operator sites) are
/// given as published. Longer regulatory regions are given as <i>representative</i> sequences that
/// carry the correct functional motifs but are not guaranteed base-for-base identical to the Registry
/// entry. All coding sequences in this catalogue are deliberately <b>truncated placeholders</b>: the
/// authentic open reading frames are 600-1200 bp and are not reproduced here, so each CDS part holds
/// only the first few authentic codons followed by a stop codon, and says so in its description.
/// Use this catalogue for composition, annotation and modelling exercises; fetch the real sequence
/// from the Registry or GenBank before ordering DNA.
/// </para>
/// </summary>
public static class PartsRegistry
{
    private static readonly BioPart[] Catalogue = BuildCatalogue();
    private static readonly Dictionary<string, BioPart> ByName = BuildIndex(Catalogue);

    /// <summary>Every part in the catalogue, in catalogue order (promoters, RBSs, CDSs, terminators, operators, scar).</summary>
    public static IReadOnlyList<BioPart> All => Catalogue;

    /// <summary>
    /// Looks up a part by short name (e.g. <c>J23119</c>) or by Registry identifier (e.g. <c>BBa_J23119</c>),
    /// case-insensitively.
    /// </summary>
    /// <param name="name">Short name or Registry identifier.</param>
    /// <exception cref="KeyNotFoundException">Thrown when no part matches, with a list of the closest names.</exception>
    public static BioPart Get(string name)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        if (ByName.TryGetValue(name, out var part)) return part;

        var suggestions = new List<string>();
        foreach (var candidate in Catalogue)
        {
            if (candidate.Name.Contains(name, StringComparison.OrdinalIgnoreCase) ||
                name.Contains(candidate.Name, StringComparison.OrdinalIgnoreCase))
            {
                suggestions.Add(candidate.Name);
            }
        }
        string hint = suggestions.Count > 0
            ? $" Did you mean: {string.Join(", ", suggestions)}?"
            : $" The catalogue holds {Catalogue.Length} parts; enumerate PartsRegistry.All to see them.";
        throw new KeyNotFoundException($"No part named '{name}' in the built-in parts registry.{hint}");
    }

    /// <summary>Attempts to look up a part by short name or Registry identifier, case-insensitively.</summary>
    /// <param name="name">Short name or Registry identifier.</param>
    /// <param name="part">The matching part when the method returns true.</param>
    public static bool TryGet(string name, out BioPart part)
    {
        if (string.IsNullOrEmpty(name))
        {
            part = null!;
            return false;
        }
        return ByName.TryGetValue(name, out part!);
    }

    /// <summary>All catalogue parts of a given functional class.</summary>
    /// <param name="type">The part class to filter by.</param>
    public static IEnumerable<BioPart> OfType(PartType type)
    {
        foreach (var part in Catalogue)
        {
            if (part.Type == type) yield return part;
        }
    }

    /// <summary>
    /// Free-text search over part name, Registry identifier, description and regulator name,
    /// case-insensitive substring matching.
    /// </summary>
    /// <param name="text">The text to look for; an empty query returns everything.</param>
    public static IEnumerable<BioPart> Search(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
        {
            foreach (var part in Catalogue) yield return part;
            yield break;
        }
        foreach (var part in Catalogue)
        {
            if (part.Name.Contains(text, StringComparison.OrdinalIgnoreCase) ||
                (part.RegistryId?.Contains(text, StringComparison.OrdinalIgnoreCase) ?? false) ||
                part.Description.Contains(text, StringComparison.OrdinalIgnoreCase) ||
                (part.RegulatedBy?.Contains(text, StringComparison.OrdinalIgnoreCase) ?? false) ||
                (part.Inducer?.Contains(text, StringComparison.OrdinalIgnoreCase) ?? false))
            {
                yield return part;
            }
        }
    }

    /// <summary>
    /// Composes a <see cref="TranscriptionalUnit"/> from catalogue part names, in 5' to 3' order,
    /// for example <c>BuildUnit("gfp_device", "J23119", "B0034", "GFP", "B0015")</c>.
    /// </summary>
    /// <param name="name">Name for the assembled device.</param>
    /// <param name="partNames">Catalogue part names or Registry identifiers, in order.</param>
    /// <exception cref="KeyNotFoundException">Thrown when any name is not in the catalogue.</exception>
    public static TranscriptionalUnit BuildUnit(string name, params string[] partNames)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentNullException.ThrowIfNull(partNames);
        var parts = new BioPart[partNames.Length];
        for (int i = 0; i < partNames.Length; i++) parts[i] = Get(partNames[i]);
        return new TranscriptionalUnit(name, parts);
    }

    private static Dictionary<string, BioPart> BuildIndex(BioPart[] parts)
    {
        var index = new Dictionary<string, BioPart>(StringComparer.OrdinalIgnoreCase);
        foreach (var part in parts)
        {
            index[part.Name] = part;
            if (part.RegistryId is not null) index[part.RegistryId] = part;
        }
        return index;
    }

    private static BioPart[] BuildCatalogue() =>
    [
        // ---------------- Promoters ----------------

        // The Anderson constitutive promoter family (BBa_J23100..J23119) is a set of 35 bp sigma-70
        // promoters differing only in the -35 and -10 hexamers and the spacer. The sequences below are
        // the published family members; relative strengths are the Registry's fluorescence measurements
        // scaled to J23100 = 1.0 and vary by roughly a factor of two between laboratories.
        new BioPart("J23100", PartType.Promoter, "ttgacggctagctcagtcctaggtacagtgctagc",
            "Anderson family constitutive sigma-70 promoter; the family reference, conventionally scored 1.00.")
        { RelativeStrength = 1.00, RegistryId = "BBa_J23100" },

        new BioPart("J23106", PartType.Promoter, "tttacggctagctcagtcctaggtatagtgctagc",
            "Anderson family constitutive sigma-70 promoter of intermediate strength (~0.47 of J23100).")
        { RelativeStrength = 0.47, RegistryId = "BBa_J23106" },

        new BioPart("J23119", PartType.Promoter, "ttgacagctagctcagtcctaggtataatgctagc",
            "Anderson family 'consensus' promoter, carrying the canonical -35 TTGACA and -10 TATAAT hexamers; the strongest member of the family.")
        { RelativeStrength = 1.20, RegistryId = "BBa_J23119" },

        new BioPart("R0010", PartType.Promoter, "caatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacaca",
            "approximate/representative sequence: the lac promoter region (pLac) with the lacO1 operator. Carries the authentic lacO1 site ggaattgtgagcggataacaattt and the CAP site, but the full-length Registry entry may differ in the flanking bases.")
        { RelativeStrength = 0.30, RegulatedBy = "LacI", IsRepressible = true, IsInducible = true, Inducer = "IPTG", RegistryId = "BBa_R0010" },

        new BioPart("R0011", PartType.Promoter, "aattgtgagcggataacaattgacattgtgagcggataacaagatactgagcaca",
            "approximate/representative sequence: a hybrid lacI-repressed promoter with two lac operator sites flanking the -35/-10 region, weaker and tighter than R0010.")
        { RelativeStrength = 0.15, RegulatedBy = "LacI", IsRepressible = true, IsInducible = true, Inducer = "IPTG", RegistryId = "BBa_R0011" },

        new BioPart("R0040", PartType.Promoter, "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac",
            "The tetR-repressed promoter pTet (PLtetO-1 style), with two tetO2 operators. Widely reproduced 54 bp sequence; relieved by anhydrotetracycline.")
        { RelativeStrength = 0.60, RegulatedBy = "TetR", IsRepressible = true, IsInducible = true, Inducer = "aTc", RegistryId = "BBa_R0040" },

        new BioPart("R0051", PartType.Promoter, "taacaccgtgcgtgttgactattttacctctggcggtgataatggttgc",
            "approximate/representative sequence: the bacteriophage lambda PR promoter with the OR1/OR2 operators, repressed by the lambda cI repressor.")
        { RelativeStrength = 0.55, RegulatedBy = "cI", IsRepressible = true, RegistryId = "BBa_R0051" },

        // ---------------- Ribosome binding sites ----------------

        // Relative strength ordering in the Registry's community RBS family is
        // B0030 > B0034 > B0032 > B0033, spanning roughly two orders of magnitude.
        new BioPart("B0030", PartType.RibosomeBindingSite, "attaaagaggagaaa",
            "approximate/representative sequence: strong community RBS, an extended variant of the B0034 Shine-Dalgarno core. Strongest of the B003x family.")
        { RelativeStrength = 1.60, RegistryId = "BBa_B0030" },

        new BioPart("B0032", PartType.RibosomeBindingSite, "tcacacaggaaag",
            "approximate/representative sequence: medium-weak community RBS; weaker than B0034, stronger than B0033.")
        { RelativeStrength = 0.30, RegistryId = "BBa_B0032" },

        new BioPart("B0033", PartType.RibosomeBindingSite, "tcacacaggac",
            "approximate/representative sequence: weak community RBS with a degenerate Shine-Dalgarno; the weakest of the B003x family.")
        { RelativeStrength = 0.01, RegistryId = "BBa_B0033" },

        new BioPart("B0034", PartType.RibosomeBindingSite, "aaagaggagaaa",
            "The canonical 12 bp BioBrick RBS, a near-consensus Shine-Dalgarno (AGGAGG) element; the de facto standard medium-strong RBS.")
        { RelativeStrength = 1.00, RegistryId = "BBa_B0034" },

        // ---------------- Coding sequences (truncated placeholders) ----------------

        new BioPart("GFP", PartType.CodingSequence, "atgcgtaaaggagaagaacttttcactggataa",
            "TRUNCATED PLACEHOLDER, not a functional CDS: the first ten codons of the GFPmut3 open reading frame followed by a stop codon. The authentic BBa_E0040 coding sequence is 720 bp and is not reproduced here.")
        { RegistryId = "BBa_E0040" },

        new BioPart("mCherry", PartType.CodingSequence, "atggtgagcaagggcgaggaggataacatgtaa",
            "TRUNCATED PLACEHOLDER, not a functional CDS: the first ten codons of the mCherry open reading frame followed by a stop codon. The authentic mCherry CDS is 711 bp.")
        { RegistryId = "BBa_J06504" },

        new BioPart("LacI", PartType.CodingSequence, "atgaaaccagtaacgttatacgatgtcgcagagtatgcctaa",
            "TRUNCATED PLACEHOLDER, not a functional CDS: the first thirteen codons of lacI followed by a stop codon. The authentic BBa_C0012 coding sequence is ~1083 bp.")
        { RegulatedBy = "IPTG", RegistryId = "BBa_C0012" },

        new BioPart("TetR", PartType.CodingSequence, "atgtccagattagataaaagtaaagtgattaacagcgcataa",
            "TRUNCATED PLACEHOLDER, not a functional CDS: the first thirteen codons of tetR followed by a stop codon. The authentic BBa_C0040 coding sequence is ~621 bp.")
        { RegulatedBy = "aTc", RegistryId = "BBa_C0040" },

        new BioPart("cI", PartType.CodingSequence, "atgagcacaaaaaagaaaccattaacacaagagcagcttgaataa",
            "TRUNCATED PLACEHOLDER, not a functional CDS: the first fourteen codons of the bacteriophage lambda cI repressor followed by a stop codon. The authentic BBa_C0051 coding sequence is ~714 bp.")
        { RegistryId = "BBa_C0051" },

        // ---------------- Terminators ----------------

        new BioPart("B0010", PartType.Terminator,
            "ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctc",
            "Intrinsic (Rho-independent) terminator derived from the E. coli rrnB T1 terminator: a GC-rich hairpin followed by a U-rich tract.")
        { RegistryId = "BBa_B0010" },

        new BioPart("B0012", PartType.Terminator,
            "tcacactggctcaccttcgggtgggcctttctgcgtttata",
            "Intrinsic terminator TE derived from coliphage T7; commonly used as the second half of the B0015 double terminator.")
        { RegistryId = "BBa_B0012" },

        new BioPart("B0015", PartType.Terminator,
            "ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctcacactggctcaccttcgggtgggcctttctgcgtttata",
            "The standard BioBrick double terminator: B0010 followed by B0012, giving near-complete transcriptional termination.")
        { RegistryId = "BBa_B0015" },

        // ---------------- Operators and assembly parts ----------------

        new BioPart("lacO", PartType.Operator, "aattgtgagcggataacaatt",
            "The symmetric lac operator O1 bound by the LacI tetramer; binding is relieved by IPTG.")
        { RegulatedBy = "LacI", IsRepressible = true, IsInducible = true, Inducer = "IPTG", RegistryId = "BBa_R0063" },

        new BioPart("tetO", PartType.Operator, "tccctatcagtgatagaga",
            "The tetO2 operator bound by the TetR dimer; binding is relieved by anhydrotetracycline.")
        { RegulatedBy = "TetR", IsRepressible = true, IsInducible = true, Inducer = "aTc", RegistryId = "BBa_R0061" },

        new BioPart("scar", PartType.Scar, "tactagag",
            "The 8 bp mixed site left between two parts after BioBrick RFC[10] assembly (EcoRI/XbaI/SpeI/PstI prefix-suffix ligation).")
        { RelativeStrength = 0.0 },

        new BioPart("spacer20", PartType.Spacer, "gatccgtacgatcgtaagct",
            "approximate/representative sequence: a neutral 20 bp spacer with balanced GC content and no internal BioBrick restriction sites, used to insulate adjacent parts.")
        { RelativeStrength = 0.0 },
    ];
}
