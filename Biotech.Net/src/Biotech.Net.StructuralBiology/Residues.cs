namespace Biotech.Net.StructuralBiology;

/// <summary>Residue-name dictionaries for amino acids, nucleotides and water.</summary>
public static class ResidueNames
{
    private static readonly Dictionary<string, char> ThreeToOne = new(StringComparer.OrdinalIgnoreCase)
    {
        ["ALA"] = 'A', ["ARG"] = 'R', ["ASN"] = 'N', ["ASP"] = 'D', ["CYS"] = 'C', ["GLN"] = 'Q', ["GLU"] = 'E',
        ["GLY"] = 'G', ["HIS"] = 'H', ["ILE"] = 'I', ["LEU"] = 'L', ["LYS"] = 'K', ["MET"] = 'M', ["PHE"] = 'F',
        ["PRO"] = 'P', ["SER"] = 'S', ["THR"] = 'T', ["TRP"] = 'W', ["TYR"] = 'Y', ["VAL"] = 'V',
        ["SEC"] = 'U', ["PYL"] = 'O', ["MSE"] = 'M', ["HID"] = 'H', ["HIE"] = 'H', ["HIP"] = 'H', ["HSD"] = 'H',
        ["HSE"] = 'H', ["HSP"] = 'H', ["CYX"] = 'C', ["ASX"] = 'B', ["GLX"] = 'Z', ["UNK"] = 'X',
    };

    private static readonly HashSet<string> Nucleotides = new(StringComparer.OrdinalIgnoreCase)
    { "A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU", "I", "DI" };

    private static readonly HashSet<string> Waters = new(StringComparer.OrdinalIgnoreCase) { "HOH", "WAT", "H2O", "DOD", "SOL", "TIP3" };

    /// <summary>
    /// Maximum solvent accessibility (Å²) of each amino acid in a Gly-X-Gly tripeptide, the theoretical values of
    /// Tien et al. 2013 (PLoS ONE 8:e80635), used to compute relative solvent accessibility.
    /// </summary>
    public static readonly IReadOnlyDictionary<char, double> MaxAsaTien2013 = new Dictionary<char, double>
    {
        ['A'] = 129, ['R'] = 274, ['N'] = 195, ['D'] = 193, ['C'] = 167, ['Q'] = 225, ['E'] = 223, ['G'] = 104,
        ['H'] = 224, ['I'] = 197, ['L'] = 201, ['K'] = 236, ['M'] = 224, ['F'] = 240, ['P'] = 159, ['S'] = 155,
        ['T'] = 172, ['W'] = 285, ['Y'] = 263, ['V'] = 174,
    };

    /// <summary>One-letter code of an amino-acid residue name, or <c>'X'</c> when unknown.</summary>
    public static char ToOneLetter(string residueName) => ThreeToOne.TryGetValue(residueName.Trim(), out char c) ? c : 'X';

    /// <summary>True for the standard (and common modified) amino-acid names.</summary>
    public static bool IsAminoAcid(string residueName) => ThreeToOne.ContainsKey(residueName.Trim()) && !residueName.Trim().Equals("UNK", StringComparison.OrdinalIgnoreCase);

    /// <summary>True for nucleotide residue names (A, C, G, U, DA, DC, DG, DT, ...).</summary>
    public static bool IsNucleotide(string residueName) => Nucleotides.Contains(residueName.Trim());

    /// <summary>True for water.</summary>
    public static bool IsWater(string residueName) => Waters.Contains(residueName.Trim());
}
