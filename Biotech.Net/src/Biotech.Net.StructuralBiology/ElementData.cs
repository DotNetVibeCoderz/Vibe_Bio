namespace Biotech.Net.StructuralBiology;

/// <summary>
/// Per-element constants used by structure analysis: average atomic mass (IUPAC), van der Waals radius used for
/// solvent accessibility (the set Biopython's <c>Bio.PDB.SASA</c> uses, after Bondi 1964 / Rowland &amp; Taylor 1996 /
/// Mantina 2009) and covalent radius (Cordero et al. 2008) used for bond perception.
/// </summary>
public static class ElementData
{
    private static readonly Dictionary<string, (double Mass, double Vdw, double Covalent)> Table = new(StringComparer.OrdinalIgnoreCase)
    {
        ["H"] = (1.008, 1.20, 0.31), ["D"] = (2.014, 1.20, 0.31), ["HE"] = (4.0026, 1.40, 0.28),
        ["C"] = (12.0107, 1.70, 0.76), ["N"] = (14.0067, 1.55, 0.71), ["O"] = (15.9994, 1.52, 0.66),
        ["F"] = (18.9984, 1.47, 0.57), ["NA"] = (22.9898, 2.27, 1.66), ["MG"] = (24.305, 1.73, 1.41),
        ["P"] = (30.9738, 1.80, 1.07), ["S"] = (32.065, 1.80, 1.05), ["CL"] = (35.453, 1.75, 1.02),
        ["K"] = (39.0983, 2.75, 2.03), ["CA"] = (40.078, 2.31, 1.76), ["MN"] = (54.938, 1.61, 1.39),
        ["FE"] = (55.845, 1.56, 1.32), ["CO"] = (58.933, 1.52, 1.26), ["NI"] = (58.6934, 1.63, 1.24),
        ["CU"] = (63.546, 1.40, 1.32), ["ZN"] = (65.38, 1.39, 1.22), ["SE"] = (78.96, 1.90, 1.20),
        ["BR"] = (79.904, 1.85, 1.20), ["I"] = (126.904, 1.98, 1.39), ["CD"] = (112.411, 1.58, 1.44),
        ["HG"] = (200.59, 1.55, 1.32), ["B"] = (10.811, 1.92, 0.84), ["SI"] = (28.0855, 2.10, 1.11),
    };

    /// <summary>Whether the element symbol is known to this table.</summary>
    public static bool IsKnown(string element) => Table.ContainsKey(element);

    /// <summary>Average atomic mass in daltons (12.0107 for unknown elements, i.e. treated as carbon).</summary>
    public static double Mass(string element) => Table.TryGetValue(element, out var e) ? e.Mass : 12.0107;

    /// <summary>Van der Waals radius in Å (1.80 for unknown elements).</summary>
    public static double VdwRadius(string element) => Table.TryGetValue(element, out var e) ? e.Vdw : 1.80;

    /// <summary>Covalent radius in Å (0.77 for unknown elements).</summary>
    public static double CovalentRadius(string element) => Table.TryGetValue(element, out var e) ? e.Covalent : 0.77;

    /// <summary>
    /// Guesses an element from a PDB atom name when the element columns (77–78) are blank. Follows the PDB
    /// convention that a two-letter element starts in column 13 while one-letter elements start in column 14.
    /// </summary>
    /// <param name="atomNameField">The raw 4-character name field (columns 13–16).</param>
    /// <param name="isHetero">HETATM records may hold metal ions such as <c>CA</c> or <c>ZN</c>.</param>
    public static string GuessFromAtomName(string atomNameField, bool isHetero)
    {
        string field = atomNameField.PadRight(4);
        string trimmed = field.Trim();
        if (trimmed.Length == 0) return "X";
        if (field[0] != ' ' && !char.IsDigit(field[0]))
        {
            // Name starts in column 13: two-letter element for HETATM (FE, ZN, CL, ...), hydrogens like HG21 otherwise.
            string two = field[..2].ToUpperInvariant();
            if (field[0] is 'H' or 'h' && !isHetero) return "H";
            if (Table.ContainsKey(two)) return two;
            return field[0].ToString().ToUpperInvariant();
        }
        foreach (char c in trimmed)
            if (char.IsLetter(c)) return char.ToUpperInvariant(c).ToString();
        return "X";
    }

    /// <summary>Normalises an element symbol to upper case without surrounding blanks (<c>" Fe"</c> → <c>"FE"</c>).</summary>
    public static string Normalize(string symbol) => symbol.Trim().ToUpperInvariant();
}
