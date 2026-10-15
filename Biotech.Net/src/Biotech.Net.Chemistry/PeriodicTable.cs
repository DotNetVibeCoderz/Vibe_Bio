using System.Collections.Frozen;

namespace Biotech.Net.Chemistry;

/// <summary>Physical and chemical data for one chemical element.</summary>
/// <param name="AtomicNumber">Number of protons (Z).</param>
/// <param name="Symbol">IUPAC symbol, e.g. <c>Cl</c>.</param>
/// <param name="AverageMass">Standard atomic weight in daltons (natural isotopic abundance).</param>
/// <param name="MonoisotopicMass">Mass of the most abundant isotope in daltons.</param>
/// <param name="MostCommonIsotope">Mass number of the most abundant isotope.</param>
/// <param name="Valences">Allowed (neutral) valences, lowest first; used for implicit-hydrogen computation.</param>
/// <param name="Electronegativity">Pauling electronegativity (0 when undefined, e.g. noble gases).</param>
/// <param name="CovalentRadius">Single-bond covalent radius in ångström (Cordero et al. 2008).</param>
/// <param name="OuterElectrons">Number of valence-shell electrons.</param>
public sealed record Element(int AtomicNumber, string Symbol, double AverageMass, double MonoisotopicMass,
    int MostCommonIsotope, int[] Valences, double Electronegativity, double CovalentRadius, int OuterElectrons)
{
    /// <summary>The lowest allowed valence, or 0 when the element has none (noble gases, metals without a default).</summary>
    public int DefaultValence => Valences.Length > 0 ? Valences[0] : 0;

    /// <summary>True for B, C, N, O, P, S, F, Cl, Br, I — the SMILES "organic subset" that may be written without brackets.</summary>
    public bool IsOrganicSubset => AtomicNumber is 5 or 6 or 7 or 8 or 9 or 15 or 16 or 17 or 35 or 53;

    /// <summary>True for metals (alkali, alkaline earth, transition, post-transition, lanthanides, actinides).</summary>
    public bool IsMetal => AtomicNumber is 3 or 4 or (>= 11 and <= 13) or (>= 19 and <= 31) or (>= 37 and <= 50)
        or (>= 55 and <= 84) or >= 87;

    public override string ToString() => $"{Symbol} (Z={AtomicNumber}, {AverageMass:F3} Da)";
}

/// <summary>
/// The periodic table for elements 1–92. Masses and covalent radii follow the RDKit/Blue Obelisk data set
/// (IUPAC standard atomic weights; Cordero 2008 radii); electronegativities are Pauling values.
/// </summary>
public static class PeriodicTable
{
    private static readonly Element[] _elements =
    [
        new(0, "*", 0, 0, 0, [], 0, 0, 0),
        new(1, "H", 1.008, 1.007825032, 1, [1], 2.2, 0.31, 1),
        new(2, "He", 4.003, 4.002603254, 4, [], 0, 0.28, 2),
        new(3, "Li", 6.941, 7.01600455, 7, [1], 0.98, 1.28, 1),
        new(4, "Be", 9.012, 9.0121822, 9, [2], 1.57, 0.96, 2),
        new(5, "B", 10.812, 11.0093054, 11, [3], 2.04, 0.84, 3),
        new(6, "C", 12.011, 12.0, 12, [4], 2.55, 0.76, 4),
        new(7, "N", 14.007, 14.003074, 14, [3, 5], 3.04, 0.71, 5),
        new(8, "O", 15.999, 15.99491462, 16, [2], 3.44, 0.66, 6),
        new(9, "F", 18.998, 18.99840322, 19, [1], 3.98, 0.57, 7),
        new(10, "Ne", 20.18, 19.99244018, 20, [], 0, 0.58, 8),
        new(11, "Na", 22.99, 22.98976928, 23, [1], 0.93, 1.66, 1),
        new(12, "Mg", 24.305, 23.9850417, 24, [2], 1.31, 1.41, 2),
        new(13, "Al", 26.982, 26.98153863, 27, [3], 1.61, 1.21, 3),
        new(14, "Si", 28.086, 27.97692653, 28, [4], 1.9, 1.11, 4),
        new(15, "P", 30.974, 30.97376163, 31, [3, 5], 2.19, 1.07, 5),
        new(16, "S", 32.067, 31.972071, 32, [2, 4, 6], 2.58, 1.05, 6),
        new(17, "Cl", 35.453, 34.96885268, 35, [1], 3.16, 1.02, 7),
        new(18, "Ar", 39.948, 39.96238312, 40, [], 0, 1.06, 8),
        new(19, "K", 39.098, 38.96370668, 39, [1], 0.82, 2.03, 1),
        new(20, "Ca", 40.078, 39.96259098, 40, [2], 1.0, 1.76, 2),
        new(21, "Sc", 44.956, 44.9559119, 45, [], 1.36, 1.7, 3),
        new(22, "Ti", 47.867, 47.9479463, 48, [], 1.54, 1.6, 4),
        new(23, "V", 50.944, 50.9439595, 51, [], 1.63, 1.52, 5),
        new(24, "Cr", 51.996, 51.9405075, 52, [], 1.66, 1.39, 6),
        new(25, "Mn", 54.938, 54.9380451, 55, [], 1.55, 1.39, 7),
        new(26, "Fe", 55.845, 55.9349375, 56, [], 1.83, 1.32, 8),
        new(27, "Co", 58.933, 58.933195, 59, [], 1.88, 1.26, 9),
        new(28, "Ni", 58.693, 57.9353429, 58, [], 1.91, 1.24, 10),
        new(29, "Cu", 63.546, 62.9295975, 63, [], 1.9, 1.32, 11),
        new(30, "Zn", 65.39, 63.9291422, 64, [], 1.65, 1.22, 2),
        new(31, "Ga", 69.723, 68.9255736, 69, [3], 1.81, 1.22, 3),
        new(32, "Ge", 72.61, 73.9211778, 74, [4], 2.01, 1.2, 4),
        new(33, "As", 74.922, 74.9215965, 75, [3, 5], 2.18, 1.19, 5),
        new(34, "Se", 78.96, 79.9165213, 80, [2, 4, 6], 2.55, 1.2, 6),
        new(35, "Br", 79.904, 78.9183371, 79, [1], 2.96, 1.2, 7),
        new(36, "Kr", 83.8, 83.911507, 84, [], 3.0, 1.16, 8),
        new(37, "Rb", 85.468, 84.91178974, 85, [1], 0.82, 2.2, 1),
        new(38, "Sr", 87.62, 87.9056121, 88, [2], 0.95, 1.95, 2),
        new(39, "Y", 88.906, 88.9058483, 89, [], 1.22, 1.9, 3),
        new(40, "Zr", 91.224, 89.9047044, 90, [], 1.33, 1.75, 4),
        new(41, "Nb", 92.906, 92.9063781, 93, [], 1.6, 1.64, 5),
        new(42, "Mo", 95.94, 97.9054082, 98, [], 2.16, 1.54, 6),
        new(43, "Tc", 98.0, 96.906365, 97, [], 1.9, 1.47, 7),
        new(44, "Ru", 101.07, 101.9043493, 102, [], 2.2, 1.46, 8),
        new(45, "Rh", 102.906, 102.905504, 103, [], 2.28, 1.42, 9),
        new(46, "Pd", 106.42, 105.903486, 106, [], 2.2, 1.39, 10),
        new(47, "Ag", 107.868, 106.905097, 107, [], 1.93, 1.45, 11),
        new(48, "Cd", 112.412, 113.9033585, 114, [], 1.69, 1.44, 2),
        new(49, "In", 114.818, 114.903878, 115, [3], 1.78, 1.42, 3),
        new(50, "Sn", 118.711, 119.9021947, 120, [2, 4], 1.96, 1.39, 4),
        new(51, "Sb", 121.76, 120.9038157, 121, [3, 5], 2.05, 1.39, 5),
        new(52, "Te", 127.6, 129.9062244, 130, [2, 4, 6], 2.1, 1.38, 6),
        new(53, "I", 126.904, 126.904473, 127, [1, 3, 5], 2.66, 1.39, 7),
        new(54, "Xe", 131.29, 131.9041535, 132, [0, 2, 4, 6], 2.6, 1.4, 8),
        new(55, "Cs", 132.905, 132.9054519, 133, [1], 0.79, 2.44, 1),
        new(56, "Ba", 137.328, 137.9052472, 138, [2], 0.89, 2.15, 2),
        new(57, "La", 138.906, 138.9063533, 139, [], 1.1, 2.07, 3),
        new(58, "Ce", 140.116, 139.9054387, 140, [], 1.12, 2.04, 4),
        new(59, "Pr", 140.908, 140.9076528, 141, [], 1.13, 2.03, 3),
        new(60, "Nd", 144.24, 141.9077233, 142, [], 1.14, 2.01, 4),
        new(61, "Pm", 145.0, 144.912749, 145, [], 1.13, 1.99, 5),
        new(62, "Sm", 150.36, 151.9197324, 152, [], 1.17, 1.98, 6),
        new(63, "Eu", 151.964, 152.9212303, 153, [], 1.2, 1.98, 7),
        new(64, "Gd", 157.25, 157.9241039, 158, [], 1.2, 1.96, 8),
        new(65, "Tb", 158.925, 158.9253468, 159, [], 1.1, 1.94, 9),
        new(66, "Dy", 162.5, 163.9291748, 164, [], 1.22, 1.92, 10),
        new(67, "Ho", 164.93, 164.9303221, 165, [], 1.23, 1.92, 11),
        new(68, "Er", 167.26, 165.9302931, 166, [], 1.24, 1.89, 12),
        new(69, "Tm", 168.934, 168.9342133, 169, [], 1.25, 1.9, 13),
        new(70, "Yb", 173.04, 173.9388621, 174, [], 1.1, 1.87, 14),
        new(71, "Lu", 174.967, 174.9407718, 175, [], 1.27, 1.87, 15),
        new(72, "Hf", 178.49, 179.94655, 180, [], 1.3, 1.75, 4),
        new(73, "Ta", 180.948, 180.9479958, 181, [], 1.5, 1.7, 5),
        new(74, "W", 183.84, 183.9509312, 184, [], 2.36, 1.62, 6),
        new(75, "Re", 186.207, 186.9557531, 187, [], 1.9, 1.51, 7),
        new(76, "Os", 190.23, 191.9614807, 192, [], 2.2, 1.44, 8),
        new(77, "Ir", 192.217, 192.9629264, 193, [], 2.2, 1.41, 9),
        new(78, "Pt", 195.078, 194.9647911, 195, [], 2.28, 1.36, 10),
        new(79, "Au", 196.967, 196.9665687, 197, [], 2.54, 1.36, 11),
        new(80, "Hg", 200.59, 201.970643, 202, [], 2.0, 1.32, 2),
        new(81, "Tl", 204.383, 204.9744275, 205, [], 1.62, 1.45, 3),
        new(82, "Pb", 207.2, 207.9766521, 208, [2, 4], 2.33, 1.46, 4),
        new(83, "Bi", 208.98, 208.9803987, 209, [3, 5], 2.02, 1.48, 5),
        new(84, "Po", 209.0, 208.9824304, 209, [2, 4, 6], 2.0, 1.4, 6),
        new(85, "At", 210.0, 209.987148, 210, [1, 3, 5], 2.2, 1.5, 7),
        new(86, "Rn", 222.0, 222.0175706, 222, [], 2.2, 1.5, 8),
        new(87, "Fr", 223.0, 223.0197359, 223, [1], 0.7, 2.6, 1),
        new(88, "Ra", 226.0, 226.0254026, 226, [2], 0.9, 2.2, 2),
        new(89, "Ac", 227.0, 227.0277521, 227, [], 1.1, 2.15, 3),
        new(90, "Th", 232.038, 232.0380553, 232, [], 1.3, 2.06, 4),
        new(91, "Pa", 231.036, 231.035884, 231, [], 1.5, 2.0, 3),
        new(92, "U", 238.029, 238.0507882, 238, [], 1.38, 1.96, 4),
    ];

    private static readonly FrozenDictionary<string, Element> _bySymbol =
        _elements.Skip(1).ToFrozenDictionary(e => e.Symbol, StringComparer.Ordinal);

    private static readonly FrozenDictionary<string, Element> _bySymbolIgnoreCase =
        _elements.Skip(1).ToFrozenDictionary(e => e.Symbol, StringComparer.OrdinalIgnoreCase);

    // Exact masses of isotopes that commonly appear in labelled compounds (AME2016).
    private static readonly FrozenDictionary<(int Z, int A), double> _isotopeMasses = new Dictionary<(int, int), double>
    {
        [(1, 1)] = 1.00782503207, [(1, 2)] = 2.0141017778, [(1, 3)] = 3.0160492777,
        [(6, 11)] = 11.0114336, [(6, 12)] = 12.0, [(6, 13)] = 13.0033548378, [(6, 14)] = 14.003241989,
        [(7, 14)] = 14.0030740048, [(7, 15)] = 15.0001088982,
        [(8, 16)] = 15.99491461956, [(8, 17)] = 16.99913170, [(8, 18)] = 17.9991610,
        [(9, 18)] = 18.0009380, [(9, 19)] = 18.99840322,
        [(15, 32)] = 31.97390727, [(16, 32)] = 31.97207100, [(16, 34)] = 33.96786690, [(16, 35)] = 34.96903216,
        [(17, 35)] = 34.96885268, [(17, 37)] = 36.96590259, [(35, 79)] = 78.9183371, [(35, 81)] = 80.9162906,
        [(53, 123)] = 122.905589, [(53, 125)] = 124.9046302, [(53, 127)] = 126.904473, [(53, 131)] = 130.9061246,
    }.ToFrozenDictionary();

    /// <summary>Number of real elements in the table.</summary>
    public static int Count => _elements.Length - 1;

    /// <summary>The element with atomic number <paramref name="atomicNumber"/> (0 is the dummy/wildcard atom <c>*</c>).</summary>
    public static Element Get(int atomicNumber)
    {
        if ((uint)atomicNumber >= (uint)_elements.Length) throw new ArgumentOutOfRangeException(nameof(atomicNumber));
        return _elements[atomicNumber];
    }

    /// <summary>The element with symbol <paramref name="symbol"/> (case-sensitive, e.g. <c>Cl</c>).</summary>
    public static Element Get(string symbol)
        => TryGet(symbol, out var e) ? e : throw new KeyNotFoundException($"Unknown element symbol '{symbol}'.");

    /// <summary>Looks up an element by its case-sensitive symbol; <c>*</c> gives the dummy atom.</summary>
    public static bool TryGet(string symbol, out Element element)
    {
        if (symbol == "*") { element = _elements[0]; return true; }
        return _bySymbol.TryGetValue(symbol, out element!);
    }

    /// <summary>Looks up an element ignoring case (useful for PDB/MOL columns written in upper case).</summary>
    public static bool TryGetIgnoreCase(string symbol, out Element element)
        => _bySymbolIgnoreCase.TryGetValue(symbol.Trim(), out element!);

    /// <summary>Atomic number for <paramref name="symbol"/>, or −1 if unknown.</summary>
    public static int AtomicNumberOf(string symbol) => TryGet(symbol, out var e) ? e.AtomicNumber : -1;

    /// <summary>Exact mass of isotope <paramref name="massNumber"/> of element <paramref name="atomicNumber"/>.</summary>
    /// <remarks>Uses a table of common labelling isotopes; other isotopes fall back to the mass number itself.</remarks>
    public static double IsotopeMass(int atomicNumber, int massNumber)
    {
        if (_isotopeMasses.TryGetValue((atomicNumber, massNumber), out double m)) return m;
        var e = Get(atomicNumber);
        return massNumber == e.MostCommonIsotope ? e.MonoisotopicMass : massNumber;
    }

    /// <summary>All elements (Z = 1…92).</summary>
    public static IEnumerable<Element> All => _elements.Skip(1);
}
