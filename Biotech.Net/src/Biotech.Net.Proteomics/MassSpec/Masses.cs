using System.Globalization;
using System.Text;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>Physical constants for peptide mass spectrometry (CODATA / IUPAC values).</summary>
public static class MassConstants
{
    /// <summary>Proton mass, 1.007276466621 Da.</summary>
    public const double Proton = 1.007276466621;

    /// <summary>Hydrogen atom (¹H), 1.00782503207 Da.</summary>
    public const double Hydrogen = 1.00782503207;

    /// <summary>Monoisotopic water, 18.0105646863 Da.</summary>
    public const double WaterMonoisotopic = 18.0105646863;

    /// <summary>Average water, 18.01528 Da.</summary>
    public const double WaterAverage = 18.01528;

    /// <summary>Monoisotopic ammonia, 17.0265491015 Da.</summary>
    public const double Ammonia = 17.0265491015;

    /// <summary>Monoisotopic carbon monoxide, 27.9949146221 Da.</summary>
    public const double CarbonMonoxide = 27.9949146221;

    /// <summary>¹³C − ¹²C mass difference, used as the isotope spacing, 1.0033548378 Da.</summary>
    public const double C13Spacing = 1.0033548378;

    /// <summary>m/z of an ion of neutral mass <paramref name="neutralMass"/> carrying <paramref name="charge"/> protons.</summary>
    public static double ToMz(double neutralMass, int charge)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(charge, 1);
        return (neutralMass + charge * Proton) / charge;
    }

    /// <summary>Neutral mass of a protonated ion observed at <paramref name="mz"/> with <paramref name="charge"/>.</summary>
    public static double ToNeutralMass(double mz, int charge)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(charge, 1);
        return mz * charge - charge * Proton;
    }

    /// <summary>Mass error of <paramref name="observed"/> relative to <paramref name="theoretical"/> in ppm.</summary>
    public static double PpmError(double observed, double theoretical) => (observed - theoretical) / theoretical * 1e6;
}

/// <summary>
/// Residue (in-chain, water removed) masses of the amino acids, monoisotopic (Unimod) and average, plus their
/// elemental compositions. Tables are indexed by ASCII byte and accept either case.
/// </summary>
public static class ResidueMasses
{
    private static readonly double[] Mono = new double[256];
    private static readonly double[] Avg = new double[256];
    private static readonly string?[] Formulas = new string?[256];

    static ResidueMasses()
    {
        Array.Fill(Mono, double.NaN);
        Array.Fill(Avg, double.NaN);
        (char, double, double, string)[] t =
        [
            ('G', 57.021463720, 57.05132, "C2H3NO"), ('A', 71.037113785, 71.0779, "C3H5NO"),
            ('S', 87.032028405, 87.0773, "C3H5NO2"), ('P', 97.052763850, 97.11518, "C5H7NO"),
            ('V', 99.068413915, 99.13106, "C5H9NO"), ('T', 101.047678470, 101.10388, "C4H7NO2"),
            ('C', 103.009184505, 103.1429, "C3H5NOS"), ('L', 113.084063980, 113.15764, "C6H11NO"),
            ('I', 113.084063980, 113.15764, "C6H11NO"), ('N', 114.042927470, 114.10264, "C4H6N2O2"),
            ('D', 115.026943065, 115.0874, "C4H5NO3"), ('Q', 128.058577540, 128.12922, "C5H8N2O2"),
            ('K', 128.094963050, 128.17228, "C6H12N2O"), ('E', 129.042593135, 129.11398, "C5H7NO3"),
            ('M', 131.040484645, 131.19606, "C5H9NOS"), ('H', 137.058911875, 137.13928, "C6H7N3O"),
            ('F', 147.068413980, 147.17386, "C9H9NO"), ('R', 156.101111050, 156.18568, "C6H12N4O"),
            ('Y', 163.063328600, 163.17326, "C9H9NO2"), ('W', 186.079312980, 186.2099, "C11H10N2O"),
            ('U', 150.953633405, 150.0379, "C3H5NOSe"), ('O', 237.147726925, 237.29816, "C12H19N3O2"),
        ];
        foreach (var (c, m, a, f) in t)
        {
            Mono[c] = m; Mono[char.ToLowerInvariant(c)] = m;
            Avg[c] = a; Avg[char.ToLowerInvariant(c)] = a;
            Formulas[c] = f; Formulas[char.ToLowerInvariant(c)] = f;
        }
    }

    /// <summary>Monoisotopic residue mass (NaN for unknown residues).</summary>
    public static double Monoisotopic(byte residue) => Mono[residue];

    /// <summary>Monoisotopic residue mass (NaN for unknown residues).</summary>
    public static double Monoisotopic(char residue) => residue < 256 ? Mono[residue] : double.NaN;

    /// <summary>Average residue mass (NaN for unknown residues).</summary>
    public static double Average(byte residue) => Avg[residue];

    /// <summary>Average residue mass (NaN for unknown residues).</summary>
    public static double Average(char residue) => residue < 256 ? Avg[residue] : double.NaN;

    /// <summary>Elemental composition of the residue, e.g. <c>C3H5NO</c> for alanine; null when unknown.</summary>
    public static string? Formula(char residue) => residue < 256 ? Formulas[residue] : null;

    /// <summary>Neutral mass of an unmodified peptide (residues + water). Throws on non-standard residues.</summary>
    public static double PeptideMass(ReadOnlySpan<byte> residues, bool monoisotopic = true)
    {
        var table = monoisotopic ? Mono : Avg;
        double m = monoisotopic ? MassConstants.WaterMonoisotopic : MassConstants.WaterAverage;
        foreach (byte b in residues)
        {
            double r = table[b];
            if (double.IsNaN(r)) throw new ArgumentException($"Residue '{(char)b}' has no defined mass.");
            m += r;
        }
        return m;
    }

    /// <summary>Neutral mass of an unmodified peptide given as a string.</summary>
    public static double PeptideMass(string residues, bool monoisotopic = true) => PeptideMass(Encoding.ASCII.GetBytes(residues), monoisotopic);
}

/// <summary>An elemental composition such as <c>C34H53N7O15</c>.</summary>
public sealed class ChemicalFormula : IEquatable<ChemicalFormula>
{
    private readonly SortedDictionary<string, int> _counts = new(StringComparer.Ordinal);

    /// <summary>Creates an empty formula.</summary>
    public ChemicalFormula() { }

    /// <summary>Creates a formula from element counts.</summary>
    public ChemicalFormula(IEnumerable<KeyValuePair<string, int>> counts)
    {
        foreach (var (e, n) in counts) Add(e, n);
    }

    /// <summary>Element → atom count (zero counts are removed).</summary>
    public IReadOnlyDictionary<string, int> Counts => _counts;

    /// <summary>Atoms of one element.</summary>
    public int this[string element] => _counts.TryGetValue(element, out int n) ? n : 0;

    /// <summary>Adds (or with a negative count removes) atoms in place.</summary>
    public ChemicalFormula Add(string element, int count)
    {
        if (!Isotopes.Known(element)) throw new ArgumentException($"Unknown element '{element}'.");
        int n = this[element] + count;
        if (n == 0) _counts.Remove(element); else _counts[element] = n;
        return this;
    }

    /// <summary>Adds another formula in place.</summary>
    public ChemicalFormula Add(ChemicalFormula other)
    {
        foreach (var (e, n) in other._counts) Add(e, n);
        return this;
    }

    /// <summary>Parses Hill-style formulas with optional counts, e.g. <c>C2H5OH</c> or <c>C3H5NOSe</c>. Negative counts allowed (<c>H-1</c>).</summary>
    public static ChemicalFormula Parse(string formula)
    {
        ArgumentNullException.ThrowIfNull(formula);
        var f = new ChemicalFormula();
        int i = 0;
        while (i < formula.Length)
        {
            if (!char.IsUpper(formula[i])) throw new FormatException($"Bad formula '{formula}' at position {i}.");
            int s = i++;
            while (i < formula.Length && char.IsLower(formula[i])) i++;
            string element = formula[s..i];
            int ns = i;
            if (i < formula.Length && formula[i] == '-') i++;
            while (i < formula.Length && char.IsDigit(formula[i])) i++;
            int count = i == ns ? 1 : int.Parse(formula.AsSpan(ns, i - ns), CultureInfo.InvariantCulture);
            f.Add(element, count);
        }
        return f;
    }

    /// <summary>Composition of an unmodified peptide: the residues plus one water.</summary>
    public static ChemicalFormula FromPeptide(string residues)
    {
        ArgumentNullException.ThrowIfNull(residues);
        var f = Parse("H2O");
        foreach (char c in residues)
            f.Add(Parse(ResidueMasses.Formula(c) ?? throw new ArgumentException($"Residue '{c}' has no defined composition.")));
        return f;
    }

    /// <summary>
    /// The averagine composition (Senko et al., JASMS 6:229, 1995: C4.9384 H7.7583 N1.3577 O1.4773 S0.0417 per
    /// 111.1254 Da) scaled to <paramref name="averageMass"/>, rounded, with hydrogens adjusted to absorb the
    /// rounding error.
    /// </summary>
    public static ChemicalFormula Averagine(double averageMass)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(averageMass);
        double units = averageMass / 111.1254;
        var f = new ChemicalFormula()
            .Add("C", (int)Math.Round(4.9384 * units))
            .Add("N", (int)Math.Round(1.3577 * units))
            .Add("O", (int)Math.Round(1.4773 * units));
        int s = (int)Math.Round(0.0417 * units);
        if (s > 0) f.Add("S", s);
        double partial = f.AverageMass;
        int h = (int)Math.Round((averageMass - partial) / Isotopes.AverageMass("H"));
        if (h > 0) f.Add("H", h);
        return f;
    }

    /// <summary>Monoisotopic mass (every element at its most abundant isotope, which for C/H/N/O/S/P is also the lightest; Se uses ⁸⁰Se as Unimod does).</summary>
    public double MonoisotopicMass
    {
        get
        {
            double m = 0;
            foreach (var (e, n) in _counts) m += n * Isotopes.MonoisotopicMass(e);
            return m;
        }
    }

    /// <summary>Average (natural-abundance weighted) mass.</summary>
    public double AverageMass
    {
        get
        {
            double m = 0;
            foreach (var (e, n) in _counts) m += n * Isotopes.AverageMass(e);
            return m;
        }
    }

    /// <summary>Hill notation: C, H, then the other elements alphabetically.</summary>
    public override string ToString()
    {
        var sb = new StringBuilder();
        void Emit(string e) { int n = this[e]; if (n == 0) return; sb.Append(e); if (n != 1) sb.Append(n.ToString(CultureInfo.InvariantCulture)); }
        Emit("C"); Emit("H");
        foreach (var e in _counts.Keys) if (e is not "C" and not "H") Emit(e);
        return sb.ToString();
    }

    /// <inheritdoc/>
    public bool Equals(ChemicalFormula? other)
    {
        if (other is null || other._counts.Count != _counts.Count) return false;
        foreach (var (e, n) in _counts) if (other[e] != n) return false;
        return true;
    }

    /// <inheritdoc/>
    public override bool Equals(object? obj) => obj is ChemicalFormula f && Equals(f);

    /// <inheritdoc/>
    public override int GetHashCode() => ToString().GetHashCode(StringComparison.Ordinal);
}

/// <summary>Stable isotopes (mass, natural abundance) of the elements found in peptides and common modifications (IUPAC 2013).</summary>
public static class Isotopes
{
    private static readonly Dictionary<string, (double Mass, double Abundance)[]> Table = new(StringComparer.Ordinal)
    {
        ["H"] = [(1.00782503207, 0.999885), (2.0141017778, 0.000115)],
        ["C"] = [(12.0, 0.9893), (13.0033548378, 0.0107)],
        ["N"] = [(14.0030740048, 0.99636), (15.0001088982, 0.00364)],
        ["O"] = [(15.99491461956, 0.99757), (16.99913170, 0.00038), (17.9991610, 0.00205)],
        ["S"] = [(31.97207100, 0.9499), (32.97145876, 0.0075), (33.96786690, 0.0425), (35.96708076, 0.0001)],
        ["P"] = [(30.97376163, 1.0)],
        ["Se"] = [(73.9224764, 0.0089), (75.9192136, 0.0937), (76.9199140, 0.0763), (77.9173091, 0.2377), (79.9165213, 0.4961), (81.9166994, 0.0873)],
        ["Na"] = [(22.9897692809, 1.0)],
        ["K"] = [(38.96370668, 0.932581), (39.96399848, 0.000117), (40.96182576, 0.067302)],
        ["Cl"] = [(34.96885268, 0.7576), (36.96590259, 0.2424)],
        ["Br"] = [(78.9183371, 0.5069), (80.9162906, 0.4931)],
        ["F"] = [(18.99840322, 1.0)],
        ["I"] = [(126.904473, 1.0)],
        ["Fe"] = [(53.9396105, 0.05845), (55.9349375, 0.91754), (56.9353940, 0.02119), (57.9332756, 0.00282)],
    };

    /// <summary>True when isotope data exists for <paramref name="element"/>.</summary>
    public static bool Known(string element) => Table.ContainsKey(element);

    /// <summary>The (mass, abundance) pairs of an element, lightest first.</summary>
    public static IReadOnlyList<(double Mass, double Abundance)> Of(string element)
        => Table.TryGetValue(element, out var v) ? v : throw new ArgumentException($"Unknown element '{element}'.");

    /// <summary>Mass of the most abundant isotope (the monoisotopic contribution by proteomics convention).</summary>
    public static double MonoisotopicMass(string element)
    {
        var iso = Of(element);
        var best = iso[0];
        for (int i = 1; i < iso.Count; i++) if (iso[i].Abundance > best.Abundance) best = iso[i];
        return best.Mass;
    }

    /// <summary>Abundance-weighted mass.</summary>
    public static double AverageMass(string element)
    {
        double m = 0, a = 0;
        foreach (var (mass, ab) in Of(element)) { m += mass * ab; a += ab; }
        return m / a;
    }
}
