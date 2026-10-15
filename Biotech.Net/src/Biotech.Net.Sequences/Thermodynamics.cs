using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>Buffer conditions for melting-temperature calculations (concentrations in mM unless stated).</summary>
public sealed record MeltingConditions
{
    /// <summary>Na⁺ (mM).</summary>
    public double Sodium { get; init; } = 50;

    /// <summary>K⁺ (mM).</summary>
    public double Potassium { get; init; }

    /// <summary>Tris buffer (mM); half of it counts as monovalent cation.</summary>
    public double Tris { get; init; }

    /// <summary>Mg²⁺ (mM).</summary>
    public double Magnesium { get; init; }

    /// <summary>dNTPs (mM); they chelate Mg²⁺.</summary>
    public double Dntps { get; init; }

    /// <summary>Concentration of the higher-concentration strand (nM).</summary>
    public double PrimerConcentration { get; init; } = 25;

    /// <summary>Concentration of the lower-concentration strand (nM).</summary>
    public double TargetConcentration { get; init; } = 25;

    /// <summary>
    /// Sodium-equivalent monovalent concentration in molar: Na + K + Tris/2, plus 120·√(Mg − dNTPs) when Mg²⁺
    /// exceeds dNTPs (von Ahsen et al. 2001).
    /// </summary>
    public double MonovalentEquivalentMolar
    {
        get
        {
            double mon = Sodium + Potassium + Tris / 2.0;
            if (Potassium + Magnesium + Tris + Dntps > 0 && Dntps < Magnesium) mon += 120 * Math.Sqrt(Magnesium - Dntps);
            return mon * 1e-3;
        }
    }
}

/// <summary>
/// Oligonucleotide melting temperature (°C) by the Wallace rule, GC-content formulas and the nearest-neighbour
/// model with the SantaLucia (1998) unified parameters. Defaults reproduce Biopython's <c>MeltingTemp</c>.
/// </summary>
public static class MeltingTemperature
{
    private const double R = 1.987; // cal/(K·mol)

    /// <summary>Wallace rule: 2·(A+T) + 4·(G+C). Only meaningful for primers of 14–20 nt.</summary>
    public static double Wallace(Sequence oligo)
    {
        var c = oligo.NucleotideCounts();
        return 2.0 * (c.A + c.T) + 4.0 * (c.G + c.C);
    }

    /// <summary>Basic GC formula: 64.9 + 41·(G+C − 16.4)/N (Marmur &amp; Doty 1962, as used by OligoCalc).</summary>
    public static double BasicGc(Sequence oligo)
    {
        var c = oligo.NucleotideCounts();
        if (oligo.Length == 0) return double.NaN;
        return 64.9 + 41.0 * (c.G + c.C - 16.4) / oligo.Length;
    }

    /// <summary>
    /// Salt-adjusted GC formula 81.5 + 0.41·%GC − 600/N + 16.6·log10[Na⁺] (Biopython <c>Tm_GC</c> default value set 7).
    /// </summary>
    public static double SaltAdjustedGc(Sequence oligo, MeltingConditions? conditions = null)
    {
        conditions ??= new MeltingConditions();
        var c = oligo.NucleotideCounts();
        double pgc = 100.0 * (c.G + c.C) / oligo.Length;
        return 81.5 + 0.41 * pgc - 600.0 / oligo.Length + 16.6 * Math.Log10(conditions.MonovalentEquivalentMolar);
    }

    // SantaLucia (1998) unified nearest-neighbour parameters: ΔH (kcal/mol), ΔS (cal/K·mol), keyed by the top-strand dinucleotide.
    private static readonly Dictionary<string, (double H, double S)> NearestNeighbour = new()
    {
        ["AA"] = (-7.9, -22.2), ["TT"] = (-7.9, -22.2),
        ["AT"] = (-7.2, -20.4),
        ["TA"] = (-7.2, -21.3),
        ["CA"] = (-8.5, -22.7), ["TG"] = (-8.5, -22.7),
        ["GT"] = (-8.4, -22.4), ["AC"] = (-8.4, -22.4),
        ["CT"] = (-7.8, -21.0), ["AG"] = (-7.8, -21.0),
        ["GA"] = (-8.2, -22.2), ["TC"] = (-8.2, -22.2),
        ["CG"] = (-10.6, -27.2),
        ["GC"] = (-9.8, -24.4),
        ["GG"] = (-8.0, -19.9), ["CC"] = (-8.0, -19.9),
    };

    /// <summary>
    /// Nearest-neighbour Tm of a perfectly matched duplex, SantaLucia (1998) unified parameters, terminal A·T / G·C
    /// initiation, symmetry correction for self-complementary oligos and the SantaLucia entropy salt correction
    /// ΔS += 0.368·(N−1)·ln[Na⁺eq]. Tm = 1000·ΔH / (ΔS + R·ln(Ct)) − 273.15, with Ct = C1 − C2/2 (or C1 when
    /// self-complementary). Matches Biopython <c>Tm_NN</c> (DNA_NN3, saltcorr=5).
    /// </summary>
    public static double NearestNeighbor(Sequence oligo, MeltingConditions? conditions = null, bool? selfComplementary = null)
    {
        conditions ??= new MeltingConditions();
        string s = oligo.ToString().ToUpperInvariant().Replace('U', 'T');
        if (s.Length < 2) throw new ArgumentException("Nearest-neighbour Tm needs at least two bases.", nameof(oligo));
        foreach (char ch in s)
            if (ch is not ('A' or 'C' or 'G' or 'T')) throw new ArgumentException($"Base '{ch}' is not supported by the nearest-neighbour model.", nameof(oligo));
        double dh = 0, ds = 0;
        // Initiation: one term per terminal base pair.
        foreach (char end in new[] { s[0], s[^1] })
        {
            if (end is 'A' or 'T') { dh += 2.3; ds += 4.1; }
            else { dh += 0.1; ds += -2.8; }
        }
        for (int i = 0; i + 1 < s.Length; i++)
        {
            var (h, e) = NearestNeighbour[s.Substring(i, 2)];
            dh += h; ds += e;
        }
        bool self = selfComplementary ?? s == oligo.ReverseComplement().ToString().ToUpperInvariant().Replace('U', 'T');
        double ct = (conditions.PrimerConcentration - conditions.TargetConcentration / 2.0) * 1e-9;
        if (self)
        {
            ct = conditions.PrimerConcentration * 1e-9;
            ds += -1.4;
        }
        ds += 0.368 * (s.Length - 1) * Math.Log(conditions.MonovalentEquivalentMolar);
        return 1000 * dh / (ds + R * Math.Log(ct)) - 273.15;
    }

    /// <summary>ΔG°37 (kcal/mol) of duplex formation by nearest neighbours (no salt correction).</summary>
    public static double FreeEnergy37(Sequence oligo)
    {
        string s = oligo.ToString().ToUpperInvariant().Replace('U', 'T');
        double dh = 0, ds = 0;
        foreach (char end in new[] { s[0], s[^1] })
        {
            if (end is 'A' or 'T') { dh += 2.3; ds += 4.1; } else { dh += 0.1; ds += -2.8; }
        }
        for (int i = 0; i + 1 < s.Length; i++) { var (h, e) = NearestNeighbour[s.Substring(i, 2)]; dh += h; ds += e; }
        return dh - 310.15 * ds / 1000.0;
    }
}

/// <summary>Molecular weight of nucleic acids (average masses, as Biopython <c>SeqUtils.molecular_weight</c>).</summary>
public static class NucleicAcidMass
{
    // Average masses of nucleoside monophosphates (Da) from Biopython IUPACData; one water is lost per phosphodiester bond.
    private const double Water = 18.0153;
    private static readonly Dictionary<char, double> DnaWeights = new() { ['A'] = 331.2218, ['C'] = 307.1971, ['G'] = 347.2212, ['T'] = 322.2085 };
    private static readonly Dictionary<char, double> RnaWeights = new() { ['A'] = 347.2212, ['C'] = 323.1965, ['G'] = 363.2206, ['U'] = 324.1813 };

    /// <summary>
    /// Molecular weight (Da) of a DNA or RNA strand: Σ monophosphate masses − (N−1)·H₂O (one more water for circular
    /// molecules). With <paramref name="doubleStranded"/> the complementary strand is added.
    /// </summary>
    /// <exception cref="ArgumentException">The sequence contains ambiguous bases.</exception>
    public static double MolecularWeight(Sequence sequence, bool doubleStranded = false, bool circular = false)
    {
        if (!sequence.Alphabet.IsNucleotide) throw new ArgumentException("Use Protein.MolecularWeight for proteins.", nameof(sequence));
        bool rna = sequence.Type == SequenceType.Rna;
        var table = rna ? RnaWeights : DnaWeights;
        double w = Strand(sequence, table, circular);
        if (doubleStranded) w += Strand(sequence.ReverseComplement(), table, circular);
        return w;
    }

    private static double Strand(Sequence s, Dictionary<char, double> table, bool circular)
    {
        if (s.Length == 0) return 0;
        double w = 0;
        foreach (byte b in s.Span)
        {
            char c = char.ToUpperInvariant((char)b);
            if (!table.TryGetValue(c, out double m)) throw new ArgumentException($"'{c}' has no defined mass; molecular weight needs unambiguous bases.");
            w += m;
        }
        w -= (s.Length - 1) * Water;
        if (circular) w -= Water;
        return w;
    }
}
