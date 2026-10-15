using System.Globalization;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>Backbone fragment ion series (Roepstorff–Fohlman–Biemann nomenclature).</summary>
public enum IonType
{
    /// <summary>a = b − CO (N-terminal).</summary>
    A,

    /// <summary>b: N-terminal acylium ion, Σ residues + H⁺ (CID/HCD).</summary>
    B,

    /// <summary>c = b + NH₃ (N-terminal, ETD/ECD).</summary>
    C,

    /// <summary>x = y + CO − 2H (C-terminal).</summary>
    X,

    /// <summary>y: C-terminal, Σ residues + H₂O + H⁺ (CID/HCD).</summary>
    Y,

    /// <summary>z = y − NH₃ (even-electron z, as in Mascot's "z").</summary>
    Z,

    /// <summary>z• = y − NH₃ + H (radical z-dot ion observed in ETD/ECD).</summary>
    ZDot,
}

/// <summary>A theoretical fragment ion.</summary>
/// <param name="Type">Ion series.</param>
/// <param name="Number">Number of residues in the fragment (b3 has 3).</param>
/// <param name="Charge">Charge state.</param>
/// <param name="Mz">Monoisotopic m/z.</param>
public readonly record struct FragmentIon(IonType Type, int Number, int Charge, double Mz)
{
    /// <summary>True for a/b/c ions.</summary>
    public bool IsNTerminal => Type is IonType.A or IonType.B or IonType.C;

    /// <summary>Label such as <c>b3</c> or <c>y7++</c>.</summary>
    public string Label => $"{Type switch { IonType.ZDot => "z•", _ => Type.ToString().ToLowerInvariant() }}{Number.ToString(CultureInfo.InvariantCulture)}{new string('+', Charge > 1 ? Charge : 0)}";

    /// <inheritdoc/>
    public override string ToString() => $"{Label} {Mz.ToString("F5", CultureInfo.InvariantCulture)}";
}

/// <summary>Theoretical fragment ladders for (modified) peptides.</summary>
public static class FragmentCalculator
{
    /// <summary>Neutral mass offset of each series relative to the summed residue masses of the fragment.</summary>
    public static double SeriesOffset(IonType type) => type switch
    {
        IonType.A => -MassConstants.CarbonMonoxide,
        IonType.B => 0,
        IonType.C => MassConstants.Ammonia,
        IonType.X => MassConstants.WaterMonoisotopic + MassConstants.CarbonMonoxide - 2 * MassConstants.Hydrogen,
        IonType.Y => MassConstants.WaterMonoisotopic,
        IonType.Z => MassConstants.WaterMonoisotopic - MassConstants.Ammonia,
        IonType.ZDot => MassConstants.WaterMonoisotopic - MassConstants.Ammonia + MassConstants.Hydrogen,
        _ => throw new ArgumentOutOfRangeException(nameof(type)),
    };

    /// <summary>
    /// The m/z ladder of one series: element <c>i</c> is the ion with <c>i + 1</c> residues, for
    /// <c>i = 0 .. Length − 2</c> (the full-length "fragment" is the precursor and is omitted).
    /// </summary>
    public static double[] Ladder(Peptide peptide, IonType type, int charge = 1)
    {
        ArgumentNullException.ThrowIfNull(peptide);
        var output = new double[Math.Max(0, peptide.Length - 1)];
        Ladder(peptide, type, charge, output);
        return output;
    }

    /// <summary>Writes the ladder of one series into <paramref name="output"/> (length ≥ <c>Length − 1</c>).</summary>
    public static void Ladder(Peptide peptide, IonType type, int charge, Span<double> output)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(charge, 1);
        var res = peptide.ResidueMassesMonoisotopic;
        int n = res.Length;
        if (output.Length < n - 1) throw new ArgumentException("Output span too short.", nameof(output));
        bool nTerminal = type is IonType.A or IonType.B or IonType.C;
        double sum = SeriesOffset(type) + (nTerminal
            ? peptide.NTermModification?.MonoisotopicDelta ?? 0
            : peptide.CTermModification?.MonoisotopicDelta ?? 0);
        for (int i = 0; i < n - 1; i++)
        {
            sum += nTerminal ? res[i] : res[n - 1 - i];
            output[i] = (sum + charge * MassConstants.Proton) / charge;
        }
    }

    /// <summary>All fragments of the requested series and charges 1..<paramref name="maxCharge"/>.</summary>
    public static List<FragmentIon> Fragments(Peptide peptide, IReadOnlyList<IonType>? types = null, int maxCharge = 1)
    {
        types ??= [IonType.B, IonType.Y];
        var list = new List<FragmentIon>(types.Count * maxCharge * peptide.Length);
        var buf = new double[Math.Max(0, peptide.Length - 1)];
        foreach (var t in types)
            for (int z = 1; z <= maxCharge; z++)
            {
                Ladder(peptide, t, z, buf);
                for (int i = 0; i < buf.Length; i++) list.Add(new(t, i + 1, z, buf[i]));
            }
        return list;
    }

    /// <summary>Fragments of an unmodified peptide.</summary>
    public static List<FragmentIon> Fragments(string sequence, IReadOnlyList<IonType>? types = null, int maxCharge = 1)
        => Fragments(new Peptide(sequence), types, maxCharge);
}
