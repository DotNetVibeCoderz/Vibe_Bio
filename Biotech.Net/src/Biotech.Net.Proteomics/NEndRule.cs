namespace Biotech.Net.Proteomics;

/// <summary>Experimental system an N-end-rule half-life refers to (the three ExPASy ProtParam reports).</summary>
public enum HalfLifeSystem
{
    /// <summary>Mammalian reticulocytes, in vitro.</summary>
    MammalianReticulocytes,

    /// <summary>Yeast, in vivo.</summary>
    Yeast,

    /// <summary>Escherichia coli, in vivo.</summary>
    EscherichiaColi,
}

/// <summary>An N-end-rule half-life estimate.</summary>
/// <param name="Hours">Half-life in hours (NaN when unknown).</param>
/// <param name="IsLowerBound">True for values reported as "&gt;" (e.g. "&gt;20 hours").</param>
/// <param name="Text">ProtParam-style text, e.g. "30 hours" or "&gt;10 hours".</param>
public readonly record struct HalfLife(double Hours, bool IsLowerBound, string Text);

/// <summary>
/// The N-end rule (Bachmair, Finley &amp; Varshavsky 1986; Gonda et al. 1989; Tobias et al. 1991): half-life as a
/// function of the N-terminal residue, with the values reported by ExPASy ProtParam.
/// </summary>
public static class NEndRule
{
    private static readonly Dictionary<char, (string Mammal, string Yeast, string Ecoli)> Table = new()
    {
        ['A'] = ("4.4 hours", ">20 hours", ">10 hours"),
        ['R'] = ("1 hour", "2 min", "2 min"),
        ['N'] = ("1.4 hours", "3 min", ">10 hours"),
        ['D'] = ("1.1 hours", "3 min", ">10 hours"),
        ['C'] = ("1.2 hours", ">20 hours", ">10 hours"),
        ['Q'] = ("0.8 hours", "10 min", ">10 hours"),
        ['E'] = ("1 hour", "30 min", ">10 hours"),
        ['G'] = ("30 hours", ">20 hours", ">10 hours"),
        ['H'] = ("3.5 hours", "10 min", ">10 hours"),
        ['I'] = ("20 hours", "30 min", ">10 hours"),
        ['L'] = ("5.5 hours", "3 min", "2 min"),
        ['K'] = ("1.3 hours", "3 min", "2 min"),
        ['M'] = ("30 hours", ">20 hours", ">10 hours"),
        ['F'] = ("1.1 hours", "3 min", "2 min"),
        ['P'] = (">20 hours", ">20 hours", "?"),
        ['S'] = ("1.9 hours", ">20 hours", ">10 hours"),
        ['T'] = ("7.2 hours", ">20 hours", ">10 hours"),
        ['W'] = ("2.8 hours", "3 min", "2 min"),
        ['Y'] = ("2.8 hours", "10 min", "2 min"),
        ['V'] = ("100 hours", ">20 hours", ">10 hours"),
    };

    /// <summary>Half-life of a protein whose N-terminal residue is <paramref name="nTerminal"/>.</summary>
    public static HalfLife HalfLife(char nTerminal, HalfLifeSystem system = HalfLifeSystem.MammalianReticulocytes)
    {
        if (!Table.TryGetValue(char.ToUpperInvariant(nTerminal), out var row)) return new(double.NaN, false, "unknown");
        string text = system switch
        {
            HalfLifeSystem.MammalianReticulocytes => row.Mammal,
            HalfLifeSystem.Yeast => row.Yeast,
            _ => row.Ecoli,
        };
        return Parse(text);
    }

    private static HalfLife Parse(string text)
    {
        if (text == "?") return new(double.NaN, false, "unknown");
        bool lower = text.StartsWith('>');
        var parts = text.TrimStart('>').Split(' ');
        double v = double.Parse(parts[0], System.Globalization.CultureInfo.InvariantCulture);
        double hours = parts[1].StartsWith("min", StringComparison.Ordinal) ? v / 60.0 : v;
        return new(hours, lower, text);
    }
}
