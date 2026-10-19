using System.Globalization;
using System.Reflection;

namespace Biotech.Net.Chemistry;

/// <summary>One experimental torsion preference: a SMARTS with mapped atoms 1–4 and a six-term cosine series.</summary>
internal sealed class TorsionPreference
{
    public required SmartsPattern Pattern { get; init; }
    public required int[] MappedAtoms { get; init; }   // query-atom indices of map numbers 1..4
    public required double[] V { get; init; }          // V1..V6 (kcal/mol)
    public required int[] Signs { get; init; }         // s1..s6
    public required string Smarts { get; init; }
}

/// <summary>
/// The ETKDG experimental torsion library (RDKit's version-2 patterns plus the ETKDGv3 macrocycle patterns), parsed
/// once from embedded data files. Patterns our SMARTS engine cannot parse are skipped and counted.
/// </summary>
internal static class EtkdgTorsionLibrary
{
    private static readonly Lazy<(IReadOnlyList<TorsionPreference> Regular, IReadOnlyList<TorsionPreference> Macrocycle, int Skipped)> Data = new(Load);

    public static IReadOnlyList<TorsionPreference> Regular => Data.Value.Regular;
    public static IReadOnlyList<TorsionPreference> Macrocycle => Data.Value.Macrocycle;
    public static int Skipped => Data.Value.Skipped;

    private static (IReadOnlyList<TorsionPreference>, IReadOnlyList<TorsionPreference>, int) Load()
    {
        int skipped = 0;
        var regular = Parse("etkdg_torsions_v2.txt", ref skipped);
        var macro = Parse("etkdg_torsions_macrocycles.txt", ref skipped);
        return (regular, macro, skipped);
    }

    private static List<TorsionPreference> Parse(string resource, ref int skipped)
    {
        var asm = typeof(EtkdgTorsionLibrary).Assembly;
        string name = asm.GetManifestResourceNames().First(n => n.EndsWith(resource, StringComparison.Ordinal));
        using var stream = asm.GetManifestResourceStream(name)!;
        using var reader = new StreamReader(stream);
        var list = new List<TorsionPreference>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var tok = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            if (tok.Length != 13) { skipped++; continue; }
            var pattern = SmartsPattern.TryParse(tok[0], mergeHydrogens: false);
            if (pattern is null) { skipped++; continue; }
            var mapped = new int[4];
            Array.Fill(mapped, -1);
            for (int q = 0; q < pattern.AtomCount; q++)
            {
                int m = pattern.AtomMap(q);
                if (m is >= 1 and <= 4) mapped[m - 1] = q;
            }
            if (mapped.Any(x => x < 0)) { skipped++; continue; }
            var signs = new int[6];
            var v = new double[6];
            for (int k = 0; k < 6; k++)
            {
                signs[k] = int.Parse(tok[1 + 2 * k], CultureInfo.InvariantCulture);
                v[k] = double.Parse(tok[2 + 2 * k], CultureInfo.InvariantCulture);
            }
            list.Add(new TorsionPreference { Pattern = pattern, MappedAtoms = mapped, V = v, Signs = signs, Smarts = tok[0] });
        }
        return list;
    }
}
