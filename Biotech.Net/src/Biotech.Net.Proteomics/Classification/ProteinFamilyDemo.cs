using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Proteomics.Classification;

/// <summary>
/// A small, reproducible <b>synthetic</b> protein-family dataset for demos, tutorials and tests. Sequences are generated
/// from a seeded RNG with family-specific architectures; they are not real proteins and carry no biological claims.
/// </summary>
/// <remarks>
/// Families: <c>Globin-like</c> (helix-former-rich, His-containing), <c>Collagen-like</c> (Gly-Pro-X repeats),
/// <c>ZincFinger-like</c> (C2H2 repeats C-x2-C-x12-H-x3-H), <c>Transmembrane-like</c> (19–23-residue hydrophobic segments
/// separated by charged loops). Every residue is additionally mutated to a random residue with probability
/// <c>noise</c>, so the task is not trivial.
/// </remarks>
public static class ProteinFamilyDemo
{
    /// <summary>The family names.</summary>
    public static IReadOnlyList<string> Families { get; } = ["Collagen-like", "Globin-like", "Transmembrane-like", "ZincFinger-like"];

    /// <summary>Generates <paramref name="perFamily"/> proteins per family.</summary>
    public static List<(Protein Protein, string Family)> Generate(int perFamily = 20, int seed = 42, double noise = 0.15)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(perFamily, 1);
        if (noise is < 0 or > 1) throw new ArgumentOutOfRangeException(nameof(noise));
        var rng = new Random(seed);
        var list = new List<(Protein, string)>(perFamily * Families.Count);
        foreach (var family in Families)
            for (int i = 0; i < perFamily; i++)
            {
                var sb = new StringBuilder("M");
                int length = rng.Next(90, 180);
                switch (family)
                {
                    case "Globin-like":
                        while (sb.Length < length) sb.Append(Pick(rng, "AAALLLLEEEKKKHVVQMDNG"));
                        break;
                    case "Collagen-like":
                        while (sb.Length < length) sb.Append('G').Append(rng.NextDouble() < 0.6 ? 'P' : Pick(rng, "AEKRSQ")).Append(Pick(rng, "PPPAEKRSDQ"));
                        break;
                    case "ZincFinger-like":
                        while (sb.Length < length)
                        {
                            sb.Append(Random(rng, 3, "STKRGQEP")).Append('C').Append(Random(rng, 2, "PEKSA")).Append('C')
                              .Append(Random(rng, 12, "GKSFRQSSNLTRHLE")).Append('H').Append(Random(rng, 3, "TRGEK")).Append('H');
                        }
                        break;
                    default:
                        while (sb.Length < length)
                        {
                            sb.Append(Random(rng, rng.Next(19, 24), "LLLIIIVVVFFAAMGW"));
                            sb.Append(Random(rng, rng.Next(5, 12), "KKRRDDEENSTGQP"));
                        }
                        break;
                }
                var chars = sb.ToString(0, Math.Min(sb.Length, length)).ToCharArray();
                for (int j = 1; j < chars.Length; j++)
                    if (rng.NextDouble() < noise) chars[j] = AminoAcids.Standard[rng.Next(20)];
                string id = $"{family.Split('-')[0].ToUpperInvariant()}_{i + 1:D3}";
                list.Add((Protein.FromString(id, new string(chars), $"Synthetic {family} protein {i + 1}"), family));
            }
        return list;
    }

    private static char Pick(Random rng, string alphabet) => alphabet[rng.Next(alphabet.Length)];

    private static string Random(Random rng, int n, string alphabet)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(c);
    }
}
