using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Proteomics;

/// <summary>
/// The 20 standard amino acids and fast byte → index lookups shared by the analysis, mass-spectrometry
/// and classification code. Indices follow alphabetical one-letter order <c>ACDEFGHIKLMNPQRSTVWY</c>.
/// </summary>
public static class AminoAcids
{
    /// <summary>The 20 standard one-letter codes in alphabetical order (the index order used by feature vectors).</summary>
    public const string Standard = "ACDEFGHIKLMNPQRSTVWY";

    private static readonly sbyte[] IndexTable = BuildIndex();

    private static sbyte[] BuildIndex()
    {
        var t = new sbyte[256];
        Array.Fill(t, (sbyte)-1);
        for (int i = 0; i < Standard.Length; i++)
        {
            t[Standard[i]] = (sbyte)i;
            t[char.ToLowerInvariant(Standard[i])] = (sbyte)i;
        }
        return t;
    }

    /// <summary>Index (0–19) of a standard residue in <see cref="Standard"/>, or −1 for anything else. Case-insensitive.</summary>
    public static int IndexOf(byte residue) => IndexTable[residue];

    /// <summary>Index (0–19) of a standard residue, or −1.</summary>
    public static int IndexOf(char residue) => residue < 256 ? IndexTable[residue] : -1;

    /// <summary>True for the 20 standard residues (either case).</summary>
    public static bool IsStandard(byte residue) => IndexTable[residue] >= 0;

    /// <summary>Upper-cased ASCII copy of a residue string.</summary>
    internal static byte[] ToUpperBytes(ReadOnlySpan<byte> residues)
    {
        var b = new byte[residues.Length];
        for (int i = 0; i < b.Length; i++)
        {
            byte c = residues[i];
            b[i] = c is >= (byte)'a' and <= (byte)'z' ? (byte)(c - 32) : c;
        }
        return b;
    }

    /// <summary>Upper-cased ASCII bytes of a string.</summary>
    internal static byte[] ToUpperBytes(string residues) => ToUpperBytes(Encoding.ASCII.GetBytes(residues));

    /// <summary>Residue bytes of a Core sequence after checking that it is a protein.</summary>
    internal static ReadOnlySpan<byte> RequireProtein(Sequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Type != SequenceType.Protein)
            throw new ArgumentException($"Expected a protein sequence, got {sequence.Type}.", nameof(sequence));
        return sequence.Span;
    }
}
