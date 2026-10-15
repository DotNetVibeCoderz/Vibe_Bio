using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>
/// A proteolytic enzyme described by its cleavage rule: cut C-terminal to <see cref="CleaveAfter"/> residues and/or
/// N-terminal to <see cref="CleaveBefore"/> residues, unless the neighbouring residue is in <see cref="BlockedBy"/>.
/// </summary>
/// <param name="Name">Display name.</param>
/// <param name="CleaveAfter">Residues at P1 after which the enzyme cuts.</param>
/// <param name="CleaveBefore">Residues at P1' before which the enzyme cuts.</param>
/// <param name="BlockedBy">Residues at P1' that block a C-terminal cut (P1 that block an N-terminal cut), e.g. P for trypsin.</param>
public sealed record Enzyme(string Name, string CleaveAfter, string CleaveBefore = "", string BlockedBy = "")
{
    /// <summary>Trypsin with the Keil proline rule: after K/R, not before P.</summary>
    public static Enzyme Trypsin { get; } = new("Trypsin", "KR", "", "P");

    /// <summary>Trypsin/P: after K/R, also before P.</summary>
    public static Enzyme TrypsinP { get; } = new("Trypsin/P", "KR");

    /// <summary>Lys-C: after K.</summary>
    public static Enzyme LysC { get; } = new("Lys-C", "K");

    /// <summary>Arg-C: after R, not before P.</summary>
    public static Enzyme ArgC { get; } = new("Arg-C", "R", "", "P");

    /// <summary>Glu-C (V8, bicarbonate buffer): after E, not before P.</summary>
    public static Enzyme GluC { get; } = new("Glu-C", "E", "", "P");

    /// <summary>Asp-N: before D.</summary>
    public static Enzyme AspN { get; } = new("Asp-N", "", "D");

    /// <summary>Chymotrypsin, high specificity: after F/Y/W, not before P.</summary>
    public static Enzyme Chymotrypsin { get; } = new("Chymotrypsin", "FYW", "", "P");

    /// <summary>Chymotrypsin, low specificity: after F/Y/W/L/M, not before P.</summary>
    public static Enzyme ChymotrypsinLow { get; } = new("Chymotrypsin (low)", "FYWLM", "", "P");

    /// <summary>Pepsin (pH &gt; 2), simplified search-engine rule: after F/L/W/Y, not before P.</summary>
    public static Enzyme Pepsin { get; } = new("Pepsin", "FLWY", "", "P");

    /// <summary>The built-in enzymes.</summary>
    public static IReadOnlyList<Enzyme> All { get; } = [Trypsin, TrypsinP, LysC, ArgC, GluC, AspN, Chymotrypsin, ChymotrypsinLow, Pepsin];

    /// <summary>Looks up a built-in enzyme by name (case-insensitive; "trypsin/p", "lysc", "lys-c", ...).</summary>
    public static Enzyme Parse(string name)
    {
        static string Norm(string s) => s.Replace("-", "", StringComparison.Ordinal).Replace(" ", "", StringComparison.Ordinal).ToUpperInvariant();
        string key = Norm(name);
        foreach (var e in All) if (Norm(e.Name) == key) return e;
        throw new ArgumentException($"Unknown enzyme '{name}'.", nameof(name));
    }

    /// <summary>True when the enzyme cuts between <c>residues[site − 1]</c> and <c>residues[site]</c>.</summary>
    public bool Cleaves(ReadOnlySpan<byte> residues, int site)
    {
        if (site <= 0 || site >= residues.Length) return false;
        char p1 = char.ToUpperInvariant((char)residues[site - 1]), p1p = char.ToUpperInvariant((char)residues[site]);
        if (CleaveAfter.Contains(p1) && !BlockedBy.Contains(p1p)) return true;
        return CleaveBefore.Contains(p1p) && !BlockedBy.Contains(p1);
    }

    /// <summary>All cleavage positions including 0 and <c>Length</c>, ascending.</summary>
    public int[] CleavageSites(ReadOnlySpan<byte> residues)
    {
        var sites = new List<int>(residues.Length / 8 + 2) { 0 };
        for (int i = 1; i < residues.Length; i++) if (Cleaves(residues, i)) sites.Add(i);
        if (residues.Length > 0) sites.Add(residues.Length);
        return sites.ToArray();
    }

    /// <inheritdoc/>
    public override string ToString() => Name;
}

/// <summary>Digestion settings.</summary>
/// <param name="MissedCleavages">Maximum missed cleavages per peptide.</param>
/// <param name="MinLength">Minimum peptide length.</param>
/// <param name="MaxLength">Maximum peptide length.</param>
/// <param name="ClipInitiatorMethionine">Also emit protein N-terminal peptides without the initiator Met (as most search engines do).</param>
public sealed record DigestionOptions(int MissedCleavages = 2, int MinLength = 7, int MaxLength = 50, bool ClipInitiatorMethionine = true)
{
    /// <summary>Defaults: 2 missed cleavages, 7–50 residues, initiator-Met clipping.</summary>
    public static DigestionOptions Default { get; } = new();
}

/// <summary>A peptide produced by in-silico digestion.</summary>
/// <param name="Sequence">Residues (upper case).</param>
/// <param name="ProteinId">Source protein.</param>
/// <param name="Start">0-based start in the protein.</param>
/// <param name="End">0-based exclusive end.</param>
/// <param name="MissedCleavages">Internal cleavage sites skipped.</param>
public readonly record struct DigestedPeptide(string Sequence, string ProteinId, int Start, int End, int MissedCleavages)
{
    /// <summary>Peptide length.</summary>
    public int Length => End - Start;

    /// <summary>Neutral monoisotopic mass of the unmodified peptide.</summary>
    public double MonoisotopicMass => ResidueMasses.PeptideMass(Sequence);
}

/// <summary>In-silico proteolytic digestion.</summary>
public static class Digestion
{
    /// <summary>Digests ASCII residues.</summary>
    public static List<DigestedPeptide> Digest(ReadOnlySpan<byte> residues, Enzyme enzyme, DigestionOptions? options = null, string proteinId = "protein")
    {
        ArgumentNullException.ThrowIfNull(enzyme);
        options ??= DigestionOptions.Default;
        var upper = AminoAcids.ToUpperBytes(residues);
        var sites = enzyme.CleavageSites(upper);
        var result = new List<DigestedPeptide>();
        string text = Encoding.ASCII.GetString(upper);
        for (int a = 0; a < sites.Length - 1; a++)
            for (int mc = 0; mc <= options.MissedCleavages && a + mc + 1 < sites.Length; mc++)
            {
                int start = sites[a], end = sites[a + mc + 1];
                Add(start, end, mc);
                if (start == 0 && options.ClipInitiatorMethionine && upper.Length > 1 && upper[0] == 'M' && sites.Length > 1)
                    Add(1, end, mc);
            }
        return result;

        void Add(int start, int end, int mc)
        {
            int len = end - start;
            if (len < options.MinLength || len > options.MaxLength) return;
            result.Add(new(text.Substring(start, len), proteinId, start, end, mc));
        }
    }

    /// <summary>Digests a one-letter string.</summary>
    public static List<DigestedPeptide> Digest(string residues, Enzyme enzyme, DigestionOptions? options = null, string proteinId = "protein")
        => Digest(Encoding.ASCII.GetBytes(residues), enzyme, options, proteinId);

    /// <summary>Digests a Core protein.</summary>
    public static List<DigestedPeptide> Digest(Protein protein, Enzyme enzyme, DigestionOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(protein);
        return Digest(protein.Sequence.Span, enzyme, options, protein.Id);
    }

    /// <summary>
    /// Number of fully cleaved peptides (no missed cleavages) of length 6–30 — the "observable peptide" count used as
    /// the iBAQ denominator (Schwanhäusser et al., Nature 473:337, 2011).
    /// </summary>
    public static int ObservablePeptideCount(ReadOnlySpan<byte> residues, Enzyme? enzyme = null, int minLength = 6, int maxLength = 30)
    {
        var sites = (enzyme ?? Enzyme.Trypsin).CleavageSites(AminoAcids.ToUpperBytes(residues));
        int n = 0;
        for (int i = 0; i + 1 < sites.Length; i++)
        {
            int len = sites[i + 1] - sites[i];
            if (len >= minLength && len <= maxLength) n++;
        }
        return n;
    }
}
