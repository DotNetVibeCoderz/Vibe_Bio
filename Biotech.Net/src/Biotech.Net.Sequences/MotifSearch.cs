using System.Text;
using System.Text.RegularExpressions;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>A motif occurrence. <see cref="Start"/>/<see cref="End"/> are 0-based half-open on the forward strand.</summary>
/// <param name="Matched">The matched residues read on <see cref="Strand"/> (reverse-complemented for reverse hits).</param>
public readonly record struct MotifMatch(int Start, int End, Strand Strand, string Matched)
{
    /// <summary>Length of the match.</summary>
    public int Length => End - Start;

    public override string ToString() => $"{Matched} at {Start + 1}..{End} ({(Strand == Strand.Reverse ? "-" : "+")})";
}

/// <summary>
/// IUPAC nucleotide codes as 4-bit masks (A=1, C=2, G=4, T/U=8). A pattern position matches a sequence base
/// when the base's mask is a subset of the pattern's mask, so <c>N</c> in a pattern matches anything while an
/// <c>N</c> in the sequence only matches a pattern <c>N</c>.
/// </summary>
public static class Iupac
{
    private static readonly byte[] Masks = BuildMasks();

    private static byte[] BuildMasks()
    {
        var m = new byte[256];
        void Set(char c, int v) { m[c] = (byte)v; m[char.ToLowerInvariant(c)] = (byte)v; }
        Set('A', 1); Set('C', 2); Set('G', 4); Set('T', 8); Set('U', 8);
        Set('R', 1 | 4); Set('Y', 2 | 8); Set('S', 2 | 4); Set('W', 1 | 8); Set('K', 4 | 8); Set('M', 1 | 2);
        Set('B', 2 | 4 | 8); Set('D', 1 | 4 | 8); Set('H', 1 | 2 | 8); Set('V', 1 | 2 | 4); Set('N', 15);
        return m;
    }

    /// <summary>The 4-bit mask of an IUPAC symbol (0 for gaps and non-nucleotide symbols).</summary>
    public static byte Mask(byte symbol) => Masks[symbol];

    /// <summary>Whether sequence base <paramref name="baseSymbol"/> is compatible with pattern symbol <paramref name="patternSymbol"/>.</summary>
    public static bool Matches(byte patternSymbol, byte baseSymbol)
    {
        int b = Masks[baseSymbol], p = Masks[patternSymbol];
        return b != 0 && (b & ~p) == 0;
    }

    /// <summary>
    /// A regular expression equivalent to an IUPAC pattern (e.g. <c>GANTC</c> → <c>GA[ACGTU][TU]C</c>). The classes
    /// match <see cref="Matches"/>, where T and U share one mask, so U is accepted wherever T is.
    /// </summary>
    public static string ToRegex(string pattern)
    {
        var sb = new StringBuilder();
        foreach (char c in pattern.ToUpperInvariant())
            sb.Append(c switch
            {
                'A' or 'C' or 'G' => c.ToString(),
                'T' or 'U' => "[TU]",
                'R' => "[AG]", 'Y' => "[CTU]", 'S' => "[CG]", 'W' => "[ATU]", 'K' => "[GTU]", 'M' => "[AC]",
                'B' => "[CGTU]", 'D' => "[AGTU]", 'H' => "[ACTU]", 'V' => "[ACG]", 'N' => "[ACGTU]",
                _ => throw new ArgumentException($"'{c}' is not an IUPAC nucleotide code.", nameof(pattern)),
            });
        return sb.ToString();
    }

    /// <summary>The reverse complement of an IUPAC pattern.</summary>
    public static string ReverseComplement(string pattern)
    {
        var src = Encoding.ASCII.GetBytes(pattern);
        var dst = new byte[src.Length];
        Native.BioKernels.ReverseComplement(src, dst);
        return Encoding.ASCII.GetString(dst);
    }
}

/// <summary>
/// Finds every occurrence of a motif on one or both strands: IUPAC degenerate nucleotide patterns (with optional
/// mismatches), .NET regular expressions and PROSITE patterns for proteins. Overlapping matches are all reported.
/// </summary>
public static class MotifSearch
{
    /// <summary>
    /// All occurrences of the IUPAC <paramref name="pattern"/> (e.g. <c>GAATTC</c>, <c>CANNTG</c>). The reverse strand
    /// is searched by matching the reverse-complemented pattern on the forward strand (no copy of the sequence);
    /// for a palindromic pattern each site is reported once, on the forward strand.
    /// </summary>
    /// <param name="maxMismatches">Number of positions allowed not to match.</param>
    public static IReadOnlyList<MotifMatch> FindIupac(Sequence sequence, string pattern, bool bothStrands = true, int maxMismatches = 0)
    {
        ArgumentException.ThrowIfNullOrEmpty(pattern);
        var fwd = Encoding.ASCII.GetBytes(pattern.ToUpperInvariant());
        foreach (byte b in fwd)
            if (Iupac.Mask(b) == 0) throw new ArgumentException($"'{(char)b}' is not an IUPAC nucleotide code.", nameof(pattern));
        var result = new List<MotifMatch>();
        var s = sequence.Span;
        Scan(s, fwd, maxMismatches, Strand.Forward, result);
        if (bothStrands && sequence.Alphabet.IsNucleotide)
        {
            var rev = Encoding.ASCII.GetBytes(Iupac.ReverseComplement(pattern.ToUpperInvariant()));
            if (!rev.AsSpan().SequenceEqual(fwd)) Scan(s, rev, maxMismatches, Strand.Reverse, result);
        }
        result.Sort((a, b) => a.Start != b.Start ? a.Start.CompareTo(b.Start) : a.Strand.CompareTo(b.Strand));
        return result;
    }

    private static void Scan(ReadOnlySpan<byte> s, byte[] pattern, int maxMismatches, Strand strand, List<MotifMatch> result)
    {
        int m = pattern.Length;
        Span<byte> masks = m <= 256 ? stackalloc byte[m] : new byte[m];
        for (int j = 0; j < m; j++) masks[j] = Iupac.Mask(pattern[j]);
        for (int i = 0; i + m <= s.Length; i++)
        {
            int mm = 0;
            for (int j = 0; j < m; j++)
            {
                int b = Iupac.Mask(s[i + j]);
                if (b == 0 || (b & ~masks[j]) != 0) { if (++mm > maxMismatches) break; }
            }
            if (mm > maxMismatches) continue;
            string text = Encoding.ASCII.GetString(s.Slice(i, m));
            if (strand == Strand.Reverse) text = Iupac.ReverseComplement(text);
            result.Add(new MotifMatch(i, i + m, strand, text));
        }
    }

    /// <summary>
    /// All (overlapping) matches of a regular expression. Matching is case-insensitive. Reverse-strand hits come
    /// from running the expression over the reverse complement and mapping the coordinates back.
    /// </summary>
    public static IReadOnlyList<MotifMatch> FindRegex(Sequence sequence, string regex, bool bothStrands = false)
    {
        // A lookahead around the expression makes overlapping matches visible.
        var re = new Regex($"(?=({regex}))", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
        var result = new List<MotifMatch>();
        string fwd = sequence.ToString();
        foreach (Match m in re.Matches(fwd))
        {
            var g = m.Groups[1];
            if (g.Length == 0) continue;
            result.Add(new MotifMatch(g.Index, g.Index + g.Length, Strand.Forward, g.Value));
        }
        if (bothStrands && sequence.Alphabet.IsNucleotide)
        {
            string rc = sequence.ReverseComplement().ToString();
            int n = rc.Length;
            foreach (Match m in re.Matches(rc))
            {
                var g = m.Groups[1];
                if (g.Length == 0) continue;
                result.Add(new MotifMatch(n - g.Index - g.Length, n - g.Index, Strand.Reverse, g.Value));
            }
        }
        result.Sort((a, b) => a.Start != b.Start ? a.Start.CompareTo(b.Start) : a.Strand.CompareTo(b.Strand));
        return result;
    }

    /// <summary>
    /// Converts a PROSITE pattern to a .NET regular expression: <c>-</c> separators, <c>x</c> any residue,
    /// <c>[...]</c> allowed, <c>{...}</c> forbidden, <c>(n)</c>/<c>(n,m)</c> repeats, <c>&lt;</c>/<c>&gt;</c> anchors.
    /// Example: <c>C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H</c> (zinc finger C2H2, PS00028).
    /// </summary>
    public static string PrositeToRegex(string prosite)
    {
        var p = prosite.Trim().TrimEnd('.');
        var sb = new StringBuilder();
        foreach (var raw in p.Split('-', StringSplitOptions.RemoveEmptyEntries))
        {
            var el = raw.Trim();
            if (el.StartsWith('<')) { sb.Append('^'); el = el[1..]; }
            bool anchorEnd = el.EndsWith('>');
            if (anchorEnd) el = el[..^1];
            string repeat = "";
            int paren = el.IndexOf('(');
            if (paren >= 0)
            {
                repeat = "{" + el[(paren + 1)..el.IndexOf(')')] + "}";
                el = el[..paren];
            }
            string core = el switch
            {
                "x" or "X" => ".",
                _ when el.StartsWith('[') => el.Contains('>') ? $"(?:[{el[1..^1].Replace(">", "")}]|$)" : el,
                _ when el.StartsWith('{') => $"[^{el[1..^1]}]",
                _ => el,
            };
            sb.Append(core).Append(repeat);
            if (anchorEnd) sb.Append('$');
        }
        return sb.ToString();
    }

    /// <summary>All matches of a PROSITE pattern in a protein sequence.</summary>
    public static IReadOnlyList<MotifMatch> FindProsite(Sequence protein, string prosite)
        => FindRegex(protein, PrositeToRegex(prosite), bothStrands: false);

    /// <summary>Number of occurrences of an exact (case-insensitive) word, overlapping, on the forward strand.</summary>
    public static int Count(Sequence sequence, string word)
    {
        int count = 0, i = 0;
        while ((i = sequence.IndexOf(word, i)) >= 0) { count++; i++; }
        return count;
    }
}
