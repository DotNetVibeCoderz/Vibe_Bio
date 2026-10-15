using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.AI.Preprocessing;

/// <summary>
/// A window of a longer sequence, produced by <see cref="SequencePreprocessor.Chunk"/>.
/// </summary>
/// <param name="Start">0-based start of the window in the original sequence.</param>
/// <param name="Text">The residues in the window.</param>
public readonly record struct SequenceChunk(int Start, string Text)
{
    /// <summary>0-based exclusive end of the window.</summary>
    public int End => Start + Text.Length;

    /// <summary>Number of residues in the window.</summary>
    public int Length => Text.Length;
}

/// <summary>
/// Turns raw nucleotide and protein text into the exact string forms biological language models were trained on:
/// cleaned and validated sequences, overlapping k-mers (DNABERT), space-separated residues (ProtBERT),
/// raw upper-case residues (ESM) and overlapping windows for sequences longer than a model's context.
/// </summary>
/// <remarks>
/// Cleaning accepts FASTA input (header lines starting with <c>&gt;</c> or <c>;</c> are dropped), and ignores
/// whitespace and digits, so a GenBank <c>ORIGIN</c> block or a numbered alignment row can be pasted as-is.
/// </remarks>
public static class SequencePreprocessor
{
    private const string IupacNucleotides = "ACGTUNRYSWKMBDHV";
    private const string ProteinLetters = "ACDEFGHIKLMNPQRSTVWYXBZUOJ";

    /// <summary>
    /// Cleans a nucleotide sequence: strips FASTA headers, whitespace and digits, upper-cases, converts U to T and
    /// (by default) maps IUPAC ambiguity codes other than N to N.
    /// </summary>
    /// <param name="raw">The raw text.</param>
    /// <param name="ambiguityToN">Map R, Y, S, W, K, M, B, D, H, V to N. When false they are kept.</param>
    /// <returns>An upper-case DNA string over A, C, G, T, N (plus IUPAC codes when <paramref name="ambiguityToN"/> is false).</returns>
    /// <exception cref="FormatException">The text contains a character that is not an IUPAC nucleotide code.</exception>
    public static string CleanDna(string raw, bool ambiguityToN = true)
    {
        ArgumentNullException.ThrowIfNull(raw);
        var sb = new StringBuilder(raw.Length);
        int position = 0;
        foreach (var line in EnumerateBodyLines(raw))
        {
            foreach (char ch in line)
            {
                if (char.IsWhiteSpace(ch) || char.IsAsciiDigit(ch)) continue;
                char c = char.ToUpperInvariant(ch);
                if (IupacNucleotides.IndexOf(c) < 0)
                    throw new FormatException($"Character '{ch}' at residue {position + 1} is not a nucleotide (IUPAC) code.");
                if (c == 'U') c = 'T';
                else if (ambiguityToN && c is not ('A' or 'C' or 'G' or 'T' or 'N')) c = 'N';
                sb.Append(c);
                position++;
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Cleans a protein sequence: strips FASTA headers, whitespace, digits and one trailing stop (<c>*</c>),
    /// upper-cases, and maps <c>J</c> (Leu/Ile) to <c>X</c>.
    /// </summary>
    /// <param name="raw">The raw text.</param>
    /// <param name="removeGaps">Drop alignment gap characters (<c>-</c> and <c>.</c>) instead of rejecting them.</param>
    /// <returns>An upper-case one-letter amino-acid string.</returns>
    /// <exception cref="FormatException">The text contains an internal stop or a character that is not an amino-acid code.</exception>
    public static string CleanProtein(string raw, bool removeGaps = false)
    {
        ArgumentNullException.ThrowIfNull(raw);
        var sb = new StringBuilder(raw.Length);
        int position = 0;
        bool stopSeen = false;
        foreach (var line in EnumerateBodyLines(raw))
        {
            foreach (char ch in line)
            {
                if (char.IsWhiteSpace(ch) || char.IsAsciiDigit(ch)) continue;
                if (removeGaps && ch is '-' or '.') continue;
                if (ch == '*')
                {
                    stopSeen = true;
                    continue;
                }

                if (stopSeen)
                    throw new FormatException($"Internal stop codon '*' before residue {position + 1}.");
                char c = char.ToUpperInvariant(ch);
                if (ProteinLetters.IndexOf(c) < 0)
                    throw new FormatException($"Character '{ch}' at residue {position + 1} is not an amino-acid code.");
                sb.Append(c == 'J' ? 'X' : c);
                position++;
            }
        }

        return sb.ToString();
    }

    /// <summary>Returns whether <paramref name="raw"/> cleans to a non-empty nucleotide sequence.</summary>
    public static bool IsValidDna(string? raw)
    {
        if (string.IsNullOrWhiteSpace(raw)) return false;
        try { return CleanDna(raw).Length > 0; }
        catch (FormatException) { return false; }
    }

    /// <summary>Returns whether <paramref name="raw"/> cleans to a non-empty protein sequence.</summary>
    public static bool IsValidProtein(string? raw)
    {
        if (string.IsNullOrWhiteSpace(raw)) return false;
        try { return CleanProtein(raw).Length > 0; }
        catch (FormatException) { return false; }
    }

    /// <summary>
    /// Splits a cleaned DNA sequence into k-mers. With <paramref name="stride"/> 1 the k-mers overlap (DNABERT);
    /// with <paramref name="stride"/> equal to <paramref name="k"/> they tile the sequence (Nucleotide Transformer).
    /// </summary>
    /// <param name="dna">The sequence; it is cleaned first.</param>
    /// <param name="k">The k-mer length (1–64).</param>
    /// <param name="stride">Step between consecutive k-mers.</param>
    /// <returns>The k-mers in order. Empty when the sequence is shorter than <paramref name="k"/>.</returns>
    public static IReadOnlyList<string> Kmers(string dna, int k, int stride = 1)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(k, 64);
        ArgumentOutOfRangeException.ThrowIfLessThan(stride, 1);
        string s = CleanDna(dna);
        if (s.Length < k) return [];
        var list = new List<string>((s.Length - k) / stride + 1);
        for (int i = 0; i + k <= s.Length; i += stride) list.Add(s.Substring(i, k));
        return list;
    }

    /// <summary>
    /// The DNABERT input form: overlapping k-mers separated by single spaces, for example
    /// <c>ATGCAT TGCATG GCATGC</c> for <c>ATGCATGC</c> with k = 6.
    /// </summary>
    /// <param name="raw">The raw DNA.</param>
    /// <param name="k">The k-mer length the model was trained with (3–6 for DNABERT).</param>
    /// <exception cref="ArgumentException">The cleaned sequence is shorter than <paramref name="k"/>.</exception>
    public static string KmerTokenize(string raw, int k)
    {
        var kmers = Kmers(raw, k);
        if (kmers.Count == 0)
            throw new ArgumentException($"A sequence of at least {k} nucleotides is needed for {k}-mer tokenisation.", nameof(raw));
        return string.Join(' ', kmers);
    }

    /// <summary>The DNABERT input form for a Biotech.Net <see cref="Sequence"/>.</summary>
    public static string KmerTokenize(Sequence sequence, int k)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return KmerTokenize(sequence.ToString(), k);
    }

    /// <summary>Reassembles a DNA sequence from overlapping (stride 1) k-mers.</summary>
    public static string FromOverlappingKmers(IReadOnlyList<string> kmers)
    {
        ArgumentNullException.ThrowIfNull(kmers);
        if (kmers.Count == 0) return string.Empty;
        var sb = new StringBuilder(kmers[0]);
        for (int i = 1; i < kmers.Count; i++) sb.Append(kmers[i][^1]);
        return sb.ToString();
    }

    /// <summary>The input form for nucleotide models that tokenise their own input: the cleaned upper-case sequence.</summary>
    public static string DnaInput(string raw) => CleanDna(raw);

    /// <summary>
    /// The ProtBERT / ProtT5 input form: residues separated by single spaces, with the rare amino acids
    /// U, Z, O and B mapped to X as the Rostlab models expect.
    /// </summary>
    public static string SpaceResidues(string raw)
    {
        string s = CleanProtein(raw);
        var sb = new StringBuilder(s.Length * 2);
        for (int i = 0; i < s.Length; i++)
        {
            if (i > 0) sb.Append(' ');
            char c = s[i];
            sb.Append(c is 'U' or 'Z' or 'O' or 'B' ? 'X' : c);
        }

        return sb.ToString();
    }

    /// <summary>The ESM input form: the cleaned upper-case one-letter sequence (the ESM tokenizer is per character).</summary>
    public static string EsmInput(string raw) => CleanProtein(raw);

    /// <summary>The ESM input form for a Biotech.Net <see cref="Sequence"/> or <see cref="Protein"/> sequence.</summary>
    public static string EsmInput(Sequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CleanProtein(sequence.ToString());
    }

    /// <summary>
    /// Splits a sequence into windows of at most <paramref name="chunkLength"/> residues that overlap by
    /// <paramref name="overlap"/> residues, so a model with a limited context can see all of it. The final window is
    /// aligned to the end of the sequence, so every window except possibly the only one is full length.
    /// </summary>
    /// <param name="sequence">The (already cleaned) sequence.</param>
    /// <param name="chunkLength">Maximum window length in residues.</param>
    /// <param name="overlap">Residues shared by consecutive windows; must be smaller than <paramref name="chunkLength"/>.</param>
    /// <returns>The windows in order. A sequence no longer than <paramref name="chunkLength"/> gives one window.</returns>
    public static IReadOnlyList<SequenceChunk> Chunk(string sequence, int chunkLength, int overlap = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(chunkLength, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(overlap);
        if (overlap >= chunkLength)
            throw new ArgumentOutOfRangeException(nameof(overlap), "The overlap must be smaller than the chunk length.");
        if (sequence.Length <= chunkLength) return [new SequenceChunk(0, sequence)];

        var chunks = new List<SequenceChunk>();
        int step = chunkLength - overlap;
        int start = 0;
        while (true)
        {
            if (start + chunkLength >= sequence.Length)
            {
                int last = sequence.Length - chunkLength;
                chunks.Add(new SequenceChunk(last, sequence.Substring(last, chunkLength)));
                break;
            }

            chunks.Add(new SequenceChunk(start, sequence.Substring(start, chunkLength)));
            start += step;
        }

        return chunks;
    }

    private static IEnumerable<string> EnumerateBodyLines(string raw)
    {
        foreach (var line in raw.Split('\n'))
        {
            var trimmed = line.TrimStart();
            if (trimmed.StartsWith('>') || trimmed.StartsWith(';')) continue;
            yield return line;
        }
    }
}
