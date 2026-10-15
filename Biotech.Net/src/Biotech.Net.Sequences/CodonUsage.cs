using System.Globalization;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>
/// Codon counts of a set of coding sequences, with derived measures: frequency per thousand, RSCU (relative
/// synonymous codon usage) and Sharp &amp; Li relative adaptiveness (w). Codons are indexed in kernel order
/// (A=0, C=1, G=2, T=3, index = 16·b1 + 4·b2 + b3), like <see cref="GeneticCode"/>.
/// </summary>
public sealed class CodonUsageTable
{
    private readonly long[] _counts = new long[64];

    /// <summary>Creates an empty table for <paramref name="code"/>.</summary>
    public CodonUsageTable(GeneticCode? code = null) => Code = code ?? GeneticCode.Standard;

    /// <summary>Builds a table from in-frame coding sequences (codons with ambiguous bases are skipped).</summary>
    public static CodonUsageTable FromSequences(IEnumerable<Sequence> codingSequences, GeneticCode? code = null)
    {
        var t = new CodonUsageTable(code);
        foreach (var s in codingSequences) t.Add(s);
        return t;
    }

    /// <summary>The genetic code used to group synonymous codons.</summary>
    public GeneticCode Code { get; }

    /// <summary>Total codons counted.</summary>
    public long Total { get; private set; }

    /// <summary>Adds the codons of an in-frame coding sequence.</summary>
    public CodonUsageTable Add(Sequence codingSequence)
    {
        var s = codingSequence.Span;
        for (int i = 0; i + 3 <= s.Length; i += 3)
        {
            int k = GeneticCode.Index(s.Slice(i, 3));
            if (k < 0) continue;
            _counts[k]++;
            Total++;
        }
        return this;
    }

    /// <summary>Count of a codon (DNA or RNA, any case).</summary>
    public long this[string codon] => _counts[IndexOf(codon)];

    /// <summary>Count at kernel index <paramref name="index"/>.</summary>
    public long CountAt(int index) => _counts[index];

    /// <summary>Occurrences per thousand codons.</summary>
    public double PerThousand(string codon) => Total == 0 ? 0 : 1000.0 * this[codon] / Total;

    /// <summary>
    /// RSCU: observed count divided by the count expected if all synonymous codons were used equally
    /// (1 = no bias, &gt;1 preferred). NaN when the amino acid never occurs.
    /// </summary>
    public double Rscu(string codon)
    {
        int k = IndexOf(codon);
        long sum = 0; int n = 0;
        foreach (int j in Synonyms(k)) { sum += _counts[j]; n++; }
        return sum == 0 ? double.NaN : (double)_counts[k] * n / sum;
    }

    /// <summary>
    /// Relative adaptiveness w = count / count of the most used synonymous codon (Sharp &amp; Li 1987). Codons never
    /// seen get a count of 0.5, as in the original paper and Biopython.
    /// </summary>
    public double RelativeAdaptiveness(string codon) => RelativeAdaptiveness(IndexOf(codon));

    private double RelativeAdaptiveness(int k)
    {
        double max = 0;
        foreach (int j in Synonyms(k)) max = Math.Max(max, Adjusted(j));
        return Adjusted(k) / max;
    }

    private double Adjusted(int k) => _counts[k] == 0 ? 0.5 : _counts[k];

    private IEnumerable<int> Synonyms(int k)
    {
        char aa = Code.AminoAcidAt(k);
        for (int j = 0; j < 64; j++) if (Code.AminoAcidAt(j) == aa) yield return j;
    }

    /// <summary>
    /// The codon adaptation index of a coding sequence: the geometric mean of w over its codons, excluding ATG and
    /// TGG (whose w is always 1) and codons with ambiguous bases. Matches Biopython's
    /// <c>CodonAdaptationIndex.calculate</c> when this table was built from the same reference genes.
    /// </summary>
    /// <remarks>
    /// Stop codons are included. Biopython builds the index over all 64 codons — the stop codons form their own
    /// synonymous family — and <c>calculate</c> skips only ATG and TGG, falling back to skipping a stop codon just
    /// when the supplied index has no entry for it. Dropping stops here would make a trailing TAA invisible and
    /// shift the index (0.664 instead of 0.699 on the Biopython test case).
    /// </remarks>
    public double Cai(Sequence codingSequence)
    {
        var w = new double[64];
        for (int k = 0; k < 64; k++) w[k] = Math.Log(RelativeAdaptiveness(k));
        double sum = 0; int n = 0;
        var s = codingSequence.Span;
        for (int i = 0; i + 3 <= s.Length; i += 3)
        {
            int k = GeneticCode.Index(s.Slice(i, 3));
            if (k < 0 || k == Atg || k == Tgg) continue;
            sum += w[k]; n++;
        }
        return n == 0 ? double.NaN : Math.Exp(sum / n);
    }

    private const int Atg = 0b00_11_10, Tgg = 0b11_10_10;

    /// <summary>
    /// Back-translates <paramref name="protein"/> using the most adaptive codon for each amino acid (simple codon
    /// optimisation for expression in the organism this table describes). Ties pick the first codon in ACGT order.
    /// </summary>
    public Sequence Optimize(Sequence protein)
    {
        var best = new Dictionary<char, int>();
        for (int k = 0; k < 64; k++)
        {
            char aa = Code.AminoAcidAt(k);
            if (!best.TryGetValue(aa, out int b) || Adjusted(k) > Adjusted(b)) best[aa] = k;
        }
        var dst = new byte[protein.Length * 3];
        for (int i = 0; i < protein.Length; i++)
        {
            char aa = char.ToUpperInvariant((char)protein[i]);
            if (!best.TryGetValue(aa, out int k)) throw new ArgumentException($"Amino acid '{aa}' has no codon in {Code.Name}.", nameof(protein));
            Encoding.ASCII.GetBytes(GeneticCode.CodonAt(k), dst.AsSpan(i * 3, 3));
        }
        return Sequence.FromBuffer(protein.Id, dst, Alphabet.Dna, protein.Description, validate: false);
    }

    /// <summary>A tab-separated table: codon, amino acid, count, per-thousand, RSCU, w.</summary>
    public string ToTable()
    {
        var sb = new StringBuilder("codon\taa\tcount\tper_thousand\trscu\tw\n");
        for (int k = 0; k < 64; k++)
        {
            string c = GeneticCode.CodonAt(k);
            sb.Append(CultureInfo.InvariantCulture, $"{c}\t{Code.AminoAcidAt(k)}\t{_counts[k]}\t{PerThousand(c):F2}\t{Rscu(c):F3}\t{RelativeAdaptiveness(k):F3}\n");
        }
        return sb.ToString();
    }

    private static int IndexOf(string codon)
    {
        int k = GeneticCode.Index(codon.AsSpan());
        if (codon.Length != 3 || k < 0) throw new ArgumentException($"'{codon}' is not an unambiguous codon.", nameof(codon));
        return k;
    }
}
