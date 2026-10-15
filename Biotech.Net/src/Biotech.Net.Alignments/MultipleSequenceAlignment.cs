using System.Text;
using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Core;

namespace Biotech.Net.Alignments;

/// <summary>
/// A multiple sequence alignment: equal-length gapped rows (gap symbol <c>-</c>) with identifiers. Provides consensus,
/// column conservation, pairwise identity and FASTA / Clustal output.
/// </summary>
public sealed class MultipleSequenceAlignment
{
    private readonly string[] _ids;
    private readonly string[] _rows;

    /// <summary>Creates an alignment from identifiers and gapped rows of equal length.</summary>
    public MultipleSequenceAlignment(IEnumerable<string> ids, IEnumerable<string> rows, PhylogeneticTree? guideTree = null)
    {
        ArgumentNullException.ThrowIfNull(ids);
        ArgumentNullException.ThrowIfNull(rows);
        _ids = ids.ToArray();
        _rows = rows.ToArray();
        if (_ids.Length != _rows.Length) throw new ArgumentException("One identifier per row is required.");
        if (_rows.Length > 0 && _rows.Any(r => r.Length != _rows[0].Length)) throw new ArgumentException("All rows must have the same length.");
        GuideTree = guideTree;
    }

    /// <summary>Row identifiers.</summary>
    public IReadOnlyList<string> Ids => _ids;

    /// <summary>Gapped rows.</summary>
    public IReadOnlyList<string> Rows => _rows;

    /// <summary>Number of rows.</summary>
    public int Count => _rows.Length;

    /// <summary>Number of columns.</summary>
    public int Length => _rows.Length == 0 ? 0 : _rows[0].Length;

    /// <summary>The gapped row at <paramref name="index"/>.</summary>
    public string this[int index] => _rows[index];

    /// <summary>The guide tree used to build the alignment (progressive alignment only).</summary>
    public PhylogeneticTree? GuideTree { get; }

    /// <summary>True when the residues look like nucleotides (see <see cref="Alphabet.Detect"/>).</summary>
    public bool IsNucleotide => Alphabet.Detect(Encoding.ASCII.GetBytes(string.Concat(_rows).Replace("-", ""))).IsNucleotide;

    /// <summary>Column <paramref name="index"/> top to bottom.</summary>
    public string Column(int index)
    {
        var c = new char[Count];
        for (int r = 0; r < Count; r++) c[r] = _rows[r][index];
        return new string(c);
    }

    /// <summary>The rows without gaps, as sequences.</summary>
    public IReadOnlyList<Sequence> ToSequences()
        => _rows.Select((r, i) => Sequence.Parse(r.Replace("-", "").Replace(".", ""), _ids[i])).ToList();

    private static bool IsGap(char c) => c is '-' or '.';

    /// <summary>
    /// Majority consensus: for each column the most frequent residue (case-insensitive, gaps ignored) if its share of
    /// <em>all</em> rows is at least <paramref name="threshold"/>, otherwise <paramref name="ambiguous"/>; columns that are
    /// at least half gaps give <c>-</c>. Ties go to the alphabetically first residue.
    /// </summary>
    public string Consensus(double threshold = 0.5, char ambiguous = 'X')
    {
        var sb = new StringBuilder(Length);
        for (int col = 0; col < Length; col++)
        {
            var (residue, count, gaps) = MostCommon(col);
            if (gaps * 2 >= Count) sb.Append('-');
            else if (count > 0 && (double)count / Count >= threshold) sb.Append(residue);
            else sb.Append(ambiguous);
        }
        return sb.ToString();
    }

    private (char Residue, int Count, int Gaps) MostCommon(int col)
    {
        Span<int> counts = stackalloc int[128];
        counts.Clear();
        int gaps = 0;
        for (int r = 0; r < Count; r++)
        {
            char c = char.ToUpperInvariant(_rows[r][col]);
            if (IsGap(c)) gaps++;
            else if (c < 128) counts[c]++;
        }
        char best = '\0';
        int bestCount = 0;
        for (int c = 0; c < 128; c++) if (counts[c] > bestCount) { bestCount = counts[c]; best = (char)c; }
        return (best, bestCount, gaps);
    }

    /// <summary>For each column, the fraction of rows carrying its most frequent residue (gaps count against it).</summary>
    public double[] ColumnConservation()
    {
        var result = new double[Length];
        for (int col = 0; col < Length; col++) result[col] = Count == 0 ? 0 : (double)MostCommon(col).Count / Count;
        return result;
    }

    /// <summary>Shannon entropy (bits) of the residue distribution in each column, ignoring gaps (0 = fully conserved).</summary>
    public double[] ColumnEntropy()
    {
        var result = new double[Length];
        for (int col = 0; col < Length; col++)
        {
            var groups = Column(col).Where(c => !IsGap(c)).GroupBy(char.ToUpperInvariant).Select(g => g.Count()).ToList();
            double total = groups.Sum();
            result[col] = total == 0 ? 0 : -groups.Sum(g => g / total * Math.Log2(g / total));
        }
        return result;
    }

    /// <summary>Identities over columns where neither row has a gap (0 when there are none).</summary>
    public double PercentIdentity(int i, int j)
    {
        int same = 0, compared = 0;
        string a = _rows[i], b = _rows[j];
        for (int k = 0; k < a.Length; k++)
        {
            if (IsGap(a[k]) || IsGap(b[k])) continue;
            compared++;
            if (char.ToUpperInvariant(a[k]) == char.ToUpperInvariant(b[k])) same++;
        }
        return compared == 0 ? 0 : 100.0 * same / compared;
    }

    /// <summary>Sum-of-pairs score with a substitution matrix and linear gap score; gap–gap pairs score 0.</summary>
    public long SumOfPairsScore(SubstitutionMatrix matrix, int gapScore = -1)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        long total = 0;
        for (int i = 0; i < Count; i++)
            for (int j = i + 1; j < Count; j++)
                for (int k = 0; k < Length; k++)
                {
                    bool ga = IsGap(_rows[i][k]), gb = IsGap(_rows[j][k]);
                    if (ga && gb) continue;
                    total += ga || gb ? gapScore : matrix.Score(_rows[i][k], _rows[j][k]);
                }
        return total;
    }

    // ---- Output -----------------------------------------------------------------------------------

    /// <summary>Aligned FASTA (gapped rows), wrapped at <paramref name="lineWidth"/> (0 = no wrapping).</summary>
    public string ToFasta(int lineWidth = 60)
    {
        var sb = new StringBuilder();
        for (int r = 0; r < Count; r++)
        {
            sb.Append('>').AppendLine(_ids[r]);
            string row = _rows[r];
            if (lineWidth <= 0 || row.Length == 0) { sb.AppendLine(row); continue; }
            for (int k = 0; k < row.Length; k += lineWidth) sb.AppendLine(row.Substring(k, Math.Min(lineWidth, row.Length - k)));
        }
        return sb.ToString();
    }

    /// <summary>Parses aligned FASTA (all records must have the same length).</summary>
    public static MultipleSequenceAlignment FromFasta(string fasta)
    {
        ArgumentNullException.ThrowIfNull(fasta);
        var ids = new List<string>();
        var rows = new List<StringBuilder>();
        foreach (string raw in fasta.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (line.StartsWith('>'))
            {
                string header = line[1..].Trim();
                int ws = header.IndexOfAny([' ', '\t']);
                ids.Add(ws < 0 ? header : header[..ws]);
                rows.Add(new StringBuilder());
            }
            else if (rows.Count > 0) rows[^1].Append(line.Trim());
        }
        return new MultipleSequenceAlignment(ids, rows.Select(r => r.ToString()));
    }

    private static readonly string[] StrongGroups = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"];
    private static readonly string[] WeakGroups = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"];

    /// <summary>
    /// Clustal conservation line: <c>*</c> fully conserved, <c>:</c> one of the strong (Gonnet &gt; 0.5) groups, <c>.</c>
    /// one of the weak groups (proteins only), space otherwise or when the column has a gap.
    /// </summary>
    public string ConservationLine()
    {
        bool protein = !IsNucleotide;
        var sb = new StringBuilder(Length);
        for (int col = 0; col < Length; col++)
        {
            string c = Column(col).ToUpperInvariant();
            if (c.Length == 0 || c.Any(IsGap)) { sb.Append(' '); continue; }
            if (c.All(x => x == c[0])) sb.Append('*');
            else if (protein && StrongGroups.Any(g => c.All(g.Contains))) sb.Append(':');
            else if (protein && WeakGroups.Any(g => c.All(g.Contains))) sb.Append('.');
            else sb.Append(' ');
        }
        return sb.ToString();
    }

    /// <summary>Clustal (ALN) format with 60-column blocks, cumulative residue counts and the conservation line.</summary>
    public string ToClustal(int blockWidth = 60)
    {
        if (blockWidth < 1) throw new ArgumentOutOfRangeException(nameof(blockWidth));
        var sb = new StringBuilder();
        sb.AppendLine("CLUSTAL W (1.83) multiple sequence alignment - Biotech.Net").AppendLine().AppendLine();
        int nameWidth = Math.Max(10, (_ids.Length == 0 ? 0 : _ids.Max(i => i.Length)) + 6);
        string cons = ConservationLine();
        var counts = new int[Count];
        for (int start = 0; start < Length; start += blockWidth)
        {
            int len = Math.Min(blockWidth, Length - start);
            for (int r = 0; r < Count; r++)
            {
                string seg = _rows[r].Substring(start, len);
                counts[r] += seg.Count(c => !IsGap(c));
                sb.Append(_ids[r].PadRight(nameWidth)).Append(seg).Append(' ').Append(counts[r]).AppendLine();
            }
            sb.Append(new string(' ', nameWidth)).AppendLine(cons.Substring(start, len)).AppendLine();
        }
        return sb.ToString();
    }

    /// <inheritdoc />
    public override string ToString() => ToClustal();
}
