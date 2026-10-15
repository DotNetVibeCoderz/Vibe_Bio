using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments;

/// <summary>
/// The result of <see cref="PairwiseAligner.Align(Sequence, Sequence)"/>: score, gapped strings, coordinates,
/// statistics, CIGAR and a Biopython-style pretty print. <see cref="Target"/> is the reference, <see cref="Query"/> the
/// sequence aligned to it; coordinates are 0-based half-open.
/// </summary>
public sealed class PairwiseAlignment
{
    private readonly SubstitutionMatrix? _matrix;
    private string? _alignedTarget, _alignedQuery;

    internal PairwiseAlignment(Sequence target, Sequence query, AlignmentSummary summary, string operations, AlignmentMode mode, SubstitutionMatrix? matrix)
    {
        Target = target;
        Query = query;
        Score = summary.Score;
        TargetStart = summary.AStart;
        TargetEnd = summary.AEnd;
        QueryStart = summary.BStart;
        QueryEnd = summary.BEnd;
        Operations = operations;
        Mode = mode;
        _matrix = matrix;
        int ti = TargetStart, qi = QueryStart;
        foreach (char op in operations)
        {
            switch (op)
            {
                case 'M':
                    if (char.ToUpperInvariant((char)target[ti]) == char.ToUpperInvariant((char)query[qi])) Identities++; else Mismatches++;
                    if (matrix is not null && matrix.Score(target[ti], query[qi]) > 0) Positives++;
                    ti++; qi++;
                    break;
                case 'I': qi++; Gaps++; break;
                default: ti++; Gaps++; break;
            }
        }
        if (matrix is null) Positives = Identities;
        for (int k = 0; k < operations.Length; k++)
            if (operations[k] != 'M' && (k == 0 || operations[k - 1] != operations[k])) GapOpenings++;
    }

    /// <summary>The reference (first) sequence.</summary>
    public Sequence Target { get; }

    /// <summary>The query (second) sequence.</summary>
    public Sequence Query { get; }

    /// <summary>Alignment score.</summary>
    public int Score { get; }

    /// <summary>Mode the alignment was computed in.</summary>
    public AlignmentMode Mode { get; }

    /// <summary>Edit operations: <c>M</c> aligned pair, <c>I</c> query residue against a gap, <c>D</c> target residue against a gap.</summary>
    public string Operations { get; }

    /// <summary>Start of the aligned region in <see cref="Target"/> (0 for global/semi-global).</summary>
    public int TargetStart { get; }

    /// <summary>End (exclusive) of the aligned region in <see cref="Target"/>.</summary>
    public int TargetEnd { get; }

    /// <summary>Start of the aligned region in <see cref="Query"/>.</summary>
    public int QueryStart { get; }

    /// <summary>End (exclusive) of the aligned region in <see cref="Query"/>.</summary>
    public int QueryEnd { get; }

    /// <summary>Number of alignment columns.</summary>
    public int Length => Operations.Length;

    /// <summary>Columns with identical residues (case-insensitive).</summary>
    public int Identities { get; }

    /// <summary>Columns with different residues.</summary>
    public int Mismatches { get; }

    /// <summary>Columns with a positive substitution score (equals <see cref="Identities"/> without a matrix).</summary>
    public int Positives { get; }

    /// <summary>Gap columns.</summary>
    public int Gaps { get; }

    /// <summary>Number of gaps (runs of gap columns in the same sequence).</summary>
    public int GapOpenings { get; }

    /// <summary>Identities over alignment length (0 for an empty alignment), as reported by EMBOSS needle/water.</summary>
    public double Identity => Length == 0 ? 0 : (double)Identities / Length;

    /// <summary>Target row with <c>-</c> for gaps.</summary>
    public string AlignedTarget => _alignedTarget ??= BuildRow(true);

    /// <summary>Query row with <c>-</c> for gaps.</summary>
    public string AlignedQuery => _alignedQuery ??= BuildRow(false);

    private string BuildRow(bool target)
    {
        var sb = new StringBuilder(Operations.Length);
        int ti = TargetStart, qi = QueryStart;
        foreach (char op in Operations)
        {
            if (op == 'M') { sb.Append(target ? (char)Target[ti] : (char)Query[qi]); ti++; qi++; }
            else if (op == 'I') { sb.Append(target ? '-' : (char)Query[qi]); qi++; }
            else { sb.Append(target ? (char)Target[ti] : '-'); ti++; }
        }
        return sb.ToString();
    }

    /// <summary>
    /// Path vertices in the style of Biopython's <c>alignment.coordinates</c>: row 0 holds target positions and row 1
    /// query positions at the start, at every change between aligned/gapped segments, and at the end.
    /// </summary>
    public int[][] Coordinates
    {
        get
        {
            var t = new List<int> { TargetStart };
            var q = new List<int> { QueryStart };
            int ti = TargetStart, qi = QueryStart;
            for (int k = 0; k < Operations.Length; k++)
            {
                char op = Operations[k];
                if (op != 'I') ti++;
                if (op != 'D') qi++;
                if (k == Operations.Length - 1 || Operations[k + 1] != op) { t.Add(ti); q.Add(qi); }
            }
            return [t.ToArray(), q.ToArray()];
        }
    }

    /// <summary>Aligned (gap-free) blocks as pairs of half-open ranges, like Biopython's <c>alignment.aligned</c>.</summary>
    public IReadOnlyList<(int TargetStart, int TargetEnd, int QueryStart, int QueryEnd)> AlignedBlocks
    {
        get
        {
            var blocks = new List<(int, int, int, int)>();
            int ti = TargetStart, qi = QueryStart, k = 0;
            while (k < Operations.Length)
            {
                char op = Operations[k];
                int run = 0;
                while (k < Operations.Length && Operations[k] == op) { run++; k++; }
                if (op == 'M') { blocks.Add((ti, ti + run, qi, qi + run)); ti += run; qi += run; }
                else if (op == 'I') qi += run;
                else ti += run;
            }
            return blocks;
        }
    }

    /// <summary>
    /// SAM CIGAR of the query against the target: <c>M</c> (or <c>=</c>/<c>X</c> when <paramref name="extended"/>),
    /// <c>I</c> for query insertions, <c>D</c> for target residues missing from the query. With
    /// <paramref name="softClip"/>, unaligned query ends of a local alignment are reported as <c>S</c>.
    /// </summary>
    public string Cigar(bool extended = false, bool softClip = false)
    {
        var sb = new StringBuilder();
        if (softClip && QueryStart > 0) sb.Append(QueryStart).Append('S');
        char prev = '\0';
        int run = 0, ti = TargetStart, qi = QueryStart;
        foreach (char op in Operations)
        {
            char c = op;
            if (op == 'M')
            {
                if (extended) c = char.ToUpperInvariant((char)Target[ti]) == char.ToUpperInvariant((char)Query[qi]) ? '=' : 'X';
                ti++; qi++;
            }
            else if (op == 'I') qi++;
            else ti++;
            if (c != prev && run > 0) { sb.Append(run).Append(prev); run = 0; }
            prev = c; run++;
        }
        if (run > 0) sb.Append(run).Append(prev);
        if (softClip && QueryEnd < Query.Length) sb.Append(Query.Length - QueryEnd).Append('S');
        return sb.ToString();
    }

    /// <summary>
    /// Biopython-style text: blocks of three lines (target, match line, query) with start/end coordinates, where the
    /// match line shows <c>|</c> for identities, <c>.</c> for mismatches and <c>-</c> for gaps.
    /// </summary>
    public string Format(int width = 60)
    {
        if (width < 1) throw new ArgumentOutOfRangeException(nameof(width));
        string tName = Label(Target.Id, "target"), qName = Label(Query.Id, "query");
        string at = AlignedTarget, aq = AlignedQuery;
        var match = new StringBuilder(at.Length);
        for (int k = 0; k < at.Length; k++)
            match.Append(at[k] == '-' || aq[k] == '-' ? '-' : char.ToUpperInvariant(at[k]) == char.ToUpperInvariant(aq[k]) ? '|' : '.');
        var sb = new StringBuilder();
        int tp = TargetStart, qp = QueryStart;
        for (int start = 0; start < at.Length || start == 0; start += width)
        {
            int len = Math.Min(width, at.Length - start);
            string ts = at.Substring(start, len), qs = aq.Substring(start, len), ms = match.ToString(start, len);
            int tEnd = tp + ts.Count(c => c != '-'), qEnd = qp + qs.Count(c => c != '-');
            if (start > 0) sb.AppendLine();
            sb.AppendLine(Line(tName, tp, ts, tEnd));
            sb.AppendLine(Line(string.Empty, start, ms, start + len));
            sb.AppendLine(Line(qName, qp, qs, qEnd));
            tp = tEnd; qp = qEnd;
            if (at.Length == 0) break;
        }
        return sb.ToString();
    }

    private static string Label(string id, string fallback)
    {
        string s = string.IsNullOrWhiteSpace(id) ? fallback : id;
        return s.Length > 9 ? s[..9] : s;
    }

    private static string Line(string name, int start, string segment, int end)
        => $"{name,-10}{start.ToString(CultureInfo.InvariantCulture),9} {segment} {end.ToString(CultureInfo.InvariantCulture)}";

    /// <inheritdoc />
    public override string ToString() => Format();
}
