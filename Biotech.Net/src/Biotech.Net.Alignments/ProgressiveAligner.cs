using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments;

/// <summary>How <see cref="ProgressiveAligner"/> estimates the distances that drive its guide tree.</summary>
public enum GuideDistance
{
    /// <summary>Alignment-free k-mer distance <c>1 − Σ min(cₐ, c_b) / (min(Lₐ, L_b) − k + 1)</c> (Edgar 2004). Fast.</summary>
    Kmer,

    /// <summary><c>1 − identity</c> of a global pairwise alignment of each pair. Slower, more accurate on divergent sets.</summary>
    Alignment,
}

/// <summary>Clustering method for the guide tree.</summary>
public enum GuideTreeMethod
{
    /// <summary>UPGMA (rooted, ultrametric).</summary>
    Upgma,

    /// <summary>Neighbor-Joining.</summary>
    NeighborJoining,
}

/// <summary>
/// Progressive multiple sequence alignment in the Clustal/MUSCLE tradition: pairwise distances (k-mer or alignment) →
/// guide tree (UPGMA or Neighbor-Joining) → profile–profile alignment up the tree, with affine gaps and a
/// sum-of-pairs column score.
/// </summary>
/// <remarks>Rows of the result are in input order; the guide tree is available as <see cref="MultipleSequenceAlignment.GuideTree"/>.</remarks>
public sealed class ProgressiveAligner
{
    /// <summary>Substitution matrix; null picks NUC.4.4 for nucleotides and BLOSUM62 for proteins.</summary>
    public SubstitutionMatrix? SubstitutionMatrix { get; set; }

    /// <summary>Score of the first gap position (≤ 0). Default −10.</summary>
    public double OpenGapScore { get; set; } = -10;

    /// <summary>Score of each further gap position (≤ 0). Default −1.</summary>
    public double ExtendGapScore { get; set; } = -1;

    /// <summary>Guide-tree distance. Default <see cref="GuideDistance.Kmer"/>.</summary>
    public GuideDistance Distance { get; set; } = GuideDistance.Kmer;

    /// <summary>Guide-tree method. Default <see cref="GuideTreeMethod.Upgma"/>.</summary>
    public GuideTreeMethod TreeMethod { get; set; } = GuideTreeMethod.Upgma;

    /// <summary>k for <see cref="GuideDistance.Kmer"/>; null = 3 for proteins, 4 for nucleotides.</summary>
    public int? KmerSize { get; set; }

    /// <summary>Aligns the sequences. Identifiers that repeat are suffixed (<c>_2</c>, ...) in the guide tree only.</summary>
    public MultipleSequenceAlignment Align(IEnumerable<Sequence> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var seqs = sequences.ToList();
        if (seqs.Count == 0) return new MultipleSequenceAlignment([], []);
        if (OpenGapScore > 0 || ExtendGapScore > 0) throw new InvalidOperationException("Gap scores are zero or negative.");
        var matrix = ResolveMatrix(seqs);
        if (seqs.Count == 1) return new MultipleSequenceAlignment([seqs[0].Id], [seqs[0].ToString()]);

        var dm = ComputeDistances(seqs, matrix);
        var tree = TreeMethod == GuideTreeMethod.Upgma ? TreeBuilder.Upgma(dm) : TreeBuilder.NeighborJoining(dm, clampNegative: true);
        var indexOf = dm.Names.Select((n, i) => (n, i)).ToDictionary(x => x.n, x => x.i);

        var scorer = new ProfileScorer(seqs, matrix);
        var profiles = new Dictionary<TreeNode, Profile>();
        foreach (var node in tree.Root.PostOrder())
        {
            if (node.IsLeaf)
            {
                int idx = indexOf[node.Name!];
                profiles[node] = new Profile([idx], [seqs[idx].ToString().ToCharArray()]);
                continue;
            }
            Profile acc = profiles[node.Children[0]];
            for (int c = 1; c < node.Children.Count; c++) acc = AlignProfiles(acc, profiles[node.Children[c]], scorer);
            profiles[node] = acc;
            foreach (var child in node.Children) profiles.Remove(child);
        }
        var final = profiles[tree.Root];
        var rows = new string[seqs.Count];
        for (int k = 0; k < final.Members.Count; k++) rows[final.Members[k]] = new string(final.Rows[k]);
        return new MultipleSequenceAlignment(seqs.Select(s => s.Id), rows, tree);
    }

    /// <summary>The guide-tree distance matrix for <paramref name="sequences"/> under the current settings.</summary>
    public DistanceMatrix ComputeDistances(IReadOnlyList<Sequence> sequences)
        => ComputeDistances(sequences, ResolveMatrix(sequences));

    private SubstitutionMatrix ResolveMatrix(IReadOnlyList<Sequence> seqs)
        => SubstitutionMatrix ?? (seqs.All(s => s.Alphabet.IsNucleotide) ? SubstitutionMatrix.Nuc44 : SubstitutionMatrix.Blosum62);

    private DistanceMatrix ComputeDistances(IReadOnlyList<Sequence> seqs, SubstitutionMatrix matrix)
    {
        var names = new List<string>();
        var seen = new Dictionary<string, int>();
        foreach (var s in seqs)
        {
            int k = seen.TryGetValue(s.Id, out int c) ? c + 1 : 1;
            seen[s.Id] = k;
            names.Add(k == 1 ? s.Id : $"{s.Id}_{k}");
        }
        var dm = new DistanceMatrix(names);
        if (Distance == GuideDistance.Kmer)
        {
            int k = KmerSize ?? (seqs.All(s => s.Alphabet.IsNucleotide) ? 4 : 3);
            if (k < 1) throw new InvalidOperationException("KmerSize must be positive.");
            var counts = seqs.Select(s => KmerCounts(s, k)).ToList();
            for (int i = 0; i < seqs.Count; i++)
                for (int j = 0; j < i; j++)
                {
                    int denom = Math.Min(seqs[i].Length, seqs[j].Length) - k + 1;
                    if (denom <= 0) { dm[i, j] = 1; continue; }
                    int shared = 0;
                    var (small, large) = counts[i].Count <= counts[j].Count ? (counts[i], counts[j]) : (counts[j], counts[i]);
                    foreach (var (kmer, c) in small) if (large.TryGetValue(kmer, out int c2)) shared += Math.Min(c, c2);
                    dm[i, j] = Math.Max(0, 1 - (double)shared / denom);
                }
        }
        else
        {
            var aligner = new PairwiseAligner
            {
                SubstitutionMatrix = matrix,
                OpenGapScore = (int)Math.Round(OpenGapScore),
                ExtendGapScore = (int)Math.Round(ExtendGapScore),
            };
            for (int i = 0; i < seqs.Count; i++)
                for (int j = 0; j < i; j++)
                {
                    var aln = aligner.Align(seqs[i], seqs[j]);
                    int pairs = aln.Identities + aln.Mismatches;
                    dm[i, j] = pairs == 0 ? 1 : 1 - (double)aln.Identities / pairs;
                }
        }
        return dm;
    }

    private static Dictionary<string, int> KmerCounts(Sequence s, int k)
    {
        var d = new Dictionary<string, int>();
        string text = s.ToString().ToUpperInvariant();
        for (int i = 0; i + k <= text.Length; i++)
        {
            string kmer = text.Substring(i, k);
            d[kmer] = d.TryGetValue(kmer, out int c) ? c + 1 : 1;
        }
        return d;
    }

    // ---- Profile alignment ----------------------------------------------------------------------

    private sealed record Profile(List<int> Members, List<char[]> Rows)
    {
        public int Length => Rows[0].Length;
    }

    /// <summary>Maps residues to a compact alphabet and holds the score table over it.</summary>
    private sealed class ProfileScorer
    {
        public readonly int[] Code = new int[128];
        public readonly double[,] Score;
        public readonly int Size;

        public ProfileScorer(IEnumerable<Sequence> seqs, SubstitutionMatrix matrix)
        {
            Array.Fill(Code, -1);
            var symbols = new List<char>();
            foreach (var s in seqs)
                foreach (byte b in s.Span)
                {
                    char c = char.ToUpperInvariant((char)(b & 0x7F));
                    if (c is '-' or '.' || Code[c] >= 0) continue;
                    Code[c] = symbols.Count;
                    Code[char.ToLowerInvariant(c)] = symbols.Count;
                    symbols.Add(c);
                }
            Size = symbols.Count;
            Score = new double[Size, Size];
            for (int a = 0; a < Size; a++)
                for (int b = 0; b < Size; b++) Score[a, b] = matrix.Score(symbols[a], symbols[b]);
        }

        /// <summary>Per column: residue frequencies (count / rows).</summary>
        public double[][] Frequencies(Profile p)
        {
            var f = new double[p.Length][];
            double inv = 1.0 / p.Rows.Count;
            for (int col = 0; col < p.Length; col++)
            {
                var v = new double[Size];
                foreach (var row in p.Rows)
                {
                    char c = row[col];
                    int code = c < 128 ? Code[c] : -1;
                    if (code >= 0) v[code] += inv;
                }
                f[col] = v;
            }
            return f;
        }
    }

    private Profile AlignProfiles(Profile p, Profile q, ProfileScorer scorer)
    {
        int n = p.Length, m = q.Length, A = scorer.Size;
        double o = -OpenGapScore, e = -ExtendGapScore;
        var fp = scorer.Frequencies(p);
        var fq = scorer.Frequencies(q);
        // w[i][b] = Σa fP_i(a)·S(a,b), so a column pair scores Σb w[i][b]·fQ_j(b).
        var w = new double[n][];
        for (int i = 0; i < n; i++)
        {
            var v = new double[A];
            for (int a = 0; a < A; a++)
            {
                double fa = fp[i][a];
                if (fa == 0) continue;
                for (int b = 0; b < A; b++) v[b] += fa * scorer.Score[a, b];
            }
            w[i] = v;
        }
        var sparseQ = new (int Code, double F)[m][];
        for (int j = 0; j < m; j++)
            sparseQ[j] = Enumerable.Range(0, A).Where(b => fq[j][b] != 0).Select(b => (b, fq[j][b])).ToArray();

        const double neg = double.NegativeInfinity;
        var tb = new byte[(n + 1) * (m + 1)];
        var hprev = new double[m + 1];
        var hcur = new double[m + 1];
        var f = new double[m + 1];
        hprev[0] = 0;
        for (int j = 1; j <= m; j++) { hprev[j] = -(o + (j - 1) * e); tb[j] = 1 | 4; f[j] = neg; }
        tb[1] = 1;
        for (int i = 1; i <= n; i++)
        {
            double left = -(o + (i - 1) * e);
            hcur[0] = left;
            tb[i * (m + 1)] = (byte)(i == 1 ? 2 : 2 | 8);
            double eRun = neg, hLeft = left, diag = hprev[0];
            var wi = w[i - 1];
            for (int j = 1; j <= m; j++)
            {
                double s = 0;
                foreach (var (code, fr) in sparseQ[j - 1]) s += wi[code] * fr;
                double eOpen = hLeft - o, eExt = eRun - e;
                double ev; byte eb;
                if (eExt > eOpen) { ev = eExt; eb = 4; } else { ev = eOpen; eb = 0; }
                double up = hprev[j];
                double fOpen = up - o, fExt = f[j] - e;
                double fv; byte fb;
                if (fExt > fOpen) { fv = fExt; fb = 8; } else { fv = fOpen; fb = 0; }
                double d = diag + s;
                double h; byte src;
                if (d >= ev && d >= fv) { h = d; src = 0; }
                else if (ev >= fv) { h = ev; src = 1; }
                else { h = fv; src = 2; }
                tb[i * (m + 1) + j] = (byte)(src | eb | fb);
                diag = up;
                hcur[j] = h;
                f[j] = fv;
                eRun = ev;
                hLeft = h;
            }
            (hprev, hcur) = (hcur, hprev);
        }

        // Traceback from (n, m).
        var ops = new List<char>(n + m);
        int ti = n, tj = m, state = 0;
        while (ti > 0 || tj > 0)
        {
            if (ti == 0) { ops.Add('I'); tj--; continue; }
            if (tj == 0) { ops.Add('D'); ti--; continue; }
            byte t = tb[ti * (m + 1) + tj];
            if (state == 0)
            {
                int src = t & 3;
                if (src == 0) { ops.Add('M'); ti--; tj--; }
                else state = src;
            }
            else if (state == 1) { ops.Add('I'); tj--; state = (t & 4) != 0 ? 1 : 0; }
            else { ops.Add('D'); ti--; state = (t & 8) != 0 ? 2 : 0; }
        }
        ops.Reverse();

        var rows = new List<char[]>();
        foreach (var _ in p.Rows) rows.Add(new char[ops.Count]);
        foreach (var _ in q.Rows) rows.Add(new char[ops.Count]);
        int pi = 0, qj = 0;
        for (int k = 0; k < ops.Count; k++)
        {
            char op = ops[k];
            for (int r = 0; r < p.Rows.Count; r++) rows[r][k] = op == 'I' ? '-' : p.Rows[r][pi];
            for (int r = 0; r < q.Rows.Count; r++) rows[p.Rows.Count + r][k] = op == 'D' ? '-' : q.Rows[r][qj];
            if (op != 'I') pi++;
            if (op != 'D') qj++;
        }
        return new Profile([.. p.Members, .. q.Members], rows);
    }
}
