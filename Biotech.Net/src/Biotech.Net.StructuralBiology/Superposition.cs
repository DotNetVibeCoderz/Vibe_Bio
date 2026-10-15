namespace Biotech.Net.StructuralBiology;

/// <summary>
/// A rigid-body transform <c>x → R·x + t</c> that superposes a mobile point set onto a target, with the RMSD
/// of the fitted pairs.
/// </summary>
public readonly record struct SuperpositionResult(Matrix3 Rotation, Vec3 Translation, double Rmsd)
{
    /// <summary>The identity transform.</summary>
    public static readonly SuperpositionResult Identity = new(Matrix3.Identity, Vec3.Zero, 0);

    /// <summary>Transforms one point.</summary>
    public Vec3 Apply(Vec3 p) => Rotation * p + Translation;

    /// <summary>Transforms points in place.</summary>
    public void Apply(Span<Vec3> points)
    {
        for (int i = 0; i < points.Length; i++) points[i] = Rotation * points[i] + Translation;
    }

    /// <summary>Transforms every atom of a model in place.</summary>
    public void Apply(Model model) => model.Transform(Rotation, Translation);

    /// <summary>Transforms the given atoms in place.</summary>
    public void Apply(IEnumerable<Atom> atoms)
    {
        foreach (var a in atoms) a.Position = Apply(a.Position);
    }
}

/// <summary>
/// Optimal rigid-body superposition (Kabsch 1976 least-squares fit, solved with Horn's 1987 unit-quaternion
/// eigenproblem, which never returns a reflection) and RMSD. Managed implementation.
/// </summary>
public static class Superposition
{
    /// <summary>RMSD between corresponding points without any fitting.</summary>
    public static double Rmsd(ReadOnlySpan<Vec3> a, ReadOnlySpan<Vec3> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Point sets differ in length.");
        if (a.IsEmpty) return 0;
        double s = 0;
        for (int i = 0; i < a.Length; i++) s += Vec3.DistanceSquared(a[i], b[i]);
        return Math.Sqrt(s / a.Length);
    }

    /// <summary>
    /// Finds the rotation and translation that minimise the (weighted) RMSD of <paramref name="mobile"/> onto
    /// <paramref name="target"/>; the result maps mobile coordinates into the target frame.
    /// </summary>
    /// <param name="weights">Optional non-negative weight per pair (empty = all 1).</param>
    public static SuperpositionResult Kabsch(ReadOnlySpan<Vec3> mobile, ReadOnlySpan<Vec3> target, ReadOnlySpan<double> weights = default)
    {
        int n = mobile.Length;
        if (n != target.Length) throw new ArgumentException("Point sets differ in length.");
        if (!weights.IsEmpty && weights.Length != n) throw new ArgumentException("One weight per point is required.", nameof(weights));
        if (n == 0) return SuperpositionResult.Identity;

        double wsum = 0;
        double mx = 0, my = 0, mz = 0, tx = 0, ty = 0, tz = 0;
        for (int i = 0; i < n; i++)
        {
            double w = weights.IsEmpty ? 1 : weights[i];
            wsum += w;
            mx += w * mobile[i].X; my += w * mobile[i].Y; mz += w * mobile[i].Z;
            tx += w * target[i].X; ty += w * target[i].Y; tz += w * target[i].Z;
        }
        if (wsum <= 0) throw new ArgumentException("Weights sum to zero.", nameof(weights));
        var cm = new Vec3(mx / wsum, my / wsum, mz / wsum);
        var ct = new Vec3(tx / wsum, ty / wsum, tz / wsum);

        double sxx = 0, sxy = 0, sxz = 0, syx = 0, syy = 0, syz = 0, szx = 0, szy = 0, szz = 0, g = 0;
        for (int i = 0; i < n; i++)
        {
            double w = weights.IsEmpty ? 1 : weights[i];
            var x = mobile[i] - cm; var y = target[i] - ct;
            sxx += w * x.X * y.X; sxy += w * x.X * y.Y; sxz += w * x.X * y.Z;
            syx += w * x.Y * y.X; syy += w * x.Y * y.Y; syz += w * x.Y * y.Z;
            szx += w * x.Z * y.X; szy += w * x.Z * y.Y; szz += w * x.Z * y.Z;
            g += w * (x.LengthSquared + y.LengthSquared);
        }

        // Horn's symmetric 4×4 key matrix; its top eigenvector is the optimal rotation quaternion.
        Span<double> k = stackalloc double[16];
        k[0] = sxx + syy + szz; k[1] = syz - szy; k[2] = szx - sxz; k[3] = sxy - syx;
        k[5] = sxx - syy - szz; k[6] = sxy + syx; k[7] = szx + sxz;
        k[10] = -sxx + syy - szz; k[11] = syz + szy;
        k[15] = -sxx - syy + szz;
        k[4] = k[1]; k[8] = k[2]; k[12] = k[3]; k[9] = k[6]; k[13] = k[7]; k[14] = k[11];

        Span<double> vectors = stackalloc double[16];
        Span<double> values = stackalloc double[4];
        JacobiEigen4(k, values, vectors);
        int best = 0;
        for (int i = 1; i < 4; i++) if (values[i] > values[best]) best = i;
        var r = Matrix3.FromQuaternion(vectors[best], vectors[4 + best], vectors[8 + best], vectors[12 + best]);
        double msd = Math.Max(0, (g - 2 * values[best]) / wsum);
        return new SuperpositionResult(r, ct - r * cm, Math.Sqrt(msd));
    }

    /// <summary>Superposes atoms pairwise (same count, paired by index).</summary>
    public static SuperpositionResult Kabsch(IReadOnlyList<Atom> mobile, IReadOnlyList<Atom> target)
        => Kabsch(Positions(mobile), Positions(target));

    /// <summary>RMSD after optimal superposition, without changing the inputs.</summary>
    public static double FittedRmsd(ReadOnlySpan<Vec3> mobile, ReadOnlySpan<Vec3> target) => Kabsch(mobile, target).Rmsd;

    /// <summary>
    /// Sequence-independent CA superposition of two chains with the same number of CA atoms: residues are paired
    /// by order, not by name or number, so mutants, homologues or renumbered chains of equal length can be fitted.
    /// </summary>
    public static SuperpositionResult SuperposeCa(Chain mobile, Chain target)
    {
        var a = CaTrace(mobile); var b = CaTrace(target);
        if (a.Length != b.Length) throw new ArgumentException($"CA counts differ ({a.Length} vs {b.Length}); use StructuralAlignment.Align.");
        return Kabsch(a, b);
    }

    /// <summary>The CA coordinates of a chain's amino acids in order.</summary>
    public static Vec3[] CaTrace(Chain chain) => chain.AminoAcids().Select(r => r["CA"].Position).ToArray();

    internal static Vec3[] Positions(IReadOnlyList<Atom> atoms)
    {
        var p = new Vec3[atoms.Count];
        for (int i = 0; i < p.Length; i++) p[i] = atoms[i].Position;
        return p;
    }

    /// <summary>Cyclic Jacobi eigen-decomposition of a symmetric 4×4 matrix (row-major); eigenvectors are columns.</summary>
    private static void JacobiEigen4(Span<double> a, Span<double> d, Span<double> v)
    {
        const int n = 4;
        v.Clear();
        for (int i = 0; i < n; i++) v[i * n + i] = 1;
        for (int sweep = 0; sweep < 100; sweep++)
        {
            double off = 0;
            for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) off += a[p * n + q] * a[p * n + q];
            if (off < 1e-22) break;
            for (int p = 0; p < n; p++)
                for (int q = p + 1; q < n; q++)
                {
                    double apq = a[p * n + q];
                    if (Math.Abs(apq) < 1e-300) continue;
                    double theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
                    double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    if (theta == 0) t = 1;
                    double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
                    for (int k = 0; k < n; k++)
                    {
                        double akp = a[k * n + p], akq = a[k * n + q];
                        a[k * n + p] = c * akp - s * akq;
                        a[k * n + q] = s * akp + c * akq;
                    }
                    for (int k = 0; k < n; k++)
                    {
                        double apk = a[p * n + k], aqk = a[q * n + k];
                        a[p * n + k] = c * apk - s * aqk;
                        a[q * n + k] = s * apk + c * aqk;
                    }
                    for (int k = 0; k < n; k++)
                    {
                        double vkp = v[k * n + p], vkq = v[k * n + q];
                        v[k * n + p] = c * vkp - s * vkq;
                        v[k * n + q] = s * vkp + c * vkq;
                    }
                }
        }
        for (int i = 0; i < n; i++) d[i] = a[i * n + i];
    }
}

/// <summary>Result of a TM-score evaluation or a structural alignment.</summary>
/// <param name="TmScore">TM-score (0–1], normalised by <paramref name="NormalizationLength"/>.</param>
/// <param name="Rmsd">RMSD of the aligned pairs under <paramref name="Transform"/>.</param>
/// <param name="D0">Distance scale d₀ = 1.24·∛(L − 15) − 1.8 (at least 0.5).</param>
/// <param name="AlignedPairs">Pairs (mobile index, target index) that were aligned.</param>
/// <param name="Transform">Superposition of the mobile set onto the target that maximises the TM-score.</param>
public sealed record TmScoreResult(double TmScore, double Rmsd, double D0, int NormalizationLength,
    IReadOnlyList<(int Mobile, int Target)> AlignedPairs, SuperpositionResult Transform)
{
    /// <summary>Aligned pairs within 5 Å after superposition.</summary>
    public int PairsWithin5A { get; init; }

    /// <summary>Sequence identity is not used; this is the fraction of target residues aligned.</summary>
    public double Coverage => NormalizationLength == 0 ? 0 : (double)AlignedPairs.Count / NormalizationLength;
}

/// <summary>
/// TM-score (Zhang &amp; Skolnick 2004, Proteins 57:702) and a TM-align-style sequence-independent structural
/// alignment (Zhang &amp; Skolnick 2005, NAR 33:2302) on CA traces. The superposition search follows the TM-score
/// program: seed fragments of decreasing length, then iteratively re-fit on pairs closer than d₀-search.
/// </summary>
public static class StructuralAlignment
{
    /// <summary>The TM-score distance scale for a target of <paramref name="length"/> residues.</summary>
    public static double D0(int length) => length <= 21 ? 0.5 : Math.Max(0.5, 1.24 * Math.Cbrt(length - 15) - 1.8);

    /// <summary>
    /// TM-score of <paramref name="model"/> against <paramref name="native"/> with residues paired by index (equal
    /// lengths), maximised over rigid superpositions. Normalised by the native length unless
    /// <paramref name="normalizationLength"/> is given.
    /// </summary>
    public static TmScoreResult TmScore(ReadOnlySpan<Vec3> model, ReadOnlySpan<Vec3> native, int? normalizationLength = null)
    {
        if (model.Length != native.Length) throw new ArgumentException("TM-score with fixed correspondence needs equal lengths.");
        var pairs = new (int, int)[model.Length];
        for (int i = 0; i < pairs.Length; i++) pairs[i] = (i, i);
        return Score(model, native, pairs, normalizationLength ?? native.Length);
    }

    /// <summary>TM-score of two chains of equal CA count, paired by order.</summary>
    public static TmScoreResult TmScore(Chain model, Chain native) => TmScore(Superposition.CaTrace(model), Superposition.CaTrace(native));

    /// <summary>
    /// Sequence-independent structural alignment of two CA traces of any lengths (TM-align-like: gapless threading
    /// seeds, then dynamic programming on the TM-score similarity matrix iterated with re-superposition).
    /// The TM-score is normalised by the target length.
    /// </summary>
    public static TmScoreResult Align(ReadOnlySpan<Vec3> mobile, ReadOnlySpan<Vec3> target)
    {
        int la = mobile.Length, lb = target.Length;
        if (la < 3 || lb < 3) throw new ArgumentException("Need at least three residues on each side.");
        double d0 = D0(lb);
        var a = mobile.ToArray(); var b = target.ToArray();

        // 1. Gapless threading seeds: every offset with an overlap of at least half the shorter chain.
        TmScoreResult? best = null;
        int minOverlap = Math.Max(3, Math.Min(la, lb) / 2);
        for (int shift = -(la - minOverlap); shift <= lb - minOverlap; shift++)
        {
            var pairs = new List<(int, int)>();
            for (int i = 0; i < la; i++)
            {
                int j = i + shift;
                if (j >= 0 && j < lb) pairs.Add((i, j));
            }
            if (pairs.Count < minOverlap) continue;
            var r = Score(a, b, pairs, lb, quick: true);
            if (best is null || r.TmScore > best.TmScore) best = r;
        }
        best ??= Score(a, b, [(0, 0), (1, 1), (2, 2)], lb);

        // 2. Iterate: superpose → DP on similarity 1/(1+(d/d0)²) with gap-open −0.6 → re-score.
        var transform = best.Transform;
        var lastPairs = best.AlignedPairs;
        for (int iter = 0; iter < 30; iter++)
        {
            var pairs = DynamicProgramming(a, b, transform, d0, gapOpen: -0.6);
            if (pairs.Count < 3) break;
            var r = Score(a, b, pairs, lb);
            if (r.TmScore > best.TmScore + 1e-9) best = r;
            if (pairs.SequenceEqual(lastPairs)) break;
            lastPairs = pairs;
            transform = r.Transform;
        }
        return best;
    }

    /// <summary>Structural alignment of the CA traces of two chains.</summary>
    public static TmScoreResult Align(Chain mobile, Chain target) => Align(Superposition.CaTrace(mobile), Superposition.CaTrace(target));

    private static List<(int, int)> DynamicProgramming(Vec3[] a, Vec3[] b, SuperpositionResult t, double d0, double gapOpen)
    {
        int la = a.Length, lb = b.Length;
        var moved = new Vec3[la];
        for (int i = 0; i < la; i++) moved[i] = t.Apply(a[i]);
        var h = new double[(la + 1) * (lb + 1)];
        var dir = new byte[(la + 1) * (lb + 1)]; // 0 diag, 1 up (gap in b), 2 left (gap in a)
        var open = new bool[(la + 1) * (lb + 1)];
        int w = lb + 1;
        double d02 = d0 * d0;
        for (int i = 1; i <= la; i++)
        {
            for (int j = 1; j <= lb; j++)
            {
                double s = 1 / (1 + Vec3.DistanceSquared(moved[i - 1], b[j - 1]) / d02);
                double diag = h[(i - 1) * w + j - 1] + s;
                double up = h[(i - 1) * w + j] + (open[(i - 1) * w + j] ? gapOpen : 0);
                double left = h[i * w + j - 1] + (open[i * w + j - 1] ? gapOpen : 0);
                int idx = i * w + j;
                if (diag >= up && diag >= left) { h[idx] = diag; dir[idx] = 0; open[idx] = true; }
                else if (up >= left) { h[idx] = up; dir[idx] = 1; open[idx] = false; }
                else { h[idx] = left; dir[idx] = 2; open[idx] = false; }
            }
        }
        var pairs = new List<(int, int)>();
        int x = la, y = lb;
        while (x > 0 && y > 0)
        {
            byte d = dir[x * w + y];
            if (d == 0) { pairs.Add((x - 1, y - 1)); x--; y--; }
            else if (d == 1) x--;
            else y--;
        }
        pairs.Reverse();
        return pairs;
    }

    /// <summary>Maximises the TM-score over superpositions for a fixed list of aligned pairs.</summary>
    private static TmScoreResult Score(ReadOnlySpan<Vec3> a, ReadOnlySpan<Vec3> b, IReadOnlyList<(int Mobile, int Target)> pairs, int lnorm, bool quick = false)
    {
        int n = pairs.Count;
        double d0 = D0(lnorm);
        double dSearch = Math.Clamp(d0, 4.5, 8.0);
        var pa = new Vec3[n]; var pb = new Vec3[n];
        for (int i = 0; i < n; i++) { pa[i] = a[pairs[i].Mobile]; pb[i] = b[pairs[i].Target]; }

        double bestTm = -1;
        SuperpositionResult bestT = SuperpositionResult.Identity;
        var dist = new double[n];
        var selected = new List<int>(n);
        var subA = new Vec3[n]; var subB = new Vec3[n];
        int minFrag = Math.Min(n, 4);
        int maxIter = quick ? 2 : 20;
        for (int frag = n; frag >= minFrag; frag /= 2)
        {
            int step = quick ? Math.Max(1, frag) : Math.Max(1, frag / 2);
            for (int start = 0; start + frag <= n; start += step)
            {
                selected.Clear();
                for (int i = start; i < start + frag; i++) selected.Add(i);
                for (int iter = 0; iter < maxIter; iter++)
                {
                    for (int s = 0; s < selected.Count; s++) { subA[s] = pa[selected[s]]; subB[s] = pb[selected[s]]; }
                    var t = Superposition.Kabsch(subA.AsSpan(0, selected.Count), subB.AsSpan(0, selected.Count));
                    double tm = 0;
                    for (int i = 0; i < n; i++)
                    {
                        dist[i] = Vec3.Distance(t.Apply(pa[i]), pb[i]);
                        tm += 1 / (1 + dist[i] * dist[i] / (d0 * d0));
                    }
                    tm /= lnorm;
                    if (tm > bestTm) { bestTm = tm; bestT = t; }
                    double cut = dSearch;
                    var next = new List<int>();
                    while (true)
                    {
                        next.Clear();
                        for (int i = 0; i < n; i++) if (dist[i] < cut) next.Add(i);
                        if (next.Count >= Math.Min(3, n)) break;
                        cut += 0.5;
                    }
                    if (next.SequenceEqual(selected)) break;
                    selected.Clear();
                    selected.AddRange(next);
                }
            }
            if (frag == minFrag) break;
            if (frag / 2 < minFrag) frag = minFrag * 2;
        }

        double ss = 0; int within5 = 0;
        for (int i = 0; i < n; i++)
        {
            double d2 = Vec3.DistanceSquared(bestT.Apply(pa[i]), pb[i]);
            ss += d2;
            if (d2 <= 25) within5++;
        }
        var final = bestT with { Rmsd = n == 0 ? 0 : Math.Sqrt(ss / n) };
        return new TmScoreResult(bestTm, final.Rmsd, d0, lnorm, pairs.ToArray(), final) { PairsWithin5A = within5 };
    }
}
