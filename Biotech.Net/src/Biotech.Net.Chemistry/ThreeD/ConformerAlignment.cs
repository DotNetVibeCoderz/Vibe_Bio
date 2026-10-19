namespace Biotech.Net.Chemistry;

/// <summary>
/// RMSD and optimal superposition of 3-D coordinate sets (quaternion method of Horn, J. Opt. Soc. Am. A 1987, 4, 629),
/// with an optional symmetry-aware variant that tries every graph automorphism of the heavy atoms.
/// </summary>
public static class ConformerAlignment
{
    /// <summary>
    /// RMSD between two equally long point sets. With <paramref name="align"/> the second set is first optimally
    /// rotated and translated onto the first; otherwise the coordinates are compared in place.
    /// </summary>
    public static double Rmsd(ReadOnlySpan<Point3> a, ReadOnlySpan<Point3> b, bool align = true)
    {
        if (a.Length != b.Length) throw new ArgumentException("Point sets must have the same length.");
        int n = a.Length;
        if (n == 0) return 0;
        if (!align)
        {
            double s = 0;
            for (int i = 0; i < n; i++) { var d = a[i] - b[i]; s += Point3.Dot(d, d); }
            return Math.Sqrt(s / n);
        }
        return Math.Sqrt(Math.Max(0, SuperposeSquared(a, b, out _, out _, out _) / n));
    }

    /// <summary>
    /// Superposes <paramref name="mobile"/> onto <paramref name="reference"/> and returns the transformed copy of
    /// <paramref name="mobile"/> (rotation about its centroid, then translation onto the reference centroid).
    /// </summary>
    public static Point3[] Superpose(ReadOnlySpan<Point3> reference, ReadOnlySpan<Point3> mobile)
    {
        SuperposeSquared(reference, mobile, out var rot, out var cRef, out var cMob);
        var result = new Point3[mobile.Length];
        for (int i = 0; i < mobile.Length; i++) result[i] = Rotate(rot, mobile[i] - cMob) + cRef;
        return result;
    }

    /// <summary>
    /// Heavy-atom RMSD between two conformations of the same molecule, minimised over the molecule's heavy-atom
    /// graph automorphisms (so symmetric groups such as phenyl rings or carboxylates do not inflate it).
    /// </summary>
    /// <param name="mol">The molecule (its atom order defines the coordinate order).</param>
    /// <param name="a">First conformation (one point per atom of <paramref name="mol"/>).</param>
    /// <param name="b">Second conformation.</param>
    /// <param name="align">Superpose before measuring (conformer comparison) or compare in place (docking poses).</param>
    /// <param name="maxAutomorphisms">Cap on the number of symmetry mappings tried.</param>
    public static double BestRmsd(Molecule mol, ReadOnlySpan<Point3> a, ReadOnlySpan<Point3> b, bool align = true, int maxAutomorphisms = 1000)
    {
        var maps = HeavyAtomAutomorphisms(mol, maxAutomorphisms);
        double best = double.MaxValue;
        var pa = new Point3[maps[0].Length];
        var pb = new Point3[maps[0].Length];
        foreach (var map in maps)
        {
            for (int k = 0; k < map.Length; k++) { pa[k] = a[maps[0][k]]; pb[k] = b[map[k]]; }
            double r = Rmsd(pa, pb, align);
            if (r < best) best = r;
        }
        return best;
    }

    /// <summary>
    /// Heavy-atom automorphisms of <paramref name="mol"/>: each array maps heavy-atom position k (in index order) to
    /// the atom index it may be exchanged with. The identity comes first.
    /// </summary>
    public static IReadOnlyList<int[]> HeavyAtomAutomorphisms(Molecule mol, int maxAutomorphisms = 1000)
    {
        var heavy = Enumerable.Range(0, mol.AtomCount).Where(i => mol.Atoms[i].AtomicNumber != 1).ToArray();
        var identity = heavy.ToArray();
        if (heavy.Length == 0) return [identity];
        var stripped = mol.AtomCount == heavy.Length ? mol : mol.RemoveHydrogens();
        if (stripped.AtomCount != heavy.Length) return [identity];
        var query = Substructure.FromMolecule(stripped);
        var result = new List<int[]> { identity };
        var seen = new HashSet<string> { string.Join(",", identity) };
        foreach (var m in query.Matches(stripped, unique: false, maxMatches: maxAutomorphisms))
        {
            var mapped = m.Select(x => heavy[x]).ToArray();
            if (seen.Add(string.Join(",", mapped))) result.Add(mapped);
        }
        return result;
    }

    // ---- quaternion superposition -------------------------------------------------------------------

    private static double SuperposeSquared(ReadOnlySpan<Point3> a, ReadOnlySpan<Point3> b, out double[] rot, out Point3 ca, out Point3 cb)
    {
        int n = a.Length;
        ca = default; cb = default;
        for (int i = 0; i < n; i++) { ca += a[i]; cb += b[i]; }
        ca *= 1.0 / n; cb *= 1.0 / n;
        double sxx = 0, sxy = 0, sxz = 0, syx = 0, syy = 0, syz = 0, szx = 0, szy = 0, szz = 0, ga = 0, gb = 0;
        for (int i = 0; i < n; i++)
        {
            var p = a[i] - ca;
            var q = b[i] - cb;
            ga += Point3.Dot(p, p);
            gb += Point3.Dot(q, q);
            sxx += q.X * p.X; sxy += q.X * p.Y; sxz += q.X * p.Z;
            syx += q.Y * p.X; syy += q.Y * p.Y; syz += q.Y * p.Z;
            szx += q.Z * p.X; szy += q.Z * p.Y; szz += q.Z * p.Z;
        }
        double[,] k =
        {
            { sxx + syy + szz, syz - szy, szx - sxz, sxy - syx },
            { syz - szy, sxx - syy - szz, sxy + syx, szx + sxz },
            { szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy },
            { sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz },
        };
        var (lambda, v) = LargestEigen4(k);
        rot = QuaternionToMatrix(v[0], v[1], v[2], v[3]);
        return ga + gb - 2 * lambda;
    }

    private static Point3 Rotate(double[] r, Point3 p)
        => new(r[0] * p.X + r[1] * p.Y + r[2] * p.Z, r[3] * p.X + r[4] * p.Y + r[5] * p.Z, r[6] * p.X + r[7] * p.Y + r[8] * p.Z);

    private static double[] QuaternionToMatrix(double w, double x, double y, double z)
    {
        double n = Math.Sqrt(w * w + x * x + y * y + z * z);
        if (n < 1e-12) return [1, 0, 0, 0, 1, 0, 0, 0, 1];
        w /= n; x /= n; y /= n; z /= n;
        return
        [
            w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y),
            2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x),
            2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z,
        ];
    }

    /// <summary>Largest eigenvalue/eigenvector of a symmetric 4×4 matrix (cyclic Jacobi).</summary>
    private static (double Value, double[] Vector) LargestEigen4(double[,] m)
    {
        var a = (double[,])m.Clone();
        var v = new double[4, 4];
        for (int i = 0; i < 4; i++) v[i, i] = 1;
        for (int sweep = 0; sweep < 60; sweep++)
        {
            double off = 0;
            for (int p = 0; p < 4; p++) for (int q = p + 1; q < 4; q++) off += a[p, q] * a[p, q];
            if (off < 1e-24) break;
            for (int p = 0; p < 4; p++)
                for (int q = p + 1; q < 4; q++)
                {
                    if (Math.Abs(a[p, q]) < 1e-300) continue;
                    double theta = (a[q, q] - a[p, p]) / (2 * a[p, q]);
                    double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    if (theta == 0) t = 1;
                    double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
                    for (int k = 0; k < 4; k++)
                    {
                        double akp = a[k, p], akq = a[k, q];
                        a[k, p] = c * akp - s * akq;
                        a[k, q] = s * akp + c * akq;
                    }
                    for (int k = 0; k < 4; k++)
                    {
                        double apk = a[p, k], aqk = a[q, k];
                        a[p, k] = c * apk - s * aqk;
                        a[q, k] = s * apk + c * aqk;
                    }
                    for (int k = 0; k < 4; k++)
                    {
                        double vkp = v[k, p], vkq = v[k, q];
                        v[k, p] = c * vkp - s * vkq;
                        v[k, q] = s * vkp + c * vkq;
                    }
                }
        }
        int best = 0;
        for (int i = 1; i < 4; i++) if (a[i, i] > a[best, best]) best = i;
        return (a[best, best], [v[0, best], v[1, best], v[2, best], v[3, best]]);
    }
}
