namespace Biotech.Net.StructuralBiology;

/// <summary>Solvent accessibility of one residue.</summary>
/// <param name="Absolute">Accessible surface area in Å².</param>
/// <param name="Relative">Absolute / maximum ASA of the residue type (Tien et al. 2013), NaN for non-amino acids.</param>
public readonly record struct ResidueSasa(Residue Residue, double Absolute, double Relative)
{
    /// <summary>Conventional buried/exposed call: relative accessibility below 0.25 (25%).</summary>
    public bool IsBuried => Relative < 0.25;
}

/// <summary>
/// Solvent accessible surface area by the Shrake–Rupley (1973) dot method: each atom is a sphere of radius
/// r(vdW) + probe carrying <c>n</c> golden-spiral test points; a point is accessible when it lies outside every
/// neighbouring sphere. Neighbours come from a cell list, so the cost is O(n·points·k) for k local neighbours.
/// The point set and the default radii are those of Biopython's <c>Bio.PDB.SASA.ShrakeRupley</c>, so results agree
/// with it to within rounding.
/// </summary>
public static class ShrakeRupley
{
    /// <summary>Unit-sphere test points on a golden spiral (Biopython's construction).</summary>
    public static Vec3[] SpherePoints(int n)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(n, 1);
        var pts = new Vec3[n];
        double dl = Math.PI * (3 - Math.Sqrt(5)), dz = 2.0 / n, lon = 0, z = 1 - dz / 2;
        for (int k = 0; k < n; k++)
        {
            double r = Math.Sqrt(1 - z * z);
            pts[k] = new Vec3(Math.Cos(lon) * r, Math.Sin(lon) * r, z);
            z -= dz;
            lon += dl;
        }
        return pts;
    }

    /// <summary>Per-atom accessible area (Å²) for explicit centres and radii.</summary>
    public static double[] Compute(ReadOnlySpan<Vec3> centers, ReadOnlySpan<double> radii, double probe = 1.4, int points = 100)
    {
        if (centers.Length != radii.Length) throw new ArgumentException("One radius per centre is required.");
        int n = centers.Length;
        var result = new double[n];
        if (n == 0) return result;
        var sphere = SpherePoints(points);
        var c = centers.ToArray();
        var rp = new double[n];
        double maxR = 0;
        for (int i = 0; i < n; i++) { rp[i] = radii[i] + probe; maxR = Math.Max(maxR, rp[i]); }
        var ns = new NeighborSearch(c, 2 * maxR);
        Parallel.For(0, n, () => (new List<int>(64), new List<int>(64)), (i, _, buffers) =>
        {
            var (found, neighbours) = buffers;
            found.Clear(); neighbours.Clear();
            ns.Query(c[i], rp[i] + maxR, found);
            foreach (int j in found)
                if (j != i && Vec3.Distance(c[i], c[j]) < rp[i] + rp[j]) neighbours.Add(j);
            int accessible = 0;
            int last = 0; // index into neighbours of the last occluder: neighbouring test points are often hidden by the same atom
            foreach (var p in sphere)
            {
                var q = c[i] + p * rp[i];
                bool buried = false;
                if (neighbours.Count > 0)
                {
                    int j0 = neighbours[last];
                    if (Vec3.DistanceSquared(q, c[j0]) <= rp[j0] * rp[j0]) buried = true;
                    else
                        for (int k = 0; k < neighbours.Count; k++)
                        {
                            int j = neighbours[k];
                            if (Vec3.DistanceSquared(q, c[j]) <= rp[j] * rp[j]) { buried = true; last = k; break; }
                        }
                }
                if (!buried) accessible++;
            }
            result[i] = 4 * Math.PI * rp[i] * rp[i] * accessible / sphere.Length;
            return buffers;
        }, _ => { });
        return result;
    }

    /// <summary>Per-atom accessible area (Å²) of the given atoms (radii from <see cref="ElementData.VdwRadius"/>).</summary>
    public static double[] Compute(IReadOnlyList<Atom> atoms, double probe = 1.4, int points = 100)
    {
        var centers = Superposition.Positions(atoms);
        var radii = new double[atoms.Count];
        for (int i = 0; i < radii.Length; i++) radii[i] = atoms[i].VdwRadius;
        return Compute(centers, radii, probe, points);
    }

    /// <summary>Total accessible area of a model (all atoms, including waters and ligands unless removed).</summary>
    public static double Total(Model model, double probe = 1.4, int points = 100) => Compute(model.GetAtoms(), probe, points).Sum();

    /// <summary>Per-residue absolute and relative accessibility of a model.</summary>
    /// <param name="excludeWater">Drop water molecules before computing (default true).</param>
    public static IReadOnlyList<ResidueSasa> PerResidue(Model model, double probe = 1.4, int points = 100, bool excludeWater = true)
    {
        var atoms = model.Atoms.Where(a => !excludeWater || a.Residue is not { IsWater: true }).ToList();
        var sasa = Compute(atoms, probe, points);
        var sums = new Dictionary<Residue, double>();
        var order = new List<Residue>();
        for (int i = 0; i < atoms.Count; i++)
        {
            var r = atoms[i].Residue!;
            if (!sums.ContainsKey(r)) { sums[r] = 0; order.Add(r); }
            sums[r] += sasa[i];
        }
        return order.Select(r => new ResidueSasa(r, sums[r],
            r.IsAminoAcid && ResidueNames.MaxAsaTien2013.TryGetValue(r.OneLetterCode, out double max) ? sums[r] / max : double.NaN)).ToList();
    }
}
