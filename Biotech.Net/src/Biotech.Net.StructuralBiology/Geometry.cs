namespace Biotech.Net.StructuralBiology;

/// <summary>Broad Ramachandran regions (rectangular approximations of the populated φ/ψ basins).</summary>
public enum RamachandranRegion
{
    /// <summary>φ or ψ undefined (chain termini or breaks).</summary>
    Undefined,

    /// <summary>Right-handed α-helix basin: −160° ≤ φ ≤ −20°, −120° ≤ ψ ≤ 50°.</summary>
    AlphaRight,

    /// <summary>β-sheet basin: −180° ≤ φ ≤ −45°, ψ ≥ 90° or ψ ≤ −150° (excluding the PPII box).</summary>
    Beta,

    /// <summary>Polyproline II: −95° ≤ φ ≤ −50°, 110° ≤ ψ ≤ 180°.</summary>
    PolyprolineII,

    /// <summary>Left-handed helix basin (common for Gly, Asn, Asp): 20° ≤ φ ≤ 120°, −20° ≤ ψ ≤ 90°.</summary>
    AlphaLeft,

    /// <summary>Everything else (rare, strained or glycine-only regions).</summary>
    Other,
}

/// <summary>Backbone torsions of one residue in degrees (NaN when undefined), IUPAC definitions:
/// φᵢ = C(i−1)–N–CA–C, ψᵢ = N–CA–C–N(i+1), ωᵢ = CA–C–N(i+1)–CA(i+1).</summary>
public readonly record struct BackboneTorsion(Residue Residue, double Phi, double Psi, double Omega)
{
    /// <summary>The Ramachandran region of (φ, ψ).</summary>
    public RamachandranRegion Region => Ramachandran.Classify(Phi, Psi);

    /// <summary>True when ω is within 30° of 0° (a cis peptide bond to the next residue).</summary>
    public bool IsCis => !double.IsNaN(Omega) && Math.Abs(Omega) < 30;
}

/// <summary>Ramachandran classification helpers.</summary>
public static class Ramachandran
{
    /// <summary>Classifies a (φ, ψ) pair in degrees into a broad region.</summary>
    public static RamachandranRegion Classify(double phi, double psi)
    {
        if (double.IsNaN(phi) || double.IsNaN(psi)) return RamachandranRegion.Undefined;
        if (phi is >= -95 and <= -50 && psi is >= 110 and <= 180) return RamachandranRegion.PolyprolineII;
        if (phi is >= -160 and <= -20 && psi is >= -120 and <= 50) return RamachandranRegion.AlphaRight;
        if (phi is >= -180 and <= -45 && (psi >= 90 || psi <= -150)) return RamachandranRegion.Beta;
        if (phi is >= 20 and <= 120 && psi is >= -20 and <= 90) return RamachandranRegion.AlphaLeft;
        return RamachandranRegion.Other;
    }

    /// <summary>Counts residues per region for a chain (termini are <see cref="RamachandranRegion.Undefined"/>).</summary>
    public static Dictionary<RamachandranRegion, int> Summarize(Chain chain)
    {
        var d = new Dictionary<RamachandranRegion, int>();
        foreach (var t in Geometry.BackboneTorsions(chain)) d[t.Region] = d.GetValueOrDefault(t.Region) + 1;
        return d;
    }
}

/// <summary>Distances, angles, torsions and whole-structure shape descriptors.</summary>
public static class Geometry
{
    /// <summary>Peptide C–N distance above which two consecutive residues are considered not bonded.</summary>
    public const double PeptideBondCutoff = 2.0;

    /// <summary>Distance between two atoms in Å.</summary>
    public static double Distance(Atom a, Atom b) => Vec3.Distance(a.Position, b.Position);

    /// <summary>Angle a–b–c in degrees (b is the vertex).</summary>
    public static double Angle(Vec3 a, Vec3 b, Vec3 c) => Vec3.AngleBetween(a - b, c - b);

    /// <summary>Angle between three atoms in degrees.</summary>
    public static double Angle(Atom a, Atom b, Atom c) => Angle(a.Position, b.Position, c.Position);

    /// <summary>
    /// Dihedral (torsion) angle a–b–c–d in degrees in (−180°, 180°], IUPAC sign convention (clockwise positive
    /// looking down b→c), the same convention as Biopython's <c>calc_dihedral</c>.
    /// </summary>
    public static double Dihedral(Vec3 a, Vec3 b, Vec3 c, Vec3 d)
    {
        var b1 = b - a; var b2 = c - b; var b3 = d - c;
        var n1 = Vec3.Cross(b1, b2);
        var n2 = Vec3.Cross(b2, b3);
        double x = Vec3.Dot(n1, n2);
        double y = Vec3.Dot(Vec3.Cross(n1, n2), b2.Normalized());
        return Math.Atan2(y, x) * (180 / Math.PI);
    }

    /// <summary>Dihedral angle between four atoms in degrees.</summary>
    public static double Dihedral(Atom a, Atom b, Atom c, Atom d) => Dihedral(a.Position, b.Position, c.Position, d.Position);

    /// <summary>Unweighted centroid of the points.</summary>
    public static Vec3 Centroid(ReadOnlySpan<Vec3> points)
    {
        if (points.IsEmpty) throw new ArgumentException("No points.", nameof(points));
        double x = 0, y = 0, z = 0;
        foreach (var p in points) { x += p.X; y += p.Y; z += p.Z; }
        return new Vec3(x / points.Length, y / points.Length, z / points.Length);
    }

    /// <summary>Mass-weighted centre (IUPAC average masses) of the atoms.</summary>
    public static Vec3 CenterOfMass(IEnumerable<Atom> atoms)
    {
        double x = 0, y = 0, z = 0, m = 0;
        foreach (var a in atoms)
        {
            double w = a.Mass;
            x += w * a.Position.X; y += w * a.Position.Y; z += w * a.Position.Z; m += w;
        }
        if (m == 0) throw new ArgumentException("No atoms.", nameof(atoms));
        return new Vec3(x / m, y / m, z / m);
    }

    /// <summary>Total mass of the atoms in daltons.</summary>
    public static double TotalMass(IEnumerable<Atom> atoms) => atoms.Sum(a => a.Mass);

    /// <summary>Radius of gyration in Å, mass-weighted by default (geometric when <paramref name="massWeighted"/> is false).</summary>
    public static double RadiusOfGyration(IEnumerable<Atom> atoms, bool massWeighted = true)
    {
        var list = atoms as IReadOnlyList<Atom> ?? atoms.ToList();
        if (list.Count == 0) throw new ArgumentException("No atoms.", nameof(atoms));
        Vec3 c;
        if (massWeighted) c = CenterOfMass(list);
        else
        {
            var pts = new Vec3[list.Count];
            for (int i = 0; i < pts.Length; i++) pts[i] = list[i].Position;
            c = Centroid(pts);
        }
        double s = 0, m = 0;
        foreach (var a in list)
        {
            double w = massWeighted ? a.Mass : 1;
            s += w * Vec3.DistanceSquared(a.Position, c);
            m += w;
        }
        return Math.Sqrt(s / m);
    }

    /// <summary>Radius of gyration of bare points (unweighted).</summary>
    public static double RadiusOfGyration(ReadOnlySpan<Vec3> points)
    {
        var c = Centroid(points);
        double s = 0;
        foreach (var p in points) s += Vec3.DistanceSquared(p, c);
        return Math.Sqrt(s / points.Length);
    }

    /// <summary>Axis-aligned bounding box of the points.</summary>
    public static (Vec3 Min, Vec3 Max) BoundingBox(ReadOnlySpan<Vec3> points)
    {
        if (points.IsEmpty) throw new ArgumentException("No points.", nameof(points));
        Vec3 min = points[0], max = points[0];
        foreach (var p in points) { min = Vec3.Min(min, p); max = Vec3.Max(max, p); }
        return (min, max);
    }

    /// <summary>Axis-aligned bounding box of the atoms.</summary>
    public static (Vec3 Min, Vec3 Max) BoundingBox(IEnumerable<Atom> atoms) => BoundingBox(atoms.Select(a => a.Position).ToArray());

    /// <summary>
    /// φ/ψ/ω for every amino-acid residue of the chain that has N, CA and C atoms. Torsions spanning a chain break
    /// (C–N longer than <see cref="PeptideBondCutoff"/>) are NaN.
    /// </summary>
    public static IReadOnlyList<BackboneTorsion> BackboneTorsions(Chain chain)
    {
        var res = chain.Residues.Where(r => r.IsAminoAcid && r.HasAtom("N") && r.HasAtom("CA") && r.HasAtom("C")).ToList();
        var result = new List<BackboneTorsion>(res.Count);
        for (int i = 0; i < res.Count; i++)
        {
            var r = res[i];
            Vec3 n = r["N"].Position, ca = r["CA"].Position, c = r["C"].Position;
            double phi = double.NaN, psi = double.NaN, omega = double.NaN;
            if (i > 0 && Bonded(res[i - 1], r))
                phi = Dihedral(res[i - 1]["C"].Position, n, ca, c);
            if (i + 1 < res.Count && Bonded(r, res[i + 1]))
            {
                var next = res[i + 1];
                psi = Dihedral(n, ca, c, next["N"].Position);
                omega = Dihedral(ca, c, next["N"].Position, next["CA"].Position);
            }
            result.Add(new BackboneTorsion(r, phi, psi, omega));
        }
        return result;
    }

    /// <summary>True when residue <paramref name="a"/>'s C is within peptide-bond distance of <paramref name="b"/>'s N.</summary>
    public static bool Bonded(Residue a, Residue b)
    {
        var c = a.GetAtom("C"); var n = b.GetAtom("N");
        return c is not null && n is not null && Vec3.Distance(c.Position, n.Position) <= PeptideBondCutoff;
    }

    /// <summary>Side-chain χ1 (N–CA–CB–XG) in degrees, or NaN for Gly/Ala or missing atoms.</summary>
    public static double Chi1(Residue r)
    {
        string? g = r.Name switch
        {
            "SER" => "OG", "CYS" => "SG", "THR" => "OG1", "VAL" => "CG1", "ILE" => "CG1",
            "GLY" or "ALA" => null, _ => "CG",
        };
        if (g is null) return double.NaN;
        var n = r.GetAtom("N"); var ca = r.GetAtom("CA"); var cb = r.GetAtom("CB"); var x = r.GetAtom(g);
        return n is null || ca is null || cb is null || x is null ? double.NaN : Dihedral(n, ca, cb, x);
    }
}
