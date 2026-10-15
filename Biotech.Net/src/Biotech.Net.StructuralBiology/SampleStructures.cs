namespace Biotech.Net.StructuralBiology;

/// <summary>
/// Small structures embedded in the library for demos, notebooks and tests (no network needed).
/// </summary>
public static class SampleStructures
{
    private static string? _crambinText;

    /// <summary>
    /// Crambin (PDB 1CRN; Teeter &amp; Hendrickson, 0.945 Å): 46 residues, 327 atoms, three disulfides
    /// (Cys3–Cys40, Cys4–Cys32, Cys16–Cys26), two helices and a small antiparallel β-sheet. Returns a fresh,
    /// independent copy on every call so callers may transform it.
    /// </summary>
    public static Structure Crambin() => PdbReader.Parse(CrambinPdb, "1CRN");

    /// <summary>The embedded 1CRN coordinate file (HEADER, TITLE, CRYST1, HELIX, SHEET, SSBOND, ATOM records).</summary>
    public static string CrambinPdb => _crambinText ??= LoadResource("Biotech.Net.StructuralBiology.Data.1crn.pdb");

    /// <summary>
    /// An ideal right-handed α-helical poly-alanine backbone of <paramref name="length"/> residues (N, CA, C, O, CB)
    /// built from standard geometry with φ = −57°, ψ = −47°, ω = 180°. Handy for geometry and DSSP demos.
    /// </summary>
    public static Structure AlanineHelix(int length = 12) => PeptideBuilder.Build("ALA_HELIX", length, -57, -47);

    /// <summary>An extended poly-alanine β-strand (φ = −120°, ψ = +130°).</summary>
    public static Structure AlanineStrand(int length = 8) => PeptideBuilder.Build("ALA_STRAND", length, -120, 130);

    private static string LoadResource(string name)
    {
        using var stream = typeof(SampleStructures).Assembly.GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded resource '{name}' is missing.");
        using var reader = new StreamReader(stream);
        return reader.ReadToEnd();
    }
}

/// <summary>Builds peptide backbones from internal coordinates (NeRF: natural extension reference frame).</summary>
public static class PeptideBuilder
{
    // Engh & Huber 1991 ideal geometry.
    private const double NCa = 1.458, CaC = 1.525, CN = 1.329, CO = 1.231, CaCb = 1.530;
    private const double AngNCaC = 111.2, AngCaCN = 116.2, AngCNCa = 121.7, AngCaCO = 120.5;

    /// <summary>
    /// Places atom D given A, B, C, the bond length |CD|, the angle BCD and the dihedral ABCD (degrees).
    /// </summary>
    public static Vec3 PlaceAtom(Vec3 a, Vec3 b, Vec3 c, double bond, double angle, double torsion)
    {
        double th = angle * Math.PI / 180, ph = torsion * Math.PI / 180;
        var d2 = new Vec3(-bond * Math.Cos(th), bond * Math.Sin(th) * Math.Cos(ph), bond * Math.Sin(th) * Math.Sin(ph));
        var bc = (c - b).Normalized();
        var n = Vec3.Cross(b - a, bc).Normalized();
        var m = Vec3.Cross(n, bc);
        return c + bc * d2.X + m * d2.Y + n * d2.Z;
    }

    /// <summary>Builds a poly-alanine chain A with constant φ/ψ (ω = 180°).</summary>
    public static Structure Build(string id, int length, double phi, double psi, double omega = 180)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
        var structure = new Structure(id);
        var model = new Model(1);
        structure.AddModel(model);
        var chain = new Chain("A");
        model.AddChain(chain);
        var n = new Vec3(0, 0, 0);
        var ca = new Vec3(NCa, 0, 0);
        var c = ca + new Vec3(-CaC * Math.Cos(AngNCaC * Math.PI / 180), CaC * Math.Sin(AngNCaC * Math.PI / 180), 0);
        int serial = 0;
        for (int i = 0; i < length; i++)
        {
            if (i > 0)
            {
                var prevN = n; var prevCa = ca; var prevC = c;
                n = PlaceAtom(prevN, prevCa, prevC, CN, AngCaCN, psi);
                ca = PlaceAtom(prevCa, prevC, n, NCa, AngCNCa, omega);
                c = PlaceAtom(prevC, n, ca, CaC, AngNCaC, phi);
            }
            // Carbonyl O lies trans to the next N (ψ + 180); CB from the L-amino-acid improper.
            var o = PlaceAtom(n, ca, c, CO, AngCaCO, psi + 180);
            var cb = PlaceAtom(c, n, ca, CaCb, 110.5, -122.6);
            var res = new Residue("ALA", i + 1);
            res.AddAtom(new Atom("N", "N", n) { Serial = ++serial });
            res.AddAtom(new Atom("CA", "C", ca) { Serial = ++serial });
            res.AddAtom(new Atom("C", "C", c) { Serial = ++serial });
            res.AddAtom(new Atom("O", "O", o) { Serial = ++serial });
            res.AddAtom(new Atom("CB", "C", cb) { Serial = ++serial });
            chain.AddResidue(res);
        }
        return structure;
    }
}
