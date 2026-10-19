using System.Text.Json;
using Biotech.Net.Chemistry;

namespace Biotech.Net.Chemistry.Tests;

/// <summary>
/// ETKDG-style conformer generation compared with RDKit 2026.09.1 ETKDGv3 (Data/conformer_reference.json):
/// bond lengths and angles against an RDKit conformer, stereochemistry, ring planarity, clashes and coverage of the
/// RDKit ensemble relative to RDKit's own seed-to-seed variability.
/// </summary>
public class ConformerTests
{
    public static TheoryData<string> Names() => ConformerReference.Theory();

    private static readonly Dictionary<(string, bool), ConformerEnsemble> Cache = [];

    private static ConformerEnsemble Ensemble(string name, bool cleanup)
    {
        lock (Cache)
        {
            if (Cache.TryGetValue((name, cleanup), out var e)) return e;
            var (_, r) = ConformerReference.Load(name);
            var mol = Molecule.FromSmiles(r.GetProperty("smiles").GetString()!);
            e = ConformerGenerator.Generate(mol, 10, new ConformerOptions { ForceFieldCleanup = cleanup });
            Cache[(name, cleanup)] = e;
            return e;
        }
    }

    private static double Angle(Point3[] c, int i, int j, int k)
    {
        var u = c[i] - c[j];
        var v = c[k] - c[j];
        return Math.Acos(Math.Clamp(Point3.Dot(u, v) / (u.Length * v.Length), -1, 1)) * 180 / Math.PI;
    }

    [Fact]
    public void Torsion_library_parses_completely()
    {
        Assert.Equal(0, EtkdgTorsionLibrary.Skipped);
        Assert.True(EtkdgTorsionLibrary.Regular.Count >= 360);
        Assert.True(EtkdgTorsionLibrary.Macrocycle.Count >= 360);
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Bond_lengths_and_angles_agree_with_rdkit(string name)
    {
        var ens = Ensemble(name, cleanup: true);
        Assert.True(ens.Count >= 1, $"{name}: no conformer");
        var (_, r) = ConformerReference.Load(name);
        var rd = ConformerReference.Coords(r.GetProperty("xyz"));
        var m = ens.Molecule;
        double sumB = 0, maxB = 0, sumA = 0, maxA = 0;
        int nB = 0, nA = 0;
        foreach (var c in ens.Conformers)
        {
            foreach (var b in m.Bonds)
            {
                double d = Math.Abs(Point3.Distance(c[b.Begin], c[b.End]) - Point3.Distance(rd[b.Begin], rd[b.End]));
                sumB += d; maxB = Math.Max(maxB, d); nB++;
            }
            for (int j = 0; j < m.AtomCount; j++)
            {
                var nb = m.Neighbors(j).ToArray();
                for (int x = 0; x < nb.Length; x++)
                    for (int y = x + 1; y < nb.Length; y++)
                    {
                        double d = Math.Abs(Angle(c, nb[x], j, nb[y]) - Angle(rd, nb[x], j, nb[y]));
                        sumA += d; maxA = Math.Max(maxA, d); nA++;
                    }
            }
        }
        Assert.True(sumB / nB < 0.03 && maxB < 0.15, $"{name}: bond deviation mean {sumB / nB:F3} max {maxB:F3} Å");
        Assert.True(sumA / nA < 4.0 && maxA < 30, $"{name}: angle deviation mean {sumA / nA:F2} max {maxA:F1}°");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Stereocentres_match_the_rdkit_conformer(string name)
    {
        var (_, r) = ConformerReference.Load(name);
        var rd = ConformerReference.Coords(r.GetProperty("xyz"));
        var ens = Ensemble(name, cleanup: true);
        var m = ens.Molecule;
        foreach (var tag in r.GetProperty("chiral_tags").EnumerateArray())
        {
            int c = tag[0].GetInt32();
            var nb = m.Neighbors(c).ToArray().Where(x => m.Atoms[x].AtomicNumber != 1).OrderBy(x => x).Take(3).ToArray();
            double Vol(Point3[] p) => Point3.Dot(p[nb[0]] - p[c], Point3.Cross(p[nb[1]] - p[c], p[nb[2]] - p[c]));
            double reference = Vol(rd);
            foreach (var conf in ens.Conformers)
                Assert.True(Math.Sign(Vol(conf)) == Math.Sign(reference), $"{name}: centre {c} inverted");
        }
    }

    [Fact]
    public void Double_bond_geometry_is_preserved()
    {
        var (_, r) = ConformerReference.Load("tamoxifen_like_EZ");
        var rd = ConformerReference.Coords(r.GetProperty("xyz"));
        var ens = Ensemble("tamoxifen_like_EZ", cleanup: true);
        var m = ens.Molecule;
        var stereo = m.Bonds.ToArray().Where(b => b.Stereo is BondStereo.Cis or BondStereo.Trans).ToList();
        Assert.NotEmpty(stereo);
        foreach (var b in stereo)
        {
            double Dih(Point3[] p)
            {
                var b1 = p[b.Begin] - p[b.StereoAtomBegin];
                var b2 = p[b.End] - p[b.Begin];
                var b3 = p[b.StereoAtomEnd] - p[b.End];
                var n1 = Point3.Cross(b1, b2);
                var n2 = Point3.Cross(b2, b3);
                return Math.Abs(Math.Atan2(Point3.Dot(Point3.Cross(n1, n2), b2) / b2.Length, Point3.Dot(n1, n2))) * 180 / Math.PI;
            }
            bool rdTrans = Dih(rd) > 90;
            Assert.Equal(b.Stereo == BondStereo.Trans, rdTrans);
            foreach (var c in ens.Conformers) Assert.Equal(rdTrans, Dih(c) > 90);
        }
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Aromatic_rings_are_planar_and_there_are_no_clashes(string name)
    {
        var ens = Ensemble(name, cleanup: true);
        var m = ens.Molecule;
        var dist = m.TopologicalDistances();
        foreach (var c in ens.Conformers)
        {
            foreach (var ring in m.Rings.AtomRings.Where(rr => rr.All(a => m.Atoms[a].IsAromatic)))
            {
                double dev = PlaneDeviation(ring.Select(a => c[a]).ToArray());
                // UFF puckers strained fused aromatics slightly (RDKit ETKDG + UFF gives 0.07 Å for morphine too)
                Assert.True(dev < 0.1, $"{name}: aromatic ring not planar ({dev:F3} Å)");
            }
            for (int i = 0; i < m.AtomCount; i++)
                for (int j = i + 1; j < m.AtomCount; j++)
                {
                    if (m.Atoms[i].AtomicNumber == 1 || m.Atoms[j].AtomicNumber == 1 || dist[i, j] is >= 0 and < 4) continue;
                    Assert.True(Point3.Distance(c[i], c[j]) > 2.2, $"{name}: heavy atoms {i} and {j} clash");
                }
        }
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Ensemble_covers_rdkit_conformers_about_as_well_as_rdkit_itself(string name)
    {
        var (mol, r) = ConformerReference.Load(name);
        var heavyMol = Molecule.FromSmiles(r.GetProperty("smiles").GetString()!);
        if (heavyMol.Fragments().Count > 1) return; // fragments are laid side by side here, overlaid by RDKit
        var rdEns = r.GetProperty("ensemble").EnumerateArray().Select(ConformerReference.Coords).ToList();
        double selfMean = r.GetProperty("rdkit_self_coverage")[0].GetDouble();
        foreach (bool cleanup in new[] { false, true })
        {
            var ens = Ensemble(name, cleanup);
            var heavy = Enumerable.Range(0, ens.Molecule.AtomCount).Where(i => ens.Molecule.Atoms[i].AtomicNumber != 1).ToArray();
            double mean = rdEns.Average(rc => ens.Conformers.Min(c => ConformerAlignment.BestRmsd(heavyMol, rc, heavy.Select(i => c[i]).ToArray())));
            double slack = cleanup ? 0.5 : 0.4;
            Assert.True(mean <= selfMean + slack, $"{name} (cleanup={cleanup}): mean best RMSD {mean:F2} Å vs RDKit self-coverage {selfMean:F2} Å");
        }
        _ = mol;
    }

    [Fact]
    public void Generation_is_deterministic_for_a_seed()
    {
        var mol = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O");
        var a = ConformerGenerator.Generate(mol, 4, new ConformerOptions { RandomSeed = 11 });
        var b = ConformerGenerator.Generate(mol, 4, new ConformerOptions { RandomSeed = 11, MaxDegreeOfParallelism = 1 });
        var c = ConformerGenerator.Generate(mol, 4, new ConformerOptions { RandomSeed = 12 });
        Assert.Equal(a.Count, b.Count);
        for (int k = 0; k < a.Count; k++) Assert.Equal(a.Conformers[k], b.Conformers[k]);
        Assert.NotEqual(a.Conformers[0], c.Conformers[0]);
    }

    [Fact]
    public void Pruning_keeps_conformers_apart()
    {
        var mol = Molecule.FromSmiles("CCCCCCOc1ccccc1C(=O)NCCN(C)C");
        var ens = ConformerGenerator.Generate(mol, 20, new ConformerOptions { PruneRmsThreshold = 1.0 });
        Assert.InRange(ens.Count, 2, 20);
        var heavy = Enumerable.Range(0, ens.Molecule.AtomCount).Where(i => ens.Molecule.Atoms[i].AtomicNumber != 1).ToArray();
        for (int i = 0; i < ens.Count; i++)
            for (int j = i + 1; j < ens.Count; j++)
                Assert.True(ConformerAlignment.BestRmsd(mol, heavy.Select(x => ens.Conformers[i][x]).ToArray(), heavy.Select(x => ens.Conformers[j][x]).ToArray()) >= 1.0);
        Assert.All(ens.Energies, e => Assert.False(double.IsNaN(e)));
    }

    [Fact]
    public void Embed_returns_a_3d_molecule_with_hydrogens_that_round_trips_through_sdf()
    {
        var m = ConformerGenerator.Embed(Molecule.FromSmiles("OC(=O)c1ccccc1OC(C)=O", "aspirin"));
        Assert.True(m.Has3DCoordinates);
        Assert.Equal(21, m.AtomCount);
        Assert.Contains("UFF_ENERGY", m.Properties.Keys);
        var sdf = SdfWriter.ToSdf([m]);
        var back = SdfReader.Parse(sdf, new MolReaderOptions { RemoveHydrogens = false }).Single();
        Assert.Equal(m.AtomCount, back.AtomCount);
        for (int i = 0; i < m.AtomCount; i++) Assert.True(Point3.Distance(m.Coordinates![i], back.Coordinates![i]) < 1e-3);
    }

    [Fact]
    public void Salts_are_embedded_fragment_by_fragment()
    {
        var ens = ConformerGenerator.Generate(Molecule.FromSmiles("C[NH3+].CC(=O)[O-]"), 2);
        Assert.True(ens.Count >= 1);
        var c = ens.Conformers[0];
        var frags = ens.Molecule.Fragments();
        Assert.Equal(2, frags.Count);
        double minD = frags[0].SelectMany(i => frags[1].Select(j => Point3.Distance(c[i], c[j]))).Min();
        Assert.True(minD > 2.0, $"fragments overlap ({minD:F2} Å)");
    }

    [Fact]
    public void Rmsd_and_superposition_handle_rigid_motion_and_symmetry()
    {
        var mol = Molecule.FromSmiles("Cc1ccccc1");
        var ens = ConformerGenerator.Generate(mol, 1);
        var c = ens.Conformers[0];
        var rot = c.Select(p => new Point3(p.Y + 3, -p.X - 1, p.Z + 2)).ToArray(); // 90° about z + shift
        Assert.True(ConformerAlignment.Rmsd(c, rot) < 1e-6);
        Assert.True(ConformerAlignment.Rmsd(c, rot, align: false) > 1);
        var back = ConformerAlignment.Superpose(c, rot);
        for (int i = 0; i < c.Length; i++) Assert.True(Point3.Distance(c[i], back[i]) < 1e-6);
        // swapping the two ortho carbons is a symmetry of toluene
        var heavy = Enumerable.Range(0, ens.Molecule.AtomCount).Where(i => ens.Molecule.Atoms[i].AtomicNumber != 1).ToArray();
        var a = heavy.Select(i => c[i]).ToArray();
        var b = (Point3[])a.Clone();
        (b[2], b[6]) = (b[6], b[2]);
        (b[3], b[5]) = (b[5], b[3]);
        Assert.True(ConformerAlignment.Rmsd(a, b, align: false) > 0.5);
        Assert.True(ConformerAlignment.BestRmsd(mol, a, b, align: false) < 1e-6);
    }

    private static double PlaneDeviation(Point3[] pts)
    {
        var c = new Point3(pts.Average(p => p.X), pts.Average(p => p.Y), pts.Average(p => p.Z));
        // normal from the cross products of consecutive vectors (adequate for near-planar rings)
        var n = new Point3(0, 0, 0);
        for (int i = 0; i < pts.Length; i++) n += Point3.Cross(pts[i] - c, pts[(i + 1) % pts.Length] - c);
        n *= 1 / n.Length;
        return pts.Max(p => Math.Abs(Point3.Dot(p - c, n)));
    }
}
