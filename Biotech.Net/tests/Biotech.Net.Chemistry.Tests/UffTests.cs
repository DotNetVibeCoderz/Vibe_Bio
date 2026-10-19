using System.Text.Json;
using Biotech.Net.Chemistry;

namespace Biotech.Net.Chemistry.Tests;

/// <summary>Loads Data/conformer_reference.json (RDKit 2026.09.1, produced by Data/conformer_reference_gen.py).</summary>
internal static class ConformerReference
{
    private static readonly Lazy<JsonDocument> Doc = new(() =>
        JsonDocument.Parse(File.ReadAllText(Path.Combine(AppContext.BaseDirectory, "Data", "conformer_reference.json"))));

    public static IEnumerable<string> Names => Doc.Value.RootElement.GetProperty("mols").EnumerateObject().Select(p => p.Name);

    public static TheoryData<string> Theory()
    {
        var d = new TheoryData<string>();
        foreach (var n in Names) d.Add(n);
        return d;
    }

    /// <summary>The molecule with explicit hydrogens (same atom order as RDKit's AddHs) and the reference record.</summary>
    public static (Molecule Mol, JsonElement Ref) Load(string name)
    {
        var r = Doc.Value.RootElement.GetProperty("mols").GetProperty(name);
        var mol = Molecule.FromSmiles(r.GetProperty("smiles").GetString()!).AddHydrogens();
        return (mol, r);
    }

    public static Point3[] Coords(JsonElement arr)
        => arr.EnumerateArray().Select(p => new Point3(p[0].GetDouble(), p[1].GetDouble(), p[2].GetDouble())).ToArray();
}

public class UffTests
{
    public static TheoryData<string> Names() => ConformerReference.Theory();

    [Theory]
    [MemberData(nameof(Names))]
    public void Explicit_hydrogen_graph_matches_rdkit(string name)
    {
        var (m, r) = ConformerReference.Load(name);
        var atoms = r.GetProperty("atoms").EnumerateArray().Select(a => a.GetInt32()).ToArray();
        Assert.Equal(atoms, m.Atoms.ToArray().Select(a => a.AtomicNumber).ToArray());
        foreach (var b in r.GetProperty("bonds").EnumerateArray())
        {
            int bi = m.BondBetween(b[0].GetInt32(), b[1].GetInt32());
            Assert.True(bi >= 0, $"{name}: missing bond {b[0]}-{b[1]}");
            Assert.Equal(b[2].GetDouble(), UffForceField.BondOrderValue(m.Bonds[bi]));
        }
        Assert.Equal(r.GetProperty("bonds").GetArrayLength(), m.BondCount);
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Hybridization_and_conjugation_match_rdkit(string name)
    {
        var (m, r) = ConformerReference.Load(name);
        var conj = AtomHybridization.ConjugatedBonds(m);
        var hyb = AtomHybridization.Compute(m, conj);
        var expected = r.GetProperty("hybridization").EnumerateArray().Select(h => h.GetString()).ToArray();
        for (int i = 0; i < m.AtomCount; i++)
        {
            if (m.Atoms[i].AtomicNumber == 1 && expected[i] == "UNSPECIFIED") continue; // RDKit's AddHs leaves new H unset
            Assert.True(expected[i] == hyb[i].ToString().ToUpperInvariant(), $"{name} atom {i}: RDKit {expected[i]}, ours {hyb[i]}");
        }
        var rc = r.GetProperty("conjugated").EnumerateArray().Select(c => c.GetBoolean()).ToArray();
        var rb = r.GetProperty("bonds").EnumerateArray().ToArray();
        for (int k = 0; k < rb.Length; k++)
        {
            int bi = m.BondBetween(rb[k][0].GetInt32(), rb[k][1].GetInt32());
            Assert.True(rc[k] == conj[bi], $"{name} bond {rb[k][0]}-{rb[k][1]}: RDKit conjugated={rc[k]}");
        }
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Bond_and_vdw_parameters_match_rdkit(string name)
    {
        var (m, r) = ConformerReference.Load(name);
        var types = UffForceField.AssignAtomTypes(m);
        Assert.All(types, t => Assert.NotNull(t));
        var vdw = r.GetProperty("uff_vdw").EnumerateArray().ToArray();
        for (int i = 0; i < m.AtomCount; i++)
        {
            Assert.True(UffParameters.TryGet(types[i]!, out var p));
            Assert.True(Math.Abs(vdw[i][0].GetDouble() - p.X1) < 1e-6 && Math.Abs(vdw[i][1].GetDouble() - p.D1) < 1e-6,
                $"{name} atom {i} type {types[i]}: vdW x/D differ from RDKit");
        }
        var bonds = r.GetProperty("bonds").EnumerateArray().ToArray();
        var bp = r.GetProperty("uff_bonds").EnumerateArray().ToArray();
        for (int k = 0; k < bonds.Length; k++)
        {
            int a = bonds[k][0].GetInt32(), b = bonds[k][1].GetInt32();
            UffParameters.TryGet(types[a]!, out var pa);
            UffParameters.TryGet(types[b]!, out var pb);
            double r0 = UffForceField.RestLength(bonds[k][2].GetDouble(), pa, pb);
            Assert.True(Math.Abs(bp[k][1].GetDouble() - r0) < 1e-6, $"{name} bond {a}-{b} ({types[a]}-{types[b]}): r0 {r0} vs RDKit {bp[k][1]}");
        }
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Energy_and_gradient_match_rdkit(string name)
    {
        var (m, r) = ConformerReference.Load(name);
        var xyz = ConformerReference.Coords(r.GetProperty("xyz"));
        var ff = UffForceField.Create(m);
        var grad = new Point3[m.AtomCount];
        double e = ff.Energy(xyz, grad);
        double expected = r.GetProperty("uff_energy").GetDouble();
        // coordinates are stored with 4 decimals, so allow a small absolute + relative slack
        Assert.True(Math.Abs(e - expected) < 1e-4 + 1e-6 * Math.Abs(expected), $"{name}: energy {e:F4} vs RDKit {expected:F4}");
        var rg = r.GetProperty("uff_grad").EnumerateArray().Select(v => v.GetDouble()).ToArray();
        double maxDiff = 0, maxRef = 0;
        for (int i = 0; i < m.AtomCount; i++)
        {
            maxDiff = Math.Max(maxDiff, Math.Abs(grad[i].X - rg[3 * i]));
            maxDiff = Math.Max(maxDiff, Math.Abs(grad[i].Y - rg[3 * i + 1]));
            maxDiff = Math.Max(maxDiff, Math.Abs(grad[i].Z - rg[3 * i + 2]));
            maxRef = Math.Max(maxRef, Math.Abs(rg[3 * i]));
        }
        Assert.True(maxDiff < 2e-3 + 1e-4 * maxRef, $"{name}: max gradient deviation {maxDiff:F4} (max |g| {maxRef:F2})");
    }

    [Theory]
    [MemberData(nameof(Names))]
    public void Rdkit_minimum_is_a_minimum_of_our_field(string name)
    {
        var (m, r) = ConformerReference.Load(name);
        var xyz = ConformerReference.Coords(r.GetProperty("uff_min_xyz"));
        var ff = UffForceField.Create(m);
        double expected = r.GetProperty("uff_min_energy").GetDouble();
        double atMin = ff.Energy(xyz);
        Assert.True(Math.Abs(atMin - expected) < 1e-3 + 1e-5 * Math.Abs(expected), $"{name}: energy at RDKit minimum {atMin:F5} vs {expected:F5}");
        var res = ff.Minimize(xyz, 10000, 1e-6);
        Assert.True(Math.Abs(res.Energy - expected) < 0.01, $"{name}: re-minimised {res.Energy:F5} vs RDKit {expected:F5}");
    }

    [Fact]
    public void Minimisation_from_the_same_start_reaches_rdkit_energies()
    {
        int close = 0, total = 0;
        var report = new List<string>();
        foreach (var name in ConformerReference.Names)
        {
            var (m, r) = ConformerReference.Load(name);
            var xyz = ConformerReference.Coords(r.GetProperty("xyz"));
            var res = UffForceField.Create(m).Minimize(xyz, 10000, 1e-6);
            double expected = r.GetProperty("uff_min_energy").GetDouble();
            total++;
            if (Math.Abs(res.Energy - expected) < 0.05) close++;
            else report.Add($"{name}: {res.Energy:F3} vs {expected:F3}");
            Assert.True(res.Energy < expected + 3.0, $"{name}: {res.Energy:F3} far above RDKit {expected:F3}");
        }
        // different optimisers may settle in different nearby minima; most must coincide
        Assert.True(close >= total * 0.8, $"only {close}/{total} reach RDKit's minimum: " + string.Join("; ", report));
    }

    [Fact]
    public void Analytic_gradient_matches_finite_differences()
    {
        foreach (var name in new[] { "atorvastatin", "nitrobenzene", "dimethylsulfone_phosphate", "phenylacetylene_nitrile", "cyclopropyl_cyclobutene", "biotin" })
        {
            var (m, r) = ConformerReference.Load(name);
            var xyz = ConformerReference.Coords(r.GetProperty("xyz"));
            var rnd = new Random(1);
            for (int i = 0; i < xyz.Length; i++) xyz[i] += new Point3(rnd.NextDouble() * 0.2 - 0.1, rnd.NextDouble() * 0.2 - 0.1, rnd.NextDouble() * 0.2 - 0.1);
            var ff = UffForceField.Create(m);
            var grad = new Point3[m.AtomCount];
            ff.Energy(xyz, grad);
            const double h = 1e-5;
            for (int i = 0; i < xyz.Length; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    var plus = (Point3[])xyz.Clone();
                    var minus = (Point3[])xyz.Clone();
                    var d = k == 0 ? new Point3(h, 0, 0) : k == 1 ? new Point3(0, h, 0) : new Point3(0, 0, h);
                    plus[i] += d; minus[i] -= d;
                    double fd = (ff.Energy(plus) - ff.Energy(minus)) / (2 * h);
                    double an = k == 0 ? grad[i].X : k == 1 ? grad[i].Y : grad[i].Z;
                    Assert.True(Math.Abs(fd - an) < 1e-3 * Math.Max(1, Math.Abs(fd)), $"{name} atom {i} axis {k}: analytic {an} vs numeric {fd}");
                }
            }
        }
    }

    [Fact]
    public void Optimize_lowers_the_energy_of_a_distorted_molecule()
    {
        var (m, r) = ConformerReference.Load("ibuprofen");
        var xyz = ConformerReference.Coords(r.GetProperty("xyz"));
        var rnd = new Random(3);
        for (int i = 0; i < xyz.Length; i++) xyz[i] += new Point3(rnd.NextDouble() * 0.4 - 0.2, rnd.NextDouble() * 0.4 - 0.2, rnd.NextDouble() * 0.4 - 0.2);
        m.Coordinates = xyz;
        double before = UffForceField.Create(m).Energy(xyz);
        var res = UffForceField.Optimize(m);
        Assert.True(res.Energy < before);
        Assert.True(res.Converged);
        Assert.Equal(res.Energy, UffForceField.Create(m).Energy(m.Coordinates!), 6);
    }
}
