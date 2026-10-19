using System.Text.Json;
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.DrugDiscovery.Tests;

/// <summary>
/// Vina-style scoring and flexible docking checked against AutoDock Vina 1.2.7 (Data/vina_reference.json, produced by
/// Data/docking_reference_gen.py) on two public complexes: trypsin–benzamidine (3PTB) and streptavidin–biotin (1STP).
/// </summary>
public class FlexibleDockingTests
{
    private static readonly string DataDir = Path.Combine(AppContext.BaseDirectory, "Data");
    private static readonly Lazy<JsonDocument> Reference = new(() => JsonDocument.Parse(File.ReadAllText(Path.Combine(DataDir, "vina_reference.json"))));
    private static readonly Dictionary<string, VinaReceptor> Receptors = [];

    internal static VinaReceptor Receptor(string id)
    {
        lock (Receptors)
        {
            if (!Receptors.TryGetValue(id, out var r)) Receptors[id] = r = VinaReceptor.FromPdb(Path.Combine(DataDir, id + "_receptor.pdb"));
            return r;
        }
    }

    internal static Molecule Crystal(string id) => SdfReader.ReadFile(Path.Combine(DataDir, id + "_ligand.sdf")).First();

    private static JsonElement Ref(string id) => Reference.Value.RootElement.GetProperty("complexes").GetProperty(id);

    private static Vec3 Centre(string id)
    {
        var c = Ref(id).GetProperty("center");
        return new Vec3(c[0].GetDouble(), c[1].GetDouble(), c[2].GetDouble());
    }

    /// <summary>
    /// In-place heavy-atom RMSD to the crystal ligand, minimised over graph symmetries ignoring bond orders (so the two
    /// amidine nitrogens of benzamidine, or carboxylate oxygens, are interchangeable as in the crystal).
    /// </summary>
    internal static double RmsdToCrystal(Molecule crystal, Molecule pose)
    {
        var skeleton = crystal.Clone();
        for (int b = 0; b < skeleton.BondCount; b++) skeleton.BondAt(b).Order = BondOrder.Single;
        for (int a = 0; a < skeleton.AtomCount; a++) { skeleton.AtomAt(a).IsAromatic = false; skeleton.AtomAt(a).FormalCharge = 0; }
        skeleton.InvalidateTopology();
        var heavy = Enumerable.Range(0, pose.AtomCount).Where(i => pose.Atoms[i].AtomicNumber != 1).ToArray();
        return ConformerAlignment.BestRmsd(skeleton, crystal.Coordinates!, heavy.Select(i => pose.Coordinates![i]).ToArray(), align: false);
    }

    private static readonly Dictionary<(string, bool), DockingResult> Docked = [];

    private static DockingResult Dock(string id, bool fromCrystalGeometry)
    {
        lock (Docked)
        {
            if (Docked.TryGetValue((id, fromCrystalGeometry), out var r)) return r;
            var crystal = Crystal(id);
            var start = fromCrystalGeometry ? crystal : ConformerGenerator.Embed(crystal, new ConformerOptions { RandomSeed = 3 });
            r = VinaDocking.Dock(start, Receptor(id), new DockingOptions { Center = Centre(id) });
            Docked[(id, fromCrystalGeometry)] = r;
            return r;
        }
    }

    [Theory]
    [InlineData("3ptb", 1.0, 0.3)]
    [InlineData("1stp", 1.5, 0.5)]
    public void Redocking_from_the_crystal_geometry_recovers_the_pose(string id, double maxRmsd, double maxEnergyGap)
    {
        var res = Dock(id, fromCrystalGeometry: true);
        var vina = Ref(id).GetProperty("redock")[0];
        double rmsd = RmsdToCrystal(Crystal(id), res.Best.Molecule);
        Assert.True(rmsd < maxRmsd, $"{id}: top pose {rmsd:F2} Å from the crystal ligand (Vina: {vina.GetProperty("rmsd_to_crystal").GetDouble():F2} Å)");
        Assert.True(Math.Abs(res.Best.Affinity - vina.GetProperty("affinity").GetDouble()) < maxEnergyGap,
            $"{id}: best affinity {res.Best.Affinity:F3} vs Vina {vina.GetProperty("affinity").GetDouble():F3}");
        Assert.True(res.Poses.Zip(res.Poses.Skip(1)).All(p => p.First.Affinity <= p.Second.Affinity + 1e-9));
        Assert.All(res.Poses, p => Assert.True(p.Affinity <= res.Best.Affinity + 3.0 + 1e-9));
    }

    [Fact]
    public void Redocking_biotin_from_a_generated_conformer_finds_the_site()
    {
        var res = Dock("1stp", fromCrystalGeometry: false);
        double best = RmsdToCrystal(Crystal("1stp"), res.Best.Molecule);
        double top3 = res.Poses.Take(3).Min(p => RmsdToCrystal(Crystal("1stp"), p.Molecule));
        Assert.True(top3 < 2.0, $"best of top-3 poses {top3:F2} Å (top pose {best:F2} Å)");
        Assert.True(res.Best.Affinity < -6.0);
    }

    [Fact]
    public void Docking_is_deterministic_for_a_seed()
    {
        var opts = new DockingOptions { Center = Centre("3ptb"), MonteCarloSteps = 200, Exhaustiveness = 3, Seed = 9 };
        var a = VinaDocking.Dock(Crystal("3ptb"), Receptor("3ptb"), opts);
        var b = VinaDocking.Dock(Crystal("3ptb"), Receptor("3ptb"), opts with { MaxDegreeOfParallelism = 1 });
        Assert.Equal(a.Poses.Count, b.Poses.Count);
        for (int k = 0; k < a.Poses.Count; k++) Assert.Equal(a.Poses[k].Molecule.Coordinates, b.Poses[k].Molecule.Coordinates);
    }

    [Fact]
    public void Poses_export_to_sdf_and_pdbqt()
    {
        var res = Dock("1stp", fromCrystalGeometry: true);
        var back = SdfReader.Parse(res.ToSdf());
        Assert.Equal(res.Poses.Count, back.Count);
        Assert.Equal(res.Best.Affinity.ToString("F3", System.Globalization.CultureInfo.InvariantCulture), back[0].Properties["VINA_AFFINITY"]);
        string pdbqt = res.ToPdbqt();
        var lines = pdbqt.Split('\n');
        Assert.Equal(res.Poses.Count, lines.Count(l => l.StartsWith("MODEL")));
        Assert.Equal(lines.Count(l => l.StartsWith("BRANCH")), lines.Count(l => l.StartsWith("ENDBRANCH")));
        Assert.Equal(5 * res.Poses.Count, lines.Count(l => l.StartsWith("BRANCH")));
        Assert.Contains("TORSDOF 6", pdbqt);
        Assert.Contains("REMARK VINA RESULT:", pdbqt);
        // the crystal ligand has implicit hydrogens only, so each model lists its 16 heavy atoms
        Assert.Equal(16 * res.Poses.Count, lines.Count(l => l.StartsWith("ATOM")));
    }

    [Theory]
    [InlineData("3ptb")]
    [InlineData("1stp")]
    public void Crystal_pose_score_matches_vina(string id)
    {
        var score = VinaDocking.Score(Crystal(id), Receptor(id));
        var v = Ref(id).GetProperty("crystal_score");
        Assert.Equal(v.GetProperty("affinity").GetDouble(), score.Affinity, 0.02);
        Assert.Equal(v.GetProperty("intermolecular").GetDouble(), score.Intermolecular, 0.02);
        Assert.Equal(v.GetProperty("torsional").GetDouble(), score.Torsional, 0.02);
        // Vina also counts pairs whose distance is fixed (atoms on a rotor axis), so its internal energy differs slightly
        Assert.Equal(v.GetProperty("intramolecular").GetDouble(), score.Intramolecular, 0.06);
    }

    [Fact]
    public void Ligand_typing_follows_vina_conventions()
    {
        var types = VinaAtomTypes.AssignLigand(Molecule.FromSmiles("NC(=O)c1ccc(O)cc1Cl"));
        Assert.Equal(XsAtomType.N_D, types[0]);      // amide N: donor only
        Assert.Equal(XsAtomType.C_P, types[1]);
        Assert.Equal(XsAtomType.O_A, types[2]);
        Assert.Equal(XsAtomType.C_H, types[4]);
        Assert.Equal(XsAtomType.O_DA, types[7]);     // phenol
        Assert.Equal(XsAtomType.Cl_H, types[10]);
        Assert.Equal(XsAtomType.N_A, VinaAtomTypes.AssignLigand(Molecule.FromSmiles("c1ccncc1"))[3]);
        Assert.Equal(XsAtomType.N_D, VinaAtomTypes.AssignLigand(Molecule.FromSmiles("C[NH3+]"))[1]);
    }

    [Fact]
    public void Receptor_without_hydrogens_uses_residue_templates()
    {
        var s = PdbReader.Read(Path.Combine(DataDir, "3ptb_receptor.pdb"));
        var noH = VinaReceptor.FromAtoms(s.FirstModel.Atoms.Where(a => a.Element != "H").ToList());
        var withH = Receptor("3ptb");
        Assert.False(noH.DonorsFromHydrogens);
        Assert.True(withH.DonorsFromHydrogens);
        Assert.Equal(withH.Count, noH.Count);
        int agree = Enumerable.Range(0, noH.Count).Count(i => noH.Types[i] == withH.Types[i]);
        Assert.True(agree >= 0.97 * noH.Count, $"{agree}/{noH.Count} types agree");
    }

    [Fact]
    public void Conformation_gradient_matches_finite_differences()
    {
        var lig = new FlexibleLigand(Crystal("1stp"));
        var grid = VinaGrid.Build(Receptor("1stp"), Centre("1stp"), new Vec3(22.5, 22.5, 22.5), lig.Scored.Select(i => lig.Types[i]!.Value));
        var x = new Vec3[lig.AtomCount];
        var ag = new Vec3[lig.AtomCount];
        var rng = new Random(5);
        double E(LigandConf c, Span<double> g)
        {
            lig.Build(c, x);
            Array.Clear(ag);
            double e = 0;
            foreach (int i in lig.Scored) { e += grid.Evaluate(lig.Types[i]!.Value, x[i], 1000, out var gi); ag[i] += gi; }
            e += VinaDocking.Intra(lig, x, VinaScoreOptions.Default, 1000, ag);
            if (!g.IsEmpty) lig.ConfGradient(c, x, ag, g);
            return e;
        }
        for (int trial = 0; trial < 5; trial++)
        {
            var c = lig.InputConf();
            c.Position += new Vec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
            for (int t = 0; t < c.Torsions.Length; t++) c.Torsions[t] = (rng.NextDouble() - 0.5) * 0.6;
            int n = 6 + c.Torsions.Length;
            var g = new double[n];
            E(c, g);
            for (int k = 0; k < n; k++)
            {
                const double h = 1e-5;
                var step = new double[n];
                step[k] = 1;
                var plus = c.Clone(); plus.Increment(step, h);
                var minus = c.Clone(); minus.Increment(step, -h);
                double fd = (E(plus, []) - E(minus, [])) / (2 * h);
                Assert.True(Math.Abs(fd - g[k]) < 1e-3 * Math.Max(1, Math.Abs(fd)), $"dof {k}: analytic {g[k]:F5} numeric {fd:F5}");
            }
        }
    }

    [Fact]
    public void Torsion_tree_counts_match_meeko_and_vina()
    {
        var lig = new FlexibleLigand(Crystal("1stp"));
        Assert.Equal(Ref("1stp").GetProperty("torsdof").GetInt32(), lig.RotatableBonds.Length);  // Meeko TORSDOF 6
        Assert.Equal(5, lig.Torsions.Length);   // the carboxylic O–H rotation moves only a hydrogen
        Assert.Equal(5.5, lig.NumRotors);       // Vina's num_tors (half credit for the OH rotor)
        Assert.Equal(1.0, new FlexibleLigand(Crystal("3ptb")).NumRotors);
    }
}
