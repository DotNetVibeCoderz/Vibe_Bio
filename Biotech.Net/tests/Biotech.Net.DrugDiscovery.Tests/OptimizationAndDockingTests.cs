using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.DrugDiscovery.Tests;

public class BioisostereTests
{
    [Fact]
    public void Catalogue_entries_all_compile()
        => Assert.All(Bioisosteres.Default, t =>
        {
            Assert.True(t.Pattern.AtomCount > 0, t.Name);
            Assert.True(t.Replacement.AtomCount > 0, t.Name);
        });

    [Fact]
    public void Acid_to_tetrazole_produces_the_expected_analogue()
    {
        var transform = Bioisosteres.Default.First(t => t.Name.Contains("tetrazole"));
        var products = LeadOptimizer.Apply(Drugs.Parse(Drugs.Ibuprofen), transform);
        var smiles = products.Select(p => p.ToSmiles(includeStereo: false)).ToList();
        Assert.Contains(Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)c1nnn[nH]1").ToSmiles(includeStereo: false), smiles);
    }

    [Fact]
    public void Ring_transforms_keep_substituents_attached()
    {
        var aza = Bioisosteres.Default.First(t => t.Name.Contains("pyridine"));
        foreach (var p in LeadOptimizer.Apply(Drugs.Parse(Drugs.Ibuprofen), aza))
        {
            Assert.Single(p.Fragments());
            Assert.Equal(1, p.Atoms.ToArray().Count(a => a.AtomicNumber == 7));
            Assert.Equal(MolecularDescriptors.HeavyAtomCount(Drugs.Parse(Drugs.Ibuprofen)), MolecularDescriptors.HeavyAtomCount(p));
        }
    }

    [Fact]
    public void Transformations_that_do_not_apply_produce_nothing()
    {
        var acidTransform = Bioisosteres.Default.First(t => t.Name.Contains("tetrazole"));
        Assert.False(acidTransform.Applies(Drugs.Parse(Drugs.Caffeine)));
        Assert.Empty(LeadOptimizer.Apply(Drugs.Parse(Drugs.Caffeine), acidTransform));
    }

    [Fact]
    public void Applicable_lists_only_matching_transforms()
    {
        var applicable = Bioisosteres.Applicable(Drugs.Parse(Drugs.Ibuprofen));
        Assert.Contains(applicable, t => t.Name.Contains("tetrazole"));
        Assert.All(applicable, t => Assert.True(t.Applies(Drugs.Parse(Drugs.Ibuprofen))));
    }
}

public class LeadOptimizerTests
{
    [Fact]
    public void Enumerates_distinct_connected_analogues()
    {
        var lead = Drugs.Parse(Drugs.Ibuprofen, "ibuprofen");
        var analogues = LeadOptimizer.Enumerate(lead);
        Assert.NotEmpty(analogues);
        Assert.Equal(analogues.Count, analogues.Select(a => a.Smiles).Distinct().Count());
        Assert.DoesNotContain(lead.ToSmiles(includeStereo: false), analogues.Select(a => a.Smiles));
        Assert.All(analogues, a =>
        {
            Assert.Single(a.Molecule.Fragments());
            Assert.Equal(lead.ToSmiles(includeStereo: false), a.ParentSmiles);
            Assert.Equal(a.Smiles, a.Molecule.ToSmiles(includeStereo: false));
        });
    }

    [Fact]
    public void Depth_two_reaches_further_and_respects_the_cap()
    {
        var lead = Drugs.Parse(Drugs.Paracetamol, "paracetamol");
        int one = LeadOptimizer.Enumerate(lead, depth: 1).Count;
        int two = LeadOptimizer.Enumerate(lead, depth: 2).Count;
        Assert.True(two >= one);
        Assert.Equal(3, LeadOptimizer.Enumerate(lead, depth: 2, maxAnalogues: 3).Count);
        Assert.Throws<ArgumentOutOfRangeException>(() => LeadOptimizer.Enumerate(lead, depth: 0));
    }

    [Fact]
    public void Optimize_ranks_the_lead_together_with_its_analogues()
    {
        var lead = Drugs.Parse(Drugs.Ibuprofen, "ibuprofen");
        var ranked = LeadOptimizer.Optimize(lead);
        Assert.Contains(ranked, r => r.Analogue is null); // the lead is included
        Assert.True(ranked.Zip(ranked.Skip(1)).All(p => p.First.Score >= p.Second.Score));
        Assert.All(ranked, r => Assert.InRange(r.Score, 0, 1));
    }
}

public class MpoScoreTests
{
    [Fact]
    public void Desirability_is_one_inside_the_window_and_zero_beyond_the_limit()
    {
        Assert.Equal(1.0, MpoScore.Desirability(300, 400, 550));
        Assert.Equal(0.0, MpoScore.Desirability(600, 400, 550));
        Assert.Equal(0.5, MpoScore.Desirability(475, 400, 550), 9);
        // an inverted window (higher is better)
        Assert.Equal(1.0, MpoScore.Desirability(5, 3, 1));
        Assert.Equal(0.0, MpoScore.Desirability(0, 3, 1));
    }

    [Fact]
    public void Score_is_bounded_and_prefers_drug_like_molecules()
    {
        double ibuprofen = MpoScore.Compute(Drugs.Parse(Drugs.Ibuprofen));
        double vanco = MpoScore.Compute(Drugs.Parse(Drugs.Vancomycin));
        Assert.InRange(ibuprofen, 0, 1);
        Assert.InRange(vanco, 0, 1);
        Assert.True(ibuprofen > vanco);
    }

    [Fact]
    public void Alerts_lower_the_score()
    {
        var quinone = Drugs.Parse("O=C1C=CC(=O)C=C1");
        double withPenalty = MpoScore.Compute(quinone);
        double withoutPenalty = MpoScore.Compute(quinone, MpoProfile.Oral with { AlertPenalty = 0 });
        Assert.True(withoutPenalty > withPenalty);
    }

    [Fact]
    public void Cns_profile_is_stricter_about_polarity()
    {
        var polar = Drugs.Parse("OC(=O)CC(O)(CC(=O)O)C(=O)O"); // citric acid
        Assert.True(MpoScore.Compute(polar, MpoProfile.Cns) <= MpoScore.Compute(polar, MpoProfile.Oral));
        var components = MpoScore.Components(Drugs.Parse(Drugs.Diazepam), MpoProfile.Cns);
        Assert.Equal(6, components.Count);
        Assert.All(components.Values, v => Assert.InRange(v, 0, 1));
    }
}

public class TargetPredictionTests
{
    [Fact]
    public void Demo_reference_set_parses_and_covers_several_targets()
    {
        Assert.All(DemoReferenceSet.Ligands, l => Assert.True(Molecule.FromSmiles(l.Smiles).AtomCount > 0, l.LigandId));
        Assert.True(DemoReferenceSet.Targets.Count >= 10);
        Assert.Contains("COX", DemoReferenceSet.Targets);
    }

    [Theory]
    [InlineData(Drugs.Aspirin, "COX")]
    [InlineData(Drugs.Naproxen, "COX")]
    [InlineData(Drugs.Caffeine, "AdenosineReceptor")]
    [InlineData(Drugs.Propranolol, "BetaAdrenergic")]
    [InlineData(Drugs.Diazepam, "GabaA")]
    [InlineData(Drugs.Morphine, "OpioidReceptor")]
    [InlineData(Drugs.Imatinib, "TyrosineKinase")]
    [InlineData("COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", "DHFR")]
    public void Known_ligands_recover_their_own_target(string smiles, string target)
    {
        var predictions = TargetPredictor.Demo().Predict(Drugs.Parse(smiles));
        Assert.Equal(target, predictions[0].Target);
        Assert.True(predictions[0].ZScore > 3);
        Assert.True(predictions[0].MaxSimilarity > 0.5);
    }

    [Fact]
    public void An_unrelated_molecule_has_no_significant_target()
    {
        var predictions = TargetPredictor.Demo().PredictSignificant(Drugs.Parse("CCO"));
        Assert.Empty(predictions);
    }

    [Fact]
    public void A_close_analogue_is_matched_even_when_it_is_not_in_the_set()
    {
        // ketorolac-like NSAID: not a reference ligand, but close to the COX series
        var predictions = TargetPredictor.Demo().Predict(Drugs.Parse("CC(C(=O)O)c1ccc(cc1)C(C)C"));
        Assert.Equal("COX", predictions[0].Target);
        Assert.True(predictions[0].MaxSimilarity < 1.0);
    }

    [Fact]
    public void Custom_reference_sets_and_metadata()
    {
        var predictor = TargetPredictor.FromReferences(
        [
            new ReferenceLigand("T1", "a", "c1ccccc1O"),
            new ReferenceLigand("T1", "b", "c1ccccc1OC"),
            new ReferenceLigand("T2", "c", "CCCCCCCC"),
        ]);
        Assert.Equal(["T1", "T2"], predictor.Targets.Order().ToArray());
        Assert.Equal(["a", "b"], predictor.LigandsOf("T1"));
        Assert.Empty(predictor.LigandsOf("nope"));
        var top = predictor.Predict(Drugs.Parse("c1ccccc1O"))[0];
        Assert.Equal("T1", top.Target);
        Assert.Equal(2, top.ReferenceCount);
        Assert.Equal("a", top.NearestLigand);
    }

    [Fact]
    public void Extreme_value_p_value_is_monotonic_and_bounded()
    {
        Assert.InRange(TargetPredictor.ExtremeValuePValue(0), 0, 1);
        Assert.True(TargetPredictor.ExtremeValuePValue(10) < TargetPredictor.ExtremeValuePValue(1));
        Assert.True(TargetPredictor.ExtremeValuePValue(-5) > 0.9);
    }
}

public class DockingTests
{
    private static Structure Receptor() => SampleStructures.Crambin();

    private static LigandPose SmallPose(Vec3 centre) =>
        new("probe",
        [
            new PoseAtom(centre, "C", true, false, false),
            new PoseAtom(centre + new Vec3(1.5, 0, 0), "C", true, false, false),
            new PoseAtom(centre + new Vec3(0, 1.4, 0), "O", false, true, true),
        ]);

    [Fact]
    public void Term_shapes_follow_the_published_definitions()
    {
        Assert.Equal(1.0, VinaScore.Hydrophobic(0.0));
        Assert.Equal(0.0, VinaScore.Hydrophobic(2.0));
        Assert.Equal(0.5, VinaScore.Hydrophobic(1.0), 9);
        Assert.Equal(1.0, VinaScore.HydrogenBond(-1.0));
        Assert.Equal(0.0, VinaScore.HydrogenBond(0.5));
        Assert.Equal(0.5, VinaScore.HydrogenBond(-0.35), 9);
    }

    /// <summary>A probe placed in the best pocket and locally optimised, i.e. a favourable (negative-score) pose.</summary>
    private static (LigandPose Pose, IReadOnlyList<Biotech.Net.StructuralBiology.Atom> Atoms) FavourablePose()
    {
        var receptor = Receptor();
        var atoms = receptor.FirstModel.Atoms.Where(a => a.Element != "H").ToList();
        var pocket = PocketDetector.Detect(receptor, new PocketOptions { GridSpacing = 1.2, MinPocketPoints = 5 })[0];
        var (pose, _) = VinaScore.Optimize(SmallPose(pocket.Center), atoms, iterations: 200, seed: 7);
        return (pose, atoms);
    }

    [Fact]
    public void A_pocket_bound_pose_scores_better_than_a_distant_one()
    {
        var (pose, atoms) = FavourablePose();
        var bound = VinaScore.Score(pose, atoms);
        var faraway = VinaScore.Score(pose.Translated(new Vec3(500, 0, 0)), atoms);
        Assert.Equal(0, faraway.ContactCount);
        Assert.Equal(0, faraway.Affinity);
        Assert.True(bound.ContactCount > 0);
        Assert.True(bound.Affinity < 0, $"bound affinity {bound.Affinity}");
        Assert.True(bound.Affinity < faraway.Affinity);
    }

    [Fact]
    public void The_torsional_term_penalises_flexible_ligands()
    {
        var (pose, atoms) = FavourablePose();
        var rigid = VinaScore.Score(pose, atoms, rotatableBonds: 0);
        var flexible = VinaScore.Score(pose, atoms, rotatableBonds: 10);
        Assert.Equal(rigid.Intermolecular, flexible.Intermolecular, 9);
        Assert.True(rigid.Affinity < 0);
        Assert.True(flexible.Affinity > rigid.Affinity); // a favourable score is diluted by the entropy term
        Assert.Equal(10, flexible.RotatableBonds);
    }

    [Fact]
    public void Clashing_poses_are_punished_by_the_repulsion_term()
    {
        var atoms = Receptor().FirstModel.Atoms.Where(a => a.Element != "H").ToList();
        var onTopOfAnAtom = new LigandPose("clash", [new PoseAtom(atoms[0].Position, "C", true, false, false)]);
        var score = VinaScore.Score(onTopOfAnAtom, atoms);
        Assert.True(score.Repulsion > 0);
        Assert.True(score.Affinity > 0); // a clash makes the score unfavourable
    }

    [Fact]
    public void Optimize_never_worsens_the_starting_pose_and_is_reproducible()
    {
        var atoms = Receptor().FirstModel.Atoms.Where(a => a.Element != "H").ToList();
        var pocket = PocketDetector.Detect(Receptor(), new PocketOptions { GridSpacing = 1.2, MinPocketPoints = 5 })[0];
        var start = SmallPose(pocket.Center);
        var before = VinaScore.Score(start, atoms);
        var (pose, after) = VinaScore.Optimize(start, atoms, iterations: 120, seed: 99);
        Assert.True(after.Affinity <= before.Affinity + 1e-9);
        Assert.Equal(start.Atoms.Count, pose.Atoms.Count);
        var (_, again) = VinaScore.Optimize(start, atoms, iterations: 120, seed: 99);
        Assert.Equal(after.Affinity, again.Affinity, 12);
    }

    [Fact]
    public void Poses_can_be_built_from_a_molecule_and_ranked()
    {
        var mol = Molecule.FromSmiles("CCO", "ethanol");
        Assert.Throws<InvalidOperationException>(() => LigandPose.FromMolecule(mol));
        Depiction.Compute2DCoordinates(mol);
        var pose = LigandPose.FromMolecule(mol);
        Assert.Equal(3, pose.Atoms.Count);
        Assert.Contains(pose.Atoms, a => a is { Element: "O", IsAcceptor: true, IsDonor: true });
        var (bound, atoms) = FavourablePose();
        var ranked = VinaScore.Rank([pose.Translated(bound.Center - pose.Center), pose.Translated(new Vec3(500, 0, 0) - pose.Center)], atoms);
        Assert.Equal(2, ranked.Count);
        Assert.True(ranked[0].Score.Affinity <= ranked[1].Score.Affinity);
    }

    [Fact]
    public void Pose_transformations_preserve_geometry()
    {
        var pose = SmallPose(new Vec3(1, 2, 3));
        var moved = pose.Translated(new Vec3(10, 0, 0));
        Assert.Equal(pose.Center + new Vec3(10, 0, 0), moved.Center);
        var rotated = pose.Rotated(Matrix3.FromAxisAngle(Vec3.UnitZ, 90), pose.Center);
        Assert.Equal(
            Vec3.Distance(pose.Atoms[0].Position, pose.Atoms[1].Position),
            Vec3.Distance(rotated.Atoms[0].Position, rotated.Atoms[1].Position), 9);
    }

    [Fact]
    public void Scoring_in_a_detected_pocket_improves_the_placement()
    {
        var receptor = Receptor();
        var pockets = PocketDetector.Detect(receptor, new PocketOptions { GridSpacing = 1.2, MinPocketPoints = 5 });
        Assert.NotEmpty(pockets);
        var pose = SmallPose(new Vec3(0, 0, 0));
        var (refined, score) = VinaScore.ScoreInPocket(pose, receptor, pockets[0], iterations: 80, seed: 7);
        Assert.True(score.ContactCount > 0);
        Assert.True(Vec3.Distance(refined.Center, pockets[0].Center) < 10);
    }
}
