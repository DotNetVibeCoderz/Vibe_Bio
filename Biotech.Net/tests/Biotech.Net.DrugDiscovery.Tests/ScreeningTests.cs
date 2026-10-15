using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

namespace Biotech.Net.DrugDiscovery.Tests;

public class SimilarityScreenTests
{
    private static IReadOnlyList<LibraryCompound> Library() =>
    [
        LibraryCompound.FromSmiles("aspirin", Drugs.Aspirin),
        LibraryCompound.FromSmiles("salicylic", "O=C(O)c1ccccc1O"),
        LibraryCompound.FromSmiles("ibuprofen", Drugs.Ibuprofen),
        LibraryCompound.FromSmiles("naproxen", Drugs.Naproxen),
        LibraryCompound.FromSmiles("caffeine", Drugs.Caffeine),
        LibraryCompound.FromSmiles("morphine", Drugs.Morphine),
        LibraryCompound.FromSmiles("ethanol", "CCO"),
        LibraryCompound.FromSmiles("glucose", "OCC1OC(O)C(O)C(O)C1O"),
    ];

    [Fact]
    public void The_query_itself_ranks_first()
    {
        var hits = SimilarityScreen.Run(Drugs.Parse(Drugs.Aspirin), Library());
        Assert.Equal("aspirin", hits[0].Id);
        Assert.Equal(1, hits[0].Rank);
        Assert.Equal(Library().Count, hits.Count);
        Assert.Equal(1.0, hits[0].PerFingerprint[FingerprintType.Ecfp4], 9);
        // ranks are 1..n and scores are non-increasing
        Assert.Equal(Enumerable.Range(1, hits.Count), hits.Select(h => h.Rank));
        Assert.True(hits.Zip(hits.Skip(1)).All(p => p.First.Score >= p.Second.Score));
    }

    [Fact]
    public void Close_analogues_outrank_unrelated_compounds()
    {
        var hits = SimilarityScreen.Run(Drugs.Parse(Drugs.Aspirin), Library());
        int salicylic = hits.ToList().FindIndex(h => h.Id == "salicylic");
        int glucose = hits.ToList().FindIndex(h => h.Id == "glucose");
        Assert.True(salicylic < glucose);
        Assert.True(salicylic <= 1);
    }

    [Fact]
    public void Group_fusion_over_several_queries_finds_both_series()
    {
        var queries = new[] { Drugs.Parse(Drugs.Ibuprofen), Drugs.Parse(Drugs.Caffeine) };
        var hits = SimilarityScreen.Run(queries, Library(), new SimilarityScreenOptions { QueryFusion = FusionRule.Max });
        var top3 = hits.Take(3).Select(h => h.Id).ToHashSet();
        Assert.Contains("ibuprofen", top3);
        Assert.Contains("caffeine", top3);
    }

    [Fact]
    public void Top_k_limits_the_result()
        => Assert.Equal(3, SimilarityScreen.Run(Drugs.Parse(Drugs.Aspirin), Library(), new SimilarityScreenOptions { TopK = 3 }).Count);

    [Theory]
    [InlineData(FusionRule.Max)]
    [InlineData(FusionRule.Mean)]
    [InlineData(FusionRule.Min)]
    [InlineData(FusionRule.RankSum)]
    [InlineData(FusionRule.ReciprocalRank)]
    public void Every_fusion_rule_keeps_the_exact_match_on_top(FusionRule rule)
    {
        var hits = SimilarityScreen.Run(Drugs.Parse(Drugs.Aspirin), Library(),
            new SimilarityScreenOptions { FingerprintFusion = rule, QueryFusion = rule });
        Assert.Equal("aspirin", hits[0].Id);
    }

    [Fact]
    public void Fuse_computes_the_documented_combinations()
    {
        double[][] lists = [[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]];
        Assert.Equal([1.0, 1.0, 0.5], SimilarityScreen.Fuse(lists, FusionRule.Max));
        Assert.Equal([0.0, 0.0, 0.5], SimilarityScreen.Fuse(lists, FusionRule.Min));
        Assert.Equal([0.5, 0.5, 0.5], SimilarityScreen.Fuse(lists, FusionRule.Mean));
    }

    [Fact]
    public void Ranks_handle_ties_by_averaging()
        => Assert.Equal([1.0, 2.5, 2.5, 4.0], SimilarityScreen.Ranks([0.9, 0.5, 0.5, 0.1]));

    [Fact]
    public void Enrichment_and_roc_auc_describe_a_ranking()
    {
        var library = Library();
        var actives = new HashSet<string> { "salicylic", "ibuprofen", "naproxen", "aspirin" };
        var hits = SimilarityScreen.Run(Drugs.Parse(Drugs.Aspirin), library);
        double auc = SimilarityScreen.RocAuc(hits, actives);
        Assert.InRange(auc, 0.7, 1.0);
        Assert.True(SimilarityScreen.EnrichmentFactor(hits, actives, 0.5) >= 1.0);
        // a perfect ranking scores 1, a reversed one 0
        var perfect = actives.Select((a, i) => new ScreeningHit(a, i + 1, 1, new Dictionary<FingerprintType, double>()))
            .Concat(library.Where(c => !actives.Contains(c.Id)).Select((c, i) => new ScreeningHit(c.Id, actives.Count + i + 1, 0, new Dictionary<FingerprintType, double>()))).ToList();
        Assert.Equal(1.0, SimilarityScreen.RocAuc(perfect, actives));
        Assert.Equal(0.0, SimilarityScreen.RocAuc(perfect.AsEnumerable().Reverse().ToList(), actives));
    }

    [Fact]
    public void An_empty_query_list_is_rejected()
        => Assert.Throws<ArgumentException>(() => SimilarityScreen.Run([], Library()));
}

public class AlertTests
{
    [Fact]
    public void Catalogues_parse_and_are_non_empty()
    {
        Assert.NotEmpty(StructuralAlerts.Pains);
        Assert.NotEmpty(StructuralAlerts.Brenk);
        Assert.Equal(StructuralAlerts.Pains.Count + StructuralAlerts.Brenk.Count, StructuralAlerts.All.Count);
        Assert.All(StructuralAlerts.All, a => Assert.True(a.Pattern.AtomCount > 0, a.Name));
    }

    [Theory]
    [InlineData("O=C1C=CC(=O)C=C1", "quinone_A")]
    [InlineData("Oc1ccccc1O", "catechol_A")]
    [InlineData("c1ccc(cc1)N=Nc1ccccc1", "azo_compound")]
    [InlineData("CCOC(=O)Cl", "acyl_halide")]
    [InlineData("C1OC1C", "epoxide_aziridine")]
    [InlineData("CCSSCC", "disulfide")]
    [InlineData("CC(=O)OO", "peroxide")]
    [InlineData("c1ccccc1[N+](=O)[O-]", "nitro")]
    public void Known_liabilities_fire_the_right_alert(string smiles, string alertName)
    {
        var hits = StructuralAlerts.Match(Drugs.Parse(smiles));
        Assert.Contains(hits, a => a.Name == alertName);
    }

    [Theory]
    [InlineData(Drugs.Aspirin)]
    [InlineData(Drugs.Ibuprofen)]
    [InlineData(Drugs.Caffeine)]
    public void Clean_drugs_trigger_no_alert(string smiles) => Assert.True(StructuralAlerts.IsClean(Drugs.Parse(smiles)));

    [Fact]
    public void Families_can_be_selected_by_name()
    {
        Assert.Same(StructuralAlerts.Pains, StructuralAlerts.Family("pains"));
        Assert.Same(StructuralAlerts.Brenk, StructuralAlerts.Family("BRENK"));
        Assert.Throws<ArgumentException>(() => StructuralAlerts.Family("nope"));
        // family filtering limits the hits
        var quinone = Drugs.Parse("O=C1C=CC(=O)C=C1");
        Assert.All(StructuralAlerts.Match(quinone, StructuralAlerts.Pains), a => Assert.Equal("PAINS", a.Family));
    }
}

public class PharmacophoreTests
{
    [Fact]
    public void Feature_counts_on_aspirin()
    {
        var m = Drugs.Parse(Drugs.Aspirin);
        Assert.Equal(1, Pharmacophore.Count(m, PharmacophoreFeature.Donor));           // the COOH
        Assert.Equal(1, Pharmacophore.Count(m, PharmacophoreFeature.NegativeIonisable));
        Assert.Equal(1, Pharmacophore.Count(m, PharmacophoreFeature.Aromatic));
        Assert.Equal(0, Pharmacophore.Count(m, PharmacophoreFeature.PositiveIonisable));
        Assert.Equal(0, Pharmacophore.Count(m, PharmacophoreFeature.HalogenDonor));
    }

    [Fact]
    public void Basic_amines_and_halogens_are_typed()
    {
        Assert.Equal(1, Pharmacophore.Count(Drugs.Parse(Drugs.Propranolol), PharmacophoreFeature.PositiveIonisable));
        Assert.Equal(1, Pharmacophore.Count(Drugs.Parse("Clc1ccccc1"), PharmacophoreFeature.HalogenDonor));
        Assert.Equal(2, Pharmacophore.Count(Drugs.Parse("c1ccc2ccccc2c1"), PharmacophoreFeature.Aromatic));
    }

    [Fact]
    public void Counts_vector_and_profile_similarity()
    {
        var all = Pharmacophore.Counts(Drugs.Parse(Drugs.Ibuprofen));
        Assert.Equal(Enum.GetValues<PharmacophoreFeature>().Length, all.Count);
        Assert.Equal(all.Values.Select(v => (double)v), Pharmacophore.ToVector(Drugs.Parse(Drugs.Ibuprofen)));
        var a = Drugs.Parse(Drugs.Ibuprofen);
        Assert.Equal(1.0, Pharmacophore.ProfileSimilarity(a, a));
        double related = Pharmacophore.ProfileSimilarity(Drugs.Parse(Drugs.Ibuprofen), Drugs.Parse(Drugs.Naproxen));
        double unrelated = Pharmacophore.ProfileSimilarity(Drugs.Parse(Drugs.Ibuprofen), Drugs.Parse("OCC1OC(O)C(O)C(O)C1O"));
        Assert.InRange(related, 0, 1);
        Assert.True(related > unrelated);
    }

    [Fact]
    public void Definitions_are_exposed_and_parse()
        => Assert.All(Enum.GetValues<PharmacophoreFeature>(), f => Assert.NotNull(SmartsPattern.TryParse(Pharmacophore.Definition(f))));
}

public class CompoundFilterTests
{
    private static IReadOnlyList<LibraryCompound> Library() =>
    [
        LibraryCompound.FromSmiles("aspirin", Drugs.Aspirin),
        LibraryCompound.FromSmiles("ibuprofen", Drugs.Ibuprofen),
        LibraryCompound.FromSmiles("vancomycin", Drugs.Vancomycin),
        LibraryCompound.FromSmiles("quinone", "O=C1C=CC(=O)C=C1"),
        LibraryCompound.FromSmiles("ethanol", "CCO"),
    ];

    [Fact]
    public void Property_windows_reject_outliers()
    {
        var filter = new CompoundFilter().WithPropertyRange("MW", MolecularDescriptors.MolecularWeight, 150, 500);
        var kept = filter.Apply(Library()).Select(c => c.Id).ToList();
        Assert.Contains("aspirin", kept);
        Assert.DoesNotContain("vancomycin", kept);
        Assert.DoesNotContain("ethanol", kept);
        var r = filter.Evaluate(Drugs.Parse("CCO"), "ethanol");
        Assert.False(r.Passed);
        Assert.Contains("MW", r.Reasons[0]);
    }

    [Fact]
    public void Substructure_requirements_and_exclusions()
    {
        var aromatic = new CompoundFilter().RequireSubstructure("c1ccccc1", "benzene");
        Assert.DoesNotContain("ethanol", aromatic.Apply(Library()).Select(c => c.Id));
        var noAcid = new CompoundFilter().ExcludeSubstructure("C(=O)[OH]", "carboxylic acid");
        Assert.DoesNotContain("aspirin", noAcid.Apply(Library()).Select(c => c.Id));
        Assert.Contains("ethanol", noAcid.Apply(Library()).Select(c => c.Id));
    }

    [Fact]
    public void Alerts_can_be_excluded()
    {
        var kept = new CompoundFilter().ExcludeAlerts(StructuralAlerts.Pains).Apply(Library()).Select(c => c.Id).ToList();
        Assert.DoesNotContain("quinone", kept);
        Assert.Contains("aspirin", kept);
    }

    [Fact]
    public void Rules_can_be_required_and_criteria_accumulate()
    {
        var filter = new CompoundFilter()
            .RequireLipinski()
            .WithPropertyRange("logP", Crippen.LogP, max: 4)
            .ExcludeAlerts();
        Assert.Equal(3, filter.Criteria.Count);
        Assert.True(filter.Passes(Drugs.Parse(Drugs.Aspirin)));
        Assert.False(filter.Passes(Drugs.Parse(Drugs.Vancomycin)));
        // every failed criterion is reported, not just the first
        var vanco = filter.Evaluate(Drugs.Parse(Drugs.Vancomycin), "vancomycin");
        Assert.False(vanco.Passed);
        Assert.NotEmpty(vanco.Reasons);
        Assert.Contains("reject", vanco.ToString());
    }

    [Fact]
    public void Preset_filters_behave_as_documented()
    {
        Assert.True(CompoundFilter.FragmentLike().Passes(Drugs.Parse("c1ccc(O)cc1")));
        Assert.False(CompoundFilter.FragmentLike().Passes(Drugs.Parse(Drugs.Imatinib)));
        Assert.True(CompoundFilter.LeadLike().Passes(Drugs.Parse(Drugs.Naproxen)));
        Assert.False(CompoundFilter.LeadLike().Passes(Drugs.Parse(Drugs.Atorvastatin)));
    }

    [Fact]
    public void Evaluate_keeps_the_library_order()
    {
        var results = new CompoundFilter().RequireLipinski().Evaluate(Library());
        Assert.Equal(Library().Select(c => c.Id), results.Select(r => r.Id));
    }
}
