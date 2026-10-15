using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

namespace Biotech.Net.DrugDiscovery.Tests;

/// <summary>Molecules used across the DrugDiscovery tests.</summary>
public static class Drugs
{
    public const string Aspirin = "CC(=O)Oc1ccccc1C(=O)O";
    public const string Caffeine = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
    public const string Ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
    public const string Paracetamol = "CC(=O)Nc1ccc(O)cc1";
    public const string Atorvastatin = "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O";
    public const string Imatinib = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1";
    public const string Celecoxib = "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1";
    public const string Diazepam = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21";
    public const string Morphine = "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5";
    public const string Propranolol = "CC(C)NCC(O)COc1cccc2ccccc12";
    public const string Naproxen = "COc1ccc2cc(ccc2c1)C(C)C(=O)O";
    public const string Vancomycin = "CNC(CC(C)C)C(=O)NC1C(O)c2ccc(Oc3cc4cc(Oc5ccc(cc5Cl)C(O)C5NC(=O)C(NC(=O)C4NC(=O)C(CC(N)=O)NC1=O)c1ccc(O)c(c1)-c1c(O)cc(O)cc1C(NC5=O)C(=O)O)c3O)c(Cl)c2";

    public static Molecule Parse(string smiles, string name = "") => Molecule.FromSmiles(smiles, name);
}

public class EsolTests
{
    /// <summary>
    /// Delaney's ESOL equation evaluated with Crippen logP; the expected values were computed from the published
    /// coefficients with RDKit descriptors, so these pin the implementation of the equation itself.
    /// </summary>
    [Theory]
    [InlineData(Drugs.Aspirin, -1.992)]
    [InlineData(Drugs.Ibuprofen, -3.087)]
    [InlineData("c1ccccc1", -2.127)]
    [InlineData("CCO", -0.125)]
    [InlineData("c1ccc2ccccc2c1", -3.164)]
    [InlineData(Drugs.Caffeine, -0.871)]
    public void Esol_matches_the_published_equation(string smiles, double expected)
        => Assert.Equal(expected, Admet.EsolLogS(Drugs.Parse(smiles)), 3);

    [Fact]
    public void Esol_orders_molecules_by_known_solubility()
    {
        // ethanol is far more soluble than naphthalene, which beats cholesterol
        double ethanol = Admet.EsolLogS(Drugs.Parse("CCO"));
        double naphthalene = Admet.EsolLogS(Drugs.Parse("c1ccc2ccccc2c1"));
        double cholesterol = Admet.EsolLogS(Drugs.Parse("CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C"));
        Assert.True(ethanol > naphthalene);
        Assert.True(naphthalene > cholesterol);
    }

    [Fact]
    public void Aromatic_proportion_and_solubility_class()
    {
        Assert.Equal(1.0, Admet.AromaticProportion(Drugs.Parse("c1ccccc1")));
        Assert.Equal(0.0, Admet.AromaticProportion(Drugs.Parse("CCO")));
        Assert.Equal(6.0 / 9, Admet.AromaticProportion(Drugs.Parse("Cc1ccccc1CC")), 6);
        Assert.Equal("very soluble", Admet.SolubilityClass(-1));
        Assert.Equal("poorly soluble", Admet.SolubilityClass(-7));
    }

    [Fact]
    public void Solubility_in_mass_units_is_consistent_with_log_s()
    {
        var m = Drugs.Parse(Drugs.Aspirin);
        double expected = Math.Pow(10, Admet.EsolLogS(m)) * MolecularDescriptors.MolecularWeight(m) * 1000;
        Assert.Equal(expected, Admet.EsolSolubilityMgPerL(m), 6);
    }
}

public class RuleTests
{
    [Theory]
    [InlineData(Drugs.Aspirin, true)]
    [InlineData(Drugs.Caffeine, true)]
    [InlineData(Drugs.Ibuprofen, true)]
    [InlineData(Drugs.Paracetamol, true)]
    [InlineData(Drugs.Imatinib, true)]
    [InlineData(Drugs.Diazepam, true)]
    [InlineData(Drugs.Morphine, true)]
    [InlineData(Drugs.Vancomycin, false)] // the textbook rule-of-five exception
    public void Lipinski_on_known_drugs(string smiles, bool passes)
        => Assert.Equal(passes, DrugLikeness.Lipinski(MolecularDescriptors.Compute(Drugs.Parse(smiles))).Passed);

    [Fact]
    public void Lipinski_reports_which_criteria_failed()
    {
        var r = DrugLikeness.Lipinski(MolecularDescriptors.Compute(Drugs.Parse(Drugs.Vancomycin)), allowedViolations: 0);
        Assert.False(r.Passed);
        Assert.Contains(r.Violations, v => v.StartsWith("MW"));
        Assert.Contains(r.Violations, v => v.StartsWith("HBD"));
    }

    [Fact]
    public void Veber_egan_ghose_and_leadlike()
    {
        var ibuprofen = MolecularDescriptors.Compute(Drugs.Parse(Drugs.Ibuprofen));
        Assert.True(DrugLikeness.Veber(ibuprofen).Passed);
        Assert.True(DrugLikeness.Egan(ibuprofen).Passed);
        Assert.True(DrugLikeness.Ghose(ibuprofen).Passed);
        // atorvastatin is too big and too greasy for lead-likeness and for Ghose
        var atorva = MolecularDescriptors.Compute(Drugs.Parse(Drugs.Atorvastatin));
        Assert.False(DrugLikeness.LeadLike(atorva).Passed);
        Assert.False(DrugLikeness.Ghose(atorva).Passed);
        Assert.False(DrugLikeness.Veber(atorva).Passed);
    }

    [Fact]
    public void All_rules_returns_every_filter()
    {
        var rules = Admet.AllRules(Drugs.Parse(Drugs.Aspirin));
        Assert.Equal(["Lipinski", "Veber", "Egan", "Ghose", "Muegge", "LeadLike"], rules.Select(r => r.Rule).ToArray());
    }
}

public class AdmetHeuristicTests
{
    [Fact]
    public void Boiled_egg_flags_cns_drugs_and_rejects_polar_ones()
    {
        Assert.True(Admet.BoiledEgg(Drugs.Parse(Drugs.Diazepam)).Bbb);
        Assert.True(Admet.BoiledEgg(Drugs.Parse(Drugs.Propranolol)).Bbb);
        // a large polar sugar neither crosses nor is absorbed passively
        var sucrose = Drugs.Parse("OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O");
        Assert.False(Admet.BoiledEgg(sucrose).Bbb);
        Assert.False(Admet.BoiledEgg(sucrose).GiAbsorption);
        Assert.True(Admet.BoiledEgg(Drugs.Parse(Drugs.Aspirin)).GiAbsorption);
    }

    [Fact]
    public void Bbb_score_prefers_cns_like_molecules()
    {
        double diazepam = Admet.BbbScore(Drugs.Parse(Drugs.Diazepam));
        double vanco = Admet.BbbScore(Drugs.Parse(Drugs.Vancomycin));
        Assert.InRange(diazepam, 0, 1);
        Assert.InRange(vanco, 0, 1);
        Assert.True(diazepam > vanco);
    }

    [Fact]
    public void Herg_heuristic_flags_basic_lipophilic_amines()
    {
        var (risk, reasons) = Admet.HergRisk(Drugs.Parse("c1ccc(cc1)C(CCN(C)C)c1ccccc1")); // basic, greasy, two aryls
        Assert.True(risk > 0.5);
        Assert.Contains(reasons, r => r.Contains("basic"));
        var (lowRisk, _) = Admet.HergRisk(Drugs.Parse("CCO"));
        Assert.Equal(0, lowRisk);
    }

    [Fact]
    public void Cyp_flags_cover_the_five_isoforms()
    {
        var flags = Admet.CypInhibitionFlags(Drugs.Parse(Drugs.Celecoxib));
        Assert.Equal(["CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4"], flags.Keys.Order(StringComparer.Ordinal).ToArray());
        Assert.True(flags["CYP2C9"]); // lipophilic, non-basic, aromatic
        Assert.All(Admet.CypInhibitionFlags(Drugs.Parse("CCO")).Values, v => Assert.False(v));
    }

    [Fact]
    public void Pgp_and_complexity_heuristics()
    {
        Assert.True(Admet.LikelyPgpSubstrate(Drugs.Parse(Drugs.Vancomycin)));
        Assert.False(Admet.LikelyPgpSubstrate(Drugs.Parse(Drugs.Aspirin)));
        Assert.True(Admet.ComplexityScore(Drugs.Parse(Drugs.Morphine)) > Admet.ComplexityScore(Drugs.Parse(Drugs.Aspirin)));
        Assert.InRange(Admet.ComplexityScore(Drugs.Parse("CCO")), 1, 10);
        Assert.InRange(Admet.ComplexityScore(Drugs.Parse(Drugs.Vancomycin)), 1, 10);
    }

    [Fact]
    public void Profile_collects_everything_consistently()
    {
        var m = Drugs.Parse(Drugs.Ibuprofen, "ibuprofen");
        var p = Admet.Profile(m);
        Assert.Equal(Admet.EsolLogS(m), p.EsolLogS, 9);
        Assert.Equal(Qed.Compute(m), p.Qed, 9);
        Assert.Equal(6, p.Rules.Count);
        Assert.Contains("Lipinski", p.PassedRules);
        Assert.Equal(5, p.CypFlags.Count);
        Assert.Contains("logP", p.ToString(), StringComparison.OrdinalIgnoreCase);
    }
}
