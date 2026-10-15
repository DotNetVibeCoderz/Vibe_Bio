using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>
/// Rule-based ADMET estimates. The physicochemical rules (Lipinski, Veber, Egan, Ghose, Muegge, lead-likeness) and
/// the Delaney ESOL solubility model are published equations reproduced faithfully; everything else here
/// (BBB, hERG, CYP, P-gp) is a <b>transparent heuristic</b>, not a trained model, and is labelled as such on every
/// result. Use it for triage and as a worked example, never as a substitute for measured data or a validated
/// QSAR model.
/// </summary>
public static class Admet
{
    /// <summary>Every rule-based drug-likeness filter, evaluated at once.</summary>
    public static IReadOnlyList<RuleResult> AllRules(Molecule mol)
    {
        var d = MolecularDescriptors.Compute(mol);
        return
        [
            DrugLikeness.Lipinski(d),
            DrugLikeness.Veber(d),
            DrugLikeness.Egan(d),
            DrugLikeness.Ghose(d),
            DrugLikeness.Muegge(mol, d),
            DrugLikeness.LeadLike(d),
        ];
    }

    /// <summary>
    /// Delaney's ESOL (J. Chem. Inf. Comput. Sci. 2004, 44, 1000): log S (mol/L) =
    /// 0.16 − 0.63·clogP − 0.0062·MW + 0.066·RB − 0.74·AP, with the aromatic proportion AP = aromatic heavy
    /// atoms / heavy atoms. Crippen logP stands in for clogP, as in the common RDKit implementation.
    /// </summary>
    public static double EsolLogS(Molecule mol)
    {
        double logP = Crippen.LogP(mol);
        double mw = MolecularDescriptors.MolecularWeight(mol);
        int rb = MolecularDescriptors.RotatableBonds(mol);
        double ap = AromaticProportion(mol);
        return 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * rb - 0.74 * ap;
    }

    /// <summary>Solubility in mg/L from <see cref="EsolLogS"/>.</summary>
    public static double EsolSolubilityMgPerL(Molecule mol) => Math.Pow(10, EsolLogS(mol)) * MolecularDescriptors.MolecularWeight(mol) * 1000;

    /// <summary>Fraction of heavy atoms that are aromatic (Delaney's AP term).</summary>
    public static double AromaticProportion(Molecule mol)
    {
        int heavy = 0, aromatic = 0;
        for (int i = 0; i < mol.AtomCount; i++)
        {
            if (mol.Atoms[i].AtomicNumber <= 1) continue;
            heavy++;
            if (mol.Atoms[i].IsAromatic) aromatic++;
        }
        return heavy == 0 ? 0 : (double)aromatic / heavy;
    }

    /// <summary>
    /// Solubility classes used for ESOL reporting (log S in mol/L).
    /// </summary>
    public static string SolubilityClass(double logS) => logS switch
    {
        > 0 => "highly soluble",
        > -2 => "very soluble",
        > -4 => "soluble",
        > -6 => "moderately soluble",
        > -10 => "poorly soluble",
        _ => "insoluble",
    };

    /// <summary>
    /// Blood–brain-barrier penetration heuristic from the BOILED-Egg physicochemical space (Daina &amp; Zoete,
    /// ChemMedChem 2016): the ellipse TPSA ≤ 79 Å² and −0.4 ≤ WLOGP ≤ 6.0 marks likely brain penetration, the wider
    /// region TPSA ≤ 142 and logP ≤ 5.88 likely gastro-intestinal absorption. <b>Heuristic, not a trained model.</b>
    /// </summary>
    public static (bool Bbb, bool GiAbsorption) BoiledEgg(Molecule mol)
    {
        double tpsa = Tpsa.Compute(mol), logp = Crippen.LogP(mol);
        bool bbb = tpsa <= 79 && logp is >= -0.4 and <= 6.0;
        bool gi = tpsa <= 142 && logp <= 5.88;
        return (bbb, gi);
    }

    /// <summary>
    /// A continuous BBB score in the spirit of Clark's CNS-MPO-like sums: high for small, lipophilic, low-TPSA,
    /// low-HBD molecules. <b>Heuristic, not calibrated against measured logBB.</b> Range 0–1.
    /// </summary>
    public static double BbbScore(Molecule mol)
    {
        var d = MolecularDescriptors.Compute(mol);
        double Ramp(double x, double good, double bad) => Math.Clamp((bad - x) / (bad - good), 0, 1);
        double s = 0;
        s += Ramp(d.Tpsa, 60, 120);
        s += Ramp(d.MolecularWeight, 300, 500);
        s += Ramp(d.HBondDonors, 1, 4);
        s += Ramp(d.RotatableBonds, 4, 10);
        s += 1 - Math.Abs(Math.Clamp(d.LogP, 0, 6) - 2.5) / 3.5;
        return s / 5;
    }

    private static readonly Lazy<SmartsPattern> _basicNitrogen = new(() =>
        SmartsPattern.Parse("[$([NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[S,P]=O);!$(N#*);!$([N;R0]a)]),$([NX4+;!$(N~[O-])])]"));

    private static readonly Lazy<SmartsPattern> _aromaticRing = new(() => SmartsPattern.Parse("a1aaaaa1"));

    /// <summary>
    /// hERG liability heuristic (Aronov-style pharmacophore counting): a basic amine together with lipophilic
    /// aromatic bulk raises the flag. Returns a 0–1 risk score and the reasons. <b>Heuristic, not a trained
    /// classifier; do not use it for safety decisions.</b>
    /// </summary>
    public static (double Risk, IReadOnlyList<string> Reasons) HergRisk(Molecule mol)
    {
        var reasons = new List<string>();
        double score = 0;
        double logp = Crippen.LogP(mol);
        bool basic = _basicNitrogen.Value.IsMatch(mol);
        int aromatic = MolecularDescriptors.AromaticRingCount(mol);
        double mw = MolecularDescriptors.MolecularWeight(mol);
        if (basic) { score += 0.35; reasons.Add("basic nitrogen (likely protonated at pH 7.4)"); }
        if (logp > 3.7) { score += 0.3; reasons.Add($"logP {logp:F2} > 3.7"); }
        else if (logp > 3.0) { score += 0.15; reasons.Add($"logP {logp:F2} > 3.0"); }
        if (aromatic >= 2) { score += 0.2; reasons.Add($"{aromatic} aromatic rings"); }
        if (mw > 250 && basic && aromatic >= 2) { score += 0.15; reasons.Add("basic amine with aromatic bulk (classic hERG pharmacophore)"); }
        return (Math.Clamp(score, 0, 1), reasons);
    }

    /// <summary>
    /// CYP-inhibition heuristics for the five main isoforms. Each flag combines lipophilicity, size and a few
    /// isoform-typical motifs. <b>Heuristics, not trained models.</b>
    /// </summary>
    public static IReadOnlyDictionary<string, bool> CypInhibitionFlags(Molecule mol)
    {
        double logp = Crippen.LogP(mol), mw = MolecularDescriptors.MolecularWeight(mol), tpsa = Tpsa.Compute(mol);
        int aromatic = MolecularDescriptors.AromaticRingCount(mol);
        bool basic = _basicNitrogen.Value.IsMatch(mol);
        bool planarAromatic = aromatic >= 2 && MolecularDescriptors.FractionCsp3(mol) < 0.3;
        return new Dictionary<string, bool>(StringComparer.OrdinalIgnoreCase)
        {
            ["CYP1A2"] = planarAromatic && mw < 400 && tpsa < 80,
            ["CYP2C9"] = logp > 3 && tpsa < 90 && aromatic >= 1 && !basic,
            ["CYP2C19"] = logp > 2.5 && aromatic >= 1 && mw < 500,
            ["CYP2D6"] = basic && logp > 2.5 && aromatic >= 1,
            ["CYP3A4"] = mw > 350 && logp > 3 && aromatic >= 2,
        };
    }

    /// <summary>
    /// P-glycoprotein substrate heuristic (large, polar, flexible molecules are typical substrates).
    /// <b>Heuristic, not a trained model.</b>
    /// </summary>
    public static bool LikelyPgpSubstrate(Molecule mol)
    {
        var d = MolecularDescriptors.Compute(mol);
        return d.MolecularWeight > 400 && d.Tpsa > 75 && (d.HBondAcceptors >= 6 || d.RotatableBonds >= 8);
    }

    /// <summary>
    /// Synthetic-accessibility-like complexity penalty in 1 (easy) to 10 (hard), from ring complexity, stereo
    /// centres, macrocycles, spiro atoms and size. <b>Heuristic</b>, inspired by Ertl &amp; Schuffenhauer's SA score
    /// but without their fragment-frequency term.
    /// </summary>
    public static double ComplexityScore(Molecule mol)
    {
        var rings = mol.Rings;
        int heavy = MolecularDescriptors.HeavyAtomCount(mol);
        if (heavy == 0) return 1;
        int stereo = 0;
        for (int i = 0; i < mol.AtomCount; i++) if (mol.Atoms[i].Chirality != ChiralTag.None) stereo++;
        int macro = rings.AtomRings.Count(r => r.Length > 8);
        int bridged = rings.RingSystems.Count(s => s.Length >= 3);
        int spiro = Enumerable.Range(0, mol.AtomCount).Count(i => rings.AtomRingCount(i) >= 2 && rings.RingBondCount(i) == 4);
        double score = 1
            + 0.8 * Math.Log(1 + heavy / 10.0)
            + 0.5 * stereo
            + 1.5 * macro
            + 0.7 * bridged
            + 0.6 * spiro
            + 0.5 * rings.Count;
        return Math.Clamp(score, 1, 10);
    }

    /// <summary>Computes the full <see cref="AdmetProfile"/>.</summary>
    public static AdmetProfile Profile(Molecule mol)
    {
        var d = MolecularDescriptors.Compute(mol);
        double logS = EsolLogS(mol);
        var (bbb, gi) = BoiledEgg(mol);
        var (herg, hergReasons) = HergRisk(mol);
        return new AdmetProfile
        {
            Descriptors = d,
            Rules = AllRules(mol),
            EsolLogS = logS,
            EsolSolubilityMgPerL = Math.Pow(10, logS) * d.MolecularWeight * 1000,
            SolubilityClass = SolubilityClass(logS),
            BloodBrainBarrier = bbb,
            GiAbsorption = gi,
            BbbScore = BbbScore(mol),
            HergRisk = herg,
            HergReasons = hergReasons,
            CypFlags = CypInhibitionFlags(mol),
            LikelyPgpSubstrate = LikelyPgpSubstrate(mol),
            Complexity = ComplexityScore(mol),
            Qed = Qed.Compute(mol),
        };
    }
}

/// <summary>
/// A complete ADMET profile. The rule results and <see cref="EsolLogS"/> implement published models; the
/// <see cref="BbbScore"/>, <see cref="HergRisk"/>, <see cref="CypFlags"/>, <see cref="LikelyPgpSubstrate"/> and
/// <see cref="Complexity"/> fields are documented heuristics, not trained predictors.
/// </summary>
public sealed record AdmetProfile
{
    /// <summary>Standard 2-D descriptors.</summary>
    public required DescriptorSet Descriptors { get; init; }

    /// <summary>Lipinski, Veber, Egan, Ghose, Muegge and lead-likeness results.</summary>
    public required IReadOnlyList<RuleResult> Rules { get; init; }

    /// <summary>Delaney ESOL log S (mol/L).</summary>
    public required double EsolLogS { get; init; }

    /// <summary>Predicted aqueous solubility in mg/L.</summary>
    public required double EsolSolubilityMgPerL { get; init; }

    /// <summary>Qualitative solubility class.</summary>
    public required string SolubilityClass { get; init; }

    /// <summary>BOILED-Egg yolk: likely brain penetration (heuristic).</summary>
    public required bool BloodBrainBarrier { get; init; }

    /// <summary>BOILED-Egg white: likely passive GI absorption (heuristic).</summary>
    public required bool GiAbsorption { get; init; }

    /// <summary>Continuous 0–1 BBB desirability (heuristic).</summary>
    public required double BbbScore { get; init; }

    /// <summary>0–1 hERG liability score (heuristic).</summary>
    public required double HergRisk { get; init; }

    /// <summary>Why the hERG score was raised.</summary>
    public required IReadOnlyList<string> HergReasons { get; init; }

    /// <summary>Per-isoform CYP inhibition flags (heuristic).</summary>
    public required IReadOnlyDictionary<string, bool> CypFlags { get; init; }

    /// <summary>P-glycoprotein substrate flag (heuristic).</summary>
    public required bool LikelyPgpSubstrate { get; init; }

    /// <summary>1–10 synthetic complexity (heuristic).</summary>
    public required double Complexity { get; init; }

    /// <summary>Quantitative estimate of drug-likeness.</summary>
    public required double Qed { get; init; }

    /// <summary>Rules that passed.</summary>
    public IEnumerable<string> PassedRules => Rules.Where(r => r.Passed).Select(r => r.Rule);

    /// <summary>Rules that failed, with their violations.</summary>
    public IEnumerable<RuleResult> FailedRules => Rules.Where(r => !r.Passed);

    /// <inheritdoc/>
    public override string ToString()
        => $"MW {Descriptors.MolecularWeight:F1}, logP {Descriptors.LogP:F2}, TPSA {Descriptors.Tpsa:F1}, QED {Qed:F2}, "
         + $"log S {EsolLogS:F2} ({SolubilityClass}); passes {string.Join("/", PassedRules)}";
}
