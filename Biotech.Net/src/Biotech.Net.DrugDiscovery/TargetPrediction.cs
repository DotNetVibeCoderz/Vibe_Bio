using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>A known ligand of a target, used as reference data for <see cref="TargetPredictor"/>.</summary>
/// <param name="Target">Target name.</param>
/// <param name="LigandId">Ligand name.</param>
/// <param name="Smiles">Ligand structure.</param>
public sealed record ReferenceLigand(string Target, string LigandId, string Smiles);

/// <summary>One predicted target.</summary>
/// <param name="Target">Target name.</param>
/// <param name="RawScore">Sum of the similarities above the threshold.</param>
/// <param name="ZScore">Raw score standardised against the random-background model.</param>
/// <param name="PValue">Probability of seeing this raw score by chance (extreme-value model).</param>
/// <param name="MaxSimilarity">Similarity to the nearest reference ligand of the target.</param>
/// <param name="NearestLigand">The nearest reference ligand.</param>
/// <param name="ReferenceCount">How many reference ligands the target has.</param>
public sealed record TargetPrediction(string Target, double RawScore, double ZScore, double PValue,
    double MaxSimilarity, string NearestLigand, int ReferenceCount)
{
    /// <inheritdoc/>
    public override string ToString() => $"{Target}: max Tc {MaxSimilarity:F2} to {NearestLigand}, z {ZScore:F1}, p {PValue:G3}";
}

/// <summary>
/// A small reference set of ligand→target pairs for demonstrating <see cref="TargetPredictor"/>.
/// </summary>
/// <remarks>
/// <b>This is a demonstration data set, not a pharmacological database.</b> It holds a few well-known drugs per
/// target, hand-written from common knowledge, and is far too small and too biased for real target prediction.
/// Supply your own reference set (e.g. ChEMBL actives at a defined activity cut-off) through
/// <see cref="TargetPredictor.FromReferences"/> for anything beyond a worked example.
/// </remarks>
public static class DemoReferenceSet
{
    /// <summary>The demo ligand→target pairs.</summary>
    public static IReadOnlyList<ReferenceLigand> Ligands { get; } =
    [
        // Cyclooxygenase inhibitors (NSAIDs)
        new("COX", "aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        new("COX", "salicylic acid", "O=C(O)c1ccccc1O"),
        new("COX", "ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        new("COX", "naproxen", "COc1ccc2cc(ccc2c1)C(C)C(=O)O"),
        new("COX", "ketoprofen", "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1"),
        new("COX", "flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1"),
        new("COX", "diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
        new("COX2", "celecoxib", "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"),
        new("COX2", "rofecoxib", "O=C1OCC(=C1c1ccccc1)c1ccc(cc1)S(C)(=O)=O"),
        // Adenosine-receptor / xanthine ligands
        new("AdenosineReceptor", "caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
        new("AdenosineReceptor", "theophylline", "Cn1c2nc[nH]c2c(=O)n(C)c1=O"),
        new("AdenosineReceptor", "theobromine", "Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
        new("AdenosineReceptor", "paraxanthine", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"),
        // Beta-adrenergic blockers
        new("BetaAdrenergic", "propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
        new("BetaAdrenergic", "atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
        new("BetaAdrenergic", "metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1"),
        new("BetaAdrenergic", "pindolol", "CC(C)NCC(O)COc1cccc2[nH]ccc12"),
        // Benzodiazepine (GABA-A) ligands
        new("GabaA", "diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
        new("GabaA", "lorazepam", "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O"),
        new("GabaA", "oxazepam", "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O"),
        // Opioid receptor ligands
        new("OpioidReceptor", "morphine", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
        new("OpioidReceptor", "codeine", "COc1ccc2CC3N(C)CCC45C(Oc1c24)C(O)C=CC35"),
        new("OpioidReceptor", "naloxone", "C=CCN1CCC23c4c5ccc(O)c4OC2C(=O)CCC3(O)C1C5"),
        // Dihydrofolate reductase
        new("DHFR", "methotrexate", "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1"),
        new("DHFR", "trimethoprim", "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC"),
        new("DHFR", "pyrimethamine", "CCc1nc(N)nc(N)c1-c1ccc(Cl)cc1"),
        // Beta-lactam antibiotics (PBP)
        new("PenicillinBindingProtein", "penicillin G", "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"),
        new("PenicillinBindingProtein", "amoxicillin", "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O"),
        new("PenicillinBindingProtein", "ampicillin", "CC1(C)SC2C(NC(=O)C(N)c3ccccc3)C(=O)N2C1C(=O)O"),
        // HMG-CoA reductase (statins)
        new("HmgCoaReductase", "atorvastatin", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"),
        new("HmgCoaReductase", "fluvastatin", "CC(C)n1c(C=CC(O)CC(O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21"),
        // Kinase inhibitors
        new("TyrosineKinase", "imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
        new("TyrosineKinase", "nilotinib", "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1"),
        new("TyrosineKinase", "gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
        new("TyrosineKinase", "erlotinib", "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"),
        // Proton-pump inhibitors
        new("ProtonPump", "omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
        new("ProtonPump", "lansoprazole", "Cc1c(OCC(F)(F)F)ccnc1CS(=O)c1nc2ccccc2[nH]1"),
        // Carbonic anhydrase
        new("CarbonicAnhydrase", "acetazolamide", "CC(=O)Nc1nnc(S(N)(=O)=O)s1"),
        new("CarbonicAnhydrase", "sulfanilamide", "NS(=O)(=O)c1ccc(N)cc1"),
        new("CarbonicAnhydrase", "hydrochlorothiazide", "NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O"),
    ];

    /// <summary>Distinct target names.</summary>
    public static IReadOnlyList<string> Targets { get; } = Ligands.Select(l => l.Target).Distinct().ToList();
}

/// <summary>
/// Ligand-based target prediction in the style of the Similarity Ensemble Approach (Keiser et al., Nat. Biotechnol.
/// 2007, 25, 197): the query is compared with every reference ligand of a target, similarities above a threshold
/// are summed into a raw set score, and that score is standardised against a random background so that targets with
/// many reference ligands do not win automatically. The significance uses an extreme-value model, as in SEA.
/// </summary>
/// <remarks>
/// With <see cref="Demo"/> the reference data is <see cref="DemoReferenceSet"/>, a tiny illustrative set — the
/// predictions are a worked example, not pharmacology. Build a predictor over real actives with
/// <see cref="FromReferences"/> for genuine use.
/// </remarks>
public sealed class TargetPredictor
{
    private readonly Dictionary<string, List<(string Id, BitVector Fp)>> _byTarget = [];
    private readonly Func<Molecule, BitVector> _fingerprinter;
    private double _backgroundMean, _backgroundStd;

    private TargetPredictor(Func<Molecule, BitVector> fingerprinter) => _fingerprinter = fingerprinter;

    /// <summary>Similarity below which a pair contributes nothing to the raw score (SEA's threshold).</summary>
    public double SimilarityThreshold { get; init; } = 0.28;

    /// <summary>The targets covered.</summary>
    public IReadOnlyCollection<string> Targets => _byTarget.Keys;

    /// <summary>Reference ligands of <paramref name="target"/>.</summary>
    public IReadOnlyList<string> LigandsOf(string target) => _byTarget.TryGetValue(target, out var l) ? l.Select(x => x.Id).ToList() : [];

    /// <summary>Builds a predictor from ligand→target pairs.</summary>
    public static TargetPredictor FromReferences(IEnumerable<ReferenceLigand> references, Func<Molecule, BitVector>? fingerprinter = null,
        double similarityThreshold = 0.28)
    {
        var p = new TargetPredictor(fingerprinter ?? (m => MorganFingerprint.Ecfp4(m))) { SimilarityThreshold = similarityThreshold };
        foreach (var r in references)
        {
            var mol = Molecule.FromSmiles(r.Smiles, r.LigandId);
            if (!p._byTarget.TryGetValue(r.Target, out var list)) p._byTarget[r.Target] = list = [];
            list.Add((r.LigandId, p._fingerprinter(mol)));
        }
        p.CalibrateBackground();
        return p;
    }

    /// <summary>A predictor over the embedded demo reference set. <b>Demonstration data only.</b></summary>
    public static TargetPredictor Demo(Func<Molecule, BitVector>? fingerprinter = null)
        => FromReferences(DemoReferenceSet.Ligands, fingerprinter);

    /// <summary>
    /// Mean and standard deviation of the pairwise similarity of the reference ligands to each other, used as the
    /// random background for the z-score.
    /// </summary>
    private void CalibrateBackground()
    {
        var all = _byTarget.SelectMany(kv => kv.Value.Select(v => (kv.Key, v.Fp))).ToList();
        var sims = new List<double>();
        for (int i = 0; i < all.Count; i++)
            for (int j = i + 1; j < all.Count; j++)
            {
                if (all[i].Key == all[j].Key) continue; // different targets only: an unrelated-pair background
                double s = Similarity.Tanimoto(all[i].Fp, all[j].Fp);
                sims.Add(s < SimilarityThreshold ? 0 : s);
            }
        if (sims.Count == 0) { _backgroundMean = 0; _backgroundStd = 1; return; }
        _backgroundMean = sims.Average();
        double v = sims.Sum(s => (s - _backgroundMean) * (s - _backgroundMean)) / Math.Max(1, sims.Count - 1);
        _backgroundStd = Math.Max(1e-6, Math.Sqrt(v));
    }

    /// <summary>Predicts targets for <paramref name="query"/>, most significant first.</summary>
    public IReadOnlyList<TargetPrediction> Predict(Molecule query, int topK = int.MaxValue)
    {
        var qfp = _fingerprinter(query);
        var results = new List<TargetPrediction>();
        foreach (var (target, ligands) in _byTarget)
        {
            double raw = 0, max = 0;
            string nearest = "";
            foreach (var (id, fp) in ligands)
            {
                double s = Similarity.Tanimoto(qfp, fp);
                if (s > max) { max = s; nearest = id; }
                if (s >= SimilarityThreshold) raw += s;
            }
            // expected raw score if the ligands were unrelated to the query
            double expected = _backgroundMean * ligands.Count;
            double sd = _backgroundStd * Math.Sqrt(ligands.Count);
            double z = sd < 1e-9 ? 0 : (raw - expected) / sd;
            results.Add(new TargetPrediction(target, raw, z, ExtremeValuePValue(z), max, nearest, ligands.Count));
        }
        return results.OrderByDescending(r => r.ZScore).ThenByDescending(r => r.MaxSimilarity).Take(topK).ToList();
    }

    /// <summary>Predicted targets whose z-score is at least <paramref name="minZ"/>.</summary>
    public IReadOnlyList<TargetPrediction> PredictSignificant(Molecule query, double minZ = 3.0)
        => Predict(query).Where(r => r.ZScore >= minZ).ToList();

    /// <summary>
    /// SEA's extreme-value survival function P(Z ≥ z) = 1 − exp(−e<sup>−x</sup>) with x = −(z·π/√6 + Γ), the
    /// standardised Gumbel form.
    /// </summary>
    public static double ExtremeValuePValue(double z)
    {
        const double euler = 0.577215664901533;
        double x = -(z * Math.PI / Math.Sqrt(6) + euler);
        double e = Math.Exp(x);
        return e > 30 ? 1 : 1 - Math.Exp(-Math.Exp(x));
    }
}
