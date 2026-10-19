# Tutorial: Drug discovery

Package: `Gravicode.Biotech.Net.DrugDiscovery` (builds on Chemistry and StructuralBiology)

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Virtual screening (multi-fingerprint similarity with data fusion), target prediction (SEA-style), ADMET triage, structural alerts (PAINS, Brenk), compound filters, pharmacophore counts, lead optimisation (bioisosteres plus MPO scoring), QSAR (ridge, random forest) and a Vina-like pose scoring function.

## 1. Virtual screening

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

var library = new List<LibraryCompound>
{
    LibraryCompound.FromSmiles("salicylic",   "OC(=O)c1ccccc1O"),
    LibraryCompound.FromSmiles("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    LibraryCompound.FromSmiles("ibuprofen",   "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    LibraryCompound.FromSmiles("naproxen",    "COc1ccc2cc(ccc2c1)C(C)C(=O)O"),
    LibraryCompound.FromSmiles("caffeine",    "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
};
var query = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin");

var hits = SimilarityScreen.Run(query, library, new SimilarityScreenOptions { TopK = 3 });
foreach (var h in hits) Console.WriteLine($"{h.Rank}. {h.Id,-12} {h.Score:F3}");

// Benchmark a screen against known actives
var actives = new HashSet<string> { "salicylic" };
Console.WriteLine($"ROC AUC {SimilarityScreen.RocAuc(SimilarityScreen.Run(query, library), actives):F2}");
```

## 2. ADMET triage, alerts and filters

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

var mol = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
AdmetProfile admet = Admet.Profile(mol);
Console.WriteLine(admet);                       // summary text
Console.WriteLine($"ESOL logS {admet.EsolLogS:F2} ({admet.SolubilityClass}), BBB {admet.BloodBrainBarrier}, GI {admet.GiAbsorption}");
Console.WriteLine($"hERG risk {admet.HergRisk:F2}, QED {admet.Qed:F2}");

var rhodanine = Molecule.FromSmiles("O=C1NC(=S)SC1=Cc1ccccc1", "rhodanine-like");
foreach (var alert in StructuralAlerts.Match(rhodanine)) Console.WriteLine($"{alert.Family}: {alert.Name}");

var filter = new CompoundFilter().RequireLipinski().ExcludeAlerts(StructuralAlerts.Pains);
Console.WriteLine(filter.Evaluate(rhodanine));
Console.WriteLine(filter.Evaluate(mol));
```

## 3. Target prediction

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

// Demo reference set: a few well-known ligands per target (see Limitations)
var predictor = TargetPredictor.Demo();
var query = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
foreach (var p in predictor.Predict(query, topK: 3)) Console.WriteLine(p);

// For real use, build from your own actives (e.g. ChEMBL at a pChEMBL cut-off):
//   TargetPredictor.FromReferences(refs)   where refs is IEnumerable<ReferenceLigand>(target, id, smiles)
```

## 4. Lead optimisation and QSAR

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

var lead = Molecule.FromSmiles("OC(=O)c1ccc(Cl)cc1", "lead");
foreach (var (analogue, molecule, score) in LeadOptimizer.Optimize(lead, MpoProfile.Oral).Take(5))
    Console.WriteLine($"{score:F3}  {molecule.ToSmiles()}  via {analogue?.Transform.Name ?? "(parent)"}");

// QSAR: descriptors + fingerprints -> random forest, with k-fold cross-validation
string[] smiles = ["CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO", "CC(C)O", "CC(C)CO", "OCCO", "OCC(O)CO", "CCCCCCCO"];
double[] logP = [-0.31, 0.25, 0.88, 1.51, 2.03, 0.05, 0.76, -1.36, -1.76, 2.62];
var compounds = smiles.Select((s, i) => LibraryCompound.FromSmiles($"c{i}", s)).ToList();
var model = QsarModel.Build(compounds, logP, folds: 5);
Console.WriteLine(model.CrossValidation);
Console.WriteLine($"predicted logP of octanol: {model.Predict(Molecule.FromSmiles("CCCCCCCCO")):F2}");
```

## 5. Pose scoring

`VinaScore` implements an AutoDock Vina-like empirical scoring function. It needs a ligand with 3-D coordinates, for example a HETATM ligand from a PDB complex:

```csharp
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

// With a protein–ligand complex loaded from disk:
//   var complex = PdbReader.Read("complex.pdb");
//   var ligand  = complex.Ligands.First();
//   var pose    = LigandPose.FromResidue(ligand);
//   var score   = VinaScore.Score(pose, complex);              // kcal/mol-like affinity
//   var pocket  = PocketDetector.Detect(complex)[0];
//   var (best, s) = VinaScore.ScoreInPocket(pose, complex, pocket);   // rigid-body local search
var receptor = SampleStructures.Crambin();
Console.WriteLine($"{PocketDetector.Detect(receptor).Count} candidate pockets on {receptor.Id}");
```

## 6. Flexible docking

`VinaDocking` reimplements the AutoDock Vina 1.2 method: X-Score atom typing, grid maps with trilinear interpolation, a torsion-tree ligand, Monte-Carlo search with BFGS local optimisation, RMSD clustering, and a final refinement and rescoring of every pose with exact receptor interactions. The receptor is rigid; give it hydrogens (for example from PDBFixer or Reduce) for the most faithful donor typing, otherwise residue templates are used.

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

// Streptavidin pocket and the crystal biotin from tests/Biotech.Net.DrugDiscovery.Tests/Data (PDB 1STP)
VinaReceptor receptor = VinaReceptor.FromPdb("1stp_receptor.pdb");
Molecule crystal = SdfReader.ReadFile("1stp_ligand.sdf").First();
VinaPoseScore score = VinaDocking.Score(crystal, receptor);           // like `vina --score_only`
Console.WriteLine($"crystal pose: {score}");

// Dock biotin starting from a generated 3-D conformer
var biotin = Molecule.FromSmiles("OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12", "biotin");
var options = new DockingOptions { Center = new Vec3(11.1, 1.7, -10.8), Size = new Vec3(22.5, 22.5, 22.5), Exhaustiveness = 8 };
DockingResult result = VinaDocking.Dock(ConformerGenerator.Embed(biotin), receptor, options);
foreach (DockedPose pose in result.Poses) Console.WriteLine(pose);   // #1 -7.4 kcal/mol ...
File.WriteAllText("biotin_docked.sdf", result.ToSdf());               // poses with VINA_AFFINITY etc.
File.WriteAllText("biotin_docked.pdbqt", result.ToPdbqt());           // Vina-style multi-model PDBQT
```

Without a known site, `PocketDetector.Detect(structure)[0].Center` gives a box centre.

Validation against AutoDock Vina 1.2.7 on two public complexes (fixtures and generator script in `tests/Biotech.Net.DrugDiscovery.Tests/Data`; receptors protonated with PDBFixer and cropped to the pocket):

| Complex | Crystal-pose affinity, Biotech.Net / Vina | Redock from crystal geometry: top pose RMSD, best affinity (Vina) |
| --- | --- | --- |
| 3PTB trypsin–benzamidine | −5.806 / −5.806 kcal/mol | < 1 Å, −6.20 (Vina 0.38 Å, −6.20) |
| 1STP streptavidin–biotin | −6.582 / −6.576 kcal/mol | 1.0 Å, −7.02 (Vina 0.57 Å, −7.09) |

Vina's own best poses rescored by Biotech.Net agree to 0.001 kcal/mol. Redocking biotin from a freshly generated conformer places the top pose about 1.3 Å from the crystal ligand. With exhaustiveness 8, grid maps take about 0.3 s and the search 2 s (benzamidine) to 9 s (biotin) on a desktop CPU.

## Limitations

- **ADMET is heuristic.** Lipinski, Veber, Egan, Ghose, Muegge, lead-likeness and the Delaney ESOL equation are published rules reproduced faithfully. BBB, hERG, CYP inhibition and P-gp are transparent rule-of-thumb heuristics, not trained models. Use them for triage only, never in place of measured data.
- **The target-prediction reference set is a demo.** `DemoReferenceSet` holds a handful of hand-picked drugs per target and is far too small and biased for real predictions. Supply ChEMBL-derived references through `TargetPredictor.FromReferences`.
- **`VinaScore` is a quick scorer; `VinaDocking` is the Vina method.** `VinaScore` uses simplified element typing, so its values are indicative. `VinaDocking` follows Vina closely but is not bit-identical: map pair energies come from 0.005 Å tables (Vina samples r²), random streams differ, and flexible receptor side chains, hydrated docking, macrocycle closure and the AD4/Vinardo scoring functions are not implemented. The intramolecular energy differs from Vina's by a few hundredths of a kcal/mol because pairs whose distance cannot change are skipped. The PDBQT writer omits branches that only rotate hydrogens, so Vina re-scoring such a file counts slightly fewer rotors. Like Vina, scores are a ranking aid, not binding free energies.
- The QSAR example uses toy data. Cross-validated R² on ten molecules says nothing about real predictive power.

---
Previous: [Chemistry](chemistry.md) · Next: [Synthetic biology & CRISPR](synthetic-biology.md) · [Tutorial index](index.md)
