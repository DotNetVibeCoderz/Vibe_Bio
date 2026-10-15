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

## Limitations

- **ADMET is heuristic.** Lipinski, Veber, Egan, Ghose, Muegge, lead-likeness and the Delaney ESOL equation are published rules reproduced faithfully. BBB, hERG, CYP inhibition and P-gp are transparent rule-of-thumb heuristics, not trained models. Use them for triage only, never in place of measured data.
- **The target-prediction reference set is a demo.** `DemoReferenceSet` holds a handful of hand-picked drugs per target and is far too small and biased for real predictions. Supply ChEMBL-derived references through `TargetPredictor.FromReferences`.
- **Docking is pose scoring, not docking.** `VinaScore` scores a given pose and offers a rigid-body local optimisation. It does no conformer search, flexible docking or 3-D embedding (3-D conformers are on the [roadmap](../roadmap.md)). Atom typing is simplified, so absolute values are indicative; ranking poses of the same ligand is the intended use.
- The QSAR example uses toy data. Cross-validated R² on ten molecules says nothing about real predictive power.

---
Previous: [Chemistry](chemistry.md) · Next: [Synthetic biology & CRISPR](synthetic-biology.md) · [Tutorial index](index.md)
