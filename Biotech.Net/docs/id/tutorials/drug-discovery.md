# Tutorial: Penemuan obat

Paket: `Gravicode.Biotech.Net.DrugDiscovery` (berdiri di atas Chemistry dan StructuralBiology)

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Penapisan virtual (kemiripan multi-sidik-jari dengan fusi data), prediksi target (gaya SEA), triase ADMET, peringatan struktural (PAINS, Brenk), filter senyawa, hitungan farmakofor, optimasi lead (bioisoster plus skor MPO), QSAR (ridge, random forest), serta fungsi penilaian pose mirip Vina.

## 1. Penapisan virtual

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

// Mengukur mutu penapisan terhadap senyawa aktif yang diketahui
var actives = new HashSet<string> { "salicylic" };
Console.WriteLine($"ROC AUC {SimilarityScreen.RocAuc(SimilarityScreen.Run(query, library), actives):F2}");
```

## 2. Triase ADMET, peringatan, dan filter

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

var mol = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
AdmetProfile admet = Admet.Profile(mol);
Console.WriteLine(admet);                       // ringkasan teks
Console.WriteLine($"ESOL logS {admet.EsolLogS:F2} ({admet.SolubilityClass}), BBB {admet.BloodBrainBarrier}, GI {admet.GiAbsorption}");
Console.WriteLine($"risiko hERG {admet.HergRisk:F2}, QED {admet.Qed:F2}");

var rhodanine = Molecule.FromSmiles("O=C1NC(=S)SC1=Cc1ccccc1", "rhodanine-like");
foreach (var alert in StructuralAlerts.Match(rhodanine)) Console.WriteLine($"{alert.Family}: {alert.Name}");

var filter = new CompoundFilter().RequireLipinski().ExcludeAlerts(StructuralAlerts.Pains);
Console.WriteLine(filter.Evaluate(rhodanine));
Console.WriteLine(filter.Evaluate(mol));
```

## 3. Prediksi target

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

// Set acuan demo: beberapa ligan terkenal per target (lihat bagian Batasan)
var predictor = TargetPredictor.Demo();
var query = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
foreach (var p in predictor.Predict(query, topK: 3)) Console.WriteLine(p);

// Untuk penggunaan nyata, bangun dari senyawa aktif Anda sendiri (misalnya ChEMBL pada ambang pChEMBL):
//   TargetPredictor.FromReferences(refs)   dengan refs berupa IEnumerable<ReferenceLigand>(target, id, smiles)
```

## 4. Optimasi lead dan QSAR

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

var lead = Molecule.FromSmiles("OC(=O)c1ccc(Cl)cc1", "lead");
foreach (var (analogue, molecule, score) in LeadOptimizer.Optimize(lead, MpoProfile.Oral).Take(5))
    Console.WriteLine($"{score:F3}  {molecule.ToSmiles()}  via {analogue?.Transform.Name ?? "(induk)"}");

// QSAR: deskriptor + sidik jari -> random forest, dengan validasi silang k-fold
string[] smiles = ["CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO", "CC(C)O", "CC(C)CO", "OCCO", "OCC(O)CO", "CCCCCCCO"];
double[] logP = [-0.31, 0.25, 0.88, 1.51, 2.03, 0.05, 0.76, -1.36, -1.76, 2.62];
var compounds = smiles.Select((s, i) => LibraryCompound.FromSmiles($"c{i}", s)).ToList();
var model = QsarModel.Build(compounds, logP, folds: 5);
Console.WriteLine(model.CrossValidation);
Console.WriteLine($"prediksi logP oktanol: {model.Predict(Molecule.FromSmiles("CCCCCCCCO")):F2}");
```

## 5. Penilaian pose

`VinaScore` mengimplementasikan fungsi penilaian empiris mirip AutoDock Vina. Ia memerlukan ligan dengan koordinat 3-D, misalnya ligan HETATM dari kompleks PDB:

```csharp
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

// Dengan kompleks protein-ligan dari disk:
//   var complex = PdbReader.Read("complex.pdb");
//   var ligand  = complex.Ligands.First();
//   var pose    = LigandPose.FromResidue(ligand);
//   var score   = VinaScore.Score(pose, complex);              // afinitas berskala kcal/mol
//   var pocket  = PocketDetector.Detect(complex)[0];
//   var (best, s) = VinaScore.ScoreInPocket(pose, complex, pocket);   // pencarian lokal benda tegar
var receptor = SampleStructures.Crambin();
Console.WriteLine($"{PocketDetector.Detect(receptor).Count} kandidat kantong pada {receptor.Id}");
```

## Batasan

- **ADMET bersifat heuristik.** Lipinski, Veber, Egan, Ghose, Muegge, lead-likeness, dan persamaan ESOL Delaney adalah aturan terpublikasi yang direproduksi setia. BBB, hERG, inhibisi CYP, dan P-gp hanyalah heuristik transparan, bukan model terlatih. Gunakan untuk triase saja, jangan menggantikan data terukur.
- **Set acuan prediksi target hanyalah demo.** `DemoReferenceSet` berisi segelintir obat terkenal per target dan terlalu kecil serta bias untuk prediksi nyata. Pakailah acuan turunan ChEMBL lewat `TargetPredictor.FromReferences`.
- **Ini penilaian pose, bukan docking.** `VinaScore` menilai pose yang diberikan dan menyediakan optimasi lokal benda tegar. Tidak ada pencarian konformer, docking fleksibel, maupun penanaman 3-D (konformer 3-D ada di [roadmap](../roadmap.md)). Penentuan tipe atomnya disederhanakan, jadi nilai absolutnya bersifat indikatif; peruntukannya adalah memeringkat pose dari ligan yang sama.
- Contoh QSAR memakai data mainan. R² hasil validasi silang pada sepuluh molekul tidak berarti apa-apa untuk daya prediksi nyata.

---
Sebelumnya: [Kimia](chemistry.md) · Selanjutnya: [Biologi sintetik](synthetic-biology.md) · [Indeks tutorial](index.md)
