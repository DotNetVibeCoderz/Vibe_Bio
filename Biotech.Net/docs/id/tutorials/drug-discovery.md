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

## 6. Docking fleksibel

`VinaDocking` mengimplementasikan ulang metode AutoDock Vina 1.2: penentuan tipe atom X-Score, peta grid dengan interpolasi trilinear, ligan berbentuk pohon torsi, pencarian Monte-Carlo dengan optimasi lokal BFGS, pengelompokan RMSD, serta penyempurnaan dan penilaian ulang akhir setiap pose dengan interaksi reseptor yang eksak. Reseptornya kaku; beri hidrogen (misalnya dari PDBFixer atau Reduce) agar penentuan donor paling setia, jika tidak dipakai templat residu.

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;

// Kantong streptavidin dan biotin kristal dari tests/Biotech.Net.DrugDiscovery.Tests/Data (PDB 1STP)
VinaReceptor receptor = VinaReceptor.FromPdb("1stp_receptor.pdb");
Molecule crystal = SdfReader.ReadFile("1stp_ligand.sdf").First();
VinaPoseScore score = VinaDocking.Score(crystal, receptor);           // seperti `vina --score_only`
Console.WriteLine($"pose kristal: {score}");

// Docking biotin mulai dari konformer 3-D hasil bangkitan
var biotin = Molecule.FromSmiles("OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12", "biotin");
var options = new DockingOptions { Center = new Vec3(11.1, 1.7, -10.8), Size = new Vec3(22.5, 22.5, 22.5), Exhaustiveness = 8 };
DockingResult result = VinaDocking.Dock(ConformerGenerator.Embed(biotin), receptor, options);
foreach (DockedPose pose in result.Poses) Console.WriteLine(pose);   // #1 -7.4 kcal/mol ...
File.WriteAllText("biotin_docked.sdf", result.ToSdf());               // pose dengan VINA_AFFINITY dll.
File.WriteAllText("biotin_docked.pdbqt", result.ToPdbqt());           // PDBQT multi-model gaya Vina
```

Tanpa situs yang diketahui, `PocketDetector.Detect(structure)[0].Center` memberi pusat kotak.

Validasi terhadap AutoDock Vina 1.2.7 pada dua kompleks publik (fixture dan skrip pembangkitnya di `tests/Biotech.Net.DrugDiscovery.Tests/Data`; reseptor diprotonasi dengan PDBFixer dan dipangkas ke kantong):

| Kompleks | Afinitas pose kristal, Biotech.Net / Vina | Redocking dari geometri kristal: RMSD pose teratas, afinitas terbaik (Vina) |
| --- | --- | --- |
| 3PTB tripsin–benzamidin | −5,806 / −5,806 kkal/mol | < 1 Å, −6,20 (Vina 0,38 Å, −6,20) |
| 1STP streptavidin–biotin | −6,582 / −6,576 kkal/mol | 1,0 Å, −7,02 (Vina 0,57 Å, −7,09) |

Pose terbaik Vina sendiri yang dinilai ulang oleh Biotech.Net sesuai hingga 0,001 kkal/mol. Redocking biotin dari konformer yang baru dibangkitkan menempatkan pose teratas sekitar 1,3 Å dari ligan kristal. Dengan exhaustiveness 8, peta grid butuh sekitar 0,3 detik dan pencarian 2 detik (benzamidin) hingga 9 detik (biotin) pada CPU desktop.

## Batasan

- **ADMET bersifat heuristik.** Lipinski, Veber, Egan, Ghose, Muegge, lead-likeness, dan persamaan ESOL Delaney adalah aturan terpublikasi yang direproduksi setia. BBB, hERG, inhibisi CYP, dan P-gp hanyalah heuristik transparan, bukan model terlatih. Gunakan untuk triase saja, jangan menggantikan data terukur.
- **Set acuan prediksi target hanyalah demo.** `DemoReferenceSet` berisi segelintir obat terkenal per target dan terlalu kecil serta bias untuk prediksi nyata. Pakailah acuan turunan ChEMBL lewat `TargetPredictor.FromReferences`.
- **`VinaScore` adalah penilai cepat; `VinaDocking` adalah metode Vina.** `VinaScore` memakai penentuan tipe unsur yang disederhanakan, jadi nilainya indikatif. `VinaDocking` mengikuti Vina dengan dekat tetapi tidak identik bit demi bit: energi pasangan untuk peta berasal dari tabel 0,005 Å (Vina menyampel r²), aliran bilangan acaknya berbeda, dan rantai samping reseptor yang fleksibel, docking terhidrasi, penutupan makrosiklus, serta fungsi skor AD4/Vinardo belum diimplementasikan. Energi intramolekulnya berbeda dari Vina beberapa perseratus kkal/mol karena pasangan yang jaraknya tidak bisa berubah dilewati. Penulis PDBQT menghilangkan cabang yang hanya memutar hidrogen, sehingga penilaian ulang berkas itu oleh Vina menghitung sedikit lebih sedikit rotor. Seperti Vina, skor ini alat bantu pemeringkatan, bukan energi bebas pengikatan.
- Contoh QSAR memakai data mainan. R² hasil validasi silang pada sepuluh molekul tidak berarti apa-apa untuk daya prediksi nyata.

---
Sebelumnya: [Kimia](chemistry.md) · Selanjutnya: [Biologi sintetik](synthetic-biology.md) · [Indeks tutorial](index.md)
