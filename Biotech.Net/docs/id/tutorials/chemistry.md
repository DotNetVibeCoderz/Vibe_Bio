# Tutorial: Kemoinformatika

Paket: `Gravicode.Biotech.Net.Chemistry`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Input dan output SMILES serta MOL/SDF, SMILES kanonik, persepsi aromatisitas dan cincin, deskriptor (logP Crippen, TPSA, Lipinski dan kawan-kawan, QED), sidik jari molekul (Morgan/ECFP, FCFP, MACCS, path, atom-pair, topological torsion), pencarian kemiripan, pencarian substruktur SMARTS, scaffold Murcko, dan penggambaran 2-D. Seluruhnya C# murni, tanpa ketergantungan pada RDKit.

## 1. Molekul dan deskriptor

```csharp
using Biotech.Net.Chemistry;

var aspirin = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin");
Console.WriteLine($"{aspirin.Name}: {aspirin.Formula}, {aspirin.AtomCount} atom berat");
Console.WriteLine($"SMILES kanonik: {aspirin.ToSmiles()}");

DescriptorSet d = MolecularDescriptors.Compute(aspirin);
Console.WriteLine($"BM {d.MolecularWeight:F2}  logP {d.LogP:F2}  TPSA {d.Tpsa:F1}  HBD {d.HBondDonors}  HBA {d.HBondAcceptors}  RotB {d.RotatableBonds}");
Console.WriteLine(DrugLikeness.Lipinski(d));
Console.WriteLine(DrugLikeness.Veber(d));
Console.WriteLine($"QED {Qed.Compute(aspirin):F3}");

Console.WriteLine($"Scaffold Murcko: {MurckoScaffold.GetScaffoldSmiles(aspirin)}");
```

## 2. Sidik jari dan kemiripan

```csharp
using Biotech.Net.Chemistry;

var library = new[]
{
    ("aspirin",     "CC(=O)Oc1ccccc1C(=O)O"),
    ("salicylic",   "OC(=O)c1ccccc1O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen",   "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("caffeine",    "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
}.Select(x => Molecule.FromSmiles(x.Item2, x.Item1)).ToList();

BitVector ecfp4 = MorganFingerprint.Ecfp4(library[0]);          // radius 2, 2048 bit
BitVector maccs = MaccsKeys.Compute(library[0]);                 // 166 kunci
Console.WriteLine($"bit ECFP4 menyala: {ecfp4.Count}, kunci MACCS menyala: {maccs.Count}");

var index = new FingerprintIndex(m => MorganFingerprint.Ecfp4(m));
index.AddRange(library);
foreach (var hit in index.TopK(library[0], k: 3))
    Console.WriteLine($"{hit.Id,-12} Tanimoto {hit.Score:F3}");

// Metrik lain: Dice, Cosine, Tversky, Sokal...
double dice = Similarity.Dice(ecfp4, MorganFingerprint.Ecfp4(library[1]));
```

## 3. Pencarian substruktur dan format berkas

```csharp
using Biotech.Net.Chemistry;

var ibuprofen = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
Console.WriteLine(Substructure.HasMatch(ibuprofen, "C(=O)[OH]"));     // asam karboksilat: True
var pattern = SmartsPattern.Parse("c1ccccc1");
Console.WriteLine($"kecocokan benzena: {pattern.CountMatches(ibuprofen)}");

// Bolak-balik MOL / SDF (koordinat 2-D dibuat saat diperlukan)
string molBlock = MolFile.Write(ibuprofen);
Molecule back = MolFile.Parse(molBlock);
string sdf = SdfWriter.ToSdf([ibuprofen, back]);
Console.WriteLine($"{SdfReader.Parse(sdf).Count} rekaman dalam SDF");

// Tata letak 2-D untuk penggambaran
Point3[] xy = Depiction.Compute2DCoordinates(ibuprofen);
```

`MoleculeViewer` di `Biotech.Net.Visualization` merender gambar 2-D ke PNG atau SVG:

![Gambar molekul](../../images/molecule.png)

## 4. Konformer 3-D dan medan gaya UFF

`ConformerGenerator` menanamkan molekul dalam 3-D dengan alur geometri jarak bergaya ETKDG (pengaturan ETKDGv3 RDKit: matriks batas, penanaman 4-D, pemeriksaan kiralitas, preferensi torsi eksperimental, suku planaritas), lalu merapikan setiap konformer dengan UFF dan membuang duplikat yang hampir sama berdasarkan RMSD yang memperhitungkan simetri. Hidrogen ditambahkan otomatis dan hasilnya dapat direproduksi untuk seed yang sama.

```csharp
using Biotech.Net.Chemistry;

var ibuprofen = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "ibuprofen");
ConformerEnsemble ensemble = ConformerGenerator.Generate(ibuprofen, 10, new ConformerOptions { RandomSeed = 42 });
Console.WriteLine($"{ensemble.Count} dari {ensemble.Attempted} konformer tersisa setelah pemangkasan RMSD");
for (int i = 0; i < ensemble.Count; i++)
    Console.WriteLine($"konformer {i}: energi UFF {ensemble.Energies[i]:F2} kkal/mol");
SdfWriter.WriteFile("ibuprofen_conformers.sdf", ensemble.ToMolecules());   // SDF 3-D dengan H eksplisit
double rmsd = ConformerAlignment.BestRmsd(ensemble.Molecule, ensemble.Conformers[0], ensemble.Conformers[1]);
Console.WriteLine($"RMSD atom berat antara konformer 0 dan 1: {rmsd:F2} Å");

// Satu konformer, lalu Universal Force Field secara langsung
Molecule phenol = ConformerGenerator.Embed(Molecule.FromSmiles("Oc1ccccc1", "phenol"));
var uff = UffForceField.Create(phenol);                  // perlu hidrogen eksplisit (Embed menambahkannya)
Console.WriteLine(string.Join(" ", uff.AtomTypes));      // O_R C_R C_R C_R C_R C_R C_R H_ ...
Console.WriteLine(uff.EnergyTerms(phenol.Coordinates!)); // ikatan, sudut, torsi, inversi, vdW
MinimizationResult min = UffForceField.Optimize(phenol); // L-BFGS, di tempat
Console.WriteLine($"minimum UFF {min.Energy:F3} kkal/mol setelah {min.Iterations} iterasi");
```

Akurasi terhadap RDKit 2026.09.1 (30 molekul, lihat `tests/Biotech.Net.Chemistry.Tests/Data/conformer_reference_gen.py`):

| Pemeriksaan | Hasil |
| --- | --- |
| Tipe atom UFF, hibridisasi, konjugasi | identik pada ke-30 molekul |
| Energi UFF pada geometri RDKit | selisih di bawah 1×10⁻⁴ kkal/mol |
| Gradien UFF | selisih di bawah 2×10⁻³ kkal/mol/Å |
| Minimum UFF RDKit | juga minimum medan gaya kami (minimisasi ulang mengubah E < 0,01 kkal/mol) |
| Panjang ikatan vs konformer ETKDGv3 RDKit | deviasi rata-rata 0,015 Å (maks 0,10 Å) |
| Sudut ikatan vs konformer ETKDGv3 RDKit | deviasi rata-rata 2,4° (maks 25°) |
| Stereo tetrahedral dan E/Z | selalu dipertahankan |
| Cakupan ansambel RDKit (rata-rata RMSD terbaik) | dalam 0,13 Å dari cakupan RDKit sendiri antar-seed (0,07 Å tanpa perapian UFF) |

Membangkitkan sepuluh konformer molekul seukuran obat butuh sekitar 0,1–0,9 detik pada satu inti; secara bawaan konformer dikerjakan paralel.

## Batasan

- Pembangkitan konformer mengikuti algoritme ETKDG RDKit tetapi tidak identik bit demi bit: batas jarak 1-5, batas 1-4 khusus makrosiklus, stereo atropisomer dan non-tetrahedral, penanaman terkekang (peta koordinat), serta peminimal "all-in-one" RDKit yang lebih baru belum diimplementasikan, dan kecocokan pola torsi yang simetris bisa memilih atom berbeda. Perapian UFF (aktif secara bawaan, `ForceFieldCleanup`) sedikit menggeser konformer dari torsi mirip kristal hasil ETKDG; matikan agar paling dekat dengan keluaran RDKit.
- Garam dan molekul multi-fragmen ditanamkan per fragmen lalu dijajarkan berjarak 4 Å (RDKit menumpuknya).
- Medan gayanya UFF, bukan MMFF94; tanpa elektrostatik (bawaan UFF di RDKit). Pusat trigonal-bipiramidal dan sebagian besar keadaan muatan logam belum diparameterisasi.
- Stereokimia diurai dan ditulis dalam SMILES, tetapi peringkat kanoniknya tidak dijamin sama persis dengan string SMILES kanonik RDKit. Bandingkan molekul lewat sidik jari atau lewat SMILES kanonik yang dihasilkan Biotech.Net sendiri.
- Model aromatisitas mirip RDKit. Sistem heteroaromatik yang tidak lazim bisa dipersepsi berbeda dari RDKit atau Open Babel, sehingga sidik jari dan kecocokan SMARTS ikut berubah.

---
Sebelumnya: [Biologi struktural](structural-biology.md) · Selanjutnya: [Penemuan obat](drug-discovery.md) · [Indeks tutorial](index.md)
