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

## Batasan

- Koordinat hanya 2-D. Pembangkitan konformer 3-D (gaya ETKDG) masih di [roadmap](../roadmap.md).
- Stereokimia diurai dan ditulis dalam SMILES, tetapi peringkat kanoniknya tidak dijamin sama persis dengan string SMILES kanonik RDKit. Bandingkan molekul lewat sidik jari atau lewat SMILES kanonik yang dihasilkan Biotech.Net sendiri.
- Model aromatisitas mirip RDKit. Sistem heteroaromatik yang tidak lazim bisa dipersepsi berbeda dari RDKit atau Open Babel, sehingga sidik jari dan kecocokan SMARTS ikut berubah.

---
Sebelumnya: [Biologi struktural](structural-biology.md) · Selanjutnya: [Penemuan obat](drug-discovery.md) · [Indeks tutorial](index.md)
