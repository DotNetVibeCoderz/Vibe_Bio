# Tutorial: Biologi struktural

Paket: `Gravicode.Biotech.Net.StructuralBiology`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Membaca dan menulis PDB serta mmCIF, analisis geometri dan Ramachandran, superposisi Kabsch dan TM-score, struktur sekunder gaya DSSP, SASA Shrake–Rupley, kontak dan ikatan disulfida, deteksi kantong ikatan, serta sidik jari interaksi protein–ligan.

## 1. Memuat struktur

```csharp
using Biotech.Net.StructuralBiology;

// Crambin (PDB 1CRN) disertakan untuk contoh dan pengujian
Structure crambin = SampleStructures.Crambin();
Console.WriteLine(crambin);
Chain chainA = crambin["A"];
Console.WriteLine(chainA.ToSequence());   // TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN

// Dari disk: PdbReader.Read("1abc.pdb"), MmCifReader.Read("1abc.cif")
// Menulis:   PdbWriter.ToPdb(structure), MmCifReader.ToMmCif(structure)
```

## 2. Geometri, struktur sekunder, dan permukaan

```csharp
using Biotech.Net.StructuralBiology;

var structure = SampleStructures.Crambin();
var chain = structure["A"];

double rg = Geometry.RadiusOfGyration(structure.Atoms);
Console.WriteLine($"Rg = {rg:F2} Å");

foreach (var t in Geometry.BackboneTorsions(chain).Take(3))
    Console.WriteLine($"{t.Residue.Id} phi {t.Phi:F1} psi {t.Psi:F1} {t.Region}");

DsspResult dssp = Dssp.Assign(structure);
Console.WriteLine(dssp.ToSecondaryStructureString());

var sasa = ShrakeRupley.PerResidue(structure.FirstModel);
Console.WriteLine($"residu terkubur: {sasa.Count(r => r.IsBuried)}/{sasa.Count}");

foreach (var ss in Disulfides.Find(structure))
    Console.WriteLine($"{ss.CysteineA.Id} - {ss.CysteineB.Id} {ss.Distance:F2} Å");

var contacts = ContactMap.Compute(chain, cutoff: 8.0);
Console.WriteLine($"{contacts.ContactCount} kontak Cα, RCO {contacts.RelativeContactOrder():F3}");
```

## 3. Superposisi dan TM-score

```csharp
using Biotech.Net.StructuralBiology;

var native = SampleStructures.Crambin();
var moved = SampleStructures.Crambin();
// Putar dan geser salinannya, lalu pulihkan transformasinya
moved.FirstModel.Transform(Matrix3.FromAxisAngle(new Vec3(1, 1, 0).Normalized(), 35), new Vec3(5, -3, 2));

var fit = Superposition.SuperposeCa(moved["A"], native["A"]);
Console.WriteLine($"RMSD setelah fit {fit.Rmsd:F4} Å, rotasi {fit.Rotation.RotationAngle:F1}°");

var tm = StructuralAlignment.TmScore(moved["A"], native["A"]);
Console.WriteLine($"TM-score {tm.TmScore:F3}");
```

`Superposition.Kabsch` menyelesaikan fit kuadrat terkecil dengan metode kuaternion Horn, sehingga tidak pernah menghasilkan pencerminan. Implementasinya managed. Fit yang sama juga tersedia sebagai kernel native `BioKernels.KabschRmsd` yang bekerja pada array koordinat datar.

## 4. Kantong ikatan dan interaksi

```csharp
using Biotech.Net.StructuralBiology;

var structure = SampleStructures.Crambin();
foreach (var pocket in PocketDetector.Detect(structure).Take(3))
    Console.WriteLine(pocket);

// Untuk kompleks dengan ligan (rekaman HETATM):
//   foreach (var profile in InteractionFingerprint.AnalyzeAll(complex))
//       Console.WriteLine(profile);       // ikatan hidrogen, hidrofobik, ionik, pi-stacking per residu
```

Visualisasi: `ProteinViewer.CaTrace(structure)` menggambar jejak Cα statis, `ProteinViewer.ContactMap(structure)` menggambar peta kontak, dan `ProteinViewer.To3DHtml(structure)` menghasilkan halaman 3Dmol.js interaktif.

## Batasan

- Struktur sekunder mengikuti algoritma DSSP (energi ikatan hidrogen Kabsch & Sander, heliks, ladder, dan bridge), tetapi bisa berbeda dari program acuan `mkdssp` pada kasus tepi seperti prioritas heliks-π dan putusnya rantai.
- Deteksi kantong memakai metode grid gaya LIGSITE yang menghitung kejadian protein–pelarut–protein. Metode ini tidak memprediksi druggability.
- Tidak ada prediksi struktur dari sekuens. Embedding ESM dibahas di [tutorial AI](ai.md), dan model AlphaFold bisa dimuat sebagai berkas PDB atau mmCIF.

---
Sebelumnya: [Proteomika](proteomics.md) · Selanjutnya: [Kimia](chemistry.md) · [Indeks tutorial](index.md)
