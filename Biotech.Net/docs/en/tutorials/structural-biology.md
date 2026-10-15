# Tutorial: Structural biology

Package: `Gravicode.Biotech.Net.StructuralBiology`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

PDB and mmCIF reading and writing, geometry and Ramachandran analysis, Kabsch superposition and TM-score, DSSP-style secondary structure, Shrake–Rupley SASA, contacts and disulfides, pocket detection and protein–ligand interaction fingerprints.

## 1. Load a structure

```csharp
using Biotech.Net.StructuralBiology;

// Crambin (PDB 1CRN) is embedded for examples and tests
Structure crambin = SampleStructures.Crambin();
Console.WriteLine(crambin);
Chain chainA = crambin["A"];
Console.WriteLine(chainA.ToSequence());   // TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN

// From disk: PdbReader.Read("1abc.pdb"), MmCifReader.Read("1abc.cif")
// Writing:   PdbWriter.ToPdb(structure), MmCifReader.ToMmCif(structure)
```

## 2. Geometry, secondary structure and surface

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
Console.WriteLine($"buried residues: {sasa.Count(r => r.IsBuried)}/{sasa.Count}");

foreach (var ss in Disulfides.Find(structure))
    Console.WriteLine($"{ss.CysteineA.Id} - {ss.CysteineB.Id} {ss.Distance:F2} Å");

var contacts = ContactMap.Compute(chain, cutoff: 8.0);
Console.WriteLine($"{contacts.ContactCount} Cα contacts, RCO {contacts.RelativeContactOrder():F3}");
```

## 3. Superposition and TM-score

```csharp
using Biotech.Net.StructuralBiology;

var native = SampleStructures.Crambin();
var moved = SampleStructures.Crambin();
// Rotate and translate a copy, then recover the transform
moved.FirstModel.Transform(Matrix3.FromAxisAngle(new Vec3(1, 1, 0).Normalized(), 35), new Vec3(5, -3, 2));

var fit = Superposition.SuperposeCa(moved["A"], native["A"]);
Console.WriteLine($"RMSD after fit {fit.Rmsd:F4} Å, rotation {fit.Rotation.RotationAngle:F1}°");

var tm = StructuralAlignment.TmScore(moved["A"], native["A"]);
Console.WriteLine($"TM-score {tm.TmScore:F3}");
```

`Superposition.Kabsch` solves the least-squares fit with Horn's quaternion method, so it never returns a reflection. It is managed code. The same fit is also available as the native kernel `BioKernels.KabschRmsd`, which works on flat coordinate arrays.

## 4. Pockets and interactions

```csharp
using Biotech.Net.StructuralBiology;

var structure = SampleStructures.Crambin();
foreach (var pocket in PocketDetector.Detect(structure).Take(3))
    Console.WriteLine(pocket);

// For a complex with a ligand (HETATM records):
//   foreach (var profile in InteractionFingerprint.AnalyzeAll(complex))
//       Console.WriteLine(profile);       // H-bonds, hydrophobic, ionic, pi-stacking per residue
```

Visualisation: `ProteinViewer.CaTrace(structure)` draws a static Cα trace, `ProteinViewer.ContactMap(structure)` a contact map, and `ProteinViewer.To3DHtml(structure)` an interactive 3Dmol.js page.

## Limitations

- The secondary structure follows the DSSP algorithm (Kabsch & Sander H-bond energy, helices, ladders and bridges), but it can differ from the reference `mkdssp` program at edge cases such as π-helix priority and chain breaks.
- Pocket detection is a LIGSITE-style grid method that counts protein–solvent–protein events. It does not predict druggability.
- There is no structure prediction from sequence. ESM embeddings are covered in the [AI tutorial](ai.md), and AlphaFold models can be loaded as PDB or mmCIF files.

---
Previous: [Proteomics](proteomics.md) · Next: [Chemistry](chemistry.md) · [Tutorial index](index.md)
