# Tutorial: Cheminformatics

Package: `Gravicode.Biotech.Net.Chemistry`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

SMILES and MOL/SDF input and output, canonical SMILES, aromaticity and ring perception, descriptors (Crippen logP, TPSA, Lipinski and friends, QED), fingerprints (Morgan/ECFP, FCFP, MACCS, path, atom-pair, topological torsion), similarity search, SMARTS substructure search, Murcko scaffolds and 2-D depiction. It is written in pure C#, with no RDKit dependency.

## 1. Molecules and descriptors

```csharp
using Biotech.Net.Chemistry;

var aspirin = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin");
Console.WriteLine($"{aspirin.Name}: {aspirin.Formula}, {aspirin.AtomCount} heavy atoms");
Console.WriteLine($"canonical SMILES: {aspirin.ToSmiles()}");

DescriptorSet d = MolecularDescriptors.Compute(aspirin);
Console.WriteLine($"MW {d.MolecularWeight:F2}  logP {d.LogP:F2}  TPSA {d.Tpsa:F1}  HBD {d.HBondDonors}  HBA {d.HBondAcceptors}  RotB {d.RotatableBonds}");
Console.WriteLine(DrugLikeness.Lipinski(d));
Console.WriteLine(DrugLikeness.Veber(d));
Console.WriteLine($"QED {Qed.Compute(aspirin):F3}");

Console.WriteLine($"Murcko scaffold: {MurckoScaffold.GetScaffoldSmiles(aspirin)}");
```

## 2. Fingerprints and similarity

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

BitVector ecfp4 = MorganFingerprint.Ecfp4(library[0]);          // radius 2, 2048 bits
BitVector maccs = MaccsKeys.Compute(library[0]);                 // 166 keys
Console.WriteLine($"ECFP4 bits set: {ecfp4.Count}, MACCS keys set: {maccs.Count}");

var index = new FingerprintIndex(m => MorganFingerprint.Ecfp4(m));
index.AddRange(library);
foreach (var hit in index.TopK(library[0], k: 3))
    Console.WriteLine($"{hit.Id,-12} Tanimoto {hit.Score:F3}");

// Other metrics: Dice, Cosine, Tversky, Sokal...
double dice = Similarity.Dice(ecfp4, MorganFingerprint.Ecfp4(library[1]));
```

## 3. Substructure search and file formats

```csharp
using Biotech.Net.Chemistry;

var ibuprofen = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
Console.WriteLine(Substructure.HasMatch(ibuprofen, "C(=O)[OH]"));     // carboxylic acid: True
var pattern = SmartsPattern.Parse("c1ccccc1");
Console.WriteLine($"benzene matches: {pattern.CountMatches(ibuprofen)}");

// MOL / SDF round trip (2-D coordinates are generated on demand)
string molBlock = MolFile.Write(ibuprofen);
Molecule back = MolFile.Parse(molBlock);
string sdf = SdfWriter.ToSdf([ibuprofen, back]);
Console.WriteLine($"{SdfReader.Parse(sdf).Count} records in SDF");

// 2-D layout for drawing
Point3[] xy = Depiction.Compute2DCoordinates(ibuprofen);
```

`MoleculeViewer` in `Biotech.Net.Visualization` renders a 2-D depiction to PNG or SVG:

![Molecule depiction](../../images/molecule.png)

## Limitations

- Coordinates are 2-D only. 3-D conformer generation (ETKDG-like) is on the [roadmap](../roadmap.md).
- Stereochemistry is parsed and written in SMILES, but the canonical ranking is not guaranteed to match RDKit's canonical SMILES string for string. Compare molecules by fingerprint or canonical SMILES produced by Biotech.Net itself.
- The aromaticity model is RDKit-like. Unusual heteroaromatic systems can be perceived differently from RDKit or Open Babel, which changes fingerprints and SMARTS matches.

---
Previous: [Structural biology](structural-biology.md) · Next: [Drug discovery](drug-discovery.md) · [Tutorial index](index.md)
