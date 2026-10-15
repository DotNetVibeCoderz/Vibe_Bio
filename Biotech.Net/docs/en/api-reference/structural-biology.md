# API reference: Structural biology

Package `Gravicode.Biotech.Net.StructuralBiology` · namespace `Biotech.Net.StructuralBiology`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Object model

`Structure` → `Model` → `Chain` → `Residue` → `Atom`.

| Type | Key members |
|---|---|
| `Structure(id)` | `Title`, `Classification`, `Cell`, `Models`, `FirstModel`, `Chains`, `Atoms`, `Residues`, `Ligands`, `this[chainId]`, `Helices`, `Sheets`, `Disulfides`, `Connections`, `AtomCount` |
| `Model(number)` | `Chains`, `this[chainId]`, `GetChain`, `Residues`, `Atoms`, `GetAtoms()`, `GetCoordinates()`, `Ligands`, `Transform(rotation, translation)` |
| `Chain(id)` | `Residues`, `AminoAcids()`, `FindResidue(number, insertionCode)`, `Atoms`, `ToSequence(id)` |
| `Residue(name, number, insertionCode, isHetero)` | `Atoms`, `GetAtom(name)`, `this[atomName]`, `HasAtom`, `IsAminoAcid`, `IsWater`, `IsLigand`, `OneLetterCode`, `Id` |
| `Atom(name, element, position)` | `Position`, `Serial`, `AltLoc`, `Occupancy`, `BFactor`, `Charge`, `IsHetero`, `Alternates`, `Mass`, `VdwRadius`, `IsBackbone`, `ResidueName`, `ChainId` |

Maths: `Vec3` (`Dot`, `Cross`, `Distance`, `AngleBetween`, `Normalized`) and `Matrix3` (`Identity`, `FromAxisAngle`, `FromQuaternion`, `Transpose`, `Determinant`, `RotationAngle`).

## I/O

- `PdbReader.Read(path)`, `Parse(text, id)`, `Read(TextReader, id)`.
- `PdbWriter.ToPdb(structure, options)`, `Write(...)`, `FormatAtom`; `PdbWriteOptions` (`IncludeHeader`, `IncludeConnections`, `IncludeAlternates`, `Renumber`).
- `MmCifReader.Read`, `Parse`, `ParseCif`, `ToMmCif`; `CifBlock`, `CifLoop`.
- `SampleStructures.Crambin()` (PDB 1CRN, embedded), `CrambinPdb`, `AlanineHelix(length)`, `AlanineStrand(length)`, `PeptideBuilder.Build`.
- `ResidueNames.ToOneLetter`, `IsAminoAcid`, `IsNucleotide`, `IsWater`, `MaxAsaTien2013`; `ElementData.Mass`, `VdwRadius`, `CovalentRadius`, `GuessFromAtomName`.

## Geometry and secondary structure

- `Geometry.Distance`, `Angle`, `Dihedral`, `Centroid`, `CenterOfMass`, `TotalMass`, `RadiusOfGyration`, `BoundingBox`, `BackboneTorsions(chain)`, `Bonded`, `Chi1`.
- `Ramachandran.Classify(phi, psi)`, `Summarize(chain)`; `BackboneTorsion`, `RamachandranRegion`.
- `Dssp.Assign(structure or model, options)` returning `DsspResult` with `Residues`, `HydrogenBonds`, `ToSecondaryStructureString()`, `ToThreeStateString()`, `Composition()`; `SecondaryStructure`, `DsspResidue`, `BackboneHBond`, `Dssp.HBondEnergy`.
- `ShrakeRupley.Compute`, `Total`, `PerResidue(model, probe, points)`, `SpherePoints`; `ResidueSasa(Residue, Absolute, Relative)` with `IsBuried`.

## Contacts, superposition, pockets, interactions

- `ContactMap.Compute(chain or residues, cutoff, mode)`: `Pairs`, `ContactCount`, `RelativeContactOrder`, `ToMatrix`, `DistanceMatrix(chain)`; `ContactMode`.
- `Disulfides.Find(model or structure, minDistance, maxDistance)` returning `Disulfide(CysteineA, CysteineB, Distance)` with `Chi3`.
- `Superposition.Kabsch(mobile, target, weights)`, `Rmsd`, `FittedRmsd`, `SuperposeCa(mobile, target)`, `CaTrace(chain)`; `SuperpositionResult(Rotation, Translation, Rmsd)` with `Apply(...)`.
- `StructuralAlignment.TmScore(...)`, `Align(...)`, `D0(length)`; `TmScoreResult(TmScore, Rmsd, D0, NormalizationLength, AlignedPairs, Transform)` with `Coverage`.
- `PocketDetector.Detect(structure, model or atoms, options)` returning `Pocket(Rank, Points, Center, Volume, MeanBuriedness, LiningResidues)`; `PocketOptions` (`GridSpacing`, `MinBuriedness`, `ScanDistance`, `MinPocketPoints`, `LiningDistance`, `IgnoreHetero`).
- `InteractionFingerprint.AnalyzeAll(structure, options)`, `Analyze(model, ligand, options)`; `InteractionProfile` (`Interactions`, `Residues`, `OfType`, `Has`, `ToFingerprint`), `Interaction(Type, LigandAtom, ProteinAtom, Distance, Detail)`, `InteractionType`, `InteractionOptions`, `InteractionBits.Tanimoto`.
- `NeighborSearch(points, cellSize)` / `ForAtoms(atoms)`: `Query(center, radius)`, `Any`, `Nearest`, `AllPairs(radius)`.

---
[API reference index](index.md) · [Structural biology tutorial](../tutorials/structural-biology.md)
