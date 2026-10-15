# API reference: Chemistry

Package `Gravicode.Biotech.Net.Chemistry` · namespace `Biotech.Net.Chemistry`. Pure C#, no RDKit dependency.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Molecule graph

| Type | Key members |
|---|---|
| `Molecule(name, atomCapacity)` | `FromSmiles(smiles, name)`, `AtomCount`, `BondCount`, `Atoms`, `Bonds`, `AtomAt(i)`, `BondAt(i)`, `AddAtom`, `AddBond`, `RemoveAtoms`, `RemoveBond`, `Clone`, `Coordinates`, `Has3DCoordinates`, `HeavyAtomCount`, `Neighbors`, `NeighborBonds`, `Degree`, `HeavyDegree`, `TotalHydrogens`, `BondBetween`, `Rings`, `KekuleOrder`, `ExplicitValence`, `TotalValence`, `Fragments()`, `TopologicalDistances()`, `AddHydrogens`, `RemoveHydrogens`, `UpdateImplicitHydrogens`, `ToSmiles(canonical, includeStereo)`, `Formula`, `Properties` |
| `Atom` | `AtomicNumber`, `FormalCharge`, `Isotope`, `ImplicitHydrogens`, `ExplicitHydrogens`, `IsAromatic`, `Chirality`, `MapNumber`, `Element`, `Symbol` |
| `Bond` | `Begin`, `End`, `Order`, `Stereo`, `IsAromatic`, `Other(atom)`, `Valence` |
| `RingInfo` (`mol.Rings`) | `AtomRings`, `BondRings`, `RingSystems`, `IsAtomInRing`, `IsBondInRing`, `AtomRingCount`, `SmallestRingSize`, `IsAtomInRingOfSize` |
| `PeriodicTable` / `Element` | `Get(symbol or Z)`, `TryGet`, `IsotopeMass`, `All`; `AverageMass`, `MonoisotopicMass`, `Valences`, `IsOrganicSubset`, `IsMetal` |

Enums: `BondOrder`, `BondStereo`, `ChiralTag`; `Point3` for coordinates.

## Parsing, writing and depiction

- `SmilesParser.Parse/TryParse` with `SmilesParserOptions`; `SmilesWriter.Write(mol, canonical, includeStereo)`; `CanonicalRanking.Compute`, `SymmetryClasses`.
- `MolFile.Parse/Read/Write`, `SdfReader.Read/ReadFile/Parse`, `SdfWriter.Write/WriteFile/WriteRecord/ToSdf`, `MolReaderOptions`.
- `Aromaticity.Perceive`, `Kekulize`.
- `Depiction.Compute2DCoordinates(mol, bondLength)`, `CountClashes`.

## Descriptors and drug-likeness

- `MolecularDescriptors.Compute(mol)` returning `DescriptorSet` (`Formula`, `MolecularWeight`, `ExactMass`, `HeavyAtoms`, `Heteroatoms`, `LipinskiHbd/Hba`, `HBondDonors/Acceptors`, `RotatableBonds`, `Tpsa`, `LogP`, `MolarRefractivity`, `RingCount`, `AromaticRingCount`, `FractionCsp3`, `FormalCharge`, `TotalAtoms`, `ToVector()`, `Names`), plus the individual methods.
- `Crippen.LogP`, `MolarRefractivity`, `Compute`; `Tpsa.Compute`, `Contribution`.
- `DrugLikeness.Lipinski`, `Veber`, `Ghose`, `Egan`, `Muegge`, `LeadLike`, each returning `RuleResult(Rule, Passed, Violations)`.
- `Qed.Compute(mol, weights)`, `Properties`, `MatchedAlerts`, `Desirability`, `MeanWeights`, `MaxWeights`.
- `MurckoScaffold.GetScaffold`, `GetScaffoldSmiles`, `GetGenericFramework`, `ScaffoldAtoms`, `Extract`; `MolecularFormula.Parse/FromMolecule` with `AverageMass`, `MonoisotopicMass`, `MassToCharge`.

## Fingerprints and similarity

- `Fingerprint.Compute(mol, type)`, `For(type)`; `FingerprintType`.
- `MorganFingerprint.Compute(mol, radius, nBits, useFeatures, useBondTypes)`, `Ecfp4`, `Fcfp4`, `ComputeCounts`, `ComputeHashedCounts`, `BitInfo`, `ConnectivityInvariants`, `FeatureInvariants`.
- `MaccsKeys.Compute`, `Definition(key)`; `PathFingerprint.Compute(mol, options)`; `AtomPairFingerprint.Compute/ComputeHashed`; `TopologicalTorsionFingerprint.Compute/ComputeHashed`; `Subgraphs.Enumerate`.
- `BitVector` (`Count`, `OnBits`, `Fold`, `ToBitString`, `FromBitString`, `ToBase64`, `FromBase64`, operators) and `CountVector`.
- `Similarity.Tanimoto/Dice/Cosine/Tversky/Sokal/Kulczynski/Russel`, `Compute(metric, …)`, `Bulk`, `BulkTanimoto`, `Matrix`; `SimilarityMetric`.
- `FingerprintIndex(fingerprinter)`: `Add`, `AddRange`, `Search(query, threshold, metric)`, `TopK(query, k, metric)`, returning `SimilarityHit(Index, Id, Score)`.

## Substructure search

`SmartsPattern.Parse/TryParse`, `IsMatch`, `Matches(mol, unique, maxMatches)`, `CountMatches`, `MatchesAtom`, `MatchingAtoms`, `AtomMap`; `Substructure.HasMatch`, `FindAll`, `FromMolecule(query)`.

---
[API reference index](index.md) · [Chemistry tutorial](../tutorials/chemistry.md)
