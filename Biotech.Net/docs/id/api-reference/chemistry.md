# Referensi API: Chemistry

Paket `Gravicode.Biotech.Net.Chemistry` · namespace `Biotech.Net.Chemistry`. C# murni, tanpa ketergantungan RDKit.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Graf molekul

| Tipe | Anggota utama |
|---|---|
| `Molecule(name, atomCapacity)` | `FromSmiles(smiles, name)`, `AtomCount`, `BondCount`, `Atoms`, `Bonds`, `AtomAt(i)`, `BondAt(i)`, `AddAtom`, `AddBond`, `RemoveAtoms`, `RemoveBond`, `Clone`, `Coordinates`, `Has3DCoordinates`, `HeavyAtomCount`, `Neighbors`, `NeighborBonds`, `Degree`, `HeavyDegree`, `TotalHydrogens`, `BondBetween`, `Rings`, `KekuleOrder`, `ExplicitValence`, `TotalValence`, `Fragments()`, `TopologicalDistances()`, `AddHydrogens`, `RemoveHydrogens`, `UpdateImplicitHydrogens`, `ToSmiles(canonical, includeStereo)`, `Formula`, `Properties` |
| `Atom` | `AtomicNumber`, `FormalCharge`, `Isotope`, `ImplicitHydrogens`, `ExplicitHydrogens`, `IsAromatic`, `Chirality`, `MapNumber`, `Element`, `Symbol` |
| `Bond` | `Begin`, `End`, `Order`, `Stereo`, `IsAromatic`, `Other(atom)`, `Valence` |
| `RingInfo` (`mol.Rings`) | `AtomRings`, `BondRings`, `RingSystems`, `IsAtomInRing`, `IsBondInRing`, `AtomRingCount`, `SmallestRingSize`, `IsAtomInRingOfSize` |
| `PeriodicTable` / `Element` | `Get(simbol atau Z)`, `TryGet`, `IsotopeMass`, `All`; `AverageMass`, `MonoisotopicMass`, `Valences`, `IsOrganicSubset`, `IsMetal` |

Enum: `BondOrder`, `BondStereo`, `ChiralTag`; `Point3` untuk koordinat.

## Parsing, penulisan, dan penggambaran

- `SmilesParser.Parse/TryParse` dengan `SmilesParserOptions`; `SmilesWriter.Write(mol, canonical, includeStereo)`; `CanonicalRanking.Compute`, `SymmetryClasses`.
- `MolFile.Parse/Read/Write`, `SdfReader.Read/ReadFile/Parse`, `SdfWriter.Write/WriteFile/WriteRecord/ToSdf`, `MolReaderOptions`.
- `Aromaticity.Perceive`, `Kekulize`.
- `Depiction.Compute2DCoordinates(mol, bondLength)`, `CountClashes`.

## Deskriptor dan kemiripan obat

- `MolecularDescriptors.Compute(mol)` mengembalikan `DescriptorSet` (`Formula`, `MolecularWeight`, `ExactMass`, `HeavyAtoms`, `Heteroatoms`, `LipinskiHbd/Hba`, `HBondDonors/Acceptors`, `RotatableBonds`, `Tpsa`, `LogP`, `MolarRefractivity`, `RingCount`, `AromaticRingCount`, `FractionCsp3`, `FormalCharge`, `TotalAtoms`, `ToVector()`, `Names`), beserta metode individualnya.
- `Crippen.LogP`, `MolarRefractivity`, `Compute`; `Tpsa.Compute`, `Contribution`.
- `DrugLikeness.Lipinski`, `Veber`, `Ghose`, `Egan`, `Muegge`, `LeadLike`, masing-masing mengembalikan `RuleResult(Rule, Passed, Violations)`.
- `Qed.Compute(mol, weights)`, `Properties`, `MatchedAlerts`, `Desirability`, `MeanWeights`, `MaxWeights`.
- `MurckoScaffold.GetScaffold`, `GetScaffoldSmiles`, `GetGenericFramework`, `ScaffoldAtoms`, `Extract`; `MolecularFormula.Parse/FromMolecule` dengan `AverageMass`, `MonoisotopicMass`, `MassToCharge`.

## Sidik jari dan kemiripan

- `Fingerprint.Compute(mol, type)`, `For(type)`; `FingerprintType`.
- `MorganFingerprint.Compute(mol, radius, nBits, useFeatures, useBondTypes)`, `Ecfp4`, `Fcfp4`, `ComputeCounts`, `ComputeHashedCounts`, `BitInfo`, `ConnectivityInvariants`, `FeatureInvariants`.
- `MaccsKeys.Compute`, `Definition(key)`; `PathFingerprint.Compute(mol, options)`; `AtomPairFingerprint.Compute/ComputeHashed`; `TopologicalTorsionFingerprint.Compute/ComputeHashed`; `Subgraphs.Enumerate`.
- `BitVector` (`Count`, `OnBits`, `Fold`, `ToBitString`, `FromBitString`, `ToBase64`, `FromBase64`, operator) dan `CountVector`.
- `Similarity.Tanimoto/Dice/Cosine/Tversky/Sokal/Kulczynski/Russel`, `Compute(metric, …)`, `Bulk`, `BulkTanimoto`, `Matrix`; `SimilarityMetric`.
- `FingerprintIndex(fingerprinter)`: `Add`, `AddRange`, `Search(query, threshold, metric)`, `TopK(query, k, metric)`, mengembalikan `SimilarityHit(Index, Id, Score)`.

## 3-D: konformer, medan gaya, penjajaran

- `ConformerGenerator.Generate(mol, count, options)` yang mengembalikan `ConformerEnsemble` (`Molecule` dengan H eksplisit, `Conformers`, `Energies`, `Attempted`, `Count`, `ToMolecule(i)`, `ToMolecules()`); `ConformerGenerator.Embed(mol, options)` untuk satu molekul 3-D.
- `ConformerOptions`: `RandomSeed`, `MaxAttempts`, `UseExperimentalTorsions`, `UseBasicKnowledge`, `UseMacrocycleTorsions`, `ForceTransAmides`, `EnforceChirality`, `PruneRmsThreshold`, `ForceFieldCleanup`, `CleanupMaxIterations`, `MaxDegreeOfParallelism`, `Default` (mirip ETKDGv3).
- `UffForceField.Create/TryCreate(mol)`, `AssignAtomTypes(mol)`, `AtomTypes`, `TermCounts`, `Energy(coords)`, `Energy(coords, gradient)`, `EnergyTerms(coords)` → `UffEnergyTerms(BondStretch, AngleBend, Torsion, Inversion, VanDerWaals, Total)`, `Minimize(coords, maxIterations, gradientTolerance)`, `Optimize(mol)` → `MinimizationResult(Energy, Iterations, Converged)`; `UffParameters.TryGet(label, out UffAtomParameters)`, `Labels`.
- `AtomHybridization.Compute(mol)`, `ConjugatedBonds(mol)`, `HasConjugatedBond`; enum `Hybridization` (aturan RDKit).
- `ConformerAlignment.Rmsd(a, b, align)`, `Superpose(reference, mobile)`, `BestRmsd(mol, a, b, align, maxAutomorphisms)` (RMSD atom berat yang memperhitungkan simetri), `HeavyAtomAutomorphisms(mol)`.

## Pencarian substruktur

`SmartsPattern.Parse/TryParse`, `IsMatch`, `Matches(mol, unique, maxMatches)`, `CountMatches`, `MatchesAtom`, `MatchingAtoms`, `AtomMap`; `Substructure.HasMatch`, `FindAll`, `FromMolecule(query)`. Selain SMARTS Daylight, ekstensi RDKit `^n` (hibridisasi) dan `r{a-b}` (rentang ukuran cincin) juga dipahami.

---
[Indeks referensi API](index.md) · [Tutorial kemoinformatika](../tutorials/chemistry.md)
