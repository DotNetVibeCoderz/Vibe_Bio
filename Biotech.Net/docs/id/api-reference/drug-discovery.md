# Referensi API: DrugDiscovery

Paket `Gravicode.Biotech.Net.DrugDiscovery` · namespace `Biotech.Net.DrugDiscovery` · berdiri di atas Chemistry dan StructuralBiology.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Penapisan

- `LibraryCompound(Id, Molecule)`, `LibraryCompound.FromSmiles(id, smiles)`.
- `SimilarityScreen.Run(query atau queries, library, options)` mengembalikan `ScreeningHit(Id, Rank, Score, PerFingerprint)`; `Fuse(lists, rule)`, `EnrichmentFactor(ranked, actives, fraction)`, `RocAuc(ranked, actives)`.
- `SimilarityScreenOptions`: `Fingerprints`, `Metric`, `FingerprintFusion`, `QueryFusion`, `TopK`; `FusionRule`.

## Prediksi target

- `TargetPredictor.FromReferences(references, fingerprinter, …)`, `TargetPredictor.Demo()`, `Predict(query, topK)`, `PredictSignificant(query, minZ)`, `Targets`, `LigandsOf`, `SimilarityThreshold`, `ExtremeValuePValue`.
- `ReferenceLigand(Target, LigandId, Smiles)`, `TargetPrediction(Target, RawScore, ZScore, PValue, MaxSimilarity, NearestLigand, ReferenceCount)`.
- `DemoReferenceSet.Ligands`, `Targets` — **set demonstrasi saja** (lihat bagian Batasan pada tutorial).

## ADMET dan peringatan

- `Admet.Profile(mol)` mengembalikan `AdmetProfile` (`Descriptors`, `Rules`, `EsolLogS`, `EsolSolubilityMgPerL`, `SolubilityClass`, `BloodBrainBarrier`, `GiAbsorption`, `BbbScore`, `HergRisk`, `HergReasons`, `CypFlags`, `LikelyPgpSubstrate`, `Complexity`, `Qed`, `PassedRules`, `FailedRules`).
- Estimasi individual: `AllRules`, `EsolLogS`, `EsolSolubilityMgPerL`, `AromaticProportion`, `SolubilityClass`, `BoiledEgg`, `BbbScore`, `HergRisk`, `CypInhibitionFlags`, `LikelyPgpSubstrate`, `ComplexityScore`.
- `StructuralAlerts.Pains`, `Brenk`, `All`, `Family(name)`, `Match(mol, catalogue)`, `IsClean`; `StructuralAlert(Family, Name, Smarts)`.
- `Pharmacophore.FindAtoms`, `Count`, `Counts`, `ToVector`, `ProfileSimilarity`, `Definition`; `PharmacophoreFeature`.
- `CompoundFilter`: `Where`, `WithPropertyRange`, `RequireSubstructure`, `ExcludeSubstructure`, `ExcludeAlerts`, `RequireRule`, `RequireLipinski`, `LeadLike()`, `FragmentLike()`, `Evaluate`, `Passes`, `Apply`; `FilterResult(Id, Passed, Reasons)`.

## Optimasi lead

- `Bioisosteres.Default`, `Applicable(mol, catalogue)`; `BioisostereTransform(Name, FromSmarts, ToSmiles, Rationale)`.
- `LeadOptimizer.Enumerate(lead, catalogue, depth, maxAnalogues)`, `Apply(mol, transform)`, `Optimize(lead, profile, catalogue, depth, maxAnalogues)`; `Analogue(Smiles, Molecule, Transform, ParentSmiles)`.
- `MpoScore.Compute(mol, profile)`, `Components`, `Desirability`; `MpoProfile.Oral`, `MpoProfile.Cns` dengan `MolecularWeight`, `LogP`, `Tpsa`, `HBondDonors`, `RotatableBonds`, `QedWeight`, `PropertyWeight`, `AlertPenalty`.

## QSAR

- `DescriptorMatrixBuilder.Standard`, `Fingerprints`, `Combined`, `Custom`; `DescriptorMatrix` dengan `RemoveLowVariance`, `Standardize`, `ToCsv`.
- `IRegressionModel` dengan `RidgeRegression(Alpha)` dan `RandomForestRegressor(TreeCount, MaxDepth, MinSamplesSplit, MaxFeatures, Seed, OutOfBagR2, FeatureImportances)`.
- `CrossValidation.KFold`, `LeaveOneOut`; `CrossValidationResult(R2, Rmse, Mae, FoldR2, Predictions)`; `RegressionMetrics.R2/Rmse/Mae/Pearson/Spearman`.
- `QsarModel.Build(compounds, response, fit, folds)`, `Predict(mol)`, `Model`, `Matrix`, `CrossValidation`.

## Penilaian pose

- `LigandPose.FromMolecule(mol)` (perlu koordinat 3-D), `FromResidue(residue)`, `Translated`, `Rotated`, `Center`, `HeavyAtoms`; `PoseAtom(Position, Element, IsHydrophobic, IsDonor, IsAcceptor)`.
- `VinaScore.Score(pose, receptor atau structure, rotatableBonds, options)`, `Score(molecule, receptor, options)`, `Optimize(...)`, `ScoreInPocket(pose, receptor, pocket, …)`, `Rank(poses, receptor, …)`, `Hydrophobic`, `HydrogenBond`.
- `VinaScoreTerms(Total, Gauss1, Gauss2, Repulsion, Hydrophobic, HydrogenBond, RotatableBonds, ContactCount)` dengan `Intermolecular` dan `Affinity`; `VinaScoreOptions` untuk bobot dan `CutOff`.

Penilaian mengevaluasi pose yang diberikan dan menyediakan pencarian lokal benda tegar; tidak ada pembangkitan konformer maupun docking fleksibel.

---
[Indeks referensi API](index.md) · [Tutorial penemuan obat](../tutorials/drug-discovery.md)
