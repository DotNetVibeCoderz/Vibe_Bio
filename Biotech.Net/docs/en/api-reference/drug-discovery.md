# API reference: Drug discovery

Package `Gravicode.Biotech.Net.DrugDiscovery` · namespace `Biotech.Net.DrugDiscovery` · builds on Chemistry and StructuralBiology.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Screening

- `LibraryCompound(Id, Molecule)`, `LibraryCompound.FromSmiles(id, smiles)`.
- `SimilarityScreen.Run(query or queries, library, options)` returning `ScreeningHit(Id, Rank, Score, PerFingerprint)`; `Fuse(lists, rule)`, `EnrichmentFactor(ranked, actives, fraction)`, `RocAuc(ranked, actives)`.
- `SimilarityScreenOptions`: `Fingerprints`, `Metric`, `FingerprintFusion`, `QueryFusion`, `TopK`; `FusionRule`.

## Target prediction

- `TargetPredictor.FromReferences(references, fingerprinter, …)`, `TargetPredictor.Demo()`, `Predict(query, topK)`, `PredictSignificant(query, minZ)`, `Targets`, `LigandsOf`, `SimilarityThreshold`, `ExtremeValuePValue`.
- `ReferenceLigand(Target, LigandId, Smiles)`, `TargetPrediction(Target, RawScore, ZScore, PValue, MaxSimilarity, NearestLigand, ReferenceCount)`.
- `DemoReferenceSet.Ligands`, `Targets` — **a demonstration set only** (see the tutorial's limitations).

## ADMET and alerts

- `Admet.Profile(mol)` returning `AdmetProfile` (`Descriptors`, `Rules`, `EsolLogS`, `EsolSolubilityMgPerL`, `SolubilityClass`, `BloodBrainBarrier`, `GiAbsorption`, `BbbScore`, `HergRisk`, `HergReasons`, `CypFlags`, `LikelyPgpSubstrate`, `Complexity`, `Qed`, `PassedRules`, `FailedRules`).
- Individual estimates: `AllRules`, `EsolLogS`, `EsolSolubilityMgPerL`, `AromaticProportion`, `SolubilityClass`, `BoiledEgg`, `BbbScore`, `HergRisk`, `CypInhibitionFlags`, `LikelyPgpSubstrate`, `ComplexityScore`.
- `StructuralAlerts.Pains`, `Brenk`, `All`, `Family(name)`, `Match(mol, catalogue)`, `IsClean`; `StructuralAlert(Family, Name, Smarts)`.
- `Pharmacophore.FindAtoms`, `Count`, `Counts`, `ToVector`, `ProfileSimilarity`, `Definition`; `PharmacophoreFeature`.
- `CompoundFilter`: `Where`, `WithPropertyRange`, `RequireSubstructure`, `ExcludeSubstructure`, `ExcludeAlerts`, `RequireRule`, `RequireLipinski`, `LeadLike()`, `FragmentLike()`, `Evaluate`, `Passes`, `Apply`; `FilterResult(Id, Passed, Reasons)`.

## Lead optimisation

- `Bioisosteres.Default`, `Applicable(mol, catalogue)`; `BioisostereTransform(Name, FromSmarts, ToSmiles, Rationale)`.
- `LeadOptimizer.Enumerate(lead, catalogue, depth, maxAnalogues)`, `Apply(mol, transform)`, `Optimize(lead, profile, catalogue, depth, maxAnalogues)`; `Analogue(Smiles, Molecule, Transform, ParentSmiles)`.
- `MpoScore.Compute(mol, profile)`, `Components`, `Desirability`; `MpoProfile.Oral`, `MpoProfile.Cns` with `MolecularWeight`, `LogP`, `Tpsa`, `HBondDonors`, `RotatableBonds`, `QedWeight`, `PropertyWeight`, `AlertPenalty`.

## QSAR

- `DescriptorMatrixBuilder.Standard`, `Fingerprints`, `Combined`, `Custom`; `DescriptorMatrix` with `RemoveLowVariance`, `Standardize`, `ToCsv`.
- `IRegressionModel` with `RidgeRegression(Alpha)` and `RandomForestRegressor(TreeCount, MaxDepth, MinSamplesSplit, MaxFeatures, Seed, OutOfBagR2, FeatureImportances)`.
- `CrossValidation.KFold`, `LeaveOneOut`; `CrossValidationResult(R2, Rmse, Mae, FoldR2, Predictions)`; `RegressionMetrics.R2/Rmse/Mae/Pearson/Spearman`.
- `QsarModel.Build(compounds, response, fit, folds)`, `Predict(mol)`, `Model`, `Matrix`, `CrossValidation`.

## Pose scoring

- `LigandPose.FromMolecule(mol)` (needs 3-D coordinates), `FromResidue(residue)`, `Translated`, `Rotated`, `Center`, `HeavyAtoms`; `PoseAtom(Position, Element, IsHydrophobic, IsDonor, IsAcceptor)`.
- `VinaScore.Score(pose, receptor or structure, rotatableBonds, options)`, `Score(molecule, receptor, options)`, `Optimize(...)`, `ScoreInPocket(pose, receptor, pocket, …)`, `Rank(poses, receptor, …)`, `Hydrophobic`, `HydrogenBond`.
- `VinaScoreTerms(Total, Gauss1, Gauss2, Repulsion, Hydrophobic, HydrogenBond, RotatableBonds, ContactCount)` with `Intermolecular` and `Affinity`; `VinaScoreOptions` for the weights and `CutOff`.

Scoring evaluates a given pose and offers a rigid-body local search; there is no conformer generation or flexible docking.

---
[API reference index](index.md) · [Drug discovery tutorial](../tutorials/drug-discovery.md)
