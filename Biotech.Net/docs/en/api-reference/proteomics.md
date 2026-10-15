# API reference: Proteomics

Package `Gravicode.Biotech.Net.Proteomics` · namespaces `Biotech.Net.Proteomics`, `.MassSpec`, `.Classification`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Protein analysis

- `ProteinAnalysis(residues / Sequence / Protein)`: `MolecularWeight(monoisotopic)`, `IsoelectricPoint(pKaSet)`, `NetCharge(pH, pKaSet)`, `ChargeCurve`, `ExtinctionCoefficient()`, `Absorbance01Percent`, `InstabilityIndex()`, `IsStable`, `AliphaticIndex()`, `Gravy(scale)`, `Aromaticity()`, `SecondaryStructureFraction()`, `Flexibility()`, `Profile(scale, window, edge)`, `AminoAcidCounts()`, `AminoAcidFractions()`, `HalfLife(system)`, `Report()`. `PKaSet` selects the pKa table.
- `AminoAcidScale`: `KyteDoolittle`, `HoppWoods`, `Eisenberg`, `Flexibility`, `EminiSurface`, `Janin`; `Average`, `Profile`.
- `AminoAcids.Standard`, `IndexOf`, `IsStandard`; `NEndRule.HalfLife(nTerminal, system)`.
- `ChouFasman.Predict`, `TransmembranePredictor.Predict(residues, window, threshold)` returning `TransmembraneSegment`, `SignalPeptidePredictor.Predict` returning `SignalPeptidePrediction`.

## Mass spectrometry (`.MassSpec`)

| Area | Types and members |
|---|---|
| Masses | `MassConstants` (`Proton`, `WaterMonoisotopic`, `C13Spacing`, `ToMz`, `ToNeutralMass`, `PpmError`), `ResidueMasses.Monoisotopic/Average/PeptideMass`, `ChemicalFormula.Parse/FromPeptide/Averagine`, `Isotopes.Of` |
| Digestion | `Enzyme.Trypsin/TrypsinP/LysC/ArgC/GluC/AspN/Chymotrypsin/Pepsin`, `All`, `Parse`, `CleavageSites`; `Digestion.Digest(...)`, `DigestionOptions(MissedCleavages, MinLength, MaxLength, ClipInitiatorMethionine)`, `DigestedPeptide` |
| Modifications | `Modification.Carbamidomethyl/Oxidation/Phospho/Acetyl/Deamidation/Methyl/GlyGly/PyroGluQ/Tmt6plex`, `ModificationPosition`; `Peptide(sequence, residueModifications, nTerm, cTerm)` with `MonoisotopicMass`, `AverageMass`, `Mz(charge)`, `ModifiedSequence`; `PeptideModifier.Expand` |
| Fragments | `IonType`, `FragmentIon(Type, Number, Charge, Mz)`, `FragmentCalculator.Ladder/Fragments/SeriesOffset` |
| Isotopes | `IsotopeDistribution.Compute/FromPeptide/Averagine/ToMz`, `IsotopePeak`, `IsotopeNormalization` |
| Spectra | `Spectrum(id, mz, intensity)` with `MsLevel`, `RetentionTime`, `PrecursorMz`, `PrecursorCharge`, `TotalIonCurrent`, `FindPeak`, `TopPeaks`; `MgfReader.Read/ReadFile`, `MzMLReader.Read/ReadFile/Decode/Encode` |
| Search | `PeptideSpectrumMatcher(proteins, parameters)` with `Search(spectra)`, `Match(spectrum)`, `IndexedPeptides`; `SearchParameters` (`Enzyme`, `Digestion`, `PrecursorTolerance`, `FragmentTolerance`, `FixedModifications`, `VariableModifications`, `MaxVariableModifications`, `IonTypes`, `MaxFragmentCharge`, `DefaultCharges`, `GenerateDecoys`, `DecoyPrefix`, `TopPeaks`); `Tolerance.Ppm/Da`; `PeptideSpectrumMatch` (hyperscore, matched ions, `PrecursorErrorPpm`, `QValue`); `TargetDecoy.AssignQValues/Filter` |
| Inference and quantification | `ProteinInference.Parsimony`, `ProteinGroup`; `LabelFreeQuantification.Quantify`, `Nsaf`, `ProteinQuantity(ProteinId, Length, SpectralCount, Nsaf, Intensity, ObservablePeptides, Ibaq)` |

## Classification (`.Classification`)

- `ProteinFeatures.AminoAcidComposition`, `DipeptideComposition`, `KmerSpectrum`, `Physicochemical`, `Combined`, `Extract`.
- `SequenceSimilarity.KmerJaccard`, `KmerCosine`, `CosineMatrix`.
- `IProteinClassifier` with `KNearestNeighborsClassifier(k, metric)` and `LogisticRegressionClassifier(l2, learningRate, iterations)`.
- `ProteinFamilyClassifier(model, featureExtractor)`: `Fit`, `Predict`, `PredictProbabilities`, `Classes`, `CrossValidate`.
- `Classifiers.StratifiedFolds`, `CrossValidate`; `CrossValidationResult(FoldAccuracies, Accuracy, ConfusionMatrix, Classes)` with `MacroF1`.
- `ProteinFamilyDemo.Generate(perFamily, seed, noise)` — synthetic families for examples and tests only.

---
[API reference index](index.md) · [Proteomics tutorial](../tutorials/proteomics.md)
