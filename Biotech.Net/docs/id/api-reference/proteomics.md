# Referensi API: Proteomics

Paket `Gravicode.Biotech.Net.Proteomics` · namespace `Biotech.Net.Proteomics`, `.MassSpec`, `.Classification`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Analisis protein

- `ProteinAnalysis(residues / Sequence / Protein)`: `MolecularWeight(monoisotopic)`, `IsoelectricPoint(pKaSet)`, `NetCharge(pH, pKaSet)`, `ChargeCurve`, `ExtinctionCoefficient()`, `Absorbance01Percent`, `InstabilityIndex()`, `IsStable`, `AliphaticIndex()`, `Gravy(scale)`, `Aromaticity()`, `SecondaryStructureFraction()`, `Flexibility()`, `Profile(scale, window, edge)`, `AminoAcidCounts()`, `AminoAcidFractions()`, `HalfLife(system)`, `Report()`. `PKaSet` memilih tabel pKa.
- `AminoAcidScale`: `KyteDoolittle`, `HoppWoods`, `Eisenberg`, `Flexibility`, `EminiSurface`, `Janin`; `Average`, `Profile`.
- `AminoAcids.Standard`, `IndexOf`, `IsStandard`; `NEndRule.HalfLife(nTerminal, system)`.
- `ChouFasman.Predict`, `TransmembranePredictor.Predict(residues, window, threshold)` yang mengembalikan `TransmembraneSegment`, `SignalPeptidePredictor.Predict` yang mengembalikan `SignalPeptidePrediction`.

## Spektrometri massa (`.MassSpec`)

| Bidang | Tipe dan anggota |
|---|---|
| Massa | `MassConstants` (`Proton`, `WaterMonoisotopic`, `C13Spacing`, `ToMz`, `ToNeutralMass`, `PpmError`), `ResidueMasses.Monoisotopic/Average/PeptideMass`, `ChemicalFormula.Parse/FromPeptide/Averagine`, `Isotopes.Of` |
| Digesti | `Enzyme.Trypsin/TrypsinP/LysC/ArgC/GluC/AspN/Chymotrypsin/Pepsin`, `All`, `Parse`, `CleavageSites`; `Digestion.Digest(...)`, `DigestionOptions(MissedCleavages, MinLength, MaxLength, ClipInitiatorMethionine)`, `DigestedPeptide` |
| Modifikasi | `Modification.Carbamidomethyl/Oxidation/Phospho/Acetyl/Deamidation/Methyl/GlyGly/PyroGluQ/Tmt6plex`, `ModificationPosition`; `Peptide(sequence, residueModifications, nTerm, cTerm)` dengan `MonoisotopicMass`, `AverageMass`, `Mz(charge)`, `ModifiedSequence`; `PeptideModifier.Expand` |
| Fragmen | `IonType`, `FragmentIon(Type, Number, Charge, Mz)`, `FragmentCalculator.Ladder/Fragments/SeriesOffset` |
| Isotop | `IsotopeDistribution.Compute/FromPeptide/Averagine/ToMz`, `IsotopePeak`, `IsotopeNormalization` |
| Spektrum | `Spectrum(id, mz, intensity)` dengan `MsLevel`, `RetentionTime`, `PrecursorMz`, `PrecursorCharge`, `TotalIonCurrent`, `FindPeak`, `TopPeaks`; `MgfReader.Read/ReadFile`, `MzMLReader.Read/ReadFile/Decode/Encode` |
| Pencarian | `PeptideSpectrumMatcher(proteins, parameters)` dengan `Search(spectra)`, `Match(spectrum)`, `IndexedPeptides`; `SearchParameters` (`Enzyme`, `Digestion`, `PrecursorTolerance`, `FragmentTolerance`, `FixedModifications`, `VariableModifications`, `MaxVariableModifications`, `IonTypes`, `MaxFragmentCharge`, `DefaultCharges`, `GenerateDecoys`, `DecoyPrefix`, `TopPeaks`); `Tolerance.Ppm/Da`; `PeptideSpectrumMatch` (hyperscore, ion tercocok, `PrecursorErrorPpm`, `QValue`); `TargetDecoy.AssignQValues/Filter` |
| Inferensi dan kuantifikasi | `ProteinInference.Parsimony`, `ProteinGroup`; `LabelFreeQuantification.Quantify`, `Nsaf`, `ProteinQuantity(ProteinId, Length, SpectralCount, Nsaf, Intensity, ObservablePeptides, Ibaq)` |

## Klasifikasi (`.Classification`)

- `ProteinFeatures.AminoAcidComposition`, `DipeptideComposition`, `KmerSpectrum`, `Physicochemical`, `Combined`, `Extract`.
- `SequenceSimilarity.KmerJaccard`, `KmerCosine`, `CosineMatrix`.
- `IProteinClassifier` dengan `KNearestNeighborsClassifier(k, metric)` dan `LogisticRegressionClassifier(l2, learningRate, iterations)`.
- `ProteinFamilyClassifier(model, featureExtractor)`: `Fit`, `Predict`, `PredictProbabilities`, `Classes`, `CrossValidate`.
- `Classifiers.StratifiedFolds`, `CrossValidate`; `CrossValidationResult(FoldAccuracies, Accuracy, ConfusionMatrix, Classes)` dengan `MacroF1`.
- `ProteinFamilyDemo.Generate(perFamily, seed, noise)` — famili sintetis khusus contoh dan pengujian.

---
[Indeks referensi API](index.md) · [Tutorial proteomika](../tutorials/proteomics.md)
