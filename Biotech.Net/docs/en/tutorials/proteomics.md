# Tutorial: Proteomics

Package: `Gravicode.Biotech.Net.Proteomics`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

ProtParam-style protein analysis, sequence-based structure predictors, enzymatic digestion, fragment ions and isotope patterns, peptide–spectrum matching with target–decoy FDR, protein inference, label-free quantification and protein family classification.

## 1. Protein properties (ProtParam)

```csharp
using Biotech.Net.Proteomics;

// Human insulin B chain
var analysis = new ProteinAnalysis("FVNQHLCGSHLVEALYLVCGERGFFYTPKT");
Console.WriteLine($"MW        {analysis.MolecularWeight():F2} Da");
Console.WriteLine($"pI        {analysis.IsoelectricPoint():F2}");
Console.WriteLine($"GRAVY     {analysis.Gravy():F3}");
Console.WriteLine($"Instab.   {analysis.InstabilityIndex():F2} (stable: {analysis.IsStable})");
Console.WriteLine($"Aliphatic {analysis.AliphaticIndex():F2}");
Console.WriteLine($"Ext. coef {analysis.ExtinctionCoefficient().Reduced} M^-1 cm^-1");
Console.WriteLine($"Half-life {analysis.HalfLife().Text}");

// Hydropathy profile (Kyte-Doolittle, window 9) and a full text report
double[] profile = analysis.Profile(AminoAcidScale.KyteDoolittle, window: 9);
Console.WriteLine(analysis.Report());

// Sequence-based predictors
Console.WriteLine(ChouFasman.Predict("FVNQHLCGSHLVEALYLVCGERGFFYTPKT"));
var tm = TransmembranePredictor.Predict("MKTLLILAVVAAALAGSLLLVVFLLAGIAALLGYKRRSS");
var signal = SignalPeptidePredictor.Predict("MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQ");
Console.WriteLine($"{tm.Count} TM segment(s); signal peptide: {signal.HasSignalPeptide} (cleavage after {signal.CleavageSite})");
```

`MolecularWeight` and `IsoelectricPoint` go through the native protein kernel, and the `Protein` type in Core exposes the same values.

## 2. Digestion, fragments and isotopes

```csharp
using Biotech.Net.Proteomics.MassSpec;

var peptides = Digestion.Digest("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF",
    Enzyme.Trypsin, new DigestionOptions(MissedCleavages: 1, MinLength: 5));
foreach (var p in peptides.Take(4))
    Console.WriteLine($"{p.Sequence,-20} {p.MonoisotopicMass:F4} Da  missed={p.MissedCleavages}");

var peptide = new Peptide("SAMPLER");
Console.WriteLine($"[M+2H]2+ = {peptide.Mz(2):F4}");
foreach (var ion in FragmentCalculator.Fragments(peptide, [IonType.B, IonType.Y]).Take(6))
    Console.WriteLine(ion);

foreach (var peak in IsotopeDistribution.FromPeptide("SAMPLER", maxPeaks: 4))
    Console.WriteLine($"{peak.Mass:F4}  {peak.Abundance:F3}");
```

## 3. Peptide–spectrum matching with target–decoy FDR

```csharp
using Biotech.Net.Core;
using Biotech.Net.Proteomics.MassSpec;

var proteins = new[]
{
    Protein.FromString("P1", "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFEDHVK"),
    Protein.FromString("P2", "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFELRRQACGRTRAPPPPPLRSGC"),
};

// Build a synthetic MS/MS spectrum for DLGEEHFK (2+) from its b/y ions
var target = new Peptide("DLGEEHFK");
var ions = FragmentCalculator.Fragments(target, [IonType.B, IonType.Y]);
var spectrum = new Spectrum("scan1", ions.Select(i => i.Mz).ToArray(), ions.Select(_ => 100.0).ToArray())
{
    PrecursorMz = target.Mz(2),
    PrecursorCharge = 2,
};

var matcher = new PeptideSpectrumMatcher(proteins, new SearchParameters { VariableModifications = [] });
var psms = matcher.Search([spectrum]);
TargetDecoy.AssignQValues(psms);
foreach (var psm in psms) Console.WriteLine(psm);

// Real data: MgfReader.ReadFile("run.mgf") or MzMLReader.ReadFile("run.mzML")
```

`ProteinInference.Parsimony(psms)` groups proteins by shared peptides, and `LabelFreeQuantification.Quantify(psms, proteins)` reports spectral counts, NSAF and iBAQ.

## 4. Protein family classification

```csharp
using Biotech.Net.Proteomics.Classification;

// A synthetic four-family data set for demonstrations and tests
var examples = ProteinFamilyDemo.Generate(perFamily: 20, seed: 1);
var cv = ProteinFamilyClassifier.CrossValidate(examples, () => new LogisticRegressionClassifier(), folds: 5);
Console.WriteLine($"5-fold accuracy {cv.Accuracy:P1}, macro-F1 {cv.MacroF1:F3}");

var classifier = new ProteinFamilyClassifier(new KNearestNeighborsClassifier(k: 3)).Fit(examples);
Console.WriteLine(classifier.Predict(examples[0].Protein));
```

## Limitations

- `ProteinFamilyDemo` generates synthetic families for demonstrations. Train on real labelled data (for example Pfam families) before drawing biological conclusions.
- Chou–Fasman, the transmembrane predictor and the signal-peptide predictor are classical heuristics, far less accurate than SignalP, TMHMM/DeepTMHMM or AlphaFold-derived annotations. For protein language-model features, see the [AI tutorial](ai.md).
- `PeptideSpectrumMatcher` is a compact search engine with reversed-sequence decoys and an X!Tandem-style hyperscore. It is not a replacement for MSFragger, Comet or Mascot in production proteomics.

---
Previous: [Transcriptomics](transcriptomics.md) · Next: [Structural biology](structural-biology.md) · [Tutorial index](index.md)
