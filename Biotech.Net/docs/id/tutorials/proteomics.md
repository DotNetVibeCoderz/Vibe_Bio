# Tutorial: Proteomika

Paket: `Gravicode.Biotech.Net.Proteomics`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Analisis protein gaya ProtParam, prediktor struktur berbasis sekuens, digesti enzimatik, ion fragmen dan pola isotop, pencocokan peptida–spektrum dengan FDR target–decoy, inferensi protein, kuantifikasi label-free, serta klasifikasi famili protein.

## 1. Sifat protein (ProtParam)

```csharp
using Biotech.Net.Proteomics;

// Rantai B insulin manusia
var analysis = new ProteinAnalysis("FVNQHLCGSHLVEALYLVCGERGFFYTPKT");
Console.WriteLine($"BM        {analysis.MolecularWeight():F2} Da");
Console.WriteLine($"pI        {analysis.IsoelectricPoint():F2}");
Console.WriteLine($"GRAVY     {analysis.Gravy():F3}");
Console.WriteLine($"Instab.   {analysis.InstabilityIndex():F2} (stabil: {analysis.IsStable})");
Console.WriteLine($"Alifatik  {analysis.AliphaticIndex():F2}");
Console.WriteLine($"Koef. ext {analysis.ExtinctionCoefficient().Reduced} M^-1 cm^-1");
Console.WriteLine($"Waktu paruh {analysis.HalfLife().Text}");

// Profil hidropati (Kyte-Doolittle, jendela 9) dan laporan teks lengkap
double[] profile = analysis.Profile(AminoAcidScale.KyteDoolittle, window: 9);
Console.WriteLine(analysis.Report());

// Prediktor berbasis sekuens
Console.WriteLine(ChouFasman.Predict("FVNQHLCGSHLVEALYLVCGERGFFYTPKT"));
var tm = TransmembranePredictor.Predict("MKTLLILAVVAAALAGSLLLVVFLLAGIAALLGYKRRSS");
var signal = SignalPeptidePredictor.Predict("MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQ");
Console.WriteLine($"{tm.Count} segmen TM; peptida sinyal: {signal.HasSignalPeptide} (pemotongan setelah {signal.CleavageSite})");
```

`MolecularWeight` dan `IsoelectricPoint` berjalan lewat kernel protein native, dan tipe `Protein` di Core memberikan nilai yang sama.

## 2. Digesti, fragmen, dan isotop

```csharp
using Biotech.Net.Proteomics.MassSpec;

var peptides = Digestion.Digest("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF",
    Enzyme.Trypsin, new DigestionOptions(MissedCleavages: 1, MinLength: 5));
foreach (var p in peptides.Take(4))
    Console.WriteLine($"{p.Sequence,-20} {p.MonoisotopicMass:F4} Da  terlewat={p.MissedCleavages}");

var peptide = new Peptide("SAMPLER");
Console.WriteLine($"[M+2H]2+ = {peptide.Mz(2):F4}");
foreach (var ion in FragmentCalculator.Fragments(peptide, [IonType.B, IonType.Y]).Take(6))
    Console.WriteLine(ion);

foreach (var peak in IsotopeDistribution.FromPeptide("SAMPLER", maxPeaks: 4))
    Console.WriteLine($"{peak.Mass:F4}  {peak.Abundance:F3}");
```

## 3. Pencocokan peptida–spektrum dengan FDR target–decoy

```csharp
using Biotech.Net.Core;
using Biotech.Net.Proteomics.MassSpec;

var proteins = new[]
{
    Protein.FromString("P1", "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFEDHVK"),
    Protein.FromString("P2", "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEFELRRQACGRTRAPPPPPLRSGC"),
};

// Spektrum MS/MS sintetis untuk DLGEEHFK (2+) dari ion b/y-nya
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

// Data nyata: MgfReader.ReadFile("run.mgf") atau MzMLReader.ReadFile("run.mzML")
```

`ProteinInference.Parsimony(psms)` mengelompokkan protein berdasarkan peptida bersama, dan `LabelFreeQuantification.Quantify(psms, proteins)` melaporkan spectral count, NSAF, dan iBAQ.

## 4. Klasifikasi famili protein

```csharp
using Biotech.Net.Proteomics.Classification;

// Data sintetis empat famili untuk demonstrasi dan pengujian
var examples = ProteinFamilyDemo.Generate(perFamily: 20, seed: 1);
var cv = ProteinFamilyClassifier.CrossValidate(examples, () => new LogisticRegressionClassifier(), folds: 5);
Console.WriteLine($"akurasi 5-fold {cv.Accuracy:P1}, macro-F1 {cv.MacroF1:F3}");

var classifier = new ProteinFamilyClassifier(new KNearestNeighborsClassifier(k: 3)).Fit(examples);
Console.WriteLine(classifier.Predict(examples[0].Protein));
```

## Batasan

- `ProteinFamilyDemo` menghasilkan famili sintetis untuk demonstrasi. Latihlah dengan data berlabel nyata (misalnya famili Pfam) sebelum menarik kesimpulan biologis.
- Chou–Fasman, prediktor transmembran, dan prediktor peptida sinyal adalah heuristik klasik yang jauh kurang akurat dibanding SignalP, TMHMM/DeepTMHMM, atau anotasi turunan AlphaFold. Untuk fitur dari model bahasa protein, lihat [tutorial AI](ai.md).
- `PeptideSpectrumMatcher` adalah mesin pencarian ringkas dengan decoy sekuens terbalik dan hyperscore gaya X!Tandem. Bukan pengganti MSFragger, Comet, atau Mascot di proteomika produksi.

---
Sebelumnya: [Transkriptomika](transcriptomics.md) · Selanjutnya: [Biologi struktural](structural-biology.md) · [Indeks tutorial](index.md)
