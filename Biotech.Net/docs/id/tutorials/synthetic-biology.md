# Tutorial: Biologi sintetik

Paket: `Gravicode.Biotech.Net.SyntheticBiology`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Optimasi kodon dengan kendala kelayakan sintesis, desain primer (Tm nearest-neighbour, dimer, hairpin), simulasi perakitan Gibson dan Golden Gate, penanganan plasmid, model kekuatan ribosome binding site, registri parts, serta simulasi ODE sirkuit genetik. Modul ini juga memuat tipe analisis CRISPR yang diringkas di akhir halaman.

## 1. Optimasi kodon

`CodonUsageTable` menyediakan tabel penggunaan kodon untuk *E. coli* K-12, *S. cerevisiae*, manusia, CHO (*C. griseus*), dan *P. pastoris*. `SequenceConstraints` menggambarkan apa yang bisa diterima vendor sintesis, dan pengoptimal memperbaiki draf pertama sampai semua kendala terpenuhi.

```csharp
using Biotech.Net.SyntheticBiology.Codons;

var options = new CodonOptimizationOptions
{
    Target = CodonUsageTable.EColiK12,
    Method = CodonOptimizationMethod.MostFrequent,
    Constraints = SequenceConstraints.ForGeneSynthesis,
};
var result = CodonOptimizer.Optimize("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", options);
Console.WriteLine($"CAI {result.Cai:F3}, GC {result.GcContent:P1}, sisa pelanggaran {result.RemainingViolations.Count}");
Console.WriteLine(result.Dna);

// Mengkode ulang gen yang ada untuk inang lain, dengan protein yang identik
var humanOptions = new CodonOptimizationOptions
{
    Target = CodonUsageTable.HomoSapiens,
    Constraints = SequenceConstraints.ForGeneSynthesis,
};
var recoded = CodonOptimizer.OptimizeCoding(result.Dna, humanOptions);
Console.WriteLine($"kodon berubah: {recoded.CodonsChanged}");
```

## 2. Desain primer

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;
using Biotech.Net.SyntheticBiology.Primers;

string template = RandomSequence.Dna(1_200, seed: 8, gcContent: 0.5).ToString();
var pairs = PrimerDesigner.DesignPairs(template, new PrimerDesignOptions
{
    OptTm = 60,
    ProductSizeMin = 300,
    ProductSizeMax = 600,
    MaxPairsReturned = 3,
});
foreach (var p in pairs)
    Console.WriteLine($"{p.Forward.Sequence} / {p.Reverse.Sequence}  {p.ProductSize} bp  Tm {p.Forward.Tm:F1}/{p.Reverse.Tm:F1}  penalti {p.Penalty:F2}");

// Mengevaluasi oligo yang sudah ada
var primer = PrimerDesigner.Evaluate("AGCGGATAACAATTTCACACAGG");
Console.WriteLine($"Tm {primer.Tm:F1} °C, GC {primer.GcFraction:P0}, hairpin dG {primer.HairpinDeltaG:F2} kkal/mol");
Console.WriteLine($"Tm nearest-neighbour: {NearestNeighbourTm.Tm("AGCGGATAACAATTTCACACAGG"):F1} °C");
```

## 3. Simulasi perakitan

```csharp
using Biotech.Net.SyntheticBiology.Cloning;

// Gibson: rancang tumpang-tindih antar fragmen, lalu simulasikan penyambungannya
string a = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTC";
string b = "AGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACA";
var design = GibsonAssembly.DesignOverlaps([a, b]);
Console.WriteLine($"{design.Junctions.Count} sambungan, ideal: {design.IsIdeal}");
var gibson = GibsonAssembly.Simulate(design.TailedFragments);
Console.WriteLine(gibson);

// Golden Gate dengan enzim Tipe IIS
var bsaI = RestrictionEnzyme.FromName("BsaI");
string part1 = GoldenGate.AddFlanks(a, bsaI, "AATG", "GCTT");
string part2 = GoldenGate.AddFlanks(b, bsaI, "GCTT", "AATG");
var gg = GoldenGate.Simulate([part1, part2], bsaI);
Console.WriteLine(gg);

// Plasmid: pemotong tunggal, digesti, dan ringkasan gaya GenBank
var plasmid = Plasmid.DemonstrationBackbone;
Console.WriteLine($"{plasmid}  pemotong tunggal: {string.Join(", ", plasmid.UniqueCutters().Take(5).Select(e => e.Name))}");
```

## 4. Kekuatan RBS dan sirkuit genetik

```csharp
using Biotech.Net.SyntheticBiology.Circuits;

var rbs = RbsStrengthModel.Predict("TAAGGAGGTAAAAA", "ATGAAAGAAGCAGAAGCAAAAGAAGCATAA");
Console.WriteLine($"laju inisiasi {rbs.TranslationInitiationRate:F0} (a.u.), spasi {rbs.SpacingNucleotides} nt, kecocokan SD {rbs.ShineDalgarnoMatch}");

// Registri parts dan unit transkripsi
var unit = PartsRegistry.BuildUnit("reporter", "J23100", "B0034", "GFP", "B0015");
Console.WriteLine($"{unit}  GC {unit.ToSequence().GcContent():P1}");

// Repressilator Elowitz, diintegrasikan dengan RK4
var solution = CircuitModels.SimulateRepressilator(tEnd: 300);
double period = CircuitModels.EstimateOscillationPeriod(solution, stateIndex: 3);
Console.WriteLine($"{solution.StepCount} langkah keluaran, periode ≈ {period:F1} satuan waktu");

// Toggle switch dan ekspresi terinduksi juga tersedia
var toggle = CircuitModels.SimulateToggleSwitch(alpha1: 50, alpha2: 16, beta: 2.5, gamma: 1,
    initial: [0.1, 0.2], tEnd: 50);
Console.WriteLine($"keadaan akhir: {string.Join(", ", toggle.Final.Select(v => v.ToString("F3")))}");
```

## 5. Tipe analisis CRISPR

Namespace `Biotech.Net.SyntheticBiology.Crispr` memodelkan skema penilaian terpublikasi agar kumpulan guide dapat dianalisis dan dibandingkan di .NET:

| Tipe | Kegunaan |
|---|---|
| `CasNuclease` | PAM, panjang guide, dan geometri pemotongan untuk SpCas9, SaCas9, Cas12a, SpG, SpRY, NmeCas9, dan varian rekayasa (`CasNuclease.All`) |
| `GuideFinder`, `Guide` | mencacah protospacer bersebelahan PAM, dengan pemeriksaan GC, homopolimer, dan poli-T |
| `OnTargetScorer` | Rule Set 1 Doench dkk. 2014 pada konteks 30 nt |
| `CrisprSiteIndex`, `OffTargetSearch` | pencarian ketidakcocokan berumpan pada sekuens terindeks |
| `CfdScorer`, `MitSpecificityScorer` | agregasi spesifisitas CFD (Doench dkk. 2016) dan MIT |
| `GuideDesigner`, `GuideDesign` | menggabungkan semuanya menjadi tabel berperingkat |
| `BaseEditing`, `BaseEditor` | model jendela pengeditan untuk base editor sitosin dan adenin terpublikasi |
| `PrimeEditor`, `PegRna` | tata letak komponen pegRNA (spacer, scaffold, PBS, RTT) dan kandidat nicking guide |

Lihat [referensi API](../api-reference/synthetic-biology.md) untuk daftar anggotanya, dan halaman CRISPR Designer di [Gallery](../samples.md) untuk tampilan interaktifnya.

## Batasan

- Penilaian on-target memakai Rule Set 1 (model 2014). Rule Set 2 / Azimuth dan model gaya DeepCRISPR tidak disertakan, dan aksesibilitas kromatin sama sekali tidak dimodelkan.
- Pencarian off-target berbasis sekuens pada data yang Anda indeks sendiri. Ini bukan layanan pencarian genom utuh, dan `MaxBulges > 0` lebih lambat karena penjajaran bergelembung dicacah secara eksplisit.
- Model base editor dan prime editor adalah model bobot jendela dan tata letak dari literatur primer, bukan prediktor efisiensi terlatih; gunakan untuk membandingkan kandidat, bukan memperkirakan laju pengeditan absolut.
- Model RBS adalah model energi bebas sederhana dalam semangat RBS Calculator dan dilaporkan dalam satuan arbitrer. Belum dikalibrasi terhadap data inang tertentu.
- Simulasi sirkuit mencakup model ODE klasik (toggle switch, repressilator, ekspresi terinduksi) dengan RK4/RK45. Belum ada simulasi stokastik (Gillespie).

---
Sebelumnya: [Penemuan obat](drug-discovery.md) · Selanjutnya: [AI](ai.md) · [Indeks tutorial](index.md)
