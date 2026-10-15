# Pendahuluan

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Biotech.Net adalah platform bioteknologi dan bioinformatika open-source untuk .NET 10 — bayangkan Bioconductor + Biopython, dengan kernel yang dipercepat Rust dan fallback managed yang selalu tersedia.

## Instalasi singkat

```bash
dotnet new console -o BioDemo
cd BioDemo
dotnet add package Gravicode.Biotech.Net.Core
dotnet add package Gravicode.Biotech.Net.Sequences
```

Petunjuk lengkap, daftar paket, dan catatan offline/CI ada di [Instalasi](installation.md).

## Lima menit pertama

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;
using Biotech.Net.Sequences.IO;

Console.WriteLine(BiotechInfo.Banner);     // versi, backend kernel, dan kredit

// 1. Sequence menyimpan huruf apa adanya; kernel dan perbandingan mengabaikan besar-kecil huruf
var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
Console.WriteLine($"{dna.Length} bp, GC {dna.GcContent():P1}");
Console.WriteLine(dna.ReverseComplement());
Console.WriteLine(dna.Translate(toStop: true));

// 2. Membaca dan menulis format umum (streaming; berkas .gz ditangani otomatis)
foreach (var record in FastaReader.Parse(">a\nACGTACGT\n>b\nGGGCCC\n"))
    Console.WriteLine($"{record.Id}: {record.Length} bp");

// 3. Analisis sehari-hari
Console.WriteLine($"{OrfFinder.Find(dna, new OrfOptions { MinLength = 30 }).Count} ORF");
Console.WriteLine($"Tm {MeltingTemperature.NearestNeighbor(Sequence.Dna("AGCGGATAACAATTTCACACAGG")):F1} °C");
```

## Peta modul

| Paket | Isinya |
|---|---|
| `Gravicode.Biotech.Net.Native` | kernel Rust (sekuens, penjajaran, statistik, protein, struktur) dengan fallback managed |
| `Gravicode.Biotech.Net.Core` | `Sequence`, `Alphabet`, `GeneticCode`, `GenomicInterval`, `Gene`/`Transcript`/`Genome`, `Protein`, `Variant`, `Dataset` |
| `Gravicode.Biotech.Net.Sequences` | I/O FASTA/FASTQ/GenBank/EMBL/PDB, ORF, motif, enzim restriksi, k-mer, Tm, penggunaan kodon |
| `Gravicode.Biotech.Net.Alignments` | matriks substitusi, penjajaran berpasangan dan ganda, filogenetika, pencarian basis data |
| `Gravicode.Biotech.Net.Genomics` | VCF, efek varian, pemetaan bacaan, pemanggilan SNV, perakitan, pencarian gen, genetika populasi |
| `Gravicode.Biotech.Net.Transcriptomics` | normalisasi, DESeq/limma-voom, pengayaan, sel tunggal |
| `Gravicode.Biotech.Net.Proteomics` | ProtParam, digesti, pencarian MS/MS, kuantifikasi, klasifikasi |
| `Gravicode.Biotech.Net.StructuralBiology` | PDB/mmCIF, geometri, DSSP, SASA, superposisi, kantong ikatan |
| `Gravicode.Biotech.Net.Chemistry` | SMILES/SDF, deskriptor, sidik jari molekul, pencarian substruktur |
| `Gravicode.Biotech.Net.DrugDiscovery` | penapisan, prediksi target, ADMET, optimasi lead, QSAR |
| `Gravicode.Biotech.Net.SyntheticBiology` | optimasi kodon, primer, perakitan DNA, sirkuit, analisis CRISPR |
| `Gravicode.Biotech.Net.AI` | model zoo HF.Net, embedding, skor varian ESM, kopilot LLM |
| `Gravicode.Biotech.Net.Visualization` | genome browser, heatmap, volcano, PCA/UMAP, penampil protein dan molekul |
| `Gravicode.Biotech.Net.Notebooks` | formatter kaya untuk .NET Interactive / Polyglot Notebooks |

Arah ketergantungan satu arah: Native ← Core ← modul lain, dengan Visualization dan Notebooks di lapisan paling atas.

## Konvensi yang perlu diketahui

- **Koordinat.** `GenomicInterval` dan `FeatureLocation` berbasis 0 dan setengah terbuka. `Variant.Position` berbasis 1 seperti VCF, dan `ToString()` mencetak bentuk berbasis 1.
- **Sekuens.** `Sequence` menyimpan byte ASCII dan mempertahankan besar-kecil huruf; perbandingan dan kernel mengabaikannya. Parser membangun `byte[]` lalu memanggil `Sequence.FromBuffer` tanpa penyalinan.
- **Kernel.** `BioKernels` memilih Rust bila `biotech_native` berhasil dimuat, dan C# bila tidak. `BioKernels.Backend` memberitahu yang dipakai, `BioKernels.ForceManaged = true` mengunci jalur managed, dan variabel lingkungan `BIOTECH_NET_DISABLE_NATIVE=1` melakukan hal yang sama dari luar. Hasil kedua backend sama; yang berbeda hanya kecepatan ([benchmark](benchmarks.md)).
- **Streaming.** Pembaca berkas mengembalikan `IEnumerable<T>`, jadi berkas besar tidak pernah dimuat seluruhnya.

## Langkah berikutnya

- [Tutorial](tutorials/index.md) — satu per modul, semua cuplikan dikompilasi dan dijalankan terhadap API nyata
- [Referensi API](api-reference/index.md) — tipe dan metode utama per modul
- [Contoh](samples.md) — Gallery Avalonia dan notebook contoh
- [Benchmark](benchmarks.md) — angka Rust vs managed hasil pengukuran
- [Roadmap](roadmap.md)

English: [docs/en](../en/index.md).
