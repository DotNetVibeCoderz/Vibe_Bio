# Biotech.Net

**Ultimate Biotechnology & Bioinformatics Platform for .NET** — the Bioconductor + Biopython of .NET 10, with Rust-accelerated kernels.

*Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil*

---

## English

Biotech.Net is an open-source platform for biotechnology, bioinformatics, computational biology, genomics, proteomics, drug discovery, synthetic biology and biotech AI, built for .NET 10. Hot paths run on Rust kernels through a C ABI; every kernel also has a managed C# implementation, so the library works even without the native binary.

### Quick start

```bash
dotnet add package Gravicode.Biotech.Net.Core
dotnet add package Gravicode.Biotech.Net.Sequences
```

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
Console.WriteLine($"{dna.Length} bp, GC {dna.GcContent():P1}");
Console.WriteLine(dna.ReverseComplement());
Console.WriteLine(dna.Translate(toStop: true));              // MAIVMGR
Console.WriteLine($"{OrfFinder.Find(dna, new OrfOptions { MinLength = 30 }).Count} ORFs");
```

### Modules

| Package | What it does |
|---|---|
| `Gravicode.Biotech.Net.Native` | Rust kernels (sequence, alignment, statistics, protein, structure) with a managed fallback |
| `Gravicode.Biotech.Net.Core` | `Sequence`, `Alphabet`, `GeneticCode`, `GenomicInterval`, `Gene`/`Transcript`/`Genome`, `Protein`, `Variant`, `Dataset` |
| `Gravicode.Biotech.Net.Sequences` | FASTA/FASTQ/GenBank/EMBL/PDB I/O, ORFs, motifs, restriction, k-mers, Tm, codon usage |
| `Gravicode.Biotech.Net.Alignments` | substitution matrices, pairwise and multiple alignment, phylogenetics, database search |
| `Gravicode.Biotech.Net.Genomics` | VCF, variant effects, read mapping, SNV calling, assembly, gene finding, population genetics |
| `Gravicode.Biotech.Net.Transcriptomics` | normalisation, DESeq2-style and limma-voom-style DE, ORA/GSEA, single cell |
| `Gravicode.Biotech.Net.Proteomics` | ProtParam, digestion, MS/MS search with target–decoy FDR, quantification, classification |
| `Gravicode.Biotech.Net.StructuralBiology` | PDB/mmCIF, geometry, DSSP, SASA, superposition, pockets, interactions |
| `Gravicode.Biotech.Net.Chemistry` | SMILES/SDF, descriptors, fingerprints, substructure search, scaffolds |
| `Gravicode.Biotech.Net.DrugDiscovery` | virtual screening, target prediction, ADMET, lead optimisation, QSAR |
| `Gravicode.Biotech.Net.SyntheticBiology` | codon optimisation, primers, Gibson/Golden Gate, circuits, CRISPR analysis |
| `Gravicode.Biotech.Net.AI` | HF.Net model zoo, embeddings, ESM variant scoring, LLM copilot |
| `Gravicode.Biotech.Net.Visualization` | genome browser, heatmap, volcano, PCA/UMAP, protein and molecule viewers |
| `Gravicode.Biotech.Net.Notebooks` | rich formatters for .NET Interactive / Polyglot Notebooks |

### Documentation

- English: [`docs/en`](docs/en/index.md) — getting started, installation, tutorials, API reference, samples, benchmarks, roadmap
- Bahasa Indonesia: [`docs/id`](docs/id/index.md)

Every C# snippet in the documentation is compiled and executed against the real API before publication, and every figure is rendered by `Biotech.Net.Visualization` itself.

### Performance

Measured on an Intel Core i7-8650U (BenchmarkDotNet, short job). Ratio is Rust time ÷ managed time, so lower is faster:

| Kernel | Managed | Rust | Ratio |
|---|---:|---:|---:|
| Local score batch (300 aa × 500 targets) | 295,240 µs | 6,270 µs | 0.02 |
| Reverse complement (1 Mb) | 625.8 µs | 80.7 µs | 0.13 |
| Kabsch RMSD (10,000 atoms) | 246.5 µs | 96.1 µs | 0.39 |
| k-mer counting, k = 8 (1 Mb) | 4,698 µs | 2,014 µs | 0.44 |

Full results, including the cases where the managed path currently wins: [`docs/en/benchmarks.md`](docs/en/benchmarks.md).

---

## Bahasa Indonesia

Biotech.Net adalah platform open-source untuk bioteknologi, bioinformatika, biologi komputasi, genomika, proteomika, penemuan obat, biologi sintetik, dan AI bioteknologi, dibangun untuk .NET 10. Jalur kritis berjalan di kernel Rust lewat C ABI; setiap kernel juga punya implementasi C# managed, jadi pustaka ini tetap bekerja tanpa biner native.

### Mulai cepat

```bash
dotnet add package Gravicode.Biotech.Net.Core
dotnet add package Gravicode.Biotech.Net.Sequences
```

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
Console.WriteLine($"{dna.Length} bp, GC {dna.GcContent():P1}");
Console.WriteLine(dna.ReverseComplement());
Console.WriteLine(dna.Translate(toStop: true));              // MAIVMGR
Console.WriteLine($"{OrfFinder.Find(dna, new OrfOptions { MinLength = 30 }).Count} ORF");
```

### Modul

Daftar paket sama seperti tabel di bagian Inggris di atas: Native, Core, Sequences, Alignments, Genomics, Transcriptomics, Proteomics, StructuralBiology, Chemistry, DrugDiscovery, SyntheticBiology, AI, Visualization, dan Notebooks — semuanya diawali `Gravicode.Biotech.Net.`.

### Dokumentasi

- Bahasa Indonesia: [`docs/id`](docs/id/index.md) — pendahuluan, instalasi, tutorial, referensi API, contoh, benchmark, roadmap
- English: [`docs/en`](docs/en/index.md)

Setiap cuplikan C# di dokumentasi dikompilasi dan dijalankan terhadap API yang sebenarnya sebelum diterbitkan, dan setiap gambar dirender oleh `Biotech.Net.Visualization` sendiri.

---

## Screenshots — Biotech.Net Gallery

The Avalonia showcase (`samples/Biotech.Net.Gallery`): six modules, each running real library code with a live code tab and PNG/SVG export.

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-dna-explorer.png" alt="DNA Explorer" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-genome-browser.png" alt="Genome Browser" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-protein-explorer.png" alt="Protein Explorer" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-drug-discovery-lab.png" alt="Drug Discovery Lab" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-ai-playground.png" alt="AI Playground" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/screenshots/gallery-crispr-designer.png" alt="CRISPR Designer" width="720" />

Example figures produced by `Biotech.Net.Visualization`:

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/images/genome-browser.png" alt="Genome browser plot" width="720" />

<img src="https://raw.githubusercontent.com/DotNetVibeCoderz/Vibe_Bio/main/Biotech.Net/docs/images/volcano.png" alt="Volcano plot" width="560" />

---

## Building from source

```bash
dotnet build Biotech.Net.slnx                  # builds the Rust crate too when cargo is on PATH
dotnet build -p:SkipRustBuild=true             # reuse an existing native binary
dotnet test Biotech.Net.slnx
cd native/biotech-native && cargo test --release
```

Requires the .NET 10 SDK. The Rust toolchain is optional — without it every kernel falls back to managed C#.

## Credits & license

Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil · Created by Gravicode Studios, led by Kang Fadhil.

Licensed under the [MIT License](LICENSE). Copyright © Gravicode Studios.
