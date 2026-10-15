# Getting started

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Biotech.Net is an open-source biotechnology and bioinformatics platform for .NET 10 — think Bioconductor + Biopython, with Rust-accelerated kernels and a managed fallback that is always available.

## Install

```bash
dotnet new console -o BioDemo
cd BioDemo
dotnet add package Gravicode.Biotech.Net.Core
dotnet add package Gravicode.Biotech.Net.Sequences
```

Full instructions, the package list and offline/CI notes are in [Installation](installation.md).

## Five minutes

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;
using Biotech.Net.Sequences.IO;

Console.WriteLine(BiotechInfo.Banner);     // version, kernel backend and credits

// 1. A sequence keeps the case you give it; kernels and comparisons ignore case
var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
Console.WriteLine($"{dna.Length} bp, GC {dna.GcContent():P1}");
Console.WriteLine(dna.ReverseComplement());
Console.WriteLine(dna.Translate(toStop: true));

// 2. Read and write the usual formats (streaming; .gz handled transparently)
foreach (var record in FastaReader.Parse(">a\nACGTACGT\n>b\nGGGCCC\n"))
    Console.WriteLine($"{record.Id}: {record.Length} bp");

// 3. Everyday analyses
Console.WriteLine($"{OrfFinder.Find(dna, new OrfOptions { MinLength = 30 }).Count} ORFs");
Console.WriteLine($"Tm {MeltingTemperature.NearestNeighbor(Sequence.Dna("AGCGGATAACAATTTCACACAGG")):F1} °C");
```

## Where things live

| Package | What it does |
|---|---|
| `Gravicode.Biotech.Net.Native` | Rust kernels (sequence, alignment, statistics, protein, structure) with a managed fallback |
| `Gravicode.Biotech.Net.Core` | `Sequence`, `Alphabet`, `GeneticCode`, `GenomicInterval`, `Gene`/`Transcript`/`Genome`, `Protein`, `Variant`, `Dataset` |
| `Gravicode.Biotech.Net.Sequences` | FASTA/FASTQ/GenBank/EMBL/PDB I/O, ORFs, motifs, restriction, k-mers, Tm, codon usage |
| `Gravicode.Biotech.Net.Alignments` | matrices, pairwise and multiple alignment, phylogenetics, database search |
| `Gravicode.Biotech.Net.Genomics` | VCF, variant effects, read mapping, SNV calling, assembly, gene finding, population genetics |
| `Gravicode.Biotech.Net.Transcriptomics` | normalisation, DESeq/limma-voom, enrichment, single-cell |
| `Gravicode.Biotech.Net.Proteomics` | ProtParam, digestion, MS/MS search, quantification, classification |
| `Gravicode.Biotech.Net.StructuralBiology` | PDB/mmCIF, geometry, DSSP, SASA, superposition, pockets |
| `Gravicode.Biotech.Net.Chemistry` | SMILES/SDF, descriptors, fingerprints, substructure search |
| `Gravicode.Biotech.Net.DrugDiscovery` | screening, target prediction, ADMET, lead optimisation, QSAR |
| `Gravicode.Biotech.Net.SyntheticBiology` | codon optimisation, primers, assembly, circuits, CRISPR analysis |
| `Gravicode.Biotech.Net.AI` | HF.Net model zoo, embeddings, ESM variant scoring, LLM copilot |
| `Gravicode.Biotech.Net.Visualization` | genome browser, heatmap, volcano, PCA/UMAP, protein and molecule viewers |
| `Gravicode.Biotech.Net.Notebooks` | rich formatters for .NET Interactive / Polyglot Notebooks |

Dependencies point one way: Native ← Core ← everything else, with Visualization and Notebooks at the top.

## Conventions worth knowing

- **Coordinates.** `GenomicInterval` and `FeatureLocation` are 0-based half-open. `Variant.Position` is 1-based like VCF, and `ToString()` prints the 1-based form.
- **Sequences.** `Sequence` stores ASCII bytes and preserves case; comparisons and kernels are case-insensitive. Parsers build a `byte[]` and call `Sequence.FromBuffer` with no copy.
- **Kernels.** `BioKernels` dispatches to Rust when `biotech_native` loads and to C# otherwise. `BioKernels.Backend` tells you which, `BioKernels.ForceManaged = true` pins the managed path, and the environment variable `BIOTECH_NET_DISABLE_NATIVE=1` does the same from outside. Results agree between back-ends; only speed differs ([benchmarks](benchmarks.md)).
- **Streaming.** Readers return `IEnumerable<T>`, so large files are never fully materialised.

## Next steps

- [Tutorials](tutorials/index.md) — one per module, every snippet compiled and run against the real API
- [API reference](api-reference/index.md) — key types and methods per module
- [Samples](samples.md) — the Avalonia Gallery and the sample notebooks
- [Benchmarks](benchmarks.md) — measured Rust-vs-managed numbers
- [Roadmap](roadmap.md)

Bahasa Indonesia: [docs/id](../id/index.md).
