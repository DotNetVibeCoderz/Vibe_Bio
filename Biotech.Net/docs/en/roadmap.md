# Roadmap

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

The authoritative roadmap is [`PLAN.md`](../../PLAN.md) at the repository root, and day-to-day status is tracked in [`Progress.md`](../../Progress.md). This page summarises where the project stands.

## Shipped

| Milestone | Contents |
|---|---|
| **v1 — Foundation (0.1.x)** | Native Rust kernels with a managed fallback; Core; Sequences (FASTA/faidx, FASTQ, GenBank, EMBL, PDB, GFF/BED, ORFs, motifs, restriction, Tm, codon usage); Alignments (matrices, pairwise, MSA, phylogenetics, database search); Visualization; Notebooks; the Avalonia Gallery |
| **v2 — Omics (0.2.x)** | Genomics (VCF, variant effects, population genetics, de Bruijn assembly, gene finding, read mapping, SNV calling); Transcriptomics (normalisation, DESeq2-style and limma-voom-style DE, ORA/GSEA, single cell); Proteomics (ProtParam, digestion, MS/MS with target–decoy FDR, quantification, classification); StructuralBiology (PDB/mmCIF, geometry, superposition, DSSP, SASA, pockets) |
| **v3 — Discovery (0.3.x)** | Chemistry (SMILES/SDF, descriptors, fingerprints, substructure search); DrugDiscovery (screening, target prediction, ADMET, lead optimisation, QSAR); SyntheticBiology (guide design and scoring, base and prime editing models, codon optimisation, primers, Gibson and Golden Gate, circuit simulation); AI (HF.Net model zoo, offline preprocessing, embeddings, ESM variant scoring, LLM copilot) |

## Added after 0.1.0 (next release)

| Area | Contents |
|---|---|
| **Statistics hardening** | DESeq2 pipeline (size factors incl. poscounts, dispersions, Cook's replacement, Wald/LRT, independent filtering) and apeglm / normal LFC shrinkage reproduce R DESeq2 1.52.0 and apeglm 1.34.0 on pasilla and airway — see the [Transcriptomics tutorial](tutorials/transcriptomics.md) |
| **BAM/CRAM** | BGZF, SAM/BAM read/write with BAI/CSI, CRAM 3.0/3.1 reading with CRAI, tabix, depth — verified against samtools/htslib 1.24 ([Genomics tutorial](tutorials/genomics.md)) |
| **3-D chemistry** | ETKDG-style conformers, UFF, Vina-style flexible docking with scores matching AutoDock Vina 1.2.7 ([Chemistry](tutorials/chemistry.md), [Drug discovery](tutorials/drug-discovery.md)) |
| **GPU kernels** | `Biotech.Net.Gpu` (ILGPU: CUDA / OpenCL / CPU) for all-vs-all alignment scores and batched docking, with CPU fallback ([GPU tutorial](tutorials/gpu.md)) |
| **Workflow engine** | `Biotech.Net.Workflows`: Snakemake-style rules, DAG scheduling, incremental re-runs ([Workflows tutorial](tutorials/workflows.md)) |
| **Cloud data connectors** | `Biotech.Net.Data`: NCBI E-utilities, UniProt, RCSB PDB, AlphaFold DB, Ensembl REST with caching and offline mode ([Data tutorial](tutorials/data.md)) |
| **Single cell at atlas scale** | AnnData h5ad read/write with backed mode, out-of-core sparse matrices, QC / normalisation / HVG / PCA matching scanpy 1.12 |
| **AI** | Biotech.Net.AI now runs ESM and Geneformer on HF.Net 0.11 instead of its own copies |

## Still open

- Web gallery (Blazor WebAssembly) and a hosted playground — in progress.
- CRAM writing and LZMA-compressed CRAM blocks.
- MMFF94, flexible receptor side chains, Vinardo/AD4 scoring.
- Vulkan backend; GPU traceback and global alignment.
- Approximate nearest neighbours for million-cell kNN graphs; out-of-core markers.
- Cluster/container executors for workflows.

## Principles the roadmap is held to

| | |
|---|---|
| Native | idiomatic .NET 10 APIs, spans and streaming everywhere |
| Fast | Rust kernels for hot paths, with the managed fallback always available |
| Cross-platform | win-x64, linux-x64, osx-arm64 |
| Notebook-first | Polyglot / .NET Interactive / Jupyter formatters |
| AI-ready | HF.Net bio models plus an LLM copilot |
| Honest | every algorithm is tested against a published reference, and every page states its limitations |

## Releases

1. The solution moves into the `Vibe_Bio` monorepo as the `Biotech.Net/` subfolder.
2. CI builds the Rust kernels for all three RIDs, runs every test, and packs every library with the native binaries inside.
3. Pushing a `biotechnet-v<version>` tag publishes the `Gravicode.Biotech.Net.*` packages to nuget.org.

---
[Documentation index](index.md) · Bahasa Indonesia: [roadmap](../id/roadmap.md)
