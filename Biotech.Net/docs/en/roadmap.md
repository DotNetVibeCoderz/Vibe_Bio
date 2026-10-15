# Roadmap

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

The authoritative roadmap is [`PLAN.md`](../../PLAN.md) at the repository root, and day-to-day status is tracked in [`Progress.md`](../../Progress.md). This page summarises where the project stands.

## Shipped

| Milestone | Contents |
|---|---|
| **v1 — Foundation (0.1.x)** | Native Rust kernels with a managed fallback; Core; Sequences (FASTA/faidx, FASTQ, GenBank, EMBL, PDB, GFF/BED, ORFs, motifs, restriction, Tm, codon usage); Alignments (matrices, pairwise, MSA, phylogenetics, database search); Visualization; Notebooks; the Avalonia Gallery |
| **v2 — Omics (0.2.x)** | Genomics (VCF, variant effects, population genetics, de Bruijn assembly, gene finding, read mapping, SNV calling); Transcriptomics (normalisation, DESeq2-style and limma-voom-style DE, ORA/GSEA, single cell); Proteomics (ProtParam, digestion, MS/MS with target–decoy FDR, quantification, classification); StructuralBiology (PDB/mmCIF, geometry, superposition, DSSP, SASA, pockets) |
| **v3 — Discovery (0.3.x)** | Chemistry (SMILES/SDF, descriptors, fingerprints, substructure search); DrugDiscovery (screening, target prediction, ADMET, lead optimisation, QSAR); SyntheticBiology (guide design and scoring, base and prime editing models, codon optimisation, primers, Gibson and Golden Gate, circuit simulation); AI (HF.Net model zoo, offline preprocessing, embeddings, ESM variant scoring, LLM copilot) |

## In progress / next

- **Statistics hardening** — exact DESeq2 parity on public data sets, apeglm-style LFC shrinkage.
- **BAM/CRAM** — BGZF reading with an index, plus tabix for VCF and GFF.
- **GPU kernels** — ILGPU or Vulkan for all-vs-all alignment and docking.
- **3-D chemistry** — ETKDG-like conformer generation and flexible docking.

## Beyond

- A Blazor web gallery and a hosted playground.
- A workflow engine: a Snakemake/Nextflow-like DAG in C#.
- Cloud data connectors (NCBI E-utilities, UniProt, PDB, Ensembl REST) with caching.
- Single cell at atlas scale: out-of-core sparse matrices and an AnnData/h5ad reader.

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
