# Biotech.Net — Roadmap / Rencana Pengembangan

> Visi: kombinasi **Bioconductor + Biopython** untuk platform .NET.
> Vision: the **Bioconductor + Biopython** of .NET.

Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.

Status per item is tracked in [Progress.md](Progress.md). The architecture is described in [solution-design.md](solution-design.md).

## Principles / Prinsip

| | English | Bahasa Indonesia |
|---|---|---|
| Native | Idiomatic .NET 10 APIs, spans and streaming everywhere | API .NET 10 yang idiomatis, span dan streaming di mana-mana |
| Fast | Rust kernels for hot paths, managed fallback always available | Kernel Rust untuk jalur kritis, fallback managed selalu ada |
| Cross-platform | win-x64, linux-x64, osx-arm64 | win-x64, linux-x64, osx-arm64 |
| Notebook-first | Polyglot / .NET Interactive / Jupyter formatters | Formatter Polyglot / .NET Interactive / Jupyter |
| AI-ready | HF.Net bio models (BioBERT, ESM, DNABERT...) + LLM copilot | Model bio HF.Net + kopilot LLM |
| Honest | Every algorithm is tested against a published reference | Setiap algoritma diuji terhadap referensi terpublikasi |

## v1 — Foundation (0.1.x)

- [x] **Native** — Rust crate `biotech-native` with a C ABI: sequence engine (reverse complement, translation, composition, k-mers, Hamming), statistics engine (Welch t, BH-FDR), protein engine (mass, charge, pI, hydropathy), and alignment engine (Gotoh global/local/semi-global, score-only search, Kabsch RMSD).
- [x] **Core** — Sequence, Alphabet, GeneticCode (NCBI 1/2/4/11), SequenceRead, GenomicInterval, SequenceFeature, Gene/Transcript/Chromosome/Genome, Protein, Variant/Genotype, Dataset.
- [x] **Sequences** — FASTA (+faidx), FASTQ, GenBank, EMBL, PDB sequences, GFF/BED; ORFs, motifs, restriction digestion, Tm, codon usage.
- [x] **Alignments** — substitution matrices, pairwise and multiple alignment, phylogenetic trees, database search.
- [x] **Visualization** — genome browser, heatmap, volcano, PCA, UMAP, protein and molecule viewers (Visual.Net).
- [x] **Notebooks** — rich formatters plus sample `.dib` notebooks.
- [x] **Gallery** — Avalonia showcase: DNA Explorer, Genome Browser, Protein Explorer, Drug Discovery Lab, AI Playground, CRISPR Designer.

## v2 — Omics (0.2.x)

- [x] **Genomics** — VCF, variant effects, population genetics, de Bruijn assembly, gene finding, read mapping, SNV calling.
- [x] **Transcriptomics** — normalisation, DESeq2-style / limma-voom-style differential expression, enrichment (ORA, GSEA), single-cell (QC, HVG, PCA, Leiden/Louvain, UMAP, t-SNE, markers).
- [x] **Proteomics** — ProtParam, digestion, MS/MS matching with target-decoy FDR, quantification, classification.
- [x] **StructuralBiology** — PDB/mmCIF, geometry, superposition, DSSP-style secondary structure, SASA, pockets.
- [x] **Statistics hardening** — DESeq2 pipeline and apeglm / normal LFC shrinkage reproduce R (DESeq2 1.52.0, apeglm 1.34.0) on pasilla and airway; Cook's replacement, LRT, independent filtering, poscounts.
- [x] **BAM/CRAM** — BGZF, SAM/BAM read+write, BAI/CSI, CRAM 3.0/3.1 reading (all codecs except LZMA) + CRAI, tabix (.tbi/.csi) with indexed VCF queries, depth; verified against samtools/htslib 1.24 and the hts-specs conformance files.

## v3 — Discovery (0.3.x)

- [x] **Chemistry** — SMILES/SDF, descriptors, fingerprints, substructure search.
- [x] **DrugDiscovery** — virtual screening, target prediction, ADMET, lead optimisation, QSAR.
- [x] **SyntheticBiology / CRISPR** — guide design with on-/off-target scoring, base and prime editing, codon optimisation, primers, Gibson and Golden Gate assembly, circuit simulation.
- [x] **AI Platform** — HF.Net model zoo with offline preprocessing, embeddings, ESM variant scoring, and an LLM copilot.
- [x] **GPU kernels** (`Biotech.Net.Gpu`, ILGPU: CUDA / OpenCL / CPU) — all-vs-all local alignment scores (bit-identical to the CPU path) and batched Vina pose scoring, BFGS and Monte-Carlo docking, with automatic CPU fallback. Vulkan not done.
- [x] **3-D conformers and flexible docking** — ETKDG-style embedding with the ETKDG v2 + macrocycle torsion library, UFF (matches RDKit to 1e-4 kcal/mol), Vina-style Monte Carlo docking with torsion trees (scores match AutoDock Vina 1.2.7; redock RMSD ≤ ~1 Å on 3PTB and 1STP). MMFF94 and flexible side chains not yet.

## Beyond / Selanjutnya

- [x] **Web gallery** (`samples/Biotech.Net.WebGallery`, Blazor WebAssembly) — seven pages computed in the browser plus a C# playground compiled by Roslyn in WASM; GitHub Pages workflow ready in `build/workflows/biotechnet-webgallery.yml` (not enabled yet).
- [x] **Workflow engine** (`Biotech.Net.Workflows`) — Snakemake-style rules with wildcards and sample sheets, DAG with cycle detection, resource-aware parallel scheduling, hash-based incremental re-runs, atomic outputs, shell and C# steps, reports and DOT/Mermaid export.
- [x] **Cloud data connectors** (`Biotech.Net.Data`) — NCBI E-utilities, UniProt, RCSB PDB + AlphaFold DB, Ensembl REST (incl. VEP); disk/memory cache with TTL and offline mode, per-host rate limiting, retries.
- [x] **Single-cell at atlas scale** — AnnData h5ad read/write (PureHDF) with backed mode, memory-mapped out-of-core sparse matrices, out-of-core QC, normalisation, HVG (Seurat v3 / Seurat) and PCA matching scanpy 1.12; 1M cells × 30k genes in under 1 GB for QC/HVG. Approximate kNN not yet.

## Upstream / Hulu

- [x] **HF.Net 0.11.0** (`Vibe_ML/HF.Net`, published to nuget.org on 2026-10-15) — ESM / ESM-2 encoder and tokenizer, NumPy scalars in the pickle reader, models without a tokenizer, and a new package icon. `Biotech.Net.AI` now references 0.11.0.
- [x] `Biotech.Net.AI` dropped its own ESM implementation, restricted pickle reader and Geneformer `vocab.txt` placeholder in favour of HF.Net 0.11.0 (breaking: `EsmModel`/`EsmTokenizer`/`EsmConfig` removed).
- [ ] **HF.Net 0.11.1** (fixed locally, unpublished) — ESM rotary embeddings use the fp16-rounded `inv_freq` stored in the checkpoint, as transformers 5 does.

## Release / Rilis

1. Move the solution into the monorepo `Vibe_Bio` as the `Biotech.Net/` subfolder (Pindahkan solusi ke monorepo `Vibe_Bio`).
2. CI builds the Rust kernels for win-x64, linux-x64 and osx-arm64, runs every test, and packs every library with the native binaries.
3. Pushing a `biotechnet-v<version>` tag publishes the `Gravicode.Biotech.Net.*` packages to nuget.org (secret `NUGET_API_KEY`).
