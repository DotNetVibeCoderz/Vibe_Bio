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
- [ ] Statistics hardening: exact DESeq2 parity on public datasets, apeglm-style LFC shrinkage.
- [ ] BAM/CRAM reading (BGZF + index), tabix.

## v3 — Discovery (0.3.x)

- [x] **Chemistry** — SMILES/SDF, descriptors, fingerprints, substructure search.
- [x] **DrugDiscovery** — virtual screening, target prediction, ADMET, lead optimisation, QSAR.
- [x] **SyntheticBiology / CRISPR** — guide design with on-/off-target scoring, base and prime editing, codon optimisation, primers, Gibson and Golden Gate assembly, circuit simulation.
- [x] **AI Platform** — HF.Net model zoo with offline preprocessing, embeddings, ESM variant scoring, and an LLM copilot.
- [ ] GPU kernels (ILGPU / Vulkan) for all-vs-all alignment and docking.
- [ ] 3-D conformer generation (ETKDG-like) and flexible docking.

## Beyond / Selanjutnya

- [ ] Blazor web gallery and a hosted playground.
- [ ] Workflow engine (Snakemake/Nextflow-like DAG in C#).
- [ ] Cloud data connectors (NCBI E-utilities, UniProt, PDB, Ensembl REST) with caching.
- [ ] Single-cell at atlas scale (out-of-core sparse matrices, AnnData/h5ad reader).

## Upstream / Hulu

- [x] **HF.Net 0.11.0** (`Vibe_ML/HF.Net`, published to nuget.org on 2026-10-15) — ESM / ESM-2 encoder and tokenizer, NumPy scalars in the pickle reader, models without a tokenizer, and a new package icon. `Biotech.Net.AI` now references 0.11.0.
- [ ] Let `Biotech.Net.AI` drop its own ESM implementation, restricted pickle reader and Geneformer `vocab.txt` placeholder in favour of HF.Net 0.11.0.

## Release / Rilis

1. Move the solution into the monorepo `Vibe_Bio` as the `Biotech.Net/` subfolder (Pindahkan solusi ke monorepo `Vibe_Bio`).
2. CI builds the Rust kernels for win-x64, linux-x64 and osx-arm64, runs every test, and packs every library with the native binaries.
3. Pushing a `biotechnet-v<version>` tag publishes the `Gravicode.Biotech.Net.*` packages to nuget.org (secret `NUGET_API_KEY`).
