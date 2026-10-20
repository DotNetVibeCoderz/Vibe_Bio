# Biotech.Net — Progress / Kemajuan

Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.

Tracks development against [PLAN.md](PLAN.md). Status: ✅ done · 🚧 in progress · ⏳ not started.

## Modules / Modul

| Module | Source | Tests | Notes |
|---|---|---|---|
| Native (Rust + managed fallback) | ✅ | ✅ | Rust kernels + real P/Invoke bridge; managed fallback for every kernel (parity-tested); `BIOTECH_NET_DISABLE_NATIVE=1` forces managed |
| Core | ✅ | ✅ |  |
| Sequences | ✅ | ✅ | FASTA/faidx, FASTQ, GenBank, EMBL, PDB, ORFs, motifs, restriction, Tm, codon usage |
| Alignments | ✅ | ✅ | Substitution matrices, pairwise (Gotoh), progressive MSA, UPGMA/NJ/Newick, BLAST-like search |
| Genomics | ✅ | ✅ | VCF, variant effects + HGVS, popgen, De Bruijn assembly, gene finding, read mapping, SNV calling; BGZF, SAM/BAM, BAI/CSI, CRAM 3.0/3.1 reading, tabix (vs samtools 1.24) |
| Transcriptomics | ✅ | ✅ | Normalisation, DESeq2 + apeglm (parity with R DESeq2 1.52.0), limma-voom-style DE, ORA/GSEA, single-cell incl. h5ad and out-of-core (vs scanpy 1.12) |
| Proteomics | ✅ | ✅ | ProtParam, digestion, MS/MS with target-decoy FDR, quantification, classification |
| StructuralBiology | ✅ | ✅ | PDB/mmCIF, geometry, superposition, DSSP-style, SASA, pockets |
| Chemistry | ✅ | ✅ | Verified against RDKit 2026.09.1; fingerprints, SDF/MOL, 2-D depiction, scaffolds, ETKDG-style 3-D conformers, UFF |
| DrugDiscovery | ✅ | ✅ | Screening, Vina-style flexible docking (vs AutoDock Vina 1.2.7), SEA-style target prediction (demo set), ADMET, QSAR |
| SyntheticBiology | ✅ | ✅ | CRISPR design + off-target, base/prime editing, primers, Gibson, Golden Gate |
| AI | ✅ | ✅ | HF.Net model zoo; ESM and Geneformer on HF.Net 0.11, LLM copilot |
| Visualization | ✅ | ✅ | Genome browser, heatmap, volcano/MA, PCA/UMAP, sequence logo, dot plot, protein and molecule viewers |
| Notebooks | ✅ | ✅ | Formatters + `extension.dib`; six sample `.dib` notebooks |
| Data | ✅ | ✅ | NCBI E-utilities, UniProt, RCSB PDB, AlphaFold DB, Ensembl REST; cache + offline mode |
| Workflows | ✅ | ✅ | Snakemake-style DAG engine, incremental re-runs |
| Gpu | ✅ | ✅ | ILGPU all-vs-all alignment and docking, CPU fallback |
| Gallery (Avalonia) | ✅ | — | 7 pages + live code; `--screenshots docs/screenshots`; desktop window verified |
| Web Gallery (Blazor WASM) | ✅ | — | 7 pages + Roslyn-in-WASM C# playground; Pages workflow ready, not enabled |
| Docs (EN/ID) | ✅ | — | `docs/en` + `docs/id` mirrored; every C# snippet compiled and run; figures rendered by Visualization |
| Benchmarks | ✅ | — | `benchmarks/Biotech.Net.Benchmarks` (BenchmarkDotNet); results in `docs/en/benchmarks.md` |
| CI / publish workflows | ✅ | — | `build/workflows/`, copied to the monorepo root at release |

## Log

- **2026-10-13** — Solution scaffolded; Rust crate; Core, Sequences, Proteomics, StructuralBiology, Transcriptomics, SyntheticBiology and Chemistry written (session interrupted).
- **2026-10-14** — Build fixes (Chemistry ref-local); remaining modules and test suites under way; CI/publish workflows written.
- **2026-10-15** — Alignments, Genomics, Chemistry, DrugDiscovery, SyntheticBiology, AI, Visualization, Notebooks and the Avalonia Gallery completed; real P/Invoke bridge to the Rust kernels; upstream ESM support, NumPy-scalar pickles and tokenizer-less models contributed to HF.Net (0.11.0, unpublished).
- **2026-10-15** — Benchmark suite (kernels, parsers, chemistry, normalisation) run and recorded; bilingual documentation (EN/ID) written and verified; root README; docs figures rendered.
- **2026-10-15** — Sequences verified against Biopython 1.88: GenBank/EMBL readers now upper-case the sequence (as Biopython does), EMBL sequence lines fixed to 80 columns, CAI now counts stop codons; full suite green (1242 passed, 1 skipped).
- **2026-10-15** — HF.Net 0.11.0 published to nuget.org; `Biotech.Net.AI` moved to it. Solution moved into the `Vibe_Bio` monorepo (`Biotech.Net/`), workflows installed at the monorepo root, first release `biotechnet-v0.1.0` published to nuget.org (14 packages). CI on Linux caught two gaps fixed before tagging: `*.pdb` ignore rule hid Protein Data Bank files, and sequence logos drew rectangles without system fonts (now embedded A–Z outlines).
- **2026-10-18 – 2026-10-20** — Remaining PLAN items done: DESeq2/apeglm parity with R, BAM/CRAM/tabix, 3-D conformers + flexible docking, cloud data connectors, workflow engine, GPU kernels, atlas-scale single cell (h5ad), Blazor web gallery + playground; AI moved onto HF.Net 0.11. Foundation packages raised to Visual.Net 1.2.0 / Gravicode.Science 1.3.0. Full suite green (1897 passed, 5 skipped). HF.Net 0.11.1 published with the ESM rotary fix (fp16 `inv_freq`); `biotechnet-v0.2.0` released (17 packages); web gallery deployed to GitHub Pages.
