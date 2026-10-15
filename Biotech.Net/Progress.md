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
| Genomics | ✅ | ✅ | VCF, variant effects + HGVS, popgen, De Bruijn assembly, gene finding, read mapping, SNV calling |
| Transcriptomics | ✅ | ✅ | Normalisation, DESeq2-style and limma-voom-style DE, ORA/GSEA, single-cell |
| Proteomics | ✅ | ✅ | ProtParam, digestion, MS/MS with target-decoy FDR, quantification, classification |
| StructuralBiology | ✅ | ✅ | PDB/mmCIF, geometry, superposition, DSSP-style, SASA, pockets |
| Chemistry | ✅ | ✅ | Verified against RDKit 2026.09.1; fingerprints, SDF/MOL, 2-D depiction, scaffolds |
| DrugDiscovery | ✅ | ✅ | Screening, Vina-like scoring, SEA-style target prediction (demo set), ADMET, QSAR |
| SyntheticBiology | ✅ | ✅ | CRISPR design + off-target, base/prime editing, primers, Gibson, Golden Gate |
| AI | ✅ | ✅ | HF.Net model zoo; managed ESM implementation, Geneformer, LLM copilot |
| Visualization | ✅ | ✅ | Genome browser, heatmap, volcano/MA, PCA/UMAP, sequence logo, dot plot, protein and molecule viewers |
| Notebooks | ✅ | ✅ | Formatters + `extension.dib`; six sample `.dib` notebooks |
| Gallery (Avalonia) | ✅ | — | 7 pages + live code; `--screenshots docs/screenshots`; desktop window verified |
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
