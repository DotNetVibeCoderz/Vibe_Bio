# Roadmap

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Roadmap resmi ada di [`PLAN.md`](../../PLAN.md) pada akar repositori, dan status harian dilacak di [`Progress.md`](../../Progress.md). Halaman ini merangkum posisi proyek saat ini.

## Sudah dirilis

| Tonggak | Isi |
|---|---|
| **v1 — Fondasi (0.1.x)** | kernel Rust native dengan fallback managed; Core; Sequences (FASTA/faidx, FASTQ, GenBank, EMBL, PDB, GFF/BED, ORF, motif, restriksi, Tm, penggunaan kodon); Alignments (matriks, berpasangan, MSA, filogenetika, pencarian basis data); Visualization; Notebooks; Gallery Avalonia |
| **v2 — Omics (0.2.x)** | Genomics (VCF, efek varian, genetika populasi, perakitan de Bruijn, pencarian gen, pemetaan bacaan, pemanggilan SNV); Transcriptomics (normalisasi, DE gaya DESeq2 dan limma-voom, ORA/GSEA, sel tunggal); Proteomics (ProtParam, digesti, MS/MS dengan FDR target–decoy, kuantifikasi, klasifikasi); StructuralBiology (PDB/mmCIF, geometri, superposisi, DSSP, SASA, kantong ikatan) |
| **v3 — Penemuan (0.3.x)** | Chemistry (SMILES/SDF, deskriptor, sidik jari, pencarian substruktur); DrugDiscovery (penapisan, prediksi target, ADMET, optimasi lead, QSAR); SyntheticBiology (desain dan penilaian guide, model base dan prime editing, optimasi kodon, primer, Gibson dan Golden Gate, simulasi sirkuit); AI (model zoo HF.Net, praproses offline, embedding, penilaian varian ESM, kopilot LLM) |

## Sedang dikerjakan / berikutnya

- **Pengerasan statistik** — kesetaraan persis dengan DESeq2 pada data publik, penyusutan LFC gaya apeglm.
- **BAM/CRAM** — pembacaan BGZF berindeks, plus tabix untuk VCF dan GFF.
- **Kernel GPU** — ILGPU atau Vulkan untuk penjajaran all-vs-all dan docking.
- **Kimia 3-D** — pembangkitan konformer gaya ETKDG dan docking fleksibel.

## Selanjutnya

- Gallery web Blazor dan playground terhosting.
- Mesin alur kerja: DAG gaya Snakemake/Nextflow dalam C#.
- Konektor data awan (NCBI E-utilities, UniProt, PDB, Ensembl REST) dengan cache.
- Sel tunggal skala atlas: matriks sparse out-of-core dan pembaca AnnData/h5ad.

## Prinsip yang dipegang roadmap ini

| | |
|---|---|
| Native | API .NET 10 yang idiomatis, span dan streaming di mana-mana |
| Cepat | kernel Rust untuk jalur kritis, dengan fallback managed yang selalu tersedia |
| Lintas platform | win-x64, linux-x64, osx-arm64 |
| Notebook-first | formatter Polyglot / .NET Interactive / Jupyter |
| Siap AI | model bio HF.Net plus kopilot LLM |
| Jujur | setiap algoritma diuji terhadap referensi terpublikasi, dan setiap halaman menyebutkan batasannya |

## Rilis

1. Solusi dipindahkan ke monorepo `Vibe_Bio` sebagai subfolder `Biotech.Net/`.
2. CI membangun kernel Rust untuk ketiga RID, menjalankan semua uji, dan memaketkan setiap pustaka beserta biner native di dalamnya.
3. Mendorong tag `biotechnet-v<versi>` menerbitkan paket `Gravicode.Biotech.Net.*` ke nuget.org.

---
[Indeks dokumentasi](index.md) · English: [roadmap](../en/roadmap.md)
