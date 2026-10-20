# Roadmap

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Roadmap resmi ada di [`PLAN.md`](../../PLAN.md) pada akar repositori, dan status harian dilacak di [`Progress.md`](../../Progress.md). Halaman ini merangkum posisi proyek saat ini.

## Sudah dirilis

| Tonggak | Isi |
|---|---|
| **v1 — Fondasi (0.1.x)** | kernel Rust native dengan fallback managed; Core; Sequences (FASTA/faidx, FASTQ, GenBank, EMBL, PDB, GFF/BED, ORF, motif, restriksi, Tm, penggunaan kodon); Alignments (matriks, berpasangan, MSA, filogenetika, pencarian basis data); Visualization; Notebooks; Gallery Avalonia |
| **v2 — Omics (0.2.x)** | Genomics (VCF, efek varian, genetika populasi, perakitan de Bruijn, pencarian gen, pemetaan bacaan, pemanggilan SNV); Transcriptomics (normalisasi, DE gaya DESeq2 dan limma-voom, ORA/GSEA, sel tunggal); Proteomics (ProtParam, digesti, MS/MS dengan FDR target–decoy, kuantifikasi, klasifikasi); StructuralBiology (PDB/mmCIF, geometri, superposisi, DSSP, SASA, kantong ikatan) |
| **v3 — Penemuan (0.3.x)** | Chemistry (SMILES/SDF, deskriptor, sidik jari, pencarian substruktur); DrugDiscovery (penapisan, prediksi target, ADMET, optimasi lead, QSAR); SyntheticBiology (desain dan penilaian guide, model base dan prime editing, optimasi kodon, primer, Gibson dan Golden Gate, simulasi sirkuit); AI (model zoo HF.Net, praproses offline, embedding, penilaian varian ESM, kopilot LLM) |

## Ditambahkan di 0.2.0

| Bidang | Isi |
|---|---|
| **Pengerasan statistik** | pipeline DESeq2 (faktor ukuran termasuk poscounts, dispersi, penggantian Cook, Wald/LRT, penyaringan independen) dan penyusutan LFC apeglm / normal mereproduksi R DESeq2 1.52.0 dan apeglm 1.34.0 pada pasilla dan airway — lihat [tutorial Transkriptomika](tutorials/transcriptomics.md) |
| **BAM/CRAM** | BGZF, baca/tulis SAM/BAM dengan BAI/CSI, pembacaan CRAM 3.0/3.1 dengan CRAI, tabix, kedalaman baca — diverifikasi terhadap samtools/htslib 1.24 ([tutorial Genomika](tutorials/genomics.md)) |
| **Kimia 3-D** | konformer gaya ETKDG, UFF, docking fleksibel gaya Vina dengan skor yang cocok dengan AutoDock Vina 1.2.7 ([Kimia](tutorials/chemistry.md), [Penemuan obat](tutorials/drug-discovery.md)) |
| **Kernel GPU** | `Biotech.Net.Gpu` (ILGPU: CUDA / OpenCL / CPU) untuk skor penjajaran all-vs-all dan docking berkelompok, dengan fallback CPU ([tutorial GPU](tutorials/gpu.md)) |
| **Mesin alur kerja** | `Biotech.Net.Workflows`: aturan gaya Snakemake, penjadwalan DAG, jalankan ulang inkremental ([tutorial Workflows](tutorials/workflows.md)) |
| **Konektor data awan** | `Biotech.Net.Data`: NCBI E-utilities, UniProt, RCSB PDB, AlphaFold DB, Ensembl REST dengan cache dan mode offline ([tutorial Data](tutorials/data.md)) |
| **Sel tunggal skala atlas** | baca/tulis h5ad AnnData dengan mode backed, matriks sparse out-of-core, QC / normalisasi / HVG / PCA yang cocok dengan scanpy 1.12 |
| **Gallery web** | gallery Blazor WebAssembly dan playground C# di browser ([Web Gallery](web-gallery.md)) |
| **AI** | Biotech.Net.AI kini menjalankan ESM dan Geneformer di atas HF.Net 0.11, bukan salinannya sendiri |

## Masih terbuka

- Penulisan CRAM dan blok CRAM berkompresi LZMA.
- MMFF94, rantai samping reseptor yang fleksibel, penilaian Vinardo/AD4.
- Backend Vulkan; traceback GPU dan penjajaran global.
- Tetangga terdekat aproksimasi untuk graf kNN sejuta sel; penanda (marker) out-of-core.
- Eksekutor klaster/kontainer untuk alur kerja.

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
