# Dokumentasi Biotech.Net (Bahasa Indonesia)

**Ultimate Biotechnology & Bioinformatics Platform for .NET** — Bioconductor + Biopython-nya .NET 10, dengan kernel yang dipercepat Rust.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

English: [documentation](../en/index.md)

## Daftar isi

| Halaman | Isi |
|---|---|
| [Pendahuluan](getting-started.md) | apa itu Biotech.Net, tur lima menit, konvensi yang dipakai |
| [Instalasi](installation.md) | kebutuhan, daftar paket, membangun dari sumber, kernel native, CI |
| [Tutorial](tutorials/index.md) | satu per modul, setiap cuplikan kode dikompilasi dan dijalankan |
| [Referensi API](api-reference/index.md) | tipe dan metode utama per modul |
| [Contoh](samples.md) | aplikasi Gallery (Avalonia) dan notebook contoh |
| [Benchmark](benchmarks.md) | angka nyata Rust vs managed, parser, kimia, normalisasi |
| [Roadmap](roadmap.md) | yang sudah selesai dan yang berikutnya |

Nama berkas di `docs/id/` sengaja sama persis dengan `docs/en/` supaya kedua pohon dokumentasi bercermin satu lawan satu (`getting-started.md` = Pendahuluan, `installation.md` = Instalasi, `tutorials/` = Tutorial, `api-reference/` = Referensi API, `samples.md` = Contoh, `benchmarks.md` = Benchmark).

## Tutorial sekilas

[Sekuens](tutorials/sequences.md) · [Penjajaran](tutorials/alignments.md) · [Genomika](tutorials/genomics.md) · [Transkriptomika](tutorials/transcriptomics.md) · [Proteomika](tutorials/proteomics.md) · [Biologi struktural](tutorials/structural-biology.md) · [Kimia](tutorials/chemistry.md) · [Penemuan obat](tutorials/drug-discovery.md) · [Biologi sintetik](tutorials/synthetic-biology.md) · [AI](tutorials/ai.md) · [Visualisasi](tutorials/visualization.md) · [Notebook](tutorials/notebooks.md)

## Aturan main dokumentasi ini

- Setiap cuplikan C# diekstrak ke proyek sementara, dikompilasi, lalu dijalankan terhadap API yang sebenarnya sebelum diterbitkan. Kalau sebuah cuplikan mencetak nilai, nilai itu berasal dari eksekusi nyata.
- Setiap gambar dihasilkan oleh `Biotech.Net.Visualization` sendiri.
- Setiap halaman diakhiri bagian **Batasan** yang jujur. Heuristik disebut heuristik, data demo disebut demo, dan model yang tidak didukung dinyatakan tidak didukung.
- Pohon `docs/en/` dan `docs/id/` bercermin berkas per berkas.

## Lisensi

MIT. Hak cipta © Gravicode Studios.
