# Referensi API

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Ikhtisar tipe dan metode utama tiap modul — apa yang tersedia dan di mana letaknya, bukan dump seluruh anggota. Contoh yang bisa dijalankan ada di [tutorial](../tutorials/index.md); untuk tanda tangan persis gunakan IntelliSense atau dokumentasi XML yang ikut dalam setiap paket.

| Modul | Akar namespace | Halaman |
|---|---|---|
| Native | `Biotech.Net.Native` | [native.md](native.md) |
| Core | `Biotech.Net.Core` | [core.md](core.md) |
| Workflows | `Biotech.Net.Workflows` | [workflows.md](workflows.md) |
| Sequences | `Biotech.Net.Sequences` | [sequences.md](sequences.md) |
| Alignments | `Biotech.Net.Alignments` | [alignments.md](alignments.md) |
| Genomics | `Biotech.Net.Genomics` | [genomics.md](genomics.md) |
| Transcriptomics | `Biotech.Net.Transcriptomics` | [transcriptomics.md](transcriptomics.md) |
| Proteomics | `Biotech.Net.Proteomics` | [proteomics.md](proteomics.md) |
| StructuralBiology | `Biotech.Net.StructuralBiology` | [structural-biology.md](structural-biology.md) |
| Chemistry | `Biotech.Net.Chemistry` | [chemistry.md](chemistry.md) |
| DrugDiscovery | `Biotech.Net.DrugDiscovery` | [drug-discovery.md](drug-discovery.md) |
| SyntheticBiology | `Biotech.Net.SyntheticBiology` | [synthetic-biology.md](synthetic-biology.md) |
| Data | `Biotech.Net.Data` | [data.md](data.md) |
| GPU | `Biotech.Net.Gpu` | [gpu.md](gpu.md) |
| AI | `Biotech.Net.AI` | [ai.md](ai.md) |
| Visualization | `Biotech.Net.Visualization` | [visualization.md](visualization.md) |
| Notebooks | `Biotech.Net.Notebooks` | [notebooks.md](notebooks.md) |

## Konvensi yang berlaku di semua modul

- **Koordinat.** `GenomicInterval` dan `FeatureLocation` berbasis 0 dan setengah terbuka; `Variant.Position` berbasis 1 seperti VCF.
- **Span.** Jalur kritis menerima `ReadOnlySpan<byte>`/`Span<T>` dan menulis ke buffer milik pemanggil.
- **Streaming.** Pembaca berkas mengembalikan `IEnumerable<T>` dan tidak pernah memuat seluruh berkas.
- **Objek opsi.** Daftar parameter panjang dibungkus tipe `record` dengan properti `init` dan nilai bawaan yang masuk akal (`OrfOptions`, `DESeqOptions`, `SearchParameters`, `MapperOptions`, …). Sebagian besar menerima `null` yang berarti "pakai bawaan".
- **Besar-kecil huruf.** Data sekuens dipertahankan apa adanya; perbandingan dan kernel mengabaikannya.

---
[Indeks dokumentasi](../index.md) · English: [api-reference](../../en/api-reference/index.md)
