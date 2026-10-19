# Benchmark

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Pengukuran nyata dari `benchmarks/Biotech.Net.Benchmarks` (BenchmarkDotNet). Tidak ada angka perkiraan di halaman ini.

## Cara mengulanginya

```bash
cd benchmarks/Biotech.Net.Benchmarks
dotnet run -c Release -p:SkipRustBuild=true -- --filter '*'          # job singkat, ~4 menit
dotnet run -c Release -- --filter '*Kernel*' --full                  # job bawaan BenchmarkDotNet, jauh lebih lama
```

Benchmark kernel memanggil `NativeBridge` (Rust) dan `ManagedKernels` (C#) langsung dalam satu proses, sehingga kedua backend diukur berdampingan pada data yang sama. Metode managed di tiap pasangan menjadi baseline, jadi **Ratio < 1 berarti Rust lebih cepat**.

## Mesin uji

```text
BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2)
Intel Core i7-8650U CPU 1.90GHz (maks 2.11GHz) (Kaby Lake R), 1 CPU, 8 logis dan 4 inti fisik, RAM 16 GB
.NET SDK 10.0.401, [Host] .NET 10.0.12, X64 RyuJIT x86-64-v3
Job = ShortRun-InProcess, IterationCount = 3, WarmupCount = 3, LaunchCount = 1
Kernel Rust: biotech-native 0.1.0 (cargo 1.97.1, release: opt-level 3, fat LTO)
Dijalankan 2026-10-15; total waktu 3 menit 54 detik untuk 31 benchmark.
```

Ini CPU ultrabook tahun 2017. Anggap angka absolutnya sebagai batas bawah, dan rasionya sebagai bagian yang menarik. Tiga iterasi membuat eksekusinya singkat sehingga rentang galatnya lebar — perbedaan di bawah sekitar 30 % tidak berarti di sini.

## Kernel Rust vs fallback managed

| Benchmark | Managed | Rust | Rasio | Kesimpulan |
|---|---:|---:|---:|---|
| Skor lokal berkelompok, kueri 300 aa × 500 target | 295.240 µs | 6.270 µs | **0,02** | Rust ~47× lebih cepat |
| Komplemen balik, 1 Mb | 625,8 µs | 80,7 µs | **0,13** | Rust ~7,8× lebih cepat |
| Kabsch RMSD, 10.000 atom | 246,5 µs | 96,1 µs | **0,39** | Rust ~2,6× lebih cepat |
| Hitung k-mer (k = 8), 1 Mb | 4.698 µs | 2.014 µs | **0,44** | Rust ~2,3× lebih cepat |
| Translasi, 1 Mb | 985,4 µs | 718,3 µs | **0,73** | Rust ~1,4× lebih cepat |
| Titik isoelektrik, 10.000 aa | 13,1 µs | 11,4 µs | 0,87 | setara |
| Massa protein, 10.000 aa | 10,7 µs | 9,8 µs | 0,92 | setara |
| Benjamini–Hochberg, 100.000 nilai p | 10.417 µs | 11.864 µs | 1,14 | managed sedikit unggul |
| Gotoh global, 2 kb × 2 kb | 42.743 µs | 51.906 µs | 1,22 | managed sedikit unggul |
| Uji t Welch, 20.000 × 12 | 4.610 µs | 5.739 µs | 1,24 | managed sedikit unggul |

Alokasi memori menceritakan hal yang sama dari sudut lain: jalur Rust nyaris tidak mengalokasikan apa pun (0–458 B per operasi), sementara jalur managed mengalokasikan 1,2 MB untuk BH atas 100 ribu nilai p dan 4,0 MB untuk matriks Gotoh 2 kb × 2 kb.

Yang bisa disimpulkan:

- **Penjajaran lokal berkelompok adalah sorotan utamanya.** Pencarian basis data (`SequenceDatabase.Search`, `BioKernels.LocalScoreBatch`) adalah tempat kernel native membuktikan diri — hampir lima puluh kali lebih cepat di mesin ini.
- **Pekerjaan sekuens per byte** (komplemen balik, penghitungan k-mer) beberapa kali lebih cepat di Rust.
- **Kernel numerik kecil** (massa protein, pI) sudah terbatas memori di C#; overhead pemanggilan FFI menghapus sebagian besar keuntungannya.
- **Welch, BH, dan satu penjajaran Gotoh besar** saat ini justru sedikit *lebih lambat* lewat jembatan di CPU ini: kode managed divektorisasi dengan baik oleh RyuJIT, dan pola marshalling serta alokasinya berbeda. Ketiganya tetap benar dan bebas alokasi di sisi native, dan menjadi kandidat optimasi. Tidak ada yang memaksa Anda memakainya — `BioKernels.ForceManaged = true` memilih jalur tercepat per beban kerja.

## Parser

| Benchmark | Rerata | Alokasi | Throughput |
|---|---:|---:|---|
| FASTA, 2.000 rekaman × 5 kb (10 Mb sekuens) | 61,1 ms | 33,0 MB | ~164 Mb/d |
| FASTQ, 50.000 bacaan × 150 bp (7,5 Mb basa) | 64,1 ms | 73,6 MB | ~780.000 bacaan/d |
| VCF, 20.000 rekaman × 10 sampel | 25,5 ms | 27,8 MB | ~784.000 rekaman/d |

Masukannya teks di memori, jadi angka ini mengukur parsing, bukan disk. Throughput diturunkan dari waktu rerata dan ukuran masukan.

Dekode BAM (`BamBenchmarks`, job ShortRun in-process): BAM di memori berisi 200.000 bacaan berpasangan × 150 bp (kualitas acak, 34,5 MB terkompresi).

| Benchmark | Rerata | Alokasi | Throughput |
|---|---:|---:|---|
| Inflate BGZF, seluruh berkas (CRC32 diperiksa) | 230,8 ms | 2,4 MB | ~150 MB/d terkompresi |
| Dekode BAM menjadi `AlignmentRecord` | 289,1 ms | 111,5 MB | ~690.000 rekaman/d |

Sebagian besar waktu BAM habis untuk dekompresi deflate; dekode rekaman hanya menambah sekitar 0,3 µs per bacaan.

## Kemoinformatika

| Benchmark | Rerata | Alokasi |
|---|---:|---:|
| Parsing 20 SMILES obat | 429,4 µs | 410 KB |
| Morgan (ECFP4, r = 2, 2048 bit) × 20 molekul | 175,4 µs | 134 KB |
| Tanimoto top-10 atas 10.000 sidik jari | 633,3 µs | 278 KB |
| Tanimoto massal atas 10.000 sidik jari | 238,9 µs | 159 KB |

Sekitar 21 µs per parsing SMILES dan 8,8 µs per sidik jari ECFP4; Tanimoto massal atas 10.000 sidik jari memakan 0,24 ms, jadi membandingkan satu kueri terhadap sejuta sidik jari yang sudah dihitung memerlukan sekitar 24 ms satu utas — membangun sidik jarinya (sekitar 8,8 detik per juta molekul) justru biaya sesungguhnya.

## Normalisasi RNA-seq (20.000 gen × 12 sampel)

| Benchmark | Rerata | Alokasi |
|---|---:|---:|
| Size factor median-of-ratios (DESeq2) | 25,1 ms | 313 KB |
| Normalized counts | 1,3 ms | 2,4 MB |
| Faktor TMM (edgeR) | 62,9 ms | 19,9 MB |
| log-CPM | 5,2 ms | 2,4 MB |

`DESeq.Run` penuh jauh lebih mahal daripada normalisasi saja karena pencocokan GLM per gen; angka di atas hanya mengisolasi tahap normalisasi.

## Kernel GPU (Biotech.Net.Gpu, ILGPU)

Diukur pada 2026-10-20 di laptop yang sama, yang GPU satu-satunya adalah **Intel UHD Graphics 620** terintegrasi (OpenCL, 24 execution unit, memori berbagi). Tidak ada GPU NVIDIA, jadi tidak ada angka CUDA. `ILGPU-CPU` adalah akselerator CPU milik ILGPU yang menjalankan kernel GPU di thread CPU untuk pengujian; ia memang tidak dimaksudkan cepat. Kelas `AllVsAllBenchmarks`, `DockingPoseBenchmarks`, dan `DockingSearchBenchmarks`; setiap backend adalah parameter `Backend`, dan baris OpenCL/CUDA hanya muncul bila perangkatnya ada.

```bash
dotnet run -c Release -p:SkipRustBuild=true -- --filter '*AllVsAll*' '*DockingPose*' '*DockingSearch*'
```

**Skor Smith–Waterman all-vs-all**, 64 query × 128 target sepanjang 100–300 aa, BLOSUM62, gap 11/1 (≈ 330 juta sel DP). Baris CPU menjalankan satu query per inti.

| Backend | Rerata | ≈ GCUPS |
|---|---:|---:|
| Rust (AVX2 striped, 8 thread) | 33,7 ms | 9,7 |
| ILGPU-OpenCL (UHD 620) | 113,0 ms | 2,9 |
| C# managed (8 thread) | 1.008 ms | 0,33 |
| Akselerator ILGPU-CPU | 3.362 ms | 0,10 |

**Objektif Vina + gradien untuk 4.096 pose** biotin di streptavidin (1STP, 5 torsi):

| Backend | Rerata |
|---|---:|
| ILGPU-OpenCL (UHD 620) | 5,5 ms |
| Managed, semua inti | 8,3 ms |
| Managed, satu thread | 22,6 ms |
| Akselerator ILGPU-CPU | 55,0 ms |

**Satu proses docking utuh**, benzamidin ke tripsin (3PTB), anggaran bawaan, termasuk peta grid dan pascapemrosesan di CPU:

| Backend | Rerata | Afinitas terbaik |
|---|---:|---:|
| Pencarian Vina CPU (exhaustiveness 8, semua inti) | 1,55 s | −6,19 kkal/mol |
| Pencarian GPU, 1.024 rantai × 106 langkah (UHD 620) | 1,55 s | −6,19 kkal/mol |

Bacaan jujurnya:

- **Pada GPU terintegrasi ini kernel Rust memenangkan benchmark penjajaran 3,4×.** Empat inti Kaby Lake dengan AVX2 striped memang lebih cepat daripada 24 execution unit berclock rendah yang diumpan dari DDR4 bersama. Kernel GPU tetap 9× lebih cepat daripada fallback managed. DP bilangan bulat pada GPU diskret semestinya satu sampai dua orde besaran lebih cepat daripada iGPU ini, tetapi itu belum diukur di sini.
- **Evaluasi pose satu-satunya tempat iGPU menang** (1,5× atas kode managed di semua inti), karena setiap pose adalah kerja float independen tanpa percabangan pada jalur datanya.
- **Docking impas.** Pencarian GPU memakai jumlah evaluasi energi yang mirip dengan pencarian CPU dan mendarat di pose terbaik yang sama; pada 1STP ia butuh 7,8 s berbanding 9,5 s di CPU untuk tiga seed (afinitas terbaik −7,01 s.d. −7,21 vs −7,11 s.d. −7,31 kkal/mol). Eksekusi GPU pertama dalam satu proses juga membayar kompilasi kernel sekali jalan selama beberapa detik.
- Tidak ada kernel docking Rust, jadi tabel docking tidak punya baris Rust; Vulkan belum diimplementasikan.

## Laporan mentah

Markdown BenchmarkDotNet yang belum disunting disimpan di repositori pada `benchmarks/Biotech.Net.Benchmarks/results/2026-10-15-win-x64-i7-8650U/`, satu berkas per kelas benchmark, lengkap dengan header lingkungannya. Hasil GPU ada di `results/2026-10-20-win-x64-i7-8650U-uhd620-gpu/`.

---
[Indeks dokumentasi](index.md) · [API Native](api-reference/native.md) · English: [benchmarks](../en/benchmarks.md)
