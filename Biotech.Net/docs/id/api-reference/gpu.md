# Referensi API: Gpu

Paket `Gravicode.Biotech.Net.Gpu` · namespace `Biotech.Net.Gpu`, `Biotech.Net.Gpu.Alignment`, `Biotech.Net.Gpu.Docking`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Kernel GPU yang ditulis dengan [ILGPU](https://ilgpu.net) 1.5 (CUDA, OpenCL, dan akselerator CPU milik ILGPU). Bergantung pada Alignments dan DrugDiscovery; tidak ada modul lain di solusi yang bergantung padanya, jadi semua modul lain tetap berjalan tanpa GPU. Setiap API menerima `GpuDevice` yang boleh `null`: `null` menjalankan kode CPU Rust/managed yang sudah ada, dan kegagalan kernel juga jatuh kembali ke jalur itu.

## Perangkat: `Biotech.Net.Gpu`

| Tipe | Anggota |
|---|---|
| `GpuBackend` | `Auto` (CUDA, lalu OpenCL; tidak pernah akselerator CPU), `Cuda`, `OpenCL`, `Cpu` (akselerator CPU ILGPU, deterministik), `Off` |
| `GpuDeviceInfo` (record) | `Backend`, `Name`, `MemoryBytes`, `MaxThreadsPerGroup` |
| `GpuDevice` (`IDisposable`) | `Enumerate()`, `IsGpuAvailable`, `TryCreate(backend = Auto)` (null bila tidak ada, dimatikan, atau driver gagal), `CreateCpu()`, `EnvironmentVariable` = `BIOTECH_NET_GPU`, `EnvironmentBackend`, `ParseBackend(string?)`; instans `Backend`, `Info`, `Name`, `IsGpu`, `Dispose()` |

`BIOTECH_NET_GPU` (`cuda`, `opencl`, `cpu`, `off`, `auto`) berlaku saat pemanggil meminta `Auto`; backend yang disebut eksplisit selalu menang. Kernel dikompilasi sekali per perangkat lalu disimpan; peluncuran pada satu perangkat diserialkan, sehingga satu perangkat aman dipakai bersama oleh banyak thread.

## Penjajaran all-vs-all: `Biotech.Net.Gpu.Alignment`

| Tipe | Anggota |
|---|---|
| `GpuLocalAligner` (mengimplementasikan `ILocalScoreBackend`) | `(device, matrix = BLOSUM62, gapOpen = 11, gapExtend = 1)`, `FromPairwiseAligner(aligner, device)` (mode lokal, tanpa band); `ScoreAll(queries, targets)`, `ScoreAll(sequences)`, `TopHits(queries, targets, k)`, `Score(queries, queryOffsets, targets, targetOffsets, matrix128, open, extend, scores)` mentah, `LocalScoreBatch(...)`; `Device`, `Matrix`, `GapOpen`, `GapExtend`, `Name`, `LastBackend`, `LastFallbackError` |
| `LocalScoreMatrix` | `QueryCount`, `TargetCount`, `Backend`, `this[q, t]`, `Row(q)`, `Values` (baris demi baris), `ToArray()`, `Top(q, k)`, `TopAll(k)` |
| `LocalScoreHit` (record struct) | `Target`, `Score` |

Biaya gap bernilai positif (konvensi BLAST: gap sepanjang k berbiaya `open + (k − 1)·extend`). Query mengindeks baris matriks, sama seperti `SequenceDatabase`. Skornya adalah bilangan bulat yang persis sama dengan `BioKernels.LocalScore` dan `PairwiseAligner.Score` mode lokal. Bila skor top-k seri, indeks target yang lebih kecil didahulukan.

**Kernel.** Satu thread per pasangan (query, target): target diurutkan dari yang terpanjang dan dikemas dalam satu buffer, residu dikodekan ulang ke alfabet ringkas dengan matriks K×K, baris DP sepanjang query disimpan di memori global yang diselang-seling per lane (akses tergabung/coalesced), empat residu target per sapuan. Rentang pasangan dipecah menjadi beberapa peluncuran yang muat dalam anggaran scratch 128 MB (baris H dan F).

**Pencarian basis data.** `SearchOptions.ScoreBackend` (baru di Alignments, bertipe `ILocalScoreBackend`) menerima penilai apa pun; `new SearchOptions { ScoreBackend = new GpuLocalAligner(device) }` memindahkan tahap Smith–Waterman milik `SequenceDatabase.Search` ke GPU dengan hit yang identik.

## Docking: `Biotech.Net.Gpu.Docking`

| Tipe | Anggota |
|---|---|
| `GpuVinaScorer` (`IDisposable`) | `(ligand, grid, device, weights = null)`, `Create(ligand, receptor, center, size, device, spacing = 0.375, weights = null)`; `PoseLength` (7 + T), `GradientLength` (6 + T), `TorsionCount`, `Grid`, `Device`; `InputPose()`, `RandomPoses(count, seed, spread = 3)`, `Coordinates(pose)`, `Evaluate(poses, cap = 1000, gradients = default)`, `EvaluateCpu(...)` (acuan presisi ganda), `Optimize(poses, maxSteps = 20, cap = 1000)` (BFGS batch di tempat), `OptimizeCpu(...)`, `LastBackend`, `LastFallbackError` |
| `GpuVinaDocker` | `(device)`; `Dock(ligand, receptor, DockingOptions? = null, GpuDockingOptions? = null)` → `DockingResult`; `LastBackend`, `LastFallbackError`, `LastLaunches`, `LastStepsPerChain` |
| `GpuDockingOptions` (record) | `Chains` (1024), `StepsPerChain` (0 = anggaran CPU `langkah × exhaustiveness / chains`, minimal 64), `LaunchBudget` (250 ms per peluncuran) |

Sebuah pose terdiri atas posisi (Å, absolut), kuaternion (w, x, y, z), dan sudut torsi (radian) relatif terhadap geometri masukan. `Evaluate` mengembalikan fungsi objektif pencarian — energi antarmolekul dari grid beserta penalti di luar kotak ditambah energi intramolekul, keduanya dengan *curl cap* — dan gradien konformasinya (gaya, torsi putar, dE/dsudut), bukan afinitas yang dilaporkan.

**Presisi.** Perangkat menghitung dalam presisi tunggal (float) pada kerangka yang berpusat di tengah kotak. Dibandingkan mesin CPU (double) pada 3PTB dan 1STP, masing-masing 257 pose: selisih energi terbesar 7e-5 kkal/mol pada akselerator CPU dan 1,2e-4 kkal/mol pada Intel UHD 620 (OpenCL) untuk pose dengan |E| < 100, galat relatif di bawah 6e-6, gradien dalam 1,1e-4 relatif. Pengujian menegaskan `|ΔE| ≤ 1e-4 + 1e-5·|E|`.

**Pencarian.** Setiap thread menjalankan satu rantai Monte-Carlo (mutasi satu derajat kebebasan → BFGS dengan cap perburuan → Metropolis pada T = 1,2 → penghalusan setiap rekor baru dengan cap penuh), yaitu algoritme mesin CPU dalam float. Rantai berjalan dalam irisan beberapa langkah per peluncuran, ukurannya disesuaikan dengan kecepatan terukur agar satu peluncuran tetap di bawah watchdog driver layar Windows. Pose terbaik setiap rantai kemudian melewati pascapemrosesan `VinaDocking` di CPU (pengelompokan RMSD, penghalusan presisi ganda dengan interaksi reseptor eksak, penilaian ulang Vina), sehingga afinitas yang dilaporkan adalah nilai CPU untuk pose-pose itu.

## Infrastruktur kernel

`LocalScoreKernel`, `LocalScoreArgs`, `VinaKernels`, `VinaKernelData`, `VinaScratch`, `VinaEvalArgs`, `VinaMcArgs`, dan `VinaLayout` publik hanya karena akselerator CPU ILGPU menulis peluncurnya ke assembly dinamis yang tidak dapat melihat tipe internal. Tipe-tipe ini disembunyikan dari IntelliSense dan bukan API yang didukung.

## Batasan

- **Tidak ada backend Vulkan.** ILGPU menargetkan CUDA dan OpenCL (ditambah akselerator CPU-nya); jalur Vulkan/SPIR-V belum diimplementasikan. Metal (GPU Apple) juga tidak tersedia.
- **GPU terintegrasi lambat.** Pada Intel UHD 620 yang dipakai saat pengembangan, kernel all-vs-all lebih lambat daripada kernel Rust AVX2 di empat inti CPU; lihat [benchmark](../benchmarks.md). Kernel ditulis untuk GPU diskret dan belum diukur pada GPU semacam itu.
- **Hanya skor lokal.** Tidak ada traceback di GPU, mode global atau semi-global, maupun band; jajarkan hit teratas dengan `PairwiseAligner` di CPU.
- **Trajektori docking khas perangkat.** Hasil dapat diulang untuk satu seed pada satu perangkat, tetapi pembulatan float berbeda antara akselerator CPU, driver OpenCL, dan CUDA, sehingga rantai — dan kadang pose peringkat bawah — berbeda.
- **Hanya pose terbaik per rantai.** Setiap rantai melaporkan satu pose terbaik; keragaman pose berasal dari banyaknya rantai.

---
[Indeks referensi API](index.md) · [Tutorial GPU](../tutorials/gpu.md)
