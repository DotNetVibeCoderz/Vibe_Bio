# Tutorial: Kernel GPU

Paket: `Gravicode.Biotech.Net.Gpu`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Ada dua beban kerja yang cukup mudah diparalelkan untuk dipindah ke GPU: menilai banyak sekuens terhadap banyak sekuens, dan mengevaluasi atau mencari banyak pose docking sekaligus. `Biotech.Net.Gpu` menjalankan keduanya di atas [ILGPU](https://ilgpu.net) (CUDA atau OpenCL) dan memberi jawaban yang sama dengan kode CPU — skor penjajaran sama persis bit demi bit, energi docking sampai sekitar 1e-4 kkal/mol. Tanpa GPU, setiap pemanggilan jatuh kembali ke implementasi Rust/managed yang sudah ada, sehingga program yang sama berjalan di mana saja.

Cuplikan di bawah dikompilasi dan dijalankan pada laptop yang GPU satu-satunya adalah Intel UHD Graphics 620 (OpenCL); keluaran yang ditampilkan berasal dari eksekusi itu.

## 1. Persiapan dan perangkat

```text
dotnet add package Gravicode.Biotech.Net.Gpu
```

```csharp
using Biotech.Net.Gpu;

foreach (GpuDeviceInfo info in GpuDevice.Enumerate())
    Console.WriteLine(info);

using GpuDevice? device = GpuDevice.TryCreate();   // CUDA, else OpenCL; null without a GPU or with BIOTECH_NET_GPU=off
Console.WriteLine(device is null ? "no GPU - CPU path" : $"using {device}");
```

```text
Cpu: CPUAccelerator
OpenCL: Intel(R) UHD Graphics 620 (6.4 GiB)
using OpenCL: Intel(R) UHD Graphics 620 (6.4 GiB)
```

`TryCreate()` tidak pernah memilih akselerator CPU ILGPU dengan sendirinya: tanpa GPU hasilnya `null`, dan setiap API lalu menjalankan kode CPU. Mintalah secara eksplisit dengan `GpuDevice.CreateCpu()` (atau `BIOTECH_NET_GPU=cpu`) untuk menjalankan kernel GPU di thread CPU — begitulah rangkaian uji memeriksanya di mesin tanpa GPU. `BIOTECH_NET_GPU=off` mematikan jalur GPU tanpa mengubah kode. Buang (dispose) perangkat setelah selesai; ia memiliki kernel yang sudah dikompilasi.

## 2. Skor penjajaran all-vs-all

`GpuLocalAligner` menilai setiap query terhadap setiap target dengan Smith–Waterman dan gap afin, satu thread GPU per pasangan:

```csharp
using Biotech.Net.Alignments;
using Biotech.Net.Core;
using Biotech.Net.Gpu.Alignment;
using Biotech.Net.Native;

Sequence[] proteins =
[
    Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ", "p1"),
    Sequence.Protein("MKTAYIAKQRQISFVKAHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ", "p2"),
    Sequence.Protein("MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGK", "p3"),
    Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQLEERLG", "p4"),
];
var aligner = new GpuLocalAligner(device, SubstitutionMatrix.Blosum62, gapOpen: 11, gapExtend: 1);
LocalScoreMatrix scores = aligner.ScoreAll(proteins);          // 4 x 4, one GPU thread per pair
Console.WriteLine(scores.Backend);
for (int q = 0; q < scores.QueryCount; q++)
    Console.WriteLine(string.Join(" ", scores.Row(q).ToArray().Select(s => s.ToString().PadLeft(4))));
foreach (LocalScoreHit hit in scores.Top(query: 0, k: 2))
    Console.WriteLine($"p1 -> {proteins[hit.Target].Id}: {hit.Score}");

var cpu = new PairwiseAligner { Mode = AlignmentMode.Local, SubstitutionMatrix = SubstitutionMatrix.Blosum62, OpenGapScore = -11, ExtendGapScore = -1 };
Console.WriteLine(cpu.Score(proteins[0], proteins[2]) == scores[0, 2]);   // True: same integers as the CPU kernel
```

```text
OpenCL: Intel(R) UHD Graphics 620
 321  318   23  138
 318  321   23  135
  23   23  329   18
 138  135   18  138
p1 -> p1: 321
p1 -> p2: 318
True
```

Biaya gap di sini positif, gaya BLAST (`PairwiseAligner` memakai *skor* gap negatif; `GpuLocalAligner.FromPairwiseAligner` mengonversinya). Skornya bukan sekadar mirip dengan skor CPU — bilangan bulatnya sama, karena kernel menjalankan rekurens yang sama dengan aritmetika bilangan bulat. `scores.Values` memberi seluruh matriks baris demi baris dan `TopAll(k)` daftar top-k setiap query. Skor adalah titik awal; jajarkan hit yang Anda pertahankan dengan `PairwiseAligner` untuk mendapatkan penjajarannya.

## 3. Backend GPU untuk pencarian basis data

`SequenceDatabase.Search` menilai subjek hasil seeding dengan Smith–Waterman. `SearchOptions.ScoreBackend` menerima `ILocalScoreBackend` apa pun, dan `GpuLocalAligner` adalah salah satunya:

```csharp
using Biotech.Net.Alignments.Search;

var db = new SequenceDatabase(proteins);
var hits = db.Search(proteins[3], new SearchOptions { UseSeeds = false, ScoreBackend = aligner });
foreach (var h in hits) Console.WriteLine(h);
```

```text
p1	score=138	bits=57.8	E=2.53E-14
p4	score=138	bits=57.8	E=2.53E-14
p2	score=135	bits=56.6	E=5.63E-14
p3	score=18	bits=11.5	E=2.08
```

Daftar hit identik dengan pencarian CPU, termasuk nilai E-nya. Satu query terhadap basis data kecil tidak sepadan dengan biaya peluncuran GPU; backend ini baru menguntungkan dengan ribuan subjek atau dengan `UseSeeds = false`.

## 4. Menilai ribuan pose docking

`GpuVinaScorer` mengevaluasi fungsi objektif pencarian Vina — energi reseptor dari peta grid ditambah energi internal ligan — beserta gradiennya untuk satu batch pose. Sebuah pose terdiri atas posisi akar, kuaternion orientasi, dan satu sudut per torsi:

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.Gpu.Docking;
using Biotech.Net.StructuralBiology;

VinaReceptor receptor = VinaReceptor.FromPdb("1stp_receptor.pdb");
Molecule biotin = SdfReader.ReadFile("1stp_ligand.sdf").First();
var center = new Vec3(11.12, 1.68, -10.75);                 // pocket centre
using var scorer = GpuVinaScorer.Create(biotin, receptor, center, new Vec3(22.5, 22.5, 22.5), device);
double[] poses = scorer.RandomPoses(count: 4096, seed: 1, spread: 3.0);   // 4096 x (7 + 5 torsions)
double[] gradients = new double[4096 * scorer.GradientLength];
double[] energies = scorer.Evaluate(poses, cap: 1000, gradients);
double[] reference = scorer.EvaluateCpu(poses);
Console.WriteLine($"{scorer.LastBackend}: best {energies.Min():F3}, crystal {scorer.Evaluate(scorer.InputPose())[0]:F3} kcal/mol");
Console.WriteLine($"max |GPU - CPU| = {energies.Zip(reference).Where(p => Math.Abs(p.Second) < 100).Max(p => Math.Abs(p.First - p.Second)):E1}");
```

```text
OpenCL: Intel(R) UHD Graphics 620: best -4.922, crystal -8.180 kcal/mol
max |GPU - CPU| = 1.1E-004
```

Berkas reseptor dan ligan adalah fixture 1STP di `tests/Biotech.Net.DrugDiscovery.Tests/Data`. GPU menghitung dalam presisi tunggal, sehingga hasilnya sesuai dengan mesin CPU presisi ganda (`EvaluateCpu`) sampai sekitar 1e-4 kkal/mol; pose acak kebanyakan bertabrakan dengan reseptor, itulah sebabnya pose kristal jauh lebih baik daripada yang terbaik dari 4.096 pose acak. `Optimize` menjalankan BFGS milik mesin docking pada semua pose sekaligus:

```csharp
double[] start = scorer.RandomPoses(256, seed: 2);
double[] before = scorer.Evaluate(start);
double[] after = scorer.Optimize(start, maxSteps: 20);       // in place, one BFGS per thread
Console.WriteLine($"mean {before.Average():F1} -> {after.Average():F1} kcal/mol");
```

```text
mean 111.2 -> 10.9 kcal/mol
```

## 5. Docking di GPU

`GpuVinaDocker` mengganti pencarian Monte-Carlo milik `VinaDocking.Dock` dengan beberapa ribu rantai pendek, satu per thread GPU, lalu menyerahkan pose terbaik setiap rantai ke pascapemrosesan CPU (pengelompokan, penghalusan eksak, penilaian ulang Vina). Ia menerima `DockingOptions` yang sama dan mengembalikan `DockingResult` yang sama:

```csharp
Molecule benzamidine = SdfReader.ReadFile("3ptb_ligand.sdf").First();
VinaReceptor trypsin = VinaReceptor.FromPdb("3ptb_receptor.pdb");
var docker = new GpuVinaDocker(device);
DockingResult result = docker.Dock(benzamidine, trypsin, new DockingOptions { Seed = 1 }, new GpuDockingOptions { Chains = 1024 });
Console.WriteLine($"{docker.LastBackend}, {docker.LastStepsPerChain} steps x 1024 chains in {docker.LastLaunches} launches");
foreach (DockedPose pose in result.Poses.Take(3)) Console.WriteLine(pose);
Console.WriteLine($"search {result.SearchTime.TotalSeconds:F1} s");
```

```text
OpenCL: Intel(R) UHD Graphics 620, 106 steps x 1024 chains in 7 launches
#1 -6.19 kcal/mol (rmsd l.b. 0.00, u.b. 0.00)
#2 -6.18 kcal/mol (rmsd l.b. 0.02, u.b. 1.60)
#3 -6.09 kcal/mol (rmsd l.b. 1.07, u.b. 1.93)
search 5.3 s
```

Pencarian CPU juga menemukan −6,19 kkal/mol untuk kompleks ini (AutoDock Vina: −6,20). Secara bawaan rantai-rantai berbagi anggaran pencarian CPU (heuristik langkah Vina × exhaustiveness); `StepsPerChain` menimpanya. Pencarian berjalan dalam irisan beberapa langkah per peluncuran kernel supaya tidak ada peluncuran yang memicu watchdog driver layar Windows. Eksekusi pertama dalam proses ini sudah termasuk beberapa detik kompilasi kernel; docking kedua untuk kompleks yang sama butuh sekitar 1,5 s di sini, sama dengan pencarian CPU.

## 6. Seberapa cepat?

Pada laptop pengembangan (i7-8650U, Intel UHD 620), dari [halaman benchmark](../benchmarks.md):

| Beban kerja | Jalur CPU terbaik | GPU (UHD 620) |
|---|---:|---:|
| Skor protein 64 × 128 (≈ 330 juta sel) | 34 ms (Rust, 8 thread) | 113 ms |
| Energi + gradien 4.096 pose | 8,3 ms (managed, semua inti) | 5,5 ms |
| Docking 3PTB, anggaran bawaan | 1,55 s | 1,55 s |

GPU terintegrasi tidak mengalahkan CPU empat inti yang tervektorisasi dalam pemrograman dinamis bilangan bulat. Gunakan paket ini bila tersedia GPU diskret, atau untuk menjalankan batch yang sangat besar tanpa membebani CPU.

## Batasan

- **Vulkan belum diimplementasikan**; ILGPU menyediakan CUDA dan OpenCL. Tidak ada backend Metal untuk GPU Apple.
- **Hanya diukur pada GPU terintegrasi.** Kernel belum di-benchmark pada kartu NVIDIA atau AMD diskret; percepatan yang biasa dikutip untuk Smith–Waterman GPU dan Vina-GPU tidak dibuktikan di sini.
- **Skor, bukan penjajaran.** GPU hanya menghitung skor lokal (tanpa traceback, tanpa mode global atau band).
- **Docking dalam float.** Energi sesuai dengan CPU sampai sekitar 1e-4 kkal/mol untuk pose yang sama, tetapi trajektori pencarian berbeda antarperangkat, sehingga docking dengan seed yang sama memberi hasil yang sama hanya pada perangkat yang sama. Setiap rantai menyimpan satu pose terbaik.
- **Keanehan ILGPU.** Akselerator CPU ILGPU membutuhkan tipe argumen kernel yang publik (ditandai sebagai infrastruktur), dan backend OpenCL ILGPU 1.5 salah mengompilasi sebagian konversi int ke float, yang dihindari kernel dengan literal float.

---
Sebelumnya: [Workflow](workflows.md) · [Indeks tutorial](index.md)
