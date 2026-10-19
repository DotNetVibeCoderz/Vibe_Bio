# Benchmarks

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Real measurements from `benchmarks/Biotech.Net.Benchmarks` (BenchmarkDotNet). Nothing on this page is estimated.

## How to reproduce

```bash
cd benchmarks/Biotech.Net.Benchmarks
dotnet run -c Release -p:SkipRustBuild=true -- --filter '*'          # short job, ~4 minutes
dotnet run -c Release -- --filter '*Kernel*' --full                  # BenchmarkDotNet default job, much longer
```

The kernel benchmarks call `NativeBridge` (Rust) and `ManagedKernels` (C#) directly in one process, so both back-ends are measured side by side on the same data. The managed method of each pair is the baseline, so **Ratio < 1 means Rust is faster**.

## Machine

```text
BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2)
Intel Core i7-8650U CPU 1.90GHz (Max 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores, 16 GB RAM
.NET SDK 10.0.401, [Host] .NET 10.0.12, X64 RyuJIT x86-64-v3
Job = ShortRun-InProcess, IterationCount = 3, WarmupCount = 3, LaunchCount = 1
Rust kernels: biotech-native 0.1.0 (cargo 1.97.1, release: opt-level 3, fat LTO)
Run on 2026-10-15; total run time 3 min 54 s for 31 benchmarks.
```

This is a 2017 ultrabook CPU. Treat the absolute numbers as a floor and the ratios as the interesting part. Three iterations keep the run short, so the error bars are wide — differences under about 30 % are not meaningful here.

## Rust kernels vs the managed fallback

| Benchmark | Managed | Rust | Ratio | Verdict |
|---|---:|---:|---:|---|
| Local score batch, 300 aa query × 500 targets | 295,240 µs | 6,270 µs | **0.02** | Rust ~47× faster |
| Reverse complement, 1 Mb | 625.8 µs | 80.7 µs | **0.13** | Rust ~7.8× faster |
| Kabsch RMSD, 10,000 atoms | 246.5 µs | 96.1 µs | **0.39** | Rust ~2.6× faster |
| Count k-mers (k = 8), 1 Mb | 4,698 µs | 2,014 µs | **0.44** | Rust ~2.3× faster |
| Translate, 1 Mb | 985.4 µs | 718.3 µs | **0.73** | Rust ~1.4× faster |
| Isoelectric point, 10,000 aa | 13.1 µs | 11.4 µs | 0.87 | about equal |
| Protein mass, 10,000 aa | 10.7 µs | 9.8 µs | 0.92 | about equal |
| Benjamini–Hochberg, 100,000 p-values | 10,417 µs | 11,864 µs | 1.14 | managed slightly ahead |
| Gotoh global, 2 kb × 2 kb | 42,743 µs | 51,906 µs | 1.22 | managed slightly ahead |
| Welch t-test, 20,000 × 12 | 4,610 µs | 5,739 µs | 1.24 | managed slightly ahead |

Allocation tells the same story from another angle: the Rust path allocates essentially nothing (0–458 B per operation) while the managed path allocates 1.2 MB for BH over 100 k p-values and 4.0 MB for a 2 kb × 2 kb Gotoh matrix.

What to take from this:

- **Batched local alignment is the headline.** Searching a database (`SequenceDatabase.Search`, `BioKernels.LocalScoreBatch`) is where the native kernel earns its place — nearly fifty times faster on this machine.
- **Byte-wise sequence work** (reverse complement, k-mer counting) is several times faster in Rust.
- **Small numeric kernels** (protein mass, pI) are already memory-bound in C#; the FFI call overhead cancels most of the gain.
- **Welch, BH and single large Gotoh alignments** currently come out slightly *slower* through the bridge on this CPU: the managed code is vectorised well by RyuJIT, and the marshalling and allocation pattern differ. They remain correct and allocation-free on the native side, and they are candidates for optimisation. Nothing forces you to use them — `BioKernels.ForceManaged = true` picks the faster path per workload.

## Parsers

| Benchmark | Mean | Allocated | Throughput |
|---|---:|---:|---|
| FASTA, 2,000 records × 5 kb (10 Mb of sequence) | 61.1 ms | 33.0 MB | ~164 Mb/s |
| FASTQ, 50,000 reads × 150 bp (7.5 Mb of bases) | 64.1 ms | 73.6 MB | ~780,000 reads/s |
| VCF, 20,000 records × 10 samples | 25.5 ms | 27.8 MB | ~784,000 records/s |

Input is in-memory text, so the numbers measure parsing, not disk. Throughput is derived from the mean time and the input size.

BAM decoding (`BamBenchmarks`, ShortRun in-process job): an in-memory BAM of 200,000 paired reads × 150 bp (random qualities, 34.5 MB compressed).

| Benchmark | Mean | Allocated | Throughput |
|---|---:|---:|---|
| BGZF inflate, whole file (CRC32 checked) | 230.8 ms | 2.4 MB | ~150 MB/s compressed |
| BAM decode to `AlignmentRecord` | 289.1 ms | 111.5 MB | ~690,000 records/s |

Most of the BAM time is deflate decompression; record decoding adds about 0.3 µs per read.

## Cheminformatics

| Benchmark | Mean | Allocated |
|---|---:|---:|
| Parse 20 drug SMILES | 429.4 µs | 410 KB |
| Morgan (ECFP4, r = 2, 2048 bits) × 20 molecules | 175.4 µs | 134 KB |
| Tanimoto top-10 over 10,000 fingerprints | 633.3 µs | 278 KB |
| Bulk Tanimoto over 10,000 fingerprints | 238.9 µs | 159 KB |

About 21 µs per SMILES parse and 8.8 µs per ECFP4 fingerprint; bulk Tanimoto over 10,000 fingerprints takes 0.24 ms, so comparing a query against a million precomputed fingerprints is about 24 ms single-threaded — building those fingerprints (roughly 8.8 s per million molecules) is the real cost.

## RNA-seq normalisation (20,000 genes × 12 samples)

| Benchmark | Mean | Allocated |
|---|---:|---:|
| Median-of-ratios size factors (DESeq2) | 25.1 ms | 313 KB |
| Normalized counts | 1.3 ms | 2.4 MB |
| TMM factors (edgeR) | 62.9 ms | 19.9 MB |
| log-CPM | 5.2 ms | 2.4 MB |

A full `DESeq.Run` costs considerably more than normalisation alone because of the per-gene GLM fits; the numbers above isolate the normalisation step.

## GPU kernels (Biotech.Net.Gpu, ILGPU)

Measured on 2026-10-20 on the same laptop, whose only GPU is the integrated **Intel UHD Graphics 620** (OpenCL, 24 execution units, shared memory). No NVIDIA GPU was available, so there is no CUDA number. `ILGPU-CPU` is ILGPU's CPU accelerator, which runs the GPU kernels on CPU threads for testing; it is not meant to be fast. Classes `AllVsAllBenchmarks`, `DockingPoseBenchmarks` and `DockingSearchBenchmarks`; each backend is a `Backend` parameter, and OpenCL/CUDA rows appear only when such a device is present.

```bash
dotnet run -c Release -p:SkipRustBuild=true -- --filter '*AllVsAll*' '*DockingPose*' '*DockingSearch*'
```

**All-vs-all Smith–Waterman scores**, 64 queries × 128 targets of 100–300 aa, BLOSUM62, gaps 11/1 (≈ 330 M DP cells). The CPU rows run one query per core.

| Backend | Mean | ≈ GCUPS |
|---|---:|---:|
| Rust (AVX2 striped, 8 threads) | 33.7 ms | 9.7 |
| ILGPU-OpenCL (UHD 620) | 113.0 ms | 2.9 |
| Managed C# (8 threads) | 1,008 ms | 0.33 |
| ILGPU-CPU accelerator | 3,362 ms | 0.10 |

**Vina objective + gradient for 4,096 poses** of biotin in streptavidin (1STP, 5 torsions):

| Backend | Mean |
|---|---:|
| ILGPU-OpenCL (UHD 620) | 5.5 ms |
| Managed, all cores | 8.3 ms |
| Managed, one thread | 22.6 ms |
| ILGPU-CPU accelerator | 55.0 ms |

**Whole docking run**, benzamidine into trypsin (3PTB), default budget, including grid maps and the CPU post-processing:

| Backend | Mean | Best affinity |
|---|---:|---:|
| CPU Vina search (exhaustiveness 8, all cores) | 1.55 s | −6.19 kcal/mol |
| GPU search, 1,024 chains × 106 steps (UHD 620) | 1.55 s | −6.19 kcal/mol |

Honest reading:

- **On this integrated GPU the Rust kernel wins the alignment benchmark by 3.4×.** Four Kaby Lake cores with AVX2 striping are simply faster than 24 low-clock execution units fed from shared DDR4. The GPU kernel is still 9× faster than the managed fallback. Integer DP on a discrete GPU should be one to two orders of magnitude faster than this iGPU, but that has not been measured here.
- **Pose evaluation is the one place the iGPU wins** (1.5× over managed code on all cores), because each pose is independent float work with no branches on the data path.
- **Docking breaks even.** The GPU search spends a similar number of energy evaluations as the CPU search and lands on the same best pose; on 1STP it took 7.8 s against 9.5 s on the CPU over three seeds (best affinities −7.01 to −7.21 vs −7.11 to −7.31 kcal/mol). The first GPU run in a process also pays a one-off kernel compilation of a few seconds.
- No Rust docking kernel exists, so the docking tables have no Rust row; Vulkan is not implemented.

## Raw reports

The unedited BenchmarkDotNet markdown is kept in the repository at `benchmarks/Biotech.Net.Benchmarks/results/2026-10-15-win-x64-i7-8650U/`, one file per benchmark class, including the full environment header. The GPU runs are in `results/2026-10-20-win-x64-i7-8650U-uhd620-gpu/`.

---
[Documentation index](index.md) · [Native API](api-reference/native.md) · Bahasa Indonesia: [benchmark](../id/benchmarks.md)
