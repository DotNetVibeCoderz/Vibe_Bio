# API reference: Gpu

Package `Gravicode.Biotech.Net.Gpu` · namespaces `Biotech.Net.Gpu`, `Biotech.Net.Gpu.Alignment`, `Biotech.Net.Gpu.Docking`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

GPU kernels written with [ILGPU](https://ilgpu.net) 1.5 (CUDA, OpenCL and ILGPU's CPU accelerator). Depends on Alignments and DrugDiscovery; nothing in the solution depends on it, so every other module keeps working without a GPU. Every API takes a nullable `GpuDevice`: `null` runs the existing Rust/managed CPU code, and a kernel failure falls back to it as well.

## Devices: `Biotech.Net.Gpu`

| Type | Members |
|---|---|
| `GpuBackend` | `Auto` (CUDA, then OpenCL; never the CPU accelerator), `Cuda`, `OpenCL`, `Cpu` (ILGPU's CPU accelerator, deterministic), `Off` |
| `GpuDeviceInfo` (record) | `Backend`, `Name`, `MemoryBytes`, `MaxThreadsPerGroup` |
| `GpuDevice` (`IDisposable`) | `Enumerate()`, `IsGpuAvailable`, `TryCreate(backend = Auto)` (null when absent, off or the driver fails), `CreateCpu()`, `EnvironmentVariable` = `BIOTECH_NET_GPU`, `EnvironmentBackend`, `ParseBackend(string?)`; instance `Backend`, `Info`, `Name`, `IsGpu`, `Dispose()` |

`BIOTECH_NET_GPU` (`cuda`, `opencl`, `cpu`, `off`, `auto`) applies when the caller asks for `Auto`; an explicit backend wins. Kernels are compiled once per device and cached; launches on one device are serialised, so a device can be shared between threads.

## All-vs-all alignment: `Biotech.Net.Gpu.Alignment`

| Type | Members |
|---|---|
| `GpuLocalAligner` (implements `ILocalScoreBackend`) | `(device, matrix = BLOSUM62, gapOpen = 11, gapExtend = 1)`, `FromPairwiseAligner(aligner, device)` (local mode, no band); `ScoreAll(queries, targets)`, `ScoreAll(sequences)`, `TopHits(queries, targets, k)`, raw `Score(queries, queryOffsets, targets, targetOffsets, matrix128, open, extend, scores)`, `LocalScoreBatch(...)`; `Device`, `Matrix`, `GapOpen`, `GapExtend`, `Name`, `LastBackend`, `LastFallbackError` |
| `LocalScoreMatrix` | `QueryCount`, `TargetCount`, `Backend`, `this[q, t]`, `Row(q)`, `Values` (row-major), `ToArray()`, `Top(q, k)`, `TopAll(k)` |
| `LocalScoreHit` (record struct) | `Target`, `Score` |

Gap costs are positive (BLAST convention: a gap of length k costs `open + (k − 1)·extend`). The query indexes the matrix rows, as in `SequenceDatabase`. Scores are the exact integers of `BioKernels.LocalScore` and `PairwiseAligner.Score` in local mode. Top-k ties go to the lower target index.

**Kernel.** One thread per (query, target) pair: targets sorted longest first and packed in one buffer, residues re-coded to a compact alphabet with a K×K matrix, the DP row along the query kept in lane-interleaved global scratch (coalesced), four target residues per sweep. The pair range is cut into launches that fit a 128 MB scratch budget (H and F rows).

**Database search.** `SearchOptions.ScoreBackend` (new in Alignments, type `ILocalScoreBackend`) takes any scorer; `new SearchOptions { ScoreBackend = new GpuLocalAligner(device) }` moves the Smith–Waterman stage of `SequenceDatabase.Search` to the GPU with identical hits.

## Docking: `Biotech.Net.Gpu.Docking`

| Type | Members |
|---|---|
| `GpuVinaScorer` (`IDisposable`) | `(ligand, grid, device, weights = null)`, `Create(ligand, receptor, center, size, device, spacing = 0.375, weights = null)`; `PoseLength` (7 + T), `GradientLength` (6 + T), `TorsionCount`, `Grid`, `Device`; `InputPose()`, `RandomPoses(count, seed, spread = 3)`, `Coordinates(pose)`, `Evaluate(poses, cap = 1000, gradients = default)`, `EvaluateCpu(...)` (double-precision reference), `Optimize(poses, maxSteps = 20, cap = 1000)` (batched BFGS in place), `OptimizeCpu(...)`, `LastBackend`, `LastFallbackError` |
| `GpuVinaDocker` | `(device)`; `Dock(ligand, receptor, DockingOptions? = null, GpuDockingOptions? = null)` → `DockingResult`; `LastBackend`, `LastFallbackError`, `LastLaunches`, `LastStepsPerChain` |
| `GpuDockingOptions` (record) | `Chains` (1024), `StepsPerChain` (0 = CPU budget `steps × exhaustiveness / chains`, at least 64), `LaunchBudget` (250 ms per launch) |

A pose is position (Å, absolute), quaternion (w, x, y, z) and torsion angles (radians) relative to the input geometry. `Evaluate` returns the search objective — grid intermolecular energy with the out-of-box penalty plus the intramolecular energy, both curl-capped — and its conformation gradient (force, torque, dE/dangle), not the reported affinity.

**Precision.** The device works in single precision in a frame centred on the box. Against the CPU engine (double) on 3PTB and 1STP, for 257 poses each: largest energy difference 7e-5 kcal/mol on the CPU accelerator and 1.2e-4 kcal/mol on an Intel UHD 620 (OpenCL) for poses with |E| < 100, relative error below 6e-6, gradients within 1.1e-4 relative. The tests assert `|ΔE| ≤ 1e-4 + 1e-5·|E|`.

**Search.** Each thread runs one Monte-Carlo chain (mutate one degree of freedom → BFGS with the hunting cap → Metropolis at T = 1.2 → refine each new best with the full cap), the CPU engine's algorithm in float. Chains run in slices of a few steps per launch, sized from the measured speed so a launch stays under the Windows display-driver watchdog. The best pose of every chain then goes through `VinaDocking`'s post-processing on the CPU (RMSD clustering, double-precision refinement with exact receptor interactions, Vina rescoring), so reported affinities are CPU values for those poses.

## Kernel infrastructure

`LocalScoreKernel`, `LocalScoreArgs`, `VinaKernels`, `VinaKernelData`, `VinaScratch`, `VinaEvalArgs`, `VinaMcArgs` and `VinaLayout` are public only because ILGPU's CPU accelerator emits its launchers into a dynamic assembly that cannot see internal types. They are hidden from IntelliSense and not a supported API.

## Limitations

- **No Vulkan backend.** ILGPU targets CUDA and OpenCL (plus its CPU accelerator); a Vulkan/SPIR-V path is not implemented. Metal (Apple GPUs) is not available either.
- **Integrated GPUs are slow.** On the Intel UHD 620 used for development the all-vs-all kernel is slower than the Rust AVX2 kernel on four CPU cores; see [benchmarks](../benchmarks.md). The kernels are written for discrete GPUs and have not been measured on one.
- **Local scores only.** No GPU traceback, global or semi-global mode, or banding; align the reported top hits with `PairwiseAligner` on the CPU.
- **Docking trajectories are device-specific.** Results are reproducible for a seed on one device, but float rounding differs between the CPU accelerator, OpenCL drivers and CUDA, so the chains, and sometimes the lower-ranked poses, differ.
- **Per-chain best only.** Each chain reports a single best pose; pose diversity comes from the number of chains.

---
[API reference index](index.md) · [GPU tutorial](../tutorials/gpu.md)
