# Tutorial: GPU kernels

Package: `Gravicode.Biotech.Net.Gpu`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Two workloads parallelise well enough to put on a GPU: scoring many sequences against many sequences, and evaluating or searching many docking poses at once. `Biotech.Net.Gpu` runs both on [ILGPU](https://ilgpu.net) (CUDA or OpenCL) and gives the same answers as the CPU code — the alignment scores bit for bit, the docking energies to about 1e-4 kcal/mol. Without a GPU every call falls back to the existing Rust/managed implementation, so the same program runs everywhere.

The snippets below were compiled and run on a laptop whose only GPU is an Intel UHD Graphics 620 (OpenCL); the outputs are from that run.

## 1. Set-up and devices

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

`TryCreate()` never picks ILGPU's CPU accelerator by itself: no GPU means `null`, and every API then runs the CPU code. Ask for it explicitly with `GpuDevice.CreateCpu()` (or `BIOTECH_NET_GPU=cpu`) to run the GPU kernels on CPU threads, which is how the test suite checks them on machines without a GPU. `BIOTECH_NET_GPU=off` turns the GPU path off without touching code. Dispose the device when you are done; it owns the compiled kernels.

## 2. All-vs-all alignment scores

`GpuLocalAligner` scores every query against every target with Smith–Waterman and affine gaps, one GPU thread per pair:

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

Gap costs are positive here, BLAST style (`PairwiseAligner` uses negative gap *scores*; `GpuLocalAligner.FromPairwiseAligner` converts). The scores are not approximately equal to the CPU ones — they are the same integers, because the kernel runs the same recurrence in integer arithmetic. `scores.Values` gives the whole matrix row-major and `TopAll(k)` the top-k list of every query. The scores are the starting point; align the hits you keep with `PairwiseAligner` to get the alignments themselves.

## 3. A GPU backend for database search

`SequenceDatabase.Search` scores seeded subjects with Smith–Waterman. `SearchOptions.ScoreBackend` takes any `ILocalScoreBackend`, and `GpuLocalAligner` is one:

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

The hit list is identical to the CPU search, E-values included. A single query against a small database is not worth a GPU launch; the backend pays off with thousands of subjects or with `UseSeeds = false`.

## 4. Scoring thousands of docking poses

`GpuVinaScorer` evaluates the Vina search objective — grid-map receptor energy plus the ligand's internal energy — and its gradient for a whole batch of poses. A pose is the root position, an orientation quaternion and one angle per torsion:

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

The receptor and ligand files are the 1STP fixtures in `tests/Biotech.Net.DrugDiscovery.Tests/Data`. The GPU computes in single precision, so it agrees with the double-precision CPU engine (`EvaluateCpu`) to about 1e-4 kcal/mol; random poses mostly clash, which is why the crystal pose scores far better than the best of 4,096 random ones. `Optimize` runs the docking engine's BFGS on every pose at once:

```csharp
double[] start = scorer.RandomPoses(256, seed: 2);
double[] before = scorer.Evaluate(start);
double[] after = scorer.Optimize(start, maxSteps: 20);       // in place, one BFGS per thread
Console.WriteLine($"mean {before.Average():F1} -> {after.Average():F1} kcal/mol");
```

```text
mean 111.2 -> 10.9 kcal/mol
```

## 5. Docking on the GPU

`GpuVinaDocker` replaces the Monte-Carlo search of `VinaDocking.Dock` with a few thousand short chains, one per GPU thread, and then hands the best pose of every chain to the CPU post-processing (clustering, exact refinement, Vina rescoring). It takes the same `DockingOptions` and returns the same `DockingResult`:

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

The CPU search finds −6.19 kcal/mol for this complex too (AutoDock Vina: −6.20). By default the chains share the CPU search budget (Vina's step heuristic × exhaustiveness); `StepsPerChain` overrides it. The search runs in slices of a few steps per kernel launch so that no launch trips the Windows display-driver watchdog. This first run in the process includes a few seconds of kernel compilation; a second docking of the same complex takes about 1.5 s here, the same as the CPU search.

## 6. How fast is it?

On the development laptop (i7-8650U, Intel UHD 620), from [the benchmarks page](../benchmarks.md):

| Workload | Best CPU path | GPU (UHD 620) |
|---|---:|---:|
| 64 × 128 protein scores (≈ 330 M cells) | 34 ms (Rust, 8 threads) | 113 ms |
| 4,096 pose energies + gradients | 8.3 ms (managed, all cores) | 5.5 ms |
| Docking 3PTB, default budget | 1.55 s | 1.55 s |

An integrated GPU does not beat a vectorised four-core CPU at integer dynamic programming. Use this package where a discrete GPU is available, or to run very large batches without tying up the CPU.

## Limitations

- **Vulkan is not implemented**; ILGPU provides CUDA and OpenCL. There is no Metal backend for Apple GPUs.
- **Measured on an integrated GPU only.** The kernels have not been benchmarked on a discrete NVIDIA or AMD card; the speed-ups usually quoted for GPU Smith–Waterman and Vina-GPU are not demonstrated here.
- **Scores, not alignments.** The GPU computes local scores only (no traceback, no global or banded mode).
- **Float docking.** Energies match the CPU to about 1e-4 kcal/mol for a given pose, but search trajectories differ between devices, so docking with the same seed gives the same result only on the same device. Each chain keeps one best pose.
- **ILGPU quirks.** ILGPU's CPU accelerator needs the kernel argument types to be public (they are marked as infrastructure), and ILGPU 1.5's OpenCL backend miscompiles some int-to-float conversions, which the kernels avoid with float literals.

---
Previous: [Workflows](workflows.md) · [Tutorial index](index.md)
