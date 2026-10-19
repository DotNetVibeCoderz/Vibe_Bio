using System.Diagnostics;
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;
using ILGPU;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu.Docking;

/// <summary>Search settings for <see cref="GpuVinaDocker"/>; scoring, box, pose count and clustering come from <see cref="DockingOptions"/>.</summary>
public sealed record GpuDockingOptions
{
    /// <summary>Independent Monte-Carlo chains, one GPU thread each. Default 1024.</summary>
    public int Chains { get; init; } = 1024;

    /// <summary>
    /// Monte-Carlo steps per chain; 0 (default) spreads the CPU search budget (Vina's step heuristic × exhaustiveness)
    /// over the chains, with at least 64 steps per chain.
    /// </summary>
    public int StepsPerChain { get; init; }

    /// <summary>
    /// Target wall time of one kernel launch. The search runs in slices of a few steps so a launch stays well under the
    /// display-driver watchdog (2 s on Windows); the slice length adapts to the measured speed. Default 250 ms.
    /// </summary>
    public TimeSpan LaunchBudget { get; init; } = TimeSpan.FromMilliseconds(250);

    /// <summary>Default settings.</summary>
    public static GpuDockingOptions Default { get; } = new();
}

/// <summary>
/// Vina-style flexible-ligand docking with the Monte-Carlo search on a GPU (the Vina-GPU / QuickVina-GPU idea): thousands of
/// short independent chains, one per thread, each running mutation → BFGS (hunting cap) → Metropolis → refinement entirely
/// on the device in single precision against the same grid maps as <see cref="VinaDocking"/>. The best pose of every
/// chain then goes through the CPU post-processing of <see cref="VinaDocking.Dock"/> (RMSD clustering, double-precision
/// BFGS refinement with exact receptor interactions, Vina rescoring), so reported affinities are CPU-exact for those poses.
/// </summary>
/// <remarks>
/// With a null device, or if the device fails, the CPU <see cref="VinaDocking.Dock"/> runs instead
/// (<see cref="LastFallbackError"/> records why). Results are reproducible for a seed on a given device; different devices
/// (or the CPU path) follow different trajectories because of float rounding.
/// </remarks>
public sealed class GpuVinaDocker
{
    /// <summary>Docker on a device (null = CPU docking).</summary>
    public GpuVinaDocker(GpuDevice? device) => Device = device;

    /// <summary>The device, or null.</summary>
    public GpuDevice? Device { get; }

    /// <summary>Backend of the most recent <see cref="Dock"/> call.</summary>
    public string LastBackend { get; private set; } = "";

    /// <summary>The exception that made the most recent call fall back to the CPU, or null.</summary>
    public Exception? LastFallbackError { get; private set; }

    /// <summary>Kernel launches of the most recent GPU search.</summary>
    public int LastLaunches { get; private set; }

    /// <summary>Monte-Carlo steps per chain of the most recent GPU search.</summary>
    public int LastStepsPerChain { get; private set; }

    /// <summary>Test hook: makes the GPU path throw so the fallback can be exercised.</summary>
    internal bool InjectFailure { get; set; }

    /// <summary>Docks a ligand (3-D coordinates; embedded when missing) into a rigid receptor.</summary>
    public DockingResult Dock(Molecule ligand, VinaReceptor receptor, DockingOptions? options = null, GpuDockingOptions? gpu = null)
    {
        ArgumentNullException.ThrowIfNull(ligand);
        ArgumentNullException.ThrowIfNull(receptor);
        options ??= DockingOptions.Default;
        gpu ??= GpuDockingOptions.Default;
        ArgumentOutOfRangeException.ThrowIfLessThan(gpu.Chains, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(gpu.StepsPerChain);
        LastFallbackError = null;
        if (Device is { } device)
        {
            try
            {
                if (InjectFailure) throw new InvalidOperationException("Injected kernel failure (test hook).");
                var r = DockOnDevice(device, ligand, receptor, options, gpu);
                LastBackend = $"{device.Backend}: {device.Name}";
                return r;
            }
            catch (Exception ex) when (ex is not ObjectDisposedException and not ArgumentException)
            {
                LastFallbackError = ex;
            }
        }
        LastBackend = "CPU (managed)";
        return VinaDocking.Dock(ligand, receptor, options);
    }

    private DockingResult DockOnDevice(GpuDevice device, Molecule ligand, VinaReceptor receptor, DockingOptions options, GpuDockingOptions gpu)
    {
        if (ligand.Coordinates is null || !ligand.Has3DCoordinates) ligand = ConformerGenerator.Embed(ligand);
        var lig = new FlexibleLigand(ligand);
        if (lig.Scored.Length == 0) throw new ArgumentException("The ligand has no atoms Vina can score.", nameof(ligand));
        var center = options.Center ?? lig.ReferenceOrigin;
        var sw = Stopwatch.StartNew();
        var grid = VinaGrid.Build(receptor, center, options.Size, lig.Scored.Select(i => lig.Types[i]!.Value), options.GridSpacing, options.Weights);
        var gridTime = sw.Elapsed;
        sw.Restart();

        var (cpuSteps, localSteps) = VinaDocking.SearchSteps(ligand, lig, options);
        int chains = gpu.Chains;
        int steps = gpu.StepsPerChain > 0
            ? gpu.StepsPerChain
            : (int)Math.Max(64, Math.Ceiling((double)cpuSteps * Math.Max(1, options.Exhaustiveness) / chains));
        LastStepsPerChain = steps;

        using var model = new VinaGpuModel(device, lig, grid, options.Weights);
        int nc = model.ConfLength;
        var bestE = new float[chains];
        var bestPose = new float[chains * nc];
        var evals = new int[chains];
        int launches = 0;
        lock (device.Sync)
        {
            device.ThrowIfDisposed();
            var acc = device.Accelerator;
            var kernel = device.Kernel("vina-mc", a => a.LoadAutoGroupedStreamKernel<Index1D, VinaKernelData, VinaScratch, VinaMcArgs>(VinaKernels.MonteCarlo));
            using var dS = acc.Allocate1D<float>((long)chains * model.ScratchPerLane);
            using var dEv = acc.Allocate1D<int>(chains);
            using var dRng = acc.Allocate1D<ulong>(chains);
            using var dBestE = acc.Allocate1D<float>(chains);
            using var dBest = acc.Allocate1D<float>((long)chains * nc);
            var scratch = new VinaScratch { Data = dS.View, Evaluations = dEv.View, Lanes = chains };
            var args = new VinaMcArgs
            {
                Rng = dRng.View, BestEnergy = dBestE.View, BestPose = dBest.View,
                Seed = Scramble((ulong)(uint)options.Seed), LocalSteps = localSteps,
            };
            int slice = 1;
            var clock = new Stopwatch();
            for (int first = 0; first < steps;)
            {
                int n = Math.Min(slice, steps - first);
                args.FirstStep = first;
                args.Steps = n;
                clock.Restart();
                kernel(chains, model.Data, scratch, args);
                acc.Synchronize();
                launches++;
                first += n;
                // grow (or shrink) the slice towards the launch budget; at most ×4 per launch
                double perStep = clock.Elapsed.TotalMilliseconds / n;
                int target = perStep > 0 ? (int)(gpu.LaunchBudget.TotalMilliseconds / perStep) : slice * 4;
                slice = Math.Clamp(target, 1, Math.Max(1, slice * 4));
            }
            dBestE.CopyToCPU(bestE);
            dBest.CopyToCPU(bestPose);
            dEv.CopyToCPU(evals);
        }
        LastLaunches = launches;

        // best pose of every chain -> CPU post-processing (clustering, exact refinement, rescoring)
        var candidates = new List<(LigandConf Conf, double E, Vec3[] Heavy)>(chains);
        var x = new Vec3[lig.AtomCount];
        for (int c = 0; c < chains; c++)
        {
            if (!float.IsFinite(bestE[c])) continue;
            var conf = model.FromDevice(bestPose.AsSpan(c * nc, nc));
            lig.Build(conf, x);
            candidates.Add((conf, bestE[c], lig.Scored.Select(i => x[i]).ToArray()));
        }
        // keep the best distinct poses only, so the exact refinement stays cheap with thousands of chains
        candidates.Sort((a, b) => a.E.CompareTo(b.E));
        var distinct = new List<(LigandConf Conf, double E, Vec3[] Heavy)>();
        int keep = Math.Max(3 * options.NumPoses, 20);
        foreach (var cand in candidates)
        {
            if (distinct.Count >= keep) break;
            if (distinct.All(m => VinaDocking.Rmsd(m.Heavy, cand.Heavy) >= options.MinRmsd)) distinct.Add(cand);
        }
        if (distinct.Count == 0) throw new InvalidOperationException("The GPU search produced no finite pose.");
        long total = evals.Sum(e => (long)e);
        return VinaDocking.Finish(ligand, lig, grid, receptor, options, distinct, localSteps, (int)Math.Min(int.MaxValue, total), gridTime, sw);
    }

    /// <summary>SplitMix64 finaliser, so neighbouring seeds give unrelated chain streams.</summary>
    private static ulong Scramble(ulong z)
    {
        z += 0x9E3779B97F4A7C15UL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        return z ^ (z >> 31);
    }
}
