using BenchmarkDotNet.Attributes;
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.Gpu;
using Biotech.Net.Gpu.Alignment;
using Biotech.Net.Gpu.Docking;
using Biotech.Net.Native;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.Benchmarks;

/// <summary>Backends present on this machine: the CPU paths always, ILGPU's CPU accelerator always, OpenCL/CUDA when found.</summary>
internal static class GpuBackends
{
    public static IEnumerable<string> Alignment()
    {
        yield return "Managed";
        if (BioKernels.IsNativeAvailable) yield return "Rust";
        yield return "ILGPU-CPU";
        foreach (var b in Gpus()) yield return b;
    }

    public static IEnumerable<string> Docking()
    {
        yield return "Managed-1T";
        yield return "Managed-Parallel";
        yield return "ILGPU-CPU";
        foreach (var b in Gpus()) yield return b;
    }

    public static IEnumerable<string> Gpus()
    {
        var devices = GpuDevice.Enumerate();
        if (devices.Any(d => d.Backend == GpuBackend.Cuda)) yield return "ILGPU-CUDA";
        if (devices.Any(d => d.Backend == GpuBackend.OpenCL)) yield return "ILGPU-OpenCL";
    }

    public static GpuDevice? Open(string backend) => backend switch
    {
        "ILGPU-CPU" => GpuDevice.CreateCpu(),
        "ILGPU-CUDA" => GpuDevice.TryCreate(GpuBackend.Cuda),
        "ILGPU-OpenCL" => GpuDevice.TryCreate(GpuBackend.OpenCL),
        _ => null,
    };

    public static readonly string DataDir = Path.Combine(AppContext.BaseDirectory, "Data");
}

/// <summary>
/// All-vs-all Smith–Waterman scores (BLOSUM62, gaps 11/1): 64 queries × 128 targets of 100–300 aa (≈ 330 M cells).
/// The CPU backends run one query per core (Parallel.For over queries, as <see cref="GpuLocalAligner"/> does without a device).
/// </summary>
public class AllVsAllBenchmarks
{
    public const int Queries = 64, Targets = 128;

    private byte[] _q = [], _t = [];
    private int[] _qo = [], _to = [], _scores = [], _matrix = [];
    private GpuDevice? _device;
    private GpuLocalAligner? _aligner;

    [ParamsSource(nameof(Backends))]
    public string Backend { get; set; } = "Managed";

    public static IEnumerable<string> Backends => GpuBackends.Alignment();

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(42);
        const string aa = "ACDEFGHIKLMNPQRSTVWY";
        (byte[], int[]) Make(int n)
        {
            var off = new int[n + 1];
            var buf = new List<byte>();
            for (int i = 0; i < n; i++)
            {
                int len = rng.Next(100, 301);
                for (int j = 0; j < len; j++) buf.Add((byte)aa[rng.Next(aa.Length)]);
                off[i + 1] = buf.Count;
            }
            return ([.. buf], off);
        }
        (_q, _qo) = Make(Queries);
        (_t, _to) = Make(Targets);
        _scores = new int[Queries * Targets];
        _matrix = Biotech.Net.Alignments.SubstitutionMatrix.Blosum62.KernelMatrix.ToArray();
        _device = GpuBackends.Open(Backend);
        if (_device is not null) _aligner = new GpuLocalAligner(_device);
    }

    [GlobalCleanup]
    public void Cleanup() => _device?.Dispose();

    [BenchmarkCategory("All-vs-all SW 64 x 128 (100-300 aa)"), Benchmark]
    public void AllVsAll()
    {
        switch (Backend)
        {
            case "Managed":
                Parallel.For(0, Queries, i => ManagedKernels.LocalScoreBatch(_q.AsSpan(_qo[i], _qo[i + 1] - _qo[i]), _t, _to, _matrix, 11, 1, _scores.AsSpan(i * Targets, Targets)));
                break;
            case "Rust":
                Parallel.For(0, Queries, i => NativeBridge.LocalScoreBatch(_q.AsSpan(_qo[i], _qo[i + 1] - _qo[i]), _t, _to, _matrix, 11, 1, _scores.AsSpan(i * Targets, Targets)));
                break;
            default:
                _aligner!.Score(_q, _qo, _t, _to, _matrix, 11, 1, _scores);
                if (_aligner.LastFallbackError is { } e) throw new InvalidOperationException("GPU path fell back", e);
                break;
        }
    }
}

/// <summary>Vina objective + gradient for 4096 random poses of biotin in streptavidin (1STP, 5 torsions, 16 heavy atoms).</summary>
public class DockingPoseBenchmarks
{
    public const int Poses = 4096;

    private GpuVinaScorer? _scorer;
    private GpuDevice? _device;
    private double[] _poses = [], _gradients = [];

    [ParamsSource(nameof(Backends))]
    public string Backend { get; set; } = "Managed-1T";

    public static IEnumerable<string> Backends => GpuBackends.Docking();

    [GlobalSetup]
    public void Setup()
    {
        var ligand = SdfReader.ReadFile(Path.Combine(GpuBackends.DataDir, "1stp_ligand.sdf")).First();
        var receptor = VinaReceptor.FromPdb(Path.Combine(GpuBackends.DataDir, "1stp_receptor.pdb"));
        var center = ligand.Coordinates!.Aggregate(Vec3.Zero, (s, p) => s + new Vec3(p.X, p.Y, p.Z)) / ligand.AtomCount;
        _device = GpuBackends.Open(Backend);
        _scorer = GpuVinaScorer.Create(ligand, receptor, center, new Vec3(22.5, 22.5, 22.5), _device);
        _poses = _scorer.RandomPoses(Poses, seed: 1, spread: 4);
        _gradients = new double[Poses * _scorer.GradientLength];
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        _scorer?.Dispose();
        _device?.Dispose();
    }

    [BenchmarkCategory("Vina pose energy+gradient x4096 (1STP)"), Benchmark]
    public double Evaluate()
    {
        var s = _scorer!;
        switch (Backend)
        {
            case "Managed-1T":
                return s.EvaluateCpu(_poses, 1000, _gradients)[0];
            case "Managed-Parallel":
                {
                    int chunk = Poses / Environment.ProcessorCount + 1, pl = s.PoseLength, gl = s.GradientLength;
                    var results = new double[Poses];
                    Parallel.For(0, (Poses + chunk - 1) / chunk, c =>
                    {
                        int a = c * chunk, n = Math.Min(chunk, Poses - a);
                        s.EvaluateCpu(_poses.AsSpan(a * pl, n * pl), 1000, _gradients.AsSpan(a * gl, n * gl)).CopyTo(results, a);
                    });
                    return results[0];
                }
            default:
                var e = s.Evaluate(_poses, 1000, _gradients);
                if (s.LastFallbackError is { } ex) throw new InvalidOperationException("GPU path fell back", ex);
                return e[0];
        }
    }
}

/// <summary>
/// Whole docking runs of benzamidine into trypsin (3PTB): CPU Vina search (exhaustiveness 8, all cores) against the GPU
/// search (1024 chains, same total step budget). ILGPU-CPU is left out: it emulates thousands of threads and takes minutes.
/// </summary>
public class DockingSearchBenchmarks
{
    private Molecule _ligand = null!;
    private VinaReceptor _receptor = null!;
    private DockingOptions _options = null!;
    private GpuDevice? _device;

    [ParamsSource(nameof(Backends))]
    public string Backend { get; set; } = "CPU-Vina";

    public static IEnumerable<string> Backends => GpuBackends.Gpus().Prepend("CPU-Vina");

    [GlobalSetup]
    public void Setup()
    {
        _ligand = SdfReader.ReadFile(Path.Combine(GpuBackends.DataDir, "3ptb_ligand.sdf")).First();
        _receptor = VinaReceptor.FromPdb(Path.Combine(GpuBackends.DataDir, "3ptb_receptor.pdb"));
        _options = new DockingOptions { Seed = 1 };
        _device = GpuBackends.Open(Backend);
    }

    [GlobalCleanup]
    public void Cleanup() => _device?.Dispose();

    [BenchmarkCategory("Vina docking 3PTB (default budget)"), Benchmark]
    public double Dock()
    {
        if (_device is null) return VinaDocking.Dock(_ligand, _receptor, _options).Best.Affinity;
        var docker = new GpuVinaDocker(_device);
        var r = docker.Dock(_ligand, _receptor, _options);
        if (docker.LastFallbackError is { } ex) throw new InvalidOperationException("GPU path fell back", ex);
        return r.Best.Affinity;
    }
}
