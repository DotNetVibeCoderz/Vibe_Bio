using System.Text.Json;
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.Gpu.Docking;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.Gpu.Tests;

/// <summary>GPU Vina scoring and search on trypsin–benzamidine (3PTB) and streptavidin–biotin (1STP).</summary>
public class DockingTests(ITestOutputHelper output)
{
    private static readonly string DataDir = Path.Combine(AppContext.BaseDirectory, "Data");
    private static readonly Lazy<JsonDocument> Reference = new(() => JsonDocument.Parse(File.ReadAllText(Path.Combine(DataDir, "vina_reference.json"))));
    private static readonly Dictionary<string, VinaReceptor> Receptors = [];
    private static readonly Dictionary<string, VinaGrid> Grids = [];

    internal static VinaReceptor Receptor(string id)
    {
        lock (Receptors)
        {
            if (!Receptors.TryGetValue(id, out var r)) Receptors[id] = r = VinaReceptor.FromPdb(Path.Combine(DataDir, id + "_receptor.pdb"));
            return r;
        }
    }

    internal static Molecule Crystal(string id) => SdfReader.ReadFile(Path.Combine(DataDir, id + "_ligand.sdf")).First();

    private static JsonElement Ref(string id) => Reference.Value.RootElement.GetProperty("complexes").GetProperty(id);

    internal static Vec3 Centre(string id)
    {
        var c = Ref(id).GetProperty("center");
        return new Vec3(c[0].GetDouble(), c[1].GetDouble(), c[2].GetDouble());
    }

    internal static VinaGrid Grid(string id)
    {
        lock (Grids)
        {
            if (!Grids.TryGetValue(id, out var g))
            {
                var types = VinaAtomTypes.AssignLigand(Crystal(id)).Where(t => t is not null).Select(t => t!.Value);
                Grids[id] = g = VinaGrid.Build(Receptor(id), Centre(id), new Vec3(22.5, 22.5, 22.5), types);
            }
            return g;
        }
    }

    /// <summary>Heavy-atom RMSD to the crystal ligand, minimised over graph symmetries (as in the DrugDiscovery tests).</summary>
    private static double RmsdToCrystal(Molecule crystal, Molecule pose)
    {
        var skeleton = crystal.Clone();
        for (int b = 0; b < skeleton.BondCount; b++) skeleton.BondAt(b).Order = BondOrder.Single;
        for (int a = 0; a < skeleton.AtomCount; a++) { skeleton.AtomAt(a).IsAromatic = false; skeleton.AtomAt(a).FormalCharge = 0; }
        skeleton.InvalidateTopology();
        var heavy = Enumerable.Range(0, pose.AtomCount).Where(i => pose.Atoms[i].AtomicNumber != 1).ToArray();
        return ConformerAlignment.BestRmsd(skeleton, crystal.Coordinates!, heavy.Select(i => pose.Coordinates![i]).ToArray(), align: false);
    }

    [Theory]
    [InlineData("3ptb", 1000.0)]
    [InlineData("1stp", 1000.0)]
    [InlineData("1stp", 10.0)]
    public void Pose_energies_and_gradients_match_the_cpu_engine(string id, double cap)
    {
        foreach (var device in GpuTestDevices.All)
        {
            using var scorer = new GpuVinaScorer(Crystal(id), Grid(id), device);
            var poses = scorer.RandomPoses(256, seed: 7, spread: 2.5).Concat(scorer.InputPose()).ToArray();
            int count = poses.Length / scorer.PoseLength;
            var g = new double[count * scorer.GradientLength];
            var gRef = new double[count * scorer.GradientLength];
            var e = scorer.Evaluate(poses, cap, g);
            Assert.Null(scorer.LastFallbackError);
            Assert.DoesNotContain("CPU (", scorer.LastBackend);
            var eRef = scorer.EvaluateCpu(poses, cap, gRef);
            double maxAbs = 0, maxRel = 0, maxGrad = 0;
            for (int k = 0; k < count; k++)
            {
                double diff = Math.Abs(e[k] - eRef[k]);
                maxAbs = Math.Max(maxAbs, Math.Abs(eRef[k]) < 100 ? diff : 0);
                maxRel = Math.Max(maxRel, diff / Math.Max(1, Math.Abs(eRef[k])));
                // single precision: 1e-4 kcal/mol absolute (plus 1e-5 relative for clashing poses with large energies)
                Assert.True(diff <= 1e-4 + 1e-5 * Math.Abs(eRef[k]), $"{device.Backend} pose {k}: {e[k]} vs {eRef[k]}");
                for (int j = 0; j < scorer.GradientLength; j++)
                {
                    double a = g[k * scorer.GradientLength + j], b = gRef[k * scorer.GradientLength + j];
                    double gd = Math.Abs(a - b) / Math.Max(1, Math.Abs(b));
                    maxGrad = Math.Max(maxGrad, gd);
                    Assert.True(gd < 1e-3, $"{device.Backend} pose {k} dof {j}: {a} vs {b}");
                }
            }
            output.WriteLine($"{id} cap={cap} {device}: max |ΔE| {maxAbs:E2} kcal/mol (|E|<100), max rel {maxRel:E2}, max grad rel {maxGrad:E2}");
        }
    }

    [Fact]
    public void Batched_bfgs_lowers_energies_like_the_cpu_engine()
    {
        foreach (var device in GpuTestDevices.All)
        {
            using var scorer = new GpuVinaScorer(Crystal("1stp"), Grid("1stp"), device);
            var start = scorer.RandomPoses(64, seed: 3, spread: 3.0);
            var e0 = scorer.EvaluateCpu(start);
            var gpu = (double[])start.Clone();
            var cpu = (double[])start.Clone();
            var eg = scorer.Optimize(gpu, maxSteps: 14);
            Assert.Null(scorer.LastFallbackError);
            Assert.DoesNotContain("CPU (", scorer.LastBackend);
            var ec = scorer.OptimizeCpu(cpu, maxSteps: 14);
            // float and double trajectories drift apart on these clashing random starts, but every pose must improve,
            // a fair share must land where the CPU lands, and the averages must be comparable
            for (int k = 0; k < e0.Length; k++) Assert.True(eg[k] <= e0[k] + 1e-3, $"{device.Backend} pose {k}: {e0[k]} -> {eg[k]}");
            int close = eg.Zip(ec).Count(p => Math.Abs(p.First - p.Second) < 0.05);
            output.WriteLine($"{device}: mean start {e0.Average():F2}, GPU {eg.Average():F2}, CPU {ec.Average():F2}, {close}/64 within 0.05");
            Assert.True(close >= 16, $"{device.Backend}: only {close}/64 optimised poses agree with the CPU engine");
            Assert.True(Math.Abs(eg.Average() - ec.Average()) < 0.1 * e0.Average(), $"{device.Backend}: mean {eg.Average():F2} vs CPU {ec.Average():F2}");
            // the returned poses carry the reported energies
            var re = scorer.EvaluateCpu(gpu);
            for (int k = 0; k < re.Length; k++) Assert.True(Math.Abs(re[k] - eg[k]) < 1e-3 + 1e-5 * Math.Abs(re[k]));
        }
    }

    [Fact]
    public void Crystal_pose_scores_like_the_cpu_and_null_device_uses_the_cpu()
    {
        using var gpu = new GpuVinaScorer(Crystal("1stp"), Grid("1stp"), GpuTestDevices.Cpu);
        using var cpu = new GpuVinaScorer(Crystal("1stp"), Grid("1stp"), null);
        var pose = gpu.InputPose();
        double eg = gpu.Evaluate(pose)[0], ec = cpu.Evaluate(pose)[0];
        Assert.Equal("CPU (managed)", cpu.LastBackend);
        Assert.Equal(cpu.EvaluateCpu(pose)[0], ec);
        Assert.Equal(ec, eg, 1e-4);
        Assert.True(ec < -8, $"crystal biotin should bind strongly: {ec}");
        // the input pose reproduces the input coordinates
        var xyz = gpu.Coordinates(pose);
        var crystal = Crystal("1stp").Coordinates!;
        for (int i = 0; i < xyz.Length; i++) Assert.True(Math.Abs(xyz[i].X - crystal[i].X) + Math.Abs(xyz[i].Y - crystal[i].Y) + Math.Abs(xyz[i].Z - crystal[i].Z) < 1e-9);
        Assert.Throws<ArgumentException>(() => gpu.Evaluate(new double[gpu.PoseLength + 1]));
        Assert.Empty(gpu.Evaluate([]));
    }

    [Theory]
    [InlineData("3ptb", 1.0, 0.3)]
    [InlineData("1stp", 1.5, 0.5)]
    public void Gpu_search_redocks_the_crystal_ligand(string id, double maxRmsd, double maxEnergyGap)
    {
        var docker = new GpuVinaDocker(GpuTestDevices.Cpu);
        var options = new DockingOptions { Center = Centre(id), Seed = 3 };
        var res = docker.Dock(Crystal(id), Receptor(id), options, new GpuDockingOptions { Chains = 256, StepsPerChain = 120 });
        Assert.Null(docker.LastFallbackError);
        Assert.StartsWith("Cpu:", docker.LastBackend);
        var vina = Ref(id).GetProperty("redock")[0];
        double rmsd = RmsdToCrystal(Crystal(id), res.Best.Molecule);
        output.WriteLine($"{id}: best {res.Best.Affinity:F3} kcal/mol (Vina {vina.GetProperty("affinity").GetDouble():F3}), RMSD {rmsd:F2} Å, " +
            $"{res.Poses.Count} poses, {res.Evaluations} evaluations, {docker.LastLaunches} launches, search {res.SearchTime.TotalSeconds:F1} s");
        Assert.True(rmsd < maxRmsd, $"{id}: top pose {rmsd:F2} Å from the crystal ligand");
        Assert.True(Math.Abs(res.Best.Affinity - vina.GetProperty("affinity").GetDouble()) < maxEnergyGap,
            $"{id}: best affinity {res.Best.Affinity:F3} vs Vina {vina.GetProperty("affinity").GetDouble():F3}");
        Assert.True(res.Poses.Zip(res.Poses.Skip(1)).All(p => p.First.Affinity <= p.Second.Affinity + 1e-9));
        Assert.True(docker.LastLaunches >= 1);
    }

    [Fact]
    public void Gpu_search_is_deterministic_on_the_cpu_accelerator()
    {
        var options = new DockingOptions { Center = Centre("3ptb"), Seed = 11, NumPoses = 3 };
        var gpuOptions = new GpuDockingOptions { Chains = 32, StepsPerChain = 20 };
        var a = new GpuVinaDocker(GpuTestDevices.Cpu).Dock(Crystal("3ptb"), Receptor("3ptb"), options, gpuOptions);
        var b = new GpuVinaDocker(GpuTestDevices.Cpu).Dock(Crystal("3ptb"), Receptor("3ptb"), options, gpuOptions with { LaunchBudget = TimeSpan.FromTicks(1) });
        Assert.Equal(a.Poses.Count, b.Poses.Count);
        for (int k = 0; k < a.Poses.Count; k++) Assert.Equal(a.Poses[k].Molecule.Coordinates, b.Poses[k].Molecule.Coordinates);
    }

    [Fact]
    public void Null_device_and_kernel_failure_fall_back_to_cpu_docking()
    {
        var options = new DockingOptions { Center = Centre("3ptb"), MonteCarloSteps = 100, Exhaustiveness = 2, Seed = 9 };
        var cpu = VinaDocking.Dock(Crystal("3ptb"), Receptor("3ptb"), options);
        var none = new GpuVinaDocker(null);
        var r1 = none.Dock(Crystal("3ptb"), Receptor("3ptb"), options);
        Assert.Equal("CPU (managed)", none.LastBackend);
        var failing = new GpuVinaDocker(GpuTestDevices.Cpu) { InjectFailure = true };
        var r2 = failing.Dock(Crystal("3ptb"), Receptor("3ptb"), options);
        Assert.IsType<InvalidOperationException>(failing.LastFallbackError);
        Assert.Equal(cpu.Best.Molecule.Coordinates, r1.Best.Molecule.Coordinates);
        Assert.Equal(cpu.Best.Molecule.Coordinates, r2.Best.Molecule.Coordinates);
    }

    [Fact]
    public void Real_gpu_docks_when_present()
    {
        if (GpuTestDevices.Gpu is not { } gpu) return;   // GPU-gated
        var docker = new GpuVinaDocker(gpu);
        var res = docker.Dock(Crystal("3ptb"), Receptor("3ptb"), new DockingOptions { Center = Centre("3ptb") }, new GpuDockingOptions { Chains = 512, StepsPerChain = 100 });
        Assert.Null(docker.LastFallbackError);
        double rmsd = RmsdToCrystal(Crystal("3ptb"), res.Best.Molecule);
        output.WriteLine($"{gpu}: best {res.Best.Affinity:F3} kcal/mol, RMSD {rmsd:F2} Å, {docker.LastLaunches} launches, search {res.SearchTime.TotalSeconds:F2} s");
        Assert.True(rmsd < 1.0);
        Assert.True(Math.Abs(res.Best.Affinity - Ref("3ptb").GetProperty("redock")[0].GetProperty("affinity").GetDouble()) < 0.3);
    }
}
