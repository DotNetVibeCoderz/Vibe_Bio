using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;
using ILGPU;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu.Docking;

/// <summary>
/// Batched evaluation of the Vina search objective — grid-map intermolecular energy (trilinear interpolation, out-of-box
/// penalty) plus the intramolecular pair energy, both curl-capped — and its conformation-space gradient, for many ligand
/// poses at once. A pose is <see cref="PoseLength"/> doubles: root position (Å, absolute), orientation quaternion
/// (w, x, y, z) and one angle (radians) per torsion relative to the input geometry.
/// </summary>
/// <remarks>
/// The device computes in single precision in a frame centred on the box; the CPU reference (<see cref="EvaluateCpu"/>,
/// the docking engine's own double-precision code) agrees to about 1e-4 kcal/mol for poses inside the box. Energies here
/// are the search objective, not the reported affinity (which divides by 1 + w·N<sub>rot</sub> and uses exact receptor
/// interactions). With a null device <see cref="Evaluate"/> runs the CPU code.
/// </remarks>
public sealed class GpuVinaScorer : IDisposable
{
    private readonly VinaGpuModel? _model;
    private readonly VinaScoreOptions _weights;

    /// <summary>Scorer for a ligand (3-D coordinates required) and prebuilt grid maps.</summary>
    public GpuVinaScorer(Molecule ligand, VinaGrid grid, GpuDevice? device, VinaScoreOptions? weights = null)
    {
        ArgumentNullException.ThrowIfNull(ligand);
        ArgumentNullException.ThrowIfNull(grid);
        if (ligand.Coordinates is null || !ligand.Has3DCoordinates) throw new ArgumentException("The ligand needs 3-D coordinates.", nameof(ligand));
        _weights = weights ?? VinaScoreOptions.Default;
        Ligand = new FlexibleLigand(ligand);
        if (Ligand.Scored.Length == 0) throw new ArgumentException("The ligand has no atoms Vina can score.", nameof(ligand));
        foreach (int i in Ligand.Scored)
            if (!grid.HasType(Ligand.Types[i]!.Value)) throw new ArgumentException($"The grid has no map for {Ligand.Types[i]}.", nameof(grid));
        Grid = grid;
        Device = device;
        if (device is not null) _model = new VinaGpuModel(device, Ligand, grid, _weights);
    }

    /// <summary>Builds the grid maps for the ligand's atom types and creates a scorer.</summary>
    public static GpuVinaScorer Create(Molecule ligand, VinaReceptor receptor, Vec3 center, Vec3 size, GpuDevice? device,
        double spacing = 0.375, VinaScoreOptions? weights = null)
    {
        ArgumentNullException.ThrowIfNull(ligand);
        var types = VinaAtomTypes.AssignLigand(ligand).Where(t => t is not null).Select(t => t!.Value);
        var grid = VinaGrid.Build(receptor, center, size, types, spacing, weights);
        return new GpuVinaScorer(ligand, grid, device, weights);
    }

    internal FlexibleLigand Ligand { get; }

    /// <summary>The grid maps.</summary>
    public VinaGrid Grid { get; }

    /// <summary>The device, or null for the CPU path.</summary>
    public GpuDevice? Device { get; }

    /// <summary>Torsional degrees of freedom.</summary>
    public int TorsionCount => Ligand.Torsions.Length;

    /// <summary>Values per pose: 3 position + 4 quaternion + torsions.</summary>
    public int PoseLength => 7 + TorsionCount;

    /// <summary>Values per gradient: 3 force + 3 torque + torsions (dE/dposition, rotation-vector derivative, dE/dangle).</summary>
    public int GradientLength => 6 + TorsionCount;

    /// <summary>The exception that made the most recent <see cref="Evaluate"/> fall back to the CPU, or null.</summary>
    public Exception? LastFallbackError { get; private set; }

    /// <summary>Backend of the most recent <see cref="Evaluate"/> call.</summary>
    public string LastBackend { get; private set; } = "";

    /// <summary>The pose that reproduces the input coordinates.</summary>
    public double[] InputPose() => ToPose(Ligand.InputConf());

    /// <summary>
    /// <paramref name="count"/> poses: random root positions within <paramref name="spread"/> Å of the input root (inside
    /// the box), uniform random orientations and torsions. Deterministic for a seed.
    /// </summary>
    public double[] RandomPoses(int count, int seed, double spread = 3.0)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        var rng = new Random(seed);
        var poses = new double[count * PoseLength];
        var lo = Grid.Origin;
        var hi = Grid.Origin + Grid.Size;
        var c = Ligand.InputConf().Position;
        for (int k = 0; k < count; k++)
        {
            var conf = new LigandConf(
                new Vec3(Math.Clamp(c.X + (rng.NextDouble() * 2 - 1) * spread, lo.X, hi.X),
                         Math.Clamp(c.Y + (rng.NextDouble() * 2 - 1) * spread, lo.Y, hi.Y),
                         Math.Clamp(c.Z + (rng.NextDouble() * 2 - 1) * spread, lo.Z, hi.Z)),
                Quat.Random(rng),
                Enumerable.Range(0, TorsionCount).Select(_ => (rng.NextDouble() * 2 - 1) * Math.PI).ToArray());
            ToPose(conf).CopyTo(poses.AsSpan(k * PoseLength));
        }
        return poses;
    }

    /// <summary>Heavy and hydrogen atom coordinates of a pose (input atom order).</summary>
    public Point3[] Coordinates(ReadOnlySpan<double> pose)
    {
        var x = new Vec3[Ligand.AtomCount];
        Ligand.Build(ToConf(pose), x);
        return x.Select(v => new Point3(v.X, v.Y, v.Z)).ToArray();
    }

    /// <summary>
    /// Energies of every pose in <paramref name="poses"/> (count × <see cref="PoseLength"/>); gradients
    /// (count × <see cref="GradientLength"/>) are written when <paramref name="gradients"/> is non-empty.
    /// <paramref name="cap"/> is Vina's curl cap (1000 = none; the search uses 10 while hunting).
    /// </summary>
    public double[] Evaluate(ReadOnlySpan<double> poses, double cap = 1000, Span<double> gradients = default)
    {
        int count = CheckPoses(poses, gradients);
        if (_model is null)
        {
            LastBackend = "CPU (managed)";
            return EvaluateCpu(poses, cap, gradients);
        }
        var energies = new double[count];
        if (count == 0) { LastBackend = $"{_model.Device.Backend}: {_model.Device.Name}"; return energies; }
        int nc = _model.ConfLength, ng = _model.GradLength;
        var conf = new float[count * nc];
        for (int k = 0; k < count; k++) _model.ToDevice(poses.Slice(k * nc, nc), conf.AsSpan(k * nc, nc));
        var e = new float[count];
        var g = new float[gradients.IsEmpty ? 1 : count * ng];
        var device = _model.Device;
        LastFallbackError = null;
        try
        {
            Launch(device, conf, count, cap, e, g, !gradients.IsEmpty);
        }
        catch (Exception ex) when (ex is not ObjectDisposedException)
        {
            LastFallbackError = ex;
            LastBackend = "CPU (managed)";
            return EvaluateCpu(poses, cap, gradients);
        }
        for (int k = 0; k < count; k++) energies[k] = e[k];
        if (!gradients.IsEmpty) for (int k = 0; k < count * ng; k++) gradients[k] = g[k];
        LastBackend = $"{device.Backend}: {device.Name}";
        return energies;
    }

    private void Launch(GpuDevice device, float[] conf, int count, double cap, float[] e, float[] g, bool wantGradients)
    {
        lock (device.Sync)
        {
            device.ThrowIfDisposed();
            var acc = device.Accelerator;
            var kernel = device.Kernel("vina-eval", a => a.LoadAutoGroupedStreamKernel<Index1D, VinaKernelData, VinaScratch, VinaEvalArgs>(VinaKernels.EvaluatePoses));
            using var dPoses = acc.Allocate1D(conf);
            using var dE = acc.Allocate1D<float>(count);
            using var dG = acc.Allocate1D<float>(g.Length);
            using var dS = acc.Allocate1D<float>((long)count * _model!.ScratchPerLane);
            using var dEv = acc.Allocate1D<int>(count);
            kernel(count, _model.Data, new VinaScratch { Data = dS.View, Evaluations = dEv.View, Lanes = count },
                new VinaEvalArgs { Poses = dPoses.View, Energies = dE.View, Gradients = dG.View, Cap = (float)cap, WantGradients = wantGradients ? 1 : 0 });
            acc.Synchronize();
            dE.CopyToCPU(e);
            if (wantGradients) dG.CopyToCPU(g);
        }
    }

    /// <summary>
    /// BFGS local optimisation (the docking engine's algorithm, at most <paramref name="maxSteps"/> iterations) of every
    /// pose, in place; returns the final energies. With a null device, or on a kernel failure, the CPU engine runs.
    /// </summary>
    public double[] Optimize(Span<double> poses, int maxSteps = 20, double cap = 1000)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxSteps);
        int count = CheckPoses(poses, default);
        if (_model is not null && count > 0)
        {
            int nc = _model.ConfLength;
            var conf = new float[count * nc];
            for (int k = 0; k < count; k++) _model.ToDevice(poses.Slice(k * nc, nc), conf.AsSpan(k * nc, nc));
            var e = new float[count];
            var device = _model.Device;
            LastFallbackError = null;
            try
            {
                lock (device.Sync)
                {
                    device.ThrowIfDisposed();
                    var acc = device.Accelerator;
                    var kernel = device.Kernel("vina-opt", a => a.LoadAutoGroupedStreamKernel<Index1D, VinaKernelData, VinaScratch, VinaEvalArgs, int>(VinaKernels.OptimizePoses));
                    using var dPoses = acc.Allocate1D(conf);
                    using var dE = acc.Allocate1D<float>(count);
                    using var dS = acc.Allocate1D<float>((long)count * _model.ScratchPerLane);
                    using var dEv = acc.Allocate1D<int>(count);
                    kernel(count, _model.Data, new VinaScratch { Data = dS.View, Evaluations = dEv.View, Lanes = count },
                        new VinaEvalArgs { Poses = dPoses.View, Energies = dE.View, Gradients = dE.View, Cap = (float)cap }, maxSteps);
                    acc.Synchronize();
                    dE.CopyToCPU(e);
                    dPoses.CopyToCPU(conf);
                }
                for (int k = 0; k < count; k++)
                {
                    ToPose(_model.FromDevice(conf.AsSpan(k * nc, nc))).CopyTo(poses.Slice(k * nc, nc));
                }
                LastBackend = $"{device.Backend}: {device.Name}";
                return e.Select(v => (double)v).ToArray();
            }
            catch (Exception ex) when (ex is not ObjectDisposedException)
            {
                LastFallbackError = ex;
            }
        }
        LastBackend = "CPU (managed)";
        return OptimizeCpu(poses, maxSteps, cap);
    }

    /// <summary>The CPU engine's BFGS on every pose, in place (double precision reference for <see cref="Optimize"/>).</summary>
    public double[] OptimizeCpu(Span<double> poses, int maxSteps = 20, double cap = 1000)
    {
        int count = CheckPoses(poses, default);
        var engine = new VinaDocking.Engine(Ligand, Grid, _weights);
        var energies = new double[count];
        for (int k = 0; k < count; k++)
        {
            var conf = ToConf(poses.Slice(k * PoseLength, PoseLength));
            energies[k] = engine.Bfgs(conf, maxSteps, cap);
            ToPose(conf).CopyTo(poses.Slice(k * PoseLength, PoseLength));
        }
        return energies;
    }

    /// <summary>The same quantities from the CPU docking engine in double precision (the reference for parity checks).</summary>
    public double[] EvaluateCpu(ReadOnlySpan<double> poses, double cap = 1000, Span<double> gradients = default)
    {
        int count = CheckPoses(poses, gradients);
        var energies = new double[count];
        var engine = new VinaDocking.Engine(Ligand, Grid, _weights);
        var g = new double[GradientLength];
        for (int k = 0; k < count; k++)
        {
            energies[k] = engine.Evaluate(ToConf(poses.Slice(k * PoseLength, PoseLength)), cap, g);
            if (!gradients.IsEmpty) g.CopyTo(gradients.Slice(k * GradientLength, GradientLength));
        }
        return energies;
    }

    private int CheckPoses(ReadOnlySpan<double> poses, Span<double> gradients)
    {
        if (poses.Length % PoseLength != 0) throw new ArgumentException($"Poses come in blocks of {PoseLength} values.", nameof(poses));
        int count = poses.Length / PoseLength;
        if (!gradients.IsEmpty && gradients.Length < count * GradientLength)
            throw new ArgumentException($"The gradient buffer needs {count * GradientLength} values.", nameof(gradients));
        return count;
    }

    internal LigandConf ToConf(ReadOnlySpan<double> pose)
    {
        var tors = pose.Slice(7, TorsionCount).ToArray();
        return new LigandConf(new Vec3(pose[0], pose[1], pose[2]), new Quat(pose[3], pose[4], pose[5], pose[6]), tors);
    }

    internal double[] ToPose(LigandConf c)
    {
        var p = new double[PoseLength];
        p[0] = c.Position.X; p[1] = c.Position.Y; p[2] = c.Position.Z;
        p[3] = c.Orientation.W; p[4] = c.Orientation.X; p[5] = c.Orientation.Y; p[6] = c.Orientation.Z;
        c.Torsions.CopyTo(p, 7);
        return p;
    }

    /// <summary>Releases the device buffers (the device itself stays open).</summary>
    public void Dispose() => _model?.Dispose();
}
