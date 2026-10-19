using Biotech.Net.DrugDiscovery;
using Biotech.Net.StructuralBiology;
using ILGPU;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu.Docking;

/// <summary>A ligand torsion tree plus grid maps uploaded to one device (single precision, box-centred frame).</summary>
internal sealed class VinaGpuModel : IDisposable
{
    private readonly List<MemoryBuffer> _buffers = [];

    public VinaGpuModel(GpuDevice device, FlexibleLigand lig, VinaGrid grid, VinaScoreOptions w)
    {
        Device = device;
        Ligand = lig;
        Grid = grid;
        Frame = grid.Center;
        var acc = device.Accelerator;
        int n = lig.AtomCount;

        var types = lig.Scored.Select(i => lig.Types[i]!.Value).Distinct().ToList();
        var (nx, ny, nz) = grid.Dimensions;
        long cells = (long)nx * ny * nz;
        var maps = new float[types.Count * cells];
        for (int k = 0; k < types.Count; k++)
        {
            var m = grid.MapFor(types[k]) ?? throw new InvalidOperationException($"The grid has no map for {types[k]}.");
            Array.Copy(m, 0, maps, k * cells, cells);
        }
        var local = new float[3 * n];
        var lc = lig.LocalCoordinates;
        for (int i = 0; i < n; i++) { local[3 * i] = (float)lc[i].X; local[3 * i + 1] = (float)lc[i].Y; local[3 * i + 2] = (float)lc[i].Z; }
        var scoredMap = lig.Scored.Select(i => types.IndexOf(lig.Types[i]!.Value)).ToArray();
        var tor = lig.Torsions;
        var moveOff = new int[tor.Length + 1];
        for (int t = 0; t < tor.Length; t++) moveOff[t + 1] = moveOff[t] + tor[t].Moving.Length;
        var moving = tor.SelectMany(t => t.Moving).ToArray();
        var pairs = lig.IntraPairs;
        var flags = pairs.Select(p =>
        {
            XsAtomType a = lig.Types[p.I]!.Value, b = lig.Types[p.J]!.Value;
            return (VinaAtomTypes.IsHydrophobic(a) && VinaAtomTypes.IsHydrophobic(b) ? 1 : 0) | (VinaAtomTypes.HBondPossible(a, b) ? 2 : 0);
        }).ToArray();
        var radius = pairs.Select(p => (float)(VinaAtomTypes.Radius(lig.Types[p.I]!.Value) + VinaAtomTypes.Radius(lig.Types[p.J]!.Value))).ToArray();
        var heavyR = lig.GyrationRadius;

        var origin = grid.Origin - Frame;
        Data = new VinaKernelData
        {
            Maps = Upload(acc, maps),
            Local = Upload(acc, local),
            ScoredAtom = Upload(acc, lig.Scored),
            ScoredMap = Upload(acc, scoredMap),
            TorFrom = Upload(acc, tor.Select(t => t.From).ToArray()),
            TorTo = Upload(acc, tor.Select(t => t.To).ToArray()),
            MoveOffsets = Upload(acc, moveOff),
            MoveAtoms = Upload(acc, moving),
            PairI = Upload(acc, pairs.Select(p => p.I).ToArray()),
            PairJ = Upload(acc, pairs.Select(p => p.J).ToArray()),
            PairFlags = Upload(acc, flags),
            PairRadius = Upload(acc, radius),
            Atoms = n,
            Scored = lig.Scored.Length,
            Torsions = tor.Length,
            Pairs = pairs.Length,
            Nx = nx, Ny = ny, Nz = nz,
            Ox = (float)origin.X, Oy = (float)origin.Y, Oz = (float)origin.Z,
            Spacing = (float)grid.Spacing,
            BoxX = (float)grid.Size.X, BoxY = (float)grid.Size.Y, BoxZ = (float)grid.Size.Z,
            W1 = (float)w.Gauss1Weight, W2 = (float)w.Gauss2Weight, WRep = (float)w.RepulsionWeight,
            WHyd = (float)w.HydrophobicWeight, WHb = (float)w.HydrogenBondWeight, CutOff = (float)w.CutOff,
            GyrationRadius = (float)heavyR,
        };
    }

    public GpuDevice Device { get; }

    public FlexibleLigand Ligand { get; }

    public VinaGrid Grid { get; }

    /// <summary>Origin of the device frame (the box centre).</summary>
    public Vec3 Frame { get; }

    public VinaKernelData Data { get; }

    public int ConfLength => VinaLayout.Conf(Ligand.Torsions.Length);

    public int GradLength => VinaLayout.Grad(Ligand.Torsions.Length);

    /// <summary>Floats of scratch per lane (layout plus the chain's two energies).</summary>
    public int ScratchPerLane => VinaLayout.Size(Ligand.AtomCount, Ligand.Torsions.Length) + 2;

    /// <summary>Absolute pose (double) to a device conformation (float, box frame).</summary>
    public void ToDevice(ReadOnlySpan<double> pose, Span<float> conf)
    {
        conf[0] = (float)(pose[0] - Frame.X);
        conf[1] = (float)(pose[1] - Frame.Y);
        conf[2] = (float)(pose[2] - Frame.Z);
        for (int k = 3; k < ConfLength; k++) conf[k] = (float)pose[k];
    }

    /// <summary>Device conformation back to a CPU conformation.</summary>
    public LigandConf FromDevice(ReadOnlySpan<float> conf)
    {
        var tors = new double[Ligand.Torsions.Length];
        for (int t = 0; t < tors.Length; t++) tors[t] = conf[7 + t];
        var q = new Quat(conf[3], conf[4], conf[5], conf[6]);
        return new LigandConf(new Vec3(conf[0] + Frame.X, conf[1] + Frame.Y, conf[2] + Frame.Z), q.Normalized(), tors);
    }

    private ArrayView<T> Upload<T>(Accelerator acc, T[] data) where T : unmanaged
    {
        var buf = acc.Allocate1D(data.Length == 0 ? new T[1] : data);
        _buffers.Add(buf);
        return buf.View;
    }

    public void Dispose()
    {
        foreach (var b in _buffers) b.Dispose();
        _buffers.Clear();
    }
}
