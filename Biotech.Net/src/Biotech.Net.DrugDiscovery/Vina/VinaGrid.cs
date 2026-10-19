using Biotech.Net.StructuralBiology;

namespace Biotech.Net.DrugDiscovery;

/// <summary>Unweighted Vina term values for one atom pair (or a sum over pairs).</summary>
/// <param name="Gauss1">Σ exp(−(d/0.5)²).</param>
/// <param name="Gauss2">Σ exp(−((d−3)/2)²).</param>
/// <param name="Repulsion">Σ d² for d &lt; 0.</param>
/// <param name="Hydrophobic">Σ hydrophobic step terms.</param>
/// <param name="HydrogenBond">Σ hydrogen-bond step terms.</param>
public readonly record struct VinaTermSums(double Gauss1, double Gauss2, double Repulsion, double Hydrophobic, double HydrogenBond)
{
    /// <summary>Weighted sum with the given weights.</summary>
    public double Weighted(VinaScoreOptions w)
        => w.Gauss1Weight * Gauss1 + w.Gauss2Weight * Gauss2 + w.RepulsionWeight * Repulsion + w.HydrophobicWeight * Hydrophobic + w.HydrogenBondWeight * HydrogenBond;

    /// <summary>Element-wise sum.</summary>
    public static VinaTermSums operator +(VinaTermSums a, VinaTermSums b)
        => new(a.Gauss1 + b.Gauss1, a.Gauss2 + b.Gauss2, a.Repulsion + b.Repulsion, a.Hydrophobic + b.Hydrophobic, a.HydrogenBond + b.HydrogenBond);
}

/// <summary>The five distance-dependent Vina terms with analytic derivatives, on X-Score types.</summary>
internal static class VinaPotential
{
    public static VinaTermSums Terms(XsAtomType a, XsAtomType b, double r)
    {
        double d = r - VinaAtomTypes.Radius(a) - VinaAtomTypes.Radius(b);
        double g1 = Math.Exp(-(d / 0.5) * (d / 0.5));
        double g2 = Math.Exp(-((d - 3.0) / 2.0) * ((d - 3.0) / 2.0));
        double rep = d < 0 ? d * d : 0;
        double phob = VinaAtomTypes.IsHydrophobic(a) && VinaAtomTypes.IsHydrophobic(b) ? (d < 0.5 ? 1 : d > 1.5 ? 0 : 1.5 - d) : 0;
        double hb = VinaAtomTypes.HBondPossible(a, b) ? (d < -0.7 ? 1 : d > 0 ? 0 : -d / 0.7) : 0;
        return new VinaTermSums(g1, g2, rep, phob, hb);
    }

    /// <summary>Weighted pair energy and dE/dr (zero beyond the cut-off).</summary>
    public static double Energy(XsAtomType a, XsAtomType b, double r, VinaScoreOptions w, out double dEdr)
    {
        dEdr = 0;
        if (r >= w.CutOff) return 0;
        double d = r - VinaAtomTypes.Radius(a) - VinaAtomTypes.Radius(b);
        double g1 = Math.Exp(-(d / 0.5) * (d / 0.5));
        double g2 = Math.Exp(-((d - 3.0) / 2.0) * ((d - 3.0) / 2.0));
        double e = w.Gauss1Weight * g1 + w.Gauss2Weight * g2;
        dEdr = w.Gauss1Weight * g1 * (-2 * d / 0.25) + w.Gauss2Weight * g2 * (-2 * (d - 3.0) / 4.0);
        if (d < 0) { e += w.RepulsionWeight * d * d; dEdr += w.RepulsionWeight * 2 * d; }
        if (VinaAtomTypes.IsHydrophobic(a) && VinaAtomTypes.IsHydrophobic(b))
        {
            if (d < 0.5) e += w.HydrophobicWeight;
            else if (d <= 1.5) { e += w.HydrophobicWeight * (1.5 - d); dEdr -= w.HydrophobicWeight; }
        }
        if (VinaAtomTypes.HBondPossible(a, b))
        {
            if (d < -0.7) e += w.HydrogenBondWeight;
            else if (d <= 0) { e += w.HydrogenBondWeight * (-d / 0.7); dEdr -= w.HydrogenBondWeight / 0.7; }
        }
        return e;
    }

    /// <summary>Vina's "curl": caps positive energies smoothly at <paramref name="v"/> (no-op for v ≥ 1000).</summary>
    public static double Curl(double e, double v, ref double deriv)
    {
        if (e > 0 && v < 1000)
        {
            double t = v / (v + e);
            e *= t;
            deriv *= t * t;
        }
        return e;
    }
}

/// <summary>
/// Precomputed Vina grid maps: for each ligand X-Score type, the receptor interaction energy on a regular grid
/// (default spacing 0.375 Å, as in Vina), evaluated with trilinear interpolation. Atoms outside the box pay a steep
/// linear penalty that pulls them back.
/// </summary>
public sealed class VinaGrid
{
    internal const double OutsideSlope = 1e6;
    private readonly float[]?[] _maps = new float[VinaAtomTypes.Count][];
    private readonly int _nx, _ny, _nz;

    /// <summary>Lower corner of the box.</summary>
    public Vec3 Origin { get; }

    /// <summary>Box edge lengths (Å).</summary>
    public Vec3 Size { get; }

    /// <summary>Grid spacing (Å).</summary>
    public double Spacing { get; }

    /// <summary>Box centre.</summary>
    public Vec3 Center => Origin + Size * 0.5;

    /// <summary>Grid points along each axis.</summary>
    public (int X, int Y, int Z) Dimensions => (_nx, _ny, _nz);

    private VinaGrid(Vec3 origin, Vec3 size, double spacing)
    {
        Origin = origin;
        Size = size;
        Spacing = spacing;
        _nx = (int)Math.Ceiling(size.X / spacing) + 1;
        _ny = (int)Math.Ceiling(size.Y / spacing) + 1;
        _nz = (int)Math.Ceiling(size.Z / spacing) + 1;
    }

    /// <summary>Builds the maps for the given ligand types (parallel over grid slabs).</summary>
    public static VinaGrid Build(VinaReceptor receptor, Vec3 center, Vec3 size, IEnumerable<XsAtomType> ligandTypes,
        double spacing = 0.375, VinaScoreOptions? weights = null)
    {
        weights ??= VinaScoreOptions.Default;
        var g = new VinaGrid(center - size * 0.5, size, spacing);
        var types = ligandTypes.Distinct().ToArray();
        // tabulate pair energies over r for every (receptor type, ligand type), 0.005 Å resolution
        const double dr = 0.005;
        int nr = (int)(weights.CutOff / dr) + 2;
        var table = new double[VinaAtomTypes.Count][][];
        var recTypes = receptor.Types.Distinct().ToArray();
        foreach (var rt in recTypes)
        {
            table[(int)rt] = new double[types.Length][];
            for (int k = 0; k < types.Length; k++)
            {
                var row = new double[nr];
                for (int s = 0; s < nr; s++) row[s] = VinaPotential.Energy(rt, types[k], s * dr, weights, out _);
                table[(int)rt][k] = row;
            }
        }
        foreach (var t in types) g._maps[(int)t] = new float[g._nx * g._ny * g._nz];
        // receptor atoms near the box only
        double reach = weights.CutOff;
        var lo = g.Origin - new Vec3(reach, reach, reach);
        var hi = g.Origin + size + new Vec3(reach, reach, reach);
        var idx = new List<int>();
        for (int i = 0; i < receptor.Count; i++)
        {
            var p = receptor.Positions[i];
            if (p.X >= lo.X && p.X <= hi.X && p.Y >= lo.Y && p.Y <= hi.Y && p.Z >= lo.Z && p.Z <= hi.Z) idx.Add(i);
        }
        if (idx.Count == 0) return g;
        var pts = idx.Select(i => receptor.Positions[i]).ToArray();
        var rtype = idx.Select(i => receptor.Types[i]).ToArray();
        var search = new NeighborSearch(pts, reach);
        double cut2 = weights.CutOff * weights.CutOff;
        Parallel.For(0, g._nx, ix =>
        {
            var found = new List<int>();
            var acc = new double[types.Length];
            for (int iy = 0; iy < g._ny; iy++)
                for (int iz = 0; iz < g._nz; iz++)
                {
                    var p = g.Origin + new Vec3(ix * spacing, iy * spacing, iz * spacing);
                    found.Clear();
                    search.Query(p, reach, found);
                    Array.Clear(acc);
                    foreach (int j in found)
                    {
                        double r2 = Vec3.DistanceSquared(p, pts[j]);
                        if (r2 >= cut2) continue;
                        double x = Math.Sqrt(r2) / dr;
                        int s = (int)x;
                        double f = x - s;
                        var rows = table[(int)rtype[j]];
                        for (int k = 0; k < types.Length; k++) acc[k] += rows[k][s] * (1 - f) + rows[k][s + 1] * f;
                    }
                    int cell = (ix * g._ny + iy) * g._nz + iz;
                    for (int k = 0; k < types.Length; k++) g._maps[(int)types[k]]![cell] = (float)acc[k];
                }
        });
        return g;
    }

    /// <summary>The raw map of a type (x-major: <c>(ix·ny + iy)·nz + iz</c>), or null; used by the GPU backend.</summary>
    internal float[]? MapFor(XsAtomType t) => _maps[(int)t];

    /// <summary>Whether a map exists for the type.</summary>
    public bool HasType(XsAtomType t) => _maps[(int)t] is not null;

    /// <summary>Interpolated energy (capped by <paramref name="cap"/>) plus the outside-box penalty, and its gradient.</summary>
    public double Evaluate(XsAtomType type, Vec3 p, double cap, out Vec3 gradient)
    {
        var map = _maps[(int)type] ?? throw new InvalidOperationException($"No grid map for type {type}.");
        double inv = 1.0 / Spacing;
        Span<double> s = [(p.X - Origin.X) * inv, (p.Y - Origin.Y) * inv, (p.Z - Origin.Z) * inv];
        Span<int> dims = [_nx, _ny, _nz];
        Span<int> a = stackalloc int[3];
        Span<int> region = stackalloc int[3];
        double penalty = 0;
        for (int k = 0; k < 3; k++)
        {
            if (s[k] < 0) { penalty += -s[k] * Spacing; region[k] = -1; a[k] = 0; s[k] = 0; }
            else if (s[k] >= dims[k] - 1) { penalty += (s[k] - (dims[k] - 1)) * Spacing; region[k] = 1; a[k] = dims[k] - 2; s[k] = 1; }
            else { region[k] = 0; a[k] = (int)s[k]; s[k] -= a[k]; }
        }
        penalty *= OutsideSlope;
        int Id(int x, int y, int z) => (x * _ny + y) * _nz + z;
        double f000 = map[Id(a[0], a[1], a[2])], f100 = map[Id(a[0] + 1, a[1], a[2])];
        double f010 = map[Id(a[0], a[1] + 1, a[2])], f110 = map[Id(a[0] + 1, a[1] + 1, a[2])];
        double f001 = map[Id(a[0], a[1], a[2] + 1)], f101 = map[Id(a[0] + 1, a[1], a[2] + 1)];
        double f011 = map[Id(a[0], a[1] + 1, a[2] + 1)], f111 = map[Id(a[0] + 1, a[1] + 1, a[2] + 1)];
        double x = s[0], y = s[1], z = s[2], mx = 1 - x, my = 1 - y, mz = 1 - z;
        double f = f000 * mx * my * mz + f100 * x * my * mz + f010 * mx * y * mz + f110 * x * y * mz
                 + f001 * mx * my * z + f101 * x * my * z + f011 * mx * y * z + f111 * x * y * z;
        double gx = (f100 - f000) * my * mz + (f110 - f010) * y * mz + (f101 - f001) * my * z + (f111 - f011) * y * z;
        double gy = (f010 - f000) * mx * mz + (f110 - f100) * x * mz + (f011 - f001) * mx * z + (f111 - f101) * x * z;
        double gz = (f001 - f000) * mx * my + (f101 - f100) * x * my + (f011 - f010) * mx * y + (f111 - f110) * x * y;
        double scale = 1.0;
        f = VinaPotential.Curl(f, cap, ref scale);
        gx *= scale * inv; gy *= scale * inv; gz *= scale * inv;
        gradient = new Vec3(region[0] == 0 ? gx : OutsideSlope * region[0], region[1] == 0 ? gy : OutsideSlope * region[1], region[2] == 0 ? gz : OutsideSlope * region[2]);
        return f + penalty;
    }
}
