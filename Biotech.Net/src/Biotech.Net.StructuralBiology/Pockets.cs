namespace Biotech.Net.StructuralBiology;

/// <summary>Options for <see cref="PocketDetector"/>.</summary>
public sealed record PocketOptions
{
    /// <summary>Grid spacing in Å (default 1.0).</summary>
    public double GridSpacing { get; init; } = 1.0;

    /// <summary>Minimum number of the 7 scan directions that must hit protein on both sides (default 5).</summary>
    public int MinBuriedness { get; init; } = 5;

    /// <summary>How far each direction is scanned, in Å (default 10).</summary>
    public double ScanDistance { get; init; } = 10.0;

    /// <summary>Clusters smaller than this many grid points are discarded (default 10).</summary>
    public int MinPocketPoints { get; init; } = 10;

    /// <summary>Extra distance added to van der Waals radii when marking protein grid points (default 0.0 Å).</summary>
    public double RadiusPadding { get; init; }

    /// <summary>Atoms within this distance of a pocket point line the pocket (default 4.0 Å).</summary>
    public double LiningDistance { get; init; } = 4.0;

    /// <summary>Ignore HETATM atoms (ligands, waters, ions) so pockets of holo structures can be found (default true).</summary>
    public bool IgnoreHetero { get; init; } = true;
}

/// <summary>A detected pocket: a connected cluster of buried solvent grid points.</summary>
public sealed record Pocket(int Rank, IReadOnlyList<Vec3> Points, Vec3 Center, double Volume, double MeanBuriedness,
    IReadOnlyList<Residue> LiningResidues)
{
    /// <summary>Whether any pocket point is within <paramref name="distance"/> Å of <paramref name="p"/>.</summary>
    public bool Contains(Vec3 p, double distance = 1.5)
    {
        double d2 = distance * distance;
        foreach (var q in Points) if (Vec3.DistanceSquared(p, q) <= d2) return true;
        return false;
    }

    public override string ToString() => $"Pocket {Rank}: {Volume:F0} Å³, {LiningResidues.Count} lining residues, centre {Center}";
}

/// <summary>
/// Grid-based binding-site detection in the spirit of LIGSITE (Hendlich, Rippmann &amp; Barnickel 1997): the protein
/// is mapped onto a grid, and every solvent grid point is scanned along the three axes and four cube diagonals.
/// A direction counts as a protein–solvent–protein (PSP) event when protein is met on both sides within
/// <see cref="PocketOptions.ScanDistance"/>. Points with at least <see cref="PocketOptions.MinBuriedness"/> PSP
/// events are pocket points; they are clustered by 26-connectivity and ranked by volume.
/// </summary>
public static class PocketDetector
{
    private static readonly (int X, int Y, int Z)[] Directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)];

    /// <summary>Detects pockets in the first model.</summary>
    public static IReadOnlyList<Pocket> Detect(Structure structure, PocketOptions? options = null) => Detect(structure.FirstModel, options);

    /// <summary>Detects pockets in a model.</summary>
    public static IReadOnlyList<Pocket> Detect(Model model, PocketOptions? options = null)
    {
        options ??= new PocketOptions();
        var atoms = model.Atoms.Where(a => a.Element != "H" && !(options.IgnoreHetero && a.IsHetero)).ToList();
        return Detect(atoms, options);
    }

    /// <summary>Detects pockets formed by the given atoms.</summary>
    public static IReadOnlyList<Pocket> Detect(IReadOnlyList<Atom> atoms, PocketOptions? options = null)
    {
        options ??= new PocketOptions();
        if (atoms.Count == 0) return [];
        double s = options.GridSpacing;
        var (min, max) = Geometry.BoundingBox(atoms);
        var origin = min - new Vec3(2, 2, 2);
        int nx = (int)Math.Ceiling((max.X - origin.X + 2) / s) + 1;
        int ny = (int)Math.Ceiling((max.Y - origin.Y + 2) / s) + 1;
        int nz = (int)Math.Ceiling((max.Z - origin.Z + 2) / s) + 1;
        if ((long)nx * ny * nz > 64_000_000) throw new InvalidOperationException("Grid too large; increase GridSpacing.");
        int Index(int x, int y, int z) => (x * ny + y) * nz + z;
        Vec3 Point(int x, int y, int z) => origin + new Vec3(x * s, y * s, z * s);

        // 1. Stamp atom spheres onto the grid.
        var protein = new bool[nx * ny * nz];
        foreach (var a in atoms)
        {
            double r = a.VdwRadius + options.RadiusPadding, r2 = r * r;
            var p = a.Position;
            int x0 = Math.Max(0, (int)Math.Floor((p.X - r - origin.X) / s)), x1 = Math.Min(nx - 1, (int)Math.Ceiling((p.X + r - origin.X) / s));
            int y0 = Math.Max(0, (int)Math.Floor((p.Y - r - origin.Y) / s)), y1 = Math.Min(ny - 1, (int)Math.Ceiling((p.Y + r - origin.Y) / s));
            int z0 = Math.Max(0, (int)Math.Floor((p.Z - r - origin.Z) / s)), z1 = Math.Min(nz - 1, (int)Math.Ceiling((p.Z + r - origin.Z) / s));
            for (int x = x0; x <= x1; x++)
                for (int y = y0; y <= y1; y++)
                    for (int z = z0; z <= z1; z++)
                        if (Vec3.DistanceSquared(Point(x, y, z), p) <= r2) protein[Index(x, y, z)] = true;
        }

        // 2. PSP scan of every solvent point (parallel over x slices).
        var buried = new byte[nx * ny * nz];
        Parallel.For(0, nx, x =>
        {
            for (int y = 0; y < ny; y++)
                for (int z = 0; z < nz; z++)
                {
                    if (protein[Index(x, y, z)]) continue;
                    int count = 0;
                    foreach (var d in Directions)
                    {
                        double stepLength = s * Math.Sqrt(d.X * d.X + d.Y * d.Y + d.Z * d.Z);
                        int maxSteps = (int)(options.ScanDistance / stepLength);
                        if (Hit(x, y, z, d.X, d.Y, d.Z, maxSteps) && Hit(x, y, z, -d.X, -d.Y, -d.Z, maxSteps)) count++;
                    }
                    buried[Index(x, y, z)] = (byte)count;
                }
        });
        bool Hit(int x, int y, int z, int dx, int dy, int dz, int steps)
        {
            for (int k = 1; k <= steps; k++)
            {
                x += dx; y += dy; z += dz;
                if ((uint)x >= (uint)nx || (uint)y >= (uint)ny || (uint)z >= (uint)nz) return false;
                if (protein[Index(x, y, z)]) return true;
            }
            return false;
        }

        // 3. Cluster pocket points (26-connectivity flood fill).
        var visited = new bool[buried.Length];
        var clusters = new List<List<int>>();
        var stack = new Stack<int>();
        for (int i = 0; i < buried.Length; i++)
        {
            if (visited[i] || buried[i] < options.MinBuriedness) continue;
            var cluster = new List<int>();
            visited[i] = true;
            stack.Push(i);
            while (stack.Count > 0)
            {
                int cur = stack.Pop();
                cluster.Add(cur);
                int cx = cur / (ny * nz), cy = cur / nz % ny, cz = cur % nz;
                for (int dx = -1; dx <= 1; dx++)
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            int x = cx + dx, y = cy + dy, z = cz + dz;
                            if ((uint)x >= (uint)nx || (uint)y >= (uint)ny || (uint)z >= (uint)nz) continue;
                            int j = Index(x, y, z);
                            if (visited[j] || buried[j] < options.MinBuriedness) continue;
                            visited[j] = true;
                            stack.Push(j);
                        }
            }
            if (cluster.Count >= options.MinPocketPoints) clusters.Add(cluster);
        }

        // 4. Describe clusters, largest first.
        var ns = NeighborSearch.ForAtoms(atoms, options.LiningDistance);
        var pockets = new List<Pocket>();
        var found = new List<int>();
        int rank = 0;
        foreach (var cluster in clusters.OrderByDescending(c => c.Count))
        {
            var pts = cluster.Select(i => Point(i / (ny * nz), i / nz % ny, i % nz)).ToArray();
            var lining = new HashSet<Residue>();
            foreach (var p in pts)
            {
                found.Clear();
                ns.Query(p, options.LiningDistance, found);
                foreach (int k in found) if (atoms[k].Residue is { } r) lining.Add(r);
            }
            var liningList = lining.OrderBy(r => r.Chain?.Id).ThenBy(r => r.Number).ToList();
            double meanBuried = cluster.Average(i => (double)buried[i]);
            pockets.Add(new Pocket(++rank, pts, Geometry.Centroid(pts), pts.Length * s * s * s, meanBuried, liningList));
        }
        return pockets;
    }
}
