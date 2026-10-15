namespace Biotech.Net.StructuralBiology;

/// <summary>
/// Fixed-radius neighbour search over a static point set using a cell list (uniform grid). Points are bucketed by
/// counting sort into contiguous arrays, so building is O(n) and a query visits only the 27 cells around the query
/// point. Thread-safe for concurrent queries after construction.
/// </summary>
public sealed class NeighborSearch
{
    private const int MaxCells = 1 << 22;
    private readonly Vec3[] _points;
    private readonly int[] _cellStart;   // length cells + 1
    private readonly int[] _sorted;      // point indices ordered by cell
    private readonly Vec3 _origin;
    private readonly double _cell;
    private readonly int _nx, _ny, _nz;

    /// <summary>Indexes <paramref name="points"/> (copied) with a cell edge of <paramref name="cellSize"/> Å.</summary>
    /// <param name="cellSize">Use the largest radius you will query with (larger radii still work, just slower).</param>
    public NeighborSearch(ReadOnlySpan<Vec3> points, double cellSize = 5.0)
    {
        if (!(cellSize > 0)) throw new ArgumentOutOfRangeException(nameof(cellSize));
        _points = points.ToArray();
        if (_points.Length == 0)
        {
            _origin = Vec3.Zero; _cell = cellSize; _nx = _ny = _nz = 1;
            _cellStart = new int[2]; _sorted = [];
            return;
        }
        var (min, max) = Geometry.BoundingBox(_points);
        var ext = max - min;
        // Grow the cell when a sparse, wide point set would need too many cells.
        while ((long)((int)(ext.X / cellSize) + 1) * ((int)(ext.Y / cellSize) + 1) * ((int)(ext.Z / cellSize) + 1) > MaxCells) cellSize *= 1.5;
        _cell = cellSize;
        _origin = min;
        _nx = (int)(ext.X / cellSize) + 1; _ny = (int)(ext.Y / cellSize) + 1; _nz = (int)(ext.Z / cellSize) + 1;
        int cells = _nx * _ny * _nz;
        var cellOf = new int[_points.Length];
        _cellStart = new int[cells + 1];
        for (int i = 0; i < _points.Length; i++)
        {
            int c = CellIndex(_points[i]);
            cellOf[i] = c;
            _cellStart[c + 1]++;
        }
        for (int c = 0; c < cells; c++) _cellStart[c + 1] += _cellStart[c];
        var fill = new int[cells];
        _sorted = new int[_points.Length];
        for (int i = 0; i < _points.Length; i++)
        {
            int c = cellOf[i];
            _sorted[_cellStart[c] + fill[c]++] = i;
        }
    }

    /// <summary>Indexes the positions of <paramref name="atoms"/>; result indices refer to that list.</summary>
    public static NeighborSearch ForAtoms(IReadOnlyList<Atom> atoms, double cellSize = 5.0)
        => new(Superposition.Positions(atoms), cellSize);

    /// <summary>Number of indexed points.</summary>
    public int Count => _points.Length;

    /// <summary>The indexed points.</summary>
    public ReadOnlySpan<Vec3> Points => _points;

    /// <summary>Appends to <paramref name="results"/> the indices of points within <paramref name="radius"/> of <paramref name="center"/>.</summary>
    public void Query(Vec3 center, double radius, List<int> results)
    {
        if (_points.Length == 0) return;
        double r2 = radius * radius;
        int reach = (int)Math.Ceiling(radius / _cell);
        Cell(center, out int cx, out int cy, out int cz);
        int x0 = Math.Max(0, cx - reach), x1 = Math.Min(_nx - 1, cx + reach);
        int y0 = Math.Max(0, cy - reach), y1 = Math.Min(_ny - 1, cy + reach);
        int z0 = Math.Max(0, cz - reach), z1 = Math.Min(_nz - 1, cz + reach);
        for (int x = x0; x <= x1; x++)
            for (int y = y0; y <= y1; y++)
                for (int z = z0; z <= z1; z++)
                {
                    int c = (x * _ny + y) * _nz + z;
                    for (int k = _cellStart[c]; k < _cellStart[c + 1]; k++)
                    {
                        int i = _sorted[k];
                        if (Vec3.DistanceSquared(_points[i], center) <= r2) results.Add(i);
                    }
                }
    }

    /// <summary>Indices of points within <paramref name="radius"/> of <paramref name="center"/> (ascending).</summary>
    public int[] Query(Vec3 center, double radius)
    {
        var list = new List<int>();
        Query(center, radius, list);
        list.Sort();
        return list.ToArray();
    }

    /// <summary>Whether any point lies within <paramref name="radius"/> of <paramref name="center"/>.</summary>
    public bool Any(Vec3 center, double radius)
    {
        if (_points.Length == 0) return false;
        double r2 = radius * radius;
        int reach = (int)Math.Ceiling(radius / _cell);
        Cell(center, out int cx, out int cy, out int cz);
        for (int x = Math.Max(0, cx - reach); x <= Math.Min(_nx - 1, cx + reach); x++)
            for (int y = Math.Max(0, cy - reach); y <= Math.Min(_ny - 1, cy + reach); y++)
                for (int z = Math.Max(0, cz - reach); z <= Math.Min(_nz - 1, cz + reach); z++)
                {
                    int c = (x * _ny + y) * _nz + z;
                    for (int k = _cellStart[c]; k < _cellStart[c + 1]; k++)
                        if (Vec3.DistanceSquared(_points[_sorted[k]], center) <= r2) return true;
                }
        return false;
    }

    /// <summary>Index of the nearest point (−1 when empty) and its distance, searching up to <paramref name="maxRadius"/>.</summary>
    public (int Index, double Distance) Nearest(Vec3 center, double maxRadius = double.PositiveInfinity)
    {
        int best = -1; double bestD2 = double.PositiveInfinity;
        if (_points.Length == 0) return (-1, double.PositiveInfinity);
        double limit = double.IsPositiveInfinity(maxRadius) ? (_nx + _ny + _nz) * _cell + 1e9 : maxRadius;
        var list = new List<int>();
        for (double r = _cell; ; r *= 2)
        {
            list.Clear();
            Query(center, Math.Min(r, limit), list);
            foreach (int i in list)
            {
                double d2 = Vec3.DistanceSquared(_points[i], center);
                if (d2 < bestD2) { bestD2 = d2; best = i; }
            }
            if (best >= 0 || r >= limit) break;
            if (double.IsPositiveInfinity(maxRadius) && r > 1e7) break;
        }
        return (best, Math.Sqrt(bestD2));
    }

    /// <summary>All pairs (i &lt; j) closer than <paramref name="radius"/>, found cell by cell (no O(n²) scan).</summary>
    public List<(int I, int J)> AllPairs(double radius)
    {
        var pairs = new List<(int, int)>();
        double r2 = radius * radius;
        int reach = (int)Math.Ceiling(radius / _cell);
        for (int x = 0; x < _nx; x++)
            for (int y = 0; y < _ny; y++)
                for (int z = 0; z < _nz; z++)
                {
                    int c = (x * _ny + y) * _nz + z;
                    if (_cellStart[c] == _cellStart[c + 1]) continue;
                    for (int ox = x; ox <= Math.Min(_nx - 1, x + reach); ox++)
                        for (int oy = Math.Max(0, y - reach); oy <= Math.Min(_ny - 1, y + reach); oy++)
                            for (int oz = Math.Max(0, z - reach); oz <= Math.Min(_nz - 1, z + reach); oz++)
                            {
                                int o = (ox * _ny + oy) * _nz + oz;
                                // Visit each unordered cell pair once: lexicographically (ox, oy, oz) >= (x, y, z).
                                if (o < c) continue;
                                for (int k = _cellStart[c]; k < _cellStart[c + 1]; k++)
                                {
                                    int i = _sorted[k];
                                    int l0 = o == c ? k + 1 : _cellStart[o];
                                    for (int l = l0; l < _cellStart[o + 1]; l++)
                                    {
                                        int j = _sorted[l];
                                        if (Vec3.DistanceSquared(_points[i], _points[j]) <= r2)
                                            pairs.Add(i < j ? (i, j) : (j, i));
                                    }
                                }
                            }
                }
        return pairs;
    }

    private int CellIndex(Vec3 p)
    {
        Cell(p, out int x, out int y, out int z);
        return (x * _ny + y) * _nz + z;
    }

    private void Cell(Vec3 p, out int x, out int y, out int z)
    {
        x = Math.Clamp((int)Math.Floor((p.X - _origin.X) / _cell), 0, _nx - 1);
        y = Math.Clamp((int)Math.Floor((p.Y - _origin.Y) / _cell), 0, _ny - 1);
        z = Math.Clamp((int)Math.Floor((p.Z - _origin.Z) / _cell), 0, _nz - 1);
    }
}
