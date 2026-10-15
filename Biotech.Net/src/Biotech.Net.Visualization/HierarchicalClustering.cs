namespace Biotech.Net.Visualization;

/// <summary>Distance between two observation vectors.</summary>
public enum DistanceMetric
{
    /// <summary>Euclidean distance.</summary>
    Euclidean,

    /// <summary>1 − Pearson correlation (the usual choice for expression profiles).</summary>
    Correlation,

    /// <summary>Manhattan (city-block) distance.</summary>
    Manhattan,

    /// <summary>1 − cosine similarity.</summary>
    Cosine,
}

/// <summary>Agglomerative linkage criterion.</summary>
public enum LinkageMethod
{
    Single,
    Complete,
    Average,

    /// <summary>Ward's minimum-variance criterion (Lance–Williams form, as in SciPy).</summary>
    Ward,
}

/// <summary>One merge of a dendrogram (SciPy linkage-matrix convention).</summary>
/// <param name="Left">Cluster id: leaves are 0..n−1, the cluster formed at merge k is n+k.</param>
/// <param name="Right">Second cluster id.</param>
/// <param name="Height">Linkage distance.</param>
/// <param name="Size">Number of leaves in the merged cluster.</param>
public readonly record struct DendrogramMerge(int Left, int Right, double Height, int Size);

/// <summary>The result of hierarchical clustering: merges in increasing height and the leaf display order.</summary>
public sealed class Dendrogram
{
    internal Dendrogram(int leafCount, DendrogramMerge[] merges)
    {
        LeafCount = leafCount;
        Merges = merges;
        Order = ComputeOrder();
    }

    /// <summary>Number of observations.</summary>
    public int LeafCount { get; }

    /// <summary>n − 1 merges sorted by height.</summary>
    public IReadOnlyList<DendrogramMerge> Merges { get; }

    /// <summary>Leaves in left-to-right display order.</summary>
    public int[] Order { get; }

    /// <summary>Height of the root (0 for a single leaf).</summary>
    public double RootHeight => Merges.Count == 0 ? 0 : Merges[^1].Height;

    /// <summary>Cuts the tree into <paramref name="k"/> clusters; returns a 0-based cluster label per leaf.</summary>
    public int[] Cut(int k)
    {
        k = Math.Clamp(k, 1, Math.Max(1, LeafCount));
        var parent = Enumerable.Range(0, LeafCount + Merges.Count).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        for (int m = 0; m < Merges.Count - (k - 1); m++)
        {
            int id = LeafCount + m;
            parent[Find(Merges[m].Left)] = id;
            parent[Find(Merges[m].Right)] = id;
        }
        var labels = new int[LeafCount];
        var map = new Dictionary<int, int>();
        foreach (int leaf in Order)
        {
            int root = Find(leaf);
            if (!map.TryGetValue(root, out int l)) map[root] = l = map.Count;
            labels[leaf] = l;
        }
        return labels;
    }

    /// <summary>
    /// Line segments that draw the tree: x in leaf units (leaf i of <see cref="Order"/> sits at i + 0.5), y = height.
    /// Each merge yields three segments (two legs and the bar).
    /// </summary>
    public IReadOnlyList<(double X1, double Y1, double X2, double Y2)> Segments()
    {
        var pos = new double[LeafCount + Merges.Count];
        var height = new double[LeafCount + Merges.Count];
        for (int i = 0; i < Order.Length; i++) pos[Order[i]] = i + 0.5;
        var segs = new List<(double, double, double, double)>(Merges.Count * 3);
        for (int m = 0; m < Merges.Count; m++)
        {
            var mg = Merges[m];
            int id = LeafCount + m;
            double xl = pos[mg.Left], xr = pos[mg.Right], h = mg.Height;
            segs.Add((xl, height[mg.Left], xl, h));
            segs.Add((xr, height[mg.Right], xr, h));
            segs.Add((xl, h, xr, h));
            pos[id] = (xl + xr) / 2;
            height[id] = h;
        }
        return segs;
    }

    private int[] ComputeOrder()
    {
        if (LeafCount == 0) return [];
        if (Merges.Count == 0) return Enumerable.Range(0, LeafCount).ToArray();
        var order = new List<int>(LeafCount);
        var stack = new Stack<int>();
        stack.Push(LeafCount + Merges.Count - 1);
        while (stack.Count > 0)
        {
            int id = stack.Pop();
            if (id < LeafCount) { order.Add(id); continue; }
            var m = Merges[id - LeafCount];
            stack.Push(m.Right);
            stack.Push(m.Left);
        }
        // Disconnected leaves (never happens for a full tree) keep input order.
        if (order.Count < LeafCount) order.AddRange(Enumerable.Range(0, LeafCount).Except(order));
        return order.ToArray();
    }
}

/// <summary>
/// Agglomerative hierarchical clustering with the nearest-neighbour-chain algorithm (O(n²) time, Müllner 2011)
/// and Lance–Williams updates — the same results as SciPy's <c>linkage</c>.
/// </summary>
public static class HierarchicalClustering
{
    /// <summary>Clusters the rows of <paramref name="observations"/>.</summary>
    public static Dendrogram Cluster(IReadOnlyList<double[]> observations, DistanceMetric metric = DistanceMetric.Euclidean,
        LinkageMethod linkage = LinkageMethod.Average)
        => Cluster(DistanceMatrix(observations, metric), linkage);

    /// <summary>Full symmetric distance matrix between rows.</summary>
    public static double[,] DistanceMatrix(IReadOnlyList<double[]> observations, DistanceMetric metric = DistanceMetric.Euclidean)
    {
        int n = observations.Count;
        var d = new double[n, n];
        double[][] prepared = metric switch
        {
            DistanceMetric.Correlation => observations.Select(Standardise).ToArray(),
            DistanceMetric.Cosine => observations.Select(Normalise).ToArray(),
            _ => observations.ToArray(),
        };
        Parallel.For(0, n, i =>
        {
            for (int j = i + 1; j < n; j++)
            {
                double v = Distance(prepared[i], prepared[j], metric);
                d[i, j] = v;
                d[j, i] = v;
            }
        });
        return d;
    }

    /// <summary>Clusters from a precomputed symmetric distance matrix (not modified).</summary>
    public static Dendrogram Cluster(double[,] distances, LinkageMethod linkage = LinkageMethod.Average)
    {
        int n = distances.GetLength(0);
        if (distances.GetLength(1) != n) throw new ArgumentException("The distance matrix must be square.");
        if (n <= 1) return new Dendrogram(n, []);
        var d = (double[,])distances.Clone();
        var size = new int[n];
        Array.Fill(size, 1);
        var active = new bool[n];
        Array.Fill(active, true);
        var raw = new List<(int A, int B, double H)>(n - 1);
        var chain = new List<int>(n);

        for (int step = 0; step < n - 1; step++)
        {
            if (chain.Count == 0)
                for (int i = 0; i < n; i++) if (active[i]) { chain.Add(i); break; }
            int a, b;
            while (true)
            {
                a = chain[^1];
                int prev = chain.Count > 1 ? chain[^2] : -1;
                // Nearest active neighbour of a; prefer the previous chain element on ties (guarantees termination).
                b = prev;
                double best = prev >= 0 ? d[a, prev] : double.PositiveInfinity;
                for (int j = 0; j < n; j++)
                {
                    if (!active[j] || j == a) continue;
                    double v = d[a, j];
                    if (v < best) { best = v; b = j; }
                }
                if (b == prev && prev >= 0) break;
                chain.Add(b);
            }
            chain.RemoveAt(chain.Count - 1);
            chain.RemoveAt(chain.Count - 1);
            double h = d[a, b];
            if (a > b) (a, b) = (b, a);
            raw.Add((a, b, h));

            // Merge b into a (a keeps representing the new cluster).
            int na = size[a], nb = size[b];
            for (int k = 0; k < n; k++)
            {
                if (!active[k] || k == a || k == b) continue;
                double dka = d[k, a], dkb = d[k, b], nk = size[k];
                double v = linkage switch
                {
                    LinkageMethod.Single => Math.Min(dka, dkb),
                    LinkageMethod.Complete => Math.Max(dka, dkb),
                    LinkageMethod.Average => (na * dka + nb * dkb) / (na + nb),
                    _ => Math.Sqrt(Math.Max(0, ((na + nk) * dka * dka + (nb + nk) * dkb * dkb - nk * h * h) / (na + nb + nk))),
                };
                d[k, a] = v;
                d[a, k] = v;
            }
            active[b] = false;
            size[a] = na + nb;
        }

        // Sort merges by height (stable) and relabel with union-find, SciPy style.
        var sorted = raw.Select((m, i) => (m, i)).OrderBy(t => t.m.H).ThenBy(t => t.i).Select(t => t.m).ToList();
        var parent = Enumerable.Range(0, 2 * n - 1).ToArray();
        var clusterSize = new int[2 * n - 1];
        for (int i = 0; i < n; i++) clusterSize[i] = 1;
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        var merges = new DendrogramMerge[n - 1];
        for (int k = 0; k < sorted.Count; k++)
        {
            int ra = Find(sorted[k].A), rb = Find(sorted[k].B);
            int id = n + k;
            if (ra > rb) (ra, rb) = (rb, ra);
            parent[ra] = id;
            parent[rb] = id;
            clusterSize[id] = clusterSize[ra] + clusterSize[rb];
            merges[k] = new DendrogramMerge(ra, rb, sorted[k].H, clusterSize[id]);
        }
        return new Dendrogram(n, merges);
    }

    private static double Distance(double[] x, double[] y, DistanceMetric metric)
    {
        int len = Math.Min(x.Length, y.Length);
        double s = 0;
        switch (metric)
        {
            case DistanceMetric.Manhattan:
                for (int i = 0; i < len; i++) s += Math.Abs(x[i] - y[i]);
                return s;
            case DistanceMetric.Correlation:
            case DistanceMetric.Cosine:
                for (int i = 0; i < len; i++) s += x[i] * y[i];
                return Math.Max(0, 1 - s);
            default:
                for (int i = 0; i < len; i++) { double t = x[i] - y[i]; s += t * t; }
                return Math.Sqrt(s);
        }
    }

    private static double[] Standardise(double[] x)
    {
        double mean = x.Average();
        var c = x.Select(v => v - mean).ToArray();
        return Normalise(c);
    }

    private static double[] Normalise(double[] x)
    {
        double norm = Math.Sqrt(x.Sum(v => v * v));
        return norm > 0 ? x.Select(v => v / norm).ToArray() : new double[x.Length];
    }
}
