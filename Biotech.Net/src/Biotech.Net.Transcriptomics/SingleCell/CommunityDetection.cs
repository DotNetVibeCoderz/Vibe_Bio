namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Graph clustering algorithm.</summary>
public enum CommunityAlgorithm
{
    /// <summary>Leiden (Traag, Waltman &amp; van Eck 2019): local moving, refinement, aggregation — guarantees connected communities.</summary>
    Leiden,

    /// <summary>Louvain (Blondel et al. 2008): local moving and aggregation.</summary>
    Louvain,
}

/// <summary>
/// Modularity-based community detection on a weighted undirected graph (e.g. <see cref="NeighborGraph.Connectivities"/>)
/// with a resolution parameter γ, maximising Q = (1/2m) Σ [A_ij − γ k_i k_j / 2m] δ(c_i, c_j) — the quality
/// function scanpy's <c>tl.leiden</c>/<c>tl.louvain</c> use by default (RBConfiguration).
/// </summary>
/// <remarks>
/// Node visiting order and Leiden's randomised refinement (θ = 0.01) use a seeded RNG, so results are reproducible
/// for a given seed. Cluster labels are renumbered by decreasing size (0 = largest), as scanpy does.
/// </remarks>
public static class CommunityDetection
{
    /// <summary>Clusters the nodes of a symmetric weighted adjacency matrix. Returns one label per node.</summary>
    public static int[] Cluster(SparseMatrix adjacency, double resolution = 1.0, CommunityAlgorithm algorithm = CommunityAlgorithm.Leiden,
        int seed = 0, int maxLevels = 50)
    {
        ArgumentNullException.ThrowIfNull(adjacency);
        if (adjacency.Rows != adjacency.Columns) throw new ArgumentException("Adjacency must be square.");
        var rng = new Random(seed);
        var g = Graph.FromSparse(adjacency);
        int n0 = g.N;
        var membership = new int[n0];          // original node -> current aggregate node
        for (int i = 0; i < n0; i++) membership[i] = i;
        var partition = new int[g.N];
        for (int i = 0; i < g.N; i++) partition[i] = i;

        bool finished = false;
        for (int level = 0; level < maxLevels; level++)
        {
            bool moved = algorithm == CommunityAlgorithm.Leiden
                ? FastLocalMove(g, partition, resolution, rng)
                : LocalMove(g, partition, resolution, rng);
            int nComm = Renumber(partition);
            if (!moved || nComm == g.N)
            {
                for (int i = 0; i < n0; i++) membership[i] = partition[membership[i]];
                finished = true;
                break;
            }
            int[] aggregateOf;
            int[] nextPartition;
            if (algorithm == CommunityAlgorithm.Leiden)
            {
                var refined = Refine(g, partition, resolution, rng);
                int nRef = Renumber(refined);
                aggregateOf = refined;
                // Initial partition of the aggregate graph: the non-refined community of each refined community.
                nextPartition = new int[nRef];
                for (int v = 0; v < g.N; v++) nextPartition[refined[v]] = partition[v];
                if (nRef == g.N)
                {
                    // Refinement could not merge anything; fall back to aggregating by the partition.
                    aggregateOf = partition;
                    nextPartition = Enumerable.Range(0, nComm).ToArray();
                }
            }
            else
            {
                aggregateOf = partition;
                nextPartition = Enumerable.Range(0, nComm).ToArray();
            }
            for (int i = 0; i < n0; i++) membership[i] = aggregateOf[membership[i]];
            g = g.Aggregate(aggregateOf, nextPartition.Length);
            partition = nextPartition;
        }
        // Level budget exhausted: membership points at aggregate nodes, whose (Leiden) partition is not the identity.
        if (!finished) for (int i = 0; i < n0; i++) membership[i] = partition[membership[i]];
        return RenumberBySize(membership);
    }

    /// <summary>Clusters <see cref="SingleCellExperiment.Neighbors"/> and stores labels in <c>CellLabels[key]</c>.</summary>
    public static int[] Cluster(SingleCellExperiment sc, double resolution = 1.0, CommunityAlgorithm algorithm = CommunityAlgorithm.Leiden,
        int seed = 0, string? key = null)
    {
        ArgumentNullException.ThrowIfNull(sc);
        var g = sc.Neighbors ?? throw new InvalidOperationException("Compute the neighbour graph first.");
        var labels = Cluster(g.Connectivities, resolution, algorithm, seed);
        sc.CellLabels[key ?? (algorithm == CommunityAlgorithm.Leiden ? "leiden" : "louvain")] =
            labels.Select(l => l.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        return labels;
    }

    /// <summary>Modularity of a labelling (with resolution γ) on a symmetric weighted adjacency.</summary>
    public static double Modularity(SparseMatrix adjacency, ReadOnlySpan<int> labels, double resolution = 1.0)
    {
        ArgumentNullException.ThrowIfNull(adjacency);
        int n = adjacency.Rows;
        double m2 = 0;
        var k = new double[n];
        for (int i = 0; i < n; i++) { foreach (float v in adjacency.RowValues(i)) k[i] += v; m2 += k[i]; }
        if (m2 == 0) return 0;
        int c = 0;
        foreach (int l in labels) c = Math.Max(c, l + 1);
        var tot = new double[c];
        double inside = 0;
        for (int i = 0; i < n; i++)
        {
            tot[labels[i]] += k[i];
            var idx = adjacency.RowIndices(i); var v = adjacency.RowValues(i);
            for (int t = 0; t < idx.Length; t++) if (labels[idx[t]] == labels[i]) inside += v[t];
        }
        double q = inside / m2;
        foreach (double t in tot) q -= resolution * (t / m2) * (t / m2);
        return q;
    }

    // ---- internals -------------------------------------------------------------------------------

    private sealed class Graph
    {
        public int N;
        public int[] Ptr = [];
        public int[] Adj = [];
        public double[] W = [];
        public double[] Strength = [];   // k_i (self-loops counted twice)
        public double[] NodeSize = [];   // number of original nodes
        public double TotalWeight;       // 2m

        public static Graph FromSparse(SparseMatrix a)
        {
            int n = a.Rows;
            var g = new Graph { N = n, Ptr = new int[n + 1], Strength = new double[n], NodeSize = new double[n] };
            var adj = new List<int>(a.NonZeroCount); var w = new List<double>(a.NonZeroCount);
            for (int i = 0; i < n; i++)
            {
                var idx = a.RowIndices(i); var v = a.RowValues(i);
                for (int t = 0; t < idx.Length; t++)
                {
                    g.Strength[i] += idx[t] == i ? 2 * v[t] : v[t];
                    if (idx[t] == i || v[t] == 0) continue;
                    adj.Add(idx[t]); w.Add(v[t]);
                }
                g.Ptr[i + 1] = adj.Count;
                g.NodeSize[i] = 1;
                g.TotalWeight += g.Strength[i];
            }
            g.Adj = adj.ToArray(); g.W = w.ToArray();
            return g;
        }

        public Graph Aggregate(int[] map, int count)
        {
            var g = new Graph { N = count, Ptr = new int[count + 1], Strength = new double[count], NodeSize = new double[count], TotalWeight = TotalWeight };
            var members = new List<int>[count];
            for (int c = 0; c < count; c++) members[c] = [];
            for (int v = 0; v < N; v++) { members[map[v]].Add(v); g.Strength[map[v]] += Strength[v]; g.NodeSize[map[v]] += NodeSize[v]; }
            var adj = new List<int>(); var w = new List<double>();
            var acc = new Dictionary<int, double>();
            for (int c = 0; c < count; c++)
            {
                acc.Clear();
                foreach (int v in members[c])
                    for (int t = Ptr[v]; t < Ptr[v + 1]; t++)
                    {
                        int d = map[Adj[t]];
                        if (d == c) continue;
                        acc[d] = acc.GetValueOrDefault(d) + W[t];
                    }
                foreach (var kv in acc.OrderBy(kv => kv.Key)) { adj.Add(kv.Key); w.Add(kv.Value); }
                g.Ptr[c + 1] = adj.Count;
            }
            g.Adj = adj.ToArray(); g.W = w.ToArray();
            return g;
        }
    }

    /// <summary>Louvain local moving: sweeps nodes in random order until no move improves modularity.</summary>
    private static bool LocalMove(Graph g, int[] part, double gamma, Random rng)
    {
        var tot = new double[g.N];
        for (int v = 0; v < g.N; v++) tot[part[v]] += g.Strength[v];
        var (size, empty) = Occupancy(part, g.N);
        var neighW = new double[g.N];
        var touched = new List<int>();
        bool anyMove = false, improved = true;
        int sweeps = 0;
        while (improved && sweeps++ < 100)
        {
            improved = false;
            foreach (int v in Shuffled(g.N, rng))
            {
                if (MoveNode(g, part, tot, size, empty, neighW, touched, v, gamma)) { improved = true; anyMove = true; }
            }
        }
        return anyMove;
    }

    /// <summary>Leiden fast local moving with a queue of nodes whose neighbourhood changed.</summary>
    private static bool FastLocalMove(Graph g, int[] part, double gamma, Random rng)
    {
        var tot = new double[g.N];
        for (int v = 0; v < g.N; v++) tot[part[v]] += g.Strength[v];
        var (size, empty) = Occupancy(part, g.N);
        var neighW = new double[g.N];
        var touched = new List<int>();
        var queue = new Queue<int>(Shuffled(g.N, rng));
        var inQueue = new bool[g.N];
        Array.Fill(inQueue, true);
        bool anyMove = false;
        while (queue.Count > 0)
        {
            int v = queue.Dequeue();
            inQueue[v] = false;
            if (!MoveNode(g, part, tot, size, empty, neighW, touched, v, gamma)) continue;
            anyMove = true;
            for (int t = g.Ptr[v]; t < g.Ptr[v + 1]; t++)
            {
                int u = g.Adj[t];
                if (!inQueue[u] && part[u] != part[v]) { inQueue[u] = true; queue.Enqueue(u); }
            }
        }
        return anyMove;
    }

    private static (int[] Size, Stack<int> Empty) Occupancy(int[] part, int n)
    {
        var size = new int[n];
        foreach (int c in part) size[c]++;
        var empty = new Stack<int>();
        for (int c = n - 1; c >= 0; c--) if (size[c] == 0) empty.Push(c);
        return (size, empty);
    }

    private static bool MoveNode(Graph g, int[] part, double[] tot, int[] size, Stack<int> empty, double[] neighW, List<int> touched, int v, double gamma)
    {
        int old = part[v];
        double kv = g.Strength[v], m2 = g.TotalWeight;
        touched.Clear();
        for (int t = g.Ptr[v]; t < g.Ptr[v + 1]; t++)
        {
            int c = part[g.Adj[t]];
            if (neighW[c] == 0) touched.Add(c);
            neighW[c] += g.W[t];
        }
        tot[old] -= kv;
        size[old]--;
        double bestGain = neighW[old] - gamma * kv * tot[old] / m2;
        int best = old;
        foreach (int c in touched)
        {
            if (c == old) continue;
            double gain = neighW[c] - gamma * kv * tot[c] / m2;
            if (gain > bestGain + 1e-12) { bestGain = gain; best = c; }
        }
        // An empty community has gain 0.
        if (bestGain < -1e-12 && size[old] > 0 && empty.Count > 0) best = empty.Pop();
        if (size[old] == 0 && best != old) empty.Push(old);
        tot[best] += kv;
        size[best]++;
        part[v] = best;
        foreach (int c in touched) neighW[c] = 0;
        return best != old;
    }

    /// <summary>
    /// Leiden refinement: inside each community, starting from singletons, well-connected nodes merge into
    /// well-connected sub-communities chosen at random with probability ∝ exp(gain / θ), θ = 0.01.
    /// </summary>
    private static int[] Refine(Graph g, int[] part, double gamma, Random rng)
    {
        const double theta = 0.01;
        int n = g.N;
        double m2 = g.TotalWeight;
        var refined = new int[n];
        for (int v = 0; v < n; v++) refined[v] = v;
        var commTot = new double[n];
        for (int v = 0; v < n; v++) commTot[part[v]] += g.Strength[v];
        var refTot = new double[n];
        var refExt = new double[n];   // weight from refined community to the rest of its community
        var singleton = new bool[n];
        for (int v = 0; v < n; v++)
        {
            refTot[v] = g.Strength[v];
            singleton[v] = true;
            for (int t = g.Ptr[v]; t < g.Ptr[v + 1]; t++) if (part[g.Adj[t]] == part[v]) refExt[v] += g.W[t];
        }
        var neighW = new double[n];
        var touched = new List<int>();
        var weights = new List<double>();
        foreach (int v in Shuffled(n, rng))
        {
            if (!singleton[v]) continue;
            int c = part[v];
            double kv = g.Strength[v];
            // Node must be well connected to its community.
            if (refExt[v] < gamma * kv * (commTot[c] - kv) / m2) continue;
            touched.Clear();
            for (int t = g.Ptr[v]; t < g.Ptr[v + 1]; t++)
            {
                int u = g.Adj[t];
                if (part[u] != c) continue;
                int r = refined[u];
                if (neighW[r] == 0) touched.Add(r);
                neighW[r] += g.W[t];
            }
            refTot[refined[v]] -= kv;
            int own = refined[v];
            weights.Clear();
            var cands = new List<int>();
            double maxGain = double.NegativeInfinity;
            foreach (int r in touched)
            {
                if (r == own) continue;
                // Target must be well connected to the rest of its community.
                if (refExt[r] < gamma * refTot[r] * (commTot[c] - refTot[r]) / m2) continue;
                double gain = neighW[r] - gamma * kv * refTot[r] / m2;
                if (gain < 0) continue;
                cands.Add(r); weights.Add(gain);
                maxGain = Math.Max(maxGain, gain);
            }
            int chosen = own;
            if (cands.Count > 0)
            {
                double sum = 0;
                for (int t = 0; t < weights.Count; t++) { weights[t] = Math.Exp((weights[t] - maxGain) / theta); sum += weights[t]; }
                double u = rng.NextDouble() * sum;
                chosen = cands[^1];
                for (int t = 0; t < cands.Count; t++) { u -= weights[t]; if (u <= 0) { chosen = cands[t]; break; } }
            }
            if (chosen != own)
            {
                // Update the external weight of the target: edges between v and chosen become internal,
                // the rest of v's in-community edges become external edges of the merged community.
                double wvComm = refExt[v];
                refExt[chosen] = refExt[chosen] + wvComm - 2 * neighW[chosen];
                refined[v] = chosen;
                singleton[v] = false;
                singleton[chosen] = false;
            }
            refTot[refined[v]] += kv;
            foreach (int r in touched) neighW[r] = 0;
        }
        return refined;
    }

    private static int[] Shuffled(int n, Random rng)
    {
        var a = new int[n];
        for (int i = 0; i < n; i++) a[i] = i;
        rng.Shuffle(a);
        return a;
    }

    private static int Renumber(int[] labels)
    {
        var map = new Dictionary<int, int>();
        for (int i = 0; i < labels.Length; i++)
        {
            if (!map.TryGetValue(labels[i], out int m)) { m = map.Count; map[labels[i]] = m; }
            labels[i] = m;
        }
        return map.Count;
    }

    private static int[] RenumberBySize(int[] labels)
    {
        var size = new Dictionary<int, int>();
        var first = new Dictionary<int, int>();
        for (int i = 0; i < labels.Length; i++) { size[labels[i]] = size.GetValueOrDefault(labels[i]) + 1; first.TryAdd(labels[i], i); }
        var order = size.Keys.OrderByDescending(k => size[k]).ThenBy(k => first[k]).ToArray();
        var map = new Dictionary<int, int>();
        for (int r = 0; r < order.Length; r++) map[order[r]] = r;
        var outp = new int[labels.Length];
        for (int i = 0; i < labels.Length; i++) outp[i] = map[labels[i]];
        return outp;
    }
}
