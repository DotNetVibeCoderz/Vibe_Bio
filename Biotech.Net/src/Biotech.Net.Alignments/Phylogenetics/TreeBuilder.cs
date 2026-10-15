namespace Biotech.Net.Alignments.Phylogenetics;

/// <summary>Distance-based tree construction: UPGMA and Neighbor-Joining.</summary>
public static class TreeBuilder
{
    /// <summary>
    /// UPGMA (Sokal &amp; Michener 1958): repeatedly joins the closest pair of clusters and averages distances weighted by
    /// cluster size. The result is rooted and ultrametric; branch lengths are differences of node heights (half the
    /// joining distance). Ties are resolved by the lowest row indices.
    /// </summary>
    public static PhylogeneticTree Upgma(DistanceMatrix distances)
    {
        ArgumentNullException.ThrowIfNull(distances);
        int n = distances.Count;
        if (n == 0) throw new ArgumentException("The distance matrix is empty.");
        var nodes = new List<TreeNode>();
        var heights = new List<double>();
        var sizes = new List<int>();
        var d = new List<List<double>>();
        for (int i = 0; i < n; i++)
        {
            nodes.Add(new TreeNode(distances.Names[i]));
            heights.Add(0);
            sizes.Add(1);
            d.Add(Enumerable.Range(0, n).Select(j => distances[i, j]).ToList());
        }
        while (nodes.Count > 1)
        {
            int bi = 0, bj = 1;
            double best = double.PositiveInfinity;
            for (int i = 0; i < nodes.Count; i++)
                for (int j = i + 1; j < nodes.Count; j++)
                    if (d[i][j] < best) { best = d[i][j]; bi = i; bj = j; }
            if (double.IsPositiveInfinity(best)) best = d[bi][bj];
            double h = best / 2;
            var u = new TreeNode();
            nodes[bi].BranchLength = h - heights[bi];
            nodes[bj].BranchLength = h - heights[bj];
            u.AddChild(nodes[bi]);
            u.AddChild(nodes[bj]);
            int si = sizes[bi], sj = sizes[bj];
            var row = new List<double>();
            for (int k = 0; k < nodes.Count; k++)
                row.Add(k == bi || k == bj ? 0 : (d[bi][k] * si + d[bj][k] * sj) / (si + sj));
            // Replace bi with u, remove bj.
            nodes[bi] = u; heights[bi] = h; sizes[bi] = si + sj;
            for (int k = 0; k < nodes.Count; k++) { d[bi][k] = row[k]; d[k][bi] = row[k]; }
            d[bi][bi] = 0;
            Remove(bj, nodes, heights, sizes, d);
        }
        return new PhylogeneticTree(nodes[0]);
    }

    /// <summary>
    /// Neighbor-Joining (Saitou &amp; Nei 1987; Studier &amp; Keppler 1988 formulation): repeatedly joins the pair minimising
    /// <c>Q(i,j) = (r − 2)·d(i,j) − Σₖd(i,k) − Σₖd(j,k)</c>. The returned tree is unrooted in meaning, with the last three
    /// nodes joined at a trifurcating root; use <see cref="PhylogeneticTree.RootAtMidpoint"/> to root it. Negative branch
    /// lengths are kept as computed unless <paramref name="clampNegative"/> is set.
    /// </summary>
    public static PhylogeneticTree NeighborJoining(DistanceMatrix distances, bool clampNegative = false)
    {
        ArgumentNullException.ThrowIfNull(distances);
        int n = distances.Count;
        if (n == 0) throw new ArgumentException("The distance matrix is empty.");
        var nodes = Enumerable.Range(0, n).Select(i => new TreeNode(distances.Names[i])).ToList();
        if (n == 1) return new PhylogeneticTree(nodes[0]);
        var d = Enumerable.Range(0, n).Select(i => Enumerable.Range(0, n).Select(j => distances[i, j]).ToList()).ToList();
        double Fix(double x) => clampNegative && x < 0 ? 0 : x;
        if (n == 2)
        {
            var root = new TreeNode();
            nodes[0].BranchLength = Fix(d[0][1] / 2); nodes[1].BranchLength = Fix(d[0][1] / 2);
            root.AddChild(nodes[0]); root.AddChild(nodes[1]);
            return new PhylogeneticTree(root);
        }
        while (nodes.Count > 3)
        {
            int r = nodes.Count;
            var sum = new double[r];
            for (int i = 0; i < r; i++) for (int k = 0; k < r; k++) sum[i] += d[i][k];
            int bi = 0, bj = 1;
            double best = double.PositiveInfinity;
            for (int i = 0; i < r; i++)
                for (int j = i + 1; j < r; j++)
                {
                    double q = (r - 2) * d[i][j] - sum[i] - sum[j];
                    if (q < best) { best = q; bi = i; bj = j; }
                }
            double dij = d[bi][bj];
            double li = dij / 2 + (sum[bi] - sum[bj]) / (2.0 * (r - 2));
            double lj = dij - li;
            var u = new TreeNode();
            nodes[bi].BranchLength = Fix(li);
            nodes[bj].BranchLength = Fix(lj);
            u.AddChild(nodes[bi]);
            u.AddChild(nodes[bj]);
            var row = new double[r];
            for (int k = 0; k < r; k++) row[k] = k == bi || k == bj ? 0 : (d[bi][k] + d[bj][k] - dij) / 2;
            nodes[bi] = u;
            for (int k = 0; k < r; k++) { d[bi][k] = row[k]; d[k][bi] = row[k]; }
            d[bi][bi] = 0;
            var dummyH = new List<double>(new double[r]);
            var dummyS = new List<int>(new int[r]);
            Remove(bj, nodes, dummyH, dummyS, d);
        }
        // Three nodes left: join them at a central (root) node.
        var center = new TreeNode();
        double d01 = d[0][1], d02 = d[0][2], d12 = d[1][2];
        nodes[0].BranchLength = Fix((d01 + d02 - d12) / 2);
        nodes[1].BranchLength = Fix((d01 + d12 - d02) / 2);
        nodes[2].BranchLength = Fix((d02 + d12 - d01) / 2);
        foreach (var node in nodes) center.AddChild(node);
        return new PhylogeneticTree(center);
    }

    private static void Remove(int j, List<TreeNode> nodes, List<double> heights, List<int> sizes, List<List<double>> d)
    {
        nodes.RemoveAt(j); heights.RemoveAt(j); sizes.RemoveAt(j);
        d.RemoveAt(j);
        foreach (var row in d) row.RemoveAt(j);
    }
}
