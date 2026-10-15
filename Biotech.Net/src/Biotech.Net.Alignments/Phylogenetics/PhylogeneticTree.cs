using System.Globalization;
using System.Text;

namespace Biotech.Net.Alignments.Phylogenetics;

/// <summary>A node of a <see cref="PhylogeneticTree"/>: a leaf (taxon) or an internal node (clade).</summary>
public sealed class TreeNode
{
    private readonly List<TreeNode> _children = [];

    /// <summary>Creates a node.</summary>
    public TreeNode(string? name = null, double? branchLength = null)
    {
        Name = name;
        BranchLength = branchLength;
    }

    /// <summary>Taxon name for leaves; optional label (often a support value) for internal nodes.</summary>
    public string? Name { get; set; }

    /// <summary>Length of the branch to the parent; null when unknown.</summary>
    public double? BranchLength { get; set; }

    /// <summary>The parent, or null for the root.</summary>
    public TreeNode? Parent { get; private set; }

    /// <summary>Child nodes in order.</summary>
    public IReadOnlyList<TreeNode> Children => _children;

    /// <summary>True when the node has no children.</summary>
    public bool IsLeaf => _children.Count == 0;

    /// <summary>Appends <paramref name="child"/> (detaching it from any previous parent) and returns it.</summary>
    public TreeNode AddChild(TreeNode child)
    {
        ArgumentNullException.ThrowIfNull(child);
        child.Parent?._children.Remove(child);
        child.Parent = this;
        _children.Add(child);
        return child;
    }

    /// <summary>Removes <paramref name="child"/> from this node.</summary>
    public bool RemoveChild(TreeNode child)
    {
        if (!_children.Remove(child)) return false;
        child.Parent = null;
        return true;
    }

    /// <summary>This node and all descendants in pre-order.</summary>
    public IEnumerable<TreeNode> PreOrder()
    {
        var stack = new Stack<TreeNode>();
        stack.Push(this);
        while (stack.Count > 0)
        {
            var n = stack.Pop();
            yield return n;
            for (int i = n._children.Count - 1; i >= 0; i--) stack.Push(n._children[i]);
        }
    }

    /// <summary>This node and all descendants in post-order (children before parents).</summary>
    public IEnumerable<TreeNode> PostOrder()
    {
        var stack = new Stack<(TreeNode Node, int Next)>();
        stack.Push((this, 0));
        while (stack.Count > 0)
        {
            var (n, next) = stack.Pop();
            if (next < n._children.Count)
            {
                stack.Push((n, next + 1));
                stack.Push((n._children[next], 0));
            }
            else yield return n;
        }
    }

    /// <summary>Leaves under this node, left to right.</summary>
    public IEnumerable<TreeNode> Leaves() => PreOrder().Where(n => n.IsLeaf);

    /// <inheritdoc />
    public override string ToString() => Name ?? (IsLeaf ? "(leaf)" : $"(clade of {Leaves().Count()})");
}

/// <summary>
/// A rooted phylogenetic tree with Newick input/output, midpoint rooting, patristic distances and the
/// Robinson–Foulds distance. Trees from <see cref="TreeBuilder.NeighborJoining"/> are unrooted in meaning; their root
/// is the final trifurcation.
/// </summary>
public sealed class PhylogeneticTree
{
    /// <summary>Wraps a root node.</summary>
    public PhylogeneticTree(TreeNode root)
    {
        ArgumentNullException.ThrowIfNull(root);
        Root = root;
    }

    /// <summary>The root node.</summary>
    public TreeNode Root { get; private set; }

    /// <summary>Leaves, left to right.</summary>
    public IReadOnlyList<TreeNode> Leaves => Root.Leaves().ToList();

    /// <summary>Leaf names, left to right.</summary>
    public IReadOnlyList<string> LeafNames => Root.Leaves().Select(l => l.Name ?? string.Empty).ToList();

    /// <summary>Sum of all branch lengths (unknown lengths count as 0).</summary>
    public double TotalBranchLength => Root.PreOrder().Where(n => n != Root).Sum(n => n.BranchLength ?? 0);

    /// <summary>The leaf named <paramref name="name"/>, or null.</summary>
    public TreeNode? FindLeaf(string name) => Root.Leaves().FirstOrDefault(l => l.Name == name);

    /// <summary>Path length between two nodes (sum of branch lengths along the path).</summary>
    public static double Distance(TreeNode a, TreeNode b)
    {
        var depth = new Dictionary<TreeNode, double>();
        double d = 0;
        for (var n = a; n is not null; n = n.Parent) { depth[n] = d; d += n.BranchLength ?? 0; }
        d = 0;
        for (var n = b; n is not null; n = n.Parent)
        {
            if (depth.TryGetValue(n, out double da)) return da + d;
            d += n.BranchLength ?? 0;
        }
        throw new ArgumentException("The nodes are not in the same tree.");
    }

    /// <summary>Path length between two named leaves.</summary>
    public double Distance(string leafA, string leafB)
        => Distance(FindLeaf(leafA) ?? throw new KeyNotFoundException(leafA), FindLeaf(leafB) ?? throw new KeyNotFoundException(leafB));

    // ---- Newick -----------------------------------------------------------------------------------

    /// <summary>Parses a Newick string (quoted labels, <c>[comments]</c>, internal labels and branch lengths supported).</summary>
    public static PhylogeneticTree ParseNewick(string newick)
    {
        ArgumentNullException.ThrowIfNull(newick);
        var p = new NewickParser(newick);
        var root = p.ParseSubtree();
        p.SkipSpace();
        if (p.Pos < newick.Length && newick[p.Pos] == ';') p.Pos++;
        p.SkipSpace();
        if (p.Pos != newick.Length) throw new FormatException($"Unexpected text after the tree at position {p.Pos}.");
        return new PhylogeneticTree(root);
    }

    /// <summary>
    /// Writes Newick, e.g. <c>((A:0.1,B:0.2):0.05,C:0.3);</c>. Lengths use the shortest round-trip form unless
    /// <paramref name="decimals"/> is given.
    /// </summary>
    public string ToNewick(bool branchLengths = true, int? decimals = null)
    {
        var sb = new StringBuilder();
        Write(Root, sb, branchLengths, decimals);
        return sb.Append(';').ToString();
    }

    private static void Write(TreeNode n, StringBuilder sb, bool lengths, int? decimals)
    {
        if (!n.IsLeaf)
        {
            sb.Append('(');
            for (int i = 0; i < n.Children.Count; i++)
            {
                if (i > 0) sb.Append(',');
                Write(n.Children[i], sb, lengths, decimals);
            }
            sb.Append(')');
        }
        if (n.Name is not null) sb.Append(QuoteLabel(n.Name));
        if (lengths && n.BranchLength is double bl)
            sb.Append(':').Append(decimals is int d ? Math.Round(bl, d).ToString("F" + d, CultureInfo.InvariantCulture) : bl.ToString("R", CultureInfo.InvariantCulture));
    }

    private static string QuoteLabel(string name)
    {
        if (name.Length > 0 && name.IndexOfAny([' ', '(', ')', '[', ']', '\'', ':', ';', ',', '\t', '_']) < 0) return name;
        return "'" + name.Replace("'", "''") + "'";
    }

    private sealed class NewickParser(string s)
    {
        public int Pos;

        public void SkipSpace()
        {
            while (Pos < s.Length)
            {
                if (char.IsWhiteSpace(s[Pos])) Pos++;
                else if (s[Pos] == '[')
                {
                    int end = s.IndexOf(']', Pos);
                    if (end < 0) throw new FormatException("Unterminated [comment].");
                    Pos = end + 1;
                }
                else break;
            }
        }

        public TreeNode ParseSubtree()
        {
            SkipSpace();
            var node = new TreeNode();
            if (Pos < s.Length && s[Pos] == '(')
            {
                Pos++;
                while (true)
                {
                    node.AddChild(ParseSubtree());
                    SkipSpace();
                    if (Pos >= s.Length) throw new FormatException("Unbalanced parentheses.");
                    if (s[Pos] == ',') { Pos++; continue; }
                    if (s[Pos] == ')') { Pos++; break; }
                    throw new FormatException($"Expected ',' or ')' at position {Pos}.");
                }
            }
            SkipSpace();
            node.Name = ParseLabel();
            SkipSpace();
            if (Pos < s.Length && s[Pos] == ':')
            {
                Pos++;
                SkipSpace();
                int start = Pos;
                while (Pos < s.Length && (char.IsDigit(s[Pos]) || s[Pos] is '.' or '-' or '+' or 'e' or 'E')) Pos++;
                if (!double.TryParse(s.AsSpan(start, Pos - start), NumberStyles.Float, CultureInfo.InvariantCulture, out double bl))
                    throw new FormatException($"Invalid branch length at position {start}.");
                node.BranchLength = bl;
            }
            return node;
        }

        private string? ParseLabel()
        {
            if (Pos >= s.Length) return null;
            if (s[Pos] == '\'')
            {
                var sb = new StringBuilder();
                Pos++;
                while (true)
                {
                    if (Pos >= s.Length) throw new FormatException("Unterminated quoted label.");
                    if (s[Pos] == '\'')
                    {
                        if (Pos + 1 < s.Length && s[Pos + 1] == '\'') { sb.Append('\''); Pos += 2; continue; }
                        Pos++;
                        break;
                    }
                    sb.Append(s[Pos++]);
                }
                return sb.ToString();
            }
            int start = Pos;
            while (Pos < s.Length && s[Pos] is not ('(' or ')' or ',' or ':' or ';' or '[') && !char.IsWhiteSpace(s[Pos])) Pos++;
            // Unquoted underscores stand for blanks in the Newick standard.
            return Pos > start ? s[start..Pos].Replace('_', ' ') : null;
        }
    }

    // ---- Rooting ----------------------------------------------------------------------------------

    /// <summary>
    /// Re-roots the tree at the midpoint of the longest leaf-to-leaf path (in place) and returns it. A former
    /// bifurcating root left with a single child is removed and its two branches merged.
    /// </summary>
    public PhylogeneticTree RootAtMidpoint()
    {
        var leaves = Root.Leaves().ToList();
        if (leaves.Count < 2) return this;
        var adj = BuildAdjacency();
        var (a, _) = Farthest(adj, leaves[0]);
        var (b, dist) = Farthest(adj, a);
        double half = dist / 2;

        // Walk from a to b and find the edge containing the midpoint.
        var path = PathBetween(adj, a, b);
        double acc = 0;
        for (int k = 0; k + 1 < path.Count; k++)
        {
            var (u, v) = (path[k], path[k + 1]);
            double len = adj[u].First(e => e.Node == v).Length;
            if (acc + len >= half)
            {
                double fromU = half - acc;
                var newRoot = new TreeNode();
                // Split the edge u–v by inserting newRoot.
                adj[u].RemoveAll(e => e.Node == v);
                adj[v].RemoveAll(e => e.Node == u);
                adj[newRoot] = [(u, fromU), (v, len - fromU)];
                adj[u].Add((newRoot, fromU));
                adj[v].Add((newRoot, len - fromU));
                Rebuild(adj, newRoot);
                return this;
            }
            acc += len;
        }
        return this;
    }

    private Dictionary<TreeNode, List<(TreeNode Node, double Length)>> BuildAdjacency()
    {
        var adj = new Dictionary<TreeNode, List<(TreeNode, double)>>();
        foreach (var n in Root.PreOrder()) adj[n] = [];
        foreach (var n in Root.PreOrder())
            if (n.Parent is not null)
            {
                double l = n.BranchLength ?? 0;
                adj[n].Add((n.Parent, l));
                adj[n.Parent].Add((n, l));
            }
        // A bifurcating root is just a point on an edge: remove it so it does not survive re-rooting.
        if (Root.Children.Count == 2)
        {
            var (c1, l1) = adj[Root][0];
            var (c2, l2) = adj[Root][1];
            adj[c1].RemoveAll(e => e.Item1 == Root);
            adj[c2].RemoveAll(e => e.Item1 == Root);
            adj[c1].Add((c2, l1 + l2));
            adj[c2].Add((c1, l1 + l2));
            adj.Remove(Root);
        }
        return adj;
    }

    private static (TreeNode Node, double Distance) Farthest(Dictionary<TreeNode, List<(TreeNode Node, double Length)>> adj, TreeNode from)
    {
        var dist = new Dictionary<TreeNode, double> { [from] = 0 };
        var stack = new Stack<TreeNode>();
        stack.Push(from);
        TreeNode best = from;
        while (stack.Count > 0)
        {
            var u = stack.Pop();
            foreach (var (v, l) in adj[u])
            {
                if (dist.ContainsKey(v)) continue;
                dist[v] = dist[u] + l;
                if (v.IsLeaf && dist[v] > dist[best]) best = v;
                stack.Push(v);
            }
        }
        return (best, dist[best]);
    }

    private static List<TreeNode> PathBetween(Dictionary<TreeNode, List<(TreeNode Node, double Length)>> adj, TreeNode a, TreeNode b)
    {
        var prev = new Dictionary<TreeNode, TreeNode?> { [a] = null };
        var stack = new Stack<TreeNode>();
        stack.Push(a);
        while (stack.Count > 0)
        {
            var u = stack.Pop();
            if (u == b) break;
            foreach (var (v, _) in adj[u])
                if (prev.TryAdd(v, u)) stack.Push(v);
        }
        var path = new List<TreeNode>();
        for (TreeNode? n = b; n is not null; n = prev[n]) path.Add(n);
        path.Reverse();
        return path;
    }

    private void Rebuild(Dictionary<TreeNode, List<(TreeNode Node, double Length)>> adj, TreeNode newRoot)
    {
        // Detach everything, then re-link following the adjacency outward from the new root.
        foreach (var n in adj.Keys.Append(Root))
            foreach (var c in n.Children.ToList()) n.RemoveChild(c);
        newRoot.BranchLength = null;
        var seen = new HashSet<TreeNode> { newRoot };
        var queue = new Queue<TreeNode>();
        queue.Enqueue(newRoot);
        while (queue.Count > 0)
        {
            var u = queue.Dequeue();
            foreach (var (v, l) in adj[u])
            {
                if (!seen.Add(v)) continue;
                v.BranchLength = l;
                u.AddChild(v);
                queue.Enqueue(v);
            }
        }
        Root = newRoot;
    }

    // ---- Comparison -------------------------------------------------------------------------------

    /// <summary>
    /// Robinson–Foulds distance: the number of non-trivial bipartitions (splits) present in exactly one of the two trees,
    /// both treated as unrooted. With <paramref name="normalized"/> it is divided by the maximum, 2(n − 3).
    /// </summary>
    /// <exception cref="ArgumentException">The trees do not have the same leaf names.</exception>
    public double RobinsonFoulds(PhylogeneticTree other, bool normalized = false)
    {
        ArgumentNullException.ThrowIfNull(other);
        var names = LeafNames.OrderBy(x => x, StringComparer.Ordinal).ToList();
        var otherNames = other.LeafNames.OrderBy(x => x, StringComparer.Ordinal).ToList();
        if (!names.SequenceEqual(otherNames) || names.Distinct().Count() != names.Count)
            throw new ArgumentException("Robinson–Foulds needs two trees with the same, unique leaf names.");
        var index = names.Select((n, i) => (n, i)).ToDictionary(x => x.n, x => x.i);
        var a = Splits(index);
        var b = other.Splits(index);
        int rf = a.Count(s => !b.Contains(s)) + b.Count(s => !a.Contains(s));
        int n = names.Count;
        if (!normalized) return rf;
        return n <= 3 ? 0 : rf / (2.0 * (n - 3));
    }

    /// <summary>Non-trivial splits, each canonicalised as the side that excludes leaf 0, encoded as a sorted index string.</summary>
    internal HashSet<string> Splits(Dictionary<string, int> index)
    {
        int n = index.Count;
        var below = new Dictionary<TreeNode, List<int>>();
        var result = new HashSet<string>();
        foreach (var node in Root.PostOrder())
        {
            List<int> set = node.IsLeaf ? [index[node.Name ?? string.Empty]] : node.Children.SelectMany(c => below[c]).ToList();
            below[node] = set;
            if (node == Root) continue;
            var side = set.Contains(0) ? Enumerable.Range(0, n).Except(set).ToList() : [.. set];
            if (side.Count < 2 || side.Count > n - 2) continue;
            side.Sort();
            result.Add(string.Join(',', side));
        }
        return result;
    }

    /// <inheritdoc />
    public override string ToString() => ToNewick();
}
