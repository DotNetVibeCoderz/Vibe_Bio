namespace Biotech.Net.Chemistry;

/// <summary>
/// Ring perception for a molecule: ring atoms/bonds (exact, via bridge detection) and the smallest set of
/// smallest rings (SSSR), computed from Horton's candidate cycles with Gaussian elimination over GF(2).
/// </summary>
/// <remarks>
/// The SSSR has exactly <c>bonds − atoms + components</c> rings. Each ring lists its atoms in cyclic order,
/// and the matching bonds. Ring systems group SSSR rings that share at least one atom (fused, bridged or spiro).
/// </remarks>
public sealed class RingInfo
{
    private readonly int[][] _atomRings;
    private readonly int[][] _bondRings;
    private readonly bool[] _ringBond;
    private readonly int[] _atomRingCount;
    private readonly int[] _smallestRing;
    private readonly int[] _ringBondCount;

    private RingInfo(int[][] atomRings, int[][] bondRings, bool[] ringBond, int atomCount, Molecule mol)
    {
        _atomRings = atomRings;
        _bondRings = bondRings;
        _ringBond = ringBond;
        _atomRingCount = new int[atomCount];
        _smallestRing = new int[atomCount];
        _ringBondCount = new int[atomCount];
        foreach (var ring in atomRings)
            foreach (int a in ring)
            {
                _atomRingCount[a]++;
                if (_smallestRing[a] == 0 || ring.Length < _smallestRing[a]) _smallestRing[a] = ring.Length;
            }
        var bonds = mol.Bonds;
        for (int b = 0; b < bonds.Length; b++)
            if (ringBond[b]) { _ringBondCount[bonds[b].Begin]++; _ringBondCount[bonds[b].End]++; }
        RingSystems = BuildSystems(atomRings);
    }

    /// <summary>SSSR rings as atom indices in ring order.</summary>
    public IReadOnlyList<int[]> AtomRings => _atomRings;

    /// <summary>SSSR rings as bond indices (bond i joins atoms i and i+1 of the matching atom ring).</summary>
    public IReadOnlyList<int[]> BondRings => _bondRings;

    /// <summary>Number of SSSR rings.</summary>
    public int Count => _atomRings.Length;

    /// <summary>Groups of SSSR ring indices sharing atoms (fused/spiro/bridged systems).</summary>
    public IReadOnlyList<int[]> RingSystems { get; }

    /// <summary>Whether the atom belongs to any ring.</summary>
    public bool IsAtomInRing(int atom) => _ringBondCount[atom] > 0;

    /// <summary>Whether the bond belongs to any ring (is not a bridge of the graph).</summary>
    public bool IsBondInRing(int bond) => _ringBond[bond];

    /// <summary>Number of SSSR rings containing the atom (SMARTS <c>R&lt;n&gt;</c>).</summary>
    public int AtomRingCount(int atom) => _atomRingCount[atom];

    /// <summary>Number of ring bonds at the atom (SMARTS <c>x&lt;n&gt;</c>).</summary>
    public int RingBondCount(int atom) => _ringBondCount[atom];

    /// <summary>Size of the smallest SSSR ring containing the atom, or 0.</summary>
    public int SmallestRingSize(int atom) => _smallestRing[atom];

    /// <summary>Whether the atom is in an SSSR ring of exactly <paramref name="size"/> atoms (SMARTS <c>r&lt;n&gt;</c>).</summary>
    public bool IsAtomInRingOfSize(int atom, int size)
    {
        foreach (var r in _atomRings)
            if (r.Length == size && Array.IndexOf(r, atom) >= 0) return true;
        return false;
    }

    /// <summary>Whether the bond is in an SSSR ring of exactly <paramref name="size"/> atoms.</summary>
    public bool IsBondInRingOfSize(int bond, int size)
    {
        foreach (var r in _bondRings)
            if (r.Length == size && Array.IndexOf(r, bond) >= 0) return true;
        return false;
    }

    internal static RingInfo Compute(Molecule mol)
    {
        int n = mol.AtomCount, m = mol.BondCount;
        var ringBond = FindRingBonds(mol);
        int components = mol.Fragments().Count;
        int expected = m - n + components;
        if (expected <= 0) return new RingInfo([], [], ringBond, n, mol);

        // Horton candidates restricted to the ring subgraph.
        var candidates = new List<int[]>(); // atom cycles
        var seen = new HashSet<string>();
        var dist = new int[n];
        var parent = new int[n];
        var queue = new int[n];
        var onPath = new int[n];
        int stamp = 0;
        var bonds = mol.Bonds;
        for (int r = 0; r < n; r++)
        {
            if (!IsRingAtom(mol, r, ringBond)) continue;
            Array.Fill(dist, -1);
            dist[r] = 0; parent[r] = -1;
            int head = 0, tail = 0;
            queue[tail++] = r;
            while (head < tail)
            {
                int u = queue[head++];
                var nb = mol.Neighbors(u);
                var nbb = mol.NeighborBonds(u);
                for (int k = 0; k < nb.Length; k++)
                {
                    if (!ringBond[nbb[k]]) continue;
                    int v = nb[k];
                    if (dist[v] >= 0) continue;
                    dist[v] = dist[u] + 1;
                    parent[v] = u;
                    queue[tail++] = v;
                }
            }
            for (int b = 0; b < m; b++)
            {
                if (!ringBond[b]) continue;
                int x = bonds[b].Begin, y = bonds[b].End;
                if (dist[x] < 0 || dist[y] < 0) continue;
                if (parent[x] == y || parent[y] == x) continue; // tree edge
                // paths r..x and r..y must only share r
                stamp++;
                for (int p = x; p != -1; p = parent[p]) onPath[p] = stamp;
                bool disjoint = true;
                for (int p = y; p != r; p = parent[p]) if (onPath[p] == stamp) { disjoint = false; break; }
                if (!disjoint) continue;
                var cycle = new List<int>(dist[x] + dist[y] + 1);
                var px = new List<int>();
                for (int p = x; p != -1; p = parent[p]) px.Add(p);
                px.Reverse(); // r .. x
                cycle.AddRange(px);
                for (int p = y; p != r; p = parent[p]) cycle.Add(p); // y .. child of r
                var arr = cycle.ToArray();
                var key = string.Join(",", arr.OrderBy(a => a));
                if (seen.Add(key)) candidates.Add(arr);
            }
        }
        candidates.Sort((a, b) => a.Length.CompareTo(b.Length));

        // Gaussian elimination over GF(2) on bond incidence vectors.
        int words = (m + 63) / 64;
        var basis = new List<ulong[]>();
        var pivots = new List<int>();
        var atomRings = new List<int[]>();
        var bondRings = new List<int[]>();
        foreach (var cyc in candidates)
        {
            if (atomRings.Count == expected) break;
            var bondList = new int[cyc.Length];
            var vec = new ulong[words];
            for (int i = 0; i < cyc.Length; i++)
            {
                int bi = mol.BondBetween(cyc[i], cyc[(i + 1) % cyc.Length]);
                bondList[i] = bi;
                vec[bi >> 6] |= 1UL << (bi & 63);
            }
            var reduced = (ulong[])vec.Clone();
            for (int k = 0; k < basis.Count; k++)
            {
                int p = pivots[k];
                if ((reduced[p >> 6] & (1UL << (p & 63))) != 0)
                    for (int w = 0; w < words; w++) reduced[w] ^= basis[k][w];
            }
            int pivot = -1;
            for (int w = 0; w < words && pivot < 0; w++)
                if (reduced[w] != 0) pivot = w * 64 + System.Numerics.BitOperations.TrailingZeroCount(reduced[w]);
            if (pivot < 0) continue;
            basis.Add(reduced);
            pivots.Add(pivot);
            atomRings.Add(cyc);
            bondRings.Add(bondList);
        }
        return new RingInfo(atomRings.ToArray(), bondRings.ToArray(), ringBond, n, mol);
    }

    private static bool IsRingAtom(Molecule mol, int atom, bool[] ringBond)
    {
        foreach (int b in mol.NeighborBonds(atom)) if (ringBond[b]) return true;
        return false;
    }

    /// <summary>Marks bonds that are not bridges (iterative Tarjan low-link).</summary>
    private static bool[] FindRingBonds(Molecule mol)
    {
        int n = mol.AtomCount;
        var ring = new bool[mol.BondCount];
        Array.Fill(ring, true);
        var disc = new int[n];
        var low = new int[n];
        Array.Fill(disc, -1);
        int time = 0;
        var stack = new Stack<(int Atom, int ParentBond, int Next)>();
        for (int s = 0; s < n; s++)
        {
            if (disc[s] >= 0) continue;
            stack.Push((s, -1, 0));
            disc[s] = low[s] = time++;
            while (stack.Count > 0)
            {
                var (u, pb, next) = stack.Pop();
                var nb = mol.Neighbors(u);
                var nbb = mol.NeighborBonds(u);
                if (next < nb.Length)
                {
                    stack.Push((u, pb, next + 1));
                    int v = nb[next], b = nbb[next];
                    if (b == pb) continue;
                    if (disc[v] < 0)
                    {
                        disc[v] = low[v] = time++;
                        stack.Push((v, b, 0));
                    }
                    else low[u] = Math.Min(low[u], disc[v]);
                }
                else if (pb >= 0)
                {
                    int p = mol.Bonds[pb].Other(u);
                    low[p] = Math.Min(low[p], low[u]);
                    if (low[u] > disc[p]) ring[pb] = false;
                }
            }
        }
        return ring;
    }

    private static int[][] BuildSystems(int[][] rings)
    {
        int k = rings.Length;
        var parent = Enumerable.Range(0, k).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        for (int i = 0; i < k; i++)
            for (int j = i + 1; j < k; j++)
                if (rings[i].Intersect(rings[j]).Any()) parent[Find(i)] = Find(j);
        return Enumerable.Range(0, k).GroupBy(Find).Select(g => g.ToArray()).ToArray();
    }
}
