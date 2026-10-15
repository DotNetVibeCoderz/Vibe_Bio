using System.Text;

namespace Biotech.Net.Chemistry;

/// <summary>
/// Canonical atom ranking: atom invariants (degree, element, isotope, charge, H count, aromaticity, ring
/// membership) refined iteratively with sorted neighbour ranks and bond orders (extended-connectivity /
/// Morgan-style), with deterministic tie-breaking until every atom has a distinct rank.
/// </summary>
public static class CanonicalRanking
{
    /// <summary>Canonical ranks 0…n−1; symmetry-equivalent atoms are separated by tie-breaking.</summary>
    public static int[] Compute(Molecule mol) => Compute(mol, breakTies: true);

    /// <summary>Symmetry classes (refined ranks without tie-breaking): equal values mean topologically equivalent atoms.</summary>
    public static int[] SymmetryClasses(Molecule mol) => Compute(mol, breakTies: false);

    private static int[] Compute(Molecule mol, bool breakTies)
    {
        int n = mol.AtomCount;
        if (n == 0) return [];
        var atoms = mol.Atoms;
        var rings = mol.Rings;
        var inv = new long[n];
        for (int i = 0; i < n; i++)
        {
            var a = atoms[i];
            long v = mol.Degree(i);
            v = v * 128 + a.AtomicNumber;
            v = v * 512 + Math.Min(a.Isotope, 511);
            v = v * 16 + (a.FormalCharge + 8);
            v = v * 8 + Math.Min(mol.TotalHydrogens(i), 7);
            v = v * 2 + (a.IsAromatic ? 1 : 0);
            v = v * 8 + Math.Min(rings.AtomRingCount(i), 7);
            v = v * 64 + Math.Min(a.MapNumber, 63);
            inv[i] = v;
        }
        var rank = DenseRank(inv);
        rank = Refine(mol, rank);
        if (!breakTies) return rank;
        while (true)
        {
            int classes = rank.Max() + 1;
            if (classes == n) break;
            // smallest tied rank
            var counts = new int[classes];
            foreach (int r in rank) counts[r]++;
            int tied = Array.FindIndex(counts, c => c > 1);
            int chosen = Array.IndexOf(rank, tied);
            var next = new long[n];
            for (int i = 0; i < n; i++) next[i] = rank[i] * 2L + (i == chosen || rank[i] != tied ? 0 : 1);
            rank = Refine(mol, DenseRank(next));
        }
        return rank;
    }

    private static int[] Refine(Molecule mol, int[] rank)
    {
        int n = rank.Length;
        int classes = rank.Max() + 1;
        var bonds = mol.Bonds;
        while (true)
        {
            var keys = new long[n][];
            for (int i = 0; i < n; i++)
            {
                var nb = mol.Neighbors(i);
                var nbb = mol.NeighborBonds(i);
                var k = new long[nb.Length + 1];
                k[0] = rank[i];
                for (int j = 0; j < nb.Length; j++) k[j + 1] = rank[nb[j]] * 8L + (int)bonds[nbb[j]].Order;
                Array.Sort(k, 1, nb.Length);
                keys[i] = k;
            }
            var idx = Enumerable.Range(0, n).ToArray();
            Array.Sort(idx, (a, b) => CompareKeys(keys[a], keys[b]));
            var nr = new int[n];
            int r = 0;
            for (int i = 0; i < n; i++)
            {
                if (i > 0 && CompareKeys(keys[idx[i]], keys[idx[i - 1]]) != 0) r++;
                nr[idx[i]] = r;
            }
            int newClasses = r + 1;
            rank = nr;
            if (newClasses == classes) return rank;
            classes = newClasses;
        }
    }

    private static int CompareKeys(long[] a, long[] b)
    {
        int len = Math.Min(a.Length, b.Length);
        for (int i = 0; i < len; i++) { int c = a[i].CompareTo(b[i]); if (c != 0) return c; }
        return a.Length.CompareTo(b.Length);
    }

    private static int[] DenseRank(long[] values)
    {
        var distinct = values.Distinct().OrderBy(v => v).ToArray();
        var map = new Dictionary<long, int>(distinct.Length);
        for (int i = 0; i < distinct.Length; i++) map[distinct[i]] = i;
        return values.Select(v => map[v]).ToArray();
    }
}

/// <summary>
/// Writes SMILES. In canonical mode atoms are ranked with <see cref="CanonicalRanking"/> and the depth-first
/// traversal visits neighbours in rank order, so equal molecules give identical strings. Tetrahedral
/// stereo is re-expressed for the output neighbour order (permutation parity); cis/trans double bonds are
/// written with <c>/</c> and <c>\</c>.
/// </summary>
public static class SmilesWriter
{
    /// <summary>SMILES for <paramref name="mol"/>.</summary>
    public static string Write(Molecule mol, bool canonical = true, bool includeStereo = true)
    {
        int n = mol.AtomCount;
        if (n == 0) return string.Empty;
        var rank = canonical ? CanonicalRanking.Compute(mol) : Enumerable.Range(0, n).ToArray();
        return new Writer(mol, rank, includeStereo).Run();
    }

    /// <summary>
    /// The implicit hydrogen count a SMILES reader would infer for an organic-subset atom written without
    /// brackets, or −1 if the atom cannot be written that way.
    /// </summary>
    internal static int OrganicSubsetHydrogens(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        if (!PeriodicTable.Get(a.AtomicNumber).IsOrganicSubset && a.AtomicNumber != 0) return -1;
        if (a.AtomicNumber == 0) return 0;
        int sigma = 0;
        foreach (int b in mol.NeighborBonds(atom))
        {
            var bond = mol.Bonds[b];
            sigma += bond.Order == BondOrder.Aromatic ? 1 : (int)bond.Order;
        }
        var allowed = Molecule.AllowedValences(a.AtomicNumber, 0);
        int t = -1;
        foreach (int v in allowed) if (v >= sigma) { t = v; break; }
        if (t < 0) return a.IsAromatic ? -1 : 0;
        if (a.IsAromatic) return t - sigma >= 1 ? t - sigma - 1 : 0;
        return t - sigma;
    }

    private sealed class Writer(Molecule mol, int[] rank, bool stereo)
    {
        private readonly StringBuilder _sb = new();
        private readonly int[] _pos = Enumerable.Repeat(-1, mol.AtomCount).ToArray();
        private readonly int[] _parent = new int[mol.AtomCount];
        private readonly int[] _parentBond = new int[mol.AtomCount];
        private readonly List<int>[] _children = Enumerable.Range(0, mol.AtomCount).Select(_ => new List<int>()).ToArray();
        private readonly List<int>[] _closures = Enumerable.Range(0, mol.AtomCount).Select(_ => new List<int>()).ToArray();
        private readonly HashSet<int> _closureBonds = [];
        private readonly Dictionary<int, char> _dir = [];
        private readonly Dictionary<int, int> _digitOf = [];
        private readonly SortedSet<int> _freeDigits = new(Enumerable.Range(1, 99));
        private int _counter;

        public string Run()
        {
            int n = mol.AtomCount;
            var fragments = mol.Fragments()
                .Select(f => f.OrderBy(a => rank[a]).First())
                .OrderBy(a => rank[a]).ToList();
            foreach (int start in fragments)
            {
                _parent[start] = -1;
                _parentBond[start] = -1;
                Dfs(start);
            }
            if (stereo) AssignDoubleBondDirections();
            bool first = true;
            foreach (int start in fragments)
            {
                if (!first) _sb.Append('.');
                first = false;
                Emit(start);
            }
            _ = n;
            return _sb.ToString();
        }

        private void Dfs(int u)
        {
            _pos[u] = _counter++;
            var nb = mol.Neighbors(u).ToArray();
            var nbb = mol.NeighborBonds(u).ToArray();
            var order = Enumerable.Range(0, nb.Length).OrderBy(k => rank[nb[k]]).ToArray();
            foreach (int k in order)
            {
                int v = nb[k], b = nbb[k];
                if (b == _parentBond[u]) continue;
                if (_pos[v] >= 0)
                {
                    if (_closureBonds.Add(b)) { _closures[u].Add(b); _closures[v].Add(b); }
                    continue;
                }
                _parent[v] = u;
                _parentBond[v] = b;
                _children[u].Add(v);
                Dfs(v);
            }
        }

        private void Emit(int u)
        {
            // closures at u: closing ones (partner already emitted → has a digit) first, then openings by partner position.
            var closing = _closures[u].Where(b => _digitOf.ContainsKey(b)).OrderBy(b => _digitOf[b]).ToList();
            var opening = _closures[u].Where(b => !_digitOf.ContainsKey(b)).OrderBy(b => _pos[mol.Bonds[b].Other(u)]).ToList();
            var outOrder = new List<int>();
            if (_parent[u] >= 0) outOrder.Add(_parent[u]);
            bool hasH = mol.Atoms[u].HydrogenCount > 0;
            if (hasH) outOrder.Add(-1);
            foreach (int b in closing) outOrder.Add(mol.Bonds[b].Other(u));
            foreach (int b in opening) outOrder.Add(mol.Bonds[b].Other(u));
            outOrder.AddRange(_children[u]);

            if (_parentBond[u] >= 0) _sb.Append(BondSymbol(_parentBond[u], _parent[u], u));
            WriteAtom(u, outOrder);
            foreach (int b in closing)
            {
                int d = _digitOf[b];
                _digitOf.Remove(b);
                AppendDigit(d);
                _freeDigits.Add(d);
            }
            foreach (int b in opening)
            {
                int d = _freeDigits.Min;
                _freeDigits.Remove(d);
                _digitOf[b] = d;
                _sb.Append(BondSymbol(b, u, mol.Bonds[b].Other(u)));
                AppendDigit(d);
            }
            var ch = _children[u];
            for (int i = 0; i < ch.Count; i++)
            {
                if (i < ch.Count - 1) { _sb.Append('('); Emit(ch[i]); _sb.Append(')'); }
                else Emit(ch[i]);
            }
        }

        private void AppendDigit(int d)
        {
            if (d < 10) _sb.Append((char)('0' + d));
            else _sb.Append('%').Append(d.ToString("00"));
        }

        private string BondSymbol(int b, int from, int to)
        {
            var bond = mol.Bonds[b];
            var aFrom = mol.Atoms[from];
            var aTo = mol.Atoms[to];
            switch (bond.Order)
            {
                case BondOrder.Double: return "=";
                case BondOrder.Triple: return "#";
                case BondOrder.Aromatic: return aFrom.IsAromatic && aTo.IsAromatic ? "" : ":";
                default:
                    if (_dir.TryGetValue(b, out char c)) return c.ToString();
                    return aFrom.IsAromatic && aTo.IsAromatic ? "-" : "";
            }
        }

        private void WriteAtom(int u, List<int> outOrder)
        {
            var a = mol.Atoms[u];
            var chir = ChiralTag.None;
            if (stereo && a.Chirality != ChiralTag.None)
            {
                var reference = mol.GetChiralReference(u);
                if (reference is not null && SameSet(reference, outOrder))
                {
                    bool odd = PermutationParityOdd(reference, outOrder);
                    chir = odd ? Flip(a.Chirality) : a.Chirality;
                }
            }
            int totalH = a.HydrogenCount;
            string symbol = a.AtomicNumber == 0 ? "*" : a.Symbol;
            if (a.IsAromatic) symbol = symbol.ToLowerInvariant();
            bool bracket = a.FormalCharge != 0 || a.Isotope != 0 || chir != ChiralTag.None || a.MapNumber != 0 || a.AtomicNumber == 1;
            if (!bracket)
            {
                int def = OrganicSubsetHydrogens(mol, u);
                if (def != totalH) bracket = true;
                if (a.IsAromatic && a.AtomicNumber is not (5 or 6 or 7 or 8 or 15 or 16 or 0)) bracket = true;
            }
            if (!bracket) { _sb.Append(symbol); return; }
            _sb.Append('[');
            if (a.Isotope != 0) _sb.Append(a.Isotope);
            _sb.Append(symbol);
            if (chir == ChiralTag.CounterClockwise) _sb.Append('@');
            else if (chir == ChiralTag.Clockwise) _sb.Append("@@");
            if (totalH > 0 && a.AtomicNumber != 1) { _sb.Append('H'); if (totalH > 1) _sb.Append(totalH); }
            else if (a.AtomicNumber == 1 && totalH > 0) { _sb.Append('H'); if (totalH > 1) _sb.Append(totalH); }
            if (a.FormalCharge != 0)
            {
                _sb.Append(a.FormalCharge > 0 ? '+' : '-');
                if (Math.Abs(a.FormalCharge) > 1) _sb.Append(Math.Abs(a.FormalCharge));
            }
            if (a.MapNumber != 0) _sb.Append(':').Append(a.MapNumber);
            _sb.Append(']');
        }

        private static ChiralTag Flip(ChiralTag t) => t == ChiralTag.Clockwise ? ChiralTag.CounterClockwise : ChiralTag.Clockwise;

        private static bool SameSet(int[] a, List<int> b) => a.Length == b.Count && a.All(b.Contains);

        private static bool PermutationParityOdd(int[] reference, List<int> output)
        {
            var perm = output.Select(x => Array.IndexOf(reference, x)).ToArray();
            int swaps = 0;
            for (int i = 0; i < perm.Length; i++)
                for (int j = i + 1; j < perm.Length; j++)
                    if (perm[i] > perm[j]) swaps++;
            return (swaps & 1) == 1;
        }

        // ---- cis/trans ---------------------------------------------------------------------------

        private bool? GetUp(int a, int x, int bond)
        {
            if (!_dir.TryGetValue(bond, out char c)) return null;
            return (c == '/') == (_pos[a] < _pos[x]);
        }

        private void SetUp(int a, int x, int bond, bool up) => _dir[bond] = up == (_pos[a] < _pos[x]) ? '/' : '\\';

        private void AssignDoubleBondDirections()
        {
            var bonds = mol.Bonds;
            for (int bi = 0; bi < bonds.Length; bi++)
            {
                var db = bonds[bi];
                if (db.Order != BondOrder.Double || db.Stereo is not (BondStereo.Cis or BondStereo.Trans)) continue;
                int a = db.Begin, b = db.End;
                var sideA = Candidates(a, b);
                var sideB = Candidates(b, a);
                if (sideA.Count == 0 || sideB.Count == 0) continue;
                bool done = false;
                foreach (var (x, bx) in sideA)
                {
                    foreach (var (y, by) in sideB)
                    {
                        bool flip = (x != db.StereoAtomBegin) ^ (y != db.StereoAtomEnd);
                        bool trans = (db.Stereo == BondStereo.Trans) ^ flip;
                        bool upA = GetUp(a, x, bx) ?? false;
                        bool needB = trans ? !upA : upA;
                        bool? haveB = GetUp(b, y, by);
                        if (haveB.HasValue && haveB.Value != needB) continue;
                        if (!_dir.ContainsKey(bx)) SetUp(a, x, bx, upA);
                        if (!haveB.HasValue) SetUp(b, y, by, needB);
                        done = true;
                        break;
                    }
                    if (done) break;
                }
            }
        }

        private List<(int Atom, int Bond)> Candidates(int a, int other)
        {
            var list = new List<(int, int, int)>();
            var nb = mol.Neighbors(a);
            var nbb = mol.NeighborBonds(a);
            for (int k = 0; k < nb.Length; k++)
            {
                if (nb[k] == other) continue;
                var bond = mol.Bonds[nbb[k]];
                if (bond.Order != BondOrder.Single) continue;
                int priority = _dir.ContainsKey(nbb[k]) ? 0 : _closureBonds.Contains(nbb[k]) ? 2 : 1;
                list.Add((nb[k], nbb[k], priority));
            }
            return list.OrderBy(t => t.Item3).ThenBy(t => _pos[t.Item1]).Select(t => (t.Item1, t.Item2)).ToList();
        }
    }
}
