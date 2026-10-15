using System.Runtime.CompilerServices;

namespace Biotech.Net.Chemistry;

/// <summary>
/// A compiled SMARTS query. Supported: atom primitives <c>*</c> <c>a</c> <c>A</c>, element symbols (aliphatic
/// upper case, aromatic lower case), <c>#n</c>, isotope, <c>D</c> <c>X</c> <c>H</c> <c>h</c> <c>v</c>
/// <c>R</c> <c>r</c> <c>x</c> with optional counts, charges (<c>+</c>, <c>-2</c>, <c>++</c>), recursive
/// <c>$(...)</c>, logical <c>!</c> <c>&amp;</c> <c>,</c> <c>;</c> (with SMARTS precedence); bond primitives
/// <c>- = # : ~ @</c> (plus <c>/</c> <c>\</c> as single) with the same logic; ring closures, branches and
/// <c>.</c>. Chirality and atom-map labels are accepted and ignored.
/// </summary>
/// <remarks>
/// Matching is a VF2-style depth-first search: query atoms are ordered so each one (after the first of its
/// component) is adjacent to an already-mapped atom, so candidates come from the neighbours of the mapped
/// atom. With <c>mergeHydrogens</c> (default) explicit <c>[#1]</c> query atoms are folded into an
/// "at least n hydrogens" constraint on their neighbour (like RDKit's MergeQueryHs), so patterns written with
/// explicit hydrogens (e.g. PAINS) match hydrogen-suppressed molecules.
/// </remarks>
public sealed class SmartsPattern
{
    private readonly List<AtomExpr> _atoms = [];
    private readonly List<int> _maps = [];
    private readonly List<(int A, int B, BondExpr E)> _bonds = [];
    private int[] _order = [];
    private int[] _parentQ = [];
    private int[] _parentBond = [];
    private List<(int Other, int Bond)>[] _adj = [];

    private SmartsPattern(string text) => Text = text;

    /// <summary>The source SMARTS.</summary>
    public string Text { get; }

    /// <summary>Number of query atoms.</summary>
    public int AtomCount => _atoms.Count;

    /// <summary>Atom-map number of query atom <paramref name="queryAtom"/> (<c>[C:1]</c> gives 1), or 0 when unmapped.</summary>
    public int AtomMap(int queryAtom) => _maps[queryAtom];

    /// <summary>Parses SMARTS.</summary>
    /// <exception cref="FormatException">The SMARTS is malformed or uses unsupported syntax.</exception>
    public static SmartsPattern Parse(string smarts, bool mergeHydrogens = true)
    {
        ArgumentNullException.ThrowIfNull(smarts);
        var p = new SmartsPattern(smarts);
        new SmartsParser(smarts, p).Run();
        if (mergeHydrogens) p.MergeHydrogens();
        p.Compile();
        return p;
    }

    /// <summary>Parses SMARTS, returning null when it is malformed.</summary>
    public static SmartsPattern? TryParse(string smarts, bool mergeHydrogens = true)
    {
        try { return Parse(smarts, mergeHydrogens); }
        catch (FormatException) { return null; }
    }

    /// <summary>Whether the pattern occurs in <paramref name="mol"/>.</summary>
    public bool IsMatch(Molecule mol) => FindMatches(mol, unique: false, max: 1, new MatchContext(mol), fixedRoot: -1).Count > 0;

    /// <summary>All matches (query atom i → molecule atom) — unique by atom set when <paramref name="unique"/>.</summary>
    public IReadOnlyList<int[]> Matches(Molecule mol, bool unique = true, int maxMatches = 10_000)
        => FindMatches(mol, unique, maxMatches, new MatchContext(mol), fixedRoot: -1);

    /// <summary>Number of unique matches.</summary>
    public int CountMatches(Molecule mol) => Matches(mol).Count;

    /// <summary>Whether query atom 0 can be mapped onto <paramref name="atom"/> in a full match.</summary>
    public bool MatchesAtom(Molecule mol, int atom) => FindMatches(mol, false, 1, new MatchContext(mol), atom).Count > 0;

    /// <summary>Atoms of <paramref name="mol"/> that query atom 0 can map to (recursive-SMARTS style atom typing).</summary>
    public IEnumerable<int> MatchingAtoms(Molecule mol)
    {
        var ctx = new MatchContext(mol);
        for (int i = 0; i < mol.AtomCount; i++)
            if (FindMatches(mol, false, 1, ctx, i).Count > 0) yield return i;
    }

    public override string ToString() => Text;

    // ---- matching ----------------------------------------------------------------------------------

    internal sealed class MatchContext(Molecule mol)
    {
        public readonly Molecule Mol = mol;
        public readonly RingInfo Rings = mol.Rings;
        public Dictionary<(SmartsPattern, int), bool>? Recursive;
    }

    internal List<int[]> FindMatches(Molecule mol, bool unique, int max, MatchContext ctx, int fixedRoot)
    {
        var results = new List<int[]>();
        int qn = _atoms.Count;
        if (qn == 0 || mol.AtomCount == 0) return results;
        var map = new int[qn];
        Array.Fill(map, -1);
        var used = new bool[mol.AtomCount];
        HashSet<string>? seen = unique ? [] : null;
        Search(0);
        return results;

        bool Search(int depth)
        {
            if (depth == qn)
            {
                var m = (int[])map.Clone();
                if (seen is not null)
                {
                    var key = string.Join(",", m.OrderBy(x => x));
                    if (!seen.Add(key)) return false;
                }
                results.Add(m);
                return results.Count >= max;
            }
            int q = _order[depth];
            int pq = _parentQ[q];
            if (pq < 0)
            {
                if (depth == 0 && fixedRoot >= 0)
                {
                    if (TryAssign(q, fixedRoot) && Search(depth + 1)) return true;
                    if (map[q] >= 0) Unassign(q);
                    return false;
                }
                for (int t = 0; t < mol.AtomCount; t++)
                {
                    if (used[t]) continue;
                    if (TryAssign(q, t))
                    {
                        if (Search(depth + 1)) return true;
                        Unassign(q);
                    }
                }
                return false;
            }
            int pt = map[pq];
            var nb = mol.Neighbors(pt);
            var nbb = mol.NeighborBonds(pt);
            var be = _bonds[_parentBond[q]].E;
            for (int k = 0; k < nb.Length; k++)
            {
                int t = nb[k];
                if (used[t]) continue;
                if (!be.Match(mol, nbb[k])) continue;
                if (TryAssign(q, t))
                {
                    if (Search(depth + 1)) return true;
                    Unassign(q);
                }
                nb = mol.Neighbors(pt);
                nbb = mol.NeighborBonds(pt);
            }
            return false;
        }

        bool TryAssign(int q, int t)
        {
            if (used[t]) return false;
            if (!_atoms[q].Match(t, ctx)) return false;
            // check bonds to other mapped query atoms (excluding the parent bond, already checked)
            foreach (var (other, bi) in _adj[q])
            {
                int ot = map[other];
                if (ot < 0) continue;
                if (bi == _parentBond[q]) continue;
                int tb = mol.BondBetween(t, ot);
                if (tb < 0 || !_bonds[bi].E.Match(mol, tb)) return false;
            }
            map[q] = t;
            used[t] = true;
            return true;
        }

        void Unassign(int q)
        {
            used[map[q]] = false;
            map[q] = -1;
        }
    }

    private void Compile()
    {
        int n = _atoms.Count;
        _adj = Enumerable.Range(0, n).Select(_ => new List<(int, int)>()).ToArray();
        for (int i = 0; i < _bonds.Count; i++)
        {
            _adj[_bonds[i].A].Add((_bonds[i].B, i));
            _adj[_bonds[i].B].Add((_bonds[i].A, i));
        }
        _order = new int[n];
        _parentQ = Enumerable.Repeat(-1, n).ToArray();
        _parentBond = Enumerable.Repeat(-1, n).ToArray();
        var visited = new bool[n];
        int k = 0;
        for (int s = 0; s < n; s++)
        {
            if (visited[s]) continue;
            var stack = new Stack<int>();
            stack.Push(s);
            visited[s] = true;
            while (stack.Count > 0)
            {
                int u = stack.Pop();
                _order[k++] = u;
                for (int j = _adj[u].Count - 1; j >= 0; j--)
                {
                    var (v, b) = _adj[u][j];
                    if (visited[v]) continue;
                    visited[v] = true;
                    _parentQ[v] = u;
                    _parentBond[v] = b;
                    stack.Push(v);
                }
            }
        }
    }

    private void MergeHydrogens()
    {
        int n = _atoms.Count;
        if (n <= 1) return;
        var degree = new int[n];
        foreach (var b in _bonds) { degree[b.A]++; degree[b.B]++; }
        var remove = new HashSet<int>();
        var extraH = new int[n];
        for (int i = 0; i < n; i++)
        {
            if (degree[i] != 1 || !_atoms[i].IsPlainHydrogen) continue;
            var bond = _bonds.First(b => b.A == i || b.B == i);
            int other = bond.A == i ? bond.B : bond.A;
            if (_atoms[other].IsPlainHydrogen) continue;
            remove.Add(i);
            extraH[other]++;
        }
        if (remove.Count == 0) return;
        var map = new int[n];
        var atoms = new List<AtomExpr>();
        var maps = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (remove.Contains(i)) { map[i] = -1; continue; }
            map[i] = atoms.Count;
            var e = _atoms[i];
            if (extraH[i] > 0) e = AtomExpr.And(e, new AtomExpr(AtomExpr.K.MinTotalH, extraH[i]));
            atoms.Add(e);
            maps.Add(_maps[i]);
        }
        var bonds = _bonds.Where(b => !remove.Contains(b.A) && !remove.Contains(b.B)).Select(b => (map[b.A], map[b.B], b.E)).ToList();
        _atoms.Clear(); _atoms.AddRange(atoms);
        _maps.Clear(); _maps.AddRange(maps);
        _bonds.Clear(); _bonds.AddRange(bonds);
    }

    // ---- expressions -------------------------------------------------------------------------------

    internal sealed class AtomExpr
    {
        internal enum K { True, AtomicNum, Aromatic, Aliphatic, Degree, TotalH, MinTotalH, ImplicitH, MinImplicitH, Connectivity, Valence, RingCount, InRing, RingSize, RingConn, MinRingConn, Charge, Isotope, Recursive, Not, And, Or }

        public readonly K Kind;
        public readonly int Value;
        public AtomExpr? L, R;
        public SmartsPattern? Rec;

        public AtomExpr(K kind, int value = 0) { Kind = kind; Value = value; }

        public static AtomExpr And(AtomExpr a, AtomExpr b) => new(K.And) { L = a, R = b };
        public static AtomExpr Or(AtomExpr a, AtomExpr b) => new(K.Or) { L = a, R = b };
        public static AtomExpr Not(AtomExpr a) => new(K.Not) { L = a };

        public bool IsPlainHydrogen => Kind == K.AtomicNum && Value == 1;

        public bool Match(int atom, MatchContext ctx)
        {
            var mol = ctx.Mol;
            switch (Kind)
            {
                case K.True: return true;
                case K.AtomicNum: return mol.Atoms[atom].AtomicNumber == Value;
                case K.Aromatic: return mol.Atoms[atom].IsAromatic;
                case K.Aliphatic: return !mol.Atoms[atom].IsAromatic;
                case K.Degree: return mol.Degree(atom) == Value;
                case K.TotalH: return mol.TotalHydrogens(atom) == Value;
                case K.MinTotalH: return mol.TotalHydrogens(atom) >= Value;
                case K.ImplicitH: return mol.Atoms[atom].HydrogenCount == Value;
                case K.MinImplicitH: return mol.Atoms[atom].HydrogenCount >= Value;
                case K.Connectivity: return mol.Degree(atom) + mol.Atoms[atom].HydrogenCount == Value;
                case K.Valence: return mol.TotalValence(atom) == Value;
                case K.RingCount: return ctx.Rings.AtomRingCount(atom) == Value;
                case K.InRing: return ctx.Rings.IsAtomInRing(atom);
                case K.RingSize: return ctx.Rings.IsAtomInRingOfSize(atom, Value);
                case K.RingConn: return ctx.Rings.RingBondCount(atom) == Value;
                case K.MinRingConn: return ctx.Rings.RingBondCount(atom) >= Value;
                case K.Charge: return mol.Atoms[atom].FormalCharge == Value;
                case K.Isotope: return mol.Atoms[atom].Isotope == Value;
                case K.Recursive:
                    {
                        ctx.Recursive ??= [];
                        if (ctx.Recursive.TryGetValue((Rec!, atom), out bool cached)) return cached;
                        bool r = Rec!.FindMatches(mol, false, 1, ctx, atom).Count > 0;
                        ctx.Recursive[(Rec, atom)] = r;
                        return r;
                    }
                case K.Not: return !L!.Match(atom, ctx);
                case K.And: return L!.Match(atom, ctx) && R!.Match(atom, ctx);
                case K.Or: return L!.Match(atom, ctx) || R!.Match(atom, ctx);
                default: return false;
            }
        }
    }

    internal sealed class BondExpr
    {
        internal enum K { Default, Single, Double, Triple, Aromatic, Any, Ring, Not, And, Or }

        public readonly K Kind;
        public BondExpr? L, R;

        public BondExpr(K kind) => Kind = kind;

        public static readonly BondExpr DefaultBond = new(K.Default);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool Match(Molecule mol, int bond)
        {
            var order = mol.Bonds[bond].Order;
            return Kind switch
            {
                K.Default => order is BondOrder.Single or BondOrder.Aromatic,
                K.Single => order == BondOrder.Single,
                K.Double => order == BondOrder.Double,
                K.Triple => order == BondOrder.Triple,
                K.Aromatic => order == BondOrder.Aromatic,
                K.Any => true,
                K.Ring => mol.Rings.IsBondInRing(bond),
                K.Not => !L!.Match(mol, bond),
                K.And => L!.Match(mol, bond) && R!.Match(mol, bond),
                K.Or => L!.Match(mol, bond) || R!.Match(mol, bond),
                _ => false,
            };
        }
    }

    // ---- parser ------------------------------------------------------------------------------------

    private sealed class SmartsParser(string s, SmartsPattern target)
    {
        private int _pos;

        public void Run()
        {
            int prev = -1;
            BondExpr? bond = null;
            var branch = new Stack<int>();
            var rings = new Dictionary<int, (int Atom, BondExpr? Bond)>();
            while (_pos < s.Length)
            {
                char c = s[_pos];
                if (c == '(') { if (prev < 0) Fail("branch without atom"); branch.Push(prev); _pos++; continue; }
                if (c == ')') { if (branch.Count == 0) Fail("unbalanced ')'"); prev = branch.Pop(); _pos++; continue; }
                if (c == '.') { prev = -1; _pos++; continue; }
                if (IsBondStart(c)) { bond = ParseBondLowAnd(); continue; }
                if (char.IsAsciiDigit(c) || c == '%')
                {
                    if (prev < 0) Fail("ring closure without atom");
                    int num;
                    if (c == '%') { num = int.Parse(s.AsSpan(_pos + 1, 2)); _pos += 3; }
                    else { num = c - '0'; _pos++; }
                    if (rings.Remove(num, out var open))
                        target._bonds.Add((open.Atom, prev, bond ?? open.Bond ?? BondExpr.DefaultBond));
                    else rings[num] = (prev, bond);
                    bond = null;
                    continue;
                }
                _currentMap = 0;
                var expr = ParseAtom();
                target._atoms.Add(expr);
                target._maps.Add(_currentMap);
                int idx = target._atoms.Count - 1;
                if (prev >= 0) target._bonds.Add((prev, idx, bond ?? BondExpr.DefaultBond));
                bond = null;
                prev = idx;
            }
            if (branch.Count > 0) Fail("unclosed branch");
            if (rings.Count > 0) Fail("unclosed ring");
        }

        private static bool IsBondStart(char c) => c is '-' or '=' or '#' or ':' or '~' or '@' or '!' or '/' or '\\' or '&' or ',' or ';';

        private BondExpr ParseBondLowAnd()
        {
            var e = ParseBondOr();
            while (Peek() == ';') { _pos++; e = new BondExpr(BondExpr.K.And) { L = e, R = ParseBondOr() }; }
            return e;
        }

        private BondExpr ParseBondOr()
        {
            var e = ParseBondHighAnd();
            while (Peek() == ',') { _pos++; e = new BondExpr(BondExpr.K.Or) { L = e, R = ParseBondHighAnd() }; }
            return e;
        }

        private BondExpr ParseBondHighAnd()
        {
            var e = ParseBondUnary();
            while (true)
            {
                char c = Peek();
                if (c == '&') { _pos++; e = new BondExpr(BondExpr.K.And) { L = e, R = ParseBondUnary() }; continue; }
                if (c is '-' or '=' or '#' or ':' or '~' or '@' or '!' or '/' or '\\') { e = new BondExpr(BondExpr.K.And) { L = e, R = ParseBondUnary() }; continue; }
                return e;
            }
        }

        private BondExpr ParseBondUnary()
        {
            char c = Peek();
            _pos++;
            return c switch
            {
                '!' => new BondExpr(BondExpr.K.Not) { L = ParseBondUnary() },
                '-' or '/' or '\\' => new BondExpr(BondExpr.K.Single),
                '=' => new BondExpr(BondExpr.K.Double),
                '#' => new BondExpr(BondExpr.K.Triple),
                ':' => new BondExpr(BondExpr.K.Aromatic),
                '~' => new BondExpr(BondExpr.K.Any),
                '@' => new BondExpr(BondExpr.K.Ring),
                _ => throw Fail($"unexpected bond character '{c}'"),
            };
        }

        private AtomExpr ParseAtom()
        {
            char c = s[_pos];
            if (c == '[')
            {
                _pos++;
                _bracketStart = _pos;
                var e = ParseLowAnd();
                if (Peek() != ']') Fail($"expected ']' at {_pos}");
                _pos++;
                return e;
            }
            if (c == '*') { _pos++; return new AtomExpr(AtomExpr.K.True); }
            if (c == 'a') { _pos++; return new AtomExpr(AtomExpr.K.Aromatic); }
            if (c == 'A') { _pos++; return new AtomExpr(AtomExpr.K.Aliphatic); }
            if (c == 'C' && Peek(1) == 'l') { _pos += 2; return Element(17, false); }
            if (c == 'B' && Peek(1) == 'r') { _pos += 2; return Element(35, false); }
            int z = c switch { 'B' => 5, 'C' => 6, 'N' => 7, 'O' => 8, 'P' => 15, 'S' => 16, 'F' => 9, 'I' => 53, 'b' => 5, 'c' => 6, 'n' => 7, 'o' => 8, 'p' => 15, 's' => 16, _ => -1 };
            if (z < 0) Fail($"unexpected character '{c}' at {_pos}");
            _pos++;
            return Element(z, char.IsLower(c));
        }

        private static AtomExpr Element(int z, bool aromatic)
            => AtomExpr.And(new AtomExpr(AtomExpr.K.AtomicNum, z), new AtomExpr(aromatic ? AtomExpr.K.Aromatic : AtomExpr.K.Aliphatic));

        private int _bracketStart;

        private AtomExpr ParseLowAnd()
        {
            var e = ParseOr();
            while (Peek() == ';') { _pos++; e = AtomExpr.And(e, ParseOr()); }
            return e;
        }

        private AtomExpr ParseOr()
        {
            var e = ParseHighAnd();
            while (Peek() == ',') { _pos++; e = AtomExpr.Or(e, ParseHighAnd()); }
            return e;
        }

        private AtomExpr ParseHighAnd()
        {
            var e = ParseUnary();
            while (true)
            {
                char c = Peek();
                if (c is ']' or ';' or ',' or ')' or '\0') return e;
                if (c == '&') { _pos++; }
                if (Peek() == ':') { SkipMap(); continue; }
                e = AtomExpr.And(e, ParseUnary());
            }
        }

        private int _currentMap;

        private void SkipMap()
        {
            _pos++;
            int v = 0;
            while (char.IsAsciiDigit(Peek())) v = v * 10 + (s[_pos++] - '0');
            _currentMap = v;
        }

        private AtomExpr ParseUnary()
        {
            if (Peek() == '!') { _pos++; return AtomExpr.Not(ParseUnary()); }
            return ParsePrimitive();
        }

        private AtomExpr ParsePrimitive()
        {
            bool atStart = _pos == _bracketStart;
            char c = Peek();
            if (c == '\0') throw Fail("unexpected end");
            if (char.IsAsciiDigit(c))
            {
                int iso = ReadInt(0);
                _bracketStart = _pos; // allow "[2H]"
                return new AtomExpr(AtomExpr.K.Isotope, iso);
            }
            switch (c)
            {
                case '*': _pos++; return new AtomExpr(AtomExpr.K.True);
                case '#': _pos++; return new AtomExpr(AtomExpr.K.AtomicNum, ReadInt(-1) is var z and >= 0 ? z : throw Fail("'#' needs a number"));
                case '$':
                    {
                        _pos++;
                        if (Peek() != '(') throw Fail("expected '(' after '$'");
                        int depth = 0, start = _pos + 1;
                        do
                        {
                            if (Peek() == '(') depth++;
                            else if (Peek() == ')') depth--;
                            else if (Peek() == '\0') throw Fail("unterminated recursive SMARTS");
                            _pos++;
                        } while (depth > 0);
                        var inner = Parse(s[start..(_pos - 1)], mergeHydrogens: true);
                        return new AtomExpr(AtomExpr.K.Recursive) { Rec = inner };
                    }
                case '+' or '-':
                    {
                        int sign = c == '+' ? 1 : -1;
                        _pos++;
                        int v = ReadInt(-1);
                        if (v < 0) { v = 1; while (Peek() == c) { v++; _pos++; } }
                        return new AtomExpr(AtomExpr.K.Charge, sign * v);
                    }
                case '@':
                    _pos++;
                    if (Peek() == '@') _pos++;
                    if (Peek() == '?') _pos++;
                    return new AtomExpr(AtomExpr.K.True);
            }
            if (char.IsAsciiLetterUpper(c))
            {
                // Hydrogen atom: "[H]", "[H+]", "[2H]".
                if (c == 'H' && atStart && Peek(1) is ']' or '+' or '-' or ';' or ',' or '&')
                {
                    _pos++;
                    return new AtomExpr(AtomExpr.K.AtomicNum, 1);
                }
                char lower = Peek(1);
                if (char.IsAsciiLetterLower(lower) && PeriodicTable.TryGet(new string([c, lower]), out var e2))
                {
                    _pos += 2;
                    return Element(e2.AtomicNumber, false);
                }
                switch (c)
                {
                    case 'A': _pos++; return new AtomExpr(AtomExpr.K.Aliphatic);
                    case 'D': _pos++; return new AtomExpr(AtomExpr.K.Degree, ReadInt(1));
                    case 'X': _pos++; return new AtomExpr(AtomExpr.K.Connectivity, ReadInt(1));
                    case 'H': _pos++; return new AtomExpr(AtomExpr.K.TotalH, ReadInt(1));
                    case 'R':
                        {
                            _pos++;
                            int v = ReadInt(-1);
                            if (v < 0) return new AtomExpr(AtomExpr.K.InRing);
                            return v == 0 ? AtomExpr.Not(new AtomExpr(AtomExpr.K.InRing)) : new AtomExpr(AtomExpr.K.RingCount, v);
                        }
                }
                if (PeriodicTable.TryGet(c.ToString(), out var e1)) { _pos++; return Element(e1.AtomicNumber, false); }
                throw Fail($"unknown primitive '{c}'");
            }
            if (char.IsAsciiLetterLower(c))
            {
                string two = new([c, Peek(1)]);
                if (two is "se" or "as" or "te") { _pos += 2; return Element(two == "se" ? 34 : two == "as" ? 33 : 52, true); }
                switch (c)
                {
                    case 'a': _pos++; return new AtomExpr(AtomExpr.K.Aromatic);
                    case 'c': _pos++; return Element(6, true);
                    case 'n': _pos++; return Element(7, true);
                    case 'o': _pos++; return Element(8, true);
                    case 's': _pos++; return Element(16, true);
                    case 'p': _pos++; return Element(15, true);
                    case 'b': _pos++; return Element(5, true);
                    case 'h':
                        {
                            _pos++;
                            int v = ReadInt(-1);
                            return v < 0 ? new AtomExpr(AtomExpr.K.MinImplicitH, 1) : new AtomExpr(AtomExpr.K.ImplicitH, v);
                        }
                    case 'v': _pos++; return new AtomExpr(AtomExpr.K.Valence, ReadInt(1));
                    case 'r':
                        {
                            _pos++;
                            int v = ReadInt(-1);
                            if (v < 0) return new AtomExpr(AtomExpr.K.InRing);
                            return v == 0 ? AtomExpr.Not(new AtomExpr(AtomExpr.K.InRing)) : new AtomExpr(AtomExpr.K.RingSize, v);
                        }
                    case 'x':
                        {
                            _pos++;
                            int v = ReadInt(-1);
                            return v < 0 ? new AtomExpr(AtomExpr.K.MinRingConn, 1) : new AtomExpr(AtomExpr.K.RingConn, v);
                        }
                }
            }
            throw Fail($"unexpected '{c}' at {_pos}");
        }

        private int ReadInt(int fallback)
        {
            if (!char.IsAsciiDigit(Peek())) return fallback;
            int v = 0;
            while (char.IsAsciiDigit(Peek())) v = v * 10 + (s[_pos++] - '0');
            return v;
        }

        private char Peek(int offset = 0) => _pos + offset < s.Length ? s[_pos + offset] : '\0';

        private FormatException Fail(string message) => throw new FormatException($"Invalid SMARTS '{s}': {message}.");
    }
}

/// <summary>Convenience substructure search API.</summary>
public static class Substructure
{
    /// <summary>Whether <paramref name="mol"/> contains the SMARTS <paramref name="smarts"/>.</summary>
    public static bool HasMatch(Molecule mol, string smarts) => SmartsPattern.Parse(smarts).IsMatch(mol);

    /// <summary>Whether <paramref name="mol"/> contains <paramref name="query"/> (a molecule used as a query: element, aromaticity and bond order must agree; charges must agree when the query atom is charged).</summary>
    public static bool HasMatch(Molecule mol, Molecule query) => FromMolecule(query).IsMatch(mol);

    /// <summary>Unique matches of a SMARTS pattern.</summary>
    public static IReadOnlyList<int[]> FindAll(Molecule mol, string smarts) => SmartsPattern.Parse(smarts).Matches(mol);

    /// <summary>Converts a molecule into an equivalent SMARTS query (via its canonical SMILES written as SMARTS-compatible text).</summary>
    public static SmartsPattern FromMolecule(Molecule query)
    {
        var sb = new System.Text.StringBuilder();
        // Build an explicit SMARTS: each atom as [#Z;a/A;+c], each bond explicit.
        var q = query;
        var atomText = new string[q.AtomCount];
        for (int i = 0; i < q.AtomCount; i++)
        {
            var a = q.Atoms[i];
            var t = a.AtomicNumber == 0 ? "*" : $"#{a.AtomicNumber};{(a.IsAromatic ? "a" : "A")}";
            if (a.FormalCharge != 0) t += $";{(a.FormalCharge > 0 ? "+" : "-")}{Math.Abs(a.FormalCharge)}";
            atomText[i] = "[" + t + "]";
        }
        // Emit as a DFS with ring closures.
        var visited = new bool[q.AtomCount];
        var ringDigits = new Dictionary<int, int>();
        int nextDigit = 1;
        var closures = new List<int>[q.AtomCount];
        for (int i = 0; i < closures.Length; i++) closures[i] = [];
        var treeBonds = new HashSet<int>();
        var order = new List<int>();
        for (int s = 0; s < q.AtomCount; s++)
        {
            if (visited[s]) continue;
            Visit(s, -1);
        }
        void Visit(int u, int pb)
        {
            visited[u] = true;
            var nb = q.Neighbors(u).ToArray();
            var nbb = q.NeighborBonds(u).ToArray();
            for (int k = 0; k < nb.Length; k++)
            {
                if (nbb[k] == pb) continue;
                if (visited[nb[k]]) { if (!treeBonds.Contains(nbb[k]) && !closures[u].Contains(nbb[k])) { closures[u].Add(nbb[k]); closures[nb[k]].Add(nbb[k]); } continue; }
                treeBonds.Add(nbb[k]);
                Visit(nb[k], nbb[k]);
            }
        }
        var written = new bool[q.AtomCount];
        bool first = true;
        for (int s = 0; s < q.AtomCount; s++)
        {
            if (written[s]) continue;
            if (!first) sb.Append('.');
            first = false;
            Emit(s, -1);
        }
        return SmartsPattern.Parse(sb.ToString(), mergeHydrogens: false);

        string BondText(int b) => q.Bonds[b].Order switch
        {
            BondOrder.Double => "=",
            BondOrder.Triple => "#",
            BondOrder.Aromatic => ":",
            _ => "-",
        };

        void Emit(int u, int pb)
        {
            written[u] = true;
            sb.Append(atomText[u]);
            foreach (int b in closures[u])
            {
                if (ringDigits.TryGetValue(b, out int d)) { sb.Append(BondText(b)).Append('%').Append(d.ToString("00")); ringDigits.Remove(b); }
                else { d = nextDigit++; ringDigits[b] = d; sb.Append('%').Append(d.ToString("00")); }
            }
            var kids = new List<(int, int)>();
            var nb = q.Neighbors(u).ToArray();
            var nbb = q.NeighborBonds(u).ToArray();
            for (int k = 0; k < nb.Length; k++)
                if (nbb[k] != pb && treeBonds.Contains(nbb[k]) && !written[nb[k]]) kids.Add((nb[k], nbb[k]));
            for (int i = 0; i < kids.Count; i++)
            {
                bool branch = i < kids.Count - 1;
                if (branch) sb.Append('(');
                sb.Append(BondText(kids[i].Item2));
                Emit(kids[i].Item1, kids[i].Item2);
                if (branch) sb.Append(')');
            }
        }
    }
}
