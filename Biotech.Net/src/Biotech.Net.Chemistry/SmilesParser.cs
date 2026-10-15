namespace Biotech.Net.Chemistry;

/// <summary>Options for <see cref="SmilesParser"/>.</summary>
public sealed record SmilesParserOptions
{
    /// <summary>Re-derive aromaticity with <see cref="Aromaticity.Perceive"/> so Kekulé and aromatic inputs agree (default true).</summary>
    public bool PerceiveAromaticity { get; init; } = true;

    /// <summary>Default options.</summary>
    public static SmilesParserOptions Default { get; } = new();
}

/// <summary>
/// OpenSMILES parser: organic-subset and bracket atoms (isotope, chirality, H count, charge, atom class),
/// bonds <c>- = # $ : / \</c>, branches, ring closures (<c>1</c>…<c>9</c>, <c>%nn</c>), dot-disconnected
/// components and lower-case aromatic atoms. Implicit hydrogens are computed after kekulization; tetrahedral
/// <c>@</c>/<c>@@</c> are stored with their neighbour order and <c>/</c> <c>\</c> are converted to cis/trans
/// annotations on the double bond.
/// </summary>
public static class SmilesParser
{
    /// <summary>Parses a SMILES string. Text after the first whitespace is used as the name when <paramref name="name"/> is null.</summary>
    /// <exception cref="FormatException">The SMILES is malformed.</exception>
    public static Molecule Parse(string smiles, string? name = null, SmilesParserOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(smiles);
        options ??= SmilesParserOptions.Default;
        smiles = smiles.Trim();
        int ws = smiles.IndexOfAny([' ', '\t']);
        if (ws >= 0)
        {
            name ??= smiles[(ws + 1)..].Trim();
            smiles = smiles[..ws];
        }
        var p = new State(smiles, name ?? string.Empty);
        p.Run();
        var mol = p.Mol;
        mol.UpdateImplicitHydrogens();
        if (options.PerceiveAromaticity) Aromaticity.Perceive(mol);
        return mol;
    }

    /// <summary>Parses without throwing; returns null on malformed input.</summary>
    public static Molecule? TryParse(string smiles, string? name = null)
    {
        try { return Parse(smiles, name); }
        catch (FormatException) { return null; }
        catch (InvalidOperationException) { return null; }
    }

    private sealed class State(string s, string name)
    {
        public readonly Molecule Mol = new(name, Math.Max(4, s.Length));
        private int _pos;
        private readonly List<List<int>> _order = []; // neighbour order per atom (-1 = implicit H)
        private readonly Dictionary<int, (int Atom, char Bond, int Slot)> _rings = [];
        private readonly List<(int Bond, int From, char Dir)> _directional = [];

        public void Run()
        {
            int prev = -1;
            char bond = '\0';
            var branchStack = new Stack<int>();
            bool expectAtom = true;
            while (_pos < s.Length)
            {
                char c = s[_pos];
                switch (c)
                {
                    case '(':
                        if (prev < 0) Fail("branch without a preceding atom");
                        branchStack.Push(prev);
                        _pos++;
                        break;
                    case ')':
                        if (branchStack.Count == 0) Fail("unbalanced ')'");
                        if (bond != '\0') Fail("bond before ')'");
                        prev = branchStack.Pop();
                        _pos++;
                        break;
                    case '-' or '=' or '#' or '$' or ':' or '/' or '\\':
                        if (bond != '\0') Fail("two consecutive bonds");
                        bond = c;
                        _pos++;
                        break;
                    case '.':
                        if (bond != '\0') Fail("bond before '.'");
                        prev = -1;
                        _pos++;
                        break;
                    case '%':
                    case >= '0' and <= '9':
                        {
                            if (prev < 0) Fail("ring closure without an atom");
                            int num;
                            if (c == '%')
                            {
                                if (_pos + 2 >= s.Length || !char.IsAsciiDigit(s[_pos + 1]) || !char.IsAsciiDigit(s[_pos + 2])) Fail("bad %nn ring closure");
                                num = (s[_pos + 1] - '0') * 10 + (s[_pos + 2] - '0');
                                _pos += 3;
                            }
                            else { num = c - '0'; _pos++; }
                            RingClosure(prev, num, bond);
                            bond = '\0';
                            break;
                        }
                    default:
                        {
                            int atom = ReadAtom();
                            if (prev >= 0) { Connect(prev, atom, bond); _hasPreceding.Add(atom); }
                            else if (bond != '\0') Fail("bond without a preceding atom");
                            bond = '\0';
                            prev = atom;
                            expectAtom = false;
                            break;
                        }
                }
            }
            _ = expectAtom;
            if (bond != '\0') Fail("dangling bond at end");
            if (branchStack.Count > 0) Fail("unclosed branch");
            if (_rings.Count > 0) Fail($"unclosed ring {_rings.Keys.First()}");
            FinishStereo();
        }

        private void Connect(int a, int b, char bond)
        {
            var order = BondFor(a, b, bond);
            if (Mol.BondBetween(a, b) >= 0) Fail("duplicate bond");
            int bi = Mol.AddBond(a, b, order);
            if (bond is '/' or '\\') _directional.Add((bi, a, bond));
            _order[a].Add(b);
            _order[b].Add(a);
        }

        private BondOrder BondFor(int a, int b, char bond) => bond switch
        {
            '=' => BondOrder.Double,
            '#' => BondOrder.Triple,
            '$' => BondOrder.Triple, // quadruple not modelled
            ':' => BondOrder.Aromatic,
            '-' or '/' or '\\' => BondOrder.Single,
            _ => Mol.Atoms[a].IsAromatic && Mol.Atoms[b].IsAromatic ? BondOrder.Aromatic : BondOrder.Single,
        };

        private void RingClosure(int atom, int num, char bond)
        {
            if (_rings.TryGetValue(num, out var open))
            {
                _rings.Remove(num);
                if (open.Atom == atom) Fail("ring closure to itself");
                char use = bond != '\0' ? bond : open.Bond;
                if (bond != '\0' && open.Bond != '\0' && bond != open.Bond && !(IsDir(bond) && IsDir(open.Bond))) Fail("conflicting ring-closure bonds");
                if (Mol.BondBetween(open.Atom, atom) >= 0) Fail("duplicate ring-closure bond");
                int bi = Mol.AddBond(open.Atom, atom, BondFor(open.Atom, atom, use));
                if (bond is '/' or '\\') _directional.Add((bi, atom, bond));
                else if (open.Bond is '/' or '\\') _directional.Add((bi, open.Atom, open.Bond));
                _order[open.Atom][open.Slot] = atom;
                _order[atom].Add(open.Atom);
            }
            else
            {
                _order[atom].Add(int.MinValue); // placeholder
                _rings[num] = (atom, bond, _order[atom].Count - 1);
            }
        }

        private static bool IsDir(char c) => c is '/' or '\\';

        private int ReadAtom()
        {
            char c = s[_pos];
            if (c == '[') return ReadBracket();
            int z;
            bool aromatic = false;
            if (c == 'C' && Peek(1) == 'l') { z = 17; _pos += 2; }
            else if (c == 'B' && Peek(1) == 'r') { z = 35; _pos += 2; }
            else
            {
                z = c switch
                {
                    'B' => 5, 'C' => 6, 'N' => 7, 'O' => 8, 'P' => 15, 'S' => 16, 'F' => 9, 'I' => 53, '*' => 0,
                    'b' => 5, 'c' => 6, 'n' => 7, 'o' => 8, 'p' => 15, 's' => 16,
                    _ => -1,
                };
                if (z < 0) Fail($"unexpected character '{c}'");
                aromatic = char.IsLower(c);
                _pos++;
            }
            var atom = new Atom(z) { IsAromatic = aromatic, NoImplicit = z == 0 };
            _order.Add([]);
            return Mol.AddAtom(atom);
        }

        private int ReadBracket()
        {
            _pos++; // [
            int isotope = 0;
            while (_pos < s.Length && char.IsAsciiDigit(s[_pos])) isotope = isotope * 10 + (s[_pos++] - '0');
            if (_pos >= s.Length) Fail("unterminated bracket atom");
            int z;
            bool aromatic = false;
            char c = s[_pos];
            if (c == '*') { z = 0; _pos++; }
            else if (char.IsLower(c))
            {
                string two = _pos + 1 < s.Length ? s.Substring(_pos, 2) : "";
                if (two is "se" or "as" or "te") { z = two == "se" ? 34 : two == "as" ? 33 : 52; _pos += 2; }
                else
                {
                    z = c switch { 'b' => 5, 'c' => 6, 'n' => 7, 'o' => 8, 'p' => 15, 's' => 16, _ => -1 };
                    if (z < 0) Fail($"unknown aromatic symbol '{c}'");
                    _pos++;
                }
                aromatic = true;
            }
            else if (char.IsUpper(c))
            {
                z = -1;
                if (_pos + 1 < s.Length && char.IsLower(s[_pos + 1]) && PeriodicTable.TryGet(s.Substring(_pos, 2), out var e2))
                { z = e2.AtomicNumber; _pos += 2; }
                else if (PeriodicTable.TryGet(c.ToString(), out var e1)) { z = e1.AtomicNumber; _pos++; }
                if (z < 0) Fail($"unknown element at {_pos}");
            }
            else { Fail($"bad bracket atom at {_pos}"); return -1; }

            var chir = ChiralTag.None;
            if (Peek(0) == '@')
            {
                _pos++;
                chir = ChiralTag.CounterClockwise;
                if (Peek(0) == '@') { _pos++; chir = ChiralTag.Clockwise; }
                else if (Peek(0) == 'T' && Peek(1) == 'H') { _pos += 2; if (Peek(0) == '2') chir = ChiralTag.Clockwise; _pos++; }
                else if (char.IsAsciiLetterUpper(Peek(0)) && Peek(0) != 'H') { while (_pos < s.Length && (char.IsAsciiLetterUpper(s[_pos]) || char.IsAsciiDigit(s[_pos]))) _pos++; chir = ChiralTag.None; }
            }
            int h = 0;
            if (Peek(0) == 'H')
            {
                _pos++;
                h = 1;
                if (char.IsAsciiDigit(Peek(0))) { h = Peek(0) - '0'; _pos++; }
            }
            int charge = 0;
            if (Peek(0) is '+' or '-')
            {
                int sign = Peek(0) == '+' ? 1 : -1;
                _pos++;
                if (char.IsAsciiDigit(Peek(0)))
                {
                    int v = 0;
                    while (char.IsAsciiDigit(Peek(0))) v = v * 10 + (s[_pos++] - '0');
                    charge = sign * v;
                }
                else
                {
                    charge = sign;
                    while (Peek(0) == (sign > 0 ? '+' : '-')) { charge += sign; _pos++; }
                }
            }
            int map = 0;
            if (Peek(0) == ':')
            {
                _pos++;
                while (char.IsAsciiDigit(Peek(0))) map = map * 10 + (s[_pos++] - '0');
            }
            if (Peek(0) != ']') Fail($"expected ']' at {_pos}");
            _pos++;
            var atom = new Atom(z, charge)
            {
                Isotope = isotope,
                IsAromatic = aromatic,
                NoImplicit = true,
                ExplicitHydrogens = h,
                Chirality = chir,
                MapNumber = map,
            };
            _order.Add([]);
            int idx = Mol.AddAtom(atom);
            if (chir != ChiralTag.None && h > 0) _pendingHydrogenChiral.Add(idx);
            return idx;
        }

        private readonly HashSet<int> _pendingHydrogenChiral = [];
        private readonly HashSet<int> _hasPreceding = [];

        private void FinishStereo()
        {
            // Tetrahedral: neighbour order = [preceding atom], [implicit H], ring closures, branches/chain.
            for (int i = 0; i < Mol.AtomCount; i++)
            {
                var a = Mol.Atoms[i];
                if (a.Chirality == ChiralTag.None) continue;
                var order = new List<int>(_order[i]);
                if (_pendingHydrogenChiral.Contains(i))
                {
                    // H goes after the preceding atom when there is one: the preceding atom is the first neighbour
                    // added before any ring-closure/branch of this atom, i.e. an atom with a smaller index bonded via chain.
                    order.Insert(_hasPreceding.Contains(i) ? 1 : 0, -1);
                }
                int heavyAndH = order.Count;
                if (heavyAndH is < 3 or > 4) { Mol.AtomAt(i).Chirality = ChiralTag.None; continue; }
                Mol.SetChiralReference(i, order.ToArray());
            }

            // Double-bond stereo from directional single bonds.
            if (_directional.Count == 0) return;
            var dirOf = new Dictionary<int, (int From, char Dir)>();
            foreach (var (b, from, d) in _directional) dirOf[b] = (from, d);
            for (int bi = 0; bi < Mol.BondCount; bi++)
            {
                var bond = Mol.Bonds[bi];
                if (bond.Order != BondOrder.Double) continue;
                if (!SideReference(bond.Begin, bond.End, dirOf, out int x, out bool upX)) continue;
                if (!SideReference(bond.End, bond.Begin, dirOf, out int y, out bool upY)) continue;
                ref var db = ref Mol.BondAt(bi);
                db.Stereo = upX != upY ? BondStereo.Trans : BondStereo.Cis;
                db.StereoAtomBegin = x;
                db.StereoAtomEnd = y;
            }
        }

        /// <summary>For double-bond atom <paramref name="a"/> (partner <paramref name="other"/>), the directional neighbour and whether it is "up" relative to <paramref name="a"/>.</summary>
        private bool SideReference(int a, int other, Dictionary<int, (int From, char Dir)> dirOf, out int nbr, out bool up)
        {
            var nb = Mol.Neighbors(a);
            var nbb = Mol.NeighborBonds(a);
            for (int k = 0; k < nb.Length; k++)
            {
                if (nb[k] == other) continue;
                if (!dirOf.TryGetValue(nbb[k], out var d)) continue;
                nbr = nb[k];
                // '/' written as "p/q" means q is above p.
                bool slash = d.Dir == '/';
                up = d.From == a ? slash : !slash;
                return true;
            }
            nbr = -1; up = false;
            return false;
        }

        private char Peek(int offset) => _pos + offset < s.Length ? s[_pos + offset] : '\0';

        private void Fail(string message) => throw new FormatException($"Invalid SMILES '{s}': {message}.");
    }
}
