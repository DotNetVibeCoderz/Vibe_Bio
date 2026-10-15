using System.Globalization;
using System.Text;

namespace Biotech.Net.Chemistry;

/// <summary>
/// Bemis–Murcko scaffolds (J. Med. Chem. 1996, 39, 2887): ring systems plus the linkers between them, with side
/// chains removed. As in RDKit's <c>MurckoScaffold.GetScaffoldForMol</c>, atoms attached by a double bond to a
/// scaffold atom (e.g. carbonyl O) are kept by default. Acyclic molecules have an empty scaffold.
/// </summary>
public static class MurckoScaffold
{
    /// <summary>The scaffold as a new molecule (empty when <paramref name="mol"/> has no rings).</summary>
    public static Molecule GetScaffold(Molecule mol, bool keepExocyclicDoubleBonds = true)
    {
        var keep = ScaffoldAtoms(mol, keepExocyclicDoubleBonds);
        return Extract(mol, keep);
    }

    /// <summary>Canonical SMILES of the scaffold ("" for acyclic molecules).</summary>
    public static string GetScaffoldSmiles(Molecule mol, bool keepExocyclicDoubleBonds = true)
        => GetScaffold(mol, keepExocyclicDoubleBonds).ToSmiles(includeStereo: false);

    /// <summary>
    /// The generic framework (graph framework): the scaffold with every atom turned into an sp3 carbon and every bond
    /// into a single bond (RDKit <c>MakeScaffoldGeneric</c> applied to the Murcko scaffold, without exocyclic atoms).
    /// </summary>
    public static Molecule GetGenericFramework(Molecule mol)
    {
        var scaffold = GetScaffold(mol, keepExocyclicDoubleBonds: false);
        var g = new Molecule(mol.Name, scaffold.AtomCount);
        for (int i = 0; i < scaffold.AtomCount; i++) g.AddAtom(6);
        foreach (var b in scaffold.Bonds) g.AddBond(b.Begin, b.End);
        g.UpdateImplicitHydrogens();
        return g;
    }

    /// <summary>Canonical SMILES of the generic framework.</summary>
    public static string GetGenericFrameworkSmiles(Molecule mol) => GetGenericFramework(mol).ToSmiles(includeStereo: false);

    /// <summary>Atom indices of <paramref name="mol"/> that belong to the scaffold.</summary>
    public static IReadOnlyList<int> ScaffoldAtoms(Molecule mol, bool keepExocyclicDoubleBonds = true)
    {
        int n = mol.AtomCount;
        var rings = mol.Rings;
        if (rings.Count == 0) return [];
        var alive = new bool[n];
        var deg = new int[n];
        for (int i = 0; i < n; i++) { alive[i] = true; deg[i] = mol.Degree(i); }
        // prune terminal non-ring atoms repeatedly
        var stack = new Stack<int>();
        for (int i = 0; i < n; i++) if (!rings.IsAtomInRing(i) && deg[i] <= 1) stack.Push(i);
        while (stack.Count > 0)
        {
            int u = stack.Pop();
            if (!alive[u]) continue;
            alive[u] = false;
            foreach (int v in mol.Neighbors(u))
            {
                if (!alive[v]) continue;
                deg[v]--;
                if (!rings.IsAtomInRing(v) && deg[v] <= 1) stack.Push(v);
            }
        }
        var keep = new List<int>();
        for (int i = 0; i < n; i++) if (alive[i]) keep.Add(i);
        if (keepExocyclicDoubleBonds)
        {
            var core = keep.ToArray();
            foreach (int a in core)
            {
                var nb = mol.Neighbors(a);
                var nbb = mol.NeighborBonds(a);
                for (int k = 0; k < nb.Length; k++)
                    if (!alive[nb[k]] && mol.Bonds[nbb[k]].Order == BondOrder.Double && mol.Degree(nb[k]) == 1)
                    {
                        alive[nb[k]] = true;
                        keep.Add(nb[k]);
                    }
            }
            keep.Sort();
        }
        return keep;
    }

    /// <summary>
    /// A new molecule made of the <paramref name="atoms"/> of <paramref name="mol"/>; atoms that lose bonds receive
    /// hydrogens for the lost valence (so an N-substituted aromatic n becomes [nH]).
    /// </summary>
    public static Molecule Extract(Molecule mol, IReadOnlyList<int> atoms)
    {
        var keep = new HashSet<int>(atoms);
        var m = mol.Clone();
        m.Coordinates = null;
        for (int i = 0; i < m.AtomCount; i++)
        {
            if (!keep.Contains(i)) continue;
            int lost = 0;
            foreach (int b in mol.NeighborBonds(i))
                if (!keep.Contains(mol.Bonds[b].Other(i))) lost += mol.KekuleOrder(b);
            if (lost == 0) continue;
            ref var a = ref m.AtomAt(i);
            a.ExplicitHydrogens += a.ImplicitHydrogens + lost;
            a.ImplicitHydrogens = 0;
            a.NoImplicit = true;
        }
        m.RemoveAtoms(Enumerable.Range(0, mol.AtomCount).Where(i => !keep.Contains(i)));
        for (int i = 0; i < m.AtomCount; i++) m.SetChiralReference(i, null);
        for (int i = 0; i < m.AtomCount; i++) m.AtomAt(i).Chirality = ChiralTag.None;
        if (m.AtomCount > 0) Aromaticity.Perceive(m);
        return m;
    }
}

/// <summary>A parsed molecular formula (element counts and net charge) with mass calculations.</summary>
public sealed class MolecularFormula
{
    private readonly SortedDictionary<string, int> _counts;

    private MolecularFormula(SortedDictionary<string, int> counts, int charge)
    {
        _counts = counts;
        Charge = charge;
    }

    /// <summary>Element symbol → count (D and T are reported as separate pseudo-elements).</summary>
    public IReadOnlyDictionary<string, int> Counts => _counts;

    /// <summary>Net charge.</summary>
    public int Charge { get; }

    /// <summary>The formula of a molecule (implicit hydrogens included).</summary>
    public static MolecularFormula FromMolecule(Molecule mol) => Parse(MolecularDescriptors.Formula(mol));

    /// <summary>
    /// Parses formulas such as <c>C9H8O4</c>, <c>CH3COO-</c>, <c>Ca(OH)2</c>, <c>CuSO4·5H2O</c> or <c>C6H5D</c>.
    /// </summary>
    /// <exception cref="FormatException">The text is not a valid formula.</exception>
    public static MolecularFormula Parse(string formula)
    {
        ArgumentNullException.ThrowIfNull(formula);
        string s = formula.Replace(" ", "");
        int charge = 0;
        // trailing charge: "+", "-", "++", "+2", "-3"
        var m = System.Text.RegularExpressions.Regex.Match(s, @"(?:([+-])(\d+)|(\++|-+))$");
        if (m.Success)
        {
            if (m.Groups[1].Success)
                charge = (m.Groups[1].Value == "-" ? -1 : 1) * int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture);
            else
                charge = (m.Groups[3].Value[0] == '-' ? -1 : 1) * m.Groups[3].Value.Length;
            s = s[..m.Index];
        }
        var counts = new SortedDictionary<string, int>(StringComparer.Ordinal);
        foreach (var part in s.Split('·', '.', '*'))
        {
            if (part.Length == 0) continue;
            int pos = 0;
            int mult = ReadNumber(part, ref pos, 1);
            var c = ParseGroup(part, ref pos, formula);
            if (pos != part.Length) throw new FormatException($"Invalid formula '{formula}'.");
            foreach (var (k, v) in c) counts[k] = counts.GetValueOrDefault(k) + v * mult;
        }
        if (counts.Count == 0) throw new FormatException($"Empty formula '{formula}'.");
        return new MolecularFormula(counts, charge);
    }

    private static Dictionary<string, int> ParseGroup(string s, ref int pos, string original)
    {
        var counts = new Dictionary<string, int>();
        while (pos < s.Length && s[pos] != ')')
        {
            if (s[pos] is '(' or '[')
            {
                pos++;
                var inner = ParseGroup(s, ref pos, original);
                if (pos >= s.Length || s[pos] is not (')' or ']')) throw new FormatException($"Unbalanced parentheses in '{original}'.");
                pos++;
                int m = ReadNumber(s, ref pos, 1);
                foreach (var (k, v) in inner) counts[k] = counts.GetValueOrDefault(k) + v * m;
                continue;
            }
            if (s[pos] == ']') break;
            if (!char.IsAsciiLetterUpper(s[pos])) throw new FormatException($"Unexpected '{s[pos]}' in '{original}'.");
            string sym = s[pos].ToString();
            if (pos + 1 < s.Length && char.IsAsciiLetterLower(s[pos + 1]) && PeriodicTable.TryGet(s.Substring(pos, 2), out _)) sym = s.Substring(pos, 2);
            if (sym is not ("D" or "T") && !PeriodicTable.TryGet(sym, out _)) throw new FormatException($"Unknown element '{sym}' in '{original}'.");
            pos += sym.Length;
            int n = ReadNumber(s, ref pos, 1);
            counts[sym] = counts.GetValueOrDefault(sym) + n;
        }
        return counts;
    }

    private static int ReadNumber(string s, ref int pos, int fallback)
    {
        int start = pos;
        while (pos < s.Length && char.IsAsciiDigit(s[pos])) pos++;
        return pos == start ? fallback : int.Parse(s.AsSpan(start, pos - start), CultureInfo.InvariantCulture);
    }

    /// <summary>Average (standard atomic weight) mass in Da.</summary>
    public double AverageMass => _counts.Sum(kv => kv.Value * Mass(kv.Key, average: true));

    /// <summary>Monoisotopic mass in Da (most abundant isotopes), neutral.</summary>
    public double MonoisotopicMass => _counts.Sum(kv => kv.Value * Mass(kv.Key, average: false));

    /// <summary>m/z of the ion with the formula's charge (electron mass 0.000548580 Da accounted for).</summary>
    public double MassToCharge => Charge == 0 ? MonoisotopicMass : (MonoisotopicMass - Charge * 0.00054857990946) / Math.Abs(Charge);

    private static double Mass(string sym, bool average) => sym switch
    {
        "D" => PeriodicTable.IsotopeMass(1, 2),
        "T" => PeriodicTable.IsotopeMass(1, 3),
        _ => average ? PeriodicTable.Get(sym).AverageMass : PeriodicTable.Get(sym).MonoisotopicMass,
    };

    /// <summary>Hill-order formula string with charge suffix.</summary>
    public override string ToString()
    {
        var sb = new StringBuilder();
        var keys = _counts.Keys.ToList();
        void Emit(string k) { if (_counts.TryGetValue(k, out int c) && keys.Remove(k)) { sb.Append(k); if (c != 1) sb.Append(c); } }
        if (_counts.ContainsKey("C")) { Emit("C"); Emit("H"); }
        foreach (var k in keys.ToList()) Emit(k);
        if (Charge != 0) sb.Append(Charge > 0 ? '+' : '-').Append(Math.Abs(Charge) > 1 ? Math.Abs(Charge).ToString(CultureInfo.InvariantCulture) : "");
        return sb.ToString();
    }
}
