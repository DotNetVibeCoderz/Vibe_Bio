using System.Text;

namespace Biotech.Net.StructuralBiology;

/// <summary>DSSP secondary-structure classes (Kabsch &amp; Sander 1983).</summary>
public enum SecondaryStructure
{
    /// <summary>No regular structure (DSSP blank, shown as <c>-</c> or space).</summary>
    Loop,

    /// <summary>α-helix (<c>H</c>).</summary>
    AlphaHelix,

    /// <summary>Isolated β-bridge (<c>B</c>).</summary>
    Bridge,

    /// <summary>Extended strand in a ladder (<c>E</c>).</summary>
    Strand,

    /// <summary>3₁₀-helix (<c>G</c>).</summary>
    Helix310,

    /// <summary>π-helix (<c>I</c>).</summary>
    PiHelix,

    /// <summary>H-bonded turn (<c>T</c>).</summary>
    Turn,

    /// <summary>Bend (<c>S</c>).</summary>
    Bend,
}

/// <summary>Conversions for <see cref="SecondaryStructure"/>.</summary>
public static class SecondaryStructureExtensions
{
    /// <summary>The DSSP one-letter code (<c>' '</c> for loop).</summary>
    public static char ToCode(this SecondaryStructure s) => s switch
    {
        SecondaryStructure.AlphaHelix => 'H', SecondaryStructure.Bridge => 'B', SecondaryStructure.Strand => 'E',
        SecondaryStructure.Helix310 => 'G', SecondaryStructure.PiHelix => 'I', SecondaryStructure.Turn => 'T',
        SecondaryStructure.Bend => 'S', _ => ' ',
    };

    /// <summary>Three-state reduction: H/G/I → <c>H</c>, E/B → <c>E</c>, everything else → <c>C</c>.</summary>
    public static char ToThreeState(this SecondaryStructure s) => s switch
    {
        SecondaryStructure.AlphaHelix or SecondaryStructure.Helix310 or SecondaryStructure.PiHelix => 'H',
        SecondaryStructure.Strand or SecondaryStructure.Bridge => 'E',
        _ => 'C',
    };
}

/// <summary>A backbone hydrogen bond N–H(donor) ··· O=C(acceptor) with its DSSP electrostatic energy.</summary>
public readonly record struct BackboneHBond(Residue Donor, Residue Acceptor, double Energy);

/// <summary>Per-residue DSSP output.</summary>
public sealed class DsspResidue
{
    internal DsspResidue(Residue residue, int index) { Residue = residue; Index = index; }

    /// <summary>The residue.</summary>
    public Residue Residue { get; }

    /// <summary>Position in <see cref="DsspResult.Residues"/>.</summary>
    public int Index { get; }

    /// <summary>Assigned secondary structure.</summary>
    public SecondaryStructure Structure { get; internal set; }

    /// <summary>DSSP one-letter code.</summary>
    public char Code => Structure.ToCode();

    /// <summary>Virtual bond angle κ (CA i−2, i, i+2) in degrees, NaN at the ends; κ &gt; 70° is a bend.</summary>
    public double Kappa { get; internal set; } = double.NaN;

    /// <summary>The computed backbone amide hydrogen position.</summary>
    public Vec3 HydrogenPosition { get; internal set; }

    /// <summary>The two strongest acceptors of this residue's N–H, as (offset, energy) with offset = acceptor − this.</summary>
    public (int Offset, double Energy)[] NHToO { get; } = [(0, 0), (0, 0)];

    /// <summary>The two strongest donors to this residue's C=O, as (offset, energy).</summary>
    public (int Offset, double Energy)[] OToHN { get; } = [(0, 0), (0, 0)];

    /// <summary>Bridge partner indices (−1 when none), up to two.</summary>
    public int[] BridgePartners { get; } = [-1, -1];

    /// <summary>True when a chain break precedes this residue.</summary>
    public bool BreakBefore { get; internal set; }

    public override string ToString() => $"{Residue.Id} {Code}";
}

/// <summary>Result of <see cref="Dssp.Assign(Model, DsspOptions?)"/>.</summary>
public sealed class DsspResult
{
    internal DsspResult(List<DsspResidue> residues, List<BackboneHBond> hbonds)
    {
        Residues = residues;
        HydrogenBonds = hbonds;
    }

    /// <summary>Residues with complete backbones, in model order.</summary>
    public IReadOnlyList<DsspResidue> Residues { get; }

    /// <summary>Every backbone H-bond with energy below −0.5 kcal/mol.</summary>
    public IReadOnlyList<BackboneHBond> HydrogenBonds { get; }

    /// <summary>The DSSP string (one character per residue, loop as <paramref name="loopChar"/>).</summary>
    public string ToSecondaryStructureString(char loopChar = '-')
    {
        var sb = new StringBuilder(Residues.Count);
        foreach (var r in Residues) sb.Append(r.Structure == SecondaryStructure.Loop ? loopChar : r.Code);
        return sb.ToString();
    }

    /// <summary>Three-state (H/E/C) string.</summary>
    public string ToThreeStateString() => string.Concat(Residues.Select(r => r.Structure.ToThreeState()));

    /// <summary>The assignment for one residue, or null when it was not assigned (incomplete backbone).</summary>
    public DsspResidue? this[Residue residue] => Residues.FirstOrDefault(r => ReferenceEquals(r.Residue, residue));

    /// <summary>Counts per class.</summary>
    public Dictionary<SecondaryStructure, int> Composition() =>
        Residues.GroupBy(r => r.Structure).ToDictionary(g => g.Key, g => g.Count());

    public override string ToString() => ToSecondaryStructureString();
}

/// <summary>Options for <see cref="Dssp"/>.</summary>
public sealed record DsspOptions
{
    /// <summary>π-helices take precedence over α-helices (DSSP ≥ 2.1 default, true).</summary>
    public bool PreferPiHelices { get; init; } = true;
}

/// <summary>
/// DSSP secondary-structure assignment (Kabsch &amp; Sander 1983, Biopolymers 22:2577), following the reference
/// implementation (DSSP 2.x/4.x): electrostatic H-bond energy
/// E = 0.084·332·(1/r(ON) + 1/r(CH) − 1/r(OH) − 1/r(CN)) kcal/mol with a −0.5 cut-off, the two best partners per
/// N–H and C=O, n-turns (n = 3, 4, 5), minimal helices, bridges/ladders with β-bulges, and bends (κ &gt; 70°).
/// The polyproline (P) class of DSSP 4 is not assigned.
/// </summary>
public static class Dssp
{
    private const double Coupling = -27.888; // −0.084 × 332
    private const double MinEnergy = -9.9;
    private const double MaxHBondEnergy = -0.5;
    private const double MinCaDistance = 9.0;

    /// <summary>
    /// DSSP electrostatic hydrogen-bond energy (kcal/mol) between donor N/H and acceptor C=O, rounded to 0.001
    /// as DSSP does and floored at −9.9.
    /// </summary>
    public static double HBondEnergy(Vec3 donorN, Vec3 donorH, Vec3 acceptorC, Vec3 acceptorO)
    {
        double dHO = Vec3.Distance(donorH, acceptorO), dHC = Vec3.Distance(donorH, acceptorC);
        double dNC = Vec3.Distance(donorN, acceptorC), dNO = Vec3.Distance(donorN, acceptorO);
        if (dHO < 0.5 || dHC < 0.5 || dNC < 0.5 || dNO < 0.5) return MinEnergy;
        double e = Coupling / dHO - Coupling / dHC + Coupling / dNC - Coupling / dNO;
        e = Math.Round(e * 1000) / 1000;
        return Math.Max(e, MinEnergy);
    }

    /// <summary>Assigns secondary structure to the first model of a structure.</summary>
    public static DsspResult Assign(Structure structure, DsspOptions? options = null) => Assign(structure.FirstModel, options);

    /// <summary>Assigns secondary structure to every amino acid with N, CA, C and O in the model.</summary>
    public static DsspResult Assign(Model model, DsspOptions? options = null)
    {
        options ??= new DsspOptions();
        var res = new List<DsspResidue>();
        foreach (var chain in model.Chains)
            foreach (var r in chain.Residues)
                if (r.IsAminoAcid && r.HasAtom("N") && r.HasAtom("CA") && r.HasAtom("C") && r.HasAtom("O"))
                    res.Add(new DsspResidue(r, res.Count));
        int n = res.Count;
        var N = new Vec3[n]; var CA = new Vec3[n]; var C = new Vec3[n]; var O = new Vec3[n]; var H = new Vec3[n];
        var isPro = new bool[n];
        for (int i = 0; i < n; i++)
        {
            var r = res[i].Residue;
            N[i] = r["N"].Position; CA[i] = r["CA"].Position; C[i] = r["C"].Position; O[i] = r["O"].Position;
            isPro[i] = r.Name == "PRO";
        }
        // Chain breaks: different chain, or C(i−1)–N(i) longer than 2.5 Å.
        var breakBefore = new bool[n];
        for (int i = 0; i < n; i++)
        {
            breakBefore[i] = i == 0 || !ReferenceEquals(res[i].Residue.Chain, res[i - 1].Residue.Chain) || Vec3.Distance(C[i - 1], N[i]) > 2.5;
            res[i].BreakBefore = i > 0 && breakBefore[i];
            H[i] = N[i];
            if (i > 0 && !breakBefore[i]) H[i] = N[i] + (C[i - 1] - O[i - 1]).Normalized();
            res[i].HydrogenPosition = H[i];
        }
        // Prefix count of breaks for O(1) NoChainBreak(i, j).
        var breaks = new int[n + 1];
        for (int i = 0; i < n; i++) breaks[i + 1] = breaks[i] + (i > 0 && breakBefore[i] ? 1 : 0);
        bool NoBreak(int i, int j)
        {
            if (i > j) (i, j) = (j, i);
            if (i < 0 || j >= n) return false;
            return breaks[j + 1] - breaks[i + 1] == 0;
        }

        // H-bond energies between residues with CA closer than 9 Å (cell list).
        var acc = new (int Partner, double E)[n, 2];   // best acceptors of residue i's NH
        var don = new (int Partner, double E)[n, 2];   // best donors to residue i's CO
        for (int i = 0; i < n; i++) { acc[i, 0] = acc[i, 1] = don[i, 0] = don[i, 1] = (-1, 0); }
        var hbonds = new List<BackboneHBond>();
        var ns = new NeighborSearch(CA, MinCaDistance);
        foreach (var (i, j) in ns.AllPairs(MinCaDistance).OrderBy(p => p.I).ThenBy(p => p.J))
        {
            if (Vec3.Distance(CA[i], CA[j]) >= MinCaDistance) continue;
            Calc(i, j);
            if (j != i + 1) Calc(j, i);
        }
        void Calc(int d, int a)
        {
            if (isPro[d]) return;
            double e = HBondEnergy(N[d], H[d], C[a], O[a]);
            if (e < MaxHBondEnergy) hbonds.Add(new BackboneHBond(res[d].Residue, res[a].Residue, e));
            if (e < acc[d, 0].E) { acc[d, 1] = acc[d, 0]; acc[d, 0] = (a, e); }
            else if (e < acc[d, 1].E) acc[d, 1] = (a, e);
            if (e < don[a, 0].E) { don[a, 1] = don[a, 0]; don[a, 0] = (d, e); }
            else if (e < don[a, 1].E) don[a, 1] = (d, e);
        }
        for (int i = 0; i < n; i++)
            for (int k = 0; k < 2; k++)
            {
                res[i].NHToO[k] = acc[i, k].Partner < 0 ? (0, 0) : (acc[i, k].Partner - i, acc[i, k].E);
                res[i].OToHN[k] = don[i, k].Partner < 0 ? (0, 0) : (don[i, k].Partner - i, don[i, k].E);
            }

        // TestBond(a, b): N-H of a is H-bonded to C=O of b (among a's two best acceptors).
        bool Bond(int a, int b) =>
            a >= 0 && b >= 0 && a < n && b < n &&
            ((acc[a, 0].Partner == b && acc[a, 0].E < MaxHBondEnergy) || (acc[a, 1].Partner == b && acc[a, 1].E < MaxHBondEnergy));

        var ss = new SecondaryStructure[n];

        // ---- β-bridges and ladders ----
        var ladders = new List<Ladder>();
        for (int i = 1; i + 4 < n; i++)
            for (int j = i + 3; j + 1 < n; j++)
            {
                var type = TestBridge(i, j);
                if (type == 0) continue;
                bool found = false;
                foreach (var l in ladders)
                {
                    if (l.Type != type || i != l.I[^1] + 1) continue;
                    if (type == 1 && l.J[^1] + 1 == j) { l.I.Add(i); l.J.Add(j); l.Pairs.Add((i, j)); found = true; break; }
                    if (type == 2 && l.J[0] - 1 == j) { l.I.Add(i); l.J.Insert(0, j); l.Pairs.Add((i, j)); found = true; break; }
                }
                if (!found) ladders.Add(new Ladder(type, i, j));
            }
        int TestBridge(int i, int j)
        {
            int a = i - 1, b = i, c = i + 1, d = j - 1, e = j, f = j + 1;
            if (!NoBreak(a, c) || !NoBreak(d, f)) return 0;
            if ((Bond(c, e) && Bond(e, a)) || (Bond(f, b) && Bond(b, d))) return 1;          // parallel
            if ((Bond(c, d) && Bond(f, a)) || (Bond(e, b) && Bond(b, e))) return 2;          // antiparallel
            return 0;
        }
        // β-bulges: merge ladders of the same type separated by small gaps.
        for (int x = 0; x < ladders.Count; x++)
            for (int y = x + 1; y < ladders.Count; y++)
            {
                var li = ladders[x]; var lj = ladders[y];
                int ibi = li.I[0], iei = li.I[^1], jbi = li.J[0], jei = li.J[^1];
                int ibj = lj.I[0], iej = lj.I[^1], jbj = lj.J[0], jej = lj.J[^1];
                if (li.Type != lj.Type || !NoBreak(Math.Min(ibi, ibj), Math.Max(iei, iej)) || !NoBreak(Math.Min(jbi, jbj), Math.Max(jei, jej))
                    || ibj - iei >= 6 || (iei >= ibj && ibi <= iej)) continue;
                bool bulge = li.Type == 1
                    ? (jbj - jei < 6 && ibj - iei < 3) || jbj - jei < 3
                    : (jbi - jej < 6 && ibj - iei < 3) || jbi - jej < 3;
                if (!bulge) continue;
                li.I.AddRange(lj.I);
                if (li.Type == 1) li.J.AddRange(lj.J); else li.J.InsertRange(0, lj.J);
                foreach (var p in lj.Pairs) li.Pairs.Add(p);
                ladders.RemoveAt(y);
                y--;
            }
        foreach (var l in ladders)
        {
            var type = l.I.Count > 1 ? SecondaryStructure.Strand : SecondaryStructure.Bridge;
            for (int i = l.I.Min(); i <= l.I.Max(); i++) if (ss[i] != SecondaryStructure.Strand) ss[i] = type;
            for (int j = l.J.Min(); j <= l.J.Max(); j++) if (ss[j] != SecondaryStructure.Strand) ss[j] = type;
            foreach (var (pi, pj) in l.Pairs) { AddPartner(pi, pj); AddPartner(pj, pi); }
        }
        void AddPartner(int i, int j)
        {
            var bp = res[i].BridgePartners;
            if (bp[0] == j || bp[1] == j) return;
            if (bp[0] < 0) bp[0] = j; else if (bp[1] < 0) bp[1] = j;
        }

        // ---- n-turns and helices ----
        var start = new bool[3, n]; // stride 3,4,5 → index 0,1,2: residue i starts an n-turn
        for (int s = 3; s <= 5; s++)
            for (int i = 0; i + s < n; i++)
                if (NoBreak(i, i + s) && Bond(i + s, i)) start[s - 3, i] = true;
        for (int i = 1; i + 4 < n; i++)
            if (start[1, i] && start[1, i - 1])
                for (int j = i; j <= i + 3; j++) ss[j] = SecondaryStructure.AlphaHelix;
        for (int i = 1; i + 3 < n; i++)
            if (start[0, i] && start[0, i - 1])
            {
                bool empty = true;
                for (int j = i; empty && j <= i + 2; j++) empty = ss[j] is SecondaryStructure.Loop or SecondaryStructure.Helix310;
                if (empty) for (int j = i; j <= i + 2; j++) ss[j] = SecondaryStructure.Helix310;
            }
        for (int i = 1; i + 5 < n; i++)
            if (start[2, i] && start[2, i - 1])
            {
                bool empty = true;
                for (int j = i; empty && j <= i + 4; j++)
                    empty = ss[j] is SecondaryStructure.Loop or SecondaryStructure.PiHelix
                        || (options.PreferPiHelices && ss[j] == SecondaryStructure.AlphaHelix);
                if (empty) for (int j = i; j <= i + 4; j++) ss[j] = SecondaryStructure.PiHelix;
            }

        // ---- bends (κ) and turns ----
        for (int i = 2; i + 2 < n; i++)
            if (NoBreak(i - 2, i + 2))
                res[i].Kappa = Vec3.AngleBetween(CA[i] - CA[i - 2], CA[i + 2] - CA[i]);
        for (int i = 1; i + 1 < n; i++)
        {
            if (ss[i] != SecondaryStructure.Loop) continue;
            bool turn = false;
            for (int s = 3; s <= 5 && !turn; s++)
                for (int k = 1; k < s && !turn; k++)
                    turn = i >= k && start[s - 3, i - k];
            if (turn) ss[i] = SecondaryStructure.Turn;
            else if (res[i].Kappa > 70) ss[i] = SecondaryStructure.Bend;
        }
        for (int i = 0; i < n; i++) res[i].Structure = ss[i];
        return new DsspResult(res, hbonds);
    }

    private sealed class Ladder
    {
        public Ladder(int type, int i, int j) { Type = type; I.Add(i); J.Add(j); Pairs.Add((i, j)); }
        public int Type { get; }          // 1 parallel, 2 antiparallel
        public List<int> I { get; } = [];
        public List<int> J { get; } = [];
        public List<(int, int)> Pairs { get; } = [];
    }
}
