using System.Globalization;
using System.Text;

namespace Biotech.Net.Chemistry;

/// <summary>
/// Wildman–Crippen atom-contribution logP and molar refractivity (J. Chem. Inf. Comput. Sci. 1999, 39, 868).
/// All 110 published atom-type SMARTS (C1–C27, CS, H1–H4, HS, N1–N14, NS, O1–O12, OS, halogens, P, S1–S3,
/// Me1, Me2) are implemented, so values reproduce RDKit's <c>MolLogP</c>/<c>MolMR</c>.
/// </summary>
public static class Crippen
{
    private static readonly Lazy<(string Type, SmartsPattern Pattern, double LogP, double MR)[]> _types = new(() =>
        CrippenData.Types.Select(t => (t.Type, SmartsPattern.Parse(t.Smarts, mergeHydrogens: false), t.LogP, t.MR)).ToArray());

    /// <summary>Wildman–Crippen logP (a.k.a. WLOGP/ALOGP-style).</summary>
    public static double LogP(Molecule mol) => Compute(mol).LogP;

    /// <summary>Wildman–Crippen molar refractivity.</summary>
    public static double MolarRefractivity(Molecule mol) => Compute(mol).MR;

    /// <summary>Both values plus the type assigned to every atom of the hydrogen-explicit molecule.</summary>
    public static (double LogP, double MR, string[] AtomTypes) Compute(Molecule mol)
    {
        var h = mol.AddHydrogens();
        int n = h.AtomCount;
        var assigned = new string?[n];
        double logp = 0, mr = 0;
        var types = _types.Value;
        var ctx = new SmartsPattern.MatchContext(h);
        foreach (var (type, pattern, lp, r) in types)
        {
            for (int i = 0; i < n; i++)
            {
                if (assigned[i] is not null) continue;
                if (pattern.FindMatches(h, false, 1, ctx, i).Count == 0) continue;
                assigned[i] = type;
                logp += lp;
                mr += r;
            }
        }
        return (logp, mr, assigned.Select(a => a ?? "?").ToArray());
    }
}

/// <summary>
/// Topological polar surface area (Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714) from N and O
/// fragment contributions (S and P are excluded, as in the original paper and RDKit's default).
/// </summary>
public static class Tpsa
{
    /// <summary>TPSA in Å².</summary>
    public static double Compute(Molecule mol)
    {
        double total = 0;
        for (int i = 0; i < mol.AtomCount; i++) total += Contribution(mol, i);
        return total;
    }

    /// <summary>The Ertl contribution of one atom.</summary>
    public static double Contribution(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        if (a.AtomicNumber is not (7 or 8)) return 0;
        int nSing = 0, nDoub = 0, nTrip = 0, nArom = 0, nNbrs = 0;
        var nb = mol.Neighbors(atom);
        var nbb = mol.NeighborBonds(atom);
        for (int k = 0; k < nb.Length; k++)
        {
            if (mol.Atoms[nb[k]].AtomicNumber == 1) continue;
            nNbrs++;
            switch (mol.Bonds[nbb[k]].Order)
            {
                case BondOrder.Single: nSing++; break;
                case BondOrder.Double: nDoub++; break;
                case BondOrder.Triple: nTrip++; break;
                case BondOrder.Aromatic: nArom++; break;
            }
        }
        int nHs = mol.TotalHydrogens(atom);
        int chg = a.FormalCharge;
        bool in3 = mol.Rings.IsAtomInRingOfSize(atom, 3);
        double t = -1;
        if (a.AtomicNumber == 7)
        {
            switch (nNbrs)
            {
                case 1:
                    if (nHs == 0 && chg == 0 && nTrip == 1) t = 23.79;
                    else if (nHs == 1 && chg == 0 && nDoub == 1) t = 23.85;
                    else if (nHs == 2 && chg == 0 && nSing == 1) t = 26.02;
                    else if (nHs == 2 && chg == 1 && nDoub == 1) t = 25.59;
                    else if (nHs == 3 && chg == 1 && nSing == 1) t = 27.64;
                    break;
                case 2:
                    if (nHs == 0 && chg == 0 && nSing == 1 && nDoub == 1) t = 12.36;
                    else if (nHs == 0 && chg == 0 && nTrip == 1 && nDoub == 1) t = 13.60;
                    else if (nHs == 1 && chg == 0 && nSing == 2 && in3) t = 21.94;
                    else if (nHs == 1 && chg == 0 && nSing == 2 && !in3) t = 12.03;
                    else if (nHs == 0 && chg == 1 && nTrip == 1 && nSing == 1) t = 4.36;
                    else if (nHs == 1 && chg == 1 && nDoub == 1 && nSing == 1) t = 13.97;
                    else if (nHs == 2 && chg == 1 && nSing == 2) t = 16.61;
                    else if (nHs == 0 && chg == 0 && nArom == 2) t = 12.89;
                    else if (nHs == 1 && chg == 0 && nArom == 2) t = 15.79;
                    else if (nHs == 1 && chg == 1 && nArom == 2) t = 14.14;
                    break;
                case 3:
                    if (nHs == 0 && chg == 0 && nSing == 3 && in3) t = 3.01;
                    else if (nHs == 0 && chg == 0 && nSing == 3 && !in3) t = 3.24;
                    else if (nHs == 0 && chg == 0 && nSing == 1 && nDoub == 2) t = 11.68;
                    else if (nHs == 0 && chg == 1 && nSing == 2 && nDoub == 1) t = 3.01;
                    else if (nHs == 1 && chg == 1 && nSing == 3) t = 4.44;
                    else if (nHs == 0 && chg == 0 && nArom == 3) t = 4.41;
                    else if (nHs == 0 && chg == 0 && nSing == 1 && nArom == 2) t = 4.93;
                    else if (nHs == 0 && chg == 0 && nDoub == 1 && nArom == 2) t = 8.39;
                    else if (nHs == 0 && chg == 1 && nArom == 3) t = 4.10;
                    else if (nHs == 0 && chg == 1 && nSing == 1 && nArom == 2) t = 3.88;
                    break;
                case 4:
                    if (nHs == 0 && nSing == 4 && chg == 1) t = 0.00;
                    break;
            }
            if (t < 0) t = Math.Max(0, 30.5 - nNbrs * 8.2 + nHs * 1.5);
        }
        else
        {
            switch (nNbrs)
            {
                case 1:
                    if (nHs == 0 && chg == 0 && nDoub == 1) t = 17.07;
                    else if (nHs == 1 && chg == 0 && nSing == 1) t = 20.23;
                    else if (nHs == 0 && chg == -1 && nSing == 1) t = 23.06;
                    break;
                case 2:
                    if (nHs == 0 && chg == 0 && nSing == 2 && in3) t = 12.53;
                    else if (nHs == 0 && chg == 0 && nSing == 2 && !in3) t = 9.23;
                    else if (nHs == 0 && chg == 0 && nArom == 2) t = 13.14;
                    break;
            }
            if (t < 0) t = Math.Max(0, 28.5 - nNbrs * 8.6 + nHs * 1.5);
        }
        return t;
    }
}

/// <summary>A table of standard 2-D descriptors for one molecule.</summary>
public sealed record DescriptorSet
{
    public required string Formula { get; init; }
    public required double MolecularWeight { get; init; }
    public required double ExactMass { get; init; }
    public required int HeavyAtoms { get; init; }
    public required int Heteroatoms { get; init; }
    /// <summary>Lipinski donors: number of N–H and O–H hydrogens (RDKit <c>NHOHCount</c>).</summary>
    public required int LipinskiHbd { get; init; }
    /// <summary>Lipinski acceptors: number of N and O atoms (RDKit <c>NOCount</c>).</summary>
    public required int LipinskiHba { get; init; }
    public required int HBondDonors { get; init; }
    public required int HBondAcceptors { get; init; }
    public required int RotatableBonds { get; init; }
    public required double Tpsa { get; init; }
    public required double LogP { get; init; }
    public required double MolarRefractivity { get; init; }
    public required int RingCount { get; init; }
    public required int AromaticRingCount { get; init; }
    public required double FractionCsp3 { get; init; }
    public required int FormalCharge { get; init; }
    public required int TotalAtoms { get; init; }

    /// <summary>Descriptor values in a fixed order (see <see cref="Names"/>), for feature matrices.</summary>
    public double[] ToVector() =>
    [
        MolecularWeight, HeavyAtoms, Heteroatoms, LipinskiHbd, LipinskiHba, HBondDonors, HBondAcceptors,
        RotatableBonds, Tpsa, LogP, MolarRefractivity, RingCount, AromaticRingCount, FractionCsp3, FormalCharge,
    ];

    /// <summary>Names matching <see cref="ToVector"/>.</summary>
    public static IReadOnlyList<string> Names { get; } =
    [
        "MW", "HeavyAtoms", "Heteroatoms", "NHOH", "NO", "HBD", "HBA", "RotB", "TPSA", "LogP", "MR", "Rings",
        "AromaticRings", "FractionCSP3", "Charge",
    ];
}

/// <summary>
/// Molecular descriptors. Definitions follow RDKit so results can be compared directly:
/// H-bond donors/acceptors use the Lipinski module SMARTS, rotatable bonds the "strict" definition
/// (excluding amide/ester C–X bonds and CX3/t-butyl rotors), aromatic rings are SSSR rings whose bonds are all aromatic.
/// </summary>
public static class MolecularDescriptors
{
    private static readonly Lazy<SmartsPattern> _hbd = new(() => SmartsPattern.Parse("[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]"));
    private static readonly Lazy<SmartsPattern> _hba = new(() => SmartsPattern.Parse(
        "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&X2&H0&+0,$([o,s;+0])]"));
    private static readonly Lazy<SmartsPattern> _rotStrict = new(() => SmartsPattern.Parse(
        "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]"));
    private static readonly Lazy<SmartsPattern> _rotSimple = new(() => SmartsPattern.Parse("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"));

    /// <summary>
    /// Molecular formula in Hill order (C, H, then alphabetical; alphabetical when there is no carbon), with net charge
    /// suffix. Deuterium/tritium are counted as H unless <paramref name="separateIsotopes"/> (RDKit <c>CalcMolFormula</c> default).
    /// </summary>
    public static string Formula(Molecule mol, bool separateIsotopes = false)
    {
        var counts = new SortedDictionary<string, int>(StringComparer.Ordinal);
        int charge = 0;
        for (int i = 0; i < mol.AtomCount; i++)
        {
            var a = mol.Atoms[i];
            charge += a.FormalCharge;
            if (a.AtomicNumber == 0) continue;
            string sym = a.Symbol;
            if (separateIsotopes && a.AtomicNumber == 1 && a.Isotope == 2) sym = "D";
            else if (separateIsotopes && a.AtomicNumber == 1 && a.Isotope == 3) sym = "T";
            counts[sym] = counts.GetValueOrDefault(sym) + 1;
            if (a.HydrogenCount > 0) counts["H"] = counts.GetValueOrDefault("H") + a.HydrogenCount;
        }
        var sb = new StringBuilder();
        void Emit(string s) { if (counts.Remove(s, out int c)) { sb.Append(s); if (c > 1) sb.Append(c); } }
        if (counts.ContainsKey("C")) { Emit("C"); Emit("H"); }
        foreach (var s in counts.Keys.ToList()) Emit(s);
        if (charge != 0) sb.Append(charge > 0 ? '+' : '-').Append(Math.Abs(charge) > 1 ? Math.Abs(charge).ToString(CultureInfo.InvariantCulture) : "");
        return sb.ToString();
    }

    /// <summary>Average molecular weight (standard atomic weights; labelled atoms use isotope masses).</summary>
    public static double MolecularWeight(Molecule mol)
    {
        double m = 0, h = PeriodicTable.Get(1).AverageMass;
        foreach (ref readonly var a in mol.Atoms)
        {
            if (a.AtomicNumber == 0) continue;
            m += a.Isotope != 0 ? PeriodicTable.IsotopeMass(a.AtomicNumber, a.Isotope) : a.Element.AverageMass;
            m += a.HydrogenCount * h;
        }
        return m;
    }

    /// <summary>Monoisotopic (exact) mass.</summary>
    public static double ExactMass(Molecule mol)
    {
        double m = 0, h = PeriodicTable.Get(1).MonoisotopicMass;
        foreach (ref readonly var a in mol.Atoms)
        {
            if (a.AtomicNumber == 0) continue;
            m += a.Isotope != 0 ? PeriodicTable.IsotopeMass(a.AtomicNumber, a.Isotope) : a.Element.MonoisotopicMass;
            m += a.HydrogenCount * h;
        }
        return m;
    }

    /// <summary>Heavy atoms (not H, not dummy).</summary>
    public static int HeavyAtomCount(Molecule mol)
    {
        int n = 0;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber > 1) n++;
        return n;
    }

    /// <summary>Atoms other than C and H.</summary>
    public static int HeteroatomCount(Molecule mol)
    {
        int n = 0;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber is > 1 and not 6) n++;
        return n;
    }

    /// <summary>Lipinski donor count: hydrogens on N and O (each H counts).</summary>
    public static int LipinskiDonors(Molecule mol)
    {
        int n = 0;
        for (int i = 0; i < mol.AtomCount; i++)
            if (mol.Atoms[i].AtomicNumber is 7 or 8) n += mol.TotalHydrogens(i);
        return n;
    }

    /// <summary>Lipinski acceptor count: N + O atoms.</summary>
    public static int LipinskiAcceptors(Molecule mol)
    {
        int n = 0;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber is 7 or 8) n++;
        return n;
    }

    /// <summary>H-bond donors (RDKit Lipinski.NumHDonors definition).</summary>
    public static int HBondDonors(Molecule mol) => _hbd.Value.CountMatches(mol);

    /// <summary>H-bond acceptors (RDKit <c>CalcNumHBA</c> / <c>Lipinski.NumHAcceptors</c> definition; aromatic n only when two-connected).</summary>
    public static int HBondAcceptors(Molecule mol) => _hba.Value.CountMatches(mol);

    /// <summary>Rotatable bonds; <paramref name="strict"/> (default) excludes amide/ester bonds and CF3/CCl3/CBr3/t-butyl rotors.</summary>
    public static int RotatableBonds(Molecule mol, bool strict = true) => (strict ? _rotStrict : _rotSimple).Value.CountMatches(mol);

    /// <summary>Number of SSSR rings.</summary>
    public static int RingCount(Molecule mol) => mol.Rings.Count;

    /// <summary>SSSR rings whose bonds are all aromatic.</summary>
    public static int AromaticRingCount(Molecule mol)
    {
        int n = 0;
        foreach (var r in mol.Rings.BondRings) if (r.All(b => mol.Bonds[b].IsAromatic)) n++;
        return n;
    }

    /// <summary>Fraction of carbons that are sp3 (only single bonds, not aromatic).</summary>
    public static double FractionCsp3(Molecule mol)
    {
        int c = 0, sp3 = 0;
        for (int i = 0; i < mol.AtomCount; i++)
        {
            if (mol.Atoms[i].AtomicNumber != 6) continue;
            c++;
            bool sat = !mol.Atoms[i].IsAromatic;
            foreach (int b in mol.NeighborBonds(i)) if (mol.Bonds[b].Order != BondOrder.Single) { sat = false; break; }
            if (sat) sp3++;
        }
        return c == 0 ? 0 : (double)sp3 / c;
    }

    /// <summary>Net formal charge.</summary>
    public static int NetCharge(Molecule mol)
    {
        int q = 0;
        foreach (ref readonly var a in mol.Atoms) q += a.FormalCharge;
        return q;
    }

    /// <summary>Total atoms including hydrogens.</summary>
    public static int TotalAtomCount(Molecule mol)
    {
        int n = 0;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber > 0) n += 1 + a.HydrogenCount;
        return n;
    }

    /// <summary>All descriptors at once.</summary>
    public static DescriptorSet Compute(Molecule mol)
    {
        var (logp, mr, _) = Crippen.Compute(mol);
        return new DescriptorSet
        {
            Formula = Formula(mol),
            MolecularWeight = MolecularWeight(mol),
            ExactMass = ExactMass(mol),
            HeavyAtoms = HeavyAtomCount(mol),
            Heteroatoms = HeteroatomCount(mol),
            LipinskiHbd = LipinskiDonors(mol),
            LipinskiHba = LipinskiAcceptors(mol),
            HBondDonors = HBondDonors(mol),
            HBondAcceptors = HBondAcceptors(mol),
            RotatableBonds = RotatableBonds(mol),
            Tpsa = Tpsa.Compute(mol),
            LogP = logp,
            MolarRefractivity = mr,
            RingCount = RingCount(mol),
            AromaticRingCount = AromaticRingCount(mol),
            FractionCsp3 = FractionCsp3(mol),
            FormalCharge = NetCharge(mol),
            TotalAtoms = TotalAtomCount(mol),
        };
    }
}

/// <summary>Outcome of a rule-based filter: pass/fail plus which criteria were violated.</summary>
public sealed record RuleResult(string Rule, bool Passed, IReadOnlyList<string> Violations)
{
    public override string ToString() => Passed ? $"{Rule}: pass" : $"{Rule}: fail ({string.Join(", ", Violations)})";
}

/// <summary>Drug-likeness and lead-likeness rules.</summary>
public static class DrugLikeness
{
    /// <summary>
    /// Lipinski's rule of five (Adv. Drug Deliv. Rev. 1997): MW ≤ 500, logP ≤ 5, NH/OH ≤ 5, N+O ≤ 10. Passes with at most
    /// <paramref name="allowedViolations"/> violations (1 by the original formulation).
    /// </summary>
    public static RuleResult Lipinski(DescriptorSet d, int allowedViolations = 1)
    {
        var v = new List<string>();
        if (d.MolecularWeight > 500) v.Add($"MW {d.MolecularWeight:F1} > 500");
        if (d.LogP > 5) v.Add($"logP {d.LogP:F2} > 5");
        if (d.LipinskiHbd > 5) v.Add($"HBD {d.LipinskiHbd} > 5");
        if (d.LipinskiHba > 10) v.Add($"HBA {d.LipinskiHba} > 10");
        return new("Lipinski", v.Count <= allowedViolations, v);
    }

    /// <summary>Veber (J. Med. Chem. 2002): rotatable bonds ≤ 10 and TPSA ≤ 140 Å².</summary>
    public static RuleResult Veber(DescriptorSet d)
    {
        var v = new List<string>();
        if (d.RotatableBonds > 10) v.Add($"RotB {d.RotatableBonds} > 10");
        if (d.Tpsa > 140) v.Add($"TPSA {d.Tpsa:F1} > 140");
        return new("Veber", v.Count == 0, v);
    }

    /// <summary>Ghose (J. Comb. Chem. 1999): 160 ≤ MW ≤ 480, −0.4 ≤ logP ≤ 5.6, 40 ≤ MR ≤ 130, 20 ≤ atoms ≤ 70.</summary>
    public static RuleResult Ghose(DescriptorSet d)
    {
        var v = new List<string>();
        if (d.MolecularWeight is < 160 or > 480) v.Add($"MW {d.MolecularWeight:F1} outside 160–480");
        if (d.LogP is < -0.4 or > 5.6) v.Add($"logP {d.LogP:F2} outside −0.4–5.6");
        if (d.MolarRefractivity is < 40 or > 130) v.Add($"MR {d.MolarRefractivity:F1} outside 40–130");
        if (d.TotalAtoms is < 20 or > 70) v.Add($"atoms {d.TotalAtoms} outside 20–70");
        return new("Ghose", v.Count == 0, v);
    }

    /// <summary>Egan (J. Med. Chem. 2000): TPSA ≤ 131.6 and WLOGP ≤ 5.88.</summary>
    public static RuleResult Egan(DescriptorSet d)
    {
        var v = new List<string>();
        if (d.Tpsa > 131.6) v.Add($"TPSA {d.Tpsa:F1} > 131.6");
        if (d.LogP > 5.88) v.Add($"logP {d.LogP:F2} > 5.88");
        return new("Egan", v.Count == 0, v);
    }

    /// <summary>
    /// Muegge (J. Med. Chem. 2001): 200 ≤ MW ≤ 600, −2 ≤ logP ≤ 5, TPSA ≤ 150, rings ≤ 7, C &gt; 4, heteroatoms &gt; 1,
    /// rotatable bonds ≤ 15, HBA ≤ 10, HBD ≤ 5.
    /// </summary>
    public static RuleResult Muegge(Molecule mol, DescriptorSet d)
    {
        var v = new List<string>();
        int carbons = 0;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber == 6) carbons++;
        if (d.MolecularWeight is < 200 or > 600) v.Add("MW");
        if (d.LogP is < -2 or > 5) v.Add("logP");
        if (d.Tpsa > 150) v.Add("TPSA");
        if (d.RingCount > 7) v.Add("rings");
        if (carbons <= 4) v.Add("carbons");
        if (d.Heteroatoms <= 1) v.Add("heteroatoms");
        if (d.RotatableBonds > 15) v.Add("RotB");
        if (d.HBondAcceptors > 10) v.Add("HBA");
        if (d.HBondDonors > 5) v.Add("HBD");
        return new("Muegge", v.Count == 0, v);
    }

    /// <summary>Lead-likeness (Teague/Oprea; SwissADME form): 250 ≤ MW ≤ 350, logP ≤ 3.5, rotatable bonds ≤ 7.</summary>
    public static RuleResult LeadLike(DescriptorSet d)
    {
        var v = new List<string>();
        if (d.MolecularWeight is < 250 or > 350) v.Add($"MW {d.MolecularWeight:F1} outside 250–350");
        if (d.LogP > 3.5) v.Add($"logP {d.LogP:F2} > 3.5");
        if (d.RotatableBonds > 7) v.Add($"RotB {d.RotatableBonds} > 7");
        return new("LeadLike", v.Count == 0, v);
    }
}

/// <summary>The eight QED properties (Bickerton et al., Nat. Chem. 2012).</summary>
public sealed record QedProperties(double MW, double ALogP, int HBA, int HBD, double PSA, int RotB, int Arom, int Alerts);

/// <summary>
/// Quantitative Estimate of Drug-likeness (Bickerton, Paolini, Besnard, Muresan, Hopkins, Nat. Chem. 2012, 4, 90),
/// using the published asymmetric-double-sigmoid parameters, mean weights and the 116 unwanted-group alerts,
/// with property definitions matching RDKit's <c>QED</c> module.
/// </summary>
public static class Qed
{
    private static readonly (double A, double B, double C, double D, double E, double F, double DMax)[] Ads =
    [
        (2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561),
        (3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604),
        (2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046),
        (1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616),
        (1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167),
        (0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403),
        (3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610),
        (0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140),
    ];

    /// <summary>Mean weights (QEDw,mo).</summary>
    public static IReadOnlyList<double> MeanWeights { get; } = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];

    /// <summary>Max weights (QEDw,max).</summary>
    public static IReadOnlyList<double> MaxWeights { get; } = [0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00];

    private static readonly string[] AcceptorSmarts =
    [
        "[oH0;X2]", "[OH1;X2;v2]", "[OH0;X2;v2]", "[OH0;X1;v2]", "[O-;X1]", "[SH0;X2;v2]", "[SH0;X1;v2]",
        "[S-;X1]", "[nH0;X2]", "[NH0;X1;v3]", "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
    ];

    /// <summary>The 116 structural alerts of the QED paper (as SMARTS).</summary>
    public static IReadOnlyList<string> AlertSmarts { get; } =
    [
        "*1[O,S,N]*1", "[S,C](=[O,S])[F,Br,Cl,I]", "[CX4][Cl,Br,I]", "[#6]S(=O)(=O)O[#6]",
        "[$([CH]),$(CC)]#CC(=O)[#6]", "[$([CH]),$(CC)]#CC(=O)O[#6]", "n[OH]",
        "[$([CH]),$(CC)]#CS(=O)(=O)[#6]", "C=C(C=O)C=O", "n1c([F,Cl,Br,I])cccc1", "[CH1](=O)", "[#8][#8]",
        "[C;!R]=[N;!R]", "[N!R]=[N!R]", "[#6](=O)[#6](=O)", "[#16][#16]", "[#7][NH2]", "C(=O)N[NH2]",
        "[#6]=S", "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
        "C1(=[O,N])C=CC(=[O,N])C=C1", "C1(=[O,N])C(=[O,N])C=CC=C1", "a21aa3a(aa1aaaa2)aaaa3",
        "a31a(a2a(aa1)aaaa2)aaaa3", "a1aa2a3a(a1)A=AA=A3=AA=A2", "c1cc([NH2])ccc1",
        "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
        "I", "OS(=O)(=O)[O-]", "[N+](=O)[O-]", "C(=O)N[OH]",
        "C1NC(=O)NC(=O)1", "[SH]", "[S-]", "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
        "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]", "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
        "[CR1]1[CR1][CR1]cc[CR1][CR1]1", "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
        "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1", "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
        "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1", "C#C",
        "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]", "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
        "[#6]=N[OH]", "[#6]=NOC=O", "[#6](=O)[CX4,CR0X3,O][#6](=O)", "c1ccc2c(c1)ccc(=O)o2",
        "[O+,o+,S+,s+]", "N=C=O", "[NX3,NX4][F,Cl,Br,I]", "c1ccccc1OC(=O)[#6]", "[CR0]=[CR0][CR0]=[CR0]",
        "[C+,c+,C-,c-]", "N=[N+]=[N-]", "C12C(NC(N1)=O)CSC2", "c1c([OH])c([OH,NH2,NH])ccc1", "P",
        "[N,O,S]C#N", "C=C=O", "[Si][F,Cl,Br,I]", "[SX2]O", "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
        "O1CCCCC1OC2CCC3CCCCC3C2", "N=[CR0][N,n,O,S]",
        "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
        "C=[C!r]C#N", "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
        "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
        "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])", "[OH]c1ccc([OH,NH2,NH])cc1",
        "c1ccccc1OC(=O)O", "[SX2H0][N]", "c12ccccc1(SC(S)=N2)", "c12ccccc1(SC(=S)N2)", "c1nnnn1C=O",
        "s1c(S)nnc1NC=O", "S1C=CSC1=S", "C(=O)Onnn", "OS(=O)(=O)C(F)(F)F", "N#CC[OH]", "N#CC(=O)",
        "S(=O)(=O)C#N", "N[CH2]C#N", "C1(=O)NCC1", "S(=O)(=O)[O-,OH]", "NC[F,Cl,Br,I]", "C=[C!r]O",
        "[NX2+0]=[O+0]", "[OR0,NR0][OR0,NR0]", "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]", "[CX2R0][NX3R0]",
        "c1ccccc1[C;!R]=[C;!R]c2ccccc2", "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
        "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
        "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]", "[*]=[N+]=[*]", "[SX3](=O)[O-,OH]",
        "N#N", "F.F.F.F", "[R0;D2][R0;D2][R0;D2][R0;D2]", "[cR,CR]~C(=O)NC(=O)~[cR,CR]", "C=!@CC=[O,S]",
        "[#6,#8,#16][#6](=O)O[#6]", "c[C;R0](=[O,S])[#6]", "c[SX2][C;!R]", "C=C=C",
        "c1nc([F,Cl,Br,I,S])ncc1", "c1ncnc([F,Cl,Br,I,S])c1", "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
        "[#6]S(=O)(=O)c1ccc(cc1)F", "[15N]", "[13C]", "[18O]", "[34S]",
    ];

    private static readonly Lazy<SmartsPattern[]> _acceptors = new(() => AcceptorSmarts.Select(s => SmartsPattern.Parse(s)).ToArray());
    private static readonly Lazy<SmartsPattern[]> _alerts = new(() => AlertSmarts.Select(s => SmartsPattern.Parse(s)).ToArray());
    private static readonly Lazy<SmartsPattern> _aliphaticRing = new(() => SmartsPattern.Parse("[$([A;R][!a])]"));

    /// <summary>The eight QED input properties.</summary>
    public static QedProperties Properties(Molecule mol)
    {
        int hba = 0;
        foreach (var p in _acceptors.Value) hba += p.CountMatches(mol);
        int alerts = 0;
        foreach (var p in _alerts.Value) if (p.IsMatch(mol)) alerts++;
        return new QedProperties(
            MolecularDescriptors.MolecularWeight(mol),
            Crippen.LogP(mol),
            hba,
            MolecularDescriptors.HBondDonors(mol),
            Tpsa.Compute(mol),
            MolecularDescriptors.RotatableBonds(mol, strict: true),
            AromaticRingsAfterRemovingAliphatic(mol),
            alerts);
    }

    /// <summary>Number of QED structural alerts matched.</summary>
    public static IReadOnlyList<string> MatchedAlerts(Molecule mol)
        => _alerts.Value.Where(p => p.IsMatch(mol)).Select(p => p.Text).ToList();

    /// <summary>The asymmetric double sigmoidal desirability of property <paramref name="index"/> (0 = MW … 7 = alerts).</summary>
    public static double Desirability(int index, double x)
    {
        var p = Ads[index];
        double exp1 = 1 + Math.Exp(-(x - p.C + p.D / 2) / p.E);
        double exp2 = 1 + Math.Exp(-(x - p.C - p.D / 2) / p.F);
        double dx = p.A + p.B / exp1 * (1 - 1 / exp2);
        return dx / p.DMax;
    }

    /// <summary>QED with mean weights (or <paramref name="weights"/>), in (0, 1].</summary>
    public static double Compute(Molecule mol, IReadOnlyList<double>? weights = null) => Compute(Properties(mol), weights);

    /// <summary>QED from precomputed properties.</summary>
    public static double Compute(QedProperties p, IReadOnlyList<double>? weights = null)
    {
        weights ??= MeanWeights;
        double[] x = [p.MW, p.ALogP, p.HBA, p.HBD, p.PSA, p.RotB, p.Arom, p.Alerts];
        double t = 0, sw = 0;
        for (int i = 0; i < 8; i++) { t += weights[i] * Math.Log(Desirability(i, x[i])); sw += weights[i]; }
        return Math.Exp(t / sw);
    }

    private static int AromaticRingsAfterRemovingAliphatic(Molecule mol)
    {
        var remove = _aliphaticRing.Value.MatchingAtoms(mol).ToList();
        var m = mol.Clone();
        m.RemoveAtoms(remove);
        return m.BondCount - m.AtomCount + m.Fragments().Count;
    }
}
