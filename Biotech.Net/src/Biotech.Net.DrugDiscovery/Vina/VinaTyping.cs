using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;
using SbAtom = Biotech.Net.StructuralBiology.Atom;

namespace Biotech.Net.DrugDiscovery;

/// <summary>X-Score atom types used by the AutoDock Vina scoring function (Trott &amp; Olson 2010).</summary>
public enum XsAtomType : byte
{
    /// <summary>Carbon not bonded to a heteroatom (hydrophobic).</summary>
    C_H,
    /// <summary>Carbon bonded to a heteroatom.</summary>
    C_P,
    /// <summary>Nitrogen, neither donor nor acceptor.</summary>
    N_P,
    /// <summary>Nitrogen donor.</summary>
    N_D,
    /// <summary>Nitrogen acceptor.</summary>
    N_A,
    /// <summary>Nitrogen donor and acceptor.</summary>
    N_DA,
    /// <summary>Oxygen, neither donor nor acceptor.</summary>
    O_P,
    /// <summary>Oxygen donor.</summary>
    O_D,
    /// <summary>Oxygen acceptor.</summary>
    O_A,
    /// <summary>Oxygen donor and acceptor.</summary>
    O_DA,
    /// <summary>Sulfur.</summary>
    S_P,
    /// <summary>Phosphorus.</summary>
    P_P,
    /// <summary>Fluorine (hydrophobic).</summary>
    F_H,
    /// <summary>Chlorine (hydrophobic).</summary>
    Cl_H,
    /// <summary>Bromine (hydrophobic).</summary>
    Br_H,
    /// <summary>Iodine (hydrophobic).</summary>
    I_H,
    /// <summary>Silicon.</summary>
    Si,
    /// <summary>Astatine.</summary>
    At,
    /// <summary>Metal (treated as a donor).</summary>
    Met_D,
}

/// <summary>Vina atom-type properties and typing rules for ligands (Meeko AD4 rules) and receptors.</summary>
public static class VinaAtomTypes
{
    internal const int Count = 19;
    private static readonly double[] Radii = [1.9, 1.9, 1.8, 1.8, 1.8, 1.8, 1.7, 1.7, 1.7, 1.7, 2.0, 2.1, 1.5, 1.8, 2.0, 2.2, 2.2, 2.3, 1.2];

    /// <summary>X-Score van der Waals radius (Å) used for the surface distance d = r − R₁ − R₂.</summary>
    public static double Radius(XsAtomType t) => Radii[(int)t];

    /// <summary>Hydrophobic types: C_H and the halogens.</summary>
    public static bool IsHydrophobic(XsAtomType t) => t is XsAtomType.C_H or XsAtomType.F_H or XsAtomType.Cl_H or XsAtomType.Br_H or XsAtomType.I_H;

    /// <summary>H-bond donor types.</summary>
    public static bool IsDonor(XsAtomType t) => t is XsAtomType.N_D or XsAtomType.N_DA or XsAtomType.O_D or XsAtomType.O_DA or XsAtomType.Met_D;

    /// <summary>H-bond acceptor types.</summary>
    public static bool IsAcceptor(XsAtomType t) => t is XsAtomType.N_A or XsAtomType.N_DA or XsAtomType.O_A or XsAtomType.O_DA;

    /// <summary>Whether a donor–acceptor pair is possible between the two types.</summary>
    public static bool HBondPossible(XsAtomType a, XsAtomType b) => (IsDonor(a) && IsAcceptor(b)) || (IsDonor(b) && IsAcceptor(a));

    private static readonly Lazy<SmartsPattern[]> NonAcceptorN = new(() =>
    [
        SmartsPattern.Parse("[#7X3v3][a]"),          // pyrrole, aniline
        SmartsPattern.Parse("[#7X3v3][#6X3v4]"),     // amide, amidine NH2
        SmartsPattern.Parse("[#7X3v3][NX2]=[*]"),    // triazene
        SmartsPattern.Parse("[#7+1]"),               // ammonium, pyridinium
    ]);

    /// <summary>
    /// Types every atom of a ligand; hydrogens (and atoms Vina does not score) get <c>null</c>. Rules follow Meeko's
    /// AD4 typing as read by Vina: N is an acceptor unless it is an aniline/pyrrole/amide/triazene N or cationic;
    /// O is always an acceptor; N/O carrying hydrogens are donors; C is polar when bonded to any heteroatom.
    /// </summary>
    public static XsAtomType?[] AssignLigand(Molecule mol)
    {
        int n = mol.AtomCount;
        var result = new XsAtomType?[n];
        var nonAcc = new bool[n];
        foreach (var p in NonAcceptorN.Value)
            foreach (var m in p.Matches(mol, unique: false)) nonAcc[m[0]] = true;
        for (int i = 0; i < n; i++)
        {
            int z = mol.Atoms[i].AtomicNumber;
            bool donor = mol.TotalHydrogens(i) > 0;
            result[i] = z switch
            {
                6 => BondedToHetero(mol, i) ? XsAtomType.C_P : XsAtomType.C_H,
                7 => (!nonAcc[i], donor) switch { (true, true) => XsAtomType.N_DA, (true, false) => XsAtomType.N_A, (false, true) => XsAtomType.N_D, _ => XsAtomType.N_P },
                8 => donor ? XsAtomType.O_DA : XsAtomType.O_A,
                16 or 34 => XsAtomType.S_P,
                15 => XsAtomType.P_P,
                9 => XsAtomType.F_H,
                17 => XsAtomType.Cl_H,
                35 => XsAtomType.Br_H,
                53 => XsAtomType.I_H,
                14 => XsAtomType.Si,
                85 => XsAtomType.At,
                0 or 1 or 5 => null,
                _ => PeriodicTable.Get(z).IsMetal ? XsAtomType.Met_D : null,
            };
        }
        return result;
    }

    private static bool BondedToHetero(Molecule mol, int atom)
    {
        foreach (int nb in mol.Neighbors(atom)) if (mol.Atoms[nb].AtomicNumber is not (1 or 6)) return true;
        return false;
    }

    /// <summary>AutoDock covalent radius (Å) used by Vina's distance-based bond perception.</summary>
    internal static double CovalentRadius(string element) => element switch
    {
        "C" => 0.77, "N" => 0.75, "O" => 0.73, "P" => 1.06, "S" => 1.02, "H" => 0.37, "F" => 0.71, "I" => 1.33,
        "CL" => 0.99, "BR" => 1.14, "SI" => 1.11, "MG" => 1.30, "MN" => 1.39, "ZN" => 1.31, "CA" => 1.74, "FE" => 1.25,
        "SE" => 1.02,
        _ => 1.75,
    };
}

/// <summary>
/// A rigid receptor typed for Vina scoring: heavy atoms with X-Score types. Donors come from explicit hydrogens when
/// the structure has them (e.g. after protonation), otherwise from standard amino-acid templates (histidine as HIE).
/// </summary>
public sealed class VinaReceptor
{
    private static readonly Dictionary<string, string[]> TemplateDonors = new()
    {
        ["ARG"] = ["NE", "NH1", "NH2"], ["ASN"] = ["ND2"], ["GLN"] = ["NE2"], ["HIS"] = ["NE2"], ["LYS"] = ["NZ"],
        ["TRP"] = ["NE1"], ["SER"] = ["OG"], ["THR"] = ["OG1"], ["TYR"] = ["OH"], ["HOH"] = ["O"], ["WAT"] = ["O"],
    };

    private VinaReceptor(Vec3[] positions, XsAtomType[] types, SbAtom[] atoms, bool fromHydrogens)
    {
        Positions = positions;
        Types = types;
        Atoms = atoms;
        DonorsFromHydrogens = fromHydrogens;
    }

    /// <summary>Heavy-atom positions (Å).</summary>
    public IReadOnlyList<Vec3> Positions { get; }

    /// <summary>X-Score type of each heavy atom.</summary>
    public IReadOnlyList<XsAtomType> Types { get; }

    /// <summary>The source atoms (parallel to <see cref="Positions"/>).</summary>
    public IReadOnlyList<SbAtom> Atoms { get; }

    /// <summary>True when donors were derived from explicit hydrogens rather than residue templates.</summary>
    public bool DonorsFromHydrogens { get; }

    /// <summary>Number of scored (heavy) atoms.</summary>
    public int Count => Positions.Count;

    /// <summary>Builds a receptor from the first model of <paramref name="structure"/>; ligands and waters are excluded by default.</summary>
    public static VinaReceptor FromStructure(Structure structure, bool includeHetero = false, bool includeWater = false)
        => FromAtoms(structure.FirstModel.Atoms.Where(a =>
        {
            var r = a.Residue;
            if (r is null) return true;
            if (r.IsWater) return includeWater;
            return includeHetero || !a.IsHetero || r.IsAminoAcid;
        }).ToList());

    /// <summary>Reads a PDB file and builds the receptor (protein atoms only).</summary>
    public static VinaReceptor FromPdb(string path) => FromStructure(PdbReader.Read(path));

    /// <summary>Builds a receptor from explicit atoms (hydrogens, if present, are used for donor typing only).</summary>
    public static VinaReceptor FromAtoms(IReadOnlyList<SbAtom> atoms)
    {
        int n = atoms.Count;
        bool hasH = atoms.Any(a => a.Element == "H");
        var pos = atoms.Select(a => a.Position).ToArray();
        var search = new NeighborSearch(pos, 3.0);
        var neighbors = new List<int>[n];
        var found = new List<int>();
        for (int i = 0; i < n; i++)
        {
            neighbors[i] = [];
            double ri = VinaAtomTypes.CovalentRadius(atoms[i].Element);
            found.Clear();
            search.Query(pos[i], 1.1 * (ri + 1.75), found);
            foreach (int j in found)
            {
                if (j == i) continue;
                double rj = VinaAtomTypes.CovalentRadius(atoms[j].Element);
                if (Vec3.Distance(pos[i], pos[j]) < 1.1 * (ri + rj)) neighbors[i].Add(j);
            }
        }
        var heavy = new List<int>();
        var types = new List<XsAtomType>();
        for (int i = 0; i < n; i++)
        {
            var a = atoms[i];
            string el = a.Element;
            if (el == "H") continue;
            XsAtomType? t;
            if (el == "C") t = neighbors[i].Any(j => atoms[j].Element is not ("C" or "H")) ? XsAtomType.C_P : XsAtomType.C_H;
            else if (el is "N" or "O")
            {
                bool donor, acceptor;
                int heavyDeg = neighbors[i].Count(j => atoms[j].Element != "H");
                if (hasH)
                {
                    bool hasHn = neighbors[i].Any(j => atoms[j].Element == "H");
                    donor = hasHn;
                    acceptor = el == "O" || (!hasHn && heavyDeg <= 2);
                }
                else
                {
                    string res = a.ResidueName;
                    bool aa = a.Residue?.IsAminoAcid ?? false;
                    donor = (aa && a.Name == "N" && res != "PRO") || (TemplateDonors.TryGetValue(res, out var d) && d.Contains(a.Name));
                    acceptor = el == "O" || (res == "HIS" && a.Name == "ND1") || (!aa && heavyDeg <= 2 && !donor);
                }
                t = el == "N"
                    ? (acceptor, donor) switch { (true, true) => XsAtomType.N_DA, (true, false) => XsAtomType.N_A, (false, true) => XsAtomType.N_D, _ => XsAtomType.N_P }
                    : (acceptor, donor) switch { (true, true) => XsAtomType.O_DA, (true, false) => XsAtomType.O_A, (false, true) => XsAtomType.O_D, _ => XsAtomType.O_P };
            }
            else t = el switch
            {
                "S" or "SE" => XsAtomType.S_P,
                "P" => XsAtomType.P_P,
                "F" => XsAtomType.F_H,
                "CL" => XsAtomType.Cl_H,
                "BR" => XsAtomType.Br_H,
                "I" => XsAtomType.I_H,
                "SI" => XsAtomType.Si,
                "ZN" or "MG" or "CA" or "FE" or "MN" or "CU" or "NI" or "CO" or "NA" or "K" or "CD" or "HG" => XsAtomType.Met_D,
                _ => null,
            };
            if (t is null) continue;
            heavy.Add(i);
            types.Add(t.Value);
        }
        return new VinaReceptor(heavy.Select(i => pos[i]).ToArray(), [.. types], heavy.Select(i => atoms[i]).ToArray(), hasH);
    }
}
