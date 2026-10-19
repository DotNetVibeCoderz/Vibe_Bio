using System.Globalization;
using System.Text;
using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>
/// Writes docked ligands in PDBQT layout (ROOT/BRANCH torsion tree, AutoDock atom types, zero partial charges;
/// non-polar hydrogens omitted as in AutoDock-prepared ligands). Readable by AutoDock Vina and Meeko.
/// </summary>
public static class PdbqtWriter
{
    /// <summary>All poses of a docking result as MODEL records with Vina's <c>REMARK VINA RESULT</c> lines.</summary>
    public static string Write(DockingResult result)
    {
        var lig = result.Ligand ?? throw new InvalidOperationException("The result carries no ligand model.");
        var sb = new StringBuilder();
        foreach (var p in result.Poses)
        {
            sb.Append(CultureInfo.InvariantCulture, $"MODEL {p.Rank}\n");
            sb.Append(CultureInfo.InvariantCulture, $"REMARK VINA RESULT: {p.Affinity,9:F3} {p.RmsdLowerBound,10:F3} {p.RmsdUpperBound,10:F3}\n");
            WriteLigand(sb, lig, p.Molecule);
            sb.Append("ENDMDL\n");
        }
        return sb.ToString();
    }

    /// <summary>A single ligand (with 3-D coordinates) as PDBQT, building its torsion tree.</summary>
    public static string Write(Molecule ligand)
    {
        var sb = new StringBuilder();
        WriteLigand(sb, new FlexibleLigand(ligand), ligand);
        return sb.ToString();
    }

    private static void WriteLigand(StringBuilder sb, FlexibleLigand lig, Molecule mol)
    {
        int n = mol.AtomCount;
        var coords = mol.Coordinates!;
        bool Keep(int i)
        {
            if (mol.Atoms[i].AtomicNumber != 1) return true;
            foreach (int nb in mol.Neighbors(i)) if (mol.Atoms[nb].AtomicNumber is 7 or 8) return true;
            return false;
        }
        var torsions = lig.Torsions;
        var inSub = torsions.Select(t => new HashSet<int>(t.Moving)).ToArray();
        // own atoms of each torsion node: its moving set minus those of nested torsions
        var own = new List<int>[torsions.Length];
        for (int t = 0; t < torsions.Length; t++)
        {
            own[t] = torsions[t].Moving.Where(a => !Enumerable.Range(0, torsions.Length).Any(u => u != t && inSub[u].Count < inSub[t].Count && inSub[u].Contains(a) && inSub[t].IsSupersetOf(inSub[u]))).ToList();
        }
        var rootAtoms = Enumerable.Range(0, n).Where(a => !inSub.Any(s => s.Contains(a))).ToList();
        var serial = new int[n];
        int next = 1;
        sb.Append(CultureInfo.InvariantCulture, $"REMARK SMILES {mol.ToSmiles(includeStereo: true)}\n");
        sb.Append("ROOT\n");
        foreach (int a in rootAtoms) if (Keep(a)) Atom(sb, mol, a, coords[a], serial, ref next);
        sb.Append("ENDROOT\n");
        void Branches(IEnumerable<int> parentAtoms)
        {
            var parents = new HashSet<int>(parentAtoms);
            for (int t = torsions.Length - 1; t >= 0; t--)
            {
                if (!parents.Contains(torsions[t].From)) continue;
                int from = serial[torsions[t].From];
                int startSerial = next;
                sb.Append(CultureInfo.InvariantCulture, $"BRANCH {from,3} {startSerial,3}\n");
                // write the axis atom first, as AutoDock expects
                var atoms = own[t].OrderBy(a => a == torsions[t].To ? 0 : 1).ToList();
                foreach (int a in atoms) if (Keep(a)) Atom(sb, mol, a, coords[a], serial, ref next);
                Branches(own[t]);
                sb.Append(CultureInfo.InvariantCulture, $"ENDBRANCH {from,3} {startSerial,3}\n");
            }
        }
        Branches(rootAtoms);
        sb.Append(CultureInfo.InvariantCulture, $"TORSDOF {lig.RotatableBonds.Length}\n");
    }

    private static void Atom(StringBuilder sb, Molecule mol, int i, Point3 p, int[] serial, ref int next)
    {
        var a = mol.Atoms[i];
        string type = AdType(mol, i);
        string name = (a.Symbol + (i + 1).ToString(CultureInfo.InvariantCulture));
        if (name.Length > 4) name = name[..4];
        serial[i] = next;
        sb.Append(string.Format(CultureInfo.InvariantCulture, "ATOM  {0,5} {1,-4} UNL     1    {2,8:F3}{3,8:F3}{4,8:F3}  1.00  0.00    {5,6:F3} {6,-2}\n",
            next, name, p.X, p.Y, p.Z, 0.0, type));
        next++;
    }

    private static string AdType(Molecule mol, int i)
    {
        var a = mol.Atoms[i];
        return a.AtomicNumber switch
        {
            1 => "HD",
            6 => a.IsAromatic ? "A" : "C",
            7 => VinaAtomTypes.AssignLigand(mol)[i] is XsAtomType.N_A or XsAtomType.N_DA ? "NA" : "N",
            8 => "OA",
            16 => "SA",
            _ => a.Symbol.Length == 1 ? a.Symbol : a.Symbol,
        };
    }
}
