namespace Biotech.Net.Chemistry;

/// <summary>Orbital hybridisation of an atom, as perceived by <see cref="AtomHybridization"/>.</summary>
public enum Hybridization : byte
{
    /// <summary>Not determined (dummy atoms, unusual valence states).</summary>
    Unspecified,
    /// <summary>s (hydrogen, halide ions).</summary>
    S,
    /// <summary>sp (linear).</summary>
    SP,
    /// <summary>sp² (trigonal planar).</summary>
    SP2,
    /// <summary>sp³ (tetrahedral).</summary>
    SP3,
    /// <summary>sp³d (trigonal bipyramidal).</summary>
    SP3D,
    /// <summary>sp³d² (octahedral).</summary>
    SP3D2,
}

/// <summary>
/// Conjugation and hybridisation perception following the RDKit rules (<c>MolOps::setConjugation</c> and
/// <c>MolOps::setHybridization</c>), so force-field atom types and distance-geometry bounds agree with RDKit.
/// </summary>
/// <remarks>
/// A bond is conjugated when it is aromatic, or when a multiple bond and a neighbouring bond meet at a "conjugation
/// candidate" atom (an atom with π electrons or a lone pair to share, two or three substituents, and normal valence).
/// Hybridisation comes from the number of σ bonds plus lone pairs; four such orbitals give sp² instead of sp³ when the
/// atom carries a conjugated bond (amide N, ester O, aniline N, pyrrole N).
/// </remarks>
public static class AtomHybridization
{
    /// <summary>Per-bond conjugation flags (index = bond index).</summary>
    public static bool[] ConjugatedBonds(Molecule mol)
    {
        int n = mol.AtomCount;
        var bonds = mol.Bonds;
        var conj = new bool[bonds.Length];
        for (int b = 0; b < bonds.Length; b++) conj[b] = bonds[b].Order == BondOrder.Aromatic;
        var cand = new bool[n];
        var subst = new int[n];
        for (int i = 0; i < n; i++)
        {
            cand[i] = IsConjugationCandidate(mol, i);
            subst[i] = cand[i] ? TotalDegree(mol, i) : 0;
        }
        for (int at = 0; at < n; at++)
        {
            if (!cand[at] || subst[at] < 2 || subst[at] > 3) continue;
            var nb = mol.Neighbors(at);
            var nbb = mol.NeighborBonds(at);
            for (int x = 0; x < nb.Length; x++)
            {
                var o = bonds[nbb[x]].Order;
                if (o == BondOrder.Single || !cand[nb[x]]) continue; // valence contribution must be ≥ 1.5
                for (int y = 0; y < nb.Length; y++)
                {
                    if (y == x) continue;
                    int other = nb[y];
                    if (subst[other] > 3) continue;
                    if (cand[other]) { conj[nbb[x]] = true; conj[nbb[y]] = true; }
                }
            }
        }
        return conj;
    }

    /// <summary>Hybridisation of every atom.</summary>
    public static Hybridization[] Compute(Molecule mol) => Compute(mol, ConjugatedBonds(mol));

    /// <summary>Hybridisation of every atom, given precomputed bond conjugation flags.</summary>
    public static Hybridization[] Compute(Molecule mol, bool[] conjugated)
    {
        int n = mol.AtomCount;
        var result = new Hybridization[n];
        for (int i = 0; i < n; i++)
        {
            var a = mol.Atoms[i];
            if (a.AtomicNumber == 0) { result[i] = Hybridization.Unspecified; continue; }
            int totalDegree = TotalDegree(mol, i);
            if (a.Chirality != ChiralTag.None && totalDegree == 4) { result[i] = Hybridization.SP3; continue; }
            int norbs = a.AtomicNumber < 89 ? BondsPlusLonePairs(mol, i) : totalDegree;
            result[i] = norbs switch
            {
                0 or 1 => Hybridization.S,
                2 => Hybridization.SP,
                3 => Hybridization.SP2,
                4 => totalDegree > 3 || !HasConjugatedBond(mol, i, conjugated) ? Hybridization.SP3 : Hybridization.SP2,
                5 => Hybridization.SP3D,
                6 => Hybridization.SP3D2,
                _ => Hybridization.Unspecified,
            };
        }
        return result;
    }

    /// <summary>Whether any bond of <paramref name="atom"/> is conjugated.</summary>
    public static bool HasConjugatedBond(Molecule mol, int atom, bool[] conjugated)
    {
        foreach (int b in mol.NeighborBonds(atom)) if (conjugated[b]) return true;
        return false;
    }

    /// <summary>Explicit neighbours plus implicit/explicit hydrogen counts (RDKit <c>getTotalDegree</c>).</summary>
    internal static int TotalDegree(Molecule mol, int atom) => mol.Degree(atom) + mol.Atoms[atom].HydrogenCount;

    /// <summary>RDKit default valence: the first allowed valence, −1 for elements without one (metals), 0 for noble gases.</summary>
    internal static int DefaultValence(int atomicNumber)
    {
        if (atomicNumber <= 0) return -1;
        var e = PeriodicTable.Get(atomicNumber);
        if (e.Valences.Length > 0) return e.Valences[0];
        return atomicNumber is 2 or 10 or 18 or 36 or 54 or 86 ? 0 : -1;
    }

    private static bool IsConjugationCandidate(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        int z = a.AtomicNumber;
        if (z <= 0) return false;
        int dv = DefaultValence(z);
        int totalValence = mol.TotalValence(atom);
        if (a.FormalCharge == 0 && dv >= 0 && totalValence > dv) return false;
        int nouter = PeriodicTable.Get(z).OuterElectrons;
        bool ok = z <= 10 || (nouter != 5 && nouter != 6) || (nouter == 6 && TotalDegree(mol, atom) < 2);
        return ok && CountAtomElectrons(mol, atom) > 0;
    }

    // Sign-equivalent port of RDKit countAtomElec: only "> 0" matters for conjugation candidates.
    private static int CountAtomElectrons(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        int dv = DefaultValence(a.AtomicNumber);
        if (dv <= 1) return -1;
        int degree = TotalDegree(mol, atom);
        if (degree > 3) return -1;
        int nlp = PeriodicTable.Get(a.AtomicNumber).OuterElectrons - dv;
        nlp = Math.Max(nlp - a.FormalCharge, 0);
        return dv - degree + nlp;
    }

    private static int BondsPlusLonePairs(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        int deg = TotalDegree(mol, atom);
        if (a.AtomicNumber <= 1) return deg;
        int nouter = PeriodicTable.Get(a.AtomicNumber).OuterElectrons;
        int totalValence = mol.TotalValence(atom);
        int chg = a.FormalCharge;
        int free = nouter - (totalValence + chg);
        return deg + free / 2;
    }
}
