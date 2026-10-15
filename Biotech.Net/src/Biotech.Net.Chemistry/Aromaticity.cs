namespace Biotech.Net.Chemistry;

/// <summary>
/// Kekulization (assigning alternating single/double bonds to aromatic systems) and Hückel aromaticity perception.
/// </summary>
/// <remarks>
/// <para><b>Kekulization</b>: aromatic atoms that still need a π bond (e.g. <c>c</c> with two ring bonds, pyridine
/// <c>n</c>) must be paired through aromatic bonds; a perfect matching is found by most-constrained-first
/// backtracking. When no matching exists, aromatic nitrogens left unpaired are given a hydrogen (the common
/// "c1ccnc1" pyrrole typo) before giving up.</para>
/// <para><b>Perception</b> follows the RDKit-style model: each SSSR ring, and each pair of fused rings, is aromatic
/// when every atom can contribute π electrons and the total is 4n+2. Contributions: atom in an endocyclic double
/// bond 1; exocyclic double bond to N/O/S 0; lone-pair donors (pyrrole N, furan O, thiophene S, carbanion) 2;
/// carbocation/boron 0.</para>
/// </remarks>
public static class Aromaticity
{
    /// <summary>A copy of <paramref name="molecule"/> in which aromatic bonds are replaced by a Kekulé structure and aromatic flags are cleared.</summary>
    /// <exception cref="InvalidOperationException">No Kekulé structure exists.</exception>
    public static Molecule Kekulize(Molecule molecule)
    {
        var orders = ComputeKekuleOrders(molecule, throwOnFailure: true);
        var m = molecule.Clone();
        for (int i = 0; i < m.BondCount; i++)
        {
            ref var b = ref m.BondAt(i);
            if (b.Order == BondOrder.Aromatic) b.Order = (BondOrder)orders[i];
        }
        for (int i = 0; i < m.AtomCount; i++) m.AtomAt(i).IsAromatic = false;
        m.InvalidateTopology();
        return m;
    }

    /// <summary>Whether an aromatic atom needs a double bond within the aromatic system.</summary>
    internal static bool NeedsPiBond(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        if (!a.IsAromatic && !HasAromaticBond(mol, atom)) return false;
        int sigma = a.NoImplicit ? a.ExplicitHydrogens : 0;
        foreach (int b in mol.NeighborBonds(atom))
        {
            var bond = mol.Bonds[b];
            sigma += bond.Order == BondOrder.Aromatic ? 1 : (int)bond.Order;
        }
        var allowed = Molecule.AllowedValences(a.AtomicNumber, a.FormalCharge);
        foreach (int v in allowed)
            if (v >= sigma) return v - sigma >= 1;
        return false;
    }

    private static bool HasAromaticBond(Molecule mol, int atom)
    {
        foreach (int b in mol.NeighborBonds(atom)) if (mol.Bonds[b].Order == BondOrder.Aromatic) return true;
        return false;
    }

    /// <summary>Kekulé bond orders (1/2/3) for every bond; aromatic bonds receive 1 or 2.</summary>
    internal static byte[] ComputeKekuleOrders(Molecule mol, bool throwOnFailure)
    {
        int m = mol.BondCount, n = mol.AtomCount;
        var orders = new byte[m];
        var bonds = mol.Bonds;
        bool any = false;
        for (int i = 0; i < m; i++)
        {
            if (bonds[i].Order == BondOrder.Aromatic) { orders[i] = 1; any = true; }
            else orders[i] = (byte)bonds[i].Order;
        }
        if (!any) return orders;

        var need = new bool[n];
        for (int i = 0; i < n; i++) need[i] = NeedsPiBond(mol, i);
        if (TryMatch(mol, need, orders)) return orders;

        // Repair: unmatched organic-subset aromatic nitrogens become [nH].
        var trial = (byte[])orders.Clone();
        var need2 = (bool[])need.Clone();
        bool changed = false;
        for (int i = 0; i < n; i++)
        {
            var a = mol.Atoms[i];
            if (need2[i] && a.AtomicNumber == 7 && a.IsAromatic && !a.NoImplicit && a.FormalCharge == 0)
            {
                need2[i] = false;
                if (TryMatch(mol, need2, trial))
                {
                    mol.AtomAt(i).ExplicitHydrogens = 1;
                    mol.AtomAt(i).NoImplicit = true;
                    return trial;
                }
                need2[i] = true;
                trial = (byte[])orders.Clone();
                changed = true;
            }
        }
        _ = changed;
        if (throwOnFailure) throw new InvalidOperationException($"Cannot kekulize molecule '{mol.Name}'.");
        // Best effort: leave unmatched bonds single.
        TryMatchPartial(mol, need, orders);
        return orders;
    }

    private static bool TryMatch(Molecule mol, bool[] need, byte[] orders)
    {
        int n = mol.AtomCount;
        var matched = new bool[n];
        var bonds = mol.Bonds.ToArray();
        var pending = new List<int>();
        for (int i = 0; i < n; i++) if (need[i]) pending.Add(i);
        if (pending.Count == 0) return true;
        if (pending.Count % 2 == 1) { /* still try: odd counts cannot match */ return false; }
        long budget = 200_000;
        return Solve();

        bool Solve()
        {
            if (--budget < 0) return false;
            // most constrained unmatched atom
            int best = -1, bestOptions = int.MaxValue;
            foreach (int a in pending)
            {
                if (matched[a]) continue;
                int options = 0;
                var nb = mol.Neighbors(a);
                var nbb = mol.NeighborBonds(a);
                for (int k = 0; k < nb.Length; k++)
                    if (bonds[nbb[k]].Order == BondOrder.Aromatic && need[nb[k]] && !matched[nb[k]]) options++;
                if (options < bestOptions) { bestOptions = options; best = a; if (options == 0) break; }
            }
            if (best < 0) return true;
            if (bestOptions == 0) return false;
            var nbs = mol.Neighbors(best).ToArray();
            var nbbs = mol.NeighborBonds(best).ToArray();
            for (int k = 0; k < nbs.Length; k++)
            {
                int v = nbs[k], b = nbbs[k];
                if (bonds[b].Order != BondOrder.Aromatic || !need[v] || matched[v]) continue;
                matched[best] = matched[v] = true;
                orders[b] = 2;
                if (Solve()) return true;
                orders[b] = 1;
                matched[best] = matched[v] = false;
            }
            return false;
        }
    }

    private static void TryMatchPartial(Molecule mol, bool[] need, byte[] orders)
    {
        var matched = new bool[mol.AtomCount];
        var bonds = mol.Bonds;
        for (int i = 0; i < bonds.Length; i++)
        {
            if (bonds[i].Order != BondOrder.Aromatic) continue;
            int a = bonds[i].Begin, b = bonds[i].End;
            if (need[a] && need[b] && !matched[a] && !matched[b]) { orders[i] = 2; matched[a] = matched[b] = true; }
        }
    }

    /// <summary>
    /// Clears aromatic flags and re-derives them with the Hückel rule (see type remarks). Aromatic bonds are first
    /// kekulized so that Kekulé and aromatic inputs give identical results. Implicit hydrogens are preserved.
    /// </summary>
    public static void Perceive(Molecule mol)
    {
        int n = mol.AtomCount, m = mol.BondCount;
        var kek = ComputeKekuleOrders(mol, throwOnFailure: false);
        // Freeze hydrogen counts: they were computed against the Kekulé form.
        for (int i = 0; i < m; i++)
        {
            ref var b = ref mol.BondAt(i);
            if (b.Order == BondOrder.Aromatic) b.Order = (BondOrder)kek[i];
        }
        for (int i = 0; i < n; i++) mol.AtomAt(i).IsAromatic = false;
        mol.InvalidateTopology();
        var rings = mol.Rings;
        if (rings.Count == 0) return;

        var aromaticAtom = new bool[n];
        var aromaticBond = new bool[m];
        var ringList = rings.AtomRings;
        var bondRings = rings.BondRings;

        var contrib = new int[n];
        for (int i = 0; i < n; i++) contrib[i] = PiElectrons(mol, i, rings);

        bool TryRingSet(IReadOnlyList<int> ringIdx)
        {
            var atoms = new HashSet<int>();
            var bondUse = new Dictionary<int, int>();
            foreach (int r in ringIdx)
            {
                atoms.UnionWith(ringList[r]);
                foreach (int b in bondRings[r]) bondUse[b] = bondUse.GetValueOrDefault(b) + 1;
            }
            int electrons = 0;
            foreach (int a in atoms)
            {
                if (contrib[a] < 0) return false;
                electrons += contrib[a];
            }
            if (electrons % 4 != 2) return false;
            foreach (int a in atoms) aromaticAtom[a] = true;
            // fused systems that are aromatic only as a whole (azulene) keep their fusion bonds non-aromatic, as RDKit does
            foreach (var (b, uses) in bondUse) if (uses == 1) aromaticBond[b] = true;
            return true;
        }

        for (int r = 0; r < ringList.Count; r++) TryRingSet([r]);
        // Fused pairs that fail individually (e.g. azulene, some heteroaromatics).
        for (int r = 0; r < ringList.Count; r++)
            for (int s = r + 1; s < ringList.Count; s++)
            {
                if (bondRings[r].All(b => aromaticBond[b]) && bondRings[s].All(b => aromaticBond[b])) continue;
                if (!bondRings[r].Intersect(bondRings[s]).Any()) continue;
                TryRingSet([r, s]);
            }

        for (int i = 0; i < m; i++)
        {
            if (!aromaticBond[i]) continue;
            ref var b = ref mol.BondAt(i);
            b.Order = BondOrder.Aromatic;
        }
        for (int i = 0; i < n; i++) if (aromaticAtom[i]) mol.AtomAt(i).IsAromatic = true;
        mol.InvalidateTopology();
        // Keep the Kekulé assignment we started from (consistent with the hydrogen counts).
        var cache = new byte[m];
        for (int i = 0; i < m; i++) cache[i] = kek[i];
        mol.SetKekuleCache(cache);
    }

    /// <summary>
    /// π-electron contribution of <paramref name="atom"/> to any ring it is in, or −1 if it cannot be aromatic. As in
    /// RDKit, a double bond that lies in any ring counts as one electron, whichever ring it belongs to.
    /// </summary>
    private static int PiElectrons(Molecule mol, int atom, RingInfo rings)
    {
        var a = mol.Atoms[atom];
        if (!rings.IsAtomInRing(atom)) return -1;
        int z = a.AtomicNumber, charge = a.FormalCharge;
        int heavy = 0, doubleIn = 0, doubleOutNeg = 0, doubleOutOther = 0, triple = 0;
        var nb = mol.Neighbors(atom);
        var nbb = mol.NeighborBonds(atom);
        for (int k = 0; k < nb.Length; k++)
        {
            int order = (int)mol.Bonds[nbb[k]].Order;
            heavy++;
            if (order == 2)
            {
                if (rings.IsBondInRing(nbb[k])) doubleIn++;
                else if (mol.Atoms[nb[k]].AtomicNumber is 7 or 8 or 16 or 34) doubleOutNeg++;
                else doubleOutOther++;
            }
            else if (order == 3) triple++;
        }
        int h = a.HydrogenCount;
        int connections = heavy + h;
        if (triple > 0 || doubleOutOther > 0) return -1;
        if (doubleIn > 1) return -1;
        switch (z)
        {
            case 6:
                if (doubleIn == 1) return 1;
                if (doubleOutNeg == 1) return 0;
                if (charge == -1 && connections == 3) return 2;
                if (charge == 1 && connections == 3) return 0;
                return -1;
            case 7:
            case 15:
                if (doubleIn == 1) return 1;
                if (doubleOutNeg == 1 && charge == 1) return 0;
                if (charge == 0 && connections == 3) return 2;
                if (charge == -1 && connections == 2) return 2;
                return -1;
            case 8:
            case 16:
            case 34:
            case 52:
                if (doubleIn == 1 && charge == 1) return 1;
                if (doubleIn == 1 && z != 8 && connections == 2) return 1;
                if (charge == 0 && connections == 2 && doubleOutNeg == 0) return 2;
                return -1;
            case 5:
                if (doubleIn == 1) return 1;
                if (charge == 0 && connections == 3) return 0;
                if (charge == -1 && connections == 3 && doubleIn == 1) return 1;
                return -1;
            case 33:
                if (doubleIn == 1) return 1;
                if (charge == 0 && connections == 3) return 2;
                return -1;
            default:
                return -1;
        }
    }
}
