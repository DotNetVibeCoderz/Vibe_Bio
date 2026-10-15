using System.Numerics;

namespace Biotech.Net.StructuralBiology;

/// <summary>Kinds of protein–ligand interaction recognised by <see cref="InteractionFingerprint"/>.</summary>
public enum InteractionType
{
    /// <summary>Polar heavy-atom pair (N/O ··· N/O) within H-bond distance, with compatible donor/acceptor roles.</summary>
    HydrogenBond = 0,

    /// <summary>Carbon/halogen contact between apolar ligand atoms and apolar side-chain atoms.</summary>
    Hydrophobic = 1,

    /// <summary>Salt bridge between oppositely charged groups.</summary>
    Ionic = 2,

    /// <summary>Aromatic ring stacking (face-to-face or edge-to-face/T-shaped).</summary>
    PiStacking = 3,
}

/// <summary>One detected interaction. For π-stacking the atoms are the first atoms of each ring and
/// <see cref="Distance"/> is the ring-centroid distance.</summary>
public sealed record Interaction(InteractionType Type, Atom LigandAtom, Atom ProteinAtom, double Distance, string Detail = "")
{
    /// <summary>The protein residue involved.</summary>
    public Residue ProteinResidue => ProteinAtom.Residue!;

    public override string ToString() => $"{Type} {LigandAtom.Name} – {ProteinResidue.Id}:{ProteinAtom.Name} {Distance:F2} Å {Detail}".TrimEnd();
}

/// <summary>Distance/angle criteria for <see cref="InteractionFingerprint"/> (PLIP/ProLIF-like defaults).</summary>
public sealed record InteractionOptions
{
    /// <summary>Maximum donor–acceptor heavy-atom distance (3.5 Å).</summary>
    public double HBondMaxDistance { get; init; } = 3.5;

    /// <summary>Minimum donor–acceptor distance (2.4 Å) to reject clashes and covalent contacts.</summary>
    public double HBondMinDistance { get; init; } = 2.4;

    /// <summary>Maximum apolar contact distance (4.0 Å).</summary>
    public double HydrophobicMaxDistance { get; init; } = 4.0;

    /// <summary>Maximum distance between charged atoms (4.0 Å).</summary>
    public double IonicMaxDistance { get; init; } = 4.0;

    /// <summary>Maximum centroid distance for face-to-face stacking (5.5 Å).</summary>
    public double PiParallelMaxDistance { get; init; } = 5.5;

    /// <summary>Maximum centroid distance for T-shaped stacking (6.5 Å).</summary>
    public double PiTShapedMaxDistance { get; init; } = 6.5;

    /// <summary>Maximum angle between ring normals for parallel stacking (30°); T-shaped needs 60–90°.</summary>
    public double PiParallelMaxAngle { get; init; } = 30;

    /// <summary>Treat His side-chain nitrogens as positively charged (default false).</summary>
    public bool HistidineCharged { get; init; }
}

/// <summary>The interactions of one ligand with its protein environment.</summary>
public sealed class InteractionProfile
{
    internal InteractionProfile(Residue ligand, List<Interaction> interactions)
    {
        Ligand = ligand;
        Interactions = interactions;
    }

    /// <summary>The ligand residue.</summary>
    public Residue Ligand { get; }

    /// <summary>All interactions found.</summary>
    public IReadOnlyList<Interaction> Interactions { get; }

    /// <summary>Distinct interacting protein residues, ordered by chain and number.</summary>
    public IReadOnlyList<Residue> Residues => Interactions.Select(i => i.ProteinResidue).Distinct()
        .OrderBy(r => r.Chain?.Id).ThenBy(r => r.Number).ToList();

    /// <summary>Interactions of one type.</summary>
    public IEnumerable<Interaction> OfType(InteractionType type) => Interactions.Where(i => i.Type == type);

    /// <summary>Whether <paramref name="residue"/> makes an interaction of <paramref name="type"/>.</summary>
    public bool Has(Residue residue, InteractionType type) => Interactions.Any(i => i.Type == type && ReferenceEquals(i.ProteinResidue, residue));

    /// <summary>
    /// The interaction fingerprint over a fixed residue list: 4 bits per residue in the order HydrogenBond,
    /// Hydrophobic, Ionic, PiStacking (residues are matched by chain id and number, so poses of the same target are
    /// comparable).
    /// </summary>
    public InteractionBits ToFingerprint(IReadOnlyList<Residue> residues)
    {
        var bits = new InteractionBits(residues.Count * 4);
        for (int r = 0; r < residues.Count; r++)
        {
            var res = residues[r];
            foreach (var i in Interactions)
                if (i.ProteinResidue.Number == res.Number && i.ProteinResidue.Chain?.Id == res.Chain?.Id && i.ProteinResidue.InsertionCode == res.InsertionCode)
                    bits.Set(r * 4 + (int)i.Type);
        }
        return bits;
    }

    /// <summary>Fingerprint over the residues this profile itself interacts with.</summary>
    public InteractionBits ToFingerprint() => ToFingerprint(Residues);

    /// <summary>Short text, e.g. <c>A:ASP189 Ionic,HydrogenBond</c> per line.</summary>
    public override string ToString() => string.Join(Environment.NewLine, Residues.Select(r =>
        $"{r.Id} {string.Join(",", Interactions.Where(i => i.ProteinResidue == r).Select(i => i.Type).Distinct())}"));
}

/// <summary>A bit-packed interaction fingerprint.</summary>
public sealed class InteractionBits
{
    private readonly ulong[] _words;

    /// <summary>Creates an all-zero fingerprint of <paramref name="length"/> bits.</summary>
    public InteractionBits(int length)
    {
        Length = length;
        _words = new ulong[(length + 63) / 64];
    }

    /// <summary>Number of bits.</summary>
    public int Length { get; }

    /// <summary>Gets a bit.</summary>
    public bool this[int index] => (_words[index >> 6] & (1UL << (index & 63))) != 0;

    /// <summary>Sets a bit.</summary>
    public void Set(int index) => _words[index >> 6] |= 1UL << (index & 63);

    /// <summary>Number of set bits.</summary>
    public int PopCount() => _words.Sum(w => BitOperations.PopCount(w));

    /// <summary>Tanimoto (Jaccard) similarity |A∧B| / |A∨B| (1 when both are empty).</summary>
    public static double Tanimoto(InteractionBits a, InteractionBits b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Fingerprints differ in length.");
        int both = 0, either = 0;
        for (int i = 0; i < a._words.Length; i++)
        {
            both += BitOperations.PopCount(a._words[i] & b._words[i]);
            either += BitOperations.PopCount(a._words[i] | b._words[i]);
        }
        return either == 0 ? 1 : (double)both / either;
    }

    /// <summary>The bits as a 0/1 string.</summary>
    public override string ToString() => string.Concat(Enumerable.Range(0, Length).Select(i => this[i] ? '1' : '0'));
}

/// <summary>
/// Rule-based protein–ligand interaction fingerprint (in the spirit of PLIP, Salentin 2015, and ProLIF):
/// hydrogen bonds from heavy-atom geometry and donor/acceptor roles, hydrophobic contacts, salt bridges between
/// charged groups (ligand charges perceived from carboxylate/phosphate/sulfonate, amidine/guanidine and aliphatic
/// amine patterns, or explicit formal charges), and π-stacking between planar 5/6-membered rings. Ligand bonds are
/// perceived from covalent radii. Hydrogens are not required.
/// </summary>
public static class InteractionFingerprint
{
    private static readonly Dictionary<string, string[]> AromaticRings = new()
    {
        ["PHE"] = ["CG,CD1,CD2,CE1,CE2,CZ"], ["TYR"] = ["CG,CD1,CD2,CE1,CE2,CZ"],
        ["TRP"] = ["CD2,CE2,CE3,CZ2,CZ3,CH2", "CG,CD1,NE1,CE2,CD2"], ["HIS"] = ["CG,ND1,CD2,CE1,NE2"],
    };

    private static readonly HashSet<string> HydrophobicSideChain =
    [
        "ALA.CB", "VAL.CB", "VAL.CG1", "VAL.CG2", "LEU.CB", "LEU.CG", "LEU.CD1", "LEU.CD2", "ILE.CB", "ILE.CG1", "ILE.CG2", "ILE.CD1",
        "MET.CB", "MET.CG", "MET.SD", "MET.CE", "PHE.CB", "PHE.CG", "PHE.CD1", "PHE.CD2", "PHE.CE1", "PHE.CE2", "PHE.CZ",
        "TRP.CB", "TRP.CG", "TRP.CD2", "TRP.CE3", "TRP.CZ2", "TRP.CZ3", "TRP.CH2", "TYR.CB", "TYR.CG", "TYR.CD1", "TYR.CD2",
        "TYR.CE1", "TYR.CE2", "PRO.CB", "PRO.CG", "THR.CG2", "LYS.CB", "LYS.CG", "LYS.CD", "ARG.CB", "ARG.CG", "CYS.SG",
        "GLU.CB", "GLN.CB", "ILE.CA", "VAL.CA", "LEU.CA",
    ];

    /// <summary>Profiles every ligand of the first model.</summary>
    public static IReadOnlyList<InteractionProfile> AnalyzeAll(Structure structure, InteractionOptions? options = null)
        => structure.Ligands.Where(l => l.Atoms.Count > 1).Select(l => Analyze(structure.FirstModel, l, options)).ToList();

    /// <summary>Finds the interactions between <paramref name="ligand"/> and the protein atoms of <paramref name="model"/>.</summary>
    public static InteractionProfile Analyze(Model model, Residue ligand, InteractionOptions? options = null)
    {
        options ??= new InteractionOptions();
        var protein = model.Atoms.Where(a => a.Residue is { IsAminoAcid: true } && a.Element != "H").ToList();
        var lig = ligand.Atoms.Where(a => a.Element != "H").ToList();
        var result = new List<Interaction>();
        if (lig.Count == 0 || protein.Count == 0) return new InteractionProfile(ligand, result);

        double reach = Math.Max(Math.Max(options.HBondMaxDistance, options.HydrophobicMaxDistance), options.IonicMaxDistance);
        var ns = NeighborSearch.ForAtoms(protein, reach);
        var neighbours = LigandBonds(lig);
        var charge = LigandCharges(lig, neighbours);
        var found = new List<int>();

        for (int li = 0; li < lig.Count; li++)
        {
            var la = lig[li];
            found.Clear();
            ns.Query(la.Position, reach, found);
            bool ligPolar = la.Element is "N" or "O";
            bool ligApolar = la.Element is "C" or "CL" or "BR" or "I" or "S" && neighbours[li].All(n => lig[n].Element is "C" or "H" or "CL" or "BR" or "I" or "S" or "F");
            Atom? bestHyd = null; double bestHydD = double.MaxValue;
            foreach (int pi in found.OrderBy(i => i))
            {
                var pa = protein[pi];
                double d = Vec3.Distance(la.Position, pa.Position);
                var (pDonor, pAcceptor) = PolarRole(pa);
                if (ligPolar && (pDonor || pAcceptor) && d >= options.HBondMinDistance && d <= options.HBondMaxDistance)
                {
                    // Ligand O and N may donate or accept; reject N(ligand, no H possible: 3 heavy neighbours, sp3/aromatic-tertiary) with protein acceptor only.
                    bool ligCanDonate = la.Element == "O" ? neighbours[li].Count <= 1 : neighbours[li].Count <= 2;
                    bool ligCanAccept = la.Element == "O" || neighbours[li].Count <= 2;
                    if ((pDonor && ligCanAccept) || (pAcceptor && ligCanDonate))
                        result.Add(new Interaction(InteractionType.HydrogenBond, la, pa, d, pDonor && ligCanAccept ? "protein donor" : "ligand donor"));
                }
                if (ligApolar && d <= options.HydrophobicMaxDistance && HydrophobicSideChain.Contains($"{pa.ResidueName}.{pa.Name}") && d < bestHydD)
                {
                    bestHyd = pa; bestHydD = d;
                }
                if (charge[li] != 0 && d <= options.IonicMaxDistance)
                {
                    int pc = ProteinCharge(pa, options.HistidineCharged);
                    if (pc != 0 && pc != charge[li])
                        result.Add(new Interaction(InteractionType.Ionic, la, pa, d, charge[li] > 0 ? "ligand cation" : "ligand anion"));
                }
            }
            if (bestHyd is not null) result.Add(new Interaction(InteractionType.Hydrophobic, la, bestHyd, bestHydD));
        }

        // π-stacking: ligand rings vs aromatic side chains.
        var ligRings = Rings(lig, neighbours).Where(r => IsPlanar(r.Select(i => lig[i].Position).ToArray())).ToList();
        if (ligRings.Count > 0)
        {
            var protRings = new List<Atom[]>();
            var nearRes = new HashSet<Residue>();
            foreach (var la in lig)
            {
                found.Clear();
                ns.Query(la.Position, options.PiTShapedMaxDistance + 1.5, found);
                foreach (int pi in found) nearRes.Add(protein[pi].Residue!);
            }
            foreach (var r in nearRes)
                if (AromaticRings.TryGetValue(r.Name, out var defs))
                    foreach (var def in defs)
                    {
                        var ring = def.Split(',').Select(r.GetAtom).ToArray();
                        if (ring.All(a => a is not null)) protRings.Add(ring!);
                    }
            foreach (var lr in ligRings)
            {
                var lp = lr.Select(i => lig[i].Position).ToArray();
                var (lc, ln) = CentroidNormal(lp);
                foreach (var pr in protRings)
                {
                    var (pc, pn) = CentroidNormal(pr.Select(a => a.Position).ToArray());
                    double d = Vec3.Distance(lc, pc);
                    double angle = Vec3.AngleBetween(ln, pn);
                    if (angle > 90) angle = 180 - angle;
                    string? kind = null;
                    if (d <= options.PiParallelMaxDistance && angle <= options.PiParallelMaxAngle)
                    {
                        // Offset check: the projection of the centroid offset onto the ring plane must be small (≤ 2.0 Å as in PLIP).
                        var v = pc - lc;
                        double offset = (v - ln * Vec3.Dot(v, ln)).Length;
                        if (offset <= 2.0 || d <= 4.4) kind = "parallel";
                    }
                    else if (d <= options.PiTShapedMaxDistance && angle >= 60) kind = "T-shaped";
                    if (kind is not null) result.Add(new Interaction(InteractionType.PiStacking, lig[lr[0]], pr[0], d, kind));
                }
            }
        }
        return new InteractionProfile(ligand, result);
    }

    private static (bool Donor, bool Acceptor) PolarRole(Atom a)
    {
        if (a.Element is not ("N" or "O")) return (false, false);
        return (a.ResidueName, a.Name) switch
        {
            (_, "N") => (a.ResidueName != "PRO", false),
            (_, "O" or "OXT") => (false, true),
            ("SER", "OG") or ("THR", "OG1") or ("TYR", "OH") => (true, true),
            ("ASP", _) or ("GLU", _) => (false, true),
            ("ASN", "OD1") or ("GLN", "OE1") => (false, true),
            ("ASN", "ND2") or ("GLN", "NE2") => (true, false),
            ("HIS", _) => (true, true),
            ("LYS", "NZ") or ("ARG", _) or ("TRP", "NE1") => (true, false),
            _ => (true, true),
        };
    }

    private static int ProteinCharge(Atom a, bool his) => (a.ResidueName, a.Name) switch
    {
        ("ASP", "OD1" or "OD2") or ("GLU", "OE1" or "OE2") => -1,
        ("LYS", "NZ") or ("ARG", "NE" or "NH1" or "NH2") => 1,
        ("HIS", "ND1" or "NE2") when his => 1,
        (_, "OXT") => -1,
        _ => 0,
    };

    /// <summary>Ligand bond graph from covalent radii (bond when d ≤ r₁ + r₂ + 0.45 Å).</summary>
    internal static List<int>[] LigandBonds(IReadOnlyList<Atom> atoms)
    {
        var nb = new List<int>[atoms.Count];
        for (int i = 0; i < nb.Length; i++) nb[i] = [];
        var ns = NeighborSearch.ForAtoms(atoms, 2.5);
        foreach (var (i, j) in ns.AllPairs(2.5))
        {
            double max = atoms[i].Element.Length > 0 ? ElementData.CovalentRadius(atoms[i].Element) + ElementData.CovalentRadius(atoms[j].Element) + 0.45 : 0;
            if (Vec3.Distance(atoms[i].Position, atoms[j].Position) <= max) { nb[i].Add(j); nb[j].Add(i); }
        }
        return nb;
    }

    /// <summary>Formal charges of ligand atoms from explicit charges or common ionisable groups at pH 7.</summary>
    internal static int[] LigandCharges(IReadOnlyList<Atom> atoms, List<int>[] nb)
    {
        var q = new int[atoms.Count];
        var rings = Rings(atoms, nb);
        var inRing = new bool[atoms.Count];
        foreach (var r in rings) foreach (int i in r) inRing[i] = true;
        for (int i = 0; i < atoms.Count; i++)
        {
            if (atoms[i].Charge != 0) { q[i] = Math.Sign(atoms[i].Charge); continue; }
            var a = atoms[i];
            var oTerminal = nb[i].Where(j => atoms[j].Element == "O" && nb[j].Count == 1).ToList();
            if (a.Element == "C" && oTerminal.Count == 2 && nb[i].Count == 3) foreach (int o in oTerminal) q[o] = -1;      // carboxylate
            if (a.Element is "P" or "S" && oTerminal.Count >= 3) foreach (int o in oTerminal) q[o] = -1;                 // phosphate/sulfonate
            if (a.Element == "C" && !inRing[i])
            {
                // Amidine / guanidine carbon: ≥ 2 acyclic nitrogens, at least one terminal, no oxygen (excludes ureas).
                var n = nb[i].Where(j => atoms[j].Element == "N" && !inRing[j]).ToList();
                if (n.Count >= 2 && !nb[i].Any(j => atoms[j].Element == "O") && n.Count(j => nb[j].Count == 1) >= 1)
                    foreach (int j in n) if (nb[j].Count <= 2) q[j] = 1;
            }
        }
        for (int i = 0; i < atoms.Count; i++)
        {
            // Aliphatic amine: N bonded only to sp3-looking carbons (no C=O/aromatic neighbours), not in ring-aromatic systems.
            var a = atoms[i];
            if (a.Element != "N" || q[i] != 0 || inRing[i] || nb[i].Count == 0 || nb[i].Count > 3) continue;
            bool aliphatic = nb[i].All(j => atoms[j].Element == "C" && !inRing[j] && nb[j].All(k => atoms[k].Element is "C" or "N" or "H") && nb[j].Count(k => atoms[k].Element == "N") == 1);
            if (aliphatic) q[i] = 1;
        }
        return q;
    }

    /// <summary>Simple cycles of 5 or 6 atoms in a small molecular graph.</summary>
    internal static List<int[]> Rings(IReadOnlyList<Atom> atoms, List<int>[] nb)
    {
        var rings = new List<int[]>();
        var seen = new HashSet<string>();
        var path = new List<int>();
        for (int s = 0; s < atoms.Count; s++)
        {
            path.Clear();
            path.Add(s);
            Dfs(s, s);
        }
        return rings;

        void Dfs(int start, int cur)
        {
            foreach (int next in nb[cur])
            {
                if (next == start && path.Count >= 5)
                {
                    var key = string.Join(",", path.OrderBy(x => x));
                    if (seen.Add(key)) rings.Add(path.ToArray());
                    continue;
                }
                if (next <= start || path.Contains(next) || path.Count >= 6) continue;
                path.Add(next);
                Dfs(start, next);
                path.RemoveAt(path.Count - 1);
            }
        }
    }

    private static (Vec3 Centroid, Vec3 Normal) CentroidNormal(Vec3[] ring)
    {
        var c = Geometry.Centroid(ring);
        var n = Vec3.Zero;
        for (int i = 0; i < ring.Length; i++)
            n += Vec3.Cross(ring[i] - c, ring[(i + 1) % ring.Length] - c);
        return (c, n.Normalized());
    }

    private static bool IsPlanar(Vec3[] ring, double tolerance = 0.25)
    {
        var (c, n) = CentroidNormal(ring);
        return ring.All(p => Math.Abs(Vec3.Dot(p - c, n)) <= tolerance);
    }
}
