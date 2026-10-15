namespace Biotech.Net.Chemistry;

/// <summary>Fingerprint families available through <see cref="Fingerprint.Compute"/>.</summary>
public enum FingerprintType
{
    /// <summary>Morgan radius 2 with connectivity invariants (ECFP4-like), 2048 bits.</summary>
    Ecfp4,
    /// <summary>Morgan radius 3 with connectivity invariants (ECFP6-like), 2048 bits.</summary>
    Ecfp6,
    /// <summary>Morgan radius 2 with pharmacophoric feature invariants (FCFP4-like), 2048 bits.</summary>
    Fcfp4,
    /// <summary>MACCS 166 public keys (167-bit vector, bit 0 unused, as in RDKit).</summary>
    Maccs,
    /// <summary>RDKit-style topological (branched path / subgraph) fingerprint, 1–7 bonds, 2048 bits.</summary>
    Path,
    /// <summary>Hashed atom-pair fingerprint, 2048 bits.</summary>
    AtomPair,
    /// <summary>Hashed topological-torsion fingerprint, 2048 bits.</summary>
    TopologicalTorsion,
}

/// <summary>One-stop access to the fingerprint families with their default parameters.</summary>
public static class Fingerprint
{
    /// <summary>Computes a fingerprint of the given <paramref name="type"/> with default settings.</summary>
    public static BitVector Compute(Molecule mol, FingerprintType type) => type switch
    {
        FingerprintType.Ecfp4 => MorganFingerprint.Compute(mol, 2),
        FingerprintType.Ecfp6 => MorganFingerprint.Compute(mol, 3),
        FingerprintType.Fcfp4 => MorganFingerprint.Compute(mol, 2, useFeatures: true),
        FingerprintType.Maccs => MaccsKeys.Compute(mol),
        FingerprintType.Path => PathFingerprint.Compute(mol),
        FingerprintType.AtomPair => AtomPairFingerprint.ComputeHashed(mol),
        FingerprintType.TopologicalTorsion => TopologicalTorsionFingerprint.ComputeHashed(mol),
        _ => throw new ArgumentOutOfRangeException(nameof(type)),
    };

    /// <summary>A fingerprint function for <paramref name="type"/>, e.g. for <see cref="FingerprintIndex"/>.</summary>
    public static Func<Molecule, BitVector> For(FingerprintType type) => m => Compute(m, type);
}

/// <summary>Where a Morgan feature came from: its centre atom and radius.</summary>
/// <param name="Atom">Centre atom index.</param>
/// <param name="Radius">Radius (0 = the atom alone).</param>
public readonly record struct MorganEnvironment(int Atom, int Radius);

/// <summary>
/// Morgan / extended-connectivity fingerprints (Rogers &amp; Hahn 2010) following RDKit's algorithm: the same
/// connectivity invariants (atomic number, total degree, total H, formal charge, isotope mass delta, ring
/// membership) or feature invariants (donor, acceptor, aromatic, halogen, basic, acidic), the same 32-bit
/// boost-style hash combination, bond types (single 1, double 2, triple 3, aromatic 12) and duplicate
/// environment removal. With default options the bit positions are intended to match RDKit's
/// <c>GetMorganFingerprintAsBitVect(mol, radius, nBits)</c>; chirality is not used.
/// </summary>
public static class MorganFingerprint
{
    private static readonly Lazy<SmartsPattern[]> _featurePatterns = new(() =>
    [
        SmartsPattern.Parse("[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]"),
        SmartsPattern.Parse("[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"),
        SmartsPattern.Parse("[a]"),
        SmartsPattern.Parse("[F,Cl,Br,I]"),
        SmartsPattern.Parse("[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]"),
        SmartsPattern.Parse("[$([C,S](=[O,S,P])-[O;H1,-1])]"),
    ]);

    /// <summary>Morgan bit vector (default: radius 2, 2048 bits = ECFP4-like).</summary>
    public static BitVector Compute(Molecule mol, int radius = 2, int nBits = 2048, bool useFeatures = false, bool useBondTypes = true)
    {
        var fp = new BitVector(nBits);
        Run(mol, radius, useFeatures, useBondTypes, (id, _, _) => fp.Set((int)(id % (uint)nBits)));
        return fp;
    }

    /// <summary>ECFP4-like fingerprint (radius 2).</summary>
    public static BitVector Ecfp4(Molecule mol, int nBits = 2048) => Compute(mol, 2, nBits);

    /// <summary>FCFP4-like fingerprint (radius 2, feature invariants).</summary>
    public static BitVector Fcfp4(Molecule mol, int nBits = 2048) => Compute(mol, 2, nBits, useFeatures: true);

    /// <summary>Unfolded Morgan count fingerprint keyed by the 32-bit environment identifiers (RDKit <c>GetMorganFingerprint</c>).</summary>
    public static CountVector ComputeCounts(Molecule mol, int radius = 2, bool useFeatures = false, bool useBondTypes = true)
    {
        var v = new CountVector();
        Run(mol, radius, useFeatures, useBondTypes, (id, _, _) => v.Add(id));
        return v;
    }

    /// <summary>Folded count fingerprint: environment <c>id</c> is counted at <c>id mod nBits</c>.</summary>
    public static CountVector ComputeHashedCounts(Molecule mol, int radius = 2, int nBits = 2048, bool useFeatures = false)
    {
        var v = new CountVector(nBits);
        Run(mol, radius, useFeatures, true, (id, _, _) => v.Add(id % (uint)nBits));
        return v;
    }

    /// <summary>The environments that set each bit of the folded fingerprint (RDKit <c>bitInfo</c>).</summary>
    public static IReadOnlyDictionary<int, IReadOnlyList<MorganEnvironment>> BitInfo(Molecule mol, int radius = 2, int nBits = 2048, bool useFeatures = false)
    {
        var info = new SortedDictionary<int, List<MorganEnvironment>>();
        Run(mol, radius, useFeatures, true, (id, atom, r) =>
        {
            int bit = (int)(id % (uint)nBits);
            if (!info.TryGetValue(bit, out var list)) info[bit] = list = [];
            list.Add(new MorganEnvironment(atom, r));
        });
        return info.ToDictionary(kv => kv.Key, kv => (IReadOnlyList<MorganEnvironment>)kv.Value);
    }

    /// <summary>RDKit connectivity invariants (ECFP atom identifiers before iteration).</summary>
    public static uint[] ConnectivityInvariants(Molecule mol, bool includeRingMembership = true)
    {
        int n = mol.AtomCount;
        var inv = new uint[n];
        var rings = mol.Rings;
        Span<uint> comp = stackalloc uint[6];
        for (int i = 0; i < n; i++)
        {
            var a = mol.Atoms[i];
            int k = 0;
            comp[k++] = (uint)a.AtomicNumber;
            comp[k++] = (uint)(mol.Degree(i) + a.HydrogenCount);
            comp[k++] = (uint)mol.TotalHydrogens(i);
            comp[k++] = unchecked((uint)a.FormalCharge);
            int deltaMass = 0;
            if (a.Isotope != 0 && a.AtomicNumber > 0)
                deltaMass = (int)(PeriodicTable.IsotopeMass(a.AtomicNumber, a.Isotope) - a.Element.AverageMass);
            comp[k++] = unchecked((uint)deltaMass);
            if (includeRingMembership && rings.AtomRingCount(i) > 0) comp[k++] = 1;
            uint seed = 0;
            for (int j = 0; j < k; j++) HashCombine(ref seed, comp[j]);
            inv[i] = seed;
        }
        return inv;
    }

    /// <summary>RDKit feature invariants (FCFP): bit 0 donor, 1 acceptor, 2 aromatic, 3 halogen, 4 basic, 5 acidic.</summary>
    public static uint[] FeatureInvariants(Molecule mol)
    {
        var inv = new uint[mol.AtomCount];
        var pats = _featurePatterns.Value;
        for (int j = 0; j < pats.Length; j++)
            foreach (int a in pats[j].MatchingAtoms(mol)) inv[a] |= 1u << j;
        return inv;
    }

    /// <summary>boost::hash_combine with 32-bit seeds (as used by RDKit's <c>gboost</c>).</summary>
    internal static void HashCombine(ref uint seed, uint value)
        => seed ^= unchecked(value + 0x9e3779b9u + (seed << 6) + (seed >> 2));

    private static void Run(Molecule mol, int radius, bool useFeatures, bool useBondTypes, Action<uint, int, int> emit)
    {
        int n = mol.AtomCount, m = mol.BondCount;
        if (n == 0) return;
        var invariants = useFeatures ? FeatureInvariants(mol) : ConnectivityInvariants(mol);
        var current = (uint[])invariants.Clone();
        for (int i = 0; i < n; i++) emit(invariants[i], i, 0);

        int words = Math.Max(1, (m + 63) >> 6);
        var atomNbhd = new ulong[n][];
        for (int i = 0; i < n; i++) atomNbhd[i] = new ulong[words];
        var seen = new HashSet<NeighborhoodKey>();
        var dead = new bool[n];
        var bonds = mol.Bonds;
        var nbrs = new List<(int Bt, uint Inv)>();

        for (int layer = 0; layer < radius; layer++)
        {
            var roundInv = new uint[n];
            var roundNbhd = new ulong[n][];
            for (int i = 0; i < n; i++) roundNbhd[i] = (ulong[])atomNbhd[i].Clone();
            var thisRound = new List<(ulong[] Nbhd, uint Inv, int Atom)>();
            for (int i = 0; i < n; i++)
            {
                if (dead[i]) continue;
                if (mol.Degree(i) == 0) { dead[i] = true; continue; }
                nbrs.Clear();
                var nb = mol.Neighbors(i);
                var nbb = mol.NeighborBonds(i);
                for (int k = 0; k < nb.Length; k++)
                {
                    int b = nbb[k], o = nb[k];
                    roundNbhd[i][b >> 6] |= 1UL << (b & 63);
                    var on = atomNbhd[o];
                    for (int w = 0; w < words; w++) roundNbhd[i][w] |= on[w];
                    int bt = useBondTypes ? bonds[b].Order switch
                    {
                        BondOrder.Double => 2,
                        BondOrder.Triple => 3,
                        BondOrder.Aromatic => 12,
                        _ => 1,
                    } : 1;
                    nbrs.Add((bt, current[o]));
                }
                nbrs.Sort((x, y) => x.Bt != y.Bt ? x.Bt.CompareTo(y.Bt) : x.Inv.CompareTo(y.Inv));
                uint invar = (uint)layer;
                HashCombine(ref invar, current[i]);
                foreach (var (bt, ni) in nbrs)
                {
                    uint pair = 0;
                    HashCombine(ref pair, unchecked((uint)bt));
                    HashCombine(ref pair, ni);
                    HashCombine(ref invar, pair);
                }
                roundInv[i] = invar;
                thisRound.Add((roundNbhd[i], invar, i));
            }
            thisRound.Sort((x, y) =>
            {
                int c = CompareBitsets(x.Nbhd, y.Nbhd);
                if (c != 0) return c;
                c = x.Inv.CompareTo(y.Inv);
                return c != 0 ? c : x.Atom.CompareTo(y.Atom);
            });
            foreach (var (nbhd, inv, atom) in thisRound)
            {
                if (seen.Add(new NeighborhoodKey(nbhd))) emit(inv, atom, layer + 1);
                else dead[atom] = true;
            }
            current = roundInv;
            atomNbhd = roundNbhd;
        }
    }

    // boost::dynamic_bitset operator< compares from the most significant bit.
    private static int CompareBitsets(ulong[] a, ulong[] b)
    {
        for (int w = a.Length - 1; w >= 0; w--)
            if (a[w] != b[w]) return a[w].CompareTo(b[w]);
        return 0;
    }

    private readonly struct NeighborhoodKey(ulong[] bits) : IEquatable<NeighborhoodKey>
    {
        private readonly ulong[] _bits = bits;

        public bool Equals(NeighborhoodKey other) => _bits.AsSpan().SequenceEqual(other._bits);

        public override bool Equals(object? obj) => obj is NeighborhoodKey k && Equals(k);

        public override int GetHashCode()
        {
            var h = new HashCode();
            foreach (ulong w in _bits) h.Add(w);
            return h.ToHashCode();
        }
    }
}

/// <summary>
/// MACCS 166 structural keys (Durant et al. 2002) using the SMARTS definitions distributed with RDKit. The
/// result has 167 bits with bit 0 unused, so bit <c>k</c> is key <c>k</c>, matching RDKit's
/// <c>MACCSkeys.GenMACCSKeys</c>.
/// </summary>
public static class MaccsKeys
{
    private static readonly Lazy<(int Key, SmartsPattern? Pattern, int Count)[]> _keys = new(() =>
        MaccsData.Keys.Select(k => (k.Key, k.Smarts == "?" ? null : SmartsPattern.Parse(k.Smarts), k.Count)).ToArray());

    /// <summary>The SMARTS definition of key <paramref name="key"/> (1–166), "?" for keys computed specially.</summary>
    public static string Definition(int key) => MaccsData.Keys[key - 1].Smarts;

    /// <summary>Computes the 167-bit MACCS fingerprint.</summary>
    public static BitVector Compute(Molecule mol)
    {
        var fp = new BitVector(167);
        var ctx = new SmartsPattern.MatchContext(mol);
        foreach (var (key, pattern, count) in _keys.Value)
        {
            if (pattern is null) continue;
            bool on = count == 0
                ? pattern.FindMatches(mol, false, 1, ctx, -1).Count > 0
                : pattern.FindMatches(mol, true, count + 1, ctx, -1).Count > count;
            if (on) fp.Set(key);
        }
        // 125: more than one aromatic ring; 166: more than one fragment.
        if (MolecularDescriptors.AromaticRingCount(mol) > 1) fp.Set(125);
        if (mol.Fragments().Count > 1) fp.Set(166);
        return fp;
    }
}

/// <summary>
/// RDKit-style topological fingerprint (Daylight-like): every connected subgraph of
/// <see cref="Options.MinPath"/>–<see cref="Options.MaxPath"/> bonds (or only linear paths) is hashed from its
/// bond descriptors (bond order, end-atom element/aromaticity, and the end atoms' degrees within the subgraph) and
/// sets <see cref="Options.BitsPerHash"/> pseudo-random bits. It follows the design of RDKit's
/// <c>RDKFingerprint</c> but is not bit-identical to it.
/// </summary>
public static class PathFingerprint
{
    /// <summary>Path fingerprint parameters (RDKit defaults).</summary>
    public sealed record Options
    {
        /// <summary>Smallest subgraph size in bonds (default 1).</summary>
        public int MinPath { get; init; } = 1;
        /// <summary>Largest subgraph size in bonds (default 7).</summary>
        public int MaxPath { get; init; } = 7;
        /// <summary>Fingerprint length (default 2048).</summary>
        public int Size { get; init; } = 2048;
        /// <summary>Bits set per subgraph (default 2).</summary>
        public int BitsPerHash { get; init; } = 2;
        /// <summary>Enumerate branched subgraphs (true, default) or only linear paths.</summary>
        public bool BranchedPaths { get; init; } = true;
        /// <summary>Include bond orders in the hash (default true).</summary>
        public bool UseBondOrder { get; init; } = true;
        /// <summary>Safety cap on the number of subgraphs enumerated (default 200 000).</summary>
        public int MaxSubgraphs { get; init; } = 200_000;
    }

    /// <summary>Computes the fingerprint.</summary>
    public static BitVector Compute(Molecule mol, Options? options = null)
    {
        options ??= new Options();
        var fp = new BitVector(options.Size);
        var bonds = mol.Bonds;
        var atomInv = new uint[mol.AtomCount];
        for (int i = 0; i < atomInv.Length; i++)
        {
            var a = mol.Atoms[i];
            atomInv[i] = (uint)(a.AtomicNumber % 128) | (a.IsAromatic ? 128u : 0u);
        }
        var degree = new int[mol.AtomCount];
        var descriptors = new List<ulong>();
        foreach (var sub in Subgraphs.Enumerate(mol, options.MinPath, options.MaxPath, options.BranchedPaths, options.MaxSubgraphs))
        {
            foreach (int b in sub) { degree[bonds[b].Begin]++; degree[bonds[b].End]++; }
            descriptors.Clear();
            foreach (int b in sub)
            {
                var bond = bonds[b];
                uint x = atomInv[bond.Begin] | ((uint)degree[bond.Begin] << 8);
                uint y = atomInv[bond.End] | ((uint)degree[bond.End] << 8);
                if (x > y) (x, y) = (y, x);
                uint order = options.UseBondOrder ? (uint)bond.Order : 1u;
                descriptors.Add(((ulong)order << 48) | ((ulong)x << 24) | y);
            }
            descriptors.Sort();
            uint seed = (uint)sub.Count;
            foreach (ulong d in descriptors)
            {
                MorganFingerprint.HashCombine(ref seed, (uint)d);
                MorganFingerprint.HashCombine(ref seed, (uint)(d >> 32));
            }
            foreach (int b in sub) { degree[bonds[b].Begin]--; degree[bonds[b].End]--; }
            // minstd_rand-style generator seeded by the hash picks the bits.
            ulong state = seed == 0 ? 1u : seed;
            for (int k = 0; k < options.BitsPerHash; k++)
            {
                state = state * 48271 % 2147483647;
                fp.Set((int)(state % (ulong)options.Size));
            }
        }
        return fp;
    }
}

/// <summary>Enumeration of connected bond subgraphs (each subgraph exactly once, ESU-style extension).</summary>
public static class Subgraphs
{
    /// <summary>
    /// All connected subgraphs with <paramref name="minBonds"/>…<paramref name="maxBonds"/> bonds, as bond-index lists.
    /// With <paramref name="branched"/> false only linear paths (every atom of degree ≤ 2 in the subgraph, no cycles) are returned.
    /// </summary>
    public static IEnumerable<IReadOnlyList<int>> Enumerate(Molecule mol, int minBonds, int maxBonds, bool branched = true, int limit = int.MaxValue)
    {
        var result = new List<IReadOnlyList<int>>();
        int m = mol.BondCount;
        var bonds = mol.Bonds.ToArray();
        // bond adjacency
        var adj = new List<int>[m];
        for (int b = 0; b < m; b++)
        {
            var set = new List<int>();
            foreach (int e in mol.NeighborBonds(bonds[b].Begin)) if (e != b) set.Add(e);
            foreach (int e in mol.NeighborBonds(bonds[b].End)) if (e != b && !set.Contains(e)) set.Add(e);
            adj[b] = set;
        }
        var sub = new List<int>();
        var inSub = new bool[m];
        var inNbhd = new int[m];
        var atomDeg = new int[mol.AtomCount];
        int count = 0;
        for (int start = 0; start < m && count < limit; start++)
        {
            sub.Add(start);
            inSub[start] = true;
            atomDeg[bonds[start].Begin]++; atomDeg[bonds[start].End]++;
            var ext = adj[start].Where(e => e > start).ToList();
            Mark(start, +1);
            Extend(ext, start);
            Mark(start, -1);
            atomDeg[bonds[start].Begin]--; atomDeg[bonds[start].End]--;
            inSub[start] = false;
            sub.RemoveAt(0);
        }
        return result;

        void Mark(int b, int delta)
        {
            foreach (int e in adj[b]) inNbhd[e] += delta;
        }

        void Extend(List<int> ext, int root)
        {
            if (sub.Count >= minBonds && (branched || IsLinear()))
            {
                if (count++ >= limit) return;
                result.Add(sub.ToArray());
            }
            if (sub.Count == maxBonds) return;
            var work = new List<int>(ext);
            while (work.Count > 0)
            {
                int w = work[^1];
                work.RemoveAt(work.Count - 1);
                if (!branched)
                {
                    // linear: w must attach at an atom of subgraph-degree 1 and not close a cycle
                    var bw = bonds[w];
                    int da = atomDeg[bw.Begin], db = atomDeg[bw.End];
                    if (da >= 2 || db >= 2 || (da >= 1 && db >= 1)) continue;
                }
                // exclusive neighbourhood of w
                var newExt = new List<int>(work);
                foreach (int u in adj[w])
                    if (u > root && !inSub[u] && inNbhd[u] == 0 && !newExt.Contains(u)) newExt.Add(u);
                sub.Add(w);
                inSub[w] = true;
                atomDeg[bonds[w].Begin]++; atomDeg[bonds[w].End]++;
                Mark(w, +1);
                Extend(newExt, root);
                Mark(w, -1);
                atomDeg[bonds[w].Begin]--; atomDeg[bonds[w].End]--;
                inSub[w] = false;
                sub.RemoveAt(sub.Count - 1);
                if (count >= limit) return;
            }
        }

        bool IsLinear()
        {
            // a tree in which every atom has degree ≤ 2 → path
            int atoms = 0;
            foreach (int b in sub)
            {
                if (atomDeg[bonds[b].Begin] > 2 || atomDeg[bonds[b].End] > 2) return false;
            }
            var seen = new HashSet<int>();
            foreach (int b in sub) { seen.Add(bonds[b].Begin); seen.Add(bonds[b].End); }
            atoms = seen.Count;
            return atoms == sub.Count + 1;
        }
    }
}

/// <summary>
/// Atom-pair fingerprints (Carhart, Smith &amp; Venkataraghavan 1985) with RDKit's atom codes: number of heavy
/// neighbours (3 bits), π electrons (2 bits) and an element class (4 bits), combined with the topological
/// distance (5 bits) into a 23-bit pair code.
/// </summary>
public static class AtomPairFingerprint
{
    private static readonly int[] AtomNumberTypes = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 43, 0];
    internal const int CodeSize = 9;
    private const int PathBits = 5;

    /// <summary>The RDKit atom code of an atom (<paramref name="branchSubtract"/> is subtracted from the degree).</summary>
    public static uint AtomCode(Molecule mol, int atom, int branchSubtract = 0)
    {
        int degree = mol.Degree(atom);
        int branches = degree > branchSubtract ? degree - branchSubtract : 0;
        uint code = (uint)(branches % 7);
        uint pi = (uint)(PiElectrons(mol, atom) % 3);
        code |= pi << 3;
        int z = mol.Atoms[atom].AtomicNumber;
        int typeIdx = 0;
        int nTypes = AtomNumberTypes.Length;
        while (typeIdx < nTypes)
        {
            if (AtomNumberTypes[typeIdx] == z) break;
            if (AtomNumberTypes[typeIdx] > z) { typeIdx = nTypes; break; }
            typeIdx++;
        }
        if (typeIdx == nTypes) typeIdx--;
        code |= (uint)typeIdx << 5;
        return code;
    }

    internal static int PiElectrons(Molecule mol, int atom)
    {
        if (mol.Atoms[atom].IsAromatic) return 1;
        if (mol.Degree(atom) + mol.Atoms[atom].HydrogenCount >= 4) return 0; // sp3 (e.g. sulfonyl S, phosphate P)
        int pi = 0;
        foreach (int b in mol.NeighborBonds(atom))
        {
            int o = (int)mol.Bonds[b].Order;
            if (o is 2 or 3) pi += o - 1;
        }
        return pi;
    }

    /// <summary>The 23-bit pair code for two atom codes at topological distance <paramref name="distance"/>.</summary>
    public static uint PairCode(uint codeI, uint codeJ, int distance)
    {
        uint d = (uint)Math.Min(distance, (1 << PathBits) - 1);
        return d | (Math.Min(codeI, codeJ) << PathBits) | (Math.Max(codeI, codeJ) << (PathBits + CodeSize));
    }

    /// <summary>Unfolded atom-pair count fingerprint keyed by pair code (all graph atoms, distances 1–30; matches RDKit <c>GetAtomPairFingerprint</c>).</summary>
    public static CountVector Compute(Molecule mol, int minDistance = 1, int maxDistance = 30)
    {
        var v = new CountVector(1L << 23);
        var heavy = Enumerable.Range(0, mol.AtomCount).ToArray();
        var codes = new uint[mol.AtomCount];
        foreach (int i in heavy) codes[i] = AtomCode(mol, i);
        var d = mol.TopologicalDistances();
        for (int x = 0; x < heavy.Length; x++)
            for (int y = x + 1; y < heavy.Length; y++)
            {
                int i = heavy[x], j = heavy[y];
                int dist = d[i, j];
                if (dist < minDistance || dist > maxDistance) continue;
                v.Add(PairCode(codes[i], codes[j], dist));
            }
        return v;
    }

    /// <summary>
    /// Hashed atom-pair bit vector with count simulation: each pair hashes to a block of 4 bits, and bit k of the
    /// block is set when the pair occurs more than {0, 1, 3, 7}[k] times (RDKit's default <c>nBitsPerEntry = 4</c>).
    /// </summary>
    public static BitVector ComputeHashed(Molecule mol, int nBits = 2048, int minDistance = 1, int maxDistance = 30)
        => CountSimulation(HashedCounts(Compute(mol, minDistance, maxDistance), nBits / 4), nBits);

    internal static Dictionary<int, int> HashedCounts(CountVector v, int blocks)
    {
        var r = new Dictionary<int, int>();
        foreach (var (id, c) in v.Counts)
        {
            uint seed = 0;
            MorganFingerprint.HashCombine(ref seed, (uint)id);
            MorganFingerprint.HashCombine(ref seed, (uint)(id >> 32));
            int block = (int)(seed % (uint)blocks);
            r[block] = r.GetValueOrDefault(block) + c;
        }
        return r;
    }

    internal static BitVector CountSimulation(Dictionary<int, int> blockCounts, int nBits)
    {
        int[] thresholds = [1, 2, 4, 8];
        var fp = new BitVector(nBits);
        foreach (var (block, c) in blockCounts)
            for (int k = 0; k < 4; k++)
                if (c >= thresholds[k] && block * 4 + k < nBits) fp.Set(block * 4 + k);
        return fp;
    }
}

/// <summary>
/// Topological-torsion fingerprints (Nilakantan et al. 1987): linear paths of four heavy atoms encoded by RDKit-style
/// atom codes (end atoms with one branch subtracted, inner atoms with two), canonicalised by direction.
/// </summary>
public static class TopologicalTorsionFingerprint
{
    /// <summary>Unfolded torsion count fingerprint.</summary>
    public static CountVector Compute(Molecule mol, int pathLength = 4)
    {
        var v = new CountVector(long.MaxValue);
        int n = mol.AtomCount;
        var path = new List<int>();
        var onPath = new bool[n];
        var heavy = new bool[n];
        for (int i = 0; i < n; i++) heavy[i] = true;
        for (int s = 0; s < n; s++)
        {
            if (!heavy[s]) continue;
            path.Add(s); onPath[s] = true;
            Walk();
            path.Clear(); onPath[s] = false;
        }
        return v;

        void Walk()
        {
            if (path.Count == pathLength)
            {
                if (path[0] > path[^1]) return; // each path once
                var codes = new uint[pathLength];
                for (int i = 0; i < pathLength; i++)
                    codes[i] = AtomPairFingerprint.AtomCode(mol, path[i], i == 0 || i == pathLength - 1 ? 1 : 2);
                bool reverse = false;
                for (int i = 0; i < pathLength / 2; i++)
                {
                    uint x = codes[i], y = codes[pathLength - 1 - i];
                    if (x != y) { reverse = x > y; break; }
                }
                if (reverse) Array.Reverse(codes);
                long code = 0;
                for (int i = 0; i < pathLength; i++) code |= (long)codes[i] << (i * AtomPairFingerprint.CodeSize);
                v.Add(code);
                return;
            }
            int last = path[^1];
            foreach (int nb in mol.Neighbors(last))
            {
                if (onPath[nb] || !heavy[nb]) continue;
                path.Add(nb); onPath[nb] = true;
                Walk();
                path.RemoveAt(path.Count - 1); onPath[nb] = false;
            }
        }
    }

    /// <summary>Hashed torsion bit vector with 4-bit count simulation.</summary>
    public static BitVector ComputeHashed(Molecule mol, int nBits = 2048)
        => AtomPairFingerprint.CountSimulation(AtomPairFingerprint.HashedCounts(Compute(mol), nBits / 4), nBits);
}
