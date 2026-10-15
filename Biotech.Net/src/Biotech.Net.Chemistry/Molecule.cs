using System.Runtime.CompilerServices;

namespace Biotech.Net.Chemistry;

/// <summary>Bond multiplicity. Aromatic bonds keep a Kekulé assignment internally (see <see cref="Molecule.KekuleOrder"/>).</summary>
public enum BondOrder : byte
{
    /// <summary>Single bond.</summary>
    Single = 1,
    /// <summary>Double bond.</summary>
    Double = 2,
    /// <summary>Triple bond.</summary>
    Triple = 3,
    /// <summary>Aromatic (delocalised) bond.</summary>
    Aromatic = 4,
}

/// <summary>Stereo annotation on a bond.</summary>
public enum BondStereo : byte
{
    /// <summary>No stereo information.</summary>
    None,
    /// <summary>Wedge (MOL "up"), relative to the begin atom.</summary>
    Up,
    /// <summary>Hash (MOL "down"), relative to the begin atom.</summary>
    Down,
    /// <summary>Unknown/either configuration.</summary>
    Either,
    /// <summary>Double bond whose reference atoms (<see cref="Bond.StereoAtomBegin"/>, <see cref="Bond.StereoAtomEnd"/>) are cis (Z-like).</summary>
    Cis,
    /// <summary>Double bond whose reference atoms are trans (E-like).</summary>
    Trans,
}

/// <summary>SMILES tetrahedral chirality: the order of the remaining neighbours seen from the first neighbour.</summary>
public enum ChiralTag : byte
{
    /// <summary>Not a stereocentre (or unspecified).</summary>
    None,
    /// <summary><c>@</c>: anticlockwise.</summary>
    CounterClockwise,
    /// <summary><c>@@</c>: clockwise.</summary>
    Clockwise,
}

/// <summary>A 3-D point (2-D depictions use <see cref="Z"/> = 0).</summary>
public readonly record struct Point3(double X, double Y, double Z = 0)
{
    public static Point3 operator +(Point3 a, Point3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    public static Point3 operator -(Point3 a, Point3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    public static Point3 operator *(Point3 a, double s) => new(a.X * s, a.Y * s, a.Z * s);
    public static Point3 operator *(double s, Point3 a) => new(a.X * s, a.Y * s, a.Z * s);

    /// <summary>Euclidean length.</summary>
    public double Length => Math.Sqrt(X * X + Y * Y + Z * Z);

    /// <summary>Distance between two points.</summary>
    public static double Distance(Point3 a, Point3 b) => (a - b).Length;

    /// <summary>Dot product.</summary>
    public static double Dot(Point3 a, Point3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

    /// <summary>Cross product.</summary>
    public static Point3 Cross(Point3 a, Point3 b) => new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
}

/// <summary>
/// One atom of a <see cref="Molecule"/>. Atoms are value types stored contiguously; mutate them through
/// <see cref="Molecule.AtomAt"/> (which returns a reference).
/// </summary>
public struct Atom
{
    /// <summary>Creates an atom of element <paramref name="atomicNumber"/> (0 = dummy/attachment point <c>*</c>).</summary>
    public Atom(int atomicNumber, int formalCharge = 0)
    {
        AtomicNumber = atomicNumber;
        FormalCharge = formalCharge;
    }

    /// <summary>Atomic number (0 = dummy atom <c>*</c>).</summary>
    public int AtomicNumber { get; set; }

    /// <summary>Formal charge.</summary>
    public int FormalCharge { get; set; }

    /// <summary>Mass number for isotopically labelled atoms; 0 means natural abundance.</summary>
    public int Isotope { get; set; }

    /// <summary>Hydrogens implied by valence rules (recomputed by <see cref="Molecule.UpdateImplicitHydrogens"/>).</summary>
    public int ImplicitHydrogens { get; set; }

    /// <summary>Hydrogens stated explicitly but not stored as atoms (bracket <c>[NH4+]</c>, MOL H counts).</summary>
    public int ExplicitHydrogens { get; set; }

    /// <summary>Aromatic flag.</summary>
    public bool IsAromatic { get; set; }

    /// <summary>When true, implicit hydrogens are never added (SMILES bracket atoms).</summary>
    public bool NoImplicit { get; set; }

    /// <summary>Tetrahedral chirality as written (relative to <see cref="Molecule.GetChiralReference"/>).</summary>
    public ChiralTag Chirality { get; set; }

    /// <summary>Atom-map / R-group label (<c>[*:1]</c>), 0 when unset.</summary>
    public int MapNumber { get; set; }

    /// <summary>Hydrogens carried by the atom itself (implicit + explicit counts; excludes H neighbour atoms).</summary>
    public readonly int HydrogenCount => ImplicitHydrogens + ExplicitHydrogens;

    /// <summary>The element.</summary>
    public readonly Element Element => PeriodicTable.Get(AtomicNumber);

    /// <summary>Element symbol (<c>*</c> for a dummy atom).</summary>
    public readonly string Symbol => PeriodicTable.Get(AtomicNumber).Symbol;

    public override readonly string ToString()
        => $"{(IsAromatic ? Symbol.ToLowerInvariant() : Symbol)}{(HydrogenCount > 0 ? $"H{HydrogenCount}" : "")}{(FormalCharge > 0 ? "+" + FormalCharge : FormalCharge < 0 ? FormalCharge.ToString() : "")}";
}

/// <summary>A bond between two atoms (by index).</summary>
public struct Bond
{
    /// <summary>Creates a bond.</summary>
    public Bond(int begin, int end, BondOrder order)
    {
        Begin = begin;
        End = end;
        Order = order;
        StereoAtomBegin = -1;
        StereoAtomEnd = -1;
    }

    /// <summary>First atom index.</summary>
    public int Begin { get; set; }

    /// <summary>Second atom index.</summary>
    public int End { get; set; }

    /// <summary>Bond order.</summary>
    public BondOrder Order { get; set; }

    /// <summary>Stereo annotation.</summary>
    public BondStereo Stereo { get; set; }

    /// <summary>For cis/trans double bonds: the reference neighbour of <see cref="Begin"/>.</summary>
    public int StereoAtomBegin { get; set; }

    /// <summary>For cis/trans double bonds: the reference neighbour of <see cref="End"/>.</summary>
    public int StereoAtomEnd { get; set; }

    /// <summary>True for aromatic bonds.</summary>
    public readonly bool IsAromatic => Order == BondOrder.Aromatic;

    /// <summary>The atom at the other end from <paramref name="atom"/>.</summary>
    public readonly int Other(int atom) => atom == Begin ? End : Begin;

    /// <summary>Whether the bond touches <paramref name="atom"/>.</summary>
    public readonly bool Contains(int atom) => Begin == atom || End == atom;

    /// <summary>Bond order as a number (aromatic = 1.5).</summary>
    public readonly double Valence => Order == BondOrder.Aromatic ? 1.5 : (int)Order;

    public override readonly string ToString() => $"{Begin}{Order switch { BondOrder.Double => "=", BondOrder.Triple => "#", BondOrder.Aromatic => ":", _ => "-" }}{End}";
}

/// <summary>
/// A molecular graph: atoms (struct array), bonds (struct array) and a lazily built compressed adjacency
/// index. Ring information and Kekulé assignments are derived on demand and cached until the graph changes.
/// </summary>
/// <remarks>
/// Hydrogens are normally implicit (counts on heavy atoms). Use <see cref="AddHydrogens"/> to make them explicit
/// atoms (needed for 3-D work and for H-typed atom contributions such as Crippen logP).
/// </remarks>
public sealed class Molecule
{
    private Atom[] _atoms;
    private Bond[] _bonds;
    private int _atomCount;
    private int _bondCount;
    private int[]? _adjStart, _adjAtom, _adjBond;
    private RingInfo? _rings;
    private byte[]? _kekule;
    private Dictionary<int, int[]>? _chiralRefs;
    private Point3[]? _coords;
    private Dictionary<string, string>? _properties;

    /// <summary>Creates an empty molecule.</summary>
    public Molecule(string? name = null, int atomCapacity = 16)
    {
        _atoms = new Atom[Math.Max(4, atomCapacity)];
        _bonds = new Bond[Math.Max(4, atomCapacity + 2)];
        Name = name ?? string.Empty;
    }

    /// <summary>Parses SMILES (shortcut for <see cref="SmilesParser.Parse"/>).</summary>
    public static Molecule FromSmiles(string smiles, string? name = null) => SmilesParser.Parse(smiles, name);

    /// <summary>Molecule title (SDF first line, SMILES name field, ...).</summary>
    public string Name { get; set; }

    /// <summary>Data fields (SDF <c>&gt; &lt;tag&gt;</c> items) and other annotations.</summary>
    public Dictionary<string, string> Properties => _properties ??= new(StringComparer.OrdinalIgnoreCase);

    /// <summary>Number of atoms (including explicit hydrogen atoms).</summary>
    public int AtomCount => _atomCount;

    /// <summary>Number of bonds.</summary>
    public int BondCount => _bondCount;

    /// <summary>All atoms, without copying.</summary>
    public ReadOnlySpan<Atom> Atoms => _atoms.AsSpan(0, _atomCount);

    /// <summary>All bonds, without copying.</summary>
    public ReadOnlySpan<Bond> Bonds => _bonds.AsSpan(0, _bondCount);

    /// <summary>A mutable reference to atom <paramref name="index"/>.</summary>
    public ref Atom AtomAt(int index)
    {
        if ((uint)index >= (uint)_atomCount) throw new ArgumentOutOfRangeException(nameof(index));
        return ref _atoms[index];
    }

    /// <summary>A mutable reference to bond <paramref name="index"/>. Changing endpoints requires <see cref="InvalidateTopology"/>.</summary>
    public ref Bond BondAt(int index)
    {
        if ((uint)index >= (uint)_bondCount) throw new ArgumentOutOfRangeException(nameof(index));
        return ref _bonds[index];
    }

    /// <summary>Atom coordinates (2-D or 3-D), or null when the molecule has none.</summary>
    public Point3[]? Coordinates
    {
        get => _coords;
        set
        {
            if (value is not null && value.Length != _atomCount)
                throw new ArgumentException($"Expected {_atomCount} coordinates, got {value.Length}.", nameof(value));
            _coords = value;
        }
    }

    /// <summary>True when coordinates exist and at least one atom has a non-zero Z.</summary>
    public bool Has3DCoordinates => _coords is not null && _coords.Any(p => Math.Abs(p.Z) > 1e-4);

    /// <summary>Number of heavy (non-hydrogen) atoms.</summary>
    public int HeavyAtomCount
    {
        get
        {
            int n = 0;
            foreach (ref readonly var a in Atoms) if (a.AtomicNumber != 1) n++;
            return n;
        }
    }

    // ---- Editing ---------------------------------------------------------------------------------

    /// <summary>Appends an atom and returns its index.</summary>
    public int AddAtom(Atom atom)
    {
        if (_atomCount == _atoms.Length) Array.Resize(ref _atoms, _atoms.Length * 2);
        _atoms[_atomCount] = atom;
        if (_coords is not null) Array.Resize(ref _coords, _atomCount + 1);
        InvalidateTopology();
        return _atomCount++;
    }

    /// <summary>Appends an atom of element <paramref name="atomicNumber"/> and returns its index.</summary>
    public int AddAtom(int atomicNumber, int formalCharge = 0, bool aromatic = false)
        => AddAtom(new Atom(atomicNumber, formalCharge) { IsAromatic = aromatic });

    /// <summary>Appends an atom by symbol (e.g. <c>"Cl"</c>) and returns its index.</summary>
    public int AddAtom(string symbol, int formalCharge = 0) => AddAtom(PeriodicTable.Get(symbol).AtomicNumber, formalCharge);

    /// <summary>Adds a bond and returns its index.</summary>
    /// <exception cref="ArgumentException">The atoms are the same, out of range, or already bonded.</exception>
    public int AddBond(int begin, int end, BondOrder order = BondOrder.Single)
    {
        if ((uint)begin >= (uint)_atomCount || (uint)end >= (uint)_atomCount) throw new ArgumentOutOfRangeException(nameof(begin));
        if (begin == end) throw new ArgumentException("An atom cannot bond to itself.");
        if (BondBetween(begin, end) >= 0) throw new ArgumentException($"Atoms {begin} and {end} are already bonded.");
        if (_bondCount == _bonds.Length) Array.Resize(ref _bonds, _bonds.Length * 2);
        _bonds[_bondCount] = new Bond(begin, end, order);
        InvalidateTopology();
        return _bondCount++;
    }

    /// <summary>Removes atoms (and their bonds), renumbering the rest. Returns the old→new index map (−1 for removed).</summary>
    public int[] RemoveAtoms(IEnumerable<int> atoms)
    {
        var remove = new bool[_atomCount];
        foreach (int a in atoms) remove[a] = true;
        var map = new int[_atomCount];
        int n = 0;
        for (int i = 0; i < _atomCount; i++)
        {
            if (remove[i]) { map[i] = -1; continue; }
            map[i] = n;
            _atoms[n] = _atoms[i];
            if (_coords is not null) _coords[n] = _coords[i];
            n++;
        }
        int m = 0;
        for (int i = 0; i < _bondCount; i++)
        {
            var b = _bonds[i];
            if (remove[b.Begin] || remove[b.End]) continue;
            b.Begin = map[b.Begin];
            b.End = map[b.End];
            if (b.StereoAtomBegin >= 0) b.StereoAtomBegin = map[b.StereoAtomBegin];
            if (b.StereoAtomEnd >= 0) b.StereoAtomEnd = map[b.StereoAtomEnd];
            if (b.StereoAtomBegin < 0 || b.StereoAtomEnd < 0) { if (b.Stereo is BondStereo.Cis or BondStereo.Trans) b.Stereo = BondStereo.None; b.StereoAtomBegin = b.StereoAtomEnd = -1; }
            _bonds[m++] = b;
        }
        if (_chiralRefs is not null)
        {
            var refs = new Dictionary<int, int[]>();
            foreach (var (atom, nb) in _chiralRefs)
            {
                if (map[atom] < 0) continue;
                if (nb.Any(x => x >= 0 && map[x] < 0)) { _atoms[map[atom]].Chirality = ChiralTag.None; continue; }
                refs[map[atom]] = nb.Select(x => x < 0 ? x : map[x]).ToArray();
            }
            _chiralRefs = refs;
        }
        _atomCount = n;
        _bondCount = m;
        if (_coords is not null) Array.Resize(ref _coords, n);
        InvalidateTopology();
        return map;
    }

    /// <summary>Removes one bond.</summary>
    public void RemoveBond(int bondIndex)
    {
        if ((uint)bondIndex >= (uint)_bondCount) throw new ArgumentOutOfRangeException(nameof(bondIndex));
        Array.Copy(_bonds, bondIndex + 1, _bonds, bondIndex, _bondCount - bondIndex - 1);
        _bondCount--;
        InvalidateTopology();
    }

    /// <summary>Discards cached adjacency, ring and Kekulé data. Call after mutating bonds through <see cref="BondAt"/>.</summary>
    public void InvalidateTopology()
    {
        _adjStart = _adjAtom = _adjBond = null;
        _rings = null;
        _kekule = null;
    }

    /// <summary>A deep copy.</summary>
    public Molecule Clone()
    {
        var m = new Molecule(Name, _atomCount)
        {
            _atoms = (Atom[])_atoms.Clone(),
            _bonds = (Bond[])_bonds.Clone(),
            _atomCount = _atomCount,
            _bondCount = _bondCount,
            _coords = (Point3[]?)_coords?.Clone(),
            _kekule = (byte[]?)_kekule?.Clone(),
        };
        if (_chiralRefs is not null) m._chiralRefs = _chiralRefs.ToDictionary(kv => kv.Key, kv => (int[])kv.Value.Clone());
        if (_properties is not null) foreach (var kv in _properties) m.Properties[kv.Key] = kv.Value;
        return m;
    }

    /// <summary>
    /// The neighbour order that <see cref="Atom.Chirality"/> refers to (−1 stands for the atom's own implicit
    /// hydrogen), or null when none was recorded.
    /// </summary>
    public int[]? GetChiralReference(int atom) => _chiralRefs is not null && _chiralRefs.TryGetValue(atom, out var r) ? r : null;

    /// <summary>Records the neighbour order that the atom's <see cref="Atom.Chirality"/> refers to.</summary>
    public void SetChiralReference(int atom, int[]? neighbours)
    {
        if (neighbours is null) { _chiralRefs?.Remove(atom); return; }
        (_chiralRefs ??= [])[atom] = neighbours;
    }

    // ---- Topology --------------------------------------------------------------------------------

    /// <summary>Indices of the atoms bonded to <paramref name="atom"/>.</summary>
    public ReadOnlySpan<int> Neighbors(int atom)
    {
        EnsureAdjacency();
        return _adjAtom.AsSpan(_adjStart![atom], _adjStart[atom + 1] - _adjStart[atom]);
    }

    /// <summary>Indices of the bonds of <paramref name="atom"/>, parallel to <see cref="Neighbors"/>.</summary>
    public ReadOnlySpan<int> NeighborBonds(int atom)
    {
        EnsureAdjacency();
        return _adjBond.AsSpan(_adjStart![atom], _adjStart[atom + 1] - _adjStart[atom]);
    }

    /// <summary>Number of explicit neighbours (SMARTS <c>D</c>).</summary>
    public int Degree(int atom)
    {
        EnsureAdjacency();
        return _adjStart![atom + 1] - _adjStart[atom];
    }

    /// <summary>Number of non-hydrogen neighbours.</summary>
    public int HeavyDegree(int atom)
    {
        int n = 0;
        foreach (int j in Neighbors(atom)) if (_atoms[j].AtomicNumber != 1) n++;
        return n;
    }

    /// <summary>Total hydrogens on an atom: implicit + explicit counts + bonded H atoms.</summary>
    public int TotalHydrogens(int atom)
    {
        int n = _atoms[atom].HydrogenCount;
        foreach (int j in Neighbors(atom)) if (_atoms[j].AtomicNumber == 1) n++;
        return n;
    }

    /// <summary>Index of the bond between two atoms, or −1.</summary>
    public int BondBetween(int a, int b)
    {
        if (_adjStart is not null)
        {
            var nb = Neighbors(a);
            for (int i = 0; i < nb.Length; i++) if (nb[i] == b) return _adjBond![_adjStart[a] + i];
            return -1;
        }
        for (int i = 0; i < _bondCount; i++)
        {
            ref var bond = ref _bonds[i];
            if ((bond.Begin == a && bond.End == b) || (bond.Begin == b && bond.End == a)) return i;
        }
        return -1;
    }

    /// <summary>Ring perception results (SSSR), computed on first use.</summary>
    public RingInfo Rings => _rings ??= RingInfo.Compute(this);

    /// <summary>
    /// The Kekulé order (1, 2 or 3) of bond <paramref name="bondIndex"/>; aromatic bonds get their single/double
    /// assignment from <see cref="Aromaticity.Kekulize"/>.
    /// </summary>
    public int KekuleOrder(int bondIndex)
    {
        var b = _bonds[bondIndex];
        if (b.Order != BondOrder.Aromatic) return (int)b.Order;
        _kekule ??= Aromaticity.ComputeKekuleOrders(this, throwOnFailure: false);
        return _kekule[bondIndex];
    }

    internal void SetKekuleCache(byte[]? orders) => _kekule = orders;

    /// <summary>Sum of Kekulé bond orders to explicit neighbours (H atoms included).</summary>
    public int ExplicitValence(int atom)
    {
        int v = 0;
        foreach (int b in NeighborBonds(atom)) v += KekuleOrder(b);
        return v;
    }

    /// <summary>Total valence: explicit valence plus the atom's own hydrogen counts (SMARTS <c>v</c>).</summary>
    public int TotalValence(int atom) => ExplicitValence(atom) + _atoms[atom].HydrogenCount;

    /// <summary>Connected components as arrays of atom indices.</summary>
    public List<int[]> Fragments()
    {
        var seen = new bool[_atomCount];
        var result = new List<int[]>();
        var stack = new Stack<int>();
        for (int s = 0; s < _atomCount; s++)
        {
            if (seen[s]) continue;
            var comp = new List<int>();
            stack.Push(s);
            seen[s] = true;
            while (stack.Count > 0)
            {
                int u = stack.Pop();
                comp.Add(u);
                foreach (int v in Neighbors(u)) if (!seen[v]) { seen[v] = true; stack.Push(v); }
            }
            comp.Sort();
            result.Add(comp.ToArray());
        }
        return result;
    }

    /// <summary>Shortest-path (bond count) distance matrix by BFS; unreachable pairs are −1.</summary>
    public int[,] TopologicalDistances()
    {
        int n = _atomCount;
        var d = new int[n, n];
        var queue = new int[n];
        for (int s = 0; s < n; s++)
        {
            for (int j = 0; j < n; j++) d[s, j] = -1;
            d[s, s] = 0;
            int head = 0, tail = 0;
            queue[tail++] = s;
            while (head < tail)
            {
                int u = queue[head++];
                foreach (int v in Neighbors(u))
                {
                    if (d[s, v] >= 0) continue;
                    d[s, v] = d[s, u] + 1;
                    queue[tail++] = v;
                }
            }
        }
        return d;
    }

    // ---- Hydrogens -------------------------------------------------------------------------------

    /// <summary>
    /// The valences an atom may adopt given its charge. p-block atoms use the valences of the isoelectronic
    /// neutral element (N⁺ behaves like C, O⁻ like F, C⁻ like N); metals get none.
    /// </summary>
    public static ReadOnlySpan<int> AllowedValences(int atomicNumber, int charge)
    {
        if (atomicNumber <= 0) return [];
        var e = PeriodicTable.Get(atomicNumber);
        if (charge == 0) return e.IsMetal ? [] : e.Valences;
        if (atomicNumber is 1) return charge == -1 || charge == 1 ? [0] : [];
        bool pBlock = atomicNumber is (>= 5 and <= 9) or (>= 13 and <= 17) or (>= 31 and <= 35) or (>= 49 and <= 53);
        if (!pBlock) return [];
        int eff = atomicNumber - charge;
        int row(int z) => z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : z <= 36 ? 4 : 5;
        if (eff < 1 || eff > 92 || row(eff) != row(atomicNumber)) return [];
        var ee = PeriodicTable.Get(eff);
        if (ee.IsMetal && eff != 5) return [];
        return ee.Valences;
    }

    /// <summary>
    /// Recomputes <see cref="Atom.ImplicitHydrogens"/> for every atom that allows implicit hydrogens, from the
    /// lowest allowed valence not below the current bonding (aromatic bonds use the Kekulé assignment).
    /// </summary>
    public void UpdateImplicitHydrogens()
    {
        _kekule = null;
        bool anyAromatic = false;
        for (int i = 0; i < _bondCount; i++) if (_bonds[i].Order == BondOrder.Aromatic) { anyAromatic = true; break; }
        if (anyAromatic) _kekule = Aromaticity.ComputeKekuleOrders(this, throwOnFailure: false);
        for (int i = 0; i < _atomCount; i++)
        {
            ref var a = ref _atoms[i];
            if (a.NoImplicit) { a.ImplicitHydrogens = 0; continue; }
            int v = a.ExplicitHydrogens;
            foreach (int b in NeighborBonds(i)) v += KekuleOrder(b);
            var allowed = AllowedValences(a.AtomicNumber, a.FormalCharge);
            int h = 0;
            foreach (int val in allowed)
                if (val >= v) { h = val - v; break; }
            a.ImplicitHydrogens = h;
        }
    }

    /// <summary>
    /// A copy in which every hydrogen is an explicit atom bonded to its parent. When the molecule has coordinates, the
    /// new hydrogens are placed 1.0 Å (3-D) or 0.7 × mean bond length (2-D) from their parent, spread into the free directions.
    /// </summary>
    public Molecule AddHydrogens()
    {
        var m = Clone();
        var oldCoords = _coords;
        m._coords = null;
        int n = m._atomCount;
        var parents = new List<int>();
        for (int i = 0; i < n; i++)
        {
            ref var a = ref m._atoms[i];
            int h = a.HydrogenCount;
            if (h == 0) continue;
            a.ImplicitHydrogens = 0;
            a.ExplicitHydrogens = 0;
            a.NoImplicit = true;
            var chiral = m.GetChiralReference(i);
            int firstH = -1;
            for (int k = 0; k < h; k++)
            {
                int hi = m.AddAtom(new Atom(1) { NoImplicit = true });
                m.AddBond(i, hi);
                parents.Add(i);
                if (firstH < 0) firstH = hi;
            }
            if (chiral is not null && firstH >= 0)
                m.SetChiralReference(i, chiral.Select(x => x < 0 ? firstH : x).ToArray());
        }
        for (int i = 0; i < m._atomCount; i++) m._atoms[i].NoImplicit = true;
        if (oldCoords is not null) m._coords = PlaceHydrogens(m, oldCoords, n);
        return m;
    }

    private static Point3[] PlaceHydrogens(Molecule m, Point3[] heavy, int firstH)
    {
        var c = new Point3[m._atomCount];
        Array.Copy(heavy, c, firstH);
        bool is3D = heavy.Any(p => Math.Abs(p.Z) > 1e-4);
        double meanBond = 1.5;
        if (!is3D)
        {
            var lens = new List<double>();
            for (int b = 0; b < m._bondCount; b++)
                if (m._bonds[b].Begin < firstH && m._bonds[b].End < firstH) lens.Add(Point3.Distance(heavy[m._bonds[b].Begin], heavy[m._bonds[b].End]));
            if (lens.Count > 0) meanBond = lens.Average();
        }
        double hLen = is3D ? 1.0 : 0.7 * meanBond;
        for (int p = 0; p < firstH; p++)
        {
            var hs = new List<int>();
            var placed = new List<Point3>();
            foreach (int v in m.Neighbors(p)) { if (v >= firstH) hs.Add(v); else placed.Add(c[v] - c[p]); }
            if (hs.Count == 0) continue;
            if (!is3D)
            {
                var angs = placed.Select(d => Math.Atan2(d.Y, d.X)).Select(a => a < 0 ? a + 2 * Math.PI : a).OrderBy(a => a).ToList();
                double start = 0, gap = 2 * Math.PI;
                if (angs.Count == 1) { start = angs[0]; gap = 2 * Math.PI; }
                else if (angs.Count > 1)
                {
                    gap = -1;
                    for (int i = 0; i < angs.Count; i++)
                    {
                        double a0 = angs[i], a1 = i + 1 < angs.Count ? angs[i + 1] : angs[0] + 2 * Math.PI;
                        if (a1 - a0 > gap) { gap = a1 - a0; start = a0; }
                    }
                }
                for (int k = 0; k < hs.Count; k++)
                {
                    double a = angs.Count == 0 ? 2 * Math.PI * k / hs.Count : start + gap * (k + 1) / (hs.Count + 1);
                    c[hs[k]] = c[p] + new Point3(Math.Cos(a), Math.Sin(a)) * hLen;
                }
            }
            else
            {
                var sum = new Point3(0, 0, 0);
                foreach (var d in placed) sum += d * (1 / Math.Max(1e-9, d.Length));
                var axis = sum.Length < 1e-6 ? new Point3(1, 0, 0) : sum * (-1 / sum.Length);
                var perp = Point3.Cross(axis, Math.Abs(axis.X) < 0.9 ? new Point3(1, 0, 0) : new Point3(0, 1, 0));
                perp = perp * (1 / perp.Length);
                var perp2 = Point3.Cross(axis, perp);
                // tetrahedral-ish cone around the free axis (angle from the axis: 70.5° for XH3, 54.7° for X2H2)
                if (placed.Count == 2)
                {
                    var normal = Point3.Cross(placed[0], placed[1]);
                    if (normal.Length > 1e-6) { perp = normal * (1 / normal.Length); perp2 = Point3.Cross(axis, perp); }
                }
                double theta = hs.Count == 1 ? 0 : placed.Count == 0 ? 1.9106 : placed.Count == 1 ? 1.2310 : 0.9553;
                for (int k = 0; k < hs.Count; k++)
                {
                    if (placed.Count == 0 && k == 0) { c[hs[k]] = c[p] + axis * hLen; continue; }
                    int idx = placed.Count == 0 ? k - 1 : k;
                    int cnt = placed.Count == 0 ? hs.Count - 1 : hs.Count;
                    double phi = 2 * Math.PI * idx / Math.Max(1, cnt);
                    var dir = axis * Math.Cos(theta) + (perp * Math.Cos(phi) + perp2 * Math.Sin(phi)) * Math.Sin(theta);
                    c[hs[k]] = c[p] + dir * (hLen / dir.Length);
                }
            }
        }
        return c;
    }

    /// <summary>A copy with hydrogen atoms folded back into implicit counts (isotopic and bridging H are kept).</summary>
    public Molecule RemoveHydrogens()
    {
        var m = Clone();
        var remove = new List<int>();
        for (int i = 0; i < m._atomCount; i++)
        {
            var a = m._atoms[i];
            if (a.AtomicNumber != 1 || a.Isotope != 0 || a.FormalCharge != 0 || m.Degree(i) != 1) continue;
            int parent = m.Neighbors(i)[0];
            if (m._atoms[parent].AtomicNumber == 1) continue;
            m._atoms[parent].ExplicitHydrogens++;
            remove.Add(i);
        }
        m.RemoveAtoms(remove);
        return m;
    }

    // ---- Misc ------------------------------------------------------------------------------------

    /// <summary>Canonical SMILES (shortcut for <see cref="SmilesWriter.Write"/>).</summary>
    public string ToSmiles(bool canonical = true, bool includeStereo = true) => SmilesWriter.Write(this, canonical, includeStereo);

    /// <summary>Molecular formula in Hill order, e.g. <c>C9H8O4</c>.</summary>
    public string Formula => MolecularDescriptors.Formula(this);

    public override string ToString() => string.IsNullOrEmpty(Name) ? ToSmiles() : $"{Name}: {ToSmiles()}";

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void EnsureAdjacency()
    {
        if (_adjStart is null) BuildAdjacency();
    }

    private void BuildAdjacency()
    {
        var start = new int[_atomCount + 1];
        for (int i = 0; i < _bondCount; i++) { start[_bonds[i].Begin + 1]++; start[_bonds[i].End + 1]++; }
        for (int i = 0; i < _atomCount; i++) start[i + 1] += start[i];
        var fill = new int[_atomCount];
        var atom = new int[2 * _bondCount];
        var bond = new int[2 * _bondCount];
        for (int i = 0; i < _bondCount; i++)
        {
            int a = _bonds[i].Begin, b = _bonds[i].End;
            int pa = start[a] + fill[a]++, pb = start[b] + fill[b]++;
            atom[pa] = b; bond[pa] = i;
            atom[pb] = a; bond[pb] = i;
        }
        _adjAtom = atom;
        _adjBond = bond;
        _adjStart = start;
    }
}
