namespace Biotech.Net.StructuralBiology;

/// <summary>Which atoms represent a residue in a contact map.</summary>
public enum ContactMode
{
    /// <summary>Cα atoms (default cut-off 8 Å).</summary>
    CAlpha,

    /// <summary>Cβ atoms (Cα for glycine), the CASP convention.</summary>
    CBeta,

    /// <summary>Minimum distance over all heavy atoms (default cut-off 4.5 Å).</summary>
    AnyHeavyAtom,
}

/// <summary>
/// A symmetric residue–residue contact map stored as a bit matrix, built with a cell-list neighbour search.
/// </summary>
public sealed class ContactMap
{
    private readonly ulong[] _bits;
    private readonly int _words;

    private ContactMap(IReadOnlyList<Residue> residues, double cutoff, ContactMode mode)
    {
        Residues = residues;
        Cutoff = cutoff;
        Mode = mode;
        _words = (residues.Count + 63) / 64;
        _bits = new ulong[residues.Count * _words];
    }

    /// <summary>Residues in map order.</summary>
    public IReadOnlyList<Residue> Residues { get; }

    /// <summary>Distance cut-off in Å.</summary>
    public double Cutoff { get; }

    /// <summary>Representation used.</summary>
    public ContactMode Mode { get; }

    /// <summary>Number of residues (rows/columns).</summary>
    public int Size => Residues.Count;

    /// <summary>True when residues <paramref name="i"/> and <paramref name="j"/> are in contact (the diagonal is false).</summary>
    public bool this[int i, int j] => (_bits[i * _words + (j >> 6)] & (1UL << (j & 63))) != 0;

    /// <summary>Number of contacting pairs (i &lt; j).</summary>
    public int ContactCount
    {
        get
        {
            long total = 0;
            foreach (ulong w in _bits) total += System.Numerics.BitOperations.PopCount(w);
            return (int)(total / 2);
        }
    }

    /// <summary>Contacting pairs (i &lt; j) with sequence separation at least <paramref name="minSeparation"/>.</summary>
    public IEnumerable<(int I, int J)> Pairs(int minSeparation = 1)
    {
        for (int i = 0; i < Size; i++)
            for (int j = i + Math.Max(1, minSeparation); j < Size; j++)
                if (this[i, j]) yield return (i, j);
    }

    /// <summary>Contact order: mean sequence separation of contacts divided by chain length (Plaxco 1998).</summary>
    public double RelativeContactOrder(int minSeparation = 1)
    {
        long sum = 0; int n = 0;
        foreach (var (i, j) in Pairs(minSeparation)) { sum += j - i; n++; }
        return n == 0 ? 0 : (double)sum / n / Size;
    }

    /// <summary>The map as a dense boolean matrix.</summary>
    public bool[,] ToMatrix()
    {
        var m = new bool[Size, Size];
        for (int i = 0; i < Size; i++) for (int j = 0; j < Size; j++) m[i, j] = this[i, j];
        return m;
    }

    private void Set(int i, int j)
    {
        if (i == j) return;
        _bits[i * _words + (j >> 6)] |= 1UL << (j & 63);
        _bits[j * _words + (i >> 6)] |= 1UL << (i & 63);
    }

    /// <summary>Contact map of the amino acids of a chain.</summary>
    public static ContactMap Compute(Chain chain, double? cutoff = null, ContactMode mode = ContactMode.CAlpha)
        => Compute(chain.AminoAcids(), cutoff, mode);

    /// <summary>Contact map of an arbitrary residue list.</summary>
    /// <param name="cutoff">Å; defaults to 8 for Cα/Cβ and 4.5 for heavy atoms.</param>
    public static ContactMap Compute(IReadOnlyList<Residue> residues, double? cutoff = null, ContactMode mode = ContactMode.CAlpha)
    {
        double cut = cutoff ?? (mode == ContactMode.AnyHeavyAtom ? 4.5 : 8.0);
        var map = new ContactMap(residues, cut, mode);
        var points = new List<Vec3>();
        var owner = new List<int>();
        for (int r = 0; r < residues.Count; r++)
        {
            var res = residues[r];
            switch (mode)
            {
                case ContactMode.CAlpha:
                    if (res.GetAtom("CA") is { } ca) { points.Add(ca.Position); owner.Add(r); }
                    break;
                case ContactMode.CBeta:
                    if ((res.GetAtom("CB") ?? res.GetAtom("CA")) is { } cb) { points.Add(cb.Position); owner.Add(r); }
                    break;
                default:
                    foreach (var a in res.Atoms)
                        if (a.Element != "H") { points.Add(a.Position); owner.Add(r); }
                    break;
            }
        }
        var ns = new NeighborSearch(System.Runtime.InteropServices.CollectionsMarshal.AsSpan(points), cut);
        foreach (var (i, j) in ns.AllPairs(cut)) map.Set(owner[i], owner[j]);
        return map;
    }

    /// <summary>Dense Cα distance matrix of a chain.</summary>
    public static double[,] DistanceMatrix(Chain chain)
    {
        var ca = Superposition.CaTrace(chain);
        var d = new double[ca.Length, ca.Length];
        for (int i = 0; i < ca.Length; i++)
            for (int j = i + 1; j < ca.Length; j++)
                d[i, j] = d[j, i] = Vec3.Distance(ca[i], ca[j]);
        return d;
    }
}

/// <summary>A disulfide bridge detected from geometry.</summary>
public sealed record Disulfide(Residue CysteineA, Residue CysteineB, double Distance)
{
    /// <summary>The Cβ–Sγ–Sγ–Cβ torsion (χ3) in degrees; ±90° is ideal.</summary>
    public double Chi3 => CysteineA.GetAtom("CB") is { } a && CysteineB.GetAtom("CB") is { } b
        ? Geometry.Dihedral(a, CysteineA["SG"], CysteineB["SG"], b) : double.NaN;
}

/// <summary>Detects disulfide bridges from Sγ–Sγ distances.</summary>
public static class Disulfides
{
    /// <summary>
    /// Cysteine pairs whose Sγ atoms are between <paramref name="minDistance"/> and <paramref name="maxDistance"/> Å
    /// apart (an ideal S–S bond is 2.04 Å), ordered by chain and residue.
    /// </summary>
    public static IReadOnlyList<Disulfide> Find(Model model, double minDistance = 1.6, double maxDistance = 2.5)
    {
        var sg = model.Residues.Where(r => r.Name is "CYS" or "CYX" && r.HasAtom("SG")).Select(r => r["SG"]).ToList();
        var ns = NeighborSearch.ForAtoms(sg, maxDistance);
        var result = new List<Disulfide>();
        foreach (var (i, j) in ns.AllPairs(maxDistance).OrderBy(p => p.I).ThenBy(p => p.J))
        {
            double d = Geometry.Distance(sg[i], sg[j]);
            if (d >= minDistance) result.Add(new Disulfide(sg[i].Residue!, sg[j].Residue!, d));
        }
        return result;
    }

    /// <summary>Disulfides of the first model.</summary>
    public static IReadOnlyList<Disulfide> Find(Structure structure) => Find(structure.FirstModel);
}
