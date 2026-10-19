using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.DrugDiscovery;

/// <summary>A unit quaternion orientation (w, x, y, z).</summary>
internal readonly record struct Quat(double W, double X, double Y, double Z)
{
    public static readonly Quat Identity = new(1, 0, 0, 0);

    public static Quat FromRotationVector(Vec3 v)
    {
        double angle = v.Length;
        if (angle < 1e-12) return Identity;
        double s = Math.Sin(angle / 2) / angle;
        return new Quat(Math.Cos(angle / 2), v.X * s, v.Y * s, v.Z * s);
    }

    public static Quat operator *(Quat a, Quat b) => new(
        a.W * b.W - a.X * b.X - a.Y * b.Y - a.Z * b.Z,
        a.W * b.X + a.X * b.W + a.Y * b.Z - a.Z * b.Y,
        a.W * b.Y - a.X * b.Z + a.Y * b.W + a.Z * b.X,
        a.W * b.Z + a.X * b.Y - a.Y * b.X + a.Z * b.W);

    public Quat Normalized()
    {
        double n = Math.Sqrt(W * W + X * X + Y * Y + Z * Z);
        return new Quat(W / n, X / n, Y / n, Z / n);
    }

    public Matrix3 ToMatrix()
    {
        var q = Normalized();
        double w = q.W, x = q.X, y = q.Y, z = q.Z;
        return new Matrix3(
            1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
            2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
            2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y));
    }

    public static Quat Random(Random rng)
    {
        // uniform random rotation (Shoemake)
        double u1 = rng.NextDouble(), u2 = rng.NextDouble() * 2 * Math.PI, u3 = rng.NextDouble() * 2 * Math.PI;
        double a = Math.Sqrt(1 - u1), b = Math.Sqrt(u1);
        return new Quat(b * Math.Cos(u3), a * Math.Sin(u2), a * Math.Cos(u2), b * Math.Sin(u3));
    }
}

/// <summary>A ligand conformation: root position, orientation and torsion angles (relative to the input geometry).</summary>
internal sealed class LigandConf(Vec3 position, Quat orientation, double[] torsions)
{
    public Vec3 Position = position;
    public Quat Orientation = orientation;
    public readonly double[] Torsions = torsions;

    public LigandConf Clone() => new(Position, Orientation, (double[])Torsions.Clone());

    public void Increment(ReadOnlySpan<double> p, double factor)
    {
        Position += new Vec3(p[0], p[1], p[2]) * factor;
        Orientation = (Quat.FromRotationVector(new Vec3(p[3], p[4], p[5]) * factor) * Orientation).Normalized();
        for (int t = 0; t < Torsions.Length; t++) Torsions[t] = NormalizeAngle(Torsions[t] + factor * p[6 + t]);
    }

    public static double NormalizeAngle(double a)
    {
        if (!double.IsFinite(a)) return 0;
        return Math.IEEERemainder(a, 2 * Math.PI);
    }
}

/// <summary>
/// The ligand as a torsion tree (Vina/AutoDock style): a rigid root fragment plus branches across rotatable bonds.
/// Rotatable bonds follow Meeko's rules (acyclic single bonds that are not amide C–N or next to a triple bond);
/// rotations that only move hydrogens are not degrees of freedom, but still count towards Vina's torsional term.
/// </summary>
internal sealed class FlexibleLigand
{
    public sealed record Torsion(int From, int To, int[] Moving, int Depth);

    public Molecule Molecule { get; }
    public int AtomCount { get; }
    public XsAtomType?[] Types { get; }
    public int[] Scored { get; }               // atoms with an XS type
    public Torsion[] Torsions { get; }         // ordered deepest first
    public (int I, int J)[] IntraPairs { get; }
    public double NumRotors { get; }           // Vina's num_tors for the conf-independent term
    public int RootAtom { get; }
    public int[] RotatableBonds { get; }       // bond indices of Meeko-rotatable bonds
    public int[] DofBonds { get; }             // bond index per torsion (same order as Torsions)
    private readonly Vec3[] _local;            // reference coordinates relative to the root origin
    public Vec3 ReferenceOrigin { get; }
    public double GyrationRadius { get; }

    private static readonly Lazy<SmartsPattern> Amide = new(() => SmartsPattern.Parse("[NX3]-[CX3]=[O,N,S]"));

    public FlexibleLigand(Molecule mol)
    {
        var coords = mol.Coordinates ?? throw new InvalidOperationException("The ligand needs 3-D coordinates.");
        Molecule = mol;
        int n = AtomCount = mol.AtomCount;
        Types = VinaAtomTypes.AssignLigand(mol);
        Scored = Enumerable.Range(0, n).Where(i => Types[i] is not null).ToArray();
        var bonds = mol.Bonds;
        var rings = mol.Rings;
        var amide = new HashSet<int>();
        foreach (var m in Amide.Value.Matches(mol, unique: false)) amide.Add(mol.BondBetween(m[0], m[1]));
        var nearTriple = new HashSet<int>();
        for (int b = 0; b < bonds.Length; b++)
        {
            if (bonds[b].Order != BondOrder.Triple) continue;
            foreach (int end in new[] { bonds[b].Begin, bonds[b].End })
                foreach (int nb in mol.NeighborBonds(end)) if (nb != b) nearTriple.Add(nb);
        }
        bool IsHeavy(int a) => mol.Atoms[a].AtomicNumber != 1;
        int HeavyDeg(int a) { int c = 0; foreach (int x in mol.Neighbors(a)) if (IsHeavy(x)) c++; return c; }
        bool HasPolarH(int a) => mol.Atoms[a].AtomicNumber is 7 or 8 or 16 && mol.TotalHydrogens(a) > 0;
        var rotatable = new List<int>();
        var dof = new List<int>();
        double numTors = 0;
        for (int b = 0; b < bonds.Length; b++)
        {
            var bd = bonds[b];
            if (bd.Order != BondOrder.Single || rings.IsBondInRing(b) || amide.Contains(b) || nearTriple.Contains(b)) continue;
            int x = bd.Begin, y = bd.End;
            if (!IsHeavy(x) || !IsHeavy(y)) continue;
            int hx = HeavyDeg(x), hy = HeavyDeg(y);
            bool xMoves = hx > 1 || HasPolarH(x), yMoves = hy > 1 || HasPolarH(y);
            if (!xMoves || !yMoves) continue; // rotates nothing scored or explicit
            rotatable.Add(b);
            numTors += (hy > 1 ? 0.5 : 0) + (hx > 1 ? 0.5 : 0);
            if (hx > 1 && hy > 1) dof.Add(b);
        }
        RotatableBonds = [.. rotatable];
        NumRotors = numTors;

        // rigid fragments = components after cutting DOF bonds
        var cut = new HashSet<int>(dof);
        var frag = Enumerable.Repeat(-1, n).ToArray();
        int nf = 0;
        for (int s = 0; s < n; s++)
        {
            if (frag[s] >= 0) continue;
            var stack = new Stack<int>();
            stack.Push(s); frag[s] = nf;
            while (stack.Count > 0)
            {
                int u = stack.Pop();
                var nb = mol.Neighbors(u);
                var nbb = mol.NeighborBonds(u);
                for (int k = 0; k < nb.Length; k++)
                    if (!cut.Contains(nbb[k]) && frag[nb[k]] < 0) { frag[nb[k]] = nf; stack.Push(nb[k]); }
            }
            nf++;
        }
        // fragment graph; choose the root minimising the maximum branch depth (ties: more heavy atoms)
        var adj = Enumerable.Range(0, nf).Select(_ => new List<(int F, int Bond)>()).ToArray();
        foreach (int b in dof) { int fa = frag[bonds[b].Begin], fb = frag[bonds[b].End]; adj[fa].Add((fb, b)); adj[fb].Add((fa, b)); }
        int[] Depths(int root)
        {
            var d = Enumerable.Repeat(-1, nf).ToArray();
            var q = new Queue<int>(); q.Enqueue(root); d[root] = 0;
            while (q.Count > 0) { int u = q.Dequeue(); foreach (var (v, _) in adj[u]) if (d[v] < 0) { d[v] = d[u] + 1; q.Enqueue(v); } }
            return d;
        }
        int bestRoot = 0, bestDepth = int.MaxValue, bestSize = -1;
        for (int f = 0; f < nf; f++)
        {
            int depth = Depths(f).Max();
            int size = Enumerable.Range(0, n).Count(i => frag[i] == f && IsHeavy(i));
            if (depth < bestDepth || (depth == bestDepth && size > bestSize)) { bestRoot = f; bestDepth = depth; bestSize = size; }
        }
        var fd = Depths(bestRoot);
        var torsions = new List<(Torsion T, int Bond)>();
        foreach (int b in dof)
        {
            int a1 = bonds[b].Begin, a2 = bonds[b].End;
            int parent = fd[frag[a1]] < fd[frag[a2]] ? a1 : a2;
            int child = parent == a1 ? a2 : a1;
            // moving atoms: everything reachable from child without crossing this bond
            var seen = new bool[n];
            seen[parent] = true; seen[child] = true;
            var moving = new List<int> { child };
            var st = new Stack<int>(); st.Push(child);
            while (st.Count > 0)
            {
                int u = st.Pop();
                foreach (int v in mol.Neighbors(u)) if (!seen[v]) { seen[v] = true; moving.Add(v); st.Push(v); }
            }
            torsions.Add((new Torsion(parent, child, [.. moving], fd[frag[child]]), b));
        }
        torsions.Sort((p, q) => q.T.Depth.CompareTo(p.T.Depth));
        Torsions = torsions.Select(t => t.T).ToArray();
        DofBonds = torsions.Select(t => t.Bond).ToArray();

        // intramolecular pairs: different rigid fragments, more than three bonds apart
        var topo = mol.TopologicalDistances();
        var pairs = new List<(int, int)>();
        for (int x = 0; x < Scored.Length; x++)
            for (int y = x + 1; y < Scored.Length; y++)
            {
                int i = Scored[x], j = Scored[y];
                if (frag[i] != frag[j] && topo[i, j] > 3) pairs.Add((i, j));
            }
        IntraPairs = [.. pairs];

        var rootAtoms = Enumerable.Range(0, n).Where(i => frag[i] == bestRoot && IsHeavy(i)).ToArray();
        if (rootAtoms.Length == 0) rootAtoms = [0];
        RootAtom = rootAtoms[0];
        var origin = rootAtoms.Aggregate(Vec3.Zero, (s, i) => s + ToVec(coords[i])) / rootAtoms.Length;
        ReferenceOrigin = origin;
        _local = coords.Select(p => ToVec(p) - origin).ToArray();
        var heavyAll = Enumerable.Range(0, n).Where(IsHeavy).ToArray();
        GyrationRadius = Math.Sqrt(heavyAll.Average(i => _local[i].LengthSquared));
    }

    internal static Vec3 ToVec(Point3 p) => new(p.X, p.Y, p.Z);

    /// <summary>Reference coordinates relative to <see cref="ReferenceOrigin"/> (used by the GPU backend).</summary>
    internal ReadOnlySpan<Vec3> LocalCoordinates => _local;

    /// <summary>The conformation that reproduces the input coordinates.</summary>
    public LigandConf InputConf() => new(ReferenceOrigin, Quat.Identity, new double[Torsions.Length]);

    /// <summary>Cartesian coordinates of a conformation.</summary>
    public void Build(LigandConf c, Vec3[] x)
    {
        Array.Copy(_local, x, AtomCount);
        for (int t = 0; t < Torsions.Length; t++)
        {
            double angle = c.Torsions[t];
            if (angle == 0) continue;
            var tor = Torsions[t];
            var a = x[tor.From];
            var axis = (x[tor.To] - a).Normalized();
            var rot = Matrix3.FromAxisAngle(axis, angle * 180 / Math.PI);
            foreach (int i in tor.Moving) x[i] = a + rot * (x[i] - a);
        }
        var r = c.Orientation.ToMatrix();
        for (int i = 0; i < AtomCount; i++) x[i] = c.Position + r * x[i];
    }

    /// <summary>Maps atomic gradients (dE/dx) to the conformation-space gradient (position, rotation vector, torsions).</summary>
    public void ConfGradient(LigandConf c, Vec3[] x, Vec3[] atomGrad, Span<double> g)
    {
        var sum = Vec3.Zero;
        var torque = Vec3.Zero;
        for (int i = 0; i < AtomCount; i++)
        {
            sum += atomGrad[i];
            torque += Vec3.Cross(x[i] - c.Position, atomGrad[i]);
        }
        g[0] = sum.X; g[1] = sum.Y; g[2] = sum.Z;
        g[3] = torque.X; g[4] = torque.Y; g[5] = torque.Z;
        for (int t = 0; t < Torsions.Length; t++)
        {
            var tor = Torsions[t];
            var a = x[tor.From];
            var u = (x[tor.To] - a).Normalized();
            var tq = Vec3.Zero;
            foreach (int i in tor.Moving) tq += Vec3.Cross(x[i] - a, atomGrad[i]);
            g[6 + t] = Vec3.Dot(tq, u);
        }
    }
}
