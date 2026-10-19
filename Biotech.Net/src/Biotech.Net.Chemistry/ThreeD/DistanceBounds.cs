namespace Biotech.Net.Chemistry;

/// <summary>
/// Square matrix of interatomic distance bounds: upper bounds in the upper triangle (i &lt; j), lower bounds in the
/// lower triangle, as in RDKit's <c>DistGeom::BoundsMatrix</c>.
/// </summary>
internal sealed class BoundsMatrix
{
    public readonly int N;
    private readonly double[] _d;

    public BoundsMatrix(int n, double lower = 0.0, double upper = 1000.0)
    {
        N = n;
        _d = new double[n * n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                _d[i * n + j] = i < j ? upper : i > j ? lower : 0.0;
    }

    public double Upper(int i, int j) => i < j ? _d[i * N + j] : _d[j * N + i];
    public double Lower(int i, int j) => i < j ? _d[j * N + i] : _d[i * N + j];
    public void SetUpper(int i, int j, double v) { if (i < j) _d[i * N + j] = v; else _d[j * N + i] = v; }
    public void SetLower(int i, int j, double v) { if (i < j) _d[j * N + i] = v; else _d[i * N + j] = v; }

    /// <summary>Triangle-inequality smoothing (Floyd–Warshall style); false when bounds become inconsistent.</summary>
    public bool TriangleSmooth(double tol = 0.0)
    {
        int n = N;
        var d = _d;
        for (int k = 0; k < n; k++)
        {
            for (int i = 0; i < n - 1; i++)
            {
                if (i == k) continue;
                int ii = Math.Min(i, k), ik = Math.Max(i, k);
                double uik = d[ii * n + ik], lik = d[ik * n + ii];
                for (int j = i + 1; j < n; j++)
                {
                    if (j == k) continue;
                    int jj = Math.Min(j, k), jk = Math.Max(j, k);
                    double ukj = d[jj * n + jk];
                    double sum = uik + ukj;
                    if (d[i * n + j] > sum) d[i * n + j] = sum;
                    double diff1 = lik - ukj;
                    double diff2 = d[jk * n + jj] - uik;
                    if (d[j * n + i] < diff1) d[j * n + i] = diff1;
                    else if (d[j * n + i] < diff2) d[j * n + i] = diff2;
                    double lb = d[j * n + i], ub = d[i * n + j];
                    if (tol > 0 && (lb - ub) / lb > 0 && (lb - ub) / lb < tol) d[i * n + j] = lb;
                    else if (lb - ub > 0) return false;
                }
            }
        }
        return true;
    }
}

/// <summary>Torsion class of a 1-4 path used when setting its distance bounds.</summary>
internal enum TorsionClass { Cis, Trans, CisTrans, Flexible, Custom }

/// <summary>A 1-4 path with its torsion class (amide/ester paths may be "forced" cis/trans for ETKDG).</summary>
internal readonly record struct Path14(int A1, int A2, int A3, int A4, TorsionClass Class, bool Forced);

/// <summary>
/// Builds the distance-geometry bounds matrix from the molecular graph, following RDKit's
/// <c>DGeomHelpers::setTopolBounds</c>: 1-2 bounds from UFF rest lengths (±0.01 Å), 1-3 bounds from ideal angles
/// (with ring special cases) via the law of cosines, 1-4 bounds from the cis/trans range of each torsion class,
/// and van der Waals lower bounds for everything else (scaled for 1-5 / 1-6 pairs and H-bonded pairs).
/// </summary>
/// <remarks>RDKit's optional 1-5 bounds and macrocycle-specific 1-4 heuristics are not implemented.</remarks>
internal sealed class DistanceBoundsBuilder
{
    private const double Dist12Delta = 0.01;
    private const double AngleDelta = 0.035;
    private const double GenDistTol = 0.06;
    private const double VdwScale15 = 0.7;
    private const double HBondLength = 1.8;
    private const double MaxUpper = 1000.0;

    // RDKit periodic-table van der Waals radii (Z = 0..92).
    private static readonly double[] Rvdw =
    [
        0.0, 1.2, 1.4, 2.2, 1.9, 1.8, 1.7, 1.6, 1.55, 1.5, 1.54, 2.4, 2.2, 2.1, 2.1, 1.95, 1.8, 1.8, 1.88, 2.8, 2.4, 2.3, 2.15,
        2.05, 2.05, 2.05, 2.05, 2.0, 2.0, 2.0, 2.1, 2.1, 2.1, 2.05, 1.9, 1.9, 2.02, 2.9, 2.55, 2.4, 2.3, 2.15, 2.1, 2.05, 2.05, 2.0,
        2.05, 2.1, 2.2, 2.2, 2.25, 2.2, 2.1, 2.1, 2.16, 3.0, 2.7, 2.5, 2.48, 2.47, 2.45, 2.43, 2.42, 2.4, 2.38, 2.37, 2.35, 2.33,
        2.32, 2.3, 2.28, 2.27, 2.25, 2.2, 2.1, 2.05, 2.0, 2.0, 2.05, 2.1, 2.05, 2.2, 2.3, 2.3, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.4, 2.0, 2.3,
    ];

    internal static double VdwRadius(int z) => (uint)z < (uint)Rvdw.Length ? Rvdw[z] : 2.0;

    private readonly Molecule _mol;
    private readonly Hybridization[] _hyb;
    private readonly bool[] _conj;
    private readonly int _n, _nb;
    private readonly int[,] _dmat;
    private readonly double[] _bondLength;
    private readonly double[] _bondAngle; // nb*nb, -1 unset
    private readonly int[] _bondAdj;      // nb*nb, -1 unset
    private readonly bool[] _v12, _v13, _v14;
    private readonly bool _forceTransAmides;

    public List<Path14> Paths14 { get; } = [];

    public DistanceBoundsBuilder(Molecule mol, Hybridization[] hyb, bool[] conjugated, bool forceTransAmides)
    {
        _mol = mol;
        _hyb = hyb;
        _conj = conjugated;
        _n = mol.AtomCount;
        _nb = mol.BondCount;
        _dmat = mol.TopologicalDistances();
        _bondLength = new double[_nb];
        _bondAngle = new double[_nb * _nb];
        Array.Fill(_bondAngle, -1.0);
        _bondAdj = new int[_nb * _nb];
        Array.Fill(_bondAdj, -1);
        _v12 = new bool[_n * _n];
        _v13 = new bool[_n * _n];
        _v14 = new bool[_n * _n];
        _forceTransAmides = forceTransAmides;
    }

    public BoundsMatrix Build()
    {
        var m = new BoundsMatrix(_n);
        Set12(m);
        Set13(m);
        Set14(m);
        SetLowerVdw(m);
        return m;
    }

    private int Pid(int a, int b) => a < b ? a * _n + b : b * _n + a;
    private int Bp(int a, int b) => a < b ? a * _nb + b : b * _nb + a;
    private double BondAngle(int b1, int b2) => _bondAngle[Bp(b1, b2)];
    private int BondAdj(int b1, int b2) => _bondAdj[Bp(b1, b2)];

    private static void CheckAndSet(BoundsMatrix m, int i, int j, double lb, double ub)
    {
        double clb = m.Lower(i, j), cub = m.Upper(i, j);
        if (ub <= lb) ub = lb + 1e-3;
        if (clb <= Dist12Delta) m.SetLower(i, j, lb);
        else if (lb < clb && lb > Dist12Delta) m.SetLower(i, j, lb);
        if (cub >= MaxUpper) m.SetUpper(i, j, ub);
        else if (ub > cub && ub < MaxUpper) m.SetUpper(i, j, ub);
    }

    // ---- 1-2 -----------------------------------------------------------------------------------------

    private void Set12(BoundsMatrix m)
    {
        var types = UffForceField.AssignAtomTypes(_mol);
        var squish = new bool[_n];
        var rings = _mol.Rings;
        var bonds = _mol.Bonds;
        for (int b = 0; b < _nb; b++)
        {
            var bd = bonds[b];
            if (_conj[b] && (_mol.Atoms[bd.Begin].AtomicNumber > 10 || _mol.Atoms[bd.End].AtomicNumber > 10) && rings.IsBondInRingOfSize(b, 5))
                squish[bd.Begin] = squish[bd.End] = true;
        }
        for (int b = 0; b < _nb; b++)
        {
            var bd = bonds[b];
            int i = bd.Begin, j = bd.End;
            double bo = UffForceField.BondOrderValue(bd);
            if (types[i] is { } ti && types[j] is { } tj && UffParameters.TryGet(ti, out var pi) && UffParameters.TryGet(tj, out var pj))
            {
                double bl = UffForceField.RestLength(bo, pi, pj);
                double extra = squish[i] && squish[j] ? 0.2 : 0.0;
                _bondLength[b] = bl;
                m.SetUpper(i, j, bl + extra + Dist12Delta);
                m.SetLower(i, j, bl - extra - Dist12Delta);
            }
            else
            {
                double bl = PeriodicTable.Get(_mol.Atoms[i].AtomicNumber).CovalentRadius + PeriodicTable.Get(_mol.Atoms[j].AtomicNumber).CovalentRadius;
                bl -= 0.1332 * Math.Log(bo) * bl;
                _bondLength[b] = bl;
                m.SetUpper(i, j, 1.1 * bl);
                m.SetLower(i, j, 0.9 * bl);
            }
            _v12[Pid(i, j)] = true;
        }
    }

    // ---- 1-3 -----------------------------------------------------------------------------------------

    private bool IsLargerSp2(int a) => _mol.Atoms[a].AtomicNumber > 13 && _hyb[a] == Hybridization.SP2 && _mol.Rings.IsAtomInRing(a);

    private double ScaledTolerance(int a1, int a2, int a3)
    {
        double t = AngleDelta;
        if (IsLargerSp2(a1)) t *= 2;
        if (IsLargerSp2(a2)) t *= 2;
        if (IsLargerSp2(a3)) t *= 2;
        return t;
    }

    private static (double L, double U) AngleRange(double angle, double tol)
    {
        if (angle + tol >= Math.PI) return (Math.PI - tol, Math.PI);
        if (angle - tol <= 0) return (0, tol);
        return (angle - tol, angle + tol);
    }

    private static double Dist13(double d1, double d2, double angle) => Math.Sqrt(d1 * d1 + d2 * d2 - 2 * d1 * d2 * Math.Cos(angle));

    private void Set13Helper(BoundsMatrix m, int a1, int a, int a3, double angle)
    {
        var (lo, hi) = AngleRange(angle, ScaledTolerance(a1, a, a3));
        double du = Dist13(m.Upper(a1, a), m.Upper(a, a3), hi);
        double dl = Dist13(m.Lower(a1, a), m.Lower(a, a3), lo);
        CheckAndSet(m, a1, a3, dl, du);
    }

    private double RingAngle(int atom, int ringSize)
    {
        var h = _hyb[atom];
        if ((h == Hybridization.SP2 && ringSize <= 8) || ringSize == 3 || ringSize == 4) return Math.PI * (1.0 - 2.0 / ringSize);
        if (h == Hybridization.SP3) return (ringSize == 5 ? 104 : 109.5) * FFMath.Deg2Rad;
        if (h == Hybridization.SP3D) return 105.0 * FFMath.Deg2Rad;
        if (h == Hybridization.SP3D2) return 90.0 * FFMath.Deg2Rad;
        return 120 * FFMath.Deg2Rad;
    }

    private void Set13(BoundsMatrix m)
    {
        var rings = _mol.Rings;
        var atomRings = rings.AtomRings.OrderBy(r => r.Length).ToList();
        var visited = new int[_n];
        var angleTaken = new double[_n];
        var donePaths = new HashSet<int>();
        foreach (var ring in atomRings)
        {
            int r = ring.Length;
            int a1 = ring[r - 1];
            for (int i = 0; i < r; i++)
            {
                int a2 = ring[i], a3 = ring[(i + 1) % r];
                int b1 = _mol.BondBetween(a1, a2), b2 = _mol.BondBetween(a2, a3);
                int key = Bp(b1, b2);
                if (donePaths.Add(key))
                {
                    double angle = RingAngle(a2, r);
                    int pid = Pid(a1, a3);
                    if (!_v12[pid])
                    {
                        Set13Helper(m, a1, a2, a3, angle);
                        _v13[pid] = true;
                    }
                    _bondAngle[key] = angle;
                    _bondAdj[key] = a2;
                    visited[a2]++;
                    angleTaken[a2] += angle;
                }
                a1 = a2;
            }
        }
        for (int a2 = 0; a2 < _n; a2++)
        {
            int deg = _mol.Degree(a2);
            int n13 = deg * (deg - 1) / 2;
            if (n13 == visited[a2]) continue;
            var h = _hyb[a2];
            var nb = _mol.Neighbors(a2).ToArray();
            var nbb = _mol.NeighborBonds(a2).ToArray();
            for (int x = 0; x < nb.Length; x++)
            {
                for (int y = 0; y < x; y++)
                {
                    int b1 = nbb[x], b2 = nbb[y], a1 = nb[x], a3 = nb[y];
                    int key = Bp(b1, b2);
                    double angle;
                    if (visited[a2] >= 1)
                    {
                        if (_bondAngle[key] >= 0) continue;
                        if (h == Hybridization.SP2) angle = (2 * Math.PI - angleTaken[a2]) / (n13 - visited[a2]);
                        else if (h == Hybridization.SP3)
                        {
                            angle = 109.5 * FFMath.Deg2Rad;
                            if (rings.IsAtomInRingOfSize(a2, 3)) angle = 116.0 * FFMath.Deg2Rad;
                            else if (rings.IsAtomInRingOfSize(a2, 4)) angle = 112.0 * FFMath.Deg2Rad;
                        }
                        else angle = (deg == 5 ? 105.0 : deg == 6 ? 135.0 : 120.0) * FFMath.Deg2Rad;
                        int pid = Pid(a1, a3);
                        if (!_v12[pid]) { Set13Helper(m, a1, a2, a3, angle); _v13[pid] = true; }
                    }
                    else
                    {
                        angle = h switch
                        {
                            Hybridization.SP => Math.PI,
                            Hybridization.SP2 => 2 * Math.PI / 3,
                            Hybridization.SP3 => 109.5 * FFMath.Deg2Rad,
                            Hybridization.SP3D => 105.0 * FFMath.Deg2Rad,
                            Hybridization.SP3D2 => 135.0 * FFMath.Deg2Rad,
                            _ => 120.0 * FFMath.Deg2Rad,
                        };
                        int pid = Pid(a1, a3);
                        if (!_v12[pid])
                        {
                            if (deg <= 4) Set13Helper(m, a1, a2, a3, angle);
                            else CheckAndSet(m, a1, a3, 1.0, 1.2 * (_bondLength[b1] + _bondLength[b2]));
                            _v13[pid] = true;
                        }
                    }
                    _bondAngle[key] = angle;
                    _bondAdj[key] = a2;
                    angleTaken[a2] += angle;
                    visited[a2]++;
                }
            }
        }
    }

    // ---- 1-4 -----------------------------------------------------------------------------------------

    private enum Type14 { InChain, InRing, ShareRingBond, TwoInDiffRing, TwoInSameRing }

    private BondStereo StereoAlongPath(int bond, int a1, int a2, int a3, int a4)
    {
        var b = _mol.Bonds[bond];
        if (b.Stereo is not (BondStereo.Cis or BondStereo.Trans) || b.StereoAtomBegin < 0 || b.StereoAtomEnd < 0) return BondStereo.None;
        int s2 = b.Begin == a2 ? b.StereoAtomBegin : b.StereoAtomEnd;
        int s3 = b.Begin == a2 ? b.StereoAtomEnd : b.StereoAtomBegin;
        int flips = (s2 != a1 ? 1 : 0) + (s3 != a4 ? 1 : 0);
        if (flips % 2 == 0) return b.Stereo;
        return b.Stereo == BondStereo.Cis ? BondStereo.Trans : BondStereo.Cis;
    }

    private TorsionClass InRingType(int bond2, int a1, int a2, int a3, int a4, int ringSize)
    {
        var st = StereoAlongPath(bond2, a1, a2, a3, a4);
        bool sp2 = _hyb[a2] == Hybridization.SP2 && _hyb[a3] == Hybridization.SP2;
        if ((ringSize <= 8 && sp2 && st != BondStereo.Trans) || st == BondStereo.Cis) return TorsionClass.Cis;
        if (st == BondStereo.Trans) return TorsionClass.Trans;
        return sp2 ? TorsionClass.CisTrans : TorsionClass.Flexible;
    }

    private TorsionClass TwoInSameRingType(int bond2, int a1, int a2, int a3, int a4, bool preferTrans)
    {
        var st = StereoAlongPath(bond2, a1, a2, a3, a4);
        bool sp2 = _hyb[a2] == Hybridization.SP2 && _hyb[a3] == Hybridization.SP2;
        if ((preferTrans && sp2 && st != BondStereo.Cis) || st == BondStereo.Trans) return TorsionClass.Trans;
        if (st == BondStereo.Cis) return TorsionClass.Cis;
        return sp2 ? TorsionClass.CisTrans : TorsionClass.Flexible;
    }

    private bool AmideEster14(int b1, int b3, int a2, int a3, int a4)
    {
        var bonds = _mol.Bonds;
        int z2 = _mol.Atoms[a2].AtomicNumber, z3 = _mol.Atoms[a3].AtomicNumber, z4 = _mol.Atoms[a4].AtomicNumber;
        return z3 == 6 && bonds[b3].Order == BondOrder.Double && (z4 is 8 or 7) && bonds[b1].Order == BondOrder.Single
            && (z2 == 8 || (z2 == 7 && _mol.TotalHydrogens(a2) == 1));
    }

    private bool IsCarbonyl(int a)
    {
        if (_mol.Atoms[a].AtomicNumber != 6 || _mol.Degree(a) <= 2) return false;
        var nb = _mol.Neighbors(a);
        var nbb = _mol.NeighborBonds(a);
        for (int k = 0; k < nb.Length; k++)
            if (_mol.Atoms[nb[k]].AtomicNumber is 8 or 7 && _mol.Bonds[nbb[k]].Order == BondOrder.Double) return true;
        return false;
    }

    private bool AmideEster15(int b1, int b3, int a2, int a3)
    {
        int z2 = _mol.Atoms[a2].AtomicNumber;
        return (z2 == 8 || (z2 == 7 && _mol.TotalHydrogens(a2) == 1)) && _mol.Bonds[b1].Order == BondOrder.Single
            && _mol.Atoms[a3].AtomicNumber == 6 && _mol.Bonds[b3].Order == BondOrder.Single && IsCarbonyl(a3);
    }

    private bool SecondaryAmideH(int aH, int aN)
        => _mol.Atoms[aH].AtomicNumber == 1 && _mol.Atoms[aN].AtomicNumber == 7 && _mol.Degree(aN) == 3 && _mol.TotalHydrogens(aN) == 1;

    private (TorsionClass Class, bool Forced) ChainType(int b1, int b2, int b3, int a1, int a2, int a3, int a4)
    {
        var bonds = _mol.Bonds;
        switch (bonds[b2].Order)
        {
            case BondOrder.Double:
                if (bonds[b1].Order == BondOrder.Double || bonds[b3].Order == BondOrder.Double) return (TorsionClass.Cis, false);
                var st = StereoAlongPath(b2, a1, a2, a3, a4);
                if (st == BondStereo.Cis) return (TorsionClass.Cis, false);
                if (st == BondStereo.Trans) return (TorsionClass.Trans, false);
                if (_hyb[a2] == Hybridization.SP2 && _hyb[a3] == Hybridization.SP2) return (TorsionClass.CisTrans, false);
                return (TorsionClass.Flexible, false);
            case BondOrder.Single:
                if (_mol.Atoms[a2].AtomicNumber == 16 && _mol.Atoms[a3].AtomicNumber == 16 && _mol.Degree(a2) == 2 && _mol.Degree(a3) == 2)
                    return (TorsionClass.Custom, false);
                if (AmideEster14(b1, b3, a2, a3, a4) || AmideEster14(b3, b1, a3, a2, a1))
                {
                    if (!_forceTransAmides) return (TorsionClass.CisTrans, false);
                    return SecondaryAmideH(a1, a2) || SecondaryAmideH(a4, a3) ? (TorsionClass.Trans, true) : (TorsionClass.Cis, true);
                }
                if (AmideEster15(b1, b3, a2, a3) || AmideEster15(b3, b1, a3, a2))
                {
                    if (!_forceTransAmides) return (TorsionClass.CisTrans, false);
                    return SecondaryAmideH(a1, a2) || SecondaryAmideH(a4, a3) ? (TorsionClass.Cis, true) : (TorsionClass.Trans, true);
                }
                return (TorsionClass.Flexible, false);
            default:
                return (TorsionClass.Flexible, false);
        }
    }

    private static double Dist14Cis(double d1, double d2, double d3, double a12, double a23)
    {
        double dx = d2 - d3 * Math.Cos(a23) - d1 * Math.Cos(a12);
        double dy = d3 * Math.Sin(a23) - d1 * Math.Sin(a12);
        return Math.Sqrt(dx * dx + dy * dy);
    }

    private static double Dist14Trans(double d1, double d2, double d3, double a12, double a23)
    {
        double dx = d2 - d3 * Math.Cos(a23) - d1 * Math.Cos(a12);
        double dy = d3 * Math.Sin(a23) + d1 * Math.Sin(a12);
        return Math.Sqrt(dx * dx + dy * dy);
    }

    private static double Dist14(double d1, double d2, double d3, double a12, double a23, double tor)
    {
        double p1x = d1 * Math.Cos(a12), p1y = d1 * Math.Sin(a12);
        double p4x = d2 - d3 * Math.Cos(a23), p4y = d3 * Math.Sin(a23) * Math.Cos(tor), p4z = d3 * Math.Sin(a23) * Math.Sin(tor);
        double dx = p4x - p1x, dy = p4y - p1y;
        return Math.Sqrt(dx * dx + dy * dy + p4z * p4z);
    }

    private void Collect14(BoundsMatrix m, int b1, int b2, int b3, Type14 type, int ringSize, bool preferTrans,
        Dictionary<int, List<(double L, double U, int A1, int A4)>> collected)
    {
        int a2 = BondAdj(b1, b2), a3 = BondAdj(b2, b3);
        if (a2 < 0 || a3 < 0) return;
        int a1 = _mol.Bonds[b1].Other(a2), a4 = _mol.Bonds[b3].Other(a3);
        int pid = Pid(a1, a4);
        if (_v12[pid] || _v13[pid] || _dmat[a1, a4] < 3) return;
        double ba12 = BondAngle(b1, b2), ba23 = BondAngle(b2, b3);
        if (ba12 <= 0 || ba23 <= 0) return;
        if (type != Type14.InChain && (_mol.BondBetween(a1, a3) >= 0 || _mol.BondBetween(a4, a2) >= 0)) return;
        double blU1 = m.Upper(a1, a2), blU2 = m.Upper(a2, a3), blU3 = m.Upper(a3, a4);
        double blL1 = m.Lower(a1, a2), blL2 = m.Lower(a2, a3), blL3 = m.Lower(a3, a4);
        (TorsionClass cls, bool forced) = type switch
        {
            Type14.InChain => ChainType(b1, b2, b3, a1, a2, a3, a4),
            Type14.InRing => (InRingType(b2, a1, a2, a3, a4, ringSize), false),
            Type14.TwoInSameRing => (TwoInSameRingType(b2, a1, a2, a3, a4, preferTrans), false),
            _ => (InRingType(b2, a1, a2, a3, a4, 0), false),
        };
        var (l12, u12) = AngleRange(ba12, ScaledTolerance(a1, a2, a3));
        var (l23, u23) = AngleRange(ba23, ScaledTolerance(a2, a3, a4));
        double dl, du;
        switch (cls)
        {
            case TorsionClass.Cis:
                dl = Dist14Cis(blL1, blL2, blL3, l12, l23);
                du = Dist14Cis(blU1, blU2, blU3, u12, u23);
                break;
            case TorsionClass.Trans:
                dl = Dist14Trans(blL1, blL2, blL3, l12, l23);
                du = Dist14Trans(blU1, blU2, blU3, u12, u23);
                break;
            case TorsionClass.Custom:
                dl = Dist14(blL1, blL2, blL3, l12, l23, Math.PI / 2);
                du = Dist14(blU1, blU2, blU3, u12, u23, Math.PI / 2);
                break;
            default:
                dl = Dist14Cis(blL1, blL2, blL3, l12, l23);
                du = Dist14Trans(blU1, blU2, blU3, u12, u23);
                if (du < dl) (dl, du) = (du, dl);
                if (Math.Abs(du - dl) < Dist12Delta) { dl -= GenDistTol; du += GenDistTol; }
                break;
        }
        if (Math.Abs(du - dl) < Dist12Delta) { dl -= GenDistTol; du += GenDistTol; }
        Paths14.Add(new Path14(a1, a2, a3, a4, cls, forced));
        _v14[pid] = true;
        if (!collected.TryGetValue(pid, out var list)) collected[pid] = list = [];
        list.Add((dl, du, a1, a4));
    }

    private void Set14(BoundsMatrix m)
    {
        var rings = _mol.Rings;
        var bondRings = rings.BondRings.OrderByDescending(r => r.Length).ToList();
        var ringBondPairs = new HashSet<int>();
        var cisRingBondPairs = new HashSet<int>();
        var donePaths = new HashSet<long>();
        var collected = new Dictionary<int, List<(double L, double U, int A1, int A4)>>();
        long Triple(int x, int y, int z) => x < z ? ((long)x * _nb + y) * _nb + z : ((long)z * _nb + y) * _nb + x;
        foreach (var ring in bondRings)
        {
            int r = ring.Length;
            if (r < 3) continue;
            int b1 = ring[r - 1];
            for (int i = 0; i < r; i++)
            {
                int b2 = ring[i], b3 = ring[(i + 1) % r];
                int pid = Bp(b1, b2);
                ringBondPairs.Add(pid);
                donePaths.Add(Triple(b1, b2, b3));
                if (r > 5)
                {
                    Collect14(m, b1, b2, b3, Type14.InRing, r, false, collected);
                    if (r <= 8) cisRingBondPairs.Add(pid);
                }
                else
                {
                    RecordRingPath(b1, b2, b3);
                    cisRingBondPairs.Add(pid);
                }
                b1 = b2;
            }
        }
        var bonds = _mol.Bonds;
        for (int b2 = 0; b2 < _nb; b2++)
        {
            int a2 = bonds[b2].Begin, a3 = bonds[b2].End;
            foreach (int b1 in _mol.NeighborBonds(a2))
            {
                if (b1 == b2) continue;
                foreach (int b3 in _mol.NeighborBonds(a3))
                {
                    if (b3 == b2) continue;
                    if (donePaths.Contains(Triple(b1, b2, b3))) continue;
                    int p1 = Bp(b1, b2), p2 = Bp(b2, b3);
                    if (ringBondPairs.Contains(p1) || ringBondPairs.Contains(p2))
                        Collect14(m, b1, b2, b3, Type14.TwoInSameRing, 0, cisRingBondPairs.Contains(p1) || cisRingBondPairs.Contains(p2), collected);
                    else if ((rings.IsBondInRing(b1) && rings.IsBondInRing(b2)) || (rings.IsBondInRing(b2) && rings.IsBondInRing(b3)))
                        Collect14(m, b1, b2, b3, Type14.TwoInDiffRing, 0, false, collected);
                    else if (rings.IsBondInRing(b2))
                        Collect14(m, b1, b2, b3, Type14.ShareRingBond, 0, false, collected);
                    else
                        Collect14(m, b1, b2, b3, Type14.InChain, 0, false, collected);
                }
            }
        }
        foreach (var (_, list) in collected)
        {
            var (lo, hi) = Merge(list);
            CheckAndSet(m, list[0].A1, list[0].A4, lo, hi);
        }
    }

    private void RecordRingPath(int b1, int b2, int b3)
    {
        int a2 = BondAdj(b1, b2), a3 = BondAdj(b2, b3);
        if (a2 < 0 || a3 < 0) return;
        int a1 = _mol.Bonds[b1].Other(a2), a4 = _mol.Bonds[b3].Other(a3);
        var cls = _hyb[a2] == Hybridization.SP2 && _hyb[a3] == Hybridization.SP2 ? TorsionClass.Cis : TorsionClass.Flexible;
        Paths14.Add(new Path14(a1, a2, a3, a4, cls, false));
    }

    /// <summary>RDKit's merge of several 1-4 ranges for the same pair (union of maximal intersections).</summary>
    private static (double Lower, double Upper) Merge(List<(double L, double U, int A1, int A4)> bounds)
    {
        var sorted = bounds.OrderBy(b => b.L).ToList();
        double curL = sorted[0].L, curU = sorted[0].U, componentUpper = curU;
        double? resultLower = null;
        for (int k = 1; k < sorted.Count; k++)
        {
            var b = sorted[k];
            if (b.L <= curU)
            {
                curL = b.L;
                curU = Math.Min(curU, b.U);
            }
            else
            {
                resultLower ??= curL;
                curL = b.L;
                curU = b.L <= componentUpper ? Math.Min(componentUpper, b.U) : b.U;
            }
            componentUpper = Math.Max(componentUpper, b.U);
        }
        return (resultLower ?? curL, curU);
    }

    // ---- van der Waals lower bounds ----------------------------------------------------------------------

    private void SetLowerVdw(BoundsMatrix m)
    {
        var hDonorH = new bool[_n];
        var acceptor = new bool[_n];
        for (int i = 0; i < _n; i++)
        {
            int z = _mol.Atoms[i].AtomicNumber;
            acceptor[i] = z is 7 or 8;
            if (z == 1) foreach (int nb in _mol.Neighbors(i)) if (_mol.Atoms[nb].AtomicNumber is 7 or 8) hDonorH[i] = true;
        }
        for (int i = 1; i < _n; i++)
        {
            double v1 = VdwRadius(_mol.Atoms[i].AtomicNumber);
            for (int j = 0; j < i; j++)
            {
                if (m.Lower(i, j) >= Dist12Delta) continue;
                double v2 = VdwRadius(_mol.Atoms[j].AtomicNumber);
                int d = _dmat[i, j];
                if ((hDonorH[i] && acceptor[j]) || (acceptor[i] && hDonorH[j])) m.SetLower(i, j, HBondLength);
                else if (d == 4) m.SetLower(i, j, VdwScale15 * (v1 + v2));
                else if (d == 5) m.SetLower(i, j, (VdwScale15 + 0.5 * (1.0 - VdwScale15)) * (v1 + v2));
                else m.SetLower(i, j, v1 + v2);
            }
        }
    }
}
