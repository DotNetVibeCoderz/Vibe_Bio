namespace Biotech.Net.Chemistry;

/// <summary>Energy contributions of a UFF evaluation (kcal/mol).</summary>
/// <param name="BondStretch">Harmonic bond stretching.</param>
/// <param name="AngleBend">Cosine-Fourier angle bending (plus RDKit's small-angle penalty).</param>
/// <param name="Torsion">Torsional terms.</param>
/// <param name="Inversion">Out-of-plane (inversion) terms.</param>
/// <param name="VanDerWaals">Lennard-Jones 12-6 van der Waals for atom pairs separated by three or more bonds.</param>
public readonly record struct UffEnergyTerms(double BondStretch, double AngleBend, double Torsion, double Inversion, double VanDerWaals)
{
    /// <summary>Sum of all terms.</summary>
    public double Total => BondStretch + AngleBend + Torsion + Inversion + VanDerWaals;
}

/// <summary>
/// The Universal Force Field (Rappé et al., J. Am. Chem. Soc. 1992, 114, 10024) as implemented in RDKit: atom typing
/// from hybridisation/aromaticity/conjugation, harmonic bonds with the Pauling and electronegativity corrections,
/// cosine-Fourier angles (with RDKit's special cases for linear, trigonal, square-planar and 3-/4-ring sp² centres),
/// torsions, inversions and a Lennard-Jones van der Waals term. Electrostatics are not included (RDKit's default).
/// </summary>
/// <remarks>
/// Build the field from a molecule with <b>explicit hydrogens</b> (<see cref="Molecule.AddHydrogens"/>). Energies and
/// gradients are verified against RDKit 2026.09 in the test suite (see the documentation for measured deviations).
/// Not supported: trigonal-bipyramidal (sp³d, 5-coordinate) angle special-casing and charge-state flags for most metals.
/// </remarks>
public sealed class UffForceField
{
    private const double G = 332.06;
    private const double Lambda = 0.1332;

    private readonly int _n;
    private readonly string[] _types;
    private readonly (int I, int J, double K, double R0)[] _bonds;
    private readonly AngleTerm[] _angles;
    private readonly (int I, int J, int K, int L, double V, int Order, double CosTerm)[] _torsions;
    private readonly (int I, int J, int K, int L, double K0, double C0, double C1, double C2)[] _inversions;
    private readonly (int I, int J, double X, double D, double Cut2)[] _vdw;

    private readonly record struct AngleTerm(int I, int J, int K, double Ka, double Theta0, int Order, double C0, double C1, double C2);

    private UffForceField(int n, string[] types, List<(int, int, double, double)> bonds, List<AngleTerm> angles,
        List<(int, int, int, int, double, int, double)> torsions, List<(int, int, int, int, double, double, double, double)> inversions,
        List<(int, int, double, double, double)> vdw)
    {
        _n = n;
        _types = types;
        _bonds = [.. bonds];
        _angles = [.. angles];
        _torsions = [.. torsions];
        _inversions = [.. inversions];
        _vdw = [.. vdw];
    }

    /// <summary>Number of atoms the field was built for.</summary>
    public int AtomCount => _n;

    /// <summary>UFF atom type of every atom (e.g. <c>C_R</c>, <c>N_3</c>, <c>H_</c>).</summary>
    public IReadOnlyList<string> AtomTypes => _types;

    /// <summary>Number of bond, angle, torsion, inversion and van der Waals terms.</summary>
    public (int Bonds, int Angles, int Torsions, int Inversions, int VanDerWaals) TermCounts
        => (_bonds.Length, _angles.Length, _torsions.Length, _inversions.Length, _vdw.Length);

    /// <summary>Assigns UFF atom types (RDKit rules). Atoms that cannot be typed get <c>null</c>.</summary>
    public static string?[] AssignAtomTypes(Molecule mol)
    {
        var conj = AtomHybridization.ConjugatedBonds(mol);
        var hyb = AtomHybridization.Compute(mol, conj);
        var labels = new string?[mol.AtomCount];
        for (int i = 0; i < mol.AtomCount; i++)
        {
            string label = Label(mol, i, hyb[i], AtomHybridization.HasConjugatedBond(mol, i, conj));
            labels[i] = UffParameters.TryGet(label, out _) ? label : null;
        }
        return labels;
    }

    /// <summary>Builds the force field, or returns false (with a message) when an atom has no UFF type.</summary>
    public static bool TryCreate(Molecule mol, out UffForceField? forceField, out string? error)
    {
        var conj = AtomHybridization.ConjugatedBonds(mol);
        var hyb = AtomHybridization.Compute(mol, conj);
        int n = mol.AtomCount;
        var types = new string[n];
        var p = new UffAtomParameters[n];
        for (int i = 0; i < n; i++)
        {
            types[i] = Label(mol, i, hyb[i], AtomHybridization.HasConjugatedBond(mol, i, conj));
            if (!UffParameters.TryGet(types[i], out p[i]))
            {
                forceField = null;
                error = $"No UFF parameters for atom {i} ({mol.Atoms[i].Symbol}, type '{types[i]}').";
                return false;
            }
        }
        var bonds = mol.Bonds;
        var bondTerms = new List<(int, int, double, double)>();
        foreach (var b in bonds)
        {
            double r0 = RestLength(BondOrderValue(b), p[b.Begin], p[b.End]);
            bondTerms.Add((b.Begin, b.End, 2.0 * G * p[b.Begin].Z1 * p[b.End].Z1 / (r0 * r0 * r0), r0));
        }

        var rings = mol.Rings;
        var angleTerms = new List<AngleTerm>();
        for (int j = 0; j < n; j++)
        {
            var nb = mol.Neighbors(j);
            var nbb = mol.NeighborBonds(j);
            if (nb.Length < 2) continue;
            if (hyb[j] == Hybridization.SP3D && nb.Length == 5) continue; // trigonal bipyramid: not supported
            for (int x = 0; x < nb.Length; x++)
                for (int y = x + 1; y < nb.Length; y++)
                {
                    int i = nb[x], k = nb[y];
                    int order = hyb[j] switch
                    {
                        Hybridization.SP => 1,
                        Hybridization.SP2 => 3,
                        Hybridization.SP3D2 => 4,
                        _ => 0,
                    };
                    double theta0 = p[j].Theta0 * FFMath.Deg2Rad;
                    if (hyb[j] == Hybridization.SP2)
                    {
                        if (rings.IsAtomInRingOfSize(j, 3))
                        {
                            bool ri = rings.IsAtomInRingOfSize(i, 3), rk = rings.IsAtomInRingOfSize(k, 3);
                            if (ri != rk) { order = 0; theta0 = 150 * FFMath.Deg2Rad; }
                            else if (ri) { order = 0; theta0 = 60 * FFMath.Deg2Rad; }
                        }
                        else if (rings.IsAtomInRingOfSize(j, 4))
                        {
                            bool ri = rings.IsAtomInRingOfSize(i, 4), rk = rings.IsAtomInRingOfSize(k, 4);
                            if (ri != rk) { order = 0; theta0 = 135 * FFMath.Deg2Rad; }
                            else if (ri) { order = 0; theta0 = 90 * FFMath.Deg2Rad; }
                        }
                    }
                    double ka = AngleForceConstant(theta0, BondOrderValue(bonds[nbb[x]]), BondOrderValue(bonds[nbb[y]]), p[i], p[j], p[k]);
                    double c0 = 0, c1 = 0, c2 = 0;
                    if (order == 0)
                    {
                        double s = Math.Sin(theta0), c = Math.Cos(theta0);
                        c2 = 1.0 / (4.0 * Math.Max(s * s, 1e-8));
                        c1 = -4.0 * c2 * c;
                        c0 = c2 * (2.0 * c * c + 1.0);
                    }
                    angleTerms.Add(new AngleTerm(i, j, k, ka, theta0, order, c0, c1, c2));
                }
        }

        var torsionTerms = new List<(int, int, int, int, double, int, double)>();
        for (int bi = 0; bi < bonds.Length; bi++)
        {
            var b = bonds[bi];
            int j = b.Begin, k = b.End;
            if (mol.Degree(j) < 2 || mol.Degree(k) < 2 || HasTripleBond(mol, j) || HasTripleBond(mol, k)) continue;
            if (hyb[j] is not (Hybridization.SP2 or Hybridization.SP3) || hyb[k] is not (Hybridization.SP2 or Hybridization.SP3)) continue;
            int start = torsionTerms.Count;
            foreach (int i in mol.Neighbors(j))
            {
                if (i == k) continue;
                foreach (int l in mol.Neighbors(k))
                {
                    if (l == j || l == i) continue;
                    bool endSp2 = hyb[i] == Hybridization.SP2 || hyb[l] == Hybridization.SP2;
                    var (v, order, cosTerm) = TorsionParams(BondOrderValue(b), mol.Atoms[j].AtomicNumber, mol.Atoms[k].AtomicNumber, hyb[j], hyb[k], p[j], p[k], endSp2);
                    torsionTerms.Add((i, j, k, l, v, order, cosTerm));
                }
            }
            int count = torsionTerms.Count - start;
            for (int t = start; t < torsionTerms.Count; t++)
            {
                var tt = torsionTerms[t];
                torsionTerms[t] = tt with { Item5 = tt.Item5 / count };
            }
        }

        var inversionTerms = new List<(int, int, int, int, double, double, double, double)>();
        for (int j = 0; j < n; j++)
        {
            int z = mol.Atoms[j].AtomicNumber;
            if (z is not (6 or 7 or 8 or 15 or 33 or 51 or 83) || mol.Degree(j) != 3) continue;
            if (z is 6 or 7 or 8 && hyb[j] != Hybridization.SP2) continue;
            var nb = mol.Neighbors(j);
            bool boundToSp2O = false;
            if (z == 6)
                foreach (int x in nb) if (mol.Atoms[x].AtomicNumber == 8 && hyb[x] == Hybridization.SP2) boundToSp2O = true;
            var (k0, c0, c1, c2) = InversionParams(z, boundToSp2O);
            inversionTerms.Add((nb[0], j, nb[1], nb[2], k0, c0, c1, c2));
            inversionTerms.Add((nb[0], j, nb[2], nb[1], k0, c0, c1, c2));
            inversionTerms.Add((nb[1], j, nb[2], nb[0], k0, c0, c1, c2));
        }

        // van der Waals: pairs separated by >= 3 bonds, same fragment (RDKit: ignoreInterfragInteractions).
        var frag = new int[n];
        var frags = mol.Fragments();
        for (int f = 0; f < frags.Count; f++) foreach (int a in frags[f]) frag[a] = f;
        var close = new HashSet<long>();
        foreach (var b in bonds) close.Add(PairKey(b.Begin, b.End, n));
        for (int j = 0; j < n; j++)
        {
            var nb = mol.Neighbors(j);
            for (int x = 0; x < nb.Length; x++)
                for (int y = x + 1; y < nb.Length; y++) close.Add(PairKey(nb[x], nb[y], n));
        }
        var vdwTerms = new List<(int, int, double, double, double)>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                if (frag[i] != frag[j] || close.Contains(PairKey(i, j, n))) continue;
                double xij = Math.Sqrt(p[i].X1 * p[j].X1), dij = Math.Sqrt(p[i].D1 * p[j].D1);
                double cut = 10.0 * xij;
                vdwTerms.Add((i, j, xij, dij, cut * cut));
            }

        forceField = new UffForceField(n, types, bondTerms, angleTerms, torsionTerms, inversionTerms, vdwTerms);
        error = null;
        return true;
    }

    /// <summary>Builds the force field.</summary>
    /// <exception cref="InvalidOperationException">An atom has no UFF type.</exception>
    public static UffForceField Create(Molecule mol)
        => TryCreate(mol, out var ff, out var error) ? ff! : throw new InvalidOperationException(error);

    /// <summary>Total energy (kcal/mol) at <paramref name="coordinates"/>.</summary>
    public double Energy(ReadOnlySpan<Point3> coordinates) => Evaluate(Flatten(coordinates), [], null);

    /// <summary>Total energy and its gradient (kcal/mol/Å) at <paramref name="coordinates"/>.</summary>
    public double Energy(ReadOnlySpan<Point3> coordinates, Span<Point3> gradient)
    {
        if (gradient.Length != _n) throw new ArgumentException("Gradient length must equal the atom count.", nameof(gradient));
        var g = new double[3 * _n];
        double e = Evaluate(Flatten(coordinates), g, null);
        for (int i = 0; i < _n; i++) gradient[i] = new Point3(g[3 * i], g[3 * i + 1], g[3 * i + 2]);
        return e;
    }

    /// <summary>Energy split into its terms.</summary>
    public UffEnergyTerms EnergyTerms(ReadOnlySpan<Point3> coordinates)
    {
        var parts = new double[5];
        Evaluate(Flatten(coordinates), [], parts);
        return new UffEnergyTerms(parts[0], parts[1], parts[2], parts[3], parts[4]);
    }

    /// <summary>Minimises the energy in place (L-BFGS); returns the final energy.</summary>
    public MinimizationResult Minimize(Point3[] coordinates, int maxIterations = 2000, double gradientTolerance = 1e-4)
    {
        if (coordinates.Length != _n) throw new ArgumentException("Coordinate count must equal the atom count.", nameof(coordinates));
        var x = Flatten(coordinates);
        var r = Lbfgs.Minimize(x, (xs, g) => { g.Clear(); return Evaluate(xs, g, null); }, maxIterations, gradientTolerance);
        for (int i = 0; i < _n; i++) coordinates[i] = new Point3(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
        return r;
    }

    /// <summary>
    /// Convenience: UFF-minimises the molecule's own 3-D coordinates in place. The molecule should carry explicit
    /// hydrogens; returns the minimisation outcome.
    /// </summary>
    /// <exception cref="InvalidOperationException">No coordinates, or an atom has no UFF type.</exception>
    public static MinimizationResult Optimize(Molecule mol, int maxIterations = 2000)
    {
        var coords = mol.Coordinates ?? throw new InvalidOperationException("The molecule has no coordinates.");
        var ff = Create(mol);
        var copy = (Point3[])coords.Clone();
        var r = ff.Minimize(copy, maxIterations);
        mol.Coordinates = copy;
        return r;
    }

    internal double Evaluate(ReadOnlySpan<double> x, Span<double> g) => Evaluate(x, g, null);

    private double Evaluate(ReadOnlySpan<double> x, Span<double> g, double[]? parts)
    {
        bool grad = !g.IsEmpty;
        double eb = 0, ea = 0, et = 0, ei = 0, ev = 0;
        foreach (var (i, j, k, r0) in _bonds)
        {
            double r = Math.Sqrt(FFMath.Distance2(x, i, j, 3));
            double dr = r - r0;
            eb += 0.5 * k * dr * dr;
            if (grad) FFMath.AddDistanceGradient(x, g, i, j, 3, r, k * dr);
        }
        foreach (var a in _angles)
        {
            double c = FFMath.CosAngle(x, a.I, a.J, a.K, 3);
            double e, dEdc;
            if (a.Order == 0)
            {
                e = a.Ka * (a.C0 + a.C1 * c + a.C2 * (2 * c * c - 1));
                dEdc = a.Ka * (a.C1 + 4 * a.C2 * c);
            }
            else
            {
                if (a.Order == 1)
                {
                    // linear centre: E = k (1 + cos θ), minimum at 180°
                    e = a.Ka * (1 + c);
                    dEdc = a.Ka;
                }
                else
                {
                    // E = k (1 − cos nθ) / n²
                    double cosN = a.Order switch { 2 => 2 * c * c - 1, 3 => c * (4 * c * c - 3), _ => 8 * c * c * c * c - 8 * c * c + 1 };
                    double sinRatio = a.Order switch { 2 => 2 * c, 3 => 4 * c * c - 1, _ => 4 * c * (2 * c * c - 1) }; // sin(nθ)/sin θ
                    e = a.Ka * (1 - cosN) / (a.Order * a.Order);
                    dEdc = -a.Ka / a.Order * sinRatio;
                }
                if (c > 0.8660)
                {
                    double theta = Math.Acos(c);
                    double pen = Math.Exp(-20.0 * (theta - a.Theta0 + 0.25));
                    e += pen;
                    double sinT = Math.Max(Math.Sqrt(1 - c * c), 1e-8);
                    dEdc += -20.0 * pen * (-1.0 / sinT);
                }
            }
            ea += e;
            if (grad) FFMath.AddCosAngleGradient(x, g, a.I, a.J, a.K, 3, dEdc);
        }
        foreach (var (i, j, k, l, v, order, cosTerm) in _torsions)
        {
            double phi = FFMath.Dihedral(x, i, j, k, l, 3);
            et += v / 2.0 * (1.0 - cosTerm * Math.Cos(order * phi));
            if (grad) FFMath.AddDihedralGradient(x, g, i, j, k, l, 3, v / 2.0 * cosTerm * order * Math.Sin(order * phi));
        }
        foreach (var (i, j, k, l, k0, c0, c1, c2) in _inversions)
        {
            double cosY = FFMath.InversionCosY(x, i, j, k, l, 3);
            double sinY2 = 1 - cosY * cosY;
            double sinY = sinY2 > 0 ? Math.Sqrt(sinY2) : 0;
            double cos2W = 2 * sinY * sinY - 1;
            ei += k0 * (c0 + c1 * sinY + c2 * cos2W);
            if (grad) FFMath.AddInversionGradient(x, g, i, j, k, l, 3, k0, c1, c2);
        }
        foreach (var (i, j, xij, dij, cut2) in _vdw)
        {
            double r2 = FFMath.Distance2(x, i, j, 3);
            if (r2 > cut2 || r2 <= 0) continue;
            double u2 = xij * xij / r2;
            double u6 = u2 * u2 * u2;
            double u12 = u6 * u6;
            ev += dij * (u12 - 2 * u6);
            if (grad)
            {
                double r = Math.Sqrt(r2);
                FFMath.AddDistanceGradient(x, g, i, j, 3, r, -12.0 * dij / r * (u12 - u6));
            }
        }
        if (parts is not null) { parts[0] = eb; parts[1] = ea; parts[2] = et; parts[3] = ei; parts[4] = ev; }
        return eb + ea + et + ei + ev;
    }

    // ---- parameters ---------------------------------------------------------------------------------

    internal static double BondOrderValue(Bond b) => b.Order switch
    {
        BondOrder.Aromatic => 1.5,
        BondOrder.Double => 2,
        BondOrder.Triple => 3,
        _ => 1,
    };

    internal static double RestLength(double bondOrder, in UffAtomParameters a, in UffAtomParameters b)
    {
        double ri = a.R1, rj = b.R1;
        double rBO = -Lambda * (ri + rj) * Math.Log(bondOrder);
        double xi = a.Xi, xj = b.Xi;
        double d = Math.Sqrt(xi) - Math.Sqrt(xj);
        double rEN = ri * rj * d * d / (xi * ri + xj * rj);
        return ri + rj + rBO - rEN;
    }

    private static double AngleForceConstant(double theta0, double bo12, double bo23, in UffAtomParameters p1, in UffAtomParameters p2, in UffAtomParameters p3)
    {
        double cosT = Math.Cos(theta0);
        double r12 = RestLength(bo12, p1, p2), r23 = RestLength(bo23, p2, p3);
        double r13 = Math.Sqrt(r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * cosT);
        double beta = 2.0 * G / (r12 * r23);
        double pre = beta * p1.Z1 * p3.Z1 / Math.Pow(r13, 5);
        double rTerm = r12 * r23;
        double inner = 3.0 * rTerm * (1.0 - cosT * cosT) - r13 * r13 * cosT;
        return pre * rTerm * inner;
    }

    private static bool InGroup6(int z) => z is 8 or 16 or 34 or 52 or 84;

    private static double Equation17(double bo, in UffAtomParameters a, in UffAtomParameters b)
        => 5.0 * Math.Sqrt(a.U1 * b.U1) * (1.0 + 4.18 * Math.Log(bo));

    private static (double V, int Order, double CosTerm) TorsionParams(double bo, int z2, int z3, Hybridization h2, Hybridization h3,
        in UffAtomParameters p2, in UffAtomParameters p3, bool endAtomIsSp2)
    {
        if (h2 == Hybridization.SP3 && h3 == Hybridization.SP3)
        {
            if (bo == 1.0 && InGroup6(z2) && InGroup6(z3))
            {
                double v2 = z2 == 8 ? 2.0 : 6.8, v3 = z3 == 8 ? 2.0 : 6.8;
                return (Math.Sqrt(v2 * v3), 2, -1);
            }
            return (Math.Sqrt(p2.V1 * p3.V1), 3, -1);
        }
        if (h2 == Hybridization.SP2 && h3 == Hybridization.SP2)
            return (Equation17(bo, p2, p3), 2, 1);
        // sp2–sp3
        if (bo == 1.0)
        {
            if ((h2 == Hybridization.SP3 && InGroup6(z2) && !InGroup6(z3)) || (h3 == Hybridization.SP3 && InGroup6(z3) && !InGroup6(z2)))
                return (Equation17(bo, p2, p3), 2, -1);
            if (endAtomIsSp2) return (2.0, 3, -1);
        }
        return (1.0, 6, 1);
    }

    private static (double K, double C0, double C1, double C2) InversionParams(int z, bool boundToSp2O)
    {
        double c0, c1, c2, k;
        if (z is 6 or 7 or 8)
        {
            c0 = 1; c1 = -1; c2 = 0;
            k = boundToSp2O ? 50.0 : 6.0;
        }
        else
        {
            double w0 = FFMath.Deg2Rad * z switch { 15 => 84.4339, 33 => 86.9735, 51 => 87.7047, _ => 90.0 };
            c2 = 1.0;
            c1 = -4.0 * Math.Cos(w0);
            c0 = -(c1 * Math.Cos(w0) + c2 * Math.Cos(2.0 * w0));
            k = 22.0 / (c0 + c1 + c2);
        }
        return (k / 3.0, c0, c1, c2);
    }

    private static bool HasTripleBond(Molecule mol, int atom)
    {
        foreach (int b in mol.NeighborBonds(atom)) if (mol.Bonds[b].Order == BondOrder.Triple) return true;
        return false;
    }

    private static long PairKey(int a, int b, int n) => a < b ? (long)a * n + b : (long)b * n + a;

    /// <summary>The UFF label of one atom (RDKit <c>UFF::Tools::getAtomLabel</c> rules).</summary>
    internal static string Label(Molecule mol, int atom, Hybridization hyb, bool hasConjugatedBond)
    {
        var a = mol.Atoms[atom];
        int z = a.AtomicNumber;
        string key = a.Symbol;
        if (key.Length == 1) key += "_";
        if (z <= 0) return key;
        int dv = AtomHybridization.DefaultValence(z);
        int nouter = PeriodicTable.Get(z).OuterElectrons;
        if (dv == -1 || (nouter != 1 && nouter != 7))
        {
            switch (z)
            {
                case 12 or 13 or 14 or 15 or 50 or 51 or 52 or 81 or 82 or 83 or 84: key += '3'; break;
                case 80: key += '1'; break;
                default:
                    key += hyb switch
                    {
                        Hybridization.SP => "1",
                        Hybridization.SP2 => (a.IsAromatic || hasConjugatedBond) && z is 6 or 7 or 8 or 16 ? "R" : "2",
                        Hybridization.SP3 => "3",
                        Hybridization.SP3D => "5",
                        Hybridization.SP3D2 => "6",
                        _ => "",
                    };
                    break;
            }
        }
        int tv = mol.TotalValence(atom);
        switch (z)
        {
            case 15: key += tv == 3 ? "+3" : "+5"; break;
            case 16:
                if (hyb != Hybridization.SP2) key += tv switch { 2 => "+2", 4 => "+4", _ => "+6" };
                break;
            case 34 or 52 or 84: key += "+2"; break;
            case 33 or 31 or 49 or 51 or 81 or 82 or 83: key += "+3"; break;
            case 29 or 47: key += "+1"; break;
            case 4 or 12 or 20 or 25 or 26 or 28 or 30 or 46 or 48 or 78 or 80: key += "+2"; break;
            case 21 or 24 or 27 or 79: key += "+3"; break;
            case 22: key += "+4"; break;
            case 23 or 41 or 43 or 73: key += "+5"; break;
            case 42: key += "+6"; break;
        }
        return key;
    }

    private static double[] Flatten(ReadOnlySpan<Point3> c)
    {
        var x = new double[3 * c.Length];
        for (int i = 0; i < c.Length; i++) { x[3 * i] = c[i].X; x[3 * i + 1] = c[i].Y; x[3 * i + 2] = c[i].Z; }
        return x;
    }
}
