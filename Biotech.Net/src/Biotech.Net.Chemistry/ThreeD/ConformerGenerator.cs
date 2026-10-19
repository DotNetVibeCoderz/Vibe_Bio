namespace Biotech.Net.Chemistry;

/// <summary>Options for <see cref="ConformerGenerator"/>. The defaults correspond to RDKit's ETKDGv3 settings.</summary>
public sealed record ConformerOptions
{
    /// <summary>Random seed; conformer <c>k</c> uses <c>(k + 1) × seed</c>, so results are reproducible.</summary>
    public int RandomSeed { get; init; } = 42;

    /// <summary>Embedding attempts per conformer before giving up (0 = 10 × atom count, as in RDKit).</summary>
    public int MaxAttempts { get; init; }

    /// <summary>Use the experimental torsion-angle preferences ("ET" in ETKDG).</summary>
    public bool UseExperimentalTorsions { get; init; } = true;

    /// <summary>Use "basic knowledge" terms: planar sp² centres, flat aromatic rings, linear sp centres ("K").</summary>
    public bool UseBasicKnowledge { get; init; } = true;

    /// <summary>Add the ETKDGv3 macrocycle torsion patterns (only affect rings of 9+ atoms).</summary>
    public bool UseMacrocycleTorsions { get; init; } = true;

    /// <summary>Keep chain amides and esters trans (RDKit <c>forceTransAmides</c>).</summary>
    public bool ForceTransAmides { get; init; } = true;

    /// <summary>Reject conformers whose tetrahedral stereocentres or E/Z bonds disagree with the input.</summary>
    public bool EnforceChirality { get; init; } = true;

    /// <summary>Discard a conformer whose symmetry-aware heavy-atom RMSD to an already accepted one is below this (Å); ≤ 0 disables.</summary>
    public double PruneRmsThreshold { get; init; } = 0.5;

    /// <summary>Finish each conformer with a UFF minimisation (<see cref="UffForceField"/>).</summary>
    public bool ForceFieldCleanup { get; init; } = true;

    /// <summary>Iteration cap for the UFF clean-up.</summary>
    public int CleanupMaxIterations { get; init; } = 500;

    /// <summary>Degree of parallelism across conformers (−1 = all cores). Results do not depend on it.</summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;

    /// <summary>ETKDGv3-like defaults.</summary>
    public static ConformerOptions Default { get; } = new();
}

/// <summary>A set of 3-D conformers of one molecule (with explicit hydrogens).</summary>
public sealed class ConformerEnsemble
{
    internal ConformerEnsemble(Molecule molecule, List<Point3[]> conformers, List<double> energies, int attempted)
    {
        Molecule = molecule;
        Conformers = conformers;
        Energies = energies;
        Attempted = attempted;
    }

    /// <summary>The molecule with explicit hydrogens; its atom order matches every conformer.</summary>
    public Molecule Molecule { get; }

    /// <summary>Conformer coordinates (Å), in generation order.</summary>
    public IReadOnlyList<Point3[]> Conformers { get; }

    /// <summary>UFF energy (kcal/mol) of each conformer, or NaN when no clean-up was done / UFF could not type the molecule.</summary>
    public IReadOnlyList<double> Energies { get; }

    /// <summary>Number of conformers requested (before failures and RMSD pruning).</summary>
    public int Attempted { get; }

    /// <summary>Number of conformers kept.</summary>
    public int Count => Conformers.Count;

    /// <summary>A copy of the molecule carrying conformer <paramref name="index"/> as its coordinates (energy in property <c>UFF_ENERGY</c>).</summary>
    public Molecule ToMolecule(int index)
    {
        var m = Molecule.Clone();
        m.Coordinates = (Point3[])Conformers[index].Clone();
        m.Name = string.IsNullOrEmpty(Molecule.Name) ? $"conf{index}" : $"{Molecule.Name}_conf{index}";
        if (!double.IsNaN(Energies[index])) m.Properties["UFF_ENERGY"] = Energies[index].ToString("F4", System.Globalization.CultureInfo.InvariantCulture);
        return m;
    }

    /// <summary>All conformers as separate molecules (e.g. for <see cref="SdfWriter"/>).</summary>
    public IEnumerable<Molecule> ToMolecules() => Enumerable.Range(0, Count).Select(ToMolecule);
}

/// <summary>
/// 3-D conformer generation by distance geometry with experimental-torsion knowledge — an ETKDG-style pipeline
/// (Riniker &amp; Landrum, J. Chem. Inf. Model. 2015, 55, 2562; Wang et al., J. Chem. Inf. Model. 2020, 60, 2044)
/// following RDKit's implementation:
/// <list type="number">
/// <item>bounds matrix from the graph (1-2 UFF bond lengths, 1-3 ideal angles, 1-4 torsion ranges, vdW lower bounds) and triangle smoothing;</item>
/// <item>random distances within the bounds, metric-matrix embedding in 4-D (3-D without stereocentres) from the top eigenvectors;</item>
/// <item>minimisation of the distance-violation + chiral-volume error function, then squeezing out the fourth dimension;</item>
/// <item>tetrahedral/chirality checks, minimisation with the ETK force field (experimental torsions, planarity, linear centres);</item>
/// <item>optional UFF clean-up, and symmetry-aware RMSD pruning of near-duplicates.</item>
/// </list>
/// </summary>
/// <remarks>
/// Hydrogens are added automatically. Not implemented compared with RDKit: 1-5 distance bounds, macrocycle-specific
/// 1-4 bounds, atropisomer and non-tetrahedral stereo, coordinate maps (constrained embedding) and the newer
/// all-in-one (AIO) minimiser. Torsion-pattern matching order can differ from RDKit's for symmetric matches.
/// </remarks>
public static class ConformerGenerator
{
    /// <summary>Generates one conformer and returns a copy of the molecule (with explicit hydrogens) carrying it.</summary>
    /// <exception cref="InvalidOperationException">Embedding failed within the attempt budget.</exception>
    public static Molecule Embed(Molecule mol, ConformerOptions? options = null)
    {
        options ??= ConformerOptions.Default;
        var ens = Generate(mol, 1, options with { PruneRmsThreshold = -1 });
        if (ens.Count == 0) throw new InvalidOperationException($"Could not embed '{mol.Name}' in 3-D.");
        return ens.ToMolecule(0);
    }

    /// <summary>
    /// Generates up to <paramref name="count"/> conformers. Conformers that fail to embed or fall within
    /// <see cref="ConformerOptions.PruneRmsThreshold"/> of an earlier one are dropped, so fewer may be returned.
    /// </summary>
    public static ConformerEnsemble Generate(Molecule mol, int count, ConformerOptions? options = null)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(count, 1);
        options ??= ConformerOptions.Default;
        var withH = mol.AddHydrogens();
        withH.Coordinates = null;
        int n = withH.AtomCount;
        if (n == 0) throw new ArgumentException("The molecule has no atoms.", nameof(mol));

        var frags = withH.Fragments();
        var contexts = frags.Select(f => EmbedContext.Create(withH, f, options)).ToList();
        UffForceField? uff = null;
        if (options.ForceFieldCleanup) UffForceField.TryCreate(withH, out uff, out _);

        var results = new Point3[]?[count];
        var energies = new double[count];
        var po = new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism };
        Parallel.For(0, count, po, ci =>
        {
            int seed = unchecked((ci + 1) * options.RandomSeed);
            var rng = new Random(seed);
            var coords = new Point3[n];
            double cursor = 0;
            bool first = true;
            foreach (var ctx in contexts)
            {
                var c = ctx.TryEmbed(rng, options);
                if (c is null) return;
                // centre each fragment and lay fragments side by side along x (4 Å apart)
                var cen = new Point3(c.Average(p => p.X), c.Average(p => p.Y), c.Average(p => p.Z));
                double minX = c.Min(p => p.X) - cen.X, maxX = c.Max(p => p.X) - cen.X;
                double shift = first ? 0 : cursor - minX + 4.0;
                for (int k = 0; k < c.Length; k++) coords[ctx.AtomMap[k]] = c[k] - cen + new Point3(shift, 0, 0);
                cursor = shift + maxX;
                first = false;
            }
            double e = double.NaN;
            if (uff is not null)
            {
                var r = uff.Minimize(coords, options.CleanupMaxIterations, 1e-3);
                e = r.Energy;
            }
            results[ci] = coords;
            energies[ci] = e;
        });

        var kept = new List<Point3[]>();
        var keptE = new List<double>();
        IReadOnlyList<int[]>? autos = null;
        for (int ci = 0; ci < count; ci++)
        {
            var c = results[ci];
            if (c is null) continue;
            if (options.PruneRmsThreshold > 0 && kept.Count > 0)
            {
                autos ??= ConformerAlignment.HeavyAtomAutomorphisms(withH);
                bool far = true;
                foreach (var k in kept)
                    if (MinRmsd(autos, k, c) < options.PruneRmsThreshold) { far = false; break; }
                if (!far) continue;
            }
            kept.Add(c);
            keptE.Add(energies[ci]);
        }
        withH.Coordinates = kept.Count > 0 ? (Point3[])kept[0].Clone() : null;
        return new ConformerEnsemble(withH, kept, keptE, count);
    }

    private static double MinRmsd(IReadOnlyList<int[]> maps, Point3[] a, Point3[] b)
    {
        double best = double.MaxValue;
        int m = maps[0].Length;
        var pa = new Point3[m];
        var pb = new Point3[m];
        foreach (var map in maps)
        {
            for (int k = 0; k < m; k++) { pa[k] = a[maps[0][k]]; pb[k] = b[map[k]]; }
            best = Math.Min(best, ConformerAlignment.Rmsd(pa, pb));
        }
        return best;
    }

    // ---- per-fragment embedding ---------------------------------------------------------------------

    private sealed class EmbedContext
    {
        private const double BasinThresh = 5.0;
        public required Molecule Mol { get; init; }
        public required int[] AtomMap { get; init; }
        public required BoundsMatrix Bounds { get; init; }
        public required List<ChiralSet> Chiral { get; init; }
        public required List<ChiralSet> Tetrahedral { get; init; }
        public required List<(int I, int J, int K, int L, double[] V, int[] S)> Torsions { get; init; }
        public required List<(int N0, int C, int N1, int N2, bool BoundToO)> Impropers { get; init; }
        public required List<(int I, int J, int K, bool Linear)> AngleList { get; init; }
        public required List<(int Nbr, int Atom, int Other)> DoubleBondEnds { get; init; }
        public required List<(int A, int B, int C, int D, int Sign)> StereoDoubleBonds { get; init; }

        public static EmbedContext Create(Molecule full, int[] fragment, ConformerOptions o)
        {
            Molecule mol;
            if (fragment.Length == full.AtomCount) mol = full;
            else
            {
                var keep = new HashSet<int>(fragment);
                mol = full.Clone();
                mol.RemoveAtoms(Enumerable.Range(0, full.AtomCount).Where(i => !keep.Contains(i)));
            }
            int n = mol.AtomCount;
            var conj = AtomHybridization.ConjugatedBonds(mol);
            var hyb = AtomHybridization.Compute(mol, conj);
            var builder = new DistanceBoundsBuilder(mol, hyb, conj, o.ForceTransAmides);
            var bounds = builder.Build();
            if (!bounds.TriangleSmooth())
            {
                // retry without forced amide geometry; then give up on smoothing and keep raw bounds
                builder = new DistanceBoundsBuilder(mol, hyb, conj, false);
                var b2 = builder.Build();
                if (b2.TriangleSmooth()) bounds = b2;
                else { bounds = new DistanceBoundsBuilder(mol, hyb, conj, false).Build(); }
            }
            var (chiral, tetra) = ChiralSets(mol);
            var ctx = new EmbedContext
            {
                Mol = mol,
                AtomMap = fragment,
                Bounds = bounds,
                Chiral = chiral,
                Tetrahedral = tetra,
                Torsions = [],
                Impropers = [],
                AngleList = [],
                DoubleBondEnds = [],
                StereoDoubleBonds = [],
            };
            ctx.CollectGeometry(hyb);
            if (o.UseExperimentalTorsions || o.UseBasicKnowledge) ctx.CollectTorsions(hyb, builder.Paths14, o);
            return ctx;
        }

        private static (List<ChiralSet>, List<ChiralSet>) ChiralSets(Molecule mol)
        {
            var chiral = new List<ChiralSet>();
            var tetra = new List<ChiralSet>();
            var rings = mol.Rings;
            for (int a = 0; a < mol.AtomCount; a++)
            {
                var at = mol.Atoms[a];
                if (at.AtomicNumber == 1) continue;
                bool tagged = at.Chirality != ChiralTag.None;
                if (!tagged && !(at.AtomicNumber is 6 or 7 && mol.Degree(a) == 4)) continue;
                List<int> nbrs;
                if (tagged)
                {
                    var r = mol.GetChiralReference(a);
                    nbrs = r is not null ? r.Select(x => x < 0 ? a : x).ToList() : [.. mol.Neighbors(a).ToArray()];
                }
                else nbrs = [.. mol.Neighbors(a).ToArray()];
                if (nbrs.Count < 3) continue;
                double lower = 5.0, upper = 100.0;
                if (nbrs.Count < 4) { lower = 2.0; nbrs.Add(a); }
                if (nbrs.Count != 4) continue;
                int small = rings.AtomRings.Count(rr => rr.Length < 5 && Array.IndexOf(rr, a) >= 0);
                bool fused = small > 1;
                if (tagged)
                {
                    if (at.Chirality == ChiralTag.CounterClockwise) chiral.Add(new ChiralSet(a, nbrs[0], nbrs[1], nbrs[2], nbrs[3], lower, upper, fused));
                    else chiral.Add(new ChiralSet(a, nbrs[0], nbrs[1], nbrs[2], nbrs[3], -upper, -lower, fused));
                }
                else if (rings.AtomRingCount(a) >= 2 && !rings.IsAtomInRingOfSize(a, 3))
                    tetra.Add(new ChiralSet(a, nbrs[0], nbrs[1], nbrs[2], nbrs[3], 0, 0, fused));
            }
            return (chiral, tetra);
        }

        private void CollectGeometry(Hybridization[] hyb)
        {
            var bonds = Mol.Bonds;
            for (int b1 = 0; b1 < bonds.Length; b1++)
                for (int b2 = b1 + 1; b2 < bonds.Length; b2++)
                {
                    var x = bonds[b1];
                    var y = bonds[b2];
                    int center = x.Begin == y.Begin || x.Begin == y.End ? x.Begin : x.End == y.Begin || x.End == y.End ? x.End : -1;
                    if (center < 0) continue;
                    int i = x.Other(center), k = y.Other(center);
                    bool linear = x.Order == BondOrder.Triple || y.Order == BondOrder.Triple
                        || (x.Order == BondOrder.Double && y.Order == BondOrder.Double && Mol.Degree(center) == 2);
                    AngleList.Add((i, center, k, linear));
                }
            for (int a = 0; a < Mol.AtomCount; a++)
            {
                int z = Mol.Atoms[a].AtomicNumber;
                if (z is not (6 or 7 or 8) || hyb[a] != Hybridization.SP2 || Mol.Degree(a) != 3) continue;
                var nb = Mol.Neighbors(a);
                bool toO = false;
                if (z == 6) foreach (int x in nb) if (Mol.Atoms[x].AtomicNumber == 8 && hyb[x] == Hybridization.SP2) toO = true;
                Impropers.Add((nb[0], a, nb[1], nb[2], toO));
            }
            for (int b = 0; b < bonds.Length; b++)
            {
                var bd = bonds[b];
                if (bd.Order != BondOrder.Double) continue;
                foreach (int atom in new[] { bd.Begin, bd.End })
                {
                    if (Mol.Degree(atom) < 2) continue;
                    int other = bd.Other(atom);
                    var nb = Mol.Neighbors(atom);
                    var nbb = Mol.NeighborBonds(atom);
                    for (int k = 0; k < nb.Length; k++)
                    {
                        if (nb[k] == other) continue;
                        if (bonds[nbb[k]].Order != BondOrder.Single && Mol.Degree(atom) == 2) continue;
                        DoubleBondEnds.Add((nb[k], atom, other));
                    }
                }
                if (bd.Stereo is BondStereo.Cis or BondStereo.Trans && bd.StereoAtomBegin >= 0 && bd.StereoAtomEnd >= 0)
                    StereoDoubleBonds.Add((bd.StereoAtomBegin, bd.Begin, bd.End, bd.StereoAtomEnd, bd.Stereo == BondStereo.Cis ? -1 : 1));
            }
        }

        private void CollectTorsions(Hybridization[] hyb, List<Path14> paths, ConformerOptions o)
        {
            var bonds = Mol.Bonds;
            int nb = bonds.Length;
            var rings = Mol.Rings;
            var excluded = new bool[nb];
            var bondRings = rings.BondRings;
            var bondRingCount = new int[nb];
            foreach (var r in bondRings) foreach (int b in r) bondRingCount[b]++;
            for (int x = 0; x < bondRings.Count; x++)
                for (int y = x + 1; y < bondRings.Count; y++)
                {
                    if (bondRings[x].Length >= 9 && bondRings[y].Length >= 9) continue;
                    if (bondRings[x].Intersect(bondRings[y]).Count() > 1)
                    {
                        if (bondRings[x].Length < 9) foreach (int b in bondRings[x]) excluded[b] = true;
                        if (bondRings[y].Length < 9) foreach (int b in bondRings[y]) excluded[b] = true;
                    }
                }
            var done = new bool[nb];
            bool Skip(int bond)
            {
                if (excluded[bond] || bondRingCount[bond] > 3) done[bond] = true;
                return done[bond];
            }
            if (o.UseBasicKnowledge)
            {
                foreach (var p in paths)
                {
                    if (!p.Forced || p.Class is not (TorsionClass.Cis or TorsionClass.Trans)) continue;
                    int bond = Mol.BondBetween(p.A2, p.A3);
                    if (bond < 0 || Skip(bond)) continue;
                    var v = new double[6];
                    var s = new int[] { 1, 1, 1, 1, 1, 1 };
                    v[0] = 75.0;
                    if (p.Class == TorsionClass.Cis) s[0] = -1;
                    Torsions.Add((p.A1, p.A2, p.A3, p.A4, v, s));
                    done[bond] = true;
                }
            }
            if (o.UseExperimentalTorsions)
            {
                var library = o.UseMacrocycleTorsions ? EtkdgTorsionLibrary.Regular.Concat(EtkdgTorsionLibrary.Macrocycle) : EtkdgTorsionLibrary.Regular;
                foreach (var t in library)
                {
                    foreach (var match in t.Pattern.Matches(Mol, unique: false))
                    {
                        int a1 = match[t.MappedAtoms[0]], a2 = match[t.MappedAtoms[1]], a3 = match[t.MappedAtoms[2]], a4 = match[t.MappedAtoms[3]];
                        int bond = Mol.BondBetween(a2, a3);
                        if (bond < 0 || Skip(bond)) continue;
                        Torsions.Add((a1, a2, a3, a4, t.V, t.Signs));
                        done[bond] = true;
                    }
                }
            }
            if (o.UseBasicKnowledge)
            {
                foreach (var ring in rings.AtomRings)
                {
                    int r = ring.Length;
                    if (r < 4 || r > 6) continue;
                    for (int i = 0; i < r; i++)
                    {
                        int a1 = ring[i], a2 = ring[(i + 1) % r], a3 = ring[(i + 2) % r], a4 = ring[(i + 3) % r];
                        int bond = Mol.BondBetween(a2, a3);
                        if (done[bond]) continue;
                        if (hyb[a1] != Hybridization.SP2 || hyb[a2] != Hybridization.SP2 || hyb[a3] != Hybridization.SP2 || hyb[a4] != Hybridization.SP2) continue;
                        done[bond] = true;
                        Torsions.Add((a1, a2, a3, a4, [0, 100.0, 0, 0, 0, 0], [1, -1, 1, 1, 1, 1]));
                    }
                }
            }
        }

        public Point3[]? TryEmbed(Random rng, ConformerOptions o)
        {
            int n = Mol.AtomCount;
            if (n == 1) return [new Point3(0, 0, 0)];
            int maxAttempts = o.MaxAttempts > 0 ? o.MaxAttempts : 10 * n;
            bool fourD = Chiral.Count > 0;
            int dim = fourD ? 4 : 3;
            var dist = new double[n * n];
            var x = new double[n * dim];
            for (int attempt = 0; attempt < maxAttempts; attempt++)
            {
                for (int i = 1; i < n; i++)
                    for (int j = 0; j < i; j++)
                    {
                        double ub = Bounds.Upper(i, j), lb = Bounds.Lower(i, j);
                        double d = lb + rng.NextDouble() * (ub - lb);
                        dist[i * n + j] = dist[j * n + i] = d;
                    }
                if (!MetricMatrixEmbedding.Embed(dist, n, dim, rng, x)) continue;
                // first minimisation: distance + chiral violations, small 4th-dimension weight
                var dg1 = BuildDg(dim, 1.0, 0.1);
                Lbfgs.Minimize(x, dg1.Evaluate, 4000, 1e-3);
                double e1 = dg1.Evaluate(x, []);
                if (e1 / n >= 0.05) continue;
                if (!CheckTetrahedral(x, dim) || (o.EnforceChirality && !CheckChiral(x, dim))) continue;
                if (fourD)
                {
                    var dg2 = BuildDg(dim, 0.2, 1.0);
                    Lbfgs.Minimize(x, dg2.Evaluate, 4000, 1e-3);
                }
                var x3 = new double[3 * n];
                for (int i = 0; i < n; i++) for (int k = 0; k < 3; k++) x3[3 * i + k] = x[dim * i + k];
                if (o.UseExperimentalTorsions || o.UseBasicKnowledge)
                {
                    var etk = BuildEtk(x3, o.UseBasicKnowledge);
                    Lbfgs.Minimize(x3, etk.Evaluate, 300, 1e-3);
                    if (o.UseBasicKnowledge && !PlanarityOk(x3)) continue;
                }
                if (!DoubleBondGeometryOk(x3)) continue;
                if (o.EnforceChirality)
                {
                    if (Chiral.Count > 0 && !(CheckChiral(x3, 3) && ChiralBoundsOk(x3) && CentersInVolume(x3))) continue;
                    if (!DoubleBondStereoOk(x3)) continue;
                }
                var result = new Point3[n];
                for (int i = 0; i < n; i++) result[i] = new Point3(x3[3 * i], x3[3 * i + 1], x3[3 * i + 2]);
                return result;
            }
            return null;
        }

        private DgForceField BuildDg(int dim, double chiralWeight, double fourthDimWeight)
        {
            int n = Mol.AtomCount;
            var terms = new List<(int, int, double, double, double)>();
            for (int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    double l = Bounds.Lower(i, j), u = Bounds.Upper(i, j);
                    if (u - l <= BasinThresh) terms.Add((i, j, l * l, u * u, 1.0));
                }
            var chiral = new List<ChiralSet>(Chiral);
            return new DgForceField(dim, terms, chiral, chiralWeight, fourthDimWeight);
        }

        private EtkForceField BuildEtk(double[] x, bool basicKnowledge)
        {
            int n = Mol.AtomCount;
            var ff = new EtkForceField();
            var pairs = new bool[n * n];
            void Mark(int a, int b) { pairs[Math.Min(a, b) * n + Math.Max(a, b)] = true; }
            var doneBond = new HashSet<long>();
            foreach (var t in Torsions)
            {
                long key = (long)Math.Min(t.J, t.K) * n + Math.Max(t.J, t.K);
                if (!doneBond.Add(key)) continue;
                Mark(t.I, t.L);
                ff.Torsions.Add(t);
            }
            var improperCentre = new bool[n];
            if (basicKnowledge)
            {
                foreach (var (n0, c, n1, n2, toO) in Impropers)
                {
                    double kf = 10.0 * (toO ? 50.0 : 6.0) / 3.0;
                    ff.Impropers.Add((n0, c, n1, n2, kf));
                    ff.Impropers.Add((n0, c, n2, n1, kf));
                    ff.Impropers.Add((n1, c, n2, n0, kf));
                    improperCentre[c] = true;
                }
            }
            foreach (var b in Mol.Bonds)
            {
                Mark(b.Begin, b.End);
                double d = Math.Sqrt(FFMath.Distance2(x, b.Begin, b.End, 3));
                ff.Distances.Add((b.Begin, b.End, d - 0.01, d + 0.01, 100.0));
            }
            foreach (var (i, j, k, linear) in AngleList)
            {
                Mark(i, k);
                if (basicKnowledge && linear) ff.Angles.Add((i, j, k, 179.0, 180.0, 1.0));
                else if (improperCentre[j]) ff.Distances.Add((i, k, Bounds.Lower(i, k), Bounds.Upper(i, k), 100.0));
                else
                {
                    double d = Math.Sqrt(FFMath.Distance2(x, i, k, 3));
                    ff.Distances.Add((i, k, d - 0.01, d + 0.01, 100.0));
                }
            }
            for (int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                    if (!pairs[j * n + i]) ff.Distances.Add((i, j, Bounds.Lower(i, j), Bounds.Upper(i, j), 10.0));
            return ff;
        }

        private bool PlanarityOk(double[] x)
        {
            var ff = new EtkForceField();
            foreach (var (n0, c, n1, n2, toO) in Impropers)
            {
                double kf = 10.0 * (toO ? 50.0 : 6.0) / 3.0;
                ff.Impropers.Add((n0, c, n1, n2, kf));
                ff.Impropers.Add((n0, c, n2, n1, kf));
                ff.Impropers.Add((n1, c, n2, n0, kf));
            }
            foreach (var (i, j, k, linear) in AngleList) if (linear) ff.Angles.Add((i, j, k, 179.0, 180.0, 10.0));
            return ff.Evaluate(x, []) <= Impropers.Count * 0.7;
        }

        private bool CheckTetrahedral(double[] x, int dim)
        {
            foreach (var t in Tetrahedral)
                if (!VolumeTest(x, dim, t) || !CenterInVolume(x, dim, t, 0.30)) return false;
            return true;
        }

        private bool CheckChiral(double[] x, int dim)
        {
            foreach (var c in Chiral)
            {
                double vol = FFMath.ChiralVolume(x, c.N1, c.N2, c.N3, c.N4, dim);
                double lb = c.Lower, ub = c.Upper;
                if (lb > 0 && vol < lb && (vol / lb < 0.8 || Math.Sign(vol) != Math.Sign(lb))) return false;
                if (ub < 0 && vol > ub && (vol / ub < 0.8 || Math.Sign(vol) != Math.Sign(ub))) return false;
            }
            return true;
        }

        private bool ChiralBoundsOk(double[] x)
        {
            var atoms = new SortedSet<int>();
            foreach (var c in Chiral)
            {
                if (c.Center == c.N4) continue;
                atoms.Add(c.Center); atoms.Add(c.N1); atoms.Add(c.N2); atoms.Add(c.N3); atoms.Add(c.N4);
            }
            var list = atoms.ToArray();
            for (int p = 0; p < list.Length - 1; p++)
                for (int q = p + 1; q < list.Length; q++)
                {
                    int a = list[p], b = list[q];
                    double d = Math.Sqrt(FFMath.Distance2(x, a, b, 3));
                    double lb = Bounds.Lower(a, b), ub = Bounds.Upper(a, b);
                    if ((d < lb && Math.Abs(d - lb) > 0.1 * ub) || (d > ub && Math.Abs(d - ub) > 0.1 * ub)) return false;
                }
            return true;
        }

        private bool CentersInVolume(double[] x)
        {
            foreach (var c in Chiral)
                if (c.Center != c.N4 && !CenterInVolume(x, 3, c, 0.1)) return false;
            return true;
        }

        private static Point3 P(double[] x, int dim, int i) => new(x[dim * i], x[dim * i + 1], x[dim * i + 2]);

        private static bool VolumeTest(double[] x, int dim, ChiralSet c)
        {
            var p0 = P(x, dim, c.Center);
            Point3 U(int i) { var v = p0 - P(x, dim, i); double l = v.Length; return l > 0 ? v * (1 / l) : v; }
            var v1 = U(c.N1); var v2 = U(c.N2); var v3 = U(c.N3); var v4 = U(c.N4);
            double scale = c.FusedSmallRings ? 0.25 : 1.0;
            double min = scale * 0.5;
            return Math.Abs(Point3.Dot(Point3.Cross(v1, v2), v3)) >= min && Math.Abs(Point3.Dot(Point3.Cross(v1, v2), v4)) >= min
                && Math.Abs(Point3.Dot(Point3.Cross(v1, v3), v4)) >= min && Math.Abs(Point3.Dot(Point3.Cross(v2, v3), v4)) >= min;
        }

        private static bool SameSide(Point3 v1, Point3 v2, Point3 v3, Point3 v4, Point3 p0, double tol)
        {
            var normal = Point3.Cross(v2 - v1, v3 - v1);
            double d1 = Point3.Dot(normal, v4 - v1), d2 = Point3.Dot(normal, p0 - v1);
            if (Math.Abs(d1) < tol || Math.Abs(d2) < tol) return false;
            return (d1 < 0) == (d2 < 0);
        }

        private static bool CenterInVolume(double[] x, int dim, ChiralSet c, double tol)
        {
            if (c.Center == c.N4) return true;
            var p0 = P(x, dim, c.Center);
            Point3 p1 = P(x, dim, c.N1), p2 = P(x, dim, c.N2), p3 = P(x, dim, c.N3), p4 = P(x, dim, c.N4);
            return SameSide(p1, p2, p3, p4, p0, tol) && SameSide(p2, p3, p4, p1, p0, tol) && SameSide(p3, p4, p1, p2, p0, tol) && SameSide(p4, p1, p2, p3, p0, tol);
        }

        private bool DoubleBondGeometryOk(double[] x)
        {
            foreach (var (nbr, atom, other) in DoubleBondEnds)
            {
                var p0 = P(x, 3, nbr); var p1 = P(x, 3, atom); var p2 = P(x, 3, other);
                var v1 = p1 - p0; var v2 = p1 - p2;
                double l1 = v1.Length, l2 = v2.Length;
                if (l1 < 1e-8 || l2 < 1e-8) return false;
                if (Point3.Dot(v1, v2) / (l1 * l2) + 1.0 < 1e-3) return false;
            }
            return true;
        }

        private bool DoubleBondStereoOk(double[] x)
        {
            foreach (var (a, b, c, d, sign) in StereoDoubleBonds)
            {
                double phi = Math.Abs(FFMath.Dihedral(x, a, b, c, d, 3));
                if ((phi - Math.PI / 2) * sign < 0) return false;
            }
            return true;
        }
    }
}
