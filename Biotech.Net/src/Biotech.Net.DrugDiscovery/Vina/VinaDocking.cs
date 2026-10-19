using System.Diagnostics;
using System.Globalization;
using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;

namespace Biotech.Net.DrugDiscovery;

/// <summary>Options for <see cref="VinaDocking.Dock"/>. Defaults follow AutoDock Vina 1.2.</summary>
public sealed record DockingOptions
{
    /// <summary>Box centre (Å). Null = centre of the input ligand coordinates.</summary>
    public Vec3? Center { get; init; }

    /// <summary>Box edge lengths (Å).</summary>
    public Vec3 Size { get; init; } = new(22.5, 22.5, 22.5);

    /// <summary>Grid spacing for the precomputed maps (Å).</summary>
    public double GridSpacing { get; init; } = 0.375;

    /// <summary>Number of independent Monte-Carlo runs (Vina's "exhaustiveness").</summary>
    public int Exhaustiveness { get; init; } = 8;

    /// <summary>Monte-Carlo steps per run; 0 = Vina's heuristic 70·3·(50 + atoms + 10·DOF)/2.</summary>
    public int MonteCarloSteps { get; init; }

    /// <summary>Maximum number of poses returned.</summary>
    public int NumPoses { get; init; } = 9;

    /// <summary>Only poses within this many kcal/mol of the best are returned.</summary>
    public double EnergyRange { get; init; } = 3.0;

    /// <summary>Poses closer than this heavy-atom RMSD (Å) are clustered together, keeping the better one.</summary>
    public double MinRmsd { get; init; } = 1.0;

    /// <summary>Random seed; run k uses a seed derived from it, so results are reproducible.</summary>
    public int Seed { get; init; } = 42;

    /// <summary>Scoring-function weights and cut-off.</summary>
    public VinaScoreOptions Weights { get; init; } = VinaScoreOptions.Default;

    /// <summary>Degree of parallelism across runs (−1 = all cores). Results do not depend on it.</summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;

    /// <summary>Default options.</summary>
    public static DockingOptions Default { get; } = new();
}

/// <summary>Vina score of one pose.</summary>
/// <param name="Affinity">Estimated binding free energy (kcal/mol): (inter + intra − intra_ref) / (1 + w·N<sub>rot</sub>).</param>
/// <param name="Intermolecular">Ligand–receptor energy.</param>
/// <param name="Intramolecular">Ligand internal energy (pairs in different rigid fragments, more than three bonds apart).</param>
/// <param name="Terms">Unweighted intermolecular term sums.</param>
/// <param name="RotatableBonds">Vina's torsion count N<sub>rot</sub> used in the conformation-independent term.</param>
public sealed record VinaPoseScore(double Affinity, double Intermolecular, double Intramolecular, VinaTermSums Terms, double RotatableBonds)
{
    /// <summary>Vina's "torsional free energy" line: affinity − intermolecular − intramolecular + intra_ref.</summary>
    public double Torsional => Affinity - Intermolecular;

    /// <inheritdoc/>
    public override string ToString() => string.Create(CultureInfo.InvariantCulture, $"{Affinity:F3} kcal/mol (inter {Intermolecular:F3}, intra {Intramolecular:F3})");
}

/// <summary>One docked pose.</summary>
public sealed class DockedPose
{
    internal DockedPose(int rank, Molecule molecule, VinaPoseScore score, double rmsdLower, double rmsdUpper)
    {
        Rank = rank;
        Molecule = molecule;
        Score = score;
        RmsdLowerBound = rmsdLower;
        RmsdUpperBound = rmsdUpper;
    }

    /// <summary>1-based rank (1 = best).</summary>
    public int Rank { get; }

    /// <summary>The ligand with the pose's coordinates (same atoms as the input) and Vina properties for SDF export.</summary>
    public Molecule Molecule { get; }

    /// <summary>The pose score.</summary>
    public VinaPoseScore Score { get; }

    /// <summary>Affinity shortcut (kcal/mol).</summary>
    public double Affinity => Score.Affinity;

    /// <summary>Vina's RMSD lower bound to the best pose (nearest same-type atom matching).</summary>
    public double RmsdLowerBound { get; }

    /// <summary>Heavy-atom RMSD to the best pose with atoms matched by index (Vina's upper bound).</summary>
    public double RmsdUpperBound { get; }

    /// <inheritdoc/>
    public override string ToString() => string.Create(CultureInfo.InvariantCulture, $"#{Rank} {Affinity:F2} kcal/mol (rmsd l.b. {RmsdLowerBound:F2}, u.b. {RmsdUpperBound:F2})");
}

/// <summary>Result of a docking run.</summary>
public sealed class DockingResult
{
    internal DockingResult(IReadOnlyList<DockedPose> poses, VinaGrid grid, int torsions, double rotors, int evaluations, TimeSpan gridTime, TimeSpan searchTime)
    {
        Poses = poses;
        Grid = grid;
        TorsionCount = torsions;
        RotatableBonds = rotors;
        Evaluations = evaluations;
        GridTime = gridTime;
        SearchTime = searchTime;
    }

    /// <summary>Poses, best first.</summary>
    public IReadOnlyList<DockedPose> Poses { get; }

    /// <summary>The grid maps used (box, spacing).</summary>
    public VinaGrid Grid { get; }

    /// <summary>Torsional degrees of freedom searched.</summary>
    public int TorsionCount { get; }

    /// <summary>Vina's N<sub>rot</sub> for the torsional penalty.</summary>
    public double RotatableBonds { get; }

    /// <summary>Total energy evaluations.</summary>
    public int Evaluations { get; }

    /// <summary>Time spent building grid maps.</summary>
    public TimeSpan GridTime { get; }

    /// <summary>Time spent in the Monte-Carlo search.</summary>
    public TimeSpan SearchTime { get; }

    /// <summary>The best pose.</summary>
    public DockedPose Best => Poses[0];

    /// <summary>All poses as an SDF string (properties VINA_AFFINITY, VINA_INTER, VINA_INTRA, RMSD_LB, RMSD_UB).</summary>
    public string ToSdf() => SdfWriter.ToSdf(Poses.Select(p => p.Molecule));

    /// <summary>All poses as a multi-model PDBQT string (Vina output layout).</summary>
    public string ToPdbqt() => PdbqtWriter.Write(this);

    internal FlexibleLigand? Ligand { get; init; }
}

/// <summary>
/// Flexible-ligand docking into a rigid receptor with the AutoDock Vina scoring function and search
/// (Trott &amp; Olson, J. Comput. Chem. 2010, 31, 455; Eberhardt et al., J. Chem. Inf. Model. 2021, 61, 3891):
/// X-Score atom typing, precomputed grid maps with trilinear interpolation, a torsion-tree ligand model, and
/// Monte-Carlo search where every step mutates one degree of freedom (position, orientation or a torsion) and is
/// followed by a BFGS local optimisation in conformation space, with Metropolis acceptance at T = 1.2.
/// Poses are clustered by RMSD and scored as Vina does (intramolecular energy relative to the best pose,
/// divided by 1 + 0.0585·N<sub>rot</sub>).
/// </summary>
/// <remarks>
/// Receptor side chains are rigid. The pair potentials are evaluated analytically for intramolecular pairs and from
/// 0.005 Å lookup tables while building the maps (Vina uses r²-sampled tables), so energies agree with Vina to about
/// 0.01–0.1 kcal/mol rather than exactly. No hydration, macrocycle closure or AD4/Vinardo scoring.
/// </remarks>
public static class VinaDocking
{
    /// <summary>
    /// Scores a ligand pose (3-D coordinates required) against the receptor exactly (no grid), like
    /// <c>vina --score_only</c>.
    /// </summary>
    public static VinaPoseScore Score(Molecule ligand, VinaReceptor receptor, VinaScoreOptions? weights = null)
    {
        weights ??= VinaScoreOptions.Default;
        var lig = new FlexibleLigand(ligand);
        var coords = ligand.Coordinates!.Select(FlexibleLigand.ToVec).ToArray();
        var search = new NeighborSearch(receptor.Positions.ToArray(), weights.CutOff);
        var terms = default(VinaTermSums);
        var found = new List<int>();
        foreach (int i in lig.Scored)
        {
            found.Clear();
            search.Query(coords[i], weights.CutOff, found);
            foreach (int j in found)
            {
                double r = Vec3.Distance(coords[i], receptor.Positions[j]);
                if (r < weights.CutOff) terms += VinaPotential.Terms(lig.Types[i]!.Value, receptor.Types[j], r);
            }
        }
        double inter = terms.Weighted(weights);
        double intra = Intra(lig, coords, weights, 1000, null);
        return new VinaPoseScore(inter / (1 + weights.RotatableWeight * lig.NumRotors), inter, intra, terms, lig.NumRotors);
    }

    /// <summary>Docks a ligand (with 3-D coordinates; embedded with <see cref="ConformerGenerator"/> when missing).</summary>
    public static DockingResult Dock(Molecule ligand, VinaReceptor receptor, DockingOptions? options = null)
    {
        options ??= DockingOptions.Default;
        if (ligand.Coordinates is null || !ligand.Has3DCoordinates) ligand = ConformerGenerator.Embed(ligand);
        var lig = new FlexibleLigand(ligand);
        if (lig.Scored.Length == 0) throw new ArgumentException("The ligand has no atoms Vina can score.", nameof(ligand));
        var center = options.Center ?? lig.ReferenceOrigin;
        var sw = Stopwatch.StartNew();
        var grid = VinaGrid.Build(receptor, center, options.Size, lig.Scored.Select(i => lig.Types[i]!.Value), options.GridSpacing, options.Weights);
        var gridTime = sw.Elapsed;
        sw.Restart();

        var (steps, localSteps) = SearchSteps(ligand, lig, options);
        var runs = new List<(LigandConf Conf, double E, Vec3[] Heavy)>[options.Exhaustiveness];
        int evaluations = 0;
        var po = new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism };
        Parallel.For(0, options.Exhaustiveness, po, run =>
        {
            var engine = new Engine(lig, grid, options.Weights);
            var rng = new Random(unchecked(options.Seed * 7919 + run * 104729 + 17));
            runs[run] = engine.MonteCarlo(rng, steps, localSteps, options.MinRmsd, Math.Max(options.NumPoses, 9));
            Interlocked.Add(ref evaluations, engine.Evaluations);
        });
        return Finish(ligand, lig, grid, receptor, options, runs.SelectMany(r => r), localSteps, evaluations, gridTime, sw);
    }

    /// <summary>Monte-Carlo steps per run and BFGS steps per local optimisation (Vina's heuristics, or the options).</summary>
    internal static (int Steps, int LocalSteps) SearchSteps(Molecule ligand, FlexibleLigand lig, DockingOptions options)
    {
        int movable = lig.Scored.Length + Enumerable.Range(0, lig.AtomCount).Count(i => ligand.Atoms[i].AtomicNumber == 1 && ligand.Neighbors(i).ToArray().Any(nb => ligand.Atoms[nb].AtomicNumber is 7 or 8));
        int dof = 6 + lig.Torsions.Length;
        int steps = options.MonteCarloSteps > 0 ? options.MonteCarloSteps : 70 * 3 * (50 + movable + 10 * dof) / 2;
        return (steps, (25 + movable) / 3);
    }

    /// <summary>
    /// Search post-processing shared by the CPU and GPU searches: merges candidate poses by RMSD, refines each with
    /// exact receptor interactions, rescores and builds the result. <paramref name="sw"/> times the search phase.
    /// </summary>
    internal static DockingResult Finish(Molecule ligand, FlexibleLigand lig, VinaGrid grid, VinaReceptor receptor, DockingOptions options,
        IEnumerable<(LigandConf Conf, double E, Vec3[] Heavy)> candidates, int localSteps, int evaluations, TimeSpan gridTime, Stopwatch sw)
    {
        var po = new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism };
        var all = candidates.OrderBy(r => r.E).ToList();
        var merged = new List<(LigandConf Conf, double E, Vec3[] Heavy)>();
        foreach (var c in all)
            if (merged.All(m => Rmsd(m.Heavy, c.Heavy) >= options.MinRmsd)) merged.Add(c);

        // Vina 1.2 post-processing: refine every pose with exact (non-grid) receptor interactions, then rescore
        var exact = new ExactReceptor(receptor, grid, options.Weights);
        var refinedPoses = new (LigandConf Conf, Vec3[] X, double Inter, double Intra, Vec3[] Heavy)[merged.Count];
        Parallel.For(0, merged.Count, po, k =>
        {
            var engine = new Engine(lig, grid, options.Weights) { Exact = exact };
            var conf = merged[k].Conf.Clone();
            engine.Bfgs(conf, localSteps, 1000);
            Interlocked.Add(ref evaluations, engine.Evaluations);
            var x = new Vec3[lig.AtomCount];
            lig.Build(conf, x);
            refinedPoses[k] = (conf, x, exact.Energy(lig, x, null), Intra(lig, x, options.Weights, 1000, null), lig.Scored.Select(i => x[i]).ToArray());
        });
        var scored = refinedPoses.OrderBy(s => s.Inter + s.Intra).ToList();
        double intraRef = scored.Count > 0 ? scored[0].Intra : 0;
        double denom = 1 + options.Weights.RotatableWeight * lig.NumRotors;
        var poses = new List<DockedPose>();
        double bestAff = double.NaN;
        foreach (var s in scored)
        {
            double aff = (s.Inter + s.Intra - intraRef) / denom;
            if (poses.Count == 0) bestAff = aff;
            else if (aff > bestAff + options.EnergyRange || poses.Count >= options.NumPoses) break;
            var heavy0 = scored[0].Heavy;
            double ub = Rmsd(heavy0, s.Heavy), lb = RmsdLowerBound(lig, scored[0].X, s.X);
            var mol = ligand.Clone();
            mol.Coordinates = s.X.Select(v => new Point3(v.X, v.Y, v.Z)).ToArray();
            mol.Name = $"{(string.IsNullOrEmpty(ligand.Name) ? "ligand" : ligand.Name)}_pose{poses.Count + 1}";
            mol.Properties["VINA_AFFINITY"] = aff.ToString("F3", CultureInfo.InvariantCulture);
            mol.Properties["VINA_INTER"] = s.Inter.ToString("F3", CultureInfo.InvariantCulture);
            mol.Properties["VINA_INTRA"] = s.Intra.ToString("F3", CultureInfo.InvariantCulture);
            mol.Properties["RMSD_LB"] = lb.ToString("F3", CultureInfo.InvariantCulture);
            mol.Properties["RMSD_UB"] = ub.ToString("F3", CultureInfo.InvariantCulture);
            var score = new VinaPoseScore(aff, s.Inter, s.Intra, default, lig.NumRotors);
            poses.Add(new DockedPose(poses.Count + 1, mol, score, lb, ub));
        }
        return new DockingResult(poses, grid, lig.Torsions.Length, lig.NumRotors, evaluations, gridTime, sw.Elapsed) { Ligand = lig };
    }

    internal static double Intra(FlexibleLigand lig, Vec3[] x, VinaScoreOptions w, double cap, Vec3[]? grad)
    {
        double e = 0;
        foreach (var (i, j) in lig.IntraPairs)
        {
            var d = x[i] - x[j];
            double r = d.Length;
            double pe = VinaPotential.Energy(lig.Types[i]!.Value, lig.Types[j]!.Value, r, w, out double dEdr);
            pe = VinaPotential.Curl(pe, cap, ref dEdr);
            e += pe;
            if (grad is not null && r > 1e-9)
            {
                var gv = d * (dEdr / r);
                grad[i] += gv;
                grad[j] -= gv;
            }
        }
        return e;
    }

    internal static double Rmsd(Vec3[] a, Vec3[] b)
    {
        double s = 0;
        for (int i = 0; i < a.Length; i++) s += Vec3.DistanceSquared(a[i], b[i]);
        return Math.Sqrt(s / Math.Max(1, a.Length));
    }

    private static double RmsdLowerBound(FlexibleLigand lig, Vec3[] a, Vec3[] b)
    {
        // Vina: each atom matched to the nearest atom of the same type in the other pose, both directions
        double One(Vec3[] p, Vec3[] q)
        {
            double s = 0;
            foreach (int i in lig.Scored)
            {
                double best = double.MaxValue;
                foreach (int j in lig.Scored)
                    if (lig.Types[j] == lig.Types[i]) best = Math.Min(best, Vec3.DistanceSquared(p[i], q[j]));
                s += best;
            }
            return Math.Sqrt(s / lig.Scored.Length);
        }
        return Math.Max(One(a, b), One(b, a));
    }

    /// <summary>Per-run state: energy evaluation, BFGS and Monte Carlo.</summary>
    internal sealed class Engine
    {
        private readonly FlexibleLigand _lig;
        private readonly VinaGrid _grid;
        private readonly VinaScoreOptions _w;
        private readonly Vec3[] _x, _ag;
        private readonly int _n;
        public int Evaluations;
        public ExactReceptor? Exact { get; init; }

        public Engine(FlexibleLigand lig, VinaGrid grid, VinaScoreOptions w)
        {
            _lig = lig; _grid = grid; _w = w;
            _x = new Vec3[lig.AtomCount];
            _ag = new Vec3[lig.AtomCount];
            _n = 6 + lig.Torsions.Length;
        }

        /// <summary>Total energy (inter + intra, capped) and conformation gradient.</summary>
        public double Evaluate(LigandConf c, double cap, Span<double> g)
        {
            Evaluations++;
            _lig.Build(c, _x);
            Array.Clear(_ag);
            double e = 0;
            if (Exact is not null) e += Exact.Energy(_lig, _x, _ag);
            else
                foreach (int i in _lig.Scored)
                {
                    e += _grid.Evaluate(_lig.Types[i]!.Value, _x[i], cap, out var gi);
                    _ag[i] += gi;
                }
            e += Intra(_lig, _x, _w, cap, _ag);
            _lig.ConfGradient(c, _x, _ag, g);
            return e;
        }

        public Vec3[] HeavyCoords(LigandConf c)
        {
            _lig.Build(c, _x);
            return _lig.Scored.Select(i => _x[i]).ToArray();
        }

        /// <summary>Vina's BFGS (dense inverse Hessian, backtracking line search).</summary>
        public double Bfgs(LigandConf x, int maxSteps, double cap)
        {
            int n = _n;
            var h = new double[n * n];
            for (int i = 0; i < n; i++) h[i * n + i] = 1;
            var g = new double[n];
            var gNew = new double[n];
            var p = new double[n];
            var y = new double[n];
            double f0 = Evaluate(x, cap, g);
            double fOrig = f0;
            var xOrig = x.Clone();
            var xNew = x.Clone();
            for (int step = 0; step < maxSteps; step++)
            {
                for (int i = 0; i < n; i++)
                {
                    double s = 0;
                    for (int j = 0; j < n; j++) s += h[i * n + j] * g[j];
                    p[i] = -s;
                }
                double pg = 0;
                for (int i = 0; i < n; i++) pg += p[i] * g[i];
                double alpha = 1, f1 = 0;
                for (int trial = 0; trial < 10; trial++)
                {
                    Copy(x, xNew);
                    xNew.Increment(p, alpha);
                    f1 = Evaluate(xNew, cap, gNew);
                    if (f1 - f0 < 1e-4 * alpha * pg) break;
                    alpha *= 0.5;
                }
                for (int i = 0; i < n; i++) y[i] = gNew[i] - g[i];
                f0 = f1;
                Copy(xNew, x);
                double gnorm = 0;
                for (int i = 0; i < n; i++) gnorm += g[i] * g[i];
                if (!(Math.Sqrt(gnorm) >= 1e-5)) break;
                Array.Copy(gNew, g, n);
                if (step == 0)
                {
                    double yy = 0, yp = 0;
                    for (int i = 0; i < n; i++) { yy += y[i] * y[i]; yp += y[i] * p[i]; }
                    if (Math.Abs(yy) > 1e-12) { Array.Clear(h); for (int i = 0; i < n; i++) h[i * n + i] = alpha * yp / yy; }
                }
                // BFGS update of the inverse Hessian with s = alpha p
                double ypd = 0;
                for (int i = 0; i < n; i++) ypd += y[i] * p[i];
                if (alpha * ypd < 1e-12) continue;
                var minusHy = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double s = 0;
                    for (int j = 0; j < n; j++) s += h[i * n + j] * y[j];
                    minusHy[i] = -s;
                }
                double yhy = 0;
                for (int i = 0; i < n; i++) yhy -= y[i] * minusHy[i];
                double r = 1 / (alpha * ypd);
                for (int i = 0; i < n; i++)
                    for (int j = i; j < n; j++)
                    {
                        double v = alpha * r * (minusHy[i] * p[j] + minusHy[j] * p[i]) + alpha * alpha * (r * r * yhy + r) * p[i] * p[j];
                        h[i * n + j] += v;
                        if (j != i) h[j * n + i] += v;
                    }
            }
            if (!(f0 <= fOrig)) { Copy(xOrig, x); f0 = fOrig; }
            return f0;
        }

        private static void Copy(LigandConf from, LigandConf to)
        {
            to.Position = from.Position;
            to.Orientation = from.Orientation;
            Array.Copy(from.Torsions, to.Torsions, from.Torsions.Length);
        }

        private static Vec3 RandomInSphere(Random rng)
        {
            while (true)
            {
                var v = new Vec3(rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);
                if (v.LengthSquared < 1) return v;
            }
        }

        public List<(LigandConf Conf, double E, Vec3[] Heavy)> MonteCarlo(Random rng, int steps, int localSteps, double minRmsd, int keep)
        {
            const double huntCap = 10, authentic = 1000, temperature = 1.2, amplitude = 2.0;
            var lo = _grid.Origin;
            var size = _grid.Size;
            var current = new LigandConf(lo + new Vec3(rng.NextDouble() * size.X, rng.NextDouble() * size.Y, rng.NextDouble() * size.Z),
                Quat.Random(rng), Enumerable.Range(0, _lig.Torsions.Length).Select(_ => (rng.NextDouble() * 2 - 1) * Math.PI).ToArray());
            double currentE = double.MaxValue, bestE = double.MaxValue;
            var output = new List<(LigandConf Conf, double E, Vec3[] Heavy)>();
            int entities = 2 + _lig.Torsions.Length;
            for (int step = 0; step < steps; step++)
            {
                var cand = current.Clone();
                int which = rng.Next(entities);
                if (which == 0) cand.Position += RandomInSphere(rng) * amplitude;
                else if (which == 1)
                {
                    if (_lig.GyrationRadius > 1e-6)
                        cand.Orientation = (Quat.FromRotationVector(RandomInSphere(rng) * (amplitude / _lig.GyrationRadius)) * cand.Orientation).Normalized();
                }
                else cand.Torsions[which - 2] = (rng.NextDouble() * 2 - 1) * Math.PI;
                double e = Bfgs(cand, localSteps, huntCap);
                bool accept = step == 0 || e < currentE || rng.NextDouble() < Math.Exp((currentE - e) / temperature);
                if (!accept) continue;
                current = cand;
                currentE = e;
                if (e < bestE || output.Count < keep)
                {
                    var refined = current.Clone();
                    double re = Bfgs(refined, localSteps, authentic);
                    current = refined;
                    currentE = re;
                    AddToOutput(output, (refined.Clone(), re, HeavyCoords(refined)), minRmsd, keep);
                    if (re < bestE) bestE = re;
                }
            }
            return output;
        }

        private static void AddToOutput(List<(LigandConf Conf, double E, Vec3[] Heavy)> output, (LigandConf Conf, double E, Vec3[] Heavy) item, double minRmsd, int max)
        {
            int close = -1;
            double closest = double.MaxValue;
            for (int i = 0; i < output.Count; i++)
            {
                double r = Rmsd(output[i].Heavy, item.Heavy);
                if (r < minRmsd && r < closest) { closest = r; close = i; }
            }
            if (close >= 0)
            {
                if (item.E < output[close].E) output[close] = item;
            }
            else output.Add(item);
            output.Sort((a, b) => a.E.CompareTo(b.E));
            if (output.Count > max) output.RemoveRange(max, output.Count - max);
        }
    }

    /// <summary>Exact ligand–receptor energy with gradient (Vina's "non-cache"), with the same out-of-box penalty as the grid.</summary>
    internal sealed class ExactReceptor(VinaReceptor receptor, VinaGrid grid, VinaScoreOptions w)
    {
        private readonly NeighborSearch _search = new(receptor.Positions.ToArray(), w.CutOff);

        public double Energy(FlexibleLigand lig, Vec3[] x, Vec3[]? grad)
        {
            double e = 0;
            var found = new List<int>();
            var lo = grid.Origin;
            var hi = grid.Origin + grid.Size;
            foreach (int i in lig.Scored)
            {
                var p = x[i];
                // out-of-box penalty (slope 1e6 per Å), as in Vina
                double px = Math.Clamp(p.X, lo.X, hi.X), py = Math.Clamp(p.Y, lo.Y, hi.Y), pz = Math.Clamp(p.Z, lo.Z, hi.Z);
                e += 1e6 * (Math.Abs(p.X - px) + Math.Abs(p.Y - py) + Math.Abs(p.Z - pz));
                if (grad is not null) grad[i] += new Vec3(Math.Sign(p.X - px), Math.Sign(p.Y - py), Math.Sign(p.Z - pz)) * 1e6;
                found.Clear();
                _search.Query(p, w.CutOff, found);
                var t = lig.Types[i]!.Value;
                foreach (int j in found)
                {
                    var d = p - receptor.Positions[j];
                    double r = d.Length;
                    if (r >= w.CutOff) continue;
                    e += VinaPotential.Energy(t, receptor.Types[j], r, w, out double dEdr);
                    if (grad is not null && r > 1e-9) grad[i] += d * (dEdr / r);
                }
            }
            return e;
        }
    }
}
