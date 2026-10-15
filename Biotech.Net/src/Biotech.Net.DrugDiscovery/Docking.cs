using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;
using SbAtom = Biotech.Net.StructuralBiology.Atom;

namespace Biotech.Net.DrugDiscovery;

/// <summary>
/// An atom of a ligand pose: a 3-D position plus the typing the scoring function needs.
/// </summary>
/// <param name="Position">Coordinates in Å.</param>
/// <param name="Element">Upper-case element symbol.</param>
/// <param name="IsHydrophobic">Carbon or halogen not bonded to a heteroatom.</param>
/// <param name="IsDonor">Can donate a hydrogen bond (N/O carrying H).</param>
/// <param name="IsAcceptor">Can accept a hydrogen bond (N/O).</param>
public readonly record struct PoseAtom(Vec3 Position, string Element, bool IsHydrophobic, bool IsDonor, bool IsAcceptor)
{
    /// <summary>Van der Waals radius in Å.</summary>
    public double VdwRadius => ElementData.VdwRadius(Element);

    /// <summary>Whether the atom is a hydrogen.</summary>
    public bool IsHydrogen => Element == "H";
}

/// <summary>A ligand conformation to be scored: typed atoms in 3-D.</summary>
public sealed class LigandPose
{
    /// <summary>Creates a pose from typed atoms.</summary>
    public LigandPose(string name, IReadOnlyList<PoseAtom> atoms)
    {
        Name = name;
        Atoms = atoms;
    }

    /// <summary>Pose name.</summary>
    public string Name { get; }

    /// <summary>The pose atoms (hydrogens may be present; the scoring function ignores them).</summary>
    public IReadOnlyList<PoseAtom> Atoms { get; }

    /// <summary>Heavy atoms only.</summary>
    public IEnumerable<PoseAtom> HeavyAtoms => Atoms.Where(a => !a.IsHydrogen);

    /// <summary>Geometric centre of the heavy atoms.</summary>
    public Vec3 Center
    {
        get
        {
            var heavy = HeavyAtoms.ToList();
            if (heavy.Count == 0) return Vec3.Zero;
            var s = heavy.Aggregate(Vec3.Zero, (acc, a) => acc + a.Position);
            return s / heavy.Count;
        }
    }

    /// <summary>A copy of the pose translated by <paramref name="delta"/>.</summary>
    public LigandPose Translated(Vec3 delta) => new(Name, Atoms.Select(a => a with { Position = a.Position + delta }).ToList());

    /// <summary>A copy rotated by <paramref name="rotation"/> about <paramref name="origin"/>.</summary>
    public LigandPose Rotated(Matrix3 rotation, Vec3 origin) => new(Name, Atoms.Select(a => a with { Position = origin + rotation * (a.Position - origin) }).ToList());

    /// <summary>
    /// Builds a pose from a molecule with 3-D coordinates, typing each atom for the scoring function.
    /// </summary>
    /// <exception cref="InvalidOperationException">The molecule has no coordinates.</exception>
    public static LigandPose FromMolecule(Molecule mol)
    {
        var coords = mol.Coordinates ?? throw new InvalidOperationException($"Molecule '{mol.Name}' has no coordinates; generate or read a 3-D conformer first.");
        var atoms = new List<PoseAtom>(mol.AtomCount);
        for (int i = 0; i < mol.AtomCount; i++)
        {
            var a = mol.Atoms[i];
            if (a.AtomicNumber == 0) continue;
            string el = a.Symbol.ToUpperInvariant();
            bool hetero = a.AtomicNumber is 7 or 8 or 15 or 16;
            bool nearHetero = false;
            foreach (int nb in mol.Neighbors(i)) if (mol.Atoms[nb].AtomicNumber is 7 or 8) nearHetero = true;
            bool hydrophobic = (a.AtomicNumber == 6 && !nearHetero) || a.AtomicNumber is 9 or 17 or 35 or 53;
            bool donor = a.AtomicNumber is 7 or 8 && mol.TotalHydrogens(i) > 0;
            bool acceptor = a.AtomicNumber is 7 or 8;
            atoms.Add(new PoseAtom(new Vec3(coords[i].X, coords[i].Y, coords[i].Z), el, hydrophobic, donor, acceptor));
        }
        return new LigandPose(mol.Name, atoms);
    }

    /// <summary>Builds a pose from a structure residue (e.g. a crystallographic ligand).</summary>
    public static LigandPose FromResidue(Residue residue)
    {
        var atoms = residue.Atoms.Where(a => a.Element != "H").Select(a =>
        {
            bool hetero = a.Element is "N" or "O" or "P" or "S";
            bool hydrophobic = a.Element is "C" or "F" or "CL" or "BR" or "I";
            return new PoseAtom(a.Position, a.Element, hydrophobic && !hetero, a.Element is "N" or "O", a.Element is "N" or "O");
        }).ToList();
        return new LigandPose(residue.Name, atoms);
    }
}

/// <summary>Weights and cut-offs of the <see cref="VinaScore"/> empirical function.</summary>
public sealed record VinaScoreOptions
{
    /// <summary>Gaussian 1 weight (AutoDock Vina: −0.0356).</summary>
    public double Gauss1Weight { get; init; } = -0.035579;

    /// <summary>Gaussian 2 weight (−0.00516).</summary>
    public double Gauss2Weight { get; init; } = -0.005156;

    /// <summary>Repulsion weight (0.840).</summary>
    public double RepulsionWeight { get; init; } = 0.840245;

    /// <summary>Hydrophobic weight (−0.0351).</summary>
    public double HydrophobicWeight { get; init; } = -0.035069;

    /// <summary>Hydrogen-bond weight (−0.587).</summary>
    public double HydrogenBondWeight { get; init; } = -0.587439;

    /// <summary>Torsional-entropy weight (0.0586), multiplied by the number of rotatable bonds.</summary>
    public double RotatableWeight { get; init; } = 0.05846;

    /// <summary>Pair cut-off in Å (8).</summary>
    public double CutOff { get; init; } = 8.0;

    /// <summary>Default parameters (the published Vina values).</summary>
    public static VinaScoreOptions Default { get; } = new();
}

/// <summary>The terms of one score.</summary>
/// <param name="Total">Weighted sum (the "score" in kcal/mol-like units; more negative is better).</param>
/// <param name="Gauss1">Unweighted gauss1 sum.</param>
/// <param name="Gauss2">Unweighted gauss2 sum.</param>
/// <param name="Repulsion">Unweighted repulsion sum.</param>
/// <param name="Hydrophobic">Unweighted hydrophobic sum.</param>
/// <param name="HydrogenBond">Unweighted H-bond sum.</param>
/// <param name="RotatableBonds">Rotatable bonds used for the entropy term.</param>
/// <param name="ContactCount">Ligand–receptor atom pairs within the cut-off.</param>
public sealed record VinaScoreTerms(double Total, double Gauss1, double Gauss2, double Repulsion, double Hydrophobic,
    double HydrogenBond, int RotatableBonds, int ContactCount)
{
    /// <summary>The intermolecular score before the torsional-entropy correction.</summary>
    public double Intermolecular { get; init; }

    /// <summary>
    /// Vina's "affinity" estimate: the intermolecular score divided by 1 + w·N<sub>rot</sub>, in kcal/mol.
    /// </summary>
    public double Affinity { get; init; }

    /// <inheritdoc/>
    public override string ToString() => $"{Affinity:F2} kcal/mol (inter {Intermolecular:F2}, {ContactCount} contacts)";
}

/// <summary>
/// A rigid-receptor empirical scoring function in the form of AutoDock Vina (Trott &amp; Olson, J. Comput. Chem.
/// 2010, 31, 455): the five distance-dependent terms (two Gaussians, a repulsion, a hydrophobic and a hydrogen-bond
/// term) over the surface distance d = r − R<sub>i</sub> − R<sub>j</sub>, with the published weights, plus Vina's
/// torsional correction.
/// </summary>
/// <remarks>
/// This scores a <b>given</b> pose; it does not search conformations or perform flexible docking.
/// <see cref="Optimize"/> does a small rigid-body local search (translation and rotation only) and
/// <see cref="ScoreInPocket"/> seeds that search at a detected pocket centre. Atom typing is simplified compared
/// to Vina's X-score types, so absolute values are indicative rather than reproducing Vina exactly; relative
/// rankings of poses of the same ligand are the intended use.
/// </remarks>
public static class VinaScore
{
    /// <summary>Scores one pose against the receptor atoms.</summary>
    public static VinaScoreTerms Score(LigandPose pose, IReadOnlyList<SbAtom> receptor, int rotatableBonds = 0, VinaScoreOptions? options = null)
    {
        options ??= VinaScoreOptions.Default;
        var protein = receptor.Where(a => a.Element != "H").ToList();
        var search = NeighborSearch.ForAtoms(protein, Math.Max(5.0, options.CutOff));
        double g1 = 0, g2 = 0, rep = 0, phob = 0, hb = 0;
        int contacts = 0;
        var found = new List<int>();
        foreach (var la in pose.HeavyAtoms)
        {
            found.Clear();
            search.Query(la.Position, options.CutOff, found);
            double lr = la.VdwRadius;
            foreach (int k in found)
            {
                var pa = protein[k];
                double r = Vec3.Distance(la.Position, pa.Position);
                if (r > options.CutOff) continue;
                double d = r - lr - pa.VdwRadius;
                contacts++;
                g1 += Math.Exp(-(d / 0.5) * (d / 0.5));
                g2 += Math.Exp(-((d - 3.0) / 2.0) * ((d - 3.0) / 2.0));
                if (d < 0) rep += d * d;
                bool pHydrophobic = pa.Element is "C" or "F" or "CL" or "BR" or "I" && !pa.IsBackbone;
                if (la.IsHydrophobic && pHydrophobic) phob += Hydrophobic(d);
                bool pDonor = pa.Element is "N" or "O";
                bool pAcceptor = pa.Element is "N" or "O";
                if ((la.IsDonor && pAcceptor) || (la.IsAcceptor && pDonor)) hb += HydrogenBond(d);
            }
        }
        double inter = options.Gauss1Weight * g1 + options.Gauss2Weight * g2 + options.RepulsionWeight * rep
                     + options.HydrophobicWeight * phob + options.HydrogenBondWeight * hb;
        double affinity = inter / (1 + options.RotatableWeight * rotatableBonds);
        return new VinaScoreTerms(inter + options.RotatableWeight * rotatableBonds, g1, g2, rep, phob, hb, rotatableBonds, contacts)
        {
            Intermolecular = inter,
            Affinity = affinity,
        };
    }

    /// <summary>Scores a pose against every heavy atom of <paramref name="structure"/>'s first model.</summary>
    public static VinaScoreTerms Score(LigandPose pose, Structure structure, int rotatableBonds = 0, VinaScoreOptions? options = null)
        => Score(pose, structure.FirstModel.Atoms.Where(a => !a.IsHetero).ToList(), rotatableBonds, options);

    /// <summary>Scores a molecule with 3-D coordinates, deriving the rotatable-bond count from the graph.</summary>
    public static VinaScoreTerms Score(Molecule ligand, Structure receptor, VinaScoreOptions? options = null)
        => Score(LigandPose.FromMolecule(ligand), receptor, MolecularDescriptors.RotatableBonds(ligand), options);

    /// <summary>Vina's hydrophobic term: 1 below 0.5 Å surface distance, falling linearly to 0 at 1.5 Å.</summary>
    public static double Hydrophobic(double surfaceDistance) => surfaceDistance switch
    {
        < 0.5 => 1,
        > 1.5 => 0,
        _ => 1.5 - surfaceDistance,
    };

    /// <summary>Vina's hydrogen-bond term: 1 below −0.7 Å surface distance, falling linearly to 0 at 0 Å.</summary>
    public static double HydrogenBond(double surfaceDistance) => surfaceDistance switch
    {
        < -0.7 => 1,
        > 0 => 0,
        _ => -surfaceDistance / 0.7,
    };

    /// <summary>
    /// A rigid-body local search (Monte-Carlo with decreasing step size) that keeps the best-scoring pose. The ligand
    /// conformation is never changed, so this refines placement rather than docking flexibly.
    /// </summary>
    public static (LigandPose Pose, VinaScoreTerms Score) Optimize(LigandPose pose, IReadOnlyList<SbAtom> receptor, int rotatableBonds = 0,
        int iterations = 300, double maxTranslation = 2.0, int seed = 12345, VinaScoreOptions? options = null)
    {
        var rng = new Random(seed);
        var best = pose;
        var bestScore = Score(pose, receptor, rotatableBonds, options);
        var current = pose;
        var currentScore = bestScore;
        for (int it = 0; it < iterations; it++)
        {
            double t = 1.0 - (double)it / iterations;
            double step = maxTranslation * t;
            var delta = new Vec3((rng.NextDouble() * 2 - 1) * step, (rng.NextDouble() * 2 - 1) * step, (rng.NextDouble() * 2 - 1) * step);
            var axis = new Vec3(rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);
            if (axis.Length < 1e-6) axis = Vec3.UnitZ;
            var rot = Matrix3.FromAxisAngle(axis.Normalized(), (rng.NextDouble() * 2 - 1) * 30 * t);
            var centre = current.Center;
            var candidate = current.Rotated(rot, centre).Translated(delta);
            var score = Score(candidate, receptor, rotatableBonds, options);
            if (score.Affinity < currentScore.Affinity)
            {
                current = candidate;
                currentScore = score;
                if (score.Affinity < bestScore.Affinity) { best = candidate; bestScore = score; }
            }
        }
        return (best, bestScore);
    }

    /// <summary>
    /// Places the ligand at the centre of <paramref name="pocket"/>, runs <see cref="Optimize"/> and returns the
    /// refined pose. A convenience entry point for "score this ligand in this pocket".
    /// </summary>
    public static (LigandPose Pose, VinaScoreTerms Score) ScoreInPocket(LigandPose pose, Structure receptor, Pocket pocket,
        int rotatableBonds = 0, int iterations = 300, int seed = 12345, VinaScoreOptions? options = null)
    {
        var atoms = receptor.FirstModel.Atoms.Where(a => !a.IsHetero && a.Element != "H").ToList();
        var centred = pose.Translated(pocket.Center - pose.Center);
        return Optimize(centred, atoms, rotatableBonds, iterations, 2.0, seed, options);
    }

    /// <summary>
    /// Ranks several ligand poses against one receptor, best (most negative affinity) first.
    /// </summary>
    public static IReadOnlyList<(LigandPose Pose, VinaScoreTerms Score)> Rank(IEnumerable<LigandPose> poses, IReadOnlyList<SbAtom> receptor,
        Func<LigandPose, int>? rotatableBonds = null, VinaScoreOptions? options = null)
        => poses.Select(p => (Pose: p, Score: Score(p, receptor, rotatableBonds?.Invoke(p) ?? 0, options)))
                .OrderBy(x => x.Score.Affinity)
                .ToList();
}
