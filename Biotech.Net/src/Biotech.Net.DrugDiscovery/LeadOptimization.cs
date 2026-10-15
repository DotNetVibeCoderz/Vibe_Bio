using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>
/// A matched-molecular-pair style transformation: a SMARTS pattern whose match is replaced by a new fragment.
/// </summary>
/// <param name="Name">Transformation name, e.g. "carboxylic acid → tetrazole".</param>
/// <param name="FromSmarts">Pattern to find. Atom maps mark the attachment points kept by the replacement.</param>
/// <param name="ToSmiles">Replacement fragment in SMILES, using the same atom maps for its attachment points.</param>
/// <param name="Rationale">Why a medicinal chemist would make the change.</param>
public sealed record BioisostereTransform(string Name, string FromSmarts, string ToSmiles, string Rationale = "")
{
    private readonly Lazy<SmartsPattern> _pattern = new(() => SmartsPattern.Parse(FromSmarts, mergeHydrogens: false));
    private readonly Lazy<Molecule> _replacement = new(() => Molecule.FromSmiles(ToSmiles));

    /// <summary>The compiled query.</summary>
    public SmartsPattern Pattern => _pattern.Value;

    /// <summary>The replacement fragment as a molecule.</summary>
    public Molecule Replacement => _replacement.Value;

    /// <summary>Whether the transformation applies to <paramref name="mol"/>.</summary>
    /// <param name="mol">Molecule to test.</param>
    public bool Applies(Molecule mol) => Pattern.IsMatch(mol);
}

/// <summary>
/// Classic bioisosteric replacements as SMARTS→SMILES transformations, applied by
/// <see cref="LeadOptimizer.Enumerate"/>.
/// </summary>
/// <remarks>
/// The catalogue is a hand-picked selection of textbook replacements (Meanwell, J. Med. Chem. 2011, 54, 2529), not
/// an exhaustive MMP database mined from ChEMBL. Atom map 1 marks the attachment point that is kept.
/// </remarks>
public static class Bioisosteres
{
    /// <summary>The built-in transformations.</summary>
    public static IReadOnlyList<BioisostereTransform> Default { get; } =
    [
        new("carboxylic acid → tetrazole", "[C:1](=O)[OX2H1]", "[C:1]1=NN=NN1", "keeps acidity, raises metabolic stability and permeability"),
        new("carboxylic acid → acyl sulfonamide", "[C:1](=O)[OX2H1]", "[C:1](=O)NS(C)(=O)=O", "acidic isostere with extra vector"),
        new("carboxylic acid → hydroxamic acid", "[C:1](=O)[OX2H1]", "[C:1](=O)NO", "metal-binding isostere"),
        new("ester → amide", "[C:1](=O)O[CH3]", "[C:1](=O)NC", "removes an esterase liability"),
        new("phenyl → pyridine (3-aza)", "[c:1]1[c:2][c:3][cH:4][c:5][c:6]1", "[c:1]1[c:2][c:3][n:4][c:5][c:6]1", "lowers logP, adds an acceptor"),
        new("benzene → cyclohexane", "[c:1]1[c:2][c:3][c:4][c:5][c:6]1", "[C:1]1[C:2][C:3][C:4][C:5][C:6]1", "raises Fsp3, removes planarity"),
        new("pyridine → pyrimidine", "[c:1]1[c:2][cH:3][n:4][c:5][c:6]1", "[c:1]1[c:2][n:3][n:4][c:5][c:6]1", "adds an acceptor, lowers logP"),
        new("amide → 1,2,4-oxadiazole", "[C:1](=O)N([H])C", "[c:1]1noc(C)n1", "amide isostere, resists hydrolysis"),
        new("amide → reverse amide", "[C:1](=O)NC", "[C:1]NC(C)=O", "flips the hydrogen-bond vector"),
        new("methyl → trifluoromethyl", "[CX4H3:1]", "[C:1](F)(F)F", "blocks metabolism, raises lipophilicity"),
        new("methyl → chlorine", "[c][CX4H3:1]", "[Cl:1]", "similar size, different electronics"),
        new("hydroxyl → fluorine", "[CX4][OX2H1:1]", "[F:1]", "removes a donor, keeps size"),
        new("tert-butyl → trifluoromethylcyclopropyl", "[CX4:1]([CH3])([CH3])[CH3]", "[C:1]1(C(F)(F)F)CC1", "lowers metabolic liability"),
        new("piperazine → piperidine", "[NX3:1]1[CH2:2][CH2:3][NX3;H1:4][CH2:5][CH2:6]1", "[N:1]1[C:2][C:3][C:4][C:5][C:6]1", "removes a basic centre (hERG mitigation)"),
        new("morpholine → piperazine", "[N:1]1[C:2][C:3][O:4][C:5][C:6]1", "[N:1]1[C:2][C:3][N:4][C:5][C:6]1", "adds a basic centre and a donor vector"),
        new("primary amide → nitrile", "[C:1](=O)[NX3H2]", "[C:1]#N", "smaller acceptor"),
        new("ketone → sulfone", "[C:1](=O)C", "[S:1](=O)(=O)C", "stronger acceptor, more polar"),
        new("thioether → sulfone", "[SX2;$(S([#6])[#6]):1]", "[S:1](=O)=O", "blocks oxidative metabolism"),
        new("phenol → aniline", "[c:1][OX2H1]", "[c:1]N", "donor retained, acidity lowered"),
        new("nitro → sulfonamide", "[c:1][N+](=O)[O-]", "[c:1]S(N)(=O)=O", "removes a toxicophore"),
    ];

    /// <summary>Transformations that apply to <paramref name="mol"/>.</summary>
    public static IReadOnlyList<BioisostereTransform> Applicable(Molecule mol, IReadOnlyList<BioisostereTransform>? catalogue = null)
        => (catalogue ?? Default).Where(t => t.Applies(mol)).ToList();
}

/// <summary>One enumerated analogue.</summary>
/// <param name="Smiles">Canonical SMILES of the analogue.</param>
/// <param name="Molecule">The analogue.</param>
/// <param name="Transform">Which transformation produced it.</param>
/// <param name="ParentSmiles">Canonical SMILES of the molecule it came from.</param>
public sealed record Analogue(string Smiles, Molecule Molecule, BioisostereTransform Transform, string ParentSmiles);

/// <summary>Weights and desirability windows of a multi-parameter optimisation score.</summary>
public sealed record MpoProfile
{
    /// <summary>Desired molecular-weight window (full marks inside, tapering to 0 at the outer bound).</summary>
    public (double Ideal, double Limit) MolecularWeight { get; init; } = (400, 550);

    /// <summary>Desired logP window.</summary>
    public (double Ideal, double Limit) LogP { get; init; } = (3.0, 5.5);

    /// <summary>Desired TPSA window.</summary>
    public (double Ideal, double Limit) Tpsa { get; init; } = (90, 140);

    /// <summary>Desired hydrogen-bond-donor window.</summary>
    public (double Ideal, double Limit) HBondDonors { get; init; } = (2, 5);

    /// <summary>Desired rotatable-bond window.</summary>
    public (double Ideal, double Limit) RotatableBonds { get; init; } = (6, 11);

    /// <summary>Weight of the QED term.</summary>
    public double QedWeight { get; init; } = 1.0;

    /// <summary>Weight of each physicochemical term.</summary>
    public double PropertyWeight { get; init; } = 1.0;

    /// <summary>Penalty subtracted per matched structural alert.</summary>
    public double AlertPenalty { get; init; } = 0.1;

    /// <summary>The CNS-oriented profile: smaller, less polar, fewer donors.</summary>
    public static MpoProfile Cns { get; } = new()
    {
        MolecularWeight = (360, 450),
        LogP = (2.5, 4.5),
        Tpsa = (60, 90),
        HBondDonors = (1, 3),
        RotatableBonds = (4, 8),
    };

    /// <summary>The default oral-drug profile.</summary>
    public static MpoProfile Oral { get; } = new();
}

/// <summary>
/// Multi-parameter optimisation scoring in the spirit of Wager's CNS-MPO (ACS Chem. Neurosci. 2010): each property
/// gets a 0–1 desirability that is 1 inside the ideal window and falls linearly to 0 at the limit, and the score is
/// the weighted mean, minus a penalty per structural alert.
/// </summary>
public static class MpoScore
{
    /// <summary>Desirability of <paramref name="value"/>: 1 up to <paramref name="ideal"/>, 0 beyond <paramref name="limit"/>.</summary>
    public static double Desirability(double value, double ideal, double limit)
    {
        if (limit > ideal) return value <= ideal ? 1 : value >= limit ? 0 : (limit - value) / (limit - ideal);
        return value >= ideal ? 1 : value <= limit ? 0 : (value - limit) / (ideal - limit);
    }

    /// <summary>The per-property desirabilities.</summary>
    public static IReadOnlyDictionary<string, double> Components(Molecule mol, MpoProfile? profile = null)
    {
        profile ??= MpoProfile.Oral;
        var d = MolecularDescriptors.Compute(mol);
        return new Dictionary<string, double>
        {
            ["MW"] = Desirability(d.MolecularWeight, profile.MolecularWeight.Ideal, profile.MolecularWeight.Limit),
            ["LogP"] = Desirability(d.LogP, profile.LogP.Ideal, profile.LogP.Limit),
            ["TPSA"] = Desirability(d.Tpsa, profile.Tpsa.Ideal, profile.Tpsa.Limit),
            ["HBD"] = Desirability(d.HBondDonors, profile.HBondDonors.Ideal, profile.HBondDonors.Limit),
            ["RotB"] = Desirability(d.RotatableBonds, profile.RotatableBonds.Ideal, profile.RotatableBonds.Limit),
            ["QED"] = Qed.Compute(mol),
        };
    }

    /// <summary>The MPO score in [0, 1] (higher is better).</summary>
    public static double Compute(Molecule mol, MpoProfile? profile = null)
    {
        profile ??= MpoProfile.Oral;
        var c = Components(mol, profile);
        double sum = 0, weight = 0;
        foreach (var (name, value) in c)
        {
            double w = name == "QED" ? profile.QedWeight : profile.PropertyWeight;
            sum += w * value;
            weight += w;
        }
        double score = weight == 0 ? 0 : sum / weight;
        score -= profile.AlertPenalty * StructuralAlerts.Match(mol).Count;
        return Math.Clamp(score, 0, 1);
    }
}

/// <summary>
/// Lead optimisation: enumerates analogues with bioisosteric transformations and ranks them by a multi-parameter
/// optimisation score.
/// </summary>
public static class LeadOptimizer
{
    /// <summary>
    /// Applies every applicable transformation to <paramref name="lead"/>, optionally for several rounds
    /// (<paramref name="depth"/>), returning distinct analogues (the parent is excluded).
    /// </summary>
    public static IReadOnlyList<Analogue> Enumerate(Molecule lead, IReadOnlyList<BioisostereTransform>? catalogue = null,
        int depth = 1, int maxAnalogues = 1000)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(depth, 1);
        var transforms = catalogue ?? Bioisosteres.Default;
        var results = new List<Analogue>();
        var seen = new HashSet<string> { lead.ToSmiles(includeStereo: false) };
        var frontier = new List<Molecule> { lead };
        for (int round = 0; round < depth && results.Count < maxAnalogues; round++)
        {
            var next = new List<Molecule>();
            foreach (var parent in frontier)
            {
                string parentSmiles = parent.ToSmiles(includeStereo: false);
                foreach (var t in transforms)
                {
                    foreach (var product in Apply(parent, t))
                    {
                        string smi;
                        try { smi = product.ToSmiles(includeStereo: false); }
                        catch (InvalidOperationException) { continue; }
                        if (!seen.Add(smi)) continue;
                        results.Add(new Analogue(smi, product, t, parentSmiles));
                        next.Add(product);
                        if (results.Count >= maxAnalogues) return results;
                    }
                }
            }
            frontier = next;
        }
        return results;
    }

    /// <summary>
    /// Applies one transformation to every match, returning the products. Matched atoms that carry no atom map are
    /// deleted; the replacement fragment is attached at the mapped atoms.
    /// </summary>
    public static IReadOnlyList<Molecule> Apply(Molecule mol, BioisostereTransform transform)
    {
        var products = new List<Molecule>();
        var pattern = transform.Pattern;
        foreach (var match in pattern.Matches(mol))
        {
            var product = ApplyAt(mol, transform, match);
            if (product is not null) products.Add(product);
        }
        return products;
    }

    private static Molecule? ApplyAt(Molecule mol, BioisostereTransform transform, int[] match)
    {
        var pattern = transform.Pattern;
        // anchors: query atom → molecule atom, for mapped query atoms
        var anchors = new Dictionary<int, int>();
        for (int q = 0; q < pattern.AtomCount; q++)
        {
            int map = pattern.AtomMap(q);
            if (map > 0) anchors[map] = match[q];
        }
        var remove = new HashSet<int>(match);
        foreach (int a in anchors.Values) remove.Remove(a);
        // A removed atom may only be bonded inside the match: otherwise the transformation would cut a substituent off.
        foreach (int a in remove)
            foreach (int nb in mol.Neighbors(a))
                if (!remove.Contains(nb) && !anchors.ContainsValue(nb)) return null;

        var product = mol.Clone();
        product.Coordinates = null;
        // 1. remove the atoms replaced by the new fragment
        var keepMap = product.RemoveAtoms(remove);
        // 2. add the replacement fragment, remembering where its mapped atoms land
        var frag = transform.Replacement;
        var fragIndex = new int[frag.AtomCount];
        var fragAnchor = new Dictionary<int, int>();
        for (int i = 0; i < frag.AtomCount; i++)
        {
            var a = frag.Atoms[i];
            int map = a.MapNumber;
            if (map > 0 && anchors.ContainsKey(map))
            {
                // this fragment atom replaces the anchor atom itself
                int target = keepMap[anchors[map]];
                if (target < 0) return null;
                fragIndex[i] = target;
                fragAnchor[map] = target;
                ref var t = ref product.AtomAt(target);
                t.AtomicNumber = a.AtomicNumber;
                t.FormalCharge = a.FormalCharge;
                t.IsAromatic = a.IsAromatic;
                t.Chirality = ChiralTag.None;
                t.NoImplicit = false;
                t.ExplicitHydrogens = 0;
                product.SetChiralReference(target, null);
                continue;
            }
            var copy = new Atom(a.AtomicNumber, a.FormalCharge) { IsAromatic = a.IsAromatic, Isotope = a.Isotope };
            fragIndex[i] = product.AddAtom(copy);
        }
        foreach (var b in frag.Bonds)
        {
            int x = fragIndex[b.Begin], y = fragIndex[b.End];
            if (x == y || product.BondBetween(x, y) >= 0) continue;
            product.AddBond(x, y, b.Order);
        }
        if (fragAnchor.Count != anchors.Count) return null;
        if (product.Fragments().Count != mol.Fragments().Count) return null;
        product.UpdateImplicitHydrogens();
        try
        {
            Aromaticity.Perceive(product);
            _ = product.ToSmiles();
        }
        catch (InvalidOperationException) { return null; }
        product.Name = mol.Name.Length == 0 ? transform.Name : $"{mol.Name} ({transform.Name})";
        return product;
    }

    /// <summary>
    /// Enumerates analogues and ranks them by <see cref="MpoScore"/>, best first. The lead itself is included so
    /// improvements can be read off directly.
    /// </summary>
    public static IReadOnlyList<(Analogue? Analogue, Molecule Molecule, double Score)> Optimize(Molecule lead,
        MpoProfile? profile = null, IReadOnlyList<BioisostereTransform>? catalogue = null, int depth = 1, int maxAnalogues = 1000)
    {
        var list = new List<(Analogue?, Molecule, double)> { (null, lead, MpoScore.Compute(lead, profile)) };
        foreach (var a in Enumerate(lead, catalogue, depth, maxAnalogues))
            list.Add((a, a.Molecule, MpoScore.Compute(a.Molecule, profile)));
        return list.OrderByDescending(t => t.Item3).ToList();
    }
}
