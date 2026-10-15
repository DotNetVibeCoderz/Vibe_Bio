using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>A named SMARTS alert.</summary>
/// <param name="Family">Alert family (e.g. "PAINS", "Brenk").</param>
/// <param name="Name">Alert name.</param>
/// <param name="Smarts">Pattern.</param>
public sealed record StructuralAlert(string Family, string Name, string Smarts)
{
    private readonly Lazy<SmartsPattern> _pattern = new(() => SmartsPattern.Parse(Smarts));

    /// <summary>The compiled pattern.</summary>
    public SmartsPattern Pattern => _pattern.Value;

    /// <summary>Whether the alert fires on <paramref name="mol"/>.</summary>
    public bool Matches(Molecule mol) => Pattern.IsMatch(mol);
}

/// <summary>
/// Structural-alert catalogues.
/// </summary>
/// <remarks>
/// These are <b>curated subsets written for this library, not the complete published catalogues</b>. The PAINS
/// entries cover frequently reported interference families from Baell &amp; Holloway (J. Med. Chem. 2010, 53, 2719)
/// written as simplified SMARTS; the Brenk entries cover common reactive and otherwise unsuitable groups from
/// Brenk et al. (ChemMedChem 2008, 3, 435). Matches are indicative and a clean result does not mean a compound is
/// free of the full published alerts.
/// </remarks>
public static class StructuralAlerts
{
    /// <summary>PAINS-type interference alerts (curated subset).</summary>
    public static IReadOnlyList<StructuralAlert> Pains { get; } =
    [
        new("PAINS", "quinone_A", "O=C1C=CC(=O)C=C1"),
        new("PAINS", "quinone_naphthalene", "O=C1C=CC(=O)c2ccccc21"),
        new("PAINS", "catechol_A", "c1cc([OH])c([OH])cc1"),
        new("PAINS", "hydroquinone", "[OH]c1ccc([OH])cc1"),
        new("PAINS", "ene_rhodanine", "C(=[CH])[CX3](=[O,S])[NX3][CX3]=[S,O]"),
        new("PAINS", "rhodanine_core", "S=C1SC(=O)[NX3]C1"),
        new("PAINS", "thiazolidinone_ene", "O=C1[NX3]C(=S)SC1=C"),
        new("PAINS", "phenol_sulfonamide", "[OH]c1ccccc1S(=O)(=O)[NX3]"),
        new("PAINS", "azo_compound", "[#6]N=N[#6]"),
        new("PAINS", "hydrazone_aryl", "[c,C]C=N[NX3;H1,H2]"),
        new("PAINS", "mannich_base_phenol", "[OH]c1ccccc1C[NX3]"),
        new("PAINS", "curcumin_like", "O=C(C=Cc1ccccc1)CC(=O)C=Cc1ccccc1"),
        new("PAINS", "alkylidene_barbiturate", "O=C1[NX3]C(=O)[NX3]C(=O)C1=C"),
        new("PAINS", "anil_di_alkyl", "c1cc([NX3;H1,H2])ccc1[NX3;H1,H2]"),
        new("PAINS", "catechol_ester", "O=C(O)c1cc([OH])c([OH])cc1"),
        new("PAINS", "nitroaryl", "[$(c[N+](=O)[O-])]"),
        new("PAINS", "cyanoacrylate", "N#CC(=C)[CX3]=O"),
        new("PAINS", "thiocarbonyl_aryl", "c[CX3]=[SX1]"),
        new("PAINS", "beta_keto_anil", "O=C(C[CX3]=O)[NX3]c1ccccc1"),
        new("PAINS", "indole_3_ene", "C=Cc1c[nH]c2ccccc12"),
    ];

    /// <summary>Brenk-type reactive / undesirable-group alerts (curated subset).</summary>
    public static IReadOnlyList<StructuralAlert> Brenk { get; } =
    [
        new("Brenk", "acyl_halide", "[CX3](=O)[F,Cl,Br,I]"),
        new("Brenk", "sulfonyl_halide", "S(=O)(=O)[F,Cl,Br,I]"),
        new("Brenk", "alkyl_halide", "[CX4][Cl,Br,I]"),
        new("Brenk", "aldehyde", "[CX3H1](=O)[#6]"),
        new("Brenk", "isocyanate_isothiocyanate", "[NX2]=C=[O,S]"),
        new("Brenk", "epoxide_aziridine", "[OX2,NX3]1[CX4][CX4]1"),
        new("Brenk", "michael_acceptor", "[CX3]=[CX3][CX3]=[OX1]"),
        new("Brenk", "acyclic_imine", "[CX3;!R]=[NX2;!R][#6]"),
        new("Brenk", "azide", "[NX1]~[NX2]~[NX2,NX1]"),
        new("Brenk", "azo", "[#6][NX2]=[NX2][#6]"),
        new("Brenk", "diazo", "[#6]=[N+]=[N-]"),
        new("Brenk", "nitro", "[N+](=O)[O-]"),
        new("Brenk", "nitroso", "[NX2]=[OX1]"),
        new("Brenk", "hydrazine", "[NX3;H1,H2][NX3;H1,H2]"),
        new("Brenk", "hydroxylamine", "[NX3][OX2H1]"),
        new("Brenk", "peroxide", "[OX2][OX2]"),
        new("Brenk", "disulfide", "[SX2][SX2]"),
        new("Brenk", "thiol", "[SX2H1]"),
        new("Brenk", "phosphorus", "[#15]"),
        new("Brenk", "silicon", "[#14]"),
        new("Brenk", "triflate", "OS(=O)(=O)C(F)(F)F"),
        new("Brenk", "perhalo_ketone", "[CX3](=O)C([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]"),
        new("Brenk", "aliphatic_long_chain", "[CX4;!R][CX4;!R][CX4;!R][CX4;!R][CX4;!R][CX4;!R][CX4;!R]"),
        new("Brenk", "cyanamide", "[NX3]C#N"),
        new("Brenk", "crown_ether", "[OX2;R][CX4;R][CX4;R][OX2;R][CX4;R][CX4;R][OX2;R]"),
        new("Brenk", "heavy_metal", "[Hg,Pb,Cd,As,Sb,Se,Te,Tl,Be,Ba,Sn,Cr,Ni,Co]"),
        new("Brenk", "triple_bond_terminal", "[CX2;H1]#[CX2]"),
        new("Brenk", "oxime", "[CX3]=[NX2][OX2H1]"),
        new("Brenk", "carbodiimide", "[NX2]=[CX2]=[NX2]"),
        new("Brenk", "charged_quaternary_n", "[NX4+;!$([NX4+][O-])]"),
    ];

    /// <summary>Both catalogues.</summary>
    public static IReadOnlyList<StructuralAlert> All { get; } = [.. Pains, .. Brenk];

    /// <summary>The catalogue of the given family ("PAINS", "Brenk" or "all").</summary>
    public static IReadOnlyList<StructuralAlert> Family(string family) => family.ToUpperInvariant() switch
    {
        "PAINS" => Pains,
        "BRENK" => Brenk,
        "ALL" => All,
        _ => throw new ArgumentException($"Unknown alert family '{family}'.", nameof(family)),
    };

    /// <summary>Alerts of <paramref name="catalogue"/> (default all) that fire on <paramref name="mol"/>.</summary>
    public static IReadOnlyList<StructuralAlert> Match(Molecule mol, IReadOnlyList<StructuralAlert>? catalogue = null)
        => (catalogue ?? All).Where(a => a.Matches(mol)).ToList();

    /// <summary>Whether no alert of <paramref name="catalogue"/> fires.</summary>
    public static bool IsClean(Molecule mol, IReadOnlyList<StructuralAlert>? catalogue = null)
        => !(catalogue ?? All).Any(a => a.Matches(mol));
}

/// <summary>The pharmacophoric feature types counted by <see cref="Pharmacophore"/>.</summary>
public enum PharmacophoreFeature
{
    /// <summary>Hydrogen-bond donor (N–H, O–H).</summary>
    Donor,
    /// <summary>Hydrogen-bond acceptor.</summary>
    Acceptor,
    /// <summary>Positively ionisable centre.</summary>
    PositiveIonisable,
    /// <summary>Negatively ionisable centre.</summary>
    NegativeIonisable,
    /// <summary>Aromatic ring centre.</summary>
    Aromatic,
    /// <summary>Hydrophobic atom or group.</summary>
    Hydrophobic,
    /// <summary>Halogen-bond donor (C–X with X = Cl, Br, I).</summary>
    HalogenDonor,
}

/// <summary>
/// Pharmacophore-feature typing with the SMARTS definitions used by RDKit's simplified "BaseFeatures" set:
/// counting features, locating them, and comparing molecules by feature profile.
/// </summary>
public static class Pharmacophore
{
    private static readonly Dictionary<PharmacophoreFeature, string> Definitions = new()
    {
        [PharmacophoreFeature.Donor] = "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
        [PharmacophoreFeature.Acceptor] = "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&X2&H0&+0,$([o,s;+0])]",
        [PharmacophoreFeature.PositiveIonisable] = "[$([NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[S,P]=O);!$(N#*);!$([N;R0]a)]),$([NX4+;!$(N~[O-])]),$([#7;+;!$([N+]~[O-])]),$(NC(=N)N)]",
        [PharmacophoreFeature.NegativeIonisable] = "[$([CX3](=O)[OX2H1,OX1-]),$([SX4](=O)(=O)[OX2H1,OX1-]),$([PX4](=O)([OX2H1,OX1-])),$([nH0-]1nnnc1),$(c1nn[nH]n1)]",
        [PharmacophoreFeature.Aromatic] = "[$(a1aaaaa1),$(a1aaaa1)]",
        [PharmacophoreFeature.Hydrophobic] = "[$([CX4;!$(C[N,O,S,F,Cl,Br,I])]),$([F,Cl,Br,I;$(*[CX4])]),$([S;X2;$(S([#6])[#6])])]",
        [PharmacophoreFeature.HalogenDonor] = "[Cl,Br,I;$(*[#6])]",
    };

    private static readonly Dictionary<PharmacophoreFeature, Lazy<SmartsPattern>> Patterns =
        Definitions.ToDictionary(kv => kv.Key, kv => new Lazy<SmartsPattern>(() => SmartsPattern.Parse(kv.Value)));

    /// <summary>The SMARTS used for <paramref name="feature"/>.</summary>
    public static string Definition(PharmacophoreFeature feature) => Definitions[feature];

    /// <summary>Atoms of <paramref name="mol"/> carrying <paramref name="feature"/>.</summary>
    public static IReadOnlyList<int> FindAtoms(Molecule mol, PharmacophoreFeature feature)
        => Patterns[feature].Value.MatchingAtoms(mol).ToList();

    /// <summary>Number of occurrences of <paramref name="feature"/> (aromatic features are counted per ring).</summary>
    public static int Count(Molecule mol, PharmacophoreFeature feature)
    {
        if (feature == PharmacophoreFeature.Aromatic)
            return mol.Rings.BondRings.Count(r => r.All(b => mol.Bonds[b].IsAromatic));
        return FindAtoms(mol, feature).Count;
    }

    /// <summary>Counts of every feature type.</summary>
    public static IReadOnlyDictionary<PharmacophoreFeature, int> Counts(Molecule mol)
        => Enum.GetValues<PharmacophoreFeature>().ToDictionary(f => f, f => Count(mol, f));

    /// <summary>The feature counts as a vector in <see cref="Enum.GetValues{T}"/> order (for QSAR matrices).</summary>
    public static double[] ToVector(Molecule mol) => Enum.GetValues<PharmacophoreFeature>().Select(f => (double)Count(mol, f)).ToArray();

    /// <summary>
    /// Similarity of two molecules' feature profiles: the count-vector Tanimoto Σmin / Σmax, in [0, 1].
    /// A crude but useful 3-D-free pharmacophore comparison.
    /// </summary>
    public static double ProfileSimilarity(Molecule a, Molecule b)
    {
        var x = ToVector(a);
        var y = ToVector(b);
        double mins = 0, maxs = 0;
        for (int i = 0; i < x.Length; i++)
        {
            mins += Math.Min(x[i], y[i]);
            maxs += Math.Max(x[i], y[i]);
        }
        return maxs == 0 ? 1 : mins / maxs;
    }
}
