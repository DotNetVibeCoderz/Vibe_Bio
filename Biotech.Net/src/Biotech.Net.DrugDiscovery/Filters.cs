using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>Why a compound was rejected (empty when it passed).</summary>
/// <param name="Id">Compound id.</param>
/// <param name="Passed">Whether every criterion was met.</param>
/// <param name="Reasons">Failed criteria.</param>
public sealed record FilterResult(string Id, bool Passed, IReadOnlyList<string> Reasons)
{
    /// <inheritdoc/>
    public override string ToString() => Passed ? $"{Id}: pass" : $"{Id}: reject ({string.Join("; ", Reasons)})";
}

/// <summary>
/// A composable library filter: numeric property windows, required and forbidden substructures, and alert
/// catalogues. Criteria are added fluently and evaluated in the order they were added, so cheap property checks
/// should come first.
/// </summary>
/// <example>
/// <code>
/// var filter = new CompoundFilter()
///     .WithPropertyRange("MW", m => MolecularDescriptors.MolecularWeight(m), 150, 500)
///     .RequireSubstructure("c1ccccc1", "benzene ring")
///     .ExcludeAlerts(StructuralAlerts.Brenk);
/// var kept = filter.Apply(library);
/// </code>
/// </example>
public sealed class CompoundFilter
{
    private readonly List<(string Name, Func<Molecule, string?> Check)> _criteria = [];

    /// <summary>Names of the criteria, in evaluation order.</summary>
    public IReadOnlyList<string> Criteria => _criteria.Select(c => c.Name).ToList();

    /// <summary>Adds an arbitrary criterion; the function returns null to pass or a reason to reject.</summary>
    public CompoundFilter Where(string name, Func<Molecule, string?> check)
    {
        _criteria.Add((name, check));
        return this;
    }

    /// <summary>Requires <paramref name="property"/> to lie in [<paramref name="min"/>, <paramref name="max"/>].</summary>
    public CompoundFilter WithPropertyRange(string name, Func<Molecule, double> property, double min = double.NegativeInfinity, double max = double.PositiveInfinity)
        => Where(name, m =>
        {
            double v = property(m);
            return v < min || v > max ? $"{name} {v:F2} outside [{min:F2}, {max:F2}]" : null;
        });

    /// <summary>Requires a SMARTS match.</summary>
    public CompoundFilter RequireSubstructure(string smarts, string? name = null)
    {
        var p = SmartsPattern.Parse(smarts);
        return Where(name ?? $"requires {smarts}", m => p.IsMatch(m) ? null : $"missing {name ?? smarts}");
    }

    /// <summary>Rejects compounds that contain a SMARTS match.</summary>
    public CompoundFilter ExcludeSubstructure(string smarts, string? name = null)
    {
        var p = SmartsPattern.Parse(smarts);
        return Where(name ?? $"excludes {smarts}", m => p.IsMatch(m) ? $"contains {name ?? smarts}" : null);
    }

    /// <summary>Rejects compounds matching any alert of <paramref name="catalogue"/>.</summary>
    public CompoundFilter ExcludeAlerts(IReadOnlyList<StructuralAlert>? catalogue = null)
    {
        var list = catalogue ?? StructuralAlerts.All;
        return Where("structural alerts", m =>
        {
            var hits = StructuralAlerts.Match(m, list);
            return hits.Count == 0 ? null : $"alerts: {string.Join(", ", hits.Select(a => $"{a.Family}/{a.Name}"))}";
        });
    }

    /// <summary>Requires a drug-likeness rule to pass.</summary>
    public CompoundFilter RequireRule(Func<DescriptorSet, RuleResult> rule)
        => Where("rule", m =>
        {
            var r = rule(MolecularDescriptors.Compute(m));
            return r.Passed ? null : r.ToString();
        });

    /// <summary>Lipinski's rule of five with at most <paramref name="allowedViolations"/> violations.</summary>
    public CompoundFilter RequireLipinski(int allowedViolations = 1) => RequireRule(d => DrugLikeness.Lipinski(d, allowedViolations));

    /// <summary>
    /// The usual lead-like property window: MW 200–350, logP −1–3.5, HBD ≤ 3, HBA ≤ 7, rotatable bonds ≤ 7,
    /// TPSA ≤ 120 Å².
    /// </summary>
    public static CompoundFilter LeadLike() => new CompoundFilter()
        .WithPropertyRange("MW", MolecularDescriptors.MolecularWeight, 200, 350)
        .WithPropertyRange("logP", Crippen.LogP, -1, 3.5)
        .WithPropertyRange("HBD", m => MolecularDescriptors.HBondDonors(m), max: 3)
        .WithPropertyRange("HBA", m => MolecularDescriptors.HBondAcceptors(m), max: 7)
        .WithPropertyRange("RotB", m => MolecularDescriptors.RotatableBonds(m), max: 7)
        .WithPropertyRange("TPSA", Tpsa.Compute, max: 120);

    /// <summary>
    /// The "rule of three" fragment window (Congreve et al., Drug Discov. Today 2003): MW ≤ 300, logP ≤ 3,
    /// HBD ≤ 3, HBA ≤ 3, rotatable bonds ≤ 3.
    /// </summary>
    public static CompoundFilter FragmentLike() => new CompoundFilter()
        .WithPropertyRange("MW", MolecularDescriptors.MolecularWeight, max: 300)
        .WithPropertyRange("logP", Crippen.LogP, max: 3)
        .WithPropertyRange("HBD", m => MolecularDescriptors.HBondDonors(m), max: 3)
        .WithPropertyRange("HBA", m => MolecularDescriptors.HBondAcceptors(m), max: 3)
        .WithPropertyRange("RotB", m => MolecularDescriptors.RotatableBonds(m), max: 3);

    /// <summary>Evaluates every criterion against one molecule.</summary>
    public FilterResult Evaluate(Molecule mol, string? id = null)
    {
        var reasons = new List<string>();
        foreach (var (_, check) in _criteria)
        {
            string? why = check(mol);
            if (why is not null) reasons.Add(why);
        }
        return new FilterResult(id ?? mol.Name, reasons.Count == 0, reasons);
    }

    /// <summary>Whether a molecule passes every criterion (stops at the first failure).</summary>
    public bool Passes(Molecule mol) => _criteria.All(c => c.Check(mol) is null);

    /// <summary>Evaluates a library in parallel, keeping the input order.</summary>
    public IReadOnlyList<FilterResult> Evaluate(IReadOnlyList<LibraryCompound> library)
    {
        var r = new FilterResult[library.Count];
        Parallel.For(0, library.Count, i => r[i] = Evaluate(library[i].Molecule, library[i].Id));
        return r;
    }

    /// <summary>The compounds of <paramref name="library"/> that pass.</summary>
    public IReadOnlyList<LibraryCompound> Apply(IReadOnlyList<LibraryCompound> library)
    {
        var results = Evaluate(library);
        return Enumerable.Range(0, library.Count).Where(i => results[i].Passed).Select(i => library[i]).ToList();
    }
}
