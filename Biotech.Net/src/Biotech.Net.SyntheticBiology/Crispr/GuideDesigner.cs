using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>What a designed guide should be ranked on, and how hard each filter bites.</summary>
public sealed class GuideDesignOptions
{
    /// <summary>How candidate guides are found and which composition filters they must pass.</summary>
    public GuideSearchOptions Search { get; init; } = GuideSearchOptions.Default;

    /// <summary>
    /// The genome (or any set of sequences) to profile off-targets against. When null, no off-target search is
    /// run and the specificity terms drop out of the composite score.
    /// </summary>
    public CrisprSiteIndex? OffTargetIndex { get; init; }

    /// <summary>Limits for the off-target search, when <see cref="OffTargetIndex"/> is given.</summary>
    public OffTargetSearchOptions OffTarget { get; init; } = OffTargetSearchOptions.Default;

    /// <summary>Weight of the Doench 2014 on-target score in the composite (default 0.5).</summary>
    public double OnTargetWeight { get; init; } = 0.5;

    /// <summary>Weight of the CFD aggregate specificity in the composite (default 0.35).</summary>
    public double CfdWeight { get; init; } = 0.35;

    /// <summary>Weight of the MIT aggregate specificity in the composite (default 0.15).</summary>
    public double MitWeight { get; init; } = 0.15;

    /// <summary>Discard guides whose CFD specificity is below this (0–100). 0 disables the filter.</summary>
    public double MinCfdSpecificity { get; init; }

    /// <summary>Discard guides with any off-target of 0, 1 or 2 mismatches. Off by default.</summary>
    public bool RejectCloseOffTargets { get; init; }

    /// <summary>Maximum number of designs to return (0 = all).</summary>
    public int MaxResults { get; init; } = 20;

    /// <summary>The default options: SpCas9, no off-target index, the published weighting above.</summary>
    public static GuideDesignOptions Default { get; } = new();
}

/// <summary>
/// One designed guide with everything that went into its ranking.
/// </summary>
/// <param name="Guide">The guide itself.</param>
/// <param name="OnTargetScore">Doench 2014 Rule Set 1 (0–1), <see cref="double.NaN"/> when the model does not apply.</param>
/// <param name="Penalties">The composition penalties applied on top of the model.</param>
/// <param name="OffTargets">The off-target profile, or null when no index was supplied.</param>
/// <param name="CompositeScore">The 0–1 ranking score described on <see cref="GuideDesigner"/>.</param>
public sealed record GuideDesign(
    Guide Guide,
    double OnTargetScore,
    GuidePenalties Penalties,
    OffTargetProfile? OffTargets,
    double CompositeScore)
{
    /// <summary>CFD aggregate specificity (0–100), or <see cref="double.NaN"/> without an off-target search.</summary>
    public double CfdSpecificity => OffTargets?.CfdSpecificity ?? double.NaN;

    /// <summary>MIT aggregate specificity (0–100), or <see cref="double.NaN"/> without an off-target search.</summary>
    public double MitSpecificity => OffTargets?.MitSpecificity ?? double.NaN;

    /// <summary>Number of off-target sites found (the intended site excluded).</summary>
    public int OffTargetCount => OffTargets is null ? 0 : OffTargets.Hits.Count(static h => !h.IsOnTarget);

    /// <inheritdoc/>
    public override string ToString()
        => $"{Guide.Spacer} {Guide.Pam} score {CompositeScore:F3}"
           + (double.IsNaN(OnTargetScore) ? "" : $" on-target {OnTargetScore:F3}")
           + (OffTargets is null ? "" : $" CFD-spec {CfdSpecificity:F1} ({OffTargetCount} off-targets)");
}

/// <summary>
/// The one-call guide-design façade: find candidate spacers, score them for on-target activity (Doench 2014),
/// profile them against a genome (CFD and MIT) and rank the result.
/// </summary>
/// <remarks>
/// <para>
/// The composite score is a weighted mean of three 0–1 terms — the Rule Set 1 on-target score multiplied by the
/// composition penalties of <see cref="OnTargetScorer.Penalties(string)"/>, the CFD aggregate specificity
/// divided by 100, and the MIT aggregate specificity divided by 100 — using
/// <see cref="GuideDesignOptions.OnTargetWeight"/>, <see cref="GuideDesignOptions.CfdWeight"/> and
/// <see cref="GuideDesignOptions.MitWeight"/>. Terms that cannot be computed (no off-target index, or a
/// nuclease the Doench model was not trained on) are dropped and the remaining weights are renormalised, so the
/// score stays in 0–1 and stays comparable within one run.
/// </para>
/// <para>
/// The weighting is a convention of this library, not a published model: the three underlying scores are
/// published, the way they are blended is not. Ranking guides against each other is what it is for; the
/// absolute number means nothing.
/// </para>
/// </remarks>
public static class GuideDesigner
{
    /// <summary>Designs and ranks guides over a whole target sequence.</summary>
    /// <param name="target">The sequence to design against.</param>
    /// <param name="options">Search, scoring and ranking options.</param>
    /// <returns>Designs ordered by descending <see cref="GuideDesign.CompositeScore"/>.</returns>
    public static IReadOnlyList<GuideDesign> Design(Sequence target, GuideDesignOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        options ??= GuideDesignOptions.Default;
        return Rank(GuideFinder.Find(target, options.Search), options);
    }

    /// <summary>Designs and ranks guides over a plain DNA string.</summary>
    /// <param name="target">The sequence to design against.</param>
    /// <param name="options">Search, scoring and ranking options.</param>
    /// <param name="sequenceId">Identifier reported on the returned guides.</param>
    /// <returns>Designs ordered by descending <see cref="GuideDesign.CompositeScore"/>.</returns>
    public static IReadOnlyList<GuideDesign> Design(string target, GuideDesignOptions? options = null, string sequenceId = "target")
    {
        ArgumentNullException.ThrowIfNull(target);
        options ??= GuideDesignOptions.Default;
        return Rank(GuideFinder.Find(target, options.Search, sequenceId), options);
    }

    /// <summary>
    /// Designs guides that cut inside a region, the usual knock-out case: the region is the exon (or the
    /// functional domain) you want to disrupt.
    /// </summary>
    /// <param name="target">The sequence to design against.</param>
    /// <param name="start">0-based start of the region the cut must fall in.</param>
    /// <param name="length">Length of that region.</param>
    /// <param name="options">Search, scoring and ranking options.</param>
    /// <returns>Designs ordered by descending <see cref="GuideDesign.CompositeScore"/>.</returns>
    public static IReadOnlyList<GuideDesign> DesignInRegion(Sequence target, int start, int length, GuideDesignOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        options ??= GuideDesignOptions.Default;
        return Rank(GuideFinder.FindInRegion(target, start, length, options.Search), options);
    }

    /// <summary>Scores and ranks guides that have already been found.</summary>
    /// <param name="guides">The candidate guides.</param>
    /// <param name="options">Scoring and ranking options.</param>
    /// <returns>Designs ordered by descending <see cref="GuideDesign.CompositeScore"/>.</returns>
    public static IReadOnlyList<GuideDesign> Rank(IEnumerable<Guide> guides, GuideDesignOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(guides);
        options ??= GuideDesignOptions.Default;

        var designs = new List<GuideDesign>();
        foreach (var guide in guides)
        {
            double onTarget = double.IsNaN(guide.OnTargetScore) ? OnTargetScorer.Score(guide) : guide.OnTargetScore;
            var penalties = OnTargetScorer.Penalties(guide.Spacer);

            OffTargetProfile? profile = null;
            if (options.OffTargetIndex is { } index && guide.Spacer.Length == index.Nuclease.GuideLength)
            {
                profile = OffTargetSearch.Search(guide, index, options.OffTarget);
                if (options.MinCfdSpecificity > 0 && profile.CfdSpecificity < options.MinCfdSpecificity) continue;
                if (options.RejectCloseOffTargets && profile.OffTargets.Any(static h => h.MismatchCount <= 2)) continue;
            }

            designs.Add(new GuideDesign(guide, onTarget, penalties, profile,
                Composite(onTarget, penalties, profile, options)));
        }

        designs.Sort(static (a, b) => b.CompositeScore.CompareTo(a.CompositeScore));
        if (options.MaxResults > 0 && designs.Count > options.MaxResults)
            designs.RemoveRange(options.MaxResults, designs.Count - options.MaxResults);
        return designs;
    }

    /// <summary>
    /// The weighted mean described on the class, over whichever terms are available.
    /// </summary>
    /// <param name="onTargetScore">Rule Set 1 score, or NaN.</param>
    /// <param name="penalties">Composition penalties for the spacer.</param>
    /// <param name="profile">Off-target profile, or null.</param>
    /// <param name="options">The weights to use.</param>
    /// <returns>A 0–1 ranking score.</returns>
    public static double Composite(double onTargetScore, GuidePenalties penalties, OffTargetProfile? profile,
        GuideDesignOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        double total = 0, weight = 0;

        if (!double.IsNaN(onTargetScore) && options.OnTargetWeight > 0)
        {
            total += options.OnTargetWeight * Math.Clamp(onTargetScore * penalties.Multiplier, 0, 1);
            weight += options.OnTargetWeight;
        }
        if (profile is not null && options.CfdWeight > 0)
        {
            total += options.CfdWeight * Math.Clamp(profile.CfdSpecificity / 100.0, 0, 1);
            weight += options.CfdWeight;
        }
        if (profile is not null && options.MitWeight > 0)
        {
            total += options.MitWeight * Math.Clamp(profile.MitSpecificity / 100.0, 0, 1);
            weight += options.MitWeight;
        }

        // Nothing scoreable: fall back on the composition penalties alone, against a neutral 0.5 baseline.
        if (weight <= 0) return Math.Clamp(0.5 * penalties.Multiplier, 0, 1);
        return Math.Clamp(total / weight, 0, 1);
    }
}
