using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology;
using Biotech.Net.SyntheticBiology.Crispr;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Guide discovery, the 2-bit site index, off-target search and the ranking façade.
/// </summary>
public class GuideDesignTests
{
    // 4 nt upstream + a 20 nt spacer + TGG PAM + 3 nt downstream: exactly one forward SpCas9 site, and no
    // reverse site because the plus strand contains no CC dinucleotide.
    private const string Spacer = "ACGTACGTACGTACGTACGT";
    private const string Target = "TTTT" + Spacer + "TGG" + "AAA";

    [Fact]
    public void Find_LocatesTheSingleForwardSite()
    {
        var guides = GuideFinder.Find(Target);

        var guide = Assert.Single(guides);
        Assert.Equal(Spacer, guide.Spacer);
        Assert.Equal("TGG", guide.Pam);
        Assert.Equal(Strand.Forward, guide.Strand);
        Assert.Equal(4, guide.Start);
        Assert.Equal(24, guide.End);
        Assert.Equal(Target, guide.Context);
        Assert.Equal(30, guide.Context!.Length);
    }

    [Fact]
    public void Find_ReportsTheReverseStrandSiteInSpacerCoordinates()
    {
        // The reverse complement of the whole construct: the same guide, now on the minus strand.
        string target = DnaUtility.ReverseComplement(Target);
        var guide = Assert.Single(GuideFinder.Find(target));

        Assert.Equal(Strand.Reverse, guide.Strand);
        Assert.Equal(Spacer, guide.Spacer);
        Assert.Equal("TGG", guide.Pam);
        // The protospacer occupies the mirror window of the 30 nt construct.
        Assert.Equal(target.Length - 24, guide.Start);
        Assert.Equal(DnaUtility.ReverseComplement(target.AsSpan(guide.Start, 20).ToString()), guide.Spacer);
    }

    [Fact]
    public void CutSite_IsThreeBasesFromThePamOnBothStrands()
    {
        var forward = Assert.Single(GuideFinder.Find(Target));
        // Blunt cut 3 bp into the protospacer from the PAM end: 4 + 20 - 3.
        Assert.Equal(21, forward.CutSite);

        var reverse = Assert.Single(GuideFinder.Find(DnaUtility.ReverseComplement(Target)));
        Assert.Equal(reverse.Start + 3, reverse.CutSite);
    }

    [Fact]
    public void FullSite_And_U6Oligo_FollowTheNucleaseGeometry()
    {
        var guide = Assert.Single(GuideFinder.Find(Target));
        Assert.Equal(Spacer + "TGG", guide.FullSite);
        Assert.Equal("G" + Spacer, guide.U6CloningOligo); // the spacer starts with A, so a G is prepended
        Assert.Equal(0.5, guide.GcContent, 12);
        Assert.False(guide.HasPolyTTerminator);
    }

    [Fact]
    public void Options_FilterOnComposition()
    {
        // GC filter: the 50% spacer fails a 60% floor.
        Assert.Empty(GuideFinder.Find(Target, new GuideSearchOptions { MinGcContent = 0.60 }));
        Assert.Single(GuideFinder.Find(Target, new GuideSearchOptions { MinGcContent = 0.40, MaxGcContent = 0.60 }));

        // A leading-G requirement rejects this spacer, which starts with A.
        Assert.Empty(GuideFinder.Find(Target, new GuideSearchOptions { RequireLeadingG = true }));

        // Strand switches.
        Assert.Empty(GuideFinder.Find(Target, new GuideSearchOptions { SearchForward = false }));
        Assert.Single(GuideFinder.Find(Target, new GuideSearchOptions { SearchReverse = false }));
    }

    [Fact]
    public void Options_RejectPolyTTerminators()
    {
        const string polyT = "ACGTTTTACGTACGTACGTA";
        string target = "AAAA" + polyT + "TGG" + "AAA";

        Assert.Empty(GuideFinder.Find(target, new GuideSearchOptions { RejectPolyT = true }));
        var kept = Assert.Single(GuideFinder.Find(target, new GuideSearchOptions { RejectPolyT = false }));
        Assert.True(kept.HasPolyTTerminator);
    }

    [Fact]
    public void Find_AmbiguityCodesAreSkipped()
    {
        string target = "TTTT" + "ACGTACGTNCGTACGTACGT" + "TGG" + "AAA";
        Assert.Empty(GuideFinder.Find(target));
        Assert.Single(GuideFinder.Find(target, new GuideSearchOptions { RequireUnambiguous = false, ScoreOnTarget = false }));
    }

    [Fact]
    public void FindInRegion_KeepsOnlyCutsInsideTheRegion()
    {
        var sequence = Sequence.Dna(Target, "t");
        Assert.Single(GuideFinder.FindInRegion(sequence, 18, 6, GuideSearchOptions.Default));  // cut at 21
        Assert.Empty(GuideFinder.FindInRegion(sequence, 0, 10, GuideSearchOptions.Default));
    }

    [Fact]
    public void NagPamIsNotFoundByTheCanonicalNuclease()
    {
        string target = "TTTT" + Spacer + "TAG" + "AAA";
        Assert.Empty(GuideFinder.Find(target));
        Assert.Single(GuideFinder.Find(target, new GuideSearchOptions { Nuclease = CasNuclease.SpCas9Nag }));
    }

    [Fact]
    public void Cas12a_UsesAFivePrimePamAndItsOwnSpacerLength()
    {
        const string spacer23 = "ACGTACGTACGTACGTACGTACG";
        string target = "AAA" + "TTTA" + spacer23 + "AAA";
        var guide = Assert.Single(GuideFinder.Find(target, new GuideSearchOptions
        {
            Nuclease = CasNuclease.Cas12a,
            SearchReverse = false,
        }));

        Assert.Equal(spacer23, guide.Spacer);
        Assert.Equal("TTTA", guide.Pam);
        Assert.Equal(7, guide.Start);
        Assert.Equal(spacer23 + "", guide.Spacer);
        Assert.Equal("TTTA" + spacer23, guide.FullSite);
        // The Doench model is SpCas9-only, so there is no on-target score for Cas12a.
        Assert.True(double.IsNaN(OnTargetScorer.Score(guide)));
    }

    [Fact]
    public void CasNuclease_CatalogueIsSelfConsistent()
    {
        foreach (var nuclease in CasNuclease.All)
        {
            Assert.Equal(nuclease.GuideLength + nuclease.Pam.Length, nuclease.SiteLength);
            Assert.Same(nuclease, CasNuclease.FromName(nuclease.Name));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => CasNuclease.FromName("not a nuclease"));
        Assert.True(CasNuclease.Cas12a.IsStaggered);
        Assert.False(CasNuclease.SpCas9.IsStaggered);
    }

    // ---- site index ----

    /// <summary>A small "genome": the on-target site, a 2-mismatch site and a 4-mismatch site, all with NGG PAMs.</summary>
    private static string MiniGenome()
    {
        const string onTarget = "GAGTCCGAGCAGAAGAAGAA";          // the guide
        const string twoMismatch = "GAGTCCGAGCAGAAGAACAC";       // positions 18 and 20
        const string fourMismatch = "TAGTCCGAGCATAAGAACAC";      // four substitutions
        return "AAAAA" + onTarget + "TGG" + "TTTTTTTTTT"
               + twoMismatch + "AGG" + "TTTTTTTTTT"
               + fourMismatch + "CGG" + "AAAAA";
    }

    [Fact]
    public void Index_SeededQueryAgreesWithTheExhaustiveScan()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1", includeNagPam: false);
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";

        for (int mm = 0; mm <= 4; mm++)
        {
            var seeded = index.Query(spacer, mm).OrderBy(h => h.Start).ToList();
            var exhaustive = index.QueryExhaustive(spacer, mm).OrderBy(h => h.Start).ToList();
            Assert.Equal(exhaustive.Select(h => (h.Start, h.Strand, h.MismatchCount)),
                         seeded.Select(h => (h.Start, h.Strand, h.MismatchCount)));
        }
    }

    [Fact]
    public void Index_FindsTheOnTargetAndItsMismatchedNeighbours()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1", includeNagPam: false);
        var hits = index.Query("GAGTCCGAGCAGAAGAAGAA", 4);

        Assert.Equal([0, 2, 4], hits.Select(h => h.MismatchCount).OrderBy(x => x));
        var perfect = Assert.Single(hits, h => h.MismatchCount == 0);
        Assert.Equal(5, perfect.Start);
        Assert.Equal("TGG", perfect.Pam);
        Assert.Equal(1.0, perfect.CfdScore, 12);
        Assert.Equal(100.0, perfect.MitScore, 12);
    }

    [Fact]
    public void Index_RejectsTooManyMismatchesAndWrongLengths()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1");
        Assert.Throws<ArgumentOutOfRangeException>(() => index.Query("GAGTCCGAGCAGAAGAAGAA", 5));
        Assert.Throws<ArgumentException>(() => index.Query("ACGT", 2));
        Assert.Throws<ArgumentException>(() => index.Query("GAGTCCGAGCAGAAGAAGAN", 2));
    }

    [Fact]
    public void Index_IncludesNagSitesOnRequest()
    {
        string genome = "AAAAA" + "GAGTCCGAGCAGAAGAAGAA" + "TAG" + "AAAAA";
        Assert.Empty(CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false).Query("GAGTCCGAGCAGAAGAAGAA", 0));

        var withNag = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: true).Query("GAGTCCGAGCAGAAGAAGAA", 0);
        var hit = Assert.Single(withNag);
        Assert.Equal("TAG", hit.Pam);
        // CFD penalises the AG PAM exactly as the published table says.
        Assert.Equal(0.259259259, hit.CfdScore, 9);
    }

    [Fact]
    public void Index_FindsReverseStrandSites()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        string genome = "AAAAA" + DnaUtility.ReverseComplement(spacer + "TGG") + "AAAAA";
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);

        var hit = Assert.Single(index.Query(spacer, 0));
        Assert.Equal(Strand.Reverse, hit.Strand);
        Assert.Equal(spacer, hit.Protospacer);
        Assert.Equal("TGG", hit.Pam);
    }

    // ---- off-target search ----

    [Fact]
    public void Search_ScoresAndRanksTheWholeProfile()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1", includeNagPam: false);
        var profile = OffTargetSearch.Search("GAGTCCGAGCAGAAGAAGAA", index);

        Assert.Equal(3, profile.Hits.Count);
        Assert.Equal(2, profile.OffTargets.Count());
        Assert.True(profile.Hits[0].IsOnTarget);
        // Ranked by descending CFD after the on-target.
        Assert.True(profile.Hits[1].CfdScore >= profile.Hits[2].CfdScore);

        Assert.InRange(profile.CfdSpecificity, 0, 100);
        Assert.InRange(profile.MitSpecificity, 0, 100);

        // Specificity is exactly the published aggregate over the off-targets only.
        double sum = profile.OffTargets.Sum(h => h.CfdScore * 100);
        Assert.Equal(100.0 * 100.0 / (100 + sum), profile.CfdSpecificity, 9);

        var histogram = profile.MismatchHistogram(4);
        Assert.Equal(0, histogram[0]);
        Assert.Equal(1, histogram[2]);
        Assert.Equal(1, histogram[4]);
    }

    [Fact]
    public void Search_AGuideWithNoOffTargetsIsPerfectlySpecific()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        var index = CrisprSiteIndex.Build("AAAAA" + spacer + "TGG" + "AAAAA", "chr1", includeNagPam: false);
        var profile = OffTargetSearch.Search(spacer, index);

        Assert.Single(profile.Hits);
        Assert.Empty(profile.OffTargets);
        Assert.Equal(100.0, profile.CfdSpecificity, 9);
        Assert.Equal(100.0, profile.MitSpecificity, 9);
    }

    [Fact]
    public void SearchAt_MarksTheDesignedLocusRatherThanTheFirstPerfectMatch()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        // The same site twice: the designed one is the second copy.
        string genome = "AAAAA" + spacer + "TGG" + "TTTTTTTTTT" + spacer + "TGG" + "AAAAA";
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);

        int secondStart = 5 + 20 + 3 + 10;
        var profile = OffTargetSearch.SearchAt(spacer, index, "chr1", secondStart);

        var onTarget = Assert.Single(profile.Hits, h => h.IsOnTarget);
        Assert.Equal(secondStart, onTarget.Start);
        // The other perfect match counts fully against specificity: 100*100/(100+100).
        Assert.Equal(50.0, profile.CfdSpecificity, 9);
    }

    [Fact]
    public void Search_MaxResultsCapsTheOutput()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1", includeNagPam: false);
        var profile = OffTargetSearch.Search("GAGTCCGAGCAGAAGAAGAA", index,
            new OffTargetSearchOptions { MaxResults = 2 });
        Assert.Equal(2, profile.Hits.Count);
    }

    [Fact]
    public void Search_BulgePassFindsASiteWithAnExtraGenomicBase()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        // An extra T looped out of the DNA in the middle of an otherwise perfect site.
        string bulged = string.Concat(spacer.AsSpan(0, 10), "T", spacer.AsSpan(10));
        string genome = "AAAAA" + bulged + "TGG" + "AAAAA";
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);

        Assert.Empty(OffTargetSearch.Search(spacer, index).OffTargets);

        var withBulges = OffTargetSearch.Search(spacer, index, new OffTargetSearchOptions { MaxBulges = 1 });
        var hit = Assert.Single(withBulges.Hits, h => h.Bulge == BulgeKind.Dna);
        Assert.Equal(0, hit.MismatchCount);
        Assert.Equal(bulged, hit.Protospacer);
        // Bulge hits are penalised and kept out of the published aggregate scores.
        Assert.Equal(OffTargetSearch.BulgePenalty, hit.CfdScore, 9);
        Assert.Equal(100.0, withBulges.CfdSpecificity, 9);
    }

    [Fact]
    public void Search_RejectsMoreThanOneBulge()
    {
        var index = CrisprSiteIndex.Build(MiniGenome(), "chr1");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            OffTargetSearch.Search("GAGTCCGAGCAGAAGAAGAA", index, new OffTargetSearchOptions { MaxBulges = 2 }));
    }

    [Fact]
    public void BulgeAlignment_CollapsesTheLoopedOutBase()
    {
        const string spacer = "GAGTCCGAGCAGAAGAAGAA";
        string site = string.Concat(spacer.AsSpan(0, 5), "T", spacer.AsSpan(5));

        Assert.True(OffTargetSearch.TryBestBulgeAlignment(spacer, site, BulgeKind.Dna, 0, out int mismatches, out string collapsed));
        Assert.Equal(0, mismatches);
        Assert.Equal(spacer, collapsed);
    }

    [Fact]
    public void Alignment_MasksMatchingPositions()
    {
        var hit = new OffTarget("chr1", 0, Strand.Forward, "GAGTCCGAGCAGAAGAACAC", "TGG", 2, BulgeKind.None, 0.1, 5, false);
        // Mismatches at 0-based positions 17 and 19; everything else matches.
        Assert.Equal(".................C.C", hit.Alignment("GAGTCCGAGCAGAAGAAGAA"));
    }

    // ---- designer façade ----

    [Fact]
    public void Designer_RanksByTheCompositeScore()
    {
        var designs = GuideDesigner.Design(Target);
        var design = Assert.Single(designs);

        Assert.Equal(Spacer, design.Guide.Spacer);
        Assert.Equal(OnTargetScorer.RuleSet1(Target), design.OnTargetScore, 12);
        Assert.Null(design.OffTargets);
        Assert.True(double.IsNaN(design.CfdSpecificity));
        // Without an index only the on-target term survives, so the composite is that term alone.
        Assert.Equal(Math.Clamp(design.OnTargetScore * design.Penalties.Multiplier, 0, 1), design.CompositeScore, 12);
    }

    [Fact]
    public void Designer_AddsSpecificityWhenGivenAnIndex()
    {
        string genome = MiniGenome();
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);
        var options = new GuideDesignOptions
        {
            OffTargetIndex = index,
            Search = new GuideSearchOptions { MinGcContent = 0.0, MaxGcContent = 1.0 },
        };

        var designs = GuideDesigner.Design(genome, options, "chr1");

        Assert.NotEmpty(designs);
        Assert.All(designs, d => Assert.NotNull(d.OffTargets));
        Assert.All(designs, d => Assert.InRange(d.CompositeScore, 0, 1));
        // Sorted, descending.
        Assert.Equal(designs.Select(d => d.CompositeScore).OrderByDescending(x => x), designs.Select(d => d.CompositeScore));

        var onTargetGuide = designs.First(d => d.Guide.Spacer == "GAGTCCGAGCAGAAGAAGAA");
        Assert.Equal(2, onTargetGuide.OffTargetCount);
    }

    [Fact]
    public void Designer_CanRejectGuidesWithCloseOffTargets()
    {
        string genome = MiniGenome();
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);
        var options = new GuideDesignOptions
        {
            OffTargetIndex = index,
            RejectCloseOffTargets = true,
            Search = new GuideSearchOptions { MinGcContent = 0.0, MaxGcContent = 1.0 },
        };

        var designs = GuideDesigner.Design(genome, options, "chr1");
        Assert.DoesNotContain(designs, d => d.Guide.Spacer == "GAGTCCGAGCAGAAGAAGAA");
    }

    [Fact]
    public void Designer_MinCfdSpecificityFiltersAndMaxResultsCaps()
    {
        string genome = MiniGenome();
        var index = CrisprSiteIndex.Build(genome, "chr1", includeNagPam: false);
        var baseOptions = new GuideDesignOptions
        {
            OffTargetIndex = index,
            Search = new GuideSearchOptions { MinGcContent = 0.0, MaxGcContent = 1.0 },
        };

        var all = GuideDesigner.Design(genome, baseOptions, "chr1");
        var capped = GuideDesigner.Design(genome, new GuideDesignOptions
        {
            OffTargetIndex = index,
            MaxResults = 1,
            Search = baseOptions.Search,
        }, "chr1");

        Assert.True(all.Count >= capped.Count);
        Assert.Single(capped);

        var strict = GuideDesigner.Design(genome, new GuideDesignOptions
        {
            OffTargetIndex = index,
            MinCfdSpecificity = 99.9,
            Search = baseOptions.Search,
        }, "chr1");
        Assert.All(strict, d => Assert.True(d.CfdSpecificity >= 99.9));
    }

    [Fact]
    public void Composite_RenormalisesWhenTermsAreMissing()
    {
        var penalties = new GuidePenalties(0, 0, 0, 0);
        var options = GuideDesignOptions.Default;

        // Only the on-target term is available, so it is the whole score.
        Assert.Equal(0.8, GuideDesigner.Composite(0.8, penalties, null, options), 12);

        // Nothing at all: the neutral baseline reduced by the penalties.
        Assert.Equal(0.5, GuideDesigner.Composite(double.NaN, penalties, null, options), 12);
        Assert.Equal(0.25, GuideDesigner.Composite(double.NaN, new GuidePenalties(0.5, 0, 0, 0), null, options), 12);
    }
}
