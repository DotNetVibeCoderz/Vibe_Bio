using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology;
using Biotech.Net.SyntheticBiology.Cloning;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Restriction enzyme geometry, digestion and religation.
/// </summary>
/// <remarks>
/// Cut positions follow REBASE: a "cut position" <c>p</c> means the bond between base <c>p-1</c> and
/// <c>p</c> on the top strand. The cases below are hand-placed sites whose cuts can be counted off by eye.
/// </remarks>
public class RestrictionTests
{
    [Fact]
    public void EcoRI_CutsBetweenGAndAAttcLeavingAFourBaseOverhang()
    {
        //        0123456789...
        //        AAAAGAATTCTTT   site at index 4, G^AATTC
        var site = Assert.Single(RestrictionEnzyme.FromName("EcoRI").Find("AAAAGAATTCTTT"));

        Assert.Equal(4, site.Position);
        Assert.Equal(Strand.Forward, site.Strand);
        Assert.Equal(5, site.TopCutPosition);
        Assert.Equal(9, site.BottomCutPosition);
        Assert.Equal("AATT", site.Overhang);
        Assert.Equal(4, site.OverhangLength);
    }

    [Theory]
    [InlineData("EcoRI", "GAATTC", 1, 5, 4, "5'")]
    [InlineData("BamHI", "GGATCC", 1, 5, 4, "5'")]
    [InlineData("HindIII", "AAGCTT", 1, 5, 4, "5'")]
    [InlineData("NotI", "GCGGCCGC", 2, 6, 4, "5'")]
    [InlineData("NdeI", "CATATG", 2, 4, 2, "5'")]
    [InlineData("PstI", "CTGCAG", 5, 1, -4, "3'")]
    [InlineData("KpnI", "GGTACC", 5, 1, -4, "3'")]
    [InlineData("SmaI", "CCCGGG", 3, 3, 0, "blunt")]
    [InlineData("EcoRV", "GATATC", 3, 3, 0, "blunt")]
    public void TypeIIp_GeometryMatchesRebase(string name, string recognitionSite, int topOffset, int bottomOffset,
        int overhangLength, string overhangKind)
    {
        var enzyme = RestrictionEnzyme.FromName(name);

        Assert.Equal(recognitionSite, enzyme.RecognitionSite);
        Assert.True(enzyme.IsPalindromic);
        Assert.False(enzyme.IsTypeIis);
        Assert.Equal(topOffset, enzyme.TopCutOffsetFromSiteStart);
        Assert.Equal(bottomOffset, enzyme.BottomCutOffsetFromSiteStart);
        Assert.Equal(overhangLength, enzyme.OverhangLength);
        Assert.Equal(overhangKind, enzyme.OverhangKind);
        Assert.Equal(overhangLength == 0, enzyme.IsBlunt);

        // A palindromic site is reported once, on the forward strand.
        string dna = "AAAAA" + recognitionSite + "TTTTT";
        var site = Assert.Single(enzyme.Find(dna));
        Assert.Equal(5, site.Position);
        Assert.Equal(5 + topOffset, site.TopCutPosition);
        Assert.Equal(5 + bottomOffset, site.BottomCutPosition);
    }

    [Fact]
    public void PstI_LeavesAThreePrimeOverhang()
    {
        // CTGCA^G: the top strand is cut after 5 bases, the bottom after 1, so the overhang is on the 3' side.
        var site = Assert.Single(RestrictionEnzyme.FromName("PstI").Find("AAAAACTGCAGTTTTT"));

        Assert.Equal(5, site.Position);
        Assert.Equal(10, site.TopCutPosition);
        Assert.Equal(6, site.BottomCutPosition);
        Assert.Equal(-4, site.OverhangLength);
        Assert.Equal("TGCA", site.Overhang);
    }

    [Fact]
    public void BsaI_CutsOutsideItsSiteAndIsFoundOnBothStrands()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        Assert.True(bsaI.IsTypeIis);
        Assert.False(bsaI.IsPalindromic);
        Assert.Equal(4, bsaI.OverhangLength);
        // GGTCTC(1/5): 1 nt past the 6 bp site on top, 5 nt past on the bottom.
        Assert.Equal(1, bsaI.TopCutOffsetFromSiteEnd);
        Assert.Equal(5, bsaI.BottomCutOffsetFromSiteEnd);
        Assert.Equal(7, bsaI.TopCutOffsetFromSiteStart);
        Assert.Equal(11, bsaI.BottomCutOffsetFromSiteStart);

        //        0....5....10...
        //        AAAAAGGTCTCNCATGAAAA  site starts at 5, top cut at 5+7=12, bottom at 5+11=16
        var forward = Assert.Single(bsaI.Find("AAAAAGGTCTCACATGAAAA"));
        Assert.Equal(Strand.Forward, forward.Strand);
        Assert.Equal(5, forward.Position);
        Assert.Equal(12, forward.TopCutPosition);
        Assert.Equal(16, forward.BottomCutPosition);
        Assert.Equal("CATG", forward.Overhang);

        // The same construct read from the other strand: the enzyme now cuts on the opposite side of its site.
        string flipped = DnaUtility.ReverseComplement("AAAAAGGTCTCACATGAAAA");
        var reverse = Assert.Single(bsaI.Find(flipped));
        Assert.Equal(Strand.Reverse, reverse.Strand);
        Assert.Equal(flipped.Length - 16, reverse.TopCutPosition);
        Assert.Equal(flipped.Length - 12, reverse.BottomCutPosition);
        Assert.Equal("CATG", reverse.Overhang);
    }

    [Fact]
    public void BstXI_MatchesAnInterruptedAmbiguousSite()
    {
        var bstXI = RestrictionEnzyme.FromName("BstXI");
        Assert.Equal("CCANNNNNNTGG", bstXI.RecognitionSite);

        // CCANNNNN^NTGG: the N run is matched by anything.
        var site = Assert.Single(bstXI.Find("TTTTCCAGCTAGCTGGAAAA"));
        Assert.Equal(4, site.Position);
        Assert.Equal(4 + 8, site.TopCutPosition);
        Assert.Equal(4 + 4, site.BottomCutPosition);
        Assert.Equal(-4, site.OverhangLength);

        Assert.Empty(bstXI.Find("TTTTCCAGCTAGCTGAAAAA"));
    }

    [Fact]
    public void Find_CountsEveryOccurrenceAndNoneWhenAbsent()
    {
        var ecoRI = RestrictionEnzyme.FromName("EcoRI");
        Assert.Equal(2, ecoRI.CountSites("GAATTCAAAAAGAATTC"));
        Assert.Equal(0, ecoRI.CountSites("GAATTAAAAAAGAATTA"));
        Assert.Equal(0, ecoRI.CountSites(""));
        // Case and U are normalised.
        Assert.Equal(1, ecoRI.CountSites("aaagaattcaaa"));
    }

    [Fact]
    public void Find_DetectsOriginSpanningSitesOnlyWhenCircular()
    {
        // The site straddles the origin: the last three bases plus the first three.
        const string circle = "TTCAAAAAAAAAAGAA";
        var ecoRI = RestrictionEnzyme.FromName("EcoRI");

        Assert.Empty(ecoRI.Find(circle, circular: false));
        var site = Assert.Single(ecoRI.Find(circle, circular: true));
        Assert.Equal(13, site.Position);
        Assert.Equal("AATT", site.Overhang);
    }

    [Fact]
    public void Catalogue_IsConsistentAndLookupWorksBothWays()
    {
        foreach (var enzyme in RestrictionEnzyme.All)
        {
            Assert.True(RestrictionEnzyme.TryFromName(enzyme.Name, out var found));
            Assert.Equal(enzyme, found);
            // The two offset readings always agree.
            Assert.Equal(enzyme.TopCutOffsetFromSiteStart - enzyme.SiteLength, enzyme.TopCutOffsetFromSiteEnd);
            Assert.Equal(enzyme.BottomCutOffsetFromSiteStart - enzyme.SiteLength, enzyme.BottomCutOffsetFromSiteEnd);
            Assert.Equal(enzyme.BottomCutOffsetFromSiteStart - enzyme.TopCutOffsetFromSiteStart, enzyme.OverhangLength);
        }
        Assert.False(RestrictionEnzyme.TryFromName("NotAnEnzyme", out _));
        Assert.Throws<KeyNotFoundException>(() => RestrictionEnzyme.FromName("NotAnEnzyme"));
        // Esp3I is BsmBI under its other name.
        Assert.Equal(RestrictionEnzyme.FromName("BsmBI").RecognitionSite, RestrictionEnzyme.FromName("Esp3I").RecognitionSite);
    }

    // ---- digest and religate ----

    private static Plasmid CircleWithSites(string site, int copies, int spacerLength = 30)
    {
        // Spacer chosen so it carries no site of its own.
        string spacer = string.Concat(Enumerable.Repeat("ATATAT", (spacerLength + 5) / 6))[..spacerLength];
        return new Plasmid("test", string.Concat(Enumerable.Repeat(spacer + site, copies)));
    }

    [Theory]
    [InlineData("GAATTC", "EcoRI")]   // palindromic 5' overhang
    [InlineData("GGATCC", "BamHI")]
    [InlineData("CTGCAG", "PstI")]    // 3' overhang
    [InlineData("CCCGGG", "SmaI")]    // blunt
    public void Digest_ThenReligate_ReproducesThePlasmid(string site, string enzymeName)
    {
        var enzyme = RestrictionEnzyme.FromName(enzymeName);
        var plasmid = CircleWithSites(site, copies: 3);

        var fragments = plasmid.Fragments(enzyme);

        Assert.Equal(3, fragments.Count);
        // Every adjacent pair, including the wrap-around, must be ligatable.
        for (int i = 0; i < fragments.Count; i++)
            Assert.True(fragments[i].CanLigateTo(fragments[(i + 1) % fragments.Count]),
                $"fragment {i} cannot ligate to its neighbour");

        var joined = fragments[0];
        for (int i = 1; i < fragments.Count; i++) joined = joined.Ligate(fragments[i]);
        var circle = joined.Circularise();

        Assert.True(circle.IsCircular);
        Assert.Equal(plasmid.Length, circle.Length);
        Assert.True(IsRotationOf(plasmid.Sequence, circle.Sequence),
            $"religated circle {circle.Sequence} is not a rotation of {plasmid.Sequence}");
    }

    [Fact]
    public void Digest_WithANonPalindromicTypeIisOverhangAlsoReligates()
    {
        // Three BsaI sites, each exposing a different 4 nt overhang, all pointing the same way round the circle.
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        const string spacer = "ATATATATATATATATATATATATATAT";
        string sequence = spacer + "GGTCTCAGGAG" + spacer + "GGTCTCACCTA" + spacer + "GGTCTCATTGC";
        var plasmid = new Plasmid("typeIIs", sequence);

        var fragments = plasmid.Fragments(bsaI);

        Assert.Equal(3, fragments.Count);
        Assert.All(fragments, f => Assert.Equal(4, f.FivePrimeOverhang.Length));
        for (int i = 0; i < fragments.Count; i++)
            Assert.True(fragments[i].CanLigateTo(fragments[(i + 1) % fragments.Count]));

        var joined = fragments[0];
        for (int i = 1; i < fragments.Count; i++) joined = joined.Ligate(fragments[i]);
        Assert.True(IsRotationOf(sequence, joined.Circularise().Sequence));
    }

    [Fact]
    public void Digest_WithNoSiteReturnsTheWholeCircle()
    {
        var plasmid = new Plasmid("uncut", string.Concat(Enumerable.Repeat("ATATAT", 20)));
        var fragment = Assert.Single(plasmid.Fragments(RestrictionEnzyme.FromName("EcoRI")));
        Assert.True(fragment.IsCircular);
        Assert.Equal(plasmid.Length, fragment.Length);
    }

    [Fact]
    public void UniqueCutters_FindsTheEnzymesThatOpenTheVectorOnce()
    {
        var plasmid = CircleWithSites("GAATTC", copies: 1);
        var cutters = plasmid.UniqueCutters();

        Assert.Contains(cutters, e => e.Name == "EcoRI");
        Assert.DoesNotContain(cutters, e => e.Name == "BamHI");

        // Two copies of the site: no longer a unique cutter.
        Assert.DoesNotContain(CircleWithSites("GAATTC", copies: 2).UniqueCutters(), e => e.Name == "EcoRI");
    }

    [Fact]
    public void Fragment_LigationRejectsIncompatibleEnds()
    {
        var ecoRI = new DnaFragment("a", "AATTCGGG", "AATT", "");
        var bamHI = new DnaFragment("b", "GATCCGGG", "GATC", "");
        var blunt = new DnaFragment("c", "GGGGGG");

        Assert.False(blunt.CanLigateTo(ecoRI));           // blunt never joins sticky
        Assert.True(blunt.CanLigateTo(new DnaFragment("d", "AAAA")));
        Assert.False(ecoRI.CanLigateTo(bamHI));
        Assert.Throws<InvalidOperationException>(() => ecoRI.Ligate(bamHI));
        Assert.Throws<InvalidOperationException>(() => blunt.Circularise().Ligate(blunt));
    }

    [Fact]
    public void Fragment_ReverseComplementIsAnInvolution()
    {
        var fragment = new DnaFragment("f", "AATTCGGGTACC", "AATT", "TACC");
        var twice = fragment.ReverseComplement().ReverseComplement();

        Assert.Equal(fragment.Sequence, twice.Sequence);
        Assert.Equal(fragment.FivePrimeOverhang, twice.FivePrimeOverhang);
        Assert.Equal(fragment.ThreePrimeOverhang, twice.ThreePrimeOverhang);
        Assert.Equal("CGGG", fragment.Duplex);
    }

    private static bool IsRotationOf(string a, string b) =>
        a.Length == b.Length && (a + a).Contains(b, StringComparison.Ordinal);
}
