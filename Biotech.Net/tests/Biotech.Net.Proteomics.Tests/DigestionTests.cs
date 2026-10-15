using Biotech.Net.Core;
using Biotech.Net.Proteomics.MassSpec;

namespace Biotech.Net.Proteomics.Tests;

public class DigestionTests
{
    private const string Ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    private static readonly DigestionOptions All = new(MissedCleavages: 0, MinLength: 1, MaxLength: 1000, ClipInitiatorMethionine: false);

    [Fact]
    public void TrypticDigestOfUbiquitin()
    {
        // Fully tryptic peptides of ubiquitin (e.g. Expasy PeptideCutter / PeptideMass, 0 missed cleavages).
        string[] expected = ["MQIFVK", "TLTGK", "TITLEVEPSDTIENVK", "AK", "IQDK", "EGIPPDQQR", "LIFAGK", "QLEDGR", "TLSDYNIQK", "ESTLHLVLR", "LR", "GG"];
        var peptides = Digestion.Digest(Ubiquitin, Enzyme.Trypsin, All);
        Assert.Equal(expected, peptides.Select(p => p.Sequence));
        Assert.Equal(0, peptides[0].Start);
        Assert.Equal(Ubiquitin.Length, peptides[^1].End);
        Assert.All(peptides, p => Assert.Equal(Ubiquitin.Substring(p.Start, p.Length), p.Sequence));
    }

    [Fact]
    public void ProlineRuleBlocksTrypsinButNotTrypsinP()
    {
        const string s = "AAKPAARAAKGG";
        Assert.Equal(["AAKPAAR", "AAK", "GG"], Digestion.Digest(s, Enzyme.Trypsin, All).Select(p => p.Sequence));
        Assert.Equal(["AAK", "PAAR", "AAK", "GG"], Digestion.Digest(s, Enzyme.TrypsinP, All).Select(p => p.Sequence));
    }

    [Fact]
    public void MissedCleavagesAndLengthLimits()
    {
        var opts = new DigestionOptions(MissedCleavages: 2, MinLength: 1, MaxLength: 100, ClipInitiatorMethionine: false);
        var peps = Digestion.Digest("AAKBBKCCKDD".Replace('B', 'G').Replace('C', 'S').Replace('D', 'E'), Enzyme.Trypsin, opts);
        // 4 fragments → 4 + 3 + 2 peptides.
        Assert.Equal(9, peps.Count);
        Assert.Equal(3, peps.Count(p => p.MissedCleavages == 1));
        Assert.Equal(2, peps.Count(p => p.MissedCleavages == 2));
        var limited = Digestion.Digest(Ubiquitin, Enzyme.Trypsin, new DigestionOptions(1, 7, 20, false));
        Assert.All(limited, p => Assert.InRange(p.Length, 7, 20));
        Assert.Contains(limited, p => p.Sequence == "TITLEVEPSDTIENVK");
        Assert.Contains(limited, p => p.Sequence == "TITLEVEPSDTIENVKAK" && p.MissedCleavages == 1);
    }

    [Fact]
    public void InitiatorMethionineClipping()
    {
        var opts = new DigestionOptions(0, 1, 100, ClipInitiatorMethionine: true);
        var peps = Digestion.Digest(Ubiquitin, Enzyme.Trypsin, opts);
        Assert.Contains(peps, p => p.Sequence == "MQIFVK");
        Assert.Contains(peps, p => p.Sequence == "QIFVK" && p.Start == 1);
    }

    [Fact]
    public void OtherEnzymeSpecificities()
    {
        Assert.Equal(["AAK", "GGRAK", "P"], Digestion.Digest("AAKGGRAKP", Enzyme.LysC, All).Select(p => p.Sequence));
        Assert.Equal(["AAKGGR", "AKP"], Digestion.Digest("AAKGGRAKP", Enzyme.ArgC, All).Select(p => p.Sequence));
        Assert.Equal(["AAE", "GGEPKE", "DA"], Digestion.Digest("AAEGGEPKEDA", Enzyme.GluC, All).Select(p => p.Sequence));
        Assert.Equal(["AA", "DGG", "DK"], Digestion.Digest("AADGGDK", Enzyme.AspN, All).Select(p => p.Sequence));
        Assert.Equal(["AAF", "GGWPKLY", "S"], Digestion.Digest("AAFGGWPKLYS", Enzyme.Chymotrypsin, All).Select(p => p.Sequence));
        Assert.Equal(["AAF", "GGWPKL", "Y", "S"], Digestion.Digest("AAFGGWPKLYS", Enzyme.ChymotrypsinLow, All).Select(p => p.Sequence));
        Assert.Equal(["AAL", "GGF", "S"], Digestion.Digest("AALGGFS", Enzyme.Pepsin, All).Select(p => p.Sequence));
    }

    [Fact]
    public void EnzymeLookupAndCoreProteinInput()
    {
        Assert.Same(Enzyme.LysC, Enzyme.Parse("lysc"));
        Assert.Same(Enzyme.TrypsinP, Enzyme.Parse("Trypsin/P"));
        Assert.Throws<ArgumentException>(() => Enzyme.Parse("papain"));
        var prot = Protein.FromString("P0CG48", Ubiquitin);
        var peps = Digestion.Digest(prot, Enzyme.Trypsin);
        Assert.All(peps, p => Assert.Equal("P0CG48", p.ProteinId));
        Assert.Contains(peps, p => p.Sequence == "EGIPPDQQR");
        // Lower-case input digests identically.
        Assert.Equal(Digestion.Digest(Ubiquitin, Enzyme.Trypsin, All).Count, Digestion.Digest(Ubiquitin.ToLowerInvariant(), Enzyme.Trypsin, All).Count);
    }

    [Fact]
    public void ObservablePeptideCountForIbaq()
    {
        // Fully tryptic ubiquitin peptides of 6–30 aa: MQIFVK, TITLEVEPSDTIENVK, EGIPPDQQR, LIFAGK, QLEDGR, TLSDYNIQK, ESTLHLVLR.
        Assert.Equal(7, Digestion.ObservablePeptideCount("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"u8));
    }
}
