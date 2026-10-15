using Biotech.Net.Proteomics;

namespace Biotech.Net.Proteomics.Tests;

public class StructurePredictionTests
{
    // Human preproinsulin (UniProt P01308): signal peptide 1–24, cleaved after ...ALWGPDPAAA|FVNQ.
    private const string Preproinsulin =
        "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN";

    // Human serum albumin precursor (P02768) N-terminus: signal peptide 1–18 (MKWVTFISLLFLFSSAYS|RGVFRR).
    private const string Albumin = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL";

    // Human ubiquitin (P0CG48 unit): cytosolic, no signal peptide, no TM helix.
    private const string Ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    [Fact]
    public void SignalPeptideCleavageSitesOfKnownPrecursors()
    {
        var ins = SignalPeptidePredictor.Predict(Preproinsulin);
        Assert.True(ins.HasSignalPeptide);
        Assert.Equal(24, ins.CleavageSite);
        var alb = SignalPeptidePredictor.Predict(Albumin);
        Assert.True(alb.HasSignalPeptide);
        Assert.Equal(18, alb.CleavageSite);
    }

    [Fact]
    public void NoSignalPeptideInCytosolicProteinOrShortInput()
    {
        var ub = SignalPeptidePredictor.Predict(Ubiquitin);
        Assert.False(ub.HasSignalPeptide);
        Assert.Equal(-1, ub.CleavageSite);
        Assert.False(SignalPeptidePredictor.Predict("MKLLV").HasSignalPeptide);
    }

    [Fact]
    public void TransmembraneSegmentFoundInSyntheticSingleSpanProtein()
    {
        string loop = "MSDKERSTNGQDEKRPSG";
        string tm = "LLIVALFGVLIAVLLIFAV"; // 19 hydrophobic residues, mean KD ≈ 3
        string seq = loop + tm + "KRDEQSNGTKERDDPQSE";
        var segs = TransmembranePredictor.Predict(seq);
        var s = Assert.Single(segs);
        Assert.InRange(s.Start, loop.Length - 6, loop.Length + 1);
        Assert.InRange(s.End, loop.Length + tm.Length - 1, loop.Length + tm.Length + 6);
        Assert.True(s.MaxWindowHydropathy >= 2.5);
        Assert.Empty(TransmembranePredictor.Predict(Ubiquitin));
    }

    [Fact]
    public void TwoSeparatedHelicesGiveTwoSegments()
    {
        string loop = "DEKRSGNPQDEKRSTGNPQDEKR";
        string seq = loop + "LLIVALFGVLIAVLLIFAV" + loop + "IVLLAFGVILAVLFLIVAL" + loop;
        Assert.Equal(2, TransmembranePredictor.Predict(seq).Count);
    }

    [Fact]
    public void ChouFasmanPredictsObviousHelixAndSheet()
    {
        Assert.Equal(new string('H', 20), ChouFasman.Predict("AEELLKKAMEALLKEAMELA"));
        string strand = ChouFasman.Predict("VIVYVIVTVIYV");
        Assert.Equal(12, strand.Length);
        Assert.True(strand.Count(c => c == 'E') >= 10);
        string mixed = ChouFasman.Predict("AEELLKKAMEALLKEPGNGDPSNGVIVYVIVTVIYV");
        Assert.Contains('H', mixed);
        Assert.Contains('E', mixed);
        Assert.Contains('T', mixed);
        Assert.Equal(1.42, ChouFasman.Parameters('A').Helix, 10);   // Chou–Fasman 1978 table: Ala P(α) = 142
        Assert.Equal(1.70, ChouFasman.Parameters('V').Sheet, 10);   // Val P(β) = 170
    }
}
