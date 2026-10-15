using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Proteomics;

namespace Biotech.Net.Proteomics.Tests;

public class ProteinAnalysisTests
{
    // The sequence used by Biopython Tests/test_ProtParam.py; expected values below are copied from that test file.
    private const string BiopythonSeq =
        "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV";

    private static readonly ProteinAnalysis A = new(BiopythonSeq);

    [Fact]
    public void CountsMatchStringCounts()
    {
        var counts = A.AminoAcidCounts();
        foreach (var (aa, n) in counts) Assert.Equal(BiopythonSeq.Count(c => c == aa), n);
        Assert.Equal(152, A.Length);
        Assert.Equal(1.0, A.AminoAcidFractions().Values.Sum(), 10);
    }

    [Fact]
    public void MolecularWeightMatchesBiopython()
    {
        // Biopython test_get_molecular_weight: 17103.16 (average); test_get_monoisotopic_molecular_weight: 17092.61.
        // Average residue masses differ in the 3rd–4th decimal between tables, so allow 0.5 Da on 17 kDa.
        Assert.Equal(17103.16, A.MolecularWeight(), 0.5);
        Assert.Equal(17092.61, A.MolecularWeight(monoisotopic: true), 0.05);
    }

    [Fact]
    public void InstabilityIndexMatchesBiopython()
    {
        // Biopython test_instability_index: 41.98 (Guruprasad 1990) → unstable.
        Assert.Equal(41.98, A.InstabilityIndex(), 2);
        Assert.False(A.IsStable);
    }

    [Fact]
    public void AromaticityGravyAndExtinctionMatchBiopython()
    {
        Assert.Equal(0.10, A.Aromaticity(), 2);                                    // test_aromaticity
        Assert.Equal(-0.5974, A.Gravy(), 4);                                       // test_gravy KyteDoolitle
        Assert.Equal(-0.0435, A.Gravy(AminoAcidScale.Eisenberg), 4);               // test_gravy Eisenberg
        var e = A.ExtinctionCoefficient();
        Assert.Equal(17420, e.Reduced);                                            // test_molar_extinction_coefficient
        Assert.Equal(17545, e.Cystines);
    }

    [Fact]
    public void IsoelectricPointMatchesBiopython()
    {
        Assert.Equal(7.72, A.IsoelectricPoint(), 2);                               // test_isoelectric_point
        Assert.Equal(0.0, A.NetCharge(7.72), 1);                                   // test_charge_at_pH (|q| < 0.05)
        // Biopython IsoelectricPoint doctests: INGAR pI 9.75, charge at pH 7 = 0.76; PETER pI 4.53.
        var ingar = new ProteinAnalysis("INGAR");
        Assert.Equal(9.75, ingar.IsoelectricPoint(), 2);
        Assert.Equal(0.76, ingar.NetCharge(7.0), 2);
        Assert.Equal(4.53, new ProteinAnalysis("PETER").IsoelectricPoint(), 2);
    }

    [Fact]
    public void EmbossSetDelegatesToKernels()
    {
        Assert.Equal(BioKernels.IsoelectricPoint("ACDEFGHIKLMNPQRSTVWY"u8), new ProteinAnalysis("ACDEFGHIKLMNPQRSTVWY").IsoelectricPoint(PKaSet.Emboss), 10);
        Assert.Equal(BioKernels.NetCharge("KKKDDE"u8, 7), new ProteinAnalysis("KKKDDE").NetCharge(7, PKaSet.Emboss), 10);
    }

    [Fact]
    public void ChargeCurveIsMonotonicAndCrossesZeroAtPi()
    {
        var curve = A.ChargeCurve(0, 14, 0.25);
        Assert.Equal(57, curve.Length);
        for (int i = 1; i < curve.Length; i++) Assert.True(curve[i].Charge < curve[i - 1].Charge);
        double pi = A.IsoelectricPoint();
        Assert.True(A.NetCharge(pi - 0.1) > 0 && A.NetCharge(pi + 0.1) < 0);
        // Lehninger set: very different pKs but still a monotonic titration with a pI in range.
        double lpi = A.IsoelectricPoint(PKaSet.Lehninger);
        Assert.InRange(lpi, 6, 10);
    }

    [Fact]
    public void SecondaryStructureFractionMatchesBiopython()
    {
        var f = A.SecondaryStructureFraction();                                    // test_secondary_structure_fraction
        Assert.Equal(0.33, f.Helix, 2);
        Assert.Equal(0.29, f.Turn, 2);
        Assert.Equal(0.37, f.Sheet, 2);
    }

    [Fact]
    public void WeightedProfileMatchesBiopythonProteinScale()
    {
        // test_protein_scale: protein_scale(kd, window 9, edge 0.4), first values.
        double[] expected = [-0.0783, 0.0358, 0.1258, 0.6950, 0.8775, 0.8350, 0.2925, 0.3383, -0.1733, -0.4142,
            -0.5292, -0.6108, -0.8308, -0.8100, -0.8208, -1.0283, -1.6300, -1.8233, -2.4267, -2.2292];
        var profile = A.Profile(AminoAcidScale.KyteDoolittle, 9, 0.4);
        Assert.Equal(BiopythonSeq.Length - 8, profile.Length);
        for (int i = 0; i < expected.Length; i++) Assert.Equal(expected[i], profile[i], 4);
        // The last value of the Biopython list.
        Assert.Equal(2.3558, profile[^1], 4);
    }

    [Fact]
    public void UnweightedProfileEqualsNativeHydropathyKernel()
    {
        var profile = AminoAcidScale.KyteDoolittle.Profile(System.Text.Encoding.ASCII.GetBytes(BiopythonSeq), 7);
        var kernel = new double[BiopythonSeq.Length - 6];
        BioKernels.Hydropathy(System.Text.Encoding.ASCII.GetBytes(BiopythonSeq), 7, kernel);
        for (int i = 0; i < kernel.Length; i++) Assert.Equal(kernel[i], profile[i], 10);
    }

    [Fact]
    public void AliphaticIndexAndFlexibilityHandComputed()
    {
        // A=6, V=5, I=5, L=18 of 152: 100/152·(6 + 2.9·5 + 3.9·23) = 72.5.
        Assert.Equal(72.5, A.AliphaticIndex(), 10);
        // Vihinen window 9 on MAEGEITTF, centre E: hand-computed 1.0143214285714286 (weights 0.25..0.8125, /5.25).
        var flex = A.Flexibility();
        Assert.Equal(BiopythonSeq.Length - 8, flex.Length);
        Assert.Equal(1.0143214285714286, flex[0], 12);
    }

    [Fact]
    public void HalfLifeFollowsNEndRule()
    {
        // ExPASy ProtParam N-end rule table: Met → 30 h (mammalian), >20 h (yeast), >10 h (E. coli).
        Assert.Equal(30, A.HalfLife().Hours);
        var yeast = A.HalfLife(HalfLifeSystem.Yeast);
        Assert.True(yeast.IsLowerBound);
        Assert.Equal(20, yeast.Hours);
        Assert.Equal(">10 hours", A.HalfLife(HalfLifeSystem.EscherichiaColi).Text);
        Assert.Equal(2.0 / 60, NEndRule.HalfLife('R', HalfLifeSystem.Yeast).Hours, 10);
        Assert.True(double.IsNaN(NEndRule.HalfLife('P', HalfLifeSystem.EscherichiaColi).Hours));
        Assert.Equal(100, NEndRule.HalfLife('v').Hours);
    }

    [Fact]
    public void ConstructorsAndValidation()
    {
        var protein = Protein.FromString("P1", BiopythonSeq.ToLowerInvariant());
        Assert.Equal(A.InstabilityIndex(), new ProteinAnalysis(protein).InstabilityIndex(), 10);
        Assert.Equal(A.Gravy(), new ProteinAnalysis(protein.Sequence).Gravy(), 10);
        Assert.Throws<ArgumentException>(() => new ProteinAnalysis(Sequence.Dna("ACGT")));
        Assert.Throws<ArgumentException>(() => new ProteinAnalysis(""));
        Assert.Contains("Instability index: 41.98", A.Report());
    }

    [Fact]
    public void AbsorbanceIsExtinctionOverMass()
    {
        Assert.Equal(17545 / A.MolecularWeight(), A.Absorbance01Percent(), 12);
        Assert.Equal(17420 / A.MolecularWeight(), A.Absorbance01Percent(cystines: false), 12);
    }

    [Fact]
    public void ScaleLookupsAndCustomScale()
    {
        Assert.Equal(4.5, AminoAcidScale.KyteDoolittle['I']);
        Assert.Equal(4.5, AminoAcidScale.KyteDoolittle['i']);
        Assert.Equal(3.0, AminoAcidScale.HoppWoods['D']);
        Assert.False(AminoAcidScale.KyteDoolittle.Contains('X'));
        var custom = new AminoAcidScale("charge", new Dictionary<char, double> { ['K'] = 1, ['R'] = 1, ['D'] = -1, ['E'] = -1 });
        Assert.Equal(0.25, custom.Average("KAAA"u8));
        Assert.Empty(custom.Profile("KA"u8, 3));
    }
}
