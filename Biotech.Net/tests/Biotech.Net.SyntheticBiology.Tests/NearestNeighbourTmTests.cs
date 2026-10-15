using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Tests;

public class NearestNeighbourTmTests
{
    // ---- Hand-summed SantaLucia (1998) PNAS 95:1460, Table 1 values -------------------------------------------

    [Fact]
    public void Duplex_CGTTGA_MatchesHandSum()
    {
        // Stacks: CG (-10.6, -27.2) + GT (-8.4, -22.4) + TT=AA (-7.9, -22.2) + TG=CA (-8.5, -22.7) + GA (-8.2, -22.2)
        //       = dH -43.6, dS -116.7
        // Init:   5' C (G·C: +0.1, -2.8) + 3' A (A·T: +2.3, +4.1)  => dH -41.2, dS -115.4
        var d = NearestNeighbourTm.Duplex("CGTTGA");
        Assert.Equal(-41.2, d.EnthalpyKcalPerMol, 9);
        Assert.Equal(-115.4, d.EntropyCalPerMolK, 9);
        Assert.Equal(-41.2 + 310.15 * 115.4 / 1000.0, d.FreeEnergy37, 9);
    }

    [Fact]
    public void Duplex_GAATTC_AppliesSymmetryCorrection()
    {
        // Stacks: GA (-8.2, -22.2) + AA (-7.9, -22.2) + AT (-7.2, -20.4) + TT=AA (-7.9, -22.2) + TC=GA (-8.2, -22.2)
        //       = dH -39.4, dS -109.2
        // Init:   two G·C ends (+0.1, -2.8 each) => dH -39.2, dS -114.8
        // Symmetry (self-complementary): dS -1.4 => -116.2
        Assert.True(NearestNeighbourTm.IsSelfComplementary("GAATTC"));
        var d = NearestNeighbourTm.Duplex("GAATTC");
        Assert.Equal(-39.2, d.EnthalpyKcalPerMol, 9);
        Assert.Equal(-116.2, d.EntropyCalPerMolK, 9);
    }

    [Fact]
    public void Duplex_IsCaseInsensitive_AndReadsUAsT()
    {
        Assert.Equal(NearestNeighbourTm.Duplex("CGTTGA"), NearestNeighbourTm.Duplex("cguuga"));
    }

    [Fact]
    public void Duplex_SequenceAndItsReverseComplement_HaveTheSameThermodynamics()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        var a = NearestNeighbourTm.Duplex(s);
        var b = NearestNeighbourTm.Duplex(DnaUtility.ReverseComplement(s));
        Assert.Equal(a.EnthalpyKcalPerMol, b.EnthalpyKcalPerMol, 9);
        Assert.Equal(a.EntropyCalPerMolK, b.EntropyCalPerMolK, 9);
    }

    [Fact]
    public void Tm_MatchesClosedForm_ForHandSummedDuplex()
    {
        // Tm = 1000 dH / (dS + 0.368 (N-1) ln[Na+] + R ln(C/4)) - 273.15, SantaLucia 1998 eqs. 3 and 5.
        double ds = -115.4 + 0.368 * 5 * Math.Log(0.05) + 1.987 * Math.Log(0.25e-6 / 4);
        double expected = -41200.0 / ds - 273.15;
        Assert.Equal(expected, NearestNeighbourTm.Tm("CGTTGA"), 9);
    }

    [Fact]
    public void Tm_SelfComplementary_UsesStrandConcentrationWithoutFactorFour()
    {
        double ds = -116.2 + 0.368 * 5 * Math.Log(0.05) + 1.987 * Math.Log(0.25e-6);
        Assert.Equal(-39200.0 / ds - 273.15, NearestNeighbourTm.Tm("GAATTC"), 9);
    }

    // ---- Regression values for the SantaLucia (1998) unified parameters -------------------------------------------
    // These are the model Biopython exposes as Bio.SeqUtils.MeltingTemp.Tm_NN(seq, nn_table=DNA_NN2, saltcorr=5,
    // Na=50) — the stacks of Table 1, the G·C / A·T initiation terms and the entropic salt correction. The strand
    // term is R ln(k) with k = C_T/4 = 125 nM for two different strands and k = C_T = 250 nM for a
    // self-complementary one; SaltConditions.PrimerMolar is C_T, hence 0.5 uM and 0.25 uM below.
    // GCGTAGCTA was re-derived by hand (dH = -66.2 kcal/mol, dS = -181.7 e.u. before the salt correction) and
    // agrees with the value in the table to 1e-6, which is what anchors the rest.

    [Theory]
    [InlineData("AGCTTGCATGCCTGCAGGTC", 61.07810840337521)]
    [InlineData("ATGCGTACGTTAGCCGATCG", 58.120936961834616)]
    [InlineData("TTTTAAAATTTTAAAATTTT", 36.576253079980916)]
    [InlineData("GCGTAGCTA", 24.910304176086584)]
    [InlineData("CGTTGA", -2.9779466029789887)]
    public void Tm_MatchesBiopython_NonSelfComplementary(string sequence, double biopython)
    {
        Assert.False(NearestNeighbourTm.IsSelfComplementary(sequence));
        var conditions = new SaltConditions(MonovalentMolar: 0.050, PrimerMolar: 0.5e-6);
        Assert.Equal(biopython, NearestNeighbourTm.Tm(sequence, conditions), 6);
    }

    /// <summary>
    /// Self-complementary oligos, checked against values summed by hand from SantaLucia (1998) Table 1 rather
    /// than against another tool: ΔH and ΔS are the stacks plus one initiation term per end plus the −1.4 e.u.
    /// symmetry term, the salt correction is 0.368·(N−1)·ln[Na⁺], and the strand term is R·ln(C_T) with no
    /// factor of four because both strands of the duplex come from the same species.
    /// </summary>
    [Theory]
    // GAATTC: stacks GA+AA+AT+TT+TC, two G·C ends.
    [InlineData("GAATTC", -39200.0, -116.2)]
    // GAATTCGAATTC: the same stacks twice plus one CG, two G·C ends.
    [InlineData("GAATTCGAATTC", -89200.0, -252.6)]
    public void Tm_SelfComplementary_MatchesHandSummedThermodynamics(string sequence, double enthalpyCal, double entropy)
    {
        Assert.True(NearestNeighbourTm.IsSelfComplementary(sequence));
        var conditions = new SaltConditions(MonovalentMolar: 0.050, PrimerMolar: 0.25e-6);

        var duplex = NearestNeighbourTm.Duplex(sequence);
        Assert.Equal(enthalpyCal / 1000.0, duplex.EnthalpyKcalPerMol, 9);
        Assert.Equal(entropy, duplex.EntropyCalPerMolK, 9);

        double ds = entropy + 0.368 * (sequence.Length - 1) * Math.Log(0.05) + 1.987 * Math.Log(0.25e-6);
        Assert.Equal(enthalpyCal / ds - 273.15, NearestNeighbourTm.Tm(sequence, conditions), 9);
    }

    [Theory]
    [InlineData("ACGTACGTACGTACGTACGT")]
    [InlineData("GGGCCCGCGGCCGCGGGCCC")]
    public void Tm_SelfComplementary_AgreesWithTheExplicitThermodynamicOverload(string sequence)
    {
        Assert.True(NearestNeighbourTm.IsSelfComplementary(sequence));
        var conditions = new SaltConditions(MonovalentMolar: 0.050, PrimerMolar: 0.25e-6);

        var duplex = NearestNeighbourTm.Duplex(sequence);
        double gc = sequence.Count(c => c is 'G' or 'C') / (double)sequence.Length;
        double expected = NearestNeighbourTm.TmFromThermodynamics(duplex, conditions, sequence.Length, gc, selfComplementary: true);

        Assert.Equal(expected, NearestNeighbourTm.Tm(sequence, conditions), 9);
    }

    // ---- Behaviour -----------------------------------------------------------------------------------------------

    [Fact]
    public void Tm_RisesWithMonovalentSalt()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        double previous = double.NegativeInfinity;
        foreach (double na in new[] { 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0 })
        {
            double tm = NearestNeighbourTm.Tm(s, new SaltConditions(na));
            Assert.True(tm > previous, $"Tm at {na} M should exceed {previous}");
            previous = tm;
        }
    }

    [Fact]
    public void Tm_RisesWithMagnesium()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        double none = NearestNeighbourTm.Tm(s, new SaltConditions(0.050));
        double some = NearestNeighbourTm.Tm(s, new SaltConditions(0.050, 0.0015));
        double more = NearestNeighbourTm.Tm(s, new SaltConditions(0.050, 0.004));
        Assert.True(some > none);
        Assert.True(more > some);
    }

    [Fact]
    public void Tm_DntpsChelateMagnesium()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        double withDntp = NearestNeighbourTm.Tm(s, new SaltConditions(0.050, 0.0015, 0.0015));
        double noMg = NearestNeighbourTm.Tm(s, new SaltConditions(0.050));
        Assert.Equal(noMg, withDntp, 9); // all Mg2+ bound => behaves as Mg-free
    }

    [Fact]
    public void Tm_RisesWithPrimerConcentration()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        Assert.True(NearestNeighbourTm.Tm(s, new SaltConditions(PrimerMolar: 1e-6))
                    > NearestNeighbourTm.Tm(s, new SaltConditions(PrimerMolar: 0.1e-6)));
    }

    [Fact]
    public void Tm_RisesWithGcContent()
    {
        // Same length, increasing GC (none self-complementary).
        string[] series =
        [
            "ATTATAATTAATATTATAAT", // 0 %
            "ATTAGAATTCATATGATACT", // 20 %
            "ATGAGCATTCAGATGCTACT", // 40 %
            "AGGAGCAGTCAGCTGCGACT", // 60 %
            "GGCAGCGGTCCGCTGCGGCT", // 80 %
        ];
        double previous = double.NegativeInfinity;
        foreach (var s in series)
        {
            double tm = NearestNeighbourTm.Tm(s);
            Assert.True(tm > previous, $"{s}: {tm} should exceed {previous}");
            previous = tm;
        }
    }

    [Theory]
    [InlineData("ACGTACGT", 24.0)]   // 4 A/T x 2 + 4 G/C x 4
    [InlineData("AAAAAAAAAA", 20.0)]
    [InlineData("GGGGCCCC", 32.0)]
    [InlineData("acgu", 12.0)]
    public void WallaceTm_IsExact(string sequence, double expected)
        => Assert.Equal(expected, NearestNeighbourTm.WallaceTm(sequence));

    [Fact]
    public void GcContentTm_IsExact()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC"; // 20 nt, 12 G/C = 60 %
        double expected = 81.5 + 16.6 * Math.Log10(0.05) + 0.41 * 60.0 - 600.0 / 20;
        Assert.Equal(expected, NearestNeighbourTm.GcContentTm(s), 9);
        Assert.Equal(81.5 + 16.6 * Math.Log10(0.2) + 0.41 * 60.0 - 30.0, NearestNeighbourTm.GcContentTm(s, 0.2), 9);
    }

    [Theory]
    [InlineData("GAATTC", true)]
    [InlineData("GGATCC", true)]
    [InlineData("ACGT", true)]
    [InlineData("ACGTA", false)]
    [InlineData("AAAA", false)]
    [InlineData("GANTTC", false)]
    public void IsSelfComplementary_Works(string sequence, bool expected)
        => Assert.Equal(expected, NearestNeighbourTm.IsSelfComplementary(sequence));

    [Fact]
    public void InvalidInput_Throws()
    {
        Assert.Throws<ArgumentException>(() => NearestNeighbourTm.Duplex("A"));
        Assert.Throws<ArgumentException>(() => NearestNeighbourTm.Tm("ACGNTT"));
        Assert.Throws<ArgumentException>(() => NearestNeighbourTm.Tm("ACGT", new SaltConditions(0.0)));
    }

    [Fact]
    public void FreeEnergyAt37_EqualsFreeEnergy37()
    {
        var d = NearestNeighbourTm.Duplex("AGCTTGCATGCCTGCAGGTC");
        Assert.Equal(d.FreeEnergy37, d.FreeEnergyAt(37.0), 12);
        Assert.True(d.FreeEnergyAt(60.0) > d.FreeEnergy37); // less stable when hotter
    }
}
