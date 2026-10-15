using Biotech.Net.Native;
using Biotech.Net.Proteomics.MassSpec;

namespace Biotech.Net.Proteomics.Tests;

public class MassTests
{
    [Fact]
    public void PeptideMassesMatchReference()
    {
        // PEPTIDE: monoisotopic 799.359964 Da, [M+H]+ 800.367240, [M+2H]2+ 400.687258 (e.g. ProteoMass / Expasy PeptideMass).
        var p = new Peptide("PEPTIDE");
        Assert.Equal(799.359964, p.MonoisotopicMass, 5);
        Assert.Equal(800.367240, p.Mz(1), 5);
        Assert.Equal(400.687258, p.Mz(2), 5);
        Assert.Equal(799.359964, ResidueMasses.PeptideMass("peptide"), 5);
        // Consistent with the Core kernel to within table rounding.
        Assert.Equal(BioKernels.ProteinMass("PEPTIDE"u8, monoisotopic: true), p.MonoisotopicMass, 3);
        Assert.Equal(BioKernels.ProteinMass("PEPTIDE"u8), p.AverageMass, 1);
    }

    [Fact]
    public void MzConversionsRoundTrip()
    {
        double mz = MassConstants.ToMz(1500.0, 3);
        Assert.Equal((1500 + 3 * 1.007276466621) / 3, mz, 10);
        Assert.Equal(1500.0, MassConstants.ToNeutralMass(mz, 3), 10);
        Assert.Equal(10.0, MassConstants.PpmError(1000.01, 1000.0), 6);
        Assert.Throws<ArgumentOutOfRangeException>(() => MassConstants.ToMz(100, 0));
    }

    [Fact]
    public void BAndYLaddersOfPeptide()
    {
        // Reference b/y ions of PEPTIDE (1+), as listed by Protein Prospector MS-Product / pyteomics.
        double[] b = [98.06004, 227.10263, 324.15540, 425.20308, 538.28714, 653.31408];
        double[] y = [148.06043, 263.08738, 376.17144, 477.21912, 574.27188, 703.31448];
        var p = new Peptide("PEPTIDE");
        var lb = FragmentCalculator.Ladder(p, IonType.B);
        var ly = FragmentCalculator.Ladder(p, IonType.Y);
        for (int i = 0; i < 6; i++)
        {
            Assert.Equal(b[i], lb[i], 4);
            Assert.Equal(y[i], ly[i], 4);
        }
    }

    [Fact]
    public void ComplementaryIonsSumToPrecursor()
    {
        // b_i + y_(n−i) = M + 2H+ for singly charged ions.
        var p = new Peptide("SAMPLERK");
        var b = FragmentCalculator.Ladder(p, IonType.B);
        var y = FragmentCalculator.Ladder(p, IonType.Y);
        int n = p.Length;
        for (int i = 1; i < n; i++)
            Assert.Equal(p.MonoisotopicMass + 2 * MassConstants.Proton, b[i - 1] + y[n - i - 1], 8);
    }

    [Fact]
    public void OtherSeriesOffsetsAndCharges()
    {
        var p = new Peptide("PEPTIDE");
        var b = FragmentCalculator.Ladder(p, IonType.B);
        var a = FragmentCalculator.Ladder(p, IonType.A);
        var c = FragmentCalculator.Ladder(p, IonType.C);
        var y = FragmentCalculator.Ladder(p, IonType.Y);
        var x = FragmentCalculator.Ladder(p, IonType.X);
        var z = FragmentCalculator.Ladder(p, IonType.Z);
        var zd = FragmentCalculator.Ladder(p, IonType.ZDot);
        Assert.Equal(b[1] - 27.994915, a[1], 5);     // a = b − CO
        Assert.Equal(b[1] + 17.026549, c[1], 5);     // c = b + NH3
        Assert.Equal(y[1] + 25.979265, x[1], 5);     // x = y + CO − 2H
        Assert.Equal(y[1] - 17.026549, z[1], 5);     // z = y − NH3
        Assert.Equal(z[1] + 1.007825, zd[1], 5);     // z• = z + H
        var b2 = FragmentCalculator.Ladder(p, IonType.B, 2);
        Assert.Equal((b[3] + MassConstants.Proton) / 2, b2[3], 8);
        var all = FragmentCalculator.Fragments("PEPTIDE", [IonType.B, IonType.Y], maxCharge: 2);
        Assert.Equal(24, all.Count);
        Assert.Contains(all, f => f.Label == "y2++");
        Assert.Contains(all, f => f.Label == "b3" && Math.Abs(f.Mz - 324.15540) < 1e-4);
    }

    [Fact]
    public void ModificationsShiftMassesAndFragments()
    {
        // Carbamidomethyl C +57.021464 and Oxidation M +15.994915 (Unimod 4 and 35).
        var plain = new Peptide("ACMK");
        var mods = new Modification?[] { null, Modification.Carbamidomethyl, Modification.Oxidation, null };
        var mod = new Peptide("ACMK", mods);
        Assert.Equal(plain.MonoisotopicMass + 57.021464 + 15.994915, mod.MonoisotopicMass, 6);
        Assert.Equal("AC[+57.021]M[+15.995]K", mod.ModifiedSequence);
        var by = FragmentCalculator.Ladder(mod, IonType.B);
        var bp = FragmentCalculator.Ladder(plain, IonType.B);
        Assert.Equal(bp[0], by[0], 8);
        Assert.Equal(bp[1] + 57.021464, by[1], 6);
        Assert.Equal(2, mod.ModificationCount);
    }

    [Fact]
    public void NTerminalModificationShiftsBIonsOnly()
    {
        var ac = new Peptide("PEPTIDE", nTermModification: Modification.AcetylProteinNTerm);
        var b = FragmentCalculator.Ladder(ac, IonType.B);
        var y = FragmentCalculator.Ladder(ac, IonType.Y);
        Assert.Equal(98.06004 + 42.010565, b[0], 4);
        Assert.Equal(148.06043, y[0], 4);
        Assert.StartsWith("[+42.011]", ac.ModifiedSequence);
    }

    [Fact]
    public void VariableModificationExpansion()
    {
        Assert.Equal(4, PeptideModifier.Expand("MAMK", null, [Modification.Oxidation], 2).Count);   // none, M1, M3, both
        Assert.Equal(3, PeptideModifier.Expand("MAMK", null, [Modification.Oxidation], 1).Count);
        Assert.Equal(8, PeptideModifier.Expand("STY", null, [Modification.Phospho], 3).Count);       // 2^3
        Assert.Equal(5, PeptideModifier.Expand("STYK", null, [Modification.Phospho], 3, maxIsoforms: 5).Count);
        // Fixed modification takes its site; the variable one cannot stack on it.
        var forms = PeptideModifier.Expand("CMC", [Modification.Carbamidomethyl], [Modification.Oxidation], 2);
        Assert.Equal(2, forms.Count);
        Assert.All(forms, f => Assert.Equal(Modification.Carbamidomethyl, f.ModificationAt(0)));
        // Protein N-terminal acetylation only applies at the protein N-terminus.
        Assert.Single(PeptideModifier.Expand("PEPTIDE", null, [Modification.AcetylProteinNTerm]));
        Assert.Equal(2, PeptideModifier.Expand("PEPTIDE", null, [Modification.AcetylProteinNTerm], isProteinNTerm: true).Count);
        // Pyro-Glu needs an N-terminal Q.
        Assert.Equal(2, PeptideModifier.Expand("QPEPTIDE", null, [Modification.PyroGluQ]).Count);
        Assert.Single(PeptideModifier.Expand("PEPTIDEQ", null, [Modification.PyroGluQ]));
    }

    [Fact]
    public void ToleranceWindows()
    {
        var ppm = Tolerance.Ppm(10);
        Assert.Equal(0.01, ppm.HalfWidth(1000), 12);
        Assert.True(ppm.Contains(1000, 1000.009));
        Assert.False(ppm.Contains(1000, 1000.011));
        Assert.True(Tolerance.Da(0.5).Contains(100, 100.4));
        Assert.Equal("±10 ppm", ppm.ToString());
    }

    [Fact]
    public void UnknownResidueRejected()
    {
        Assert.Throws<ArgumentException>(() => new Peptide("PEPXIDE"));
        Assert.Throws<ArgumentException>(() => ResidueMasses.PeptideMass("PEPBIDE"));
    }
}
