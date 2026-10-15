using Biotech.Net.Proteomics.MassSpec;

namespace Biotech.Net.Proteomics.Tests;

public class IsotopeTests
{
    // IUPAC 2013 representative abundances, written out independently of the library table.
    private const double C12 = 0.9893, C13 = 0.0107, H1 = 0.999885, H2 = 0.000115, O16 = 0.99757, O17 = 0.00038, O18 = 0.00205;

    [Fact]
    public void GlucoseMonoisotopicAndFirstIsotopePeaksMatchClosedForm()
    {
        var glucose = ChemicalFormula.Parse("C6H12O6");
        Assert.Equal(180.06338810, glucose.MonoisotopicMass, 6);
        var peaks = IsotopeDistribution.Compute(glucose, maxPeaks: 10, IsotopeNormalization.Sum);

        // P(M) = C12^6 · H1^12 · O16^6 ≈ 0.9226.
        double p0 = Math.Pow(C12, 6) * Math.Pow(H1, 12) * Math.Pow(O16, 6);
        // P(M+1) = P(M) · (6·C13/C12 + 12·H2/H1 + 6·O17/O16).
        double p1 = p0 * (6 * C13 / C12 + 12 * H2 / H1 + 6 * O17 / O16);
        Assert.Equal(p0, peaks[0].Abundance, 10);
        Assert.Equal(0.92265, peaks[0].Abundance, 4);
        Assert.Equal(p1, peaks[1].Abundance, 10);
        Assert.Equal(180.06338810, peaks[0].Mass, 6);
        // M+1 centroid sits between the 2H/17O spacings and the 13C spacing, dominated by 13C.
        Assert.InRange(peaks[1].Mass - peaks[0].Mass, 1.0030, 1.0040);
        Assert.Equal(1.0, peaks.Sum(p => p.Abundance), 12);

        var basePeak = IsotopeDistribution.Compute(glucose);
        Assert.Equal(1.0, basePeak[0].Abundance, 12);
        Assert.Equal(p1 / p0, basePeak[1].Abundance, 10);
    }

    [Fact]
    public void M2PeakMatchesFullCombinatorialSum()
    {
        // Every way to gain +2 in C6H12O6: one 18O, two heavy atoms among 13C/2H/17O (same or mixed elements).
        var peaks = IsotopeDistribution.Compute(ChemicalFormula.Parse("C6H12O6"), 10, IsotopeNormalization.Sum);
        double p0 = Math.Pow(C12, 6) * Math.Pow(H1, 12) * Math.Pow(O16, 6);
        double rc = C13 / C12, rh = H2 / H1, r17 = O17 / O16, r18 = O18 / O16;
        double dominant = 6 * r18 + 15 * rc * rc;
        double exact = dominant + 66 * rh * rh + 15 * r17 * r17 + 6 * 12 * rc * rh + 6 * 6 * rc * r17 + 12 * 6 * rh * r17;
        Assert.Equal(p0 * exact, peaks[2].Abundance, 10);
        Assert.True(peaks[2].Abundance > p0 * dominant);
    }

    [Fact]
    public void DichlorineHasBinomialPatternWithOddShiftsAbsent()
    {
        var peaks = IsotopeDistribution.Compute(ChemicalFormula.Parse("Cl2"), 6, IsotopeNormalization.Sum);
        Assert.Equal(3, peaks.Count);
        const double a = 0.7576, b = 0.2424;
        Assert.Equal(a * a, peaks[0].Abundance, 12);
        Assert.Equal(2 * a * b, peaks[1].Abundance, 12);
        Assert.Equal(b * b, peaks[2].Abundance, 12);
        Assert.Equal(2 * 34.96885268, peaks[0].Mass, 8);
        Assert.Equal(34.96885268 + 36.96590259, peaks[1].Mass, 8);
        Assert.Equal(2 * 36.96590259, peaks[2].Mass, 8);
    }

    [Fact]
    public void PeptideIsotopeEnvelope()
    {
        var formula = ChemicalFormula.FromPeptide("PEPTIDE");
        Assert.Equal("C34H53N7O15", formula.ToString());
        Assert.Equal(799.359964, formula.MonoisotopicMass, 5);
        var peaks = IsotopeDistribution.FromPeptide("PEPTIDE", 8, IsotopeNormalization.Sum);
        Assert.Equal(799.359964, peaks[0].Mass, 5);
        Assert.Equal(1.0, peaks.Sum(p => p.Abundance), 12);
        // ~800 Da: monoisotopic peak is the base peak, roughly 65% of the envelope; M+1 ~ 26%.
        Assert.InRange(peaks[0].Abundance, 0.60, 0.68);
        Assert.InRange(peaks[1].Abundance, 0.23, 0.29);
        for (int i = 1; i < peaks.Count; i++) Assert.True(peaks[i].Abundance < peaks[i - 1].Abundance);

        // m/z conversion: isotope spacing shrinks to ~1.003/z.
        var mz = IsotopeDistribution.ToMz(peaks, 2);
        Assert.Equal(new Peptide("PEPTIDE").Mz(2), mz[0].Mass, 5);
        Assert.InRange(mz[1].Mass - mz[0].Mass, 0.500, 0.503);
        Assert.Equal(peaks[1].Abundance, mz[1].Abundance);
    }

    [Theory]
    [InlineData(1000.0, 0)]
    [InlineData(1500.0, 0)]
    [InlineData(2500.0, 1)]
    [InlineData(4000.0, 2)]
    public void AveragineBasePeakShiftsWithMass(double mass, int expectedBasePeak)
    {
        var formula = ChemicalFormula.Averagine(mass);
        Assert.Equal(mass, formula.AverageMass, 0.6);
        var peaks = IsotopeDistribution.Averagine(mass, 12);
        int basePeak = 0;
        for (int i = 1; i < peaks.Count; i++) if (peaks[i].Abundance > peaks[basePeak].Abundance) basePeak = i;
        Assert.Equal(expectedBasePeak, basePeak);
        Assert.Equal(1.0, peaks[basePeak].Abundance, 12);
    }

    [Fact]
    public void Averagine1000DaMonoisotopicFraction()
    {
        // Senko averagine at 1 kDa: monoisotopic peak ≈ 57% of the envelope.
        var peaks = IsotopeDistribution.Averagine(1000, 10, IsotopeNormalization.Sum);
        Assert.InRange(peaks[0].Abundance, 0.54, 0.60);
        Assert.Equal(1.0, peaks.Sum(p => p.Abundance), 12);
    }

    [Fact]
    public void InvalidInputsAreRejected()
    {
        Assert.Throws<ArgumentException>(() => IsotopeDistribution.Compute(ChemicalFormula.Parse("C2H-1")));
        Assert.Throws<ArgumentOutOfRangeException>(() => IsotopeDistribution.Compute(ChemicalFormula.Parse("C2"), 0));
        Assert.Throws<ArgumentException>(() => IsotopeDistribution.Compute(ChemicalFormula.Parse("Xx2")));
        // Empty formula: a single peak at mass 0 with abundance 1.
        var empty = IsotopeDistribution.Compute(new ChemicalFormula(), 3);
        Assert.Single(empty);
    }
}
