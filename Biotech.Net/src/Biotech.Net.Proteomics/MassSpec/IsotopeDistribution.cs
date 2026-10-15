namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>One aggregated isotope peak.</summary>
/// <param name="Mass">Abundance-weighted mean neutral mass of all isotopologues with this nominal shift.</param>
/// <param name="Abundance">Relative abundance (normalised as requested).</param>
public readonly record struct IsotopePeak(double Mass, double Abundance);

/// <summary>How <see cref="IsotopeDistribution"/> scales abundances.</summary>
public enum IsotopeNormalization
{
    /// <summary>Tallest peak = 1.</summary>
    BasePeak,

    /// <summary>Abundances sum to 1 over the returned peaks (probabilities).</summary>
    Sum,
}

/// <summary>
/// Aggregated isotope distributions computed by polynomial convolution over nominal mass shifts (the
/// "aggregated isotopic variant" approach of BRAIN / IsoSpec-style coarse calculators). Each element's isotope
/// polynomial is raised to its atom count by binary exponentiation and the results are convolved; every bin
/// keeps its probability and its probability-weighted mass, so peak masses are exact centroids.
/// </summary>
/// <remarks>Cost is O(E · log n · K²) for E elements and K retained peaks; K = 10–20 covers peptides up to ~10 kDa.</remarks>
public static class IsotopeDistribution
{
    /// <summary>Isotope pattern of an elemental composition.</summary>
    /// <param name="formula">The composition.</param>
    /// <param name="maxPeaks">Number of nominal shifts to compute (M, M+1, ...).</param>
    /// <param name="normalization">Abundance scaling.</param>
    public static IReadOnlyList<IsotopePeak> Compute(ChemicalFormula formula, int maxPeaks = 10, IsotopeNormalization normalization = IsotopeNormalization.BasePeak)
    {
        ArgumentNullException.ThrowIfNull(formula);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxPeaks, 1);
        var p = new double[maxPeaks];
        var m = new double[maxPeaks];
        p[0] = 1; // empty molecule
        foreach (var (element, count) in formula.Counts)
        {
            if (count < 0) throw new ArgumentException("Isotope distributions need non-negative atom counts.");
            var (ep, em) = ElementPolynomial(element, maxPeaks);
            var (pp, pm) = Power(ep, em, count, maxPeaks);
            (p, m) = Convolve(p, m, pp, pm, maxPeaks);
        }

        double scale = 0;
        for (int i = 0; i < maxPeaks; i++)
            scale = normalization == IsotopeNormalization.BasePeak ? Math.Max(scale, p[i]) : scale + p[i];
        var peaks = new List<IsotopePeak>(maxPeaks);
        for (int i = 0; i < maxPeaks; i++)
        {
            if (p[i] <= 0) continue;
            peaks.Add(new(m[i] / p[i], p[i] / scale));
        }
        return peaks;
    }

    /// <summary>Isotope pattern of an unmodified peptide (residues + water).</summary>
    public static IReadOnlyList<IsotopePeak> FromPeptide(string residues, int maxPeaks = 10, IsotopeNormalization normalization = IsotopeNormalization.BasePeak)
        => Compute(ChemicalFormula.FromPeptide(residues), maxPeaks, normalization);

    /// <summary>Isotope pattern of an averagine molecule of the given average mass (for unknown peptides).</summary>
    public static IReadOnlyList<IsotopePeak> Averagine(double averageMass, int maxPeaks = 10, IsotopeNormalization normalization = IsotopeNormalization.BasePeak)
        => Compute(ChemicalFormula.Averagine(averageMass), maxPeaks, normalization);

    /// <summary>Converts neutral-mass peaks to m/z for a given charge (isotope spacing shrinks to ~1.003/z).</summary>
    public static IsotopePeak[] ToMz(IReadOnlyList<IsotopePeak> peaks, int charge)
    {
        var r = new IsotopePeak[peaks.Count];
        for (int i = 0; i < r.Length; i++) r[i] = new(MassConstants.ToMz(peaks[i].Mass, charge), peaks[i].Abundance);
        return r;
    }

    private static (double[] P, double[] M) ElementPolynomial(string element, int k)
    {
        var iso = Isotopes.Of(element);
        var p = new double[k];
        var m = new double[k];
        double baseMass = iso[0].Mass;
        foreach (var (mass, ab) in iso)
        {
            int shift = (int)Math.Round(mass - baseMass);
            if (shift >= k) continue;
            p[shift] += ab;
            m[shift] += ab * mass;
        }
        return (p, m);
    }

    private static (double[] P, double[] M) Power(double[] p, double[] m, int n, int k)
    {
        var rp = new double[k];
        var rm = new double[k];
        rp[0] = 1;
        double[] bp = p, bm = m;
        while (n > 0)
        {
            if ((n & 1) != 0) (rp, rm) = Convolve(rp, rm, bp, bm, k);
            n >>= 1;
            if (n > 0) (bp, bm) = Convolve(bp, bm, bp, bm, k);
        }
        return (rp, rm);
    }

    // (p, M) are probabilities and probability-weighted mass sums per nominal shift.
    // For independent parts: P[k] = Σ pa[i] pb[j], M[k] = Σ (Ma[i] pb[j] + pa[i] Mb[j]).
    private static (double[] P, double[] M) Convolve(double[] pa, double[] ma, double[] pb, double[] mb, int k)
    {
        var p = new double[k];
        var m = new double[k];
        for (int i = 0; i < k; i++)
        {
            if (pa[i] == 0) continue;
            for (int j = 0; i + j < k; j++)
            {
                if (pb[j] == 0) continue;
                p[i + j] += pa[i] * pb[j];
                m[i + j] += ma[i] * pb[j] + pa[i] * mb[j];
            }
        }
        return (p, m);
    }
}
