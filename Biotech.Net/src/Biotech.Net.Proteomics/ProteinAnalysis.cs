using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Proteomics;

/// <summary>pKa sets for net charge and isoelectric point calculations.</summary>
public enum PKaSet
{
    /// <summary>
    /// Bjellqvist-derived set with terminus-specific values, as used by Biopython <c>IsoelectricPoint</c>
    /// (N-term 7.5 or residue-specific, C-term 3.55 / 4.55 (D) / 4.75 (E), K 10, R 12, H 5.98, D 4.05, E 4.45, C 9, Y 10).
    /// Closest to ExPASy ProtParam / Compute pI.
    /// </summary>
    Bjellqvist,

    /// <summary>EMBOSS <c>iep</c> set (the one <see cref="BioKernels.IsoelectricPoint"/> uses; Rust-accelerated).</summary>
    Emboss,

    /// <summary>Lehninger textbook values (N-term 9.69, C-term 2.34, K 10.5, R 12.4, H 6.0, D 3.86, E 4.25, C 8.33, Y 10.07).</summary>
    Lehninger,
}

/// <summary>Molar extinction coefficients at 280 nm in water, M⁻¹ cm⁻¹ (Pace et al., Protein Sci. 4:2411, 1995).</summary>
/// <param name="Reduced">All cysteines reduced: 5500·W + 1490·Y.</param>
/// <param name="Cystines">All cysteine pairs form cystines: adds 125 per pair.</param>
public readonly record struct ExtinctionCoefficient(double Reduced, double Cystines);

/// <summary>Biopython-style secondary-structure fractions from composition alone.</summary>
/// <param name="Helix">Fraction of E, M, A, L, K (helix formers).</param>
/// <param name="Turn">Fraction of N, P, G, S, D (turn formers).</param>
/// <param name="Sheet">Fraction of V, I, Y, F, W, L, T (sheet formers).</param>
public readonly record struct SecondaryStructureFraction(double Helix, double Turn, double Sheet);

/// <summary>
/// Physicochemical parameters of a protein with parity to ExPASy ProtParam and Biopython <c>ProteinAnalysis</c>:
/// molecular weight, pI, net charge, extinction coefficient, instability index, aliphatic index, GRAVY,
/// aromaticity, flexibility, composition and the N-end-rule half-life.
/// </summary>
/// <remarks>
/// The residues are copied and upper-cased once and counted once; every property is computed from those counts or
/// in a single pass. Non-standard letters (X, B, Z, U, ...) count towards the length but towards no residue class.
/// </remarks>
public sealed class ProteinAnalysis
{
    private readonly byte[] _seq;
    private readonly int[] _counts = new int[20];

    /// <summary>Analyses ASCII residues (any case).</summary>
    public ProteinAnalysis(ReadOnlySpan<byte> residues)
    {
        if (residues.IsEmpty) throw new ArgumentException("A protein needs at least one residue.", nameof(residues));
        _seq = AminoAcids.ToUpperBytes(residues);
        foreach (byte b in _seq)
        {
            int i = AminoAcids.IndexOf(b);
            if (i >= 0) _counts[i]++;
        }
    }

    /// <summary>Analyses a one-letter string.</summary>
    public ProteinAnalysis(string residues) : this(Encoding.ASCII.GetBytes(residues ?? throw new ArgumentNullException(nameof(residues)))) { }

    /// <summary>Analyses a Core protein sequence.</summary>
    public ProteinAnalysis(Sequence sequence) : this(AminoAcids.RequireProtein(sequence)) { }

    /// <summary>Analyses a Core protein.</summary>
    public ProteinAnalysis(Protein protein) : this((protein ?? throw new ArgumentNullException(nameof(protein))).Sequence.Span) { }

    /// <summary>Number of residues.</summary>
    public int Length => _seq.Length;

    /// <summary>The upper-cased residues.</summary>
    public ReadOnlySpan<byte> Residues => _seq;

    /// <summary>Count of one standard residue.</summary>
    public int Count(char residue)
    {
        int i = AminoAcids.IndexOf(residue);
        return i < 0 ? 0 : _counts[i];
    }

    /// <summary>Counts of the 20 standard residues.</summary>
    public IReadOnlyDictionary<char, int> AminoAcidCounts()
    {
        var d = new Dictionary<char, int>(20);
        for (int i = 0; i < 20; i++) d[AminoAcids.Standard[i]] = _counts[i];
        return d;
    }

    /// <summary>Fraction (0–1) of each standard residue over the full length.</summary>
    public IReadOnlyDictionary<char, double> AminoAcidFractions()
    {
        var d = new Dictionary<char, double>(20);
        for (int i = 0; i < 20; i++) d[AminoAcids.Standard[i]] = (double)_counts[i] / Length;
        return d;
    }

    /// <summary>Average (default) or monoisotopic molecular weight in Da, including one water.</summary>
    public double MolecularWeight(bool monoisotopic = false) => BioKernels.ProteinMass(_seq, monoisotopic);

    /// <summary>Net charge at <paramref name="pH"/> under the Henderson–Hasselbalch model.</summary>
    public double NetCharge(double pH, PKaSet pKaSet = PKaSet.Bjellqvist)
    {
        if (pKaSet == PKaSet.Emboss) return BioKernels.NetCharge(_seq, pH);
        var pk = Terms(pKaSet);
        static double Pos(double pK, double ph) => 1.0 / (Math.Pow(10, ph - pK) + 1.0);
        static double Neg(double pK, double ph) => 1.0 / (Math.Pow(10, pK - ph) + 1.0);
        double positive = Pos(pk.NTerm, pH) + C('K') * Pos(pk.K, pH) + C('R') * Pos(pk.R, pH) + C('H') * Pos(pk.H, pH);
        double negative = Neg(pk.CTerm, pH) + C('D') * Neg(pk.D, pH) + C('E') * Neg(pk.E, pH)
            + C('C') * Neg(pk.C, pH) + C('Y') * Neg(pk.Y, pH);
        return positive - negative;
    }

    /// <summary>Isoelectric point: the pH of zero net charge, found by bisection on [0, 14] to 1e-6.</summary>
    /// <remarks>With <see cref="PKaSet.Bjellqvist"/> this reproduces Biopython <c>isoelectric_point()</c> (e.g. INGAR 9.75, PETER 4.53).</remarks>
    public double IsoelectricPoint(PKaSet pKaSet = PKaSet.Bjellqvist)
    {
        if (pKaSet == PKaSet.Emboss) return BioKernels.IsoelectricPoint(_seq);
        double lo = 0, hi = 14;
        while (hi - lo > 1e-6)
        {
            double mid = 0.5 * (lo + hi);
            if (NetCharge(mid, pKaSet) > 0) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /// <summary>Net charge sampled from <paramref name="fromPH"/> to <paramref name="toPH"/> (inclusive) in steps of <paramref name="step"/>.</summary>
    public (double PH, double Charge)[] ChargeCurve(double fromPH = 0, double toPH = 14, double step = 0.5, PKaSet pKaSet = PKaSet.Bjellqvist)
    {
        if (step <= 0) throw new ArgumentOutOfRangeException(nameof(step));
        if (toPH < fromPH) throw new ArgumentException("toPH must be ≥ fromPH.");
        int n = (int)Math.Floor((toPH - fromPH) / step + 1e-9) + 1;
        var curve = new (double, double)[n];
        for (int i = 0; i < n; i++)
        {
            double ph = fromPH + i * step;
            curve[i] = (ph, NetCharge(ph, pKaSet));
        }
        return curve;
    }

    /// <summary>Molar extinction coefficients at 280 nm (ExPASy ProtParam, Pace 1995).</summary>
    public ExtinctionCoefficient ExtinctionCoefficient()
    {
        double reduced = C('W') * 5500.0 + C('Y') * 1490.0;
        return new(reduced, reduced + (C('C') / 2) * 125.0);
    }

    /// <summary>Absorbance of a 0.1 % (1 g/L) solution at 280 nm, assuming all cysteines form cystines (ProtParam "Abs 0.1%").</summary>
    public double Absorbance01Percent(bool cystines = true)
    {
        var e = ExtinctionCoefficient();
        return (cystines ? e.Cystines : e.Reduced) / MolecularWeight();
    }

    /// <summary>
    /// Guruprasad instability index (Protein Eng. 4:155, 1990): <c>(10 / L) · Σ DIWV(xᵢ, xᵢ₊₁)</c>. Values above 40
    /// predict an unstable protein. Dipeptides containing a non-standard residue contribute 0.
    /// </summary>
    public double InstabilityIndex()
    {
        double score = 0;
        for (int i = 0; i + 1 < _seq.Length; i++)
        {
            int a = AminoAcids.IndexOf(_seq[i]), b = AminoAcids.IndexOf(_seq[i + 1]);
            if (a >= 0 && b >= 0) score += Diwv[a * 20 + b];
        }
        return 10.0 / Length * score;
    }

    /// <summary>True when the instability index is below 40 (ProtParam's "stable" classification).</summary>
    public bool IsStable => InstabilityIndex() < 40;

    /// <summary>Aliphatic index (Ikai 1980): <c>X(A) + 2.9·X(V) + 3.9·(X(I) + X(L))</c> with mole percents.</summary>
    public double AliphaticIndex()
        => 100.0 / Length * (C('A') + 2.9 * C('V') + 3.9 * (C('I') + C('L')));

    /// <summary>Grand average of hydropathy (Kyte–Doolittle unless another scale is given).</summary>
    public double Gravy(AminoAcidScale? scale = null) => (scale ?? AminoAcidScale.KyteDoolittle).Average(_seq);

    /// <summary>Relative frequency of F + W + Y (Lobry 1994), as Biopython <c>aromaticity()</c>.</summary>
    public double Aromaticity() => (double)(C('F') + C('W') + C('Y')) / Length;

    /// <summary>Composition-based secondary-structure fractions, identical to Biopython <c>secondary_structure_fraction()</c>.</summary>
    public SecondaryStructureFraction SecondaryStructureFraction()
    {
        double l = Length;
        return new((C('E') + C('M') + C('A') + C('L') + C('K')) / l,
            (C('N') + C('P') + C('G') + C('S') + C('D')) / l,
            (C('V') + C('I') + C('Y') + C('F') + C('W') + C('L') + C('T')) / l);
    }

    /// <summary>
    /// Vihinen flexibility profile (window 9, weights 0.25 → 1 towards the centre, normalised by 5.25).
    /// Value <c>i</c> belongs to residue <c>i + 4</c>. Note: Biopython's implementation takes the centre one residue
    /// late and drops the last window; this version uses the true window centre.
    /// </summary>
    public double[] Flexibility()
    {
        ReadOnlySpan<double> w = [0.25, 0.4375, 0.625, 0.8125];
        var flex = AminoAcidScale.Flexibility;
        int n = Length - 9 + 1;
        if (n <= 0) return [];
        var output = new double[n];
        for (int i = 0; i < n; i++)
        {
            double s = flex[_seq[i + 4]];
            for (int j = 0; j < 4; j++) s += w[j] * (flex[_seq[i + j]] + flex[_seq[i + 8 - j]]);
            output[i] = s / 5.25;
        }
        return output;
    }

    /// <summary>Sliding-window profile on any scale (see <see cref="AminoAcidScale.Profile"/>).</summary>
    public double[] Profile(AminoAcidScale scale, int window = 9, double edge = 1.0)
    {
        ArgumentNullException.ThrowIfNull(scale);
        return scale.Profile(_seq, window, edge);
    }

    /// <summary>Estimated half-life from the N-terminal residue (N-end rule), as reported by ExPASy ProtParam.</summary>
    public HalfLife HalfLife(HalfLifeSystem system = HalfLifeSystem.MammalianReticulocytes) => NEndRule.HalfLife((char)_seq[0], system);

    /// <summary>A ProtParam-style multi-line text report.</summary>
    public string Report()
    {
        var ci = CultureInfo.InvariantCulture;
        var e = ExtinctionCoefficient();
        var sb = new StringBuilder();
        sb.AppendLine(ci, $"Number of amino acids: {Length}");
        sb.AppendLine(ci, $"Molecular weight: {MolecularWeight():F2}");
        sb.AppendLine(ci, $"Theoretical pI: {IsoelectricPoint():F2}");
        sb.AppendLine(ci, $"Negatively charged residues (Asp + Glu): {C('D') + C('E')}");
        sb.AppendLine(ci, $"Positively charged residues (Arg + Lys): {C('R') + C('K')}");
        sb.AppendLine(ci, $"Extinction coefficient (cystines / reduced): {e.Cystines:F0} / {e.Reduced:F0}");
        sb.AppendLine(ci, $"Estimated half-life (mammalian reticulocytes, in vitro): {HalfLife().Text}");
        double ii = InstabilityIndex();
        sb.AppendLine(ci, $"Instability index: {ii:F2} ({(ii < 40 ? "stable" : "unstable")})");
        sb.AppendLine(ci, $"Aliphatic index: {AliphaticIndex():F2}");
        sb.Append(ci, $"GRAVY: {Gravy():F3}");
        return sb.ToString();
    }

    private int C(char residue) => _counts[AminoAcids.IndexOf(residue)];

    private readonly record struct PkTable(double NTerm, double CTerm, double K, double R, double H, double D, double E, double C, double Y);

    private PkTable Terms(PKaSet set)
    {
        if (set == PKaSet.Lehninger) return new(9.69, 2.34, 10.5, 12.4, 6.0, 3.86, 4.25, 8.33, 10.07);
        // Biopython IsoelectricPoint: terminus-specific pKs.
        double nTerm = (char)_seq[0] switch
        {
            'A' => 7.59, 'M' => 7.0, 'S' => 6.93, 'P' => 8.36, 'T' => 6.82, 'V' => 7.44, 'E' => 7.7, _ => 7.5,
        };
        double cTerm = (char)_seq[^1] switch { 'D' => 4.55, 'E' => 4.75, _ => 3.55 };
        return new(nTerm, cTerm, 10.0, 12.0, 5.98, 4.05, 4.45, 9.0, 10.0);
    }

    // Guruprasad dipeptide instability weight values, DIWV[first, second], both in AminoAcids.Standard order.
    // Transcribed from Biopython Bio/SeqUtils/ProtParamData.py (same values as ExPASy ProtParam).
    private static readonly double[] Diwv =
    [
        /* A */ 1.0, 44.94, -7.49, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        /* C */ 1.0, 1.0, 20.26, 1.0, 1.0, 1.0, 33.6, 1.0, 1.0, 20.26, 33.6, 1.0, 20.26, -6.54, 1.0, 1.0, 33.6, -6.54, 24.68, 1.0,
        /* D */ 1.0, 1.0, 1.0, 1.0, -6.54, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0, -6.54, 20.26, -14.03, 1.0, 1.0, 1.0,
        /* E */ 1.0, 44.94, 20.26, 33.6, 1.0, 1.0, -6.54, 20.26, 1.0, 1.0, 1.0, 1.0, 20.26, 20.26, 1.0, 20.26, 1.0, 1.0, -14.03, 1.0,
        /* F */ 1.0, 1.0, 13.34, 1.0, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0, 1.0, 1.0, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 33.601,
        /* G */ -7.49, 1.0, 1.0, -6.54, 1.0, 13.34, 1.0, -7.49, -7.49, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, 13.34, -7.49,
        /* H */ 1.0, 1.0, 1.0, 1.0, -9.37, -9.37, 1.0, 44.94, 24.68, 1.0, 1.0, 24.68, -1.88, 1.0, 1.0, 1.0, -6.54, 1.0, -1.88, 44.94,
        /* I */ 1.0, 1.0, 1.0, 44.94, 1.0, 1.0, 13.34, 1.0, -7.49, 20.26, 1.0, 1.0, -1.88, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0,
        /* K */ 1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, -7.49, 1.0, -7.49, 33.6, 1.0, -6.54, 24.64, 33.6, 1.0, 1.0, -7.49, 1.0, 1.0,
        /* L */ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 20.26, 33.6, 20.26, 1.0, 1.0, 1.0, 24.68, 1.0,
        /* M */ 13.34, 1.0, 1.0, 1.0, 1.0, 1.0, 58.28, 1.0, 1.0, 1.0, -1.88, 1.0, 44.94, -6.54, -6.54, 44.94, -1.88, 1.0, 1.0, 24.68,
        /* N */ 1.0, -1.88, 1.0, 1.0, -14.03, -14.03, 1.0, 44.94, 24.68, 1.0, 1.0, 1.0, -1.88, -6.54, 1.0, 1.0, -7.49, 1.0, -9.37, 1.0,
        /* P */ 20.26, -6.54, -6.54, 18.38, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, -6.54, 1.0, 20.26, 20.26, -6.54, 20.26, 1.0, 20.26, -1.88, 1.0,
        /* Q */ 1.0, -6.54, 20.26, 20.26, -6.54, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 20.26, 20.26, 1.0, 44.94, 1.0, -6.54, 1.0, -6.54,
        /* R */ 1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 20.26, 1.0, 1.0, 1.0, 1.0, 13.34, 20.26, 20.26, 58.28, 44.94, 1.0, 1.0, 58.28, -6.54,
        /* S */ 1.0, 33.6, 1.0, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 44.94, 20.26, 20.26, 20.26, 1.0, 1.0, 1.0, 1.0,
        /* T */ 1.0, 1.0, 1.0, 20.26, 13.34, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0, -6.54, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0,
        /* V */ 1.0, 1.0, -14.03, 1.0, 1.0, -7.49, 1.0, 1.0, -1.88, 1.0, 1.0, 1.0, 20.26, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, -6.54,
        /* W */ -14.03, 1.0, 1.0, 1.0, 1.0, -9.37, 24.68, 1.0, 1.0, 13.34, 24.68, 13.34, 1.0, 1.0, 1.0, 1.0, -14.03, -7.49, 1.0, 1.0,
        /* Y */ 24.68, 1.0, 24.68, -6.54, 1.0, -7.49, 13.34, 1.0, 1.0, 1.0, 44.94, 1.0, 13.34, 1.0, -15.91, 1.0, -7.49, 1.0, -9.37, 13.34,
    ];
}
