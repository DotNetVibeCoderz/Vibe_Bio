using Biotech.Net.Native;

namespace Biotech.Net.Proteomics;

/// <summary>Chou–Fasman conformational parameters of one residue (×100 in the original tables).</summary>
/// <param name="Helix">P(α).</param>
/// <param name="Sheet">P(β).</param>
/// <param name="Turn">P(turn).</param>
public readonly record struct ChouFasmanParameters(double Helix, double Sheet, double Turn);

/// <summary>
/// Chou–Fasman secondary-structure prediction (Chou &amp; Fasman, Adv. Enzymol. 47:45, 1978).
/// </summary>
/// <remarks>
/// A simplified but faithful version of the classic rules: helix nuclei are 6-residue windows with at least
/// 4 helix formers (P(α) ≥ 1), sheet nuclei are 5-residue windows with at least 3 sheet formers (P(β) ≥ 1); both are
/// extended while the average of the 4-residue end window stays ≥ 1. A helix is kept when its mean P(α) ≥ 1.03 and
/// exceeds mean P(β); a sheet when mean P(β) ≥ 1.05 and exceeds mean P(α). Overlaps go to the class with the higher
/// mean over the overlap. Turns are assigned to still-unassigned 4-residue windows whose mean P(turn) exceeds
/// 1.0 and both other means (the bend-frequency product test is omitted). Accuracy is the method's historical
/// ~50–60 %; use it as a quick profile, not a structure.
/// </remarks>
public static class ChouFasman
{
    private const string Order = "ARNDCQEGHILKMFPSTWYV";
    private static readonly int[] Pa = [142, 98, 67, 101, 70, 111, 151, 57, 100, 108, 121, 114, 145, 113, 57, 77, 83, 108, 69, 106];
    private static readonly int[] Pb = [83, 93, 89, 54, 119, 110, 37, 75, 87, 160, 130, 74, 105, 138, 55, 75, 119, 137, 147, 170];
    private static readonly int[] Pt = [66, 95, 156, 146, 119, 98, 74, 156, 95, 47, 59, 101, 60, 60, 152, 143, 96, 96, 114, 50];
    private static readonly double[] TableA = Build(Pa), TableB = Build(Pb), TableT = Build(Pt);

    private static double[] Build(int[] v)
    {
        var t = new double[256];
        Array.Fill(t, 1.0); // unknown residues are neutral
        for (int i = 0; i < Order.Length; i++)
        {
            t[Order[i]] = v[i] / 100.0;
            t[char.ToLowerInvariant(Order[i])] = v[i] / 100.0;
        }
        return t;
    }

    /// <summary>Parameters of one residue (all 1.0 for non-standard residues).</summary>
    public static ChouFasmanParameters Parameters(char residue)
        => residue < 256 ? new(TableA[residue], TableB[residue], TableT[residue]) : new(1, 1, 1);

    /// <summary>
    /// Predicts one state per residue: <c>H</c> (α-helix), <c>E</c> (β-strand), <c>T</c> (turn) or <c>C</c> (coil).
    /// </summary>
    public static string Predict(ReadOnlySpan<byte> residues)
    {
        int n = residues.Length;
        var a = new double[n];
        var b = new double[n];
        var t = new double[n];
        for (int i = 0; i < n; i++) { a[i] = TableA[residues[i]]; b[i] = TableB[residues[i]]; t[i] = TableT[residues[i]]; }

        var helix = FindRegions(a, b, nucleus: 6, formers: 4, minMean: 1.03);
        var sheet = FindRegions(b, a, nucleus: 5, formers: 3, minMean: 1.05);

        var state = new char[n];
        Array.Fill(state, 'C');
        for (int i = 0; i < n; i++)
        {
            if (helix[i] && sheet[i]) continue; // resolved below
            if (helix[i]) state[i] = 'H';
            else if (sheet[i]) state[i] = 'E';
        }
        // Resolve overlaps region by region.
        for (int i = 0; i < n;)
        {
            if (!(helix[i] && sheet[i])) { i++; continue; }
            int j = i;
            double sa = 0, sb = 0;
            while (j < n && helix[j] && sheet[j]) { sa += a[j]; sb += b[j]; j++; }
            char c = sa >= sb ? 'H' : 'E';
            for (int k = i; k < j; k++) state[k] = c;
            i = j;
        }
        // Turns in unassigned stretches.
        for (int i = 0; i + 4 <= n; i++)
        {
            double mt = 0, ma = 0, mb = 0;
            bool free = true;
            for (int k = i; k < i + 4; k++) { mt += t[k]; ma += a[k]; mb += b[k]; free &= state[k] is 'C' or 'T'; }
            if (free && mt / 4 > 1.0 && mt > ma && mt > mb)
                for (int k = i; k < i + 4; k++) state[k] = 'T';
        }
        return new string(state);
    }

    /// <summary>Predicts states for a one-letter string.</summary>
    public static string Predict(string residues) => Predict(System.Text.Encoding.ASCII.GetBytes(residues));

    private static bool[] FindRegions(double[] p, double[] other, int nucleus, int formers, double minMean)
    {
        int n = p.Length;
        var mark = new bool[n];
        for (int i = 0; i + nucleus <= n; i++)
        {
            int count = 0;
            for (int k = i; k < i + nucleus; k++) if (p[k] >= 1.0) count++;
            if (count < formers) continue;
            int start = i, end = i + nucleus;
            while (end < n && Mean(p, end - 3, end + 1) >= 1.0) end++;
            while (start > 0 && Mean(p, start - 1, start + 3) >= 1.0) start--;
            double mp = Mean(p, start, end), mo = Mean(other, start, end);
            if (mp >= minMean && mp > mo)
                for (int k = start; k < end; k++) mark[k] = true;
        }
        return mark;
    }

    private static double Mean(double[] v, int from, int to)
    {
        double s = 0;
        for (int i = from; i < to; i++) s += v[i];
        return s / (to - from);
    }
}

/// <summary>A predicted transmembrane helix.</summary>
/// <param name="Start">0-based start (inclusive).</param>
/// <param name="End">0-based end (exclusive).</param>
/// <param name="MaxWindowHydropathy">Highest window-averaged Kyte–Doolittle score inside the segment.</param>
public readonly record struct TransmembraneSegment(int Start, int End, double MaxWindowHydropathy)
{
    /// <summary>Segment length.</summary>
    public int Length => End - Start;
}

/// <summary>
/// Transmembrane-helix detection by Kyte–Doolittle hydropathy (Kyte &amp; Doolittle 1982): a window of 19 residues
/// with mean hydropathy ≥ 1.6 marks a likely membrane-spanning helix. Overlapping hit windows are merged.
/// </summary>
public static class TransmembranePredictor
{
    /// <summary>Finds hydrophobic segments; uses the native <see cref="BioKernels.Hydropathy"/> kernel.</summary>
    public static IReadOnlyList<TransmembraneSegment> Predict(ReadOnlySpan<byte> residues, int window = 19, double threshold = 1.6)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(window, 1);
        int n = Math.Max(0, residues.Length - window + 1);
        var segments = new List<TransmembraneSegment>();
        if (n == 0) return segments;
        var profile = new double[n];
        BioKernels.Hydropathy(residues, window, profile);
        int i = 0;
        while (i < n)
        {
            if (profile[i] < threshold) { i++; continue; }
            int j = i;
            double max = double.NegativeInfinity;
            while (j < n && profile[j] >= threshold) { max = Math.Max(max, profile[j]); j++; }
            segments.Add(new(i, j - 1 + window, max));
            i = j;
        }
        return segments;
    }

    /// <summary>Finds hydrophobic segments in a one-letter string.</summary>
    public static IReadOnlyList<TransmembraneSegment> Predict(string residues, int window = 19, double threshold = 1.6)
        => Predict(System.Text.Encoding.ASCII.GetBytes(residues), window, threshold);
}

/// <summary>Result of the signal-peptide heuristic.</summary>
/// <param name="HasSignalPeptide">True when an n/h/c-region architecture with a valid cleavage site was found.</param>
/// <param name="CleavageSite">Length of the signal peptide (cleavage between residues <c>CleavageSite − 1</c> and <c>CleavageSite</c>, 0-based), or −1.</param>
/// <param name="HRegionStart">0-based start of the most hydrophobic 8-residue window in the first 30 residues.</param>
/// <param name="HRegionHydropathy">Mean Kyte–Doolittle hydropathy of that window.</param>
/// <param name="Score">Heuristic score (higher is more signal-like); not a probability.</param>
public readonly record struct SignalPeptidePrediction(bool HasSignalPeptide, int CleavageSite, int HRegionStart, double HRegionHydropathy, double Score);

/// <summary>
/// A rule-based Sec signal-peptide detector following von Heijne's tripartite model (J. Mol. Biol. 184:99, 1985):
/// a hydrophobic h-region (8-residue mean Kyte–Doolittle ≥ 1.5 starting within the first 22 residues, free of
/// D/E/K/R), optionally preceded by a positively charged n-region, followed 3–12 residues later by a cleavage
/// site obeying the (−3, −1) small-residue rule at positions 15–35.
/// </summary>
/// <remarks>This is a transparent heuristic for screening, not a replacement for SignalP-class predictors.</remarks>
public static class SignalPeptidePredictor
{
    /// <summary>Applies the heuristic to the N-terminus of <paramref name="residues"/>.</summary>
    public static SignalPeptidePrediction Predict(ReadOnlySpan<byte> residues)
    {
        var seq = AminoAcids.ToUpperBytes(residues.Length > 60 ? residues[..60] : residues);
        if (seq.Length < 20) return new(false, -1, -1, double.NaN, 0);
        var kd = AminoAcidScale.KyteDoolittle;

        // h-region: most hydrophobic 8-mer starting in [1, 22].
        int hStart = -1;
        double hBest = double.NegativeInfinity;
        for (int s = 1; s <= Math.Min(22, seq.Length - 8); s++)
        {
            double m = kd.Average(seq.AsSpan(s, 8));
            if (m > hBest) { hBest = m; hStart = s; }
        }
        int hEnd = hStart + 8;
        bool charged = false;
        for (int k = hStart; k < hEnd; k++) charged |= seq[k] is (byte)'D' or (byte)'E' or (byte)'K' or (byte)'R';
        bool nPositive = false;
        for (int k = 0; k < hStart; k++) nPositive |= seq[k] is (byte)'K' or (byte)'R';

        // c-region: best (−3, −1) site.
        int best = -1;
        double bestSite = 0;
        for (int c = Math.Max(15, hEnd + 3); c <= Math.Min(35, Math.Min(hEnd + 12, seq.Length - 1)); c++)
        {
            double w1 = MinusOne(seq[c - 1]), w3 = MinusThree(seq[c - 3]);
            if (w1 <= 0 || w3 <= 0) continue;
            if (seq[c] == 'P' || seq[c - 2] == 'P') continue;
            double site = w1 + w3 + (IsBulky(seq[c - 2]) ? 0.5 : 0);
            if (site > bestSite) { bestSite = site; best = c; }
        }

        double score = hBest + (nPositive ? 0.5 : 0) + bestSite - (charged ? 2 : 0);
        bool has = hBest >= 1.5 && !charged && best > 0;
        return new(has, has ? best : -1, hStart, hBest, score);
    }

    /// <summary>Applies the heuristic to a one-letter string.</summary>
    public static SignalPeptidePrediction Predict(string residues) => Predict(System.Text.Encoding.ASCII.GetBytes(residues));

    private static double MinusOne(byte r) => r switch
    {
        (byte)'A' => 3, (byte)'G' => 2, (byte)'S' => 2, (byte)'C' => 1.5, (byte)'T' => 1, _ => 0,
    };

    private static double MinusThree(byte r) => r switch
    {
        (byte)'A' => 2, (byte)'V' => 1.5, (byte)'S' or (byte)'G' or (byte)'C' or (byte)'T' => 1, (byte)'I' or (byte)'L' => 0.75, _ => 0,
    };

    private static bool IsBulky(byte r) => r is (byte)'F' or (byte)'Y' or (byte)'W' or (byte)'H' or (byte)'L' or (byte)'I' or (byte)'M' or (byte)'Q' or (byte)'E';
}
