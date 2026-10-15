namespace Biotech.Net.Proteomics;

/// <summary>
/// A numeric value per amino acid (hydropathy, hydrophilicity, flexibility, ...) with ExPASy ProtScale-style
/// sliding-window profiles. Lookups are a 256-entry table, so profiles cost one pass over the sequence.
/// </summary>
public sealed class AminoAcidScale
{
    private readonly double[] _table = new double[256];
    private readonly bool[] _known = new bool[256];

    /// <summary>Creates a scale from one value per one-letter code.</summary>
    public AminoAcidScale(string name, IReadOnlyDictionary<char, double> values)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentNullException.ThrowIfNull(values);
        Name = name;
        foreach (var (c, v) in values)
        {
            Set(char.ToUpperInvariant(c), v);
            Set(char.ToLowerInvariant(c), v);
        }
    }

    private AminoAcidScale(string name, string residues, double[] values)
    {
        Name = name;
        for (int i = 0; i < residues.Length; i++)
        {
            Set(residues[i], values[i]);
            Set(char.ToLowerInvariant(residues[i]), values[i]);
        }
    }

    private void Set(char c, double v)
    {
        if (c >= 256) throw new ArgumentException($"'{c}' is not an ASCII residue code.");
        _table[c] = v;
        _known[c] = true;
    }

    /// <summary>Scale name.</summary>
    public string Name { get; }

    /// <summary>Value for a residue (0 when the residue is not on the scale).</summary>
    public double this[char residue] => residue < 256 ? _table[residue] : 0;

    /// <summary>Value for an ASCII residue byte (0 when the residue is not on the scale).</summary>
    public double this[byte residue] => _table[residue];

    /// <summary>True when the scale defines a value for <paramref name="residue"/>.</summary>
    public bool Contains(char residue) => residue < 256 && _known[residue];

    /// <summary>
    /// Sum of the scale over the sequence divided by its full length (unknown residues contribute 0). With
    /// <see cref="KyteDoolittle"/> this is GRAVY, as in ExPASy ProtParam and Biopython <c>gravy()</c>.
    /// </summary>
    public double Average(ReadOnlySpan<byte> residues)
    {
        if (residues.IsEmpty) return 0;
        double s = 0;
        foreach (byte b in residues) s += _table[b];
        return s / residues.Length;
    }

    /// <summary>
    /// Sliding-window profile as computed by ExPASy ProtScale / Biopython <c>protein_scale</c>: the window centre has
    /// weight 1 and the weights fall linearly to <paramref name="edge"/> at both ends; each score is the weighted
    /// sum divided by the sum of the weights. Returns <c>Length − window + 1</c> values; value <c>i</c> belongs to
    /// residue <c>i + window / 2</c>.
    /// </summary>
    /// <param name="residues">Amino acids.</param>
    /// <param name="window">Odd or even window length ≥ 1.</param>
    /// <param name="edge">Relative weight of the window ends, in [0, 1]; 1 gives a plain moving average.</param>
    public double[] Profile(ReadOnlySpan<byte> residues, int window, double edge = 1.0)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(window, 1);
        if (edge is < 0 or > 1) throw new ArgumentOutOfRangeException(nameof(edge));
        int n = residues.Length - window + 1;
        if (n <= 0) return [];
        var output = new double[n];
        int half = window / 2;
        if (edge == 1.0 && window % 2 == 1)
        {
            // Plain moving average: O(n) running sum.
            double s = 0;
            for (int i = 0; i < window; i++) s += _table[residues[i]];
            output[0] = s / window;
            for (int i = 1; i < n; i++)
            {
                s += _table[residues[i + window - 1]] - _table[residues[i - 1]];
                output[i] = s / window;
            }
            return output;
        }

        // One-sided linear weights, exactly as ProtScale: w[j] = edge + j·2(1−edge)/(window−1).
        Span<double> w = stackalloc double[Math.Max(half, 1)];
        double unit = window > 1 ? 2 * (1.0 - edge) / (window - 1) : 0;
        double wsum = 1;
        for (int j = 0; j < half; j++) { w[j] = edge + unit * j; wsum += 2 * w[j]; }
        for (int i = 0; i < n; i++)
        {
            var sub = residues.Slice(i, window);
            double score = 0;
            for (int j = 0; j < half; j++) score += w[j] * (_table[sub[j]] + _table[sub[window - j - 1]]);
            score += _table[sub[half]];
            output[i] = score / wsum;
        }
        return output;
    }

    /// <inheritdoc/>
    public override string ToString() => Name;

    private const string Order = "ARNDCQEGHILKMFPSTWYV";

    /// <summary>Kyte &amp; Doolittle hydropathy (J. Mol. Biol. 157:105, 1982).</summary>
    public static AminoAcidScale KyteDoolittle { get; } = new("Kyte-Doolittle", Order,
        [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]);

    /// <summary>Hopp &amp; Woods hydrophilicity (PNAS 78:3824, 1981).</summary>
    public static AminoAcidScale HoppWoods { get; } = new("Hopp-Woods", Order,
        [-0.5, 3.0, 0.2, 3.0, -1.0, 0.2, 3.0, 0.0, -0.5, -1.8, -1.8, 3.0, -1.3, -2.5, 0.0, 0.3, -0.4, -3.4, -2.3, -1.5]);

    /// <summary>Eisenberg normalised consensus hydrophobicity (J. Mol. Biol. 179:125, 1984).</summary>
    public static AminoAcidScale Eisenberg { get; } = new("Eisenberg", Order,
        [0.62, -2.53, -0.78, -0.90, 0.29, -0.85, -0.74, 0.48, -0.40, 1.38, 1.06, -1.50, 0.64, 1.19, 0.12, -0.18, -0.05, 0.81, 0.26, 1.08]);

    /// <summary>Vihinen et al. normalised flexibility (B-values), Proteins 19:141 (1994).</summary>
    public static AminoAcidScale Flexibility { get; } = new("Vihinen flexibility", Order,
        [0.984, 1.008, 1.048, 1.068, 0.906, 1.037, 1.094, 1.031, 0.950, 0.927, 0.935, 1.102, 0.952, 0.915, 1.049, 1.046, 0.997, 0.904, 0.929, 0.931]);

    /// <summary>Emini surface fractional probability (as in ProtScale / Biopython).</summary>
    public static AminoAcidScale EminiSurface { get; } = new("Emini surface", Order,
        [0.815, 1.475, 1.296, 1.283, 0.394, 1.348, 1.445, 0.714, 1.180, 0.603, 0.603, 1.545, 0.714, 0.695, 1.236, 1.115, 1.184, 0.808, 1.089, 0.606]);

    /// <summary>Janin interior-to-surface transfer energy.</summary>
    public static AminoAcidScale Janin { get; } = new("Janin", Order,
        [0.28, -1.14, -0.55, -0.52, 0.97, -0.69, -1.01, 0.43, -0.31, 0.60, 0.60, -1.62, 0.43, 0.46, -0.42, -0.19, -0.32, 0.29, -0.15, 0.60]);
}
