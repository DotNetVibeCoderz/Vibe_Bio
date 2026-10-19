using Biotech.Net.Core;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>
/// A model (design) matrix for linear and generalised linear models: one row per sample, one column per
/// coefficient. Built from sample annotations with treatment (reference-level) coding, as R's
/// <c>model.matrix(~ batch + condition)</c> does.
/// </summary>
public sealed class DesignMatrix
{
    /// <summary>Wraps a row-major samples × coefficients matrix.</summary>
    public DesignMatrix(int samples, IReadOnlyList<string> coefficientNames, double[] values)
    {
        ArgumentNullException.ThrowIfNull(coefficientNames);
        ArgumentNullException.ThrowIfNull(values);
        if (values.Length != samples * coefficientNames.Count)
            throw new ArgumentException($"Expected {samples}×{coefficientNames.Count} values.", nameof(values));
        Samples = samples;
        CoefficientNames = coefficientNames;
        Values = values;
    }

    /// <summary>Number of samples (rows).</summary>
    public int Samples { get; }

    /// <summary>Number of coefficients (columns).</summary>
    public int Coefficients => CoefficientNames.Count;

    /// <summary>Coefficient names, e.g. <c>Intercept</c>, <c>condition_treated_vs_control</c>.</summary>
    public IReadOnlyList<string> CoefficientNames { get; }

    /// <summary>Row-major values.</summary>
    public double[] Values { get; }

    /// <summary>
    /// The categorical variables the matrix was built from (empty for matrices supplied directly). Used to resolve
    /// level-based contrasts such as <c>("condition", "treated", "untreated")</c>.
    /// </summary>
    public IReadOnlyList<DesignFactor> Factors { get; private init; } = [];

    /// <summary>Element access.</summary>
    public double this[int sample, int coefficient] => Values[sample * Coefficients + coefficient];

    /// <summary>Index of a coefficient by name.</summary>
    public int IndexOf(string coefficientName)
    {
        for (int i = 0; i < Coefficients; i++) if (CoefficientNames[i] == coefficientName) return i;
        throw new KeyNotFoundException($"No coefficient '{coefficientName}'. Available: {string.Join(", ", CoefficientNames)}.");
    }

    /// <summary>An intercept-only design.</summary>
    public static DesignMatrix Intercept(int samples)
    {
        var v = new double[samples];
        Array.Fill(v, 1.0);
        return new DesignMatrix(samples, ["Intercept"], v);
    }

    /// <summary>
    /// Two-group design over <paramref name="samples"/> columns: intercept plus an indicator of
    /// <paramref name="treatment"/>. Samples in neither group are not allowed (subset the data first).
    /// </summary>
    public static DesignMatrix TwoGroup(int samples, ReadOnlySpan<int> reference, ReadOnlySpan<int> treatment, string coefficientName = "group_B_vs_A")
    {
        var v = new double[samples * 2];
        var seen = new bool[samples];
        foreach (int c in reference) { seen[c] = true; v[c * 2] = 1; }
        foreach (int c in treatment) { seen[c] = true; v[c * 2] = 1; v[c * 2 + 1] = 1; }
        for (int i = 0; i < samples; i++)
            if (!seen[i]) throw new ArgumentException($"Sample {i} is in neither group.");
        return new DesignMatrix(samples, ["Intercept", coefficientName], v);
    }

    /// <summary>
    /// Builds <c>~ covariates + factor</c> from the column annotations of <paramref name="data"/> with treatment
    /// coding. The factor of interest uses <paramref name="reference"/> as its baseline and yields one coefficient per
    /// other level named <c>{factor}_{level}_vs_{reference}</c> (DESeq2's naming); categorical covariates use their
    /// alphabetically first level as baseline and are named <c>{covariate}_{level}</c>.
    /// </summary>
    public static DesignMatrix FromAnnotations(Dataset data, string factor, string reference, params string[] covariates)
    {
        ArgumentNullException.ThrowIfNull(data);
        int n = data.Columns;
        var names = new List<string> { "Intercept" };
        var cols = new List<double[]> { Enumerable.Repeat(1.0, n).ToArray() };
        var factors = new List<DesignFactor>();
        foreach (var cov in covariates)
        {
            var ann = Annotation(data, cov);
            var levels = ann.Distinct().Order(StringComparer.Ordinal).ToArray();
            var at = new List<int> { -1 };
            for (int l = 1; l < levels.Length; l++)
            {
                at.Add(names.Count);
                names.Add($"{cov}_{levels[l]}");
                cols.Add(ann.Select(a => a == levels[l] ? 1.0 : 0.0).ToArray());
            }
            factors.Add(new DesignFactor(cov, levels, at));
        }
        var f = Annotation(data, factor);
        if (!f.Contains(reference)) throw new ArgumentException($"Reference level '{reference}' does not occur in '{factor}'.");
        var flevels = new List<string> { reference };
        var fat = new List<int> { -1 };
        foreach (var level in f.Distinct().Order(StringComparer.Ordinal))
        {
            if (level == reference) continue;
            flevels.Add(level);
            fat.Add(names.Count);
            names.Add($"{factor}_{level}_vs_{reference}");
            cols.Add(f.Select(a => a == level ? 1.0 : 0.0).ToArray());
        }
        factors.Add(new DesignFactor(factor, flevels, fat));
        int p = names.Count;
        var v = new double[n * p];
        for (int j = 0; j < p; j++) for (int i = 0; i < n; i++) v[i * p + j] = cols[j][i];
        var d = new DesignMatrix(n, names, v) { Factors = factors };
        if (d.Rank() < p) throw new ArgumentException("The design matrix is not of full rank (confounded factors?).");
        return d;
    }

    /// <summary>Numerical rank (via Gram matrix Cholesky with pivot tolerance).</summary>
    public int Rank()
    {
        int p = Coefficients;
        var g = Gram();
        // Gaussian elimination rank.
        int rank = 0;
        var a = (double[])g.Clone();
        var used = new bool[p];
        double scale = 0;
        for (int i = 0; i < p; i++) scale = Math.Max(scale, Math.Abs(a[i * p + i]));
        for (int col = 0; col < p; col++)
        {
            int piv = -1; double best = 1e-10 * Math.Max(scale, 1e-300);
            for (int r = 0; r < p; r++) if (!used[r] && Math.Abs(a[r * p + col]) > best) { best = Math.Abs(a[r * p + col]); piv = r; }
            if (piv < 0) continue;
            used[piv] = true; rank++;
            for (int r = 0; r < p; r++)
            {
                if (r == piv) continue;
                double fct = a[r * p + col] / a[piv * p + col];
                for (int c = 0; c < p; c++) a[r * p + c] -= fct * a[piv * p + c];
            }
        }
        return rank;
    }

    /// <summary>The leverages (diagonal of the hat matrix X(XᵀX)⁻¹Xᵀ).</summary>
    public double[] Leverages()
    {
        int p = Coefficients;
        var inv = NumericUtil.InvertSpd(Gram(), p) ?? throw new InvalidOperationException("Design is not of full rank.");
        var h = new double[Samples];
        for (int i = 0; i < Samples; i++)
        {
            double s = 0;
            for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++)
                    s += this[i, a] * inv[a * p + b] * this[i, b];
            h[i] = s;
        }
        return h;
    }

    internal double[] Gram()
    {
        int p = Coefficients;
        var g = new double[p * p];
        for (int i = 0; i < Samples; i++)
            for (int a = 0; a < p; a++)
            {
                double xa = Values[i * p + a];
                if (xa == 0) continue;
                for (int b = 0; b < p; b++) g[a * p + b] += xa * Values[i * p + b];
            }
        return g;
    }

    private static string[] Annotation(Dataset data, string key)
        => data.ColumnAnnotations.TryGetValue(key, out var a) ? a : throw new KeyNotFoundException($"No column annotation '{key}'.");

    /// <inheritdoc />
    public override string ToString() => $"Design {Samples}×{Coefficients}: {string.Join(", ", CoefficientNames)}";
}

/// <summary>
/// A categorical variable of a <see cref="DesignMatrix"/> built with treatment coding.
/// </summary>
/// <param name="Name">Annotation name, e.g. <c>condition</c>.</param>
/// <param name="Levels">Levels with the reference level first.</param>
/// <param name="Columns">Coefficient index of each level (−1 for the reference level).</param>
public sealed record DesignFactor(string Name, IReadOnlyList<string> Levels, IReadOnlyList<int> Columns);
