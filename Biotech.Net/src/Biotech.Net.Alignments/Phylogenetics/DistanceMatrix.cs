using System.Globalization;
using System.Text;

namespace Biotech.Net.Alignments.Phylogenetics;

/// <summary>Evolutionary distance estimators for aligned sequences.</summary>
public enum DistanceModel
{
    /// <summary>Proportion of differing sites, p.</summary>
    PDistance,

    /// <summary>
    /// Jukes–Cantor (1969): <c>d = −b·ln(1 − p/b)</c> with <c>b = 3/4</c> for nucleotides and <c>19/20</c> for proteins
    /// (the Poisson-corrected protein distance). Saturated pairs (p ≥ b) get +∞.
    /// </summary>
    JukesCantor,

    /// <summary>
    /// Kimura two-parameter (1980), nucleotides only: <c>d = −½·ln(1 − 2P − Q) − ¼·ln(1 − 2Q)</c> with P transitions and
    /// Q transversions per site. Saturated pairs get +∞.
    /// </summary>
    Kimura2P,
}

/// <summary>A symmetric matrix of pairwise distances between named taxa.</summary>
public sealed class DistanceMatrix
{
    private readonly double[,] _d;
    private readonly string[] _names;

    /// <summary>Creates a zero matrix for <paramref name="names"/>.</summary>
    public DistanceMatrix(IEnumerable<string> names)
    {
        ArgumentNullException.ThrowIfNull(names);
        _names = names.ToArray();
        if (_names.Distinct().Count() != _names.Length) throw new ArgumentException("Taxon names must be unique.", nameof(names));
        _d = new double[_names.Length, _names.Length];
    }

    /// <summary>Creates a matrix from a square table (only the lower triangle is read and mirrored).</summary>
    public DistanceMatrix(IEnumerable<string> names, double[,] distances) : this(names)
    {
        ArgumentNullException.ThrowIfNull(distances);
        if (distances.GetLength(0) != Count || distances.GetLength(1) != Count)
            throw new ArgumentException("The table must be square with one row per name.", nameof(distances));
        for (int i = 0; i < Count; i++)
            for (int j = 0; j < i; j++) this[i, j] = distances[i, j];
    }

    /// <summary>Taxon names in row order.</summary>
    public IReadOnlyList<string> Names => _names;

    /// <summary>Number of taxa.</summary>
    public int Count => _names.Length;

    /// <summary>Distance between taxa <paramref name="i"/> and <paramref name="j"/>; setting one mirrors it.</summary>
    public double this[int i, int j]
    {
        get => _d[i, j];
        set
        {
            if (i == j && value != 0) throw new ArgumentException("The diagonal of a distance matrix is zero.");
            _d[i, j] = value;
            _d[j, i] = value;
        }
    }

    /// <summary>Distance between two named taxa.</summary>
    public double this[string a, string b]
    {
        get => _d[IndexOf(a), IndexOf(b)];
        set => this[IndexOf(a), IndexOf(b)] = value;
    }

    /// <summary>Row index of <paramref name="name"/>.</summary>
    public int IndexOf(string name)
    {
        int i = Array.IndexOf(_names, name);
        return i >= 0 ? i : throw new KeyNotFoundException($"No taxon named '{name}'.");
    }

    /// <summary>
    /// Pairwise distances between the rows of a multiple alignment. Columns where either row has a gap (or an ambiguous
    /// symbol, for nucleotides) are skipped for that pair (pairwise deletion).
    /// </summary>
    public static DistanceMatrix FromAlignment(MultipleSequenceAlignment alignment, DistanceModel model = DistanceModel.PDistance)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        bool nucleotide = alignment.IsNucleotide;
        if (model == DistanceModel.Kimura2P && !nucleotide) throw new ArgumentException("Kimura 2-parameter distances need nucleotide sequences.");
        var dm = new DistanceMatrix(UniqueNames(alignment.Ids));
        for (int i = 0; i < alignment.Count; i++)
            for (int j = 0; j < i; j++)
                dm[i, j] = PairDistance(alignment[i], alignment[j], model, nucleotide);
        return dm;
    }

    private static IEnumerable<string> UniqueNames(IReadOnlyList<string> ids)
    {
        var seen = new Dictionary<string, int>();
        foreach (string id in ids)
        {
            int k = seen.TryGetValue(id, out int c) ? c + 1 : 1;
            seen[id] = k;
            yield return k == 1 ? id : $"{id}_{k}";
        }
    }

    /// <summary>Distance between two aligned rows of equal length under <paramref name="model"/>.</summary>
    public static double PairDistance(string a, string b, DistanceModel model, bool nucleotide)
    {
        if (a.Length != b.Length) throw new ArgumentException("Aligned rows must have the same length.");
        int sites = 0, diff = 0, transitions = 0, transversions = 0;
        for (int k = 0; k < a.Length; k++)
        {
            char x = char.ToUpperInvariant(a[k]), y = char.ToUpperInvariant(b[k]);
            if (x is '-' or '.' || y is '-' or '.') continue;
            if (nucleotide)
            {
                if (x == 'U') x = 'T';
                if (y == 'U') y = 'T';
                if ("ACGT".IndexOf(x) < 0 || "ACGT".IndexOf(y) < 0) continue;
            }
            sites++;
            if (x == y) continue;
            diff++;
            bool purines = x is 'A' or 'G' && y is 'A' or 'G';
            bool pyrimidines = x is 'C' or 'T' && y is 'C' or 'T';
            if (purines || pyrimidines) transitions++; else transversions++;
        }
        if (sites == 0) return double.NaN;
        double p = (double)diff / sites;
        switch (model)
        {
            case DistanceModel.PDistance:
                return p;
            case DistanceModel.JukesCantor:
            {
                double bConst = nucleotide ? 0.75 : 19.0 / 20.0;
                double arg = 1 - p / bConst;
                return arg <= 0 ? double.PositiveInfinity : p == 0 ? 0 : -bConst * Math.Log(arg);
            }
            case DistanceModel.Kimura2P:
            {
                double P = (double)transitions / sites, Q = (double)transversions / sites;
                double a1 = 1 - 2 * P - Q, a2 = 1 - 2 * Q;
                if (a1 <= 0 || a2 <= 0) return double.PositiveInfinity;
                return p == 0 ? 0 : -0.5 * Math.Log(a1) - 0.25 * Math.Log(a2);
            }
            default:
                throw new ArgumentOutOfRangeException(nameof(model));
        }
    }

    /// <summary>Square PHYLIP format: the taxon count, then one line per taxon (name padded to 10).</summary>
    public string ToPhylip()
    {
        var sb = new StringBuilder();
        sb.Append(Count.ToString(CultureInfo.InvariantCulture)).AppendLine();
        for (int i = 0; i < Count; i++)
        {
            string name = _names[i].Length > 10 ? _names[i][..10] : _names[i];
            sb.Append(name.PadRight(10));
            for (int j = 0; j < Count; j++) sb.Append(' ').Append(_d[i, j].ToString("0.000000", CultureInfo.InvariantCulture));
            sb.AppendLine();
        }
        return sb.ToString();
    }

    /// <inheritdoc />
    public override string ToString() => ToPhylip();
}
