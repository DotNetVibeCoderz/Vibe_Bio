using System.Globalization;
using System.IO.Compression;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>
/// Loads the R reference outputs in Fixtures/DESeq2 (produced by <c>make_reference.R</c> with R 4.6.1, DESeq2 1.52.0 and
/// apeglm 1.34.0 on the Bioconductor <c>pasilla</c> data). Set <c>BIOTECH_DESEQ2_REFERENCE_DIR</c> to a directory written
/// by the script in "full" mode to compare against the whole dataset instead of the committed every-7th-gene subset.
/// </summary>
internal static class DESeq2Fixtures
{
    public static string Directory =>
        Environment.GetEnvironmentVariable("BIOTECH_DESEQ2_REFERENCE_DIR") is { Length: > 0 } d
            ? d
            : Path.Combine(AppContext.BaseDirectory, "Fixtures", "DESeq2");

    private static readonly Lazy<Dataset> PasillaLazy = new(LoadPasilla);
    private static readonly Lazy<Dataset> AirwayLazy = new(LoadAirway);
    private static readonly Dictionary<string, Dictionary<string, double[]>> Tables = new();
    private static readonly Lazy<Dictionary<(string, string), double>> ScalarsLazy = new(LoadScalars);

    /// <summary>The pasilla counts (genes × 7 samples) with column annotations <c>condition</c> and <c>type</c>.</summary>
    public static Dataset Pasilla => PasillaLazy.Value;

    /// <summary>The airway counts (genes × 8 samples) with column annotations <c>cell</c> and <c>dex</c>.</summary>
    public static Dataset Airway => AirwayLazy.Value;

    private static Dataset LoadAirway() => LoadCounts("airway", "airway_coldata.tsv", ["cell", "dex"]);

    private static Dataset LoadPasilla() => LoadCounts("pasilla", "pasilla_coldata.tsv", ["condition", "type"]);

    private static Dataset LoadCounts(string name, string coldataFile, string[] annotations)
    {
        string file = Path.Combine(Directory, name + "_subset_counts.tsv");
        if (!File.Exists(file)) file = Path.Combine(Directory, name + "_full_counts.tsv");
        var lines = File.ReadAllLines(file);
        var header = lines[0].Split('\t');
        var samples = header.Skip(1).ToArray();
        var genes = new List<string>();
        var values = new List<double>();
        foreach (var line in lines.Skip(1))
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            genes.Add(f[0]);
            for (int j = 1; j < f.Length; j++) values.Add(double.Parse(f[j], CultureInfo.InvariantCulture));
        }
        var d = new Dataset(name, genes, samples, values.ToArray());
        var coldata = File.ReadAllLines(Path.Combine(Directory, coldataFile)).Skip(1)
            .Where(l => l.Length > 0).Select(l => l.Split('\t')).ToDictionary(f => f[0], f => f);
        for (int a = 0; a < annotations.Length; a++)
            d.ColumnAnnotations[annotations[a]] = samples.Select(s => coldata[s][a + 1]).ToArray();
        return d;
    }

    /// <summary>One reference table (columns by name; R's NA becomes NaN).</summary>
    public static Dictionary<string, double[]> Table(string scenario)
    {
        lock (Tables)
        {
            if (Tables.TryGetValue(scenario, out var t)) return t;
            using var fs = File.OpenRead(Path.Combine(Directory, scenario + ".csv.gz"));
            using var gz = new GZipStream(fs, CompressionMode.Decompress);
            using var reader = new StreamReader(gz, Encoding.UTF8);
            var header = reader.ReadLine()!.Split(',');
            var cols = header.Select(_ => new List<double>()).ToArray();
            string? line;
            while ((line = reader.ReadLine()) is not null)
            {
                if (line.Length == 0) continue;
                var f = line.Split(',');
                for (int c = 1; c < f.Length; c++)
                {
                    var s = f[c].Trim();
                    cols[c].Add(s == "NA" ? double.NaN : double.Parse(s, CultureInfo.InvariantCulture));
                }
            }
            t = new Dictionary<string, double[]>();
            for (int c = 1; c < header.Length; c++) t[header[c]] = cols[c].ToArray();
            Tables[scenario] = t;
            return t;
        }
    }

    public static double Scalar(string scenario, string key) => ScalarsLazy.Value[(scenario, key)];

    private static Dictionary<(string, string), double> LoadScalars()
    {
        var res = new Dictionary<(string, string), double>();
        foreach (var line in File.ReadAllLines(Path.Combine(Directory, "scalars.csv")).Skip(1))
        {
            if (line.Length == 0) continue;
            var f = line.Split(',');
            var s = f[2].Trim();
            res[(f[0], f[1])] = s == "NA" ? double.NaN : double.Parse(s, CultureInfo.InvariantCulture);
        }
        return res;
    }
}

/// <summary>Collects per-column agreement statistics and failures for a parity check.</summary>
internal sealed class ParityReport(string scenario)
{
    private readonly List<string> _lines = [];
    private readonly List<string> _failures = [];

    public IReadOnlyList<string> Lines => _lines;

    /// <summary>
    /// Compares two columns: NaN must coincide, and |actual − expected| ≤ abs + rel·|expected| elsewhere (a looser
    /// <paramref name="outlierRel"/> / <paramref name="outlierAbs"/> may be granted to at most <paramref name="maxOutliers"/> genes).
    /// </summary>
    public void Check(string column, double[] expected, double[] actual, double rel, double abs = 0,
        int maxOutliers = 0, double outlierRel = 0, double? outlierAbs = null)
    {
        double oabs = outlierAbs ?? abs;
        if (expected.Length != actual.Length) { _failures.Add($"{column}: length {actual.Length} vs {expected.Length}"); return; }
        int nanMismatch = 0, compared = 0, outliers = 0, bad = 0, worst = -1;
        double maxAbs = 0, maxRel = 0;
        for (int i = 0; i < expected.Length; i++)
        {
            double e = expected[i], a = actual[i];
            if (double.IsNaN(e) || double.IsNaN(a))
            {
                if (double.IsNaN(e) != double.IsNaN(a)) nanMismatch++;
                continue;
            }
            compared++;
            double d = Math.Abs(a - e);
            double r = e == 0 ? (d == 0 ? 0 : double.PositiveInfinity) : d / Math.Abs(e);
            if (d > maxAbs) maxAbs = d;
            if (r > maxRel && d > abs) { maxRel = r; worst = i; }
            if (d > abs + rel * Math.Abs(e))
            {
                if (d <= oabs + outlierRel * Math.Abs(e)) outliers++;
                else bad++;
            }
        }
        string line = $"{scenario,-18} {column,-34} n={compared,6} maxAbs={maxAbs:E2} maxRel={maxRel:E2} beyondTol={bad + outliers} NaNmismatch={nanMismatch}" +
                      (worst >= 0 ? $" worst@{worst}" : "");
        _lines.Add(line);
        if (nanMismatch > 0 || bad > 0 || outliers > maxOutliers)
            _failures.Add(line);
    }

    /// <summary>Compares integer-valued / boolean columns exactly, allowing <paramref name="maxMismatch"/> differences.</summary>
    public void CheckExact(string column, double[] expected, double[] actual, int maxMismatch = 0)
    {
        int mismatch = 0, first = -1;
        for (int i = 0; i < expected.Length; i++)
        {
            bool same = (double.IsNaN(expected[i]) && double.IsNaN(actual[i])) || expected[i] == actual[i];
            if (!same) { mismatch++; if (first < 0) first = i; }
        }
        string line = $"{scenario,-18} {column,-34} exact mismatches={mismatch}" + (first >= 0 ? $" first@{first} (R {expected[first]} vs {actual[first]})" : "");
        _lines.Add(line);
        if (mismatch > maxMismatch) _failures.Add(line);
    }

    public void AssertOk(ITestOutputHelper output)
    {
        foreach (var l in _lines) output.WriteLine(l);
        var path = Environment.GetEnvironmentVariable("BIOTECH_DESEQ2_REPORT");
        if (!string.IsNullOrEmpty(path)) lock (typeof(ParityReport)) File.AppendAllLines(path, _lines);
        Assert.True(_failures.Count == 0, "Parity failures:\n" + string.Join("\n", _failures));
    }
}
