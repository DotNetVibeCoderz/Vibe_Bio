using System.Globalization;
using System.Text;

namespace Biotech.Net.Core;

/// <summary>
/// A labelled numeric matrix — features (genes, proteins, compounds) in rows, samples in columns — with
/// per-sample annotations such as condition or batch. This is the common currency between the omics
/// modules, statistics and visualisation (expression matrices, abundance tables, descriptor tables).
/// </summary>
/// <remarks>Values are stored densely in one row-major <see cref="double"/> array, so a row is a contiguous span.</remarks>
public sealed class Dataset
{
    private readonly double[] _values;
    private readonly Dictionary<string, int> _rowIndex;
    private readonly Dictionary<string, int> _columnIndex;

    public Dataset(string name, IReadOnlyList<string> rowNames, IReadOnlyList<string> columnNames, double[] values)
    {
        ArgumentNullException.ThrowIfNull(rowNames);
        ArgumentNullException.ThrowIfNull(columnNames);
        ArgumentNullException.ThrowIfNull(values);
        if (values.Length != rowNames.Count * columnNames.Count)
            throw new ArgumentException($"Expected {rowNames.Count}×{columnNames.Count} values, got {values.Length}.", nameof(values));
        Name = name;
        RowNames = rowNames;
        ColumnNames = columnNames;
        _values = values;
        _rowIndex = BuildIndex(rowNames);
        _columnIndex = BuildIndex(columnNames);
    }

    /// <summary>Creates an all-zero dataset.</summary>
    public static Dataset Zeros(string name, IReadOnlyList<string> rowNames, IReadOnlyList<string> columnNames)
        => new(name, rowNames, columnNames, new double[rowNames.Count * columnNames.Count]);

    /// <summary>Creates a dataset from a rectangular array.</summary>
    public static Dataset FromArray(string name, double[,] values, IReadOnlyList<string>? rowNames = null, IReadOnlyList<string>? columnNames = null)
    {
        int r = values.GetLength(0), c = values.GetLength(1);
        var flat = new double[r * c];
        for (int i = 0; i < r; i++) for (int j = 0; j < c; j++) flat[i * c + j] = values[i, j];
        return new(name, rowNames ?? Enumerable.Range(1, r).Select(i => $"F{i}").ToArray(),
            columnNames ?? Enumerable.Range(1, c).Select(i => $"S{i}").ToArray(), flat);
    }

    public string Name { get; }
    public IReadOnlyList<string> RowNames { get; }
    public IReadOnlyList<string> ColumnNames { get; }
    public int Rows => RowNames.Count;
    public int Columns => ColumnNames.Count;

    /// <summary>Per-sample annotations: key (e.g. <c>condition</c>) → one value per column.</summary>
    public Dictionary<string, string[]> ColumnAnnotations { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>Per-feature annotations: key (e.g. <c>symbol</c>) → one value per row.</summary>
    public Dictionary<string, string[]> RowAnnotations { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>The raw row-major values (mutable, shared).</summary>
    public Span<double> Values => _values;

    public double this[int row, int column]
    {
        get => _values[row * Columns + column];
        set => _values[row * Columns + column] = value;
    }

    public double this[string row, string column]
    {
        get => this[RowIndex(row), ColumnIndex(column)];
        set => this[RowIndex(row), ColumnIndex(column)] = value;
    }

    /// <summary>One row as a contiguous span.</summary>
    public Span<double> Row(int row) => _values.AsSpan(row * Columns, Columns);

    /// <summary>One row by name.</summary>
    public Span<double> Row(string name) => Row(RowIndex(name));

    /// <summary>A copy of one column.</summary>
    public double[] Column(int column)
    {
        var c = new double[Rows];
        for (int i = 0; i < Rows; i++) c[i] = _values[i * Columns + column];
        return c;
    }

    public int RowIndex(string name) => _rowIndex.TryGetValue(name, out int i) ? i : throw new KeyNotFoundException($"No row '{name}'.");
    public int ColumnIndex(string name) => _columnIndex.TryGetValue(name, out int i) ? i : throw new KeyNotFoundException($"No column '{name}'.");

    /// <summary>Column indices whose annotation <paramref name="key"/> equals <paramref name="value"/>.</summary>
    public int[] ColumnsWhere(string key, string value)
    {
        if (!ColumnAnnotations.TryGetValue(key, out var ann)) throw new KeyNotFoundException($"No column annotation '{key}'.");
        return Enumerable.Range(0, Columns).Where(i => ann[i] == value).ToArray();
    }

    /// <summary>Annotates samples, e.g. <c>Annotate("condition", "ctrl", "ctrl", "treated", "treated")</c>.</summary>
    public Dataset Annotate(string key, params string[] values)
    {
        if (values.Length != Columns) throw new ArgumentException($"Need {Columns} values for '{key}'.");
        ColumnAnnotations[key] = values;
        return this;
    }

    /// <summary>A new dataset with the selected rows (annotations carried over).</summary>
    public Dataset SelectRows(IReadOnlyList<int> rows, string? name = null)
    {
        var v = new double[rows.Count * Columns];
        for (int i = 0; i < rows.Count; i++) Row(rows[i]).CopyTo(v.AsSpan(i * Columns, Columns));
        var d = new Dataset(name ?? Name, rows.Select(r => RowNames[r]).ToArray(), ColumnNames, v);
        foreach (var kv in ColumnAnnotations) d.ColumnAnnotations[kv.Key] = (string[])kv.Value.Clone();
        foreach (var kv in RowAnnotations) d.RowAnnotations[kv.Key] = rows.Select(r => kv.Value[r]).ToArray();
        return d;
    }

    /// <summary>A new dataset with the selected columns (annotations carried over).</summary>
    public Dataset SelectColumns(IReadOnlyList<int> columns, string? name = null)
    {
        var v = new double[Rows * columns.Count];
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < columns.Count; j++)
                v[i * columns.Count + j] = _values[i * Columns + columns[j]];
        var d = new Dataset(name ?? Name, RowNames, columns.Select(c => ColumnNames[c]).ToArray(), v);
        foreach (var kv in ColumnAnnotations) d.ColumnAnnotations[kv.Key] = columns.Select(c => kv.Value[c]).ToArray();
        foreach (var kv in RowAnnotations) d.RowAnnotations[kv.Key] = (string[])kv.Value.Clone();
        return d;
    }

    /// <summary>Samples × features copy.</summary>
    public Dataset Transpose()
    {
        var v = new double[_values.Length];
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Columns; j++)
                v[j * Rows + i] = _values[i * Columns + j];
        var d = new Dataset(Name, ColumnNames, RowNames, v);
        foreach (var kv in ColumnAnnotations) d.RowAnnotations[kv.Key] = kv.Value;
        foreach (var kv in RowAnnotations) d.ColumnAnnotations[kv.Key] = kv.Value;
        return d;
    }

    /// <summary>A copy with every value mapped through <paramref name="f"/>.</summary>
    public Dataset Map(Func<double, double> f, string? name = null)
    {
        var v = new double[_values.Length];
        for (int i = 0; i < v.Length; i++) v[i] = f(_values[i]);
        var d = new Dataset(name ?? Name, RowNames, ColumnNames, v);
        foreach (var kv in ColumnAnnotations) d.ColumnAnnotations[kv.Key] = kv.Value;
        foreach (var kv in RowAnnotations) d.RowAnnotations[kv.Key] = kv.Value;
        return d;
    }

    /// <summary>log2(x + pseudocount), the usual transform before plotting or clustering counts.</summary>
    public Dataset Log2(double pseudocount = 1) => Map(x => Math.Log2(x + pseudocount));

    /// <summary>Mean of each row.</summary>
    public double[] RowMeans()
    {
        var m = new double[Rows];
        for (int i = 0; i < Rows; i++)
        {
            double s = 0;
            foreach (double x in Row(i)) s += x;
            m[i] = s / Columns;
        }
        return m;
    }

    /// <summary>Sample variance of each row.</summary>
    public double[] RowVariances()
    {
        var means = RowMeans();
        var v = new double[Rows];
        for (int i = 0; i < Rows; i++)
        {
            double ss = 0;
            foreach (double x in Row(i)) { double d = x - means[i]; ss += d * d; }
            v[i] = Columns > 1 ? ss / (Columns - 1) : 0;
        }
        return v;
    }

    /// <summary>Sum of each column (library sizes for count data).</summary>
    public double[] ColumnSums()
    {
        var s = new double[Columns];
        for (int i = 0; i < Rows; i++)
        {
            var row = Row(i);
            for (int j = 0; j < Columns; j++) s[j] += row[j];
        }
        return s;
    }

    /// <summary>Writes tab-separated values with a header row; column annotations are not written.</summary>
    public void WriteTsv(TextWriter writer)
    {
        writer.Write("feature");
        foreach (var c in ColumnNames) { writer.Write('\t'); writer.Write(c); }
        writer.WriteLine();
        for (int i = 0; i < Rows; i++)
        {
            writer.Write(RowNames[i]);
            foreach (double x in Row(i)) { writer.Write('\t'); writer.Write(x.ToString("R", CultureInfo.InvariantCulture)); }
            writer.WriteLine();
        }
    }

    /// <summary>Reads a matrix written as TSV/CSV: header row of sample names, first column feature names.</summary>
    public static Dataset ReadDelimited(TextReader reader, string name = "dataset", char separator = '\t')
    {
        string header = reader.ReadLine() ?? throw new InvalidDataException("Empty matrix file.");
        var columns = header.Split(separator)[1..];
        var rows = new List<string>();
        var values = new List<double>();
        string? line;
        int lineNo = 1;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (line.Length == 0) continue;
            var parts = line.Split(separator);
            if (parts.Length != columns.Length + 1)
                throw new InvalidDataException($"Line {lineNo}: expected {columns.Length + 1} fields, found {parts.Length}.");
            rows.Add(parts[0]);
            for (int j = 1; j < parts.Length; j++)
                values.Add(double.TryParse(parts[j], NumberStyles.Float, CultureInfo.InvariantCulture, out double x) ? x : double.NaN);
        }
        return new(name, rows, columns, values.ToArray());
    }

    public override string ToString()
    {
        var sb = new StringBuilder($"{Name}: {Rows} features × {Columns} samples");
        if (ColumnAnnotations.Count > 0) sb.Append($" [{string.Join(", ", ColumnAnnotations.Keys)}]");
        return sb.ToString();
    }

    private static Dictionary<string, int> BuildIndex(IReadOnlyList<string> names)
    {
        var d = new Dictionary<string, int>(names.Count, StringComparer.Ordinal);
        for (int i = 0; i < names.Count; i++) d.TryAdd(names[i], i);
        return d;
    }
}
