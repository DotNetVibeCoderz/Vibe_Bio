using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics;

/// <summary>
/// A plain row-major matrix of doubles without labels — embeddings, PCA scores and loadings, scaled
/// expression and distance matrices. Use <see cref="Dataset"/> when rows and columns need names.
/// </summary>
public sealed class DenseMatrix
{
    /// <summary>Creates a matrix, wrapping <paramref name="data"/> without copying when given.</summary>
    public DenseMatrix(int rows, int columns, double[]? data = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(rows);
        ArgumentOutOfRangeException.ThrowIfNegative(columns);
        if (data is not null && data.Length != rows * columns)
            throw new ArgumentException($"Expected {rows}×{columns} values, got {data.Length}.", nameof(data));
        Rows = rows;
        Columns = columns;
        Data = data ?? new double[rows * columns];
    }

    /// <summary>Number of rows.</summary>
    public int Rows { get; }

    /// <summary>Number of columns.</summary>
    public int Columns { get; }

    /// <summary>The row-major backing array (shared, mutable).</summary>
    public double[] Data { get; }

    /// <summary>Element access.</summary>
    public double this[int row, int column]
    {
        get => Data[row * Columns + column];
        set => Data[row * Columns + column] = value;
    }

    /// <summary>One row as a contiguous span.</summary>
    public Span<double> Row(int row) => Data.AsSpan(row * Columns, Columns);

    /// <summary>A copy of one column.</summary>
    public double[] Column(int column)
    {
        var c = new double[Rows];
        for (int i = 0; i < Rows; i++) c[i] = Data[i * Columns + column];
        return c;
    }

    /// <summary>A transposed copy.</summary>
    public DenseMatrix Transpose()
    {
        var t = new DenseMatrix(Columns, Rows);
        for (int i = 0; i < Rows; i++)
            for (int j = 0; j < Columns; j++)
                t.Data[j * Rows + i] = Data[i * Columns + j];
        return t;
    }

    /// <summary>The values of a <see cref="Dataset"/> (copied), optionally transposed to samples × features.</summary>
    public static DenseMatrix FromDataset(Dataset dataset, bool transpose = false)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        var m = new DenseMatrix(dataset.Rows, dataset.Columns, dataset.Values.ToArray());
        return transpose ? m.Transpose() : m;
    }

    /// <summary>Wraps the matrix as a labelled <see cref="Dataset"/> (data copied).</summary>
    public Dataset ToDataset(string name, IReadOnlyList<string>? rowNames = null, IReadOnlyList<string>? columnNames = null)
    {
        var rn = rowNames ?? Enumerable.Range(1, Rows).Select(i => $"R{i}").ToArray();
        var cn = columnNames ?? Enumerable.Range(1, Columns).Select(i => $"C{i}").ToArray();
        return new Dataset(name, rn, cn, (double[])Data.Clone());
    }

    /// <inheritdoc />
    public override string ToString() => $"DenseMatrix {Rows}×{Columns}";
}
