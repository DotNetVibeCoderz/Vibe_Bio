namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// A compressed sparse row (CSR) matrix of <see cref="float"/> values — the storage for single-cell count matrices
/// (cells in rows, genes in columns, as in AnnData/scanpy). Column indices are sorted within each row.
/// </summary>
/// <remarks>
/// Memory is 8 bytes per non-zero (4 for the value, 4 for the column index) plus 4 per row; a 10x run with
/// 10 000 cells × 33 000 genes at ~2 000 detected genes per cell takes about 160 MB versus 1.3 GB dense.
/// </remarks>
public sealed class SparseMatrix : ISparseMatrixSource, IRowPointers
{
    /// <summary>Wraps CSR arrays without copying. Column indices must be sorted per row.</summary>
    public SparseMatrix(int rows, int columns, int[] rowPointers, int[] columnIndices, float[] values)
    {
        ArgumentNullException.ThrowIfNull(rowPointers);
        ArgumentNullException.ThrowIfNull(columnIndices);
        ArgumentNullException.ThrowIfNull(values);
        if (rowPointers.Length != rows + 1) throw new ArgumentException("rowPointers needs rows + 1 entries.");
        if (columnIndices.Length != values.Length || rowPointers[rows] != values.Length)
            throw new ArgumentException("columnIndices and values must both have rowPointers[rows] entries.");
        Rows = rows;
        Columns = columns;
        RowPointers = rowPointers;
        ColumnIndices = columnIndices;
        Values = values;
    }

    /// <summary>Number of rows (cells).</summary>
    public int Rows { get; }

    /// <summary>Number of columns (genes).</summary>
    public int Columns { get; }

    /// <summary>Row start offsets (length Rows + 1).</summary>
    public int[] RowPointers { get; }

    /// <summary>Column index of each stored value.</summary>
    public int[] ColumnIndices { get; }

    /// <summary>Stored values.</summary>
    public float[] Values { get; }

    /// <summary>Number of stored (non-zero) entries.</summary>
    public int NonZeroCount => RowPointers[Rows];

    /// <summary>Fraction of entries that are stored.</summary>
    public double Density => Rows == 0 || Columns == 0 ? 0 : (double)NonZeroCount / ((long)Rows * Columns);

    /// <summary>Column indices of the stored entries of one row.</summary>
    public ReadOnlySpan<int> RowIndices(int row) => ColumnIndices.AsSpan(RowPointers[row], RowPointers[row + 1] - RowPointers[row]);

    /// <summary>Stored values of one row (mutable).</summary>
    public Span<float> RowValues(int row) => Values.AsSpan(RowPointers[row], RowPointers[row + 1] - RowPointers[row]);

    /// <summary>Element lookup by binary search (0 when not stored).</summary>
    public float this[int row, int column]
    {
        get
        {
            int i = RowIndices(row).BinarySearch(column);
            return i >= 0 ? Values[RowPointers[row] + i] : 0f;
        }
    }

    /// <summary>Builds a CSR matrix from (row, column, value) triplets; duplicates are summed and zeros dropped.</summary>
    public static SparseMatrix FromTriplets(int rows, int columns, ReadOnlySpan<int> rowIdx, ReadOnlySpan<int> colIdx, ReadOnlySpan<float> values)
    {
        if (rowIdx.Length != colIdx.Length || rowIdx.Length != values.Length) throw new ArgumentException("Triplet arrays differ in length.");
        var counts = new int[rows + 1];
        foreach (int r in rowIdx)
        {
            if ((uint)r >= (uint)rows) throw new ArgumentOutOfRangeException(nameof(rowIdx));
            counts[r + 1]++;
        }
        for (int r = 0; r < rows; r++) counts[r + 1] += counts[r];
        var cols = new int[rowIdx.Length];
        var vals = new float[rowIdx.Length];
        var next = (int[])counts.Clone();
        for (int k = 0; k < rowIdx.Length; k++)
        {
            if ((uint)colIdx[k] >= (uint)columns) throw new ArgumentOutOfRangeException(nameof(colIdx));
            int at = next[rowIdx[k]]++;
            cols[at] = colIdx[k]; vals[at] = values[k];
        }
        // Sort each row by column, merge duplicates, drop zeros.
        var ptr = new int[rows + 1];
        int w = 0;
        for (int r = 0; r < rows; r++)
        {
            int s = counts[r], e = counts[r + 1];
            Array.Sort(cols, vals, s, e - s);
            int rowStart = w;
            for (int k = s; k < e; k++)
            {
                if (w > rowStart && cols[w - 1] == cols[k]) vals[w - 1] += vals[k];
                else { cols[w] = cols[k]; vals[w] = vals[k]; w++; }
            }
            int keep = rowStart;
            for (int k = rowStart; k < w; k++) if (vals[k] != 0) { cols[keep] = cols[k]; vals[keep] = vals[k]; keep++; }
            w = keep;
            ptr[r + 1] = w;
        }
        return new SparseMatrix(rows, columns, ptr, cols[..w], vals[..w]);
    }

    /// <summary>A sparse copy of a dense row-major matrix.</summary>
    public static SparseMatrix FromDense(DenseMatrix dense)
    {
        ArgumentNullException.ThrowIfNull(dense);
        var ptr = new int[dense.Rows + 1];
        var cols = new List<int>(); var vals = new List<float>();
        for (int i = 0; i < dense.Rows; i++)
        {
            var row = dense.Row(i);
            for (int j = 0; j < row.Length; j++) if (row[j] != 0) { cols.Add(j); vals.Add((float)row[j]); }
            ptr[i + 1] = cols.Count;
        }
        return new SparseMatrix(dense.Rows, dense.Columns, ptr, cols.ToArray(), vals.ToArray());
    }

    /// <summary>A dense copy.</summary>
    public DenseMatrix ToDense()
    {
        var d = new DenseMatrix(Rows, Columns);
        for (int i = 0; i < Rows; i++)
        {
            var idx = RowIndices(i); var v = RowValues(i);
            for (int k = 0; k < idx.Length; k++) d.Data[i * Columns + idx[k]] = v[k];
        }
        return d;
    }

    /// <summary>A dense copy of the selected columns only (cells × selected genes).</summary>
    public DenseMatrix ToDense(ReadOnlySpan<int> columns)
    {
        var map = new int[Columns];
        Array.Fill(map, -1);
        for (int k = 0; k < columns.Length; k++) map[columns[k]] = k;
        var d = new DenseMatrix(Rows, columns.Length);
        for (int i = 0; i < Rows; i++)
        {
            var idx = RowIndices(i); var v = RowValues(i);
            for (int k = 0; k < idx.Length; k++)
            {
                int c = map[idx[k]];
                if (c >= 0) d.Data[i * columns.Length + c] = v[k];
            }
        }
        return d;
    }

    /// <summary>The transpose, also in CSR (i.e. this matrix in CSC form).</summary>
    public SparseMatrix Transpose()
    {
        var ptr = new int[Columns + 1];
        foreach (int c in ColumnIndices) ptr[c + 1]++;
        for (int c = 0; c < Columns; c++) ptr[c + 1] += ptr[c];
        var next = (int[])ptr.Clone();
        var rows = new int[NonZeroCount];
        var vals = new float[NonZeroCount];
        for (int i = 0; i < Rows; i++)
            for (int k = RowPointers[i]; k < RowPointers[i + 1]; k++)
            {
                int at = next[ColumnIndices[k]]++;
                rows[at] = i; vals[at] = Values[k];
            }
        return new SparseMatrix(Columns, Rows, ptr, rows, vals);
    }

    /// <summary>A copy with only the selected rows, in the given order.</summary>
    public SparseMatrix SelectRows(ReadOnlySpan<int> rows)
    {
        var ptr = new int[rows.Length + 1];
        for (int k = 0; k < rows.Length; k++) ptr[k + 1] = ptr[k] + RowPointers[rows[k] + 1] - RowPointers[rows[k]];
        var cols = new int[ptr[rows.Length]];
        var vals = new float[ptr[rows.Length]];
        for (int k = 0; k < rows.Length; k++)
        {
            RowIndices(rows[k]).CopyTo(cols.AsSpan(ptr[k]));
            RowValues(rows[k]).CopyTo(vals.AsSpan(ptr[k]));
        }
        return new SparseMatrix(rows.Length, Columns, ptr, cols, vals);
    }

    /// <summary>A copy with only the selected columns (renumbered in the given order, which must be ascending).</summary>
    public SparseMatrix SelectColumns(ReadOnlySpan<int> columns)
    {
        var map = new int[Columns];
        Array.Fill(map, -1);
        for (int k = 0; k < columns.Length; k++) map[columns[k]] = k;
        var ptr = new int[Rows + 1];
        var cols = new List<int>(NonZeroCount); var vals = new List<float>(NonZeroCount);
        for (int i = 0; i < Rows; i++)
        {
            for (int k = RowPointers[i]; k < RowPointers[i + 1]; k++)
            {
                int c = map[ColumnIndices[k]];
                if (c >= 0) { cols.Add(c); vals.Add(Values[k]); }
            }
            ptr[i + 1] = cols.Count;
        }
        return new SparseMatrix(Rows, columns.Length, ptr, cols.ToArray(), vals.ToArray());
    }

    /// <summary>Sum of each row.</summary>
    public double[] RowSums()
    {
        var s = new double[Rows];
        for (int i = 0; i < Rows; i++) { double t = 0; foreach (float v in RowValues(i)) t += v; s[i] = t; }
        return s;
    }

    /// <summary>Sum of each column.</summary>
    public double[] ColumnSums()
    {
        var s = new double[Columns];
        for (int k = 0; k < NonZeroCount; k++) s[ColumnIndices[k]] += Values[k];
        return s;
    }

    /// <summary>Number of non-zero entries in each row.</summary>
    public int[] RowNonZeroCounts()
    {
        var c = new int[Rows];
        for (int i = 0; i < Rows; i++) c[i] = RowPointers[i + 1] - RowPointers[i];
        return c;
    }

    /// <summary>Number of non-zero entries in each column.</summary>
    public int[] ColumnNonZeroCounts()
    {
        var c = new int[Columns];
        for (int k = 0; k < NonZeroCount; k++) if (Values[k] != 0) c[ColumnIndices[k]]++;
        return c;
    }

    /// <summary>Per-column mean and sample variance (ddof = 1), counting implicit zeros.</summary>
    public (double[] Mean, double[] Variance) ColumnMeanVariance()
    {
        var sum = new double[Columns]; var sq = new double[Columns];
        for (int k = 0; k < NonZeroCount; k++) { double v = Values[k]; sum[ColumnIndices[k]] += v; sq[ColumnIndices[k]] += v * v; }
        var mean = new double[Columns]; var variance = new double[Columns];
        for (int j = 0; j < Columns; j++)
        {
            mean[j] = sum[j] / Rows;
            variance[j] = Rows > 1 ? (sq[j] - Rows * mean[j] * mean[j]) / (Rows - 1) : 0;
            if (variance[j] < 0) variance[j] = 0;
        }
        return (mean, variance);
    }

    long ISparseMatrixSource.NonZeroCount => NonZeroCount;

    long IRowPointers.RowPointer(int row) => RowPointers[row];

    /// <summary>A copy of rows <paramref name="start"/> … <paramref name="start"/> + <paramref name="count"/> − 1.</summary>
    public SparseMatrix ReadRows(int start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        if (start + count > Rows) throw new ArgumentOutOfRangeException(nameof(count));
        int s = RowPointers[start], e = RowPointers[start + count];
        var ptr = new int[count + 1];
        for (int i = 0; i <= count; i++) ptr[i] = RowPointers[start + i] - s;
        return new SparseMatrix(count, Columns, ptr, ColumnIndices[s..e], Values[s..e]);
    }

    /// <summary>A deep copy.</summary>
    public SparseMatrix Clone()
        => new(Rows, Columns, (int[])RowPointers.Clone(), (int[])ColumnIndices.Clone(), (float[])Values.Clone());

    /// <inheritdoc />
    public override string ToString() => $"SparseMatrix {Rows}×{Columns}, {NonZeroCount} stored ({Density:P2})";
}
