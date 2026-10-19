namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>On-disk layout of a matrix element in an AnnData file.</summary>
public enum AnnMatrixFormat
{
    /// <summary>A dense 2-D array (<c>encoding-type = "array"</c>).</summary>
    Dense,

    /// <summary>Compressed sparse rows (<c>csr_matrix</c>).</summary>
    Csr,

    /// <summary>Compressed sparse columns (<c>csc_matrix</c>).</summary>
    Csc,
}

/// <summary>
/// A matrix element of an AnnData object (X, a layer, an <c>obsm</c>/<c>varm</c> entry, an <c>obsp</c>/<c>varp</c>
/// graph). Sparse data is held as a CSR <see cref="SparseMatrix"/> (float32 values) in the logical rows × columns
/// orientation whatever the on-disk format; dense data as a <see cref="DenseMatrix"/> (float64).
/// <see cref="Format"/> and <see cref="DataType"/> record the on-disk layout and element type so that
/// <see cref="H5ad.Write(string, AnnData, H5adWriteOptions?)"/> writes the element back the same way.
/// </summary>
public sealed class AnnMatrix
{
    private AnnMatrix(AnnMatrixFormat format, string dataType, SparseMatrix? sparse, DenseMatrix? dense, double[]? values64 = null)
    {
        Format = format; DataType = dataType; Sparse = sparse; Dense = dense; SparseValues64 = values64;
    }

    /// <summary>
    /// A sparse element (stored as CSR by default, or CSC). <paramref name="values64"/> optionally carries the
    /// full-precision values (aligned with <see cref="SparseMatrix.Values"/>) for float64/int64 data.
    /// </summary>
    public static AnnMatrix FromSparse(SparseMatrix matrix, AnnMatrixFormat format = AnnMatrixFormat.Csr, string dataType = "float32", double[]? values64 = null)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        if (format == AnnMatrixFormat.Dense) throw new ArgumentException("Use FromDense for dense elements.", nameof(format));
        if (values64 is not null && values64.Length != matrix.NonZeroCount) throw new ArgumentException("values64 must have one value per stored entry.", nameof(values64));
        return new AnnMatrix(format, H5Types.CheckNumericType(dataType), matrix, null, values64);
    }

    /// <summary>
    /// Full-precision values of a sparse float64/int64 element whose values do not all fit float32 exactly,
    /// aligned with <see cref="SparseMatrix.Values"/> (whose float32 copies are what the single-cell algorithms
    /// use); null otherwise.
    /// </summary>
    public double[]? SparseValues64 { get; }

    /// <summary>A dense element.</summary>
    public static AnnMatrix FromDense(DenseMatrix matrix, string dataType = "float64")
    {
        ArgumentNullException.ThrowIfNull(matrix);
        return new AnnMatrix(AnnMatrixFormat.Dense, H5Types.CheckNumericType(dataType), null, matrix);
    }

    /// <summary>On-disk layout.</summary>
    public AnnMatrixFormat Format { get; }

    /// <summary>numpy element type name used on disk: float32, float64, int8 … int64, uint8 … uint64 or bool.</summary>
    public string DataType { get; }

    /// <summary>True for CSR/CSC elements.</summary>
    public bool IsSparse => Sparse is not null;

    /// <summary>The sparse data (null for dense elements).</summary>
    public SparseMatrix? Sparse { get; }

    /// <summary>The dense data (null for sparse elements).</summary>
    public DenseMatrix? Dense { get; }

    /// <summary>Number of rows.</summary>
    public int Rows => Sparse?.Rows ?? Dense!.Rows;

    /// <summary>Number of columns.</summary>
    public int Columns => Sparse?.Columns ?? Dense!.Columns;

    /// <summary>The element as a sparse matrix (zeros dropped when converting from dense).</summary>
    public SparseMatrix ToSparse() => Sparse ?? SparseMatrix.FromDense(Dense!);

    /// <summary>The element as a dense matrix.</summary>
    public DenseMatrix ToDense()
    {
        if (Dense is not null) return Dense;
        var d = Sparse!.ToDense();
        if (SparseValues64 is not null)
            for (int i = 0; i < Sparse.Rows; i++)
            {
                var idx = Sparse.RowIndices(i);
                for (int k = 0; k < idx.Length; k++) d[i, idx[k]] = SparseValues64[Sparse.RowPointers[i] + k];
            }
        return d;
    }

    /// <summary>Element lookup.</summary>
    public double this[int row, int column]
    {
        get
        {
            if (Dense is not null) return Dense[row, column];
            int i = Sparse!.RowIndices(row).BinarySearch(column);
            if (i < 0) return 0;
            int at = Sparse.RowPointers[row] + i;
            return SparseValues64 is not null ? SparseValues64[at] : Sparse.Values[at];
        }
    }

    /// <inheritdoc />
    public override string ToString() => $"AnnMatrix {Rows}×{Columns} {Format} {DataType}";
}

/// <summary>
/// A pandas categorical column: integer <see cref="Codes"/> into <see cref="Categories"/>, −1 meaning missing.
/// Categories are usually strings (<c>string[]</c>) but may be numeric arrays.
/// </summary>
public sealed class Categorical
{
    /// <summary>Wraps codes and categories without copying.</summary>
    public Categorical(int[] codes, Array categories, bool ordered = false)
    {
        ArgumentNullException.ThrowIfNull(codes);
        ArgumentNullException.ThrowIfNull(categories);
        foreach (int c in codes)
            if (c < -1 || c >= categories.Length) throw new ArgumentException($"Code {c} is outside the {categories.Length} categories.");
        Codes = codes; Categories = categories; Ordered = ordered;
    }

    /// <summary>Builds a categorical from labels (null = missing); categories are the distinct labels in ordinal order, as pandas sorts them.</summary>
    public static Categorical FromLabels(IReadOnlyList<string?> labels, bool ordered = false, IReadOnlyList<string>? categories = null)
    {
        ArgumentNullException.ThrowIfNull(labels);
        var cats = categories?.ToArray() ?? labels.Where(l => l is not null).Distinct(StringComparer.Ordinal).Order(StringComparer.Ordinal).ToArray()!;
        var map = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < cats.Length; i++) map[cats[i]!] = i;
        var codes = new int[labels.Count];
        for (int i = 0; i < codes.Length; i++)
            codes[i] = labels[i] is null ? -1 : map.TryGetValue(labels[i]!, out int c) ? c : throw new ArgumentException($"Label '{labels[i]}' is not a category.");
        return new Categorical(codes, cats, ordered);
    }

    /// <summary>Category index of each row (−1 = missing).</summary>
    public int[] Codes { get; }

    /// <summary>The categories (a <c>string[]</c> or a numeric array).</summary>
    public Array Categories { get; }

    /// <summary>Whether the categories have an order.</summary>
    public bool Ordered { get; }

    /// <summary>Number of rows.</summary>
    public int Length => Codes.Length;

    /// <summary>The value of row <paramref name="row"/> (null when missing).</summary>
    public object? this[int row] => Codes[row] < 0 ? null : Categories.GetValue(Codes[row]);

    /// <summary>Every row as a string (null when missing; numeric categories formatted invariantly).</summary>
    public string?[] ToStrings()
    {
        var names = new string[Categories.Length];
        for (int i = 0; i < names.Length; i++) names[i] = Convert.ToString(Categories.GetValue(i), System.Globalization.CultureInfo.InvariantCulture) ?? "";
        return Codes.Select(c => c < 0 ? null : names[c]).ToArray();
    }

    /// <summary>A copy holding only the given rows.</summary>
    public Categorical Select(ReadOnlySpan<int> rows)
    {
        var codes = new int[rows.Length];
        for (int i = 0; i < rows.Length; i++) codes[i] = Codes[rows[i]];
        return new Categorical(codes, Categories, Ordered);
    }

    /// <inheritdoc />
    public override string ToString() => $"Categorical[{Length}] with {Categories.Length} categories{(Ordered ? " (ordered)" : "")}";
}

/// <summary>
/// A nullable column (pandas <c>Int64</c>/<c>boolean</c>/<c>string</c> extension arrays): <see cref="Values"/> plus
/// a <see cref="Mask"/> that is true where the value is missing.
/// </summary>
public sealed class MaskedArray
{
    /// <summary>Wraps values and mask without copying.</summary>
    public MaskedArray(Array values, bool[] mask)
    {
        ArgumentNullException.ThrowIfNull(values);
        ArgumentNullException.ThrowIfNull(mask);
        if (values.Length != mask.Length) throw new ArgumentException("Values and mask differ in length.");
        if (values is not (long[] or int[] or bool[] or string[] or double[] or float[] or short[] or sbyte[] or byte[] or ushort[] or uint[] or ulong[]))
            throw new ArgumentException($"Unsupported value array {values.GetType().Name}.", nameof(values));
        Values = values; Mask = mask;
    }

    /// <summary>The values (entries under the mask are unspecified).</summary>
    public Array Values { get; }

    /// <summary>True where the value is missing.</summary>
    public bool[] Mask { get; }

    /// <summary>
    /// For string columns: the pandas missing-value marker recorded in the file (<c>"NaN"</c> for the pandas 3
    /// <c>str</c> dtype, <c>"NA"</c> for <c>string[python]</c>); written back as the <c>na-value</c> attribute.
    /// </summary>
    public string? NaValue { get; init; }

    /// <summary>Number of rows.</summary>
    public int Length => Mask.Length;

    /// <summary>The value of row <paramref name="row"/> (null when missing).</summary>
    public object? this[int row] => Mask[row] ? null : Values.GetValue(row);

    /// <summary>A copy holding only the given rows.</summary>
    public MaskedArray Select(ReadOnlySpan<int> rows)
    {
        var v = Array.CreateInstance(Values.GetType().GetElementType()!, rows.Length);
        var m = new bool[rows.Length];
        for (int i = 0; i < rows.Length; i++) { v.SetValue(Values.GetValue(rows[i]), i); m[i] = Mask[rows[i]]; }
        return new MaskedArray(v, m) { NaValue = NaValue };
    }

    /// <inheritdoc />
    public override string ToString() => $"MaskedArray<{Values.GetType().GetElementType()!.Name}>[{Length}], {Mask.Count(x => x)} missing";
}

/// <summary>
/// An N-dimensional numeric or string array from <c>uns</c> (row-major <see cref="Data"/> with its <see cref="Shape"/>).
/// One-dimensional arrays are stored in <c>uns</c> as plain .NET arrays instead.
/// </summary>
public sealed class AnnArray
{
    /// <summary>Wraps data and shape.</summary>
    public AnnArray(Array data, long[] shape)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(shape);
        long n = 1;
        foreach (long s in shape) n *= s;
        if (n != data.Length) throw new ArgumentException("Shape does not match the data length.");
        Data = data; Shape = shape;
    }

    /// <summary>Row-major data.</summary>
    public Array Data { get; }

    /// <summary>Dimensions.</summary>
    public long[] Shape { get; }

    /// <inheritdoc />
    public override string ToString() => $"AnnArray<{Data.GetType().GetElementType()!.Name}>[{string.Join("×", Shape)}]";
}

/// <summary>An element the reader recognised but cannot represent (e.g. compound/record arrays, awkward arrays); kept so nothing is silently dropped.</summary>
public sealed record AnnUnsupported(string Path, string EncodingType, string Reason);

/// <summary>
/// A pandas DataFrame as AnnData stores it (<c>obs</c>, <c>var</c>): a string index and ordered columns. Column
/// values are <c>double[]</c>, <c>float[]</c>, <c>int[]</c>, <c>long[]</c> (and other primitive arrays),
/// <c>bool[]</c>, <c>string[]</c>, <see cref="Categorical"/> or <see cref="MaskedArray"/>.
/// </summary>
public sealed class AnnDataFrame
{
    private readonly List<string> _order = [];
    private readonly Dictionary<string, Array> _arrays = new(StringComparer.Ordinal);
    private readonly Dictionary<string, object> _columns = new(StringComparer.Ordinal);

    /// <summary>Creates an empty frame over the given index.</summary>
    public AnnDataFrame(IReadOnlyList<string> index, string indexName = "_index")
    {
        ArgumentNullException.ThrowIfNull(index);
        Index = index.ToArray();
        IndexName = string.IsNullOrEmpty(indexName) ? "_index" : indexName;
    }

    /// <summary>Row labels (cell barcodes or gene ids).</summary>
    public string[] Index { get; private set; }

    /// <summary>Name of the index (pandas <c>index.name</c>; <c>_index</c> when unnamed).</summary>
    public string IndexName { get; set; }

    /// <summary>Number of rows.</summary>
    public int RowCount => Index.Length;

    /// <summary>Column names in order.</summary>
    public IReadOnlyList<string> ColumnNames => _order;

    /// <summary>A column's values (see the class remarks for the possible types).</summary>
    public object this[string column]
    {
        get => _columns[column];
        set => Set(column, value);
    }

    /// <summary>Whether a column exists.</summary>
    public bool Contains(string column) => _columns.ContainsKey(column);

    /// <summary>Adds or replaces a column (appended at the end when new).</summary>
    public void Set(string column, object values)
    {
        ArgumentException.ThrowIfNullOrEmpty(column);
        ArgumentNullException.ThrowIfNull(values);
        int length = values switch
        {
            Categorical c => c.Length,
            MaskedArray m => m.Length,
            Array a when a.Rank == 1 && IsSupportedArray(a) => a.Length,
            _ => throw new ArgumentException($"Unsupported column type {values.GetType().Name}.", nameof(values)),
        };
        if (length != RowCount) throw new ArgumentException($"Column '{column}' has {length} rows, the frame has {RowCount}.");
        if (!_columns.ContainsKey(column)) _order.Add(column);
        _columns[column] = values;
    }

    /// <summary>Removes a column.</summary>
    public bool Remove(string column)
    {
        _order.Remove(column);
        return _columns.Remove(column);
    }

    /// <summary>A numeric column as doubles (missing values of nullable columns become NaN; bool → 0/1).</summary>
    public double[] GetDoubles(string column)
    {
        var v = _columns[column];
        if (v is MaskedArray m)
        {
            var d = ToDoubles(m.Values);
            for (int i = 0; i < d.Length; i++) if (m.Mask[i]) d[i] = double.NaN;
            return d;
        }
        if (v is Array a && a is not string[]) return ToDoubles(a);
        throw new InvalidCastException($"Column '{column}' is not numeric.");
    }

    /// <summary>A column as strings (categoricals expanded, numbers formatted invariantly, missing → null).</summary>
    public string?[] GetStrings(string column)
    {
        var v = _columns[column];
        return v switch
        {
            string[] s => s,
            Categorical c => c.ToStrings(),
            MaskedArray m => Enumerable.Range(0, m.Length).Select(i => m.Mask[i] ? null : Convert.ToString(m.Values.GetValue(i), System.Globalization.CultureInfo.InvariantCulture)).ToArray(),
            Array a => Enumerable.Range(0, a.Length).Select(i => Convert.ToString(a.GetValue(i), System.Globalization.CultureInfo.InvariantCulture)).ToArray(),
            _ => throw new InvalidCastException(),
        };
    }

    /// <summary>A copy holding only the given rows (in order).</summary>
    public AnnDataFrame SelectRows(ReadOnlySpan<int> rows)
    {
        var idx = new string[rows.Length];
        for (int i = 0; i < rows.Length; i++) idx[i] = Index[rows[i]];
        var f = new AnnDataFrame(idx, IndexName);
        foreach (var name in _order)
        {
            object v = _columns[name] switch
            {
                Categorical c => c.Select(rows),
                MaskedArray m => m.Select(rows),
                Array a => SelectArray(a, rows),
                _ => throw new InvalidOperationException(),
            };
            f.Set(name, v);
        }
        return f;
    }

    internal static Array SelectArray(Array a, ReadOnlySpan<int> rows)
    {
        var r = Array.CreateInstance(a.GetType().GetElementType()!, rows.Length);
        for (int i = 0; i < rows.Length; i++) r.SetValue(a.GetValue(rows[i]), i);
        return r;
    }

    internal static bool IsSupportedArray(Array a) => a is double[] or float[] or int[] or long[] or short[] or sbyte[] or byte[] or ushort[] or uint[] or ulong[] or bool[] or string[];

    internal static double[] ToDoubles(Array a) => a switch
    {
        double[] d => (double[])d.Clone(),
        float[] f => f.Select(x => (double)x).ToArray(),
        int[] i => i.Select(x => (double)x).ToArray(),
        long[] l => l.Select(x => (double)x).ToArray(),
        bool[] b => b.Select(x => x ? 1.0 : 0.0).ToArray(),
        _ => Enumerable.Range(0, a.Length).Select(i => Convert.ToDouble(a.GetValue(i), System.Globalization.CultureInfo.InvariantCulture)).ToArray(),
    };

    /// <inheritdoc />
    public override string ToString() => $"AnnDataFrame {RowCount} rows × [{string.Join(", ", _order)}]";
}

/// <summary>The frozen <c>raw</c> slot of an AnnData object: the unfiltered X with its own <c>var</c> and <c>varm</c>.</summary>
public sealed class AnnRaw
{
    /// <summary>Creates a raw slot.</summary>
    public AnnRaw(AnnMatrix x, AnnDataFrame var)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(var);
        if (var.RowCount != x.Columns) throw new ArgumentException("raw.var must have one row per raw.X column.");
        X = x; Var = var;
    }

    /// <summary>The raw matrix (cells × raw genes).</summary>
    public AnnMatrix X { get; }

    /// <summary>Raw gene annotations.</summary>
    public AnnDataFrame Var { get; }

    /// <summary>Raw per-gene matrices.</summary>
    public Dictionary<string, AnnMatrix> Varm { get; } = new(StringComparer.Ordinal);
}

/// <summary>
/// An in-memory AnnData object, the layout of Python's <c>anndata</c> and of <c>.h5ad</c> files: X (cells × genes),
/// <see cref="Obs"/>/<see cref="Var"/> annotations, <see cref="Obsm"/>/<see cref="Varm"/> multi-dimensional
/// annotations, <see cref="Obsp"/>/<see cref="Varp"/> pairwise graphs, <see cref="Layers"/>, unstructured
/// <see cref="Uns"/> and the optional <see cref="Raw"/> slot. Read and write it with <see cref="H5ad"/>; convert to the
/// analysis type with <see cref="ToSingleCellExperiment"/>.
/// </summary>
public sealed class AnnData
{
    /// <summary>Creates an AnnData object. <paramref name="x"/> may be null (an AnnData without X).</summary>
    public AnnData(AnnMatrix? x, AnnDataFrame obs, AnnDataFrame var)
    {
        ArgumentNullException.ThrowIfNull(obs);
        ArgumentNullException.ThrowIfNull(var);
        if (x is not null && (x.Rows != obs.RowCount || x.Columns != var.RowCount))
            throw new ArgumentException($"X is {x.Rows}×{x.Columns} but obs has {obs.RowCount} rows and var {var.RowCount}.");
        X = x; Obs = obs; Var = var;
    }

    /// <summary>The main matrix (null when the file has no X).</summary>
    public AnnMatrix? X { get; set; }

    /// <summary>Per-cell annotations.</summary>
    public AnnDataFrame Obs { get; }

    /// <summary>Per-gene annotations.</summary>
    public AnnDataFrame Var { get; }

    /// <summary>Number of cells.</summary>
    public int ObsCount => Obs.RowCount;

    /// <summary>Number of genes.</summary>
    public int VarCount => Var.RowCount;

    /// <summary>Cell-level matrices (<c>X_pca</c>, <c>X_umap</c>), each cells × k. Values are <see cref="AnnMatrix"/> or <see cref="AnnDataFrame"/>.</summary>
    public Dictionary<string, object> Obsm { get; } = new(StringComparer.Ordinal);

    /// <summary>Gene-level matrices (<c>PCs</c>). Values are <see cref="AnnMatrix"/> or <see cref="AnnDataFrame"/>.</summary>
    public Dictionary<string, object> Varm { get; } = new(StringComparer.Ordinal);

    /// <summary>Cell × cell graphs (<c>distances</c>, <c>connectivities</c>).</summary>
    public Dictionary<string, AnnMatrix> Obsp { get; } = new(StringComparer.Ordinal);

    /// <summary>Gene × gene matrices.</summary>
    public Dictionary<string, AnnMatrix> Varp { get; } = new(StringComparer.Ordinal);

    /// <summary>Matrices with the shape of X (<c>counts</c>, <c>spliced</c>).</summary>
    public Dictionary<string, AnnMatrix> Layers { get; } = new(StringComparer.Ordinal);

    /// <summary>
    /// Unstructured annotations. Values: nested <c>Dictionary&lt;string, object?&gt;</c>, <c>string</c>,
    /// <c>long</c>/<c>double</c>/<c>bool</c> scalars (other numeric scalars keep their .NET type), 1-D arrays
    /// (<c>double[]</c>, <c>long[]</c>, <c>string[]</c>, <c>bool[]</c> …), <see cref="AnnArray"/> for N-d arrays,
    /// <see cref="AnnMatrix"/>, <see cref="AnnDataFrame"/>, <see cref="Categorical"/>, <see cref="MaskedArray"/> or
    /// <see cref="AnnUnsupported"/>.
    /// </summary>
    public Dictionary<string, object?> Uns { get; } = new(StringComparer.Ordinal);

    /// <summary>The frozen raw data, if any.</summary>
    public AnnRaw? Raw { get; set; }

    /// <summary>
    /// Converts to a <see cref="SingleCellExperiment"/> (X and sparse layers, numeric/bool/categorical annotations,
    /// dense <c>obsm</c> embeddings). Gene symbols come from <paramref name="geneSymbolColumn"/> when given and present
    /// (else from a <c>gene_symbols</c>/<c>gene_name</c>/<c>feature_name</c> column, else the var index).
    /// Unsupported annotation types are skipped.
    /// </summary>
    public SingleCellExperiment ToSingleCellExperiment(string? geneSymbolColumn = null)
    {
        if (X is null) throw new InvalidOperationException("The AnnData object has no X.");
        string? symCol = geneSymbolColumn is not null && Var.Contains(geneSymbolColumn) ? geneSymbolColumn
            : new[] { "gene_symbols", "gene_symbol", "gene_name", "feature_name" }.FirstOrDefault(Var.Contains);
        var names = symCol is null ? Var.Index : Var.GetStrings(symCol).Select((s, i) => s ?? Var.Index[i]).ToArray();
        var sc = new SingleCellExperiment(X.ToSparse(), Obs.Index, Var.Index, names);
        foreach (var (k, m) in Layers) sc.Layers[k] = m.ToSparse();
        foreach (var name in Obs.ColumnNames)
        {
            switch (Obs[name])
            {
                case Categorical c: sc.CellLabels[name] = c.ToStrings().Select(s => s ?? "").ToArray(); break;
                case string[] s: sc.CellLabels[name] = s; break;
                case MaskedArray { Values: string[] } m: sc.CellLabels[name] = Obs.GetStrings(name).Select(s => s ?? "").ToArray(); break;
                case bool[] b: sc.CellMetrics[name] = b.Select(x => x ? 1.0 : 0.0).ToArray(); break;
                default: sc.CellMetrics[name] = Obs.GetDoubles(name); break;
            }
        }
        foreach (var name in Var.ColumnNames)
        {
            switch (Var[name])
            {
                case bool[] b: sc.GeneFlags[name] = b; break;
                case Categorical or string[] or MaskedArray { Values: string[] }: break;
                case MaskedArray { Values: bool[] } mb: sc.GeneFlags[name] = Enumerable.Range(0, mb.Length).Select(i => !mb.Mask[i] && ((bool[])mb.Values)[i]).ToArray(); break;
                default: sc.GeneMetrics[name] = Var.GetDoubles(name); break;
            }
        }
        foreach (var (k, v) in Obsm) if (v is AnnMatrix m) sc.Embeddings[k] = m.ToDense();
        return sc;
    }

    /// <summary>
    /// Builds an AnnData object from a <see cref="SingleCellExperiment"/>: X and layers as CSR float32, cell
    /// metrics as float64 columns, cell labels as categoricals, gene metrics/flags as var columns (plus
    /// <c>gene_symbols</c> when symbols differ from ids) and embeddings in <c>obsm</c>. The neighbour graph is
    /// written to <c>obsp["distances"]</c>/<c>obsp["connectivities"]</c> when present.
    /// </summary>
    public static AnnData FromSingleCellExperiment(SingleCellExperiment sc)
    {
        ArgumentNullException.ThrowIfNull(sc);
        var obs = new AnnDataFrame(sc.CellNames);
        foreach (var (k, v) in sc.CellMetrics) obs.Set(k, (double[])v.Clone());
        foreach (var (k, v) in sc.CellLabels) obs.Set(k, Categorical.FromLabels(v));
        var var = new AnnDataFrame(sc.GeneIds);
        if (!sc.GeneNames.SequenceEqual(sc.GeneIds)) var.Set("gene_symbols", sc.GeneNames.ToArray());
        foreach (var (k, v) in sc.GeneMetrics) var.Set(k, (double[])v.Clone());
        foreach (var (k, v) in sc.GeneFlags) var.Set(k, (bool[])v.Clone());
        var a = new AnnData(AnnMatrix.FromSparse(sc.X), obs, var);
        foreach (var (k, v) in sc.Layers) a.Layers[k] = AnnMatrix.FromSparse(v);
        foreach (var (k, v) in sc.Embeddings) a.Obsm[k] = AnnMatrix.FromDense(v);
        if (sc.Neighbors is NeighborGraph g)
        {
            var rows = new List<int>(); var cols = new List<int>(); var vals = new List<float>();
            for (int i = 0; i < g.CellCount; i++)
                for (int t = 0; t < g.K; t++) { rows.Add(i); cols.Add(g.Indices[i * g.K + t]); vals.Add((float)g.Distances[i * g.K + t]); }
            a.Obsp["distances"] = AnnMatrix.FromSparse(SparseMatrix.FromTriplets(g.CellCount, g.CellCount, rows.ToArray(), cols.ToArray(), vals.ToArray()), dataType: "float64");
            a.Obsp["connectivities"] = AnnMatrix.FromSparse(g.Connectivities, dataType: "float64");
        }
        return a;
    }

    /// <inheritdoc />
    public override string ToString()
    {
        var parts = new List<string> { $"AnnData object with n_obs × n_vars = {ObsCount} × {VarCount}" };
        if (Obs.ColumnNames.Count > 0) parts.Add($"obs: {string.Join(", ", Obs.ColumnNames)}");
        if (Var.ColumnNames.Count > 0) parts.Add($"var: {string.Join(", ", Var.ColumnNames)}");
        if (Uns.Count > 0) parts.Add($"uns: {string.Join(", ", Uns.Keys)}");
        if (Obsm.Count > 0) parts.Add($"obsm: {string.Join(", ", Obsm.Keys)}");
        if (Varm.Count > 0) parts.Add($"varm: {string.Join(", ", Varm.Keys)}");
        if (Layers.Count > 0) parts.Add($"layers: {string.Join(", ", Layers.Keys)}");
        if (Obsp.Count > 0) parts.Add($"obsp: {string.Join(", ", Obsp.Keys)}");
        if (Varp.Count > 0) parts.Add($"varp: {string.Join(", ", Varp.Keys)}");
        return string.Join(Environment.NewLine + "    ", parts);
    }
}
