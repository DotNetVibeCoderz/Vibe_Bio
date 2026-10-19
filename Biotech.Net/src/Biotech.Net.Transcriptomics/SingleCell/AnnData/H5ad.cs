using PureHDF;
using PureHDF.Filters;
using PureHDF.Selections;
using PureHDF.VOL.Native;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>Compression applied to large datasets when writing <c>.h5ad</c> files.</summary>
public enum H5adCompression
{
    /// <summary>No compression (contiguous datasets; fastest to write and read).</summary>
    None,

    /// <summary>gzip/deflate with the byte-shuffle filter, the same filters as <c>anndata.write_h5ad(compression="gzip")</c>.</summary>
    Gzip,
}

/// <summary>Options for <see cref="H5ad.Write(string, AnnData, H5adWriteOptions?)"/> and <see cref="H5ad.WriteStreaming"/>.</summary>
public sealed class H5adWriteOptions
{
    /// <summary>Compression for datasets with at least <see cref="CompressionThreshold"/> elements.</summary>
    public H5adCompression Compression { get; init; } = H5adCompression.None;

    /// <summary>Deflate level. PureHDF supports 1 (fast) and 9 (small).</summary>
    public int GzipLevel { get; init; } = 1;

    /// <summary>Apply the shuffle filter before deflate (better ratios for numeric data).</summary>
    public bool Shuffle { get; init; } = true;

    /// <summary>Chunk length (elements) of compressed 1-D datasets; 2-D datasets use row chunks of about this many elements.</summary>
    public int ChunkElements { get; init; } = 1 << 16;

    /// <summary>Datasets smaller than this many elements are written contiguous and uncompressed.</summary>
    public int CompressionThreshold { get; init; } = 4096;

    internal static H5adWriteOptions Default { get; } = new();
}

/// <summary>
/// Reads and writes AnnData <c>.h5ad</c> files (the HDF5 layout of Python's <c>anndata</c>; on-disk spec
/// element encodings <c>array</c>, <c>csr_matrix</c>, <c>csc_matrix</c>, <c>dataframe</c> 0.1.0/0.2.0,
/// <c>categorical</c>, <c>string</c>, <c>string-array</c>, <c>numeric-scalar</c>, <c>nullable-integer</c>,
/// <c>nullable-boolean</c>, <c>nullable-string-array</c>, <c>dict</c> and <c>raw</c>, plus the pre-0.7
/// <c>h5sparse</c> sparse layout). HDF5 access is pure managed code (PureHDF): no native HDF5 library is needed.
/// </summary>
/// <remarks>
/// Not supported: awkward arrays, compound (record) arrays such as the structured arrays older scanpy versions
/// wrote for <c>rank_genes_groups</c> (they are read as <see cref="AnnUnsupported"/> placeholders and not written
/// back), the pre-0.7 compound-dataset <c>obs</c>/<c>var</c> layout, zarr stores and modifying files in place.
/// Sparse values are held as float32 in memory; dense matrices as float64.
/// </remarks>
public static class H5ad
{
    // ================================================================== reading

    /// <summary>Reads a whole <c>.h5ad</c> file into memory.</summary>
    public static AnnData Read(string path)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        using var file = H5File.OpenRead(path);
        return ReadRoot(file);
    }

    /// <summary>Reads an <c>.h5ad</c> file and converts it with <see cref="AnnData.ToSingleCellExperiment"/>.</summary>
    public static SingleCellExperiment ReadSingleCell(string path, string? geneSymbolColumn = null)
        => Read(path).ToSingleCellExperiment(geneSymbolColumn);

    /// <summary>
    /// Opens a file in backed mode: <c>obs</c>, <c>var</c> and <c>uns</c> are loaded, X and the layers stay on disk
    /// and are read in row (or column) slices on demand. Dispose the result to close the file.
    /// </summary>
    public static BackedAnnData OpenBacked(string path) => new(path);

    internal static AnnData ReadRoot(IH5Group root)
    {
        var obs = root.LinkExists("obs") ? ReadDataFrame(root.Get("obs"), "obs") : null;
        var var = root.LinkExists("var") ? ReadDataFrame(root.Get("var"), "var") : null;
        AnnMatrix? x = root.LinkExists("X") ? ReadMatrix(root.Get("X"), "X") : null;
        obs ??= new AnnDataFrame(Enumerable.Range(0, x?.Rows ?? 0).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray());
        var ??= new AnnDataFrame(Enumerable.Range(0, x?.Columns ?? 0).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray());
        var a = new AnnData(x, obs, var);
        foreach (var (k, v) in ReadMapping(root, "obsm")) a.Obsm[k] = v;
        foreach (var (k, v) in ReadMapping(root, "varm")) a.Varm[k] = v;
        foreach (var (k, v) in ReadMapping(root, "obsp")) if (v is AnnMatrix m) a.Obsp[k] = m;
        foreach (var (k, v) in ReadMapping(root, "varp")) if (v is AnnMatrix m) a.Varp[k] = m;
        foreach (var (k, v) in ReadMapping(root, "layers")) if (v is AnnMatrix m) a.Layers[k] = m;
        if (root.LinkExists("uns") && root.Get("uns") is IH5Group uns)
            foreach (var (k, v) in ReadDict(uns, "uns")) a.Uns[k] = v;
        if (root.LinkExists("raw") && root.Get("raw") is IH5Group raw && raw.LinkExists("X"))
        {
            var rx = ReadMatrix(raw.Get("X"), "raw/X");
            var rv = raw.LinkExists("var") ? ReadDataFrame(raw.Get("var"), "raw/var") : new AnnDataFrame(Enumerable.Range(0, rx.Columns).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray());
            a.Raw = new AnnRaw(rx, rv);
            foreach (var (k, v) in ReadMapping(raw, "varm")) if (v is AnnMatrix m) a.Raw.Varm[k] = m;
        }
        return a;
    }

    private static IEnumerable<(string, object)> ReadMapping(IH5Group parent, string name)
    {
        if (!parent.LinkExists(name) || parent.Get(name) is not IH5Group g) yield break;
        foreach (var child in g.Children())
        {
            string path = $"{name}/{child.Name}";
            object v = Encoding(child) == "dataframe" ? ReadDataFrame(child, path) : ReadMatrix(child, path);
            yield return (child.Name, v);
        }
    }

    internal static string? Encoding(IH5Object o) => H5Types.StringAttr(o, "encoding-type");

    /// <summary>Reads a matrix element: a 2-D (or 1-D, read as a column) dataset or a CSR/CSC group.</summary>
    internal static AnnMatrix ReadMatrix(IH5Object o, string path)
    {
        if (o is IH5Dataset ds)
        {
            var dims = ds.Space.Dimensions;
            string dtype = H5Types.DTypeName(ds.Type) ?? throw new NotSupportedException($"'{path}' is not a numeric array.");
            var flat = H5Types.ToDoubles(H5Types.ReadArray(ds));
            int rows = dims.Length > 0 ? (int)dims[0] : 1;
            int cols = dims.Length > 1 ? (int)(dims.Skip(1).Aggregate(1UL, (p, q) => p * q)) : 1;
            return AnnMatrix.FromDense(new DenseMatrix(rows, cols, flat), dtype);
        }
        if (o is IH5Group g)
        {
            var (format, shape) = SparseHeader(g, path);
            var data = (IH5Dataset)g.Get("data");
            string dtype = H5Types.DTypeName(data.Type) ?? "float32";
            long nnz = (long)data.Space.Dimensions[0];
            int majorCount = format == AnnMatrixFormat.Csr ? shape.Rows : shape.Columns;
            var ptr = H5Types.ReadInt64((IH5Dataset)g.Get("indptr"), 0, (ulong)majorCount + 1);
            if (ptr[^1] > int.MaxValue) throw new NotSupportedException($"'{path}' has {ptr[^1]} stored entries; open the file with H5ad.OpenBacked and process it out of core.");
            var idx = nnz == 0 ? [] : H5Types.ReadInt32((IH5Dataset)g.Get("indices"), 0, (ulong)nnz);
            bool wide = dtype is "float64" or "int64" or "uint64" or "uint32" or "int32";
            double[]? v64 = nnz == 0 ? null : wide ? H5Types.ToDoubles(H5Types.ReadArray(data)) : null;
            var vals = v64 is not null ? Array.ConvertAll(v64, v => (float)v) : nnz == 0 ? [] : H5Types.ReadSingle(data, null);
            var p32 = Array.ConvertAll(ptr, v => (int)v);
            // Sort (and transpose CSC) a permutation so the float64 copy stays aligned.
            var perm = Permutation(vals.Length);
            SortIndices(p32, idx, perm);
            var pm = format == AnnMatrixFormat.Csr
                ? new SparseMatrix(shape.Rows, shape.Columns, p32, idx, perm)
                : new SparseMatrix(shape.Columns, shape.Rows, p32, idx, perm).Transpose();
            int count = pm.NonZeroCount;
            var outVals = new float[count];
            double[]? out64 = v64 is null ? null : new double[count];
            for (int k = 0; k < count; k++)
            {
                int src = BitConverter.SingleToInt32Bits(pm.Values[k]);
                outVals[k] = vals[src];
                if (out64 is not null) out64[k] = v64![src];
            }
            var m = new SparseMatrix(pm.Rows, pm.Columns, pm.RowPointers, pm.ColumnIndices, outVals);
            bool lossless = out64 is null || out64.All(v => (double)(float)v == v);
            return AnnMatrix.FromSparse(m, format, dtype, lossless ? null : out64);
        }
        throw new InvalidDataException($"'{path}' is not a matrix element.");
    }

    /// <summary>Entry positions 0, 1, 2 … carried through sorting/transposition as float bit patterns.</summary>
    private static float[] Permutation(int n)
    {
        var perm = new float[n];
        for (int k = 0; k < n; k++) perm[k] = BitConverter.Int32BitsToSingle(k);
        return perm;
    }

    internal static (AnnMatrixFormat Format, (int Rows, int Columns) Shape) SparseHeader(IH5Group g, string path)
    {
        string? enc = Encoding(g) ?? H5Types.StringAttr(g, "h5sparse_format") switch { "csr" => "csr_matrix", "csc" => "csc_matrix", var s => s };
        var format = enc switch
        {
            "csr_matrix" => AnnMatrixFormat.Csr,
            "csc_matrix" => AnnMatrixFormat.Csc,
            _ => throw new NotSupportedException($"'{path}' has unsupported encoding '{enc}'."),
        };
        var shape = H5Types.Int64ArrayAttr(g, "shape") ?? H5Types.Int64ArrayAttr(g, "h5sparse_shape")
            ?? throw new InvalidDataException($"Sparse matrix '{path}' has no shape attribute.");
        return (format, ((int)shape[0], (int)shape[1]));
    }

    /// <summary>scipy allows unsorted indices within a row; SparseMatrix requires them sorted.</summary>
    internal static void SortIndices(int[] ptr, int[] idx, float[] vals)
    {
        for (int r = 0; r + 1 < ptr.Length; r++)
        {
            int s = ptr[r], e = ptr[r + 1];
            for (int k = s + 1; k < e; k++)
                if (idx[k] < idx[k - 1]) { Array.Sort(idx, vals, s, e - s); break; }
        }
    }

    /// <summary>Reads a dataframe group (spec 0.1.0 or 0.2.0).</summary>
    internal static AnnDataFrame ReadDataFrame(IH5Object o, string path)
    {
        if (o is not IH5Group g)
            throw new NotSupportedException($"'{path}' is a compound dataset (anndata < 0.7 layout); re-save the file with a current anndata.");
        string indexKey = H5Types.StringAttr(g, "_index") ?? "_index";
        string[] index;
        if (g.LinkExists(indexKey)) index = AsStrings(ReadElement(g.Get(indexKey), $"{path}/{indexKey}"), $"{path}/{indexKey}");
        else index = [];
        var order = H5Types.StringArrayAttr(g, "column-order");
        if (index.Length == 0 && order.Length > 0)
        {
            int n = ElementLength(g.Get(order[0]));
            index = Enumerable.Range(0, n).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        }
        var frame = new AnnDataFrame(index, indexKey == "_index" ? "_index" : indexKey);
        IH5Group? legacyCats = g.LinkExists("__categories") ? g.Get("__categories") as IH5Group : null;
        foreach (var col in order)
        {
            if (!g.LinkExists(col)) continue;
            var child = g.Get(col);
            object? v = ReadElement(child, $"{path}/{col}");
            if (legacyCats is not null && legacyCats.LinkExists(col) && v is Array codes && codes is not string[])
            {
                var catDs = legacyCats.Get(col);
                var cats = H5Types.ReadArray((IH5Dataset)catDs);
                v = new Categorical(H5Types.ToInt32(codes), cats, H5Types.BoolAttr(catDs, "ordered"));
            }
            if (v is null or AnnUnsupported or AnnArray or Dictionary<string, object?>) continue;
            frame.Set(col, v);
        }
        return frame;
    }

    private static int ElementLength(IH5Object o) => o switch
    {
        IH5Dataset d => (int)d.Space.Dimensions[0],
        IH5Group g when g.LinkExists("codes") => (int)((IH5Dataset)g.Get("codes")).Space.Dimensions[0],
        IH5Group g when g.LinkExists("values") => (int)((IH5Dataset)g.Get("values")).Space.Dimensions[0],
        _ => 0,
    };

    private static string[] AsStrings(object? v, string path) => v switch
    {
        string[] s => s,
        MaskedArray { Values: string[] s } m => s.Select((x, i) => m.Mask[i] ? "" : x).ToArray(),
        Categorical c => c.ToStrings().Select(x => x ?? "").ToArray(),
        Array a => Enumerable.Range(0, a.Length).Select(i => Convert.ToString(a.GetValue(i), System.Globalization.CultureInfo.InvariantCulture) ?? "").ToArray(),
        _ => throw new InvalidDataException($"'{path}' is not a string index."),
    };

    private static Dictionary<string, object?> ReadDict(IH5Group g, string path)
    {
        var d = new Dictionary<string, object?>(StringComparer.Ordinal);
        foreach (var child in g.Children()) d[child.Name] = ReadElement(child, $"{path}/{child.Name}");
        return d;
    }

    /// <summary>Reads any element by its <c>encoding-type</c> (falling back to the HDF5 structure).</summary>
    internal static object? ReadElement(IH5Object o, string path)
    {
        string? enc = Encoding(o);
        if (o is IH5Group g)
        {
            switch (enc)
            {
                case "dataframe": return ReadDataFrame(g, path);
                case "csr_matrix" or "csc_matrix": return ReadMatrix(g, path);
                case "categorical":
                    {
                        var codes = H5Types.ToInt32(H5Types.ReadArray((IH5Dataset)g.Get("codes")));
                        var cats = H5Types.ReadArray((IH5Dataset)g.Get("categories"));
                        return new Categorical(codes, cats, H5Types.BoolAttr(g, "ordered"));
                    }
                case "nullable-integer" or "nullable-boolean" or "nullable-string-array":
                    {
                        var values = H5Types.ReadArray((IH5Dataset)g.Get("values"));
                        var mask = H5Types.ReadArray((IH5Dataset)g.Get("mask")) as bool[] ?? throw new InvalidDataException($"'{path}/mask' is not boolean.");
                        return new MaskedArray(values, mask) { NaValue = H5Types.StringAttr(g, "na-value") };
                    }
                case "dict" or null when H5Types.StringAttr(g, "h5sparse_format") is not null: return ReadMatrix(g, path);
                case "dict" or null or "raw": return ReadDict(g, path);
                default: return new AnnUnsupported(path, enc, "unsupported group encoding");
            }
        }
        if (o is IH5Dataset ds)
        {
            var t = ds.Type;
            if (t.Class == H5DataTypeClass.Compound) return new AnnUnsupported(path, enc ?? "", "compound (record) arrays are not supported");
            if (t.Class == H5DataTypeClass.Reference) return new AnnUnsupported(path, enc ?? "", "object references are not supported");
            if (H5Types.IsString(t))
            {
                if (ds.Space.Type == H5DataspaceType.Scalar) return ds.Read<string>();
                var s = ds.Read<string[]>();
                return ds.Space.Dimensions.Length > 1 ? new AnnArray(s, ds.Space.Dimensions.Select(x => (long)x).ToArray()) : s;
            }
            if (H5Types.DTypeName(t) is null) return new AnnUnsupported(path, enc ?? "", $"HDF5 type class {t.Class} is not supported");
            var arr = H5Types.ReadArray(ds);
            if (ds.Space.Type == H5DataspaceType.Scalar)
            {
                object v = arr.GetValue(0)!;
                return v switch { int i => (long)i, short s => (long)s, sbyte b => (long)b, float f => (double)f, _ => v };
            }
            if (ds.Space.Dimensions.Length > 1) return new AnnArray(arr, ds.Space.Dimensions.Select(x => (long)x).ToArray());
            return arr;
        }
        return null;
    }

    // ================================================================== writing

    /// <summary>
    /// Writes an AnnData object to an <c>.h5ad</c> file readable by Python <c>anndata</c> ≥ 0.8 (and by
    /// <see cref="Read"/>). Strings are UTF-8 variable-length, booleans use the numpy/h5py enum type, sparse
    /// matrices are written in their recorded <see cref="AnnMatrix.Format"/> with int64 <c>indptr</c>.
    /// </summary>
    public static void Write(string path, AnnData data, H5adWriteOptions? options = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        ArgumentNullException.ThrowIfNull(data);
        var o = options ?? H5adWriteOptions.Default;
        var file = new H5File { Attributes = new() { ["encoding-type"] = "anndata", ["encoding-version"] = "0.1.0" } };
        if (data.X is not null) file["X"] = MatrixNode(data.X, o);
        file["obs"] = DataFrameNode(data.Obs, o);
        file["var"] = DataFrameNode(data.Var, o);
        file["obsm"] = MappingNode(data.Obsm, o);
        file["varm"] = MappingNode(data.Varm, o);
        file["obsp"] = MappingNode(data.Obsp.ToDictionary(kv => kv.Key, kv => (object)kv.Value), o);
        file["varp"] = MappingNode(data.Varp.ToDictionary(kv => kv.Key, kv => (object)kv.Value), o);
        file["layers"] = MappingNode(data.Layers.ToDictionary(kv => kv.Key, kv => (object)kv.Value), o);
        file["uns"] = DictNode(data.Uns, o, "uns");
        if (data.Raw is not null)
        {
            var raw = new H5Group { Attributes = new() { ["encoding-type"] = "raw", ["encoding-version"] = "0.1.0" } };
            raw["X"] = MatrixNode(data.Raw.X, o);
            raw["var"] = DataFrameNode(data.Raw.Var, o);
            raw["varm"] = MappingNode(data.Raw.Varm.ToDictionary(kv => kv.Key, kv => (object)kv.Value), o);
            file["raw"] = raw;
        }
        var dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
        file.Write(path);
    }

    /// <summary>Writes a <see cref="SingleCellExperiment"/> (via <see cref="AnnData.FromSingleCellExperiment"/>).</summary>
    public static void Write(string path, SingleCellExperiment sc, H5adWriteOptions? options = null)
        => Write(path, AnnData.FromSingleCellExperiment(sc), options);

    /// <summary>
    /// Writes an <c>.h5ad</c> whose X (and optional layers) are streamed from out-of-core sources in row blocks, so
    /// matrices larger than memory can be exported (e.g. a <see cref="MappedSparseMatrix"/> after
    /// <see cref="OutOfCore.NormalizeTotalLog1p"/>). X and layers are written as CSR with int64 <c>indptr</c>;
    /// <paramref name="obs"/>, <paramref name="var"/>, <paramref name="obsm"/> and <paramref name="uns"/> are written
    /// from memory.
    /// </summary>
    public static void WriteStreaming(string path, ISparseMatrixSource x, AnnDataFrame obs, AnnDataFrame var,
        IReadOnlyDictionary<string, ISparseMatrixSource>? layers = null, IReadOnlyDictionary<string, AnnMatrix>? obsm = null,
        IReadOnlyDictionary<string, object?>? uns = null, H5adWriteOptions? options = null, int blockRows = SparseMatrixSource.DefaultBlockRows)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(obs);
        ArgumentNullException.ThrowIfNull(var);
        if (obs.RowCount != x.Rows || var.RowCount != x.Columns) throw new ArgumentException("obs/var do not match the shape of X.");
        var o = options ?? H5adWriteOptions.Default;
        var file = new H5File { Attributes = new() { ["encoding-type"] = "anndata", ["encoding-version"] = "0.1.0" } };
        var pending = new List<(ISparseMatrixSource Source, H5Dataset<float[]> Data, H5Dataset<int[]> Indices, H5Dataset<long[]> Indptr)>();
        H5Group Deferred(ISparseMatrixSource m)
        {
            ulong nnz = (ulong)m.NonZeroCount;
            uint chunk = (uint)Math.Max(1, Math.Min((ulong)o.ChunkElements, Math.Max(1UL, nnz)));
            uint pchunk = (uint)Math.Max(1, Math.Min(o.ChunkElements, m.Rows + 1));
            // One bounded cache per dataset: PureHDF's default writing cache keeps every chunk until the file closes.
            H5DatasetCreation Cr() => new(new SequentialWritingChunkCache(4), o.Compression == H5adCompression.Gzip ? Creation(o).Filters : []);
            var data = new H5Dataset<float[]>(fileDims: [nnz], chunks: [chunk], datasetCreation: Cr());
            var idx = new H5Dataset<int[]>(fileDims: [nnz], chunks: [chunk], datasetCreation: Cr());
            var ptr = new H5Dataset<long[]>(fileDims: [(ulong)m.Rows + 1], chunks: [pchunk], datasetCreation: Cr());
            pending.Add((m, data, idx, ptr));
            return new H5Group
            {
                Attributes = new() { ["encoding-type"] = "csr_matrix", ["encoding-version"] = "0.1.0", ["shape"] = new long[] { m.Rows, m.Columns } },
                ["data"] = data, ["indices"] = idx, ["indptr"] = ptr,
            };
        }
        file["X"] = Deferred(x);
        file["obs"] = DataFrameNode(obs, o);
        file["var"] = DataFrameNode(var, o);
        var lg = new H5Group { Attributes = Enc("dict", "0.1.0") };
        foreach (var (k, m) in layers ?? new Dictionary<string, ISparseMatrixSource>())
        {
            if (m.Rows != x.Rows || m.Columns != x.Columns) throw new ArgumentException($"Layer '{k}' does not match the shape of X.");
            lg[k] = Deferred(m);
        }
        file["layers"] = lg;
        file["obsm"] = MappingNode((obsm ?? new Dictionary<string, AnnMatrix>()).ToDictionary(kv => kv.Key, kv => (object)kv.Value), o);
        file["varm"] = MappingNode([], o);
        file["obsp"] = MappingNode([], o);
        file["varp"] = MappingNode([], o);
        file["uns"] = DictNode(uns ?? new Dictionary<string, object?>(), o, "uns");
        var dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
        using var writer = file.BeginWrite(path);
        foreach (var (src, data, idx, ptrDs) in pending)
        {
            var ptr = new long[src.Rows + 1];
            long at = 0;
            foreach (var (start, b) in src.RowBlocks(blockRows))
            {
                int nnz = b.NonZeroCount;
                for (int i = 0; i < b.Rows; i++) ptr[start + i + 1] = at + b.RowPointers[i + 1];
                if (nnz > 0)
                {
                    var sel = new HyperslabSelection((ulong)at, (ulong)nnz);
                    writer.Write(data, nnz == b.Values.Length ? b.Values : b.Values[..nnz], fileSelection: sel);
                    writer.Write(idx, nnz == b.ColumnIndices.Length ? b.ColumnIndices : b.ColumnIndices[..nnz], fileSelection: sel);
                }
                at += nnz;
            }
            writer.Write(ptrDs, ptr);
        }
    }

    internal static Dictionary<string, object> Enc(string type, string version) => new() { ["encoding-type"] = type, ["encoding-version"] = version };

    internal static H5DatasetCreation Creation(H5adWriteOptions o)
    {
        var filters = new List<H5Filter>();
        if (o.Shuffle) filters.Add(new H5Filter(ShuffleFilter.Id));
        filters.Add(new H5Filter(DeflateFilter.Id, new() { [DeflateFilter.COMPRESSION_LEVEL] = o.GzipLevel }));
        return new H5DatasetCreation(Filters: filters);
    }

    /// <summary>A dataset node, chunked and compressed when large enough.</summary>
    internal static H5Dataset DatasetNode(Array data, H5adWriteOptions o, ulong[]? dims = null, Dictionary<string, object>? attrs = null)
    {
        H5Dataset d;
        long n = data.Length;
        bool compress = o.Compression == H5adCompression.Gzip && n >= o.CompressionThreshold && data is not string[];
        if (compress)
        {
            uint[] chunks = dims is { Length: 2 }
                ? [(uint)Math.Clamp(o.ChunkElements / (long)Math.Max(1UL, dims[1]), 1, (long)dims[0]), (uint)dims[1]]
                : [(uint)Math.Min(n, o.ChunkElements)];
            d = new H5Dataset(data, chunks: chunks, fileDims: dims, datasetCreation: Creation(o));
        }
        else d = new H5Dataset(data, fileDims: dims);
        if (attrs is not null) d.Attributes = attrs;
        return d;
    }

    internal static object MatrixNode(AnnMatrix m, H5adWriteOptions o)
    {
        if (m.Dense is DenseMatrix dm)
            return DatasetNode(H5Types.ToStorage(dm.Data, m.DataType), o, [(ulong)dm.Rows, (ulong)dm.Columns], Enc("array", "0.2.0"));
        var s = m.Sparse!;
        double[]? exact = m.SparseValues64;
        if (m.Format == AnnMatrixFormat.Csc)
        {
            if (exact is not null)
            {
                var t = new SparseMatrix(s.Rows, s.Columns, s.RowPointers, s.ColumnIndices, Permutation(s.NonZeroCount)).Transpose();
                var e2 = new double[t.NonZeroCount];
                for (int k = 0; k < e2.Length; k++) e2[k] = exact[BitConverter.SingleToInt32Bits(t.Values[k])];
                exact = e2;
                s = new SparseMatrix(t.Rows, t.Columns, t.RowPointers, t.ColumnIndices, Array.ConvertAll(e2, v => (float)v));
            }
            else s = s.Transpose();
        }
        return SparseNode(s, m.Format, m.Rows, m.Columns, m.DataType, o, exact);
    }

    /// <summary><paramref name="major"/> is CSR for csr_matrix or the transposed (columns × rows) CSR for csc_matrix.</summary>
    private static H5Group SparseNode(SparseMatrix major, AnnMatrixFormat format, int rows, int columns, string dtype, H5adWriteOptions o, double[]? exact = null)
    {
        int nnz = major.NonZeroCount;
        var g = new H5Group
        {
            Attributes = new()
            {
                ["encoding-type"] = format == AnnMatrixFormat.Csc ? "csc_matrix" : "csr_matrix",
                ["encoding-version"] = "0.1.0",
                ["shape"] = new long[] { rows, columns },
            },
            ["data"] = DatasetNode(exact is not null ? H5Types.ToStorage(exact.AsSpan(0, nnz), dtype) : H5Types.ToStorage(major.Values.AsSpan(0, nnz), dtype), o),
            ["indices"] = DatasetNode(major.ColumnIndices.AsSpan(0, nnz).ToArray(), o),
            ["indptr"] = DatasetNode(Array.ConvertAll(major.RowPointers, v => (long)v), o),
        };
        return g;
    }

    internal static H5Group DataFrameNode(AnnDataFrame f, H5adWriteOptions o)
    {
        string indexKey = string.IsNullOrEmpty(f.IndexName) ? "_index" : f.IndexName;
        if (f.Contains(indexKey)) indexKey = "_index";
        var g = new H5Group
        {
            Attributes = new()
            {
                ["encoding-type"] = "dataframe",
                ["encoding-version"] = "0.2.0",
                ["_index"] = indexKey,
                ["column-order"] = f.ColumnNames.ToArray(),
            },
        };
        g[indexKey] = StringArrayNode(f.Index, o);
        foreach (var name in f.ColumnNames) g[name] = ColumnNode(f[name], o);
        return g;
    }

    private static H5Dataset StringArrayNode(string[] s, H5adWriteOptions o)
        => DatasetNode(s.Select(x => x ?? "").ToArray(), o, attrs: Enc("string-array", "0.2.0"));

    private static object ColumnNode(object v, H5adWriteOptions o) => v switch
    {
        Categorical c => new H5Group
        {
            Attributes = new() { ["encoding-type"] = "categorical", ["encoding-version"] = "0.2.0", ["ordered"] = c.Ordered ? H5Bool.TRUE : H5Bool.FALSE },
            ["codes"] = DatasetNode(CompactCodes(c), o, attrs: Enc("array", "0.2.0")),
            ["categories"] = c.Categories is string[] cs ? StringArrayNode(cs, o) : DatasetNode(H5Types.ForWrite(c.Categories), o, attrs: Enc("array", "0.2.0")),
        },
        MaskedArray m => new H5Group
        {
            Attributes = m.Values switch
            {
                bool[] => Enc("nullable-boolean", "0.1.0"),
                string[] => new Dictionary<string, object> { ["encoding-type"] = "nullable-string-array", ["encoding-version"] = "0.1.0", ["na-value"] = m.NaValue ?? "NaN" },
                _ => Enc("nullable-integer", "0.1.0"),
            },
            ["values"] = m.Values is string[] sv ? StringArrayNode(sv.Select((x, i) => m.Mask[i] ? "" : x).ToArray(), o) : DatasetNode(H5Types.ForWrite(m.Values), o, attrs: Enc("array", "0.2.0")),
            ["mask"] = DatasetNode(H5Types.ForWrite(m.Mask), o, attrs: Enc("array", "0.2.0")),
        },
        string[] s => StringArrayNode(s, o),
        Array a => DatasetNode(H5Types.ForWrite(a), o, attrs: Enc("array", "0.2.0")),
        _ => throw new NotSupportedException($"Cannot write a column of type {v.GetType().Name}."),
    };

    private static Array CompactCodes(Categorical c)
    {
        int n = c.Categories.Length;
        if (n < sbyte.MaxValue) return Array.ConvertAll(c.Codes, x => (sbyte)x);
        if (n < short.MaxValue) return Array.ConvertAll(c.Codes, x => (short)x);
        return c.Codes;
    }

    private static H5Group MappingNode(Dictionary<string, object> items, H5adWriteOptions o)
    {
        var g = new H5Group { Attributes = Enc("dict", "0.1.0") };
        foreach (var (k, v) in items)
            g[k] = v switch
            {
                AnnMatrix m => MatrixNode(m, o),
                AnnDataFrame f => DataFrameNode(f, o),
                DenseMatrix d => MatrixNode(AnnMatrix.FromDense(d), o),
                SparseMatrix s => MatrixNode(AnnMatrix.FromSparse(s), o),
                _ => throw new NotSupportedException($"Cannot write '{k}' of type {v.GetType().Name}."),
            };
        return g;
    }

    private static H5Group DictNode(IReadOnlyDictionary<string, object?> items, H5adWriteOptions o, string path)
    {
        var g = new H5Group { Attributes = Enc("dict", "0.1.0") };
        foreach (var (k, v) in items)
        {
            object? node = UnsNode(v, o, $"{path}/{k}");
            if (node is not null) g[k] = node;
        }
        return g;
    }

    private static object? UnsNode(object? v, H5adWriteOptions o, string path) => v switch
    {
        null or AnnUnsupported => null,
        string s => new H5Dataset(s) { Attributes = Enc("string", "0.2.0") },
        bool b => new H5Dataset(b ? H5Bool.TRUE : H5Bool.FALSE) { Attributes = Enc("numeric-scalar", "0.2.0") },
        double or float or long or int or short or sbyte or byte or ushort or uint or ulong => new H5Dataset(v) { Attributes = Enc("numeric-scalar", "0.2.0") },
        IReadOnlyDictionary<string, object?> d => DictNode(d, o, path),
        IDictionary<string, object> d2 => DictNode(d2.ToDictionary(kv => kv.Key, kv => (object?)kv.Value), o, path),
        AnnMatrix m => MatrixNode(m, o),
        DenseMatrix dm => MatrixNode(AnnMatrix.FromDense(dm), o),
        SparseMatrix sm => MatrixNode(AnnMatrix.FromSparse(sm), o),
        AnnDataFrame f => DataFrameNode(f, o),
        Categorical or MaskedArray => ColumnNode(v, o),
        AnnArray a => a.Data is string[] sa
            ? DatasetNode(sa, o, a.Shape.Select(x => (ulong)x).ToArray(), Enc("string-array", "0.2.0"))
            : DatasetNode(H5Types.ForWrite(a.Data), o, a.Shape.Select(x => (ulong)x).ToArray(), Enc("array", "0.2.0")),
        string[] sa => StringArrayNode(sa, o),
        IEnumerable<string> se => StringArrayNode(se.ToArray(), o),
        Array a when a.Rank == 1 && AnnDataFrame.IsSupportedArray(a) => DatasetNode(H5Types.ForWrite(a), o, attrs: Enc("array", "0.2.0")),
        _ => throw new NotSupportedException($"Cannot write uns element '{path}' of type {v.GetType().Name}."),
    };
}

/// <summary>
/// A writing chunk cache for sequential writes: keeps the most recent <c>capacity</c> chunks and writes older ones
/// out, so streamed datasets need O(capacity × chunk size) memory instead of the whole dataset (PureHDF's default
/// writing cache is unbounded). Chunks must not be revisited after eviction, which holds for the strictly
/// increasing hyperslabs <see cref="H5ad.WriteStreaming"/> issues.
/// </summary>
/// <remarks>
/// PureHDF 2.2.0's <c>H5D_Chunk.GetWriteStream</c> passes an eviction callback that ignores its index argument and
/// writes under the index of the chunk being requested, which would corrupt the file. Its <c>Flush</c> callback is
/// correct, and PureHDF calls <see cref="Flush"/> once when the dataset is created, so that delegate is kept and used
/// for evictions. Until it has been seen, nothing is evicted.
/// </remarks>
internal sealed class SequentialWritingChunkCache(int capacity) : IWritingChunkCache
{
    private readonly LinkedList<(ulong Index, Memory<byte> Chunk)> _order = new();
    private readonly Dictionary<ulong, LinkedListNode<(ulong Index, Memory<byte> Chunk)>> _map = [];
    private Action<ulong, Memory<byte>>? _writer;

    public Memory<byte> GetChunk(ulong index, Func<Memory<byte>> chunkAllocator, Action<ulong, Memory<byte>> chunkWriter)
    {
        if (_map.TryGetValue(index, out var node)) return node.Value.Chunk;
        while (_writer is not null && _order.Count >= Math.Max(1, capacity))
        {
            var oldest = _order.First!;
            _order.RemoveFirst();
            _map.Remove(oldest.Value.Index);
            _writer(oldest.Value.Index, oldest.Value.Chunk);
        }
        var chunk = chunkAllocator();
        _map[index] = _order.AddLast((index, chunk));
        return chunk;
    }

    public void Flush(Action<ulong, Memory<byte>> chunkWriter)
    {
        _writer ??= chunkWriter;
        foreach (var (i, c) in _order) chunkWriter(i, c);
        _order.Clear();
        _map.Clear();
    }
}
