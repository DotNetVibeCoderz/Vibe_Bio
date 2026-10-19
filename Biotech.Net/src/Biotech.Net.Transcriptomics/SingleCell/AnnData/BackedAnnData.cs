using PureHDF;
using PureHDF.Selections;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// An <c>.h5ad</c> file opened in backed mode (<c>anndata.read_h5ad(backed="r")</c>): <see cref="Obs"/>,
/// <see cref="Var"/> are loaded, X and the layers stay on disk as <see cref="H5adMatrix"/> sources that read row
/// or column slices on demand. Dispose to close the file.
/// </summary>
public sealed class BackedAnnData : IDisposable
{
    private readonly NativeFileHandle _handle;
    private Dictionary<string, object?>? _uns;

    internal BackedAnnData(string path)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        Path = path;
        _handle = new NativeFileHandle(H5File.OpenRead(path));
        var root = _handle.File;
        X = root.LinkExists("X") ? new H5adMatrix(_handle, "X") : null;
        Obs = root.LinkExists("obs") ? H5ad.ReadDataFrame(root.Get("obs"), "obs")
            : new AnnDataFrame(Enumerable.Range(0, X?.Rows ?? 0).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray());
        Var = root.LinkExists("var") ? H5ad.ReadDataFrame(root.Get("var"), "var")
            : new AnnDataFrame(Enumerable.Range(0, X?.Columns ?? 0).Select(i => i.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray());
        LayerNames = ChildNames(root, "layers");
        ObsmKeys = ChildNames(root, "obsm");
        if (X is not null && (X.Rows != Obs.RowCount || X.Columns != Var.RowCount))
            throw new InvalidDataException($"X is {X.Rows}×{X.Columns} but obs has {Obs.RowCount} rows and var {Var.RowCount}.");
    }

    private static string[] ChildNames(IH5Group root, string name)
        => root.LinkExists(name) && root.Get(name) is IH5Group g ? g.Children().Select(c => c.Name).ToArray() : [];

    /// <summary>The file path.</summary>
    public string Path { get; }

    /// <summary>The on-disk X (null when the file has none).</summary>
    public H5adMatrix? X { get; }

    /// <summary>Per-cell annotations (loaded).</summary>
    public AnnDataFrame Obs { get; }

    /// <summary>Per-gene annotations (loaded).</summary>
    public AnnDataFrame Var { get; }

    /// <summary>Number of cells.</summary>
    public int ObsCount => Obs.RowCount;

    /// <summary>Number of genes.</summary>
    public int VarCount => Var.RowCount;

    /// <summary>Names of the layers.</summary>
    public IReadOnlyList<string> LayerNames { get; }

    /// <summary>Names of the <c>obsm</c> entries.</summary>
    public IReadOnlyList<string> ObsmKeys { get; }

    /// <summary>Unstructured annotations (read on first access).</summary>
    public IReadOnlyDictionary<string, object?> Uns
    {
        get
        {
            if (_uns is null)
            {
                lock (_handle.Sync)
                {
                    var root = _handle.File;
                    var u = root.LinkExists("uns") && root.Get("uns") is IH5Group g ? H5ad.ReadElement(g, "uns") as Dictionary<string, object?> : null;
                    _uns = u ?? new Dictionary<string, object?>(StringComparer.Ordinal);
                }
            }
            return _uns;
        }
    }

    /// <summary>A layer as an on-disk matrix.</summary>
    public H5adMatrix Layer(string name)
    {
        if (!LayerNames.Contains(name)) throw new KeyNotFoundException($"No layer '{name}' in {Path}.");
        return new H5adMatrix(_handle, $"layers/{name}");
    }

    /// <summary>Reads one <c>obsm</c> entry into memory (an <see cref="AnnMatrix"/> or <see cref="AnnDataFrame"/>).</summary>
    public object ReadObsm(string key)
    {
        lock (_handle.Sync)
        {
            var o = _handle.File.Get($"obsm/{key}");
            return H5ad.Encoding(o) == "dataframe" ? H5ad.ReadDataFrame(o, $"obsm/{key}") : H5ad.ReadMatrix(o, $"obsm/{key}");
        }
    }

    /// <summary>Loads the whole file into memory (like <c>adata.to_memory()</c>).</summary>
    public AnnData ToMemory()
    {
        lock (_handle.Sync) return H5ad.ReadRoot(_handle.File);
    }

    /// <summary>Closes the file.</summary>
    public void Dispose() => _handle.Dispose();

    /// <inheritdoc />
    public override string ToString() => $"Backed AnnData {ObsCount} × {VarCount} ({X?.Format.ToString() ?? "no X"}) at {Path}";
}

/// <summary>An open PureHDF file plus the lock that serialises access to it (PureHDF readers are not thread-safe).</summary>
internal sealed class NativeFileHandle(PureHDF.VOL.Native.NativeFile file) : IDisposable
{
    public PureHDF.VOL.Native.NativeFile File { get; } = file;
    public object Sync { get; } = new();
    public bool Disposed { get; private set; }

    public void Dispose()
    {
        lock (Sync)
        {
            if (Disposed) return;
            Disposed = true;
            File.Dispose();
        }
    }
}

/// <summary>
/// X or a layer of a backed <c>.h5ad</c> file. CSR matrices read row blocks directly (the <c>indptr</c> array, 8 bytes
/// per row, is cached; <c>indices</c>/<c>data</c> are read by hyperslab, decompressing only the chunks touched).
/// CSC matrices read column blocks directly (<see cref="ReadColumns"/>); reading rows of a CSC matrix has to scan
/// every stored entry, so convert it once with <see cref="MappedSparseMatrix.Create(string, ISparseMatrixSource, int, IProgress{long}?)"/>
/// (which transposes CSC out of core) before row-wise streaming. Dense matrices read row hyperslabs.
/// </summary>
public sealed class H5adMatrix : ISparseMatrixSource, IRowPointers
{
    private readonly NativeFileHandle _handle;
    private readonly string _path;
    private readonly long[]? _indptr;

    internal H5adMatrix(NativeFileHandle handle, string path)
    {
        _handle = handle; _path = path;
        lock (handle.Sync)
        {
            var o = handle.File.Get(path);
            if (o is IH5Dataset ds)
            {
                Format = AnnMatrixFormat.Dense;
                var dims = ds.Space.Dimensions;
                if (dims.Length != 2) throw new InvalidDataException($"'{path}' is not a 2-D array.");
                Rows = (int)dims[0]; Columns = (int)dims[1];
                DataType = H5Types.DTypeName(ds.Type) ?? throw new InvalidDataException($"'{path}' is not numeric.");
                NonZeroCount = (long)Rows * Columns;
            }
            else
            {
                var g = (IH5Group)o;
                var (format, shape) = H5ad.SparseHeader(g, path);
                Format = format; Rows = shape.Rows; Columns = shape.Columns;
                var data = (IH5Dataset)g.Get("data");
                DataType = H5Types.DTypeName(data.Type) ?? "float32";
                int major = format == AnnMatrixFormat.Csr ? Rows : Columns;
                _indptr = H5Types.ReadInt64((IH5Dataset)g.Get("indptr"), 0, (ulong)major + 1);
                NonZeroCount = _indptr[^1];
            }
        }
    }

    long IRowPointers.RowPointer(int row) => Format == AnnMatrixFormat.Csr ? _indptr![row] : (long)row * (Format == AnnMatrixFormat.Dense ? Columns : 0);

    /// <summary>On-disk layout.</summary>
    public AnnMatrixFormat Format { get; }

    /// <summary>On-disk element type.</summary>
    public string DataType { get; }

    /// <inheritdoc />
    public int Rows { get; }

    /// <inheritdoc />
    public int Columns { get; }

    /// <inheritdoc />
    public long NonZeroCount { get; }

    /// <summary>Bytes of entries one scan step reads at a time when scanning (CSC rows, CSR columns).</summary>
    internal const int ScanEntries = 1 << 22;

    private T WithFile<T>(Func<IH5Group, T> f)
    {
        lock (_handle.Sync)
        {
            ObjectDisposedException.ThrowIf(_handle.Disposed, this);
            return f(_handle.File);
        }
    }

    /// <inheritdoc />
    public SparseMatrix ReadRows(int start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        if (start + count > Rows) throw new ArgumentOutOfRangeException(nameof(count));
        return Format switch
        {
            AnnMatrixFormat.Csr => ReadMajor(start, count, Columns),
            AnnMatrixFormat.Dense => ReadDenseRows(start, count),
            _ => ScanMinor(start, count, Rows),
        };
    }

    /// <summary>
    /// Columns <paramref name="start"/> … <paramref name="start"/> + <paramref name="count"/> − 1 as a CSR matrix of
    /// all rows × <paramref name="count"/> columns. Direct for CSC; CSR and dense matrices are scanned in full.
    /// </summary>
    public SparseMatrix ReadColumns(int start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        if (start + count > Columns) throw new ArgumentOutOfRangeException(nameof(count));
        if (Format == AnnMatrixFormat.Csc) return ReadMajor(start, count, Rows).Transpose();
        if (Format == AnnMatrixFormat.Csr) return ScanMinor(start, count, Columns).Transpose();
        var blocks = new List<SparseMatrix>();
        for (int r = 0; r < Rows; r += SparseMatrixSource.DefaultBlockRows)
        {
            var b = ReadDenseRows(r, Math.Min(SparseMatrixSource.DefaultBlockRows, Rows - r));
            blocks.Add(b.SelectColumns(Enumerable.Range(start, count).ToArray()));
        }
        return SparseMatrixSource.VStack(blocks, count);
    }

    /// <summary>Reads major-axis slices (rows of CSR, columns of CSC) as CSR over the minor axis.</summary>
    private SparseMatrix ReadMajor(int start, int count, int minorLength)
    {
        long s = _indptr![start], e = _indptr[start + count];
        long n = e - s;
        if (n > int.MaxValue) throw new InvalidOperationException("Slice holds more than 2^31 entries; read fewer rows at a time.");
        var ptr = new int[count + 1];
        for (int i = 0; i <= count; i++) ptr[i] = (int)(_indptr[start + i] - s);
        if (n == 0) return new SparseMatrix(count, minorLength, ptr, [], []);
        var (idx, vals) = WithFile(root =>
        {
            var g = (IH5Group)root.Get(_path);
            var idx = H5Types.ReadInt32((IH5Dataset)g.Get("indices"), (ulong)s, (ulong)n);
            var vals = H5Types.ReadSingle((IH5Dataset)g.Get("data"), new HyperslabSelection((ulong)s, (ulong)n));
            return (idx, vals);
        });
        H5ad.SortIndices(ptr, idx, vals);
        return new SparseMatrix(count, minorLength, ptr, idx, vals);
    }

    /// <summary>
    /// Entries whose minor index lies in [start, start + count): scans the whole matrix in chunks and returns
    /// CSR of <paramref name="count"/> minor slots × major positions.
    /// </summary>
    private SparseMatrix ScanMinor(int start, int count, int minorLength)
    {
        int major = Format == AnnMatrixFormat.Csr ? Rows : Columns;
        var rowIdx = new List<int>(); var colIdx = new List<int>(); var vals = new List<float>();
        int m = 0; // current major position
        for (long pos = 0; pos < NonZeroCount; pos += ScanEntries)
        {
            long n = Math.Min(ScanEntries, NonZeroCount - pos);
            var (idx, v) = WithFile(root =>
            {
                var g = (IH5Group)root.Get(_path);
                return (H5Types.ReadInt32((IH5Dataset)g.Get("indices"), (ulong)pos, (ulong)n),
                        H5Types.ReadSingle((IH5Dataset)g.Get("data"), new HyperslabSelection((ulong)pos, (ulong)n)));
            });
            for (int k = 0; k < idx.Length; k++)
            {
                long at = pos + k;
                while (_indptr![m + 1] <= at) m++;
                int minor = idx[k];
                if (minor >= start && minor < start + count) { rowIdx.Add(minor - start); colIdx.Add(m); vals.Add(v[k]); }
            }
        }
        return SparseMatrix.FromTriplets(count, major, rowIdx.ToArray(), colIdx.ToArray(), vals.ToArray());
    }

    private SparseMatrix ReadDenseRows(int start, int count)
    {
        if (count == 0) return new SparseMatrix(0, Columns, [0], [], []);
        var flat = WithFile(root =>
        {
            var ds = (IH5Dataset)root.Get(_path);
            return H5Types.ReadSingle(ds, new HyperslabSelection(2, [(ulong)start, 0], [(ulong)count, (ulong)Columns]));
        });
        var ptr = new int[count + 1];
        var cols = new List<int>(); var vals = new List<float>();
        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < Columns; j++)
            {
                float v = flat[(long)i * Columns + j];
                if (v != 0) { cols.Add(j); vals.Add(v); }
            }
            ptr[i + 1] = cols.Count;
        }
        return new SparseMatrix(count, Columns, ptr, cols.ToArray(), vals.ToArray());
    }

    /// <summary>The CSC column pointers (null unless CSC); used by the out-of-core CSC → CSR transpose.</summary>
    internal long[]? MajorPointers => _indptr;

    /// <summary>Reads a raw slice of indices and data in storage order (sparse only).</summary>
    internal (int[] Indices, float[] Values) ReadEntries(long pos, int count) => WithFile(root =>
    {
        var g = (IH5Group)root.Get(_path);
        return (H5Types.ReadInt32((IH5Dataset)g.Get("indices"), (ulong)pos, (ulong)count),
                H5Types.ReadSingle((IH5Dataset)g.Get("data"), new HyperslabSelection((ulong)pos, (ulong)count)));
    });

    /// <inheritdoc />
    public override string ToString() => $"H5adMatrix '{_path}' {Rows}×{Columns} {Format} {DataType}, {NonZeroCount} stored";
}
