using System.IO.MemoryMappedFiles;
using System.Runtime.InteropServices;
using System.Text.Json;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// An on-disk CSR matrix in the Biotech.Net <c>.bnsparse</c> directory format, read through memory-mapped files:
/// <c>matrix.json</c> (shape, entry count, layer names), <c>indptr.i64</c> (rows + 1 little-endian int64),
/// <c>indices.i32</c> (int32 column indices, sorted per row) and one <c>&lt;layer&gt;.f32</c> float32 value file per
/// layer. Layers share the sparsity structure, so a normalised or log-transformed layer costs 4 bytes per stored
/// entry. Use it for matrices that do not fit in memory; every out-of-core algorithm accepts it as an
/// <see cref="ISparseMatrixSource"/>.
/// </summary>
/// <remarks>
/// Reads copy one row block at a time out of the mapping, so the managed heap stays bounded by the block size.
/// The operating system caches the mapped pages it touches; they count towards the process working set but are
/// clean file-backed pages the OS reclaims under memory pressure (private memory stays bounded).
/// </remarks>
public sealed unsafe class MappedSparseMatrix : ISparseMatrixSource, IRowPointers, IDisposable
{
    /// <summary>Name of the default layer (the matrix written by <see cref="MappedSparseMatrixWriter"/>).</summary>
    public const string DefaultLayer = "X";

    private readonly Mapping _indptr, _indices, _values;

    private MappedSparseMatrix(string directory, Header header, string layer)
    {
        Directory = directory;
        Rows = header.Rows; Columns = header.Columns; NonZeroCount = header.NonZeros;
        LayerNames = header.Layers;
        Layer = layer;
        _indptr = new Mapping(System.IO.Path.Combine(directory, "indptr.i64"), (Rows + 1L) * 8);
        _indices = new Mapping(System.IO.Path.Combine(directory, "indices.i32"), NonZeroCount * 4);
        _values = new Mapping(LayerPath(directory, layer), NonZeroCount * 4);
    }

    /// <summary>Opens a store (the given layer, <c>X</c> by default).</summary>
    public static MappedSparseMatrix Open(string directory, string layer = DefaultLayer)
    {
        ArgumentException.ThrowIfNullOrEmpty(directory);
        var h = Header.Load(directory);
        if (!h.Layers.Contains(layer)) throw new KeyNotFoundException($"Layer '{layer}' not found in {directory} (layers: {string.Join(", ", h.Layers)}).");
        return new MappedSparseMatrix(System.IO.Path.GetFullPath(directory), h, layer);
    }

    /// <summary>
    /// Writes <paramref name="source"/> to a new store at <paramref name="directory"/> and opens it. CSC
    /// <see cref="H5adMatrix"/> sources are transposed out of core (two passes over the entries; memory 12 bytes per
    /// row), other sources are streamed in row blocks.
    /// </summary>
    public static MappedSparseMatrix Create(string directory, ISparseMatrixSource source, int blockRows = SparseMatrixSource.DefaultBlockRows, IProgress<long>? progress = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(directory);
        ArgumentNullException.ThrowIfNull(source);
        if (source is H5adMatrix { Format: AnnMatrixFormat.Csc } csc) return TransposeCsc(directory, csc, progress);
        using var w = MappedSparseMatrixWriter.Create(directory, source.Columns);
        long done = 0;
        foreach (var (_, block) in source.RowBlocks(blockRows))
        {
            w.AppendRows(block);
            done += block.Rows;
            progress?.Report(done);
        }
        return w.Complete();
    }

    private static MappedSparseMatrix TransposeCsc(string directory, H5adMatrix csc, IProgress<long>? progress)
    {
        int rows = csc.Rows, cols = csc.Columns;
        long nnz = csc.NonZeroCount;
        var colPtr = csc.MajorPointers!;
        // Pass 1: entries per row.
        var counts = new int[rows];
        for (long pos = 0; pos < nnz; pos += H5adMatrix.ScanEntries)
        {
            int n = (int)Math.Min(H5adMatrix.ScanEntries, nnz - pos);
            var (idx, _) = csc.ReadEntries(pos, n);
            foreach (int r in idx) counts[r]++;
        }
        var ptr = new long[rows + 1];
        for (int r = 0; r < rows; r++) ptr[r + 1] = ptr[r] + counts[r];
        counts = null!;
        System.IO.Directory.CreateDirectory(directory);
        WriteArray(System.IO.Path.Combine(directory, "indptr.i64"), ptr);
        var next = (long[])ptr.Clone();
        // Pass 2: scatter into writable mappings; columns arrive in ascending order, so each row's entries are sorted.
        using (var outIdx = new Mapping(System.IO.Path.Combine(directory, "indices.i32"), nnz * 4, writable: true))
        using (var outVal = new Mapping(LayerPath(directory, DefaultLayer), nnz * 4, writable: true))
        {
            int* pi = (int*)outIdx.Pointer; float* pv = (float*)outVal.Pointer;
            int c = 0;
            for (long pos = 0; pos < nnz; pos += H5adMatrix.ScanEntries)
            {
                int n = (int)Math.Min(H5adMatrix.ScanEntries, nnz - pos);
                var (idx, vals) = csc.ReadEntries(pos, n);
                for (int k = 0; k < n; k++)
                {
                    long at = pos + k;
                    while (colPtr[c + 1] <= at) c++;
                    long dst = next[idx[k]]++;
                    pi[dst] = c; pv[dst] = vals[k];
                }
                progress?.Report(pos + n);
            }
        }
        new Header(rows, cols, nnz, [DefaultLayer]).Save(directory);
        return Open(directory);
    }

    /// <summary>The store directory.</summary>
    public string Directory { get; }

    /// <summary>The layer this instance reads.</summary>
    public string Layer { get; }

    /// <summary>All layers in the store.</summary>
    public IReadOnlyList<string> LayerNames { get; private set; }

    /// <inheritdoc />
    public int Rows { get; }

    /// <inheritdoc />
    public int Columns { get; }

    /// <inheritdoc />
    public long NonZeroCount { get; }

    /// <summary>Opens another layer of the same store (shares the structure files).</summary>
    public MappedSparseMatrix OpenLayer(string layer) => Open(Directory, layer);

    /// <summary>Row start offset (0 ≤ row ≤ Rows).</summary>
    public long RowPointer(int row)
    {
        if ((uint)row > (uint)Rows) throw new ArgumentOutOfRangeException(nameof(row));
        return ((long*)_indptr.Pointer)[row];
    }

    /// <summary>Zero-copy view of one row's column indices (valid while this instance is open).</summary>
    public ReadOnlySpan<int> RowIndices(int row)
    {
        long s = RowPointer(row), e = RowPointer(row + 1);
        return new ReadOnlySpan<int>((int*)_indices.Pointer + s, (int)(e - s));
    }

    /// <summary>Zero-copy view of one row's values (valid while this instance is open).</summary>
    public ReadOnlySpan<float> RowValues(int row)
    {
        long s = RowPointer(row), e = RowPointer(row + 1);
        return new ReadOnlySpan<float>((float*)_values.Pointer + s, (int)(e - s));
    }

    /// <inheritdoc />
    public SparseMatrix ReadRows(int start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        if (start + count > Rows) throw new ArgumentOutOfRangeException(nameof(count));
        long* p = (long*)_indptr.Pointer;
        long s = p[start], e = p[start + count];
        if (e - s > int.MaxValue) throw new InvalidOperationException("Block holds more than 2^31 entries; read fewer rows at a time.");
        int n = (int)(e - s);
        var ptr = new int[count + 1];
        for (int i = 0; i <= count; i++) ptr[i] = (int)(p[start + i] - s);
        var idx = GC.AllocateUninitializedArray<int>(n);
        var vals = GC.AllocateUninitializedArray<float>(n);
        new ReadOnlySpan<int>((int*)_indices.Pointer + s, n).CopyTo(idx);
        new ReadOnlySpan<float>((float*)_values.Pointer + s, n).CopyTo(vals);
        return new SparseMatrix(count, Columns, ptr, idx, vals);
    }

    /// <summary>
    /// Writes a new layer whose values are this layer's values transformed block by block:
    /// <paramref name="transform"/> receives the first row of each block and the block, and modifies
    /// <see cref="SparseMatrix.Values"/> in place (the sparsity structure must not change). Returns the new layer.
    /// </summary>
    public MappedSparseMatrix WriteLayer(string name, Action<int, SparseMatrix> transform, int blockRows = SparseMatrixSource.DefaultBlockRows)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentNullException.ThrowIfNull(transform);
        if (name.IndexOfAny(System.IO.Path.GetInvalidFileNameChars()) >= 0) throw new ArgumentException("Invalid layer name.", nameof(name));
        if (name == Layer) throw new ArgumentException("Cannot overwrite the layer being read.", nameof(name));
        string tmp = LayerPath(Directory, name) + ".tmp";
        using (var fs = new FileStream(tmp, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
        {
            foreach (var (start, block) in this.RowBlocks(blockRows))
            {
                transform(start, block);
                if (block.NonZeroCount != RowPointer(start + block.Rows) - RowPointer(start))
                    throw new InvalidOperationException("The transform changed the sparsity structure.");
                fs.Write(MemoryMarshal.AsBytes(block.Values.AsSpan(0, block.NonZeroCount)));
            }
        }
        File.Move(tmp, LayerPath(Directory, name), overwrite: true);
        var h = Header.Load(Directory);
        if (!h.Layers.Contains(name)) h = h with { Layers = [.. h.Layers, name] };
        h.Save(Directory);
        LayerNames = h.Layers;
        return OpenLayer(name);
    }

    /// <summary>Closes the mappings.</summary>
    public void Dispose()
    {
        _indptr.Dispose(); _indices.Dispose(); _values.Dispose();
    }

    /// <inheritdoc />
    public override string ToString() => $"MappedSparseMatrix {Rows}×{Columns}, {NonZeroCount} stored, layer '{Layer}' at {Directory}";

    internal static string LayerPath(string directory, string layer) => System.IO.Path.Combine(directory, layer + ".f32");

    internal static void WriteArray<T>(string path, T[] data) where T : unmanaged
    {
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        fs.Write(MemoryMarshal.AsBytes(data.AsSpan()));
    }

    internal sealed record Header(int Rows, int Columns, long NonZeros, string[] Layers)
    {
        private static readonly JsonSerializerOptions Json = new() { WriteIndented = true };

        public void Save(string directory)
        {
            var doc = new Dictionary<string, object>
            {
                ["format"] = "bnsparse", ["version"] = 1, ["layout"] = "csr", ["rows"] = Rows, ["columns"] = Columns,
                ["nnz"] = NonZeros, ["indptr"] = "int64", ["indices"] = "int32", ["values"] = "float32", ["layers"] = Layers,
            };
            File.WriteAllText(System.IO.Path.Combine(directory, "matrix.json"), JsonSerializer.Serialize(doc, Json));
        }

        public static Header Load(string directory)
        {
            var p = System.IO.Path.Combine(directory, "matrix.json");
            if (!File.Exists(p)) throw new FileNotFoundException($"{directory} is not a .bnsparse store (matrix.json missing).", p);
            using var d = JsonDocument.Parse(File.ReadAllText(p));
            var r = d.RootElement;
            if (r.GetProperty("format").GetString() != "bnsparse") throw new InvalidDataException("Not a bnsparse matrix.json.");
            return new Header(r.GetProperty("rows").GetInt32(), r.GetProperty("columns").GetInt32(), r.GetProperty("nnz").GetInt64(),
                r.GetProperty("layers").EnumerateArray().Select(e => e.GetString()!).ToArray());
        }
    }

    /// <summary>A whole-file mapping with a raw pointer (null for empty files).</summary>
    internal sealed class Mapping : IDisposable
    {
        private readonly MemoryMappedFile? _file;
        private readonly MemoryMappedViewAccessor? _view;
        private bool _acquired;

        public Mapping(string path, long expectedBytes, bool writable = false)
        {
            if (writable)
            {
                using (var fs = new FileStream(path, FileMode.Create, FileAccess.ReadWrite)) fs.SetLength(expectedBytes);
            }
            else
            {
                if (!File.Exists(path)) throw new FileNotFoundException($"Missing {path}.", path);
                long len = new FileInfo(path).Length;
                if (len < expectedBytes) throw new InvalidDataException($"{path} has {len} bytes, expected {expectedBytes}.");
            }
            if (expectedBytes == 0) return;
            var access = writable ? MemoryMappedFileAccess.ReadWrite : MemoryMappedFileAccess.Read;
            _file = MemoryMappedFile.CreateFromFile(path, FileMode.Open, null, 0, access);
            _view = _file.CreateViewAccessor(0, expectedBytes, access);
            byte* p = null;
            _view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            _acquired = true;
            Pointer = p + _view.PointerOffset;
        }

        public byte* Pointer { get; }

        public void Dispose()
        {
            if (_acquired) { _view!.SafeMemoryMappedViewHandle.ReleasePointer(); _acquired = false; }
            _view?.Dispose();
            _file?.Dispose();
        }
    }
}

/// <summary>
/// Streams rows into a new <see cref="MappedSparseMatrix"/> store: append CSR row blocks in order, then
/// <see cref="Complete"/>. Indices and values go straight to disk; only the row pointers (8 bytes per row) are kept
/// in memory.
/// </summary>
public sealed class MappedSparseMatrixWriter : IDisposable
{
    private readonly string _directory;
    private readonly FileStream _indices, _values;
    private readonly List<long> _ptr = [0];
    private bool _done;

    private MappedSparseMatrixWriter(string directory, int columns)
    {
        _directory = System.IO.Path.GetFullPath(directory);
        Columns = columns;
        System.IO.Directory.CreateDirectory(_directory);
        _indices = new FileStream(System.IO.Path.Combine(_directory, "indices.i32"), FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        _values = new FileStream(MappedSparseMatrix.LayerPath(_directory, MappedSparseMatrix.DefaultLayer), FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
    }

    /// <summary>Creates (or overwrites) a store in <paramref name="directory"/>.</summary>
    public static MappedSparseMatrixWriter Create(string directory, int columns)
    {
        ArgumentException.ThrowIfNullOrEmpty(directory);
        ArgumentOutOfRangeException.ThrowIfNegative(columns);
        return new MappedSparseMatrixWriter(directory, columns);
    }

    /// <summary>Number of columns.</summary>
    public int Columns { get; }

    /// <summary>Rows written so far.</summary>
    public int Rows => _ptr.Count - 1;

    /// <summary>Entries written so far.</summary>
    public long NonZeroCount => _ptr[^1];

    /// <summary>Appends the rows of <paramref name="block"/> (column indices must be sorted per row).</summary>
    public void AppendRows(SparseMatrix block)
    {
        ArgumentNullException.ThrowIfNull(block);
        ObjectDisposedException.ThrowIf(_done, this);
        if (block.Columns != Columns) throw new ArgumentException($"Block has {block.Columns} columns, the store {Columns}.");
        int n = block.NonZeroCount;
        _indices.Write(MemoryMarshal.AsBytes(block.ColumnIndices.AsSpan(0, n)));
        _values.Write(MemoryMarshal.AsBytes(block.Values.AsSpan(0, n)));
        long baseOffset = _ptr[^1];
        for (int i = 1; i <= block.Rows; i++) _ptr.Add(baseOffset + block.RowPointers[i]);
    }

    /// <summary>Appends one row.</summary>
    public void AppendRow(ReadOnlySpan<int> columns, ReadOnlySpan<float> values)
    {
        ObjectDisposedException.ThrowIf(_done, this);
        if (columns.Length != values.Length) throw new ArgumentException("Columns and values differ in length.");
        _indices.Write(MemoryMarshal.AsBytes(columns));
        _values.Write(MemoryMarshal.AsBytes(values));
        _ptr.Add(_ptr[^1] + columns.Length);
    }

    /// <summary>Finishes the store and opens it for reading.</summary>
    public MappedSparseMatrix Complete()
    {
        ObjectDisposedException.ThrowIf(_done, this);
        _done = true;
        _indices.Dispose(); _values.Dispose();
        MappedSparseMatrix.WriteArray(System.IO.Path.Combine(_directory, "indptr.i64"), _ptr.ToArray());
        new MappedSparseMatrix.Header(Rows, Columns, NonZeroCount, [MappedSparseMatrix.DefaultLayer]).Save(_directory);
        return MappedSparseMatrix.Open(_directory);
    }

    /// <summary>Closes the files (an incomplete store is left without <c>matrix.json</c>).</summary>
    public void Dispose()
    {
        if (_done) return;
        _done = true;
        _indices.Dispose(); _values.Dispose();
    }
}
