namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// A cells × genes matrix that can be read in blocks of consecutive rows without holding the whole matrix in
/// memory: an in-memory <see cref="SparseMatrix"/>, the X (or a layer) of an <c>.h5ad</c> file opened with
/// <see cref="H5ad.OpenBacked"/>, or a memory-mapped <see cref="MappedSparseMatrix"/>. The out-of-core algorithms in
/// <see cref="OutOfCore"/> and <see cref="OutOfCorePca"/> stream over this interface.
/// </summary>
public interface ISparseMatrixSource
{
    /// <summary>Number of rows (cells).</summary>
    int Rows { get; }

    /// <summary>Number of columns (genes).</summary>
    int Columns { get; }

    /// <summary>Number of stored entries (may exceed <see cref="int.MaxValue"/>).</summary>
    long NonZeroCount { get; }

    /// <summary>
    /// Rows <paramref name="start"/> … <paramref name="start"/> + <paramref name="count"/> − 1 as an in-memory CSR
    /// block (row 0 of the block is row <paramref name="start"/> of the source). Column indices are sorted per row.
    /// The block is a copy: callers may modify it.
    /// </summary>
    SparseMatrix ReadRows(int start, int count);
}

/// <summary>Sources that can report row start offsets cheaply (used to size row blocks by entry count).</summary>
internal interface IRowPointers
{
    long RowPointer(int row);
}

/// <summary>Helpers over <see cref="ISparseMatrixSource"/>.</summary>
public static class SparseMatrixSource
{
    /// <summary>Default maximum rows per block.</summary>
    public const int DefaultBlockRows = 8192;

    /// <summary>
    /// Default maximum stored entries per block (2²¹ ≈ 2 M entries, 16 MB of indices + values). Sources that know
    /// their row pointers (<see cref="SparseMatrix"/>, <see cref="MappedSparseMatrix"/>, CSR <see cref="H5adMatrix"/>)
    /// shorten blocks to stay under it, which bounds memory however dense the rows are.
    /// </summary>
    public const int DefaultBlockEntries = 1 << 21;

    /// <summary>
    /// Streams the source as consecutive row blocks of at most <paramref name="blockRows"/> rows and (when the
    /// source exposes row pointers) about <paramref name="maxEntries"/> stored entries; a single row is never split.
    /// Only one block is alive at a time unless the caller keeps references.
    /// </summary>
    public static IEnumerable<(int Start, SparseMatrix Block)> RowBlocks(this ISparseMatrixSource source, int blockRows = DefaultBlockRows, int maxEntries = DefaultBlockEntries)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentOutOfRangeException.ThrowIfLessThan(blockRows, 1);
        var rp = source as IRowPointers;
        for (int start = 0; start < source.Rows;)
        {
            int count = Math.Min(blockRows, source.Rows - start);
            if (rp is not null && count > 1)
            {
                long limit = rp.RowPointer(start) + Math.Max(1, maxEntries);
                if (rp.RowPointer(start + count) > limit)
                {
                    // Largest count with RowPointer(start + count) <= limit (at least one row).
                    int lo = 1, hi = count;
                    while (lo < hi) { int mid = (lo + hi + 1) >> 1; if (rp.RowPointer(start + mid) <= limit) lo = mid; else hi = mid - 1; }
                    count = lo;
                }
            }
            yield return (start, source.ReadRows(start, count));
            start += count;
        }
    }

    /// <summary>
    /// Runs <paramref name="body"/> on every row block, reading blocks sequentially (I/O stays in file order) while
    /// up to <paramref name="maxParallelism"/> − 1 previously read blocks are processed on worker threads.
    /// </summary>
    internal static void ForEachBlock(ISparseMatrixSource source, int blockRows, Action<int, SparseMatrix> body, int? maxParallelism = null)
    {
        int workers = Math.Max(1, maxParallelism ?? Environment.ProcessorCount);
        if (workers == 1)
        {
            foreach (var (s, b) in source.RowBlocks(blockRows)) body(s, b);
            return;
        }
        using var gate = new SemaphoreSlim(workers);
        var tasks = new List<Task>();
        Exception? failure = null;
        foreach (var (s, b) in source.RowBlocks(blockRows))
        {
            gate.Wait();
            if (Volatile.Read(ref failure) is not null) { gate.Release(); break; }
            int start = s; var block = b;
            tasks.Add(Task.Run(() =>
            {
                try { body(start, block); }
                catch (Exception e) { Interlocked.CompareExchange(ref failure, e, null); }
                finally { gate.Release(); }
            }));
            tasks.RemoveAll(t => t.IsCompleted);
        }
        Task.WaitAll(tasks.ToArray());
        if (failure is not null) throw new AggregateException("An out-of-core block worker failed.", failure);
    }

    /// <summary>Loads the whole source into memory (only for sources that fit).</summary>
    public static SparseMatrix ToSparseMatrix(this ISparseMatrixSource source)
    {
        ArgumentNullException.ThrowIfNull(source);
        if (source is SparseMatrix m) return m.Clone();
        if (source.NonZeroCount > int.MaxValue) throw new InvalidOperationException("The matrix has more than 2^31 stored entries; it cannot be loaded into one SparseMatrix.");
        return source.ReadRows(0, source.Rows);
    }

    /// <summary>Concatenates row blocks (same column count) into one matrix.</summary>
    internal static SparseMatrix VStack(IReadOnlyList<SparseMatrix> blocks, int columns)
    {
        int rows = blocks.Sum(b => b.Rows);
        long nnz = blocks.Sum(b => (long)b.NonZeroCount);
        if (nnz > int.MaxValue) throw new InvalidOperationException("Too many stored entries for one SparseMatrix.");
        var ptr = new int[rows + 1]; var cols = new int[nnz]; var vals = new float[nnz];
        int r = 0, at = 0;
        foreach (var b in blocks)
        {
            b.ColumnIndices.AsSpan(0, b.NonZeroCount).CopyTo(cols.AsSpan(at));
            b.Values.AsSpan(0, b.NonZeroCount).CopyTo(vals.AsSpan(at));
            for (int i = 0; i < b.Rows; i++) ptr[r + i + 1] = at + b.RowPointers[i + 1];
            r += b.Rows; at += b.NonZeroCount;
        }
        return new SparseMatrix(rows, columns, ptr, cols, vals);
    }
}
