using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Search;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Gpu.Alignment;

/// <summary>One entry of a top-k list: a target index and its local score.</summary>
/// <param name="Target">Index into the target list.</param>
/// <param name="Score">Smith–Waterman score.</param>
public readonly record struct LocalScoreHit(int Target, int Score);

/// <summary>A dense query × target matrix of local alignment scores.</summary>
public sealed class LocalScoreMatrix
{
    private readonly int[] _scores;

    internal LocalScoreMatrix(int queries, int targets, int[] scores, string backend)
    {
        QueryCount = queries;
        TargetCount = targets;
        _scores = scores;
        Backend = backend;
    }

    /// <summary>Number of rows (queries).</summary>
    public int QueryCount { get; }

    /// <summary>Number of columns (targets).</summary>
    public int TargetCount { get; }

    /// <summary>Where the scores were computed, e.g. <c>OpenCL: Intel(R) UHD Graphics 620</c> or <c>CPU (Native)</c>.</summary>
    public string Backend { get; }

    /// <summary>Score of query <paramref name="query"/> against target <paramref name="target"/>.</summary>
    public int this[int query, int target]
    {
        get
        {
            ArgumentOutOfRangeException.ThrowIfNegative(query);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(query, QueryCount);
            ArgumentOutOfRangeException.ThrowIfNegative(target);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(target, TargetCount);
            return _scores[query * TargetCount + target];
        }
    }

    /// <summary>The scores of one query, in target order.</summary>
    public ReadOnlySpan<int> Row(int query) => _scores.AsSpan(query * TargetCount, TargetCount);

    /// <summary>All scores, row-major (query-major).</summary>
    public ReadOnlySpan<int> Values => _scores;

    /// <summary>A copy as a 2-D array.</summary>
    public int[,] ToArray()
    {
        var a = new int[QueryCount, TargetCount];
        for (int q = 0; q < QueryCount; q++)
            for (int t = 0; t < TargetCount; t++) a[q, t] = _scores[q * TargetCount + t];
        return a;
    }

    /// <summary>The <paramref name="k"/> best targets of a query, best first (ties: lower target index first).</summary>
    public IReadOnlyList<LocalScoreHit> Top(int query, int k)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(k);
        var row = Row(query);
        var heap = new PriorityQueue<int, (int Score, int NegIndex)>();
        for (int t = 0; t < row.Length; t++)
        {
            var key = (row[t], -t);
            if (heap.Count < k) heap.Enqueue(t, key);
            else if (k > 0 && heap.TryPeek(out _, out var min) && key.CompareTo(min) > 0)
                heap.EnqueueDequeue(t, key);
        }
        var hits = new List<LocalScoreHit>(heap.Count);
        while (heap.TryDequeue(out int t, out _)) hits.Add(new LocalScoreHit(t, row[t]));
        hits.Reverse();
        return hits;
    }

    /// <summary>Top-k lists for every query.</summary>
    public IReadOnlyList<IReadOnlyList<LocalScoreHit>> TopAll(int k)
        => Enumerable.Range(0, QueryCount).Select(q => Top(q, k)).ToList();
}

/// <summary>
/// Many-against-many Smith–Waterman (Gotoh affine gaps) scoring on a GPU: one thread per (query, target) pair over
/// length-sorted, packed sequences. Scores are bit-identical to <see cref="BioKernels.LocalScore"/> and to
/// <see cref="PairwiseAligner.Score(Sequence, Sequence)"/> in local mode (integer arithmetic, same recurrence). With a null
/// device, or if a kernel fails, the Rust/managed CPU kernel runs instead (see <see cref="LastFallbackError"/>).
/// </summary>
/// <remarks>
/// Gap costs are positive (BLAST convention): a gap of length k costs <c>GapOpen + (k − 1)·GapExtend</c>. The query
/// indexes the matrix rows, as in <see cref="SequenceDatabase"/>; for the symmetric built-in matrices orientation does not
/// matter. Top-k selection runs on the host from the full score matrix.
/// </remarks>
public sealed class GpuLocalAligner : ILocalScoreBackend
{
    private readonly int[] _matrix;

    /// <summary>Scorer with a matrix (default BLOSUM62) and BLAST gap costs (default 11/1).</summary>
    public GpuLocalAligner(GpuDevice? device, SubstitutionMatrix? matrix = null, int gapOpen = 11, int gapExtend = 1)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(gapOpen);
        ArgumentOutOfRangeException.ThrowIfNegative(gapExtend);
        Device = device;
        Matrix = matrix ?? SubstitutionMatrix.Blosum62;
        _matrix = Matrix.KernelMatrix.ToArray();
        GapOpen = gapOpen;
        GapExtend = gapExtend;
    }

    /// <summary>Takes the matrix (or match/mismatch scores) and gap scores of a local-mode <see cref="PairwiseAligner"/>.</summary>
    public static GpuLocalAligner FromPairwiseAligner(PairwiseAligner aligner, GpuDevice? device)
    {
        ArgumentNullException.ThrowIfNull(aligner);
        if (aligner.Mode != AlignmentMode.Local) throw new ArgumentException("The GPU scorer computes local (Smith–Waterman) scores; set Mode = Local.", nameof(aligner));
        if (aligner.Band >= 0) throw new ArgumentException("Banded scoring is not supported on the GPU.", nameof(aligner));
        var m = aligner.SubstitutionMatrix ?? SubstitutionMatrix.MatchMismatch(aligner.MatchScore, aligner.MismatchScore);
        return new GpuLocalAligner(device, m, -aligner.OpenGapScore, -aligner.ExtendGapScore);
    }

    /// <summary>The device, or null for the CPU path.</summary>
    public GpuDevice? Device { get; }

    /// <summary>Substitution matrix.</summary>
    public SubstitutionMatrix Matrix { get; }

    /// <summary>Gap open cost (positive).</summary>
    public int GapOpen { get; }

    /// <summary>Gap extension cost (positive).</summary>
    public int GapExtend { get; }

    /// <summary>Backend of the most recent call.</summary>
    public string LastBackend { get; private set; } = "";

    /// <summary>The exception that made the most recent call fall back to the CPU, or null.</summary>
    public Exception? LastFallbackError { get; private set; }

    /// <summary>Test hook: makes the GPU path throw so the fallback can be exercised.</summary>
    internal bool InjectFailure { get; set; }

    /// <inheritdoc/>
    public string Name => Device is { } d ? $"{d.Backend}: {d.Name}" : $"CPU ({BioKernels.Backend})";

    /// <summary>Scores every query against every target.</summary>
    public LocalScoreMatrix ScoreAll(IReadOnlyList<Sequence> queries, IReadOnlyList<Sequence> targets)
    {
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentNullException.ThrowIfNull(targets);
        var (qb, qo) = Pack(queries);
        var (tb, to) = Pack(targets);
        var scores = new int[checked(queries.Count * targets.Count)];
        string backend = Score(qb, qo, tb, to, _matrix, GapOpen, GapExtend, scores);
        return new LocalScoreMatrix(queries.Count, targets.Count, scores, backend);
    }

    /// <summary>All-vs-all within one set (the diagonal holds self scores).</summary>
    public LocalScoreMatrix ScoreAll(IReadOnlyList<Sequence> sequences) => ScoreAll(sequences, sequences);

    /// <summary>The <paramref name="k"/> best targets of each query.</summary>
    public IReadOnlyList<IReadOnlyList<LocalScoreHit>> TopHits(IReadOnlyList<Sequence> queries, IReadOnlyList<Sequence> targets, int k)
        => ScoreAll(queries, targets).TopAll(k);

    /// <summary>
    /// Raw entry point: queries and targets packed back to back with count + 1 offsets, a 128×128 kernel matrix
    /// (<see cref="SubstitutionMatrix.KernelMatrix"/>) and positive gap costs; <c>scores[q·targetCount + t]</c>.
    /// Returns the backend used.
    /// </summary>
    public string Score(ReadOnlySpan<byte> queries, ReadOnlySpan<int> queryOffsets, ReadOnlySpan<byte> targets, ReadOnlySpan<int> targetOffsets,
        ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
    {
        LocalScoreRunner.Validate(queries, queryOffsets, targets, targetOffsets, matrix, gapOpen, gapExtend, scores);
        LastFallbackError = null;
        if (Device is { } device)
        {
            try
            {
                if (InjectFailure) throw new InvalidOperationException("Injected kernel failure (test hook).");
                LocalScoreRunner.RunGpu(device, queries, queryOffsets, targets, targetOffsets, matrix, gapOpen, gapExtend, scores);
                return LastBackend = $"{device.Backend}: {device.Name}";
            }
            catch (Exception ex) when (ex is not ObjectDisposedException)
            {
                LastFallbackError = ex;
            }
        }
        LocalScoreRunner.RunCpu(queries, queryOffsets, targets, targetOffsets, matrix, gapOpen, gapExtend, scores);
        return LastBackend = $"CPU ({BioKernels.Backend})";
    }

    /// <summary><see cref="ILocalScoreBackend"/>: one query against packed targets (used by <see cref="SequenceDatabase.Search"/>).</summary>
    public void LocalScoreBatch(ReadOnlySpan<byte> query, ReadOnlySpan<byte> targets, ReadOnlySpan<int> offsets,
        ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
        => Score(query, [0, query.Length], targets, offsets, matrix, gapOpen, gapExtend, scores);

    internal static (byte[] Buffer, int[] Offsets) Pack(IReadOnlyList<Sequence> seqs)
    {
        long total = 0;
        foreach (var s in seqs) total += (s ?? throw new ArgumentException("Null sequence.", nameof(seqs))).Length;
        if (total > Array.MaxLength) throw new ArgumentException("Sequences exceed 2 GB; score in chunks.", nameof(seqs));
        var buf = new byte[total];
        var off = new int[seqs.Count + 1];
        for (int i = 0; i < seqs.Count; i++)
        {
            seqs[i].Span.CopyTo(buf.AsSpan(off[i]));
            off[i + 1] = off[i] + seqs[i].Length;
        }
        return (buf, off);
    }
}
