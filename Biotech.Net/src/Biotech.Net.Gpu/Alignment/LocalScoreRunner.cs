using Biotech.Net.Native;
using ILGPU;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu.Alignment;

/// <summary>Host side of the all-vs-all score kernel: packing, chunking, launch, and the CPU fallback.</summary>
internal static class LocalScoreRunner
{
    /// <summary>Scratch budget per launch (H and F rows of every lane), in ints.</summary>
    internal static long ScratchBudget = 16L * 1024 * 1024;

    /// <summary>
    /// Scores every query against every target: <c>scores[q·targetCount + t]</c> = local score of query q (matrix rows)
    /// against target t. Sequences are packed back to back with count + 1 offsets.
    /// </summary>
    public static void RunGpu(GpuDevice device, ReadOnlySpan<byte> queries, ReadOnlySpan<int> queryOffsets,
        ReadOnlySpan<byte> targets, ReadOnlySpan<int> targetOffsets, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
    {
        int nq = queryOffsets.Length - 1, nt = targetOffsets.Length - 1;
        long total = (long)nq * nt;
        if (total == 0) return;

        // compact alphabet: only residues that occur get a code, so the K×K matrix stays small and cache-resident
        Span<int> code = stackalloc int[128];
        code.Fill(-1);
        int k = 0;
        foreach (byte b in queries) if (code[b & 0x7F] < 0) code[b & 0x7F] = k++;
        foreach (byte b in targets) if (code[b & 0x7F] < 0) code[b & 0x7F] = k++;
        if (k == 0) k = 1;
        Span<int> residue = stackalloc int[k];
        for (int r = 0; r < 128; r++) if (code[r] >= 0) residue[code[r]] = r;
        var cm = new int[k * k];
        for (int x = 0; x < k; x++)
            for (int y = 0; y < k; y++) cm[x * k + y] = matrix[(residue[x] << 7) | residue[y]];

        var qc = new byte[queries.Length];
        for (int i = 0; i < queries.Length; i++) qc[i] = (byte)code[queries[i] & 0x7F];
        int maxQ = 0;
        for (int q = 0; q < nq; q++) maxQ = Math.Max(maxQ, queryOffsets[q + 1] - queryOffsets[q]);

        // targets sorted longest first, so the threads of a warp run similar trip counts
        var order = new int[nt];
        var lengths = new int[nt];
        for (int t = 0; t < nt; t++) { order[t] = t; lengths[t] = targetOffsets[t + 1] - targetOffsets[t]; }
        var keys = lengths.Select(l => -l).ToArray();
        Array.Sort(keys, order);
        var tc = new byte[targets.Length];
        var toff = new int[nt + 1];
        for (int s = 0; s < nt; s++)
        {
            int t = order[s], start = targetOffsets[t], len = lengths[t];
            for (int i = 0; i < len; i++) tc[toff[s] + i] = (byte)code[targets[start + i] & 0x7F];
            toff[s + 1] = toff[s] + len;
        }

        int stride = (int)Math.Clamp(ScratchBudget / Math.Max(1, maxQ), 1, 1 << 20);
        if (stride > total) stride = (int)total;
        var chunk = new int[stride];
        var qoff = queryOffsets.ToArray();

        lock (device.Sync)
        {
            device.ThrowIfDisposed();
            var acc = device.Accelerator;
            var kernel = device.Kernel("local-score", a => a.LoadAutoGroupedStreamKernel<Index1D, LocalScoreArgs>(LocalScoreKernel.Run));
            using var dq = acc.Allocate1D(qc.Length == 0 ? new byte[1] : qc);
            using var dqo = acc.Allocate1D(qoff);
            using var dt = acc.Allocate1D(tc.Length == 0 ? new byte[1] : tc);
            using var dto = acc.Allocate1D(toff);
            using var dm = acc.Allocate1D(cm);
            using var dh = acc.Allocate1D<int>((long)stride * Math.Max(1, maxQ));
            using var df = acc.Allocate1D<int>((long)stride * Math.Max(1, maxQ));
            using var ds = acc.Allocate1D<int>(stride);
            var args = new LocalScoreArgs
            {
                Queries = dq.View, QueryOffsets = dqo.View, Targets = dt.View, TargetOffsets = dto.View, Matrix = dm.View,
                H = dh.View, F = df.View, Scores = ds.View, K = k, GapOpen = gapOpen, GapExtend = gapExtend,
                TargetCount = nt, Stride = stride,
            };
            for (long first = 0; first < total; first += stride)
            {
                int count = (int)Math.Min(stride, total - first);
                args.FirstPair = first;
                kernel(count, args);
                acc.Synchronize();
                ds.View.AsContiguous().SubView(0, count).CopyToCPU(chunk.AsSpan(0, count));
                for (int i = 0; i < count; i++)
                {
                    long pair = first + i;
                    int q = (int)(pair / nt), s = (int)(pair - (long)q * nt);
                    scores[q * nt + order[s]] = chunk[i];
                }
            }
        }
    }

    /// <summary>The CPU path: <see cref="BioKernels.LocalScoreBatch"/> (Rust AVX2 or managed) per query, queries in parallel.</summary>
    public static void RunCpu(ReadOnlySpan<byte> queries, ReadOnlySpan<int> queryOffsets, ReadOnlySpan<byte> targets,
        ReadOnlySpan<int> targetOffsets, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores, int maxDegreeOfParallelism = -1)
    {
        int nq = queryOffsets.Length - 1, nt = targetOffsets.Length - 1;
        if ((long)nq * nt == 0) return;
        var q = queries.ToArray();
        var qo = queryOffsets.ToArray();
        var t = targets.ToArray();
        var to = targetOffsets.ToArray();
        var m = matrix.ToArray();
        var result = new int[nq * nt];
        Parallel.For(0, nq, new ParallelOptions { MaxDegreeOfParallelism = maxDegreeOfParallelism }, i =>
            BioKernels.LocalScoreBatch(q.AsSpan(qo[i], qo[i + 1] - qo[i]), t, to, m, gapOpen, gapExtend, result.AsSpan(i * nt, nt)));
        result.CopyTo(scores);
    }

    public static void Validate(ReadOnlySpan<byte> queries, ReadOnlySpan<int> queryOffsets, ReadOnlySpan<byte> targets,
        ReadOnlySpan<int> targetOffsets, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
    {
        if (matrix.Length < BioKernels.SubstitutionMatrixLength)
            throw new ArgumentException($"A kernel substitution matrix has {BioKernels.SubstitutionMatrixLength} entries (128×128).", nameof(matrix));
        if (gapOpen < 0) throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap penalties are non-negative costs.");
        if (gapExtend < 0) throw new ArgumentOutOfRangeException(nameof(gapExtend), "Gap penalties are non-negative costs.");
        CheckOffsets(queries.Length, queryOffsets, nameof(queryOffsets));
        CheckOffsets(targets.Length, targetOffsets, nameof(targetOffsets));
        long pairs = (long)(queryOffsets.Length - 1) * (targetOffsets.Length - 1);
        if (pairs > Array.MaxLength) throw new ArgumentException("Too many pairs for one score matrix; split the queries.");
        if (scores.Length < pairs) throw new ArgumentException($"The score buffer needs {pairs} entries.", nameof(scores));
    }

    private static void CheckOffsets(int length, ReadOnlySpan<int> offsets, string name)
    {
        if (offsets.IsEmpty) throw new ArgumentException("Offsets need at least one entry (count + 1).", name);
        if (offsets[0] < 0) throw new ArgumentException("Offsets must be ascending positions inside the buffer.", name);
        for (int i = 0; i + 1 < offsets.Length; i++)
            if (offsets[i + 1] < offsets[i] || offsets[i + 1] > length)
                throw new ArgumentException("Offsets must be ascending positions inside the buffer.", name);
    }
}
