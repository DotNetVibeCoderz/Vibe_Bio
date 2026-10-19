using System.ComponentModel;
using ILGPU;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu.Alignment;

/// <summary>Kernel arguments for <see cref="LocalScoreKernel.Run"/>. Infrastructure: public only because ILGPU's CPU accelerator
/// emits the launcher into a dynamic assembly that cannot see internal types.</summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public struct LocalScoreArgs
{
    /// <summary>Kernel field.</summary>
    public ArrayView<byte> Queries;      // compact residue codes, back to back
    public ArrayView<int> QueryOffsets;  // count + 1
    public ArrayView<byte> Targets;      // compact codes, targets sorted by length (longest first)
    public ArrayView<int> TargetOffsets; // count + 1
    public ArrayView<int> Matrix;        // K×K, row = query code
    public ArrayView<int> H;             // scratch: Stride × maxQueryLength, element (j, lane) at j·Stride + lane
    public ArrayView<int> F;
    public ArrayView<int> Scores;        // one per pair in the launch
    public int K, GapOpen, GapExtend, TargetCount, Stride;
    public long FirstPair;
}

/// <summary>
/// Inter-sequence Smith–Waterman/Gotoh score kernel: one thread per (query, target) pair, the same recurrence as
/// <c>ManagedKernels.ScoreProfile</c> (target residues outer, query inner, E/F/H clamped at zero), so scores are the
/// exact integers of <see cref="Native.BioKernels.LocalScore"/>. The DP row along the query lives in global scratch,
/// interleaved by lane so neighbouring threads touch neighbouring words; four target residues are processed per sweep
/// to cut that traffic fourfold.
/// </summary>
/// <remarks>Infrastructure (public for ILGPU's CPU accelerator); use <see cref="GpuLocalAligner"/>.</remarks>
[EditorBrowsable(EditorBrowsableState.Never)]
public static class LocalScoreKernel
{
    /// <summary>One lane = one (query, target) pair of the launch.</summary>
    public static void Run(Index1D lane, LocalScoreArgs a)
    {
        long pair = a.FirstPair + lane;
        int qi = (int)(pair / a.TargetCount);
        int ti = (int)(pair - (long)qi * a.TargetCount);
        int q0 = a.QueryOffsets[qi], qlen = a.QueryOffsets[qi + 1] - q0;
        int t0 = a.TargetOffsets[ti], tlen = a.TargetOffsets[ti + 1] - t0;
        int stride = a.Stride, o = a.GapOpen, e = a.GapExtend, k = a.K;
        if (qlen == 0 || tlen == 0) { a.Scores[lane] = 0; return; }
        for (int j = 0; j < qlen; j++)
        {
            a.H[(long)j * stride + lane] = 0;
            a.F[(long)j * stride + lane] = 0;
        }
        int best = 0;
        int c = 0;
        for (; c + 4 <= tlen; c += 4)
        {
            int r0 = a.Targets[t0 + c], r1 = a.Targets[t0 + c + 1], r2 = a.Targets[t0 + c + 2], r3 = a.Targets[t0 + c + 3];
            int d0 = 0, l0 = 0, e0 = 0, d1 = 0, l1 = 0, e1 = 0, d2 = 0, l2 = 0, e2 = 0, d3 = 0, l3 = 0, e3 = 0;
            for (int j = 0; j < qlen; j++)
            {
                long at = (long)j * stride + lane;
                int row = a.Queries[q0 + j] * k;
                int up = a.H[at], fu = a.F[at];

                e0 = Max(Max(l0 - o, e0 - e), 0);
                int f0 = Max(Max(up - o, fu - e), 0);
                int h0 = Max(Max(Max(d0 + a.Matrix[row + r0], e0), f0), 0);
                d0 = up; l0 = h0;

                e1 = Max(Max(l1 - o, e1 - e), 0);
                int f1 = Max(Max(h0 - o, f0 - e), 0);
                int h1 = Max(Max(Max(d1 + a.Matrix[row + r1], e1), f1), 0);
                d1 = h0; l1 = h1;

                e2 = Max(Max(l2 - o, e2 - e), 0);
                int f2 = Max(Max(h1 - o, f1 - e), 0);
                int h2 = Max(Max(Max(d2 + a.Matrix[row + r2], e2), f2), 0);
                d2 = h1; l2 = h2;

                e3 = Max(Max(l3 - o, e3 - e), 0);
                int f3 = Max(Max(h2 - o, f2 - e), 0);
                int h3 = Max(Max(Max(d3 + a.Matrix[row + r3], e3), f3), 0);
                d3 = h2; l3 = h3;

                best = Max(best, Max(Max(h0, h1), Max(h2, h3)));
                a.H[at] = h3;
                a.F[at] = f3;
            }
        }
        for (; c < tlen; c++)
        {
            int r = a.Targets[t0 + c];
            int diag = 0, left = 0, ev = 0;
            for (int j = 0; j < qlen; j++)
            {
                long at = (long)j * stride + lane;
                ev = Max(Max(left - o, ev - e), 0);
                int up = a.H[at];
                int fv = Max(Max(up - o, a.F[at] - e), 0);
                a.F[at] = fv;
                int hv = Max(Max(Max(diag + a.Matrix[a.Queries[q0 + j] * k + r], ev), fv), 0);
                diag = up;
                a.H[at] = hv;
                left = hv;
                best = Max(best, hv);
            }
        }
        a.Scores[lane] = best;
    }

    private static int Max(int x, int y) => x > y ? x : y;
}
