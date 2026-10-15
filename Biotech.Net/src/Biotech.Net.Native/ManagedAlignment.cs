namespace Biotech.Net.Native;

// Managed twins of native/biotech-native/src/align.rs and structure.rs. They mirror the Rust code operation for
// operation (including traceback tie-breaking and floating-point evaluation order) so both backends agree exactly.
internal static partial class ManagedKernels
{
    internal const int MatrixLength = 128 * 128;
    internal const int ErrOutOfMemory = -2, ErrBufferTooSmall = -3, ErrBadArgument = -4;
    private const int Neg = int.MinValue / 4;

    private static int Gap(int o, int e, long k) => o + (int)(k - 1) * e;

    /// <summary>
    /// Gotoh affine-gap alignment with traceback (align.rs <c>align</c>). <paramref name="result"/> receives score,
    /// a_start, a_end, b_start, b_end, ops_len. Returns 0 or a negative error code.
    /// </summary>
    public static int Align(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, ReadOnlySpan<int> matrix, int o, int e, int mode,
        long band, Span<byte> ops, Span<int> result)
    {
        if (mode is < 0 or > 2) return ErrBadArgument;
        if (result.Length < 6) return ErrBadArgument;
        int n = a.Length, m = b.Length;
        if (matrix.Length < MatrixLength || o < 0 || e < 0) return ErrBadArgument;
        if (ops.Length < (long)n + m) return ErrBufferTooSmall;
        bool global = mode == 0, local = mode == 1, semi = mode == 2;
        if (n == 0 || m == 0)
        {
            if (local) { result[..6].Clear(); return 0; }
            int k0 = n + m;
            ops[..n].Fill((byte)'D');
            ops[n..k0].Fill((byte)'I');
            result[0] = global && k0 > 0 ? -Gap(o, e, k0) : 0;
            result[1] = 0; result[2] = n; result[3] = 0; result[4] = m; result[5] = k0;
            return 0;
        }

        // Band: dmin ≤ j − i ≤ dmax, widened by |m − n|.
        bool banded = band >= 0;
        long diff = (long)m - n;
        long dmin = banded ? Math.Min(diff, 0) - band : 0;
        long dmax = banded ? Math.Max(diff, 0) + band : 0;
        bool Contains(long i, long j) { long d = j - i; return !banded || (dmin <= d && d <= dmax); }
        int Boundary(long k) => global ? -Gap(o, e, k) : 0;

        var rowLo = new int[n + 1];
        var rowOff = new long[n + 2];
        for (int i = 1; i <= n; i++)
        {
            long lo, hi;
            if (!banded) { lo = 1; hi = m; }
            else { lo = Math.Max(i + dmin, 1); hi = Math.Min(i + dmax, m); }
            rowLo[i] = (int)lo;
            rowOff[i + 1] = rowOff[i] + (hi + 1 - lo);
        }
        if (rowOff[n + 1] > Array.MaxLength) return ErrOutOfMemory;
        byte[] tb;
        int[] hprev, hcur, f;
        try
        {
            tb = new byte[rowOff[n + 1]];
            hprev = new int[m + 1]; hcur = new int[m + 1]; f = new int[m + 1];
        }
        catch (OutOfMemoryException) { return ErrOutOfMemory; }
        Array.Fill(hprev, Neg); Array.Fill(hcur, Neg); Array.Fill(f, Neg);

        hprev[0] = 0;
        for (int j = 1; j <= m; j++) hprev[j] = Contains(0, j) ? Boundary(j) : Neg;
        int best = 0, bi = 0, bj = 0;
        int colBest = hprev[m], colI = 0;

        for (int i = 1; i <= n; i++)
        {
            int lo = rowLo[i];
            int hi = lo + (int)(rowOff[i + 1] - rowOff[i]) - 1;
            int srow = (a[i - 1] & 0x7F) << 7;
            int left = lo == 1 && Contains(i, 0) ? Boundary(i) : Neg;
            hcur[lo - 1] = left;
            int hLeft = left, eRun = Neg, diag = hprev[lo - 1];
            long trow = rowOff[i];
            for (int j = lo; j <= hi; j++)
            {
                int s = matrix[srow + (b[j - 1] & 0x7F)];
                int eOpen = hLeft - o, eExt = eRun - e;
                int ev; byte eb;
                if (eExt > eOpen) { ev = eExt; eb = 4; } else { ev = eOpen; eb = 0; }
                int up = hprev[j];
                int fOpen = up - o, fExt = f[j] - e;
                int fv; byte fb;
                if (fExt > fOpen) { fv = fExt; fb = 8; } else { fv = fOpen; fb = 0; }
                int d = diag + s;
                int h; byte src;
                if (d >= ev && d >= fv) { h = d; src = 0; }
                else if (ev >= fv) { h = ev; src = 1; }
                else { h = fv; src = 2; }
                if (local && h <= 0) { h = 0; src = 3; }
                tb[trow + j - lo] = (byte)(src | eb | fb);
                diag = up;
                hcur[j] = h;
                f[j] = fv;
                eRun = ev;
                hLeft = h;
                if (local && h > best) { best = h; bi = i; bj = j; }
            }
            if (hi < m) { hcur[hi + 1] = Neg; f[hi + 1] = Neg; }
            if (semi && i < n && hi == m && hcur[m] > colBest) { colBest = hcur[m]; colI = i; }
            (hprev, hcur) = (hcur, hprev);
        }

        int score, ei, ej;
        if (global) { score = hprev[m]; ei = n; ej = m; }
        else if (local) { score = best; ei = bi; ej = bj; }
        else
        {
            score = hprev[m]; ei = n; ej = m;
            int first = rowLo[n] - 1;
            for (int j = first; j < m; j++)
                if (hprev[j] > score) { score = hprev[j]; ej = j; }
            if (colBest > score) { score = colBest; ei = colI; ej = m; }
        }

        int k = 0;
        int ti = ei, tj = ej;
        if (semi)
        {
            for (int x = tj; x < m; x++) ops[k++] = (byte)'I';
            for (int x = ti; x < n; x++) ops[k++] = (byte)'D';
        }
        int state = 0;
        while (ti > 0 && tj > 0)
        {
            int lo = rowLo[ti];
            long width = rowOff[ti + 1] - rowOff[ti];
            if (tj < lo || tj >= lo + width) break;
            byte t = tb[rowOff[ti] + tj - lo];
            if (state == 0)
            {
                int src = t & 3;
                if (src == 0) { ops[k++] = (byte)'M'; ti--; tj--; }
                else if (src == 1) state = 1;
                else if (src == 2) state = 2;
                else break;
            }
            else if (state == 1)
            {
                ops[k++] = (byte)'I'; tj--;
                state = (t & 4) != 0 ? 1 : 0;
            }
            else
            {
                ops[k++] = (byte)'D'; ti--;
                state = (t & 8) != 0 ? 2 : 0;
            }
        }
        int sa = ti, sb = tj;
        if (!local)
        {
            for (int x = 0; x < ti; x++) ops[k++] = (byte)'D';
            for (int x = 0; x < tj; x++) ops[k++] = (byte)'I';
        }
        ops[..k].Reverse();
        result[0] = score;
        if (local) { result[1] = sa; result[2] = ei; result[3] = sb; result[4] = ej; }
        else { result[1] = 0; result[2] = n; result[3] = 0; result[4] = m; }
        result[5] = k;
        return 0;
    }

    /// <summary>Linear-memory local alignment score (align.rs <c>LocalProfile::scalar_score</c>).</summary>
    public static int LocalScore(ReadOnlySpan<byte> q, ReadOnlySpan<byte> t, ReadOnlySpan<int> matrix, int o, int e)
    {
        if (matrix.Length < MatrixLength || o < 0 || e < 0) return ErrBadArgument;
        int[] profile = BuildProfile(q, matrix);
        return ScoreProfile(profile, q.Length, t, o, e, new int[q.Length], new int[q.Length]);
    }

    /// <summary>Scores <paramref name="q"/> against <c>targets[offsets[k]..offsets[k+1]]</c> (align.rs <c>local_score_batch</c>).</summary>
    public static int LocalScoreBatch(ReadOnlySpan<byte> q, ReadOnlySpan<byte> targets, ReadOnlySpan<int> offsets,
        ReadOnlySpan<int> matrix, int o, int e, Span<int> scores)
    {
        if (offsets.IsEmpty || scores.Length < offsets.Length - 1) return ErrBadArgument;
        for (int k = 0; k + 1 < offsets.Length; k++)
            if (offsets[k] < 0 || offsets[k + 1] < offsets[k] || offsets[k + 1] > targets.Length) return ErrBadArgument;
        if (matrix.Length < MatrixLength || o < 0 || e < 0) return ErrBadArgument;
        int[] profile = BuildProfile(q, matrix);
        int[] h = new int[q.Length], f = new int[q.Length];
        for (int k = 0; k + 1 < offsets.Length; k++)
            scores[k] = ScoreProfile(profile, q.Length, targets[offsets[k]..offsets[k + 1]], o, e, h, f);
        return 0;
    }

    /// <summary><c>profile[c * qlen + j] = matrix[q[j]][c]</c>: one contiguous row per target residue.</summary>
    private static int[] BuildProfile(ReadOnlySpan<byte> q, ReadOnlySpan<int> matrix)
    {
        int qlen = q.Length;
        var p = new int[128 * qlen];
        for (int c = 0; c < 128; c++)
            for (int j = 0; j < qlen; j++)
                p[c * qlen + j] = matrix[((q[j] & 0x7F) << 7) | c];
        return p;
    }

    private static int ScoreProfile(int[] profile, int qlen, ReadOnlySpan<byte> t, int o, int e, int[] h, int[] f)
    {
        if (qlen == 0 || t.IsEmpty) return 0;
        Array.Clear(h); Array.Clear(f);
        int best = 0;
        foreach (byte c in t)
        {
            var row = profile.AsSpan((c & 0x7F) * qlen, qlen);
            int diag = 0, hLeft = 0, ev = 0;
            for (int j = 0; j < qlen; j++)
            {
                ev = Math.Max(Math.Max(hLeft - o, ev - e), 0);
                int up = h[j];
                int fv = Math.Max(Math.Max(up - o, f[j] - e), 0);
                f[j] = fv;
                int hv = Math.Max(Math.Max(Math.Max(diag + row[j], ev), fv), 0);
                diag = up;
                h[j] = hv;
                hLeft = hv;
                if (hv > best) best = hv;
            }
        }
        return best;
    }

    /// <summary>
    /// Optimal-superposition RMSD (Horn's quaternion method, structure.rs <c>kabsch_rmsd</c>). <paramref name="a"/> is the
    /// target, <paramref name="b"/> the mobile set; <paramref name="rotation"/> (empty or 9 values) receives the row-major R
    /// with <c>R·(b − b̄) + ā ≈ a</c>.
    /// </summary>
    public static double KabschRmsd(ReadOnlySpan<double> a, ReadOnlySpan<double> b, int n, Span<double> rotation)
    {
        if (n == 0)
        {
            if (rotation.Length >= 9) { rotation[..9].Clear(); rotation[0] = rotation[4] = rotation[8] = 1.0; }
            return 0.0;
        }
        Span<double> ca = stackalloc double[3], cb = stackalloc double[3];
        ca.Clear(); cb.Clear();
        for (int i = 0; i < n; i++)
            for (int k = 0; k < 3; k++) { ca[k] += a[3 * i + k]; cb[k] += b[3 * i + k]; }
        double nf = n;
        for (int k = 0; k < 3; k++) { ca[k] /= nf; cb[k] /= nf; }

        Span<double> s = stackalloc double[9];
        s.Clear();
        double ga = 0.0, gb = 0.0;
        Span<double> x = stackalloc double[3], y = stackalloc double[3];
        for (int i = 0; i < n; i++)
        {
            x[0] = b[3 * i] - cb[0]; x[1] = b[3 * i + 1] - cb[1]; x[2] = b[3 * i + 2] - cb[2];
            y[0] = a[3 * i] - ca[0]; y[1] = a[3 * i + 1] - ca[1]; y[2] = a[3 * i + 2] - ca[2];
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++) s[3 * p + q] += x[p] * y[q];
            ga += y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
            gb += x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        }
        double sxx = s[0], sxy = s[1], sxz = s[2], syx = s[3], syy = s[4], syz = s[5], szx = s[6], szy = s[7], szz = s[8];
        Span<double> mm =
        [
            sxx + syy + szz, syz - szy, szx - sxz, sxy - syx,
            syz - szy, sxx - syy - szz, sxy + syx, szx + sxz,
            szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy,
            sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz,
        ];
        Span<double> v = stackalloc double[16];
        v.Clear();
        for (int k = 0; k < 4; k++) v[5 * k] = 1.0;
        Jacobi4(mm, v);
        int bestK = 0;
        for (int k = 1; k < 4; k++) if (mm[5 * k] > mm[5 * bestK]) bestK = k;
        double lambda = mm[5 * bestK];
        double msd = (ga + gb - 2.0 * lambda) / nf;
        if (rotation.Length >= 9)
        {
            double q0 = v[bestK], q1 = v[4 + bestK], q2 = v[8 + bestK], q3 = v[12 + bestK];
            double norm = Math.Sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
            rotation[0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
            rotation[1] = 2.0 * (q1 * q2 - q0 * q3);
            rotation[2] = 2.0 * (q1 * q3 + q0 * q2);
            rotation[3] = 2.0 * (q1 * q2 + q0 * q3);
            rotation[4] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
            rotation[5] = 2.0 * (q2 * q3 - q0 * q1);
            rotation[6] = 2.0 * (q1 * q3 - q0 * q2);
            rotation[7] = 2.0 * (q2 * q3 + q0 * q1);
            rotation[8] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
        }
        return msd > 0.0 ? Math.Sqrt(msd) : 0.0;
    }

    /// <summary>Cyclic Jacobi eigen-decomposition of a symmetric 4×4 (row-major) matrix; mirrors structure.rs <c>jacobi4</c>.</summary>
    private static void Jacobi4(Span<double> a, Span<double> v)
    {
        double total = 0.0;
        for (int k = 0; k < 16; k++) total += a[k] * a[k];
        if (total == 0.0) return;
        for (int sweep = 0; sweep < 64; sweep++)
        {
            double off = 0.0;
            for (int p = 0; p < 3; p++)
                for (int q = p + 1; q < 4; q++) off += a[4 * p + q] * a[4 * p + q];
            if (off <= 1e-30 * total) break;
            for (int p = 0; p < 3; p++)
            {
                for (int q = p + 1; q < 4; q++)
                {
                    double apq = a[4 * p + q];
                    if (apq == 0.0) continue;
                    double theta = (a[4 * q + q] - a[4 * p + p]) / (2.0 * apq);
                    double t;
                    if (Math.Abs(theta) > 1e150) t = 0.5 / theta;
                    else
                    {
                        double sgn = theta >= 0.0 ? 1.0 : -1.0;
                        t = sgn / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1.0));
                    }
                    double c = 1.0 / Math.Sqrt(t * t + 1.0);
                    double s = t * c;
                    for (int k = 0; k < 4; k++)
                    {
                        double akp = a[4 * k + p], akq = a[4 * k + q];
                        a[4 * k + p] = c * akp - s * akq;
                        a[4 * k + q] = s * akp + c * akq;
                    }
                    for (int k = 0; k < 4; k++)
                    {
                        double apk = a[4 * p + k], aqk = a[4 * q + k];
                        a[4 * p + k] = c * apk - s * aqk;
                        a[4 * q + k] = s * apk + c * aqk;
                    }
                    for (int k = 0; k < 4; k++)
                    {
                        double vkp = v[4 * k + p], vkq = v[4 * k + q];
                        v[4 * k + p] = c * vkp - s * vkq;
                        v[4 * k + q] = s * vkp + c * vkq;
                    }
                }
            }
        }
    }
}
