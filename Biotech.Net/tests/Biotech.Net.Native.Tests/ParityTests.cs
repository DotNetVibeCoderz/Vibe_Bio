using System.Text;

namespace Biotech.Net.Native.Tests;

/// <summary>
/// Native (Rust) vs managed (C#) parity on random inputs for every kernel. These call <see cref="NativeBridge"/> and
/// <see cref="ManagedKernels"/> directly (no global <see cref="BioKernels.ForceManaged"/> toggling) and are skipped when
/// the native library is not available.
/// </summary>
public class ParityTests
{
    private static void RequireNative()
    {
        if (!NativeBridge.IsAvailable) Assert.Skip($"Native library unavailable: {NativeBridge.UnavailableReason}");
    }

    private static byte[] RandomBytes(Random rng, int n, string alphabet)
    {
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte)alphabet[rng.Next(alphabet.Length)];
        return b;
    }

    private const string Dna = "ACGTacgtNnRYKMSWBDHVUu-*";
    private const string Aa = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwyXBZJUO*";

    [Fact]
    public void NativeLibrary_IsLoadedInTestOutput()
    {
        // The dll is copied next to the test assembly by the Native project's IncludeBiotechNative target.
        string? flag = Environment.GetEnvironmentVariable(NativeBridge.DisableEnvironmentVariable);
        if (flag == "1") Assert.Skip("Native disabled by environment.");
        if (!File.Exists(Path.Combine(AppContext.BaseDirectory, "biotech_native.dll")) && OperatingSystem.IsWindows())
            Assert.Skip("No native binary was built for this machine.");
        Assert.True(NativeBridge.IsAvailable, NativeBridge.UnavailableReason);
        Assert.StartsWith("biotech-native ", NativeBridge.Version);
        Assert.Equal(KernelBackend.Native, BioKernels.IsNativeAvailable ? KernelBackend.Native : KernelBackend.Managed);
    }

    [Fact]
    public void SequenceKernels_Match()
    {
        RequireNative();
        var rng = new Random(1);
        for (int round = 0; round < 200; round++)
        {
            int n = rng.Next(0, round < 100 ? 70 : 3000);
            var s = new byte[n];
            rng.NextBytes(s);
            if (round % 2 == 0) s = RandomBytes(rng, n, Dna);

            byte[] d1 = new byte[n], d2 = new byte[n];
            NativeBridge.ReverseComplement(s, d1); ManagedKernels.ReverseComplement(s, d2);
            Assert.Equal(d2, d1);
            NativeBridge.Complement(s, d1); ManagedKernels.Complement(s, d2);
            Assert.Equal(d2, d1);

            byte[] table = RandomBytes(rng, 64, "ACDEFGHIKLMNPQRSTVWY*");
            byte[] p1 = new byte[n / 3], p2 = new byte[n / 3];
            Assert.Equal(ManagedKernels.Translate(s, table, p2), NativeBridge.Translate(s, table, p1));
            Assert.Equal(p2, p1);

            Assert.Equal(ManagedKernels.CountNucleotides(s), NativeBridge.CountNucleotides(s));

            long[] h1 = new long[256], h2 = new long[256];
            NativeBridge.ByteHistogram(s, h1); ManagedKernels.ByteHistogram(s, h2);
            Assert.Equal(h2, h1);

            int k = 1 + round % 8;
            int[] k1 = new int[1 << (2 * k)], k2 = new int[1 << (2 * k)];
            Assert.Equal(ManagedKernels.CountKmers(s, k, k2), NativeBridge.CountKmers(s, k, k1));
            Assert.Equal(k2, k1);

            var t = (byte[])s.Clone();
            for (int i = 0; i < n; i += 1 + rng.Next(7)) t[i] = (byte)rng.Next(256);
            Assert.Equal(ManagedKernels.Hamming(s, t), NativeBridge.Hamming(s, t));
        }
    }

    [Fact]
    public void StatisticsKernels_Match()
    {
        RequireNative();
        var rng = new Random(2);
        int rows = 300, cols = 9;
        var m = new double[rows * cols];
        for (int i = 0; i < m.Length; i++) m[i] = rng.NextDouble() * 10 + (i % cols >= 5 ? rng.NextDouble() : 0);
        for (int c = 0; c < cols; c++) m[7 * cols + c] = 3; // zero-variance row
        int[] ga = [0, 1, 2, 3, 4], gb = [5, 6, 7, 8];
        double[] t1 = new double[rows], p1 = new double[rows], d1 = new double[rows];
        double[] t2 = new double[rows], p2 = new double[rows], d2 = new double[rows];
        NativeBridge.WelchTTest(m, rows, cols, ga, gb, t1, p1, d1);
        ManagedKernels.WelchTTest(m, rows, cols, ga, gb, t2, p2, d2);
        for (int r = 0; r < rows; r++)
        {
            Assert.Equal(t2[r], t1[r], 12);
            Assert.Equal(p2[r], p1[r], 12);
            Assert.Equal(d2[r], d1[r], 12);
        }

        var pv = new double[500];
        for (int i = 0; i < pv.Length; i++) pv[i] = i % 37 == 0 ? double.NaN : rng.NextDouble() * rng.NextDouble();
        double[] q1 = new double[pv.Length], q2 = new double[pv.Length];
        NativeBridge.BenjaminiHochberg(pv, q1);
        ManagedKernels.BenjaminiHochberg(pv, q2);
        for (int i = 0; i < pv.Length; i++) Assert.Equal(q2[i], q1[i], 14);
    }

    [Fact]
    public void ProteinKernels_Match()
    {
        RequireNative();
        var rng = new Random(3);
        for (int round = 0; round < 200; round++)
        {
            var s = RandomBytes(rng, rng.Next(0, 400), Aa);
            Assert.Equal(ManagedKernels.ProteinMass(s, false), NativeBridge.ProteinMass(s, false), 9);
            Assert.Equal(ManagedKernels.ProteinMass(s, true), NativeBridge.ProteinMass(s, true), 9);
            double ph = rng.NextDouble() * 14;
            Assert.Equal(ManagedKernels.NetCharge(s, ph), NativeBridge.NetCharge(s, ph), 10);
            Assert.Equal(ManagedKernels.IsoelectricPoint(s), NativeBridge.IsoelectricPoint(s), 10);
            int w = 1 + rng.Next(20);
            int len = Math.Max(0, s.Length - w + 1);
            double[] o1 = new double[len], o2 = new double[len];
            Assert.Equal(ManagedKernels.Hydropathy(s, w, o2), NativeBridge.Hydropathy(s, w, o1));
            for (int i = 0; i < len; i++) Assert.Equal(o2[i], o1[i], 10);
        }
    }

    private static int[] RandomMatrix(Random rng)
    {
        var m = new int[ManagedKernels.MatrixLength];
        for (int i = 0; i < m.Length; i++) m[i] = rng.Next(-6, 9);
        return m;
    }

    [Fact]
    public void Align_MatchesExactly_AllModes()
    {
        RequireNative();
        var rng = new Random(4);
        for (int round = 0; round < 400; round++)
        {
            int[] matrix = round % 3 == 0 ? BioKernels.MatchMismatchMatrix(rng.Next(1, 6), -rng.Next(0, 6)) : RandomMatrix(rng);
            string alphabet = round % 2 == 0 ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY";
            var a = RandomBytes(rng, rng.Next(0, 60), alphabet);
            var b = RandomBytes(rng, rng.Next(0, 60), alphabet);
            int o = rng.Next(0, 12), e = rng.Next(0, 5);
            long band = round % 4 == 0 ? rng.Next(0, 6) : -1;
            for (int mode = 0; mode < 3; mode++)
            {
                byte[] ops1 = new byte[a.Length + b.Length], ops2 = new byte[a.Length + b.Length];
                int[] r1 = new int[6], r2 = new int[6];
                Assert.Equal(0, NativeBridge.Align(a, b, matrix, o, e, mode, band, ops1, r1));
                Assert.Equal(0, ManagedKernels.Align(a, b, matrix, o, e, mode, band, ops2, r2));
                Assert.Equal(r1, r2);
                Assert.Equal(Encoding.ASCII.GetString(ops1, 0, r1[5]), Encoding.ASCII.GetString(ops2, 0, r2[5]));
            }
        }
    }

    [Fact]
    public void Align_ErrorCodes_Match()
    {
        RequireNative();
        var m = new int[ManagedKernels.MatrixLength];
        int[] r = new int[6];
        byte[] small = new byte[2];
        byte[] ac = "AC"u8.ToArray();
        Assert.Equal(ManagedKernels.ErrBufferTooSmall, NativeBridge.Align(ac, ac, m, 1, 1, 0, -1, small, r));
        Assert.Equal(ManagedKernels.ErrBufferTooSmall, ManagedKernels.Align(ac, ac, m, 1, 1, 0, -1, small, r));
        Assert.Equal(ManagedKernels.ErrBadArgument, NativeBridge.Align(ac, ac, m, 1, 1, 7, -1, new byte[4], r));
        Assert.Equal(ManagedKernels.ErrBadArgument, ManagedKernels.Align(ac, ac, m, 1, 1, 7, -1, new byte[4], r));
    }

    [Fact]
    public void LocalScore_AndBatch_Match()
    {
        RequireNative();
        var rng = new Random(5);
        for (int round = 0; round < 150; round++)
        {
            int[] matrix = RandomMatrix(rng);
            int o = rng.Next(0, 12), e = rng.Next(0, 5);
            var q = RandomBytes(rng, rng.Next(0, round < 140 ? 120 : 600), "ACDEFGHIKLMNPQRSTVWY");
            int count = rng.Next(1, 8);
            var targets = new List<byte>();
            var offsets = new int[count + 1];
            for (int k = 0; k < count; k++)
            {
                targets.AddRange(RandomBytes(rng, rng.Next(0, 150), "ACDEFGHIKLMNPQRSTVWY"));
                offsets[k + 1] = targets.Count;
            }
            var tb = targets.ToArray();
            int[] s1 = new int[count], s2 = new int[count];
            Assert.Equal(0, NativeBridge.LocalScoreBatch(q, tb, offsets, matrix, o, e, s1));
            Assert.Equal(0, ManagedKernels.LocalScoreBatch(q, tb, offsets, matrix, o, e, s2));
            Assert.Equal(s2, s1);
            for (int k = 0; k < count; k++)
            {
                var t = tb.AsSpan(offsets[k], offsets[k + 1] - offsets[k]);
                Assert.Equal(s2[k], NativeBridge.LocalScore(q, t, matrix, o, e));
                Assert.Equal(s2[k], ManagedKernels.LocalScore(q, t, matrix, o, e));
                // Score-only equals the traceback aligner's local score.
                var r = new int[6];
                ManagedKernels.Align(q, t, matrix, o, e, 1, -1, new byte[q.Length + t.Length], r);
                Assert.Equal(r[0], s2[k]);
            }
        }
    }

    [Fact]
    public void LocalScore_SaturationFallback_Matches()
    {
        RequireNative();
        var m = BioKernels.MatchMismatchMatrix(20, -5);
        var q = Enumerable.Repeat((byte)'A', 3000).ToArray();
        Assert.Equal(60000, NativeBridge.LocalScore(q, q, m, 5, 1));
        Assert.Equal(60000, ManagedKernels.LocalScore(q, q, m, 5, 1));
    }

    [Fact]
    public void KabschRmsd_MatchesBitForBit()
    {
        RequireNative();
        var rng = new Random(6);
        for (int round = 0; round < 200; round++)
        {
            int n = rng.Next(0, 60);
            var a = new double[3 * n];
            var b = new double[3 * n];
            for (int i = 0; i < a.Length; i++) { a[i] = rng.NextDouble() * 40 - 20; b[i] = a[i] * 0.7 + rng.NextDouble() * 3; }
            double[] r1 = new double[9], r2 = new double[9];
            double x1 = NativeBridge.KabschRmsd(a, b, n, r1);
            double x2 = ManagedKernels.KabschRmsd(a, b, n, r2);
            Assert.Equal(x2, x1);
            Assert.Equal(r2, r1);
            Assert.Equal(x2, NativeBridge.KabschRmsd(a, b, n, default));
        }
    }
}
