namespace Biotech.Net.Native.Tests;

/// <summary>Known-value tests. They run against the managed kernels directly so they pass with or without the Rust library.</summary>
public class KernelTests
{
    private static (int[] R, string Ops) Run(string a, string b, int[] m, int o, int e, int mode, long band = -1)
    {
        var ops = new byte[a.Length + b.Length];
        var r = new int[6];
        Assert.Equal(0, ManagedKernels.Align(System.Text.Encoding.ASCII.GetBytes(a), System.Text.Encoding.ASCII.GetBytes(b), m, o, e, mode, band, ops, r));
        return (r, System.Text.Encoding.ASCII.GetString(ops, 0, r[5]));
    }

    [Fact]
    public void Global_GattacaExample()
    {
        var m = BioKernels.MatchMismatchMatrix(1, -1);
        Assert.Equal(0, Run("GATTACA", "GCATGCU", m, 1, 1, 0).R[0]);
    }

    [Fact]
    public void Local_FindsEmbeddedMatch()
    {
        var (r, ops) = Run("TTTTACGTACGTTTTT", "GGACGTACGGG", BioKernels.MatchMismatchMatrix(2, -3), 5, 2, 1);
        Assert.Equal(14, r[0]);
        Assert.Equal("MMMMMMM", ops);
        Assert.Equal([4, 11, 2, 9], r[1..5]);
    }

    [Fact]
    public void SemiGlobal_IgnoresEndGaps()
    {
        var (r, ops) = Run("ACGTACGT", "GTAC", BioKernels.MatchMismatchMatrix(1, -2), 3, 1, 2);
        Assert.Equal(4, r[0]);
        Assert.Equal("DDMMMMDD", ops);
    }

    [Fact]
    public void Affine_PrefersOneLongGap()
    {
        var (r, ops) = Run("AAAACCCCGGGG", "AAAAGGGG", BioKernels.MatchMismatchMatrix(5, -4), 10, 1, 0);
        Assert.Equal(8 * 5 - 13, r[0]);
        Assert.Equal("MMMMDDDDMMMM", ops);
    }

    [Fact]
    public void EmptyInputs()
    {
        var m = BioKernels.MatchMismatchMatrix(1, -1);
        var (r, ops) = Run("", "ACG", m, 4, 1, 0);
        Assert.Equal(-6, r[0]);
        Assert.Equal("III", ops);
        Assert.Equal(0, Run("", "", m, 4, 1, 0).R[0]);
        (r, ops) = Run("ACG", "", m, 4, 1, 1);
        Assert.Equal(0, r[0]);
        Assert.Equal("", ops);
    }

    [Fact]
    public void Banded_EqualsFull_WhenWide_AndNeverBetterWhenNarrow()
    {
        var m = BioKernels.MatchMismatchMatrix(2, -3);
        var rng = new Random(11);
        for (int k = 0; k < 50; k++)
        {
            string a = new([.. Enumerable.Range(0, 40).Select(_ => "ACGT"[rng.Next(4)])]);
            string b = new([.. a.Where(_ => rng.Next(10) != 0).Select(c => rng.Next(10) == 0 ? 'T' : c)]);
            for (int mode = 0; mode < 3; mode++)
            {
                var full = Run(a, b, m, 5, 2, mode);
                var banded = Run(a, b, m, 5, 2, mode, 100);
                Assert.Equal(full.R, banded.R);
                Assert.Equal(full.Ops, banded.Ops);
            }
            Assert.True(Run(a, b, m, 5, 2, 0, 2).R[0] <= Run(a, b, m, 5, 2, 0).R[0]);
        }
    }

    [Fact]
    public void LocalScoreBatch_KnownValues()
    {
        var m = BioKernels.MatchMismatchMatrix(1, -1);
        int[] scores = new int[3];
        BioKernels.LocalScoreBatch("ACG"u8, "ACGTTTACG"u8, [0, 4, 4, 9], m, 2, 1, scores);
        Assert.Equal([3, 0, 3], scores);
        Assert.Throws<ArgumentException>(() => BioKernels.LocalScoreBatch("ACG"u8, "ACGTTTACG"u8, [0, 10], m, 2, 1, new int[1]));
    }

    [Fact]
    public void BioKernelsAlign_PublicApi()
    {
        var m = BioKernels.MatchMismatchMatrix(2, -3);
        var aln = BioKernels.Align("TTTTACGTACGTTTTT"u8, "ggacgtacggg"u8, m, 5, 2, AlignmentMode.Local);
        Assert.Equal(14, aln.Score); // matrix is case-insensitive
        Assert.Equal("MMMMMMM", aln.OperationString);
        Assert.Equal(new AlignmentSummary(14, 4, 11, 2, 9, 7), aln.Summary);
        Assert.Equal(14, BioKernels.LocalScore("TTTTACGTACGTTTTT"u8, "GGACGTACGGG"u8, m, 5, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => BioKernels.Align("A"u8, "A"u8, m, -1, 0));
        Assert.Throws<ArgumentException>(() => BioKernels.Align("A"u8, "A"u8, new int[10], 1, 0));
    }

    private static double[] Rotate(double[] r, double[] p)
    {
        var o = new double[p.Length];
        for (int k = 0; k < p.Length; k += 3)
            for (int i = 0; i < 3; i++) o[k + i] = r[3 * i] * p[k] + r[3 * i + 1] * p[k + 1] + r[3 * i + 2] * p[k + 2];
        return o;
    }

    [Fact]
    public void Kabsch_RotatedTranslatedCopy_HasZeroRmsd_AndRecoversInverseRotation()
    {
        double[] a = [1, 2, 3, -1, 0.5, 2, 4, -2, 0, 0, 0, 1, 2.5, 1.5, -1];
        double c = Math.Cos(0.6), s = Math.Sin(0.6);
        double[] rz = [c, -s, 0, s, c, 0, 0, 0, 1];
        double[] shift = [5, -3, 2];
        var b = Rotate(rz, a).Select((x, i) => x + shift[i % 3]).ToArray();
        var r = new double[9];
        Assert.True(ManagedKernels.KabschRmsd(a, b, 5, r) < 1e-6);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Assert.Equal(rz[3 * j + i], r[3 * i + j], 9);
        Assert.True(BioKernels.KabschRmsd(a, b, r) < 1e-6);
    }

    [Fact]
    public void Kabsch_AnalyticAndChiralCases()
    {
        Assert.Equal(1.0, ManagedKernels.KabschRmsd([-1, 0, 0, 1, 0, 0], [0, -2, 0, 0, 2, 0], 2, default), 12);
        Assert.True(ManagedKernels.KabschRmsd([0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 3], [0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, -3], 4, default) > 0.1);
        var r = new double[9];
        Assert.Equal(0.0, BioKernels.KabschRmsd([], [], r));
        Assert.Equal([1.0, 0, 0, 0, 1, 0, 0, 0, 1], r);
        Assert.Throws<ArgumentException>(() => BioKernels.KabschRmsd([1, 2, 3], [1, 2]));
    }

    [Fact]
    public void DisableEnvironmentVariable_IsDocumentedName()
        => Assert.Equal("BIOTECH_NET_DISABLE_NATIVE", NativeBridge.DisableEnvironmentVariable);
}
