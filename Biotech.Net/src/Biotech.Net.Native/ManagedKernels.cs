using System.Numerics;

namespace Biotech.Net.Native;

/// <summary>
/// Pure C# implementations of every kernel. They are the fallback when the Rust library is missing and
/// the reference the native kernels are tested against, so they favour clarity over the last few percent.
/// </summary>
internal static partial class ManagedKernels
{
    /// <summary>IUPAC complement table, case-preserving; anything unrecognised maps to itself.</summary>
    internal static readonly byte[] ComplementTable = BuildComplementTable();

    /// <summary>2-bit code for A/C/G/T(U), case-insensitive; 255 for anything else.</summary>
    internal static readonly byte[] TwoBitCode = BuildTwoBitCode();

    private static byte[] BuildComplementTable()
    {
        var t = new byte[256];
        for (int i = 0; i < 256; i++) t[i] = (byte)i;
        void Pair(char a, char b)
        {
            t[a] = (byte)b; t[b] = (byte)a;
            t[char.ToLowerInvariant(a)] = (byte)char.ToLowerInvariant(b);
            t[char.ToLowerInvariant(b)] = (byte)char.ToLowerInvariant(a);
        }
        Pair('A', 'T'); Pair('C', 'G'); Pair('R', 'Y'); Pair('K', 'M'); Pair('B', 'V'); Pair('D', 'H');
        t['U'] = (byte)'A'; t['u'] = (byte)'a';
        // S, W and N are self-complementary and already map to themselves.
        return t;
    }

    private static byte[] BuildTwoBitCode()
    {
        var t = new byte[256];
        Array.Fill(t, (byte)255);
        t['A'] = t['a'] = 0;
        t['C'] = t['c'] = 1;
        t['G'] = t['g'] = 2;
        t['T'] = t['t'] = t['U'] = t['u'] = 3;
        return t;
    }

    public static void ReverseComplement(ReadOnlySpan<byte> src, Span<byte> dst)
    {
        var table = ComplementTable;
        int n = src.Length;
        for (int i = 0; i < n; i++) dst[n - 1 - i] = table[src[i]];
    }

    public static void Complement(ReadOnlySpan<byte> src, Span<byte> dst)
    {
        var table = ComplementTable;
        for (int i = 0; i < src.Length; i++) dst[i] = table[src[i]];
    }

    public static int Translate(ReadOnlySpan<byte> nucleotides, ReadOnlySpan<byte> codonTable, Span<byte> protein)
    {
        var code = TwoBitCode;
        int codons = nucleotides.Length / 3;
        for (int i = 0, j = 0; i < codons; i++, j += 3)
        {
            int a = code[nucleotides[j]], b = code[nucleotides[j + 1]], c = code[nucleotides[j + 2]];
            protein[i] = (a | b | c) > 3 ? (byte)'X' : codonTable[(a << 4) | (b << 2) | c];
        }
        return codons;
    }

    public static NucleotideCounts CountNucleotides(ReadOnlySpan<byte> seq)
    {
        Span<long> h = stackalloc long[256];
        h.Clear();
        foreach (byte b in seq) h[b]++;
        long a = h['A'] + h['a'], c = h['C'] + h['c'], g = h['G'] + h['g'];
        long t = h['T'] + h['t'] + h['U'] + h['u'], n = h['N'] + h['n'];
        return new NucleotideCounts(a, c, g, t, n, seq.Length - a - c - g - t - n);
    }

    public static void ByteHistogram(ReadOnlySpan<byte> data, Span<long> counts)
    {
        // Four interleaved histograms break the store-to-load dependency on runs of the same byte.
        Span<long> h = stackalloc long[1024];
        h.Clear();
        int i = 0;
        for (; i + 4 <= data.Length; i += 4)
        {
            h[data[i]]++; h[256 + data[i + 1]]++; h[512 + data[i + 2]]++; h[768 + data[i + 3]]++;
        }
        for (; i < data.Length; i++) h[data[i]]++;
        for (int k = 0; k < 256; k++) counts[k] += h[k] + h[256 + k] + h[512 + k] + h[768 + k];
    }

    public static long CountKmers(ReadOnlySpan<byte> seq, int k, Span<int> counts)
    {
        var code = TwoBitCode;
        uint mask = k == 16 ? uint.MaxValue : (1u << (2 * k)) - 1;
        uint kmer = 0;
        int valid = 0;
        long total = 0;
        foreach (byte b in seq)
        {
            byte c = code[b];
            if (c > 3) { valid = 0; kmer = 0; continue; }
            kmer = ((kmer << 2) | c) & mask;
            if (++valid >= k) { counts[(int)kmer]++; total++; }
        }
        return total;
    }

    public static int Hamming(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b)
    {
        int d = 0, i = 0;
        if (Vector.IsHardwareAccelerated)
        {
            int w = Vector<byte>.Count;
            for (; i + w <= a.Length; i += w)
            {
                // Equal lanes are 0xFF: eight set bits per matching byte.
                var eq = Vector.AsVectorUInt64(Vector.Equals(new Vector<byte>(a.Slice(i, w)), new Vector<byte>(b.Slice(i, w))));
                int bits = 0;
                for (int k = 0; k < Vector<ulong>.Count; k++) bits += BitOperations.PopCount(eq[k]);
                d += w - bits / 8;
            }
        }
        for (; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    public static void WelchTTest(ReadOnlySpan<double> matrix, int rows, int cols,
        ReadOnlySpan<int> groupA, ReadOnlySpan<int> groupB, Span<double> t, Span<double> p, Span<double> meanDiff)
    {
        for (int r = 0; r < rows; r++)
        {
            var row = matrix.Slice(r * cols, cols);
            MeanVar(row, groupA, out double ma, out double va);
            MeanVar(row, groupB, out double mb, out double vb);
            double na = groupA.Length, nb = groupB.Length;
            double sa = va / na, sb = vb / nb, se2 = sa + sb;
            meanDiff[r] = mb - ma;
            if (se2 <= 0)
            {
                t[r] = ma == mb ? 0 : double.PositiveInfinity * Math.Sign(mb - ma);
                p[r] = ma == mb ? 1 : 0;
                continue;
            }
            double tv = (mb - ma) / Math.Sqrt(se2);
            double df = se2 * se2 / (sa * sa / (na - 1) + sb * sb / (nb - 1));
            t[r] = tv;
            p[r] = SpecialFunctions.StudentTTwoSided(tv, df);
        }
    }

    private static void MeanVar(ReadOnlySpan<double> row, ReadOnlySpan<int> idx, out double mean, out double variance)
    {
        double s = 0;
        foreach (int i in idx) s += row[i];
        mean = s / idx.Length;
        double ss = 0;
        foreach (int i in idx) { double d = row[i] - mean; ss += d * d; }
        variance = idx.Length > 1 ? ss / (idx.Length - 1) : 0;
    }

    public static void BenjaminiHochberg(ReadOnlySpan<double> p, Span<double> q)
    {
        int n = p.Length;
        if (n == 0) return;
        var order = new int[n];
        var keys = new double[n];
        for (int i = 0; i < n; i++) { order[i] = i; keys[i] = double.IsNaN(p[i]) ? double.PositiveInfinity : p[i]; }
        Array.Sort(keys, order);
        // NaN p-values stay NaN and do not count towards m (they sort last, after every real p-value).
        int m = 0;
        for (int i = 0; i < n; i++) { if (double.IsNaN(p[i])) q[i] = double.NaN; else m++; }
        double running = 1;
        for (int rank = m; rank >= 1; rank--)
        {
            int i = order[rank - 1];
            running = Math.Min(running, p[i] * (double)m / rank);
            q[i] = Math.Min(running, 1);
        }
    }

    // ---- Protein engine -------------------------------------------------------------------------

    internal static readonly double[] AverageResidueMass = BuildMassTable(monoisotopic: false);
    internal static readonly double[] MonoResidueMass = BuildMassTable(monoisotopic: true);
    internal const double WaterAverage = 18.01528, WaterMono = 18.0105646863;

    private static double[] BuildMassTable(bool monoisotopic)
    {
        var t = new double[256];
        // Residue (in-chain, water removed) masses: (average, monoisotopic).
        (char, double, double)[] residues =
        [
            ('A', 71.0788, 71.03711), ('R', 156.1875, 156.10111), ('N', 114.1038, 114.04293), ('D', 115.0886, 115.02694),
            ('C', 103.1388, 103.00919), ('E', 129.1155, 129.04259), ('Q', 128.1307, 128.05858), ('G', 57.0519, 57.02146),
            ('H', 137.1411, 137.05891), ('I', 113.1594, 113.08406), ('L', 113.1594, 113.08406), ('K', 128.1741, 128.09496),
            ('M', 131.1926, 131.04049), ('F', 147.1766, 147.06841), ('P', 97.1167, 97.05276), ('S', 87.0782, 87.03203),
            ('T', 101.1051, 101.04768), ('W', 186.2132, 186.07931), ('Y', 163.1760, 163.06333), ('V', 99.1326, 99.06841),
            ('U', 150.0388, 150.95364), ('O', 237.3018, 237.14773),
            // Ambiguity codes take the mean of their candidates.
            ('B', 114.5962, 114.53494), ('Z', 128.6231, 128.55059), ('J', 113.1594, 113.08406), ('X', 110.0, 110.0),
        ];
        foreach (var (c, avg, mono) in residues)
        {
            double m = monoisotopic ? mono : avg;
            t[c] = m; t[char.ToLowerInvariant(c)] = m;
        }
        return t;
    }

    public static double ProteinMass(ReadOnlySpan<byte> seq, bool monoisotopic)
    {
        if (seq.IsEmpty) return 0;
        var table = monoisotopic ? MonoResidueMass : AverageResidueMass;
        double m = monoisotopic ? WaterMono : WaterAverage;
        foreach (byte b in seq) m += table[b];
        return m;
    }

    // EMBOSS pKa set.
    private const double PkNTerm = 8.6, PkCTerm = 3.6, PkK = 10.8, PkR = 12.5, PkH = 6.5,
        PkD = 3.9, PkE = 4.1, PkC = 8.5, PkY = 10.1;

    public static double NetCharge(ReadOnlySpan<byte> seq, double pH)
    {
        Span<long> h = stackalloc long[256];
        h.Clear();
        foreach (byte b in seq) h[b & 0xDF]++; // fold to upper case
        return ChargeFromCounts(h, pH);
    }

    private static double ChargeFromCounts(ReadOnlySpan<long> h, double pH)
    {
        static double Pos(double pk, double ph) => 1 / (1 + Math.Pow(10, ph - pk));
        static double Neg(double pk, double ph) => -1 / (1 + Math.Pow(10, pk - ph));
        return Pos(PkNTerm, pH) + Neg(PkCTerm, pH)
            + h['K'] * Pos(PkK, pH) + h['R'] * Pos(PkR, pH) + h['H'] * Pos(PkH, pH)
            + h['D'] * Neg(PkD, pH) + h['E'] * Neg(PkE, pH) + h['C'] * Neg(PkC, pH) + h['Y'] * Neg(PkY, pH);
    }

    public static double IsoelectricPoint(ReadOnlySpan<byte> seq)
    {
        Span<long> h = stackalloc long[256];
        h.Clear();
        foreach (byte b in seq) h[b & 0xDF]++;
        double lo = 0, hi = 14;
        for (int i = 0; i < 100 && hi - lo > 1e-4; i++)
        {
            double mid = (lo + hi) / 2;
            if (ChargeFromCounts(h, mid) > 0) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }

    internal static readonly double[] KyteDoolittle = BuildKyteDoolittle();

    private static double[] BuildKyteDoolittle()
    {
        var t = new double[256];
        (char, double)[] v =
        [
            ('A', 1.8), ('R', -4.5), ('N', -3.5), ('D', -3.5), ('C', 2.5), ('Q', -3.5), ('E', -3.5), ('G', -0.4),
            ('H', -3.2), ('I', 4.5), ('L', 3.8), ('K', -3.9), ('M', 1.9), ('F', 2.8), ('P', -1.6), ('S', -0.8),
            ('T', -0.7), ('W', -0.9), ('Y', -1.3), ('V', 4.2),
        ];
        foreach (var (c, x) in v) { t[c] = x; t[char.ToLowerInvariant(c)] = x; }
        return t;
    }

    public static int Hydropathy(ReadOnlySpan<byte> seq, int window, Span<double> output)
    {
        int n = seq.Length - window + 1;
        if (n <= 0) return 0;
        var kd = KyteDoolittle;
        double s = 0;
        for (int i = 0; i < window; i++) s += kd[seq[i]];
        output[0] = s / window;
        for (int i = 1; i < n; i++)
        {
            s += kd[seq[i + window - 1]] - kd[seq[i - 1]];
            output[i] = s / window;
        }
        return n;
    }
}
