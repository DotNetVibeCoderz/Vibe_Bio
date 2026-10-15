namespace Biotech.Net.Native;

/// <summary>Which implementation <see cref="BioKernels"/> dispatches to.</summary>
public enum KernelBackend
{
    /// <summary>The Rust library (biotech_native) is loaded and used for every kernel that has a native version.</summary>
    Native,

    /// <summary>Pure C# implementations only.</summary>
    Managed,
}

/// <summary>
/// The performance-critical primitives every Biotech.Net module is built on: the sequence engine
/// (complement, translation, composition, k-mers), the statistics engine (row-wise t-tests, FDR) and the
/// protein engine (mass, charge, pI, hydropathy). Each call goes to the Rust library when it is present
/// and to an equivalent managed implementation otherwise; both produce identical results.
/// </summary>
/// <remarks>
/// All entry points take spans so callers can pass slices of pooled or memory-mapped buffers without
/// copying. Nucleotide and amino-acid inputs are ASCII bytes and are treated case-insensitively.
/// </remarks>
public static partial class BioKernels
{
    /// <summary>Forces the managed implementation even when the native library is available (tests, benchmarks).</summary>
    public static bool ForceManaged { get; set; }

    /// <summary>The backend calls are currently dispatched to.</summary>
    public static KernelBackend Backend => UseNative ? KernelBackend.Native : KernelBackend.Managed;

    /// <summary>True when the Rust library was found and loaded.</summary>
    public static bool IsNativeAvailable => NativeBridge.IsAvailable;

    /// <summary>Version string reported by the native library, or <c>"managed"</c>.</summary>
    public static string NativeVersion => NativeBridge.IsAvailable ? NativeBridge.Version : "managed";

    private static bool UseNative => !ForceManaged && NativeBridge.IsAvailable;

    // ---- Sequence engine ------------------------------------------------------------------------

    /// <summary>Writes the IUPAC reverse complement of <paramref name="source"/> into <paramref name="destination"/> (case preserved, U → A).</summary>
    public static void ReverseComplement(ReadOnlySpan<byte> source, Span<byte> destination)
    {
        RequireLength(destination, source.Length);
        if (UseNative) NativeBridge.ReverseComplement(source, destination);
        else ManagedKernels.ReverseComplement(source, destination);
    }

    /// <summary>Writes the IUPAC complement (not reversed) of <paramref name="source"/>.</summary>
    public static void Complement(ReadOnlySpan<byte> source, Span<byte> destination)
    {
        RequireLength(destination, source.Length);
        if (UseNative) NativeBridge.Complement(source, destination);
        else ManagedKernels.Complement(source, destination);
    }

    /// <summary>
    /// Translates whole codons of <paramref name="nucleotides"/> using a 64-entry codon table indexed by
    /// <c>16·b1 + 4·b2 + b3</c> with A=0, C=1, G=2, T/U=3. Codons containing any other symbol become <c>X</c>.
    /// Returns the number of residues written (<c>nucleotides.Length / 3</c>).
    /// </summary>
    public static int Translate(ReadOnlySpan<byte> nucleotides, ReadOnlySpan<byte> codonTable, Span<byte> protein)
    {
        if (codonTable.Length != 64) throw new ArgumentException("A codon table has exactly 64 entries.", nameof(codonTable));
        RequireLength(protein, nucleotides.Length / 3);
        return UseNative
            ? NativeBridge.Translate(nucleotides, codonTable, protein)
            : ManagedKernels.Translate(nucleotides, codonTable, protein);
    }

    /// <summary>Counts A, C, G, T/U, N and other symbols.</summary>
    public static NucleotideCounts CountNucleotides(ReadOnlySpan<byte> sequence)
        => UseNative ? NativeBridge.CountNucleotides(sequence) : ManagedKernels.CountNucleotides(sequence);

    /// <summary>Adds the occurrences of every byte value in <paramref name="data"/> to <paramref name="counts"/> (length 256).</summary>
    public static void ByteHistogram(ReadOnlySpan<byte> data, Span<long> counts)
    {
        RequireLength(counts, 256);
        if (UseNative) NativeBridge.ByteHistogram(data, counts);
        else ManagedKernels.ByteHistogram(data, counts);
    }

    /// <summary>
    /// Adds every k-mer of <paramref name="sequence"/> (1 ≤ k ≤ 12) to <paramref name="counts"/>, indexed by its
    /// 2-bit encoding (A=0, C=1, G=2, T=3, most significant base first); <paramref name="counts"/> needs 4^k entries.
    /// Windows containing a non-ACGT symbol are skipped. Returns the number of k-mers counted.
    /// </summary>
    public static long CountKmers(ReadOnlySpan<byte> sequence, int k, Span<int> counts)
    {
        if (k is < 1 or > 12) throw new ArgumentOutOfRangeException(nameof(k), "k must be between 1 and 12 for dense counting.");
        RequireLength(counts, 1 << (2 * k));
        return UseNative ? NativeBridge.CountKmers(sequence, k, counts) : ManagedKernels.CountKmers(sequence, k, counts);
    }

    /// <summary>Number of positions at which two equal-length sequences differ (byte-exact).</summary>
    public static int HammingDistance(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Hamming distance needs sequences of equal length.");
        return UseNative ? NativeBridge.Hamming(a, b) : ManagedKernels.Hamming(a, b);
    }

    // ---- Statistics engine ----------------------------------------------------------------------

    /// <summary>
    /// Welch's two-sample t-test on every row of a row-major <paramref name="rows"/>×<paramref name="cols"/> matrix,
    /// comparing the columns in <paramref name="groupA"/> with those in <paramref name="groupB"/>. Writes the t statistic,
    /// the two-sided p-value and mean(B) − mean(A) for each row.
    /// </summary>
    public static void WelchTTest(ReadOnlySpan<double> matrix, int rows, int cols, ReadOnlySpan<int> groupA,
        ReadOnlySpan<int> groupB, Span<double> t, Span<double> p, Span<double> meanDifference)
    {
        if (matrix.Length < (long)rows * cols) throw new ArgumentException("Matrix is smaller than rows × cols.", nameof(matrix));
        if (groupA.Length < 2 || groupB.Length < 2) throw new ArgumentException("Each group needs at least two samples.");
        foreach (int c in groupA) if ((uint)c >= (uint)cols) throw new ArgumentOutOfRangeException(nameof(groupA));
        foreach (int c in groupB) if ((uint)c >= (uint)cols) throw new ArgumentOutOfRangeException(nameof(groupB));
        RequireLength(t, rows); RequireLength(p, rows); RequireLength(meanDifference, rows);
        if (UseNative) NativeBridge.WelchTTest(matrix, rows, cols, groupA, groupB, t, p, meanDifference);
        else ManagedKernels.WelchTTest(matrix, rows, cols, groupA, groupB, t, p, meanDifference);
    }

    /// <summary>Benjamini–Hochberg adjusted p-values (q-values). NaN inputs stay NaN and are excluded from the ranking count.</summary>
    public static void BenjaminiHochberg(ReadOnlySpan<double> pValues, Span<double> qValues)
    {
        RequireLength(qValues, pValues.Length);
        if (UseNative) NativeBridge.BenjaminiHochberg(pValues, qValues);
        else ManagedKernels.BenjaminiHochberg(pValues, qValues);
    }

    // ---- Protein engine -------------------------------------------------------------------------

    /// <summary>Mass of a peptide/protein in daltons (average or monoisotopic), including one water.</summary>
    public static double ProteinMass(ReadOnlySpan<byte> sequence, bool monoisotopic = false)
        => UseNative ? NativeBridge.ProteinMass(sequence, monoisotopic) : ManagedKernels.ProteinMass(sequence, monoisotopic);

    /// <summary>Net charge at <paramref name="pH"/> using the EMBOSS pKa set.</summary>
    public static double NetCharge(ReadOnlySpan<byte> sequence, double pH)
        => UseNative ? NativeBridge.NetCharge(sequence, pH) : ManagedKernels.NetCharge(sequence, pH);

    /// <summary>Isoelectric point (pH of zero net charge) by bisection, EMBOSS pKa set.</summary>
    public static double IsoelectricPoint(ReadOnlySpan<byte> sequence)
        => UseNative ? NativeBridge.IsoelectricPoint(sequence) : ManagedKernels.IsoelectricPoint(sequence);

    /// <summary>
    /// Kyte–Doolittle hydropathy averaged over a sliding <paramref name="window"/>; writes
    /// <c>sequence.Length − window + 1</c> values (0 when the sequence is shorter than the window) and returns that count.
    /// </summary>
    public static int Hydropathy(ReadOnlySpan<byte> sequence, int window, Span<double> output)
    {
        if (window < 1) throw new ArgumentOutOfRangeException(nameof(window));
        RequireLength(output, Math.Max(0, sequence.Length - window + 1));
        return UseNative ? NativeBridge.Hydropathy(sequence, window, output) : ManagedKernels.Hydropathy(sequence, window, output);
    }

    private static void RequireLength<T>(Span<T> span, int length)
    {
        if (span.Length < length)
            throw new ArgumentException($"Destination needs at least {length} elements but has {span.Length}.");
    }
}
