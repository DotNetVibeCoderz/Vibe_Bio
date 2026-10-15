using System.Text;

namespace Biotech.Net.Native;

/// <summary>Which end gaps a kernel alignment scores.</summary>
public enum AlignmentMode
{
    /// <summary>Needleman–Wunsch: every gap is scored, both sequences are aligned end to end.</summary>
    Global = 0,

    /// <summary>Smith–Waterman: the best-scoring pair of substrings; scores never drop below zero.</summary>
    Local = 1,

    /// <summary>Global alignment with free leading and trailing gaps (overlap / glocal alignment).</summary>
    SemiGlobal = 2,
}

/// <summary>
/// Score and 0-based half-open coordinates of an alignment whose operations were written to a caller buffer.
/// For global and semi-global alignments the coordinates span both whole sequences.
/// </summary>
/// <param name="Score">Alignment score.</param>
/// <param name="AStart">First aligned position of sequence a.</param>
/// <param name="AEnd">End (exclusive) of the aligned region of sequence a.</param>
/// <param name="BStart">First aligned position of sequence b.</param>
/// <param name="BEnd">End (exclusive) of the aligned region of sequence b.</param>
/// <param name="OperationCount">Number of operations written.</param>
public readonly record struct AlignmentSummary(int Score, int AStart, int AEnd, int BStart, int BEnd, int OperationCount);

/// <summary>
/// A kernel alignment: the summary plus its edit operations — <c>M</c> (a[i] aligned to b[j]), <c>I</c> (b residue
/// against a gap) and <c>D</c> (a residue against a gap).
/// </summary>
public sealed class KernelAlignment
{
    internal KernelAlignment(AlignmentSummary summary, byte[] operations)
    {
        Summary = summary;
        Operations = operations;
    }

    /// <summary>Score and coordinates.</summary>
    public AlignmentSummary Summary { get; }

    /// <summary>Alignment score.</summary>
    public int Score => Summary.Score;

    /// <summary>Operations as ASCII bytes (<c>M</c>/<c>I</c>/<c>D</c>).</summary>
    public ReadOnlyMemory<byte> Operations { get; }

    /// <summary>Operations as a string, e.g. <c>MMMDDMM</c>.</summary>
    public string OperationString => Encoding.ASCII.GetString(Operations.Span);

    /// <inheritdoc />
    public override string ToString() => $"score={Score} a=[{Summary.AStart},{Summary.AEnd}) b=[{Summary.BStart},{Summary.BEnd}) {OperationString}";
}

public static partial class BioKernels
{
    /// <summary>
    /// Entries of a kernel substitution matrix: 128×128 <see cref="int"/> scores indexed by ASCII byte,
    /// <c>matrix[(a &amp; 0x7F) &lt;&lt; 7 | (b &amp; 0x7F)]</c>, rows for sequence a. Case handling is up to the matrix.
    /// </summary>
    public const int SubstitutionMatrixLength = 128 * 128;

    /// <summary>
    /// Builds a kernel matrix scoring <paramref name="match"/> for identical letters (case-insensitive) and
    /// <paramref name="mismatch"/> otherwise.
    /// </summary>
    public static int[] MatchMismatchMatrix(int match, int mismatch)
    {
        var m = new int[SubstitutionMatrixLength];
        for (int a = 0; a < 128; a++)
            for (int b = 0; b < 128; b++)
                m[(a << 7) | b] = a == b || (char.IsAsciiLetter((char)a) && (a | 0x20) == (b | 0x20)) ? match : mismatch;
        return m;
    }

    /// <summary>
    /// Affine-gap (Gotoh) alignment with traceback. A gap of length k costs <c>gapOpen + (k − 1)·gapExtend</c>
    /// (EMBOSS/Biopython convention; both non-negative). <paramref name="band"/> &lt; 0 fills the whole matrix; otherwise
    /// only cells within <paramref name="band"/> diagonals of the main band (widened by the length difference) are used.
    /// Operations are written to <paramref name="operations"/>, which needs <c>a.Length + b.Length</c> bytes.
    /// </summary>
    /// <remarks>Ties prefer the diagonal, then I, then D; gaps prefer opening over extending; local and semi-global keep the first maximum.</remarks>
    public static AlignmentSummary Align(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, ReadOnlySpan<int> matrix,
        int gapOpen, int gapExtend, AlignmentMode mode, Span<byte> operations, int band = -1)
    {
        ValidateScoring(matrix, gapOpen, gapExtend);
        if (mode is < AlignmentMode.Global or > AlignmentMode.SemiGlobal) throw new ArgumentOutOfRangeException(nameof(mode));
        RequireLength(operations, a.Length + b.Length);
        Span<int> r = stackalloc int[6];
        int status = UseNative
            ? NativeBridge.Align(a, b, matrix, gapOpen, gapExtend, (int)mode, band, operations, r)
            : ManagedKernels.Align(a, b, matrix, gapOpen, gapExtend, (int)mode, band, operations, r);
        ThrowOnError(status);
        return new AlignmentSummary(r[0], r[1], r[2], r[3], r[4], r[5]);
    }

    /// <summary>Allocating overload of <see cref="Align(ReadOnlySpan{byte}, ReadOnlySpan{byte}, ReadOnlySpan{int}, int, int, AlignmentMode, Span{byte}, int)"/>.</summary>
    public static KernelAlignment Align(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, ReadOnlySpan<int> matrix,
        int gapOpen, int gapExtend, AlignmentMode mode = AlignmentMode.Global, int band = -1)
    {
        var ops = new byte[a.Length + b.Length];
        var s = Align(a, b, matrix, gapOpen, gapExtend, mode, ops, band);
        return new KernelAlignment(s, ops.AsSpan(0, s.OperationCount).ToArray());
    }

    /// <summary>Best local (Smith–Waterman, affine gaps) score of <paramref name="query"/> against <paramref name="target"/> in linear memory.</summary>
    public static int LocalScore(ReadOnlySpan<byte> query, ReadOnlySpan<byte> target, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend)
    {
        ValidateScoring(matrix, gapOpen, gapExtend);
        int s = UseNative
            ? NativeBridge.LocalScore(query, target, matrix, gapOpen, gapExtend)
            : ManagedKernels.LocalScore(query, target, matrix, gapOpen, gapExtend);
        ThrowOnError(s);
        return s;
    }

    /// <summary>
    /// Local scores of one query against many targets stored back to back in <paramref name="targets"/>:
    /// target k is <c>targets[offsets[k]..offsets[k + 1]]</c>, so <paramref name="offsets"/> holds count + 1 ascending positions.
    /// The query profile is built once (AVX2-striped in the native backend).
    /// </summary>
    public static void LocalScoreBatch(ReadOnlySpan<byte> query, ReadOnlySpan<byte> targets, ReadOnlySpan<int> offsets,
        ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
    {
        ValidateScoring(matrix, gapOpen, gapExtend);
        if (offsets.IsEmpty) throw new ArgumentException("Offsets need at least one entry (count + 1).", nameof(offsets));
        for (int k = 0; k + 1 < offsets.Length; k++)
            if (offsets[k] < 0 || offsets[k + 1] < offsets[k] || offsets[k + 1] > targets.Length)
                throw new ArgumentException("Offsets must be ascending positions inside the target buffer.", nameof(offsets));
        RequireLength(scores, offsets.Length - 1);
        int status = UseNative
            ? NativeBridge.LocalScoreBatch(query, targets, offsets, matrix, gapOpen, gapExtend, scores)
            : ManagedKernels.LocalScoreBatch(query, targets, offsets, matrix, gapOpen, gapExtend, scores);
        ThrowOnError(status);
    }

    // ---- Structure engine -----------------------------------------------------------------------

    /// <summary>
    /// RMSD after optimally superposing <paramref name="mobile"/> onto <paramref name="target"/> (Horn's quaternion
    /// method; proper rotations only). Both are row-major xyz triples of equal length. When <paramref name="rotation"/>
    /// has 9 elements it receives the row-major 3×3 R with <c>R·(mobile − mobileCentroid) + targetCentroid ≈ target</c>.
    /// </summary>
    public static double KabschRmsd(ReadOnlySpan<double> target, ReadOnlySpan<double> mobile, Span<double> rotation = default)
    {
        if (target.Length != mobile.Length || target.Length % 3 != 0)
            throw new ArgumentException("Coordinate sets must have the same length, a multiple of 3 (x, y, z per point).");
        if (!rotation.IsEmpty && rotation.Length < 9) throw new ArgumentException("A rotation needs 9 elements.", nameof(rotation));
        if (rotation.Length > 9) rotation = rotation[..9];
        int n = target.Length / 3;
        return UseNative ? NativeBridge.KabschRmsd(target, mobile, n, rotation) : ManagedKernels.KabschRmsd(target, mobile, n, rotation);
    }

    private static void ValidateScoring(ReadOnlySpan<int> matrix, int gapOpen, int gapExtend)
    {
        if (matrix.Length < SubstitutionMatrixLength)
            throw new ArgumentException($"A kernel substitution matrix has {SubstitutionMatrixLength} entries (128×128).", nameof(matrix));
        if (gapOpen < 0) throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap penalties are non-negative costs.");
        if (gapExtend < 0) throw new ArgumentOutOfRangeException(nameof(gapExtend), "Gap penalties are non-negative costs.");
    }

    private static void ThrowOnError(int status)
    {
        switch (status)
        {
            case >= 0: return;
            case ManagedKernels.ErrOutOfMemory: throw new OutOfMemoryException("The alignment matrix does not fit in memory; use a band or a score-only kernel.");
            case ManagedKernels.ErrBufferTooSmall: throw new ArgumentException("Operation buffer too small.");
            case ManagedKernels.ErrBadArgument: throw new ArgumentException("Invalid alignment argument.");
            default: throw new InvalidOperationException($"Native kernel failed with status {status}.");
        }
    }
}
