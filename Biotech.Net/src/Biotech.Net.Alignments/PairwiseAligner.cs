using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments;

/// <summary>
/// Pairwise sequence aligner modelled on Biopython's <c>Bio.Align.PairwiseAligner</c>: global (Needleman–Wunsch),
/// local (Smith–Waterman) or semi-global alignment with affine gaps, running on the Rust/managed
/// <see cref="BioKernels.Align(ReadOnlySpan{byte}, ReadOnlySpan{byte}, ReadOnlySpan{int}, int, int, AlignmentMode, int)"/> kernel.
/// </summary>
/// <remarks>
/// Gap scores follow Biopython: a gap of length k scores <c>OpenGapScore + (k − 1)·ExtendGapScore</c>, and both are
/// zero or negative. Scores are integers. Semi-global mode is Biopython's global mode with end gap scores of zero.
/// <code>
/// var aligner = new PairwiseAligner { Mode = AlignmentMode.Local, SubstitutionMatrix = SubstitutionMatrix.Blosum62,
///                                     OpenGapScore = -11, ExtendGapScore = -1 };
/// PairwiseAlignment aln = aligner.Align(target, query);
/// Console.WriteLine(aln);   // Biopython-style three-line blocks
/// </code>
/// </remarks>
public sealed class PairwiseAligner
{
    private int _openGapScore, _extendGapScore;

    /// <summary>Global (default), local or semi-global.</summary>
    public AlignmentMode Mode { get; set; } = AlignmentMode.Global;

    /// <summary>Substitution matrix; when null, <see cref="MatchScore"/> and <see cref="MismatchScore"/> are used.</summary>
    public SubstitutionMatrix? SubstitutionMatrix { get; set; }

    /// <summary>Score for identical residues (case-insensitive) when no matrix is set. Default 1, as in Biopython.</summary>
    public int MatchScore { get; set; } = 1;

    /// <summary>Score for different residues when no matrix is set. Default 0, as in Biopython.</summary>
    public int MismatchScore { get; set; }

    /// <summary>Score of the first position of a gap (≤ 0). Default 0.</summary>
    public int OpenGapScore
    {
        get => _openGapScore;
        set => _openGapScore = value <= 0 ? value : throw new ArgumentOutOfRangeException(nameof(value), "Gap scores are zero or negative.");
    }

    /// <summary>Score of each further position of a gap (≤ 0). Default 0.</summary>
    public int ExtendGapScore
    {
        get => _extendGapScore;
        set => _extendGapScore = value <= 0 ? value : throw new ArgumentOutOfRangeException(nameof(value), "Gap scores are zero or negative.");
    }

    /// <summary>Sets <see cref="OpenGapScore"/> and <see cref="ExtendGapScore"/> to the same value (a linear gap score).</summary>
    public int GapScore
    {
        set { OpenGapScore = value; ExtendGapScore = value; }
    }

    /// <summary>Diagonal band half-width for banded alignment; negative (default) fills the whole matrix.</summary>
    public int Band { get; set; } = -1;

    private int[] KernelMatrix => SubstitutionMatrix?.KernelArray ?? BioKernels.MatchMismatchMatrix(MatchScore, MismatchScore);

    /// <summary>Aligns <paramref name="query"/> to <paramref name="target"/> and returns one optimal alignment.</summary>
    public PairwiseAlignment Align(Sequence target, Sequence query)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(query);
        var k = BioKernels.Align(target.Span, query.Span, KernelMatrix, -OpenGapScore, -ExtendGapScore, Mode, Band);
        return new PairwiseAlignment(target, query, k.Summary, k.OperationString, Mode, SubstitutionMatrix);
    }

    /// <summary>Aligns two residue strings (type detected from the residues).</summary>
    public PairwiseAlignment Align(string target, string query)
        => Align(Sequence.Parse(target, "target"), Sequence.Parse(query, "query"));

    /// <summary>Optimal score only. Local mode uses the linear-memory (AVX2 striped when native) kernel.</summary>
    public int Score(Sequence target, Sequence query)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(query);
        if (Mode == AlignmentMode.Local && Band < 0)
            return BioKernels.LocalScore(target.Span, query.Span, KernelMatrix, -OpenGapScore, -ExtendGapScore);
        return Align(target, query).Score;
    }

    /// <summary>Optimal score of two residue strings.</summary>
    public int Score(string target, string query) => Score(Sequence.Parse(target, "target"), Sequence.Parse(query, "query"));

    /// <summary>Kernel gap costs (positive open/extend) for internal use by search and MSA.</summary>
    internal (int Open, int Extend) KernelGaps => (-OpenGapScore, -ExtendGapScore);

    internal int[] KernelMatrixArray => KernelMatrix;

    /// <inheritdoc />
    public override string ToString()
        => $"PairwiseAligner(mode={Mode}, matrix={SubstitutionMatrix?.Name ?? $"match={MatchScore}/mismatch={MismatchScore}"}, open={OpenGapScore}, extend={ExtendGapScore})";
}
