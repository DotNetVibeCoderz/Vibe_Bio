namespace Biotech.Net.Alignments.Search;

/// <summary>
/// A pluggable scorer for the Smith–Waterman stage of <see cref="SequenceDatabase.Search"/>. It must return exactly the
/// scores of <see cref="Native.BioKernels.LocalScoreBatch"/> (same arguments, same integer results); the GPU backend in
/// <c>Biotech.Net.Gpu</c> implements it. Set it through <see cref="SearchOptions.ScoreBackend"/>.
/// </summary>
public interface ILocalScoreBackend
{
    /// <summary>Short name of the backend, for logs.</summary>
    string Name { get; }

    /// <summary>
    /// Local (affine-gap) scores of <paramref name="query"/> against targets stored back to back in
    /// <paramref name="targets"/> (target k is <c>targets[offsets[k]..offsets[k + 1]]</c>); the query indexes matrix rows.
    /// </summary>
    void LocalScoreBatch(ReadOnlySpan<byte> query, ReadOnlySpan<byte> targets, ReadOnlySpan<int> offsets,
        ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores);
}
