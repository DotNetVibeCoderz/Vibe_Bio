using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments.Search;

/// <summary>
/// Karlin–Altschul statistics for local alignment scores: bit score <c>S' = (λS − ln K) / ln 2</c> and
/// E-value <c>E = K·m·n·e^(−λS)</c> for a query of length m against a database of total length n.
/// </summary>
/// <param name="Lambda">λ, the scale of the score distribution.</param>
/// <param name="K">K, the search-space constant.</param>
/// <param name="H">Relative entropy per aligned pair (nats), when known.</param>
public readonly record struct KarlinAltschulParameters(double Lambda, double K, double H = double.NaN)
{
    /// <summary>BLOSUM62 with gap open 11 / extend 1, the BLASTP default (NCBI blast_stat.c): λ = 0.267, K = 0.041, H = 0.14.</summary>
    public static KarlinAltschulParameters Blosum62Gapped11_1 { get; } = new(0.267, 0.041, 0.14);

    /// <summary>BLOSUM62 ungapped (NCBI blast_stat.c): λ = 0.3176, K = 0.134, H = 0.4012.</summary>
    public static KarlinAltschulParameters Blosum62Ungapped { get; } = new(0.3176, 0.134, 0.4012);

    /// <summary>BLASTN reward 2 / penalty −3 with gap open 5 / extend 2 (NCBI blast_stat.c): λ = 0.625, K = 0.41, H = 0.78.</summary>
    public static KarlinAltschulParameters Blastn2_3Gapped5_2 { get; } = new(0.625, 0.41, 0.78);

    /// <summary>
    /// Published parameters for a matrix name and gap costs (positive open/extend as in BLAST), or null when none are
    /// tabulated. Recognises BLOSUM62 11/1, and <c>match=2,mismatch=-3</c> with 5/2.
    /// </summary>
    public static KarlinAltschulParameters? Lookup(string matrixName, int gapOpen, int gapExtend)
    {
        if (matrixName.Equals("BLOSUM62", StringComparison.OrdinalIgnoreCase) && gapOpen == 11 && gapExtend == 1) return Blosum62Gapped11_1;
        if (matrixName == "match=2,mismatch=-3" && gapOpen == 5 && gapExtend == 2) return Blastn2_3Gapped5_2;
        return null;
    }

    /// <summary>Normalised (bit) score of a raw score.</summary>
    public double BitScore(int rawScore) => (Lambda * rawScore - Math.Log(K)) / Math.Log(2);

    /// <summary>Expected number of chance hits scoring at least <paramref name="rawScore"/> in an m×n search space.</summary>
    public double EValue(int rawScore, long queryLength, long databaseLength)
        => K * queryLength * databaseLength * Math.Exp(-Lambda * rawScore);
}

/// <summary>Options for <see cref="SequenceDatabase.Search"/>.</summary>
public sealed class SearchOptions
{
    /// <summary>Scoring matrix; null = BLOSUM62 for proteins, match 2 / mismatch −3 for nucleotides (BLAST defaults).</summary>
    public SubstitutionMatrix? SubstitutionMatrix { get; set; }

    /// <summary>Gap open cost (positive, BLAST convention: a gap of length k costs open + (k − 1)·extend). Null = 11 protein / 5 nucleotide.</summary>
    public int? GapOpen { get; set; }

    /// <summary>Gap extension cost (positive). Null = 1 protein / 2 nucleotide.</summary>
    public int? GapExtend { get; set; }

    /// <summary>Report hits with an E-value at most this. Default 10 (BLAST default).</summary>
    public double MaxEValue { get; set; } = 10;

    /// <summary>Maximum number of hits returned. Default 50.</summary>
    public int MaxHits { get; set; } = 50;

    /// <summary>When true (default) only subjects sharing at least <see cref="MinSeedHits"/> exact words with the query are scored.</summary>
    public bool UseSeeds { get; set; } = true;

    /// <summary>Distinct shared words needed for a subject to be scored. Default 1.</summary>
    public int MinSeedHits { get; set; } = 1;

    /// <summary>Statistics; null looks up published parameters for the matrix and gap costs (E-values are NaN if none).</summary>
    public KarlinAltschulParameters? Statistics { get; set; }

    /// <summary>Compute a full local alignment (with traceback) for each reported hit. Default true.</summary>
    public bool IncludeAlignments { get; set; } = true;

    /// <summary>
    /// Scorer for the Smith–Waterman stage; null (default) uses <see cref="BioKernels.LocalScoreBatch"/> (Rust or managed).
    /// A GPU backend from <c>Biotech.Net.Gpu</c> gives identical scores.
    /// </summary>
    public ILocalScoreBackend? ScoreBackend { get; set; }
}

/// <summary>A database hit.</summary>
public sealed class SearchHit
{
    internal SearchHit(Sequence subject, int subjectIndex, int score, double bitScore, double eValue, PairwiseAlignment? alignment)
    {
        Subject = subject; SubjectIndex = subjectIndex; Score = score; BitScore = bitScore; EValue = eValue; Alignment = alignment;
    }

    /// <summary>The matching database sequence.</summary>
    public Sequence Subject { get; }

    /// <summary>Index of <see cref="Subject"/> in the database.</summary>
    public int SubjectIndex { get; }

    /// <summary>Raw Smith–Waterman score.</summary>
    public int Score { get; }

    /// <summary>Bit score (NaN without statistics).</summary>
    public double BitScore { get; }

    /// <summary>E-value (NaN without statistics).</summary>
    public double EValue { get; }

    /// <summary>Local alignment, target = subject and query = the search query (null unless requested).</summary>
    public PairwiseAlignment? Alignment { get; }

    /// <inheritdoc />
    public override string ToString() => $"{Subject.Id}\tscore={Score}\tbits={BitScore:F1}\tE={EValue:G3}";
}

/// <summary>
/// A searchable set of sequences: a BLAST-like pipeline with an exact k-mer word index for seeding, gapped
/// Smith–Waterman scoring of every seeded subject through <see cref="BioKernels.LocalScoreBatch"/>
/// (AVX2 striped in the native backend) and Karlin–Altschul E-values over the whole database search space.
/// </summary>
/// <remarks>
/// Unlike BLAST there is no neighbourhood-word expansion or X-drop extension: seeded subjects get an exact
/// Smith–Waterman score, and <see cref="SearchOptions.UseSeeds"/> = false scores every subject (exhaustive SSEARCH-style).
/// E-values use m = query length and n = total residues, without edge-effect correction.
/// </remarks>
public sealed class SequenceDatabase
{
    private readonly List<Sequence> _sequences;
    private readonly Dictionary<long, List<int>> _index = [];

    /// <summary>Indexes <paramref name="sequences"/> with words of <paramref name="wordSize"/> (null = 3 protein / 11 nucleotide; 1–12).</summary>
    public SequenceDatabase(IEnumerable<Sequence> sequences, int? wordSize = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        _sequences = sequences.ToList();
        IsNucleotide = _sequences.Count > 0 && _sequences.All(s => s.Alphabet.IsNucleotide);
        WordSize = wordSize ?? (IsNucleotide ? 11 : 3);
        if (WordSize is < 1 or > 12) throw new ArgumentOutOfRangeException(nameof(wordSize), "Word size must be between 1 and 12.");
        for (int i = 0; i < _sequences.Count; i++)
        {
            TotalLength += _sequences[i].Length;
            foreach (long w in Words(_sequences[i].Span, WordSize).Distinct())
            {
                if (!_index.TryGetValue(w, out var list)) _index[w] = list = [];
                list.Add(i);
            }
        }
    }

    /// <summary>The database sequences.</summary>
    public IReadOnlyList<Sequence> Sequences => _sequences;

    /// <summary>Number of sequences.</summary>
    public int Count => _sequences.Count;

    /// <summary>Total residues (the n of the search space).</summary>
    public long TotalLength { get; }

    /// <summary>Seed word length.</summary>
    public int WordSize { get; }

    /// <summary>True when every sequence is DNA/RNA.</summary>
    public bool IsNucleotide { get; }

    /// <summary>Number of distinct words in the index.</summary>
    public int IndexedWords => _index.Count;

    /// <summary>Words as 5 bits per upper-cased letter (k ≤ 12); windows with non-letters are skipped. U is folded to T.</summary>
    private static IEnumerable<long> Words(ReadOnlySpan<byte> s, int k)
    {
        var words = new List<long>();
        long mask = (1L << (5 * k)) - 1, w = 0;
        int valid = 0;
        foreach (byte b in s)
        {
            int c = b & 0xDF;
            if (c is < 'A' or > 'Z') { valid = 0; w = 0; continue; }
            if (c == 'U') c = 'T';
            w = ((w << 5) | (long)(c - 'A' + 1)) & mask;
            if (++valid >= k) words.Add(w);
        }
        return words;
    }

    /// <summary>Indices of subjects sharing at least <paramref name="minHits"/> distinct words with <paramref name="query"/>.</summary>
    public IReadOnlyList<int> Candidates(Sequence query, int minHits = 1)
    {
        ArgumentNullException.ThrowIfNull(query);
        var hits = new Dictionary<int, int>();
        foreach (long w in Words(query.Span, WordSize).Distinct())
            if (_index.TryGetValue(w, out var list))
                foreach (int i in list) hits[i] = hits.TryGetValue(i, out int c) ? c + 1 : 1;
        return hits.Where(kv => kv.Value >= minHits).Select(kv => kv.Key).Order().ToList();
    }

    /// <summary>Searches the database with <paramref name="query"/>; hits are sorted by score (best first).</summary>
    public IReadOnlyList<SearchHit> Search(Sequence query, SearchOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(query);
        options ??= new SearchOptions();
        bool nuc = IsNucleotide && query.Alphabet.IsNucleotide;
        var matrix = options.SubstitutionMatrix ?? (nuc ? SubstitutionMatrix.MatchMismatch(2, -3) : SubstitutionMatrix.Blosum62);
        int open = options.GapOpen ?? (nuc ? 5 : 11), extend = options.GapExtend ?? (nuc ? 2 : 1);
        var stats = options.Statistics ?? KarlinAltschulParameters.Lookup(matrix.Name, open, extend);

        var candidates = options.UseSeeds ? Candidates(query, Math.Max(1, options.MinSeedHits)) : Enumerable.Range(0, Count).ToList();
        if (candidates.Count == 0) return [];

        long total = candidates.Sum(i => (long)_sequences[i].Length);
        if (total > int.MaxValue) throw new InvalidOperationException("Candidate set exceeds 2 GB; search in chunks.");
        var buffer = new byte[total];
        var offsets = new int[candidates.Count + 1];
        for (int k = 0; k < candidates.Count; k++)
        {
            _sequences[candidates[k]].Span.CopyTo(buffer.AsSpan(offsets[k]));
            offsets[k + 1] = offsets[k] + _sequences[candidates[k]].Length;
        }
        var scores = new int[candidates.Count];
        // Query = matrix rows, as in the kernel's profile; matrices used here are symmetric.
        if (options.ScoreBackend is { } backend)
            backend.LocalScoreBatch(query.Span, buffer, offsets, matrix.KernelMatrix, open, extend, scores);
        else
            BioKernels.LocalScoreBatch(query.Span, buffer, offsets, matrix.KernelMatrix, open, extend, scores);

        var ranked = Enumerable.Range(0, candidates.Count)
            .Select(k => (Index: candidates[k], Score: scores[k],
                E: stats is { } p ? p.EValue(scores[k], query.Length, TotalLength) : double.NaN,
                Bits: stats is { } p2 ? p2.BitScore(scores[k]) : double.NaN))
            .Where(h => h.Score > 0 && (double.IsNaN(h.E) || h.E <= options.MaxEValue))
            .OrderByDescending(h => h.Score).ThenBy(h => h.Index)
            .Take(Math.Max(0, options.MaxHits))
            .ToList();

        var aligner = new PairwiseAligner
        {
            Mode = AlignmentMode.Local, SubstitutionMatrix = matrix, OpenGapScore = -open, ExtendGapScore = -extend,
        };
        return ranked.Select(h => new SearchHit(_sequences[h.Index], h.Index, h.Score, h.Bits, h.E,
            options.IncludeAlignments ? aligner.Align(_sequences[h.Index], query) : null)).ToList();
    }
}
