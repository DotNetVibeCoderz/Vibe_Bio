namespace Biotech.Net.AI.Preprocessing;

/// <summary>Settings for <see cref="GeneformerEncoder"/>.</summary>
public sealed record GeneformerEncoderOptions
{
    /// <summary>Settings matching Geneformer V1 (30M-gene corpus): 2,048 tokens, no special tokens.</summary>
    public static GeneformerEncoderOptions V1 { get; } = new() { MaxLength = 2048, AddSpecialTokens = false };

    /// <summary>Settings matching Geneformer V2 (104M-gene corpus): 4,096 tokens framed by &lt;cls&gt; and &lt;eos&gt;.</summary>
    public static GeneformerEncoderOptions V2 { get; } = new() { MaxLength = 4096, AddSpecialTokens = true };

    /// <summary>Maximum number of tokens per cell, special tokens included.</summary>
    public int MaxLength { get; init; } = 4096;

    /// <summary>Counts are scaled so each cell sums to this value before median scaling (Geneformer uses 10,000).</summary>
    public double TargetSum { get; init; } = 10_000;

    /// <summary>Frame the ranked genes with <c>&lt;cls&gt;</c> and <c>&lt;eos&gt;</c> (Geneformer V2).</summary>
    public bool AddSpecialTokens { get; init; } = true;
}

/// <summary>The rank-value encoding of one cell, with the intermediate values for inspection.</summary>
/// <param name="TokenIds">Model input ids, special tokens included.</param>
/// <param name="GeneIds">The genes behind the non-special tokens, highest rank first.</param>
/// <param name="NormalizedValues">Each ranked gene's median-scaled expression, in the same order as <paramref name="GeneIds"/>.</param>
public sealed record GeneformerEncoding(int[] TokenIds, IReadOnlyList<string> GeneIds, IReadOnlyList<double> NormalizedValues)
{
    /// <summary>Number of model tokens.</summary>
    public int Length => TokenIds.Length;
}

/// <summary>
/// Geneformer's rank-value encoding (Theodoris et al., <i>Nature</i> 2023): a cell's transcriptome becomes the
/// list of its expressed genes, ordered by expression relative to each gene's median across the pretraining corpus.
/// </summary>
/// <remarks>
/// Per cell: counts are divided by the cell's total and multiplied by <see cref="GeneformerEncoderOptions.TargetSum"/>,
/// then divided by the gene's non-zero median; genes with zero expression or missing from either dictionary are
/// dropped; the rest are sorted by the scaled value (descending, ties in input order), truncated and mapped to token
/// ids. Housekeeping genes, high everywhere, sink; genes that distinguish the cell state rise.
/// </remarks>
public sealed class GeneformerEncoder
{
    private readonly IReadOnlyDictionary<string, int> _tokens;
    private readonly IReadOnlyDictionary<string, double> _medians;

    /// <summary>Creates an encoder.</summary>
    /// <param name="tokenDictionary">Gene id (Ensembl) to token id, including <c>&lt;pad&gt;</c> and <c>&lt;mask&gt;</c> (and <c>&lt;cls&gt;</c>/<c>&lt;eos&gt;</c> for V2).</param>
    /// <param name="geneMedians">Gene id to its non-zero median normalised expression in the pretraining corpus.</param>
    /// <param name="options">Encoding settings; <see cref="GeneformerEncoderOptions.V2"/> by default.</param>
    public GeneformerEncoder(
        IReadOnlyDictionary<string, int> tokenDictionary,
        IReadOnlyDictionary<string, double> geneMedians,
        GeneformerEncoderOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(tokenDictionary);
        ArgumentNullException.ThrowIfNull(geneMedians);
        _tokens = tokenDictionary;
        _medians = geneMedians;
        Options = options ?? GeneformerEncoderOptions.V2;
        if (Options.MaxLength < (Options.AddSpecialTokens ? 3 : 1))
            throw new ArgumentOutOfRangeException(nameof(options), "MaxLength is too small to hold any gene.");
        if (Options.TargetSum <= 0)
            throw new ArgumentOutOfRangeException(nameof(options), "TargetSum must be positive.");
        if (Options.AddSpecialTokens && (ClsId is null || EosId is null))
            throw new ArgumentException("AddSpecialTokens needs <cls> and <eos> in the token dictionary (Geneformer V2).", nameof(tokenDictionary));
    }

    /// <summary>The settings in use.</summary>
    public GeneformerEncoderOptions Options { get; }

    /// <summary>Id of <c>&lt;pad&gt;</c> (0 in published dictionaries).</summary>
    public int PadId => _tokens.TryGetValue("<pad>", out int id) ? id : 0;

    /// <summary>Id of <c>&lt;mask&gt;</c>, when present.</summary>
    public int? MaskId => _tokens.TryGetValue("<mask>", out int id) ? id : null;

    /// <summary>Id of <c>&lt;cls&gt;</c>, when present.</summary>
    public int? ClsId => _tokens.TryGetValue("<cls>", out int id) ? id : null;

    /// <summary>Id of <c>&lt;eos&gt;</c>, when present.</summary>
    public int? EosId => _tokens.TryGetValue("<eos>", out int id) ? id : null;

    /// <summary>Number of entries in the token dictionary.</summary>
    public int VocabularySize => _tokens.Count;

    /// <summary>Encodes one cell to model token ids.</summary>
    /// <param name="geneIds">Gene ids, one per column of the expression vector.</param>
    /// <param name="counts">Raw counts for the cell, aligned with <paramref name="geneIds"/>.</param>
    /// <param name="totalCounts">The cell's total count; when null the sum of <paramref name="counts"/> is used.</param>
    public int[] Encode(IReadOnlyList<string> geneIds, ReadOnlySpan<double> counts, double? totalCounts = null)
        => EncodeDetailed(geneIds, counts, totalCounts).TokenIds;

    /// <summary>Encodes one cell given as a gene-to-count map.</summary>
    public int[] Encode(IReadOnlyDictionary<string, double> expression, double? totalCounts = null)
    {
        ArgumentNullException.ThrowIfNull(expression);
        var genes = new string[expression.Count];
        var counts = new double[expression.Count];
        int i = 0;
        foreach (var (gene, count) in expression)
        {
            genes[i] = gene;
            counts[i++] = count;
        }

        return EncodeDetailed(genes, counts, totalCounts).TokenIds;
    }

    /// <summary>Encodes one cell and returns the ranked genes and their scaled values alongside the token ids.</summary>
    /// <inheritdoc cref="Encode(IReadOnlyList{string}, ReadOnlySpan{double}, double?)"/>
    public GeneformerEncoding EncodeDetailed(IReadOnlyList<string> geneIds, ReadOnlySpan<double> counts, double? totalCounts = null)
    {
        ArgumentNullException.ThrowIfNull(geneIds);
        if (geneIds.Count != counts.Length)
            throw new ArgumentException($"{geneIds.Count} gene ids but {counts.Length} counts.", nameof(counts));

        double total = 0;
        if (totalCounts is double t)
        {
            total = t;
        }
        else
        {
            foreach (double c in counts)
            {
                if (c > 0 && double.IsFinite(c)) total += c;
            }
        }

        var ranked = new List<(int Order, string Gene, double Value, int Token)>();
        if (total > 0)
        {
            for (int i = 0; i < counts.Length; i++)
            {
                double c = counts[i];
                if (!(c > 0) || !double.IsFinite(c)) continue;
                string gene = geneIds[i];
                if (!_medians.TryGetValue(gene, out double median) || !(median > 0)) continue;
                if (!_tokens.TryGetValue(gene, out int token)) continue;
                ranked.Add((i, gene, c / total * Options.TargetSum / median, token));
            }
        }

        ranked.Sort(static (a, b) =>
        {
            int cmp = b.Value.CompareTo(a.Value);
            return cmp != 0 ? cmp : a.Order.CompareTo(b.Order);
        });

        int capacity = Options.AddSpecialTokens ? Options.MaxLength - 2 : Options.MaxLength;
        int n = Math.Min(capacity, ranked.Count);
        var ids = new List<int>(n + 2);
        var genes = new string[n];
        var values = new double[n];
        if (Options.AddSpecialTokens) ids.Add(ClsId!.Value);
        for (int i = 0; i < n; i++)
        {
            ids.Add(ranked[i].Token);
            genes[i] = ranked[i].Gene;
            values[i] = ranked[i].Value;
        }

        if (Options.AddSpecialTokens) ids.Add(EosId!.Value);
        return new GeneformerEncoding(ids.ToArray(), genes, values);
    }

    /// <summary>
    /// Encodes every row of a cell-by-gene count matrix (row-major, <paramref name="geneIds"/>.Count columns).
    /// </summary>
    public IReadOnlyList<int[]> EncodeMatrix(IReadOnlyList<string> geneIds, ReadOnlySpan<double> matrix)
    {
        ArgumentNullException.ThrowIfNull(geneIds);
        int cols = geneIds.Count;
        if (cols == 0 || matrix.Length % cols != 0)
            throw new ArgumentException("The matrix length must be a multiple of the number of genes.", nameof(matrix));
        int rows = matrix.Length / cols;
        var result = new int[rows][];
        for (int r = 0; r < rows; r++) result[r] = Encode(geneIds, matrix.Slice(r * cols, cols));
        return result;
    }
}
