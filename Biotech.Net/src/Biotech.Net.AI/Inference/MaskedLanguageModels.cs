using System.Globalization;
using System.Text.RegularExpressions;
using Biotech.Net.AI.Models;
using Biotech.Net.AI.Preprocessing;
using Gravicode.HFNet.GraviTransformers;

namespace Biotech.Net.AI.Inference;

/// <summary>One candidate for a masked position.</summary>
/// <param name="Token">The predicted token (residue, k-mer or word piece).</param>
/// <param name="Score">Its probability.</param>
/// <param name="Sequence">The input with the mask replaced by the token.</param>
public sealed record BioMaskPrediction(string Token, double Score, string Sequence);

/// <summary>
/// Fill-mask over a biological model: ESM models through <see cref="EsmModel"/>, BERT-family models
/// (ProtBERT, DNABERT, BioBERT, PubMedBERT) through HF.Net's <see cref="TransformerModel.FillMask"/>.
/// </summary>
public sealed class BioMaskedLanguageModel : IDisposable
{
    private readonly EsmModel? _esm;
    private readonly TransformerModel? _bert;

    /// <summary>Wraps an ESM model.</summary>
    public BioMaskedLanguageModel(EsmModel model) => _esm = model ?? throw new ArgumentNullException(nameof(model));

    /// <summary>Wraps an HF.Net model that has a masked-LM head.</summary>
    public BioMaskedLanguageModel(TransformerModel model) => _bert = model ?? throw new ArgumentNullException(nameof(model));

    /// <summary>Loads a model by catalogue key, repository id or directory.</summary>
    /// <exception cref="Hub.BioModelUnavailableException">Not cached and not downloadable.</exception>
    public static BioMaskedLanguageModel Load(string keyOrRepoId)
    {
        BioModelCatalog.TryGet(keyOrRepoId, out var d);
        string repo = d?.RepoId ?? keyOrRepoId;
        if (d is { Backend: BioModelBackend.NotSupported })
            throw new NotSupportedException($"{d.DisplayName} cannot run in Biotech.Net: {d.Description}");
        string dir = Hub.BioModelHub.EnsureModel(repo);
        bool esm = d?.Backend == BioModelBackend.Esm
                   || File.ReadAllText(Path.Combine(dir, "config.json")).Contains("\"esm\"", StringComparison.OrdinalIgnoreCase);
        return esm ? new BioMaskedLanguageModel(EsmModel.Open(dir, repo)) : new BioMaskedLanguageModel(TransformerModel.Open(dir));
    }

    /// <summary>The token that marks the position to predict (<c>&lt;mask&gt;</c> for ESM, usually <c>[MASK]</c> for BERT).</summary>
    public string MaskToken => _esm is not null ? _esm.Tokenizer.MaskToken
        : _bert!.Tokenizer.IdToToken(_bert.Tokenizer.MaskId);

    /// <summary>The ESM model, when this wraps one.</summary>
    public EsmModel? Esm => _esm;

    /// <summary>Predicts the token at the (first) mask in model-ready <paramref name="input"/>.</summary>
    public IReadOnlyList<BioMaskPrediction> FillMask(string input, int topK = 5)
    {
        ArgumentNullException.ThrowIfNull(input);
        ArgumentOutOfRangeException.ThrowIfLessThan(topK, 1);
        if (_bert is not null)
            return [.. _bert.FillMask(input, topK).Select(f => new BioMaskPrediction(f.Token, f.Score, f.Sequence))];

        var tok = _esm!.Tokenizer;
        var ids = tok.Encode(input);
        int pos = Array.IndexOf(ids, tok.MaskId);
        if (pos < 0) throw new ArgumentException($"The input has no {tok.MaskToken} token.", nameof(input));
        var logp = EsmModel.LogSoftmax(_esm.LogitsAt(_esm.Forward(ids), pos));
        int maskAt = input.IndexOf(tok.MaskToken, StringComparison.Ordinal);
        return [.. Enumerable.Range(0, logp.Length)
            .Where(i => !tok.IsSpecial(i))
            .OrderByDescending(i => logp[i])
            .Take(topK)
            .Select(i => new BioMaskPrediction(tok.Token(i), Math.Exp(logp[i]),
                string.Concat(input.AsSpan(0, maskAt), tok.Token(i), input.AsSpan(maskAt + tok.MaskToken.Length))))];
    }

    /// <summary>
    /// ESM only: masks residue <paramref name="position"/> (1-based) of a protein and predicts it.
    /// </summary>
    public IReadOnlyList<BioMaskPrediction> PredictResidue(string protein, int position, int topK = 5)
    {
        if (_esm is null) throw new NotSupportedException("PredictResidue needs an ESM model; use FillMask with spaced residues for ProtBERT.");
        string seq = SequencePreprocessor.CleanProtein(protein);
        if (position < 1 || position > seq.Length) throw new ArgumentOutOfRangeException(nameof(position));
        return FillMask(string.Concat(seq.AsSpan(0, position - 1), _esm.Tokenizer.MaskToken, seq.AsSpan(position)), topK);
    }

    /// <inheritdoc />
    public void Dispose()
    {
        _esm?.Dispose();
        _bert?.Dispose();
    }
}

/// <summary>A single amino-acid substitution such as <c>A24G</c> (1-based position).</summary>
/// <param name="WildType">Residue in the reference sequence.</param>
/// <param name="Position">1-based position.</param>
/// <param name="Mutant">Substituted residue.</param>
public sealed partial record ProteinMutation(char WildType, int Position, char Mutant)
{
    /// <summary>Parses <c>A24G</c> (also accepts the <c>p.</c> prefix: <c>p.A24G</c>).</summary>
    /// <exception cref="FormatException">The text is not a substitution.</exception>
    public static ProteinMutation Parse(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var m = Pattern().Match(text.Trim());
        if (!m.Success) throw new FormatException($"'{text}' is not a substitution like A24G.");
        return new ProteinMutation(char.ToUpperInvariant(m.Groups[1].Value[0]),
            int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture), char.ToUpperInvariant(m.Groups[3].Value[0]));
    }

    /// <summary>Parses a multi-mutant like <c>A24G:K30R</c> (separators <c>:</c>, <c>,</c>, <c>;</c>, <c>+</c> or space).</summary>
    public static IReadOnlyList<ProteinMutation> ParseMany(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        return [.. text.Split([':', ',', ';', '+', ' '], StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries).Select(Parse)];
    }

    /// <inheritdoc />
    public override string ToString() => $"{WildType}{Position.ToString(CultureInfo.InvariantCulture)}{Mutant}";

    [GeneratedRegex(@"^(?:p\.)?([A-Za-z])(\d+)([A-Za-z])$")]
    private static partial Regex Pattern();
}

/// <summary>A variant's predicted effect.</summary>
/// <param name="Mutations">The substitutions.</param>
/// <param name="Score">Sum of per-site log-likelihood ratios; negative means the model finds the variant less likely than wild type (likely deleterious).</param>
/// <param name="PerSite">Log-likelihood ratio of each substitution.</param>
public sealed record VariantEffectScore(IReadOnlyList<ProteinMutation> Mutations, double Score, IReadOnlyList<double> PerSite)
{
    /// <inheritdoc />
    public override string ToString() => $"{string.Join(':', Mutations)} {Score:F3}";
}

/// <summary>A deep mutational scan: the score of every amino acid at every scanned position.</summary>
/// <param name="Positions">1-based positions, in order.</param>
/// <param name="AminoAcids">Column order of <paramref name="Scores"/> (the 20 standard amino acids).</param>
/// <param name="Scores">[position, amino acid] log-likelihood ratios against the wild-type residue.</param>
public sealed record MutationalScan(IReadOnlyList<int> Positions, string AminoAcids, double[,] Scores);

/// <summary>
/// Zero-shot variant-effect prediction with an ESM model by the <b>masked-marginal</b> score of Meier et al. (2021):
/// for each mutated position i, mask it in the wild-type sequence and take
/// log P(x_i = mutant | x_{-i}) − log P(x_i = wild type | x_{-i}); a multi-mutant scores the sum over its sites.
/// </summary>
public sealed class EsmVariantScorer
{
    /// <summary>The 20 standard amino acids, in alphabetical one-letter order.</summary>
    public const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    private readonly Dictionary<(string Sequence, int Position), double[]> _cache = [];
    private readonly Lock _gate = new();

    /// <summary>Creates a scorer.</summary>
    /// <exception cref="ArgumentException">The model has no language-model head.</exception>
    public EsmVariantScorer(EsmModel model)
    {
        Model = model ?? throw new ArgumentNullException(nameof(model));
        if (!model.HasLanguageModelHead) throw new ArgumentException("Variant scoring needs a checkpoint with a masked-LM head.", nameof(model));
    }

    /// <summary>The model.</summary>
    public EsmModel Model { get; }

    /// <summary>Log-probabilities over the vocabulary at <paramref name="position"/> (1-based) with that residue masked.</summary>
    public double[] SiteLogProbabilities(string wildType, int position)
    {
        string seq = SequencePreprocessor.CleanProtein(wildType);
        if (position < 1 || position > seq.Length) throw new ArgumentOutOfRangeException(nameof(position), $"Position {position} is outside 1..{seq.Length}.");
        lock (_gate)
        {
            if (_cache.TryGetValue((seq, position), out var hit)) return hit;
        }

        // Sequences longer than the context are scored in a window centred on the site.
        int maxResidues = Math.Min(Model.Config.MaxPositions, 1024) - 2;
        int start = 0;
        string window = seq;
        if (seq.Length > maxResidues)
        {
            start = Math.Clamp(position - 1 - maxResidues / 2, 0, seq.Length - maxResidues);
            window = seq.Substring(start, maxResidues);
        }

        var tok = Model.Tokenizer;
        var ids = tok.Encode(window);
        int tokenIndex = position - start; // <cls> occupies index 0
        ids[tokenIndex] = tok.MaskId;
        var logp = EsmModel.LogSoftmax(Model.LogitsAt(Model.Forward(ids), tokenIndex));
        lock (_gate) _cache[(seq, position)] = logp;
        return logp;
    }

    /// <summary>Scores a single substitution.</summary>
    /// <exception cref="ArgumentException">The wild-type residue does not match the sequence.</exception>
    public double Score(string wildType, ProteinMutation mutation)
    {
        ArgumentNullException.ThrowIfNull(mutation);
        string seq = SequencePreprocessor.CleanProtein(wildType);
        if (mutation.Position < 1 || mutation.Position > seq.Length)
            throw new ArgumentException($"{mutation} lies outside the {seq.Length}-residue sequence.", nameof(mutation));
        if (seq[mutation.Position - 1] != mutation.WildType)
            throw new ArgumentException($"{mutation}: residue {mutation.Position} is {seq[mutation.Position - 1]}, not {mutation.WildType}.", nameof(mutation));
        int wt = Model.Tokenizer.Id(mutation.WildType.ToString());
        int mt = Model.Tokenizer.Id(mutation.Mutant.ToString());
        if (mt < 0 || wt < 0) throw new ArgumentException($"{mutation}: residue not in the model vocabulary.", nameof(mutation));
        var logp = SiteLogProbabilities(seq, mutation.Position);
        return logp[mt] - logp[wt];
    }

    /// <summary>Scores a single or multi-mutant written like <c>A24G</c> or <c>A24G:K30R</c>.</summary>
    public VariantEffectScore Score(string wildType, string mutations)
    {
        var list = ProteinMutation.ParseMany(mutations);
        if (list.Count == 0) throw new ArgumentException("No mutations given.", nameof(mutations));
        var per = list.Select(m => Score(wildType, m)).ToArray();
        return new VariantEffectScore(list, per.Sum(), per);
    }

    /// <summary>Scores many variants of one wild type; shared sites are computed once.</summary>
    public IReadOnlyList<VariantEffectScore> ScoreMany(string wildType, IEnumerable<string> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);
        return [.. variants.Select(v => Score(wildType, v))];
    }

    /// <summary>
    /// Scores all 20 amino acids at each position in [<paramref name="from"/>, <paramref name="to"/>] (1-based, inclusive;
    /// the whole sequence by default) — an in-silico deep mutational scan.
    /// </summary>
    public MutationalScan SaturationScan(string wildType, int? from = null, int? to = null)
    {
        string seq = SequencePreprocessor.CleanProtein(wildType);
        int a = from ?? 1, b = to ?? seq.Length;
        if (a < 1 || b > seq.Length || a > b) throw new ArgumentOutOfRangeException(nameof(from), "Invalid position range.");
        var positions = Enumerable.Range(a, b - a + 1).ToArray();
        var scores = new double[positions.Length, StandardAminoAcids.Length];
        var tok = Model.Tokenizer;
        for (int p = 0; p < positions.Length; p++)
        {
            var logp = SiteLogProbabilities(seq, positions[p]);
            int wt = tok.Id(seq[positions[p] - 1].ToString());
            double wtLog = wt >= 0 ? logp[wt] : 0;
            for (int j = 0; j < StandardAminoAcids.Length; j++)
                scores[p, j] = logp[tok.Id(StandardAminoAcids[j].ToString())] - wtLog;
        }

        return new MutationalScan(positions, StandardAminoAcids, scores);
    }
}

/// <summary>
/// A Geneformer checkpoint run by HF.Net's <see cref="TransformerModel"/> that embeds cells from
/// <see cref="GeneformerEncoder"/> token ids.
/// </summary>
/// <remarks>
/// Geneformer has no text tokenizer, and HF.Net refuses to open a checkpoint without one, so a placeholder
/// <c>vocab.txt</c> (<c>[PAD]</c>, <c>[MASK]</c>, <c>[gene_i]</c>…) sized to the vocabulary is written next to the weights
/// when missing. It is never used to tokenise: inputs are ids from <see cref="GeneformerEncoder"/>.
/// </remarks>
public sealed class GeneformerModel : IDisposable
{
    private GeneformerModel(TransformerModel model, BioModelDescriptor descriptor)
    {
        Model = model;
        Descriptor = descriptor;
    }

    /// <summary>The HF.Net model.</summary>
    public TransformerModel Model { get; }

    /// <summary>The catalogue entry.</summary>
    public BioModelDescriptor Descriptor { get; }

    /// <summary>Embedding width.</summary>
    public int Dimension => Model.Config.HiddenSize;

    /// <summary>
    /// Downloads (or reads from cache) and loads a Geneformer variant; <see cref="BioModelCatalog.GeneformerV1"/> by default
    /// (the 10M model, small enough for a laptop).
    /// </summary>
    public static GeneformerModel Load(BioModelDescriptor? descriptor = null)
    {
        descriptor ??= BioModelCatalog.GeneformerV1;
        if (descriptor.Backend != BioModelBackend.Geneformer) throw new ArgumentException($"{descriptor.Key} is not a Geneformer model.", nameof(descriptor));
        string prefix = descriptor.Subfolder is null ? "" : descriptor.Subfolder + "/";
        string config = Hub.BioModelHub.DownloadFile(descriptor.RepoId, prefix + "config.json");
        Hub.BioModelHub.DownloadFile(descriptor.RepoId, prefix + "model.safetensors");
        return Open(Path.GetDirectoryName(config)!, descriptor);
    }

    /// <summary>Loads a Geneformer checkpoint from a directory with <c>config.json</c> and weights.</summary>
    public static GeneformerModel Open(string directory, BioModelDescriptor? descriptor = null)
    {
        var config = PretrainedConfig.Load(Path.Combine(directory, "config.json"));
        string vocab = Path.Combine(directory, "vocab.txt");
        if (!File.Exists(vocab))
        {
            var lines = new List<string>(config.VocabularySize) { "[PAD]", "[MASK]" };
            for (int i = 2; i < config.VocabularySize; i++) lines.Add($"[gene_{i}]");
            File.WriteAllLines(vocab, lines);
        }

        return new GeneformerModel(TransformerModel.Open(directory), descriptor ?? BioModelCatalog.GeneformerV1);
    }

    /// <summary>Final hidden states for a cell's token ids, row-major <c>[tokens, Dimension]</c>.</summary>
    public float[] Hidden(int[] tokenIds)
    {
        ArgumentNullException.ThrowIfNull(tokenIds);
        if (tokenIds.Length == 0) throw new ArgumentException("The cell has no expressed genes in the dictionary.", nameof(tokenIds));
        int n = Math.Min(tokenIds.Length, Descriptor.MaxTokens);
        var ids = tokenIds[..n];
        var values = Model.Forward(ids, new int[n], Enumerable.Repeat(1, n).ToArray()).ToArray();
        var f = new float[values.Length];
        for (int i = 0; i < f.Length; i++) f[i] = (float)values[i];
        return f;
    }

    /// <summary>
    /// A cell embedding: the mean of the gene-token states (special tokens <paramref name="specialIds"/> excluded),
    /// as Geneformer's <c>EmbExtractor</c> computes in "cell" mode.
    /// </summary>
    public float[] EmbedCell(int[] tokenIds, params int[] specialIds)
    {
        var hidden = Hidden(tokenIds);
        int h = Dimension, n = hidden.Length / h;
        var sum = new double[h];
        int count = 0;
        for (int t = 0; t < n; t++)
        {
            if (specialIds.Contains(tokenIds[t])) continue;
            for (int i = 0; i < h; i++) sum[i] += hidden[t * h + i];
            count++;
        }

        var result = new float[h];
        for (int i = 0; i < h; i++) result[i] = (float)(sum[i] / Math.Max(count, 1));
        return result;
    }

    /// <inheritdoc />
    public void Dispose() => Model.Dispose();
}
