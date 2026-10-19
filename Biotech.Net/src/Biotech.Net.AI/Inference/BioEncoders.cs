using Biotech.Net.AI.Hub;
using Gravicode.HFNet.GraviTokenizers;
using Gravicode.HFNet.GraviTransformers;

namespace Biotech.Net.AI.Inference;

/// <summary>Per-token hidden states of one input.</summary>
/// <param name="Hidden">Row-major <c>[Tokens, HiddenSize]</c>.</param>
/// <param name="Tokens">Number of tokens.</param>
/// <param name="HiddenSize">Width of each row.</param>
/// <param name="IsSpecial">Whether each token is a special token (CLS, SEP, EOS, padding).</param>
public sealed record TokenStates(float[] Hidden, int Tokens, int HiddenSize, bool[] IsSpecial)
{
    /// <summary>The hidden state of token <paramref name="index"/>.</summary>
    public ReadOnlySpan<float> Row(int index) => Hidden.AsSpan(index * HiddenSize, HiddenSize);
}

/// <summary>An encoder that turns model-ready text into per-token hidden states.</summary>
public interface IBioEncoder : IDisposable
{
    /// <summary>The checkpoint id.</summary>
    string ModelId { get; }

    /// <summary>Width of the hidden states.</summary>
    int HiddenSize { get; }

    /// <summary>Maximum tokens per input, special tokens included.</summary>
    int MaxTokens { get; }

    /// <summary>Encodes one model-ready input (already preprocessed: k-mers, spaced residues, raw protein, or text).</summary>
    TokenStates Encode(string modelInput);
}

/// <summary>A BERT-family encoder (BioBERT, PubMedBERT, DNABERT, ProtBERT, Geneformer, …) run by HF.Net's <see cref="TransformerModel"/>.</summary>
public sealed class HfTransformerEncoder : IBioEncoder
{
    private readonly bool _ownsModel;

    /// <summary>Wraps a loaded model.</summary>
    public HfTransformerEncoder(TransformerModel model, int maxTokens = 512, bool ownsModel = true)
    {
        Model = model ?? throw new ArgumentNullException(nameof(model));
        _ownsModel = ownsModel;
        MaxTokens = Math.Min(maxTokens, model.Config.MaxPositions > 0 ? model.Config.MaxPositions : maxTokens);
    }

    /// <summary>Downloads (or reads from cache) and loads a BERT-family checkpoint.</summary>
    /// <exception cref="BioModelUnavailableException">The checkpoint is not cached and cannot be downloaded.</exception>
    /// <exception cref="NotSupportedException">HF.Net cannot run the architecture.</exception>
    public static HfTransformerEncoder Load(string repoIdOrDirectory, int maxTokens = 512, string revision = "main")
        => new(TransformerModel.Open(BioModelHub.EnsureModel(repoIdOrDirectory, revision)), maxTokens);

    /// <summary>The underlying HF.Net model.</summary>
    public TransformerModel Model { get; }

    /// <inheritdoc />
    public string ModelId => Model.RepoId;

    /// <inheritdoc />
    public int HiddenSize => Model.Config.HiddenSize;

    /// <inheritdoc />
    public int MaxTokens { get; }

    /// <inheritdoc />
    public TokenStates Encode(string modelInput)
    {
        ArgumentNullException.ThrowIfNull(modelInput);
        var hidden = Model.Hidden(modelInput, MaxTokens);
        var values = hidden.ToArray();
        int h = HiddenSize;
        int n = values.Length / h;
        var special = new bool[n];
        if (Model.Tokenizer is { } tok)
        {
            var enc = tok.Encode(modelInput);
            for (int i = 0; i < n && i < enc.SpecialTokensMask.Count; i++) special[i] = enc.SpecialTokensMask[i] != 0;
            // Truncation keeps the final [SEP]; mark the last row special when the encoding was cut.
            if (enc.Length > n && n > 0) special[n - 1] = true;
        }

        var f = new float[values.Length];
        for (int i = 0; i < f.Length; i++) f[i] = (float)values[i];
        return new TokenStates(f, n, h, special);
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (_ownsModel) Model.Dispose();
    }
}

/// <summary>
/// An ESM-family encoder (ESM-1b, ESM-1v, ESM-2, Nucleotide Transformer v1) run by HF.Net's <see cref="TransformerModel"/>.
/// </summary>
public sealed class EsmEncoder : IBioEncoder
{
    private readonly bool _ownsModel;

    /// <summary>Wraps a loaded ESM model.</summary>
    /// <exception cref="ArgumentException">The model is not an ESM model or has no tokenizer.</exception>
    public EsmEncoder(TransformerModel model, bool ownsModel = true)
    {
        Model = model ?? throw new ArgumentNullException(nameof(model));
        Esm.Require(model, nameof(model));
        _ownsModel = ownsModel;
    }

    /// <summary>Downloads (or reads from cache) and loads an ESM checkpoint, for example <c>facebook/esm2_t6_8M_UR50D</c>.</summary>
    /// <exception cref="BioModelUnavailableException">The checkpoint is not cached and cannot be downloaded.</exception>
    /// <exception cref="NotSupportedException">The checkpoint is not a supported ESM language model.</exception>
    public static EsmEncoder Load(string repoIdOrDirectory, string revision = "main") => new(Esm.Load(repoIdOrDirectory, revision));

    /// <summary>The underlying HF.Net model.</summary>
    public TransformerModel Model { get; }

    /// <inheritdoc />
    public string ModelId => Model.RepoId;

    /// <inheritdoc />
    public int HiddenSize => Model.Config.HiddenSize;

    /// <inheritdoc />
    public int MaxTokens => Esm.MaxTokens(Model);

    /// <inheritdoc />
    public TokenStates Encode(string modelInput)
    {
        ArgumentNullException.ThrowIfNull(modelInput);
        var ids = Esm.Encode(Model, modelInput, MaxTokens);
        var hidden = Model.Hidden(ids).ToArray();
        var f = new float[hidden.Length];
        for (int i = 0; i < f.Length; i++) f[i] = (float)hidden[i];
        var special = new bool[ids.Length];
        for (int i = 0; i < ids.Length; i++) special[i] = Esm.IsSpecial(Model.Tokenizer, ids[i]);
        return new TokenStates(f, ids.Length, HiddenSize, special);
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (_ownsModel) Model.Dispose();
    }
}

/// <summary>ESM conventions shared by the encoder, fill-mask and the variant scorer.</summary>
internal static class Esm
{
    /// <summary>Loads an ESM checkpoint through the Biotech.Net hub (cache, offline mode, clear errors).</summary>
    public static TransformerModel Load(string repoIdOrDirectory, string revision = "main")
    {
        var model = TransformerModel.Open(BioModelHub.EnsureModel(repoIdOrDirectory, revision));
        try
        {
            Require(model, nameof(repoIdOrDirectory));
            return model;
        }
        catch
        {
            model.Dispose();
            throw;
        }
    }

    public static void Require(TransformerModel model, string paramName)
    {
        if (model.Config.Esm is null)
            throw new ArgumentException($"'{model.RepoId}' is a '{model.Config.ModelType}' model, not an ESM model.", paramName);
        if (!model.HasTokenizer)
            throw new ArgumentException($"'{model.RepoId}' has no ESM vocabulary (vocab.txt).", paramName);
    }

    /// <summary>Tokens per input, special tokens included.</summary>
    public static int MaxTokens(TransformerModel model)
        => Math.Min(model.Config.MaxPositions > 0 ? model.Config.MaxPositions : 1024, 1024);

    /// <summary>
    /// Token ids with <c>&lt;cls&gt;</c> (and <c>&lt;eos&gt;</c> when the model uses one); longer inputs keep
    /// <c>&lt;cls&gt;</c>, are cut, and still end with <c>&lt;eos&gt;</c>.
    /// </summary>
    public static int[] Encode(TransformerModel model, string text, int maxTokens)
    {
        var ids = model.Tokenizer.Encode(text).ToIdArray();
        if (ids.Length <= maxTokens) return ids;
        int last = ids[^1];
        var cut = ids[..maxTokens];
        if (model.Tokenizer.IdToToken(last) == "<eos>") cut[^1] = last;
        return cut;
    }

    /// <summary>Whether a token is one of ESM's angle-bracket specials (cls, pad, eos, unk, mask, null).</summary>
    public static bool IsSpecial(HfTokenizer tokenizer, int id)
    {
        string t = tokenizer.IdToToken(id);
        return t.Length > 2 && t[0] == '<' && t[^1] == '>';
    }

    /// <summary>Id of a token that is in the vocabulary, or -1 (HF.Net maps unknown text to <c>&lt;unk&gt;</c>).</summary>
    public static int IdOf(HfTokenizer tokenizer, string token)
    {
        int id = tokenizer.TokenToId(token);
        return id >= 0 && tokenizer.IdToToken(id) == token ? id : -1;
    }

    /// <summary>Log-softmax of <paramref name="logits"/>.</summary>
    public static double[] LogSoftmax(ReadOnlySpan<double> logits)
    {
        double max = double.NegativeInfinity;
        foreach (double v in logits) max = Math.Max(max, v);
        double sum = 0;
        foreach (double v in logits) sum += Math.Exp(v - max);
        double lse = max + Math.Log(sum);
        var result = new double[logits.Length];
        for (int i = 0; i < logits.Length; i++) result[i] = logits[i] - lse;
        return result;
    }
}
