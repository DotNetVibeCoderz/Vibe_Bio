using Biotech.Net.AI.Hub;
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

/// <summary>An ESM-family encoder run by the managed <see cref="EsmModel"/>.</summary>
public sealed class EsmEncoder : IBioEncoder
{
    /// <summary>Wraps a loaded model.</summary>
    public EsmEncoder(EsmModel model) => Model = model ?? throw new ArgumentNullException(nameof(model));

    /// <summary>Downloads (or reads from cache) and loads an ESM checkpoint.</summary>
    public static EsmEncoder Load(string repoIdOrDirectory, string revision = "main") => new(EsmModel.Load(repoIdOrDirectory, revision));

    /// <summary>The underlying model.</summary>
    public EsmModel Model { get; }

    /// <inheritdoc />
    public string ModelId => Model.Id;

    /// <inheritdoc />
    public int HiddenSize => Model.HiddenSize;

    /// <inheritdoc />
    public int MaxTokens => Math.Min(Model.Config.MaxPositions, 1024);

    /// <inheritdoc />
    public TokenStates Encode(string modelInput)
    {
        var ids = Model.Tokenizer.Encode(modelInput);
        if (ids.Length > MaxTokens)
        {
            // Keep <cls>, truncate, and end with <eos> when the model uses one.
            var cut = ids[..MaxTokens];
            if (Model.Tokenizer.EosId >= 0) cut[^1] = Model.Tokenizer.EosId;
            ids = cut;
        }

        var hidden = Model.Forward(ids);
        var special = new bool[ids.Length];
        for (int i = 0; i < ids.Length; i++) special[i] = Model.Tokenizer.IsSpecial(ids[i]);
        return new TokenStates(hidden, ids.Length, HiddenSize, special);
    }

    /// <inheritdoc />
    public void Dispose() => Model.Dispose();
}
