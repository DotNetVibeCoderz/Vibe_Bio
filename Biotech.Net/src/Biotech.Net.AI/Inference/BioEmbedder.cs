using Biotech.Net.AI.Models;
using Biotech.Net.AI.Preprocessing;

namespace Biotech.Net.AI.Inference;

/// <summary>How per-token hidden states are reduced to one vector.</summary>
public enum PoolingStrategy
{
    /// <summary>Mean over the non-special tokens (residues, k-mers or word pieces). The usual choice for pretrained models.</summary>
    Mean,

    /// <summary>The first token's (<c>[CLS]</c> / <c>&lt;cls&gt;</c>) hidden state. Meaningful mainly for fine-tuned models.</summary>
    Cls,

    /// <summary>Mean over every token, special tokens included.</summary>
    MeanAllTokens,
}

/// <summary>
/// Embeds biological sequences and biomedical text with a pretrained model: preprocessing for the model's input
/// convention, windowing of sequences longer than its context, pooling and batching.
/// </summary>
/// <example>
/// <code>
/// using var embedder = BioEmbedder.Load("Esm2T6");
/// float[] v = embedder.Embed("MKTAYIAKQRQISFVKSHFSRQ");
/// </code>
/// </example>
public sealed class BioEmbedder : IDisposable
{
    /// <summary>Creates an embedder over an encoder.</summary>
    /// <param name="encoder">The encoder.</param>
    /// <param name="descriptor">How inputs are preprocessed; null passes text through unchanged.</param>
    public BioEmbedder(IBioEncoder encoder, BioModelDescriptor? descriptor = null)
    {
        Encoder = encoder ?? throw new ArgumentNullException(nameof(encoder));
        Descriptor = descriptor;
    }

    /// <summary>
    /// Loads a model by catalogue key (<c>Esm2T6</c>, <c>ProtBert</c>, <c>DnaBert6</c>, <c>PubMedBert</c>, …), Hub repository id
    /// or local directory, downloading it on first use.
    /// </summary>
    /// <exception cref="Hub.BioModelUnavailableException">Not cached and not downloadable (offline, network, gated).</exception>
    /// <exception cref="NotSupportedException">The model's architecture has no runtime here, or it is a Geneformer model (use <see cref="GeneformerModel"/>).</exception>
    public static BioEmbedder Load(string keyOrRepoId)
    {
        ArgumentException.ThrowIfNullOrEmpty(keyOrRepoId);
        BioModelCatalog.TryGet(keyOrRepoId, out var descriptor);
        string repo = descriptor?.RepoId ?? keyOrRepoId;
        var backend = descriptor?.Backend ?? InferBackend(repo);
        IBioEncoder encoder = backend switch
        {
            BioModelBackend.Esm => EsmEncoder.Load(repo),
            BioModelBackend.HfTransformer => HfTransformerEncoder.Load(repo, descriptor?.MaxTokens ?? 512),
            BioModelBackend.Geneformer => throw new NotSupportedException(
                "Geneformer consumes gene-token ids, not sequences: use GeneformerModel.Load and GeneformerEncoder."),
            _ => throw new NotSupportedException($"{descriptor?.DisplayName ?? repo} cannot run in Biotech.Net: {descriptor?.Description}"),
        };
        return new BioEmbedder(encoder, descriptor);
    }

    /// <summary>The encoder.</summary>
    public IBioEncoder Encoder { get; }

    /// <summary>The model description, when known.</summary>
    public BioModelDescriptor? Descriptor { get; }

    /// <summary>Embedding width.</summary>
    public int Dimension => Encoder.HiddenSize;

    /// <summary>How token states are pooled. Defaults to <see cref="PoolingStrategy.Mean"/>.</summary>
    public PoolingStrategy Pooling { get; init; } = PoolingStrategy.Mean;

    /// <summary>Residues shared by consecutive windows when a sequence is longer than the model's context. Defaults to 64.</summary>
    public int WindowOverlap { get; init; } = 64;

    /// <summary>Maximum parallel inputs in <see cref="EmbedBatch"/>. Defaults to 2 (each forward pass is itself parallel).</summary>
    public int MaxDegreeOfParallelism { get; init; } = 2;

    /// <summary>
    /// Embeds one sequence or text. Sequences longer than the context are split into overlapping windows whose
    /// pooled vectors are averaged, weighted by window length.
    /// </summary>
    public float[] Embed(string input)
    {
        ArgumentNullException.ThrowIfNull(input);
        var windows = Windows(input);
        if (windows.Count == 1) return Pool(Encoder.Encode(windows[0].ModelInput));
        var sum = new double[Dimension];
        double weight = 0;
        foreach (var (modelInput, length) in windows)
        {
            var v = Pool(Encoder.Encode(modelInput));
            for (int i = 0; i < v.Length; i++) sum[i] += v[i] * length;
            weight += length;
        }

        var result = new float[Dimension];
        for (int i = 0; i < result.Length; i++) result[i] = (float)(sum[i] / weight);
        return result;
    }

    /// <summary>Embeds many inputs, one vector per input, in input order.</summary>
    public float[][] EmbedBatch(IReadOnlyList<string> inputs, IProgress<int>? progress = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        var result = new float[inputs.Count][];
        int done = 0;
        Parallel.For(0, inputs.Count,
            new ParallelOptions { MaxDegreeOfParallelism = Math.Max(1, MaxDegreeOfParallelism), CancellationToken = cancellationToken },
            i =>
            {
                result[i] = Embed(inputs[i]);
                progress?.Report(Interlocked.Increment(ref done));
            });
        return result;
    }

    /// <summary>
    /// Per-residue (per non-special token) embeddings of one input that fits the context, row-major
    /// <c>[tokens, Dimension]</c>. For ESM this is one row per amino acid.
    /// </summary>
    public float[][] TokenEmbeddings(string input)
    {
        var states = Encoder.Encode(Preprocess(input));
        var rows = new List<float[]>(states.Tokens);
        for (int t = 0; t < states.Tokens; t++)
        {
            if (!states.IsSpecial[t]) rows.Add(states.Row(t).ToArray());
        }

        return [.. rows];
    }

    /// <summary>Cosine similarity of two vectors, in [-1, 1].</summary>
    public static double CosineSimilarity(ReadOnlySpan<float> a, ReadOnlySpan<float> b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Vectors differ in length.");
        double dot = 0, na = 0, nb = 0;
        for (int i = 0; i < a.Length; i++)
        {
            dot += a[i] * b[i];
            na += a[i] * a[i];
            nb += b[i] * b[i];
        }

        return na == 0 || nb == 0 ? 0 : dot / Math.Sqrt(na * nb);
    }

    /// <summary>Pools token states according to <see cref="Pooling"/>.</summary>
    public float[] Pool(TokenStates states)
    {
        ArgumentNullException.ThrowIfNull(states);
        var result = new float[states.HiddenSize];
        if (Pooling == PoolingStrategy.Cls)
        {
            states.Row(0).CopyTo(result);
            return result;
        }

        var sum = new double[states.HiddenSize];
        int count = 0;
        for (int t = 0; t < states.Tokens; t++)
        {
            if (Pooling == PoolingStrategy.Mean && states.IsSpecial[t]) continue;
            var row = states.Row(t);
            for (int i = 0; i < row.Length; i++) sum[i] += row[i];
            count++;
        }

        if (count == 0) return Pool(states with { IsSpecial = new bool[states.Tokens] });
        for (int i = 0; i < result.Length; i++) result[i] = (float)(sum[i] / count);
        return result;
    }

    /// <inheritdoc />
    public void Dispose() => Encoder.Dispose();

    private string Preprocess(string raw) => Descriptor is null ? raw : Descriptor.Preprocess(raw);

    private List<(string ModelInput, int Length)> Windows(string raw)
    {
        if (Descriptor is null || Descriptor.InputKind == BioSequenceKind.BiomedicalText)
            return [(Preprocess(raw), raw.Length)];

        string clean = Descriptor.InputKind == BioSequenceKind.Dna
            ? SequencePreprocessor.CleanDna(raw)
            : SequencePreprocessor.CleanProtein(raw);
        int budget = Math.Min(Descriptor.MaxTokens, Encoder.MaxTokens) - 2;
        int maxResidues = Descriptor.KmerSize is int k
            ? budget + k - 1
            : Descriptor.Backend == BioModelBackend.Esm && Descriptor.InputKind == BioSequenceKind.Dna ? budget * 6 : budget;
        int overlap = Math.Min(WindowOverlap, maxResidues / 2);
        return [.. SequencePreprocessor.Chunk(clean, maxResidues, overlap).Select(c => (Descriptor.Preprocess(c.Text), c.Length))];
    }

    private static BioModelBackend InferBackend(string repoOrDir)
    {
        string dir = Hub.BioModelHub.EnsureModel(repoOrDir);
        string config = File.ReadAllText(Path.Combine(dir, "config.json"));
        return config.Contains("\"esm\"", StringComparison.OrdinalIgnoreCase) ? BioModelBackend.Esm : BioModelBackend.HfTransformer;
    }
}
