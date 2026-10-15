using System.Numerics;
using System.Runtime.InteropServices;
using System.Text.Json;
using Biotech.Net.AI.Hub;
using Gravicode.HFNet.GraviTransformers;

namespace Biotech.Net.AI.Inference;

/// <summary>The parts of an ESM <c>config.json</c> the forward pass needs.</summary>
public sealed record EsmConfig
{
    /// <summary>Hidden size.</summary>
    public int HiddenSize { get; init; }

    /// <summary>Number of transformer layers.</summary>
    public int Layers { get; init; }

    /// <summary>Number of attention heads.</summary>
    public int Heads { get; init; }

    /// <summary>Feed-forward width.</summary>
    public int IntermediateSize { get; init; }

    /// <summary>Vocabulary size.</summary>
    public int VocabularySize { get; init; }

    /// <summary>Layer-norm epsilon.</summary>
    public double LayerNormEpsilon { get; init; } = 1e-5;

    /// <summary>Id of the mask token.</summary>
    public int MaskTokenId { get; init; } = 32;

    /// <summary>Id of the padding token.</summary>
    public int PadTokenId { get; init; } = 1;

    /// <summary><c>rotary</c> (ESM-2, NT v2) or <c>absolute</c> (ESM-1b, NT v1).</summary>
    public string PositionEmbeddingType { get; init; } = "rotary";

    /// <summary>Whether masked tokens are zeroed and the embeddings rescaled, as ESM was trained.</summary>
    public bool TokenDropout { get; init; }

    /// <summary>Layer norm on the embeddings before the first layer.</summary>
    public bool EmbeddingLayerNormBefore { get; init; }

    /// <summary>Maximum number of positions.</summary>
    public int MaxPositions { get; init; } = 1026;

    /// <summary>Reads <c>config.json</c>.</summary>
    /// <exception cref="NotSupportedException">The checkpoint is not a standard ESM encoder.</exception>
    public static EsmConfig Load(string path)
    {
        using var doc = JsonDocument.Parse(File.ReadAllText(path));
        var r = doc.RootElement;
        string type = Str(r, "model_type") ?? "";
        if (!type.Equals("esm", StringComparison.OrdinalIgnoreCase))
            throw new NotSupportedException($"'{path}' describes a '{type}' model, not an ESM model.");
        if (r.TryGetProperty("add_bias_fnn", out var bias) && bias.ValueKind == JsonValueKind.False)
            throw new NotSupportedException(
                "This checkpoint uses the gated (SwiGLU) feed-forward of Nucleotide Transformer v2, which the built-in ESM runner does not implement.");
        if (r.TryGetProperty("is_folding_model", out var fold) && fold.ValueKind == JsonValueKind.True)
            throw new NotSupportedException("ESMFold is a structure model; only ESM language models are supported.");
        return new EsmConfig
        {
            HiddenSize = Int(r, "hidden_size"),
            Layers = Int(r, "num_hidden_layers"),
            Heads = Int(r, "num_attention_heads"),
            IntermediateSize = Int(r, "intermediate_size"),
            VocabularySize = Int(r, "vocab_size"),
            LayerNormEpsilon = r.TryGetProperty("layer_norm_eps", out var eps) ? eps.GetDouble() : 1e-5,
            MaskTokenId = r.TryGetProperty("mask_token_id", out var m) && m.ValueKind == JsonValueKind.Number ? m.GetInt32() : 32,
            PadTokenId = r.TryGetProperty("pad_token_id", out var p) && p.ValueKind == JsonValueKind.Number ? p.GetInt32() : 1,
            PositionEmbeddingType = Str(r, "position_embedding_type") ?? "absolute",
            TokenDropout = r.TryGetProperty("token_dropout", out var td) && td.ValueKind == JsonValueKind.True,
            EmbeddingLayerNormBefore = r.TryGetProperty("emb_layer_norm_before", out var lnb) && lnb.ValueKind == JsonValueKind.True,
            MaxPositions = r.TryGetProperty("max_position_embeddings", out var mp) ? mp.GetInt32() : 1026,
        };
    }

    private static int Int(JsonElement r, string name)
        => r.TryGetProperty(name, out var v) ? v.GetInt32() : throw new FormatException($"config.json has no '{name}'.");

    private static string? Str(JsonElement r, string name)
        => r.TryGetProperty(name, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString() : null;
}

/// <summary>
/// A managed float32 implementation of the Hugging Face <c>EsmForMaskedLM</c> forward pass (ESM-1b, ESM-1v, ESM-2,
/// Nucleotide Transformer v1), reading its weights through HF.Net's <see cref="WeightStore"/>.
/// </summary>
/// <remarks>
/// HF.Net's <see cref="TransformerModel"/> runs BERT-family encoders only and refuses <c>model_type: esm</c>: ESM
/// is pre-layer-norm, uses rotary position embeddings and rescales embeddings for masked tokens. This class
/// implements exactly that architecture: embeddings (with token-dropout scaling) → N × [LN → self-attention with
/// rotary q/k → residual → LN → GELU feed-forward → residual] → final LN, then the LM head
/// (dense → GELU → LN → decoder tied to the embeddings).
/// </remarks>
public sealed class EsmModel : IDisposable
{
    private readonly float[] _wordEmbeddings;
    private readonly float[]? _positionEmbeddings;
    private readonly (float[] W, float[] B)? _embLayerNorm;
    private readonly Layer[] _layers;
    private readonly (float[] W, float[] B) _finalNorm;
    private readonly float[] _lmDenseW, _lmDenseB, _lmNormW, _lmNormB, _lmDecoderW, _lmBias;
    private readonly float[] _invFreq;

    private sealed record Layer(
        float[] AttnNormW, float[] AttnNormB,
        float[] Wq, float[] Bq, float[] Wk, float[] Bk, float[] Wv, float[] Bv,
        float[] Wo, float[] Bo,
        float[] FfnNormW, float[] FfnNormB,
        float[] W1, float[] B1, float[] W2, float[] B2);

    private EsmModel(string id, EsmConfig config, EsmTokenizer tokenizer, WeightStore weights)
    {
        Id = id;
        Config = config;
        Tokenizer = tokenizer;
        string prefix = weights.Contains("esm.embeddings.word_embeddings.weight") ? "esm." : "";
        float[] Read(string name) => weights.ReadFloats(name);
        float[]? TryRead(string name) => weights.Contains(name) ? weights.ReadFloats(name) : null;

        _wordEmbeddings = Read(prefix + "embeddings.word_embeddings.weight");
        if (config.PositionEmbeddingType == "absolute") _positionEmbeddings = Read(prefix + "embeddings.position_embeddings.weight");
        if (config.EmbeddingLayerNormBefore)
            _embLayerNorm = (Read(prefix + "embeddings.layer_norm.weight"), Read(prefix + "embeddings.layer_norm.bias"));

        _layers = new Layer[config.Layers];
        for (int i = 0; i < config.Layers; i++)
        {
            string l = $"{prefix}encoder.layer.{i}.";
            _layers[i] = new Layer(
                Read(l + "attention.LayerNorm.weight"), Read(l + "attention.LayerNorm.bias"),
                Read(l + "attention.self.query.weight"), Read(l + "attention.self.query.bias"),
                Read(l + "attention.self.key.weight"), Read(l + "attention.self.key.bias"),
                Read(l + "attention.self.value.weight"), Read(l + "attention.self.value.bias"),
                Read(l + "attention.output.dense.weight"), Read(l + "attention.output.dense.bias"),
                Read(l + "LayerNorm.weight"), Read(l + "LayerNorm.bias"),
                Read(l + "intermediate.dense.weight"), Read(l + "intermediate.dense.bias"),
                Read(l + "output.dense.weight"), Read(l + "output.dense.bias"));
        }

        _finalNorm = (Read(prefix + "encoder.emb_layer_norm_after.weight"), Read(prefix + "encoder.emb_layer_norm_after.bias"));

        HasLanguageModelHead = weights.Contains("lm_head.dense.weight");
        if (HasLanguageModelHead)
        {
            _lmDenseW = Read("lm_head.dense.weight");
            _lmDenseB = Read("lm_head.dense.bias");
            _lmNormW = Read("lm_head.layer_norm.weight");
            _lmNormB = Read("lm_head.layer_norm.bias");
            _lmDecoderW = TryRead("lm_head.decoder.weight") ?? _wordEmbeddings;
            _lmBias = TryRead("lm_head.bias") ?? TryRead("lm_head.decoder.bias") ?? new float[config.VocabularySize];
        }
        else
        {
            _lmDenseW = _lmDenseB = _lmNormW = _lmNormB = _lmDecoderW = _lmBias = [];
        }

        int headSize = config.HiddenSize / config.Heads;
        _invFreq = new float[headSize / 2];
        for (int j = 0; j < _invFreq.Length; j++) _invFreq[j] = (float)(1.0 / Math.Pow(10000, 2.0 * j / headSize));
    }

    /// <summary>The repository id or directory the model was loaded from.</summary>
    public string Id { get; }

    /// <summary>The configuration.</summary>
    public EsmConfig Config { get; }

    /// <summary>The tokenizer.</summary>
    public EsmTokenizer Tokenizer { get; }

    /// <summary>Whether the checkpoint has a masked-language-model head (needed for fill-mask and variant scoring).</summary>
    public bool HasLanguageModelHead { get; }

    /// <summary>Hidden size of the embeddings.</summary>
    public int HiddenSize => Config.HiddenSize;

    /// <summary>Downloads (or reads from cache) and loads an ESM checkpoint, for example <c>facebook/esm2_t6_8M_UR50D</c>.</summary>
    /// <exception cref="BioModelUnavailableException">The checkpoint is not cached and cannot be downloaded.</exception>
    /// <exception cref="NotSupportedException">The checkpoint is not a supported ESM language model.</exception>
    public static EsmModel Load(string repoIdOrDirectory, string revision = "main")
        => Open(BioModelHub.EnsureModel(repoIdOrDirectory, revision), repoIdOrDirectory);

    /// <summary>Loads an ESM checkpoint from a directory holding <c>config.json</c>, <c>vocab.txt</c> and weights.</summary>
    public static EsmModel Open(string directory, string? id = null)
    {
        var config = EsmConfig.Load(Path.Combine(directory, "config.json"));
        var tokenizer = EsmTokenizer.FromDirectory(directory);
        using var weights = WeightStore.Open(directory);
        return new EsmModel(id ?? directory, config, tokenizer, weights);
    }

    /// <summary>
    /// Runs the encoder on token ids and returns the final hidden states, row-major <c>[ids.Length, HiddenSize]</c>.
    /// </summary>
    public float[] Forward(ReadOnlySpan<int> ids)
    {
        int n = ids.Length, h = Config.HiddenSize;
        if (n == 0) throw new ArgumentException("No tokens.", nameof(ids));
        if (n > Config.MaxPositions) throw new ArgumentException($"{n} tokens exceed the model's {Config.MaxPositions} positions; chunk the sequence.", nameof(ids));
        var x = new float[n * h];
        int maskCount = 0, real = 0;
        for (int t = 0; t < n; t++)
        {
            int id = ids[t];
            if ((uint)id >= (uint)Config.VocabularySize) throw new ArgumentOutOfRangeException(nameof(ids), $"Token id {id} is outside the vocabulary.");
            if (id != Config.PadTokenId) real++;
            if (Config.TokenDropout && id == Config.MaskTokenId)
            {
                maskCount++;
                continue;
            }

            _wordEmbeddings.AsSpan(id * h, h).CopyTo(x.AsSpan(t * h, h));
        }

        if (Config.TokenDropout)
        {
            float scale = (float)((1 - 0.15 * 0.8) / (1 - (double)maskCount / Math.Max(real, 1)));
            for (int i = 0; i < x.Length; i++) x[i] *= scale;
        }

        if (_positionEmbeddings is not null)
        {
            for (int t = 0; t < n; t++)
            {
                int pos = Config.PadTokenId + 1 + t;
                var pe = _positionEmbeddings.AsSpan(pos * h, h);
                var row = x.AsSpan(t * h, h);
                for (int j = 0; j < h; j++) row[j] += pe[j];
            }
        }

        if (_embLayerNorm is { } eln) LayerNorm(x, n, h, eln.W, eln.B, Config.LayerNormEpsilon, x);

        var normed = new float[n * h];
        foreach (var layer in _layers)
        {
            LayerNorm(x, n, h, layer.AttnNormW, layer.AttnNormB, Config.LayerNormEpsilon, normed);
            var attn = Attention(normed, n, layer);
            var proj = Linear(attn, n, h, layer.Wo, layer.Bo, h);
            Add(x, proj);
            LayerNorm(x, n, h, layer.FfnNormW, layer.FfnNormB, Config.LayerNormEpsilon, normed);
            var inter = Linear(normed, n, h, layer.W1, layer.B1, Config.IntermediateSize);
            Gelu(inter);
            var ffn = Linear(inter, n, Config.IntermediateSize, layer.W2, layer.B2, h);
            Add(x, ffn);
        }

        LayerNorm(x, n, h, _finalNorm.W, _finalNorm.B, Config.LayerNormEpsilon, x);
        return x;
    }

    /// <summary>Hidden states for a sequence string (tokenised with special tokens), row-major <c>[tokens, HiddenSize]</c>.</summary>
    public float[] Hidden(string sequence, out int[] ids)
    {
        ids = Tokenizer.Encode(sequence);
        return Forward(ids);
    }

    /// <summary>Masked-LM logits over the vocabulary at one token position of the hidden states.</summary>
    /// <exception cref="InvalidOperationException">The checkpoint has no language-model head.</exception>
    public float[] LogitsAt(float[] hidden, int position)
    {
        if (!HasLanguageModelHead) throw new InvalidOperationException($"'{Id}' has no masked-language-model head.");
        int h = Config.HiddenSize;
        var row = hidden.AsSpan(position * h, h).ToArray();
        var dense = Linear(row, 1, h, _lmDenseW, _lmDenseB, h);
        Gelu(dense);
        LayerNorm(dense, 1, h, _lmNormW, _lmNormB, Config.LayerNormEpsilon, dense);
        return Linear(dense, 1, h, _lmDecoderW, _lmBias, Config.VocabularySize);
    }

    /// <summary>Log-softmax of <paramref name="logits"/>.</summary>
    public static double[] LogSoftmax(ReadOnlySpan<float> logits)
    {
        double max = double.NegativeInfinity;
        foreach (float v in logits) max = Math.Max(max, v);
        double sum = 0;
        foreach (float v in logits) sum += Math.Exp(v - max);
        double lse = max + Math.Log(sum);
        var result = new double[logits.Length];
        for (int i = 0; i < logits.Length; i++) result[i] = logits[i] - lse;
        return result;
    }

    /// <inheritdoc />
    public void Dispose()
    {
        // Weights are copied into managed arrays at load time; nothing is held open.
    }

    /// <inheritdoc />
    public override string ToString()
        => $"ESM {Id}: {Config.Layers} layers, hidden {Config.HiddenSize}, {Config.Heads} heads, vocab {Config.VocabularySize}";

    private float[] Attention(float[] x, int n, Layer layer)
    {
        int h = Config.HiddenSize, heads = Config.Heads, d = h / heads;
        var q = Linear(x, n, h, layer.Wq, layer.Bq, h);
        var k = Linear(x, n, h, layer.Wk, layer.Bk, h);
        var v = Linear(x, n, h, layer.Wv, layer.Bv, h);
        float scale = 1f / MathF.Sqrt(d);
        for (int i = 0; i < q.Length; i++) q[i] *= scale;
        if (_positionEmbeddings is null)
        {
            ApplyRotary(q, n, heads, d);
            ApplyRotary(k, n, heads, d);
        }

        var output = new float[n * h];
        Parallel.For(0, heads, head =>
        {
            var scores = new float[n];
            for (int i = 0; i < n; i++)
            {
                var qi = q.AsSpan(i * h + head * d, d);
                float max = float.NegativeInfinity;
                for (int j = 0; j < n; j++)
                {
                    float s = Dot(qi, k.AsSpan(j * h + head * d, d));
                    scores[j] = s;
                    if (s > max) max = s;
                }

                float sum = 0;
                for (int j = 0; j < n; j++)
                {
                    scores[j] = MathF.Exp(scores[j] - max);
                    sum += scores[j];
                }

                var oi = output.AsSpan(i * h + head * d, d);
                for (int j = 0; j < n; j++)
                {
                    float w = scores[j] / sum;
                    var vj = v.AsSpan(j * h + head * d, d);
                    for (int c = 0; c < d; c++) oi[c] += w * vj[c];
                }
            }
        });
        return output;
    }

    private void ApplyRotary(float[] x, int n, int heads, int d)
    {
        int half = d / 2;
        var cos = new float[half];
        var sin = new float[half];
        for (int t = 0; t < n; t++)
        {
            for (int j = 0; j < half; j++)
            {
                float f = t * _invFreq[j];
                cos[j] = MathF.Cos(f);
                sin[j] = MathF.Sin(f);
            }

            for (int head = 0; head < heads; head++)
            {
                var v = x.AsSpan(t * heads * d + head * d, d);
                for (int j = 0; j < half; j++)
                {
                    float a = v[j], b = v[j + half];
                    v[j] = a * cos[j] - b * sin[j];
                    v[j + half] = b * cos[j] + a * sin[j];
                }
            }
        }
    }

    private static float[] Linear(float[] x, int rows, int inDim, float[] w, float[] b, int outDim)
    {
        var y = new float[rows * outDim];
        Parallel.For(0, outDim, new ParallelOptions { MaxDegreeOfParallelism = outDim < 64 ? 1 : -1 }, o =>
        {
            var wo = w.AsSpan(o * inDim, inDim);
            float bias = b.Length > 0 ? b[o] : 0f;
            for (int r = 0; r < rows; r++) y[r * outDim + o] = Dot(x.AsSpan(r * inDim, inDim), wo) + bias;
        });
        return y;
    }

    private static float Dot(ReadOnlySpan<float> a, ReadOnlySpan<float> b)
    {
        int i = 0;
        float sum = 0;
        if (Vector.IsHardwareAccelerated && a.Length >= Vector<float>.Count)
        {
            var va = MemoryMarshal.Cast<float, Vector<float>>(a);
            var vb = MemoryMarshal.Cast<float, Vector<float>>(b);
            var acc = Vector<float>.Zero;
            for (int j = 0; j < va.Length; j++) acc += va[j] * vb[j];
            sum = Vector.Sum(acc);
            i = va.Length * Vector<float>.Count;
        }

        for (; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    private static void Add(float[] x, float[] y)
    {
        for (int i = 0; i < x.Length; i++) x[i] += y[i];
    }

    private static void LayerNorm(float[] x, int rows, int dim, float[] gamma, float[] beta, double eps, float[] output)
    {
        for (int r = 0; r < rows; r++)
        {
            var row = x.AsSpan(r * dim, dim);
            double mean = 0;
            for (int i = 0; i < dim; i++) mean += row[i];
            mean /= dim;
            double var = 0;
            for (int i = 0; i < dim; i++)
            {
                double dv = row[i] - mean;
                var += dv * dv;
            }

            float inv = (float)(1.0 / Math.Sqrt(var / dim + eps));
            float m = (float)mean;
            var outRow = output.AsSpan(r * dim, dim);
            for (int i = 0; i < dim; i++) outRow[i] = (row[i] - m) * inv * gamma[i] + beta[i];
        }
    }

    private static void Gelu(float[] x)
    {
        for (int i = 0; i < x.Length; i++)
        {
            float v = x[i];
            x[i] = 0.5f * v * (1f + Erf(v * 0.70710678118654752f));
        }
    }

    // Abramowitz & Stegun 7.1.26 refined (max abs error ~1.2e-7), adequate for float32 activations.
    private static float Erf(float x)
    {
        double ax = Math.Abs(x);
        double t = 1.0 / (1.0 + 0.5 * ax);
        double y = 1 - t * Math.Exp(-ax * ax - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
            t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
        return (float)(x >= 0 ? y : -y);
    }
}
