using System.Text;
using System.Text.Json;
using Biotech.Net.AI.Hub;
using Biotech.Net.AI.Inference;
using Biotech.Net.AI.Models;

namespace Biotech.Net.AI.Tests;

public class BioModelCatalogTests
{
    [Fact]
    public void KeysAndRepoIdsAreWellFormed()
    {
        Assert.Equal(BioModelCatalog.All.Count, BioModelCatalog.All.Select(m => m.Key.ToLowerInvariant()).Distinct().Count());
        Assert.All(BioModelCatalog.All, m => Assert.Matches(@"^[\w.-]+/[\w.-]+$", m.RepoId));
        Assert.Contains(BioModelCatalog.All, m => m.RepoId == "facebook/esm2_t6_8M_UR50D");
        Assert.Contains(BioModelCatalog.All, m => m.RepoId == "dmis-lab/biobert-base-cased-v1.2");
        Assert.Contains(BioModelCatalog.All, m => m.RepoId == "zhihan1996/DNA_bert_6");
        Assert.Contains(BioModelCatalog.All, m => m.RepoId == "Rostlab/prot_bert");
        Assert.Contains(BioModelCatalog.All, m => m.RepoId == "ctheodoris/Geneformer");
    }

    [Fact]
    public void Get_ByKeyOrRepoId_CaseInsensitive()
    {
        Assert.Same(BioModelCatalog.Esm2T6, BioModelCatalog.Get("esm2t6"));
        Assert.Same(BioModelCatalog.ProtBert, BioModelCatalog.Get("ROSTLAB/PROT_BERT"));
        Assert.Throws<KeyNotFoundException>(() => BioModelCatalog.Get("nope"));
    }

    [Fact]
    public void Smallest_PicksRunnableModels()
    {
        Assert.Same(BioModelCatalog.Esm2T6, BioModelCatalog.Smallest(BioSequenceKind.Protein));
        Assert.Same(BioModelCatalog.DnaBert6, BioModelCatalog.Smallest(BioSequenceKind.Dna));
        Assert.Same(BioModelCatalog.GeneformerV1, BioModelCatalog.Smallest(BioSequenceKind.SingleCellExpression));
    }

    [Fact]
    public void Preprocess_DispatchesOnInputConvention()
    {
        Assert.Equal("ATGCAT TGCATG", BioModelCatalog.DnaBert6.Preprocess("atgcatg"));
        Assert.Equal("M K T X", BioModelCatalog.ProtBert.Preprocess("mktu"));
        Assert.Equal("MKTU", BioModelCatalog.Esm2T6.Preprocess("mktu"));
        Assert.Equal("BRCA1 text", BioModelCatalog.PubMedBert.Preprocess("BRCA1 text"));
        Assert.Throws<NotSupportedException>(() => BioModelCatalog.GeneformerV1.Preprocess("x"));
    }

    [Fact]
    public void UnsupportedArchitectures_AreRefusedWithoutDownloading()
    {
        var ex = Assert.Throws<NotSupportedException>(() => BioEmbedder.Load("DnaBert2"));
        Assert.Contains("DNABERT-2", ex.Message);
        Assert.Throws<NotSupportedException>(() => BioEmbedder.Load("NucleotideTransformerV2_50m"));
        Assert.Throws<NotSupportedException>(() => BioEmbedder.Load("GeneformerV1"));
    }
}

public class EsmTokenizerTests
{
    internal static readonly string[] Vocab =
    [
        "<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E", "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y",
        "M", "H", "W", "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>", "<mask>",
    ];

    [Fact]
    public void Encode_AddsClsEosAndRecognisesMask()
    {
        var tok = new EsmTokenizer(Vocab);
        Assert.Equal([0, 20, 15, 32, 5, 2], tok.Encode("MK<mask>A"));
        Assert.Equal([20, 3], tok.Encode("M#", addSpecialTokens: false));
        Assert.Equal("MKA", tok.Decode(tok.Encode("MKA")));
        Assert.True(tok.IsSpecial(tok.MaskId));
    }

    [Fact]
    public void Encode_GreedyKmersLikeNucleotideTransformer()
    {
        var tok = new EsmTokenizer(["<unk>", "<pad>", "<mask>", "<cls>", "AAAAAA", "CCCCCC", "A", "C", "G", "T", "N"], appendEos: false);
        Assert.Equal([3, 4, 5, 8, 9], tok.Encode("AAAAAACCCCCCGT"));
        Assert.Equal(-1, tok.EosId);
    }
}

public class EsmModelTests
{
    // A tiny random ESM checkpoint (hidden 8, 2 heads, 1 layer) written as safetensors, to exercise the forward pass offline.
    private static string WriteTinyCheckpoint(bool rotary = true)
    {
        string dir = Path.Combine(Path.GetTempPath(), "biotechnet-esm-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        int h = 8, ff = 16, v = EsmTokenizerTests.Vocab.Length;
        File.WriteAllLines(Path.Combine(dir, "vocab.txt"), EsmTokenizerTests.Vocab);
        File.WriteAllText(Path.Combine(dir, "config.json"), JsonSerializer.Serialize(new Dictionary<string, object>
        {
            ["model_type"] = "esm", ["hidden_size"] = h, ["num_hidden_layers"] = 1, ["num_attention_heads"] = 2,
            ["intermediate_size"] = ff, ["vocab_size"] = v, ["layer_norm_eps"] = 1e-5, ["mask_token_id"] = 32, ["pad_token_id"] = 1,
            ["position_embedding_type"] = rotary ? "rotary" : "absolute", ["token_dropout"] = true,
            ["emb_layer_norm_before"] = false, ["max_position_embeddings"] = 64,
        }));
        var rng = new Random(7);
        var tensors = new List<(string Name, int[] Shape)>
        {
            ("esm.embeddings.word_embeddings.weight", [v, h]),
            ("esm.embeddings.position_embeddings.weight", [64, h]),
            ("esm.encoder.layer.0.attention.LayerNorm.weight", [h]), ("esm.encoder.layer.0.attention.LayerNorm.bias", [h]),
            ("esm.encoder.layer.0.attention.self.query.weight", [h, h]), ("esm.encoder.layer.0.attention.self.query.bias", [h]),
            ("esm.encoder.layer.0.attention.self.key.weight", [h, h]), ("esm.encoder.layer.0.attention.self.key.bias", [h]),
            ("esm.encoder.layer.0.attention.self.value.weight", [h, h]), ("esm.encoder.layer.0.attention.self.value.bias", [h]),
            ("esm.encoder.layer.0.attention.output.dense.weight", [h, h]), ("esm.encoder.layer.0.attention.output.dense.bias", [h]),
            ("esm.encoder.layer.0.LayerNorm.weight", [h]), ("esm.encoder.layer.0.LayerNorm.bias", [h]),
            ("esm.encoder.layer.0.intermediate.dense.weight", [ff, h]), ("esm.encoder.layer.0.intermediate.dense.bias", [ff]),
            ("esm.encoder.layer.0.output.dense.weight", [h, ff]), ("esm.encoder.layer.0.output.dense.bias", [h]),
            ("esm.encoder.emb_layer_norm_after.weight", [h]), ("esm.encoder.emb_layer_norm_after.bias", [h]),
            ("lm_head.dense.weight", [h, h]), ("lm_head.dense.bias", [h]),
            ("lm_head.layer_norm.weight", [h]), ("lm_head.layer_norm.bias", [h]), ("lm_head.bias", [v]),
        };
        var header = new Dictionary<string, object>();
        var data = new MemoryStream();
        foreach (var (name, shape) in tensors)
        {
            int count = shape.Aggregate(1, (a, b) => a * b);
            long begin = data.Length;
            for (int i = 0; i < count; i++)
            {
                float value = name.EndsWith("LayerNorm.weight") || name.EndsWith("layer_norm.weight") || name.EndsWith("norm_after.weight")
                    ? 1f : (float)(rng.NextDouble() - 0.5);
                data.Write(BitConverter.GetBytes(value));
            }

            header[name] = new { dtype = "F32", shape, data_offsets = new[] { begin, data.Length } };
        }

        var headerBytes = Encoding.UTF8.GetBytes(JsonSerializer.Serialize(header));
        using var file = File.Create(Path.Combine(dir, "model.safetensors"));
        file.Write(BitConverter.GetBytes((long)headerBytes.Length));
        file.Write(headerBytes);
        data.Position = 0;
        data.CopyTo(file);
        return dir;
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Forward_ProducesFiniteStatesAndNormalisedLogits(bool rotary)
    {
        string dir = WriteTinyCheckpoint(rotary);
        try
        {
            using var model = EsmModel.Open(dir);
            Assert.True(model.HasLanguageModelHead);
            var ids = model.Tokenizer.Encode("MKTAYIAK");
            var hidden = model.Forward(ids);
            Assert.Equal(ids.Length * 8, hidden.Length);
            Assert.All(hidden, x => Assert.True(float.IsFinite(x)));
            Assert.Equal(hidden, model.Forward(ids)); // deterministic

            var logp = EsmModel.LogSoftmax(model.LogitsAt(hidden, 3));
            Assert.Equal(33, logp.Length);
            Assert.Equal(1.0, logp.Sum(Math.Exp), 6);

            var masked = (int[])ids.Clone();
            masked[3] = model.Tokenizer.MaskId;
            Assert.NotEqual(hidden, model.Forward(masked));
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }

    [Fact]
    public void VariantScorer_ScoresAndValidatesWildType()
    {
        string dir = WriteTinyCheckpoint();
        try
        {
            using var model = EsmModel.Open(dir);
            var scorer = new EsmVariantScorer(model);
            const string wt = "MKTAYIAKQR";
            Assert.Equal(0, scorer.Score(wt, "K2K").Score, 12);
            var multi = scorer.Score(wt, "K2A:T3G");
            Assert.Equal(multi.PerSite.Sum(), multi.Score, 12);
            Assert.Equal(scorer.Score(wt, ProteinMutation.Parse("K2A")), multi.PerSite[0], 12);
            Assert.Throws<ArgumentException>(() => scorer.Score(wt, "A2G"));
            Assert.Throws<ArgumentException>(() => scorer.Score(wt, "A99G"));

            var scan = scorer.SaturationScan(wt, 2, 4);
            Assert.Equal([2, 3, 4], scan.Positions);
            Assert.Equal(0, scan.Scores[0, EsmVariantScorer.StandardAminoAcids.IndexOf('K')], 12);

            var fill = new BioMaskedLanguageModel(model).PredictResidue(wt, 2, 3);
            Assert.Equal(3, fill.Count);
            Assert.True(fill[0].Score >= fill[1].Score);
            Assert.Equal(10, fill[0].Sequence.Length);
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }

    [Fact]
    public void Config_RefusesNucleotideTransformerV2Block()
    {
        string path = Path.GetTempFileName();
        try
        {
            File.WriteAllText(path, """{"model_type":"esm","add_bias_fnn":false,"hidden_size":8}""");
            Assert.Throws<NotSupportedException>(() => EsmConfig.Load(path));
            File.WriteAllText(path, """{"model_type":"bert"}""");
            Assert.Throws<NotSupportedException>(() => EsmConfig.Load(path));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Theory]
    [InlineData("A24G", 'A', 24, 'G')]
    [InlineData("p.k30r", 'K', 30, 'R')]
    public void ProteinMutation_Parses(string text, char wt, int pos, char mt)
    {
        var m = ProteinMutation.Parse(text);
        Assert.Equal(new ProteinMutation(wt, pos, mt), m);
        Assert.Equal(2, ProteinMutation.ParseMany("A1G:K2R").Count);
        Assert.Throws<FormatException>(() => ProteinMutation.Parse("24G"));
    }

    [Fact]
    public void Pooling_MeanSkipsSpecialTokens()
    {
        var states = new TokenStates([1, 1, 3, 3, 5, 5, 100, 100], 4, 2, [true, false, false, true]);
        var mean = new BioEmbedder(new FakeEncoder(states)).Pool(states);
        Assert.Equal([4f, 4f], mean);
        var cls = new BioEmbedder(new FakeEncoder(states)) { Pooling = PoolingStrategy.Cls }.Pool(states);
        Assert.Equal([1f, 1f], cls);
        Assert.Equal(1.0, BioEmbedder.CosineSimilarity([1f, 2f], [2f, 4f]), 12);
    }

    [Fact]
    public void Embed_WindowsLongSequences()
    {
        var encoder = new CountingEncoder();
        var embedder = new BioEmbedder(encoder, BioModelCatalog.Esm2T6 with { MaxTokens = 12 }) { WindowOverlap = 2 };
        var v = embedder.Embed(new string('A', 25));
        Assert.True(encoder.Calls > 1);
        Assert.All(encoder.Inputs, s => Assert.True(s.Length <= 10));
        Assert.Equal(2, v.Length);
        var batch = embedder.EmbedBatch(["MKT", "MKTA"]);
        Assert.Equal(2, batch.Length);
    }

    private sealed class FakeEncoder(TokenStates states) : IBioEncoder
    {
        public string ModelId => "fake";
        public int HiddenSize => states.HiddenSize;
        public int MaxTokens => 512;
        public TokenStates Encode(string modelInput) => states;
        public void Dispose() { }
    }

    private sealed class CountingEncoder : IBioEncoder
    {
        public int Calls;
        public System.Collections.Concurrent.ConcurrentBag<string> Inputs { get; } = [];
        public string ModelId => "counting";
        public int HiddenSize => 2;
        public int MaxTokens => 1024;

        public TokenStates Encode(string modelInput)
        {
            Interlocked.Increment(ref Calls);
            Inputs.Add(modelInput);
            int n = modelInput.Length + 2;
            var special = new bool[n];
            special[0] = special[^1] = true;
            return new TokenStates(Enumerable.Repeat(1f, n * 2).ToArray(), n, 2, special);
        }

        public void Dispose() { }
    }
}

[Collection(HubCollection.Name)]
public class BioModelHubOfflineTests
{
    [Fact]
    public void OfflineMode_UncachedModel_ThrowsClearException()
    {
        string cache = Path.Combine(Path.GetTempPath(), "biotechnet-emptycache-" + Guid.NewGuid().ToString("N"));
        var previous = BioModelHub.Options;
        try
        {
            BioModelHub.Configure(new BioHubOptions { Offline = true, CacheDirectory = cache });
            Assert.False(BioModelHub.IsCached("facebook/esm2_t6_8M_UR50D"));
            var ex = Assert.Throws<BioModelUnavailableException>(() => BioModelHub.EnsureModel("facebook/esm2_t6_8M_UR50D"));
            Assert.Contains("offline", ex.Message, StringComparison.OrdinalIgnoreCase);
            Assert.Equal("facebook/esm2_t6_8M_UR50D", ex.RepoId);
            Assert.Throws<BioModelUnavailableException>(() => BioModelHub.DownloadFile("ctheodoris/Geneformer", "x.pkl"));
            Assert.Throws<BioModelUnavailableException>(() => BioEmbedder.Load("Esm2T6"));
        }
        finally
        {
            BioModelHub.Configure(previous);
        }
    }
}

[CollectionDefinition(Name, DisableParallelization = true)]
public sealed class HubCollection
{
    public const string Name = "Hub configuration";
}
