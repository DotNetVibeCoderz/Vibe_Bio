using Biotech.Net.AI.Copilot;
using Biotech.Net.AI.Hub;
using Biotech.Net.AI.Inference;
using Biotech.Net.AI.Models;
using Biotech.Net.AI.Preprocessing;

namespace Biotech.Net.AI.Tests;

/// <summary>
/// Credentials for real-model tests, read at runtime from environment variables or the developer's credential files
/// outside the repository. Keys are never printed.
/// </summary>
internal static class TestCredentials
{
    private const string CredentialDir = @"C:\Users\mifma\Documents\CodeSandbox";

    public static string? HuggingFaceToken()
    {
        if (Environment.GetEnvironmentVariable("HF_TOKEN") is { Length: > 0 } t) return t;
        string file = Environment.GetEnvironmentVariable("BIOTECHNET_HF_TOKEN_FILE") ?? Path.Combine(CredentialDir, "HFToken.txt");
        return File.Exists(file) ? File.ReadAllText(file).Trim() : null;
    }

    /// <summary>Values of the "label: value" lines in one section of the multi-provider key file.</summary>
    private static Dictionary<string, string>? Section(string header)
    {
        string file = Environment.GetEnvironmentVariable("BIOTECHNET_LLM_KEY_FILE") ?? Path.Combine(CredentialDir, "testkey.txt");
        if (!File.Exists(file)) return null;
        Dictionary<string, string>? current = null;
        foreach (var raw in File.ReadAllLines(file))
        {
            string line = raw.Trim();
            if (line.Length == 0)
            {
                if (current is not null) return current;
                continue;
            }

            int colon = line.IndexOf(':');
            bool isUrlLine = colon > 0 && line.Length > colon + 1 && line[colon + 1] == '/';
            if (colon < 0 || isUrlLine)
            {
                if (current is not null) return current;
                if (line.Equals(header, StringComparison.OrdinalIgnoreCase)) current = new(StringComparer.OrdinalIgnoreCase);
                continue;
            }

            current?.TryAdd(line[..colon].Trim().Replace(" ", ""), line[(colon + 1)..].Trim());
        }

        return current;
    }

    public static AnthropicChatClient? Anthropic()
    {
        string? key = Environment.GetEnvironmentVariable("ANTHROPIC_API_KEY");
        if (string.IsNullOrEmpty(key)) key = Section("Anthoripic")?.GetValueOrDefault("ApiKey") ?? Section("Anthropic")?.GetValueOrDefault("ApiKey");
        return ChatClientFactory.DetectProvider(key) == ChatProvider.Anthropic
            ? new AnthropicChatClient(key, Environment.GetEnvironmentVariable("ANTHROPIC_MODEL") ?? AnthropicChatClient.DefaultModel)
            : null;
    }

    public static OpenAiCompatibleChatClient? AzureOpenAI()
    {
        string? key = Environment.GetEnvironmentVariable("AZURE_OPENAI_API_KEY");
        string? endpoint = Environment.GetEnvironmentVariable("AZURE_OPENAI_ENDPOINT");
        string? deployment = Environment.GetEnvironmentVariable("AZURE_OPENAI_DEPLOYMENT");
        if (string.IsNullOrEmpty(key) || string.IsNullOrEmpty(endpoint) || string.IsNullOrEmpty(deployment))
        {
            var s = Section("Azure OpenAI");
            key = s?.GetValueOrDefault("apikey");
            endpoint = s?.GetValueOrDefault("endpoint");
            deployment = s?.GetValueOrDefault("model")?.Split(',')[0].Trim();
        }

        return string.IsNullOrEmpty(key) || string.IsNullOrEmpty(endpoint) || string.IsNullOrEmpty(deployment)
            ? null
            : OpenAiCompatibleChatClient.ForAzureOpenAI(endpoint, key, deployment);
    }
}

/// <summary>Tests against real checkpoints and APIs; each skips itself when credentials or the network are missing.</summary>
[Collection(HubCollection.Name)]
public class RealModelTests
{
    private const string Ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    private static void RequireHub()
    {
        string? token = TestCredentials.HuggingFaceToken();
        if (token is null) Assert.Skip("No Hugging Face token (HF_TOKEN or HFToken.txt).");
        BioModelHub.Configure(new BioHubOptions { Token = token });
    }

    private static T OrSkip<T>(Func<T> load)
    {
        try
        {
            return load();
        }
        catch (BioModelUnavailableException ex)
        {
            Assert.Skip("Model unavailable (offline?): " + ex.Message);
            throw;
        }
    }

    [Fact]
    public void Esm2T6_EmbedsFillsMaskAndScoresVariants()
    {
        RequireHub();
        using var embedder = OrSkip(() => BioEmbedder.Load("Esm2T6"));
        var model = ((EsmEncoder)embedder.Encoder).Model;

        // Mean-pooled embedding; the first values were checked against transformers' EsmForMaskedLM (float32, ~1e-4).
        var v = embedder.Embed(Ubiquitin);
        Assert.Equal(320, v.Length);
        Assert.Equal(-0.18984, v[0], 3);
        Assert.Equal(0.30577, v[4], 3);

        var batch = embedder.EmbedBatch([Ubiquitin, Ubiquitin[..40], "MKTAYIAKQRQISFVKSHFSRQ"]);
        Assert.Equal(v, batch[0]);
        Assert.True(BioEmbedder.CosineSimilarity(batch[0], batch[1]) > BioEmbedder.CosineSimilarity(batch[0], batch[2]));

        // Masked marginal for K11P: reference log P(K) = -1.53276, log P(P) = -3.62903 at the masked site.
        var scorer = new EsmVariantScorer(model);
        var score = scorer.Score(Ubiquitin, "K11P");
        Assert.Equal(-2.09628, score.Score, 2);

        var scan = scorer.SaturationScan(Ubiquitin, 1, 20);
        double mean = scan.Scores.Cast<double>().Average();
        Assert.True(mean < 0, $"Mean substitution score {mean} should be negative for a natural protein.");

        var fill = new BioMaskedLanguageModel(model).PredictResidue(Ubiquitin, 1, 3);
        Assert.Equal("M", fill[0].Token);
    }

    [Fact]
    public void GeneformerV1_EmbedsARankEncodedCell()
    {
        RequireHub();
        var (tokens, medians) = OrSkip(() => GeneformerDictionaries.FromHub(GeneformerDictionaryVersion.V1Gc30M));
        Assert.True(tokens.Count > 20_000);
        Assert.True(medians.Count > 20_000);
        using var model = OrSkip(() => GeneformerModel.Load(BioModelCatalog.GeneformerV1));

        var encoder = new GeneformerEncoder(tokens, medians, GeneformerEncoderOptions.V1);
        var genes = medians.Keys.Where(tokens.ContainsKey).Take(300).ToArray();
        var counts = genes.Select((_, i) => (double)(i % 17)).ToArray();
        var ids = encoder.Encode(genes, counts);
        Assert.InRange(ids.Length, 250, 300);
        var embedding = model.EmbedCell(ids);
        Assert.Equal(256, embedding.Length);
        Assert.All(embedding, x => Assert.True(float.IsFinite(x)));
    }

    private static async Task AskWithTools(IChatClient? client, string provider)
    {
        if (client is null) Assert.Skip($"No {provider} credentials (environment or key file).");
        var copilot = new BioCopilot(client!, options: new BioCopilotOptions { MaxTokens = 4096 });
        CopilotAnswer answer;
        try
        {
            answer = await copilot.AskAsync(
                "Using your tools, give the reverse complement and the GC content (percent) of ATGCGCGTTAAC.",
                TestContext.Current.CancellationToken);
        }
        catch (ChatClientException ex) when (ex.StatusCode is null)
        {
            Assert.Skip($"{provider} unreachable: {ex.Message}");
            throw;
        }
        catch (ChatClientException ex) when (ex.StatusCode is 401 or 403
            || ex.StatusCode == 400 && ex.Message.Contains("workspace", StringComparison.OrdinalIgnoreCase))
        {
            Assert.Skip($"{provider} credentials are not usable here (HTTP {ex.StatusCode}; e.g. set ANTHROPIC_WORKSPACE_ID): {ex.Message}");
            throw;
        }

        TestContext.Current.TestOutputHelper?.WriteLine($"[{provider} {client!.Model}] rounds={answer.Rounds} tools={string.Join(",", answer.ToolCalls.Select(t => t.Name))}\n{answer.Text}");
        Assert.Contains(answer.ToolCalls, t => t.Name == "reverse_complement" && !t.IsError);
        Assert.Contains(answer.ToolCalls, t => t.Name == "gc_content" && !t.IsError);
        Assert.Contains("GTTAACGCGCAT", answer.Text);
        Assert.Contains("50", answer.Text);
    }

    [Fact]
    public Task Copilot_Anthropic_UsesTools() => AskWithTools(TestCredentials.Anthropic(), "Anthropic");

    [Fact]
    public Task Copilot_AzureOpenAI_UsesTools() => AskWithTools(TestCredentials.AzureOpenAI(), "Azure OpenAI");
}
