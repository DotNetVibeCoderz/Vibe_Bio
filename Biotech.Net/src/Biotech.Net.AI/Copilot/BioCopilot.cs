using Biotech.Net.Core;

namespace Biotech.Net.AI.Copilot;

/// <summary>One tool call made while answering.</summary>
/// <param name="Name">Tool name.</param>
/// <param name="ArgumentsJson">Arguments the model passed.</param>
/// <param name="Result">What the tool returned.</param>
/// <param name="IsError">Whether the tool failed.</param>
public sealed record ToolInvocation(string Name, string ArgumentsJson, string Result, bool IsError);

/// <summary>The copilot's answer.</summary>
/// <param name="Text">The final answer text.</param>
/// <param name="ToolCalls">Every tool call made, in order.</param>
/// <param name="Rounds">Model calls made.</param>
/// <param name="StopReason">Why the last model call stopped (<c>end_turn</c>, <c>max_tokens</c>, <c>refusal</c>, <c>max_rounds</c>…).</param>
public sealed record CopilotAnswer(string Text, IReadOnlyList<ToolInvocation> ToolCalls, int Rounds, string StopReason);

/// <summary>Settings for <see cref="BioCopilot"/>.</summary>
public sealed record BioCopilotOptions
{
    /// <summary>The default system prompt.</summary>
    public const string DefaultSystemPrompt =
        "You are the Biotech.Net bioinformatics copilot. Answer questions about molecular biology, genomics, proteomics and " +
        "bioinformatics clearly and accurately. When a question involves computing something about a concrete sequence " +
        "(reverse complement, translation, GC content, protein properties, open reading frames), call the provided tools " +
        "instead of computing by hand, and base your answer on their results. Say when something is uncertain. " +
        "Use 1-based coordinates when reporting positions.";

    /// <summary>System prompt.</summary>
    public string SystemPrompt { get; init; } = DefaultSystemPrompt;

    /// <summary>Maximum model calls per question (tool rounds included).</summary>
    public int MaxRounds { get; init; } = 8;

    /// <summary>Maximum output tokens per model call.</summary>
    public int MaxTokens { get; init; } = 8192;
}

/// <summary>
/// An LLM copilot for bioinformatics that answers questions and runs Biotech.Net tools through a tool-use loop over any
/// <see cref="IChatClient"/> (Anthropic, Azure OpenAI, OpenAI, DeepSeek, …). The conversation is kept between questions.
/// </summary>
/// <example>
/// <code>
/// var copilot = new BioCopilot(new AnthropicChatClient());   // key from ANTHROPIC_API_KEY
/// var answer = await copilot.AskAsync("What is the GC content of ATGCGCGTTA?");
/// </code>
/// </example>
public sealed class BioCopilot
{
    private readonly List<ChatMessage> _history = [];

    /// <summary>Creates a copilot.</summary>
    public BioCopilot(IChatClient client, BioToolbox? toolbox = null, BioCopilotOptions? options = null)
    {
        Client = client ?? throw new ArgumentNullException(nameof(client));
        Toolbox = toolbox ?? new BioToolbox();
        Options = options ?? new BioCopilotOptions();
    }

    /// <summary>The chat model.</summary>
    public IChatClient Client { get; }

    /// <summary>The tools.</summary>
    public BioToolbox Toolbox { get; }

    /// <summary>The settings.</summary>
    public BioCopilotOptions Options { get; }

    /// <summary>The conversation so far.</summary>
    public IReadOnlyList<ChatMessage> History => _history;

    /// <summary>Credits line shown by apps that host the copilot.</summary>
    public static string Credits => BiotechInfo.CreditsId;

    /// <summary>Forgets the conversation.</summary>
    public void Reset() => _history.Clear();

    /// <summary>Asks a question, running tool calls until the model answers.</summary>
    /// <exception cref="ChatClientException">The chat API failed.</exception>
    public async Task<CopilotAnswer> AskAsync(string question, CancellationToken cancellationToken = default)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(question);
        _history.Add(ChatMessage.User(question));
        var invocations = new List<ToolInvocation>();
        var tools = Toolbox.Definitions;
        string stop = "max_rounds";
        string text = "";
        int rounds = 0;
        while (rounds < Options.MaxRounds)
        {
            rounds++;
            var completion = await Client.CompleteAsync(new ChatRequest
            {
                SystemPrompt = Options.SystemPrompt,
                Messages = _history,
                Tools = tools,
                MaxTokens = Options.MaxTokens,
            }, cancellationToken).ConfigureAwait(false);

            _history.Add(completion.Message);
            text = completion.Message.Text;
            stop = completion.StopReason;
            if (completion.Message.ToolCalls.Count == 0) break;

            foreach (var call in completion.Message.ToolCalls)
            {
                var result = Toolbox.Invoke(call.Name, call.ArgumentsJson);
                invocations.Add(new ToolInvocation(call.Name, call.ArgumentsJson, result.Content, result.IsError));
                _history.Add(ChatMessage.ToolResult(call.Id, call.Name, result.Content, result.IsError));
            }

            stop = "max_rounds";
        }

        return new CopilotAnswer(text, invocations, rounds, stop);
    }
}

/// <summary>Known chat providers.</summary>
public enum ChatProvider
{
    /// <summary>Unrecognised.</summary>
    Unknown,

    /// <summary>Anthropic Messages API.</summary>
    Anthropic,

    /// <summary>OpenAI.</summary>
    OpenAI,

    /// <summary>Azure OpenAI.</summary>
    AzureOpenAI,

    /// <summary>DeepSeek.</summary>
    DeepSeek,
}

/// <summary>Creates chat clients from the environment. Keys are only read from variables at runtime, never stored.</summary>
public static class ChatClientFactory
{
    /// <summary>
    /// Guesses a key's provider from its prefix only: <c>sk-ant-</c> Anthropic; <c>sk-proj-</c>/<c>sk-svcacct-</c> OpenAI;
    /// other <c>sk-</c> keys OpenAI-compatible (OpenAI or DeepSeek — returned as <see cref="ChatProvider.OpenAI"/>);
    /// a long key without a prefix is typical of Azure OpenAI.
    /// </summary>
    public static ChatProvider DetectProvider(string? apiKey)
    {
        if (string.IsNullOrWhiteSpace(apiKey)) return ChatProvider.Unknown;
        string k = apiKey.Trim();
        if (k.StartsWith("sk-ant-", StringComparison.Ordinal)) return ChatProvider.Anthropic;
        if (k.StartsWith("sk-", StringComparison.Ordinal)) return ChatProvider.OpenAI;
        return k.Length >= 32 && k.All(char.IsAsciiLetterOrDigit) ? ChatProvider.AzureOpenAI : ChatProvider.Unknown;
    }

    /// <summary>
    /// The first configured provider among: <c>ANTHROPIC_API_KEY</c> (+ optional <c>ANTHROPIC_MODEL</c>);
    /// <c>AZURE_OPENAI_API_KEY</c> + <c>AZURE_OPENAI_ENDPOINT</c> + <c>AZURE_OPENAI_DEPLOYMENT</c>;
    /// <c>OPENAI_API_KEY</c> (+ <c>OPENAI_MODEL</c>, default <c>gpt-5-mini</c>);
    /// <c>DEEPSEEK_API_KEY</c> (+ <c>DEEPSEEK_MODEL</c>, default <c>deepseek-chat</c>). Null when none is set.
    /// </summary>
    public static IChatClient? FromEnvironment()
    {
        static string? Env(string name) => Environment.GetEnvironmentVariable(name) is { Length: > 0 } v ? v : null;
        if (Env("ANTHROPIC_API_KEY") is { } a) return new AnthropicChatClient(a, Env("ANTHROPIC_MODEL") ?? AnthropicChatClient.DefaultModel);
        if (Env("AZURE_OPENAI_API_KEY") is { } az && Env("AZURE_OPENAI_ENDPOINT") is { } ep && Env("AZURE_OPENAI_DEPLOYMENT") is { } dep)
            return OpenAiCompatibleChatClient.ForAzureOpenAI(ep, az, dep);
        if (Env("OPENAI_API_KEY") is { } o) return OpenAiCompatibleChatClient.ForOpenAI(o, Env("OPENAI_MODEL") ?? "gpt-5-mini");
        if (Env("DEEPSEEK_API_KEY") is { } d) return OpenAiCompatibleChatClient.ForDeepSeek(d, Env("DEEPSEEK_MODEL") ?? "deepseek-chat");
        return null;
    }
}
