using System.Net.Http.Headers;
using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;

namespace Biotech.Net.AI.Copilot;

/// <summary>
/// <see cref="IChatClient"/> over the Anthropic Messages API (<c>POST /v1/messages</c>) with client-side tool use.
/// </summary>
/// <remarks>
/// The assistant's content blocks (including thinking blocks) are kept in <see cref="ChatMessage.ProviderContent"/>
/// and replayed unchanged, as the API requires within a tool-use loop. The API key is supplied by the caller or read
/// from <c>ANTHROPIC_API_KEY</c>; it is never logged.
/// </remarks>
public sealed class AnthropicChatClient : IChatClient, IDisposable
{
    /// <summary>The default model.</summary>
    public const string DefaultModel = "claude-opus-5-5";

    private readonly HttpClient _http;
    private readonly bool _ownsHttp;
    private readonly string _apiKey;

    /// <summary>Creates a client.</summary>
    /// <param name="apiKey">API key; null reads <c>ANTHROPIC_API_KEY</c>.</param>
    /// <param name="model">Model id, for example <c>claude-opus-5-5</c> or <c>claude-sonnet-5-5</c>.</param>
    /// <param name="httpClient">Optional HTTP client (for proxies or tests).</param>
    /// <exception cref="InvalidOperationException">No key was given and <c>ANTHROPIC_API_KEY</c> is unset.</exception>
    public AnthropicChatClient(string? apiKey = null, string model = DefaultModel, HttpClient? httpClient = null)
    {
        _apiKey = apiKey ?? Environment.GetEnvironmentVariable("ANTHROPIC_API_KEY")
            ?? throw new InvalidOperationException("No Anthropic API key: pass one or set ANTHROPIC_API_KEY.");
        Model = model;
        _ownsHttp = httpClient is null;
        _http = httpClient ?? new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
    }

    /// <inheritdoc />
    public string Provider => "anthropic";

    /// <inheritdoc />
    public string Model { get; }

    /// <summary>API base URL. Defaults to <c>https://api.anthropic.com</c>.</summary>
    public Uri BaseUri { get; init; } = new("https://api.anthropic.com");

    /// <summary>
    /// Workspace id sent as <c>anthropic-workspace-id</c>, required for API keys not scoped to a workspace.
    /// Defaults to the <c>ANTHROPIC_WORKSPACE_ID</c> environment variable.
    /// </summary>
    public string? WorkspaceId { get; init; } = Environment.GetEnvironmentVariable("ANTHROPIC_WORKSPACE_ID") is { Length: > 0 } w ? w : null;

    /// <summary>Optional <c>output_config.effort</c> (<c>low</c>, <c>medium</c>, <c>high</c>, <c>xhigh</c>, <c>max</c>).</summary>
    public string? Effort { get; init; }

    /// <inheritdoc />
    public async Task<ChatCompletion> CompleteAsync(ChatRequest request, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(request);
        var body = BuildBody(request);
        using var message = new HttpRequestMessage(HttpMethod.Post, new Uri(BaseUri, "/v1/messages"))
        {
            Content = new StringContent(body.ToJsonString(), Encoding.UTF8, "application/json"),
        };
        message.Headers.Add("x-api-key", _apiKey);
        message.Headers.Add("anthropic-version", "2023-06-01");
        if (WorkspaceId is not null) message.Headers.Add("anthropic-workspace-id", WorkspaceId);
        message.Headers.Accept.Add(new MediaTypeWithQualityHeaderValue("application/json"));
        using var doc = await ChatHttp.PostAsync(_http, message, Provider, cancellationToken).ConfigureAwait(false);
        return Parse(doc.RootElement);
    }

    /// <summary>Builds the JSON request body (exposed for tests).</summary>
    internal JsonObject BuildBody(ChatRequest request)
    {
        var messages = new JsonArray();
        JsonArray? pendingResults = null;
        void FlushResults()
        {
            if (pendingResults is null) return;
            messages.Add(new JsonObject { ["role"] = "user", ["content"] = pendingResults });
            pendingResults = null;
        }

        foreach (var m in request.Messages)
        {
            switch (m.Role)
            {
                case ChatRole.Tool:
                    pendingResults ??= [];
                    pendingResults.Add(new JsonObject
                    {
                        ["type"] = "tool_result",
                        ["tool_use_id"] = m.ToolCallId,
                        ["content"] = m.Text,
                        ["is_error"] = m.IsError,
                    });
                    break;
                case ChatRole.User:
                    FlushResults();
                    messages.Add(new JsonObject { ["role"] = "user", ["content"] = m.Text });
                    break;
                default:
                    FlushResults();
                    JsonNode content;
                    if (m.Provider == Provider && m.ProviderContent is { } raw)
                    {
                        content = JsonNode.Parse(raw.GetRawText())!;
                    }
                    else
                    {
                        var blocks = new JsonArray();
                        if (m.Text.Length > 0) blocks.Add(new JsonObject { ["type"] = "text", ["text"] = m.Text });
                        foreach (var c in m.ToolCalls)
                        {
                            blocks.Add(new JsonObject
                            {
                                ["type"] = "tool_use", ["id"] = c.Id, ["name"] = c.Name,
                                ["input"] = JsonNode.Parse(ChatHttp.ParseObject(c.ArgumentsJson).GetRawText()),
                            });
                        }

                        content = blocks;
                    }

                    messages.Add(new JsonObject { ["role"] = "assistant", ["content"] = content });
                    break;
            }
        }

        FlushResults();
        var body = new JsonObject
        {
            ["model"] = Model,
            ["max_tokens"] = request.MaxTokens,
            ["messages"] = messages,
        };
        if (!string.IsNullOrEmpty(request.SystemPrompt)) body["system"] = request.SystemPrompt;
        if (Effort is not null) body["output_config"] = new JsonObject { ["effort"] = Effort };
        if (request.Tools.Count > 0)
        {
            var tools = new JsonArray();
            foreach (var t in request.Tools)
            {
                tools.Add(new JsonObject
                {
                    ["name"] = t.Name,
                    ["description"] = t.Description,
                    ["input_schema"] = JsonNode.Parse(t.ParametersSchema.GetRawText()),
                });
            }

            body["tools"] = tools;
        }

        return body;
    }

    /// <summary>Parses a Messages API response (exposed for tests).</summary>
    internal ChatCompletion Parse(JsonElement root)
    {
        var text = new StringBuilder();
        var calls = new List<ChatToolCall>();
        var content = root.GetProperty("content");
        foreach (var block in content.EnumerateArray())
        {
            string type = block.GetProperty("type").GetString() ?? "";
            if (type == "text")
            {
                text.Append(block.GetProperty("text").GetString());
            }
            else if (type == "tool_use")
            {
                calls.Add(new ChatToolCall(block.GetProperty("id").GetString()!, block.GetProperty("name").GetString()!,
                    block.GetProperty("input").GetRawText()));
            }
        }

        string stop = root.TryGetProperty("stop_reason", out var sr) && sr.ValueKind == JsonValueKind.String ? sr.GetString()! : "end_turn";
        int? inTok = null, outTok = null;
        if (root.TryGetProperty("usage", out var usage))
        {
            if (usage.TryGetProperty("input_tokens", out var i)) inTok = i.GetInt32();
            if (usage.TryGetProperty("output_tokens", out var o)) outTok = o.GetInt32();
        }

        var message = new ChatMessage
        {
            Role = ChatRole.Assistant,
            Text = text.ToString(),
            ToolCalls = calls,
            ProviderContent = content.Clone(),
            Provider = Provider,
        };
        return new ChatCompletion(message, stop, inTok, outTok);
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (_ownsHttp) _http.Dispose();
    }
}
