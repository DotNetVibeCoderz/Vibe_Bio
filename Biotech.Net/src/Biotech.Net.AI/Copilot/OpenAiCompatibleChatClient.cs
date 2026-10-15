using System.Text;
using System.Text.Json;
using System.Text.Json.Nodes;

namespace Biotech.Net.AI.Copilot;

/// <summary>How the API key is sent.</summary>
public enum ApiKeyStyle
{
    /// <summary><c>Authorization: Bearer &lt;key&gt;</c> (OpenAI, DeepSeek, most compatible servers).</summary>
    Bearer,

    /// <summary><c>api-key: &lt;key&gt;</c> (Azure OpenAI).</summary>
    ApiKeyHeader,
}

/// <summary>
/// <see cref="IChatClient"/> over an OpenAI-style <c>chat/completions</c> endpoint with function calling:
/// Azure OpenAI, OpenAI, DeepSeek, and other compatible servers.
/// </summary>
public sealed class OpenAiCompatibleChatClient : IChatClient, IDisposable
{
    private readonly HttpClient _http;
    private readonly bool _ownsHttp;
    private readonly string _apiKey;

    /// <summary>Creates a client.</summary>
    /// <param name="provider">Provider label (for messages and logs).</param>
    /// <param name="chatCompletionsUri">Full URL of the chat-completions endpoint.</param>
    /// <param name="apiKey">API key.</param>
    /// <param name="model">Model or deployment name sent as <c>model</c>.</param>
    /// <param name="keyStyle">How the key is sent.</param>
    /// <param name="httpClient">Optional HTTP client.</param>
    public OpenAiCompatibleChatClient(string provider, Uri chatCompletionsUri, string apiKey, string model,
        ApiKeyStyle keyStyle = ApiKeyStyle.Bearer, HttpClient? httpClient = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(apiKey);
        Provider = provider;
        Endpoint = chatCompletionsUri ?? throw new ArgumentNullException(nameof(chatCompletionsUri));
        _apiKey = apiKey;
        Model = model;
        KeyStyle = keyStyle;
        _ownsHttp = httpClient is null;
        _http = httpClient ?? new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
    }

    /// <summary>OpenAI (<c>https://api.openai.com/v1/chat/completions</c>).</summary>
    public static OpenAiCompatibleChatClient ForOpenAI(string apiKey, string model, HttpClient? httpClient = null)
        => new("openai", new Uri("https://api.openai.com/v1/chat/completions"), apiKey, model, ApiKeyStyle.Bearer, httpClient);

    /// <summary>
    /// Azure OpenAI through the v1 API (<c>{endpoint}/openai/v1/chat/completions</c>, <c>api-key</c> header);
    /// <paramref name="deployment"/> is the deployment name.
    /// </summary>
    public static OpenAiCompatibleChatClient ForAzureOpenAI(string endpoint, string apiKey, string deployment, HttpClient? httpClient = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(endpoint);
        var baseUri = new Uri(endpoint.TrimEnd('/') + "/");
        return new("azure-openai", new Uri(baseUri, "openai/v1/chat/completions"), apiKey, deployment, ApiKeyStyle.ApiKeyHeader, httpClient);
    }

    /// <summary>DeepSeek (<c>https://api.deepseek.com/chat/completions</c> by default).</summary>
    public static OpenAiCompatibleChatClient ForDeepSeek(string apiKey, string model, string endpoint = "https://api.deepseek.com", HttpClient? httpClient = null)
        => new("deepseek", new Uri(new Uri(endpoint.TrimEnd('/') + "/"), "chat/completions"), apiKey, model, ApiKeyStyle.Bearer, httpClient)
        {
            UseMaxCompletionTokens = false,
        };

    /// <inheritdoc />
    public string Provider { get; }

    /// <inheritdoc />
    public string Model { get; }

    /// <summary>The chat-completions URL.</summary>
    public Uri Endpoint { get; }

    /// <summary>How the key is sent.</summary>
    public ApiKeyStyle KeyStyle { get; }

    /// <summary>
    /// Send <c>max_completion_tokens</c> (OpenAI reasoning models) rather than <c>max_tokens</c>. Default true.
    /// </summary>
    public bool UseMaxCompletionTokens { get; init; } = true;

    /// <inheritdoc />
    public async Task<ChatCompletion> CompleteAsync(ChatRequest request, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(request);
        using var message = new HttpRequestMessage(HttpMethod.Post, Endpoint)
        {
            Content = new StringContent(BuildBody(request).ToJsonString(), Encoding.UTF8, "application/json"),
        };
        if (KeyStyle == ApiKeyStyle.Bearer) message.Headers.Authorization = new("Bearer", _apiKey);
        else message.Headers.Add("api-key", _apiKey);
        using var doc = await ChatHttp.PostAsync(_http, message, Provider, cancellationToken).ConfigureAwait(false);
        return Parse(doc.RootElement);
    }

    /// <summary>Builds the JSON request body (exposed for tests).</summary>
    internal JsonObject BuildBody(ChatRequest request)
    {
        var messages = new JsonArray();
        if (!string.IsNullOrEmpty(request.SystemPrompt)) messages.Add(new JsonObject { ["role"] = "system", ["content"] = request.SystemPrompt });
        foreach (var m in request.Messages)
        {
            switch (m.Role)
            {
                case ChatRole.User:
                    messages.Add(new JsonObject { ["role"] = "user", ["content"] = m.Text });
                    break;
                case ChatRole.Tool:
                    messages.Add(new JsonObject { ["role"] = "tool", ["tool_call_id"] = m.ToolCallId, ["content"] = m.Text });
                    break;
                default:
                    var msg = new JsonObject { ["role"] = "assistant", ["content"] = m.Text.Length > 0 ? m.Text : null };
                    if (m.ToolCalls.Count > 0)
                    {
                        var calls = new JsonArray();
                        foreach (var c in m.ToolCalls)
                        {
                            calls.Add(new JsonObject
                            {
                                ["id"] = c.Id, ["type"] = "function",
                                ["function"] = new JsonObject { ["name"] = c.Name, ["arguments"] = c.ArgumentsJson },
                            });
                        }

                        msg["tool_calls"] = calls;
                    }

                    messages.Add(msg);
                    break;
            }
        }

        var body = new JsonObject { ["model"] = Model, ["messages"] = messages };
        body[UseMaxCompletionTokens ? "max_completion_tokens" : "max_tokens"] = request.MaxTokens;
        if (request.Tools.Count > 0)
        {
            var tools = new JsonArray();
            foreach (var t in request.Tools)
            {
                tools.Add(new JsonObject
                {
                    ["type"] = "function",
                    ["function"] = new JsonObject
                    {
                        ["name"] = t.Name,
                        ["description"] = t.Description,
                        ["parameters"] = JsonNode.Parse(t.ParametersSchema.GetRawText()),
                    },
                });
            }

            body["tools"] = tools;
        }

        return body;
    }

    /// <summary>Parses a chat-completions response (exposed for tests).</summary>
    internal ChatCompletion Parse(JsonElement root)
    {
        if (!root.TryGetProperty("choices", out var choices) || choices.GetArrayLength() == 0)
            throw new ChatClientException($"{Provider}: the response has no choices.");
        var choice = choices[0];
        var msg = choice.GetProperty("message");
        string text = msg.TryGetProperty("content", out var c) && c.ValueKind == JsonValueKind.String ? c.GetString()! : "";
        var calls = new List<ChatToolCall>();
        if (msg.TryGetProperty("tool_calls", out var tc) && tc.ValueKind == JsonValueKind.Array)
        {
            foreach (var call in tc.EnumerateArray())
            {
                var fn = call.GetProperty("function");
                string args = fn.TryGetProperty("arguments", out var a)
                    ? a.ValueKind == JsonValueKind.String ? a.GetString()! : a.GetRawText()
                    : "{}";
                calls.Add(new ChatToolCall(call.GetProperty("id").GetString()!, fn.GetProperty("name").GetString()!, args));
            }
        }

        string finish = choice.TryGetProperty("finish_reason", out var f) && f.ValueKind == JsonValueKind.String ? f.GetString()! : "stop";
        string stop = finish switch
        {
            "stop" => "end_turn",
            "tool_calls" or "function_call" => "tool_use",
            "length" => "max_tokens",
            "content_filter" => "refusal",
            _ => finish,
        };
        int? inTok = null, outTok = null;
        if (root.TryGetProperty("usage", out var usage) && usage.ValueKind == JsonValueKind.Object)
        {
            if (usage.TryGetProperty("prompt_tokens", out var i)) inTok = i.GetInt32();
            if (usage.TryGetProperty("completion_tokens", out var o)) outTok = o.GetInt32();
        }

        return new ChatCompletion(new ChatMessage { Role = ChatRole.Assistant, Text = text, ToolCalls = calls, Provider = Provider },
            stop, inTok, outTok);
    }

    /// <inheritdoc />
    public void Dispose()
    {
        if (_ownsHttp) _http.Dispose();
    }
}
