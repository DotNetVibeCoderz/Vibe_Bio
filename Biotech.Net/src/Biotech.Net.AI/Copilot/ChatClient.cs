using System.Text.Json;

namespace Biotech.Net.AI.Copilot;

/// <summary>Who wrote a chat message.</summary>
public enum ChatRole
{
    /// <summary>The user.</summary>
    User,

    /// <summary>The model.</summary>
    Assistant,

    /// <summary>The result of a tool call, sent back to the model.</summary>
    Tool,
}

/// <summary>A tool call requested by the model.</summary>
/// <param name="Id">Provider-assigned call id, echoed back with the result.</param>
/// <param name="Name">Tool name.</param>
/// <param name="ArgumentsJson">Arguments as a JSON object.</param>
public sealed record ChatToolCall(string Id, string Name, string ArgumentsJson);

/// <summary>A tool the model may call.</summary>
/// <param name="Name">Tool name (letters, digits, underscores).</param>
/// <param name="Description">What the tool does and when to use it.</param>
/// <param name="ParametersSchema">JSON Schema of the arguments object.</param>
public sealed record ChatToolDefinition(string Name, string Description, JsonElement ParametersSchema);

/// <summary>One message of a conversation, in a provider-neutral shape.</summary>
public sealed record ChatMessage
{
    /// <summary>Author.</summary>
    public required ChatRole Role { get; init; }

    /// <summary>Text content (may be empty for an assistant turn that only calls tools).</summary>
    public string Text { get; init; } = "";

    /// <summary>Tool calls made in an assistant turn.</summary>
    public IReadOnlyList<ChatToolCall> ToolCalls { get; init; } = [];

    /// <summary>For <see cref="ChatRole.Tool"/> messages: the call this answers.</summary>
    public string? ToolCallId { get; init; }

    /// <summary>For <see cref="ChatRole.Tool"/> messages: the tool's name.</summary>
    public string? ToolName { get; init; }

    /// <summary>For <see cref="ChatRole.Tool"/> messages: whether the tool failed.</summary>
    public bool IsError { get; init; }

    /// <summary>
    /// The provider's raw content for an assistant turn (for example Anthropic content blocks including thinking blocks),
    /// replayed unchanged on the next request to the same provider.
    /// </summary>
    public JsonElement? ProviderContent { get; init; }

    /// <summary>Which provider produced <see cref="ProviderContent"/>.</summary>
    public string? Provider { get; init; }

    /// <summary>A user message.</summary>
    public static ChatMessage User(string text) => new() { Role = ChatRole.User, Text = text };

    /// <summary>A tool result.</summary>
    public static ChatMessage ToolResult(string callId, string toolName, string content, bool isError = false)
        => new() { Role = ChatRole.Tool, ToolCallId = callId, ToolName = toolName, Text = content, IsError = isError };
}

/// <summary>A request to a chat model.</summary>
public sealed record ChatRequest
{
    /// <summary>System instructions.</summary>
    public string? SystemPrompt { get; init; }

    /// <summary>The conversation so far.</summary>
    public required IReadOnlyList<ChatMessage> Messages { get; init; }

    /// <summary>Tools the model may call.</summary>
    public IReadOnlyList<ChatToolDefinition> Tools { get; init; } = [];

    /// <summary>Maximum output tokens.</summary>
    public int MaxTokens { get; init; } = 8192;
}

/// <summary>The model's reply.</summary>
/// <param name="Message">The assistant message (text and/or tool calls).</param>
/// <param name="StopReason">Provider stop reason, normalised to <c>end_turn</c>, <c>tool_use</c>, <c>max_tokens</c>, <c>refusal</c> or the raw value.</param>
/// <param name="InputTokens">Prompt tokens, when reported.</param>
/// <param name="OutputTokens">Completion tokens, when reported.</param>
public sealed record ChatCompletion(ChatMessage Message, string StopReason, int? InputTokens = null, int? OutputTokens = null);

/// <summary>A chat model behind an HTTP API.</summary>
public interface IChatClient
{
    /// <summary>Provider name, for example <c>anthropic</c>, <c>azure-openai</c>, <c>openai</c>, <c>deepseek</c>.</summary>
    string Provider { get; }

    /// <summary>Model or deployment id.</summary>
    string Model { get; }

    /// <summary>Sends the conversation and returns the model's next turn.</summary>
    /// <exception cref="ChatClientException">The API returned an error or could not be reached.</exception>
    Task<ChatCompletion> CompleteAsync(ChatRequest request, CancellationToken cancellationToken = default);
}

/// <summary>A chat API call failed.</summary>
public sealed class ChatClientException : Exception
{
    /// <summary>Creates the exception.</summary>
    public ChatClientException(string message, int? statusCode = null, Exception? inner = null)
        : base(message, inner) => StatusCode = statusCode;

    /// <summary>HTTP status, when the server answered.</summary>
    public int? StatusCode { get; }
}

internal static class ChatHttp
{
    public static async Task<JsonDocument> PostAsync(HttpClient http, HttpRequestMessage request, string provider, CancellationToken ct)
    {
        HttpResponseMessage response;
        try
        {
            response = await http.SendAsync(request, ct).ConfigureAwait(false);
        }
        catch (HttpRequestException ex)
        {
            throw new ChatClientException($"{provider}: could not reach the API ({ex.Message}).", null, ex);
        }
        catch (TaskCanceledException ex) when (!ct.IsCancellationRequested)
        {
            throw new ChatClientException($"{provider}: the request timed out.", null, ex);
        }

        using (response)
        {
            string body = await response.Content.ReadAsStringAsync(ct).ConfigureAwait(false);
            if (!response.IsSuccessStatusCode)
            {
                string detail = body.Length > 600 ? body[..600] + "…" : body;
                throw new ChatClientException($"{provider}: HTTP {(int)response.StatusCode} {response.ReasonPhrase}: {detail}", (int)response.StatusCode);
            }

            try
            {
                return JsonDocument.Parse(body);
            }
            catch (JsonException ex)
            {
                throw new ChatClientException($"{provider}: the response was not JSON.", (int)response.StatusCode, ex);
            }
        }
    }

    public static JsonElement ParseObject(string json)
    {
        try
        {
            using var doc = JsonDocument.Parse(string.IsNullOrWhiteSpace(json) ? "{}" : json);
            return doc.RootElement.ValueKind == JsonValueKind.Object ? doc.RootElement.Clone() : JsonDocument.Parse("{}").RootElement.Clone();
        }
        catch (JsonException)
        {
            return JsonDocument.Parse("{}").RootElement.Clone();
        }
    }
}
