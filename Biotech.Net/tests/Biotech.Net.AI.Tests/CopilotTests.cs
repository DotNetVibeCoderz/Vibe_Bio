using System.Text.Json;
using Biotech.Net.AI.Copilot;

namespace Biotech.Net.AI.Tests;

public class BioToolboxTests
{
    private readonly BioToolbox _tools = new();

    private JsonElement Run(string name, object args)
    {
        var r = _tools.Invoke(name, JsonSerializer.Serialize(args));
        Assert.False(r.IsError, r.Content);
        return JsonDocument.Parse(r.Content).RootElement;
    }

    [Fact]
    public void BuiltInTools_AreRegisteredWithObjectSchemas()
    {
        var names = _tools.Definitions.Select(d => d.Name).ToArray();
        Assert.Equal(["reverse_complement", "transcribe", "translate", "gc_content", "protein_properties", "find_orfs"], names);
        Assert.All(_tools.Definitions, d => Assert.Equal("object", d.ParametersSchema.GetProperty("type").GetString()));
    }

    [Fact]
    public void ReverseComplement_Transcribe_Translate()
    {
        Assert.Equal("GTTAACGCGCAT", Run("reverse_complement", new { sequence = "ATGCGCGTTAAC" }).GetProperty("reverse_complement").GetString());
        Assert.Equal("AUGCGC", Run("transcribe", new { sequence = "ATGCGC" }).GetProperty("rna").GetString());
        Assert.Equal("MA*", Run("translate", new { sequence = "ATGGCCTAA" }).GetProperty("protein").GetString());
        Assert.Equal("MA", Run("translate", new { sequence = "ATGGCCTAA", to_stop = true }).GetProperty("protein").GetString());
        Assert.Equal("MA", Run("translate", new { sequence = "AATGGCCTA", frame = "1" }).GetProperty("protein").GetString());
    }

    [Fact]
    public void GcContent_AndProteinProperties()
    {
        var gc = Run("gc_content", new { sequence = "GGCCAATT" });
        Assert.Equal(0.5, gc.GetProperty("gc_fraction").GetDouble());
        Assert.Equal(8, gc.GetProperty("length").GetInt32());

        var p = Run("protein_properties", new { sequence = "MKTAYIAKQRQISFVKSHFSRQ" });
        Assert.Equal(22, p.GetProperty("length").GetInt32());
        Assert.InRange(p.GetProperty("molecular_weight_da").GetDouble(), 2600, 2700);
        Assert.InRange(p.GetProperty("isoelectric_point").GetDouble(), 10, 12);
        Assert.True(p.GetProperty("composition_percent").TryGetProperty("K", out _));
    }

    [Fact]
    public void FindOrfs_FindsBothStrands()
    {
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCT", 10)) + "TAA"; // M + 10 A + stop, 36 nt
        string dna = "CC" + orf + "GG" + Biotech.Net.Core.Sequence.Dna(orf).ReverseComplement() + "C";
        var r = Run("find_orfs", new { sequence = dna, min_protein_length = 10 });
        var orfs = r.GetProperty("orfs").EnumerateArray().ToArray();
        var plus = orfs.Single(o => o.GetProperty("Strand").GetString() == "+" && o.GetProperty("Start").GetInt32() == 3);
        Assert.Equal(38, plus.GetProperty("End").GetInt32());
        Assert.Equal("MAAAAAAAAAA", plus.GetProperty("Protein").GetString());
        var minus = orfs.Single(o => o.GetProperty("Strand").GetString() == "-" && o.GetProperty("Protein").GetString() == "MAAAAAAAAAA");
        Assert.Equal(41, minus.GetProperty("Start").GetInt32());
        Assert.Equal(76, minus.GetProperty("End").GetInt32());
    }

    [Fact]
    public void Failures_AreReturnedAsErrors()
    {
        Assert.True(_tools.Invoke("nope", "{}").IsError);
        Assert.True(_tools.Invoke("gc_content", "{}").IsError);
        Assert.True(_tools.Invoke("gc_content", """{"sequence":"HELLO WORLD"}""").IsError);
        Assert.True(_tools.Invoke("translate", "not json").IsError);
    }

    [Fact]
    public void CustomTools_CanBeAdded()
    {
        var box = new BioToolbox(includeBuiltIns: false)
            .Add("echo", "Echo", BioToolbox.Schema(("x", "string", "value", true)), a => a.GetProperty("x").GetString()!);
        Assert.Equal("hi", box.Invoke("echo", """{"x":"hi"}""").Content);
    }
}

public class BioCopilotTests
{
    private sealed class ScriptedChatClient(params ChatCompletion[] replies) : IChatClient
    {
        private int _next;
        public List<ChatRequest> Requests { get; } = [];
        public string Provider => "fake";
        public string Model => "fake-model";

        public Task<ChatCompletion> CompleteAsync(ChatRequest request, CancellationToken cancellationToken = default)
        {
            Requests.Add(request with { Messages = [.. request.Messages] });
            return Task.FromResult(replies[Math.Min(_next++, replies.Length - 1)]);
        }
    }

    private static ChatCompletion Calls(params ChatToolCall[] calls)
        => new(new ChatMessage { Role = ChatRole.Assistant, ToolCalls = calls }, "tool_use");

    private static ChatCompletion Says(string text) => new(new ChatMessage { Role = ChatRole.Assistant, Text = text }, "end_turn");

    [Fact]
    public async Task ToolLoop_RunsToolsAndReturnsFinalAnswer()
    {
        var client = new ScriptedChatClient(
            Calls(new ChatToolCall("c1", "gc_content", """{"sequence":"GGCC"}"""),
                  new ChatToolCall("c2", "reverse_complement", """{"sequence":"AAC"}""")),
            Says("GC is 100%; reverse complement is GTT."));
        var copilot = new BioCopilot(client);
        var answer = await copilot.AskAsync("GC of GGCC?", TestContext.Current.CancellationToken);

        Assert.Equal("GC is 100%; reverse complement is GTT.", answer.Text);
        Assert.Equal(2, answer.Rounds);
        Assert.Equal("end_turn", answer.StopReason);
        Assert.Equal(["gc_content", "reverse_complement"], answer.ToolCalls.Select(t => t.Name));
        Assert.Contains("\"gc_fraction\":1", answer.ToolCalls[0].Result);

        var second = client.Requests[1];
        Assert.Equal(4, second.Messages.Count); // user, assistant(tool calls), 2 tool results
        Assert.Equal(ChatRole.Tool, second.Messages[2].Role);
        Assert.Equal("c1", second.Messages[2].ToolCallId);
        Assert.Equal(6, second.Tools.Count);
        Assert.Contains("Biotech.Net", second.SystemPrompt);
    }

    [Fact]
    public async Task ToolLoop_StopsAtMaxRounds_AndKeepsHistory()
    {
        var client = new ScriptedChatClient(Calls(new ChatToolCall("c", "gc_content", """{"sequence":"A"}""")));
        var copilot = new BioCopilot(client, options: new BioCopilotOptions { MaxRounds = 3 });
        var answer = await copilot.AskAsync("loop", TestContext.Current.CancellationToken);
        Assert.Equal("max_rounds", answer.StopReason);
        Assert.Equal(3, answer.Rounds);
        Assert.Equal(1 + 3 * 2, copilot.History.Count);
        copilot.Reset();
        Assert.Empty(copilot.History);
    }

    [Fact]
    public async Task ToolErrors_AreSentBackFlagged()
    {
        var client = new ScriptedChatClient(Calls(new ChatToolCall("c", "translate", "{}")), Says("Sorry."));
        var answer = await new BioCopilot(client).AskAsync("translate", TestContext.Current.CancellationToken);
        Assert.True(answer.ToolCalls[0].IsError);
        Assert.True(client.Requests[1].Messages[2].IsError);
    }

    [Theory]
    [InlineData("sk-ant-api03-abc", ChatProvider.Anthropic)]
    [InlineData("sk-proj-abc", ChatProvider.OpenAI)]
    [InlineData("0123456789abcdef0123456789abcdef", ChatProvider.AzureOpenAI)]
    [InlineData("", ChatProvider.Unknown)]
    public void DetectProvider_UsesPrefixOnly(string key, ChatProvider expected)
        => Assert.Equal(expected, ChatClientFactory.DetectProvider(key));
}

public class ChatWireFormatTests
{
    private static readonly ChatToolDefinition Tool = new("gc_content", "GC", BioToolbox.Schema(("sequence", "string", "seq", true)));

    private static ChatRequest Conversation(ChatMessage assistant) => new()
    {
        SystemPrompt = "sys",
        Messages =
        [
            ChatMessage.User("hi"), assistant,
            ChatMessage.ToolResult("t1", "gc_content", "{\"gc\":1}"),
            ChatMessage.ToolResult("t2", "gc_content", "bad", isError: true),
        ],
        Tools = [Tool],
        MaxTokens = 100,
    };

    [Fact]
    public void Anthropic_GroupsToolResultsAndReplaysContentBlocks()
    {
        using var client = new AnthropicChatClient("test-key");
        using var doc = JsonDocument.Parse("""
            {"content":[{"type":"thinking","thinking":"","signature":"sig"},{"type":"text","text":"Let me check."},
            {"type":"tool_use","id":"t1","name":"gc_content","input":{"sequence":"GC"}}],
            "stop_reason":"tool_use","usage":{"input_tokens":10,"output_tokens":5}}
            """);
        var completion = client.Parse(doc.RootElement);
        Assert.Equal("tool_use", completion.StopReason);
        Assert.Equal("Let me check.", completion.Message.Text);
        Assert.Equal("t1", completion.Message.ToolCalls.Single().Id);
        Assert.Equal(10, completion.InputTokens);

        var body = client.BuildBody(Conversation(completion.Message));
        Assert.Equal("claude-opus-5-5", (string?)body["model"]);
        Assert.Equal("sys", (string?)body["system"]);
        var messages = body["messages"]!.AsArray();
        Assert.Equal(3, messages.Count);
        Assert.Equal("thinking", (string?)messages[1]!["content"]![0]!["type"]); // replayed unchanged
        var results = messages[2]!["content"]!.AsArray();
        Assert.Equal(2, results.Count);
        Assert.Equal("tool_result", (string?)results[0]!["type"]);
        Assert.True((bool)results[1]!["is_error"]!);
        Assert.Equal("object", (string?)body["tools"]![0]!["input_schema"]!["type"]);
    }

    [Fact]
    public void OpenAiCompatible_BuildsFunctionCallingBodyAndParses()
    {
        using var client = OpenAiCompatibleChatClient.ForAzureOpenAI("https://example.openai.azure.com/", "k", "my-deployment");
        Assert.Equal("https://example.openai.azure.com/openai/v1/chat/completions", client.Endpoint.ToString());
        using var doc = JsonDocument.Parse("""
            {"choices":[{"finish_reason":"tool_calls","message":{"role":"assistant","content":null,
            "tool_calls":[{"id":"t1","type":"function","function":{"name":"gc_content","arguments":"{\"sequence\":\"GC\"}"}}]}}],
            "usage":{"prompt_tokens":7,"completion_tokens":3}}
            """);
        var completion = client.Parse(doc.RootElement);
        Assert.Equal("tool_use", completion.StopReason);
        Assert.Equal("""{"sequence":"GC"}""", completion.Message.ToolCalls[0].ArgumentsJson);

        var body = client.BuildBody(Conversation(completion.Message));
        var messages = body["messages"]!.AsArray();
        Assert.Equal(["system", "user", "assistant", "tool", "tool"], messages.Select(m => (string?)m!["role"]));
        Assert.Equal("t2", (string?)messages[4]!["tool_call_id"]);
        Assert.Equal(100, (int)body["max_completion_tokens"]!);
        Assert.Equal("function", (string?)body["tools"]![0]!["type"]);

        using var ds = OpenAiCompatibleChatClient.ForDeepSeek("k", "deepseek-chat");
        Assert.Equal(100, (int)ds.BuildBody(Conversation(completion.Message))["max_tokens"]!);
        Assert.Equal("https://api.deepseek.com/chat/completions", ds.Endpoint.ToString());
    }
}
