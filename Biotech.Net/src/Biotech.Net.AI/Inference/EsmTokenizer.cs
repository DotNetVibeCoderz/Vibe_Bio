namespace Biotech.Net.AI.Inference;

/// <summary>
/// The ESM tokenizer: a fixed vocabulary read from <c>vocab.txt</c>, matched greedily longest-first, with
/// <c>&lt;cls&gt;</c> prepended and <c>&lt;eos&gt;</c> appended (when the vocabulary has one).
/// </summary>
/// <remarks>
/// For ESM-1/ESM-2 every token is one amino acid; for the Nucleotide Transformer (v1) the vocabulary holds 6-mers,
/// which greedy matching splits into non-overlapping 6-mers with single-nucleotide fallback — the same split
/// its Python tokenizer produces. Special tokens such as <c>&lt;mask&gt;</c> are recognised inside the text.
/// HF.Net's tokenizer cannot read these vocabularies (it falls back to unknown tokens), so this is implemented here.
/// </remarks>
public sealed class EsmTokenizer
{
    private readonly Dictionary<string, int> _ids;
    private readonly string[] _tokens;
    private readonly int _maxTokenLength;
    private readonly string[] _specials;

    /// <summary>Creates a tokenizer from the vocabulary, one token per line in id order.</summary>
    public EsmTokenizer(IReadOnlyList<string> vocabulary, bool appendEos = true)
    {
        ArgumentNullException.ThrowIfNull(vocabulary);
        if (vocabulary.Count == 0) throw new ArgumentException("The vocabulary is empty.", nameof(vocabulary));
        _tokens = [.. vocabulary];
        _ids = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < _tokens.Length; i++) _ids.TryAdd(_tokens[i], i);
        _specials = [.. _tokens.Where(t => t.Length > 2 && t[0] == '<' && t[^1] == '>').OrderByDescending(t => t.Length)];
        _maxTokenLength = _tokens.Where(t => !(t.Length > 2 && t[0] == '<')).Select(t => t.Length).DefaultIfEmpty(1).Max();
        ClsId = Id("<cls>");
        PadId = Id("<pad>");
        EosId = appendEos ? Id("<eos>") : -1;
        UnknownId = Id("<unk>");
        MaskId = Id("<mask>");
    }

    /// <summary>Reads <c>vocab.txt</c> (and <c>tokenizer_config.json</c> for <c>"eos_token": null</c>) from a model directory.</summary>
    public static EsmTokenizer FromDirectory(string directory)
    {
        var vocab = File.ReadAllLines(Path.Combine(directory, "vocab.txt")).Select(l => l.TrimEnd('\r')).Where(l => l.Length > 0).ToArray();
        bool appendEos = true;
        string cfg = Path.Combine(directory, "tokenizer_config.json");
        if (File.Exists(cfg))
        {
            using var doc = System.Text.Json.JsonDocument.Parse(File.ReadAllText(cfg));
            if (doc.RootElement.TryGetProperty("eos_token", out var eos) && eos.ValueKind == System.Text.Json.JsonValueKind.Null)
                appendEos = false;
        }

        return new EsmTokenizer(vocab, appendEos);
    }

    /// <summary>Number of tokens in the vocabulary.</summary>
    public int VocabularySize => _tokens.Length;

    /// <summary>Id of <c>&lt;cls&gt;</c>, or -1.</summary>
    public int ClsId { get; }

    /// <summary>Id of <c>&lt;pad&gt;</c>, or -1.</summary>
    public int PadId { get; }

    /// <summary>Id of <c>&lt;eos&gt;</c>, or -1 when the model does not append one.</summary>
    public int EosId { get; }

    /// <summary>Id of <c>&lt;unk&gt;</c>, or -1.</summary>
    public int UnknownId { get; }

    /// <summary>Id of <c>&lt;mask&gt;</c>, or -1.</summary>
    public int MaskId { get; }

    /// <summary>The mask token's text.</summary>
    public string MaskToken => "<mask>";

    /// <summary>Id of a token, or -1 when absent.</summary>
    public int Id(string token) => _ids.TryGetValue(token, out int id) ? id : -1;

    /// <summary>Text of a token id.</summary>
    public string Token(int id) => _tokens[id];

    /// <summary>Whether <paramref name="id"/> is a special token (cls, pad, eos, unk, mask, …).</summary>
    public bool IsSpecial(int id) => id >= 0 && id < _tokens.Length && _tokens[id].Length > 2 && _tokens[id][0] == '<';

    /// <summary>Tokenises <paramref name="text"/>, adding <c>&lt;cls&gt;</c> and <c>&lt;eos&gt;</c> when requested.</summary>
    public int[] Encode(string text, bool addSpecialTokens = true)
    {
        ArgumentNullException.ThrowIfNull(text);
        var ids = new List<int>(text.Length + 2);
        if (addSpecialTokens && ClsId >= 0) ids.Add(ClsId);
        int i = 0;
        while (i < text.Length)
        {
            char c = text[i];
            if (char.IsWhiteSpace(c))
            {
                i++;
                continue;
            }

            if (c == '<')
            {
                string? special = _specials.FirstOrDefault(s => string.CompareOrdinal(text, i, s, 0, s.Length) == 0);
                if (special is not null)
                {
                    ids.Add(_ids[special]);
                    i += special.Length;
                    continue;
                }
            }

            int matched = 0;
            for (int len = Math.Min(_maxTokenLength, text.Length - i); len >= 1; len--)
            {
                if (_ids.TryGetValue(text.Substring(i, len), out int id))
                {
                    ids.Add(id);
                    matched = len;
                    break;
                }
            }

            if (matched == 0)
            {
                ids.Add(UnknownId);
                matched = 1;
            }

            i += matched;
        }

        if (addSpecialTokens && EosId >= 0) ids.Add(EosId);
        return [.. ids];
    }

    /// <summary>Joins token texts back into a string, skipping special tokens unless asked.</summary>
    public string Decode(IEnumerable<int> ids, bool skipSpecialTokens = true)
        => string.Concat(ids.Where(id => !skipSpecialTokens || !IsSpecial(id)).Select(Token));
}
