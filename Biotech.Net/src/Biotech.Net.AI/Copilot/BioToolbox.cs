using System.Globalization;
using System.Text;
using System.Text.Json;
using Biotech.Net.AI.Preprocessing;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.AI.Copilot;

/// <summary>A tool the copilot can call: a schema for the model and a handler that runs locally.</summary>
/// <param name="Definition">Name, description and argument schema shown to the model.</param>
/// <param name="Handler">Runs the tool on the parsed arguments and returns a JSON (or plain-text) result.</param>
public sealed record BioTool(ChatToolDefinition Definition, Func<JsonElement, string> Handler)
{
    /// <summary>Tool name.</summary>
    public string Name => Definition.Name;
}

/// <summary>The result of running a tool.</summary>
/// <param name="Content">What is sent back to the model.</param>
/// <param name="IsError">Whether the tool failed (bad arguments, invalid sequence, unknown tool).</param>
public sealed record BioToolResult(string Content, bool IsError);

/// <summary>
/// The Biotech.Net tools exposed to an LLM: reverse complement, transcription, translation, GC content and
/// composition, ProtParam-style protein properties and ORF finding. All computation runs locally on Biotech.Net Core
/// and its kernels; the model only sees the results.
/// </summary>
public sealed class BioToolbox
{
    private const int MaxSequenceLength = 1_000_000;
    private readonly Dictionary<string, BioTool> _tools = new(StringComparer.Ordinal);
    private static readonly JsonSerializerOptions Json = new() { WriteIndented = false };

    /// <summary>Creates a toolbox, with the built-in tools unless <paramref name="includeBuiltIns"/> is false.</summary>
    public BioToolbox(bool includeBuiltIns = true)
    {
        if (!includeBuiltIns) return;
        Add("reverse_complement", "Reverse complement of a DNA or RNA sequence (IUPAC codes kept).",
            Schema(("sequence", "string", "Nucleotide sequence (FASTA accepted)", true)), ReverseComplement);
        Add("transcribe", "Transcribe DNA to RNA (T to U).",
            Schema(("sequence", "string", "DNA sequence", true)), args => Result(new
            {
                rna = Sequence.Dna(SequencePreprocessor.CleanDna(Str(args, "sequence"), ambiguityToN: false)).Transcribe().ToString(),
            }));
        Add("translate", "Translate a DNA/RNA sequence to protein in a reading frame with an NCBI genetic code. Stops are '*'.",
            Schema(("sequence", "string", "Nucleotide sequence", true),
                ("frame", "integer", "Reading frame 0, 1 or 2 (default 0)", false),
                ("table", "integer", "NCBI translation table id (default 1, standard)", false),
                ("to_stop", "boolean", "Stop at the first stop codon (default false)", false)), Translate);
        Add("gc_content", "GC fraction, GC skew, length and base counts of a nucleotide sequence.",
            Schema(("sequence", "string", "Nucleotide sequence", true)), GcContent);
        Add("protein_properties",
            "ProtParam-style properties of a protein: length, molecular weight (Da), isoelectric point, net charge at pH 7, GRAVY (Kyte-Doolittle) and amino-acid composition.",
            Schema(("sequence", "string", "Protein sequence in one-letter code", true)), ProteinProperties);
        Add("find_orfs", "Find open reading frames (ATG to stop) on both strands; returns 1-based coordinates, strand, frame and protein.",
            Schema(("sequence", "string", "DNA sequence", true),
                ("min_protein_length", "integer", "Minimum protein length in amino acids (default 30)", false),
                ("table", "integer", "NCBI translation table id (default 1)", false)), FindOrfs);
    }

    /// <summary>The registered tools.</summary>
    public IReadOnlyCollection<BioTool> Tools => _tools.Values;

    /// <summary>Tool definitions to send to the model.</summary>
    public IReadOnlyList<ChatToolDefinition> Definitions => [.. _tools.Values.Select(t => t.Definition)];

    /// <summary>Registers (or replaces) a tool.</summary>
    public BioToolbox Add(BioTool tool)
    {
        ArgumentNullException.ThrowIfNull(tool);
        _tools[tool.Name] = tool;
        return this;
    }

    /// <summary>Registers a tool from its parts.</summary>
    public BioToolbox Add(string name, string description, JsonElement schema, Func<JsonElement, string> handler)
        => Add(new BioTool(new ChatToolDefinition(name, description, schema), handler));

    /// <summary>Runs a tool by name. Failures are returned as error results rather than thrown, so the model can recover.</summary>
    public BioToolResult Invoke(string name, string argumentsJson)
    {
        if (!_tools.TryGetValue(name, out var tool))
            return new BioToolResult($"Unknown tool '{name}'. Available: {string.Join(", ", _tools.Keys)}.", true);
        try
        {
            using var doc = JsonDocument.Parse(string.IsNullOrWhiteSpace(argumentsJson) ? "{}" : argumentsJson);
            return new BioToolResult(tool.Handler(doc.RootElement), false);
        }
        catch (Exception ex) when (ex is JsonException or FormatException or ArgumentException or KeyNotFoundException or InvalidOperationException)
        {
            return new BioToolResult($"{name} failed: {ex.Message}", true);
        }
    }

    /// <summary>Builds a JSON Schema object for simple flat arguments.</summary>
    public static JsonElement Schema(params (string Name, string Type, string Description, bool Required)[] properties)
    {
        var props = properties.ToDictionary(p => p.Name, p => (object)new { type = p.Type, description = p.Description });
        var schema = new Dictionary<string, object>
        {
            ["type"] = "object",
            ["properties"] = props,
            ["required"] = properties.Where(p => p.Required).Select(p => p.Name).ToArray(),
        };
        return JsonSerializer.SerializeToElement(schema, Json);
    }

    private static string ReverseComplement(JsonElement args)
    {
        string s = SequencePreprocessor.CleanDna(Str(args, "sequence"), ambiguityToN: false);
        return Result(new { reverse_complement = Sequence.Dna(s).ReverseComplement().ToString(), length = s.Length });
    }

    private static string Translate(JsonElement args)
    {
        string s = SequencePreprocessor.CleanDna(Str(args, "sequence"));
        int frame = Int(args, "frame", 0);
        var code = GeneticCode.FromId(Int(args, "table", 1));
        bool toStop = args.TryGetProperty("to_stop", out var ts) && ts.ValueKind == JsonValueKind.True;
        var protein = Sequence.Dna(s).Translate(code, frame, toStop).ToString();
        return Result(new { protein, length = protein.Length, frame, table = code.Id, table_name = code.Name });
    }

    private static string GcContent(JsonElement args)
    {
        string s = SequencePreprocessor.CleanDna(Str(args, "sequence"));
        var c = Sequence.Dna(s).NucleotideCounts();
        return Result(new
        {
            length = s.Length,
            gc_fraction = Math.Round(c.GcFraction, 6),
            gc_percent = Math.Round(c.GcFraction * 100, 3),
            gc_skew = Math.Round(c.GcSkew, 6),
            counts = new { A = c.A, C = c.C, G = c.G, T = c.T, N = c.N },
        });
    }

    private static string ProteinProperties(JsonElement args)
    {
        string s = SequencePreprocessor.CleanProtein(Str(args, "sequence"));
        if (s.Length == 0) throw new ArgumentException("Empty protein sequence.");
        var bytes = Encoding.ASCII.GetBytes(s);
        var gravy = new double[1];
        BioKernels.Hydropathy(bytes, s.Length, gravy);
        var composition = s.GroupBy(ch => ch).OrderBy(g => g.Key)
            .ToDictionary(g => g.Key.ToString(), g => Math.Round(100.0 * g.Count() / s.Length, 2));
        return Result(new
        {
            length = s.Length,
            molecular_weight_da = Math.Round(BioKernels.ProteinMass(bytes), 2),
            isoelectric_point = Math.Round(BioKernels.IsoelectricPoint(bytes), 2),
            net_charge_ph7 = Math.Round(BioKernels.NetCharge(bytes, 7.0), 2),
            gravy = Math.Round(gravy[0], 4),
            composition_percent = composition,
        });
    }

    private static string FindOrfs(JsonElement args)
    {
        string s = SequencePreprocessor.CleanDna(Str(args, "sequence"));
        int minAa = Int(args, "min_protein_length", 30);
        var code = GeneticCode.FromId(Int(args, "table", 1));
        var orfs = Orfs(s, minAa, code).OrderByDescending(o => o.ProteinLength).Take(50).ToArray();
        return Result(new { count = orfs.Length, min_protein_length = minAa, orfs });
    }

    /// <summary>One open reading frame.</summary>
    /// <param name="Start">1-based first base (on the forward strand coordinates).</param>
    /// <param name="End">1-based last base, stop codon included.</param>
    /// <param name="Strand"><c>+</c> or <c>-</c>.</param>
    /// <param name="Frame">Frame 0–2 on its strand.</param>
    /// <param name="ProteinLength">Amino acids, stop excluded.</param>
    /// <param name="Protein">The translated protein.</param>
    public sealed record Orf(int Start, int End, char Strand, int Frame, int ProteinLength, string Protein);

    /// <summary>Finds ATG-initiated ORFs ending in a stop codon on both strands.</summary>
    public static IEnumerable<Orf> Orfs(string dna, int minProteinLength, GeneticCode? code = null)
    {
        code ??= GeneticCode.Standard;
        string fwd = SequencePreprocessor.CleanDna(dna);
        string rev = Sequence.Dna(fwd).ReverseComplement().ToString();
        foreach (var (strandSeq, strand) in new[] { (fwd, '+'), (rev, '-') })
        {
            for (int frame = 0; frame < 3; frame++)
            {
                string protein = Sequence.Dna(strandSeq).Translate(code, frame).ToString();
                int aa = 0;
                while (aa < protein.Length)
                {
                    int m = protein.IndexOf('M', aa);
                    if (m < 0) break;
                    int stop = protein.IndexOf('*', m);
                    if (stop < 0) break;
                    int len = stop - m;
                    if (len >= minProteinLength)
                    {
                        int ntStart = frame + 3 * m;           // 0-based on this strand
                        int ntEnd = frame + 3 * stop + 3;      // exclusive
                        (int s1, int e1) = strand == '+'
                            ? (ntStart + 1, ntEnd)
                            : (fwd.Length - ntEnd + 1, fwd.Length - ntStart);
                        yield return new Orf(s1, e1, strand, frame, len, protein.Substring(m, len));
                    }

                    aa = stop + 1;
                }
            }
        }
    }

    private static string Str(JsonElement args, string name)
    {
        if (!args.TryGetProperty(name, out var v) || v.ValueKind != JsonValueKind.String)
            throw new ArgumentException($"Missing string argument '{name}'.");
        string s = v.GetString()!;
        if (s.Length > MaxSequenceLength) throw new ArgumentException($"'{name}' is longer than {MaxSequenceLength:N0} characters.");
        return s;
    }

    private static int Int(JsonElement args, string name, int fallback)
    {
        if (!args.TryGetProperty(name, out var v)) return fallback;
        return v.ValueKind switch
        {
            JsonValueKind.Number => v.GetInt32(),
            JsonValueKind.String when int.TryParse(v.GetString(), NumberStyles.Integer, CultureInfo.InvariantCulture, out int i) => i,
            JsonValueKind.Null => fallback,
            _ => throw new ArgumentException($"'{name}' must be an integer."),
        };
    }

    private static string Result(object value) => JsonSerializer.Serialize(value, Json);
}
