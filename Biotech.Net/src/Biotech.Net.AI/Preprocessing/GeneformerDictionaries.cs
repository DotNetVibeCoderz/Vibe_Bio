using System.Buffers.Binary;
using System.Globalization;
using System.Text;
using Biotech.Net.AI.Hub;

namespace Biotech.Net.AI.Preprocessing;

/// <summary>Which published Geneformer gene dictionaries to use.</summary>
public enum GeneformerDictionaryVersion
{
    /// <summary>The dictionaries of Geneformer V1 (Genecorpus-30M, files <c>*_gc30M.pkl</c>).</summary>
    V1Gc30M,

    /// <summary>The dictionaries of Geneformer V2 (Genecorpus-104M, files <c>*_gc104M.pkl</c>).</summary>
    V2Gc104M,
}

/// <summary>
/// Loads Geneformer's token and gene-median dictionaries from the pickles in <c>ctheodoris/Geneformer</c>,
/// from tab-separated text, or straight from the Hub.
/// </summary>
/// <remarks>
/// The pickles are read by a restricted reader that understands only plain containers, strings, numbers and NumPy
/// scalars; any other callable named in the file is refused, so nothing in it is executed.
/// </remarks>
public static class GeneformerDictionaries
{
    /// <summary>The Hub repository holding the dictionaries.</summary>
    public const string RepoId = "ctheodoris/Geneformer";

    /// <summary>Repository paths of the token and median dictionaries for <paramref name="version"/>.</summary>
    public static (string TokenDictionary, string GeneMedians) Paths(GeneformerDictionaryVersion version) => version switch
    {
        GeneformerDictionaryVersion.V1Gc30M => (
            "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl",
            "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl"),
        _ => ("geneformer/token_dictionary_gc104M.pkl", "geneformer/gene_median_dictionary_gc104M.pkl"),
    };

    /// <summary>Downloads (or reads from cache) and loads both dictionaries for <paramref name="version"/>.</summary>
    /// <exception cref="BioModelUnavailableException">The files are neither cached nor downloadable.</exception>
    public static (IReadOnlyDictionary<string, int> Tokens, IReadOnlyDictionary<string, double> Medians) FromHub(
        GeneformerDictionaryVersion version = GeneformerDictionaryVersion.V2Gc104M)
    {
        var (tokens, medians) = Paths(version);
        return (LoadTokenDictionary(BioModelHub.DownloadFile(RepoId, tokens)),
            LoadGeneMedians(BioModelHub.DownloadFile(RepoId, medians)));
    }

    /// <summary>Creates a <see cref="GeneformerEncoder"/> with the published dictionaries and matching settings.</summary>
    public static GeneformerEncoder CreateEncoder(GeneformerDictionaryVersion version = GeneformerDictionaryVersion.V2Gc104M)
    {
        var (tokens, medians) = FromHub(version);
        return new GeneformerEncoder(tokens, medians,
            version == GeneformerDictionaryVersion.V1Gc30M ? GeneformerEncoderOptions.V1 : GeneformerEncoderOptions.V2);
    }

    /// <summary>Loads a gene-to-token-id dictionary from a <c>.pkl</c> file or a two-column TSV/CSV file.</summary>
    public static IReadOnlyDictionary<string, int> LoadTokenDictionary(string path)
    {
        var raw = Load(path);
        var result = new Dictionary<string, int>(raw.Count, StringComparer.Ordinal);
        foreach (var (k, v) in raw) result[k] = checked((int)v);
        return result;
    }

    /// <summary>Loads a gene-to-median dictionary from a <c>.pkl</c> file or a two-column TSV/CSV file.</summary>
    public static IReadOnlyDictionary<string, double> LoadGeneMedians(string path) => Load(path);

    /// <summary>Parses a pickled <c>dict[str, number]</c> (including NumPy scalar values).</summary>
    public static IReadOnlyDictionary<string, double> ReadPickle(Stream stream)
    {
        ArgumentNullException.ThrowIfNull(stream);
        var obj = new RestrictedUnpickler(stream).Load();
        if (obj is not Dictionary<object, object?> dict)
            throw new FormatException($"Expected a pickled dict, found {obj?.GetType().Name ?? "None"}.");
        var result = new Dictionary<string, double>(dict.Count, StringComparer.Ordinal);
        foreach (var (k, v) in dict)
        {
            if (k is not string key) throw new FormatException("Dictionary keys must be strings.");
            result[key] = v switch
            {
                long l => l,
                double d => d,
                bool b => b ? 1 : 0,
                System.Numerics.BigInteger bi => (double)bi,
                _ => throw new FormatException($"Value for '{key}' is {v?.GetType().Name ?? "None"}, not a number."),
            };
        }

        return result;
    }

    /// <summary>Parses two-column text (gene, value) separated by a tab or comma; a non-numeric first row is a header.</summary>
    public static IReadOnlyDictionary<string, double> ReadTable(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var result = new Dictionary<string, double>(StringComparer.Ordinal);
        string? line;
        int lineNo = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (line.Length == 0 || line[0] == '#') continue;
            var parts = line.Split(line.Contains('\t') ? '\t' : ',');
            if (parts.Length < 2) throw new FormatException($"Line {lineNo}: expected two columns.");
            if (!double.TryParse(parts[1].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double value))
            {
                if (lineNo == 1) continue;
                throw new FormatException($"Line {lineNo}: '{parts[1]}' is not a number.");
            }

            result[parts[0].Trim()] = value;
        }

        return result;
    }

    private static IReadOnlyDictionary<string, double> Load(string path)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        if (path.EndsWith(".pkl", StringComparison.OrdinalIgnoreCase) || path.EndsWith(".pickle", StringComparison.OrdinalIgnoreCase))
        {
            using var s = File.OpenRead(path);
            return ReadPickle(s);
        }

        using var r = new StreamReader(path);
        return ReadTable(r);
    }
}

/// <summary>
/// A pickle reader limited to data: dicts, lists, tuples, strings, bytes, numbers, booleans, None and NumPy scalars.
/// Any other global is refused.
/// </summary>
internal sealed class RestrictedUnpickler
{
    private sealed record Global(string Module, string Name);

    private sealed class NumpyDtype
    {
        public required string Code { get; init; }
        public bool BigEndian { get; set; }
    }

    private static readonly object Mark = new();
    private readonly BinaryReader _r;
    private readonly List<object?> _stack = [];
    private readonly Dictionary<long, object?> _memo = [];

    public RestrictedUnpickler(Stream stream) => _r = new BinaryReader(stream, Encoding.UTF8, leaveOpen: true);

    public object? Load()
    {
        while (true)
        {
            byte op = _r.ReadByte();
            switch (op)
            {
                case 0x80: _r.ReadByte(); break; // PROTO
                case 0x95: _r.ReadInt64(); break; // FRAME
                case (byte)'.': return Pop();
                case (byte)'(': _stack.Add(Mark); break;
                case (byte)'}': _stack.Add(new Dictionary<object, object?>()); break;
                case (byte)']': _stack.Add(new List<object?>()); break;
                case (byte)')': _stack.Add(Array.Empty<object?>()); break;
                case (byte)'N': _stack.Add(null); break;
                case 0x88: _stack.Add(true); break;
                case 0x89: _stack.Add(false); break;
                case (byte)'K': _stack.Add((long)_r.ReadByte()); break;
                case (byte)'M': _stack.Add((long)_r.ReadUInt16()); break;
                case (byte)'J': _stack.Add((long)_r.ReadInt32()); break;
                case 0x8a: _stack.Add(ReadLong(_r.ReadByte())); break;
                case 0x8b: _stack.Add(ReadLong(_r.ReadInt32())); break;
                case (byte)'G': _stack.Add(BinaryPrimitives.ReadDoubleBigEndian(_r.ReadBytes(8))); break;
                case 0x8c: _stack.Add(Encoding.UTF8.GetString(_r.ReadBytes(_r.ReadByte()))); break;
                case (byte)'X': _stack.Add(Encoding.UTF8.GetString(_r.ReadBytes(_r.ReadInt32()))); break;
                case 0x8d: _stack.Add(Encoding.UTF8.GetString(_r.ReadBytes(checked((int)_r.ReadInt64())))); break;
                case (byte)'C': _stack.Add(_r.ReadBytes(_r.ReadByte())); break;
                case (byte)'B': _stack.Add(_r.ReadBytes(_r.ReadInt32())); break;
                case 0x8e: _stack.Add(_r.ReadBytes(checked((int)_r.ReadInt64()))); break;
                case (byte)'U': _stack.Add(Encoding.Latin1.GetString(_r.ReadBytes(_r.ReadByte()))); break;
                case (byte)'T': _stack.Add(Encoding.Latin1.GetString(_r.ReadBytes(_r.ReadInt32()))); break;
                case 0x94: _memo[_memo.Count] = Peek(); break; // MEMOIZE
                case (byte)'q': _memo[_r.ReadByte()] = Peek(); break;
                case (byte)'r': _memo[_r.ReadUInt32()] = Peek(); break;
                case (byte)'h': _stack.Add(_memo[_r.ReadByte()]); break;
                case (byte)'j': _stack.Add(_memo[_r.ReadUInt32()]); break;
                case 0x85: _stack.Add(new[] { Pop() }); break;
                case 0x86: { var b = Pop(); var a = Pop(); _stack.Add(new[] { a, b }); break; }
                case 0x87: { var c = Pop(); var b = Pop(); var a = Pop(); _stack.Add(new[] { a, b, c }); break; }
                case (byte)'t': _stack.Add(PopToMark().ToArray()); break;
                case (byte)'l': _stack.Add(PopToMark()); break;
                case (byte)'a': { var v = Pop(); ((List<object?>)Peek()!).Add(v); break; }
                case (byte)'e': { var items = PopToMark(); ((List<object?>)Peek()!).AddRange(items); break; }
                case (byte)'d': { var items = PopToMark(); var d = new Dictionary<object, object?>(); AddPairs(d, items); _stack.Add(d); break; }
                case (byte)'s': { var v = Pop(); var k = Pop(); ((Dictionary<object, object?>)Peek()!)[k!] = v; break; }
                case (byte)'u': { var items = PopToMark(); AddPairs((Dictionary<object, object?>)Peek()!, items); break; }
                case 0x93: { var name = (string)Pop()!; var module = (string)Pop()!; _stack.Add(Resolve(module, name)); break; }
                case (byte)'c': { var module = ReadLine(); var name = ReadLine(); _stack.Add(Resolve(module, name)); break; }
                case (byte)'R': { var args = (object?[])Pop()!; var fn = (Global)Pop()!; _stack.Add(Call(fn, args)); break; }
                case (byte)'b': { var state = Pop(); Build(Peek(), state); break; }
                default:
                    throw new FormatException($"Unsupported pickle opcode 0x{op:X2}; only plain data dictionaries are read.");
            }
        }
    }

    private static Global Resolve(string module, string name)
    {
        string q = module + "." + name;
        return q is "numpy.core.multiarray.scalar" or "numpy._core.multiarray.scalar" or "numpy.dtype" or "_codecs.encode"
            ? new Global(module, name)
            : throw new FormatException($"The pickle refers to '{q}', which this reader will not resolve: only data and NumPy scalars are accepted.");
    }

    private static object? Call(Global fn, object?[] args)
    {
        switch (fn.Name)
        {
            case "dtype":
                return new NumpyDtype { Code = (string)args[0]! };
            case "encode":
                return Encoding.Latin1.GetBytes((string)args[0]!);
            case "scalar":
                var dtype = (NumpyDtype)args[0]!;
                var data = args[1] switch
                {
                    byte[] b => b,
                    string s => Encoding.Latin1.GetBytes(s),
                    _ => throw new FormatException("NumPy scalar without raw data."),
                };
                return DecodeScalar(dtype, data);
            default:
                throw new FormatException($"Cannot call {fn.Module}.{fn.Name}.");
        }
    }

    private static void Build(object? target, object? state)
    {
        if (target is NumpyDtype dt && state is object?[] { Length: > 1 } tuple && tuple[1] is string order)
            dt.BigEndian = order == ">";
    }

    private static object DecodeScalar(NumpyDtype dtype, byte[] data)
    {
        string code = dtype.Code.TrimStart('<', '>', '=', '|');
        bool be = dtype.BigEndian || dtype.Code.StartsWith('>');
        if (be) Array.Reverse(data);
        return code switch
        {
            "f8" => BitConverter.ToDouble(data),
            "f4" => (double)BitConverter.ToSingle(data),
            "f2" => (double)BitConverter.ToHalf(data),
            "i8" => BitConverter.ToInt64(data),
            "i4" => (long)BitConverter.ToInt32(data),
            "i2" => (long)BitConverter.ToInt16(data),
            "i1" => (long)(sbyte)data[0],
            "u8" => (long)BitConverter.ToUInt64(data),
            "u4" => (long)BitConverter.ToUInt32(data),
            "u2" => (long)BitConverter.ToUInt16(data),
            "u1" => (long)data[0],
            "b1" => data[0] != 0,
            _ => throw new FormatException($"Unsupported NumPy scalar type '{dtype.Code}'."),
        };
    }

    private object ReadLong(int length)
    {
        var bytes = _r.ReadBytes(length);
        var value = new System.Numerics.BigInteger(bytes, isUnsigned: false, isBigEndian: false);
        return value >= long.MinValue && value <= long.MaxValue ? (long)value : value;
    }

    private string ReadLine()
    {
        var sb = new StringBuilder();
        char c;
        while ((c = (char)_r.ReadByte()) != '\n') sb.Append(c);
        return sb.ToString();
    }

    private static void AddPairs(Dictionary<object, object?> d, List<object?> items)
    {
        for (int i = 0; i + 1 < items.Count; i += 2) d[items[i]!] = items[i + 1];
    }

    private object? Pop()
    {
        var v = _stack[^1];
        _stack.RemoveAt(_stack.Count - 1);
        return v;
    }

    private object? Peek() => _stack[^1];

    private List<object?> PopToMark()
    {
        int i = _stack.Count - 1;
        while (!ReferenceEquals(_stack[i], Mark)) i--;
        var items = _stack.GetRange(i + 1, _stack.Count - i - 1);
        _stack.RemoveRange(i, _stack.Count - i);
        return items;
    }
}
