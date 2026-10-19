using System.Globalization;
using Biotech.Net.AI.Hub;
using Gravicode.HFNet.GraviHub.Io;

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
/// The pickles are read by HF.Net's restricted <see cref="PickleReader"/>, which understands only plain containers,
/// strings, numbers and NumPy scalars; any other callable named in the file is refused, so nothing in it is executed.
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
    /// <remarks>
    /// Read by HF.Net's <see cref="PickleReader"/>, which admits only data and NumPy scalars: any other global named
    /// in the file is refused, so nothing in it is executed.
    /// </remarks>
    /// <exception cref="FormatException">The stream is not a pickled dictionary of numbers, or it names a disallowed global.</exception>
    public static IReadOnlyDictionary<string, double> ReadPickle(Stream stream)
    {
        ArgumentNullException.ThrowIfNull(stream);
        object? obj;
        try
        {
            obj = new PickleReader(stream).Load();
        }
        catch (Exception ex) when (ex is NotSupportedException or InvalidDataException or EndOfStreamException or InvalidCastException)
        {
            throw new FormatException("Not a pickled data dictionary: " + ex.Message, ex);
        }

        if (obj is not Dictionary<object, object?> dict)
            throw new FormatException($"Expected a pickled dict, found {obj?.GetType().Name ?? "None"}.");
        var result = new Dictionary<string, double>(dict.Count, StringComparer.Ordinal);
        foreach (var (k, v) in dict)
        {
            if (k is not string key) throw new FormatException("Dictionary keys must be strings.");
            result[key] = v switch
            {
                bool b => b ? 1 : 0,
                int i => i,
                long l => l,
                ulong u => u,
                double d => d,
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
