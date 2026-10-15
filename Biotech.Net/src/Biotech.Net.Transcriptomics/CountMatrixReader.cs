using System.Globalization;
using System.IO.Compression;
using Biotech.Net.Core;

namespace Biotech.Net.Transcriptomics;

/// <summary>
/// Imports gene-level RNA-Seq count matrices produced by featureCounts (Subread), HTSeq-count, or any
/// tab/comma-separated genes × samples table. Results are <see cref="Dataset"/>s with genes in rows and samples
/// in columns; gene lengths (when the format carries them) are stored in the <c>length</c> row annotation and
/// can be read back with <see cref="GeneLengths(Dataset)"/>.
/// </summary>
public static class CountMatrixReader
{
    /// <summary>Row-annotation key holding the gene length (bp) as an invariant-culture number.</summary>
    public const string LengthAnnotation = "length";

    /// <summary>The HTSeq-count summary rows (<c>__no_feature</c>, <c>__ambiguous</c>, ...) that are not genes.</summary>
    public static readonly IReadOnlySet<string> HtSeqSpecialCounters = new HashSet<string>(StringComparer.Ordinal)
    {
        "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique",
    };

    /// <summary>Opens a file for reading, transparently decompressing <c>.gz</c>.</summary>
    public static TextReader OpenText(string path)
    {
        Stream s = File.OpenRead(path);
        if (path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)) s = new GZipStream(s, CompressionMode.Decompress);
        return new StreamReader(s);
    }

    /// <summary>
    /// Reads featureCounts output: <c>#</c> comment lines, a header <c>Geneid Chr Start End Strand Length s1 s2 ...</c>,
    /// then one line per gene. Sample names are cleaned to the file name without directory and <c>.bam</c>/<c>.sam</c>
    /// extension when <paramref name="cleanSampleNames"/> is true. Gene lengths go to the <c>length</c> annotation.
    /// </summary>
    public static Dataset ReadFeatureCounts(TextReader reader, string name = "counts", bool cleanSampleNames = true)
    {
        ArgumentNullException.ThrowIfNull(reader);
        string? line;
        string[]? header = null;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            header = line.Split('\t');
            break;
        }
        if (header is null) throw new InvalidDataException("featureCounts file has no header line.");
        int lengthCol = Array.FindIndex(header, h => h.Equals("Length", StringComparison.OrdinalIgnoreCase));
        if (lengthCol < 0 || !header[0].Equals("Geneid", StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("Not a featureCounts table: expected 'Geneid ... Length' header.");
        int first = lengthCol + 1;
        var samples = new string[header.Length - first];
        for (int j = 0; j < samples.Length; j++) samples[j] = cleanSampleNames ? CleanSampleName(header[first + j]) : header[first + j];

        var genes = new List<string>();
        var lengths = new List<string>();
        var values = new List<double>();
        int lineNo = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (line.Length == 0 || line[0] == '#') continue;
            var parts = line.Split('\t');
            if (parts.Length != header.Length)
                throw new InvalidDataException($"featureCounts line {lineNo}: expected {header.Length} fields, found {parts.Length}.");
            genes.Add(parts[0]);
            lengths.Add(parts[lengthCol]);
            for (int j = first; j < parts.Length; j++) values.Add(ParseCount(parts[j], lineNo));
        }
        var d = new Dataset(name, genes, samples, values.ToArray());
        d.RowAnnotations[LengthAnnotation] = lengths.ToArray();
        return d;
    }

    /// <summary>
    /// Combines per-sample HTSeq-count outputs (<c>gene\tcount</c>) into one matrix. Summary counters such as
    /// <c>__no_feature</c> are dropped. Genes missing from a sample get 0; the gene order is that of the first file
    /// followed by any new genes in order of appearance.
    /// </summary>
    public static Dataset ReadHtSeq(IEnumerable<(string Sample, TextReader Reader)> samples, string name = "counts")
    {
        ArgumentNullException.ThrowIfNull(samples);
        var geneIndex = new Dictionary<string, int>(StringComparer.Ordinal);
        var genes = new List<string>();
        var columns = new List<Dictionary<int, double>>();
        var sampleNames = new List<string>();
        foreach (var (sample, reader) in samples)
        {
            var col = new Dictionary<int, double>();
            string? line;
            int lineNo = 0;
            while ((line = reader.ReadLine()) is not null)
            {
                lineNo++;
                if (line.Length == 0) continue;
                int tab = line.IndexOf('\t');
                if (tab < 0) throw new InvalidDataException($"HTSeq {sample} line {lineNo}: expected 'gene<TAB>count'.");
                string gene = line[..tab];
                if (HtSeqSpecialCounters.Contains(gene)) continue;
                if (!geneIndex.TryGetValue(gene, out int gi)) { gi = genes.Count; geneIndex[gene] = gi; genes.Add(gene); }
                col[gi] = ParseCount(line[(tab + 1)..], lineNo);
            }
            columns.Add(col);
            sampleNames.Add(sample);
        }
        var values = new double[genes.Count * columns.Count];
        for (int j = 0; j < columns.Count; j++)
            foreach (var kv in columns[j]) values[kv.Key * columns.Count + j] = kv.Value;
        return new Dataset(name, genes, sampleNames, values);
    }

    /// <summary>Reads HTSeq-count files; the sample name is the file name without extension.</summary>
    public static Dataset ReadHtSeqFiles(IEnumerable<string> paths, string name = "counts")
    {
        var readers = paths.Select(p => (Sample: CleanSampleName(p), Reader: OpenText(p))).ToList();
        try { return ReadHtSeq(readers, name); }
        finally { foreach (var r in readers) r.Reader.Dispose(); }
    }

    /// <summary>
    /// Reads a plain genes × samples table (header row of sample names, first column gene ids). Extra annotation
    /// columns can be skipped with <paramref name="annotationColumns"/> (columns 2..n that are not counts); a column
    /// named <c>length</c> among them is kept as the gene-length annotation.
    /// </summary>
    public static Dataset ReadTable(TextReader reader, string name = "counts", char separator = '\t', int annotationColumns = 0)
    {
        ArgumentNullException.ThrowIfNull(reader);
        if (annotationColumns == 0) return Dataset.ReadDelimited(reader, name, separator);
        string header = reader.ReadLine() ?? throw new InvalidDataException("Empty count table.");
        var h = header.Split(separator);
        int first = 1 + annotationColumns;
        int lengthCol = -1;
        for (int j = 1; j < first; j++) if (h[j].Equals("length", StringComparison.OrdinalIgnoreCase)) lengthCol = j;
        var genes = new List<string>();
        var lengths = new List<string>();
        var values = new List<double>();
        string? line;
        int lineNo = 1;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (line.Length == 0) continue;
            var parts = line.Split(separator);
            if (parts.Length != h.Length) throw new InvalidDataException($"Line {lineNo}: expected {h.Length} fields, found {parts.Length}.");
            genes.Add(parts[0]);
            if (lengthCol >= 0) lengths.Add(parts[lengthCol]);
            for (int j = first; j < parts.Length; j++) values.Add(ParseCount(parts[j], lineNo));
        }
        var d = new Dataset(name, genes, h[first..], values.ToArray());
        if (lengthCol >= 0) d.RowAnnotations[LengthAnnotation] = lengths.ToArray();
        return d;
    }

    /// <summary>Gene lengths from the <c>length</c> row annotation (throws if absent).</summary>
    public static double[] GeneLengths(Dataset counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!counts.RowAnnotations.TryGetValue(LengthAnnotation, out var ann))
            throw new InvalidOperationException("The dataset has no 'length' row annotation; pass gene lengths explicitly.");
        var l = new double[ann.Length];
        for (int i = 0; i < ann.Length; i++) l[i] = double.Parse(ann[i], NumberStyles.Float, CultureInfo.InvariantCulture);
        return l;
    }

    /// <summary>Stores gene lengths in the <c>length</c> row annotation.</summary>
    public static Dataset WithGeneLengths(this Dataset counts, ReadOnlySpan<double> lengths)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (lengths.Length != counts.Rows) throw new ArgumentException($"Need {counts.Rows} gene lengths.");
        var s = new string[lengths.Length];
        for (int i = 0; i < s.Length; i++) s[i] = lengths[i].ToString("R", CultureInfo.InvariantCulture);
        counts.RowAnnotations[LengthAnnotation] = s;
        return counts;
    }

    private static string CleanSampleName(string s)
    {
        string f = Path.GetFileName(s.Replace('\\', '/'));
        foreach (var ext in new[] { ".gz", ".bam", ".sam", ".cram", ".txt", ".tsv", ".counts", ".htseq" })
            if (f.EndsWith(ext, StringComparison.OrdinalIgnoreCase)) f = f[..^ext.Length];
        return f;
    }

    private static double ParseCount(string s, int lineNo)
        => double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out double x)
            ? x
            : throw new InvalidDataException($"Line {lineNo}: '{s}' is not a number.");
}
