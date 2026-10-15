using System.Buffers;
using System.Runtime.CompilerServices;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>
/// Streaming FASTA reader. Records are yielded one at a time, so a multi-gigabyte file never needs to fit in memory;
/// each record's residues are accumulated in a reusable buffer and handed to <see cref="Sequence.FromBuffer"/>.
/// </summary>
public static class FastaReader
{
    /// <summary>
    /// Reads every record of <paramref name="reader"/>. The identifier is the header up to the first whitespace and
    /// the rest is the description. With <paramref name="alphabet"/> null the type is detected per record.
    /// Lines starting with <c>;</c> (old-style comments) and blank lines are ignored.
    /// </summary>
    public static IEnumerable<Sequence> Read(TextReader reader, Alphabet? alphabet = null, bool validate = false)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var buffer = new ByteBuilder(4096);
        string? header = null;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == ';') continue;
            if (line[0] == '>')
            {
                if (header is not null) yield return Build(header, buffer, alphabet, validate);
                header = line;
                buffer.Clear();
                continue;
            }
            if (header is null) throw new InvalidDataException("FASTA data must start with a '>' header line.");
            buffer.AppendResidues(line);
        }
        if (header is not null) yield return Build(header, buffer, alphabet, validate);
    }

    /// <summary>Reads a FASTA file (optionally gzip-compressed) lazily; the file stays open while enumerating.</summary>
    public static IEnumerable<Sequence> Read(string path, Alphabet? alphabet = null, bool validate = false)
    {
        using var reader = BioFile.OpenText(path);
        foreach (var s in Read(reader, alphabet, validate)) yield return s;
    }

    /// <summary>Parses FASTA text held in memory.</summary>
    public static IReadOnlyList<Sequence> Parse(string text, Alphabet? alphabet = null)
        => Read(new StringReader(text), alphabet).ToList();

    /// <summary>Asynchronously reads records from a stream (gzip detected automatically).</summary>
    public static async IAsyncEnumerable<Sequence> ReadAsync(Stream stream, Alphabet? alphabet = null, bool validate = false,
        [EnumeratorCancellation] CancellationToken cancellationToken = default)
    {
        using var reader = new StreamReader(BioFile.WrapDecompression(stream), Encoding.UTF8, true, 1 << 16, leaveOpen: false);
        var buffer = new ByteBuilder(4096);
        string? header = null;
        string? line;
        while ((line = await reader.ReadLineAsync(cancellationToken).ConfigureAwait(false)) is not null)
        {
            if (line.Length == 0 || line[0] == ';') continue;
            if (line[0] == '>')
            {
                if (header is not null) yield return Build(header, buffer, alphabet, validate);
                header = line;
                buffer.Clear();
                continue;
            }
            if (header is null) throw new InvalidDataException("FASTA data must start with a '>' header line.");
            buffer.AppendResidues(line);
        }
        if (header is not null) yield return Build(header, buffer, alphabet, validate);
    }

    internal static (string Id, string Description) SplitHeader(ReadOnlySpan<char> header)
    {
        header = header.Trim();
        int ws = header.IndexOfAny(' ', '\t');
        if (ws < 0) return (header.ToString(), string.Empty);
        return (header[..ws].ToString(), header[(ws + 1)..].Trim().ToString());
    }

    private static Sequence Build(string header, ByteBuilder buffer, Alphabet? alphabet, bool validate)
    {
        var (id, description) = SplitHeader(header.AsSpan(1));
        var bytes = buffer.ToArray();
        var a = alphabet ?? Alphabet.Detect(bytes);
        return Sequence.FromBuffer(id, bytes, a, description, validate && alphabet is not null);
    }
}

/// <summary>FASTA writer with configurable line wrapping.</summary>
public static class FastaWriter
{
    /// <summary>Default residues per line (the NCBI convention).</summary>
    public const int DefaultLineWidth = 60;

    /// <summary>Writes one record. <paramref name="lineWidth"/> ≤ 0 writes the residues on a single line.</summary>
    public static void Write(TextWriter writer, Sequence sequence, int lineWidth = DefaultLineWidth)
    {
        writer.Write('>');
        writer.Write(sequence.Id);
        if (sequence.Description.Length > 0) { writer.Write(' '); writer.Write(sequence.Description); }
        writer.Write('\n');
        WriteResidues(writer, sequence.Span, lineWidth);
    }

    /// <summary>Writes many records.</summary>
    public static void Write(TextWriter writer, IEnumerable<Sequence> sequences, int lineWidth = DefaultLineWidth)
    {
        foreach (var s in sequences) Write(writer, s, lineWidth);
    }

    /// <summary>Writes records to <paramref name="path"/> (gzip when it ends in <c>.gz</c>).</summary>
    public static void Write(string path, IEnumerable<Sequence> sequences, int lineWidth = DefaultLineWidth)
    {
        using var w = BioFile.CreateText(path);
        Write(w, sequences, lineWidth);
    }

    /// <summary>FASTA text for the given records.</summary>
    public static string ToString(IEnumerable<Sequence> sequences, int lineWidth = DefaultLineWidth)
    {
        var sw = new StringWriter { NewLine = "\n" };
        Write(sw, sequences, lineWidth);
        return sw.ToString();
    }

    internal static void WriteResidues(TextWriter writer, ReadOnlySpan<byte> residues, int lineWidth)
    {
        if (residues.Length == 0) return;
        if (lineWidth <= 0) lineWidth = residues.Length;
        char[] chars = ArrayPool<char>.Shared.Rent(lineWidth);
        try
        {
            for (int i = 0; i < residues.Length; i += lineWidth)
            {
                var chunk = residues.Slice(i, Math.Min(lineWidth, residues.Length - i));
                int n = Encoding.ASCII.GetChars(chunk, chars);
                writer.Write(chars, 0, n);
                writer.Write('\n');
            }
        }
        finally { ArrayPool<char>.Shared.Return(chars); }
    }
}
