using System.Buffers;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>ASCII offset of FASTQ quality characters.</summary>
public enum FastqEncoding
{
    /// <summary>Detect from the first records (see <see cref="FastqReader.DetectEncoding(IEnumerable{string})"/>).</summary>
    Auto = 0,

    /// <summary>Sanger / Illumina 1.8+ (Phred + 33, <c>!</c> = Q0).</summary>
    Phred33 = 33,

    /// <summary>Illumina 1.3–1.7 (Phred + 64, <c>@</c> = Q0).</summary>
    Phred64 = 64,
}

/// <summary>
/// Streaming FASTQ reader. Tolerates multi-line records (sequence and quality wrapped over several lines):
/// quality lines are consumed until they cover the sequence length, because a quality line may itself begin
/// with <c>@</c> and cannot be told apart from a header by its first character.
/// </summary>
public static class FastqReader
{
    // How many records are buffered while the encoding is still ambiguous.
    private const int DetectionWindow = 1000;

    /// <summary>Reads every record of <paramref name="reader"/>.</summary>
    /// <exception cref="InvalidDataException">The data is not well-formed FASTQ.</exception>
    public static IEnumerable<SequenceRead> Read(TextReader reader, FastqEncoding encoding = FastqEncoding.Auto)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var raw = ReadRaw(reader);
        if (encoding != FastqEncoding.Auto)
        {
            foreach (var r in raw) yield return r.ToRead((int)encoding);
            yield break;
        }

        // Buffer records until the encoding is decided (a character below ';' proves Phred+33).
        var pending = new List<RawRecord>();
        int offset = 0;
        using var e = raw.GetEnumerator();
        int min = int.MaxValue, max = 0;
        while (e.MoveNext())
        {
            var r = e.Current;
            pending.Add(r);
            foreach (char c in r.Quality) { if (c < min) min = c; if (c > max) max = c; }
            if (min < 59) { offset = 33; break; }
            if (pending.Count >= DetectionWindow) break;
        }
        if (offset == 0) offset = Decide(min, max);
        foreach (var r in pending) yield return r.ToRead(offset);
        while (e.MoveNext()) yield return e.Current.ToRead(offset);
    }

    /// <summary>Reads a FASTQ file (optionally gzip-compressed) lazily.</summary>
    public static IEnumerable<SequenceRead> Read(string path, FastqEncoding encoding = FastqEncoding.Auto)
    {
        using var reader = BioFile.OpenText(path);
        foreach (var r in Read(reader, encoding)) yield return r;
    }

    /// <summary>Parses FASTQ text held in memory.</summary>
    public static IReadOnlyList<SequenceRead> Parse(string text, FastqEncoding encoding = FastqEncoding.Auto)
        => Read(new StringReader(text), encoding).ToList();

    /// <summary>
    /// Guesses the encoding from quality strings: any character below <c>;</c> (59) means Phred+33; if every
    /// character is at least <c>@</c> (64) and some exceed <c>J</c> (74, the Phred+33 Illumina maximum) it is
    /// Phred+64; otherwise Phred+33, the modern default.
    /// </summary>
    public static FastqEncoding DetectEncoding(IEnumerable<string> qualityStrings)
    {
        int min = int.MaxValue, max = 0;
        foreach (var q in qualityStrings)
            foreach (char c in q) { if (c < min) min = c; if (c > max) max = c; }
        return (FastqEncoding)Decide(min, max);
    }

    /// <summary>Detects the encoding of a FASTQ file from its first records.</summary>
    public static FastqEncoding DetectEncoding(string path)
    {
        using var reader = BioFile.OpenText(path);
        return DetectEncoding(ReadRaw(reader).Take(DetectionWindow).Select(r => r.Quality));
    }

    private static int Decide(int min, int max)
    {
        if (min == int.MaxValue || min < 59) return 33;
        return min >= 64 && max > 74 ? 64 : 33;
    }

    internal readonly record struct RawRecord(string Header, byte[] Residues, string Quality, long Line)
    {
        public SequenceRead ToRead(int offset)
        {
            var (id, description) = FastaReader.SplitHeader(Header.AsSpan(1));
            var q = new byte[Quality.Length];
            for (int i = 0; i < q.Length; i++)
            {
                int v = Quality[i] - offset;
                if (v < 0 || v > 93)
                    throw new InvalidDataException($"Record '{id}' (line {Line}): quality character '{Quality[i]}' is out of range for Phred+{offset}.");
                q[i] = (byte)v;
            }
            var alphabet = Alphabet.DnaIupac.IsValid(Residues) ? Alphabet.DnaIupac : Alphabet.Detect(Residues);
            return SequenceRead.FromBuffer(Sequence.FromBuffer(id, Residues, alphabet, description, validate: false), q);
        }
    }

    internal static IEnumerable<RawRecord> ReadRaw(TextReader reader)
    {
        var seq = new ByteBuilder(512);
        var qual = new System.Text.StringBuilder(512);
        string? line;
        long lineNo = 0;
        while (true)
        {
            string? header = null;
            if (header is null)
            {
                do { line = reader.ReadLine(); lineNo++; } while (line is not null && line.Length == 0);
                if (line is null) yield break;
                header = line;
            }
            long headerLine = lineNo;
            if (header[0] != '@') throw new InvalidDataException($"Line {lineNo}: expected a FASTQ header starting with '@'.");
            seq.Clear();
            while (true)
            {
                line = reader.ReadLine(); lineNo++;
                if (line is null) throw new InvalidDataException($"Record at line {headerLine} ends before its '+' separator.");
                if (line.Length > 0 && line[0] == '+') break;
                seq.AppendResidues(line);
            }
            qual.Clear();
            while (qual.Length < seq.Length)
            {
                line = reader.ReadLine(); lineNo++;
                if (line is null) break;
                qual.Append(line.AsSpan().TrimEnd());
            }
            if (qual.Length != seq.Length)
                throw new InvalidDataException($"Record at line {headerLine}: {seq.Length} bases but {qual.Length} quality values.");
            yield return new RawRecord(header, seq.ToArray(), qual.ToString(), headerLine);
        }
    }
}

/// <summary>FASTQ writer (four-line records).</summary>
public static class FastqWriter
{
    /// <summary>Writes one read with qualities encoded per <paramref name="encoding"/> (Phred+33 by default).</summary>
    public static void Write(TextWriter writer, SequenceRead read, FastqEncoding encoding = FastqEncoding.Phred33)
    {
        int offset = encoding == FastqEncoding.Phred64 ? 64 : 33;
        writer.Write('@');
        writer.Write(read.Id);
        if (read.Sequence.Description.Length > 0) { writer.Write(' '); writer.Write(read.Sequence.Description); }
        writer.Write('\n');
        FastaWriter.WriteResidues(writer, read.Sequence.Span, 0);
        if (read.Length == 0) writer.Write('\n');
        writer.Write("+\n");
        char[] q = ArrayPool<char>.Shared.Rent(Math.Max(1, read.Length));
        try
        {
            var qs = read.Quality;
            for (int i = 0; i < qs.Length; i++) q[i] = (char)(qs[i] + offset);
            writer.Write(q, 0, qs.Length);
        }
        finally { ArrayPool<char>.Shared.Return(q); }
        writer.Write('\n');
    }

    /// <summary>Writes many reads.</summary>
    public static void Write(TextWriter writer, IEnumerable<SequenceRead> reads, FastqEncoding encoding = FastqEncoding.Phred33)
    {
        foreach (var r in reads) Write(writer, r, encoding);
    }

    /// <summary>Writes reads to <paramref name="path"/> (gzip when it ends in <c>.gz</c>).</summary>
    public static void Write(string path, IEnumerable<SequenceRead> reads, FastqEncoding encoding = FastqEncoding.Phred33)
    {
        using var w = BioFile.CreateText(path);
        Write(w, reads, encoding);
    }
}
