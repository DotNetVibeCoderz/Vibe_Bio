using System.Globalization;
using Biotech.Net.Core;
using Microsoft.Win32.SafeHandles;

namespace Biotech.Net.Sequences.IO;

/// <summary>One line of a samtools <c>.fai</c> index.</summary>
/// <param name="Name">Sequence name (header up to the first whitespace).</param>
/// <param name="Length">Number of residues.</param>
/// <param name="Offset">Byte offset of the first residue in the file.</param>
/// <param name="LineBases">Residues per full line.</param>
/// <param name="LineWidth">Bytes per full line, including the line terminator.</param>
public sealed record FastaIndexEntry(string Name, long Length, long Offset, int LineBases, int LineWidth)
{
    /// <summary>File offset of the 0-based residue <paramref name="position"/>.</summary>
    public long OffsetOf(long position) => LineBases == 0 ? Offset : Offset + position / LineBases * LineWidth + position % LineBases;
}

/// <summary>
/// A samtools-compatible FASTA index (<c>.fai</c>). Building requires every line of a record except the last
/// to have the same length (the same rule samtools faidx enforces); that regularity is what makes the
/// residue → byte offset mapping a constant-time computation.
/// </summary>
public sealed class FastaIndex
{
    private readonly Dictionary<string, FastaIndexEntry> _byName;

    /// <summary>Creates an index from entries.</summary>
    public FastaIndex(IEnumerable<FastaIndexEntry> entries)
    {
        Entries = entries.ToList();
        _byName = new(Entries.Count, StringComparer.Ordinal);
        foreach (var e in Entries)
            if (!_byName.TryAdd(e.Name, e)) throw new InvalidDataException($"Duplicate sequence name '{e.Name}' in FASTA index.");
    }

    /// <summary>Entries in file order.</summary>
    public IReadOnlyList<FastaIndexEntry> Entries { get; }

    /// <summary>The entry for <paramref name="name"/>.</summary>
    public FastaIndexEntry this[string name] => _byName.TryGetValue(name, out var e) ? e
        : throw new KeyNotFoundException($"Sequence '{name}' is not in the FASTA index.");

    /// <summary>Looks up an entry.</summary>
    public bool TryGet(string name, out FastaIndexEntry entry) => _byName.TryGetValue(name, out entry!);

    /// <summary>Builds the index by scanning an uncompressed FASTA file byte by byte.</summary>
    /// <exception cref="InvalidDataException">The file is gzip-compressed or a record has uneven line lengths.</exception>
    public static FastaIndex Build(string fastaPath)
    {
        using var fs = new FileStream(fastaPath, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16, FileOptions.SequentialScan);
        return Build(fs);
    }

    /// <summary>Builds the index from a seekable or sequential stream positioned at the start of the FASTA data.</summary>
    public static FastaIndex Build(Stream stream)
    {
        var entries = new List<FastaIndexEntry>();
        var buf = new byte[1 << 16];
        long filePos = 0;
        // Per-record state.
        string? name = null;
        long length = 0, offset = 0;
        int lineBases = -1, lineWidth = -1;
        bool sawShortLine = false;
        // Per-line state.
        bool inHeader = false, atLineStart = true;
        var headerBytes = new List<byte>(128);
        int curBases = 0, curWidth = 0;
        bool firstByte = true;

        void EndSequenceLine(long lineNo)
        {
            if (name is null || curWidth == 0) return;
            if (curBases == 0) { curBases = 0; curWidth = 0; return; } // blank line inside record: tolerated only at the end
            if (sawShortLine)
                throw new InvalidDataException($"Different line length in sequence '{name}'; the FASTA file cannot be indexed.");
            if (lineBases < 0) { lineBases = curBases; lineWidth = curWidth; }
            else if (curBases != lineBases || curWidth != lineWidth)
            {
                if (curBases > lineBases) throw new InvalidDataException($"Different line length in sequence '{name}'; the FASTA file cannot be indexed.");
                sawShortLine = true;
            }
            length += curBases;
            curBases = 0; curWidth = 0;
        }

        void Finish()
        {
            if (name is null) return;
            entries.Add(new FastaIndexEntry(name, length, offset, Math.Max(lineBases, 0), Math.Max(lineWidth, 0)));
        }

        int read;
        while ((read = stream.Read(buf, 0, buf.Length)) > 0)
        {
            if (firstByte && read >= 2 && buf[0] == 0x1f && buf[1] == 0x8b)
                throw new InvalidDataException("Cannot index a gzip-compressed FASTA file; decompress it first.");
            firstByte = false;
            for (int i = 0; i < read; i++, filePos++)
            {
                byte b = buf[i];
                if (atLineStart)
                {
                    atLineStart = false;
                    if (b == (byte)'>')
                    {
                        Finish();
                        inHeader = true;
                        headerBytes.Clear();
                        continue;
                    }
                }
                if (inHeader)
                {
                    if (b == (byte)'\n')
                    {
                        inHeader = false;
                        atLineStart = true;
                        var (id, _) = FastaReader.SplitHeader(System.Text.Encoding.ASCII.GetString(headerBytes.ToArray()).TrimEnd('\r'));
                        name = id;
                        length = 0; offset = filePos + 1; lineBases = -1; lineWidth = -1; sawShortLine = false;
                        curBases = 0; curWidth = 0;
                    }
                    else headerBytes.Add(b);
                    continue;
                }
                curWidth++;
                if (b == (byte)'\n')
                {
                    EndSequenceLine(0);
                    atLineStart = true;
                }
                else if (b != (byte)'\r' && b > (byte)' ') curBases++;
            }
        }
        if (curWidth > 0)
        {
            // Last line without terminator: count its bases without imposing the width rule.
            if (name is not null && curBases > 0)
            {
                if (sawShortLine) throw new InvalidDataException($"Different line length in sequence '{name}'.");
                if (lineBases < 0) { lineBases = curBases; lineWidth = curWidth + 1; }
                length += curBases;
            }
        }
        Finish();
        return new FastaIndex(entries);
    }

    /// <summary>Writes the five tab-separated <c>.fai</c> columns.</summary>
    public void Write(TextWriter writer)
    {
        foreach (var e in Entries)
            writer.Write($"{e.Name}\t{e.Length}\t{e.Offset}\t{e.LineBases}\t{e.LineWidth}\n");
    }

    /// <summary>Saves the index to <paramref name="path"/> (conventionally the FASTA path plus <c>.fai</c>).</summary>
    public void Save(string path)
    {
        using var w = BioFile.CreateText(path);
        Write(w);
    }

    /// <summary>Reads a <c>.fai</c> file.</summary>
    public static FastaIndex Load(string path)
    {
        using var r = new StreamReader(path);
        return Read(r);
    }

    /// <summary>Parses <c>.fai</c> content.</summary>
    public static FastaIndex Read(TextReader reader)
    {
        var entries = new List<FastaIndexEntry>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            if (f.Length < 5) throw new InvalidDataException($"Malformed .fai line: '{line}'.");
            entries.Add(new FastaIndexEntry(f[0], long.Parse(f[1], CultureInfo.InvariantCulture), long.Parse(f[2], CultureInfo.InvariantCulture),
                int.Parse(f[3], CultureInfo.InvariantCulture), int.Parse(f[4], CultureInfo.InvariantCulture)));
        }
        return new FastaIndex(entries);
    }
}

/// <summary>
/// Random access to regions of a large FASTA file through its <see cref="FastaIndex"/>, like <c>samtools faidx</c>.
/// Only the bytes of the requested region are read, using positional reads, so one instance can serve
/// concurrent requests from several threads.
/// </summary>
public sealed class IndexedFastaReader : IDisposable
{
    private readonly SafeFileHandle _handle;

    private IndexedFastaReader(SafeFileHandle handle, FastaIndex index, string path)
    {
        _handle = handle;
        Index = index;
        Path = path;
    }

    /// <summary>
    /// Opens <paramref name="fastaPath"/>, loading <c>fastaPath.fai</c> when present and building it otherwise
    /// (saving it next to the FASTA when <paramref name="saveIndex"/> is true).
    /// </summary>
    public static IndexedFastaReader Open(string fastaPath, bool saveIndex = true)
    {
        string faiPath = fastaPath + ".fai";
        FastaIndex index;
        if (File.Exists(faiPath)) index = FastaIndex.Load(faiPath);
        else
        {
            index = FastaIndex.Build(fastaPath);
            if (saveIndex)
            {
                try { index.Save(faiPath); }
                catch (IOException) { } // read-only location: keep the in-memory index
                catch (UnauthorizedAccessException) { }
            }
        }
        var handle = File.OpenHandle(fastaPath, FileMode.Open, FileAccess.Read, FileShare.Read, FileOptions.RandomAccess);
        return new IndexedFastaReader(handle, index, fastaPath);
    }

    /// <summary>The FASTA file path.</summary>
    public string Path { get; }

    /// <summary>The index in use.</summary>
    public FastaIndex Index { get; }

    /// <summary>Sequence names in file order.</summary>
    public IEnumerable<string> Names => Index.Entries.Select(e => e.Name);

    /// <summary>Length of sequence <paramref name="name"/>.</summary>
    public long LengthOf(string name) => Index[name].Length;

    /// <summary>Residues <c>[start, end)</c> (0-based, half-open) of sequence <paramref name="name"/>; end is clamped to the length.</summary>
    public Sequence Fetch(string name, long start, long end)
    {
        var e = Index[name];
        if (start < 0 || start > end) throw new ArgumentOutOfRangeException(nameof(start), $"Invalid region [{start}, {end}).");
        end = Math.Min(end, e.Length);
        start = Math.Min(start, end);
        int n = checked((int)(end - start));
        var result = new byte[n];
        if (n > 0)
        {
            long first = e.OffsetOf(start), last = e.OffsetOf(end - 1);
            var raw = new byte[checked((int)(last - first + 1))];
            int got = 0;
            while (got < raw.Length)
            {
                int r = RandomAccess.Read(_handle, raw.AsSpan(got), first + got);
                if (r == 0) throw new EndOfStreamException($"FASTA file is shorter than its index claims for '{name}'.");
                got += r;
            }
            int o = 0;
            foreach (byte b in raw)
                if (b > (byte)' ') result[o++] = b;
            if (o != n) throw new InvalidDataException($"Index does not match the FASTA file for '{name}' (got {o} of {n} residues).");
        }
        string id = start == 0 && end == e.Length ? name : $"{name}:{start + 1}-{end}";
        return Sequence.FromBuffer(id, result, Alphabet.Detect(result), validate: false);
    }

    /// <summary>The bases of <paramref name="region"/>, reverse-complemented when it lies on the reverse strand.</summary>
    public Sequence Fetch(GenomicInterval region)
    {
        var s = Fetch(region.Chromosome, region.Start, region.End);
        return region.Strand == Strand.Reverse ? s.ReverseComplement() : s;
    }

    /// <summary>Parses a samtools-style region (<c>chr1</c>, <c>chr1:100-200</c>, 1-based closed) and fetches it.</summary>
    public Sequence Fetch(string region)
    {
        if (Index.TryGet(region, out var whole)) return Fetch(whole.Name, 0, whole.Length);
        int colon = region.LastIndexOf(':');
        string name = colon < 0 ? region : region[..colon];
        var iv = GenomicInterval.Parse(region, Index[name].Length);
        return Fetch(iv.Chromosome, iv.Start, Math.Min(iv.End, Index[name].Length));
    }

    /// <summary>The whole sequence <paramref name="name"/>.</summary>
    public Sequence FetchAll(string name) => Fetch(name, 0, Index[name].Length);

    /// <summary>A genome whose chromosomes carry their lengths (and, with <paramref name="loadSequences"/>, their bases).</summary>
    public Genome ToGenome(string? name = null, bool loadSequences = false)
    {
        var g = new Genome(name ?? System.IO.Path.GetFileNameWithoutExtension(Path));
        foreach (var e in Index.Entries)
            g.Add(loadSequences ? Chromosome.FromSequence(FetchAll(e.Name)) : new Chromosome(e.Name, e.Length));
        return g;
    }

    public void Dispose() => _handle.Dispose();
}
