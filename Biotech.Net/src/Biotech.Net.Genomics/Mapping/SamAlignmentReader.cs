using System.Text;
using Biotech.Net.Genomics.Cram;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>
/// SAM text reader producing <see cref="AlignmentRecord"/>s (plain or gzip/BGZF). The header is read on
/// construction. Region queries scan the file linearly (SAM has no index here); use BAM or CRAM for random access.
/// </summary>
public sealed class SamAlignmentReader : IAlignmentReader
{
    private readonly TextReader _reader;
    private string? _pending;

    /// <summary>Reads the header from <paramref name="reader"/>.</summary>
    public SamAlignmentReader(TextReader reader)
    {
        _reader = reader ?? throw new ArgumentNullException(nameof(reader));
        var sb = new StringBuilder();
        string? line;
        while ((line = _reader.ReadLine()) is not null)
        {
            if (line.Length > 0 && line[0] == '@') { sb.Append(line).Append('\n'); continue; }
            if (line.Length > 0) _pending = line;
            break;
        }
        Header = new SamHeader(sb.ToString());
    }

    /// <summary>Opens a SAM file (gzip detected from the content).</summary>
    public static SamAlignmentReader Open(string path) => new(BioFile.OpenText(path));

    /// <inheritdoc/>
    public SamHeader Header { get; }

    /// <inheritdoc/>
    public bool IsIndexed => false;

    /// <inheritdoc/>
    public IEnumerable<AlignmentRecord> Records()
    {
        if (_pending is not null)
        {
            var p = _pending;
            _pending = null;
            yield return AlignmentRecord.Parse(p, Header);
        }
        string? line;
        while ((line = _reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '@') continue;
            yield return AlignmentRecord.Parse(line, Header);
        }
    }

    /// <inheritdoc/>
    public IEnumerable<AlignmentRecord> Query(string referenceName, long start, long end)
        => Records().Where(r => r.ReferenceName == referenceName && r.Overlaps(start, end));

    /// <summary>Closes the underlying reader.</summary>
    public void Dispose() => _reader.Dispose();
}

/// <summary>Opens SAM, BAM or CRAM files by content (not by extension).</summary>
public static class AlignmentFile
{
    /// <summary>
    /// Opens <paramref name="path"/> as BAM (BGZF with the <c>BAM\1</c> magic), CRAM (<c>CRAM</c> magic) or SAM text.
    /// <paramref name="referencePath"/> is the FASTA (with <c>.fai</c>) used to decode reference-based CRAM.
    /// Indexes beside the file (<c>.bai</c>/<c>.csi</c>/<c>.crai</c>) are loaded when present.
    /// </summary>
    public static IAlignmentReader Open(string path, string? referencePath = null)
    {
        Span<byte> head = stackalloc byte[4];
        int n;
        using (var fs = File.OpenRead(path)) n = fs.ReadAtLeast(head, 4, throwOnEndOfStream: false);
        if (n == 4 && head.SequenceEqual("CRAM"u8)) return CramReader.Open(path, referencePath);
        if (n >= 2 && head[0] == 0x1f && head[1] == 0x8b && BgzfReader.IsBgzf(path))
        {
            using var bg = BgzfReader.Open(path);
            Span<byte> magic = stackalloc byte[4];
            if (bg.ReadAtLeast(magic, 4, throwOnEndOfStream: false) == 4 && magic.SequenceEqual("BAM\u0001"u8))
                return BamReader.Open(path);
        }
        return SamAlignmentReader.Open(path);
    }
}
