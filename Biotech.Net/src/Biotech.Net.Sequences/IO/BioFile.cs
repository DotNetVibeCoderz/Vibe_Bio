using System.IO.Compression;
using System.Text;

namespace Biotech.Net.Sequences.IO;

/// <summary>
/// File helpers shared by every reader and writer: transparent gzip (detected by magic bytes, not by
/// extension, so a mislabelled file still opens) and large buffers for sequential scans.
/// </summary>
public static class BioFile
{
    private const int BufferSize = 1 << 16;

    /// <summary>Opens <paramref name="path"/> for reading, decompressing gzip (including multi-member BGZF) on the fly.</summary>
    public static Stream OpenRead(string path)
    {
        var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, BufferSize, FileOptions.SequentialScan);
        return WrapDecompression(fs);
    }

    /// <summary>Wraps <paramref name="stream"/> in a <see cref="GZipStream"/> when it starts with the gzip magic bytes.</summary>
    /// <remarks>Non-seekable streams get the two peeked bytes replayed in front of them, so nothing is lost.</remarks>
    public static Stream WrapDecompression(Stream stream)
    {
        if (!stream.CanSeek) return PeekGzip(stream);
        long pos = stream.Position;
        int b1 = stream.ReadByte(), b2 = stream.ReadByte();
        stream.Position = pos;
        return b1 == 0x1f && b2 == 0x8b ? new GZipStream(stream, CompressionMode.Decompress) : stream;
    }

    private static Stream PeekGzip(Stream s)
    {
        var head = new byte[2];
        int n = s.ReadAtLeast(head, 2, throwOnEndOfStream: false);
        var joined = new PrefixedStream(head.AsMemory(0, n), s);
        return n == 2 && head[0] == 0x1f && head[1] == 0x8b ? new GZipStream(joined, CompressionMode.Decompress) : joined;
    }

    /// <summary>Opens a text reader over <paramref name="path"/> (gzip-aware, ASCII/UTF-8).</summary>
    public static TextReader OpenText(string path) => new StreamReader(OpenRead(path), Encoding.UTF8, detectEncodingFromByteOrderMarks: true, BufferSize);

    /// <summary>Creates <paramref name="path"/> for writing; a <c>.gz</c> extension produces gzip output.</summary>
    public static Stream OpenWrite(string path)
    {
        var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, BufferSize);
        return IsGzipPath(path) ? new GZipStream(fs, CompressionLevel.Optimal) : fs;
    }

    /// <summary>Creates a text writer for <paramref name="path"/> (gzip when the name ends in <c>.gz</c>), with <c>\n</c> line endings.</summary>
    public static TextWriter CreateText(string path)
        => new StreamWriter(OpenWrite(path), new UTF8Encoding(false), BufferSize) { NewLine = "\n" };

    /// <summary>Whether the path names a gzip file (<c>.gz</c>/<c>.bgz</c>).</summary>
    public static bool IsGzipPath(string path)
        => path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase) || path.EndsWith(".bgz", StringComparison.OrdinalIgnoreCase);

    /// <summary>The extension of <paramref name="path"/> ignoring a trailing <c>.gz</c>, lower case, with the dot.</summary>
    public static string InnerExtension(string path)
    {
        if (IsGzipPath(path)) path = path[..path.LastIndexOf('.')];
        return Path.GetExtension(path).ToLowerInvariant();
    }

    /// <summary>A stream that replays already-consumed bytes before continuing with the inner stream.</summary>
    private sealed class PrefixedStream(ReadOnlyMemory<byte> prefix, Stream inner) : Stream
    {
        private ReadOnlyMemory<byte> _prefix = prefix;
        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
        public override void Flush() { }
        public override int Read(byte[] buffer, int offset, int count) => Read(buffer.AsSpan(offset, count));
        public override int Read(Span<byte> buffer)
        {
            if (_prefix.Length > 0)
            {
                int n = Math.Min(buffer.Length, _prefix.Length);
                _prefix.Span[..n].CopyTo(buffer);
                _prefix = _prefix[n..];
                return n;
            }
            return inner.Read(buffer);
        }
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
        protected override void Dispose(bool disposing) { if (disposing) inner.Dispose(); base.Dispose(disposing); }
    }
}

/// <summary>A growable byte buffer used by parsers to accumulate residues before handing ownership to a <see cref="Core.Sequence"/>.</summary>
internal sealed class ByteBuilder(int capacity = 1024)
{
    private byte[] _buf = new byte[Math.Max(16, capacity)];
    public int Length { get; private set; }

    public void Clear() => Length = 0;

    public void Append(byte b)
    {
        if (Length == _buf.Length) Array.Resize(ref _buf, _buf.Length * 2);
        _buf[Length++] = b;
    }

    /// <summary>
    /// Appends the residue characters of a line, skipping whitespace and digits (GenBank/EMBL numbering).
    /// With <paramref name="upperCase"/> the residues are upper-cased, as the GenBank and EMBL readers need:
    /// those flat files store the sequence in lower case and Biopython's parsers upper-case it.
    /// </summary>
    public void AppendResidues(ReadOnlySpan<char> line, bool skipDigits = false, bool upperCase = false)
    {
        Ensure(line.Length);
        var buf = _buf;
        int n = Length;
        foreach (char c in line)
        {
            if (c <= ' ' || c > '~') continue;
            if (skipDigits && c is >= '0' and <= '9') continue;
            buf[n++] = upperCase && c is >= 'a' and <= 'z' ? (byte)(c - 32) : (byte)c;
        }
        Length = n;
    }

    /// <summary>Appends ASCII characters verbatim (no filtering).</summary>
    public void AppendAscii(ReadOnlySpan<char> text)
    {
        Ensure(text.Length);
        for (int i = 0; i < text.Length; i++) _buf[Length + i] = (byte)text[i];
        Length += text.Length;
    }

    private void Ensure(int extra)
    {
        if (Length + extra <= _buf.Length) return;
        int size = _buf.Length;
        while (size < Length + extra) size *= 2;
        Array.Resize(ref _buf, size);
    }

    /// <summary>An exact-length copy of the content (one allocation, ownership passes to the caller).</summary>
    public byte[] ToArray() => _buf.AsSpan(0, Length).ToArray();

    public ReadOnlySpan<byte> Span => _buf.AsSpan(0, Length);
}
