using System.IO.Compression;
using Biotech.Net.Genomics.IO;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

public class BgzfTests
{
    [Fact]
    public void Crc32MatchesTheStandardCheckValue()
    {
        // CRC-32/ISO-HDLC check value of "123456789"
        Assert.Equal(0xCBF43926u, Crc32.Compute("123456789"u8));
        Assert.Equal(0u, Crc32.Compute([]));
        var data = new byte[1000];
        new Random(1).NextBytes(data);
        Assert.Equal(Crc32.Compute(data), Crc32.Append(Crc32.Compute(data.AsSpan(0, 333)), data.AsSpan(333)));
    }

    [Fact]
    public void ReadsBgzipOutputLikeGzip()
    {
        // variants.vcf.gz was written by htslib's bgzip; the framework's GZipStream is the independent reference
        var raw = File.ReadAllBytes(Hts("variants.vcf.gz"));
        var expected = Gunzip(raw);
        using var bg = new BgzfReader(new MemoryStream(raw));
        using var ms = new MemoryStream();
        bg.CopyTo(ms);
        Assert.Equal(expected, ms.ToArray());
        Assert.True(expected.Length > 2 * BgzfWriter.BlockSize, "fixture should span several blocks");
        Assert.True(BgzfReader.IsBgzf(Hts("variants.vcf.gz")));
    }

    [Fact]
    public void DetectsCorruptedBlocks()
    {
        var raw = File.ReadAllBytes(Hts("features.bed.gz"));
        // flip a bit in the stored CRC32 of the first block
        int bsize = BitConverter.ToUInt16(raw, 16) + 1;
        raw[bsize - 8] ^= 1;
        using var bg = new BgzfReader(new MemoryStream(raw));
        Assert.Throws<InvalidDataException>(() => bg.CopyTo(Stream.Null));

        var plainGzip = new MemoryStream();
        using (var gz = new GZipStream(plainGzip, CompressionLevel.Fastest, leaveOpen: true)) gz.Write("hello"u8);
        plainGzip.Position = 0;
        Assert.False(BgzfReader.IsBgzf(plainGzip));
        using var notBgzf = new BgzfReader(plainGzip);
        Assert.Throws<InvalidDataException>(() => notBgzf.ReadByte());
    }

    [Fact]
    public void WriterProducesValidBgzfWithEofMarker()
    {
        var data = new byte[300_000];
        var rng = new Random(7);
        for (int i = 0; i < data.Length; i++) data[i] = (byte)"ACGT\n"[rng.Next(5)];
        var ms = new MemoryStream();
        var offsets = new List<(VirtualOffset Offset, int Index)>();
        using (var w = new BgzfWriter(ms, leaveOpen: true))
        {
            for (int i = 0; i < data.Length; i += 7919)
            {
                offsets.Add((w.VirtualPosition, i));
                w.Write(data.AsSpan(i, Math.Min(7919, data.Length - i)));
            }
        }
        var bytes = ms.ToArray();
        Assert.Equal(BgzfWriter.EofMarker.ToArray(), bytes[^28..]);
        Assert.Equal(data, Gunzip(bytes)); // any gzip reader accepts BGZF
        // every block fits the 64 KiB limit and declares its size in the BC field
        for (int p = 0; p < bytes.Length;)
        {
            Assert.Equal((byte)'B', bytes[p + 12]);
            Assert.Equal((byte)'C', bytes[p + 13]);
            int bsize = BitConverter.ToUInt16(bytes, p + 16) + 1;
            Assert.InRange(bsize, 28, 65536);
            p += bsize;
        }
        // random access through the recorded virtual offsets
        using var r = new BgzfReader(new MemoryStream(bytes));
        var buf = new byte[100];
        foreach (var (off, idx) in offsets.AsEnumerable().Reverse())
        {
            r.Seek(off);
            int n = r.ReadAtLeast(buf, Math.Min(100, data.Length - idx), throwOnEndOfStream: false);
            Assert.Equal(data.AsSpan(idx, n).ToArray(), buf[..n]);
        }
    }

    [Fact]
    public void VirtualPositionAtABlockEndPointsToTheNextBlock()
    {
        var ms = new MemoryStream();
        using (var w = new BgzfWriter(ms, leaveOpen: true))
        {
            w.Write(new byte[10]);
            w.FlushBlock();
            w.Write(new byte[5]);
        }
        using var r = new BgzfReader(new MemoryStream(ms.ToArray()));
        r.ReadExactly(new byte[10]);
        var pos = r.VirtualPosition;
        Assert.Equal(0, pos.BlockOffset);
        Assert.True(pos.BlockAddress > 0); // htslib's bgzf_tell convention
        r.ReadExactly(new byte[5]);
        Assert.Equal(-1, r.ReadByte());
        Assert.True(r.AtEnd);
    }
}
