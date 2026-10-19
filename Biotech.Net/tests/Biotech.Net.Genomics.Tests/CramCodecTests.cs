using System.Text;
using Biotech.Net.Genomics.Cram;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

/// <summary>
/// The GA4GH hts-specs codec test vectors (test/cram/codecs, produced by htscodecs): every compressed variant must
/// decode to the original data, processed as the README describes.
/// </summary>
public class CramCodecTests
{
    private static string Codec(string sub, string file) => Spec(Path.Combine("codecs", sub, file));

    /// <summary>The first whitespace-separated column of every line, concatenated (the README's perl one-liner).</summary>
    private static byte[] FirstColumns(string gz)
    {
        var sb = new StringBuilder();
        foreach (var line in ReadText(Codec("gzip", gz)).Split('\n'))
        {
            var f = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
            if (f.Length > 0) sb.Append(f[0]);
        }
        return Encoding.Latin1.GetBytes(sb.ToString());
    }

    private static byte[] Raw(string gz) => Gunzip(File.ReadAllBytes(Codec("gzip", gz)));

    [Theory]
    [InlineData("q4.0")]
    [InlineData("q4.1")]
    public void Rans4x8(string file) => Assert.Equal(FirstColumns("q4.gz"), Cram.Rans4x8.Decode(File.ReadAllBytes(Codec("rans4x8", file))));

    [Theory]
    [InlineData("q4.0")] [InlineData("q4.1")] [InlineData("q4.4")] [InlineData("q4.5")] [InlineData("q4.64")]
    [InlineData("q4.65")] [InlineData("q4.128")] [InlineData("q4.129")] [InlineData("q4.192")] [InlineData("q4.193")]
    public void RansNx16(string file) => Assert.Equal(FirstColumns("q4.gz"), Cram.RansNx16.Decode(File.ReadAllBytes(Codec("ransNx16", file))));

    [Theory]
    [InlineData("u32.1")]
    [InlineData("u32.9")]
    public void RansNx16Stripe(string file) => Assert.Equal(Raw("u32.gz"), Cram.RansNx16.Decode(File.ReadAllBytes(Codec("ransNx16", file))));

    [Theory]
    [InlineData("q4.0")] [InlineData("q4.1")] [InlineData("q4.8")] [InlineData("q4.9")] [InlineData("q4.64")]
    [InlineData("q4.65")] [InlineData("q4.128")] [InlineData("q4.129")] [InlineData("q4.192")] [InlineData("q4.193")]
    public void ArithmeticCoder(string file) => Assert.Equal(FirstColumns("q4.gz"), ArithDynamic.Decode(File.ReadAllBytes(Codec("range", file))));

    [Theory]
    [InlineData("u32.1")]
    [InlineData("u32.9")]
    public void ArithmeticCoderStripe(string file) => Assert.Equal(Raw("u32.gz"), ArithDynamic.Decode(File.ReadAllBytes(Codec("range", file))));

    [Theory]
    [InlineData("q4.0", "q4.gz")] [InlineData("q4.1", "q4.gz")] [InlineData("q4.2", "q4.gz")] [InlineData("q4.3", "q4.gz")]
    [InlineData("q40+dir.3", "q40+dir.gz")]
    public void FqzComp(string file, string original)
    {
        // CRAM stores raw Phred values; the htscodecs test tool reads and prints them as Phred+33 text
        var decoded = Cram.FqzComp.Decode(File.ReadAllBytes(Codec("fqzcomp", file))).Select(q => (byte)(q + 33)).ToArray();
        Assert.Equal(FirstColumns(original), decoded);
    }

    [Theory]
    [InlineData("01.names.1")] [InlineData("01.names.3")] [InlineData("01.names.5")] [InlineData("01.names.7")]
    [InlineData("01.names.9")] [InlineData("01.names.11")] [InlineData("01.names.13")] [InlineData("01.names.15")]
    [InlineData("01.names.17")] [InlineData("01.names.19")] [InlineData("20.names.1")] [InlineData("20.names.9")]
    [InlineData("20.names.19")]
    public void NameTokeniser(string file)
    {
        var decoded = NameTokenizer3.Decode(File.ReadAllBytes(Codec("tok3", file)));
        for (int i = 0; i < decoded.Length; i++) if (decoded[i] == 0) decoded[i] = (byte)'\n';
        Assert.Equal(Raw(file[..file.LastIndexOf('.')] + ".gz"), decoded);
    }
}
