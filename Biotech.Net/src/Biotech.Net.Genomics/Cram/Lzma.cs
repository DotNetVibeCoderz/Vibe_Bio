namespace Biotech.Net.Genomics.Cram;

/// <summary>LZMA (xz container) block decompression — placeholder until implemented.</summary>
internal static class Lzma
{
    public static byte[] DecompressXz(ReadOnlySpan<byte> data, int expectedLength)
        => throw new NotSupportedException("CRAM block compression method 3 (LZMA/xz) is not supported by Biotech.Net; re-encode the file without use_lzma.");
}
