namespace Biotech.Net.Genomics.IO;

/// <summary>
/// CRC-32 (ISO-HDLC / gzip / zlib polynomial 0x04C11DB7, reflected 0xEDB88320), as used by gzip, BGZF and CRAM.
/// Slicing-by-8 table implementation.
/// </summary>
public static class Crc32
{
    private static readonly uint[] Table = BuildTable();

    private static uint[] BuildTable()
    {
        var t = new uint[8 * 256];
        for (uint i = 0; i < 256; i++)
        {
            uint c = i;
            for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            t[i] = c;
        }
        for (int i = 0; i < 256; i++)
            for (int s = 1; s < 8; s++)
                t[s * 256 + i] = (t[(s - 1) * 256 + i] >> 8) ^ t[t[(s - 1) * 256 + i] & 0xFF];
        return t;
    }

    /// <summary>The CRC-32 of <paramref name="data"/>.</summary>
    public static uint Compute(ReadOnlySpan<byte> data) => Append(0, data);

    /// <summary>Continues a CRC-32 (<paramref name="crc"/> is the value returned for the preceding bytes, 0 to start).</summary>
    public static uint Append(uint crc, ReadOnlySpan<byte> data)
    {
        var t = Table;
        uint c = ~crc;
        int i = 0;
        while (data.Length - i >= 8)
        {
            uint a = c ^ (uint)(data[i] | data[i + 1] << 8 | data[i + 2] << 16 | data[i + 3] << 24);
            uint b = (uint)(data[i + 4] | data[i + 5] << 8 | data[i + 6] << 16 | data[i + 7] << 24);
            c = t[7 * 256 + (a & 0xFF)] ^ t[6 * 256 + ((a >> 8) & 0xFF)] ^ t[5 * 256 + ((a >> 16) & 0xFF)] ^ t[4 * 256 + (a >> 24)]
              ^ t[3 * 256 + (b & 0xFF)] ^ t[2 * 256 + ((b >> 8) & 0xFF)] ^ t[256 + ((b >> 16) & 0xFF)] ^ t[b >> 24];
            i += 8;
        }
        for (; i < data.Length; i++) c = t[(c ^ data[i]) & 0xFF] ^ (c >> 8);
        return ~c;
    }
}
