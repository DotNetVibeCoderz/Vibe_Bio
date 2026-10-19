namespace Biotech.Net.Genomics.Cram;

/// <summary>
/// A managed bzip2 decompressor (Burrows-Wheeler + MTF + Huffman), used for CRAM block method 2 and the arithmetic
/// coder's EXT mode. Block and stream CRCs are verified. Deprecated "randomised" blocks are not supported.
/// </summary>
internal static class Bzip2
{
    private static readonly uint[] CrcTable = BuildCrc();

    private static uint[] BuildCrc()
    {
        var t = new uint[256];
        for (uint i = 0; i < 256; i++)
        {
            uint c = i << 24;
            for (int k = 0; k < 8; k++) c = (c & 0x80000000) != 0 ? (c << 1) ^ 0x04C11DB7 : c << 1;
            t[i] = c;
        }
        return t;
    }

    private ref struct Bits(ReadOnlySpan<byte> data)
    {
        private readonly ReadOnlySpan<byte> _d = data;
        private int _pos;
        private ulong _buf;
        private int _n;

        public uint Read(int count)
        {
            while (_n < count)
            {
                if (_pos >= _d.Length) throw new InvalidDataException("bzip2: unexpected end of data.");
                _buf = (_buf << 8) | _d[_pos++];
                _n += 8;
            }
            _n -= count;
            return (uint)((_buf >> _n) & ((1UL << count) - 1));
        }

        public bool Bit() => Read(1) != 0;

        public bool AtEnd => _pos >= _d.Length && _n < 8;

        public void AlignToByte() => _n -= _n % 8;
    }

    public static byte[] Decompress(ReadOnlySpan<byte> data, int expectedLength = -1)
    {
        var output = new List<byte>(expectedLength > 0 ? expectedLength : 1 << 16);
        var bits = new Bits(data);
        while (!bits.AtEnd)
        {
            if (bits.Read(8) != 'B' || bits.Read(8) != 'Z' || bits.Read(8) != 'h') throw new InvalidDataException("bzip2: missing 'BZh' signature.");
            int level = (int)bits.Read(8) - '0';
            if (level is < 1 or > 9) throw new InvalidDataException("bzip2: invalid block size.");
            int maxBlock = level * 100000;
            uint combined = 0;
            while (true)
            {
                ulong magic = ((ulong)bits.Read(24) << 24) | bits.Read(24);
                uint crc = bits.Read(32);
                if (magic == 0x177245385090UL)
                {
                    if (crc != combined) throw new InvalidDataException("bzip2: stream CRC mismatch.");
                    bits.AlignToByte();
                    break;
                }
                if (magic != 0x314159265359UL) throw new InvalidDataException("bzip2: bad block signature.");
                int start = output.Count;
                DecodeBlock(ref bits, maxBlock, output);
                uint blockCrc = 0xFFFFFFFF;
                for (int i = start; i < output.Count; i++) blockCrc = (blockCrc << 8) ^ CrcTable[(blockCrc >> 24) ^ output[i]];
                blockCrc = ~blockCrc;
                if (blockCrc != crc) throw new InvalidDataException("bzip2: block CRC mismatch.");
                combined = ((combined << 1) | (combined >> 31)) ^ crc;
            }
            if (expectedLength >= 0 && output.Count >= expectedLength) break;
        }
        return output.ToArray();
    }

    private static void DecodeBlock(ref Bits bits, int maxBlock, List<byte> output)
    {
        if (bits.Bit()) throw new NotSupportedException("bzip2: randomised blocks are not supported.");
        int origPtr = (int)bits.Read(24);
        // symbol map
        var seqToUnseq = new byte[256];
        int nInUse = 0;
        uint used16 = bits.Read(16);
        for (int i = 0; i < 16; i++)
            if ((used16 & (0x8000u >> i)) != 0)
            {
                uint used = bits.Read(16);
                for (int j = 0; j < 16; j++)
                    if ((used & (0x8000u >> j)) != 0) seqToUnseq[nInUse++] = (byte)(i * 16 + j);
            }
        if (nInUse == 0) throw new InvalidDataException("bzip2: no symbols in use.");
        int alphaSize = nInUse + 2;
        int nGroups = (int)bits.Read(3);
        if (nGroups is < 2 or > 6) throw new InvalidDataException("bzip2: invalid number of Huffman groups.");
        int nSelectors = (int)bits.Read(15);
        if (nSelectors < 1) throw new InvalidDataException("bzip2: no selectors.");
        var mtfGroups = new byte[nGroups];
        for (int i = 0; i < nGroups; i++) mtfGroups[i] = (byte)i;
        var selectors = new byte[nSelectors];
        for (int i = 0; i < nSelectors; i++)
        {
            int j = 0;
            while (bits.Bit())
            {
                if (++j >= nGroups) throw new InvalidDataException("bzip2: invalid selector.");
            }
            byte v = mtfGroups[j];
            for (; j > 0; j--) mtfGroups[j] = mtfGroups[j - 1];
            mtfGroups[0] = v;
            selectors[i] = v;
        }
        var tables = new Huffman[nGroups];
        var lengths = new int[alphaSize];
        for (int g = 0; g < nGroups; g++)
        {
            int len = (int)bits.Read(5);
            for (int s = 0; s < alphaSize; s++)
            {
                while (true)
                {
                    if (len is < 1 or > 20) throw new InvalidDataException("bzip2: invalid code length.");
                    if (!bits.Bit()) break;
                    len += bits.Bit() ? -1 : 1;
                }
                lengths[s] = len;
            }
            tables[g] = new Huffman(lengths);
        }

        // Huffman + RUNA/RUNB + MTF decode into ll
        var mtf = new byte[256];
        for (int i = 0; i < 256; i++) mtf[i] = (byte)i;
        var counts = new int[256];
        var ll = new byte[maxBlock];
        int n = 0, eob = nInUse + 1, groupPos = 0, selIdx = 0;
        Huffman table = tables[selectors[0]];
        int run = 0, runWeight = 1;
        while (true)
        {
            if (groupPos == 50)
            {
                groupPos = 0;
                if (++selIdx >= nSelectors) throw new InvalidDataException("bzip2: selectors exhausted.");
                table = tables[selectors[selIdx]];
            }
            groupPos++;
            int sym = table.Decode(ref bits);
            if (sym <= 1)
            {
                run += (sym + 1) * runWeight;
                runWeight <<= 1;
                if (run > maxBlock) throw new InvalidDataException("bzip2: run exceeds the block size.");
                continue;
            }
            if (run > 0)
            {
                byte b = seqToUnseq[mtf[0]];
                if (n + run > maxBlock) throw new InvalidDataException("bzip2: block overflow.");
                ll.AsSpan(n, run).Fill(b);
                counts[b] += run;
                n += run;
                run = 0;
                runWeight = 1;
            }
            if (sym == eob) break;
            int idx = sym - 1;
            byte m = mtf[idx];
            Buffer.BlockCopy(mtf, 0, mtf, 1, idx);
            mtf[0] = m;
            byte ub = seqToUnseq[m];
            if (n >= maxBlock) throw new InvalidDataException("bzip2: block overflow.");
            ll[n++] = ub;
            counts[ub]++;
        }
        if (origPtr >= n && n > 0) throw new InvalidDataException("bzip2: invalid origin pointer.");

        // inverse BWT
        var tt = new uint[n];
        var cf = new int[256];
        for (int i = 0, s = 0; i < 256; i++) { cf[i] = s; s += counts[i]; }
        for (int i = 0; i < n; i++) tt[i] = ll[i];
        for (int i = 0; i < n; i++) tt[cf[ll[i]]++] |= (uint)i << 8;
        if (n == 0) return;
        uint pos = tt[origPtr] >> 8;
        // undo the initial run-length encoding (4 equal bytes then a count) while emitting
        int last = -1, same = 0;
        for (int k = 0; k < n; k++)
        {
            pos = tt[pos];
            byte b = (byte)pos;
            pos >>= 8;
            if (same == 4)
            {
                for (int r = 0; r < b; r++) output.Add((byte)last);
                same = 0;
                last = -1;
                continue;
            }
            if (b == last) same++;
            else { same = 1; last = b; }
            output.Add(b);
        }
    }

    private sealed class Huffman
    {
        private readonly int[] _count = new int[21];
        private readonly int[] _first = new int[21];
        private readonly int[] _offset = new int[21];
        private readonly int[] _symbols;

        public Huffman(int[] lengths)
        {
            foreach (int l in lengths) _count[l]++;
            _symbols = new int[lengths.Length];
            int code = 0, off = 0;
            for (int len = 1; len <= 20; len++)
            {
                _first[len] = code;
                _offset[len] = off;
                off += _count[len];
                code = (code + _count[len]) << 1;
            }
            var next = (int[])_offset.Clone();
            for (int s = 0; s < lengths.Length; s++) _symbols[next[lengths[s]]++] = s;
        }

        public int Decode(ref Bits bits)
        {
            int code = 0;
            for (int len = 1; len <= 20; len++)
            {
                code = (code << 1) | (int)bits.Read(1);
                int d = code - _first[len];
                if (d >= 0 && d < _count[len]) return _symbols[_offset[len] + d];
            }
            throw new InvalidDataException("bzip2: invalid Huffman code.");
        }
    }
}
