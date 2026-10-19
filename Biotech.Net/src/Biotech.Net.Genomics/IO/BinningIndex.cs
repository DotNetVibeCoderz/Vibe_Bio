using System.Buffers.Binary;
using System.IO.Compression;

namespace Biotech.Net.Genomics.IO;

/// <summary>A run of BGZF data <c>[Start, End)</c> that may hold records overlapping a bin.</summary>
public readonly record struct Chunk(VirtualOffset Start, VirtualOffset End)
{
    public override string ToString() => $"{Start}-{End}";
}

/// <summary>On-disk flavour of a binning index.</summary>
public enum BinningIndexFormat
{
    /// <summary>BAM index (<c>.bai</c>): fixed 14-bit / 5-level bins plus a 16 kb linear index; uncompressed.</summary>
    Bai,

    /// <summary>Coordinate-sorted index (<c>.csi</c>): configurable bins, per-bin linear offsets; BGZF-compressed.</summary>
    Csi,

    /// <summary>Tabix index (<c>.tbi</c>): BAI layout plus a text-format header; BGZF-compressed.</summary>
    Tbi,
}

/// <summary>Index data for one reference sequence.</summary>
public sealed class BinningIndexReference
{
    internal readonly Dictionary<uint, List<Chunk>> BinMap = [];
    internal readonly Dictionary<uint, VirtualOffset> LOffsets = [];

    /// <summary>The bins present, with their chunks (meta-data pseudo-bin excluded).</summary>
    public IReadOnlyDictionary<uint, List<Chunk>> Bins => BinMap;

    /// <summary>Linear index: the smallest record start in each 2^min_shift window (BAI/TBI; empty for CSI).</summary>
    public VirtualOffset[] Linear { get; internal set; } = [];

    /// <summary>Per-bin smallest offset (CSI <c>loff</c>; derived from the linear index for BAI/TBI).</summary>
    public IReadOnlyDictionary<uint, VirtualOffset> BinLinearOffsets => LOffsets;

    /// <summary>Virtual offsets spanning this reference's records, from the meta-data pseudo-bin, when present.</summary>
    public Chunk? Span { get; internal set; }

    /// <summary>Mapped read count from the pseudo-bin (0 when absent).</summary>
    public ulong MappedCount { get; internal set; }

    /// <summary>Placed-unmapped read count from the pseudo-bin (0 when absent).</summary>
    public ulong UnmappedCount { get; internal set; }

    /// <summary>Whether the index had a meta-data pseudo-bin for this reference.</summary>
    public bool HasMetadata { get; internal set; }
}

/// <summary>
/// The UCSC binning scheme with linear offsets shared by BAI, CSI and tabix indexes: reading, writing, region
/// queries (bins + linear index) and building, with the same algorithm as htslib (<c>hts_idx_push</c> /
/// <c>hts_idx_finish</c>), so an index built here has the same bins and chunks as <c>samtools index</c>.
/// </summary>
public sealed class BinningIndex
{
    /// <summary>The BAI/TBI shift: 2^14 = 16 kb bottom-level bins.</summary>
    public const int DefaultMinShift = 14;

    /// <summary>The BAI/TBI depth.</summary>
    public const int DefaultDepth = 5;

    internal BinningIndex(BinningIndexFormat format, int minShift, int depth, List<BinningIndexReference> references, ulong? unplaced, byte[] aux)
    {
        Format = format;
        MinShift = minShift;
        Depth = depth;
        ReferenceList = references;
        UnplacedCount = unplaced;
        Aux = aux;
    }

    internal List<BinningIndexReference> ReferenceList;

    /// <summary>Format the index was read from or built for.</summary>
    public BinningIndexFormat Format { get; internal set; }

    /// <summary>Bottom-level bin width as a power of two.</summary>
    public int MinShift { get; }

    /// <summary>Number of bin levels below the root.</summary>
    public int Depth { get; }

    /// <summary>Per-reference index data, in reference order.</summary>
    public IReadOnlyList<BinningIndexReference> References => ReferenceList;

    /// <summary>Number of unplaced unmapped records (the optional trailing <c>n_no_coor</c>), when stored.</summary>
    public ulong? UnplacedCount { get; internal set; }

    /// <summary>CSI auxiliary data (the tabix header for tabix-CSI), empty otherwise.</summary>
    public byte[] Aux { get; internal set; }

    /// <summary>Number of bins in this scheme (meta-data pseudo-bin excluded).</summary>
    public int BinCount => BinLimit(Depth);

    /// <summary>Largest position the scheme can index (exclusive).</summary>
    public long MaxPosition => 1L << (MinShift + 3 * Depth);

    internal static int BinLimit(int depth) => ((1 << (3 * depth + 3)) - 1) / 7;

    internal static uint MetaBin(int depth) => (uint)BinLimit(depth) + 1;

    /// <summary>The smallest bin containing <c>[beg, end)</c> (htslib <c>hts_reg2bin</c>).</summary>
    public static int RegionToBin(long beg, long end, int minShift = DefaultMinShift, int depth = DefaultDepth)
    {
        int s = minShift, t = ((1 << (3 * depth)) - 1) / 7;
        --end;
        for (int l = depth; l > 0; --l, s += 3, t -= 1 << (3 * l))
            if (beg >> s == end >> s) return (int)(t + (beg >> s));
        return 0;
    }

    /// <summary>Every bin that may hold features overlapping <c>[beg, end)</c>.</summary>
    public static List<int> RegionToBins(long beg, long end, int minShift = DefaultMinShift, int depth = DefaultDepth)
    {
        var bins = new List<int>();
        int s = minShift + 3 * depth;
        if (end >= 1L << s) end = 1L << s;
        if (beg < 0) beg = 0;
        if (beg >= end) return bins;
        --end;
        for (int l = 0, t = 0; l <= depth; s -= 3, t += 1 << (3 * l), ++l)
        {
            long b = t + (beg >> s), e = t + (end >> s);
            for (long i = b; i <= e; ++i) bins.Add((int)i);
        }
        return bins;
    }

    internal static int BinFirst(int level) => ((1 << (3 * level)) - 1) / 7;

    internal static uint BinParent(uint bin) => (bin - 1) >> 3;

    internal static int BinLevel(uint bin)
    {
        int l = 0;
        for (uint b = bin; b != 0; b = BinParent(b)) l++;
        return l;
    }

    internal static long BinBottom(uint bin, int depth)
    {
        int l = BinLevel(bin);
        return (long)(bin - (uint)BinFirst(l)) << (3 * (depth - l));
    }

    /// <summary>
    /// The merged, sorted chunks to scan for records of reference <paramref name="refId"/> overlapping
    /// <c>[beg, end)</c> (htslib <c>hts_itr_query</c>); empty when nothing can overlap.
    /// </summary>
    public List<Chunk> Query(int refId, long beg, long end)
    {
        var result = new List<Chunk>();
        if (refId < 0 || refId >= ReferenceList.Count) return result;
        var r = ReferenceList[refId];
        if (r.BinMap.Count == 0 || end <= beg) return result;
        if (beg < 0) beg = 0;
        long maxPos = MaxPosition;
        if (beg >= maxPos) return result;
        bool unmapped = !r.HasMetadata || r.UnmappedCount > 0;

        // min_off: smallest offset worth reading (from the nearest extant bin's loff, refined by the linear index)
        long relOff = beg >> MinShift;
        uint bin = (uint)(BinFirst(Depth) + relOff);
        bool found;
        do
        {
            found = r.BinMap.ContainsKey(bin);
            if (found) break;
            uint first = (BinParent(bin) << 3) + 1;
            if (bin > first) --bin;
            else bin = BinParent(bin);
        } while (bin != 0);
        if (bin == 0) found = r.BinMap.ContainsKey(0);
        ulong binLoff = found && r.LOffsets.TryGetValue(bin, out var lo) ? lo.Value : 0;
        ulong minOff = binLoff;
        if (r.Linear.Length > 0 && relOff < r.Linear.Length)
        {
            if (minOff < r.Linear[relOff].Value) minOff = r.Linear[relOff].Value;
            if (unmapped)
            {
                long t;
                for (t = relOff - 1; t >= 0; t--)
                    if (r.Linear[t].Value < minOff) { minOff = r.Linear[t].Value; break; }
                if (found && (minOff < binLoff || t < 0)) minOff = binLoff;
            }
        }
        else if (unmapped && found) minOff = binLoff;

        // max_off: start of the first extant bin to the right of end
        ulong maxOff = ulong.MaxValue;
        if (end <= maxPos)
        {
            long b = BinFirst(Depth) + ((end - 1) >> MinShift) + 1;
            if (b >= BinCount) b = 0;
            while (true)
            {
                while (b % 8 == 1) b = BinParent((uint)b);
                if (b == 0) { maxOff = ulong.MaxValue; break; }
                if (r.BinMap.TryGetValue((uint)b, out var l) && l.Count > 0) { maxOff = l[0].Start.Value; break; }
                b++;
            }
        }

        var off = new List<(ulong U, ulong V)>();
        foreach (int b in RegionToBins(beg, end, MinShift, Depth))
        {
            if (!r.BinMap.TryGetValue((uint)b, out var chunks)) continue;
            foreach (var c in chunks)
                if (c.End.Value > minOff && c.Start.Value < maxOff)
                    off.Add((Math.Max(minOff, c.Start.Value), Math.Min(maxOff, c.End.Value)));
        }
        if (off.Count == 0) return result;
        off.Sort((a, c) => a.U.CompareTo(c.U));
        // drop chunks contained in the previous one
        int n = 1;
        for (int i = 1; i < off.Count; i++)
            if (off[n - 1].V < off[i].V) off[n++] = off[i];
        off.RemoveRange(n, off.Count - n);
        // resolve overlaps between neighbours
        for (int i = 1; i < off.Count; i++)
            if (off[i - 1].V >= off[i].U) off[i - 1] = (off[i - 1].U, off[i].U);
        // merge chunks that meet in the same BGZF block
        n = 1;
        for (int i = 1; i < off.Count; i++)
        {
            if (off[n - 1].V >> 16 == off[i].U >> 16) off[n - 1] = (off[n - 1].U, off[i].V);
            else off[n++] = off[i];
        }
        for (int i = 0; i < n; i++) result.Add(new Chunk(new VirtualOffset(off[i].U), new VirtualOffset(off[i].V)));
        return result;
    }

    /// <summary>Offset after the last record placed on a reference (start of the unplaced unmapped records), or null.</summary>
    public VirtualOffset? UnplacedStart()
    {
        ulong off = ulong.MaxValue;
        foreach (var r in ReferenceList)
            if (r.Span is { } s && (off == ulong.MaxValue || off < s.End.Value)) off = s.End.Value;
        if (off == ulong.MaxValue) return UnplacedCount > 0 ? new VirtualOffset(0) : null;
        return new VirtualOffset(off);
    }

    // ------------------------------------------------------------------------------------------------ reading

    /// <summary>Loads a <c>.bai</c> or <c>.csi</c> file (detected from its magic bytes).</summary>
    public static BinningIndex Load(string path)
    {
        using var fs = File.OpenRead(path);
        return Read(fs);
    }

    /// <summary>Reads a BAI (raw) or CSI (BGZF) index from a stream, detected from the magic bytes.</summary>
    public static BinningIndex Read(Stream stream)
    {
        var data = ReadAll(stream);
        if (data.Length >= 2 && data[0] == 0x1f && data[1] == 0x8b) data = Decompress(data);
        if (data.Length < 4) throw new InvalidDataException("Index file is too short.");
        var r = new SpanReader(data);
        var magic = r.Bytes(4);
        if (magic.SequenceEqual("BAI\u0001"u8)) return ReadBaiBody(ref r, BinningIndexFormat.Bai, r.Int32());
        if (magic.SequenceEqual("CSI\u0001"u8))
        {
            int minShift = r.Int32(), depth = r.Int32();
            int lAux = r.Int32();
            var aux = r.Bytes(lAux).ToArray();
            int nRef = r.Int32();
            var refs = new List<BinningIndexReference>(nRef);
            uint meta = MetaBin(depth);
            for (int i = 0; i < nRef; i++)
            {
                var ri = new BinningIndexReference();
                int nBin = r.Int32();
                for (int b = 0; b < nBin; b++)
                {
                    uint bin = r.UInt32();
                    var loff = new VirtualOffset(r.UInt64());
                    int nChunk = r.Int32();
                    ReadChunks(ref r, ri, bin, nChunk, meta);
                    if (bin != meta) ri.LOffsets[bin] = loff;
                }
                refs.Add(ri);
            }
            ulong? unplaced = r.Remaining >= 8 ? r.UInt64() : null;
            return new BinningIndex(BinningIndexFormat.Csi, minShift, depth, refs, unplaced, aux);
        }
        throw new InvalidDataException("Not a BAI or CSI index (unknown magic bytes).");
    }

    /// <summary>Reads the per-reference body shared by BAI and TBI (after the header).</summary>
    internal static BinningIndex ReadBaiBody(ref SpanReader r, BinningIndexFormat format, int nRef)
    {
        var refs = new List<BinningIndexReference>(nRef);
        uint meta = MetaBin(DefaultDepth);
        for (int i = 0; i < nRef; i++)
        {
            var ri = new BinningIndexReference();
            int nBin = r.Int32();
            for (int b = 0; b < nBin; b++)
            {
                uint bin = r.UInt32();
                int nChunk = r.Int32();
                ReadChunks(ref r, ri, bin, nChunk, meta);
            }
            int nIntv = r.Int32();
            var lin = new VirtualOffset[nIntv];
            for (int k = 0; k < nIntv; k++) lin[k] = new VirtualOffset(r.UInt64());
            ri.Linear = lin;
            // htslib derives per-bin loff from the linear index on load (update_loff)
            foreach (var bin in ri.BinMap.Keys)
            {
                long bot = BinBottom(bin, DefaultDepth);
                ri.LOffsets[bin] = bot < lin.Length ? lin[bot] : new VirtualOffset(0);
            }
            refs.Add(ri);
        }
        ulong? unplaced = r.Remaining >= 8 ? r.UInt64() : null;
        return new BinningIndex(format, DefaultMinShift, DefaultDepth, refs, unplaced, []);
    }

    private static void ReadChunks(ref SpanReader r, BinningIndexReference ri, uint bin, int nChunk, uint meta)
    {
        if (bin == meta)
        {
            if (nChunk >= 1) ri.Span = new Chunk(new VirtualOffset(r.UInt64()), new VirtualOffset(r.UInt64()));
            if (nChunk >= 2) { ri.MappedCount = r.UInt64(); ri.UnmappedCount = r.UInt64(); }
            for (int c = 2; c < nChunk; c++) { r.UInt64(); r.UInt64(); }
            ri.HasMetadata = true;
            return;
        }
        var list = new List<Chunk>(nChunk);
        for (int c = 0; c < nChunk; c++) list.Add(new Chunk(new VirtualOffset(r.UInt64()), new VirtualOffset(r.UInt64())));
        ri.BinMap[bin] = list;
    }

    internal static byte[] ReadAll(Stream s)
    {
        using var ms = new MemoryStream();
        s.CopyTo(ms);
        return ms.ToArray();
    }

    internal static byte[] Decompress(byte[] gz)
    {
        using var ms = new MemoryStream();
        using (var bg = new BgzfReader(new MemoryStream(gz))) bg.CopyTo(ms);
        return ms.ToArray();
    }

    // ------------------------------------------------------------------------------------------------ writing

    /// <summary>Saves the index; the format follows <see cref="Format"/> (BAI raw, CSI BGZF-compressed).</summary>
    public void Save(string path)
    {
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None);
        Write(fs);
    }

    /// <summary>Writes the index as BAI (raw) or CSI (BGZF), per <see cref="Format"/>. TBI is written by <see cref="TabixIndex"/>.</summary>
    public void Write(Stream stream)
    {
        switch (Format)
        {
            case BinningIndexFormat.Bai:
            {
                var w = new BinaryWriter(stream, System.Text.Encoding.ASCII, leaveOpen: true);
                w.Write("BAI\u0001"u8);
                w.Write(ReferenceList.Count);
                WriteBaiBody(w);
                w.Flush();
                break;
            }
            case BinningIndexFormat.Csi:
            {
                using var bg = new BgzfWriter(stream, CompressionLevel.Optimal, leaveOpen: true);
                var w = new BinaryWriter(bg, System.Text.Encoding.ASCII, leaveOpen: true);
                w.Write("CSI\u0001"u8);
                w.Write(MinShift);
                w.Write(Depth);
                w.Write(Aux.Length);
                w.Write(Aux);
                w.Write(ReferenceList.Count);
                uint meta = MetaBin(Depth);
                foreach (var r in ReferenceList)
                {
                    int nBin = r.BinMap.Count + (r.HasMetadata ? 1 : 0);
                    w.Write(nBin);
                    foreach (var (bin, chunks) in r.BinMap.OrderBy(kv => kv.Key))
                    {
                        w.Write(bin);
                        w.Write(r.LOffsets.TryGetValue(bin, out var lo) ? lo.Value : 0UL);
                        w.Write(chunks.Count);
                        foreach (var c in chunks) { w.Write(c.Start.Value); w.Write(c.End.Value); }
                    }
                    if (r.HasMetadata) WriteMeta(w, r, meta, csi: true);
                }
                w.Write(UnplacedCount ?? 0UL);
                w.Flush();
                break;
            }
            default:
                throw new InvalidOperationException("TBI indexes are written through TabixIndex.");
        }
    }

    internal void WriteBaiBody(BinaryWriter w, bool skipEmptyReferences = false)
    {
        uint meta = MetaBin(DefaultDepth);
        foreach (var r in ReferenceList)
        {
            bool has = r.BinMap.Count > 0 || r.HasMetadata;
            if (skipEmptyReferences && !has) continue;
            w.Write(r.BinMap.Count + (r.HasMetadata ? 1 : 0));
            foreach (var (bin, chunks) in r.BinMap.OrderBy(kv => kv.Key))
            {
                w.Write(bin);
                w.Write(chunks.Count);
                foreach (var c in chunks) { w.Write(c.Start.Value); w.Write(c.End.Value); }
            }
            if (r.HasMetadata) WriteMeta(w, r, meta, csi: false);
            w.Write(r.Linear.Length);
            foreach (var v in r.Linear) w.Write(v.Value);
        }
        w.Write(UnplacedCount ?? 0UL);
    }

    private static void WriteMeta(BinaryWriter w, BinningIndexReference r, uint meta, bool csi)
    {
        w.Write(meta);
        if (csi) w.Write(0UL);
        w.Write(2);
        w.Write(r.Span?.Start.Value ?? 0UL);
        w.Write(r.Span?.End.Value ?? 0UL);
        w.Write(r.MappedCount);
        w.Write(r.UnmappedCount);
    }
}

/// <summary>
/// Builds a <see cref="BinningIndex"/> from records pushed in coordinate order, with htslib's algorithm
/// (<c>hts_idx_push</c>, <c>hts_idx_finish</c>, <c>update_loff</c>, <c>compress_binning</c>).
/// </summary>
public sealed class BinningIndexBuilder
{
    private const ulong Unset = ulong.MaxValue;
    private const ulong MinMarkerDistance = 0x10000;
    private readonly BinningIndexFormat _format;
    private readonly int _minShift, _depth;
    private readonly int _nBins;
    private readonly uint _metaBin;
    private readonly List<Dictionary<uint, List<(ulong U, ulong V)>>?> _bidx = [];
    private readonly List<List<ulong>> _lidx = [];
    private int _lastTid = -1, _saveTid = -1;
    private uint _lastBin = uint.MaxValue, _saveBin = uint.MaxValue;
    private ulong _saveOff, _lastOff, _offBeg, _offEnd;
    private long _lastCoor = 0xffffffffL;
    private ulong _nMapped, _nUnmapped, _nNoCoor;
    private bool _finished;

    /// <summary>Starts an index whose first record begins at <paramref name="firstRecord"/>.</summary>
    /// <param name="format">Target format (BAI and TBI need the default shift and depth).</param>
    /// <param name="firstRecord">Virtual offset where the records begin (after any header).</param>
    /// <param name="referenceCount">Number of references (the index grows as needed).</param>
    /// <param name="minShift">Bottom-level bin shift (14 for BAI/TBI).</param>
    /// <param name="depth">Levels (5 for BAI/TBI).</param>
    public BinningIndexBuilder(BinningIndexFormat format, VirtualOffset firstRecord, int referenceCount = 0,
        int minShift = BinningIndex.DefaultMinShift, int depth = BinningIndex.DefaultDepth)
    {
        if (format != BinningIndexFormat.Csi && (minShift != BinningIndex.DefaultMinShift || depth != BinningIndex.DefaultDepth))
            throw new ArgumentException("BAI and TBI indexes use min_shift 14 and depth 5.");
        _format = format;
        _minShift = minShift;
        _depth = depth;
        _nBins = BinningIndex.BinLimit(depth);
        _metaBin = BinningIndex.MetaBin(depth);
        _saveOff = _lastOff = _offBeg = _offEnd = firstRecord.Value;
        Grow(referenceCount);
    }

    private void Grow(int n)
    {
        while (_bidx.Count < n) { _bidx.Add(null); _lidx.Add([]); }
    }

    /// <summary>
    /// Adds a record of reference <paramref name="refId"/> (-1 for unplaced) spanning <c>[beg, end)</c> that ends at
    /// virtual offset <paramref name="endOffset"/>. Records must arrive in coordinate order.
    /// </summary>
    public void Push(int refId, long beg, long end, VirtualOffset endOffset, bool mapped)
    {
        if (_finished) throw new InvalidOperationException("The index is already finished.");
        if (refId < 0) { beg = -1; end = 0; }
        long maxPos = 1L << (_minShift + 3 * _depth);
        if (refId >= 0 && (beg > maxPos || end > maxPos))
            throw new InvalidDataException($"Region {beg}..{end} cannot be stored in a {_format} index with these parameters" +
                (_format == BinningIndexFormat.Csi ? "; use a larger min_shift or depth." : "; use a CSI index."));
        Grow(refId + 1);
        if (_lastTid != refId || (_lastTid >= 0 && refId < 0))
        {
            if (refId >= 0 && _nNoCoor > 0)
                throw new InvalidDataException("Unplaced records must come last, after every placed record.");
            if (refId >= 0 && _bidx[refId] is not null)
                throw new InvalidDataException("Records of one reference are not contiguous (input is not coordinate-sorted).");
            _lastTid = refId;
            _lastBin = uint.MaxValue;
        }
        else if (refId >= 0 && _lastCoor > beg)
            throw new InvalidDataException($"Unsorted positions on reference #{refId + 1}: {_lastCoor + 1} followed by {beg + 1}.");
        if (end < beg) throw new InvalidDataException($"Invalid record on reference #{refId + 1}: end {end} < begin {beg + 1}.");

        if (refId >= 0)
        {
            _bidx[refId] ??= [];
            if (beg < 0) beg = 0;
            if (end <= 0) end = 1;
            InsertLinear(_lidx[refId], beg, end, _lastOff);
        }
        else _nNoCoor++;

        uint bin = (uint)BinningIndex.RegionToBin(beg, end, _minShift, _depth);
        if (_lastBin != bin)
        {
            if (_saveBin != uint.MaxValue) InsertBin(_saveTid, _saveBin, _saveOff, _lastOff);
            if (_lastBin == uint.MaxValue && _saveBin != uint.MaxValue)
            {
                _offEnd = _lastOff;
                InsertBin(_saveTid, _metaBin, _offBeg, _offEnd);
                InsertBin(_saveTid, _metaBin, _nMapped, _nUnmapped);
                _nMapped = _nUnmapped = 0;
                _offBeg = _offEnd;
            }
            _saveOff = _lastOff;
            _saveBin = _lastBin = bin;
            _saveTid = refId;
        }
        if (mapped) _nMapped++; else _nUnmapped++;
        _lastOff = endOffset.Value;
        _lastCoor = beg;
    }

    private void InsertLinear(List<ulong> l, long beg, long end, ulong offset)
    {
        int b = (int)(beg >> _minShift), e = (int)((end - 1) >> _minShift);
        while (l.Count < e + 1) l.Add(Unset);
        for (int i = b; i <= e; i++)
            if (l[i] == Unset) l[i] = offset;
    }

    private void InsertBin(int tid, uint bin, ulong beg, ulong end)
    {
        if (tid < 0) return; // only reachable for unplaced records, which htslib never stores in bins
        var map = _bidx[tid]!;
        if (!map.TryGetValue(bin, out var list)) map[bin] = list = [];
        list.Add((beg, end));
    }

    /// <summary>Completes the index; <paramref name="finalOffset"/> is the virtual offset after the last record.</summary>
    public BinningIndex Finish(VirtualOffset finalOffset)
    {
        if (_finished) throw new InvalidOperationException("The index is already finished.");
        _finished = true;
        if (_saveTid >= 0)
        {
            InsertBin(_saveTid, _saveBin, _saveOff, finalOffset.Value);
            InsertBin(_saveTid, _metaBin, _offBeg, finalOffset.Value);
            InsertBin(_saveTid, _metaBin, _nMapped, _nUnmapped);
        }
        var refs = new List<BinningIndexReference>(_bidx.Count);
        for (int i = 0; i < _bidx.Count; i++)
        {
            var lidx = _lidx[i];
            for (int l = lidx.Count - 2; l >= 0; l--)
                if (lidx[l] == Unset) lidx[l] = lidx[l + 1];
            var r = new BinningIndexReference();
            var map = _bidx[i];
            if (map is not null)
            {
                Compress(map);
                foreach (var (bin, list) in map)
                {
                    if (bin == _metaBin)
                    {
                        r.HasMetadata = true;
                        r.Span = new Chunk(new VirtualOffset(list[0].U), new VirtualOffset(list[0].V));
                        if (list.Count > 1) { r.MappedCount = list[1].U; r.UnmappedCount = list[1].V; }
                        continue;
                    }
                    r.BinMap[bin] = list.Select(c => new Chunk(new VirtualOffset(c.U), new VirtualOffset(c.V))).ToList();
                    long bot = BinningIndex.BinBottom(bin, _depth);
                    r.LOffsets[bin] = new VirtualOffset(bot < lidx.Count ? lidx[(int)bot] : 0);
                }
            }
            if (_format != BinningIndexFormat.Csi) r.Linear = lidx.Select(v => new VirtualOffset(v)).ToArray();
            refs.Add(r);
        }
        return new BinningIndex(_format, _minShift, _depth, refs, _nNoCoor, []);
    }

    private void Compress(Dictionary<uint, List<(ulong U, ulong V)>> map)
    {
        // merge a bin into its parent when its chunks span less than one BGZF block's worth of file
        for (int l = _depth; l > 0; --l)
        {
            uint start = (uint)BinningIndex.BinFirst(l), stop = (uint)BinningIndex.BinFirst(l + 1);
            foreach (var bin in map.Keys.Where(b => b >= start && b < stop && b < _nBins).ToList())
            {
                var p = map[bin];
                if (l < _depth && p.Count > 1) p.Sort((a, b) => a.U.CompareTo(b.U));
                if ((p[^1].V >> 16) - (p[0].U >> 16) < MinMarkerDistance)
                {
                    if (!map.TryGetValue(BinningIndex.BinParent(bin), out var q)) continue;
                    q.AddRange(p);
                    map.Remove(bin);
                }
            }
        }
        if (map.TryGetValue(0, out var root)) root.Sort((a, b) => a.U.CompareTo(b.U));
        // merge neighbouring chunks that start in the same BGZF block
        foreach (var (bin, p) in map)
        {
            if (bin >= _nBins) continue;
            int m = 0;
            for (int l = 1; l < p.Count; l++)
            {
                if (p[m].V >> 16 >= p[l].U >> 16)
                {
                    if (p[m].V < p[l].V) p[m] = (p[m].U, p[l].V);
                }
                else p[++m] = p[l];
            }
            p.RemoveRange(m + 1, p.Count - m - 1);
        }
    }
}

/// <summary>Little-endian reader over a byte span.</summary>
internal ref struct SpanReader(ReadOnlySpan<byte> data)
{
    private readonly ReadOnlySpan<byte> _data = data;
    public int Position;

    public readonly int Remaining => _data.Length - Position;

    public ReadOnlySpan<byte> Bytes(int n)
    {
        if (n < 0 || Position + n > _data.Length) throw new InvalidDataException("Unexpected end of data.");
        var s = _data.Slice(Position, n);
        Position += n;
        return s;
    }

    public int Int32() => BinaryPrimitives.ReadInt32LittleEndian(Bytes(4));
    public uint UInt32() => BinaryPrimitives.ReadUInt32LittleEndian(Bytes(4));
    public ulong UInt64() => BinaryPrimitives.ReadUInt64LittleEndian(Bytes(8));
}
