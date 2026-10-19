using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Sequences;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>Options for <see cref="ReadMapper"/>. Scores follow BWA-MEM defaults (match 1, mismatch 4, gap open 6, extend 1).</summary>
public sealed record MapperOptions
{
    /// <summary>Minimizer k-mer length (odd, ≤ 31; default 15).</summary>
    public int K { get; init; } = 15;

    /// <summary>Minimizer window: one minimizer per <see cref="W"/> consecutive k-mers (default 10).</summary>
    public int W { get; init; } = 10;

    /// <summary>Minimizers occurring more often than this in the reference are ignored as seeds (default 200).</summary>
    public int MaxOccurrences { get; init; } = 200;

    /// <summary>Half-width of the band around the seeded diagonal explored by the gapped alignment (default 12).</summary>
    public int Band { get; init; } = 12;

    /// <summary>Accept the ungapped alignment without dynamic programming when it has at most this many mismatches (default 3).</summary>
    public int MaxUngappedMismatches { get; init; } = 3;

    /// <summary>Candidate loci (seed clusters) verified per read (default 5).</summary>
    public int MaxCandidates { get; init; } = 5;

    /// <summary>Match score.</summary>
    public int Match { get; init; } = 1;

    /// <summary>Mismatch penalty (positive).</summary>
    public int Mismatch { get; init; } = 4;

    /// <summary>Gap-open penalty (positive); a gap of length n costs open + n·extend.</summary>
    public int GapOpen { get; init; } = 6;

    /// <summary>Gap-extension penalty (positive).</summary>
    public int GapExtend { get; init; } = 1;

    /// <summary>Minimum alignment score as a fraction of the perfect score (default 0.5).</summary>
    public double MinScoreFraction { get; init; } = 0.5;
}

/// <summary>The alignment of one read (or the record of an unmapped read).</summary>
public sealed class ReadAlignment
{
    internal ReadAlignment(string readId, byte[] forwardBases, byte[]? forwardQualities)
    {
        ReadId = readId;
        Bases = forwardBases;
        Qualities = forwardQualities;
    }

    /// <summary>Read identifier.</summary>
    public string ReadId { get; }

    /// <summary>Read bases as aligned to the forward reference strand (reverse-complemented for reverse hits).</summary>
    public byte[] Bases { get; }

    /// <summary>Phred qualities in the same orientation as <see cref="Bases"/>, or null.</summary>
    public byte[]? Qualities { get; }

    /// <summary>Whether the read was placed.</summary>
    public bool IsMapped { get; internal set; }

    /// <summary>Index of the reference sequence.</summary>
    public int ReferenceIndex { get; internal set; } = -1;

    /// <summary>Reference name (<c>*</c> when unmapped).</summary>
    public string ReferenceName { get; internal set; } = "*";

    /// <summary>0-based leftmost reference position.</summary>
    public int Position { get; internal set; } = -1;

    /// <summary>Whether the read aligned to the reverse strand.</summary>
    public bool IsReverse { get; internal set; }

    /// <summary>CIGAR operations.</summary>
    public CigarElement[] Cigar { get; internal set; } = [];

    /// <summary>Alignment score.</summary>
    public int Score { get; internal set; }

    /// <summary>Best score at a different locus (0 when none).</summary>
    public int SecondBestScore { get; internal set; }

    /// <summary>Mapping quality (0–60).</summary>
    public int MapQ { get; internal set; }

    /// <summary>Edit distance to the reference (NM).</summary>
    public int EditDistance { get; internal set; }

    /// <summary>MD tag string.</summary>
    public string Md { get; internal set; } = "";

    /// <summary>0-based half-open reference interval covered.</summary>
    public GenomicInterval Interval => new(ReferenceName, Position, Position + Mapping.Cigar.ReferenceLength(Cigar), IsReverse ? Strand.Reverse : Strand.Forward);

    /// <summary>The SAM record (SEQ/QUAL in reference orientation, NM/MD/AS/XS tags).</summary>
    public SamRecord ToSam()
    {
        string seq = Encoding.ASCII.GetString(Bases);
        string qual = Qualities is null ? "*" : string.Create(Qualities.Length, Qualities, static (s, q) => { for (int i = 0; i < s.Length; i++) s[i] = (char)(q[i] + 33); });
        if (!IsMapped)
            return new SamRecord(ReadId, SamFlags.Unmapped, "*", 0, 0, "*", "*", 0, 0, seq, qual);
        return new SamRecord(ReadId, IsReverse ? SamFlags.Reverse : SamFlags.None, ReferenceName, Position + 1, MapQ,
            Mapping.Cigar.Format(Cigar), "*", 0, 0, seq, qual)
        {
            Tags = [$"NM:i:{EditDistance}", $"MD:Z:{Md}", $"AS:i:{Score}", $"XS:i:{SecondBestScore}"],
        };
    }

    /// <summary>The alignment as an <see cref="AlignmentRecord"/> (reference ids resolved through <paramref name="header"/> when given).</summary>
    public AlignmentRecord ToAlignmentRecord(SamHeader? header = null) => AlignmentRecord.FromSamRecord(ToSam(), header);

    public override string ToString() => IsMapped
        ? $"{ReadId} → {ReferenceName}:{Position + 1} {(IsReverse ? '-' : '+')} {Mapping.Cigar.Format(Cigar)} MAPQ {MapQ}"
        : $"{ReadId} unmapped";
}

/// <summary>
/// Short-read mapper: a minimizer index (canonical (k, w)-minimizers, as in minimap) of the reference; reads are
/// seeded by shared minimizers, seeds are clustered by diagonal, and each candidate locus is verified first by an
/// ungapped comparison and, when that has too many mismatches, by a banded affine-gap semi-global alignment
/// (whole read, free reference ends). MAPQ is derived from the gap between the best and second-best locus scores.
/// </summary>
public sealed class ReadMapper
{
    private readonly Sequence[] _reference;
    private readonly byte[][] _refUpper;
    private readonly Dictionary<ulong, (int Offset, int Count)> _index;
    private readonly long[] _hits; // refIndex << 33 | position << 1 | strand

    /// <summary>Indexes <paramref name="reference"/>.</summary>
    public ReadMapper(IReadOnlyList<Sequence> reference, MapperOptions? options = null)
    {
        Options = options ?? new MapperOptions();
        if (Options.K is < 5 or > 31 || Options.K % 2 == 0) throw new ArgumentOutOfRangeException(nameof(options), "K must be odd, 5–31.");
        _reference = reference.ToArray();
        _refUpper = _reference.Select(s => { var b = s.Span.ToArray(); System.Text.Ascii.ToUpperInPlace(b, out _); return b; }).ToArray();
        var entries = new List<(ulong Hash, long Value)>();
        var buf = new List<(ulong Hash, int Pos, int Strand)>();
        for (int r = 0; r < _reference.Length; r++)
        {
            buf.Clear();
            Minimizers(_refUpper[r], Options.K, Options.W, buf);
            foreach (var (h, p, st) in buf) entries.Add((h, ((long)r << 33) | ((long)p << 1) | (uint)st));
        }
        entries.Sort((a, b) => a.Hash.CompareTo(b.Hash));
        _hits = new long[entries.Count];
        _index = new Dictionary<ulong, (int, int)>();
        for (int i = 0; i < entries.Count;)
        {
            int j = i;
            while (j < entries.Count && entries[j].Hash == entries[i].Hash) { _hits[j] = entries[j].Value; j++; }
            _index[entries[i].Hash] = (i, j - i);
            i = j;
        }
    }

    /// <summary>The options in use.</summary>
    public MapperOptions Options { get; }

    /// <summary>The reference sequences.</summary>
    public IReadOnlyList<Sequence> Reference => _reference;

    /// <summary>Distinct minimizers in the index.</summary>
    public int IndexSize => _index.Count;

    private static ulong Hash64(ulong key, ulong mask)
    {
        // Thomas Wang's invertible integer hash (as used by minimap) spreads poly-A and other low-complexity k-mers.
        key = (~key + (key << 21)) & mask;
        key ^= key >> 24;
        key = (key + (key << 3) + (key << 8)) & mask;
        key ^= key >> 14;
        key = (key + (key << 2) + (key << 4)) & mask;
        key ^= key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
    }

    /// <summary>Canonical (k, w)-minimizers of a sequence: (hash, position, strand) with strand 1 when the reverse complement is smaller.</summary>
    internal static void Minimizers(ReadOnlySpan<byte> s, int k, int w, List<(ulong Hash, int Pos, int Strand)> output)
    {
        ulong mask = (1UL << (2 * k)) - 1;
        int shift = 2 * (k - 1);
        ulong fwd = 0, rev = 0;
        int valid = 0;
        var window = new (ulong Hash, int Pos, int Strand)[w];
        int filled = 0;
        int lastPos = -1;
        for (int i = 0; i < s.Length; i++)
        {
            int c = Kmer.Code(s[i]);
            if (c < 0) { valid = 0; fwd = rev = 0; filled = 0; continue; }
            fwd = ((fwd << 2) | (uint)c) & mask;
            rev = (rev >> 2) | ((ulong)(3 - c) << shift);
            if (++valid < k) continue;
            if (fwd == rev) continue; // palindrome (even k only)
            int strand = fwd < rev ? 0 : 1;
            ulong h = Hash64(strand == 0 ? fwd : rev, mask);
            window[filled % w] = (h, i - k + 1, strand);
            filled++;
            if (filled < w) continue;
            var best = window[0];
            for (int j = 1; j < w; j++)
                if (window[j].Hash < best.Hash || (window[j].Hash == best.Hash && window[j].Pos < best.Pos)) best = window[j];
            if (best.Pos != lastPos) { output.Add(best); lastPos = best.Pos; }
        }
        // Short sequences (fewer than w k-mers) still get their smallest k-mer.
        if (filled > 0 && filled < w)
        {
            var best = window[0];
            for (int j = 1; j < filled; j++) if (window[j].Hash < best.Hash) best = window[j];
            output.Add(best);
        }
    }

    /// <summary>Maps one read.</summary>
    public ReadAlignment Map(SequenceRead read) => Map(read.Id, read.Sequence.Span, read.Quality);

    /// <summary>Maps one read without qualities.</summary>
    public ReadAlignment Map(Sequence read) => Map(read.Id, read.Span, default);

    /// <summary>Maps reads in parallel, preserving input order.</summary>
    public IEnumerable<ReadAlignment> Map(IEnumerable<SequenceRead> reads)
        => reads.AsParallel().AsOrdered().Select(Map);

    /// <summary>Maps reads (no qualities) in parallel, preserving input order.</summary>
    public IEnumerable<ReadAlignment> Map(IEnumerable<Sequence> reads)
        => reads.AsParallel().AsOrdered().Select(Map);

    private ReadAlignment Map(string id, ReadOnlySpan<byte> bases, ReadOnlySpan<byte> quals)
    {
        var fwd = bases.ToArray();
        System.Text.Ascii.ToUpperInPlace(fwd, out _);
        int len = fwd.Length;
        var unmapped = new ReadAlignment(id, fwd, quals.IsEmpty ? null : quals.ToArray());
        int k = Options.K;
        if (len < k) return unmapped;

        var mins = new List<(ulong Hash, int Pos, int Strand)>();
        Minimizers(fwd, k, Options.W, mins);
        var seeds = new List<(int Ref, bool Rev, long Diag)>();
        foreach (var (h, qpos, qstrand) in mins)
        {
            if (!_index.TryGetValue(h, out var slot) || slot.Count > Options.MaxOccurrences) continue;
            for (int i = slot.Offset; i < slot.Offset + slot.Count; i++)
            {
                long v = _hits[i];
                int r = (int)(v >> 33), rpos = (int)((v >> 1) & 0xFFFFFFFF), rstrand = (int)(v & 1);
                bool rev = (qstrand ^ rstrand) == 1;
                long diag = rev ? rpos - (len - qpos - k) : rpos - qpos;
                seeds.Add((r, rev, diag));
            }
        }
        if (seeds.Count == 0) return unmapped;
        seeds.Sort();

        // Cluster seeds on the same reference/strand whose diagonals are within the band.
        var clusters = new List<(int Ref, bool Rev, long MinDiag, long MaxDiag, int Votes)>();
        int band = Options.Band;
        foreach (var s in seeds)
        {
            if (clusters.Count > 0)
            {
                var c = clusters[^1];
                if (c.Ref == s.Ref && c.Rev == s.Rev && s.Diag - c.MaxDiag <= band)
                {
                    clusters[^1] = c with { MaxDiag = s.Diag, Votes = c.Votes + 1 };
                    continue;
                }
            }
            clusters.Add((s.Ref, s.Rev, s.Diag, s.Diag, 1));
        }

        byte[]? rc = null;
        int perfect = Options.Match * len;
        int minScore = (int)Math.Ceiling(Options.MinScoreFraction * perfect);
        Candidate? best = null, second = null;
        foreach (var c in clusters.OrderByDescending(c => c.Votes).Take(Options.MaxCandidates))
        {
            byte[] q = c.Rev ? rc ??= ReverseComplement(fwd) : fwd;
            var cand = Verify(c.Ref, c.Rev, q, c.MinDiag, c.MaxDiag);
            if (cand is null) continue;
            if (best is null || cand.Score > best.Score)
            {
                if (best is not null && !SameLocus(best, cand)) second = best;
                best = cand;
            }
            else if (!SameLocus(best, cand) && (second is null || cand.Score > second.Score)) second = cand;
        }
        if (best is null || best.Score < minScore) return unmapped;

        var bq = best.Rev ? rc! : fwd;
        byte[]? q2 = quals.IsEmpty ? null : quals.ToArray();
        if (q2 is not null && best.Rev) Array.Reverse(q2);
        var aln = new ReadAlignment(id, bq, q2)
        {
            IsMapped = true,
            ReferenceIndex = best.Ref,
            ReferenceName = _reference[best.Ref].Id,
            Position = best.Pos,
            IsReverse = best.Rev,
            Cigar = best.Cigar,
            Score = best.Score,
            SecondBestScore = second?.Score ?? 0,
        };
        int s2 = Math.Max(second?.Score ?? 0, 0);
        aln.MapQ = s2 >= best.Score ? 0 : Math.Min(60, (int)Math.Round(250.0 * (best.Score - s2) / best.Score));
        (aln.EditDistance, aln.Md) = EditInfo(_refUpper[best.Ref], best.Pos, bq, best.Cigar);
        return aln;
    }

    private bool SameLocus(Candidate a, Candidate b) => a.Ref == b.Ref && a.Rev == b.Rev && Math.Abs(a.Pos - b.Pos) <= Options.Band;

    private sealed record Candidate(int Ref, bool Rev, int Pos, CigarElement[] Cigar, int Score);

    private static byte[] ReverseComplement(byte[] s)
    {
        var r = new byte[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = s[i] switch { (byte)'A' => (byte)'T', (byte)'C' => (byte)'G', (byte)'G' => (byte)'C', (byte)'T' => (byte)'A', _ => (byte)'N' };
        return r;
    }

    private Candidate? Verify(int refIdx, bool rev, byte[] q, long minDiag, long maxDiag)
    {
        var refSeq = _refUpper[refIdx];
        int len = q.Length;
        // Ungapped check on the most supported diagonals (first and last of the cluster).
        Candidate? ungapped = null;
        foreach (long d in minDiag == maxDiag ? new[] { minDiag } : new[] { minDiag, maxDiag })
        {
            if (d < 0 || d + len > refSeq.Length) continue;
            int mm = 0;
            for (int i = 0; i < len && mm <= Options.MaxUngappedMismatches; i++)
                if (q[i] != refSeq[d + i] || q[i] == (byte)'N') mm++;
            if (mm <= Options.MaxUngappedMismatches)
            {
                int score = (len - mm) * Options.Match - mm * Options.Mismatch;
                if (ungapped is null || score > ungapped.Score)
                    ungapped = new Candidate(refIdx, rev, (int)d, [new CigarElement(len, 'M')], score);
            }
        }
        if (ungapped is not null) return ungapped;
        int band = Options.Band;
        long ws = Math.Max(0, minDiag - band), we = Math.Min(refSeq.Length, maxDiag + len + band);
        if (we - ws < len / 2) return null;
        return Banded(refIdx, rev, q, refSeq, (int)ws, (int)we);
    }

    private const int NegInf = int.MinValue / 4;

    /// <summary>Affine-gap semi-global alignment (Gotoh): the whole read against any substring of ref[ws..we).</summary>
    private Candidate? Banded(int refIdx, bool rev, byte[] q, byte[] refSeq, int ws, int we)
    {
        int n = q.Length, m = we - ws;
        int match = Options.Match, mis = Options.Mismatch, open = Options.GapOpen, ext = Options.GapExtend;
        int cols = m + 1;
        var H = new int[(n + 1) * cols];
        var E = new int[(n + 1) * cols];
        var F = new int[(n + 1) * cols];
        var tb = new byte[(n + 1) * cols]; // bits 0-1: H source (0 diag, 1 E, 2 F); bit 2: E extended; bit 3: F extended
        for (int j = 0; j <= m; j++) { H[j] = 0; E[j] = NegInf; F[j] = NegInf; }
        for (int i = 1; i <= n; i++)
        {
            int row = i * cols, prev = (i - 1) * cols;
            F[row] = -open - ext * i;
            H[row] = F[row];
            E[row] = NegInf;
            tb[row] = (byte)(2 | (i > 1 ? 8 : 0));
            byte qi = q[i - 1];
            for (int j = 1; j <= m; j++)
            {
                byte t = 0;
                int eExt = E[row + j - 1] - ext, eOpen = H[row + j - 1] - open - ext;
                int e = eExt > eOpen ? eExt : eOpen;
                if (eExt > eOpen) t |= 4;
                int fExt = F[prev + j] - ext, fOpen = H[prev + j] - open - ext;
                int f = fExt > fOpen ? fExt : fOpen;
                if (fExt > fOpen) t |= 8;
                int d = H[prev + j - 1] + (qi == refSeq[ws + j - 1] && qi != (byte)'N' ? match : -mis);
                int h = d;
                if (e > h) { h = e; t |= 1; }
                if (f > h) { h = f; t = (byte)((t & ~3) | 2); }
                H[row + j] = h; E[row + j] = e; F[row + j] = f; tb[row + j] = t;
            }
        }
        int bestJ = 0, bestScore = NegInf, last = n * cols;
        for (int j = 0; j <= m; j++) if (H[last + j] > bestScore) { bestScore = H[last + j]; bestJ = j; }

        var ops = new List<char>(n + 16);
        int ci = n, cj = bestJ, state = 0;
        while (ci > 0)
        {
            byte t = tb[ci * cols + cj];
            if (state == 0)
            {
                int src = t & 3;
                if (src == 0 && cj > 0) { ops.Add('M'); ci--; cj--; }
                else if (src == 1) state = 1;
                else state = 2;
            }
            else if (state == 1)
            {
                ops.Add('D');
                bool extended = (t & 4) != 0;
                cj--;
                if (!extended) state = 0;
            }
            else
            {
                ops.Add('I');
                bool extended = (t & 8) != 0;
                ci--;
                if (!extended) state = 0;
            }
            if (cj < 0) return null;
        }
        ops.Reverse();
        var cigar = Cigar.Compress(ops);
        // Leading/trailing insertions become soft clips; leading deletions cannot occur (free reference start).
        if (cigar.Length > 0 && cigar[0].Op == 'I') cigar[0] = cigar[0] with { Op = 'S' };
        if (cigar.Length > 0 && cigar[^1].Op == 'I') cigar[^1] = cigar[^1] with { Op = 'S' };
        return new Candidate(refIdx, rev, ws + cj, cigar, bestScore);
    }

    /// <summary>Edit distance and MD string of an alignment.</summary>
    private static (int Nm, string Md) EditInfo(byte[] refSeq, int pos, byte[] q, CigarElement[] cigar)
    {
        var md = new StringBuilder();
        int nm = 0, run = 0, r = pos, qi = 0;
        foreach (var e in cigar)
        {
            switch (e.Op)
            {
                case 'M':
                case '=':
                case 'X':
                    for (int k = 0; k < e.Length; k++, r++, qi++)
                    {
                        if (q[qi] == refSeq[r]) { run++; continue; }
                        nm++;
                        md.Append(run).Append((char)refSeq[r]);
                        run = 0;
                    }
                    break;
                case 'I':
                    nm += e.Length; qi += e.Length; break;
                case 'S':
                    qi += e.Length; break;
                case 'D':
                    nm += e.Length;
                    md.Append(run).Append('^');
                    for (int k = 0; k < e.Length; k++) md.Append((char)refSeq[r++]);
                    run = 0;
                    break;
            }
        }
        md.Append(run);
        return (nm, md.ToString());
    }
}
