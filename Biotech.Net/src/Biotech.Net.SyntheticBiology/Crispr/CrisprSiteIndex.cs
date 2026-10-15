using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// A searchable index of every PAM-adjacent protospacer in a set of sequences, used for off-target discovery.
/// </summary>
/// <remarks>
/// <para>
/// Each candidate site is stored as a 2-bit packed word (two bits per base, most significant base first), so a
/// 20 nt protospacer occupies 40 bits of a <see cref="ulong"/> and comparing a guide with a site is one XOR,
/// one fold and one pop-count. The sites are additionally bucketed by <see cref="SeedBlocks"/> contiguous
/// blocks: by the pigeonhole principle, an alignment with at most <c>SeedBlocks − 1</c> mismatches must match
/// at least one block exactly, so a query only has to examine the sites sharing a block with it.
/// </para>
/// <para>
/// Memory is roughly 20 bytes per candidate site plus the hash buckets. A bacterial genome (about 400 000
/// NGG sites per megabase, both strands) indexes in well under a second; a mammalian genome is feasible but
/// wants a machine with several gigabytes free. Nothing here is persisted — build the index once per process
/// and reuse it for every guide.
/// </para>
/// </remarks>
public sealed class CrisprSiteIndex
{
    /// <summary>
    /// Number of seed blocks the protospacer is split into. Five blocks make searches with up to four
    /// mismatches exact; the index refuses queries that allow more.
    /// </summary>
    public const int SeedBlocks = 5;

    private readonly CasNuclease _nuclease;
    private readonly List<string> _ids = [];
    private readonly List<byte[]> _sequences = [];
    private readonly Dictionary<ulong, List<int>> _buckets = [];

    private ulong[] _packed = [];
    private int[] _start = [];
    private sbyte[] _strand = [];
    private int[] _sequenceIndex = [];
    private int _count;

    /// <summary>Creates an empty index for <paramref name="nuclease"/>.</summary>
    /// <exception cref="ArgumentException">The nuclease's protospacer is longer than 31 nt (the 2-bit word limit).</exception>
    public CrisprSiteIndex(CasNuclease nuclease)
    {
        ArgumentNullException.ThrowIfNull(nuclease);
        if (nuclease.GuideLength > 31)
            throw new ArgumentException("2-bit packing supports protospacers of at most 31 nucleotides.", nameof(nuclease));
        _nuclease = nuclease;
    }

    /// <summary>The nuclease whose PAM defined the indexed sites.</summary>
    public CasNuclease Nuclease => _nuclease;

    /// <summary>Number of candidate sites in the index (both strands).</summary>
    public int SiteCount => _count;

    /// <summary>Number of sequences indexed.</summary>
    public int SequenceCount => _ids.Count;

    /// <summary>The identifiers of the indexed sequences.</summary>
    public IReadOnlyList<string> SequenceIds => _ids;

    /// <summary>
    /// Builds an index over <paramref name="sequences"/>. Both the canonical PAM and, optionally, the relaxed
    /// NAG PAM are indexed for SpCas9.
    /// </summary>
    public static CrisprSiteIndex Build(IEnumerable<Sequence> sequences, CasNuclease? nuclease = null, bool includeNagPam = true)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var index = new CrisprSiteIndex(nuclease ?? CasNuclease.SpCas9);
        foreach (var sequence in sequences) index.Add(sequence, includeNagPam);
        return index;
    }

    /// <summary>Builds an index over a single sequence given as a string.</summary>
    public static CrisprSiteIndex Build(string sequence, string id = "target", CasNuclease? nuclease = null, bool includeNagPam = true)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var index = new CrisprSiteIndex(nuclease ?? CasNuclease.SpCas9);
        index.Add(id, DnaUtility.Ascii(DnaUtility.NormalizeDna(sequence)), includeNagPam);
        return index;
    }

    /// <summary>Adds every PAM-adjacent site of <paramref name="sequence"/> to the index.</summary>
    public void Add(Sequence sequence, bool includeNagPam = true)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        Add(sequence.Id, sequence.Span.ToArray(), includeNagPam);
    }

    /// <summary>Adds every PAM-adjacent site of an ASCII buffer; the buffer is retained for reporting.</summary>
    public void Add(string id, byte[] residues, bool includeNagPam = true)
    {
        ArgumentException.ThrowIfNullOrEmpty(id);
        ArgumentNullException.ThrowIfNull(residues);

        int sequenceIndex = _ids.Count;
        _ids.Add(id);
        _sequences.Add(residues);

        int guideLength = _nuclease.GuideLength;
        int pamLength = _nuclease.PamLength;

        var pams = new List<byte[]> { DnaUtility.Ascii(_nuclease.Pam) };
        if (includeNagPam && _nuclease.PamSide == PamSide.ThreePrime && _nuclease.Pam == "NGG")
            pams.Add(DnaUtility.Ascii("NAG"));

        var span = residues.AsSpan();
        Span<byte> pamRc = stackalloc byte[pamLength];
        foreach (var pam in pams)
        {
            for (int i = 0; i < pamLength; i++) pamRc[i] = DnaUtility.Complement(pam[pamLength - 1 - i]);

            for (int p = 0; p + pamLength <= span.Length; p++)
            {
                var window = span.Slice(p, pamLength);

                if (DnaUtility.MatchesIupac(window, pam))
                {
                    int protoStart = _nuclease.PamSide == PamSide.ThreePrime ? p - guideLength : p + pamLength;
                    AddSite(span, sequenceIndex, protoStart, Strand.Forward);
                }

                if (DnaUtility.MatchesIupac(window, pamRc))
                {
                    int protoStart = _nuclease.PamSide == PamSide.ThreePrime ? p + pamLength : p - guideLength;
                    AddSite(span, sequenceIndex, protoStart, Strand.Reverse);
                }
            }
        }
    }

    private void AddSite(ReadOnlySpan<byte> sequence, int sequenceIndex, int protoStart, Strand strand)
    {
        int guideLength = _nuclease.GuideLength;
        if (protoStart < 0 || protoStart + guideLength > sequence.Length) return;

        var protospacer = sequence.Slice(protoStart, guideLength);
        ulong packed;
        if (strand == Strand.Reverse)
        {
            // Pack the reverse complement directly, without materialising it.
            packed = 0;
            for (int i = guideLength - 1; i >= 0; i--)
            {
                int c = DnaUtility.Code(DnaUtility.Complement(protospacer[i]));
                if (c < 0) return;
                packed = (packed << 2) | (uint)c;
            }
        }
        else if (!DnaUtility.TryPack(protospacer, out packed))
        {
            return;
        }

        if (_count == _packed.Length) Grow();
        _packed[_count] = packed;
        _start[_count] = protoStart;
        _strand[_count] = (sbyte)strand;
        _sequenceIndex[_count] = sequenceIndex;

        for (int b = 0; b < SeedBlocks; b++)
        {
            ulong key = BlockKey(packed, b);
            if (!_buckets.TryGetValue(key, out var bucket)) _buckets[key] = bucket = [];
            bucket.Add(_count);
        }

        _count++;
    }

    private void Grow()
    {
        int capacity = _packed.Length == 0 ? 1024 : _packed.Length * 2;
        Array.Resize(ref _packed, capacity);
        Array.Resize(ref _start, capacity);
        Array.Resize(ref _strand, capacity);
        Array.Resize(ref _sequenceIndex, capacity);
    }

    private (int Offset, int Length) BlockRange(int block)
    {
        int guideLength = _nuclease.GuideLength;
        int baseSize = guideLength / SeedBlocks;
        int remainder = guideLength % SeedBlocks;
        int offset = block * baseSize + Math.Min(block, remainder);
        int length = baseSize + (block < remainder ? 1 : 0);
        return (offset, length);
    }

    private ulong BlockKey(ulong packed, int block)
    {
        var (offset, length) = BlockRange(block);
        int guideLength = _nuclease.GuideLength;
        // Base 0 is the most significant pair; shift it down to the low bits.
        int shift = 2 * (guideLength - offset - length);
        ulong mask = length >= 32 ? ulong.MaxValue : (1UL << (2 * length)) - 1;
        return ((ulong)block << 56) | ((packed >> shift) & mask);
    }

    /// <summary>
    /// Finds every indexed site within <paramref name="maxMismatches"/> mismatches of <paramref name="spacer"/>.
    /// </summary>
    /// <remarks>
    /// Only mismatches are considered here; see <see cref="OffTargetSearch"/> for the bulge-aware search. The
    /// returned sites carry per-hit CFD and MIT scores; <see cref="OffTargetSearch.Search(string, CrisprSiteIndex, OffTargetSearchOptions)"/>
    /// adds the bulge pass and the aggregate specificity scores.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is not in 0…<see cref="SeedBlocks"/> − 1.</exception>
    public IReadOnlyList<OffTarget> Query(string spacer, int maxMismatches = 4)
    {
        ArgumentException.ThrowIfNullOrEmpty(spacer);
        if (maxMismatches is < 0 or >= SeedBlocks)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches),
                $"The seed index is exact for up to {SeedBlocks - 1} mismatches.");

        string normalised = DnaUtility.NormalizeDna(spacer);
        if (normalised.Length != _nuclease.GuideLength)
            throw new ArgumentException($"{_nuclease.Name} uses {_nuclease.GuideLength} nt spacers, got {normalised.Length}.", nameof(spacer));
        if (!DnaUtility.TryPack(DnaUtility.Ascii(normalised), out ulong query))
            throw new ArgumentException("The spacer contains an ambiguity code.", nameof(spacer));

        var seen = new HashSet<int>();
        var results = new List<OffTarget>();
        int guideLength = _nuclease.GuideLength;

        for (int b = 0; b < SeedBlocks; b++)
        {
            if (!_buckets.TryGetValue(BlockKey(query, b), out var bucket)) continue;
            for (int i = 0; i < bucket.Count; i++)
            {
                int site = bucket[i];
                if (!seen.Add(site)) continue;
                int mismatches = DnaUtility.MismatchCount(query, _packed[site], guideLength);
                if (mismatches > maxMismatches) continue;
                results.Add(Materialise(site, normalised, mismatches));
            }
        }

        return results;
    }

    /// <summary>
    /// Scans every site linearly instead of using the seed buckets. Slower but unrestricted in the number of
    /// mismatches, and useful as a correctness oracle for the seeded path.
    /// </summary>
    public IReadOnlyList<OffTarget> QueryExhaustive(string spacer, int maxMismatches)
    {
        ArgumentException.ThrowIfNullOrEmpty(spacer);
        string normalised = DnaUtility.NormalizeDna(spacer);
        if (!DnaUtility.TryPack(DnaUtility.Ascii(normalised), out ulong query))
            throw new ArgumentException("The spacer contains an ambiguity code.", nameof(spacer));

        var results = new List<OffTarget>();
        int guideLength = _nuclease.GuideLength;
        for (int site = 0; site < _count; site++)
        {
            int mismatches = DnaUtility.MismatchCount(query, _packed[site], guideLength);
            if (mismatches > maxMismatches) continue;
            results.Add(Materialise(site, normalised, mismatches));
        }
        return results;
    }

    /// <summary>The protospacer of an indexed site, 5′→3′ on the targeted strand.</summary>
    public string ProtospacerAt(int site) => DnaUtility.Unpack(_packed[site], _nuclease.GuideLength);

    internal OffTarget Materialise(int site, string spacer, int mismatches)
    {
        int guideLength = _nuclease.GuideLength;
        var strand = (Strand)_strand[site];
        var sequence = _sequences[_sequenceIndex[site]].AsSpan();
        int start = _start[site];

        string protospacer = DnaUtility.Unpack(_packed[site], guideLength);
        string pam = PamAt(sequence, start, strand);

        double cfd = guideLength >= CfdScorer.TableLength && pam.Length >= 2
            ? CfdScorer.Score(spacer, protospacer, pam)
            : double.NaN;
        double mit = MitSpecificityScorer.HitScore(spacer, protospacer);

        return new OffTarget(_ids[_sequenceIndex[site]], start, strand, protospacer, pam,
            mismatches, BulgeKind.None, cfd, mit, IsOnTarget: false);
    }

    private string PamAt(ReadOnlySpan<byte> sequence, int protoStart, Strand strand)
    {
        int pamLength = _nuclease.PamLength;
        int pamStart = _nuclease.PamSide == PamSide.ThreePrime
            ? (strand == Strand.Reverse ? protoStart - pamLength : protoStart + _nuclease.GuideLength)
            : (strand == Strand.Reverse ? protoStart + _nuclease.GuideLength : protoStart - pamLength);
        if (pamStart < 0 || pamStart + pamLength > sequence.Length) return string.Empty;
        var window = sequence.Slice(pamStart, pamLength);
        return strand == Strand.Reverse
            ? DnaUtility.ReverseComplement(window)
            : Encoding.ASCII.GetString(window).ToUpperInvariant();
    }

    /// <summary>The ASCII buffer of an indexed sequence, for bulge-aware rescoring.</summary>
    internal ReadOnlySpan<byte> SequenceAt(int index) => _sequences[index];

    /// <summary>Index of the sequence a site belongs to.</summary>
    internal int SequenceIndexOf(int site) => _sequenceIndex[site];
}
