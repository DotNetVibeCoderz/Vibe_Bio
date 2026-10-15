using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Cloning;

/// <summary>
/// A double-stranded DNA fragment with optional single-stranded ends, as produced by a restriction digest or
/// by PCR with tailed primers.
/// </summary>
/// <remarks>
/// <para>
/// <b>Representation.</b> <see cref="Sequence"/> is the full top strand written 5′→3′ <i>including</i> the
/// bases that are single-stranded at either end. <see cref="FivePrimeOverhang"/> is the prefix of
/// <see cref="Sequence"/> that is unpaired, and <see cref="ThreePrimeOverhang"/> is the suffix that is
/// unpaired. Keeping the overhangs as substrings of the top strand (rather than as detached oligos) means a
/// ligation is a plain string concatenation once the ends have been checked, and it keeps a 5′ and a 3′
/// overhang in the same coordinate system.
/// </para>
/// <para>
/// A sticky-ended EcoRI fragment therefore starts with <c>AATT</c> as its <see cref="FivePrimeOverhang"/> and
/// ends with <c>AATT</c> as its <see cref="ThreePrimeOverhang"/> — those four bases are present on the top
/// strand at the 5′ end and on the bottom strand at the 3′ end respectively.
/// </para>
/// <para>
/// Both overhang strings are top-strand text of the window they occupy, even for a 3′ overhang, whose
/// single-stranded bases physically sit on the bottom strand. Compatibility is therefore plain string equality
/// (<see cref="CanLigateTo"/>): the upstream fragment's 3′ window and the downstream fragment's 5′ window must
/// read the same, which is exactly the complementarity T4 DNA ligase needs to seal both nicks. The shared
/// window is written once in each neighbour, and <see cref="Ligate"/> drops the duplicate.
/// </para>
/// </remarks>
/// <param name="Id">Identifier used in assembly reports and in the features of the product.</param>
/// <param name="Sequence">Top strand 5′→3′, overhang bases included.</param>
/// <param name="FivePrimeOverhang">Unpaired prefix of <paramref name="Sequence"/>; empty for a blunt end.</param>
/// <param name="ThreePrimeOverhang">Unpaired suffix of <paramref name="Sequence"/>; empty for a blunt end.</param>
/// <param name="IsCircular">Whether the fragment has been circularised (both ends sealed).</param>
public sealed record DnaFragment(
    string Id,
    string Sequence,
    string FivePrimeOverhang = "",
    string ThreePrimeOverhang = "",
    bool IsCircular = false)
{
    /// <summary>Total length of the top strand, overhang bases included.</summary>
    public int Length => Sequence.Length;

    /// <summary>GC fraction (0–1) of the whole top strand.</summary>
    public double GcContent => DnaUtility.GcFraction(Sequence);

    /// <summary>Whether both ends are blunt.</summary>
    public bool IsBlunt => FivePrimeOverhang.Length == 0 && ThreePrimeOverhang.Length == 0;

    /// <summary>Annotations carried along by the fragment, in fragment-local coordinates.</summary>
    public IReadOnlyList<SequenceFeature> Features { get; init; } = [];

    /// <summary>
    /// The double-stranded core: the top strand with both overhangs stripped. This is the part that is
    /// base-paired on both strands.
    /// </summary>
    public string Duplex
    {
        get
        {
            int start = FivePrimeOverhang.Length;
            int end = Sequence.Length - ThreePrimeOverhang.Length;
            return end <= start ? string.Empty : Sequence[start..end];
        }
    }

    /// <summary>
    /// Whether this fragment's 3′ end can be ligated to <paramref name="other"/>'s 5′ end.
    /// </summary>
    /// <remarks>
    /// Both overhangs are written as <b>top-strand text of the window they occupy</b> (see the type remarks),
    /// so two ends anneal exactly when those two strings are <i>equal</i>: the upstream fragment's protruding
    /// strand is the complement of its own 3′ window, and that complement base-pairs with the downstream
    /// fragment's 5′ window when the two windows read the same. Comparing with a reverse complement instead
    /// would only ever succeed for palindromic overhangs, which is wrong for a type IIS or 3′-overhang end.
    /// The comparison is case-insensitive and U is read as T. Two blunt ends are always joinable (blunt
    /// ligation is inefficient in the laboratory but it is chemically valid, so it is allowed here). A blunt
    /// end never joins a sticky one.
    /// </remarks>
    /// <param name="other">The candidate downstream fragment.</param>
    public bool CanLigateTo(DnaFragment other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (ThreePrimeOverhang.Length == 0 && other.FivePrimeOverhang.Length == 0) return true;
        if (ThreePrimeOverhang.Length != other.FivePrimeOverhang.Length) return false;
        return string.Equals(
            DnaUtility.NormalizeDna(ThreePrimeOverhang),
            DnaUtility.NormalizeDna(other.FivePrimeOverhang),
            StringComparison.OrdinalIgnoreCase);
    }

    /// <summary>
    /// The same molecule viewed from the other strand.
    /// </summary>
    /// <remarks>
    /// The whole top strand is reverse-complemented, which turns the old 3′ overhang into the new 5′ overhang
    /// and vice versa. Both overhang strings are themselves reverse-complemented, because they are written on
    /// whichever strand is now the top one. The operation is an involution: applying it twice returns the
    /// original record. Features are dropped, since their coordinates would have to be mirrored against a
    /// possibly asymmetric molecule and callers that need them should flip them explicitly.
    /// </remarks>
    public DnaFragment ReverseComplement() => new(
        Id,
        DnaUtility.ReverseComplement(Sequence),
        DnaUtility.ReverseComplement(ThreePrimeOverhang),
        DnaUtility.ReverseComplement(FivePrimeOverhang),
        IsCircular);

    /// <summary>
    /// Ligates <paramref name="other"/> downstream of this fragment.
    /// </summary>
    /// <remarks>
    /// The shared end is present once in each fragment's top strand, so the junction is formed by dropping the
    /// downstream fragment's 5′ overhang (identical, after annealing, to the bases already contributed by this
    /// fragment's 3′ overhang) and concatenating. The product keeps this fragment's 5′ end and the other
    /// fragment's 3′ end. Features of both inputs are carried over with corrected coordinates.
    /// </remarks>
    /// <param name="other">The downstream fragment.</param>
    /// <param name="id">Identifier for the product; defaults to <c>a+b</c>.</param>
    /// <exception cref="InvalidOperationException">The ends are not compatible, or an input is circular.</exception>
    public DnaFragment Ligate(DnaFragment other, string? id = null)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (IsCircular || other.IsCircular)
            throw new InvalidOperationException("A circular fragment has no free ends to ligate.");
        if (!CanLigateTo(other))
            throw new InvalidOperationException(
                $"Incompatible ends: 3' overhang '{ThreePrimeOverhang}' of '{Id}' does not anneal to 5' overhang '{other.FivePrimeOverhang}' of '{other.Id}'.");

        int drop = other.FivePrimeOverhang.Length;
        string joined = Sequence + other.Sequence[drop..];
        var features = new List<SequenceFeature>(Features.Count + other.Features.Count);
        features.AddRange(Features);
        int shift = Sequence.Length - drop;
        foreach (var feature in other.Features) features.Add(ShiftFeature(feature, shift));

        return new DnaFragment(id ?? $"{Id}+{other.Id}", joined, FivePrimeOverhang, other.ThreePrimeOverhang)
        {
            Features = features,
        };
    }

    /// <summary>
    /// Seals the fragment's own two ends into a circle.
    /// </summary>
    /// <remarks>
    /// Self-ligation requires the 3′ overhang to anneal to the fragment's own 5′ overhang. The duplicated
    /// overhang is removed from the 3′ end, so the circular sequence contains it exactly once.
    /// </remarks>
    /// <exception cref="InvalidOperationException">The two ends are not compatible.</exception>
    public DnaFragment Circularise()
    {
        if (IsCircular) return this;
        if (!CanLigateTo(this))
            throw new InvalidOperationException(
                $"Fragment '{Id}' cannot self-circularise: 3' overhang '{ThreePrimeOverhang}' does not anneal to its own 5' overhang '{FivePrimeOverhang}'.");

        int drop = ThreePrimeOverhang.Length;
        string circle = drop == 0 ? Sequence : Sequence[..^drop];
        return new DnaFragment(Id, circle, string.Empty, string.Empty, IsCircular: true) { Features = Features };
    }

    /// <summary>
    /// The fragment as a <see cref="Biotech.Net.Core.Sequence"/>, carrying its features.
    /// </summary>
    /// <remarks>
    /// The returned sequence is the full top strand including overhang bases; single-strandedness is not
    /// representable in <see cref="Biotech.Net.Core.Sequence"/>, so the overhangs are recorded as annotations
    /// instead.
    /// </remarks>
    public Biotech.Net.Core.Sequence ToSequence()
    {
        var sequence = Biotech.Net.Core.Sequence.Dna(Sequence, Id);
        sequence.Annotations["topology"] = IsCircular ? "circular" : "linear";
        if (FivePrimeOverhang.Length > 0) sequence.Annotations["five_prime_overhang"] = FivePrimeOverhang;
        if (ThreePrimeOverhang.Length > 0) sequence.Annotations["three_prime_overhang"] = ThreePrimeOverhang;
        foreach (var feature in Features) sequence.Features.Add(feature);
        return sequence;
    }

    /// <summary>Readable summary used in assembly step descriptions.</summary>
    public override string ToString()
        => $"{Id} [{FivePrimeOverhang}]{Length} bp[{ThreePrimeOverhang}]{(IsCircular ? " circular" : "")}";

    internal static SequenceFeature ShiftFeature(SequenceFeature feature, int shift)
    {
        var segments = new List<(int Start, int End)>(feature.Location.Segments.Count);
        foreach (var (start, end) in feature.Location.Segments) segments.Add((start + shift, end + shift));
        var moved = new SequenceFeature(feature.Type, new FeatureLocation(segments, feature.Location.Strand));
        foreach (var pair in feature.Qualifiers)
            foreach (var value in pair.Value)
                moved.With(pair.Key, value);
        return moved;
    }
}
