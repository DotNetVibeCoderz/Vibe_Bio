using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Cloning;

/// <summary>
/// A circular double-stranded DNA molecule with an arbitrary origin and a feature map.
/// </summary>
/// <remarks>
/// <para>
/// Coordinates are 0-based on the top strand and wrap modulo <see cref="Length"/>. Position 0 is the
/// "origin" of the numbering only; a circle has no intrinsic start, so <see cref="Rotate"/> renumbers the
/// molecule (and its features) without changing the physical DNA.
/// </para>
/// <para>
/// A digest of a circle behaves differently from a digest of a linear molecule: <i>n</i> cuts give exactly
/// <i>n</i> fragments, so a single cutter linearises the plasmid into one fragment whose length equals the
/// plasmid length, and a double cutter gives two fragments whose lengths sum to the plasmid length. That
/// invariant is what <see cref="Fragments"/> guarantees.
/// </para>
/// </remarks>
public sealed class Plasmid
{
    /// <summary>Creates a plasmid from a name and a circular top-strand sequence.</summary>
    /// <param name="name">Construct name, used in reports and in the derived sequence identifier.</param>
    /// <param name="sequence">Top strand, written from the chosen origin; case and U/T are normalised.</param>
    public Plasmid(string name, string sequence)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentException.ThrowIfNullOrEmpty(sequence);
        Name = name;
        Sequence = DnaUtility.NormalizeDna(sequence);
    }

    /// <summary>Construct name.</summary>
    public string Name { get; }

    /// <summary>Top strand from the current origin, 5′→3′.</summary>
    public string Sequence { get; }

    /// <summary>Number of base pairs in the circle.</summary>
    public int Length => Sequence.Length;

    /// <summary>Always true: a <see cref="Plasmid"/> models a closed circular molecule.</summary>
    public bool IsCircular => true;

    /// <summary>Feature map in current (post-rotation) coordinates.</summary>
    public IList<SequenceFeature> Features { get; } = [];

    /// <summary>
    /// Reads <paramref name="length"/> bases starting at <paramref name="start"/>, wrapping past the origin.
    /// </summary>
    /// <param name="start">0-based start; may be negative or beyond the end and is reduced modulo the length.</param>
    /// <param name="length">Number of bases to read; may exceed the plasmid length, repeating the circle.</param>
    public string Slice(int start, int length)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        if (length == 0) return string.Empty;
        int n = Length;
        var chars = new char[length];
        int index = Modulo(start, n);
        for (int i = 0; i < length; i++)
        {
            chars[i] = Sequence[index];
            index++;
            if (index == n) index = 0;
        }
        return new string(chars);
    }

    /// <summary>
    /// Returns the same molecule renumbered so that the base currently at <paramref name="newOrigin"/> becomes
    /// position 0. Features move with the sequence.
    /// </summary>
    /// <remarks>
    /// A feature that would be split by the new origin is kept as a two-segment join, so no annotation is lost.
    /// </remarks>
    /// <param name="newOrigin">The current position that becomes the new position 0.</param>
    public Plasmid Rotate(int newOrigin)
    {
        int n = Length;
        int shift = Modulo(newOrigin, n);
        var rotated = new Plasmid(Name, Slice(shift, n));
        foreach (var feature in Features)
        {
            var segments = new List<(int Start, int End)>();
            foreach (var (start, end) in feature.Location.Segments)
            {
                int s = Modulo(start - shift, n);
                int e = s + (end - start);
                if (e <= n) segments.Add((s, e));
                else { segments.Add((s, n)); segments.Add((0, e - n)); }
            }
            var moved = new SequenceFeature(feature.Type, new FeatureLocation(segments, feature.Location.Strand));
            foreach (var pair in feature.Qualifiers)
                foreach (var value in pair.Value)
                    moved.With(pair.Key, value);
            rotated.Features.Add(moved);
        }
        return rotated;
    }

    /// <summary>
    /// Returns a copy with one more annotation. The plasmid itself is not modified, so a map can be built up
    /// fluently.
    /// </summary>
    /// <param name="type">Feature key, e.g. <c>CDS</c>, <c>promoter</c>, <c>rep_origin</c>.</param>
    /// <param name="start">0-based inclusive start.</param>
    /// <param name="end">0-based exclusive end.</param>
    /// <param name="strand">Which strand the feature is read from.</param>
    /// <param name="label">Value stored under the <c>label</c> qualifier.</param>
    public Plasmid WithFeature(string type, int start, int end, Strand strand, string label)
    {
        var copy = new Plasmid(Name, Sequence);
        foreach (var feature in Features) copy.Features.Add(feature);
        copy.Features.Add(new SequenceFeature(type, FeatureLocation.Span(start, end, strand)).With("label", label));
        return copy;
    }

    /// <summary>All sites for one enzyme, including sites that span the origin.</summary>
    /// <param name="enzyme">The enzyme to scan for.</param>
    public IReadOnlyList<RestrictionSite> Digest(RestrictionEnzyme enzyme)
    {
        ArgumentNullException.ThrowIfNull(enzyme);
        return enzyme.Find(Sequence, circular: true);
    }

    /// <summary>
    /// Performs a real digest and returns the resulting fragments in order around the circle.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Every cut of a circle starts one fragment and ends the previous one, so <i>n</i> distinct top-strand
    /// cuts yield <i>n</i> fragments.
    /// </para>
    /// <para>
    /// <b>Overhang convention.</b> A fragment's top strand is written over the union of the extents of its two
    /// strands, that is from the leftmost of the two cuts of its own site to the rightmost of the two cuts of
    /// the next site. The single-stranded window of a cut therefore appears in <i>both</i> neighbouring
    /// fragments — as the <see cref="DnaFragment.ThreePrimeOverhang"/> of the upstream one and the
    /// <see cref="DnaFragment.FivePrimeOverhang"/> of the downstream one — and both are written as
    /// <b>top-strand text of the same window</b>, so two compatible ends carry the <i>same</i> string (see
    /// <see cref="DnaFragment.CanLigateTo"/>). This holds unchanged for a 5′-overhang enzyme such as EcoRI, for
    /// a 3′-overhang enzyme such as PstI (where the bottom-strand cut precedes the top-strand one) and for the
    /// non-palindromic overhangs of a type IIS enzyme. Because the shared window is written twice, the fragment
    /// lengths sum to the plasmid length plus one overhang length per cut; <see cref="DnaFragment.Ligate"/>
    /// removes the duplicate again, so a digest followed by a ligation round trip reproduces the circle exactly.
    /// </para>
    /// <para>
    /// With no cut site at all the plasmid is returned unchanged as a single circular fragment.
    /// </para>
    /// </remarks>
    /// <param name="enzymes">One or more enzymes applied in the same reaction (a double or triple digest).</param>
    public IReadOnlyList<DnaFragment> Fragments(params RestrictionEnzyme[] enzymes)
    {
        ArgumentNullException.ThrowIfNull(enzymes);
        var sites = new List<RestrictionSite>();
        foreach (var enzyme in enzymes) sites.AddRange(Digest(enzyme));

        // One fragment boundary per distinct top-strand cut. Two enzymes nicking the same bond is a
        // degenerate case; keep one boundary.
        var byCut = new SortedDictionary<int, RestrictionSite>();
        foreach (var site in sites) byCut.TryAdd(site.TopCutPosition, site);

        if (byCut.Count == 0)
            return [new DnaFragment(Name, Sequence, string.Empty, string.Empty, IsCircular: true)];

        var boundaries = new List<RestrictionSite>(byCut.Values);
        var fragments = new List<DnaFragment>(boundaries.Count);
        int n = Length;

        for (int i = 0; i < boundaries.Count; i++)
        {
            var from = boundaries[i];
            var to = boundaries[(i + 1) % boundaries.Count];

            // Signed overhang length comes from the enzyme, not from the difference of the two cut positions:
            // on a circle those positions are reduced modulo the length, so a site that spans the origin would
            // otherwise report a nonsensical overhang.
            int fivePrime = Math.Abs(from.Enzyme.OverhangLength);
            int threePrime = Math.Abs(to.Enzyme.OverhangLength);

            // Left edge of a cut = the leftmost of its two nicks: the top-strand nick for a 5' overhang (or a
            // blunt cut), |overhang| bases further left for a 3' overhang.
            int startBase = LeftEdge(from, n);
            int endBase = LeftEdge(to, n);

            // Distance between the two left edges; a single cut linearises the whole circle.
            int span = boundaries.Count == 1 ? n : Modulo(endBase - startBase, n);

            // The fragment's top strand runs to the right edge of the next cut, so that the single-stranded
            // window of that cut is present (it is shared with the downstream fragment).
            int topStart = startBase;
            int topLength = span + threePrime;

            string top = Slice(topStart, topLength);
            string fivePrimeSeq = fivePrime == 0 ? string.Empty : top[..fivePrime];
            string threePrimeSeq = threePrime == 0 ? string.Empty : top[^threePrime..];

            var id = $"{Name}_{from.Enzyme.Name}{i + 1}";
            var fragment = new DnaFragment(id, top, fivePrimeSeq, threePrimeSeq)
            {
                Features = CarryFeatures(topStart, topLength),
            };
            fragments.Add(fragment);
        }

        return fragments;
    }

    /// <summary>
    /// The enzymes from <paramref name="candidates"/> that cut this plasmid exactly once — the classic
    /// "where can I open the vector" query.
    /// </summary>
    /// <param name="candidates">Enzymes to test; defaults to the built-in catalogue.</param>
    public IReadOnlyList<RestrictionEnzyme> UniqueCutters(IEnumerable<RestrictionEnzyme>? candidates = null)
    {
        var result = new List<RestrictionEnzyme>();
        foreach (var enzyme in candidates ?? RestrictionEnzyme.All)
            if (Digest(enzyme).Count == 1) result.Add(enzyme);
        return result;
    }

    /// <summary>
    /// A compact, human-readable map: name, length, topology and the feature list. This is a summary for
    /// logs and console output, not a GenBank flat-file writer.
    /// </summary>
    public string ToGenBankSummary()
    {
        var builder = new StringBuilder();
        builder.Append("LOCUS       ").Append(Name).Append("  ").Append(Length).AppendLine(" bp  DNA  circular");
        builder.Append("DEFINITION  ").Append(Name).AppendLine(" (synthetic construct map)");
        builder.AppendLine("FEATURES             Location/Qualifiers");
        foreach (var feature in Features)
        {
            builder.Append("     ").Append(feature.Type.PadRight(16)).Append(feature.Location).AppendLine();
            builder.Append("                      /label=\"").Append(feature.Name).AppendLine("\"");
        }
        builder.Append("ORIGIN      ").Append(Length).AppendLine(" bp");
        return builder.ToString();
    }

    /// <summary>The plasmid as a <see cref="Biotech.Net.Core.Sequence"/> with its features and a circular topology annotation.</summary>
    /// <param name="id">Identifier for the sequence; defaults to the plasmid name.</param>
    public Biotech.Net.Core.Sequence ToSequence(string? id = null)
    {
        var sequence = Biotech.Net.Core.Sequence.Dna(Sequence, id ?? Name);
        sequence.Annotations["topology"] = "circular";
        foreach (var feature in Features) sequence.Features.Add(feature);
        return sequence;
    }

    /// <summary>
    /// A short synthetic demonstration backbone. It is <b>not</b> the sequence of any real published vector:
    /// it is an invented 240 bp circle carrying a small polylinker so that examples and tests have something
    /// to cut. Do not use it for laboratory work.
    /// </summary>
    public static Plasmid DemonstrationBackbone { get; } = BuildDemonstrationBackbone();

    private static Plasmid BuildDemonstrationBackbone()
    {
        // A polylinker with single EcoRI, BamHI, SalI, PstI and HindIII sites, padded with arbitrary filler
        // that contains none of those sites.
        const string polylinker = "GAATTCGGTACCCGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT";
        const string filler =
            "ACGTACGTTAGCATCAGCATTACGACTTAGCATGCATAATCGATCGTAAC" +
            "TTACGCATTAGCATCGATTACGCATAGCATTACGGATCATTACGCATTAG" +
            "CATCAGTTACGCATCAGCATTAGCATCAGCTTACGACTTACGCATCAGTA" +
            "ACGTTAGCATCAGCATTACGACTTAGCATGCATAATCGATCGTAACTTAC";
        var plasmid = new Plasmid("pDemoBackbone", polylinker + filler);
        return plasmid
            .WithFeature("misc_feature", 0, polylinker.Length, Strand.Forward, "MCS")
            .WithFeature("misc_feature", polylinker.Length, polylinker.Length + 100, Strand.Forward, "filler_A");
    }

    // Leftmost nick of a cut on the top-strand axis, reduced into [0, n).
    private static int LeftEdge(RestrictionSite site, int n)
        => site.Enzyme.OverhangLength < 0
            ? Modulo(site.TopCutPosition - Math.Abs(site.Enzyme.OverhangLength), n)
            : Modulo(site.TopCutPosition, n);

    internal static int Modulo(int value, int modulus)
    {
        if (modulus <= 0) return 0;
        int r = value % modulus;
        return r < 0 ? r + modulus : r;
    }

    private List<SequenceFeature> CarryFeatures(int start, int length)
    {
        var carried = new List<SequenceFeature>();
        int n = Length;
        foreach (var feature in Features)
        {
            var segments = new List<(int Start, int End)>();
            bool complete = true;
            foreach (var (s, e) in feature.Location.Segments)
            {
                int relative = Modulo(s - start, n);
                if (relative + (e - s) > length) { complete = false; break; }
                segments.Add((relative, relative + (e - s)));
            }
            if (!complete || segments.Count == 0) continue;
            var moved = new SequenceFeature(feature.Type, new FeatureLocation(segments, feature.Location.Strand));
            foreach (var pair in feature.Qualifiers)
                foreach (var value in pair.Value)
                    moved.With(pair.Key, value);
            carried.Add(moved);
        }
        return carried;
    }

    /// <summary>Readable summary, e.g. <c>pDemoBackbone (290 bp, circular, 2 features)</c>.</summary>
    public override string ToString() => $"{Name} ({Length} bp, circular, {Features.Count} features)";
}
