using System.Text;
using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Cloning;

/// <summary>
/// Design and simulation of Gibson (isothermal, exonuclease-based) assembly: sequence-independent joining of
/// DNA fragments that share identical terminal sequence.
/// </summary>
/// <remarks>
/// <para>
/// A Gibson reaction chews back the 5′ ends of every fragment, lets the exposed 3′ single strands anneal where
/// the fragments share terminal homology, then fills in and seals. The practical consequence for software is
/// simple: <b>the shared homology must be present once in the product, not twice</b>. Every routine here keeps
/// that invariant — <see cref="Simulate(IEnumerable{DnaFragment}, GibsonOptions?)"/> drops the downstream
/// copy of each overlap, and <see cref="DesignOverlaps"/> splits each overlap across the junction so that the
/// two tailed amplicons reconstruct it exactly once.
/// </para>
/// <para>
/// Overlap design follows the usual rules of thumb (Gibson DG et al., Nat Methods 2009; NEBuilder guidance):
/// 15–40 nt of homology, around 30 nt as a default, a melting temperature at or above roughly 48 °C and a
/// balanced GC content. Melting temperatures come from <see cref="NearestNeighbourTm"/>, so they are nearest
/// neighbour values, not the 2+4 rule.
/// </para>
/// </remarks>
public static class GibsonAssembly
{
    /// <summary>
    /// Designs one homology arm per junction for a set of fragments that are to be assembled in the given
    /// order, together with the tailed PCR primers that install those arms.
    /// </summary>
    /// <remarks>
    /// Each overlap straddles its junction: the 5′ half is taken from the end of the upstream fragment and the
    /// 3′ half from the start of the downstream fragment. The upstream fragment is therefore amplified with a
    /// reverse primer whose 5′ tail is the reverse complement of the start of its neighbour, and the
    /// downstream fragment with a forward primer whose 5′ tail is the end of its neighbour. Candidate overlap
    /// lengths between <see cref="GibsonOptions.MinOverlap"/> and <see cref="GibsonOptions.MaxOverlap"/> are
    /// scored on melting temperature, GC content, distance from <see cref="GibsonOptions.TargetOverlap"/>,
    /// self-complementarity and uniqueness within the whole assembly, and the best scoring length wins.
    /// </remarks>
    /// <param name="fragments">Fragment sequences in the intended assembly order; at least two.</param>
    /// <param name="options">Design window; <see langword="null"/> uses <see cref="GibsonOptions.Default"/>.</param>
    /// <returns>The junctions, the fragments with their homology arms already appended, and any warnings.</returns>
    /// <exception cref="ArgumentException">Fewer than two fragments, or a fragment too short to carry an arm.</exception>
    public static GibsonDesign DesignOverlaps(IReadOnlyList<string> fragments, GibsonOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(fragments);
        var opts = options ?? GibsonOptions.Default;
        if (fragments.Count < 2)
            throw new ArgumentException("Gibson overlap design needs at least two fragments.", nameof(fragments));

        var parts = new string[fragments.Count];
        for (int i = 0; i < fragments.Count; i++)
        {
            parts[i] = DnaUtility.NormalizeDna(fragments[i] ?? throw new ArgumentException("A fragment is null.", nameof(fragments)));
            if (parts[i].Length < opts.MaxOverlap + opts.MaxAnnealingLength)
                throw new ArgumentException(
                    $"Fragment {i} is {parts[i].Length} nt, too short for an overlap of up to {opts.MaxOverlap} nt plus an annealing region of up to {opts.MaxAnnealingLength} nt.",
                    nameof(fragments));
        }

        string whole = string.Concat(parts);
        int junctionCount = opts.Circular ? parts.Length : parts.Length - 1;
        var junctions = new List<GibsonJunction>(junctionCount);
        var warnings = new List<string>();

        for (int i = 0; i < junctionCount; i++)
        {
            int j = (i + 1) % parts.Length;
            var junction = DesignJunction(i, j, parts[i], parts[j], whole, opts);
            junctions.Add(junction);
            if (!junction.MeetsTm)
                warnings.Add($"Junction {i}->{j}: overlap Tm {junction.Tm:F1} °C is below the requested {opts.MinOverlapTm:F1} °C.");
            if (!junction.MeetsGc)
                warnings.Add($"Junction {i}->{j}: overlap GC {junction.GcContent:P0} is outside {opts.MinGc:P0}–{opts.MaxGc:P0}.");
            if (!junction.IsUnique)
                warnings.Add($"Junction {i}->{j}: overlap '{junction.Overlap}' is not unique in the assembly.");
            if (junction.IsSelfComplementary)
                warnings.Add($"Junction {i}->{j}: overlap '{junction.Overlap}' is self-complementary.");
        }

        var tailed = new string[parts.Length];
        for (int i = 0; i < parts.Length; i++)
        {
            string head = string.Empty, tail = string.Empty;
            foreach (var junction in junctions)
            {
                if (junction.DownstreamIndex == i) head = junction.UpstreamContribution;
                if (junction.UpstreamIndex == i) tail = junction.DownstreamContribution;
            }
            tailed[i] = head + parts[i] + tail;
        }

        return new GibsonDesign(junctions, tailed, opts.Circular, warnings);
    }

    private static GibsonJunction DesignJunction(int i, int j, string upstream, string downstream, string whole, GibsonOptions opts)
    {
        GibsonJunction? best = null;
        double bestScore = double.NegativeInfinity;

        for (int k = opts.MinOverlap; k <= opts.MaxOverlap; k++)
        {
            int h1 = (k + 1) / 2;
            int h2 = k - h1;
            if (h1 > upstream.Length || h2 > downstream.Length) continue;

            string upstreamPart = upstream[^h1..];
            string downstreamPart = h2 == 0 ? string.Empty : downstream[..h2];
            string overlap = upstreamPart + downstreamPart;

            double tm = NearestNeighbourTm.Tm(overlap, opts.Salt);
            double gc = DnaUtility.GcFraction(overlap);
            bool meetsTm = tm >= opts.MinOverlapTm;
            bool meetsGc = gc >= opts.MinGc && gc <= opts.MaxGc;
            bool selfComplementary = NearestNeighbourTm.IsSelfComplementary(overlap);
            bool unique = CountOccurrences(whole, overlap) <= 1 && CountOccurrences(whole, DnaUtility.ReverseComplement(overlap)) == 0;

            double score = 0;
            score += meetsTm ? 10 : -(opts.MinOverlapTm - tm);
            score += meetsGc ? 5 : -20 * (gc < opts.MinGc ? opts.MinGc - gc : gc - opts.MaxGc);
            score -= 0.15 * Math.Abs(k - opts.TargetOverlap);
            if (selfComplementary) score -= 50;
            if (!unique) score -= 50;

            if (score <= bestScore) continue;
            bestScore = score;

            int anneal = ChooseAnnealingLength(downstream, opts, h2);
            int annealUp = ChooseAnnealingLength(DnaUtility.ReverseComplement(upstream), opts, h1);
            string forward = upstreamPart + downstream[..anneal];
            string reverse = DnaUtility.ReverseComplement(upstream[^annealUp..] + downstreamPart);

            best = new GibsonJunction(
                i, j, overlap, upstreamPart, downstreamPart, tm, gc,
                reverse, forward,
                NearestNeighbourTm.Tm(downstream[..anneal], opts.Salt),
                NearestNeighbourTm.Tm(upstream[^annealUp..], opts.Salt),
                meetsTm, meetsGc, unique, selfComplementary);
        }

        return best ?? throw new ArgumentException(
            $"No overlap between {opts.MinOverlap} and {opts.MaxOverlap} nt fits between fragments {i} and {j}.");
    }

    private static int ChooseAnnealingLength(string sequence, GibsonOptions opts, int minimum)
    {
        int lower = Math.Max(opts.MinAnnealingLength, minimum);
        int upper = Math.Min(opts.MaxAnnealingLength, sequence.Length);
        if (upper < lower) return Math.Min(Math.Max(lower, 1), sequence.Length);
        for (int length = lower; length <= upper; length++)
            if (NearestNeighbourTm.Tm(sequence[..length], opts.Salt) >= opts.MinPrimerTm) return length;
        return upper;
    }

    private static int CountOccurrences(string haystack, string needle)
    {
        if (needle.Length == 0) return 0;
        int count = 0, index = 0;
        while ((index = haystack.IndexOf(needle, index, StringComparison.Ordinal)) >= 0)
        {
            count++;
            index++;
        }
        return count;
    }

    /// <summary>
    /// Simulates a Gibson reaction on fragments that already carry their homology arms.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The fragments are chained by exact terminal homology of at least <see cref="GibsonOptions.MinOverlap"/>
    /// nucleotides: the 3′ end of the growing product is matched against the 5′ end of every unused fragment
    /// (and, when <see cref="GibsonOptions.AllowReverseComplement"/> is set, against its reverse complement
    /// too). The first fragment of the input fixes the orientation and the starting point.
    /// </para>
    /// <para>
    /// Each overlap is written into the product exactly once. If the chain closes — the last fragment's end
    /// overlaps the first fragment's start — the duplicate is removed and the product is returned circular.
    /// An end that matches two different fragments is reported as a failure rather than resolved by guessing,
    /// and so is an end that matches nothing.
    /// </para>
    /// </remarks>
    /// <param name="fragments">The fragments in the reaction; order is irrelevant except that the first one anchors the product.</param>
    /// <param name="options">Overlap rules; <see langword="null"/> uses <see cref="GibsonOptions.Default"/>.</param>
    /// <returns>The product and the junctions used, or a failed result carrying the reason.</returns>
    public static GibsonAssemblyResult Simulate(IEnumerable<DnaFragment> fragments, GibsonOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(fragments);
        var opts = options ?? GibsonOptions.Default;
        var items = new List<(string Id, string Sequence)>();
        foreach (var fragment in fragments)
        {
            ArgumentNullException.ThrowIfNull(fragment);
            items.Add((fragment.Id, DnaUtility.NormalizeDna(fragment.Sequence)));
        }
        return Assemble(items, opts);
    }

    /// <summary>
    /// Simulates a Gibson reaction on plain sequences, which are given the identifiers <c>F1</c>, <c>F2</c>, …
    /// </summary>
    /// <param name="sequences">Fragment top strands; the first one anchors the product.</param>
    /// <param name="options">Overlap rules; <see langword="null"/> uses <see cref="GibsonOptions.Default"/>.</param>
    /// <returns>The product and the junctions used, or a failed result carrying the reason.</returns>
    public static GibsonAssemblyResult Simulate(IEnumerable<string> sequences, GibsonOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var items = new List<(string Id, string Sequence)>();
        int n = 1;
        foreach (var sequence in sequences)
        {
            ArgumentNullException.ThrowIfNull(sequence);
            items.Add(($"F{n++}", DnaUtility.NormalizeDna(sequence)));
        }
        return Assemble(items, options ?? GibsonOptions.Default);
    }

    private static GibsonAssemblyResult Assemble(List<(string Id, string Sequence)> items, GibsonOptions opts)
    {
        if (items.Count == 0)
            return GibsonAssemblyResult.Failed("No fragments were supplied.");
        if (items.Count == 1)
        {
            var only = items[0];
            int self = BestOverlap(only.Sequence, only.Sequence, opts.MinOverlap);
            if (self <= 0)
                return new GibsonAssemblyResult(true, new DnaFragment(only.Id, only.Sequence), [], false, null);
            var circle = new DnaFragment(only.Id, only.Sequence[..^self], IsCircular: true);
            return new GibsonAssemblyResult(true, circle, [new GibsonJoin(only.Id, only.Id, only.Sequence[..self], self)], true, null);
        }

        var used = new bool[items.Count];
        used[0] = true;
        var builder = new StringBuilder(items[0].Sequence);
        string currentId = items[0].Id;
        string currentEnd = items[0].Sequence;
        var joins = new List<GibsonJoin>();

        for (int step = 1; step < items.Count; step++)
        {
            int chosen = -1;
            bool chosenFlipped = false;
            string chosenSequence = string.Empty;
            int chosenOverlap = 0;
            int matches = 0;
            string ambiguous = string.Empty;

            for (int c = 0; c < items.Count; c++)
            {
                if (used[c]) continue;
                int forward = BestOverlap(currentEnd, items[c].Sequence, opts.MinOverlap);
                int reverse = opts.AllowReverseComplement
                    ? BestOverlap(currentEnd, DnaUtility.ReverseComplement(items[c].Sequence), opts.MinOverlap)
                    : 0;
                if (forward <= 0 && reverse <= 0) continue;

                matches++;
                ambiguous = ambiguous.Length == 0 ? items[c].Id : $"{ambiguous}, {items[c].Id}";
                bool flip = reverse > forward;
                chosen = c;
                chosenFlipped = flip;
                chosenSequence = flip ? DnaUtility.ReverseComplement(items[c].Sequence) : items[c].Sequence;
                chosenOverlap = flip ? reverse : forward;
                if (forward > 0 && reverse > 0)
                    return GibsonAssemblyResult.Failed(
                        $"Ambiguous assembly: the 3' end of '{currentId}' matches both orientations of '{items[c].Id}'.");
            }

            if (matches == 0)
                return GibsonAssemblyResult.Failed(
                    $"No fragment shares at least {opts.MinOverlap} nt of terminal homology with the 3' end of '{currentId}'.");
            if (matches > 1)
                return GibsonAssemblyResult.Failed(
                    $"Ambiguous assembly: the 3' end of '{currentId}' matches more than one fragment ({ambiguous}).");

            used[chosen] = true;
            string overlapSequence = chosenSequence[..chosenOverlap];
            builder.Append(chosenSequence, chosenOverlap, chosenSequence.Length - chosenOverlap);
            joins.Add(new GibsonJoin(currentId, items[chosen].Id, overlapSequence, chosenOverlap));
            currentId = chosenFlipped ? items[chosen].Id + "(rc)" : items[chosen].Id;
            currentEnd = chosenSequence;
        }

        string product = builder.ToString();
        int closure = BestOverlap(currentEnd, items[0].Sequence, opts.MinOverlap);
        bool circular = closure > 0 && closure < product.Length;
        if (circular)
        {
            joins.Add(new GibsonJoin(currentId, items[0].Id, product[..closure], closure));
            product = product[..^closure];
        }

        var assembled = new DnaFragment(
            string.Join('+', items.ConvertAll(static item => item.Id)),
            product,
            IsCircular: circular);
        return new GibsonAssemblyResult(true, assembled, joins, circular, null);
    }

    /// <summary>
    /// Length of the longest exact terminal homology between the 3′ end of <paramref name="upstream"/> and the
    /// 5′ end of <paramref name="downstream"/>, or 0 when it is shorter than <paramref name="minimum"/>.
    /// </summary>
    /// <param name="upstream">The fragment whose 3′ end is tested.</param>
    /// <param name="downstream">The fragment whose 5′ end is tested.</param>
    /// <param name="minimum">Shortest homology that counts as a junction.</param>
    public static int BestOverlap(string upstream, string downstream, int minimum)
    {
        ArgumentNullException.ThrowIfNull(upstream);
        ArgumentNullException.ThrowIfNull(downstream);
        if (minimum < 1) minimum = 1;
        int max = Math.Min(upstream.Length, downstream.Length);
        if (ReferenceEquals(upstream, downstream) || upstream == downstream) max--;
        for (int length = max; length >= minimum; length--)
            if (string.CompareOrdinal(upstream, upstream.Length - length, downstream, 0, length) == 0)
                return length;
        return 0;
    }
}

/// <summary>
/// Design window for <see cref="GibsonAssembly"/>: overlap length, thermodynamic targets and the rules the
/// simulator uses to chain fragments.
/// </summary>
/// <param name="MinOverlap">Shortest homology arm, and the shortest terminal homology the simulator accepts.</param>
/// <param name="MaxOverlap">Longest homology arm considered by the designer.</param>
/// <param name="TargetOverlap">Preferred homology length; candidates are penalised by their distance from it.</param>
/// <param name="MinOverlapTm">Lowest acceptable nearest-neighbour melting temperature of the overlap, in °C.</param>
/// <param name="MinGc">Lowest acceptable GC fraction of the overlap.</param>
/// <param name="MaxGc">Highest acceptable GC fraction of the overlap.</param>
/// <param name="MinAnnealingLength">Shortest template-annealing region of a designed primer.</param>
/// <param name="MaxAnnealingLength">Longest template-annealing region of a designed primer.</param>
/// <param name="MinPrimerTm">Melting temperature, in °C, the annealing region should reach if it can.</param>
/// <param name="Circular">Whether the designer should also create the junction that closes the assembly into a circle.</param>
/// <param name="AllowReverseComplement">Whether the simulator may flip a fragment to make it fit.</param>
/// <param name="Salt">Salt and oligo concentrations for melting temperatures; <see langword="null"/> uses the default.</param>
public sealed record GibsonOptions(
    int MinOverlap = 20,
    int MaxOverlap = 40,
    int TargetOverlap = 30,
    double MinOverlapTm = 48.0,
    double MinGc = 0.40,
    double MaxGc = 0.65,
    int MinAnnealingLength = 18,
    int MaxAnnealingLength = 30,
    double MinPrimerTm = 58.0,
    bool Circular = true,
    bool AllowReverseComplement = false,
    SaltConditions? Salt = null)
{
    /// <summary>The default window: 20–40 nt of homology, optimum 30 nt, overlap Tm at or above 48 °C.</summary>
    public static GibsonOptions Default { get; } = new();
}

/// <summary>
/// One designed Gibson junction: the homology arm that joins two fragments and the tailed primers that build it.
/// </summary>
/// <remarks>
/// <see cref="Overlap"/> is <see cref="UpstreamContribution"/> followed by <see cref="DownstreamContribution"/>,
/// so it is present once in the upstream amplicon (which gains the downstream half as a 3′ tail) and once in
/// the downstream amplicon (which gains the upstream half as a 5′ tail). After assembly it appears once.
/// </remarks>
/// <param name="UpstreamIndex">Index of the fragment on the 5′ side of the junction.</param>
/// <param name="DownstreamIndex">Index of the fragment on the 3′ side of the junction.</param>
/// <param name="Overlap">The full homology arm, 5′→3′ on the top strand.</param>
/// <param name="UpstreamContribution">The part of <paramref name="Overlap"/> taken from the end of the upstream fragment.</param>
/// <param name="DownstreamContribution">The part of <paramref name="Overlap"/> taken from the start of the downstream fragment.</param>
/// <param name="Tm">Nearest-neighbour melting temperature of the overlap, in °C.</param>
/// <param name="GcContent">GC fraction (0–1) of the overlap.</param>
/// <param name="UpstreamReversePrimer">Reverse primer for the upstream fragment: 5′ tail = reverse complement of the downstream contribution.</param>
/// <param name="DownstreamForwardPrimer">Forward primer for the downstream fragment: 5′ tail = the upstream contribution.</param>
/// <param name="DownstreamAnnealingTm">Melting temperature, in °C, of the template-annealing part of the forward primer.</param>
/// <param name="UpstreamAnnealingTm">Melting temperature, in °C, of the template-annealing part of the reverse primer.</param>
/// <param name="MeetsTm">Whether the overlap reaches <see cref="GibsonOptions.MinOverlapTm"/>.</param>
/// <param name="MeetsGc">Whether the overlap GC lies inside the requested bounds.</param>
/// <param name="IsUnique">Whether the overlap occurs only at this junction in the whole assembly, on either strand.</param>
/// <param name="IsSelfComplementary">Whether the overlap anneals to itself (a hairpin or primer-dimer risk).</param>
public sealed record GibsonJunction(
    int UpstreamIndex,
    int DownstreamIndex,
    string Overlap,
    string UpstreamContribution,
    string DownstreamContribution,
    double Tm,
    double GcContent,
    string UpstreamReversePrimer,
    string DownstreamForwardPrimer,
    double DownstreamAnnealingTm,
    double UpstreamAnnealingTm,
    bool MeetsTm,
    bool MeetsGc,
    bool IsUnique,
    bool IsSelfComplementary)
{
    /// <summary>Length of the homology arm in nucleotides.</summary>
    public int Length => Overlap.Length;

    /// <summary>Whether every design rule was satisfied.</summary>
    public bool IsIdeal => MeetsTm && MeetsGc && IsUnique && !IsSelfComplementary;

    /// <summary>Readable summary, e.g. <c>0->1 30 nt Tm 62.4 °C GC 50%</c>.</summary>
    public override string ToString()
        => $"{UpstreamIndex}->{DownstreamIndex} {Length} nt Tm {Tm:F1} °C GC {GcContent:P0}";
}

/// <summary>
/// The result of <see cref="GibsonAssembly.DesignOverlaps"/>: one junction per neighbouring pair plus the
/// fragments with their homology arms already attached.
/// </summary>
/// <param name="Junctions">The designed junctions, in assembly order.</param>
/// <param name="TailedFragments">
/// The input fragments with the neighbouring halves of their junctions appended, i.e. exactly the amplicons the
/// designed primers would produce. Feeding these to <see cref="GibsonAssembly.Simulate(IEnumerable{string}, GibsonOptions?)"/>
/// reproduces the intended construct.
/// </param>
/// <param name="Circular">Whether a closing junction was designed.</param>
/// <param name="Warnings">Human-readable notes about junctions that missed a design rule; empty when all are ideal.</param>
public sealed record GibsonDesign(
    IReadOnlyList<GibsonJunction> Junctions,
    IReadOnlyList<string> TailedFragments,
    bool Circular,
    IReadOnlyList<string> Warnings)
{
    /// <summary>Whether every junction satisfied every design rule.</summary>
    public bool IsIdeal
    {
        get
        {
            foreach (var junction in Junctions)
                if (!junction.IsIdeal) return false;
            return true;
        }
    }
}

/// <summary>One junction actually used by the simulator.</summary>
/// <param name="UpstreamId">Identifier of the fragment contributing the 5′ side.</param>
/// <param name="DownstreamId">Identifier of the fragment contributing the 3′ side.</param>
/// <param name="Overlap">The shared terminal sequence, written once in the product.</param>
/// <param name="Length">Length of the shared sequence in nucleotides.</param>
public sealed record GibsonJoin(string UpstreamId, string DownstreamId, string Overlap, int Length)
{
    /// <summary>Readable summary, e.g. <c>F1->F2 (24 nt)</c>.</summary>
    public override string ToString() => $"{UpstreamId}->{DownstreamId} ({Length} nt)";
}

/// <summary>Outcome of <see cref="GibsonAssembly.Simulate(IEnumerable{string}, GibsonOptions?)"/>.</summary>
/// <param name="Success">Whether a single unambiguous product was assembled.</param>
/// <param name="Product">The assembled molecule, or <see langword="null"/> when assembly failed.</param>
/// <param name="Junctions">The overlaps consumed, in the order they were used.</param>
/// <param name="IsCircular">Whether the chain closed on itself.</param>
/// <param name="FailureReason">Why assembly failed, or <see langword="null"/> on success.</param>
public sealed record GibsonAssemblyResult(
    bool Success,
    DnaFragment? Product,
    IReadOnlyList<GibsonJoin> Junctions,
    bool IsCircular,
    string? FailureReason)
{
    /// <summary>A failed result carrying only a reason.</summary>
    /// <param name="reason">Why the reaction could not be resolved to one product.</param>
    public static GibsonAssemblyResult Failed(string reason) => new(false, null, [], false, reason);

    /// <summary>Readable summary, e.g. <c>Gibson: 1_842 bp circular, 3 junctions</c>.</summary>
    public override string ToString() => Success
        ? $"Gibson: {Product!.Length} bp {(IsCircular ? "circular" : "linear")}, {Junctions.Count} junctions"
        : $"Gibson failed: {FailureReason}";
}
