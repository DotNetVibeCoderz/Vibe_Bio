using Biotech.Net.SyntheticBiology;
using Biotech.Net.SyntheticBiology.Cloning;
using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Gibson assembly and Golden Gate cloning: overlap and overhang design, and simulation of the reactions.
/// </summary>
public class AssemblyTests
{
    private static string RandomDna(int length, int seed)
    {
        var random = new Random(seed);
        return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[random.Next(4)]));
    }

    /// <summary>
    /// Random DNA with every site for <paramref name="enzyme"/> removed, on both strands — a part that has
    /// already been through domestication, which is what Golden Gate requires of its inputs.
    /// </summary>
    private static string DomesticatedDna(int length, int seed, RestrictionEnzyme enzyme)
    {
        var sequence = RandomDna(length, seed).ToCharArray();
        for (int guard = 0; guard < 100; guard++)
        {
            string current = new(sequence);
            var sites = enzyme.Find(current);
            if (sites.Count == 0) return current;
            foreach (var site in sites)
            {
                // Flip the middle base of the recognition site; A and T are never both part of GGTCTC/CGTCTC.
                int index = site.Position + enzyme.SiteLength / 2;
                sequence[index] = sequence[index] == 'A' ? 'T' : 'A';
            }
        }
        throw new InvalidOperationException("Could not domesticate the random sequence.");
    }

    // ---- Gibson: overlap design ----

    [Fact]
    public void DesignOverlaps_ProducesOneJunctionPerFragmentOfTheRequestedQuality()
    {
        var fragments = new[] { RandomDna(600, 1), RandomDna(600, 2), RandomDna(600, 3) };
        var options = GibsonOptions.Default;

        var design = GibsonAssembly.DesignOverlaps(fragments, options);

        // Circular by default: one junction per fragment, including the one closing the circle.
        Assert.Equal(3, design.Junctions.Count);
        Assert.True(design.Circular);

        foreach (var junction in design.Junctions)
        {
            Assert.InRange(junction.Length, options.MinOverlap, options.MaxOverlap);
            Assert.Equal(junction.Overlap.Length, junction.Length);
            // The arm straddles the junction: upstream end + downstream start.
            Assert.Equal(junction.UpstreamContribution + junction.DownstreamContribution, junction.Overlap);
            Assert.EndsWith(junction.UpstreamContribution, fragments[junction.UpstreamIndex], StringComparison.Ordinal);
            Assert.StartsWith(junction.DownstreamContribution, fragments[junction.DownstreamIndex], StringComparison.Ordinal);

            Assert.Equal(NearestNeighbourTm.Tm(junction.Overlap), junction.Tm, 6);
            Assert.Equal(DnaUtility.GcFraction(junction.Overlap), junction.GcContent, 12);
            Assert.Equal(junction.Tm >= options.MinOverlapTm, junction.MeetsTm);
        }
    }

    [Fact]
    public void DesignOverlaps_TailedPrimersCarryTheNeighbouringSequence()
    {
        var fragments = new[] { RandomDna(600, 11), RandomDna(600, 12) };
        var design = GibsonAssembly.DesignOverlaps(fragments, new GibsonOptions(Circular: false));

        var junction = Assert.Single(design.Junctions);

        // The downstream fragment's forward primer carries the upstream end as a 5' tail, then anneals to its own
        // fragment; the upstream fragment's reverse primer is the mirror image.
        Assert.StartsWith(junction.UpstreamContribution, junction.DownstreamForwardPrimer, StringComparison.Ordinal);
        Assert.StartsWith(DnaUtility.ReverseComplement(junction.DownstreamContribution),
                          junction.UpstreamReversePrimer, StringComparison.Ordinal);

        // The annealing region really is part of the fragment it amplifies.
        string downstreamAnnealing = junction.DownstreamForwardPrimer[junction.UpstreamContribution.Length..];
        Assert.StartsWith(downstreamAnnealing, fragments[1], StringComparison.Ordinal);

        string upstreamAnnealing = junction.UpstreamReversePrimer[junction.DownstreamContribution.Length..];
        Assert.EndsWith(DnaUtility.ReverseComplement(upstreamAnnealing), fragments[0], StringComparison.Ordinal);
    }

    [Fact]
    public void DesignOverlaps_LinearAssemblyHasOneJunctionFewer()
    {
        var fragments = new[] { RandomDna(600, 21), RandomDna(600, 22), RandomDna(600, 23) };
        var design = GibsonAssembly.DesignOverlaps(fragments, new GibsonOptions(Circular: false));
        Assert.Equal(2, design.Junctions.Count);
        Assert.False(design.Circular);
    }

    [Fact]
    public void DesignOverlaps_RejectsTooFewOrTooShortFragments()
    {
        Assert.Throws<ArgumentException>(() => GibsonAssembly.DesignOverlaps([RandomDna(600, 1)]));
        Assert.Throws<ArgumentException>(() => GibsonAssembly.DesignOverlaps([RandomDna(30, 1), RandomDna(600, 2)]));
    }

    // ---- Gibson: simulation ----

    /// <summary>Splits a sequence into <paramref name="pieces"/> pieces that share <paramref name="overlap"/> nt ends.</summary>
    private static string[] SplitWithOverlaps(string sequence, int pieces, int overlap, bool circular)
    {
        int step = sequence.Length / pieces;
        var parts = new string[pieces];
        for (int i = 0; i < pieces; i++)
        {
            int start = i * step;
            int end = i == pieces - 1 ? sequence.Length : (i + 1) * step;
            // Each piece reaches `overlap` bases into its successor.
            if (i < pieces - 1) end += overlap;
            else if (circular) { parts[i] = sequence[start..] + sequence[..overlap]; continue; }
            parts[i] = sequence[start..end];
        }
        return parts;
    }

    [Fact]
    public void Simulate_LinearAssemblyReproducesTheOriginalExactlyOnce()
    {
        string original = RandomDna(900, 31);
        var pieces = SplitWithOverlaps(original, 3, 30, circular: false);

        var result = GibsonAssembly.Simulate(pieces, new GibsonOptions(Circular: false));

        Assert.True(result.Success, result.FailureReason);
        Assert.False(result.IsCircular);
        Assert.Equal(original, result.Product!.Sequence);
        Assert.Equal(2, result.Junctions.Count);
        Assert.All(result.Junctions, j => Assert.Equal(30, j.Length));
        // The overlap appears exactly once in the product, which is the whole point.
        Assert.Equal(original.Length, result.Product.Length);
    }

    [Fact]
    public void Simulate_CircularAssemblyClosesTheRingWithoutDuplicatingTheOverlap()
    {
        string original = RandomDna(900, 41);
        var pieces = SplitWithOverlaps(original, 3, 30, circular: true);

        var result = GibsonAssembly.Simulate(pieces);

        Assert.True(result.Success, result.FailureReason);
        Assert.True(result.IsCircular);
        Assert.True(result.Product!.IsCircular);
        Assert.Equal(original.Length, result.Product.Length);
        Assert.True((original + original).Contains(result.Product.Sequence, StringComparison.Ordinal),
            "the circular product is not a rotation of the original");
        Assert.Equal(3, result.Junctions.Count);
    }

    [Fact]
    public void Simulate_OrderOfTheFragmentsDoesNotMatter()
    {
        string original = RandomDna(900, 51);
        var pieces = SplitWithOverlaps(original, 3, 30, circular: true);

        var forward = GibsonAssembly.Simulate(pieces);
        var shuffled = GibsonAssembly.Simulate([pieces[2], pieces[0], pieces[1]]);

        Assert.True(forward.Success, forward.FailureReason);
        Assert.True(shuffled.Success, shuffled.FailureReason);
        Assert.Equal(forward.Product!.Length, shuffled.Product!.Length);
        Assert.True((forward.Product.Sequence + forward.Product.Sequence)
            .Contains(shuffled.Product.Sequence, StringComparison.Ordinal));
    }

    [Fact]
    public void Simulate_ReportsInsufficientOverlapRatherThanGuessing()
    {
        string original = RandomDna(900, 61);
        // Only 8 nt of shared end, below the 20 nt minimum.
        var pieces = SplitWithOverlaps(original, 3, 8, circular: true);

        var result = GibsonAssembly.Simulate(pieces);

        Assert.False(result.Success);
        Assert.Null(result.Product);
        Assert.NotNull(result.FailureReason);
    }

    [Fact]
    public void Simulate_ReportsAmbiguityWhenTwoFragmentsShareAnEnd()
    {
        string original = RandomDna(900, 71);
        var pieces = SplitWithOverlaps(original, 3, 30, circular: true).ToList();
        // A duplicate of the second piece makes the junction ambiguous.
        pieces.Add(pieces[1]);

        var result = GibsonAssembly.Simulate(pieces);

        Assert.False(result.Success);
        Assert.NotNull(result.FailureReason);
    }

    [Fact]
    public void BestOverlap_FindsTheLongestSharedTerminalSequence()
    {
        Assert.Equal(6, GibsonAssembly.BestOverlap("AAAACCGGTT", "CCGGTTAAAA", 4));
        Assert.Equal(0, GibsonAssembly.BestOverlap("AAAACCGGTT", "GGGGGGGGGG", 4));
        // Below the minimum: not reported.
        Assert.Equal(0, GibsonAssembly.BestOverlap("AAAACCGGTT", "GTTAAAAAAA", 6));
    }

    [Fact]
    public void DesignThenSimulate_RoundTrips()
    {
        var fragments = new[] { RandomDna(600, 81), RandomDna(600, 82), RandomDna(600, 83) };
        var design = GibsonAssembly.DesignOverlaps(fragments);

        // The tailed fragments are what you would actually PCR; assembling them must give one circle whose
        // length is the sum of the inputs (every overlap counted once).
        var result = GibsonAssembly.Simulate(design.TailedFragments);

        Assert.True(result.Success, result.FailureReason);
        Assert.True(result.IsCircular);
        Assert.Equal(fragments.Sum(f => f.Length), result.Product!.Length);
    }

    // ---- Golden Gate ----

    [Fact]
    public void AddFlanks_PointsBothSitesInwardAndExposesTheRequestedOverhangs()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        string part = DomesticatedDna(300, 91, bsaI);
        Assert.True(GoldenGate.IsDomesticated(part, bsaI));

        string flanked = GoldenGate.AddFlanks(part, bsaI, "AATG", "GCTT");

        Assert.Contains(part, flanked, StringComparison.Ordinal);
        // Two sites, one on each strand, both cutting into the part.
        var sites = bsaI.Find(flanked);
        Assert.Equal(2, sites.Count);
        Assert.Contains(sites, s => s.Strand == Core.Strand.Forward);
        Assert.Contains(sites, s => s.Strand == Core.Strand.Reverse);
        Assert.Equal(["AATG", "GCTT"], sites.Select(s => s.Overhang).OrderBy(s => s, StringComparer.Ordinal).ToArray()
            .OrderBy(s => flanked.IndexOf(s, StringComparison.Ordinal)).ToArray());

        // Digesting releases the part flanked by exactly the requested overhangs.
        var released = GoldenGate.Digest(flanked, bsaI);
        var fragment = Assert.Single(released, f => f.Sequence.Contains(part, StringComparison.Ordinal));
        Assert.Equal("AATG", fragment.FivePrimeOverhang);
        Assert.Equal("GCTT", fragment.ThreePrimeOverhang);
    }

    [Fact]
    public void AddFlanks_RefusesAnUndomesticatedPartOrAWrongOverhangLength()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        string internalSite = "AAAA" + "GGTCTC" + "AAAAAAAAAA";

        Assert.False(GoldenGate.IsDomesticated(internalSite, bsaI));
        Assert.Throws<ArgumentException>(() => GoldenGate.AddFlanks(internalSite, bsaI, "AATG", "GCTT"));
        Assert.Throws<ArgumentException>(() => GoldenGate.AddFlanks(DomesticatedDna(100, 1, bsaI), bsaI, "AAT", "GCTT"));
        // A type IIP enzyme has no usable geometry here.
        Assert.Throws<ArgumentException>(() =>
            GoldenGate.AddFlanks(DomesticatedDna(100, 1, bsaI), RestrictionEnzyme.FromName("EcoRI"), "AATG", "GCTT"));
    }

    [Fact]
    public void IsDomesticated_LooksAtBothStrands()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        Assert.False(GoldenGate.IsDomesticated("AAAGGTCTCAAA", bsaI));
        Assert.False(GoldenGate.IsDomesticated("AAA" + DnaUtility.ReverseComplement("GGTCTC") + "AAA", bsaI));
        Assert.True(GoldenGate.IsDomesticated("AAAGGTCTAAAA", bsaI));
    }

    [Fact]
    public void DesignOverhangs_AreOrthogonalAndDeterministic()
    {
        var set = GoldenGate.DesignOverhangs(6);

        Assert.Equal(6, set.Count);
        Assert.Equal(6, set.Distinct(StringComparer.Ordinal).Count());
        Assert.All(set, o => Assert.Equal(4, o.Length));
        Assert.True(GoldenGate.AreOrthogonal(set, out string? reason), reason);

        // No overhang is palindromic or the reverse complement of another.
        foreach (string a in set)
        {
            Assert.NotEqual(a, DnaUtility.ReverseComplement(a));
            foreach (string b in set)
                if (!ReferenceEquals(a, b)) Assert.NotEqual(DnaUtility.ReverseComplement(a), b);
        }

        // Deterministic for a given seed, and different for another.
        Assert.Equal(set, GoldenGate.DesignOverhangs(6));
        Assert.NotEqual(set, GoldenGate.DesignOverhangs(6, seed: 7));

        // Reserved overhangs are avoided.
        var avoiding = GoldenGate.DesignOverhangs(4, reserved: set);
        Assert.Empty(avoiding.Intersect(set, StringComparer.Ordinal));
    }

    [Fact]
    public void DesignOverhangs_ValidatesItsArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => GoldenGate.DesignOverhangs(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => GoldenGate.DesignOverhangs(4, overhangLength: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => GoldenGate.DesignOverhangs(4, overhangLength: 9));
    }

    [Fact]
    public void AreOrthogonal_RejectsPalindromesDuplicatesAndReverseComplements()
    {
        Assert.False(GoldenGate.AreOrthogonal(["AATT", "GCTA"], out string? palindrome));
        Assert.NotNull(palindrome);

        Assert.False(GoldenGate.AreOrthogonal(["AATG", "AATG"], out string? duplicate));
        Assert.NotNull(duplicate);

        Assert.False(GoldenGate.AreOrthogonal(["AATG", DnaUtility.ReverseComplement("AATG")], out string? mirror));
        Assert.NotNull(mirror);
    }

    [Fact]
    public void OverhangFidelity_IsOneForAnOrthogonalSetAndLowerOtherwise()
    {
        var orthogonal = GoldenGate.DesignOverhangs(4);
        Assert.Equal(1.0, GoldenGate.OverhangFidelity(orthogonal), 9);

        // A set containing an overhang and its own reverse complement cannot ligate faithfully.
        var mirrored = new List<string>(orthogonal) { DnaUtility.ReverseComplement(orthogonal[0]) };
        Assert.True(GoldenGate.OverhangFidelity(mirrored) < 1.0);

        // Two overhangs one substitution apart are worse than two that differ everywhere.
        Assert.True(GoldenGate.OverhangFidelity(["AATG", "AATC"]) < GoldenGate.OverhangFidelity(["AATG", "GCCA"]));
    }

    [Fact]
    public void Simulate_AssemblesThreePartsIntoABackboneInOrder()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        var overhangs = GoldenGate.DesignOverhangs(4);

        var parts = new[] { DomesticatedDna(300, 101, bsaI), DomesticatedDna(300, 102, bsaI), DomesticatedDna(300, 103, bsaI) };
        var flanked = new List<string>();
        for (int i = 0; i < parts.Length; i++)
            flanked.Add(GoldenGate.AddFlanks(parts[i], bsaI, overhangs[i], overhangs[i + 1]));

        // The backbone closes the circle: it receives the last overhang and returns the first.
        string backboneBody = DomesticatedDna(1200, 104, bsaI);
        string backbone = GoldenGate.AddFlanks(backboneBody, bsaI, overhangs[3], overhangs[0]);

        var result = GoldenGate.Simulate(flanked, bsaI, backbone, backboneIsCircular: false);

        Assert.True(result.Success, result.FailureReason);
        Assert.True(result.Product!.IsCircular);
        Assert.Equal(4, result.JunctionOverhangs.Count);

        // Every part is present exactly once, in order, and no BsaI site survives.
        string doubled = result.Product.Sequence + result.Product.Sequence;
        foreach (string part in parts)
        {
            Assert.Contains(part, doubled, StringComparison.Ordinal);
            Assert.Equal(2, CountOccurrences(doubled, part));
        }
        Assert.Contains(backboneBody, doubled, StringComparison.Ordinal);
        Assert.Empty(bsaI.Find(result.Product.Sequence, circular: true));

        int indexOfFirst = doubled.IndexOf(parts[0], StringComparison.Ordinal);
        int indexOfSecond = doubled.IndexOf(parts[1], StringComparison.Ordinal);
        int indexOfThird = doubled.IndexOf(parts[2], StringComparison.Ordinal);
        Assert.True(indexOfFirst < indexOfSecond && indexOfSecond < indexOfThird,
            "the parts were not assembled in the order their overhangs imply");
    }

    [Fact]
    public void Simulate_WorksWithBsmBiToo()
    {
        var bsmBI = RestrictionEnzyme.FromName("BsmBI");
        var overhangs = GoldenGate.DesignOverhangs(3, seed: 5);

        var parts = new[] { DomesticatedDna(200, 111, bsmBI), DomesticatedDna(200, 112, bsmBI) };
        var flanked = new[]
        {
            GoldenGate.AddFlanks(parts[0], bsmBI, overhangs[0], overhangs[1]),
            GoldenGate.AddFlanks(parts[1], bsmBI, overhangs[1], overhangs[2]),
        };
        string backbone = GoldenGate.AddFlanks(DomesticatedDna(800, 113, bsmBI), bsmBI, overhangs[2], overhangs[0]);

        var result = GoldenGate.Simulate(flanked, bsmBI, backbone, backboneIsCircular: false);

        Assert.True(result.Success, result.FailureReason);
        Assert.Empty(bsmBI.Find(result.Product!.Sequence, circular: true));
    }

    [Fact]
    public void Simulate_FailsOnAmbiguousOverhangsAndOnPartsWithoutSites()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        var overhangs = GoldenGate.DesignOverhangs(3);

        // Two parts offering the same 5' overhang.
        var ambiguous = new[]
        {
            GoldenGate.AddFlanks(DomesticatedDna(200, 121, bsaI), bsaI, overhangs[0], overhangs[1]),
            GoldenGate.AddFlanks(DomesticatedDna(200, 122, bsaI), bsaI, overhangs[0], overhangs[2]),
        };
        var clash = GoldenGate.Simulate(ambiguous, bsaI);
        Assert.False(clash.Success);
        Assert.NotNull(clash.FailureReason);

        // A part with no recognition site at all releases nothing.
        var noSite = GoldenGate.Simulate([DomesticatedDna(200, 123, bsaI)], bsaI);
        Assert.False(noSite.Success);
        Assert.NotNull(noSite.FailureReason);

        // An empty part is rejected rather than silently skipped.
        Assert.False(GoldenGate.Simulate([""], bsaI).Success);
    }

    [Fact]
    public void Simulate_FailsWhenTheChainCannotClose()
    {
        var bsaI = RestrictionEnzyme.FromName("BsaI");
        var overhangs = GoldenGate.DesignOverhangs(4);

        // Two parts whose overhangs form an open chain: nothing returns to the start.
        var open = new[]
        {
            GoldenGate.AddFlanks(DomesticatedDna(200, 131, bsaI), bsaI, overhangs[0], overhangs[1]),
            GoldenGate.AddFlanks(DomesticatedDna(200, 132, bsaI), bsaI, overhangs[1], overhangs[2]),
        };

        var result = GoldenGate.Simulate(open, bsaI);
        Assert.False(result.Success);
        Assert.NotNull(result.FailureReason);
    }

    private static int CountOccurrences(string haystack, string needle)
    {
        int count = 0, index = 0;
        while ((index = haystack.IndexOf(needle, index, StringComparison.Ordinal)) >= 0) { count++; index++; }
        return count;
    }
}
