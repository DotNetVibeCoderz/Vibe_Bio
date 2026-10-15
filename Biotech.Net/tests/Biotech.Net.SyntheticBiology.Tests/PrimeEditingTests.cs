using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology;
using Biotech.Net.SyntheticBiology.Crispr;
using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// pegRNA design for prime editing. The tests assert the published geometry directly: the nick sits between
/// protospacer positions 17 and 18, the PBS anneals to the nicked 3′ flap, and the RTT is templated from the
/// <em>edited</em> strand, so reverse-complementing it has to reproduce the allele being installed.
/// </summary>
public class PrimeEditingTests
{
    private static string RandomTarget(int length, int seed)
    {
        var random = new Random(seed);
        return string.Concat(Enumerable.Range(0, length).Select(_ => "ACGT"[random.Next(4)]));
    }

    private static readonly string Target = RandomTarget(300, 20260114);

    private static PrimeEditorOptions ForwardOnly => new() { SearchReverse = false };

    // ---- the edit model ----

    [Fact]
    public void Substitution_ReplacesInPlace()
    {
        const string reference = "AAACCCGGGTTT";
        var edit = PrimeEdit.Substitution(reference, 4, "TT");

        Assert.Equal(PrimeEditKind.Substitution, edit.Kind);
        Assert.Equal(4, edit.Position);
        Assert.Equal("CC", edit.ReferenceAllele);
        Assert.Equal("TT", edit.AlternateAllele);
        Assert.Equal(0, edit.LengthDelta);
        Assert.Equal("AAACTTGGGTTT", edit.Apply(reference));
    }

    [Fact]
    public void Insertion_AddsBasesWithoutConsumingAny()
    {
        const string reference = "AAACCCGGGTTT";
        var edit = PrimeEdit.Insertion(6, "GAT");

        Assert.Equal(PrimeEditKind.Insertion, edit.Kind);
        Assert.Equal("", edit.ReferenceAllele);
        Assert.Equal("GAT", edit.AlternateAllele);
        Assert.Equal(3, edit.LengthDelta);
        Assert.Equal("AAACCCGATGGGTTT", edit.Apply(reference));
    }

    [Fact]
    public void Deletion_RemovesTheReferenceAllele()
    {
        const string reference = "AAACCCGGGTTT";
        var edit = PrimeEdit.Deletion(reference, 3, 3);

        Assert.Equal(PrimeEditKind.Deletion, edit.Kind);
        Assert.Equal("CCC", edit.ReferenceAllele);
        Assert.Equal("", edit.AlternateAllele);
        Assert.Equal(-3, edit.LengthDelta);
        Assert.Equal("AAAGGGTTT", edit.Apply(reference));
    }

    // ---- pegRNA geometry ----

    /// <summary>
    /// Checks every documented relation between a pegRNA and its target, for a forward-strand design where the
    /// nicked strand is the plus strand and the coordinates need no mirroring.
    /// </summary>
    private static void AssertForwardGeometry(string reference, string edited, PegRna peg, PrimeEditorOptions options)
    {
        Assert.Equal(Strand.Forward, peg.Strand);

        // Spacer and PAM come straight off the plus strand.
        Assert.Equal(reference.Substring(peg.ProtospacerStart, 20), peg.Spacer);
        Assert.Equal(reference.Substring(peg.ProtospacerStart + 20, 3), peg.Pam);
        Assert.Equal("GG", peg.Pam[1..]);

        // SpCas9 nicks 3 bp from the PAM: between protospacer positions 17 and 18.
        Assert.Equal(peg.ProtospacerStart + 17, peg.NickPosition);

        // The PBS anneals to the liberated 3' end, so it is the reverse complement of the bases ending at the nick.
        Assert.Equal(DnaUtility.ReverseComplement(reference[(peg.NickPosition - peg.PbsLength)..peg.NickPosition]), peg.Pbs);

        // The RTT templates the new strand: reverse-complementing it gives the EDITED sequence from the nick on.
        Assert.Equal(edited.Substring(peg.NickPosition, peg.RttLength), peg.EditedFlap);
        Assert.Equal(DnaUtility.ReverseComplement(peg.Rtt), peg.EditedFlap);

        // Oligo layout and the reported distances.
        Assert.Equal(peg.Rtt + peg.Pbs, peg.ThreePrimeExtension);
        Assert.Equal(peg.Spacer + peg.Scaffold + peg.Rtt + peg.Pbs, peg.FullSequence);
        Assert.Equal(peg.Edit.Position - peg.NickPosition, peg.NickToEditDistance);

        // Option ranges.
        Assert.InRange(peg.PbsLength, options.PbsMinLength, options.PbsMaxLength);
        Assert.InRange(peg.RttLength, options.RttMinLength, options.RttMaxLength);
        Assert.InRange(peg.NickToEditDistance, 0, options.MaxNickToEditDistance);
        Assert.True(peg.HomologyAfterEdit >= options.MinHomologyAfterEdit,
            $"RTT reaches only {peg.HomologyAfterEdit} nt past the edit");
    }

    [Fact]
    public void Design_Substitution_SatisfiesTheFullGeometry()
    {
        var options = ForwardOnly;
        var edit = PrimeEdit.Substitution(Target, 150, "A" == Target[150].ToString() ? "C" : "A");
        string edited = edit.Apply(Target);

        var pegs = PrimeEditor.Design(Target, edit, options);

        Assert.NotEmpty(pegs);
        foreach (var peg in pegs) AssertForwardGeometry(Target, edited, peg, options);

        // Ranked best first.
        Assert.Equal(pegs.Select(p => p.Score).OrderByDescending(s => s), pegs.Select(p => p.Score));
        Assert.True(pegs.Count <= options.MaxResults);
    }

    [Fact]
    public void Design_InsertionAndDeletion_RoundTripThroughTheFlap()
    {
        var options = ForwardOnly;

        var insertion = PrimeEdit.Insertion(150, "GATTACA");
        string insertedSequence = insertion.Apply(Target);
        var insertionPegs = PrimeEditor.Design(Target, insertion, options);
        Assert.NotEmpty(insertionPegs);
        foreach (var peg in insertionPegs) AssertForwardGeometry(Target, insertedSequence, peg, options);
        // The flap really does carry the inserted bases.
        Assert.All(insertionPegs, p => Assert.Contains("GATTACA", p.EditedFlap, StringComparison.Ordinal));

        var deletion = PrimeEdit.Deletion(Target, 150, 6);
        string deletedSequence = deletion.Apply(Target);
        var deletionPegs = PrimeEditor.Design(Target, deletion, options);
        Assert.NotEmpty(deletionPegs);
        foreach (var peg in deletionPegs) AssertForwardGeometry(Target, deletedSequence, peg, options);
    }

    [Fact]
    public void Design_PbsAndRttScansStayInsideTheConfiguredRanges()
    {
        var options = new PrimeEditorOptions
        {
            SearchReverse = false,
            PbsMinLength = 10,
            PbsMaxLength = 12,
            RttMinLength = 12,
            RttMaxLength = 20,
            MinHomologyAfterEdit = 6,
        };
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        string edited = edit.Apply(Target);

        var pegs = PrimeEditor.Design(Target, edit, options);

        Assert.NotEmpty(pegs);
        foreach (var peg in pegs)
        {
            AssertForwardGeometry(Target, edited, peg, options);
            Assert.InRange(peg.PbsLength, 10, 12);
            Assert.InRange(peg.RttLength, 12, 20);
            Assert.True(peg.HomologyAfterEdit >= 6);
        }
    }

    [Fact]
    public void Design_PrefersAPbsNearTheTargetMeltingTemperature()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        var pegs = PrimeEditor.Design(Target, edit, ForwardOnly);

        var best = pegs[0];
        Assert.Equal(NearestNeighbourTm.Tm(best.Pbs), best.PbsMeltingTemperature, 9);
        // 38 degrees is the published preference; the scan should land within a few degrees of it.
        Assert.InRange(best.PbsMeltingTemperature, 28.0, 48.0);
    }

    [Fact]
    public void Design_HonoursTheNickToEditDistanceCap()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");

        var tight = PrimeEditor.Design(Target, edit, new PrimeEditorOptions { SearchReverse = false, MaxNickToEditDistance = 10 });
        Assert.All(tight, p => Assert.InRange(p.NickToEditDistance, 0, 10));

        var loose = PrimeEditor.Design(Target, edit, new PrimeEditorOptions { SearchReverse = false, MaxNickToEditDistance = 30 });
        Assert.True(loose.Count >= tight.Count);
    }

    [Fact]
    public void Design_FindsNothingWhenNoPamIsInRange()
    {
        // A homopolymer carries no NGG PAM at all.
        string flat = new('A', 200);
        var edit = PrimeEdit.Substitution(flat, 100, "C");
        Assert.Empty(PrimeEditor.Design(flat, edit, ForwardOnly));
    }

    [Fact]
    public void Design_ReverseStrandPegRnasUseTheMirroredGeometry()
    {
        // Designing against the reverse complement of the target must give mirror-image pegRNAs: the same
        // molecules, found from the other strand.
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        var forward = PrimeEditor.Design(Target, edit, ForwardOnly);

        string flipped = DnaUtility.ReverseComplement(Target);
        int flippedPosition = Target.Length - 151;
        var flippedEdit = PrimeEdit.Substitution(flipped, flippedPosition,
            DnaUtility.ReverseComplement(edit.AlternateAllele));
        var reverse = PrimeEditor.Design(flipped, flippedEdit, new PrimeEditorOptions { SearchForward = false });

        Assert.NotEmpty(reverse);
        Assert.All(reverse, p => Assert.Equal(Strand.Reverse, p.Strand));
        // The same spacers, as sets.
        Assert.Equal(forward.Select(p => p.Spacer).OrderBy(s => s, StringComparer.Ordinal),
                     reverse.Select(p => p.Spacer).OrderBy(s => s, StringComparer.Ordinal));
        // And the same designed extensions, since a pegRNA is a molecule, not a coordinate system.
        Assert.Equal(forward.Select(p => p.ThreePrimeExtension).OrderBy(s => s, StringComparer.Ordinal),
                     reverse.Select(p => p.ThreePrimeExtension).OrderBy(s => s, StringComparer.Ordinal));
    }

    [Fact]
    public void PegRna_ReportsThePreferredRttEnd()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        foreach (var peg in PrimeEditor.Design(Target, edit, ForwardOnly))
        {
            Assert.Equal(peg.Rtt[0], peg.RttTerminalBase);
            Assert.Equal(peg.RttTerminalBase != 'C', peg.HasPreferredRttEnd);
        }
    }

    // ---- PE3 / PE3b nicking guides ----

    [Fact]
    public void NickingGuides_AreOnTheOppositeStrandAndInsideTheDistanceWindow()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        var peg = PrimeEditor.Design(Target, edit, ForwardOnly)[0];

        var options = new PrimeEditorOptions { NickingDistanceMin = 40, NickingDistanceMax = 90 };
        var nickers = PrimeEditor.DesignNickingGuides(Target, peg, options);

        Assert.NotEmpty(nickers);
        foreach (var nicker in nickers)
        {
            Assert.NotEqual(peg.Strand, nicker.Guide.Strand);
            Assert.Equal(nicker.Guide.CutSite, nicker.NickPosition);
            Assert.Equal(Math.Abs(nicker.NickPosition - peg.NickPosition), nicker.DistanceFromPegNick);
            Assert.InRange(nicker.DistanceFromPegNick, 40, 90);
        }
        Assert.Equal(nickers.Select(n => n.Score).OrderByDescending(s => s), nickers.Select(n => n.Score));
    }

    [Fact]
    public void NickingGuides_NarrowingTheWindowNarrowsTheResult()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        var peg = PrimeEditor.Design(Target, edit, ForwardOnly)[0];

        var wide = PrimeEditor.DesignNickingGuides(Target, peg,
            new PrimeEditorOptions { NickingDistanceMin = 20, NickingDistanceMax = 120 });
        var narrow = PrimeEditor.DesignNickingGuides(Target, peg,
            new PrimeEditorOptions { NickingDistanceMin = 55, NickingDistanceMax = 70 });

        Assert.True(wide.Count >= narrow.Count);
        Assert.All(narrow, n => Assert.InRange(n.DistanceFromPegNick, 55, 70));
    }

    [Fact]
    public void Pe3BGuidesAreFlaggedExactlyWhenTheProtospacerOverlapsTheEdit()
    {
        var edit = PrimeEdit.Substitution(Target, 150, Target[150] == 'A' ? "C" : "A");
        var peg = PrimeEditor.Design(Target, edit, ForwardOnly)[0];

        // A wide window so that guides overlapping the edit are included.
        var nickers = PrimeEditor.DesignNickingGuides(Target, peg,
            new PrimeEditorOptions { NickingDistanceMin = 0, NickingDistanceMax = 200 });

        foreach (var nicker in nickers)
        {
            bool overlaps = nicker.Guide.Start < Math.Max(edit.ReferenceEnd, edit.Position + 1)
                            && nicker.Guide.End > edit.Position;
            Assert.Equal(overlaps, nicker.IsPe3B);
        }
        Assert.Contains(nickers, n => n.IsPe3B);
        // PE3b guides score above an otherwise comparable PE3 guide, so they rank first when present.
        Assert.True(nickers[0].Score >= nickers[^1].Score);
    }
}
