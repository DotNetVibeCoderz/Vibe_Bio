using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology;
using Biotech.Net.SyntheticBiology.Crispr;

namespace Biotech.Net.SyntheticBiology.Tests;

/// <summary>
/// Cytosine and adenine base-editing outcome prediction. Protospacer positions are 1-based and counted from
/// the PAM-distal end throughout, which is the convention every editing-window paper uses and the thing most
/// worth pinning down with tests.
/// </summary>
public class BaseEditingTests
{
    // A single C at protospacer position 6, inside the BE4max window (4-8).
    private const string OneCytosineInWindow = "AAGAACAGTAGTAGTAGTAG";

    [Fact]
    public void Editor_WindowConventionIsOneBasedFromThePamDistalEnd()
    {
        var editor = BaseEditor.Be4Max;
        Assert.Equal(4, editor.EditingWindowStart);
        Assert.Equal(8, editor.EditingWindowEnd);
        Assert.Equal(5, editor.WindowLength);
        Assert.Equal(6.0, editor.WindowCentre, 12);

        Assert.False(editor.IsInWindow(3));
        Assert.True(editor.IsInWindow(4));
        Assert.True(editor.IsInWindow(8));
        Assert.False(editor.IsInWindow(9));

        Assert.Equal(0.0, editor.WeightAt(20), 12);
        Assert.True(editor.WeightAt(6) > 0);
    }

    [Fact]
    public void Editor_BasesFollowTheChemistry()
    {
        Assert.Equal('C', BaseEditor.Be4Max.TargetBase);
        Assert.Equal('T', BaseEditor.Be4Max.ProductBase);
        Assert.Equal('G', BaseEditor.Be4Max.ComplementaryTargetBase);
        Assert.Equal('A', BaseEditor.Be4Max.ComplementaryProductBase);

        Assert.Equal('A', BaseEditor.Abe710.TargetBase);
        Assert.Equal('G', BaseEditor.Abe710.ProductBase);
        Assert.Equal('T', BaseEditor.Abe710.ComplementaryTargetBase);
        Assert.Equal('C', BaseEditor.Abe710.ComplementaryProductBase);
    }

    [Fact]
    public void Catalogue_IsLookupConsistent()
    {
        foreach (var editor in BaseEditor.All)
        {
            Assert.Same(editor, BaseEditor.FromName(editor.Name));
            Assert.True(editor.EditingWindowStart <= editor.EditingWindowEnd);
            Assert.NotEqual("", editor.Reference);
        }
        Assert.Equal(BaseEditorKind.Cytosine, BaseEditor.Be3.Kind);
        Assert.Equal(BaseEditorKind.Adenine, BaseEditor.Abe8E.Kind);
    }

    [Fact]
    public void Cbe_EditsASingleCytosineInsideTheWindow()
    {
        var outcome = BaseEditing.PredictEdits(BaseEditor.Be4Max, OneCytosineInWindow);

        var site = Assert.Single(outcome.Sites);
        Assert.Equal(6, site.Position);
        Assert.Equal('C', site.From);
        Assert.Equal('T', site.To);
        Assert.True(site.IsIntended);

        Assert.True(outcome.IsUnambiguous);
        Assert.Empty(outcome.Bystanders);
        Assert.Equal("AAGAATAGTAGTAGTAGTAG", outcome.EditedProtospacer);
        Assert.Equal(outcome.EditedProtospacer, outcome.IntendedOnlyProtospacer);
    }

    [Fact]
    public void Cbe_IgnoresCytosinesOutsideTheWindow()
    {
        // C at position 15 only: outside the 4-8 window.
        const string protospacer = "AAGAAGAGTAGTAGCAGTAG";
        Assert.Equal('C', protospacer[14]);

        var outcome = BaseEditing.PredictEdits(BaseEditor.Be4Max, protospacer);

        Assert.False(outcome.HasEdit);
        Assert.Empty(outcome.Sites);
        Assert.Null(outcome.Intended);
        Assert.False(outcome.IsUnambiguous);
        Assert.Equal(protospacer, outcome.EditedProtospacer);
    }

    [Fact]
    public void Cbe_ReportsEveryBystanderInTheWindow()
    {
        // C at positions 4, 6 and 8; the caller asks for position 6.
        const string protospacer = "AAGCACACTAGTAGTAGTAG";
        Assert.Equal('C', protospacer[3]);
        Assert.Equal('C', protospacer[5]);
        Assert.Equal('C', protospacer[7]);

        var outcome = BaseEditing.PredictEdits(BaseEditor.Be4Max, protospacer, intendedPosition: 6);

        Assert.Equal([4, 6, 8], outcome.Sites.Select(s => s.Position));
        Assert.False(outcome.IsUnambiguous);
        Assert.Equal(6, outcome.Intended!.Position);
        Assert.Equal([4, 8], outcome.Bystanders.Select(s => s.Position));

        // Fully converted product versus the intended-only product.
        Assert.Equal("AAGTATATTAGTAGTAGTAG", outcome.EditedProtospacer);
        Assert.Equal("AAGCATACTAGTAGTAGTAG", outcome.IntendedOnlyProtospacer);
    }

    [Fact]
    public void Abe_EditsAdenineAndLeavesCytosineAlone()
    {
        // A at position 5 only, inside the ABE7.10 window (4-7); several C's elsewhere.
        const string protospacer = "GCGCACGCGCGCGCGCGCGC";
        Assert.Equal('A', protospacer[4]);

        var outcome = BaseEditing.PredictEdits(BaseEditor.Abe710, protospacer);

        var site = Assert.Single(outcome.Sites);
        Assert.Equal(5, site.Position);
        Assert.Equal('A', site.From);
        Assert.Equal('G', site.To);
        Assert.Equal("GCGCGCGCGCGCGCGCGCGC", outcome.EditedProtospacer);

        // The same protospacer is untouched by a cytosine editor only because its window holds no C... it does,
        // so a CBE sees several sites: the two editors are not interchangeable.
        Assert.NotEmpty(BaseEditing.PredictEdits(BaseEditor.Be4Max, protospacer).Sites);
    }

    [Fact]
    public void NarrowedWindowEditorSeesFewerSites()
    {
        const string protospacer = "AAGCACAGTAGTAGTAGTAG"; // C at 4 and 6
        var wide = BaseEditing.PredictEdits(BaseEditor.Be4Max, protospacer);     // window 4-8
        var narrow = BaseEditing.PredictEdits(BaseEditor.Ye1, protospacer);      // narrowed window

        Assert.Equal(2, wide.Sites.Count);
        Assert.True(narrow.Sites.Count <= wide.Sites.Count);
        Assert.All(narrow.Sites, s => Assert.True(BaseEditor.Ye1.IsInWindow(s.Position)));
    }

    // ---- coding consequences ----
    //
    // cds index:  0 1 2 3 4 5 6 7 8 ...
    //             A T G G G G C A A ...
    // The protospacer starts at index 3, so protospacer position 4 is cds index 6 — the first base of codon 2.

    [Theory]
    [InlineData("ATGGGGCAAAAGTGTGTGTGTGTA", "CAA", "TAA", 'Q', '*', EditConsequence.Nonsense)]
    [InlineData("ATGGGGCTGAAGTGTGTGTGTGTA", "CTG", "TTG", 'L', 'L', EditConsequence.Silent)]
    [InlineData("ATGGGGCGTAAGTGTGTGTGTGTA", "CGT", "TGT", 'R', 'C', EditConsequence.Missense)]
    public void Cbe_ProteinConsequenceIsClassified(string cds, string reference, string edited,
        char referenceAa, char editedAa, EditConsequence consequence)
    {
        var outcome = BaseEditing.PredictProteinConsequence(BaseEditor.Be4Max, cds, protospacerStartInCds: 3);

        var intended = Assert.Single(outcome.Consequences);
        Assert.Same(intended, outcome.IntendedConsequence);
        Assert.Equal(4, intended.Site.Position);
        Assert.Equal(6, intended.CdsPosition);
        Assert.Equal(2, intended.CodonIndex);
        Assert.Equal(0, intended.PositionInCodon);
        Assert.Equal(reference, intended.ReferenceCodon);
        Assert.Equal(edited, intended.EditedCodon);
        Assert.Equal(referenceAa, intended.ReferenceAminoAcid);
        Assert.Equal(editedAa, intended.EditedAminoAcid);
        Assert.Equal(consequence, intended.Consequence);

        // The combined product carries the edit and nothing else.
        Assert.Equal(cds.Length, outcome.EditedCds.Length);
        Assert.Equal(string.Concat(cds.AsSpan(0, 6), edited, cds.AsSpan(9)), outcome.EditedCds);
    }

    [Fact]
    public void Cbe_OnTheReverseStrandReadsAsGtoAOnTheCodingStrand()
    {
        // The reverse complement of the nonsense case, so the same protospacer now lies on the minus strand.
        const string cds = "ATGGGGCAAAAGTGTGTGTGTGTA";
        string flipped = DnaUtility.ReverseComplement(cds);
        // The protospacer occupied [3, 23) of cds, so it occupies [1, 21) of the flipped sequence.
        var outcome = BaseEditing.PredictProteinConsequence(BaseEditor.Be4Max, flipped, 1, Strand.Reverse);

        var intended = Assert.Single(outcome.Consequences);
        // Read on the coding strand of `flipped`, a protospacer-strand C→T is a G→A.
        Assert.Equal('G', flipped[intended.CdsPosition]);
        Assert.Equal('A', outcome.EditedCds[intended.CdsPosition]);
        Assert.Equal(flipped.Length, outcome.EditedCds.Length);
        Assert.Equal(DnaUtility.ReverseComplement(outcome.EditedCds),
                     string.Concat(cds.AsSpan(0, 6), "TAA", cds.AsSpan(9)));
    }

    [Fact]
    public void ProteinConsequence_RejectsAProtospacerOutsideTheCds() =>
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BaseEditing.PredictProteinConsequence(BaseEditor.Be4Max, "ATGGGGCAAAAG", 3));

    // ---- guide selection ----

    private const string EditingTargetSpacer = "AGTAGCAGTAGTAGTAGTAG"; // C at position 6
    private const string EditingTarget = "TTTT" + EditingTargetSpacer + "TGG" + "AAA";

    [Fact]
    public void FindEditingGuides_ReturnsGuidesWhoseWindowCoversTheRequestedBase()
    {
        const int targetPosition = 4 + 5; // the C, 0-based on the plus strand
        Assert.Equal('C', EditingTarget[targetPosition]);

        var candidates = BaseEditing.FindEditingGuides(EditingTarget, BaseEditor.Be4Max, targetPosition);

        var candidate = Assert.Single(candidates);
        Assert.Equal(EditingTargetSpacer, candidate.Guide.Spacer);
        Assert.Equal(6, candidate.ProtospacerPosition);
        Assert.True(BaseEditor.Be4Max.IsInWindow(candidate.ProtospacerPosition));
        Assert.True(candidate.Outcome.IsUnambiguous);

        // The reported position really is the requested base, recomputed from the guide.
        Assert.Equal(candidate.ProtospacerPosition,
            BaseEditing.ProtospacerPosition(candidate.Guide, targetPosition));
    }

    [Fact]
    public void FindEditingGuides_IsEmptyWhenNoWindowReachesThePosition()
    {
        // Position 0 of the target is 4 nt upstream of the only protospacer.
        Assert.Empty(BaseEditing.FindEditingGuides(EditingTarget, BaseEditor.Be4Max, 0));
        // And a position inside the protospacer but outside the window (position 15) is rejected too.
        Assert.Empty(BaseEditing.FindEditingGuides(EditingTarget, BaseEditor.Be4Max, 4 + 14));
    }

    [Fact]
    public void FindEditingGuides_RejectsAPositionOffTheSequence() =>
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BaseEditing.FindEditingGuides(EditingTarget, BaseEditor.Be4Max, EditingTarget.Length));

    [Fact]
    public void ProtospacerPosition_IsZeroOutsideTheProtospacerAndMirrorsOnTheMinusStrand()
    {
        var forward = Assert.Single(GuideFinder.Find(EditingTarget));
        Assert.Equal(0, BaseEditing.ProtospacerPosition(forward, forward.Start - 1));
        Assert.Equal(1, BaseEditing.ProtospacerPosition(forward, forward.Start));
        Assert.Equal(20, BaseEditing.ProtospacerPosition(forward, forward.End - 1));
        Assert.Equal(0, BaseEditing.ProtospacerPosition(forward, forward.End));

        string flipped = DnaUtility.ReverseComplement(EditingTarget);
        var reverse = Assert.Single(GuideFinder.Find(flipped));
        Assert.Equal(Strand.Reverse, reverse.Strand);
        // On the minus strand position 1 is the rightmost plus-strand base of the protospacer.
        Assert.Equal(1, BaseEditing.ProtospacerPosition(reverse, reverse.End - 1));
        Assert.Equal(20, BaseEditing.ProtospacerPosition(reverse, reverse.Start));
    }

    [Fact]
    public void EditingWindowInterval_IsTheWindowProjectedOntoThePlusStrand()
    {
        var forward = Assert.Single(GuideFinder.Find(EditingTarget));
        var interval = BaseEditing.EditingWindowInterval(BaseEditor.Be4Max, forward);

        // Window 4-8 of a protospacer starting at `Start`: [Start + 3, Start + 8).
        Assert.Equal(forward.Start + 3, interval.Start);
        Assert.Equal(forward.Start + 8, interval.End);
        Assert.Equal(5, interval.End - interval.Start);
        Assert.Equal(Strand.Forward, interval.Strand);

        var reverse = Assert.Single(GuideFinder.Find(DnaUtility.ReverseComplement(EditingTarget)));
        var mirrored = BaseEditing.EditingWindowInterval(BaseEditor.Be4Max, reverse);
        Assert.Equal(reverse.End - 8, mirrored.Start);
        Assert.Equal(reverse.End - 3, mirrored.End);
        Assert.Equal(5, mirrored.End - mirrored.Start);
    }
}
