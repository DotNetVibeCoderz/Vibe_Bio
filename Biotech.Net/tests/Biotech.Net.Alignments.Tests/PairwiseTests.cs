using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments.Tests;

public class SubstitutionMatrixTests
{
    [Fact]
    public void Blosum62_KnownEntries()
    {
        var m = SubstitutionMatrix.Blosum62;
        Assert.Equal(4, m['A', 'A']);
        Assert.Equal(11, m['W', 'W']);
        Assert.Equal(-4, m['*', 'A']);
        Assert.Equal(-3, m['W', 'A']);
        Assert.Equal(m['W', 'A'], m['w', 'a']); // case-insensitive
        Assert.True(m.IsSymmetric);
        Assert.Equal(-4, m.MinScore);
        Assert.Equal(11, m.MaxScore);
        Assert.Equal(m['X', 'A'], m['U', 'A']); // unknown residue falls back to X
    }

    [Fact]
    public void AllBuiltIns_Load()
    {
        foreach (string name in SubstitutionMatrix.AvailableNames)
        {
            var m = SubstitutionMatrix.Load(name);
            Assert.True(m.Symbols.Length >= 15, name);
            Assert.True(m.IsSymmetric, name);
        }
        Assert.Equal(5, SubstitutionMatrix.Load("nuc44")['A', 'A']);
        Assert.Equal(17, SubstitutionMatrix.Pam250['W', 'W']);
        Assert.Equal(5, SubstitutionMatrix.Blosum50['A', 'A']);
        Assert.Throws<KeyNotFoundException>(() => SubstitutionMatrix.Load("BLOSUM999"));
    }

    [Fact]
    public void Nuc44_MapsUracilToThymine()
    {
        var m = SubstitutionMatrix.Nuc44;
        Assert.Equal(5, m['U', 'T']);
        Assert.Equal(-4, m['A', 'C']);
        Assert.Equal(m['N', 'A'], m['X', 'A']);
    }

    [Fact]
    public void Parse_RoundTripsNcbiFormat()
    {
        string text = """
            # tiny
               A  C  G  T
            A  2 -1 -1 -1
            C -1  2 -1 -1
            G -1 -1  2 -1
            T -1 -1 -1  2
            """;
        var m = SubstitutionMatrix.Parse(text, "tiny");
        Assert.Equal("ACGT", m.Symbols);
        Assert.Equal(2, m['g', 'G']);
        var again = SubstitutionMatrix.Parse(m.ToNcbiFormat());
        Assert.Equal(m.KernelMatrix.ToArray(), again.KernelMatrix.ToArray());
        Assert.Throws<FormatException>(() => SubstitutionMatrix.Parse("   A C\nA 1 2 3\n"));
    }
}

public class PairwiseTests
{
    [Fact]
    public void Biopython_DefaultScore_GAACT_GAT_Is3()
    {
        // Biopython tutorial: PairwiseAligner() defaults (match 1, mismatch 0, gaps 0); aligner.score("GAACT", "GAT") == 3.0.
        var aligner = new PairwiseAligner();
        Assert.Equal(3, aligner.Score("GAACT", "GAT"));
        var aln = aligner.Align("GAACT", "GAT");
        Assert.Equal(3, aln.Score);
        Assert.Equal(3, aln.Identities);
        Assert.Equal("GAACT", aln.AlignedTarget.Replace("-", ""));
        Assert.Equal("GAT", aln.AlignedQuery.Replace("-", ""));
    }

    [Fact]
    public void Biopython_LocalScore_AGAACTC_GAACT_Is5()
    {
        // Biopython tutorial: local mode, aligner.score("AGAACTC", "GAACT") == 5.0.
        var aligner = new PairwiseAligner { Mode = AlignmentMode.Local };
        Assert.Equal(5, aligner.Score("AGAACTC", "GAACT"));
        var aln = aligner.Align("AGAACTC", "GAACT");
        Assert.Equal(5, aln.Score);
        Assert.Equal("GAACT", aln.AlignedTarget);
        Assert.Equal((1, 6, 0, 5), (aln.TargetStart, aln.TargetEnd, aln.QueryStart, aln.QueryEnd));
        Assert.Equal("5M", aln.Cigar());
        Assert.Equal([new[] { 1, 6 }, new[] { 0, 5 }], aln.Coordinates);
    }

    [Fact]
    public void Durbin_HEAGAWGHEE_PAWHEAE_Blosum50_Global_Is1()
    {
        // Durbin, Eddy, Krogh & Mitchison (1998), Biological Sequence Analysis, Fig. 2.5: BLOSUM50, linear gap −8 → score 1.
        var aligner = new PairwiseAligner { SubstitutionMatrix = SubstitutionMatrix.Blosum50, GapScore = -8 };
        var aln = aligner.Align(Sequence.Protein("HEAGAWGHEE"), Sequence.Protein("PAWHEAE"));
        Assert.Equal(1, aln.Score);
        Assert.Equal(1, aligner.Score(Sequence.Protein("HEAGAWGHEE"), Sequence.Protein("PAWHEAE")));
        Assert.Equal(RescoreLinear(aln, SubstitutionMatrix.Blosum50, -8), aln.Score);
    }

    [Fact]
    public void Durbin_HEAGAWGHEE_PAWHEAE_Blosum50_Local_Is28()
    {
        // Durbin et al. (1998), Fig. 2.6: Smith–Waterman gives AWGHE / AW-HE with score 28.
        var aligner = new PairwiseAligner { Mode = AlignmentMode.Local, SubstitutionMatrix = SubstitutionMatrix.Blosum50, GapScore = -8 };
        var aln = aligner.Align(Sequence.Protein("HEAGAWGHEE"), Sequence.Protein("PAWHEAE"));
        Assert.Equal(28, aln.Score);
        Assert.Equal("AWGHE", aln.AlignedTarget);
        Assert.Equal("AW-HE", aln.AlignedQuery);
        Assert.Equal("2M1D2M", aln.Cigar());
        Assert.Equal("1S2M1D2M2S", aln.Cigar(softClip: true));
        Assert.Equal("2=1D2=", aln.Cigar(extended: true));
        Assert.Equal(28, aligner.Score(Sequence.Protein("HEAGAWGHEE"), Sequence.Protein("PAWHEAE")));
    }

    private static int RescoreLinear(PairwiseAlignment aln, SubstitutionMatrix m, int gap)
    {
        int s = 0;
        for (int k = 0; k < aln.Length; k++)
            s += aln.AlignedTarget[k] == '-' || aln.AlignedQuery[k] == '-' ? gap : m[aln.AlignedTarget[k], aln.AlignedQuery[k]];
        return s;
    }

    [Fact]
    public void AffineGaps_PreferOneLongGap()
    {
        var aligner = new PairwiseAligner { MatchScore = 5, MismatchScore = -4, OpenGapScore = -10, ExtendGapScore = -1 };
        var aln = aligner.Align("AAAACCCCGGGG", "AAAAGGGG");
        Assert.Equal(40 - 13, aln.Score);
        Assert.Equal("AAAACCCCGGGG", aln.AlignedTarget);
        Assert.Equal("AAAA----GGGG", aln.AlignedQuery);
        Assert.Equal("4M4D4M", aln.Cigar());
        Assert.Equal(1, aln.GapOpenings);
        Assert.Equal(4, aln.Gaps);
        Assert.Equal(8.0 / 12, aln.Identity, 12);
        Assert.Equal([new[] { 0, 4, 8, 12 }, new[] { 0, 4, 4, 8 }], aln.Coordinates);
        Assert.Equal(2, aln.AlignedBlocks.Count);
    }

    [Fact]
    public void SemiGlobal_DoesNotPenaliseEndGaps()
    {
        var aligner = new PairwiseAligner { Mode = AlignmentMode.SemiGlobal, MatchScore = 1, MismatchScore = -2, OpenGapScore = -3, ExtendGapScore = -1 };
        var aln = aligner.Align("ACGTACGT", "GTAC");
        Assert.Equal(4, aln.Score);
        Assert.Equal("ACGTACGT", aln.AlignedTarget);
        Assert.Equal("--GTAC--", aln.AlignedQuery);
        var global = new PairwiseAligner { MatchScore = 1, MismatchScore = -2, OpenGapScore = -3, ExtendGapScore = -1 };
        Assert.Equal(4 - 4 - 4, global.Score("ACGTACGT", "GTAC"));
    }

    [Fact]
    public void Banded_MatchesFullForSimilarSequences()
    {
        var full = new PairwiseAligner { SubstitutionMatrix = SubstitutionMatrix.Nuc44, OpenGapScore = -10, ExtendGapScore = -1 };
        var banded = new PairwiseAligner { SubstitutionMatrix = SubstitutionMatrix.Nuc44, OpenGapScore = -10, ExtendGapScore = -1, Band = 5 };
        var a = Sequence.Dna("ACGTTGCAAGGCTTACGATCGATCGGATCCATGCA");
        var b = Sequence.Dna("ACGTTGCAAGCTTACGATCGATTCGGATCCATGCA");
        Assert.Equal(full.Score(a, b), banded.Score(a, b));
    }

    [Fact]
    public void Format_IsBiopythonStyle()
    {
        var aligner = new PairwiseAligner { MatchScore = 5, MismatchScore = -4, OpenGapScore = -10, ExtendGapScore = -1 };
        var aln = aligner.Align("AAAACCCCGGGG", "AAAAGGGG");
        string expected =
            "target            0 AAAACCCCGGGG 12" + Environment.NewLine +
            "                  0 ||||----|||| 12" + Environment.NewLine +
            "query             0 AAAA----GGGG 8" + Environment.NewLine;
        Assert.Equal(expected, aln.Format());
        // Wrapping keeps running coordinates.
        string wrapped = aln.Format(width: 5);
        Assert.Contains("target            5 CCCGG 10", wrapped);
        Assert.Contains("query             4 ---GG 6", wrapped);
    }

    [Fact]
    public void Mismatch_ShowsDotInFormat_AndCountsPositives()
    {
        var aligner = new PairwiseAligner { SubstitutionMatrix = SubstitutionMatrix.Blosum62, OpenGapScore = -11, ExtendGapScore = -1 };
        var aln = aligner.Align(Sequence.Protein("MKVL", "p1"), Sequence.Protein("MRVL", "p2"));
        Assert.Equal(5 + 2 + 4 + 4, aln.Score); // M/M 5, K/R 2, V/V 4, L/L 4
        Assert.Equal(1, aln.Mismatches);
        Assert.Equal(4, aln.Positives);
        Assert.Contains("|.||", aln.Format());
        Assert.StartsWith("p1 ", aln.Format());
    }

    [Fact]
    public void GapScores_MustBeNonPositive()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new PairwiseAligner { OpenGapScore = 1 });
    }

    [Fact]
    public void ManagedAndNativeBackends_GiveSameAlignment()
    {
        var aligner = new PairwiseAligner { Mode = AlignmentMode.Local, SubstitutionMatrix = SubstitutionMatrix.Blosum62, OpenGapScore = -11, ExtendGapScore = -1 };
        var a = Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAAG");
        var b = Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQ");
        var first = aligner.Align(a, b);
        bool old = BioKernels.ForceManaged;
        try
        {
            BioKernels.ForceManaged = !old;
            var second = aligner.Align(a, b);
            Assert.Equal(first.Score, second.Score);
            Assert.Equal(first.Operations, second.Operations);
        }
        finally { BioKernels.ForceManaged = old; }
    }
}
