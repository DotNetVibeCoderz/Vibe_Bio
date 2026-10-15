using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Alignments.Search;
using Biotech.Net.Core;

namespace Biotech.Net.Alignments.Tests;

public class MultipleAlignmentTests
{
    private static readonly Sequence[] Globins =
    [
        Sequence.Protein("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL", "HBA_HUMAN"),
        Sequence.Protein("VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKL", "HBB_HUMAN"),
        Sequence.Protein("VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKH", "MYG_PHYCA"),
        Sequence.Protein("VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHGKKVGDALTLAVGHLDDLPGALSNLSDLHAHKL", "HBA_HORSE"),
    ];

    [Fact]
    public void Progressive_PreservesResidues_AndInputOrder()
    {
        foreach (var dist in new[] { GuideDistance.Kmer, GuideDistance.Alignment })
            foreach (var method in new[] { GuideTreeMethod.Upgma, GuideTreeMethod.NeighborJoining })
            {
                var msa = new ProgressiveAligner { Distance = dist, TreeMethod = method }.Align(Globins);
                Assert.Equal(Globins.Select(g => g.Id), msa.Ids);
                for (int i = 0; i < Globins.Length; i++) Assert.Equal(Globins[i].ToString(), msa[i].Replace("-", ""));
                Assert.All(msa.Rows, r => Assert.Equal(msa.Length, r.Length));
                Assert.NotNull(msa.GuideTree);
                Assert.Equal(4, msa.GuideTree!.Leaves.Count);
                // The two alpha chains are by far the closest pair.
                Assert.True(msa.PercentIdentity(0, 3) > 80);
                // The famous conserved "HGKKV" motif of the haemoglobins lines up.
                int col = msa[0].IndexOf("HGKKV", StringComparison.Ordinal);
                Assert.True(col >= 0);
                Assert.Equal("HGKKV", msa[1].Substring(col, 5));
                Assert.Equal("HGKKV", msa[3].Substring(col, 5));
            }
    }

    [Fact]
    public void Progressive_NucleotideInsertion_GetsOneGapColumn()
    {
        var seqs = new[] { Sequence.Dna("ACGTACGTAC", "a"), Sequence.Dna("ACGTTACGTAC", "b"), Sequence.Dna("ACGTACGTAC", "c") };
        var msa = new ProgressiveAligner().Align(seqs);
        Assert.Equal(11, msa.Length);
        Assert.Equal(msa[0], msa[2]);
        Assert.Equal(1, msa[0].Count(c => c == '-'));
        Assert.Equal("ACGTTACGTAC", msa[1]);
        Assert.Equal(msa[0], msa.Consensus()); // the inserted column is mostly gaps
        Assert.Equal(10, msa.ColumnConservation().Count(c => c == 1.0));
    }

    [Fact]
    public void Identical_Sequences_HaveNoGaps()
    {
        var s = Sequence.Protein("MKTAYIAKQR", "x");
        var msa = new ProgressiveAligner().Align([s, s.WithId("y"), s.WithId("z")]);
        Assert.All(msa.Rows, r => Assert.Equal("MKTAYIAKQR", r));
        Assert.All(msa.ColumnEntropy(), e => Assert.Equal(0, e));
        Assert.Equal(new string('*', 10), msa.ConservationLine());
    }

    [Fact]
    public void ConsensusConservationAndOutputFormats()
    {
        var msa = new MultipleSequenceAlignment(["s1", "s2", "s3", "s4"], ["MK-LV", "MR-LI", "MKALV", "MKSLV"]);
        Assert.Equal("MK-LV", msa.Consensus());
        Assert.Equal([1.0, 0.75, 0.25, 1.0, 0.75], msa.ColumnConservation());
        Assert.Equal("*: *:", msa.ConservationLine()); // K/R strong group, gap column blank, V/I strong group
        string clustal = msa.ToClustal();
        Assert.StartsWith("CLUSTAL", clustal);
        Assert.Contains("s1        MK-LV 4", clustal);
        string fasta = msa.ToFasta();
        Assert.Equal(">s1\nMK-LV\n>s2\nMR-LI\n>s3\nMKALV\n>s4\nMKSLV\n", fasta.Replace("\r\n", "\n"));
        var parsed = MultipleSequenceAlignment.FromFasta(fasta);
        Assert.Equal(msa.Rows, parsed.Rows);
        Assert.Equal(msa.Ids, parsed.Ids);
        Assert.Equal("MKLV", msa.ToSequences()[0].ToString());
        Assert.Equal(50.0, msa.PercentIdentity(0, 1), 12); // M,L identical; K/R, V/I differ
        Assert.Throws<ArgumentException>(() => new MultipleSequenceAlignment(["a", "b"], ["AC", "A"]));
    }

    [Fact]
    public void SumOfPairs_KnownValue()
    {
        var msa = new MultipleSequenceAlignment(["a", "b"], ["AC-", "ACG"]);
        var m = SubstitutionMatrix.MatchMismatch(2, -1);
        Assert.Equal(2 + 2 - 1, msa.SumOfPairsScore(m, gapScore: -1));
    }

    [Fact]
    public void GuideTree_FromAlignmentDistances_IsNewickSerializable()
    {
        var msa = new ProgressiveAligner { TreeMethod = GuideTreeMethod.NeighborJoining }.Align(Globins);
        var tree = PhylogeneticTree.ParseNewick(msa.GuideTree!.ToNewick());
        Assert.Equal(msa.GuideTree.LeafNames.Order(), tree.LeafNames.Order());
        var dm = DistanceMatrix.FromAlignment(msa, DistanceModel.JukesCantor);
        Assert.True(dm["HBA_HUMAN", "HBA_HORSE"] < dm["HBA_HUMAN", "MYG_PHYCA"]);
    }
}

public class SearchTests
{
    private static string RandomProtein(Random rng, int n)
        => new([.. Enumerable.Range(0, n).Select(_ => "ACDEFGHIKLMNPQRSTVWY"[rng.Next(20)])]);

    [Fact]
    public void KarlinAltschul_Formulas()
    {
        var p = KarlinAltschulParameters.Blosum62Gapped11_1;
        Assert.Equal((0.267 * 100 - Math.Log(0.041)) / Math.Log(2), p.BitScore(100), 12);
        double e = p.EValue(100, 300, 1_000_000);
        Assert.Equal(0.041 * 300 * 1_000_000 * Math.Exp(-0.267 * 100), e, e * 1e-12);
        // E = m·n·2^(−bits)
        Assert.Equal(300.0 * 1_000_000 * Math.Pow(2, -p.BitScore(100)), e, e * 1e-9);
        Assert.Equal(p, KarlinAltschulParameters.Lookup("BLOSUM62", 11, 1));
        Assert.Null(KarlinAltschulParameters.Lookup("PAM30", 9, 1));
    }

    [Fact]
    public void Search_FindsHomologue_WithSignificantEValue()
    {
        var rng = new Random(42);
        string queryText = RandomProtein(rng, 150);
        var db = new List<Sequence>();
        for (int i = 0; i < 200; i++) db.Add(Sequence.Protein(RandomProtein(rng, 100 + rng.Next(200)), $"rand{i}"));
        // A homologue: 80% identity with a 3-residue deletion, embedded in random flanks.
        var chars = queryText.ToCharArray();
        for (int i = 0; i < chars.Length; i += 5) chars[i] = "ACDEFGHIKLMNPQRSTVWY"[rng.Next(20)];
        string homolog = RandomProtein(rng, 40) + new string(chars).Remove(70, 3) + RandomProtein(rng, 30);
        db.Insert(117, Sequence.Protein(homolog, "homolog"));

        var database = new SequenceDatabase(db);
        Assert.Equal(3, database.WordSize);
        var query = Sequence.Protein(queryText, "query");
        var hits = database.Search(query, new SearchOptions { MaxEValue = 1e-3 });
        Assert.NotEmpty(hits);
        Assert.Equal("homolog", hits[0].Subject.Id);
        Assert.Equal(117, hits[0].SubjectIndex);
        Assert.True(hits[0].EValue < 1e-20, hits[0].ToString());
        Assert.True(hits[0].BitScore > 100);
        Assert.Single(hits); // random sequences are not significant
        var aln = hits[0].Alignment!;
        Assert.Equal(hits[0].Score, aln.Score);
        Assert.True(aln.TargetStart >= 30);

        // Exhaustive scoring gives the same top score as seeded scoring.
        var exhaustive = database.Search(query, new SearchOptions { UseSeeds = false, MaxEValue = 1e-3, IncludeAlignments = false });
        Assert.Equal(hits[0].Score, exhaustive[0].Score);
        Assert.Null(exhaustive[0].Alignment);
    }

    [Fact]
    public void Search_Nucleotide_UsesBlastnDefaults()
    {
        var rng = new Random(7);
        string Dna(int n) => new([.. Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)])]);
        string target = Dna(60);
        var db = new SequenceDatabase([Sequence.Dna(Dna(100) + target + Dna(100), "hit"), Sequence.Dna(Dna(300), "miss")]);
        Assert.Equal(11, db.WordSize);
        var hits = db.Search(Sequence.Dna(target, "q"));
        Assert.Equal("hit", hits[0].Subject.Id);
        Assert.Equal(120, hits[0].Score); // 60 matches × 2
        Assert.False(double.IsNaN(hits[0].EValue));
        Assert.Equal(new[] { 0 }, db.Candidates(Sequence.Dna(target, "q")));
    }

    [Fact]
    public void Search_WithoutPublishedStatistics_ReportsNaN()
    {
        var db = new SequenceDatabase([Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQ", "a")]);
        var hits = db.Search(Sequence.Protein("MKTAYIAKQR", "q"), new SearchOptions { SubstitutionMatrix = SubstitutionMatrix.Pam30, GapOpen = 9, GapExtend = 1 });
        Assert.Single(hits);
        Assert.True(double.IsNaN(hits[0].EValue));
    }
}
