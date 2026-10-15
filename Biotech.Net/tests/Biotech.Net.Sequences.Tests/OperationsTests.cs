using Biotech.Net.Core;
using Biotech.Net.Native;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Sequences.Tests;

public class OrfTests
{
    [Fact]
    public void FindsPlantedCdsOfGenBankRecord()
    {
        var r = GenBankReader.Read(TestData.Path("sample.gb")).Single();
        var orfs = OrfFinder.Find(r, new OrfOptions { MinLength = 240 });
        var orf = Assert.Single(orfs, o => o.Strand == Strand.Forward && o.End == 309);
        Assert.Equal(9, orf.Start);
        Assert.Equal(1, orf.Frame);
        Assert.True(orf.HasStop);
        Assert.Equal(r.Features[2]["translation"], orf.Translate(r).ToString());
    }

    [Fact]
    public void ReverseStrandOrf_HasForwardCoordinatesAndNegativeFrame()
    {
        // "ATGAAACCCTAA" reverse-complemented and embedded: GG + TTAGGGTTTCAT + C.
        var s = Sequence.Dna("GGTTAGGGTTTCATC");
        var orf = Assert.Single(OrfFinder.Find(s, new OrfOptions { MinLength = 9 }));
        Assert.Equal((2, 14), (orf.Start, orf.End));
        Assert.Equal(Strand.Reverse, orf.Strand);
        Assert.Equal(-2, orf.Frame);
        Assert.Equal("ATGAAACCCTAA", orf.Extract(s).ToString());
        Assert.Equal("MKP", orf.Translate(s).ToString());
    }

    [Fact]
    public void StartModes_AndNestedAndPartial()
    {
        var s = Sequence.Dna("CCGTGAAAATGCCCTAGTT");
        // Standard code alternative starts include TTG/CTG; table 11 adds GTG.
        var atg = OrfFinder.Find(s, new OrfOptions { MinLength = 6, Reverse = false });
        Assert.Equal((8, 17), (atg.Single().Start, atg.Single().End));
        var gtg = OrfFinder.Find(s, new OrfOptions { MinLength = 6, Reverse = false, Code = GeneticCode.Bacterial, StartMode = OrfStartMode.AlternativeStarts });
        Assert.Equal((2, 17), (gtg.Single().Start, gtg.Single().End));
        var nested = OrfFinder.Find(s, new OrfOptions { MinLength = 6, Reverse = false, Code = GeneticCode.Bacterial, StartMode = OrfStartMode.AlternativeStarts, Nested = true });
        Assert.Equal(2, nested.Count);
        var stopToStop = OrfFinder.Find(Sequence.Dna("AAATAGCCCGGGTTTAAA"), new OrfOptions { MinLength = 6, Reverse = false, StartMode = OrfStartMode.AnyCodon, IncludePartial = true });
        Assert.Contains(stopToStop, o => o.Start == 6 && o.End == 18 && !o.HasStop);
        Assert.Contains(stopToStop, o => o.Start == 0 && o.End == 6 && o.HasStop);
        Assert.Empty(OrfFinder.Find(Sequence.Dna("ATGAAATAG"), new OrfOptions { MinLength = 10 }));
    }

    [Fact]
    public void SixFrame_ReverseFramesTranslateReverseComplement()
    {
        var s = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
        var frames = SixFrame.Translate(s);
        Assert.Equal([1, 2, 3, -1, -2, -3], frames.Select(f => f.Frame));
        Assert.Equal("MAIVMGR*KGAR*", frames[0].Protein.ToString()); // Biopython Seq.translate()
        Assert.Equal(s.ReverseComplement().Translate(frame: 1).ToString(), frames[4].Protein.ToString());
        Assert.EndsWith("_-2", frames[4].Protein.Id);
    }
}

public class MotifTests
{
    [Fact]
    public void Iupac_PalindromeReportedOncePerSite()
    {
        var s = Sequence.Dna("AAGAATCAAGACTCAA");
        var hits = MotifSearch.FindIupac(s, "GANTC");
        Assert.Equal([2, 9], hits.Select(h => h.Start));
        Assert.All(hits, h => Assert.Equal(Strand.Forward, h.Strand));
        Assert.Equal("GACTC", hits[1].Matched);
    }

    [Fact]
    public void Iupac_NonPalindrome_FindsBothStrands()
    {
        var s = Sequence.Dna("TTGGTCTCAAAAGAGACCTT");
        var hits = MotifSearch.FindIupac(s, "GGTCTC");
        Assert.Equal(2, hits.Count);
        Assert.Equal((2, Strand.Forward), (hits[0].Start, hits[0].Strand));
        Assert.Equal((12, Strand.Reverse, "GGTCTC"), (hits[1].Start, hits[1].Strand, hits[1].Matched));
        Assert.Single(MotifSearch.FindIupac(s, "GGTCTC", bothStrands: false));
    }

    [Fact]
    public void Iupac_AmbiguousSequenceBaseOnlyMatchesCompatiblePattern_AndMismatches()
    {
        Assert.Empty(MotifSearch.FindIupac(Sequence.Dna("GANTC"), "GAATC"));
        Assert.Single(MotifSearch.FindIupac(Sequence.Dna("GARTC"), "GANTC"));
        Assert.Equal(2, MotifSearch.FindIupac(Sequence.Dna("TATAAATATTAA"), "TATAAA", bothStrands: false, maxMismatches: 1).Count);
        // ToRegex has to accept exactly what Iupac.Matches accepts, and there T and U share one mask:
        // N admits U as well as T, so a literal T in the pattern must admit U too.
        Assert.Equal("GA[ACGTU][TU]C", Iupac.ToRegex("GANTC"));
    }

    [Fact]
    public void Regex_OverlappingAndReverse()
    {
        var s = Sequence.Dna("AAAAA");
        Assert.Equal(4, MotifSearch.FindRegex(s, "AA").Count);
        var r = MotifSearch.FindRegex(Sequence.Dna("CCCATGGG"), "CAT", bothStrands: true);
        Assert.Contains(r, m => m.Strand == Strand.Forward && m.Start == 2);
        Assert.Contains(r, m => m.Strand == Strand.Reverse && m.Start == 3 && m.Matched == "CAT"); // ATG on top = CAT on bottom
    }

    [Fact]
    public void Prosite_ZincFingerC2H2()
    {
        // PS00028: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
        string re = MotifSearch.PrositeToRegex("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.");
        Assert.Equal("C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H", re);
        // Zif268 finger 1 (UniProt P08046 residues ~ 338-362 region).
        var p = Sequence.Protein("MERPYACPVESCDRRFSRSDELTRHIRIHTGQK");
        var hit = Assert.Single(MotifSearch.FindProsite(p, "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"));
        Assert.Equal(6, hit.Start);
        Assert.Equal("^M[^P]", MotifSearch.PrositeToRegex("<M-{P}"));
        Assert.Equal(2, MotifSearch.Count(Sequence.Dna("ATATA"), "ATA"));
    }
}

public class RestrictionTests
{
    // Expected cut positions from Biopython 1.86 Bio.Restriction (enzyme.search returns 1-based "first base after the cut",
    // i.e. our 0-based Position + 1).
    private const string Seq =
        "GTAGGGGCAATGGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCAGTGCGTGGACACTCGCTATGAATCTCTTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCATCACCCTAAGCCGAATAATGCGTTCGCTCTATTGACTACGACGTAAGCGCAGTATGCCAAGACTATAGGCACTGTCGCATCACAAATCAGGATCGCATTGACTTGTCAGGCGGCAAATAGAGTTTCTCCGACGCCGGGAATTCAAGGTCTCAGGATCCTTTTCTGCAGCCGAGACCAAGAATTCGAGTCCCGGGCATGGCCAACGTGGCCGAAGACAAGCTCTTCAAGTTCAA";

    [Theory]
    [InlineData("EcoRI", new[] { 402, 443 })]
    [InlineData("BamHI", new[] { 417 })]
    [InlineData("BsaI", new[] { 416, 429 })]
    [InlineData("PstI", new[] { 431 })]
    [InlineData("HinfI", new[] { 68, 220, 449 })]
    [InlineData("SmaI", new[] { 455 })]
    [InlineData("BbsI", new[] { 482 })]
    [InlineData("SapI", new[] { 490 })]
    [InlineData("SfiI", new[] { 469 })]
    [InlineData("NlaIII", new[] { 462 })]
    [InlineData("Sau3AI", new[] { 354, 417 })]
    public void CutPositions_MatchBiopython(string enzyme, int[] biopythonPositions)
    {
        var cuts = RestrictionDigest.FindCuts(Sequence.Dna(Seq), RestrictionEnzymes.Get(enzyme));
        Assert.Equal(biopythonPositions, cuts.Select(c => c.Position + 1));
    }

    [Fact]
    public void Fragments_LinearAndCircular_MatchBiopythonCatalyse()
    {
        var s = Sequence.Dna(Seq);
        var linear = RestrictionDigest.Digest(s, false, "EcoRI");
        Assert.Equal([401, 41, 54], linear.Fragments.Select(f => f.Length));
        Assert.Equal([401, 54, 41], linear.FragmentSizes);
        var circular = RestrictionDigest.Digest(s, true, "EcoRI");
        Assert.Equal([41, 455], circular.Fragments.Select(f => f.Length));
        var seqs = circular.FragmentSequences(s);
        Assert.Equal(496, seqs.Sum(f => f.Length));
        Assert.StartsWith("AATTC", seqs[1].ToString());
        Assert.Equal(455, seqs[1].Length);
    }

    [Fact]
    public void EnzymeProperties()
    {
        var ecoRI = RestrictionEnzymes.Get("ecori");
        Assert.True(ecoRI.IsPalindromic);
        Assert.Equal(OverhangType.FivePrime, ecoRI.OverhangType);
        Assert.Equal(4, ecoRI.Overhang);
        Assert.Equal("EcoRI G^AATTC", ecoRI.ToString());
        Assert.Equal(OverhangType.ThreePrime, RestrictionEnzymes.Get("PstI").OverhangType);
        Assert.Equal(OverhangType.Blunt, RestrictionEnzymes.Get("EcoRV").OverhangType);
        var bsaI = RestrictionEnzymes.Get("BsaI");
        Assert.False(bsaI.IsPalindromic);
        Assert.Equal("BsaI GGTCTC(1/5)", bsaI.ToString());
        Assert.True(RestrictionEnzymes.All.Count >= 60);
        Assert.Equal(RestrictionEnzymes.All.Count, RestrictionEnzymes.All.Select(e => e.Name).Distinct().Count());
    }

    [Fact]
    public void SiteSpanningOrigin_FoundOnlyWhenCircular()
    {
        var s = Sequence.Dna("ATTCAAAAAAAAAAAAGA"); // GAATTC wraps around the origin
        Assert.Empty(RestrictionDigest.FindCuts(s, RestrictionEnzymes.Get("EcoRI")));
        var cut = Assert.Single(RestrictionDigest.FindCuts(s, RestrictionEnzymes.Get("EcoRI"), circular: true));
        Assert.Equal(17, cut.Position);
        Assert.Contains(RestrictionEnzymes.Get("EcoRI"), RestrictionEnzymes.SingleCutters(s, circular: true));
        Assert.Contains(RestrictionEnzymes.Get("NotI"), RestrictionEnzymes.NonCutters(s));
    }
}

public class KmerTests
{
    [Fact]
    public void DenseProfile_CountsAndSpectrum()
    {
        var p = KmerProfile.Count(Sequence.Dna("ACGTACGTNACG"), 2);
        Assert.Equal(9, p.Total); // windows containing N are skipped
        Assert.Equal(3, p["AC"]);
        Assert.Equal(3, p["CG"]);
        Assert.Equal(2, p["GT"]);
        Assert.Equal(2, p["TA"] + 1);
        Assert.Equal(4, p.Distinct);
        var spectrum = p.Spectrum();
        Assert.Equal(1, spectrum[1]);
        Assert.Equal(1, spectrum[2]);
        Assert.Equal(2, spectrum[3]);
        Assert.Equal(("AC", 3), p.Top(1)[0] is var t && (t.Kmer == "AC" || t.Kmer == "CG") ? ("AC", t.Count) : t);
    }

    [Fact]
    public void Canonical_MergesReverseComplements()
    {
        var s = Sequence.Dna("AAAAT");
        var p = KmerProfile.Count([s], 3, canonical: true);
        // The 3-mers are AAA, AAA, AAT. Folding each code onto the smaller of itself and its reverse complement
        // keeps AAA (< TTT) and AAT (< ATT), and neither partner occurs, so the counts stay 2 and 1.
        Assert.Equal(2, p["AAA"]);
        Assert.Equal(1, p["AAT"]);
    }

    [Fact]
    public void SparseCounter_MatchesDenseAndHandlesLargeK()
    {
        var s = RandomSequence.Dna(5000, seed: 1);
        var dense = KmerProfile.Count(s, 6);
        var sparse = new KmerCounter(6).Add(s);
        foreach (var (kmer, count) in sparse.Entries()) Assert.Equal(dense[kmer], count);
        Assert.Equal(dense.Distinct, sparse.Distinct);
        var k31 = new KmerCounter(31, canonical: true).Add(s).Add(s.ReverseComplement());
        Assert.All(k31.Codes.Values, c => Assert.True(c >= 2));
        Assert.Equal(k31[s.ToString()[..31]], k31[s.ReverseComplement().ToString()[^31..]]);
    }

    [Fact]
    public void Kmer_EncodeDecodeReverseComplement()
    {
        Assert.True(Kmer.TryEncode("ACGT", out ulong code));
        Assert.Equal(0b00_01_10_11UL, code);
        Assert.Equal("ACGT", Kmer.Decode(code, 4));
        Assert.True(Kmer.TryEncode("AACG", out ulong c2));
        Assert.Equal("CGTT", Kmer.Decode(Kmer.ReverseComplement(c2, 4), 4));
        Assert.False(Kmer.TryEncode("ACNT", out _));
    }

    [Fact]
    public void ProfileDistances()
    {
        var a = KmerProfile.Count(RandomSequence.Dna(20000, 1, 0.5), 3);
        var b = KmerProfile.Count(RandomSequence.Dna(20000, 2, 0.5), 3);
        var gc = KmerProfile.Count(RandomSequence.Dna(20000, 3, 0.8), 3);
        Assert.True(a.CosineSimilarity(b) > a.CosineSimilarity(gc));
        Assert.True(a.EuclideanDistance(b) < a.EuclideanDistance(gc));
        Assert.Equal(1.0, a.CosineSimilarity(a), 12);
    }

    [Fact]
    public void ManagedAndNativeKernelsAgree()
    {
        var s = RandomSequence.Dna(3000, 5);
        bool saved = BioKernels.ForceManaged;
        try
        {
            BioKernels.ForceManaged = true;
            var managed = KmerProfile.Count(s, 5).Counts.ToArray();
            BioKernels.ForceManaged = false;
            Assert.Equal(managed, KmerProfile.Count(s, 5).Counts.ToArray());
        }
        finally { BioKernels.ForceManaged = saved; }
    }
}

public class CodonUsageTests
{
    private static List<Sequence> Reference()
    {
        var r = GenBankReader.Read(TestData.Path("sample.gb")).Single();
        return [r.Features[2].Location.Extract(r), Sequence.Dna("ATGGCTGCTGCAAAACTGCTGCTGGAAGAAGGTAAATAA"), Sequence.Dna("ATGCGTCGTAAAGAAGAATTCGGTGCTTAA")];
    }

    [Fact]
    public void Cai_MatchesBiopython()
    {
        // Biopython 1.86: CodonAdaptationIndex(reference).calculate(...)
        var table = CodonUsageTable.FromSequences(Reference());
        Assert.Equal(0.6992303071484395, table.Cai(Sequence.Dna("ATGGCAAAACTGGAAGGCCGTTTTTAA")), 12);
        Assert.Equal(0.900867659426767, table.Cai(Reference()[1]), 12);
        Assert.Equal(1.0, table.RelativeAdaptiveness("GCT"), 12);
        Assert.Equal(0.2857142857142857, table.RelativeAdaptiveness("GCA"), 12);
        Assert.Equal(1.0, table.RelativeAdaptiveness("CTG"), 12);
    }

    [Fact]
    public void Rscu_PerThousand_Optimize()
    {
        var t = CodonUsageTable.FromSequences([Sequence.Dna("GCTGCTGCTGCA")]);
        Assert.Equal(4, t.Total);
        Assert.Equal(3.0 * 4 / 4, t.Rscu("GCT"));     // Ala has 4 codons: 3 / (4/4)
        Assert.Equal(0.0, t.Rscu("GCC"));
        Assert.Equal(750, t.PerThousand("GCT"));
        Assert.True(double.IsNaN(t.Rscu("TGG")));
        Assert.Equal("GCTGCT", t.Optimize(Sequence.Protein("AA")).ToString());
        Assert.Contains("GCT\tA\t3", t.ToTable());
    }
}

public class ThermodynamicsTests
{
    private static readonly Sequence Primer = Sequence.Dna("CGTTCCAAAGATGTGGGCATGAGCTTAC"); // Biopython MeltingTemp doctest primer

    [Fact]
    public void NearestNeighbor_MatchesBiopythonTmNN()
    {
        // Biopython 1.86 Tm_NN defaults: DNA_NN3, Na=50 mM, dnac1=dnac2=25 nM, saltcorr=5.
        Assert.Equal(60.32091949919129, MeltingTemperature.NearestNeighbor(Primer), 8);
        Assert.Equal(63.792609006046575, MeltingTemperature.NearestNeighbor(Primer, new MeltingConditions { Sodium = 100 }), 8);
        Assert.Equal(66.98544175120679, MeltingTemperature.NearestNeighbor(Primer, new MeltingConditions { Magnesium = 1.5, Dntps = 0.2 }), 8);
        Assert.Equal(29.704794847328003, MeltingTemperature.NearestNeighbor(Sequence.Dna("ATGCATGCAAAT")), 8);
        Assert.Equal(24.67181352679131, MeltingTemperature.NearestNeighbor(Sequence.Dna("GAATTCGAATTC"), selfComplementary: true), 8);
    }

    [Fact]
    public void WallaceAndGcFormulas()
    {
        Assert.Equal(84.0, MeltingTemperature.Wallace(Primer));                         // Biopython Tm_Wallace
        Assert.Equal(58.974330643406475, MeltingTemperature.SaltAdjustedGc(Primer), 8); // Biopython Tm_GC (valueset 7)
        Assert.Equal(64.9 + 41 * (14 - 16.4) / 28, MeltingTemperature.BasicGc(Primer), 10);
        Assert.True(MeltingTemperature.FreeEnergy37(Primer) < -20);
    }

    [Fact]
    public void MolecularWeight_MatchesBiopython()
    {
        // Bio.SeqUtils.molecular_weight (average masses).
        Assert.Equal(949.6095, NucleicAcidMass.MolecularWeight(Sequence.Dna("AGC")), 4);
        Assert.Equal(4960.1769, NucleicAcidMass.MolecularWeight(Sequence.Dna("ATGCATGCAAATTTCG")), 4);
        Assert.Equal(9920.3538, NucleicAcidMass.MolecularWeight(Sequence.Dna("ATGCATGCAAATTTCG"), doubleStranded: true), 4);
        Assert.Equal(2590.5169, NucleicAcidMass.MolecularWeight(Sequence.Rna("AUGCUUAG")), 4);
        Assert.Equal(2471.5748, NucleicAcidMass.MolecularWeight(Sequence.Dna("ATGCATGC"), circular: true), 4);
        Assert.Throws<ArgumentException>(() => NucleicAcidMass.MolecularWeight(Sequence.Dna("ACN")));
    }
}

public class CompositionTests
{
    [Fact]
    public void Dust_HandComputedScores()
    {
        // 10 × A: 8 AAA triplets → Σ c(c−1)/2 = 28, l − 1 = 7 → 4.0.
        Assert.Equal(4.0, Complexity.DustScore(Sequence.Dna("AAAAAAAAAA")), 12);
        // All-distinct triplets score 0.
        Assert.Equal(0.0, Complexity.DustScore(Sequence.Dna("ACGTTGCA")), 12);
    }

    [Fact]
    public void DustProfile_IncrementalEqualsDirect()
    {
        var s = RandomSequence.Dna(300, 4).Concat(Sequence.Dna(string.Concat(Enumerable.Repeat("CA", 50))));
        var profile = Complexity.DustProfile(s, 64);
        Assert.Equal(s.Length - 63, profile.Count);
        foreach (int i in new[] { 0, 17, 150, 250, profile.Count - 1 })
            Assert.Equal(Complexity.DustScore(s.Span.Slice(i, 64)), profile[i].Value, 9);
    }

    [Fact]
    public void Mask_SoftAndHardMasksRepeatOnly()
    {
        var random = RandomSequence.Dna(400, 11);
        var repeat = Sequence.Dna(string.Concat(Enumerable.Repeat("CAG", 40)));
        var s = random.Concat(repeat).Concat(RandomSequence.Dna(400, 12));
        // A 64-base window of a perfect trinucleotide repeat holds 62 triplets of three kinds (21/21/20), so its
        // DUST score is (210 + 210 + 190) / 61 = 10.0. That is well above random (~1) but below dustmasker's
        // default level of 20, which only homopolymer-like runs reach (62 identical triplets score 31), so this
        // test asks for level 8.
        const double Level = 8;
        Assert.Empty(Complexity.LowComplexityRegions(s)); // nothing reaches the default level 20
        var regions = Complexity.LowComplexityRegions(s, threshold: Level);
        var region = Assert.Single(regions);
        Assert.True(region.Start <= 420 && region.End >= 500, region.ToString());
        var soft = Complexity.Mask(s, threshold: Level);
        Assert.Equal(s, soft); // case-insensitive equality: soft-masking preserves residues
        Assert.True(char.IsLower((char)soft[460]));
        Assert.True(char.IsUpper((char)soft[100]));
        Assert.Equal((byte)'N', Complexity.Mask(s, threshold: Level, style: MaskStyle.Hard)[460]);
        Assert.True(Complexity.LinguisticComplexity(random) > Complexity.LinguisticComplexity(repeat));
    }

    [Fact]
    public void SlidingWindow_GcAndSkew()
    {
        var s = Sequence.Dna("GGGGAAAACCCC");
        var gc = SlidingWindow.GcContent(s, 4, 4);
        Assert.Equal([1.0, 0.0, 1.0], gc.Select(w => w.Value));
        var skew = SlidingWindow.GcSkew(s, 4, 4);
        Assert.Equal([1.0, 0.0, -1.0], skew.Select(w => w.Value));
        Assert.Equal([1.0, 1.0, 0.0], SlidingWindow.CumulativeGcSkew(s, 4).Select(w => w.Value));
        var step1 = SlidingWindow.GcContent(s, 5);
        Assert.Equal(8, step1.Count);
        for (int i = 0; i < step1.Count; i++) Assert.Equal(s.Subsequence(i, 5).GcContent(), step1[i].Value, 12);
        Assert.Equal(2.0, SlidingWindow.Entropy(Sequence.Dna("ACGT"), 4).Single().Value, 12);
        Assert.Equal(1.0, SlidingWindow.AtSkew(Sequence.Dna("AAAA"), 4).Single().Value);
    }

    [Fact]
    public void N50_WikipediaExample()
    {
        // Lengths 2..10: total 54, half 27; 10 + 9 + 8 = 27 → N50 = 8, L50 = 3.
        var lengths = Enumerable.Range(2, 9).Select(i => (long)i).ToList();
        Assert.Equal((8L, 3), SequenceStatistics.Nx(lengths, 50));
        Assert.Equal(8, SequenceStatistics.N50(lengths));
        Assert.Equal((4L, 7), SequenceStatistics.Nx(lengths, 90));
    }

    [Fact]
    public void Summarize_SetStatistics()
    {
        var seqs = new[] { Sequence.Dna("GGCC"), Sequence.Dna("AATTNN"), Sequence.Dna("ACGTACGTAC") };
        var sum = SequenceStatistics.Summarize(seqs);
        Assert.Equal(3, sum.Count);
        Assert.Equal(20, sum.TotalLength);
        Assert.Equal(4, sum.MinLength);
        Assert.Equal(10, sum.MaxLength);
        Assert.Equal(6, sum.MedianLength);
        Assert.Equal(10, sum.N50);
        Assert.Equal(2, sum.AmbiguousBases);
        Assert.Equal((4 + 0 + 5) / 18.0, sum.GcContent, 12);
        var hist = SequenceStatistics.LengthHistogram([1, 5, 9, 12], 5);
        Assert.Equal([(0L, 1), (5L, 2), (10L, 1)], hist);
    }

    [Fact]
    public void Random_IsSeededAndShufflePreservesComposition()
    {
        Assert.Equal(RandomSequence.Dna(100, 42), RandomSequence.Dna(100, 42));
        Assert.NotEqual(RandomSequence.Dna(100, 42), RandomSequence.Dna(100, 43));
        Assert.InRange(RandomSequence.Dna(20000, 1, 0.7).GcContent(), 0.68, 0.72);
        var s = RandomSequence.Dna(500, 9);
        var sh = RandomSequence.Shuffle(s, 1);
        Assert.NotEqual(s, sh);
        Assert.Equal(s.Composition(), sh.Composition());
        var blocks = RandomSequence.ShuffleBlocks(Sequence.Dna("AAACCCGGGTT"), 3, 5);
        Assert.Equal("TT", blocks.ToString()[^2..]);
        Assert.Equal(SequenceType.Protein, RandomSequence.Protein(50, 1).Type);
        var weighted = RandomSequence.FromComposition(1000, 3, new Dictionary<char, double> { ['A'] = 1, ['T'] = 0 });
        Assert.Equal(1000, weighted.Composition()['A']);
    }
}
