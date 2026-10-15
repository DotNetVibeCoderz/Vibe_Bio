using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Core.Tests;

public class SequenceTests
{
    [Fact]
    public void ReverseComplement_HandlesIupacAndCase()
    {
        var s = Sequence.Dna("ATGCRYacgtN");
        Assert.Equal("NacgtRYGCAT", s.ReverseComplement().ToString());
    }

    [Fact]
    public void ReverseComplement_OfRna_UsesUracil()
    {
        Assert.Equal("UACG", Sequence.Rna("CGUA").ReverseComplement().ToString());
    }

    [Fact]
    public void Translate_StandardCode()
    {
        var s = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
        Assert.Equal("MAIVMGR*KGAR*", s.Translate().ToString());
        Assert.Equal("MAIVMGR", s.Translate(toStop: true).ToString());
    }

    [Fact]
    public void Translate_AmbiguousCodonIsX_AndRnaWorks()
    {
        Assert.Equal("MX", Sequence.Dna("ATGNNN").Translate().ToString());
        Assert.Equal("M*", Sequence.Rna("AUGUAA").Translate().ToString());
    }

    [Fact]
    public void Translate_VertebrateMitochondrial_TgaIsTrp()
    {
        Assert.Equal("W", Sequence.Dna("TGA").Translate(GeneticCode.VertebrateMitochondrial).ToString());
        Assert.Equal("*", Sequence.Dna("AGA").Translate(GeneticCode.VertebrateMitochondrial).ToString());
    }

    [Fact]
    public void GeneticCode_StartAndStop()
    {
        Assert.True(GeneticCode.Standard.IsStart("ATG"));
        Assert.True(GeneticCode.Standard.IsStop("TAG"));
        Assert.Equal(6, GeneticCode.Standard.CodonsFor('L').Count);
        Assert.Equal(3, GeneticCode.Standard.CodonsFor('*').Count);
    }

    [Fact]
    public void GcContent_IgnoresAmbiguousBases()
    {
        Assert.Equal(0.5, Sequence.Dna("ACGTNNNN").GcContent(), 12);
    }

    [Fact]
    public void Invalid_Symbol_Throws()
    {
        Assert.Throws<ArgumentException>(() => Sequence.Dna("ACGTZ"));
    }

    [Fact]
    public void Detect_TellsTypesApart()
    {
        Assert.Equal(SequenceType.Dna, Sequence.Parse("ACGTACGTNN").Type);
        Assert.Equal(SequenceType.Rna, Sequence.Parse("ACGUACGU").Type);
        Assert.Equal(SequenceType.Protein, Sequence.Parse("MKVLAAGIVGLLLW").Type);
    }

    [Fact]
    public void IndexOf_IsCaseInsensitive()
    {
        var s = Sequence.Dna("ttgaattcAA");
        Assert.Equal(2, s.IndexOf("GAATTC"));
        Assert.True(s.Contains("gaattc"));
    }

    [Fact]
    public void FeatureLocation_ExtractsJoinOnReverseStrand()
    {
        var parent = Sequence.Dna("AAACCCGGGTTT");
        var loc = new FeatureLocation([(0, 3), (9, 12)], Strand.Reverse);
        Assert.Equal("AAATTT", loc.Extract(parent).ToString());
        Assert.Equal("complement(join(1..3,10..12))", loc.ToString());
    }

    [Fact]
    public void SequenceRead_TrimsLowQualityEnds()
    {
        var read = SequenceRead.FromAscii("r1", "ACGTACGTAC", "##IIIIII##");
        var trimmed = read.TrimByQuality(20, 2);
        Assert.Equal("GTACGT", trimmed.Sequence.ToString());
        Assert.Equal(40, trimmed.MeanQuality());
    }
}

public class GenomicTests
{
    [Fact]
    public void Interval_ParseAndOverlap()
    {
        var a = GenomicInterval.Parse("chr1:101-200");
        Assert.Equal(100, a.Start);
        Assert.Equal(200, a.End);
        Assert.True(a.Overlaps(GenomicInterval.Parse("chr1:200-300")));
        Assert.False(a.Overlaps(GenomicInterval.Parse("chr1:201-300")));
        Assert.Equal("chr1:101-200", a.ToString());
    }

    [Fact]
    public void Interval_Merge()
    {
        var merged = GenomicInterval.Merge([new("c", 0, 10), new("c", 5, 20), new("c", 30, 40), new("d", 0, 5)]);
        Assert.Equal(3, merged.Count);
        Assert.Equal(20, merged[0].End);
    }

    [Fact]
    public void Variant_Classification()
    {
        Assert.True(new Variant("1", 10, "A", "G").IsTransition);
        Assert.True(new Variant("1", 10, "A", "C").IsTransversion);
        Assert.Equal(VariantType.Insertion, new Variant("1", 10, "A", "ATT").Type);
        Assert.Equal(VariantType.Deletion, new Variant("1", 10, "ATT", "A").Type);
    }

    [Fact]
    public void Genotype_Parse()
    {
        var g = Genotype.Parse("0|1");
        Assert.True(g.Phased);
        Assert.True(g.IsHeterozygous);
        Assert.True(Genotype.Parse("./.").IsMissing);
        Assert.True(Genotype.Parse("1/1").IsHomozygousAlternate);
    }

    [Fact]
    public void Genome_SequenceAndN50()
    {
        var g = new Genome("toy")
            .Add(Chromosome.FromSequence(Sequence.Dna("ACGTACGTAC", "c1")))
            .Add(Chromosome.FromSequence(Sequence.Dna("GGGCCC", "c2")))
            .Add(Chromosome.FromSequence(Sequence.Dna("TT", "c3")));
        Assert.Equal(10, g.N50());
        Assert.Equal("CGT", g.GetSequence(new GenomicInterval("c1", 1, 4)).ToString());
        Assert.Equal("ACG", g.GetSequence(new GenomicInterval("c1", 1, 4, Strand.Reverse)).ToString());
    }

    [Fact]
    public void Protein_PropertiesAreSensible()
    {
        var p = Protein.FromString("INS", "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN");
        Assert.InRange(p.MolecularWeight, 11900, 12000); // preproinsulin ≈ 11,981 Da
        Assert.InRange(p.IsoelectricPoint, 5.0, 5.6);
    }
}

public class KernelTests
{
    [Fact]
    public void KmerCounting()
    {
        var counts = new int[16];
        long n = BioKernels.CountKmers("ACGTNACG"u8, 2, counts);
        Assert.Equal(5, n); // AC CG GT | AC CG
        Assert.Equal(2, counts[0b0001]); // AC
        Assert.Equal(2, counts[0b0110]); // CG
        Assert.Equal(1, counts[0b1011]); // GT
    }

    [Fact]
    public void Hamming_LongInputsUseVectorPath()
    {
        var a = new byte[1000];
        var b = new byte[1000];
        Array.Fill(a, (byte)'A');
        Array.Fill(b, (byte)'A');
        b[3] = b[500] = b[999] = (byte)'C';
        Assert.Equal(3, BioKernels.HammingDistance(a, b));
    }

    [Fact]
    public void WelchTTest_MatchesReference()
    {
        // Row 0: R t.test(c(1,2,3,4), c(3,4,5,6)) → t = 2.1909, df = 6, p = 0.07098
        double[] m = [1, 2, 3, 4, 3, 4, 5, 6];
        var t = new double[1];
        var p = new double[1];
        var d = new double[1];
        BioKernels.WelchTTest(m, 1, 8, [0, 1, 2, 3], [4, 5, 6, 7], t, p, d);
        Assert.Equal(2.19089, t[0], 4);
        Assert.Equal(0.07098, p[0], 4);
        Assert.Equal(2, d[0]);
    }

    [Fact]
    public void BenjaminiHochberg_MatchesR()
    {
        // p.adjust(c(0.01, 0.04, 0.03, 0.005), "BH") → 0.02 0.04 0.04 0.02
        double[] p = [0.01, 0.04, 0.03, 0.005];
        var q = new double[4];
        BioKernels.BenjaminiHochberg(p, q);
        Assert.Equal([0.02, 0.04, 0.04, 0.02], q.Select(x => Math.Round(x, 10)));
    }

    [Fact]
    public void SpecialFunctions_Reference()
    {
        Assert.Equal(0.975, SpecialFunctions.NormalCdf(1.959964), 6);
        Assert.Equal(1.959964, SpecialFunctions.NormalQuantile(0.975), 5);
        Assert.Equal(0.05, SpecialFunctions.ChiSquareUpper(3.841459, 1), 6);
        Assert.Equal(Math.Log(120), SpecialFunctions.LogGamma(6), 10);
    }

    [Fact]
    public void Dataset_RowStatsAndSelection()
    {
        var d = Dataset.FromArray("x", new double[,] { { 1, 2, 3 }, { 4, 5, 6 } }, ["g1", "g2"], ["a", "b", "c"])
            .Annotate("condition", "ctrl", "ctrl", "trt");
        Assert.Equal([2.0, 5.0], d.RowMeans());
        Assert.Equal([0, 1], d.ColumnsWhere("condition", "ctrl"));
        var sub = d.SelectColumns([2]);
        Assert.Equal(6, sub["g2", "c"]);
        Assert.Equal("trt", sub.ColumnAnnotations["condition"][0]);
    }
}
