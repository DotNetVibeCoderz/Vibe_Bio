using Biotech.Net.AI.Preprocessing;

namespace Biotech.Net.AI.Tests;

public class SequencePreprocessorTests
{
    [Fact]
    public void CleanDna_StripsFastaWhitespaceDigitsAndConvertsRna()
    {
        const string raw = ">seq1 some header\n  1 acgu ryNN\n 11 ACGT\n";
        Assert.Equal("ACGTNNNNACGT", SequencePreprocessor.CleanDna(raw));
        Assert.Equal("ACGTRYNNACGT", SequencePreprocessor.CleanDna(raw, ambiguityToN: false));
    }

    [Fact]
    public void CleanDna_RejectsNonNucleotides()
    {
        var ex = Assert.Throws<FormatException>(() => SequencePreprocessor.CleanDna("ACGTEQ"));
        Assert.Contains("residue 5", ex.Message);
        Assert.False(SequencePreprocessor.IsValidDna("ACGTEQ"));
        Assert.True(SequencePreprocessor.IsValidDna("acgt"));
        Assert.False(SequencePreprocessor.IsValidDna("   "));
    }

    [Fact]
    public void CleanProtein_HandlesStopsGapsAndJ()
    {
        Assert.Equal("MKTXA", SequencePreprocessor.CleanProtein(">p\nmktja*\n"));
        Assert.Throws<FormatException>(() => SequencePreprocessor.CleanProtein("MK*TA"));
        Assert.Throws<FormatException>(() => SequencePreprocessor.CleanProtein("MK-TA"));
        Assert.Equal("MKTA", SequencePreprocessor.CleanProtein("MK-T.A", removeGaps: true));
        Assert.False(SequencePreprocessor.IsValidProtein("MK1#"));
    }

    [Fact]
    public void Kmers_OverlappingAndTiled()
    {
        Assert.Equal(["ATG", "TGC", "GCA", "CAT"], SequencePreprocessor.Kmers("ATGCAT", 3));
        Assert.Equal(["ATG", "CAT"], SequencePreprocessor.Kmers("ATGCAT", 3, stride: 3));
        Assert.Empty(SequencePreprocessor.Kmers("AT", 3));
    }

    [Fact]
    public void KmerTokenize_MatchesDnabertConvention()
    {
        Assert.Equal("ATGCAT TGCATG GCATGC", SequencePreprocessor.KmerTokenize("atgcatgc", 6));
        Assert.Throws<ArgumentException>(() => SequencePreprocessor.KmerTokenize("ATG", 6));
        var kmers = SequencePreprocessor.Kmers("ATGCATGCTTAG", 6);
        Assert.Equal("ATGCATGCTTAG", SequencePreprocessor.FromOverlappingKmers(kmers));
    }

    [Fact]
    public void SpaceResidues_MapsRareAminoAcidsToX()
    {
        Assert.Equal("M K X X X X A", SequencePreprocessor.SpaceResidues("MKUZOBA"));
        Assert.Equal("MKUA", SequencePreprocessor.EsmInput("mkua"));
        Assert.Equal("MKTA", SequencePreprocessor.EsmInput(Biotech.Net.Core.Sequence.Protein("MKTA")));
    }

    [Fact]
    public void Chunk_CoversSequenceWithOverlap()
    {
        string s = new string('A', 10) + new string('C', 10) + new string('G', 5);
        var chunks = SequencePreprocessor.Chunk(s, 10, 3);
        Assert.All(chunks, c => Assert.Equal(10, c.Length));
        Assert.Equal(0, chunks[0].Start);
        Assert.Equal(s.Length, chunks[^1].End);
        for (int i = 1; i < chunks.Count; i++) Assert.True(chunks[i].Start <= chunks[i - 1].End - 3);
        foreach (var c in chunks) Assert.Equal(s.Substring(c.Start, c.Length), c.Text);

        Assert.Single(SequencePreprocessor.Chunk("ACGT", 10, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequencePreprocessor.Chunk(s, 5, 5));
    }
}

public class GeneformerEncoderTests
{
    private static readonly Dictionary<string, int> Tokens = new()
    {
        ["<pad>"] = 0, ["<mask>"] = 1, ["<cls>"] = 2, ["<eos>"] = 3,
        ["G1"] = 4, ["G2"] = 5, ["G3"] = 6, ["G4"] = 7, ["G5"] = 8,
    };

    private static readonly Dictionary<string, double> Medians = new()
    {
        ["G1"] = 1.0, ["G2"] = 10.0, ["G3"] = 0.5, ["G4"] = 2.0, ["G5"] = 1.0,
    };

    [Fact]
    public void Encode_RanksByMedianScaledExpression()
    {
        var enc = new GeneformerEncoder(Tokens, Medians, GeneformerEncoderOptions.V1);
        // counts 10,10,10,0,10 -> scaled by 1/median: 10, 1, 20, -, 10 (common factor ignored)
        var detail = enc.EncodeDetailed(["G1", "G2", "G3", "G4", "G5"], new double[] { 10, 10, 10, 0, 10 });
        Assert.Equal(["G3", "G1", "G5", "G2"], detail.GeneIds);
        Assert.Equal([6, 4, 8, 5], detail.TokenIds);
        Assert.Equal(20.0 / 40 * 10_000, detail.NormalizedValues[0], 9);
    }

    [Fact]
    public void Encode_V2AddsSpecialTokensAndTruncates()
    {
        var enc = new GeneformerEncoder(Tokens, Medians, GeneformerEncoderOptions.V2 with { MaxLength = 4 });
        var ids = enc.Encode(new Dictionary<string, double> { ["G1"] = 5, ["G2"] = 100, ["G3"] = 1, ["UNKNOWN"] = 1000 });
        Assert.Equal(4, ids.Length);
        Assert.Equal(2, ids[0]);
        Assert.Equal(3, ids[^1]);
        Assert.Equal([5, 4], ids[1..3]); // G2: 100/10=10, G1: 5/1=5, G3 (1/0.5=2) truncated; unknown gene dropped
    }

    [Fact]
    public void Encode_TiesKeepInputOrderAndEmptyCellGivesSpecialsOnly()
    {
        var enc = new GeneformerEncoder(Tokens, Medians, GeneformerEncoderOptions.V1);
        Assert.Equal([4, 8], enc.Encode(["G1", "G5"], new double[] { 3, 3 }));
        Assert.Equal([8, 4], enc.Encode(["G5", "G1"], new double[] { 3, 3 }));
        var v2 = new GeneformerEncoder(Tokens, Medians);
        Assert.Equal([2, 3], v2.Encode(["G1"], new double[] { 0 }));
    }

    [Fact]
    public void EncodeMatrix_EncodesEachRow()
    {
        var enc = new GeneformerEncoder(Tokens, Medians, GeneformerEncoderOptions.V1);
        var rows = enc.EncodeMatrix(["G1", "G2"], new double[] { 1, 0, 0, 1 });
        Assert.Equal([4], rows[0]);
        Assert.Equal([5], rows[1]);
    }

    [Fact]
    public void V2WithoutClsToken_Throws()
    {
        var noSpecials = new Dictionary<string, int> { ["<pad>"] = 0, ["G1"] = 1 };
        Assert.Throws<ArgumentException>(() => new GeneformerEncoder(noSpecials, Medians, GeneformerEncoderOptions.V2));
    }
}

public class GeneformerDictionaryTests
{
    // Protocol-4 pickle of {'<pad>': 0, 'ENSG1': 70000, 'ENSG2': numpy.float64(1.5)} built opcode by opcode.
    private static byte[] BuildPickle()
    {
        var b = new List<byte> { 0x80, 4, (byte)'}', 0x94, (byte)'(' };
        void Str(string s) { b.Add(0x8c); b.Add((byte)s.Length); b.AddRange(System.Text.Encoding.ASCII.GetBytes(s)); b.Add(0x94); }
        Str("<pad>"); b.AddRange([(byte)'K', 0]);
        Str("ENSG1"); b.Add((byte)'J'); b.AddRange(BitConverter.GetBytes(70000));
        Str("ENSG2");
        Str("numpy.core.multiarray"); Str("scalar"); b.Add(0x93); b.Add(0x94);
        Str("numpy"); Str("dtype"); b.Add(0x93); b.Add(0x94);
        Str("f8"); b.Add(0x89); b.Add(0x88); b.Add(0x87); b.Add((byte)'R'); b.Add(0x94);
        b.Add((byte)'('); b.AddRange([(byte)'K', 3]); Str("<"); b.AddRange([(byte)'N', (byte)'N', (byte)'N']);
        b.Add((byte)'J'); b.AddRange(BitConverter.GetBytes(-1)); b.Add((byte)'J'); b.AddRange(BitConverter.GetBytes(-1));
        b.AddRange([(byte)'K', 0, (byte)'t', (byte)'b']);
        b.Add((byte)'C'); b.Add(8); b.AddRange(BitConverter.GetBytes(1.5));
        b.Add(0x86); b.Add((byte)'R');
        b.AddRange([(byte)'u', (byte)'.']);
        return [.. b];
    }

    [Fact]
    public void ReadPickle_ReadsPlainValuesAndNumpyScalars()
    {
        var d = GeneformerDictionaries.ReadPickle(new MemoryStream(BuildPickle()));
        Assert.Equal(0, d["<pad>"]);
        Assert.Equal(70000, d["ENSG1"]);
        Assert.Equal(1.5, d["ENSG2"]);
    }

    [Fact]
    public void ReadPickle_RefusesArbitraryGlobals()
    {
        var bytes = new List<byte> { 0x80, 4 };
        foreach (var s in new[] { "os", "system" }) { bytes.Add(0x8c); bytes.Add((byte)s.Length); bytes.AddRange(System.Text.Encoding.ASCII.GetBytes(s)); }
        bytes.AddRange([0x93, (byte)'.']);
        var ex = Assert.Throws<FormatException>(() => GeneformerDictionaries.ReadPickle(new MemoryStream([.. bytes])));
        Assert.Contains("os.system", ex.Message);
    }

    [Fact]
    public void ReadTable_SkipsHeaderAndParsesTabsOrCommas()
    {
        var d = GeneformerDictionaries.ReadTable(new StringReader("gene\tmedian\nENSG1\t2.5\nENSG2,0.25\n"));
        Assert.Equal(2, d.Count);
        Assert.Equal(0.25, d["ENSG2"]);
    }
}
