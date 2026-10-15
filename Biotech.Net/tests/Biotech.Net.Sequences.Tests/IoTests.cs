using System.IO.Compression;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Sequences.Tests;

internal static class TestData
{
    public static string Path(string name) => System.IO.Path.Combine(AppContext.BaseDirectory, "data", name);

    public static string TempFile(string name)
    {
        string dir = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "biotechnet-seq-tests", Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        return System.IO.Path.Combine(dir, name);
    }
}

public class FastaTests
{
    private const string Multi = ">seq1 first record\nACGTACGTAC\nGTACG\n\n;comment line\n>seq2\nMKVLAAGIVG\n>empty\n>seq3 rna\nACGUACGU\n";

    [Fact]
    public void Parse_MultiRecord_SplitsIdDescriptionAndDetectsType()
    {
        var recs = FastaReader.Parse(Multi);
        Assert.Equal(4, recs.Count);
        Assert.Equal("seq1", recs[0].Id);
        Assert.Equal("first record", recs[0].Description);
        Assert.Equal("ACGTACGTACGTACG", recs[0].ToString());
        Assert.Equal(SequenceType.Dna, recs[0].Type);
        Assert.Equal(SequenceType.Protein, recs[1].Type);
        Assert.Equal(0, recs[2].Length);
        Assert.Equal(SequenceType.Rna, recs[3].Type);
    }

    [Fact]
    public void Read_DataBeforeHeader_Throws()
    {
        Assert.Throws<InvalidDataException>(() => FastaReader.Parse("ACGT\n>x\nAC"));
    }

    [Fact]
    public void Writer_WrapsAtConfiguredWidth_AndRoundTrips()
    {
        var s = Sequence.Dna(new string('A', 25) + new string('C', 10), "x", "desc here");
        string text = FastaWriter.ToString([s], lineWidth: 10);
        Assert.Equal(">x desc here\nAAAAAAAAAA\nAAAAAAAAAA\nAAAAACCCCC\nCCCCC\n", text);
        var back = FastaReader.Parse(text).Single();
        Assert.Equal(s, back);
        Assert.Equal("desc here", back.Description);
        Assert.Equal(">x desc here\n" + s + "\n", FastaWriter.ToString([s], lineWidth: 0));
    }

    [Fact]
    public void ReadsGzipTransparently_AndAsync()
    {
        string path = TestData.TempFile("t.fa.gz");
        FastaWriter.Write(path, FastaReader.Parse(Multi));
        using (var fs = File.OpenRead(path)) Assert.Equal(0x1f, fs.ReadByte()); // really compressed
        Assert.Equal(4, FastaReader.Read(path).Count());
        // gzip detection is by magic bytes, so a misnamed file still opens.
        string misnamed = TestData.TempFile("t.fa");
        File.Copy(path, misnamed);
        Assert.Equal("seq1", FastaReader.Read(misnamed).First().Id);
    }

    [Fact]
    public async Task ReadAsync_StreamsRecords()
    {
        using var ms = new MemoryStream(Encoding.ASCII.GetBytes(Multi));
        var ids = new List<string>();
        await foreach (var s in FastaReader.ReadAsync(ms)) ids.Add(s.Id);
        Assert.Equal(["seq1", "seq2", "empty", "seq3"], ids);
    }

    [Fact]
    public void ReadAsync_GzipNonSeekableStream()
    {
        var gz = new MemoryStream();
        using (var z = new GZipStream(gz, CompressionLevel.Fastest, leaveOpen: true)) z.Write(Encoding.ASCII.GetBytes(Multi));
        gz.Position = 0;
        using var stream = BioFile.WrapDecompression(new NonSeekable(gz));
        using var reader = new StreamReader(stream);
        Assert.Equal(4, FastaReader.Read(reader).Count());
    }

    private sealed class NonSeekable(Stream inner) : Stream
    {
        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
        public override void Flush() { }
        public override int Read(byte[] buffer, int offset, int count) => inner.Read(buffer, offset, count);
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
    }
}

public class FastaIndexTests
{
    // chr1: 25 bases in lines of 10; chr2: 7 bases on one line; chr3: 12 bases in lines of 5.
    private const string Fasta = ">chr1 test\nACGTACGTAA\nCCCCGGGGTT\nTTGCA\n>chr2\nGATTACA\n>chr3\nAAAAA\nCCCCC\nGG\n";

    private static string WriteFasta(string text)
    {
        string path = TestData.TempFile("genome.fa");
        File.WriteAllText(path, text);
        return path;
    }

    [Fact]
    public void Build_MatchesSamtoolsFaidxColumns()
    {
        var idx = FastaIndex.Build(WriteFasta(Fasta));
        // Expected values computed as samtools faidx does: offset = byte after the header line.
        Assert.Equal(new FastaIndexEntry("chr1", 25, 11, 10, 11), idx["chr1"]);
        Assert.Equal(new FastaIndexEntry("chr2", 7, 45, 7, 8), idx["chr2"]);
        Assert.Equal(new FastaIndexEntry("chr3", 12, 59, 5, 6), idx["chr3"]);
        var sw = new StringWriter();
        idx.Write(sw);
        Assert.StartsWith("chr1\t25\t11\t10\t11\n", sw.ToString());
        Assert.Equal(idx.Entries, FastaIndex.Read(new StringReader(sw.ToString())).Entries);
    }

    [Fact]
    public void Fetch_RegionsAcrossLineBreaks()
    {
        string path = WriteFasta(Fasta);
        using var fa = IndexedFastaReader.Open(path);
        Assert.True(File.Exists(path + ".fai"));
        Assert.Equal("ACGTACGTAACCCCGGGGTTTTGCA", fa.FetchAll("chr1").ToString());
        Assert.Equal("AACCCC", fa.Fetch("chr1", 8, 14).ToString());
        Assert.Equal("CCCCCGG", fa.Fetch("chr3", 5, 100).ToString()); // end clamped
        Assert.Equal("ATTAC", fa.Fetch("chr2:2-6").ToString());          // samtools-style 1-based region
        Assert.Equal("GTAAT", fa.Fetch(new GenomicInterval("chr2", 1, 6, Strand.Reverse)).ToString());
        Assert.Equal(3, fa.ToGenome().Chromosomes.Count);
    }

    [Fact]
    public void Fetch_IsThreadSafe()
    {
        using var fa = IndexedFastaReader.Open(WriteFasta(Fasta), saveIndex: false);
        Parallel.For(0, 200, i => Assert.Equal("GATTACA", fa.FetchAll("chr2").ToString()));
    }

    [Fact]
    public void Build_RejectsUnevenLines()
    {
        Assert.Throws<InvalidDataException>(() => FastaIndex.Build(WriteFasta(">a\nACGT\nAC\nACGT\n")));
    }

    [Fact]
    public void Build_HandlesCrLfAndMissingFinalNewline()
    {
        string path = WriteFasta(">a\r\nACGT\r\nAC\r\n>b\nGGG");
        var idx = FastaIndex.Build(path);
        Assert.Equal(new FastaIndexEntry("a", 6, 4, 4, 6), idx["a"]);
        Assert.Equal(3, idx["b"].Length);
        using var fa = IndexedFastaReader.Open(path, saveIndex: false);
        Assert.Equal("GTAC", fa.Fetch("a", 2, 6).ToString());
        Assert.Equal("GGG", fa.FetchAll("b").ToString());
    }
}

public class FastqTests
{
    [Fact]
    public void Parse_FourLineRecords_Phred33()
    {
        const string fq = "@r1 lane1\nACGTN\n+\nII#5!\n@r2\nGG\n+r2\nAB\n";
        var reads = FastqReader.Parse(fq);
        Assert.Equal(2, reads.Count);
        Assert.Equal("r1", reads[0].Id);
        Assert.Equal("lane1", reads[0].Sequence.Description);
        Assert.Equal(new byte[] { 40, 40, 2, 20, 0 }, reads[0].Quality.ToArray());
        Assert.Equal(32, reads[1].Quality[0]); // 'A' - 33
    }

    [Fact]
    public void Parse_MultiLineRecord_WithQualityStartingWithAt()
    {
        // The second quality line starts with '@', which must not be mistaken for a header.
        const string fq = "@r1\nACGT\nACGT\n+\nIIII\n@III\n@r2\nA\n+\nI\n";
        var reads = FastqReader.Parse(fq);
        Assert.Equal(2, reads.Count);
        Assert.Equal("ACGTACGT", reads[0].Sequence.ToString());
        Assert.Equal(31, reads[0].Quality[4]); // '@' - 33
    }

    [Fact]
    public void DetectsPhred64()
    {
        // Illumina 1.5: all characters ≥ '@' (64) and some above 'J' (74), e.g. 'h' = Q40.
        const string fq = "@r1\nACGT\n+\nhhhB\n";
        Assert.Equal(FastqEncoding.Phred64, FastqReader.DetectEncoding(["hhhB"]));
        var r = FastqReader.Parse(fq).Single();
        Assert.Equal(new byte[] { 40, 40, 40, 2 }, r.Quality.ToArray());
        Assert.Equal(FastqEncoding.Phred33, FastqReader.DetectEncoding(["II5#"]));
    }

    [Fact]
    public void Malformed_Throws()
    {
        Assert.Throws<InvalidDataException>(() => FastqReader.Parse("@r1\nACGT\n+\nIII\n"));
        Assert.Throws<InvalidDataException>(() => FastqReader.Parse("r1\nACGT\n+\nIIII\n"));
        Assert.Throws<InvalidDataException>(() => FastqReader.Parse("@r1\nACGT\n+\n!!!!\n", FastqEncoding.Phred64));
    }

    [Fact]
    public void Writer_RoundTrips_AndConvertsEncoding()
    {
        var reads = FastqReader.Parse("@r1 x\nACGTN\n+\nII#5!\n");
        var sw = new StringWriter();
        FastqWriter.Write(sw, reads);
        Assert.Equal("@r1 x\nACGTN\n+\nII#5!\n", sw.ToString());
        var sw64 = new StringWriter();
        FastqWriter.Write(sw64, reads, FastqEncoding.Phred64);
        Assert.Equal("@r1 x\nACGTN\n+\nhhBT@\n", sw64.ToString());
        Assert.Equal(reads[0].Quality.ToArray(), FastqReader.Parse(sw64.ToString(), FastqEncoding.Phred64)[0].Quality.ToArray());
    }

    [Fact]
    public void FileRoundTrip_Gzip()
    {
        string path = TestData.TempFile("reads.fq.gz");
        var reads = Enumerable.Range(0, 50).Select(i => SequenceRead.FromAscii($"r{i}", "ACGTACGT", "IIIIHHHH")).ToList();
        FastqWriter.Write(path, reads);
        var back = SequenceFile.ReadReads(path).ToList();
        Assert.Equal(50, back.Count);
        Assert.Equal(39, back[49].Quality[4]);
    }
}

public class GenBankTests
{
    // data/sample.gb is a synthetic record; expected values come from Biopython 1.86 SeqIO.read(..., "genbank").
    private static Sequence Sample() => GenBankReader.Read(TestData.Path("sample.gb")).Single();

    [Fact]
    public void Header_MatchesBiopython()
    {
        var r = Sample();
        Assert.Equal("TEST0001.2", r.Id);
        Assert.Equal("Escherichia coli synthetic test locus, abcA gene, complete cds; and xyzB pseudogene", r.Description);
        Assert.Equal(700, r.Length);
        Assert.Equal("TEST0001", r.Annotations[RecordKeys.Locus]);
        Assert.Equal("linear", r.Annotations[RecordKeys.Topology]);
        Assert.Equal("BCT", r.Annotations[RecordKeys.Division]);
        Assert.Equal("DNA", r.Annotations[RecordKeys.MoleculeType]);
        Assert.Equal("15-MAR-2024", r.Annotations[RecordKeys.Date]);
        Assert.Equal("Escherichia coli", r.Annotations[RecordKeys.Organism]);
        Assert.StartsWith("Bacteria; Pseudomonadota;", r.Annotations[RecordKeys.Taxonomy]);
        Assert.EndsWith("Escherichia", r.Annotations[RecordKeys.Taxonomy]);
        Assert.Equal("A synthetic GenBank record for parser tests", r.Annotations[RecordKeys.Reference(1, "title")]);
        Assert.Equal("This record is synthetic.\nSecond comment line.", r.Annotations[RecordKeys.Comment]);
        Assert.StartsWith("GTAGGGGCAA", r.ToString());
    }

    [Fact]
    public void Features_LocationsAndQualifiers()
    {
        var f = Sample().Features;
        Assert.Equal(8, f.Count);
        Assert.Equal(["source", "gene", "CDS", "CDS", "mRNA", "misc_feature", "misc_binding", "repeat_region"], f.Select(x => x.Type));
        var cds = f[2];
        Assert.Equal((9, 309), cds.Location.Segments.Single());
        Assert.Equal("1", cds["codon_start"]);
        Assert.Equal("putative ABC transporter ATP-binding protein with an unusually long product name", cds["product"]);
        Assert.Equal("contains \"quoted\" text", cds["note"]);
        Assert.Equal("MAKDNYITYTSARNLLAQCESLKGVSVMHTPLLAVSTPSDWHFYYTQKQNSGNFDRSRRGAPSSAWTLAMNLFTHSAKLQRGQFHHPKPNNAFALLTTT", cds["translation"]);
        Assert.Equal("", f[3]["pseudo"]);
        // order(...) and partial ends are recorded for faithful round-trips.
        Assert.True(InsdcLocationParser.GetDetails(f[5])!.IsOrder);
        Assert.True(InsdcLocationParser.IsPartial(f[4]));
        Assert.Equal((560, 560), f[6].Location.Segments.Single()); // 560^561 between bases
    }

    [Fact]
    public void CdsTranslation_MatchesTranslationQualifier()
    {
        var r = Sample();
        var cds = r.Features[2];
        Assert.Equal(cds["translation"], cds.Location.Extract(r).Translate(GeneticCode.Bacterial, toStop: true).ToString());
    }

    [Fact]
    public void ComplementJoin_ExtractsLikeBiopython()
    {
        var r = Sample();
        var cds2 = r.Features[3];
        Assert.Equal(Strand.Reverse, cds2.Location.Strand);
        Assert.Equal([(349, 400), (448, 520)], cds2.Location.Segments);
        var extracted = cds2.Location.Extract(r);
        Assert.StartsWith("ATGCGCTCATTCCCTTGTCG", extracted.ToString());
        Assert.Equal("MRSFPCRRVMEQGRCLRLEDRCTRPASEKLYLPPDKSMRS*", extracted.Translate().ToString());
        Assert.Equal(111, r.Features[4].Location.Extract(r).Length);
    }

    [Theory]
    [InlineData("467", 466, 467, Strand.Forward)]
    [InlineData("340..565", 339, 565, Strand.Forward)]
    [InlineData("complement(340..565)", 339, 565, Strand.Reverse)]
    [InlineData("<345..>500", 344, 500, Strand.Forward)]
    [InlineData("J00194.1:100..202,1..10", 0, 10, Strand.Forward)]
    public void LocationParser_SimpleForms(string text, int start, int end, Strand strand)
    {
        var loc = InsdcLocationParser.Parse(text.Contains(',') ? $"join({text})" : text);
        Assert.Equal(start, loc.Location.Start);
        Assert.Equal(end, loc.Location.End);
        Assert.Equal(strand, loc.Location.Strand);
    }

    [Fact]
    public void LocationParser_JoinOfComplements_EqualsComplementOfJoin()
    {
        var a = InsdcLocationParser.Parse("complement(join(1..10,20..30))");
        var b = InsdcLocationParser.Parse("join(complement(20..30),complement(1..10))");
        Assert.Equal(a.Location.Segments, b.Location.Segments);
        Assert.Equal("complement(join(1..10,20..30))", InsdcLocationParser.Format(b));
        var partial = InsdcLocationParser.Parse("complement(<1..>200)");
        Assert.True(partial.PartialStart && partial.PartialEnd);
        Assert.Equal("complement(<1..>200)", InsdcLocationParser.Format(partial));
        Assert.Throws<FormatException>(() => InsdcLocationParser.Parse("join(1..10"));
    }

    [Fact]
    public void Writer_RoundTripsRecord()
    {
        var r = Sample();
        string text = GenBankWriter.ToString(r);
        Assert.StartsWith("LOCUS       TEST0001                 700 bp    DNA     linear   BCT 15-MAR-2024\n", text);
        Assert.Contains("     CDS             complement(join(350..400,449..520))\n", text);
        Assert.Contains("     mRNA            join(<1..60,100..>150)\n", text);
        Assert.Contains("     misc_feature    order(5..8,530..540)\n", text);
        Assert.Contains("                     /codon_start=1\n", text);
        Assert.Contains("                     /pseudo\n", text);
        Assert.Contains("        1 gtaggggcaa tggctaaaga", text);
        Assert.All(text.Split('\n'), l => Assert.True(l.Length <= 80, l));
        var back = GenBankReader.Parse(text).Single();
        Assert.Equal(r, back);
        Assert.Equal(r.Id, back.Id);
        Assert.Equal(r.Description, back.Description);
        Assert.Equal(r.Features.Count, back.Features.Count);
        for (int i = 0; i < r.Features.Count; i++)
        {
            Assert.Equal(r.Features[i].Location.Segments, back.Features[i].Location.Segments);
            Assert.Equal(r.Features[i].Qualifiers.Keys, back.Features[i].Qualifiers.Keys);
            foreach (var k in r.Features[i].Qualifiers.Keys) Assert.Equal(r.Features[i][k], back.Features[i][k]);
        }
        Assert.Equal(r.Annotations[RecordKeys.Taxonomy], back.Annotations[RecordKeys.Taxonomy]);
        Assert.Equal(r.Annotations[RecordKeys.Comment], back.Annotations[RecordKeys.Comment]);
    }

    [Fact]
    public void ProteinLocus_IsProtein()
    {
        const string gp = "LOCUS       NP_000509                147 aa            linear   PRI 26-AUG-2018\nDEFINITION  hemoglobin subunit beta [Homo sapiens].\nACCESSION   NP_000509\nVERSION     NP_000509.1\nORIGIN\n        1 mvhltpeeks avtalwgkvn\n//\n";
        var r = GenBankReader.Parse(gp).Single();
        Assert.Equal(SequenceType.Protein, r.Type);
        Assert.Equal("PRI", r.Annotations[RecordKeys.Division]);
        Assert.False(r.Annotations.ContainsKey(RecordKeys.MoleculeType));
        Assert.Equal("MVHLTPEEKSAVTALWGKVN", r.ToString()); // ORIGIN is lower case in the file; Biopython upper-cases it
    }
}

public class EmblTests
{
    // data/sample.embl was written by Biopython 1.86 from data/sample.gb.
    [Fact]
    public void Read_BiopythonEmbl_MatchesGenBankRecord()
    {
        var embl = EmblReader.Read(TestData.Path("sample.embl")).Single();
        var gb = GenBankReader.Read(TestData.Path("sample.gb")).Single();
        Assert.Equal("TEST0001.2", embl.Id);
        Assert.Equal(gb, embl);
        Assert.Equal("linear", embl.Annotations[RecordKeys.Topology]);
        Assert.Equal("PRO", embl.Annotations[RecordKeys.Division]);
        Assert.Equal("Escherichia coli synthetic test locus, abcA gene, complete cds; and xyzB pseudogene", embl.Description);
        Assert.Equal(gb.Features.Count, embl.Features.Count);
        Assert.Equal(gb.Features[3].Location.Segments, embl.Features[3].Location.Segments);
        Assert.Equal(gb.Features[2]["translation"], embl.Features[2]["translation"]);
        Assert.Equal("A synthetic GenBank record for parser tests", embl.Annotations[RecordKeys.Reference(1, "title")]);
    }

    [Fact]
    public void Writer_RoundTrips_AndFormatsSequenceBlock()
    {
        var gb = GenBankReader.Read(TestData.Path("sample.gb")).Single();
        string text = EmblWriter.ToString(gb);
        Assert.StartsWith("ID   TEST0001; SV 2; linear; DNA; STD; BCT; 700 BP.\n", text);
        Assert.Contains("FT   CDS             complement(join(350..400,449..520))\n", text);
        // Same block layout as the Biopython-written file: 5 spaces, 6×10 bases, count right-aligned to column 80.
        string expectedLine = File.ReadAllLines(TestData.Path("sample.embl")).First(l => l.EndsWith("       540"));
        Assert.Contains(expectedLine + "\n", text);
        var back = EmblReader.Parse(text).Single();
        Assert.Equal(gb, back);
        Assert.Equal(gb.Features.Count, back.Features.Count);
        Assert.Equal("TEST0001.2", back.Id);
    }

    [Fact]
    public void Read_OldStyleIdLine()
    {
        const string old = "ID   AA03518    standard; DNA; FUN; 12 BP.\nAC   U03518;\nDE   Aspergillus awamori internal transcribed spacer\nSQ   Sequence 12 BP;\n     aacctgcgga aggatc                                                    12\n//\n";
        var r = EmblReader.Parse(old).Single();
        Assert.Equal("U03518", r.Id);
        Assert.Equal("aacctgcggaaggatc".Length, r.Length);
        Assert.Equal("FUN", r.Annotations[RecordKeys.Division]);
    }
}

public class PdbTests
{
    // data/sample.pdb; Biopython 1.86 gives SEQRES chain A "MAMGSKWH" and atom chain A "MAMGXXWH".
    [Fact]
    public void Seqres_MapsModifiedResidues()
    {
        var chains = PdbSequenceReader.ReadSeqres(TestData.Path("sample.pdb"));
        Assert.Equal(2, chains.Count);
        Assert.Equal("9XYZ:A", chains[0].Id);
        Assert.Equal("MAMGSKWH", chains[0].ToString());
        Assert.Equal(SequenceType.Protein, chains[0].Type);
        Assert.Equal("ACGT", chains[1].ToString());
        Assert.Equal(SequenceType.Dna, chains[1].Type);
    }

    [Fact]
    public void AtomSequences_FillGaps_SkipWaterAndLigands()
    {
        var chains = PdbSequenceReader.ReadAtomSequences(TestData.Path("sample.pdb"));
        Assert.Equal("MAMGXXWH", chains[0].ToString());
        Assert.Equal("ACGT", chains[1].ToString());
        Assert.Equal("MAMGWH", PdbSequenceReader.ReadAtomSequences(TestData.Path("sample.pdb"), fillGaps: false)[0].ToString());
    }

    [Theory]
    [InlineData("MSE", 'M')]
    [InlineData("SEP", 'S')]
    [InlineData("TPO", 'T')]
    [InlineData("PTR", 'Y')]
    [InlineData("SEC", 'U')]
    [InlineData("hyp", 'P')]
    [InlineData("XYZ", 'X')]
    public void ThreeToOne(string three, char one) => Assert.Equal(one, PdbSequenceReader.ToOneLetter(three));
}

public class SequenceFileTests
{
    [Theory]
    [InlineData("sample.gb", SequenceFormat.GenBank)]
    [InlineData("sample.embl", SequenceFormat.Embl)]
    [InlineData("sample.pdb", SequenceFormat.Pdb)]
    public void DetectsByExtension(string file, SequenceFormat expected)
        => Assert.Equal(expected, SequenceFile.DetectFormat(TestData.Path(file)));

    [Fact]
    public void DetectsByContent_WhenExtensionUnknown()
    {
        string gb = TestData.TempFile("record.txt");
        File.Copy(TestData.Path("sample.gb"), gb);
        Assert.Equal(SequenceFormat.GenBank, SequenceFile.DetectFormat(gb));
        Assert.Equal(700, SequenceFile.ReadFirst(gb).Length);
        string fq = TestData.TempFile("reads.dat");
        File.WriteAllText(fq, "@r\nAC\n+\nII\n");
        Assert.Equal(SequenceFormat.Fastq, SequenceFile.DetectFormat(fq));
        Assert.Equal("AC", SequenceFile.ReadFirst(fq).ToString());
    }

    [Fact]
    public void Write_ConvertsGenBankToFastaAndEmbl()
    {
        var recs = SequenceFile.Read(TestData.Path("sample.gb")).ToList();
        string fa = TestData.TempFile("out.fasta");
        string embl = TestData.TempFile("out.embl.gz");
        SequenceFile.Write(fa, recs);
        SequenceFile.Write(embl, recs);
        Assert.Equal(recs[0], SequenceFile.ReadFirst(fa));
        Assert.Equal(recs[0], SequenceFile.ReadFirst(embl));
        Assert.Equal(8, SequenceFile.ReadFirst(embl).Features.Count);
        Assert.Throws<NotSupportedException>(() => SequenceFile.Write(TestData.TempFile("x.pdb"), recs));
    }
}
