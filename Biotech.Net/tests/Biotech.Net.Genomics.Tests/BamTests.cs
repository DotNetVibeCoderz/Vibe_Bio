using Biotech.Net.Genomics.IO;
using Biotech.Net.Genomics.Mapping;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

public class BamTests
{
    private static void AssertSamLinesEqual(IReadOnlyList<string> expected, IReadOnlyList<string> actual)
    {
        Assert.Equal(expected.Count, actual.Count);
        for (int i = 0; i < expected.Count; i++)
        {
            if (expected[i] == actual[i]) continue;
            var e = expected[i].Split('\t');
            var a = actual[i].Split('\t');
            for (int k = 0; k < Math.Max(e.Length, a.Length); k++)
                Assert.True(k < e.Length && k < a.Length && e[k] == a[k],
                    $"record {i} ({e[0]}) column {k + 1}: expected '{(k < e.Length ? e[k] : "<none>")}' got '{(k < a.Length ? a[k] : "<none>")}'");
        }
    }

    [Fact]
    public void DecodesEveryRecordLikeSamtoolsView()
    {
        var (header, expected) = ReadSam(Hts("reads.bam.sam.gz"));
        using var reader = BamReader.Open(Hts("reads.bam"));
        Assert.Equal(string.Join('\n', header) + "\n", reader.Header.Text);
        Assert.Equal(["chr1", "chr2", "chr3"], reader.Header.References.Select(r => r.Name));
        Assert.Equal([40000L, 25000L, 8000L], reader.Header.References.Select(r => r.Length));
        Assert.Equal("coordinate", reader.Header.SortOrder);
        Assert.Equal(["grpA", "grpB"], reader.Header.ReadGroupIds);
        var actual = reader.Records().Select(r => r.ToString()).ToList();
        AssertSamLinesEqual(expected, actual);
    }

    [Fact]
    public void ExposesTypedFieldsAndFlags()
    {
        using var reader = BamReader.Open(Hts("reads.bam"));
        var all = reader.Records().ToList();
        var tagged = all.First(r => r.HasTag("XB"));
        Assert.Equal('c', tagged.GetTag("XB")!.Value.ArrayType);
        Assert.Equal(-1, tagged.GetTag("XB")!.Value.AsArray<sbyte>()[1]);
        Assert.Equal([0u, 4294967295u], tagged.GetTag("XI")!.Value.AsArray<uint>());
        Assert.Equal([1.5f, -0.25f, 3.14159f, 1e-07f, 123456789f], tagged.GetTag("Xf")!.Value.AsArray<float>());
        Assert.Empty(tagged.GetTag("Xe")!.Value.AsArray<sbyte>());
        Assert.Equal("free text with spaces", tagged.GetString("XZ"));
        Assert.Equal("1AE301", tagged.GetString("XH"));
        Assert.Equal('H', tagged.GetTag("XH")!.Value.Type);
        Assert.Equal(1, tagged.GetString("XA")!.Length);
        Assert.Equal(-32768, tagged.GetTag("Xt")!.Value.AsArray<short>()[0]);

        var spliced = all.Single(r => r.Name == "spliced1");
        Assert.Equal(1000, spliced.Position);
        Assert.Equal(1000 + 80 + 19920 + 80, spliced.End);
        Assert.Equal(160, spliced.QueryLength);
        Assert.Equal(160, spliced.AlignedPairs(matchesOnly: true).Count());

        var mate = all.First(r => r.IsPaired && r.IsFirstOfPair && r.IsProperPair);
        Assert.True(mate.IsMateReverse);
        Assert.Equal(mate.ReferenceName, mate.MateReferenceName);
        Assert.True(mate.TemplateLength > 0);

        var unplaced = all.Where(r => r.ReferenceId < 0).ToList();
        Assert.Equal(12, unplaced.Count);
        Assert.All(unplaced, r => Assert.True(r.IsUnmapped));
        Assert.Contains(unplaced, r => r.Qualities is null);
        var secondary = all.First(r => r.IsSecondary);
        Assert.Empty(secondary.Sequence);
        Assert.Equal("*", secondary.SequenceString);
        var supp = all.First(r => r.IsSupplementary);
        Assert.Equal('H', supp.Cigar[0].Op);
        Assert.NotNull(supp.GetString("SA"));
    }

    [Fact]
    public void ParsesSamTextAndEncodesTheSameBytesAsSamtools()
    {
        // reads.unsorted.bam = samtools view -b reads.sam (record order kept)
        var (_, expected) = SplitBam(Gunzip(File.ReadAllBytes(Hts("reads.unsorted.bam"))));
        using var sam = SamAlignmentReader.Open(Hts("reads.sam.gz"));
        var header = sam.Header;
        var buf = new ByteBuffer();
        var prefix = new byte[4];
        int i = 0;
        foreach (var rec in sam.Records())
        {
            buf.Clear();
            buf.Append(prefix);
            BamCodec.Encode(rec, header, buf, out _, out _);
            buf.SetInt32(0, buf.Length - 4);
            Assert.True(expected[i].AsSpan().SequenceEqual(buf.Span), $"record {i} ({rec.Name}) encodes differently from htslib");
            i++;
        }
        Assert.Equal(expected.Count, i);
    }

    [Fact]
    public void SamTextRoundTripsThroughTheRecordModel()
    {
        var lines = ReadText(Hts("reads.sam.gz")).Split('\n').Where(l => l.Length > 0 && l[0] != '@').ToList();
        using var bam = new BamReader(File.OpenRead(Hts("reads.unsorted.bam")));
        var fromBam = bam.Records().Select(r => r.ToString()).ToList();
        var header = SamAlignmentReader.Open(Hts("reads.sam.gz")).Header;
        var fromText = lines.Select(l => AlignmentRecord.Parse(l, header).ToString()).ToList();
        // the BAM side upper-cases bases (4-bit encoding); everything else must agree
        static string UpperSeq(string l) { var c = l.Split('\t'); c[9] = c[9].ToUpperInvariant(); return string.Join('\t', c); }
        AssertSamLinesEqual(fromBam, fromText.Select(UpperSeq).ToList());
        var lower = AlignmentRecord.Parse(lines.Single(l => l.StartsWith("lowercase", StringComparison.Ordinal)), header);
        Assert.True(lower.Sequence.All(b => b is >= (byte)'a' and <= (byte)'z')); // case kept, as Sequence does
    }

    [Fact]
    public void WriterReproducesTheInputBamPayload()
    {
        var original = Gunzip(File.ReadAllBytes(Hts("reads.bam")));
        var ms = new MemoryStream();
        using (var reader = BamReader.Open(Hts("reads.bam")))
        using (var writer = new BamWriter(ms, reader.Header, leaveOpen: true))
            writer.Write(reader.Records());
        var written = ms.ToArray();
        Assert.Equal(original, Gunzip(written));
        Assert.Equal(BgzfWriter.EofMarker.ToArray(), written[^28..]);

        // and reads back record by record
        using var again = new BamReader(new MemoryStream(written));
        using var orig = BamReader.Open(Hts("reads.bam"));
        Assert.Equal(orig.Records().Select(r => r.ToString()), again.Records().Select(r => r.ToString()));
    }

    private static void AssertSameIndex(BinningIndex expected, BinningIndex actual)
    {
        Assert.Equal(expected.MinShift, actual.MinShift);
        Assert.Equal(expected.Depth, actual.Depth);
        Assert.Equal(expected.References.Count, actual.References.Count);
        Assert.Equal(expected.UnplacedCount, actual.UnplacedCount);
        for (int i = 0; i < expected.References.Count; i++)
        {
            var e = expected.References[i];
            var a = actual.References[i];
            Assert.Equal(e.Bins.Keys.Order(), a.Bins.Keys.Order());
            foreach (var (bin, chunks) in e.Bins) Assert.Equal(chunks, a.Bins[bin]);
            Assert.Equal(e.Linear, a.Linear);
            Assert.Equal(e.Span, a.Span);
            Assert.Equal(e.MappedCount, a.MappedCount);
            Assert.Equal(e.UnmappedCount, a.UnmappedCount);
            if (expected.Format == BinningIndexFormat.Csi)
                foreach (var (bin, loff) in e.BinLinearOffsets) Assert.Equal(loff, a.BinLinearOffsets[bin]);
        }
    }

    [Fact]
    public void BuildsTheSameBaiAndCsiAsSamtoolsIndex()
    {
        var bai = BinningIndex.Load(Hts("reads.bam.bai"));
        Assert.Equal(BinningIndexFormat.Bai, bai.Format);
        using (var fs = File.OpenRead(Hts("reads.bam"))) AssertSameIndex(bai, BamIndexer.Build(fs));

        var csi = BinningIndex.Load(Hts("reads.bam.csi"));
        Assert.Equal(BinningIndexFormat.Csi, csi.Format);
        using (var fs = File.OpenRead(Hts("reads.bam"))) AssertSameIndex(csi, BamIndexer.Build(fs, BinningIndexFormat.Csi, csi.MinShift, csi.Depth));

        // our serialisation reads back to the same index
        foreach (var idx in new[] { bai, csi })
        {
            var ms = new MemoryStream();
            idx.Write(ms);
            ms.Position = 0;
            AssertSameIndex(idx, BinningIndex.Read(ms));
        }
    }

    [Fact]
    public void IndexMetadataMatchesSamtoolsIdxstats()
    {
        var bai = BinningIndex.Load(Hts("reads.bam.bai"));
        var rows = File.ReadAllLines(Hts("reads.idxstats.tsv")).Select(l => l.Split('\t')).ToList();
        for (int i = 0; i < bai.References.Count; i++)
        {
            Assert.Equal(ulong.Parse(rows[i][2]), bai.References[i].MappedCount);
            Assert.Equal(ulong.Parse(rows[i][3]), bai.References[i].UnmappedCount);
        }
        Assert.Equal(ulong.Parse(rows[^1][3]), bai.UnplacedCount);
    }

    [Theory]
    [InlineData("reads.bam.bai")]
    [InlineData("reads.bam.csi")]
    [InlineData(null)]
    public void RegionQueriesMatchSamtoolsView(string? indexFile)
    {
        using var reader = BamReader.Open(Hts("reads.bam"), indexFile is null ? null : Hts(indexFile));
        if (indexFile is null)
            using (var fs = File.OpenRead(Hts("reads.bam"))) reader.Index = BamIndexer.Build(fs);
        foreach (var (region, expected) in ReadRegions(Hts("reads.regions.tsv")))
        {
            var actual = reader.Query(region).Select(r => Columns(r.ToString())).ToList();
            Assert.True(expected.SequenceEqual(actual), $"region {region}: expected {expected.Count} records, got {actual.Count}");
        }
    }

    [Fact]
    public void QueriesUnplacedReads()
    {
        using var reader = BamReader.Open(Hts("reads.bam"));
        var expected = reader.Records().Where(r => r.ReferenceId < 0).Select(r => r.Name).ToList();
        Assert.Equal(expected, reader.QueryUnmapped().Select(r => r.Name));
    }

    [Fact]
    public void WriterBuildsAnIndexOnTheFly()
    {
        var dir = Directory.CreateTempSubdirectory();
        try
        {
            string path = Path.Combine(dir.FullName, "copy.bam");
            using (var reader = BamReader.Open(Hts("reads.bam")))
            using (var writer = BamWriter.Create(path, reader.Header, new BamWriterOptions { Index = BinningIndexFormat.Bai }))
                writer.Write(reader.Records());
            Assert.True(File.Exists(path + ".bai"));
            var built = BinningIndex.Load(path + ".bai");
            using (var fs = File.OpenRead(path)) AssertSameIndex(BamIndexer.Build(fs), built);
            using var copy = BamReader.Open(path);
            using var orig = BamReader.Open(Hts("reads.bam"));
            foreach (var region in new[] { "chr1:500-1500", "chr2:12000-14000", "chr1:16000-33000" })
                Assert.Equal(orig.Query(region).Select(r => r.ToString()), copy.Query(region).Select(r => r.ToString()));
        }
        finally { dir.Delete(true); }
    }

    [Fact]
    public void OpensFilesByContent()
    {
        using (var a = AlignmentFile.Open(Hts("reads.bam"))) Assert.IsType<BamReader>(a);
        using (var b = AlignmentFile.Open(Hts("reads.sam.gz"))) Assert.IsType<SamAlignmentReader>(b);
    }

    [Theory]
    [InlineData(79.3, "79.3")]
    [InlineData(123456789.0, "1.23457e+08")]
    [InlineData(1e-07, "1e-07")]
    [InlineData(0.0001, "0.0001")]
    [InlineData(0.00001234, "1.234e-05")]
    [InlineData(100000.0, "100000")]
    [InlineData(1000000.0, "1e+06")]
    [InlineData(999999.5, "1e+06")]
    [InlineData(-0.25, "-0.25")]
    [InlineData(3.14159, "3.14159")]
    [InlineData(0.0, "0")]
    [InlineData(2.5e15, "2.5e+15")]
    public void FormatsFloatsLikePrintfG(double value, string expected) => Assert.Equal(expected, SamText.FormatG(value));

    [Fact]
    public void TagEditingKeepsOrder()
    {
        var r = AlignmentRecord.Parse("r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:1\tMD:Z:5\tXX:i:-3");
        Assert.Equal(['C', 'Z', 'c'], r.Tags.Select(t => t.Type));
        r.SetTag(AuxField.Int("NM", 70000));
        Assert.Equal(["NM:i:70000", "MD:Z:5", "XX:i:-3"], r.Tags.Select(t => t.ToString()));
        Assert.True(r.RemoveTag("MD"));
        Assert.False(r.HasTag("MD"));
        r.SetTag(AuxField.String("RG", "g1"));
        Assert.Equal("r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:70000\tXX:i:-3\tRG:Z:g1", r.ToString());
        Assert.Equal(9, r.Position);
        Assert.Equal(14, r.End);
        Assert.Equal([(0, 9), (1, 10), (2, 11), (3, 12), (4, 13)], r.AlignedPairs().Select(p => (p.QueryPosition!.Value, p.ReferencePosition!.Value)));
    }
}
