using Biotech.Net.Core;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Genomics.Vcf;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

public class TabixTests
{
    public static TheoryData<string, string> Files => new()
    {
        { "variants.vcf.gz", "variants.vcf.gz.tbi" },
        { "variants.vcf.gz", "variants.vcf.gz.csi" },
        { "features.bed.gz", "features.bed.gz.tbi" },
        { "genes.gff3.gz", "genes.gff3.gz.tbi" },
    };

    private static string ExpectedFor(string file) => file.Replace(".gz", "") + ".regions.tsv";

    [Theory]
    [MemberData(nameof(Files))]
    public void QueriesMatchTheTabixCommand(string file, string index)
    {
        using var reader = TabixReader.Open(Hts(file), Hts(index));
        foreach (var (region, expected) in ReadRegions(Hts(ExpectedFor(file))))
        {
            var actual = reader.Query(region).Select(l => Columns(l)).ToList();
            Assert.True(expected.SequenceEqual(actual), $"{file} {region}: expected {expected.Count} lines, got {actual.Count}");
        }
    }

    private static void AssertSameIndex(TabixIndex e, TabixIndex a)
    {
        Assert.Equal(e.Config, a.Config);
        Assert.Equal(e.Names, a.Names);
        Assert.Equal(e.Index.MinShift, a.Index.MinShift);
        Assert.Equal(e.Index.Depth, a.Index.Depth);
        Assert.Equal(e.Index.References.Count, a.Index.References.Count);
        for (int i = 0; i < e.Index.References.Count; i++)
        {
            var x = e.Index.References[i];
            var y = a.Index.References[i];
            Assert.Equal(x.Bins.Keys.Order(), y.Bins.Keys.Order());
            foreach (var (bin, chunks) in x.Bins) Assert.Equal(chunks, y.Bins[bin]);
            Assert.Equal(x.Linear, y.Linear);
            Assert.Equal(x.Span, y.Span);
            Assert.Equal(x.MappedCount, y.MappedCount);
            if (e.Index.Format == BinningIndexFormat.Csi)
                foreach (var (bin, loff) in x.BinLinearOffsets) Assert.Equal(loff, y.BinLinearOffsets[bin]);
        }
    }

    [Fact]
    public void BuildsTheSameIndexesAsTabix()
    {
        AssertSameIndex(TabixIndex.Load(Hts("variants.vcf.gz.tbi")), TabixIndex.Build(Hts("variants.vcf.gz"), TabixConfig.VcfPreset));
        AssertSameIndex(TabixIndex.Load(Hts("genes.gff3.gz.tbi")), TabixIndex.Build(Hts("genes.gff3.gz"), TabixConfig.GffPreset));
        AssertSameIndex(TabixIndex.Load(Hts("features.bed.gz.tbi")), TabixIndex.Build(Hts("features.bed.gz"), TabixConfig.BedPreset with { SkipLines = 1 }));
        var csi = TabixIndex.Load(Hts("variants.vcf.gz.csi"));
        Assert.Equal(BinningIndexFormat.Csi, csi.Index.Format);
        AssertSameIndex(csi, TabixIndex.Build(Hts("variants.vcf.gz"), TabixConfig.VcfPreset, BinningIndexFormat.Csi));
    }

    [Fact]
    public void SavedIndexesReadBack()
    {
        foreach (var fmt in new[] { BinningIndexFormat.Tbi, BinningIndexFormat.Csi })
        {
            var built = TabixIndex.Build(Hts("variants.vcf.gz"), TabixConfig.VcfPreset, fmt);
            var ms = new MemoryStream();
            built.Write(ms);
            ms.Position = 0;
            AssertSameIndex(built, TabixIndex.Read(ms));
        }
    }

    [Fact]
    public void ParsesIntervalsLikeHtslib()
    {
        var vcf = TabixConfig.VcfPreset;
        Assert.True(vcf.TryParse("chr1\t100\t.\tACGT\tA\t50\tPASS\tDP=3"u8, out var name, out long b, out long e));
        Assert.Equal("chr1"u8.ToArray(), name.ToArray());
        Assert.Equal((99L, 103L), (b, e));
        Assert.True(vcf.TryParse("chr1\t100\t.\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=500"u8, out _, out b, out e));
        Assert.Equal((99L, 500L), (b, e));
        Assert.True(vcf.TryParse("chr1\t100\t.\tA\t<DEL>\t50\tPASS\tSVLEN=-300"u8, out _, out b, out e));
        Assert.Equal((99L, 399L), (b, e));
        Assert.True(TabixConfig.BedPreset.TryParse("chr2\t0\t10\tx"u8, out _, out b, out e));
        Assert.Equal((0L, 10L), (b, e));
        Assert.True(TabixConfig.GffPreset.TryParse("chr2\tsrc\tgene\t5\t9\t.\t+\t.\tID=a"u8, out _, out b, out e));
        Assert.Equal((4L, 9L), (b, e));
        Assert.False(TabixConfig.BedPreset.TryParse("track name=x"u8, out _, out _, out _));
    }

    [Fact]
    public void VcfReaderUsesTheTabixIndex()
    {
        using var reader = VcfReader.Open(Hts("variants.vcf.gz"));
        Assert.True(reader.IsIndexed);
        foreach (var (region, expected) in ReadRegions(Hts("variants.vcf.regions.tsv")))
        {
            var actual = reader.Query(region).Select(r => $"{r.Variant.Chromosome}\t{r.Variant.Position}").ToList();
            Assert.Equal(expected.Select(l => Columns(l, 2)), actual);
        }
        // repeatable, and GenomicInterval form (0-based half-open)
        var a = reader.Query(new GenomicInterval("chr1", 4999, 9000)).Count();
        var b = reader.Query("chr1:5000-9000").Count();
        Assert.Equal(a, b);
        Assert.True(a > 0);
    }

    [Fact]
    public void WritesBgzfVcfThatCanBeIndexed()
    {
        var dir = Directory.CreateTempSubdirectory();
        try
        {
            string path = Path.Combine(dir.FullName, "out.vcf.gz");
            List<VcfRecord> records;
            VcfHeader header;
            using (var r = VcfReader.Open(Hts("variants.vcf.gz"))) { header = r.Header; records = r.Records().ToList(); }
            using (var w = VcfWriter.Create(path, header)) foreach (var rec in records) w.Write(rec);
            Assert.True(BgzfReader.IsBgzf(path));
            TabixIndex.Build(path, TabixConfig.VcfPreset).Save(path + ".tbi");
            using var indexed = VcfReader.Open(path);
            Assert.True(indexed.IsIndexed);
            using var orig = VcfReader.Open(Hts("variants.vcf.gz"));
            foreach (var region in new[] { "chr1:5000-9000", "chr2", "chr3:7000-8000" })
                Assert.Equal(orig.Query(region).Select(x => x.Variant.Position), indexed.Query(region).Select(x => x.Variant.Position));
        }
        finally { dir.Delete(true); }
    }
}
