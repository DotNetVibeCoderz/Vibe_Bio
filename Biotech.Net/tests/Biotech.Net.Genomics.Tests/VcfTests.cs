using System.IO.Compression;
using Biotech.Net.Core;
using Biotech.Net.Genomics.Vcf;

namespace Biotech.Net.Genomics.Tests;

public class VcfTests
{
    private static string SpecPath => Path.Combine(AppContext.BaseDirectory, "data", "vcf42_spec_example.vcf");

    [Fact]
    public void ReadsTheVcf42SpecificationExample()
    {
        using var reader = VcfReader.Open(SpecPath);
        var h = reader.Header;
        Assert.Equal("VCFv4.2", h.FileFormat);
        Assert.Equal(["NA00001", "NA00002", "NA00003"], h.Samples);
        Assert.Equal(6, h.Info.Count);
        Assert.Equal("A", h.GetInfo("AF")!.Number);
        Assert.Equal("Flag", h.GetInfo("DB")!.Type);
        Assert.Equal("dbSNP membership, build 129", h.GetInfo("DB")!.Description);
        Assert.Equal(4, h.Format.Count);
        Assert.Equal(2, h.Filters.Count);
        Assert.Equal(new VcfContig("20", 62435964), h.Contigs.Single());
        Assert.Contains("phasing=partial", h.OtherLines);

        var recs = reader.Records().ToList();
        Assert.Equal(5, recs.Count);

        var r0 = recs[0];
        Assert.Equal("20", r0.Chromosome);
        Assert.Equal(14370, r0.Position);
        Assert.Equal("rs6054257", r0.Variant.Id);
        Assert.Equal(29, r0.Quality);
        Assert.True(r0.Variant.IsPass);
        Assert.Equal("", r0.Info["DB"]); // flag
        Assert.Equal("14", r0.Info["DP"]);
        Assert.True(r0.GetGenotype(0).Phased);
        Assert.True(r0.GetGenotype(0).IsHomozygousReference);
        Assert.True(r0.GetGenotype(1).IsHeterozygous);
        Assert.True(r0.GetGenotype(2).IsHomozygousAlternate);
        Assert.False(r0.GetGenotype(2).Phased);
        Assert.Equal(48, r0.GenotypeQuality(0));
        Assert.Equal(8, r0.Depth(1));
        Assert.Equal("51,51", r0.GetField(0, "HQ"));
        Assert.Equal(".,.", r0.GetField(2, "HQ")); // kept verbatim

        Assert.Equal("q10", recs[1].Filter);
        Assert.False(recs[1].Variant.IsPass);
        Assert.Null(recs[1].GetField(2, "HQ")); // trailing field dropped

        var multi = recs[2];
        Assert.Equal(["G", "T"], multi.Alternates);
        Assert.Equal([1, 2], multi.GetGenotype(0).Alleles);
        Assert.Equal("0.333,0.667", multi.Info["AF"]);

        Assert.Empty(recs[3].Alternates); // ALT "."
        Assert.Equal(VariantType.Deletion, Variant.Classify("GTC", "G"));
        Assert.Equal(VariantType.Insertion, Variant.Classify("GTC", "GTCT"));
        Assert.Equal([0, 2], recs[4].GetGenotype(1).Alleles);
    }

    [Fact]
    public void RoundTripsThroughTheWriter()
    {
        var text = File.ReadAllText(SpecPath).Replace("\r\n", "\n");
        var (header, records) = VcfReader.Parse(text);
        string written = VcfWriter.ToString(header, records);
        var (header2, records2) = VcfReader.Parse(written);
        Assert.Equal(header.Samples, header2.Samples);
        Assert.Equal(header.Info.Select(i => i.Id), header2.Info.Select(i => i.Id));
        Assert.Equal(records.Count, records2.Count);
        for (int i = 0; i < records.Count; i++)
            Assert.Equal(records[i].ToString(), records2[i].ToString());
        // Data lines are byte-identical to the specification example.
        var dataIn = text.Split('\n').Where(l => l.Length > 0 && l[0] != '#');
        var dataOut = written.Split('\n').Where(l => l.Length > 0 && l[0] != '#');
        Assert.Equal(dataIn, dataOut);
    }

    [Fact]
    public void ReadsGzipInputAndWritesGzipOutput()
    {
        string dir = Path.Combine(Path.GetTempPath(), "biotech-genomics-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            string gz = Path.Combine(dir, "spec.vcf.gz");
            using (var src = File.OpenRead(SpecPath))
            using (var dst = new GZipStream(File.Create(gz), CompressionLevel.Fastest))
                src.CopyTo(dst);
            var recs = VcfReader.ReadRecords(gz).ToList();
            Assert.Equal(5, recs.Count);

            using var r = VcfReader.Open(SpecPath);
            string outGz = Path.Combine(dir, "out.vcf.gz");
            VcfWriter.Write(outGz, r.Header, r.Records());
            Assert.Equal(0x1f, File.ReadAllBytes(outGz)[0]);
            Assert.Equal(5, VcfReader.ReadRecords(outGz).Count());
        }
        finally { Directory.Delete(dir, true); }
    }

    [Fact]
    public void FiltersByQualityPassAndRegion()
    {
        var (_, recs) = VcfReader.Parse(File.ReadAllText(SpecPath));
        Assert.Equal(4, recs.MinQuality(20).Count());
        Assert.Equal(4, recs.PassOnly().Count());
        // 0-based half-open [1110695, 1230237) covers POS 1110696 and 1230237 (0-based 1230236).
        var region = new GenomicInterval("20", 1110695, 1230237);
        Assert.Equal([1110696L, 1230237L], recs.InRegion(region).Select(r => r.Position));
        Assert.Empty(recs.InRegion(new GenomicInterval("21", 0, 2_000_000)));
        // The deletion GTC at 1234567 spans 0-based [1234566, 1234569).
        Assert.Single(recs.InRegion(new GenomicInterval("20", 1234568, 1234570)));
        Assert.Equal(3, recs.InRegions([new GenomicInterval("20", 0, 20000), new GenomicInterval("20", 1234560, 1234570)]).Count());

        using var reader = VcfReader.Open(SpecPath);
        Assert.Equal([17330L], reader.Query(new GenomicInterval("20", 17329, 17330)).Select(r => r.Position));
        Assert.Equal(2, recs.BiallelicSnvs().Count());
    }

    [Fact]
    public void ParsesGenotypeEdgeCases()
    {
        Assert.Equal([-1, -1], VcfReader.ParseGenotype("./.").Alleles);
        Assert.True(VcfReader.ParseGenotype("./.").IsMissing);
        Assert.Equal([2], VcfReader.ParseGenotype("2").Alleles);
        Assert.Equal([0, 10], VcfReader.ParseGenotype("0|10").Alleles);
        Assert.True(VcfReader.ParseGenotype("0|10").Phased);
        Assert.Equal([1, -1], VcfReader.ParseGenotype("1/.").Alleles);
        Assert.Equal([0, 1, 1], VcfReader.ParseGenotype("0/1/1").Alleles);
    }

    [Fact]
    public void BuildsRecordsAndSplitsMultiallelicSites()
    {
        var v = new Variant("chr1", 100, "A", ["C", "G"]) { Quality = 50.5 };
        v.Genotypes.Add(new Genotype([1, 2]));
        v.Genotypes.Add(new Genotype([0, 2], Phased: true));
        var rec = new VcfRecord(v);
        rec.SetField(0, "AD", [1, 5, 6]);
        rec.SetField(1, "AD", [7, 0, 8]);
        rec.SetField(0, "PL", [100, 50, 90, 40, 0, 60]); // order: 00, 01, 11, 02, 12, 22
        rec.SetField(1, "PL", [80, 70, 99, 0, 50, 60]);
        Assert.Equal("chr1\t100\t.\tA\tC,G\t50.5\tPASS\t.\tGT:AD:PL\t1/2:1,5,6:100,50,90,40,0,60\t0|2:7,0,8:80,70,99,0,50,60", rec.ToString());

        Assert.Equal(4, VcfRecord.GenotypeIndex(1, 2));
        var split = rec.SplitMultiallelic().ToList();
        Assert.Equal(2, split.Count);
        Assert.Equal("C", split[0].Alternates.Single());
        Assert.Equal([1, 0], split[0].GetGenotype(0).Alleles);
        Assert.Equal([1, 5], split[0].AlleleDepths(0)!);
        Assert.Equal([100, 50, 90], split[0].PhredLikelihoods(0)!);
        Assert.Equal([0, 1], split[1].GetGenotype(1).Alleles);
        Assert.True(split[1].GetGenotype(1).Phased);
        Assert.Equal([7, 8], split[1].AlleleDepths(1)!);
        Assert.Equal([80, 0, 60], split[1].PhredLikelihoods(1)!);
    }

    [Fact]
    public void HeaderParsesQuotedDescriptionsWithCommas()
    {
        var h = new VcfHeader();
        h.ParseMetaLine("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence, with \\\"quotes\\\", and commas\">");
        Assert.Equal("Consequence, with \"quotes\", and commas", h.GetInfo("CSQ")!.Description);
        var round = new VcfHeader();
        foreach (var line in h.ToString().Split('\n').Where(l => l.StartsWith("##INFO"))) round.ParseMetaLine(line);
        Assert.Equal(h.GetInfo("CSQ"), round.GetInfo("CSQ"));
    }
}
