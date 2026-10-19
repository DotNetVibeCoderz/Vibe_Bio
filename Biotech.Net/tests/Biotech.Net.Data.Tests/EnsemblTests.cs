using System.Text.Json;
using Biotech.Net.Core;
using Biotech.Net.Data.Ensembl;
using Biotech.Net.Data.Http;
using Biotech.Net.Genomics.Annotation;

namespace Biotech.Net.Data.Tests;

public class EnsemblTests
{
    private readonly FakeHttpHandler _handler = new();
    private readonly FakeClock _clock = new();

    public EnsemblTests()
    {
        var P = FakeHttpHandler.Path;
        _handler
            .OnFixture(P("/lookup/symbol/homo_sapiens/BRCA2"), "ensembl/lookup-symbol-BRCA2.json")
            .OnFixture(P("/lookup/id/ENSG00000244734", q => q.ContainsKey("expand")), "ensembl/lookup-id-HBB-expand.json")
            .OnFixture(P("/sequence/id/ENST00000335295", q => q["type"][0] == "protein"), "ensembl/sequence-protein.json")
            .OnFixture(P("/sequence/id/ENST00000335295", q => q["type"][0] == "cds"), "ensembl/sequence-cds.json")
            .OnFixture(P("/sequence/region/human/11:5227002..5227071:-1"), "ensembl/sequence-region.json")
            .OnFixture(P("/overlap/region/human/11:5225000-5230000"), "ensembl/overlap-region.json")
            .OnFixture(r => r.RequestUri!.AbsolutePath.StartsWith("/vep/human/hgvs/", StringComparison.Ordinal), "ensembl/vep-hgvs.json")
            .OnFixture(P("/vep/human/region/11:5227002-5227002:1/A"), "ensembl/vep-region.json")
            .OnFixture(P("/xrefs/id/ENSG00000244734"), "ensembl/xrefs-id.json")
            .OnFixture(P("/xrefs/symbol/homo_sapiens/HBB"), "ensembl/xrefs-symbol.json");
    }

    private EnsemblClient Client() => new(_handler.Client(), TestOptions.Create(_clock));

    [Fact]
    public async Task LookupSymbol_ParsesTheGene_AndConvertsCoordinates()
    {
        var brca2 = await Client().LookupSymbolAsync("homo_sapiens", "BRCA2");
        Assert.Equal("ENSG00000139618", brca2.Id);
        Assert.Equal("ENSG00000139618.19", brca2.VersionedId);
        Assert.Equal("Gene", brca2.ObjectType);
        Assert.Equal("13", brca2.SeqRegionName);
        Assert.Equal(32315086, brca2.Start);
        Assert.Equal(32400268, brca2.End);
        Assert.Equal("GRCh38", brca2.AssemblyName);
        var interval = brca2.ToInterval();
        Assert.Equal(new GenomicInterval("13", 32315085, 32400268, Strand.Forward), interval);
        var req = _handler.Requests.Single();
        Assert.Equal("application/json", req.Param("content-type"));
        Assert.Null(req.Param("expand"));
    }

    [Fact]
    public async Task LookupExpanded_BuildsACoreGeneWithTranscriptsAndCds()
    {
        var hbb = await Client().LookupIdAsync("ENSG00000244734", expand: true);
        Assert.Equal("HBB", hbb.DisplayName);
        Assert.Equal(-1, hbb.Strand);
        Assert.Equal(19, hbb.Transcripts.Count);
        var canonical = hbb.Transcripts.Single(t => t.IsCanonical);
        Assert.Equal("ENST00000335295", canonical.Id);
        Assert.Equal("ENSP00000333994", canonical.Translation!.Id);
        Assert.Equal(3, canonical.Exons.Count);

        Gene gene = hbb.ToGene();
        Assert.Equal("HBB", gene.Name);
        Assert.Equal(Strand.Reverse, gene.Strand);
        Assert.Equal("protein_coding", gene.Biotype);
        var t = gene.Transcripts.Single(x => x.Id == "ENST00000335295");
        Assert.Equal(628, t.SplicedLength); // Ensembl reports length 628
        Assert.Equal(new GenomicInterval("11", 5225597, 5227021, Strand.Reverse), t.CodingRegion);
        Assert.Equal("1", _handler.Requests.Single().Param("expand"));
    }

    [Fact]
    public async Task Sequences_CdsTranslatesToTheProtein()
    {
        var c = Client();
        var protein = await c.GetSequenceAsync("ENST00000335295", EnsemblSequenceType.Protein);
        var cds = await c.GetSequenceAsync("ENST00000335295", EnsemblSequenceType.Cds);
        Assert.Equal("ENSP00000333994", protein.Id);
        Assert.Equal(SequenceType.Protein, protein.Type);
        Assert.Equal(147, protein.Length);
        Assert.Equal(SequenceType.Dna, cds.Type);
        Assert.Equal(protein.ToString(), cds.Translate(toStop: true).ToString());
    }

    [Fact]
    public async Task RegionSequence_UsesOneBasedEnsemblNotation()
    {
        // 0-based half-open [5227001, 5227071) on the reverse strand = Ensembl 11:5227002..5227071:-1.
        var region = new GenomicInterval("11", 5227001, 5227071, Strand.Reverse);
        var s = await Client().GetRegionSequenceAsync("human", region);
        Assert.Equal(70, s.Length);
        Assert.Equal("chromosome:GRCh38:11:5227002:5227071:-1", s.Id);
        Assert.EndsWith("ATGGTGCATCTGACTCCTGA", s.ToString()); // HBB start codon, read on the minus strand
        Assert.Equal("/sequence/region/human/11:5227002..5227071:-1", _handler.Requests.Single().Path);
    }

    [Fact]
    public async Task Overlap_ReturnsGenesAndTranscripts()
    {
        var features = await Client().OverlapRegionAsync("human", GenomicInterval.Parse("11:5225000-5230000"), ["gene", "transcript"]);
        Assert.Equal(24, features.Count);
        Assert.Contains(features, f => f is { FeatureType: "gene", Id: "ENSG00000244734", ExternalName: "HBB" });
        Assert.Contains(features, f => f is { FeatureType: "transcript", Id: "ENST00000335295", IsCanonical: true });
        var req = _handler.Requests.Single();
        Assert.Equal(["gene", "transcript"], req.Query["feature"]);
    }

    [Fact]
    public async Task VepHgvs_ReportsTheSickleCellMissense()
    {
        var results = await Client().VepHgvsAsync("human", "ENST00000335295:c.20A>T");
        var r = Assert.Single(results);
        Assert.Equal("missense_variant", r.MostSevereConsequence);
        Assert.Equal(Consequence.MissenseVariant, r.MostSevere);
        Assert.Equal(5227002, r.Start);
        Assert.Equal("A/T", r.AlleleString);
        var hbb = r.TranscriptConsequences.Single(c => c.TranscriptId == "ENST00000335295");
        Assert.Equal("HBB", hbb.GeneSymbol);
        Assert.Equal("E/V", hbb.AminoAcids);
        Assert.Equal("E7V", hbb.ProteinChange);
        Assert.Equal("MODERATE", hbb.Impact);
        Assert.Equal([Consequence.MissenseVariant], hbb.Consequences);
        Assert.Contains("CM097155", r.ColocatedVariantIds);
        Assert.EndsWith("/vep/human/hgvs/ENST00000335295:c.20A%3ET", _handler.Requests.Single().Uri.AbsoluteUri.Split('?')[0]);
    }

    [Fact]
    public async Task VepRegion_FromACoreVariant()
    {
        var variant = new Variant("11", 5227002, "T", "A", "rs334");
        var r = Assert.Single(await Client().VepRegionAsync("human", variant));
        Assert.Equal("missense_variant", r.MostSevereConsequence);
        Assert.Equal("/vep/human/region/11:5227002-5227002:1/A", _handler.Requests.Single().Path);
    }

    [Theory]
    [InlineData("1", 100, "A", "G", "1:100-100:1/G")]          // SNV
    [InlineData("1", 100, "AT", "A", "1:101-101:1/-")]         // 1-bp deletion after the anchor
    [InlineData("1", 100, "A", "ACG", "1:101-100:1/CG")]       // insertion between 100 and 101
    [InlineData("X", 5, "ACG", "TTT", "X:5-7:1/TTT")]          // block substitution
    public void ToVepRegion_ConvertsVcfAnchoredAlleles(string chrom, long pos, string reference, string alt, string expected) =>
        Assert.Equal(expected, EnsemblClient.ToVepRegion(chrom, pos, reference, alt));

    [Fact]
    public async Task VepRegions_PostsVcfLinesInBatches()
    {
        var handler = new FakeHttpHandler().On(FakeHttpHandler.Path("/vep/human/region"), () => Fixtures.Response("ensembl/vep-region.json"));
        var c = new EnsemblClient(handler.Client(), TestOptions.Create(_clock));
        var variants = Enumerable.Range(0, 450).Select(i => new Variant("11", 5227002 + i, "T", "A")).ToList();
        var results = await c.VepRegionsAsync("human", variants);
        Assert.Equal(3, handler.Count); // 200 + 200 + 50
        Assert.Equal(3, results.Count); // one result per fixture response
        var first = JsonDocument.Parse(handler.Requests.First().Body!).RootElement.GetProperty("variants");
        Assert.Equal(200, first.GetArrayLength());
        Assert.Equal("11 5227002 . T A . . .", first[0].GetString());
    }

    [Fact]
    public async Task Xrefs_ByIdAndBySymbol()
    {
        var c = Client();
        var hgnc = Assert.Single(await c.GetXrefsAsync("ENSG00000244734", "HGNC"));
        Assert.Equal("HGNC:4827", hgnc.PrimaryId);
        Assert.Equal("HBB", hgnc.DisplayId);
        Assert.Contains("BETA-GLOBIN", hgnc.Synonyms);
        Assert.Equal("HGNC", _handler.Requests.Single().Param("external_db"));
        var ids = await c.XrefsBySymbolAsync("homo_sapiens", "HBB");
        Assert.Contains(new EnsemblIdentity("ENSG00000244734", "gene"), ids);
    }

    [Fact]
    public async Task ErrorsCarryTheEnsemblMessage()
    {
        var handler = new FakeHttpHandler().On(_ => true, () => Fixtures.Json("{\"error\":\"No valid lookup found for symbol NOTAGENE\"}", System.Net.HttpStatusCode.BadRequest));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() =>
            new EnsemblClient(handler.Client(), TestOptions.Create(_clock)).LookupSymbolAsync("homo_sapiens", "NOTAGENE"));
        Assert.Equal("No valid lookup found for symbol NOTAGENE", ex.ServiceMessage);
        Assert.Equal("Ensembl", ex.Service);
    }

    [Fact]
    public async Task RequestsAreSpacedAtFifteenPerSecond()
    {
        var c = new EnsemblClient(_handler.Client(), TestOptions.Create(_clock, mode: DataCacheMode.Bypass));
        for (int i = 0; i < 3; i++) await c.LookupSymbolAsync("homo_sapiens", "BRCA2");
        Assert.InRange(_clock.Delays.Sum(d => d.TotalMilliseconds), 2000.0 / 15, 2000.0 / 15 + 1); // third request starts at 2/15 s
        Assert.Equal(2, _clock.Delays.Count);
    }
}
