using Biotech.Net.Core;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.UniProt;

namespace Biotech.Net.Data.Tests;

public class UniProtTests
{
    private readonly FakeHttpHandler _handler = new();
    private readonly FakeClock _clock = new();

    public UniProtTests()
    {
        var P = FakeHttpHandler.Path;
        _handler
            .OnFixture(P("/uniprotkb/P01542", q => q["format"][0] == "json"), "uniprot/P01542.json")
            .OnFixture(P("/uniprotkb/P01542", q => q["format"][0] == "fasta"), "uniprot/P01542.fasta", "text/plain")
            .OnFixture(P("/uniprotkb/search", q => q.ContainsKey("cursor")), "uniprot/search-page2.json")
            .OnFixture(P("/uniprotkb/search"), "uniprot/search-page1.json");
    }

    private UniProtClient Client(IDataCache? cache = null, DataCacheMode mode = DataCacheMode.Default) =>
        new(_handler.Client(), TestOptions.Create(_clock, cache, mode));

    [Fact]
    public async Task Entry_ParsesTheJsonModel()
    {
        var e = await Client().GetEntryAsync("P01542");
        Assert.Equal("P01542", e.Accession);
        Assert.Equal("CRAM_CRAAB", e.EntryName);
        Assert.True(e.IsReviewed);
        Assert.Equal("Crambin", e.ProteinName);
        Assert.Equal(["THI2"], e.GeneNames);
        Assert.Equal("Crambe hispanica subsp. abyssinica", e.OrganismName);
        Assert.Equal(3721, e.TaxonId);
        Assert.Equal(46, e.Length);
        Assert.Equal("TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN", e.SequenceResidues);
        Assert.Equal(4736, e.MolecularWeight);
        Assert.Equal(2, e.SequenceVersion);
        Assert.StartsWith("The function of this hydrophobic plant seed protein is not known", e.Function);
        Assert.Contains(e.Comments, c => c.Type == "SUBCELLULAR LOCATION" && c.Text.Contains("Secreted"));
        Assert.Equal(3, e.Features.Count(f => f.Type == "Disulfide bond"));
        Assert.Contains(e.Features, f => f is { Type: "Disulfide bond", Start: 3, End: 40 });
        Assert.Contains(e.Keywords, k => k.Name == "3D-structure");
        Assert.Contains(e.CrossReferencesTo("PDB"), x => x.Id == "1CRN");
        Assert.Equal("X-ray", e.CrossReferencesTo("PDB").First(x => x.Id == "1AB1").Properties["Method"]);
        var req = _handler.Requests.Single();
        Assert.Equal("https://rest.uniprot.org/uniprotkb/P01542?format=json", req.Uri.AbsoluteUri);
    }

    [Fact]
    public async Task Entry_ConvertsToCoreProteinAndSequence()
    {
        var p = await Client().GetProteinAsync("P01542");
        Assert.Equal("P01542", p.Id);
        Assert.Equal("Crambin", p.Name);
        Assert.Equal("THI2", p.GeneName);
        Assert.Equal("Crambe hispanica subsp. abyssinica", p.Organism);
        Assert.NotNull(p.StructureId);
        Assert.Equal(46, p.Length);
        // Our own ProtParam mass agrees with UniProt's to within a dalton.
        Assert.Equal(4736, p.MolecularWeight, 0.0 + 1.0);
        Assert.Equal(SequenceType.Protein, p.Sequence.Type);
    }

    [Fact]
    public async Task Fasta_ParsesTheCanonicalSequence()
    {
        var s = await Client().GetFastaAsync("P01542");
        Assert.Equal("sp|P01542|CRAM_CRAAB", s.Id);
        Assert.StartsWith("Crambin OS=Crambe", s.Description);
        Assert.Equal(46, s.Length);
        Assert.Equal(SequenceType.Protein, s.Type);
    }

    [Fact]
    public async Task SearchPage_ReadsTotalAndNextLink()
    {
        var page = await Client().SearchPageAsync("gene:HBB AND reviewed:true", "accession,id,protein_name,gene_names,organism_name,length,sequence", 2);
        Assert.Equal(2, page.Entries.Count);
        Assert.Equal(243, page.TotalResults);
        Assert.Equal("P68871", page.Entries[0].Accession);
        Assert.Equal("Hemoglobin subunit beta", page.Entries[0].ProteinName);
        Assert.Equal(147, page.Entries[0].Length);
        Assert.NotNull(page.NextPage);
        Assert.Contains("cursor=", page.NextPage.Query);
        var req = _handler.Requests.Single();
        Assert.Equal("gene:HBB AND reviewed:true", req.Param("query"));
        Assert.Equal("2", req.Param("size"));
        Assert.Equal("json", req.Param("format"));
    }

    [Fact]
    public async Task Search_FollowsLinkHeadersAcrossPages()
    {
        var all = new List<UniProtEntry>();
        await foreach (var e in Client().SearchAsync("gene:HBB AND reviewed:true", pageSize: 2, maxResults: 4)) all.Add(e);
        Assert.Equal(4, all.Count);
        Assert.Equal(4, all.Select(e => e.Accession).Distinct().Count());
        Assert.Equal(2, _handler.Count);
        Assert.Contains("cursor=", _handler.Requests.Last().Uri.Query);
    }

    [Fact]
    public async Task Search_StopsAtMaxResults_WithoutFetchingMorePages()
    {
        var all = new List<UniProtEntry>();
        await foreach (var e in Client().SearchAsync("gene:HBB AND reviewed:true", pageSize: 2, maxResults: 2)) all.Add(e);
        Assert.Equal(2, all.Count);
        Assert.Equal(1, _handler.Count);
    }

    [Fact]
    public async Task CachedSearchPages_KeepTheirLinkHeader()
    {
        var cache = new MemoryDataCache();
        await Client(cache).SearchPageAsync("gene:HBB AND reviewed:true", pageSize: 2);
        var offline = await Client(cache, DataCacheMode.Offline).SearchPageAsync("gene:HBB AND reviewed:true", pageSize: 2);
        Assert.NotNull(offline.NextPage);
        Assert.Equal(243, offline.TotalResults);
        Assert.Equal(1, _handler.Count);
    }

    [Fact]
    public void NextLink_ParsesRelNext()
    {
        Assert.Equal(new Uri("https://rest.uniprot.org/x?cursor=a"),
            UniProtClient.NextLink("<https://rest.uniprot.org/x?cursor=a>; rel=\"next\""));
        Assert.Null(UniProtClient.NextLink("<https://rest.uniprot.org/x>; rel=\"prev\""));
        Assert.Null(UniProtClient.NextLink(null));
    }

    [Fact]
    public async Task IdMapping_SubmitsPollsAndCollectsResults()
    {
        const string job = "1axZrMTLxk";
        var handler = new FakeHttpHandler()
            .OnFixture(FakeHttpHandler.Path("/idmapping/run"), "uniprot/idmapping-run.json")
            .OnSequence(FakeHttpHandler.Path($"/idmapping/status/{job}"),
                () => Fixtures.Json("{\"jobStatus\":\"NEW\"}"),
                () => Fixtures.Json("{\"jobStatus\":\"RUNNING\"}"),
                () => Fixtures.Response("uniprot/idmapping-status-first.json"))
            .OnFixture(FakeHttpHandler.Path($"/idmapping/results/{job}"), "uniprot/idmapping-results.json");
        var cache = new MemoryDataCache();
        var client = new UniProtClient(handler.Client(), TestOptions.Create(_clock, cache)) { IdMappingPollInterval = TimeSpan.FromSeconds(1) };
        var result = await client.MapIdsAsync("UniProtKB_AC-ID", "GeneID", ["P69905", "P68871", "P01542"]);

        Assert.Equal(job, result.JobId);
        Assert.Equal(["3039", "3040"], result.BySource()["P69905"]);
        Assert.Equal(["3043"], result.BySource()["P68871"]);
        Assert.Equal(["P01542"], result.FailedIds);

        var run = handler.Requests.First();
        Assert.Equal(HttpMethod.Post, run.Method);
        Assert.Equal("from=UniProtKB_AC-ID&to=GeneID&ids=P69905%2CP68871%2CP01542", run.Body);
        Assert.Equal(3, handler.Requests.Count(r => r.Path.Contains("/status/")));
        // Poll intervals grow by half (the 100 ms gaps are the 10 requests/s spacing).
        Assert.Equal([TimeSpan.FromSeconds(1), TimeSpan.FromSeconds(1.5)], _clock.Delays.Where(d => d >= TimeSpan.FromSeconds(1)).ToArray());
        Assert.Equal(0, cache.Count); // jobs are never cached
    }

    [Fact]
    public async Task IdMapping_ResultsForUniProtKbTargets_UseThePrimaryAccession()
    {
        var handler = new FakeHttpHandler()
            .On(FakeHttpHandler.Path("/idmapping/results/J"), () => Fixtures.Json(
                "{\"results\":[{\"from\":\"HBB\",\"to\":{\"primaryAccession\":\"P68871\",\"uniProtkbId\":\"HBB_HUMAN\"}}]}"));
        var r = await new UniProtClient(handler.Client(), TestOptions.Create(_clock)).GetIdMappingResultsAsync("J");
        Assert.Equal(new IdMapping("HBB", "P68871"), r.Mappings.Single());
    }

    [Fact]
    public async Task IdMapping_ErrorStateThrows()
    {
        var handler = new FakeHttpHandler()
            .On(FakeHttpHandler.Path("/idmapping/run"), () => Fixtures.Json("{\"jobId\":\"J\"}"))
            .On(FakeHttpHandler.Path("/idmapping/status/J"), () => Fixtures.Json("{\"jobStatus\":\"ERROR\"}"));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() =>
            new UniProtClient(handler.Client(), TestOptions.Create(_clock)).MapIdsAsync("A", "B", ["x"]));
        Assert.Contains("job J ended in an error", ex.Message);
    }

    [Fact]
    public async Task UnknownAccession_IsNotFound()
    {
        var handler = new FakeHttpHandler().On(_ => true, () => Fixtures.Json(
            "{\"url\":\"http://rest.uniprot.org/uniprotkb/XXXXXX\",\"messages\":[\"The 'accession' value has invalid format. It should be a valid UniProtKB accession\"]}",
            System.Net.HttpStatusCode.BadRequest));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => new UniProtClient(handler.Client(), TestOptions.Create(_clock)).GetEntryAsync("XXXXXX"));
        Assert.Equal("UniProt", ex.Service);
        Assert.Contains("invalid format", ex.ServiceMessage);
    }
}
