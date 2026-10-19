using System.Net;
using System.Text.Json;
using Biotech.Net.Core;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.Pdb;

namespace Biotech.Net.Data.Tests;

public class PdbTests
{
    private const string Crambin = "TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN";
    private readonly FakeHttpHandler _handler = new();
    private readonly FakeClock _clock = new();

    public PdbTests()
    {
        var P = FakeHttpHandler.Path;
        _handler
            .OnFixture(P("/download/1CRN.pdb.gz"), "pdb/1CRN.pdb.gz", "application/octet-stream")
            .OnFixture(P("/download/1CRN.cif.gz"), "pdb/1CRN.cif.gz", "application/octet-stream")
            .OnFixture(P("/rest/v1/core/entry/1CRN"), "pdb/entry-1CRN.json")
            .OnFixture(P("/rest/v1/core/polymer_entity/1CRN/1"), "pdb/polymer-entity-1CRN-1.json")
            .OnFixture(P("/graphql"), "pdb/graphql-1CRN.json")
            .OnFixture(r => r.RequestUri!.Host == "search.rcsb.org" && Body(r).Contains("\"full_text\""), "pdb/search-text.json")
            .OnFixture(r => r.RequestUri!.Host == "search.rcsb.org" && Body(r).Contains("\"sequence\""), "pdb/search-sequence.json")
            .OnFixture(P("/api/prediction/P01542"), "pdb/alphafold-P01542.json")
            .OnFixture(P("/files/AF-P01542-F1-model_v6.pdb"), "pdb/AF-P01542-F1-model.pdb", "chemical/x-pdb");
    }

    private static string Body(HttpRequestMessage r) => r.Content?.ReadAsStringAsync().GetAwaiter().GetResult() ?? "";

    private PdbClient Client(IDataCache? cache = null) => new(_handler.Client(), TestOptions.Create(_clock, cache));

    [Fact]
    public async Task MmCif_IsDownloadedGzipped_AndParsed()
    {
        var s = await Client().GetStructureAsync("1crn");
        Assert.Equal("1CRN", s.Id);
        Assert.Equal(327, s.AtomCount);
        var chain = Assert.Single(s.Chains);
        Assert.Equal("A", chain.Id);
        Assert.Equal(46, chain.AminoAcids().Count);
        Assert.Equal("https://files.rcsb.org/download/1CRN.cif.gz", _handler.Requests.Single().Uri.AbsoluteUri);
    }

    [Fact]
    public async Task PdbFormat_GivesTheSameModel()
    {
        var c = Client();
        var pdb = await c.GetStructureAsync("1CRN", StructureFileFormat.Pdb);
        var cif = await c.GetStructureAsync("1CRN", StructureFileFormat.MmCif);
        Assert.Equal(cif.AtomCount, pdb.AtomCount);
        Assert.Equal(cif.Chains[0].ToSequence().ToString(), pdb.Chains[0].ToSequence().ToString());
        Assert.Equal("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", pdb.Chains[0].ToSequence().ToString()); // 1CRN is the PL/SI mixture (I25)
        Assert.Contains("CRAMBIN", pdb.Title, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public async Task CachedDownloads_StayCompressed()
    {
        var cache = new MemoryDataCache();
        await Client(cache).GetStructureTextAsync("1CRN", StructureFileFormat.Pdb);
        Assert.Equal(new FileInfo(Fixtures.File("pdb/1CRN.pdb.gz")).Length, cache.SizeBytes);
    }

    [Fact]
    public async Task MissingPdbFormat_ExplainsTheMmCifAlternative()
    {
        var handler = new FakeHttpHandler().On(_ => true, () => Fixtures.Status(HttpStatusCode.NotFound));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() =>
            new PdbClient(handler.Client(), TestOptions.Create(_clock)).GetStructureAsync("8GLV", StructureFileFormat.Pdb));
        Assert.True(ex.IsNotFound);
        Assert.Contains("StructureFileFormat.MmCif", ex.Message);
        Assert.Throws<ArgumentException>(() => PdbClient.NormalizeId("1C/N"));
    }

    [Fact]
    public async Task Entry_ParsesMethodResolutionDatesAndCitation()
    {
        var e = await Client().GetEntryAsync("1CRN");
        Assert.Equal("1CRN", e.Id);
        Assert.StartsWith("WATER STRUCTURE OF A HYDROPHOBIC PROTEIN", e.Title);
        Assert.Equal(["X-RAY DIFFRACTION"], e.ExperimentalMethods);
        Assert.Equal(1.5, e.Resolution);
        Assert.Equal(new DateOnly(1981, 4, 30), e.DepositionDate);
        Assert.Equal(new DateOnly(1981, 7, 28), e.ReleaseDate);
        Assert.Equal(327, e.DepositedAtomCount);
        Assert.Equal(1, e.PolymerEntityCount);
        Assert.Equal(4.74, e.MolecularWeightKDa);
        Assert.Equal("10.1073/pnas.81.19.6014", e.CitationDoi);
        Assert.Equal(16593516, e.CitationPubMedId);
        Assert.Equal(["Teeter, M.M."], e.CitationAuthors);
        Assert.Equal("PLANT PROTEIN", e.Keywords);
    }

    [Fact]
    public async Task PolymerEntity_ParsesSequenceChainsAndUniProtMapping()
    {
        var p = await Client().GetPolymerEntityAsync("1CRN", "1");
        Assert.Equal("1CRN", p.EntryId);
        Assert.Equal("1", p.EntityId);
        Assert.Equal("CRAMBIN", p.Description);
        Assert.Equal("Protein", p.PolymerType);
        Assert.Equal(["A"], p.ChainIds);
        Assert.Equal(["P01542"], p.UniProtAccessions);
        Assert.Equal([3721L], p.TaxonomyIds);
        var seq = p.ToSequence();
        Assert.Equal("1CRN_1", seq.Id);
        Assert.Equal(46, seq.Length);
        Assert.Equal(SequenceType.Protein, seq.Type);
    }

    [Fact]
    public async Task GraphQl_ReturnsTheDataMember()
    {
        var data = await Client().GraphQlAsync("{ entry(entry_id: \"1CRN\") { rcsb_id struct { title } exptl { method } } }");
        Assert.Equal("1CRN", data.GetProperty("entry").GetProperty("rcsb_id").GetString());
        var req = _handler.Requests.Single();
        Assert.Equal(HttpMethod.Post, req.Method);
        Assert.Equal("1CRN", JsonDocument.Parse(req.Body!).RootElement.GetProperty("query").GetString()!.Split('"')[1]);
    }

    [Fact]
    public async Task GraphQlErrors_Throw()
    {
        var handler = new FakeHttpHandler().On(_ => true, () => Fixtures.Json("{\"errors\":[{\"message\":\"Validation error: field 'nope' undefined\"}]}"));
        var ex = await Assert.ThrowsAsync<DataServiceException>(() =>
            new PdbClient(handler.Client(), TestOptions.Create(_clock)).GraphQlAsync("{ nope }"));
        Assert.Contains("field 'nope' undefined", ex.Message);
    }

    [Fact]
    public async Task FullTextSearch_BuildsTheDocumentedQuery()
    {
        var query = PdbSearchQuery.FullText("crambin").With(rows: 5);
        var r = await Client().SearchAsync(query);
        Assert.Equal(45, r.TotalCount);
        Assert.Equal("entry", r.ResultType);
        Assert.Equal(["3SZS", "2FD7", "1EJG", "2FD9", "1CCM"], r.Identifiers);
        Assert.Equal(1.0, r.Hits[0].Score);
        // The body is exactly what the recording script sent.
        Assert.Equal(
            "{\"query\":{\"type\":\"terminal\",\"service\":\"full_text\",\"parameters\":{\"value\":\"crambin\"}},\"return_type\":\"entry\",\"request_options\":{\"paginate\":{\"start\":0,\"rows\":5},\"results_content_type\":[\"experimental\"]}}",
            _handler.Requests.Single().Body);
    }

    [Fact]
    public async Task SequenceSearch_ReturnsPolymerEntities()
    {
        var query = PdbSearchQuery.Sequence(Crambin, identityCutoff: 0.9).With(rows: 5);
        Assert.Equal("polymer_entity", query.ReturnType);
        var r = await Client().SearchAsync(query);
        Assert.Equal("polymer_entity", r.ResultType);
        Assert.Equal(27, r.TotalCount);
        Assert.All(r.Identifiers, id => Assert.EndsWith("_1", id));
        var body = JsonDocument.Parse(_handler.Requests.Single().Body!).RootElement;
        var p = body.GetProperty("query").GetProperty("parameters");
        Assert.Equal(Crambin, p.GetProperty("value").GetString());
        Assert.Equal("protein", p.GetProperty("sequence_type").GetString());
        Assert.Equal(0.9, p.GetProperty("identity_cutoff").GetDouble());
    }

    [Fact]
    public void QueryBuilders_ComposeGroupsAndAttributes()
    {
        var q = PdbSearchQuery.And(PdbSearchQuery.FullText("kinase"), PdbSearchQuery.Attribute("rcsb_entry_info.resolution_combined", "less", 1.5))
            .With(rows: 3, includeComputedModels: true);
        var root = JsonDocument.Parse(q.ToJson()).RootElement;
        Assert.Equal("group", root.GetProperty("query").GetProperty("type").GetString());
        Assert.Equal("and", root.GetProperty("query").GetProperty("logical_operator").GetString());
        var attr = root.GetProperty("query").GetProperty("nodes")[1].GetProperty("parameters");
        Assert.Equal(1.5, attr.GetProperty("value").GetDouble());
        Assert.Equal(2, root.GetProperty("request_options").GetProperty("results_content_type").GetArrayLength());
    }

    [Fact]
    public async Task EmptySearch_NoContentMeansNoHits()
    {
        var handler = new FakeHttpHandler().On(_ => true, () => new HttpResponseMessage(HttpStatusCode.NoContent));
        var r = await new PdbClient(handler.Client(), TestOptions.Create(_clock)).SearchAsync(PdbSearchQuery.FullText("zzzzqqq"));
        Assert.Equal(0, r.TotalCount);
        Assert.Empty(r.Hits);
    }

    [Fact]
    public async Task AlphaFold_PredictionAndModel()
    {
        var af = new AlphaFoldClient(_handler.Client(), TestOptions.Create(_clock));
        var p = Assert.Single(await af.GetPredictionsAsync("P01542"));
        Assert.Equal("AF-P01542-F1", p.EntryId);
        Assert.Equal("P01542", p.UniProtAccession);
        Assert.Equal(Crambin, p.SequenceResidues);
        Assert.Equal(93.5, p.MeanPlddt);
        Assert.Equal(6, p.LatestVersion);
        Assert.Equal("https://alphafold.ebi.ac.uk/files/AF-P01542-F1-model_v6.pdb", p.PdbUrl!.AbsoluteUri);

        var s = await af.GetStructureAsync("P01542");
        Assert.Equal("AF-P01542-F1", s.Id);
        Assert.Equal(Crambin, s.Chains[0].ToSequence().ToString());
        // pLDDT is stored as the B-factor: every CA is confidently predicted (> 50).
        var ca = s.Atoms.Where(a => a.Name == "CA").ToList();
        Assert.Equal(46, ca.Count);
        Assert.InRange(ca.Average(a => a.BFactor), 85, 100);
    }
}
