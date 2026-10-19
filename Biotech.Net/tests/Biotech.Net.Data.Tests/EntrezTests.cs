using Biotech.Net.Core;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.Ncbi;

namespace Biotech.Net.Data.Tests;

public class EntrezTests
{
    private readonly FakeHttpHandler _handler = new();
    private readonly FakeClock _clock = new();
    private readonly MemoryDataCache _cache = new();

    public EntrezTests()
    {
        var P = FakeHttpHandler.Path;
        _handler
            .OnFixture(P("einfo.fcgi", q => q.ContainsKey("db")), "ncbi/einfo-pubmed.json")
            .OnFixture(P("einfo.fcgi"), "ncbi/einfo.json")
            .OnFixture(P("esearch.fcgi", q => Q(q, "db") == "nuccore"), "ncbi/esearch-nuccore.json")
            .OnFixture(P("esearch.fcgi", q => Q(q, "db") == "pubmed"), "ncbi/esearch-pubmed.json")
            .OnFixture(P("esearch.fcgi", q => Q(q, "db") == "protein" && q.ContainsKey("usehistory")), "ncbi/esearch-history.json")
            .OnFixture(P("esummary.fcgi", q => Q(q, "db") == "pubmed"), "ncbi/esummary-pubmed.json")
            .OnFixture(P("efetch.fcgi", q => q.ContainsKey("WebEnv") && Q(q, "retstart") == "0"), "ncbi/efetch-history-0.fasta", "text/plain")
            .OnFixture(P("efetch.fcgi", q => q.ContainsKey("WebEnv") && Q(q, "retstart") == "2"), "ncbi/efetch-history-2.fasta", "text/plain")
            .OnFixture(P("efetch.fcgi", q => q.TryGetValue("rettype", out var t) && t[0] == "fasta"), "ncbi/efetch-protein.fasta", "text/plain")
            .OnFixture(P("efetch.fcgi", q => q.TryGetValue("rettype", out var t) && t[0] == "gb"), "ncbi/efetch-nuccore.gb", "text/plain")
            .OnFixture(P("elink.fcgi"), "ncbi/elink-protein-gene.json")
            .OnFixture(r => r.Method == HttpMethod.Post && r.RequestUri!.AbsolutePath.EndsWith("esummary.fcgi"), "ncbi/esummary-pubmed.json");
    }

    private static string? Q(Dictionary<string, List<string>> q, string name) => q.TryGetValue(name, out var v) ? v[0] : null;

    private EntrezClient Client(string? apiKey = null, DataClientOptions? options = null) =>
        new(_handler.Client(), options ?? TestOptions.Create(_clock, _cache)) { ApiKey = apiKey, Email = "tests@example.org" };

    [Fact]
    public async Task EInfo_ListsDatabases_AndSendsToolAndEmail()
    {
        var dbs = await Client().EInfoDatabasesAsync();
        Assert.Contains("pubmed", dbs);
        Assert.Contains("nuccore", dbs);
        var req = _handler.Requests.Single();
        Assert.Equal("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi", req.Uri.GetLeftPart(UriPartial.Path));
        Assert.Equal("Biotech.Net", req.Param("tool"));
        Assert.Equal("tests@example.org", req.Param("email"));
        Assert.Null(req.Param("api_key"));
        Assert.Equal("json", req.Param("retmode"));
    }

    [Fact]
    public async Task EInfo_Database_ParsesFieldsAndLinks()
    {
        var info = await Client().EInfoAsync("pubmed");
        Assert.Equal("pubmed", info.Name);
        Assert.Equal("PubMed", info.MenuName);
        Assert.True(info.Count > 30_000_000);
        Assert.Contains(info.Fields, f => f.Name == "ALL" && f.FullName == "All Fields");
        Assert.Contains(info.Links, l => l.DbTo == "assembly");
    }

    [Fact]
    public async Task ESearch_ParsesIdsCountAndTranslation()
    {
        var r = await Client().ESearchAsync("nuccore", "HBB[gene] AND Homo sapiens[orgn] AND refseq[filter]", retMax: 5);
        Assert.Equal(9, r.Count);
        Assert.Equal(5, r.Ids.Count);
        Assert.Contains("28380636", r.Ids);
        Assert.Contains("\"Homo sapiens\"[Organism]", r.QueryTranslation);
        Assert.Null(r.History);
        var req = _handler.Requests.Single();
        Assert.Equal("HBB[gene] AND Homo sapiens[orgn] AND refseq[filter]", req.Param("term"));
        Assert.Equal("5", req.Param("retmax"));
        Assert.Contains("term=HBB%5Bgene%5D%20AND", req.Uri.AbsoluteUri);
    }

    [Fact]
    public async Task PubMedSummaries_ParseTitleJournalAuthorsAndDoi()
    {
        var c = Client();
        var search = await c.ESearchAsync("pubmed", "crambin structure", retMax: 3, sort: "relevance");
        Assert.Equal(["39190507", "6895315", "19562114"], search.Ids);
        var docs = await c.PubMedSummariesAsync(["6895315", "19562114"]);
        Assert.Equal(2, docs.Count);
        var teeter = docs[0];
        Assert.Equal("6895315", teeter.Pmid);
        Assert.Equal("Primary structure of the hydrophobic plant protein crambin.", teeter.Title);
        Assert.Equal("Biochemistry", teeter.Journal);
        Assert.Equal("1981 Sep 15", teeter.PubDate);
        Assert.Equal(["Teeter MM", "Mazer JA", "L'Italien JJ"], teeter.Authors);
        Assert.Equal("10.1021/bi00522a013", teeter.Doi);
        Assert.Equal("6895315,19562114", _handler.Requests.Last().Param("id"));
    }

    [Fact]
    public async Task FetchSequence_Fasta_ParsesAProtein()
    {
        var hbb = await Client().FetchSequenceAsync("protein", "NP_000509.1");
        Assert.Equal("NP_000509.1", hbb.Id);
        Assert.Equal(SequenceType.Protein, hbb.Type);
        Assert.Equal(147, hbb.Length);
        Assert.StartsWith("MVHLTPEEKSAVTALWGKV", hbb.ToString());
        Assert.Contains("hemoglobin subunit beta", hbb.Description);
        var req = _handler.Requests.Single();
        Assert.Equal("fasta", req.Param("rettype"));
        Assert.Equal("text", req.Param("retmode"));
    }

    [Fact]
    public async Task FetchSequence_GenBank_KeepsAnnotationsAndFeatures()
    {
        var mrna = await Client().FetchSequenceAsync("nuccore", "NM_000518.5", EntrezSequenceFormat.GenBank);
        Assert.Equal("NM_000518.5", mrna.Id);
        Assert.Equal(SequenceType.Dna, mrna.Type);
        Assert.Equal(628, mrna.Length);
        Assert.Equal("Homo sapiens", mrna.Annotations["organism"]);
        var cds = mrna.Features.Single(f => f.Type == "CDS");
        // The CDS translates to the HBB protein (147 aa + stop) with our own kernels.
        var protein = cds.Location.Extract(mrna).Translate(toStop: true);
        Assert.Equal(147, protein.Length);
        Assert.StartsWith("MVHLTPEEK", protein.ToString());
    }

    [Fact]
    public async Task ELink_FindsTheGeneOfAProtein()
    {
        var sets = await Client().ELinkAsync("protein", "gene", ["NP_000509.1"]);
        var link = sets.Single().Links.Single();
        Assert.Equal("protein_gene", link.LinkName);
        Assert.Equal(["3043"], link.Ids);
        Assert.Equal(["4504349"], sets[0].Ids);
    }

    [Fact]
    public async Task HistoryServer_StreamsAllRecordsInBatches()
    {
        var c = Client();
        var search = await c.ESearchAsync("protein", "(HBA1[gene] OR HBB[gene] OR HBD[gene]) AND Homo sapiens[orgn] AND refseq[filter]",
            retMax: 0, useHistory: true);
        Assert.Equal(3, search.Count);
        Assert.NotNull(search.History);
        Assert.Equal("1", search.History.QueryKey);
        Assert.StartsWith("MCID_", search.History.WebEnv);

        var all = new List<Sequence>();
        await foreach (var s in c.FetchSequencesAsync(search, batchSize: 2)) all.Add(s);
        Assert.Equal(["NP_000509.1", "NP_000510.1", "NP_000549.1"], all.Select(s => s.Id).Order().ToArray());
        var fetches = _handler.Requests.Where(r => r.Path.EndsWith("efetch.fcgi")).ToList();
        Assert.Equal(["0", "2"], fetches.Select(f => f.Param("retstart")));
        Assert.All(fetches, f =>
        {
            Assert.Equal(search.History.WebEnv, f.Param("WebEnv"));
            Assert.Equal("1", f.Param("query_key"));
            Assert.Equal("2", f.Param("retmax"));
        });
    }

    [Fact]
    public async Task ApiKey_RaisesTheRate_IsSent_ButNeverEntersTheCacheKey()
    {
        var withKey = Client(apiKey: "KEY123");
        Assert.Equal(10, withKey.RequestsPerSecond);
        Assert.Equal(3, Client().RequestsPerSecond);
        await withKey.EInfoDatabasesAsync();
        Assert.Equal("KEY123", _handler.Requests.Single().Param("api_key"));
        await Client().EInfoDatabasesAsync(); // same request without the key: served from the cache
        Assert.Equal(1, _handler.Count);
    }

    [Fact]
    public async Task RequestsAreSpacedAtThreePerSecondWithoutAKey()
    {
        var options = TestOptions.Create(_clock, _cache, DataCacheMode.Bypass);
        var c = Client(options: options);
        for (int i = 0; i < 4; i++) await c.EInfoDatabasesAsync();
        Assert.Equal(3, _clock.Delays.Count);
        double elapsed = 0;
        int k = 0;
        foreach (var d in _clock.Delays)
        {
            elapsed += d.TotalMilliseconds;
            k++;
            Assert.InRange(elapsed, k * 1000.0 / 3, k * 1000.0 / 3 + 1); // request k starts at k/3 s
        }

    }

    [Fact]
    public async Task LongIdLists_ArePosted()
    {
        var ids = Enumerable.Range(1, 250).Select(i => i.ToString()).ToList();
        await Client().PubMedSummariesAsync(ids);
        var req = _handler.Requests.Single();
        Assert.Equal(HttpMethod.Post, req.Method);
        Assert.Equal("application/x-www-form-urlencoded", req.ContentType);
        Assert.Contains("id=1%2C2%2C3", req.Body);
        Assert.Null(req.Param("id"));
        Assert.Equal("Biotech.Net", req.Param("tool"));
    }

    [Fact]
    public async Task ServiceErrors_BecomeDataServiceExceptions()
    {
        var handler = new FakeHttpHandler()
            .On(FakeHttpHandler.Path("esearch.fcgi"), () => Fixtures.Json("{\"header\":{},\"esearchresult\":{\"ERROR\":\"Invalid db name specified: nucleotid\"}}"))
            .On(FakeHttpHandler.Path("efetch.fcgi"), () => Fixtures.Status(System.Net.HttpStatusCode.BadRequest,
                "Cannot+retrieve+history+data%3A+OUT+OF+RANGE."));
        var c = new EntrezClient(handler.Client(), TestOptions.Create(_clock)) { ApiKey = "SECRET" };
        var ex = await Assert.ThrowsAsync<DataServiceException>(() => c.ESearchAsync("nucleotid", "x"));
        Assert.Contains("Invalid db name", ex.Message);
        var ex2 = await Assert.ThrowsAsync<DataServiceException>(() => c.EFetchAsync("protein", new EntrezHistory("W", "1"), "fasta", retStart: 9));
        Assert.Equal("Cannot retrieve history data: OUT OF RANGE.", ex2.ServiceMessage);
        Assert.DoesNotContain("SECRET", ex2.Message);
    }
}
