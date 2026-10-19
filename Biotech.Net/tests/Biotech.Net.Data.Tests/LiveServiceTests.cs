using Biotech.Net.Core;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Ensembl;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.Ncbi;
using Biotech.Net.Data.Pdb;
using Biotech.Net.Data.UniProt;
using Biotech.Net.Genomics.Annotation;

namespace Biotech.Net.Data.Tests;

/// <summary>
/// Calls the real services. Skipped unless <c>BIOTECH_NET_LIVE_TESTS=1</c>; NCBI calls identify themselves with
/// <c>NCBI_EMAIL</c> (and use <c>NCBI_API_KEY</c>) when those are set. Each test uses a fresh on-disk cache, so it
/// always reaches the network, and then checks the cached copy replays offline.
/// </summary>
public sealed class LiveServiceTests : IDisposable
{
    private readonly TempDir _dir = new();

    private DataClientOptions Options(DataCacheMode mode = DataCacheMode.Default) => new()
    {
        Cache = new DiskDataCache(_dir.Path),
        CacheMode = mode,
        RateLimiter = HostRateLimiter.Shared,
    };

    private static void RequireLive()
    {
        if (Environment.GetEnvironmentVariable("BIOTECH_NET_LIVE_TESTS") != "1")
            Assert.Skip("Live service tests run only with BIOTECH_NET_LIVE_TESTS=1.");
    }

    public void Dispose() => _dir.Dispose();

    [Fact]
    public async Task Ncbi_SearchFetchSummariesLinksAndHistory()
    {
        RequireLive();
        var ncbi = new EntrezClient(options: Options());
        var ct = TestContext.Current.CancellationToken;

        Assert.Contains("pubmed", await ncbi.EInfoDatabasesAsync(ct));
        var search = await ncbi.ESearchAsync("nuccore", "NM_000518[accn]", cancellationToken: ct);
        Assert.True(search.Count >= 1);

        var mrna = await ncbi.FetchSequenceAsync("nuccore", "NM_000518.5", EntrezSequenceFormat.GenBank, ct);
        Assert.Equal("NM_000518.5", mrna.Id);
        Assert.Contains(mrna.Features, f => f.Type == "CDS");

        var hbb = await ncbi.FetchSequenceAsync("protein", "NP_000509.1", cancellationToken: ct);
        Assert.Equal(147, hbb.Length);

        var paper = Assert.Single(await ncbi.PubMedSummariesAsync(["6895315"], ct));
        Assert.Contains("crambin", paper.Title, StringComparison.OrdinalIgnoreCase);

        var links = await ncbi.ELinkAsync("protein", "gene", ["NP_000509.1"], cancellationToken: ct);
        Assert.Contains("3043", links.SelectMany(l => l.Links).SelectMany(l => l.Ids));

        var history = await ncbi.ESearchAsync("protein", "(HBA1[gene] OR HBB[gene] OR HBD[gene]) AND Homo sapiens[orgn] AND refseq[filter]",
            retMax: 0, useHistory: true, cancellationToken: ct);
        var all = new List<Sequence>();
        await foreach (var s in ncbi.FetchSequencesAsync(history, batchSize: 2, cancellationToken: ct)) all.Add(s);
        Assert.Equal(history.Count, all.Count);

        // The same GenBank record replays from the disk cache with the network switched off.
        var offline = new EntrezClient(options: Options(DataCacheMode.Offline));
        Assert.Equal(mrna.ToString(), (await offline.FetchSequenceAsync("nuccore", "NM_000518.5", EntrezSequenceFormat.GenBank, ct)).ToString());
    }

    [Fact]
    public async Task UniProt_EntrySearchPagesAndIdMapping()
    {
        RequireLive();
        var uniprot = new UniProtClient(options: Options());
        var ct = TestContext.Current.CancellationToken;

        var hba = await uniprot.GetEntryAsync("P69905", ct);
        Assert.Equal("HBA_HUMAN", hba.EntryName);
        Assert.Equal(142, hba.Length);
        Assert.Contains(hba.CrossReferencesTo("PDB"), x => x.Id == "4HHB");
        Assert.Equal(142, (await uniprot.GetFastaAsync("P69905", ct)).Length);

        var entries = new List<UniProtEntry>();
        await foreach (var e in uniprot.SearchAsync("gene:HBB AND reviewed:true", "accession,length", pageSize: 2, maxResults: 5, ct)) entries.Add(e);
        Assert.Equal(5, entries.Count);
        Assert.Equal(5, entries.Select(e => e.Accession).Distinct().Count());

        var mapping = await uniprot.MapIdsAsync("UniProtKB_AC-ID", "GeneID", ["P69905", "P68871"], cancellationToken: ct);
        Assert.Contains("3043", mapping.BySource()["P68871"]);
    }

    [Fact]
    public async Task Pdb_StructureEntrySearchAndAlphaFold()
    {
        RequireLive();
        var pdb = new PdbClient(options: Options());
        var ct = TestContext.Current.CancellationToken;

        var hemoglobin = await pdb.GetStructureAsync("4HHB", cancellationToken: ct);
        Assert.Equal(4, hemoglobin.Chains.Count(c => c.AminoAcids().Count > 0));
        var pdbFormat = await pdb.GetStructureAsync("1CRN", StructureFileFormat.Pdb, ct);
        Assert.Equal(327, pdbFormat.AtomCount);

        var entry = await pdb.GetEntryAsync("4HHB", ct);
        Assert.Equal(1.74, entry.Resolution);
        var entity = await pdb.GetPolymerEntityAsync("4HHB", "1", ct);
        Assert.Contains("P69905", entity.UniProtAccessions);

        var text = await pdb.SearchAsync(PdbSearchQuery.FullText("crambin").With(rows: 5), ct);
        Assert.True(text.TotalCount > 10);
        // UniProt crambin (L25) vs 1CRN (the P/S mixture, I25): 1CRN ranks below the exact matches, so fetch every hit.
        var seq = await pdb.SearchAsync(PdbSearchQuery.Sequence("TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN", identityCutoff: 0.9).With(rows: 100), ct);
        Assert.Contains("1CRN_1", seq.Identifiers);
        var data = await pdb.GraphQlAsync("{ entry(entry_id: \"4HHB\") { rcsb_id } }", cancellationToken: ct);
        Assert.Equal("4HHB", data.GetProperty("entry").GetProperty("rcsb_id").GetString());

        var af = new AlphaFoldClient(options: Options());
        var model = await af.GetStructureAsync("P69905", cancellationToken: ct);
        Assert.Equal(142, model.Chains[0].AminoAcids().Count);
    }

    [Fact]
    public async Task Ensembl_LookupSequenceOverlapVepAndXrefs()
    {
        RequireLive();
        var ensembl = new EnsemblClient(options: Options());
        var ct = TestContext.Current.CancellationToken;

        var brca2 = await ensembl.LookupSymbolAsync("homo_sapiens", "BRCA2", cancellationToken: ct);
        Assert.Equal("ENSG00000139618", brca2.Id);
        var hbb = (await ensembl.LookupIdAsync("ENSG00000244734", expand: true, ct)).ToGene();
        Assert.Contains(hbb.Transcripts, t => t.Id == "ENST00000335295");

        var protein = await ensembl.GetSequenceAsync("ENST00000335295", EnsemblSequenceType.Protein, ct);
        var cds = await ensembl.GetSequenceAsync("ENST00000335295", EnsemblSequenceType.Cds, ct);
        Assert.Equal(protein.ToString(), cds.Translate(toStop: true).ToString());

        var overlap = await ensembl.OverlapRegionAsync("human", GenomicInterval.Parse("11:5225000-5230000"), cancellationToken: ct);
        Assert.Contains(overlap, f => f.ExternalName == "HBB");

        var vep = Assert.Single(await ensembl.VepHgvsAsync("human", "ENST00000335295:c.20A>T", ct));
        Assert.Equal(Consequence.MissenseVariant, vep.MostSevere);
        var region = Assert.Single(await ensembl.VepRegionAsync("human", new Variant("11", 5227002, "T", "A"), ct));
        Assert.Equal("missense_variant", region.MostSevereConsequence);

        Assert.Contains(await ensembl.GetXrefsAsync("ENSG00000244734", "HGNC", ct), x => x.PrimaryId == "HGNC:4827");
    }
}
