# Tutorial: Cloud data connectors

Package: `Gravicode.Biotech.Net.Data`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Fetch records straight from the public databases into Biotech.Net types: GenBank and FASTA records from NCBI, proteins from UniProt, structures from the RCSB PDB and AlphaFold DB, genes and variant effects from Ensembl. Every connector shares one engine that spaces requests to each host, retries transient failures and caches responses on disk, so a second run of a script works without the network.

Each snippet below was run against the live services; the outputs shown are from that run (database counts drift with new releases).

## 1. Set-up

```text
dotnet add package Gravicode.Biotech.Net.Data
```

`BioDataClient` bundles the five clients (`Entrez`, `UniProt`, `Pdb`, `AlphaFold`, `Ensembl`); each can also be created on its own. NCBI asks every tool to identify itself: set the `NCBI_EMAIL` environment variable to your contact address, and `NCBI_API_KEY` if you have a key (it raises the limit from 3 to 10 requests per second). Neither value is ever written to the cache or into an exception message.

## 2. A GenBank record from NCBI

```csharp
using Biotech.Net.Data;
using Biotech.Net.Data.Ncbi;

var data = new BioDataClient();   // NCBI_EMAIL / NCBI_API_KEY are picked up from the environment

var hbb = await data.Entrez.FetchSequenceAsync("nuccore", "NM_000518.5", EntrezSequenceFormat.GenBank);
Console.WriteLine($"{hbb.Id}: {hbb.Length} bp, {hbb.Annotations["organism"]}");
Console.WriteLine(hbb.Description);

var cds = hbb.Features.Single(f => f.Type == "CDS");
var protein = cds.Location.Extract(hbb).Translate(toStop: true);
Console.WriteLine($"CDS {cds.Location} -> {protein.Length} aa: {protein.Preview(20)}");
```

```text
NM_000518.5: 628 bp, Homo sapiens
Homo sapiens hemoglobin subunit beta (HBB), mRNA
CDS 51..494 -> 147 aa: MVHLTPEEKSAVTALWGKVN…VANALAHKYH (147 aa)
```

The record is parsed by the same `GenBankReader` as a local file, so annotations and the feature table come along. `EntrezSequenceFormat.Fasta` (the default) returns plain sequences.

### Searching, summaries and links

```csharp
using Biotech.Net.Data.Ncbi;

var ncbi = new EntrezClient();
var hits = await ncbi.ESearchAsync("pubmed", "crambin structure", retMax: 3, sort: "relevance");
Console.WriteLine($"{hits.Count} PubMed records; first {hits.Ids.Count}: {string.Join(", ", hits.Ids)}");
foreach (var paper in await ncbi.PubMedSummariesAsync(hits.Ids))
    Console.WriteLine($"  {paper.Pmid} {paper.PubDate[..4]} {paper.Journal}: {paper.Title}");

var links = await ncbi.ELinkAsync("protein", "gene", ["NP_000509.1"]);
Console.WriteLine($"NP_000509.1 -> gene {links[0].Links[0].Ids[0]}");
```

```text
176 PubMed records; first 3: 39190507, 6895315, 19562114
  39190507 2024 IUCrJ: Solvent organization in the ultrahigh-resolution crystal structure of crambin at room temperature.
  6895315 1981 Biochemistry: Primary structure of the hydrophobic plant protein crambin.
  19562114 2009 Mol Biosyst: Role of a salt bridge in the model protein crambin explored by chemical protein synthesis: X-ray structure of a unique protein analogue, [V15A]crambin-alpha-carboxamide.
NP_000509.1 -> gene 3043
```

`ESummaryAsync` returns every summary field of any database as JSON (`ESummaryDocument.Fields`); `PubMedSummariesAsync` types the PubMed ones. `EInfoAsync("pubmed")` lists a database's searchable fields and links.

### Large result sets: the history server

```csharp
using Biotech.Net.Data.Ncbi;

var ncbi = new EntrezClient();
var search = await ncbi.ESearchAsync("protein",
    "(HBA1[gene] OR HBB[gene] OR HBD[gene]) AND Homo sapiens[orgn] AND refseq[filter]",
    retMax: 0, useHistory: true);
Console.WriteLine($"{search.Count} records kept on the history server (query_key {search.History!.QueryKey})");

await foreach (var seq in ncbi.FetchSequencesAsync(search, batchSize: 2))
    Console.WriteLine($"  {seq.Id} {seq.Length} aa  {seq.Description}");
```

```text
3 records kept on the history server (query_key 1)
  NP_000509.1 147 aa  hemoglobin subunit beta [Homo sapiens]
  NP_000549.1 142 aa  hemoglobin subunit alpha [Homo sapiens]
  NP_000510.1 147 aa  hemoglobin subunit delta [Homo sapiens]
```

The search stays on NCBI's server (`WebEnv` + `query_key`), and `FetchSequencesAsync` pages through it with `retstart`/`retmax`, so no id list is ever sent. Id lists longer than 200 are POSTed automatically.

## 3. A protein from UniProt

```csharp
using Biotech.Net.Data.UniProt;

var uniprot = new UniProtClient();
var entry = await uniprot.GetEntryAsync("P69905");
Console.WriteLine(entry);
Console.WriteLine($"Reviewed: {entry.IsReviewed}, gene {entry.GeneNames[0]}, {entry.Features.Count} features");
Console.WriteLine($"PDB structures: {entry.CrossReferencesTo("PDB").Count()}, e.g. {entry.CrossReferencesTo("PDB").First().Id}");

var protein = entry.ToProtein();
Console.WriteLine($"{protein}  pI {protein.IsoelectricPoint:F2}");
```

```text
P69905 HBA_HUMAN Hemoglobin subunit alpha (Homo sapiens, 142 aa)
Reviewed: True, gene HBA1, 214 features
PDB structures: 355, e.g. 1A00
Hemoglobin subunit alpha (P69905, 142 aa, 15.3 kDa)  pI 9.09
```

`UniProtEntry` keeps the names, organism, sequence, comments, features (1-based, inclusive, as in UniProt), keywords and cross-references. `ToProtein()` gives a Core `Protein`, whose mass and pI are computed by Biotech.Net's own kernels. `GetFastaAsync` returns just the sequence.

### Search with pagination, and ID mapping

```csharp
using Biotech.Net.Data.UniProt;

var uniprot = new UniProtClient();
var page = await uniprot.SearchPageAsync("gene:HBB AND reviewed:true", fields: "accession,protein_name,organism_name,length", pageSize: 3);
Console.WriteLine($"{page.TotalResults} hits; next page: {(page.NextPage is null ? "none" : "yes")}");

await foreach (var e in uniprot.SearchAsync("gene:HBB AND reviewed:true", "accession,organism_name,length", pageSize: 3, maxResults: 5))
    Console.WriteLine($"  {e.Accession} {e.OrganismName} ({e.Length} aa)");

var mapping = await uniprot.MapIdsAsync("UniProtKB_AC-ID", "GeneID", ["P69905", "P68871"]);
foreach (var m in mapping.Mappings) Console.WriteLine($"  {m.From} -> GeneID {m.To}");
```

```text
243 hits; next page: yes
  P68871 Homo sapiens (147 aa)
  P02070 Bos taurus (145 aa)
  P02091 Rattus norvegicus (147 aa)
  P02088 Mus musculus (147 aa)
  P02089 Mus musculus (147 aa)
  P69905 -> GeneID 3039
  P69905 -> GeneID 3040
  P68871 -> GeneID 3043
```

Searches use the [UniProt query syntax](https://www.uniprot.org/help/query-fields) and follow the `Link: rel="next"` cursor; `maxResults` stops before a page nobody will read. `MapIdsAsync` submits an ID-mapping job, polls its status (1 s, growing to 5 s) and collects every page of results; a source id may map to several targets (here both alpha-globin genes).

## 4. A structure from the PDB (and AlphaFold DB)

```csharp
using Biotech.Net.Data.Pdb;
using Biotech.Net.StructuralBiology;

var pdb = new PdbClient();
Structure crambin = await pdb.GetStructureAsync("1CRN");          // mmCIF, gzip on the wire
Console.WriteLine($"{crambin.Id}: {crambin.AtomCount} atoms, chain A = {crambin["A"].ToSequence()}");

var info = await pdb.GetEntryAsync("1CRN");
Console.WriteLine($"{info.ExperimentalMethods[0]} at {info.Resolution} Å, released {info.ReleaseDate:yyyy-MM-dd}, doi:{info.CitationDoi}");

var hits = await pdb.SearchAsync(PdbSearchQuery.Sequence(crambin["A"].ToSequence().ToString(), identityCutoff: 0.9).With(rows: 3));
Console.WriteLine($"{hits.TotalCount} similar chains, top: {string.Join(", ", hits.Identifiers)}");

var model = await new AlphaFoldClient().GetStructureAsync("P01542");
Console.WriteLine($"{model.Id}: mean pLDDT {model.Atoms.Where(a => a.Name == "CA").Average(a => a.BFactor):F1}");
```

```text
1CRN: 327 atoms, chain A = TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN
X-RAY DIFFRACTION at 1.5 Å, released 1981-07-28, doi:10.1073/pnas.81.19.6014
27 similar chains, top: 1CRN_1, 1CCM_1, 1CCN_1
AF-P01542-F1: mean pLDDT 93.5
```

Coordinates arrive gzip-compressed (and are cached compressed) and are parsed by `MmCifReader` or, with `StructureFileFormat.Pdb`, `PdbReader`, so everything in [Structural biology](structural-biology.md) applies to the result. `PdbSearchQuery` builds full-text, sequence and attribute queries and `And`/`Or` groups; `SearchAsync(string)` takes raw Search API JSON, and `GraphQlAsync` runs Data API GraphQL. AlphaFold models keep the per-residue pLDDT in each atom's B-factor.

## 5. A gene and a variant effect from Ensembl

```csharp
using Biotech.Net.Core;
using Biotech.Net.Data.Ensembl;

var ensembl = new EnsemblClient();
var lookup = await ensembl.LookupSymbolAsync("homo_sapiens", "HBB", expand: true);
Gene hbb = lookup.ToGene();
Console.WriteLine($"{hbb.Name} {hbb.Id} {hbb.Location} ({lookup.AssemblyName}), {hbb.Transcripts.Count} transcripts, canonical {lookup.CanonicalTranscript}");

var vep = (await ensembl.VepHgvsAsync("human", "ENST00000335295:c.20A>T"))[0];
var tx = vep.TranscriptConsequences.First(c => c.TranscriptId == "ENST00000335295");
Console.WriteLine($"{vep.Input}: {vep.MostSevereConsequence} ({vep.MostSevere}), {tx.GeneSymbol} p.{tx.ProteinChange}, SIFT {tx.SiftPrediction}");

var region = (await ensembl.VepRegionAsync("human", new Variant("11", 5227002, "T", "A")))[0];
Console.WriteLine($"Same change as a VCF variant: {region.AlleleString} {region.MostSevereConsequence}");
```

```text
HBB ENSG00000244734 11:5,225,464-5,229,395(-) (GRCh38), 19 transcripts, canonical ENST00000335295.4
ENST00000335295:c.20A>T: missense_variant (MissenseVariant), HBB p.E7V, SIFT deleterious_low_confidence
Same change as a VCF variant: T/A missense_variant
```

That is the sickle-cell mutation (HbS, Glu6Val in the mature protein, E7V counting the initiator Met). `ToGene()` converts Ensembl's 1-based inclusive coordinates into Biotech.Net's 0-based half-open `GenomicInterval`s, with exons and the coding region per transcript. VEP consequence terms are mapped to the Genomics module's `Consequence` enum where it models them. `VepRegionAsync` takes a Core `Variant` (VCF-style, anchored indels are converted for VEP), and `VepRegionsAsync` POSTs many variants 200 at a time. `GetSequenceAsync` (genomic, cDNA, CDS, protein), `GetRegionSequenceAsync`, `OverlapRegionAsync` and `GetXrefsAsync` cover the rest; set `Server = EnsemblClient.Grch37Server` for GRCh37.

## 6. Caching and offline work

Every response that is safe to reuse is cached, keyed by a SHA-256 of the method, the normalised URL (sorted parameters, no `api_key`/`email`/`tool`) and the request body. By default the cache lives on disk under `%LOCALAPPDATA%\Biotech.Net\cache` on Windows (`~/.local/share/Biotech.Net/cache` on Linux); set `BIOTECH_NET_CACHE` to move it. Entries stay fresh for a per-endpoint time-to-live (sequence records and coordinate files 30 days, entries and lookups 7 days, searches 1 day, history-server searches 1 hour). ID-mapping jobs are never cached.

```csharp
using Biotech.Net.Data;
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.UniProt;

string dir = Path.Combine(Path.GetTempPath(), "biotech-cache-demo");
var online = new DataClientOptions { Cache = new DiskDataCache(dir) };
var crambin = await new UniProtClient(options: online).GetEntryAsync("P01542");   // network, then stored

var offline = online with { CacheMode = DataCacheMode.Offline };
var again = await new UniProtClient(options: offline).GetEntryAsync("P01542");    // served from disk
Console.WriteLine($"{again.Accession} {again.ProteinName} from cache; entries on disk: {new DiskDataCache(dir).Count}");

try
{
    await new UniProtClient(options: offline).GetEntryAsync("P69905");             // never fetched
}
catch (DataOfflineException ex)
{
    Console.WriteLine(ex.Message);
}
Console.WriteLine($"Default cache directory: {DiskDataCache.DefaultDirectory}");
```

```text
P01542 Crambin from cache; entries on disk: 1
Offline mode: no cached response for https://rest.uniprot.org/uniprotkb/P69905?format=json. Run the same call once online (or switch CacheMode away from Offline, or unset BIOTECH_NET_OFFLINE).
Default cache directory: C:\Users\<you>\AppData\Local\Biotech.Net\cache
```

| `DataCacheMode` | Behaviour |
|---|---|
| `Default` | serve a fresh cached copy, otherwise fetch and store |
| `Refresh` | always fetch, overwrite the cached copy |
| `Bypass` | never read or write the cache |
| `Offline` | never touch the network; serve any cached copy however old, else throw `DataOfflineException` |

`BIOTECH_NET_OFFLINE=1` makes `Offline` the default for a whole process (handy in CI once a cache has been primed). `DataClientOptions.CacheTtl` overrides every TTL. `MemoryDataCache` is an in-process LRU alternative, and any `IDataCache` implementation can be plugged in. The disk cache writes each entry to a temporary file and renames it into place, so concurrent processes never see a half-written entry, and it evicts least-recently-used entries past its size cap (2 GB by default).

## 7. Rate limits, retries and dependency injection

Requests to each host are spaced through a process-wide `HostRateLimiter`: NCBI 3/s (10/s with a key), Ensembl 15/s, the others 10/s. HTTP 408/429/500/502/503/504, timeouts and network errors are retried up to `MaxRetries` (4) times with exponential backoff and jitter; a `Retry-After` header is honoured (and pauses the whole host on 429), as is Ensembl's `X-RateLimit-Reset`. All waiting goes through `DataClientOptions.TimeProvider`, so tests can use a fake clock. Anything else ends in a `DataServiceException` carrying the service, status and the service's own error text:

```csharp
using Biotech.Net.Data.Http;
using Biotech.Net.Data.UniProt;

try
{
    await new UniProtClient().GetEntryAsync("NOT_AN_ACCESSION");
}
catch (DataServiceException ex)
{
    Console.WriteLine($"{ex.Service} {(int?)ex.StatusCode}: {ex.ServiceMessage}");
}
```

```text
UniProt 400: The 'accession' value has invalid format. It should be a valid UniProtKB accession
```

Every client takes an optional `HttpClient` and `DataClientOptions`, so it works as an `IHttpClientFactory` typed client (this snippet adds the `Microsoft.Extensions.Http` package):

```csharp
using Biotech.Net.Data.Caching;
using Biotech.Net.Data.Ensembl;
using Biotech.Net.Data.Http;
using Microsoft.Extensions.DependencyInjection;

var services = new ServiceCollection();
services.AddSingleton(new DataClientOptions { Cache = new MemoryDataCache(), MaxRetries = 2, Timeout = TimeSpan.FromSeconds(30) });
services.AddHttpClient<EnsemblClient>();   // typed client: HttpClient from the factory, options from the container
using var provider = services.BuildServiceProvider();

var ensembl = provider.GetRequiredService<EnsemblClient>();
var brca2 = await ensembl.LookupSymbolAsync("homo_sapiens", "BRCA2");
Console.WriteLine($"{brca2.DisplayName} {brca2.Id} chr{brca2.SeqRegionName}:{brca2.Start:N0}-{brca2.End:N0}");
```

```text
BRCA2 ENSG00000139618 chr13:32,315,086-32,400,268
```

## Limitations

- **Coverage.** NCBI: EInfo, ESearch, ESummary, EFetch, ELink and the history server; EPost, EGQuery, ESpell and ECitMatch are not wrapped, and EFetch output other than FASTA/GenBank is returned as raw text (`EFetchAsync`). UniProt: UniProtKB only (no UniRef, UniParc or Proteomes endpoints, no streaming endpoint). PDB: RCSB files, Data API (entry, polymer entity, GraphQL) and Search API; PDBe and PDBj are not wrapped; AlphaFold DB covers metadata and model files, not PAE matrices. Ensembl: lookup, sequence, overlap, VEP, xrefs; no comparative genomics, regulation or variation endpoints beyond VEP and overlap.
- **Freshness.** A cached response is served until its TTL expires, even if the database has changed; use `Refresh` mode or a shorter `CacheTtl` when you need the latest release. Offline mode deliberately serves stale data.
- **History-server sessions** (`WebEnv`) expire on NCBI's side after a few hours; an `ESearch` with `useHistory` is therefore cached for only one hour.
- **Rate limits are per process.** The shared limiter cannot coordinate several processes on one machine; NCBI may still answer 429, which is retried.
- **Very large structures** have no legacy PDB-format file; ask for `StructureFileFormat.MmCif` (the default).
- **The offline test suite replays recorded responses** (`tests/Biotech.Net.Data.Tests/Fixtures`, with the script that recorded them). Live tests run only with `BIOTECH_NET_LIVE_TESTS=1`.

---
Previous: [Notebooks](notebooks.md) · Next: [Workflows](workflows.md) · [Tutorial index](index.md)
