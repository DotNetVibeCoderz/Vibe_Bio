# API reference: Data

Package `Gravicode.Biotech.Net.Data` · namespaces `Biotech.Net.Data`, `.Http`, `.Caching`, `.Ncbi`, `.UniProt`, `.Pdb`, `.Ensembl`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Depends on Core, Sequences, Genomics and StructuralBiology. Every client constructor is `(HttpClient? httpClient = null, DataClientOptions? options = null)`, every call is async and takes a `CancellationToken`.

## BioDataClient

`new BioDataClient(httpClient, options)` exposes `Entrez`, `UniProt`, `Pdb`, `AlphaFold` and `Ensembl`, all sharing one `HttpClient` and one set of options.

## Engine: `Biotech.Net.Data.Http`

| Type | Members |
|---|---|
| `DataClientOptions` (record) | `Cache` (null = `DataCache.Default`), `CacheMode`, `CacheTtl` (overrides per-endpoint TTLs), `TimeProvider`, `MaxRetries` (4), `InitialBackoff` (1 s), `MaxBackoff` (30 s), `MaxRetryAfter` (2 min), `Timeout` (100 s per attempt), `RateLimiter` (null = `HostRateLimiter.Shared`), `RequestsPerSecond` (override), `Random` (jitter), `UserAgent`; `Default`; `OfflineEnvironmentVariable` = `BIOTECH_NET_OFFLINE` |
| `DataCacheMode` | `Default`, `Refresh`, `Bypass`, `Offline` |
| `HostRateLimiter` | `Shared`, `WaitAsync(host, requestsPerSecond, timeProvider, ct)`, `BlockUntil(host, until)`, `BlockedUntil(host)` |
| `DataServiceException` | `Service`, `RequestUri` (credentials removed), `StatusCode`, `ServiceMessage`, `IsNotFound` |
| `DataOfflineException` | thrown in `Offline` mode when nothing is cached (or the request is never cached) |

Retries: 408, 429, 500, 502, 503, 504, timeouts and network errors; exponential backoff with equal jitter (`[d/2, d]`, `d = min(MaxBackoff, InitialBackoff·2^attempt)`); `Retry-After` honoured, 429 pauses the host; `X-RateLimit-Remaining: 0` + `X-RateLimit-Reset` pauses the host. Default rates: NCBI 3/s (10/s with a key), Ensembl 15/s, UniProt, RCSB and AlphaFold 10/s.

## Caching: `Biotech.Net.Data.Caching`

| Type | Members |
|---|---|
| `IDataCache` | `GetAsync`, `SetAsync`, `RemoveAsync`, `ClearAsync` |
| `CacheEntry` | `Content` (raw body), `StoredAt`, `Url`, `ContentType`, `Headers` (`Link`, `X-Total-Results`, ...), `Age(now)` |
| `CacheKey` | `Compute(method, uri, body)` (SHA-256 hex), `NormalizeUrl(uri)`, `IsCredential(name)` (`api_key`, `email`, `tool`) |
| `DiskDataCache` | `(directory = DefaultDirectory, maxBytes = 2 GB)`; `DefaultDirectory` (`BIOTECH_NET_CACHE` or `LocalApplicationData/Biotech.Net/cache`), `PathFor(key)`, `Count`, `SizeBytes`, `TrimAsync()`; atomic temp-file + rename writes, LRU eviction past the cap |
| `MemoryDataCache` | `(maxEntries = 1024, maxBytes = 256 MB)`; LRU; `Count`, `SizeBytes` |
| `DataCache` | `Default`, the lazily created disk cache |

Default TTLs: sequence records, coordinate files and AlphaFold 30 days; UniProt entries, ESummary, ELink, RCSB Data API and Ensembl 7 days; searches and EInfo 1 day; history-server ESearch 1 hour; ID-mapping jobs never.

## NCBI: `EntrezClient` (`Biotech.Net.Data.Ncbi`)

Init properties `ApiKey` (`NCBI_API_KEY`), `Email` (`NCBI_EMAIL`), `Tool` (`Biotech.Net`), `BaseUrl`; `RequestsPerSecond`.

| Method | Returns |
|---|---|
| `EInfoDatabasesAsync()` | database names |
| `EInfoAsync(db)` | `EInfoDatabase` (count, last update, `Fields`, `Links`) |
| `ESearchAsync(db, term, retMax, retStart, useHistory, sort)` | `ESearchResult` (`Count`, `Ids`, `QueryTranslation`, `History`) |
| `ESummaryAsync(db, ids)` / `ESummaryAsync(db, history, retStart, retMax)` | `ESummaryDocument` (`Uid`, `Fields` JSON, `GetString`) |
| `PubMedSummariesAsync(pmids)` | `PubMedSummary` (title, journal, date, authors, DOI, volume, issue, pages) |
| `EFetchAsync(db, ids, retType, retMode)` / `EFetchAsync(db, history, ...)` | raw text |
| `FetchSequencesAsync(db, ids, EntrezSequenceFormat)` / `FetchSequenceAsync(db, id, format)` | `Sequence` (FASTA, or GenBank with features) |
| `FetchSequencesAsync(ESearchResult, format, batchSize)` | `IAsyncEnumerable<Sequence>` over a history-server search |
| `ELinkAsync(dbFrom, db, ids, linkName)` | `ELinkSet` → `ELinkResult` (`DbTo`, `LinkName`, `Ids`) |

Id lists longer than 200 are POSTed.

## UniProt: `UniProtClient` (`Biotech.Net.Data.UniProt`)

| Method | Returns |
|---|---|
| `GetEntryAsync(accession)` | `UniProtEntry` |
| `GetProteinAsync(accession)` | Core `Protein` |
| `GetFastaAsync(accession)` | `Sequence` |
| `SearchPageAsync(query, fields, pageSize)` / `NextPageAsync(page)` | `UniProtSearchPage` (`Entries`, `TotalResults`, `NextPage`) |
| `SearchAsync(query, fields, pageSize, maxResults)` | `IAsyncEnumerable<UniProtEntry>` following `Link: rel="next"` |
| `SubmitIdMappingAsync(from, to, ids, taxonId)` | job id |
| `GetIdMappingStatusAsync(jobId)` | `IdMappingJobState` |
| `GetIdMappingResultsAsync(jobId)` | `IdMappingResult` (`Mappings`, `FailedIds`, `BySource()`) |
| `MapIdsAsync(from, to, ids, taxonId)` | submit + poll + collect |

`UniProtEntry`: `Accession`, `SecondaryAccessions`, `EntryName`, `IsReviewed`, `ProteinName`, `AlternativeNames`, `GeneNames`, `OrganismName`, `OrganismCommonName`, `TaxonId`, `Lineage`, `ProteinExistence`, `AnnotationScore`, `EntryVersion`, `SequenceVersion`, `SequenceResidues`, `Length`, `MolecularWeight`, `Crc64`, `Comments`, `Features` (1-based inclusive), `Keywords`, `CrossReferences`, `Function`, `CrossReferencesTo(db)`, `ToSequence()`, `ToProtein()`, `FromJson(JsonElement)`. Settings: `IdMappingPollInterval` (1 s), `IdMappingTimeout` (5 min).

## PDB: `PdbClient`, `AlphaFoldClient` (`Biotech.Net.Data.Pdb`)

| Method | Returns |
|---|---|
| `GetStructureAsync(pdbId, StructureFileFormat)` | `Structure` (mmCIF by default; gzip download) |
| `GetStructureTextAsync(pdbId, format)` | decompressed file text |
| `GetEntryAsync(pdbId)` | `PdbEntryInfo` (title, methods, resolution, dates, counts, weight, citation) |
| `GetPolymerEntityAsync(pdbId, entityId)` | `PdbPolymerEntity` (sequence, type, chains, UniProt ids, organisms; `ToSequence()`) |
| `GraphQlAsync(query, variables)` | the `data` member as `JsonElement` |
| `SearchAsync(PdbSearchQuery)` / `SearchAsync(json)` | `PdbSearchResult` (`TotalCount`, `ResultType`, `Hits`, `Identifiers`) |
| `AlphaFoldClient.GetPredictionsAsync(uniprot)` | `AlphaFoldPrediction` (entry id, sequence, mean pLDDT, file URLs) |
| `AlphaFoldClient.GetStructureAsync(uniprot or prediction, format)` | `Structure` (pLDDT in B-factors) |

`PdbSearchQuery`: `FullText(text)`, `Sequence(residues, type, evalueCutoff, identityCutoff)` (returns polymer entities), `Attribute(attribute, operator, value)`, `And(...)`, `Or(...)`, `With(returnType, start, rows, includeComputedModels)`, `ToJson()`.

## Ensembl: `EnsemblClient` (`Biotech.Net.Data.Ensembl`)

Init property `Server` (`DefaultServer`, `Grch37Server`).

| Method | Returns |
|---|---|
| `LookupIdAsync(id, expand)` / `LookupSymbolAsync(species, symbol, expand)` | `EnsemblFeature` (`ToInterval()`, `ToGene()`, `Transcripts`, `Exons`, `Translation`) |
| `GetSequenceAsync(id, EnsemblSequenceType)` / `GetSequencesAsync(id, type)` | `Sequence` (genomic, cDNA, CDS, protein) |
| `GetRegionSequenceAsync(species, GenomicInterval)` | `Sequence` (reverse strand = reverse complement) |
| `OverlapRegionAsync(species, region, features)` | `EnsemblOverlapFeature` |
| `VepHgvsAsync(species, hgvs)` | `VepResult` |
| `VepRegionAsync(species, Variant)` / `VepRegionsAsync(species, variants)` | `VepResult` (POST in batches of 200) |
| `GetXrefsAsync(id, externalDb)` / `XrefsBySymbolAsync(species, symbol)` | `EnsemblXref` / `EnsemblIdentity` |
| `ToVepRegion(chrom, pos, ref, alt)` (static) | VEP region notation from VCF-style alleles |

`VepResult`: `Input`, `Id`, `SeqRegionName`, `Start`, `End`, `AlleleString`, `MostSevereConsequence`, `MostSevere` (Genomics `Consequence?`), `TranscriptConsequences`, `Canonical`, `ColocatedVariantIds`. `VepTranscriptConsequence`: gene and transcript ids, `ConsequenceTerms`, `Consequences`, `Impact`, `AminoAcids`, `Codons`, `ProteinChange`, protein/CDS positions, SIFT and PolyPhen, `IsCanonical`.

Ensembl coordinates in these models are 1-based inclusive, as Ensembl reports them; `ToInterval()` and `ToGene()` convert to Biotech.Net's 0-based half-open intervals.

---
[API reference index](index.md) · [Cloud data tutorial](../tutorials/data.md)
