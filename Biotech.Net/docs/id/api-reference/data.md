# Referensi API: Data

Paket `Gravicode.Biotech.Net.Data` · namespace `Biotech.Net.Data`, `.Http`, `.Caching`, `.Ncbi`, `.UniProt`, `.Pdb`, `.Ensembl`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Bergantung pada Core, Sequences, Genomics, dan StructuralBiology. Konstruktor setiap klien berbentuk `(HttpClient? httpClient = null, DataClientOptions? options = null)`. Setiap pemanggilan bersifat async dan menerima `CancellationToken`.

## BioDataClient

`new BioDataClient(httpClient, options)` menyediakan `Entrez`, `UniProt`, `Pdb`, `AlphaFold`, dan `Ensembl`. Kelimanya berbagi satu `HttpClient` dan satu set opsi.

## Mesin: `Biotech.Net.Data.Http`

| Tipe | Anggota |
|---|---|
| `DataClientOptions` (record) | `Cache` (null = `DataCache.Default`), `CacheMode`, `CacheTtl` (menimpa TTL per endpoint), `TimeProvider`, `MaxRetries` (4), `InitialBackoff` (1 detik), `MaxBackoff` (30 detik), `MaxRetryAfter` (2 menit), `Timeout` (100 detik per percobaan), `RateLimiter` (null = `HostRateLimiter.Shared`), `RequestsPerSecond` (penimpa), `Random` (jitter), `UserAgent`; `Default`; `OfflineEnvironmentVariable` = `BIOTECH_NET_OFFLINE` |
| `DataCacheMode` | `Default`, `Refresh`, `Bypass`, `Offline` |
| `HostRateLimiter` | `Shared`, `WaitAsync(host, requestsPerSecond, timeProvider, ct)`, `BlockUntil(host, until)`, `BlockedUntil(host)` |
| `DataServiceException` | `Service`, `RequestUri` (kredensial dibuang), `StatusCode`, `ServiceMessage`, `IsNotFound` |
| `DataOfflineException` | dilempar dalam mode `Offline` bila respons tidak ada di cache (atau permintaannya memang tidak pernah di-cache) |

Pengulangan dilakukan untuk 408, 429, 500, 502, 503, 504, timeout, dan galat jaringan. Backoff-nya eksponensial dengan equal jitter (`[d/2, d]`, `d = min(MaxBackoff, InitialBackoff·2^attempt)`). `Retry-After` dihormati, dan 429 menjeda seluruh host. Kombinasi `X-RateLimit-Remaining: 0` + `X-RateLimit-Reset` juga menjeda host. Laju bawaan: NCBI 3/detik (10/detik dengan kunci), Ensembl 15/detik, UniProt, RCSB, dan AlphaFold 10/detik.

## Cache: `Biotech.Net.Data.Caching`

| Tipe | Anggota |
|---|---|
| `IDataCache` | `GetAsync`, `SetAsync`, `RemoveAsync`, `ClearAsync` |
| `CacheEntry` | `Content` (isi mentah), `StoredAt`, `Url`, `ContentType`, `Headers` (`Link`, `X-Total-Results`, ...), `Age(now)` |
| `CacheKey` | `Compute(method, uri, body)` (hex SHA-256), `NormalizeUrl(uri)`, `IsCredential(name)` (`api_key`, `email`, `tool`) |
| `DiskDataCache` | `(directory = DefaultDirectory, maxBytes = 2 GB)`; `DefaultDirectory` (`BIOTECH_NET_CACHE` atau `LocalApplicationData/Biotech.Net/cache`), `PathFor(key)`, `Count`, `SizeBytes`, `TrimAsync()`. Penulisan atomik (berkas sementara + ganti nama); entri LRU dibuang bila ukuran melewati batas |
| `MemoryDataCache` | `(maxEntries = 1024, maxBytes = 256 MB)`; LRU; `Count`, `SizeBytes` |
| `DataCache` | `Default`, yaitu cache disk yang dibuat saat pertama dipakai |

TTL bawaan:

- Rekaman sekuens, berkas koordinat, dan AlphaFold: 30 hari.
- Entri UniProt, ESummary, ELink, Data API RCSB, dan Ensembl: 7 hari.
- Pencarian dan EInfo: 1 hari.
- ESearch dengan history server: 1 jam.
- Job pemetaan ID: tidak pernah di-cache.

## NCBI: `EntrezClient` (`Biotech.Net.Data.Ncbi`)

Properti init: `ApiKey` (`NCBI_API_KEY`), `Email` (`NCBI_EMAIL`), `Tool` (`Biotech.Net`), `BaseUrl`. Ada juga `RequestsPerSecond`.

| Metode | Hasil |
|---|---|
| `EInfoDatabasesAsync()` | nama-nama basis data |
| `EInfoAsync(db)` | `EInfoDatabase` (jumlah rekaman, pembaruan terakhir, `Fields`, `Links`) |
| `ESearchAsync(db, term, retMax, retStart, useHistory, sort)` | `ESearchResult` (`Count`, `Ids`, `QueryTranslation`, `History`) |
| `ESummaryAsync(db, ids)` / `ESummaryAsync(db, history, retStart, retMax)` | `ESummaryDocument` (`Uid`, `Fields` JSON, `GetString`) |
| `PubMedSummariesAsync(pmids)` | `PubMedSummary` (judul, jurnal, tanggal, penulis, DOI, volume, nomor, halaman) |
| `EFetchAsync(db, ids, retType, retMode)` / `EFetchAsync(db, history, ...)` | teks mentah |
| `FetchSequencesAsync(db, ids, EntrezSequenceFormat)` / `FetchSequenceAsync(db, id, format)` | `Sequence` (FASTA, atau GenBank lengkap dengan fitur) |
| `FetchSequencesAsync(ESearchResult, format, batchSize)` | `IAsyncEnumerable<Sequence>` atas hasil pencarian di history server |
| `ELinkAsync(dbFrom, db, ids, linkName)` | `ELinkSet` → `ELinkResult` (`DbTo`, `LinkName`, `Ids`) |

Daftar id yang lebih dari 200 dikirim lewat POST.

## UniProt: `UniProtClient` (`Biotech.Net.Data.UniProt`)

| Metode | Hasil |
|---|---|
| `GetEntryAsync(accession)` | `UniProtEntry` |
| `GetProteinAsync(accession)` | `Protein` dari Core |
| `GetFastaAsync(accession)` | `Sequence` |
| `SearchPageAsync(query, fields, pageSize)` / `NextPageAsync(page)` | `UniProtSearchPage` (`Entries`, `TotalResults`, `NextPage`) |
| `SearchAsync(query, fields, pageSize, maxResults)` | `IAsyncEnumerable<UniProtEntry>` yang mengikuti `Link: rel="next"` |
| `SubmitIdMappingAsync(from, to, ids, taxonId)` | id job |
| `GetIdMappingStatusAsync(jobId)` | `IdMappingJobState` |
| `GetIdMappingResultsAsync(jobId)` | `IdMappingResult` (`Mappings`, `FailedIds`, `BySource()`) |
| `MapIdsAsync(from, to, ids, taxonId)` | kirim job + pantau status + kumpulkan hasil |

`UniProtEntry` punya anggota berikut:

- Identitas: `Accession`, `SecondaryAccessions`, `EntryName`, `IsReviewed`, `ProteinName`, `AlternativeNames`, `GeneNames`.
- Organisme: `OrganismName`, `OrganismCommonName`, `TaxonId`, `Lineage`.
- Status entri: `ProteinExistence`, `AnnotationScore`, `EntryVersion`, `SequenceVersion`.
- Sekuens: `SequenceResidues`, `Length`, `MolecularWeight`, `Crc64`.
- Anotasi: `Comments`, `Features` (berbasis 1, inklusif), `Keywords`, `CrossReferences`, `Function`, `CrossReferencesTo(db)`.
- Konversi: `ToSequence()`, `ToProtein()`, `FromJson(JsonElement)`.

Pengaturan klien: `IdMappingPollInterval` (1 detik), `IdMappingTimeout` (5 menit).

## PDB: `PdbClient`, `AlphaFoldClient` (`Biotech.Net.Data.Pdb`)

| Metode | Hasil |
|---|---|
| `GetStructureAsync(pdbId, StructureFileFormat)` | `Structure` (bawaan mmCIF; diunduh dalam gzip) |
| `GetStructureTextAsync(pdbId, format)` | teks berkas yang sudah didekompresi |
| `GetEntryAsync(pdbId)` | `PdbEntryInfo` (judul, metode, resolusi, tanggal, jumlah, bobot, sitasi) |
| `GetPolymerEntityAsync(pdbId, entityId)` | `PdbPolymerEntity` (sekuens, tipe, rantai, id UniProt, organisme; `ToSequence()`) |
| `GraphQlAsync(query, variables)` | anggota `data` sebagai `JsonElement` |
| `SearchAsync(PdbSearchQuery)` / `SearchAsync(json)` | `PdbSearchResult` (`TotalCount`, `ResultType`, `Hits`, `Identifiers`) |
| `AlphaFoldClient.GetPredictionsAsync(uniprot)` | `AlphaFoldPrediction` (id entri, sekuens, rerata pLDDT, URL berkas) |
| `AlphaFoldClient.GetStructureAsync(uniprot atau prediction, format)` | `Structure` (pLDDT di B-factor) |

`PdbSearchQuery`:

- `FullText(text)`.
- `Sequence(residues, type, evalueCutoff, identityCutoff)`, yang mengembalikan polymer entity.
- `Attribute(attribute, operator, value)`.
- `And(...)` dan `Or(...)`.
- `With(returnType, start, rows, includeComputedModels)`.
- `ToJson()`.

## Ensembl: `EnsemblClient` (`Biotech.Net.Data.Ensembl`)

Properti init: `Server` (`DefaultServer`, `Grch37Server`).

| Metode | Hasil |
|---|---|
| `LookupIdAsync(id, expand)` / `LookupSymbolAsync(species, symbol, expand)` | `EnsemblFeature` (`ToInterval()`, `ToGene()`, `Transcripts`, `Exons`, `Translation`) |
| `GetSequenceAsync(id, EnsemblSequenceType)` / `GetSequencesAsync(id, type)` | `Sequence` (genomik, cDNA, CDS, protein) |
| `GetRegionSequenceAsync(species, GenomicInterval)` | `Sequence` (untai balik = komplemen balik) |
| `OverlapRegionAsync(species, region, features)` | `EnsemblOverlapFeature` |
| `VepHgvsAsync(species, hgvs)` | `VepResult` |
| `VepRegionAsync(species, Variant)` / `VepRegionsAsync(species, variants)` | `VepResult` (POST per 200 varian) |
| `GetXrefsAsync(id, externalDb)` / `XrefsBySymbolAsync(species, symbol)` | `EnsemblXref` / `EnsemblIdentity` |
| `ToVepRegion(chrom, pos, ref, alt)` (statis) | notasi region VEP dari alel gaya VCF |

`VepResult` punya anggota berikut:

- Masukan dan posisi: `Input`, `Id`, `SeqRegionName`, `Start`, `End`, `AlleleString`.
- Konsekuensi: `MostSevereConsequence`, `MostSevere` (`Consequence?` dari Genomics), `TranscriptConsequences`, `Canonical`.
- Varian lain di posisi yang sama: `ColocatedVariantIds`.

`VepTranscriptConsequence` punya anggota berikut:

- Id gen dan transkrip.
- Konsekuensi: `ConsequenceTerms`, `Consequences`, `Impact`.
- Perubahan protein: `AminoAcids`, `Codons`, `ProteinChange`, posisi protein dan CDS.
- Prediksi: SIFT dan PolyPhen.
- `IsCanonical`.

Koordinat Ensembl di model-model ini berbasis 1 dan inklusif, sesuai yang dilaporkan Ensembl. `ToInterval()` dan `ToGene()` mengubahnya menjadi interval Biotech.Net yang berbasis 0 dan setengah terbuka.

---
[Indeks referensi API](index.md) · [Tutorial data cloud](../tutorials/data.md)
