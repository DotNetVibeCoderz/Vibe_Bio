# Tutorial: Konektor data cloud

Paket: `Gravicode.Biotech.Net.Data`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Ambil rekaman langsung dari basis data publik ke dalam tipe Biotech.Net: rekaman GenBank dan FASTA dari NCBI, protein dari UniProt, struktur dari RCSB PDB dan AlphaFold DB, serta gen dan efek varian dari Ensembl. Semua konektor memakai satu mesin yang sama: mengatur jarak permintaan ke setiap host, mengulang kegagalan sementara, dan menyimpan respons di cache disk, sehingga skrip yang dijalankan kedua kalinya tetap jalan tanpa jaringan.

Setiap cuplikan di bawah dijalankan terhadap layanan sungguhan; keluaran yang ditampilkan berasal dari run tersebut (jumlah rekaman di basis data bisa bergeser seiring rilis baru).

## 1. Persiapan

```text
dotnet add package Gravicode.Biotech.Net.Data
```

`BioDataClient` menggabungkan lima klien (`Entrez`, `UniProt`, `Pdb`, `AlphaFold`, `Ensembl`); masing-masing juga bisa dibuat sendiri. NCBI meminta setiap perangkat lunak memperkenalkan diri: isi variabel lingkungan `NCBI_EMAIL` dengan alamat kontak Anda, dan `NCBI_API_KEY` bila Anda punya kunci (batasnya naik dari 3 menjadi 10 permintaan per detik). Kedua nilai itu tidak pernah ditulis ke cache maupun ke pesan exception.

## 2. Rekaman GenBank dari NCBI

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

Rekaman diurai oleh `GenBankReader` yang sama dengan untuk berkas lokal, jadi anotasi dan tabel fitur ikut terbawa. `EntrezSequenceFormat.Fasta` (bawaan) menghasilkan sekuens polos.

### Pencarian, ringkasan, dan tautan

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

`ESummaryAsync` mengembalikan semua kolom ringkasan dari basis data apa pun sebagai JSON (`ESummaryDocument.Fields`); `PubMedSummariesAsync` memberi tipe pada kolom PubMed. `EInfoAsync("pubmed")` mendaftar kolom yang bisa dicari dan tautan sebuah basis data.

### Hasil besar: history server

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

Hasil pencarian tetap tersimpan di server NCBI (`WebEnv` + `query_key`), dan `FetchSequencesAsync` membacanya per halaman dengan `retstart`/`retmax`, sehingga daftar id tidak pernah dikirim. Daftar id yang lebih dari 200 otomatis dikirim lewat POST.

## 3. Protein dari UniProt

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

`UniProtEntry` menyimpan nama, organisme, sekuens, komentar, fitur (berbasis 1 dan inklusif, seperti di UniProt), kata kunci, dan rujukan silang. `ToProtein()` menghasilkan `Protein` dari Core, yang massa dan pI-nya dihitung oleh kernel Biotech.Net sendiri. `GetFastaAsync` hanya mengembalikan sekuensnya.

### Pencarian berhalaman dan pemetaan ID

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

Pencarian memakai [sintaks kueri UniProt](https://www.uniprot.org/help/query-fields) dan mengikuti kursor `Link: rel="next"`; `maxResults` berhenti sebelum mengambil halaman yang tidak akan dibaca. `MapIdsAsync` mengirim job pemetaan ID, memantau statusnya (1 detik, membesar hingga 5 detik), lalu mengumpulkan semua halaman hasil; satu id sumber bisa terpetakan ke beberapa target (di sini kedua gen alfa-globin).

## 4. Struktur dari PDB (dan AlphaFold DB)

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

Koordinat dikirim dalam bentuk gzip (dan disimpan di cache tetap terkompresi), lalu diurai oleh `MmCifReader` atau, dengan `StructureFileFormat.Pdb`, oleh `PdbReader`, sehingga semua yang ada di [Biologi struktural](structural-biology.md) berlaku untuk hasilnya. `PdbSearchQuery` menyusun kueri teks bebas, sekuens, dan atribut serta grup `And`/`Or`; `SearchAsync(string)` menerima JSON Search API mentah, dan `GraphQlAsync` menjalankan GraphQL Data API. Model AlphaFold menyimpan pLDDT per residu di B-factor setiap atom.

## 5. Gen dan efek varian dari Ensembl

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

Itulah mutasi sel sabit (HbS, Glu6Val pada protein matang, E7V bila Met awal ikut dihitung). `ToGene()` mengubah koordinat Ensembl yang berbasis 1 dan inklusif menjadi `GenomicInterval` Biotech.Net yang berbasis 0 dan setengah terbuka, lengkap dengan ekson dan daerah penyandi per transkrip. Istilah konsekuensi VEP dipetakan ke enum `Consequence` milik modul Genomics bila istilah itu dimodelkan. `VepRegionAsync` menerima `Variant` dari Core (gaya VCF; indel dengan basa jangkar dikonversi untuk VEP), dan `VepRegionsAsync` mengirim banyak varian lewat POST, 200 per permintaan. `GetSequenceAsync` (genomik, cDNA, CDS, protein), `GetRegionSequenceAsync`, `OverlapRegionAsync`, dan `GetXrefsAsync` melengkapi sisanya; isi `Server = EnsemblClient.Grch37Server` untuk GRCh37.

## 6. Cache dan kerja offline

Setiap respons yang aman dipakai ulang disimpan di cache, dengan kunci SHA-256 dari metode, URL yang dinormalkan (parameter diurutkan, tanpa `api_key`/`email`/`tool`), dan isi permintaan. Secara bawaan cache berada di disk, di `%LOCALAPPDATA%\Biotech.Net\cache` pada Windows (`~/.local/share/Biotech.Net/cache` pada Linux); isi `BIOTECH_NET_CACHE` untuk memindahkannya. Entri dianggap segar selama time-to-live per endpoint (rekaman sekuens dan berkas koordinat 30 hari, entri dan lookup 7 hari, pencarian 1 hari, pencarian dengan history server 1 jam). Job pemetaan ID tidak pernah di-cache.

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

| `DataCacheMode` | Perilaku |
|---|---|
| `Default` | pakai salinan cache yang masih segar, kalau tidak ada ambil dari jaringan lalu simpan |
| `Refresh` | selalu ambil dari jaringan dan timpa salinan cache |
| `Bypass` | tidak pernah membaca atau menulis cache |
| `Offline` | tidak pernah menyentuh jaringan; pakai salinan cache apa pun (setua apa pun), kalau tidak ada lempar `DataOfflineException` |

`BIOTECH_NET_OFFLINE=1` menjadikan `Offline` bawaan untuk seluruh proses (berguna di CI setelah cache diisi). `DataClientOptions.CacheTtl` menimpa semua TTL. `MemoryDataCache` adalah alternatif LRU di dalam proses, dan implementasi `IDataCache` apa pun bisa dipasang. Cache disk menulis setiap entri ke berkas sementara lalu mengganti namanya ke tempat tujuan, sehingga proses yang berjalan bersamaan tidak pernah melihat entri setengah jadi, dan entri yang paling lama tidak dipakai dibuang bila ukurannya melewati batas (bawaan 2 GB).

## 7. Batas laju, pengulangan, dan dependency injection

Permintaan ke setiap host diberi jarak oleh `HostRateLimiter` bersama untuk seluruh proses: NCBI 3/detik (10/detik dengan kunci), Ensembl 15/detik, lainnya 10/detik. HTTP 408/429/500/502/503/504, timeout, dan galat jaringan diulang hingga `MaxRetries` (4) kali dengan backoff eksponensial plus jitter; header `Retry-After` dihormati (dan menjeda seluruh host pada 429), begitu pula `X-RateLimit-Reset` dari Ensembl. Semua penantian lewat `DataClientOptions.TimeProvider`, sehingga pengujian bisa memakai jam palsu. Selain itu, kegagalan berakhir sebagai `DataServiceException` yang membawa nama layanan, status, dan pesan galat dari layanan itu sendiri:

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

Setiap klien menerima `HttpClient` dan `DataClientOptions` opsional, sehingga bisa dipakai sebagai typed client `IHttpClientFactory` (cuplikan ini menambahkan paket `Microsoft.Extensions.Http`):

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

## Batasan

- **Cakupan.** NCBI: EInfo, ESearch, ESummary, EFetch, ELink, dan history server; EPost, EGQuery, ESpell, dan ECitMatch tidak dibungkus, dan keluaran EFetch selain FASTA/GenBank dikembalikan sebagai teks mentah (`EFetchAsync`). UniProt: hanya UniProtKB (tanpa endpoint UniRef, UniParc, Proteomes, maupun endpoint streaming). PDB: berkas RCSB, Data API (entry, polymer entity, GraphQL), dan Search API; PDBe dan PDBj tidak dibungkus; AlphaFold DB mencakup metadata dan berkas model, bukan matriks PAE. Ensembl: lookup, sequence, overlap, VEP, xrefs; tanpa endpoint genomika komparatif, regulasi, atau variasi selain VEP dan overlap.
- **Kesegaran data.** Respons di cache dipakai sampai TTL-nya habis walaupun basis datanya sudah berubah; gunakan mode `Refresh` atau `CacheTtl` yang lebih pendek bila butuh rilis terbaru. Mode offline memang sengaja menyajikan data lama.
- **Sesi history server** (`WebEnv`) kedaluwarsa di sisi NCBI setelah beberapa jam; karena itu `ESearch` dengan `useHistory` hanya di-cache satu jam.
- **Batas laju berlaku per proses.** Limiter bersama tidak bisa mengoordinasikan beberapa proses di satu mesin; NCBI tetap bisa menjawab 429, yang kemudian diulang.
- **Struktur yang sangat besar** tidak punya berkas format PDB lama; minta `StructureFileFormat.MmCif` (bawaan).
- **Rangkaian uji offline memutar ulang respons yang direkam** (`tests/Biotech.Net.Data.Tests/Fixtures`, lengkap dengan skrip perekamnya). Uji live hanya berjalan dengan `BIOTECH_NET_LIVE_TESTS=1`.

---
Sebelumnya: [Notebook](notebooks.md) · Selanjutnya: [Workflow](workflows.md) · [Indeks tutorial](index.md)
