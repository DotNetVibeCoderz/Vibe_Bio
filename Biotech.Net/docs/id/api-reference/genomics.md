# Referensi API: Genomics

Paket `Gravicode.Biotech.Net.Genomics` · namespace `Biotech.Net.Genomics.Vcf`, `.IO`, `.Mapping`, `.Cram`, `.Annotation`, `.Calling`, `.Assembly`, `.GeneFinding`, `.Population`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## VCF (`.Vcf`)

| Tipe | Anggota utama |
|---|---|
| `VcfReader` | `Open(path atau stream)`, `Open(path, indexPath)`, `Parse(text)`, `ReadRecords(path)`, `Header`, `Records()`, `Variants()`, `Query(GenomicInterval)`, `Query("chr1:100-200")`, `IsIndexed` (berkas `.vcf.gz` yang punya `.tbi`/`.csi` di sebelahnya dikueri lewat indeks) |
| `VcfWriter` | `Create(path, header)` (path `.gz` ditulis sebagai BGZF, siap diindeks tabix), `Write(record / variant / records)`, `ToString(header, records)`, `Flush` |
| `VcfHeader` | `FileFormat`, `Info`, `Format`, `Filters`, `Contigs`, `Samples`, `AddInfo`, `AddFormat`, `AddFilter`, `AddContig`, `AddStandardGenotypeFields`, `ParseMetaLine`, `Write` |
| `VcfRecord` | `Variant`, `Format`, `SampleCount`, `GetGenotype(sample)`, `GetField`, `SetField`, `AlleleDepths`, `Depth`, `GenotypeQuality`, `PhredLikelihoods`, `SplitMultiallelic()`, `GenotypeIndex(j, k)` |
| `VcfFilters` (ekstensi) | `MinQuality`, `PassOnly`, `InRegion`, `InRegions`, `BiallelicSnvs` |

## BGZF, indeks binning, dan tabix (`.IO`)

| Tipe | Anggota utama |
|---|---|
| `BgzfReader` (sebuah `Stream`) | `Open(path)`, `IsBgzf(path atau stream)`, `VirtualPosition`, `Seek(VirtualOffset)`, `ReadLine(ByteBuffer)`, `AtEnd`, `VerifyChecksums`. Memeriksa CRC32 dan ISIZE setiap blok; menyimpan cache LRU blok yang sudah didekompresi |
| `BgzfWriter` (sebuah `Stream`) | `Create(path, level)`, `VirtualPosition`, `EnsureSpace(bytes)`, `FlushBlock()`, `EofMarker`. Blok berukuran ≤ 65 280 byte, penanda EOF ditulis saat dispose |
| `VirtualOffset` | `BlockAddress`, `BlockOffset`, operator perbandingan |
| `BinningIndex` | `Load(path)` / `Read(stream)` untuk `.bai` dan `.csi`, `Save`, `Write`, `Query(refId, beg, end)` → `Chunk` yang sudah digabung, `UnplacedStart()`, `References` (`Bins`, `Linear`, `Span`, `MappedCount`, `UnmappedCount`), `UnplacedCount`, `MinShift`, `Depth`, `RegionToBin`, `RegionToBins` |
| `BinningIndexBuilder` | `Push(refId, beg, end, endOffset, mapped)`, `Finish(finalOffset)` — algoritma htslib, sehingga bin dan chunk-nya sama dengan `samtools index` |
| `TabixConfig` | `VcfPreset`, `BedPreset`, `GffPreset`, `SamPreset`, `FromPreset(name)`, `TryParse(line, out name, out beg, out end)` (aturan `tbx_parse1` htslib: panjang REF, `END=`, `SVLEN`) |
| `TabixIndex` | `Load`/`Read` (`.tbi` dan `.csi` tabix), `Build(bgzfPath, config, format)`, `Save`, `Names`, `Config`, `Index` |
| `TabixReader` | `Open(path, indexPath?)`, `Query(name, beg, end)`, `Query("chr1:100-200")`, `Lines()`, `HeaderLines()` |
| `HtsRegion` | `Parse("chr1:1,000-2,000")` (berbasis 1 inklusif, atau satu referensi utuh) → `Start`/`End` berbasis 0 setengah terbuka |
| `Crc32` | `Compute`, `Append` |

## SAM, BAM, dan rekaman penjajaran (`.Mapping`)

| Tipe | Anggota utama |
|---|---|
| `AlignmentRecord` | `Name`, `Flags`, `ReferenceId`/`ReferenceName`, `Position` (berbasis 0), `MapQ`, `Cigar`, kolom pasangan (mate), `TemplateLength`, `Sequence`, `Qualities` (Phred, null untuk `*`), `AuxData`; pembantu flag (`IsPaired`, `IsReverse`, `IsSecondary`, ...); `End`, `ReferenceLength`, `QueryLength`, `Interval`, `Overlaps`, `AlignedPairs()`; tag: `Tags`, `GetTag`, `GetInt`, `GetString`, `HasTag`, `SetTag`, `RemoveTag`, `SetTags`; `Parse(samLine, header)`, `ToString()` (baris SAM), `ToSamRecord`, `FromSamRecord`, `ToSequenceRead`, `Clone` |
| `AuxField` | kolom opsional bertipe (`A c C s S i I f Z H B`): `Int`, `String`, `Char`, `Float`, `Hex`, `Array`, `AsInt64`, `AsDouble`, `AsString`, `AsArray<T>`, `ArrayType`, `Parse` |
| `SamHeader` | `Create(references, sortOrder)`, `Text`, `References`, `GetReferenceId`, `GetReferenceName`, `Lines`, `GetLines(type)`, `SortOrder`, `ReadGroupIds`, `WithSortOrder` |
| `IAlignmentReader` | `Header`, `IsIndexed`, `Records()`, `Query(name, start, end)`; ekstensi `Query("chr1:1-100")`, `Query(GenomicInterval)` |
| `BamReader` | `Open(path, indexPath?)` (mencari `.bai`/`.csi`), `Records()`, `AllRecords()`, `Read()`, `Query(...)`, `QueryUnmapped()`, `Index`, `Position`, `Seek`, `FirstRecordOffset` |
| `BamWriter` | `Create(path, header, options)`, `Write(record / records / ReadAlignment)`, `Index`; `BamWriterOptions` (`CompressionLevel`, `Index` = BAI/CSI dibangun sambil menulis, `CsiMinShift`, `CsiDepth`) |
| `BamIndexer` | `Build(stream, format)`, `Index(bamPath, format)` |
| `SamAlignmentReader`, `AlignmentFile.Open(path, referencePath?)` | pembaca teks SAM; membuka SAM/BAM/CRAM berdasarkan isi berkas |
| `AlignmentDepth`, `DepthOptions` | `Compute(reader atau records, name, start, end, options)` — semantik `samtools depth -a` (`ExcludeFlags`, `MinMapQ`, `MinBaseQuality`, `CountDeletions` = `-J`) |
| `CoordinateComparer` | `Instance`, `Sort(records)` (stabil, seperti `samtools sort`) |
| `SamWriter` | juga `SamWriter(writer, SamHeader)`, `Create(path, SamHeader)`, `Write(AlignmentRecord)` |

## CRAM (`.Cram`)

- `CramReader.Open(path, referenceFasta?, options?)` (memuat `.crai`), `Open(stream, referenceSequences)`, `Records()`, `Query(...)`, `BuildIndex()`, `Header`, `MajorVersion`, `MinorVersion`, `Index`.
- `CramReaderOptions`: `RegenerateMdNm` (bawaan true, seperti samtools), `ValidateReferenceMd5` (bawaan true), `ReadNamePrefix`.
- `CramIndex` (`Load`, `Read`, `Entries`, `Overlapping`), `CramIndexEntry`, `CramBlockMethod`.
- Didukung: CRAM 3.0 dan 3.1 (berkas 2.1 diterima tetapi belum diuji); metode blok raw, gzip, bzip2, rANS 4x8, rANS Nx16 (dengan pack, RLE, stripe, cat), adaptive arithmetic coder, fqzcomp, name tokeniser; encoding EXTERNAL, HUFFMAN, BYTE_ARRAY_LEN, BYTE_ARRAY_STOP, BETA, SUBEXP, GAMMA (ditambah GOLOMB dan GOLOMB_RICE yang sudah usang, diimplementasikan dari spesifikasi tetapi belum diuji); referensi eksternal maupun tertanam, kontainer multi-referensi dan multi-slice, nama dan kualitas lossy. Blok LZMA (metode 3) melempar `NotSupportedException`. Penulisan CRAM belum diimplementasikan.

## Anotasi (`.Annotation`)

- `GffReader.Read`, `BuildGenes(records)`, `LoadGenes(genome, path)`; `GffRecord`, `GffWriter`.
- `VariantEffectPredictor(genome, options)`: `Predict(variant)`, `Predict(chromosome, position, reference, alternate)`, `PredictTranscript`, `Annotate(vcfRecord)`, `AddCsqHeader(header)`.
- `VariantEffectOptions`: `UpstreamDistance`, `DownstreamDistance`, `Code`, `SpliceRegionExonBases`, `SpliceRegionIntronBases`.
- `VariantEffect`: `Allele`, `Consequences`, `MostSevere`, `Impact`, `GeneId`, `GeneName`, `TranscriptId`, `HgvsC`, `HgvsP`, `CdsPosition`, `ProteinPosition`, `AminoAcids`, `Codons`, `Distance`, `ToCsq()`.
- `Consequence` (istilah Sequence Ontology) dengan `ToSoTerm`, `SoAccession`, `Impact`, `ParseSoTerm`; `VariantImpact`.

## Pemetaan dan pemanggilan (`.Mapping`, `.Calling`)

- `ReadMapper(reference, options)`: `Map(read)`, `Map(reads)`, `IndexSize`; `MapperOptions` (`K`, `W`, `MaxOccurrences`, `Band`, `MaxUngappedMismatches`, `MaxCandidates`, `Match`, `Mismatch`, `GapOpen`, `GapExtend`, `MinScoreFraction`).
- `ReadAlignment`: `IsMapped`, `ReferenceName`, `Position`, `IsReverse`, `Cigar`, `Score`, `SecondBestScore`, `MapQ`, `EditDistance`, `Md`, `Interval`, `ToSam()`, `ToAlignmentRecord(header)`.
- `Cigar.Parse`, `Format`, `ReferenceLength`, `QueryLength`, `AlignedPairs`; `CigarElement`, `SamFlags`, `SamRecord.Parse`, `SamReader.Read`, `SamReader.ReadHeader`, `SamWriter.Create`, `SamWriter.Write`.
- `Pileup(reference, options)`: `Add(alignment)`, `Add(samRecord)`, `Add(alignmentRecord)`, `AddRange` (keluaran mapper atau rekaman BAM/CRAM), `Column`, `Depth`, `VariableColumns()`; `PileupOptions` (`MinBaseQuality`, `MinMapQ`, `DefaultQuality`), `PileupColumn`.
- `SnvCaller.Call(pileup, options)` dan `Call(reference, alignments, options)`, `CreateHeader`; `SnvCallerOptions` (`Heterozygosity`, `MinAltReads`, `MinAltFraction`, `MinQuality`, `SampleName`). Likelihood genotipe diploid dengan QUAL, GQ, dan PL.

## Perakitan dan pencarian gen

- `DeBruijnAssembler.Assemble(reads, options)` mengembalikan `AssemblyResult(Contigs, Stats, SolidKmers, TipsRemoved, BubblesPopped)`; `AssemblyOptions` (`K`, `MinKmerCount`, `RemoveTips`, `PopBubbles`, `CleaningRounds`, `MinContigLength`); `AssemblyStats` (N50/L50/N90/L90, GC).
- `DeBruijnGraph.Build(reads, k, minCount)`, `Unitigs()`, `RemoveTips`, `PopBubbles`, `Spectrum()`.
- `ProkaryoticGeneFinder(options)`: `Train`, `Predict(genome)`, `Predict(contigs)`, `PredictWithModels`, `IsTrained`; `GeneFinderOptions` (`Code`, `MinGeneLength`, `MarkovOrder`, `MaxOverlap`, `TrainingIterations`, `MinScore`, `RbsConsensus`, `RbsWeight`, `StartCodonLogPrior`); `PredictedGene` dengan `Extract`, `Translate`, `ToFeature`, `Interval`.

## Genetika populasi (`.Population`)

- `GenotypeCounts.From(variant, samples)`: `AlternateFrequency`, `ObservedHeterozygosity`, `ExpectedHeterozygosity`, `InbreedingCoefficient`, `HardyWeinbergP()`.
- `HardyWeinberg.ExactTest`, `ChiSquareTest`; `AlleleFrequencies.Of`, `Counts`, `ExpectedHeterozygosity`, `ObservedHeterozygosity`.
- `Diversity.Summarize(alignment)` mengembalikan `DiversitySummary` (`Pi`, `ThetaW`, `TajimasD`); juga `WattersonTheta`, `TajimasD`, `NucleotideDiversity`, `PairwiseDifferences`, `A1`, `A2`.
- `Fst.Hudson`, `Fst.WeirCockerham`, `VarianceComponents`.
- `LinkageDisequilibrium.FromHaplotypes`, `FromCounts`, `FromPhased`, `GenotypeRSquared`, yang mengembalikan `LdResult(D, DPrime, RSquared, Haplotypes)`.
- `SiteFrequencySpectrum.Unfolded`, `Folded`, `Fold`, `FromVariants`, `ThetaPi`, `ThetaW`, `TajimasD`.

---
[Indeks referensi API](index.md) · [Tutorial genomika](../tutorials/genomics.md)
