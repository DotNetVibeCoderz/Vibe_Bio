# API reference: Genomics

Package `Gravicode.Biotech.Net.Genomics` · namespaces `Biotech.Net.Genomics.Vcf`, `.IO`, `.Mapping`, `.Cram`, `.Annotation`, `.Calling`, `.Assembly`, `.GeneFinding`, `.Population`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## VCF (`.Vcf`)

| Type | Key members |
|---|---|
| `VcfReader` | `Open(path or stream)`, `Open(path, indexPath)`, `Parse(text)`, `ReadRecords(path)`, `Header`, `Records()`, `Variants()`, `Query(GenomicInterval)`, `Query("chr1:100-200")`, `IsIndexed` (a `.vcf.gz` with a `.tbi`/`.csi` beside it is queried through the index) |
| `VcfWriter` | `Create(path, header)` (a `.gz` path is written as BGZF, ready for tabix), `Write(record / variant / records)`, `ToString(header, records)`, `Flush` |
| `VcfHeader` | `FileFormat`, `Info`, `Format`, `Filters`, `Contigs`, `Samples`, `AddInfo`, `AddFormat`, `AddFilter`, `AddContig`, `AddStandardGenotypeFields`, `ParseMetaLine`, `Write` |
| `VcfRecord` | `Variant`, `Format`, `SampleCount`, `GetGenotype(sample)`, `GetField`, `SetField`, `AlleleDepths`, `Depth`, `GenotypeQuality`, `PhredLikelihoods`, `SplitMultiallelic()`, `GenotypeIndex(j, k)` |
| `VcfFilters` (extensions) | `MinQuality`, `PassOnly`, `InRegion`, `InRegions`, `BiallelicSnvs` |

## BGZF, binning indexes and tabix (`.IO`)

| Type | Key members |
|---|---|
| `BgzfReader` (a `Stream`) | `Open(path)`, `IsBgzf(path or stream)`, `VirtualPosition`, `Seek(VirtualOffset)`, `ReadLine(ByteBuffer)`, `AtEnd`, `VerifyChecksums`. Checks every block's CRC32 and ISIZE; keeps an LRU cache of decompressed blocks |
| `BgzfWriter` (a `Stream`) | `Create(path, level)`, `VirtualPosition`, `EnsureSpace(bytes)`, `FlushBlock()`, `EofMarker`. Blocks of ≤ 65 280 bytes, EOF marker on dispose |
| `VirtualOffset` | `BlockAddress`, `BlockOffset`, ordering operators |
| `BinningIndex` | `Load(path)` / `Read(stream)` for `.bai` and `.csi`, `Save`, `Write`, `Query(refId, beg, end)` → merged `Chunk`s, `UnplacedStart()`, `References` (`Bins`, `Linear`, `Span`, `MappedCount`, `UnmappedCount`), `UnplacedCount`, `MinShift`, `Depth`, `RegionToBin`, `RegionToBins` |
| `BinningIndexBuilder` | `Push(refId, beg, end, endOffset, mapped)`, `Finish(finalOffset)` — htslib's algorithm, so bins and chunks equal `samtools index` |
| `TabixConfig` | `VcfPreset`, `BedPreset`, `GffPreset`, `SamPreset`, `FromPreset(name)`, `TryParse(line, out name, out beg, out end)` (htslib `tbx_parse1` rules: REF length, `END=`, `SVLEN`) |
| `TabixIndex` | `Load`/`Read` (`.tbi` and tabix `.csi`), `Build(bgzfPath, config, format)`, `Save`, `Names`, `Config`, `Index` |
| `TabixReader` | `Open(path, indexPath?)`, `Query(name, beg, end)`, `Query("chr1:100-200")`, `Lines()`, `HeaderLines()` |
| `HtsRegion` | `Parse("chr1:1,000-2,000")` (1-based inclusive, or a whole reference) → 0-based half-open `Start`/`End` |
| `Crc32` | `Compute`, `Append` |

## SAM, BAM and alignment records (`.Mapping`)

| Type | Key members |
|---|---|
| `AlignmentRecord` | `Name`, `Flags`, `ReferenceId`/`ReferenceName`, `Position` (0-based), `MapQ`, `Cigar`, mate fields, `TemplateLength`, `Sequence`, `Qualities` (Phred, null for `*`), `AuxData`; flag helpers (`IsPaired`, `IsReverse`, `IsSecondary`, ...); `End`, `ReferenceLength`, `QueryLength`, `Interval`, `Overlaps`, `AlignedPairs()`; tags: `Tags`, `GetTag`, `GetInt`, `GetString`, `HasTag`, `SetTag`, `RemoveTag`, `SetTags`; `Parse(samLine, header)`, `ToString()` (SAM line), `ToSamRecord`, `FromSamRecord`, `ToSequenceRead`, `Clone` |
| `AuxField` | typed optional field (`A c C s S i I f Z H B`): `Int`, `String`, `Char`, `Float`, `Hex`, `Array`, `AsInt64`, `AsDouble`, `AsString`, `AsArray<T>`, `ArrayType`, `Parse` |
| `SamHeader` | `Create(references, sortOrder)`, `Text`, `References`, `GetReferenceId`, `GetReferenceName`, `Lines`, `GetLines(type)`, `SortOrder`, `ReadGroupIds`, `WithSortOrder` |
| `IAlignmentReader` | `Header`, `IsIndexed`, `Records()`, `Query(name, start, end)`; extensions `Query("chr1:1-100")`, `Query(GenomicInterval)` |
| `BamReader` | `Open(path, indexPath?)` (finds `.bai`/`.csi`), `Records()`, `AllRecords()`, `Read()`, `Query(...)`, `QueryUnmapped()`, `Index`, `Position`, `Seek`, `FirstRecordOffset` |
| `BamWriter` | `Create(path, header, options)`, `Write(record / records / ReadAlignment)`, `Index`; `BamWriterOptions` (`CompressionLevel`, `Index` = BAI/CSI built while writing, `CsiMinShift`, `CsiDepth`) |
| `BamIndexer` | `Build(stream, format)`, `Index(bamPath, format)` |
| `SamAlignmentReader`, `AlignmentFile.Open(path, referencePath?)` | SAM text reader; open SAM/BAM/CRAM by content |
| `AlignmentDepth`, `DepthOptions` | `Compute(reader or records, name, start, end, options)` — `samtools depth -a` semantics (`ExcludeFlags`, `MinMapQ`, `MinBaseQuality`, `CountDeletions` = `-J`) |
| `CoordinateComparer` | `Instance`, `Sort(records)` (stable, as `samtools sort`) |
| `SamWriter` | also `SamWriter(writer, SamHeader)`, `Create(path, SamHeader)`, `Write(AlignmentRecord)` |

## CRAM (`.Cram`)

- `CramReader.Open(path, referenceFasta?, options?)` (loads `.crai`), `Open(stream, referenceSequences)`, `Records()`, `Query(...)`, `BuildIndex()`, `Header`, `MajorVersion`, `MinorVersion`, `Index`.
- `CramReaderOptions`: `RegenerateMdNm` (default true, like samtools), `ValidateReferenceMd5` (default true), `ReadNamePrefix`.
- `CramIndex` (`Load`, `Read`, `Entries`, `Overlapping`), `CramIndexEntry`, `CramBlockMethod`.
- Supported: CRAM 3.0 and 3.1 (2.1 files are accepted but untested); block methods raw, gzip, bzip2, rANS 4x8, rANS Nx16 (with pack, RLE, stripe, cat), adaptive arithmetic coder, fqzcomp, name tokeniser; encodings EXTERNAL, HUFFMAN, BYTE_ARRAY_LEN, BYTE_ARRAY_STOP, BETA, SUBEXP, GAMMA (plus the deprecated GOLOMB and GOLOMB_RICE, implemented from the specification but untested); external and embedded references, multi-reference and multi-slice containers, lossy names and qualities. LZMA (method 3) blocks throw `NotSupportedException`. CRAM writing is not implemented.

## Annotation (`.Annotation`)

- `GffReader.Read`, `BuildGenes(records)`, `LoadGenes(genome, path)`; `GffRecord`, `GffWriter`.
- `VariantEffectPredictor(genome, options)`: `Predict(variant)`, `Predict(chromosome, position, reference, alternate)`, `PredictTranscript`, `Annotate(vcfRecord)`, `AddCsqHeader(header)`.
- `VariantEffectOptions`: `UpstreamDistance`, `DownstreamDistance`, `Code`, `SpliceRegionExonBases`, `SpliceRegionIntronBases`.
- `VariantEffect`: `Allele`, `Consequences`, `MostSevere`, `Impact`, `GeneId`, `GeneName`, `TranscriptId`, `HgvsC`, `HgvsP`, `CdsPosition`, `ProteinPosition`, `AminoAcids`, `Codons`, `Distance`, `ToCsq()`.
- `Consequence` (Sequence Ontology terms) with `ToSoTerm`, `SoAccession`, `Impact`, `ParseSoTerm`; `VariantImpact`.

## Mapping and calling (`.Mapping`, `.Calling`)

- `ReadMapper(reference, options)`: `Map(read)`, `Map(reads)`, `IndexSize`; `MapperOptions` (`K`, `W`, `MaxOccurrences`, `Band`, `MaxUngappedMismatches`, `MaxCandidates`, `Match`, `Mismatch`, `GapOpen`, `GapExtend`, `MinScoreFraction`).
- `ReadAlignment`: `IsMapped`, `ReferenceName`, `Position`, `IsReverse`, `Cigar`, `Score`, `SecondBestScore`, `MapQ`, `EditDistance`, `Md`, `Interval`, `ToSam()`, `ToAlignmentRecord(header)`.
- `Cigar.Parse`, `Format`, `ReferenceLength`, `QueryLength`, `AlignedPairs`; `CigarElement`, `SamFlags`, `SamRecord.Parse`, `SamReader.Read`, `SamReader.ReadHeader`, `SamWriter.Create`, `SamWriter.Write`.
- `Pileup(reference, options)`: `Add(alignment)`, `Add(samRecord)`, `Add(alignmentRecord)`, `AddRange` (mapper output or BAM/CRAM records), `Column`, `Depth`, `VariableColumns()`; `PileupOptions` (`MinBaseQuality`, `MinMapQ`, `DefaultQuality`), `PileupColumn`.
- `SnvCaller.Call(pileup, options)` and `Call(reference, alignments, options)`, `CreateHeader`; `SnvCallerOptions` (`Heterozygosity`, `MinAltReads`, `MinAltFraction`, `MinQuality`, `SampleName`). Diploid genotype likelihoods with QUAL, GQ and PL.

## Assembly and gene finding

- `DeBruijnAssembler.Assemble(reads, options)` returns `AssemblyResult(Contigs, Stats, SolidKmers, TipsRemoved, BubblesPopped)`; `AssemblyOptions` (`K`, `MinKmerCount`, `RemoveTips`, `PopBubbles`, `CleaningRounds`, `MinContigLength`); `AssemblyStats` (N50/L50/N90/L90, GC).
- `DeBruijnGraph.Build(reads, k, minCount)`, `Unitigs()`, `RemoveTips`, `PopBubbles`, `Spectrum()`.
- `ProkaryoticGeneFinder(options)`: `Train`, `Predict(genome)`, `Predict(contigs)`, `PredictWithModels`, `IsTrained`; `GeneFinderOptions` (`Code`, `MinGeneLength`, `MarkovOrder`, `MaxOverlap`, `TrainingIterations`, `MinScore`, `RbsConsensus`, `RbsWeight`, `StartCodonLogPrior`); `PredictedGene` with `Extract`, `Translate`, `ToFeature`, `Interval`.

## Population genetics (`.Population`)

- `GenotypeCounts.From(variant, samples)`: `AlternateFrequency`, `ObservedHeterozygosity`, `ExpectedHeterozygosity`, `InbreedingCoefficient`, `HardyWeinbergP()`.
- `HardyWeinberg.ExactTest`, `ChiSquareTest`; `AlleleFrequencies.Of`, `Counts`, `ExpectedHeterozygosity`, `ObservedHeterozygosity`.
- `Diversity.Summarize(alignment)` returns `DiversitySummary` (`Pi`, `ThetaW`, `TajimasD`); also `WattersonTheta`, `TajimasD`, `NucleotideDiversity`, `PairwiseDifferences`, `A1`, `A2`.
- `Fst.Hudson`, `Fst.WeirCockerham`, `VarianceComponents`.
- `LinkageDisequilibrium.FromHaplotypes`, `FromCounts`, `FromPhased`, `GenotypeRSquared`, returning `LdResult(D, DPrime, RSquared, Haplotypes)`.
- `SiteFrequencySpectrum.Unfolded`, `Folded`, `Fold`, `FromVariants`, `ThetaPi`, `ThetaW`, `TajimasD`.

---
[API reference index](index.md) · [Genomics tutorial](../tutorials/genomics.md)
