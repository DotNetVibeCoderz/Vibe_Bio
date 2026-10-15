# API reference: Genomics

Package `Gravicode.Biotech.Net.Genomics` · namespaces `Biotech.Net.Genomics.Vcf`, `.Annotation`, `.Mapping`, `.Calling`, `.Assembly`, `.GeneFinding`, `.Population`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## VCF (`.Vcf`)

| Type | Key members |
|---|---|
| `VcfReader` | `Open(path or stream)`, `Parse(text)`, `ReadRecords(path)`, `Header`, `Records()`, `Variants()`, `Query(region)` |
| `VcfWriter` | `Create(path, header)`, `Write(record / variant / records)`, `ToString(header, records)`, `Flush` |
| `VcfHeader` | `FileFormat`, `Info`, `Format`, `Filters`, `Contigs`, `Samples`, `AddInfo`, `AddFormat`, `AddFilter`, `AddContig`, `AddStandardGenotypeFields`, `ParseMetaLine`, `Write` |
| `VcfRecord` | `Variant`, `Format`, `SampleCount`, `GetGenotype(sample)`, `GetField`, `SetField`, `AlleleDepths`, `Depth`, `GenotypeQuality`, `PhredLikelihoods`, `SplitMultiallelic()`, `GenotypeIndex(j, k)` |
| `VcfFilters` (extensions) | `MinQuality`, `PassOnly`, `InRegion`, `InRegions`, `BiallelicSnvs` |

## Annotation (`.Annotation`)

- `GffReader.Read`, `BuildGenes(records)`, `LoadGenes(genome, path)`; `GffRecord`, `GffWriter`.
- `VariantEffectPredictor(genome, options)`: `Predict(variant)`, `Predict(chromosome, position, reference, alternate)`, `PredictTranscript`, `Annotate(vcfRecord)`, `AddCsqHeader(header)`.
- `VariantEffectOptions`: `UpstreamDistance`, `DownstreamDistance`, `Code`, `SpliceRegionExonBases`, `SpliceRegionIntronBases`.
- `VariantEffect`: `Allele`, `Consequences`, `MostSevere`, `Impact`, `GeneId`, `GeneName`, `TranscriptId`, `HgvsC`, `HgvsP`, `CdsPosition`, `ProteinPosition`, `AminoAcids`, `Codons`, `Distance`, `ToCsq()`.
- `Consequence` (Sequence Ontology terms) with `ToSoTerm`, `SoAccession`, `Impact`, `ParseSoTerm`; `VariantImpact`.

## Mapping and calling (`.Mapping`, `.Calling`)

- `ReadMapper(reference, options)`: `Map(read)`, `Map(reads)`, `IndexSize`; `MapperOptions` (`K`, `W`, `MaxOccurrences`, `Band`, `MaxUngappedMismatches`, `MaxCandidates`, `Match`, `Mismatch`, `GapOpen`, `GapExtend`, `MinScoreFraction`).
- `ReadAlignment`: `IsMapped`, `ReferenceName`, `Position`, `IsReverse`, `Cigar`, `Score`, `SecondBestScore`, `MapQ`, `EditDistance`, `Md`, `Interval`, `ToSam()`.
- `Cigar.Parse`, `Format`, `ReferenceLength`, `QueryLength`; `CigarElement`, `SamFlags`, `SamRecord.Parse`, `SamReader.Read`, `SamReader.ReadHeader`, `SamWriter.Create`, `SamWriter.Write`.
- `Pileup(reference, options)`: `Add(alignment)`, `Add(samRecord)`, `AddRange`, `Column`, `Depth`, `VariableColumns()`; `PileupOptions` (`MinBaseQuality`, `MinMapQ`, `DefaultQuality`), `PileupColumn`.
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
