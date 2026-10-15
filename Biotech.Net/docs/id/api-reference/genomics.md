# Referensi API: Genomics

Paket `Gravicode.Biotech.Net.Genomics` · namespace `Biotech.Net.Genomics.Vcf`, `.Annotation`, `.Mapping`, `.Calling`, `.Assembly`, `.GeneFinding`, `.Population`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## VCF (`.Vcf`)

| Tipe | Anggota utama |
|---|---|
| `VcfReader` | `Open(path atau stream)`, `Parse(text)`, `ReadRecords(path)`, `Header`, `Records()`, `Variants()`, `Query(region)` |
| `VcfWriter` | `Create(path, header)`, `Write(record / variant / records)`, `ToString(header, records)`, `Flush` |
| `VcfHeader` | `FileFormat`, `Info`, `Format`, `Filters`, `Contigs`, `Samples`, `AddInfo`, `AddFormat`, `AddFilter`, `AddContig`, `AddStandardGenotypeFields`, `ParseMetaLine`, `Write` |
| `VcfRecord` | `Variant`, `Format`, `SampleCount`, `GetGenotype(sample)`, `GetField`, `SetField`, `AlleleDepths`, `Depth`, `GenotypeQuality`, `PhredLikelihoods`, `SplitMultiallelic()`, `GenotypeIndex(j, k)` |
| `VcfFilters` (ekstensi) | `MinQuality`, `PassOnly`, `InRegion`, `InRegions`, `BiallelicSnvs` |

## Anotasi (`.Annotation`)

- `GffReader.Read`, `BuildGenes(records)`, `LoadGenes(genome, path)`; `GffRecord`, `GffWriter`.
- `VariantEffectPredictor(genome, options)`: `Predict(variant)`, `Predict(chromosome, position, reference, alternate)`, `PredictTranscript`, `Annotate(vcfRecord)`, `AddCsqHeader(header)`.
- `VariantEffectOptions`: `UpstreamDistance`, `DownstreamDistance`, `Code`, `SpliceRegionExonBases`, `SpliceRegionIntronBases`.
- `VariantEffect`: `Allele`, `Consequences`, `MostSevere`, `Impact`, `GeneId`, `GeneName`, `TranscriptId`, `HgvsC`, `HgvsP`, `CdsPosition`, `ProteinPosition`, `AminoAcids`, `Codons`, `Distance`, `ToCsq()`.
- `Consequence` (istilah Sequence Ontology) dengan `ToSoTerm`, `SoAccession`, `Impact`, `ParseSoTerm`; `VariantImpact`.

## Pemetaan dan pemanggilan (`.Mapping`, `.Calling`)

- `ReadMapper(reference, options)`: `Map(read)`, `Map(reads)`, `IndexSize`; `MapperOptions` (`K`, `W`, `MaxOccurrences`, `Band`, `MaxUngappedMismatches`, `MaxCandidates`, `Match`, `Mismatch`, `GapOpen`, `GapExtend`, `MinScoreFraction`).
- `ReadAlignment`: `IsMapped`, `ReferenceName`, `Position`, `IsReverse`, `Cigar`, `Score`, `SecondBestScore`, `MapQ`, `EditDistance`, `Md`, `Interval`, `ToSam()`.
- `Cigar.Parse`, `Format`, `ReferenceLength`, `QueryLength`; `CigarElement`, `SamFlags`, `SamRecord.Parse`, `SamReader.Read`, `SamReader.ReadHeader`, `SamWriter.Create`, `SamWriter.Write`.
- `Pileup(reference, options)`: `Add(alignment)`, `Add(samRecord)`, `AddRange`, `Column`, `Depth`, `VariableColumns()`; `PileupOptions` (`MinBaseQuality`, `MinMapQ`, `DefaultQuality`), `PileupColumn`.
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
