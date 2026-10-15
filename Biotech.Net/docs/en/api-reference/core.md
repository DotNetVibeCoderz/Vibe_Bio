# API reference: Core

Package `Gravicode.Biotech.Net.Core` · namespace `Biotech.Net.Core` · depends on Native.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Sequence

ASCII byte storage, case preserved, comparisons case-insensitive.

| Member | Purpose |
|---|---|
| `Sequence.Dna/Rna/Protein(residues, id, description)` | typed constructors |
| `Sequence.Parse(residues, id)` | detects the alphabet |
| `Sequence.FromBuffer(id, byte[], alphabet, …)` | no-copy construction for parsers |
| `Span`, `Memory`, `this[int]`, `Length`, `Type`, `Alphabet` | data access |
| `Subsequence(start, length)`, `Concat(other)`, `WithId(id)`, `ToUpper()` | derivations |
| `ReverseComplement()`, `Complement()`, `Transcribe()`, `BackTranscribe()` | kernel-backed |
| `Translate(code, frame, toStop)` | NCBI tables 1, 2, 4, 11 |
| `NucleotideCounts()`, `GcContent()`, `Composition()` | composition |
| `IndexOf(motif, start)`, `Contains(motif)` | search |
| `Features`, `Annotations`, `HasFeatures` | parsed metadata |
| `Preview(head)` | short display form |

## Alphabet, GeneticCode

- `Alphabet.Dna/DnaIupac/Rna/RnaIupac/Protein/ProteinExtended`, `Detect(residues)`, `IsValid`, `FirstInvalid`, `ForType`.
- `GeneticCode.Standard/VertebrateMitochondrial/MoldMitochondrial/Bacterial`, `All`, `FromId(id)`, `Translate(codon)`, `IsStart`, `IsStop`, `CodonsFor(aminoAcid)`, `CodonTable` (kernel order), `Index(codon)`, `CodonAt(index)`.

## Coordinates and features

| Type | Notes |
|---|---|
| `Strand` | `Forward`, `Reverse`, `Unknown` |
| `GenomicInterval(Chromosome, Start, End, Strand)` | **0-based half-open**; `Parse("chr1:1,001-2,000")`, `Overlaps`, `OverlapLength`, `Contains`, `Union`, `DistanceTo`, `Expand`, `Merge`, `CompareTo` |
| `FeatureLocation(Segments, Strand)` | 0-based half-open, join-aware; `Span(start, end)`, `Extract(parent)`, `IsJoin` |
| `SequenceFeature(Type, Location)` | `Qualifiers`, `this[qualifier]`, `With(qualifier, value)`, `Name` |

## Biological entities

- `Gene(id, location, name, biotype)` with `Transcripts`, `Attributes`.
- `Transcript(id, geneId, exons)` with `CodingRegion`, `Span`, `SplicedLength`, `Introns()`, `SplicedSequence(chromosome)`.
- `Chromosome(name, length, sequence?)`, `FromSequence`, `IsCircular`.
- `Genome(name, species, assembly)`: `Add(chromosome|gene)`, `this[name]`, `TryGetChromosome`, `GetSequence(region)`, `GenesIn(region)`, `FindGene`, `N50()`.
- `Protein(id, sequence, name)` / `Protein.FromString`: `MolecularWeight`, `IsoelectricPoint`, `Organism`, `GeneName`, `Function`, `StructureId`, `Annotations`.
- `SequenceRead(sequence, phredQualities)` / `FromAscii(id, residues, quality, offset)`: `Quality`, `MeanQuality()`, `ExpectedErrors()`, `QualityString(offset)`, `Slice`, `TrimByQuality(minQuality, window)`.

## Variants

- `Variant(chromosome, position, reference, alternate(s), id)` — **1-based** `Position`; `Type`, `IsTransition`, `IsTransversion`, `Quality`, `Filter`, `IsPass`, `Info`, `Genotypes`, `Interval` (converted to 0-based), `AlternateAlleleFrequency()`, `Classify(ref, alt)`.
- `Genotype(Alleles, Phased)`: `Parse("0/1")`, `Missing`, `IsMissing`, `IsHomozygousReference`, `IsHeterozygous`, `IsHomozygousAlternate`, `AlternateCount`.
- `VariantType`: SNV, insertion, deletion, MNV and so on.

## Dataset

A named numeric matrix with row/column names and annotations — the common currency for expression data.

`Dataset(name, rowNames, columnNames, values)`, `Zeros`, `FromArray`, `ReadDelimited`, `WriteTsv`, `this[row, column]` (by index or name), `Row`, `Column`, `RowIndex`, `ColumnIndex`, `ColumnsWhere(key, value)`, `Annotate(key, values)`, `SelectRows`, `SelectColumns`, `Transpose`, `Map`, `Log2`, `RowMeans`, `RowVariances`, `ColumnSums`, `RowAnnotations`, `ColumnAnnotations`.

## BiotechInfo

`Product`, `Tagline`, `Version`, `Backend`, `Credits`, `CreditsId` ("Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil") and `Banner`.

---
[API reference index](index.md) · [Sequences tutorial](../tutorials/sequences.md)
