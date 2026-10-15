# Referensi API: Core

Paket `Gravicode.Biotech.Net.Core` · namespace `Biotech.Net.Core` · bergantung pada Native.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Sequence

Penyimpanan byte ASCII, besar-kecil huruf dipertahankan, perbandingan mengabaikannya.

| Anggota | Kegunaan |
|---|---|
| `Sequence.Dna/Rna/Protein(residues, id, description)` | konstruktor bertipe |
| `Sequence.Parse(residues, id)` | mendeteksi alfabet |
| `Sequence.FromBuffer(id, byte[], alphabet, …)` | konstruksi tanpa salin untuk parser |
| `Span`, `Memory`, `this[int]`, `Length`, `Type`, `Alphabet` | akses data |
| `Subsequence(start, length)`, `Concat(other)`, `WithId(id)`, `ToUpper()` | turunan |
| `ReverseComplement()`, `Complement()`, `Transcribe()`, `BackTranscribe()` | memakai kernel |
| `Translate(code, frame, toStop)` | tabel NCBI 1, 2, 4, 11 |
| `NucleotideCounts()`, `GcContent()`, `Composition()` | komposisi |
| `IndexOf(motif, start)`, `Contains(motif)` | pencarian |
| `Features`, `Annotations`, `HasFeatures` | metadata hasil parsing |
| `Preview(head)` | bentuk tampilan singkat |

## Alphabet, GeneticCode

- `Alphabet.Dna/DnaIupac/Rna/RnaIupac/Protein/ProteinExtended`, `Detect(residues)`, `IsValid`, `FirstInvalid`, `ForType`.
- `GeneticCode.Standard/VertebrateMitochondrial/MoldMitochondrial/Bacterial`, `All`, `FromId(id)`, `Translate(codon)`, `IsStart`, `IsStop`, `CodonsFor(aminoAcid)`, `CodonTable` (urutan kernel), `Index(codon)`, `CodonAt(index)`.

## Koordinat dan fitur

| Tipe | Catatan |
|---|---|
| `Strand` | `Forward`, `Reverse`, `Unknown` |
| `GenomicInterval(Chromosome, Start, End, Strand)` | **berbasis 0, setengah terbuka**; `Parse("chr1:1,001-2,000")`, `Overlaps`, `OverlapLength`, `Contains`, `Union`, `DistanceTo`, `Expand`, `Merge`, `CompareTo` |
| `FeatureLocation(Segments, Strand)` | berbasis 0 setengah terbuka, mendukung join; `Span(start, end)`, `Extract(parent)`, `IsJoin` |
| `SequenceFeature(Type, Location)` | `Qualifiers`, `this[qualifier]`, `With(qualifier, value)`, `Name` |

## Entitas biologis

- `Gene(id, location, name, biotype)` dengan `Transcripts`, `Attributes`.
- `Transcript(id, geneId, exons)` dengan `CodingRegion`, `Span`, `SplicedLength`, `Introns()`, `SplicedSequence(chromosome)`.
- `Chromosome(name, length, sequence?)`, `FromSequence`, `IsCircular`.
- `Genome(name, species, assembly)`: `Add(chromosome|gene)`, `this[name]`, `TryGetChromosome`, `GetSequence(region)`, `GenesIn(region)`, `FindGene`, `N50()`.
- `Protein(id, sequence, name)` / `Protein.FromString`: `MolecularWeight`, `IsoelectricPoint`, `Organism`, `GeneName`, `Function`, `StructureId`, `Annotations`.
- `SequenceRead(sequence, phredQualities)` / `FromAscii(id, residues, quality, offset)`: `Quality`, `MeanQuality()`, `ExpectedErrors()`, `QualityString(offset)`, `Slice`, `TrimByQuality(minQuality, window)`.

## Varian

- `Variant(chromosome, position, reference, alternate(s), id)` — `Position` **berbasis 1**; `Type`, `IsTransition`, `IsTransversion`, `Quality`, `Filter`, `IsPass`, `Info`, `Genotypes`, `Interval` (dikonversi ke basis 0), `AlternateAlleleFrequency()`, `Classify(ref, alt)`.
- `Genotype(Alleles, Phased)`: `Parse("0/1")`, `Missing`, `IsMissing`, `IsHomozygousReference`, `IsHeterozygous`, `IsHomozygousAlternate`, `AlternateCount`.
- `VariantType`: SNV, insersi, delesi, MNV, dan seterusnya.

## Dataset

Matriks numerik bernama dengan nama baris/kolom dan anotasi — mata uang bersama untuk data ekspresi.

`Dataset(name, rowNames, columnNames, values)`, `Zeros`, `FromArray`, `ReadDelimited`, `WriteTsv`, `this[row, column]` (indeks atau nama), `Row`, `Column`, `RowIndex`, `ColumnIndex`, `ColumnsWhere(key, value)`, `Annotate(key, values)`, `SelectRows`, `SelectColumns`, `Transpose`, `Map`, `Log2`, `RowMeans`, `RowVariances`, `ColumnSums`, `RowAnnotations`, `ColumnAnnotations`.

## BiotechInfo

`Product`, `Tagline`, `Version`, `Backend`, `Credits`, `CreditsId` ("Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil"), dan `Banner`.

---
[Indeks referensi API](index.md) · [Tutorial sekuens](../tutorials/sequences.md)
