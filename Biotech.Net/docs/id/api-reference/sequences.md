# Referensi API: Sequences

Paket `Gravicode.Biotech.Net.Sequences` · namespace `Biotech.Net.Sequences`, `Biotech.Net.Sequences.IO`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## I/O berkas (`Biotech.Net.Sequences.IO`)

| Tipe | Anggota utama |
|---|---|
| `FastaReader` / `FastaWriter` | `Read(path atau TextReader)`, `Parse(text)`, `ReadAsync(stream)`; `Write(...)`, `ToString(sequences, lineWidth)` |
| `FastqReader` / `FastqWriter` | `Read`, `Parse`, `DetectEncoding`; `FastqEncoding.Auto/Phred33/Phred64` |
| `GenBankReader` / `GenBankWriter` | `Read`, `Parse`, `Write`, `ToString`; kunci header di `RecordKeys` |
| `EmblReader` / `EmblWriter` | berbentuk sama dengan GenBank |
| `PdbSequenceReader` | `ReadSeqres`, `ReadAtomSequences`, `ResidueCodes`, `ToOneLetter` |
| `SequenceFile` | `DetectFormat`, `Read`, `ReadFirst`, `ReadReads`, `Write` — titik masuk tanpa peduli format |
| `BioFile` | `OpenRead`, `OpenText`, `OpenWrite`, `CreateText`, gzip transparan |
| `FastaIndex`, `FastaIndexEntry`, `IndexedFastaReader` | `.fai` gaya samtools: `Build`, `Save`, `Load`, `Open(path)`, `Fetch(name, start, end)`, `Fetch(region)`, `FetchAll`, `ToGenome` |
| `InsdcLocationParser` | string lokasi GenBank/EMBL: `Parse`, `Format`, detail parsial/remote |

## Analisis (`Biotech.Net.Sequences`)

| Bidang | Tipe dan anggota |
|---|---|
| ORF | `OrfFinder.Find/Longest`, `OrfOptions` (`MinLength`, `Code`, `StartMode`, `IncludePartial`, `Nested`, `Forward`, `Reverse`), `Orf.Extract/Translate/ToFeature`, `SixFrame.Translate` |
| Motif | `MotifSearch.FindIupac/FindRegex/FindProsite/Count`, `MotifMatch`, `Iupac.Mask/Matches/ToRegex/ReverseComplement`, `PrositeToRegex` |
| Restriksi | `RestrictionEnzymes.All/Get/TryGet/SingleCutters/NonCutters`, `RestrictionDigest.FindCuts/Digest`, `DigestResult.FragmentSizes/FragmentSequences`, `RestrictionEnzyme(Name, Site, CutTop, CutBottom)` |
| k-mer | `KmerProfile.Count`, `Top`, `Frequencies`, `Spectrum`, `EuclideanDistance`, `CosineSimilarity`; `KmerCounter` (berbasis hash, opsi kanonik); `Kmer.TryEncode/Decode/ReverseComplement/Canonical/ForEach` |
| Komposisi | `SlidingWindow.GcContent/GcSkew/AtSkew/CumulativeGcSkew/Entropy/Profile` → `WindowValue`; `Complexity.DustScore/DustProfile/LowComplexityRegions/Mask/LinguisticComplexity` |
| Statistik | `SequenceStatistics.Nx/N50/Summarize/LengthHistogram`, `SequenceSetSummary` |
| Data acak | `RandomSequence.Dna/Protein/FromComposition/Shuffle/ShuffleBlocks` |
| Termodinamika | `MeltingTemperature.Wallace/BasicGc/SaltAdjustedGc/NearestNeighbor/FreeEnergy37`, `MeltingConditions`, `NucleicAcidMass.MolecularWeight` |
| Penggunaan kodon | `CodonUsageTable.FromSequences/Add/PerThousand/Rscu/RelativeAdaptiveness/Cai/Optimize/ToTable` |

Contoh penggunaannya ada di [tutorial sekuens](../tutorials/sequences.md).

---
[Indeks referensi API](index.md)
