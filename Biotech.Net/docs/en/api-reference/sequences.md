# API reference: Sequences

Package `Gravicode.Biotech.Net.Sequences` · namespaces `Biotech.Net.Sequences`, `Biotech.Net.Sequences.IO`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## File I/O (`Biotech.Net.Sequences.IO`)

| Type | Key members |
|---|---|
| `FastaReader` / `FastaWriter` | `Read(path|TextReader)`, `Parse(text)`, `ReadAsync(stream)`; `Write(...)`, `ToString(sequences, lineWidth)` |
| `FastqReader` / `FastqWriter` | `Read`, `Parse`, `DetectEncoding`; `FastqEncoding.Auto/Phred33/Phred64` |
| `GenBankReader` / `GenBankWriter` | `Read`, `Parse`, `Write`, `ToString`; header keys in `RecordKeys` |
| `EmblReader` / `EmblWriter` | same shape as GenBank |
| `PdbSequenceReader` | `ReadSeqres`, `ReadAtomSequences`, `ResidueCodes`, `ToOneLetter` |
| `SequenceFile` | `DetectFormat`, `Read`, `ReadFirst`, `ReadReads`, `Write` — format-agnostic entry point |
| `BioFile` | `OpenRead`, `OpenText`, `OpenWrite`, `CreateText`, transparent gzip |
| `FastaIndex`, `FastaIndexEntry`, `IndexedFastaReader` | samtools-style `.fai`: `Build`, `Save`, `Load`, `Open(path)`, `Fetch(name, start, end)`, `Fetch(region)`, `FetchAll`, `ToGenome` |
| `InsdcLocationParser` | GenBank/EMBL location strings: `Parse`, `Format`, partial/remote details |

## Analysis (`Biotech.Net.Sequences`)

| Area | Types and members |
|---|---|
| ORFs | `OrfFinder.Find/Longest`, `OrfOptions` (`MinLength`, `Code`, `StartMode`, `IncludePartial`, `Nested`, `Forward`, `Reverse`), `Orf.Extract/Translate/ToFeature`, `SixFrame.Translate` |
| Motifs | `MotifSearch.FindIupac/FindRegex/FindProsite/Count`, `MotifMatch`, `Iupac.Mask/Matches/ToRegex/ReverseComplement`, `PrositeToRegex` |
| Restriction | `RestrictionEnzymes.All/Get/TryGet/SingleCutters/NonCutters`, `RestrictionDigest.FindCuts/Digest`, `DigestResult.FragmentSizes/FragmentSequences`, `RestrictionEnzyme(Name, Site, CutTop, CutBottom)` |
| k-mers | `KmerProfile.Count`, `Top`, `Frequencies`, `Spectrum`, `EuclideanDistance`, `CosineSimilarity`; `KmerCounter` (hash-based, canonical option); `Kmer.TryEncode/Decode/ReverseComplement/Canonical/ForEach` |
| Composition | `SlidingWindow.GcContent/GcSkew/AtSkew/CumulativeGcSkew/Entropy/Profile` → `WindowValue`; `Complexity.DustScore/DustProfile/LowComplexityRegions/Mask/LinguisticComplexity` |
| Statistics | `SequenceStatistics.Nx/N50/Summarize/LengthHistogram`, `SequenceSetSummary` |
| Random data | `RandomSequence.Dna/Protein/FromComposition/Shuffle/ShuffleBlocks` |
| Thermodynamics | `MeltingTemperature.Wallace/BasicGc/SaltAdjustedGc/NearestNeighbor/FreeEnergy37`, `MeltingConditions`, `NucleicAcidMass.MolecularWeight` |
| Codon usage | `CodonUsageTable.FromSequences/Add/PerThousand/Rscu/RelativeAdaptiveness/Cai/Optimize/ToTable` |

See the [sequences tutorial](../tutorials/sequences.md) for worked examples.

---
[API reference index](index.md)
