# API reference: Alignments

Package `Gravicode.Biotech.Net.Alignments` · namespaces `Biotech.Net.Alignments`, `.Phylogenetics`, `.Search`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Scoring

`SubstitutionMatrix` — `Blosum45/50/62/80/90`, `Pam30/70/250`, `Nuc44`, `AvailableNames`, `Load(name)`, `MatchMismatch(match, mismatch, symbols)`, `Parse(text)`, `FromFile`, `Read`, `Score(a, b)`, `KernelMatrix`, `IsNucleotide`, `IsSymmetric`, `ToNcbiFormat()`.

## Pairwise

| Type | Key members |
|---|---|
| `PairwiseAligner` | `Mode` (`Biotech.Net.Native.AlignmentMode`), `SubstitutionMatrix`, `MatchScore`, `MismatchScore`, `OpenGapScore`, `ExtendGapScore`, `GapScore` (sets both), `Band`, `Align(target, query)`, `Score(target, query)` |
| `PairwiseAlignment` | `Score`, `Operations`, `AlignedTarget`, `AlignedQuery`, `TargetStart/End`, `QueryStart/End`, `Identities`, `Mismatches`, `Positives`, `Gaps`, `GapOpenings`, `Identity`, `Coordinates`, `AlignedBlocks`, `Cigar(extended, softClip)`, `Format(width)` |

Gap scores are zero or negative (Biopython convention); the kernel underneath takes non-negative costs.

## Multiple alignment

- `ProgressiveAligner`: `SubstitutionMatrix`, `OpenGapScore`, `ExtendGapScore`, `Distance` (`GuideDistance.Kmer|Alignment`), `TreeMethod` (`GuideTreeMethod.Upgma|NeighborJoining`), `KmerSize`, `Align(sequences)`, `ComputeDistances(sequences)`.
- `MultipleSequenceAlignment`: `Ids`, `Rows`, `Count`, `Length`, `Column(i)`, `Consensus(threshold)`, `ColumnConservation()`, `ColumnEntropy()`, `PercentIdentity(i, j)`, `SumOfPairsScore(matrix)`, `ToFasta`, `FromFasta`, `ToClustal`, `ConservationLine`, `ToSequences`, `GuideTree`.

## Phylogenetics (`.Phylogenetics`)

- `DistanceMatrix(names[, distances])`, `FromAlignment(msa, model)`, `PairDistance`, `ToPhylip()`; `DistanceModel.PDistance|JukesCantor|Kimura2P`.
- `TreeBuilder.Upgma(distances)`, `TreeBuilder.NeighborJoining(distances, clampNegative)`.
- `PhylogeneticTree`: `Root`, `Leaves`, `LeafNames`, `TotalBranchLength`, `FindLeaf`, `Distance(leafA, leafB)`, `ParseNewick`, `ToNewick(branchLengths, decimals)`, `RootAtMidpoint()`, `RobinsonFoulds(other, normalized)`.
- `TreeNode`: `Name`, `BranchLength`, `Parent`, `Children`, `IsLeaf`, `AddChild`, `RemoveChild`, `PreOrder`, `PostOrder`, `Leaves`.

## Database search (`.Search`)

- `SequenceDatabase(sequences, wordSize?)`: `Count`, `TotalLength`, `WordSize`, `IsNucleotide`, `IndexedWords`, `Candidates(query, minHits)`, `Search(query, options)`.
- `SearchOptions`: `SubstitutionMatrix`, `GapOpen`, `GapExtend`, `MaxEValue`, `MaxHits`, `UseSeeds`, `MinSeedHits`, `Statistics`, `IncludeAlignments`.
- `SearchHit`: `Subject`, `SubjectIndex`, `Score`, `BitScore`, `EValue`, `Alignment`.
- `KarlinAltschulParameters`: `Blosum62Gapped11_1`, `Blosum62Ungapped`, `Blastn2_3Gapped5_2`, `Lookup(matrixName, gapOpen, gapExtend)`, `BitScore`, `EValue`. Without parameters, bit score and E-value are `NaN`.

---
[API reference index](index.md) · [Alignments tutorial](../tutorials/alignments.md)
