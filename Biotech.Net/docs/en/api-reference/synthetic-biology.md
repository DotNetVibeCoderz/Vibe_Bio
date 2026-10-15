# API reference: Synthetic biology

Package `Gravicode.Biotech.Net.SyntheticBiology` · namespaces `Biotech.Net.SyntheticBiology`, `.Codons`, `.Primers`, `.Cloning`, `.Circuits`, `.Crispr`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Codon optimisation (`.Codons`)

- `CodonUsageTable`: `EColiK12`, `SaccharomycesCerevisiae`, `HomoSapiens`, `CricetulusGriseus`, `PichiaPastoris`, `All`, `FromName`, `FromCounts`, `FromCodingSequences`; `PerThousand`, `Frequency`, `RelativeAdaptiveness`, `SynonymousCodons`, `MostFrequentCodon`, `RarestCodon`, `Cai`, `Gc3Content`, `CodonAdaptationDistance`.
- `CodonOptimizer.Optimize(protein, options)`, `OptimizeCoding(cds, options)`, `BackTranslate`, `CodonFrequencyProfile`; `CodonOptimizationOptions` (`Target`, `Source`, `Method`, `Seed`, `Constraints`, `Code`, `PreserveStop`, `MaxRepairIterations`), `CodonOptimizationMethod`, `CodonOptimizationResult(Dna, Cai, GcContent, CodonsChanged, RemainingViolations, RepairIterations)`.
- `SequenceConstraints`: `Default`, `ForGeneSynthesis`, `ForbiddenSites`, `MaxHomopolymer`, `MinGc`, `MaxGc`, `GcWindow`, `AvoidHairpins`, `AvoidRepeats`, `MaxDirectRepeat`, `Check(dna)`, `IsSatisfiedBy`; `ConstraintViolation(Kind, Start, Length, Detail)`.

## Primers (`.Primers`)

- `PrimerDesigner.DesignForward`, `DesignReverse`, `DesignPairs(template, options)`, `Evaluate(oligo, …)`.
- `PrimerDesignOptions`: lengths, `MinTm`/`OptTm`/`MaxTm`, `MaxTmDifference`, GC window, `MaxPolyX`, GC clamp, dimer and hairpin ΔG limits, product size, `SaltConditions`, result caps and the penalty weights; `TemplateRegion`.
- `Primer(Sequence, Start, Length, Strand, Tm, GcFraction, SelfDimerDeltaG, ThreePrimeSelfDimerDeltaG, HairpinDeltaG, EndStability, GcClamp, LongestHomopolymer, Penalty)`; `PrimerPair(Forward, Reverse, ProductSize, ProductSequence, TmDifference, HeteroDimerDeltaG, ThreePrimeHeteroDimerDeltaG, Penalty)`.
- `NearestNeighbourTm.Tm`, `Duplex`, `TmFromThermodynamics`, `WallaceTm`, `GcContentTm`, `IsSelfComplementary`; `SaltConditions.Default/StandardPcr/Qpcr`; `DuplexThermodynamics`.
- `OligoThermodynamics.SelfDimer`, `HeteroDimer`, `ThreePrimeDimer`, `Hairpin`, `HairpinLoopDeltaG37`; `OligoDuplex`, `OligoHairpin`.

## Cloning (`.Cloning`)

- `DnaFragment(Id, Sequence, FivePrimeOverhang, ThreePrimeOverhang, IsCircular)`: `Duplex`, `CanLigateTo`, `Ligate`, `Circularise`, `ReverseComplement`, `ToSequence`, `Features`.
- `GibsonAssembly.DesignOverlaps(fragments, options)`, `Simulate(fragments or sequences, options)`, `BestOverlap`; `GibsonOptions` (overlap length and Tm targets, GC window, annealing lengths, `Circular`), `GibsonJunction`, `GibsonDesign`, `GibsonAssemblyResult`.
- `GoldenGate.Digest`, `IsDomesticated`, `AddFlanks`, `DesignOverhangs`, `AreOrthogonal`, `OverhangFidelity`, `Simulate(parts, enzyme, backbone, backboneIsCircular)`; `GoldenGateResult`.
- `RestrictionEnzyme` (Type IIS aware) with `All`, `FromName`, `TryFromName`, `Find(dna, circular)`, `CountSites`, `OverhangLength`, `OverhangKind`, `IsPalindromic`; `RestrictionSite`.
- `Plasmid(name, sequence)`: `Slice`, `Rotate`, `WithFeature`, `Digest(enzyme)`, `Fragments(enzymes)`, `UniqueCutters`, `ToGenBankSummary`, `ToSequence`, `DemonstrationBackbone`.

## Circuits (`.Circuits`)

- `BioPart(Name, Type, Sequence, Description)` with `RelativeStrength`, `RegulatedBy`, `IsRepressible`, `IsInducible`, `Inducer`, `RegistryId`, `ToFeature`; `PartType`; `TranscriptionalUnit(Name, Parts)` with `Assemble`, `Features`, `ToSequence`.
- `PartsRegistry.All`, `Get`, `TryGet`, `OfType`, `Search`, `BuildUnit(name, partNames)` — promoters (`J23100`, `J23106`, `J23119`, `R0010`, `R0011`, `R0040`, `R0051`), RBSs (`B0030`–`B0034`), coding sequences (`GFP`, `mCherry`, `LacI`, `TetR`, `cI`), terminators (`B0010`, `B0012`, `B0015`) and operators.
- `RbsStrengthModel.Predict(utr, cds)`, `RelativeStrength`, `SpacingEnergy`, `StartCodonEnergy`, constants `AntiShineDalgarno`, `ConsensusUtr`; `RbsPrediction`.
- `CircuitModels.ToggleSwitch`, `SimulateToggleSwitch`, `Repressilator`, `SimulateRepressilator`, `EstimateOscillationPeriod`, `InducibleExpression`, `SimulateInduction`, `SteadyStateProtein`, `SteadyStateMrna`; `HillParameters`.
- `OdeIntegrator.Rk4`, `Rk45`; `OdeSystem` delegate, `OdeSolution(Times, States)` with `SeriesOf`, `InterpolateAt`, `Final`, `RejectedSteps`; `IntegratorKind`.

## CRISPR analysis (`.Crispr`)

| Type | Purpose |
|---|---|
| `CasNuclease` | `SpCas9`, `SpCas9Nag`, `SpCas9Vqr`, `SpCas9Vrer`, `XCas9`, `SpG`, `SpRy`, `SaCas9`, `SaCas9Kkh`, `Cas12a`, `LbCas12a`, `NmeCas9`, `All`, `FromName`; `Pam`, `PamSide`, `GuideLength`, `CutSite`, `MatchesPam` |
| `GuideFinder`, `GuideSearchOptions`, `Guide` | enumerate protospacers; GC, homopolymer, poly-T and leading-G filters; `CutSite`, `Interval`, `FullSite`, `U6CloningOligo` |
| `OnTargetScorer`, `GuidePenalties` | Doench 2014 Rule Set 1 (`RuleSet1`), `Score`, `Penalties`, `CombinedScore` |
| `CrisprSiteIndex` | `Build(sequences, nuclease, includeNagPam)`, `Add`, `Query(spacer, maxMismatches)`, `QueryExhaustive`, `ProtospacerAt` |
| `OffTargetSearch`, `OffTargetSearchOptions`, `OffTarget`, `OffTargetProfile` | mismatch (and optional bulge) search; `MitSpecificity`, `CfdSpecificity`, `MismatchHistogram` |
| `CfdScorer`, `MitSpecificityScorer` | CFD mismatch and PAM weights; MIT hit and guide scores |
| `GuideDesigner`, `GuideDesignOptions`, `GuideDesign` | ranked designs combining on-target score, penalties and specificity |
| `BaseEditor`, `BaseEditing` | published editor windows (`Be3`, `Be4Max`, `Ye1`, `Abe710`, `Abe8E`); `PredictEdits`, `PredictProteinConsequence`, `FindEditingGuides`, `EditingWindowInterval`; `BaseEditOutcome`, `BaseEditSite`, `ProteinConsequence`, `EditConsequence` |
| `PrimeEditor`, `PegRna`, `PrimeEdit`, `NickingGuide`, `PrimeEditorOptions` | pegRNA component layout (spacer, scaffold, PBS, RTT), `Design`, `DesignNickingGuides` |

`DnaUtility` holds the shared helpers: `ReverseComplement`, `GcFraction`, `LongestHomopolymer`, `MatchesIupac`, `TryPack`, `Unpack`, `MismatchCount`, `NormalizeDna`, `IsUnambiguous`.

---
[API reference index](index.md) · [Synthetic biology tutorial](../tutorials/synthetic-biology.md)
