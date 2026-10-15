# Referensi API: SyntheticBiology

Paket `Gravicode.Biotech.Net.SyntheticBiology` · namespace `Biotech.Net.SyntheticBiology`, `.Codons`, `.Primers`, `.Cloning`, `.Circuits`, `.Crispr`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Optimasi kodon (`.Codons`)

- `CodonUsageTable`: `EColiK12`, `SaccharomycesCerevisiae`, `HomoSapiens`, `CricetulusGriseus`, `PichiaPastoris`, `All`, `FromName`, `FromCounts`, `FromCodingSequences`; `PerThousand`, `Frequency`, `RelativeAdaptiveness`, `SynonymousCodons`, `MostFrequentCodon`, `RarestCodon`, `Cai`, `Gc3Content`, `CodonAdaptationDistance`.
- `CodonOptimizer.Optimize(protein, options)`, `OptimizeCoding(cds, options)`, `BackTranslate`, `CodonFrequencyProfile`; `CodonOptimizationOptions` (`Target`, `Source`, `Method`, `Seed`, `Constraints`, `Code`, `PreserveStop`, `MaxRepairIterations`), `CodonOptimizationMethod`, `CodonOptimizationResult(Dna, Cai, GcContent, CodonsChanged, RemainingViolations, RepairIterations)`.
- `SequenceConstraints`: `Default`, `ForGeneSynthesis`, `ForbiddenSites`, `MaxHomopolymer`, `MinGc`, `MaxGc`, `GcWindow`, `AvoidHairpins`, `AvoidRepeats`, `MaxDirectRepeat`, `Check(dna)`, `IsSatisfiedBy`; `ConstraintViolation(Kind, Start, Length, Detail)`.

## Primer (`.Primers`)

- `PrimerDesigner.DesignForward`, `DesignReverse`, `DesignPairs(template, options)`, `Evaluate(oligo, …)`.
- `PrimerDesignOptions`: panjang, `MinTm`/`OptTm`/`MaxTm`, `MaxTmDifference`, jendela GC, `MaxPolyX`, GC clamp, batas ΔG dimer dan hairpin, ukuran produk, `SaltConditions`, batas jumlah hasil, dan bobot penalti; `TemplateRegion`.
- `Primer(Sequence, Start, Length, Strand, Tm, GcFraction, SelfDimerDeltaG, ThreePrimeSelfDimerDeltaG, HairpinDeltaG, EndStability, GcClamp, LongestHomopolymer, Penalty)`; `PrimerPair(Forward, Reverse, ProductSize, ProductSequence, TmDifference, HeteroDimerDeltaG, ThreePrimeHeteroDimerDeltaG, Penalty)`.
- `NearestNeighbourTm.Tm`, `Duplex`, `TmFromThermodynamics`, `WallaceTm`, `GcContentTm`, `IsSelfComplementary`; `SaltConditions.Default/StandardPcr/Qpcr`; `DuplexThermodynamics`.
- `OligoThermodynamics.SelfDimer`, `HeteroDimer`, `ThreePrimeDimer`, `Hairpin`, `HairpinLoopDeltaG37`; `OligoDuplex`, `OligoHairpin`.

## Kloning (`.Cloning`)

- `DnaFragment(Id, Sequence, FivePrimeOverhang, ThreePrimeOverhang, IsCircular)`: `Duplex`, `CanLigateTo`, `Ligate`, `Circularise`, `ReverseComplement`, `ToSequence`, `Features`.
- `GibsonAssembly.DesignOverlaps(fragments, options)`, `Simulate(fragments atau sequences, options)`, `BestOverlap`; `GibsonOptions` (panjang dan target Tm tumpang-tindih, jendela GC, panjang anil, `Circular`), `GibsonJunction`, `GibsonDesign`, `GibsonAssemblyResult`.
- `GoldenGate.Digest`, `IsDomesticated`, `AddFlanks`, `DesignOverhangs`, `AreOrthogonal`, `OverhangFidelity`, `Simulate(parts, enzyme, backbone, backboneIsCircular)`; `GoldenGateResult`.
- `RestrictionEnzyme` (mendukung Tipe IIS) dengan `All`, `FromName`, `TryFromName`, `Find(dna, circular)`, `CountSites`, `OverhangLength`, `OverhangKind`, `IsPalindromic`; `RestrictionSite`.
- `Plasmid(name, sequence)`: `Slice`, `Rotate`, `WithFeature`, `Digest(enzyme)`, `Fragments(enzymes)`, `UniqueCutters`, `ToGenBankSummary`, `ToSequence`, `DemonstrationBackbone`.

## Sirkuit (`.Circuits`)

- `BioPart(Name, Type, Sequence, Description)` dengan `RelativeStrength`, `RegulatedBy`, `IsRepressible`, `IsInducible`, `Inducer`, `RegistryId`, `ToFeature`; `PartType`; `TranscriptionalUnit(Name, Parts)` dengan `Assemble`, `Features`, `ToSequence`.
- `PartsRegistry.All`, `Get`, `TryGet`, `OfType`, `Search`, `BuildUnit(name, partNames)` — promotor (`J23100`, `J23106`, `J23119`, `R0010`, `R0011`, `R0040`, `R0051`), RBS (`B0030`–`B0034`), sekuens pengode (`GFP`, `mCherry`, `LacI`, `TetR`, `cI`), terminator (`B0010`, `B0012`, `B0015`), dan operator.
- `RbsStrengthModel.Predict(utr, cds)`, `RelativeStrength`, `SpacingEnergy`, `StartCodonEnergy`, konstanta `AntiShineDalgarno`, `ConsensusUtr`; `RbsPrediction`.
- `CircuitModels.ToggleSwitch`, `SimulateToggleSwitch`, `Repressilator`, `SimulateRepressilator`, `EstimateOscillationPeriod`, `InducibleExpression`, `SimulateInduction`, `SteadyStateProtein`, `SteadyStateMrna`; `HillParameters`.
- `OdeIntegrator.Rk4`, `Rk45`; delegasi `OdeSystem`, `OdeSolution(Times, States)` dengan `SeriesOf`, `InterpolateAt`, `Final`, `RejectedSteps`; `IntegratorKind`.

## Analisis CRISPR (`.Crispr`)

| Tipe | Kegunaan |
|---|---|
| `CasNuclease` | `SpCas9`, `SpCas9Nag`, `SpCas9Vqr`, `SpCas9Vrer`, `XCas9`, `SpG`, `SpRy`, `SaCas9`, `SaCas9Kkh`, `Cas12a`, `LbCas12a`, `NmeCas9`, `All`, `FromName`; `Pam`, `PamSide`, `GuideLength`, `CutSite`, `MatchesPam` |
| `GuideFinder`, `GuideSearchOptions`, `Guide` | mencacah protospacer; filter GC, homopolimer, poli-T, dan G di awal; `CutSite`, `Interval`, `FullSite`, `U6CloningOligo` |
| `OnTargetScorer`, `GuidePenalties` | Rule Set 1 Doench 2014 (`RuleSet1`), `Score`, `Penalties`, `CombinedScore` |
| `CrisprSiteIndex` | `Build(sequences, nuclease, includeNagPam)`, `Add`, `Query(spacer, maxMismatches)`, `QueryExhaustive`, `ProtospacerAt` |
| `OffTargetSearch`, `OffTargetSearchOptions`, `OffTarget`, `OffTargetProfile` | pencarian ketidakcocokan (dan gelembung opsional); `MitSpecificity`, `CfdSpecificity`, `MismatchHistogram` |
| `CfdScorer`, `MitSpecificityScorer` | bobot ketidakcocokan dan PAM CFD; skor hit dan guide MIT |
| `GuideDesigner`, `GuideDesignOptions`, `GuideDesign` | desain berperingkat yang menggabungkan skor on-target, penalti, dan spesifisitas |
| `BaseEditor`, `BaseEditing` | jendela editor terpublikasi (`Be3`, `Be4Max`, `Ye1`, `Abe710`, `Abe8E`); `PredictEdits`, `PredictProteinConsequence`, `FindEditingGuides`, `EditingWindowInterval`; `BaseEditOutcome`, `BaseEditSite`, `ProteinConsequence`, `EditConsequence` |
| `PrimeEditor`, `PegRna`, `PrimeEdit`, `NickingGuide`, `PrimeEditorOptions` | tata letak komponen pegRNA (spacer, scaffold, PBS, RTT), `Design`, `DesignNickingGuides` |

`DnaUtility` memuat pembantu bersama: `ReverseComplement`, `GcFraction`, `LongestHomopolymer`, `MatchesIupac`, `TryPack`, `Unpack`, `MismatchCount`, `NormalizeDna`, `IsUnambiguous`.

---
[Indeks referensi API](index.md) · [Tutorial biologi sintetik](../tutorials/synthetic-biology.md)
