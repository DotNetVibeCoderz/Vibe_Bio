# Tutorial: Synthetic biology

Package: `Gravicode.Biotech.Net.SyntheticBiology`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Codon optimisation with manufacturability constraints, primer design (nearest-neighbour Tm, dimers, hairpins), Gibson and Golden Gate assembly simulation, plasmid handling, a ribosome-binding-site strength model, a parts registry and ODE simulation of genetic circuits. The module also contains the CRISPR analysis types summarised at the end of this page.

## 1. Codon optimisation

`CodonUsageTable` ships usage tables for *E. coli* K-12, *S. cerevisiae*, human, CHO (*C. griseus*) and *P. pastoris*. `SequenceConstraints` describes what a synthesis vendor will accept, and the optimiser repairs the first draft until the constraints hold.

```csharp
using Biotech.Net.SyntheticBiology.Codons;

var options = new CodonOptimizationOptions
{
    Target = CodonUsageTable.EColiK12,
    Method = CodonOptimizationMethod.MostFrequent,
    Constraints = SequenceConstraints.ForGeneSynthesis,
};
var result = CodonOptimizer.Optimize("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", options);
Console.WriteLine($"CAI {result.Cai:F3}, GC {result.GcContent:P1}, {result.RemainingViolations.Count} violations left");
Console.WriteLine(result.Dna);

// Re-code an existing gene for another host, keeping the protein identical
var humanOptions = new CodonOptimizationOptions
{
    Target = CodonUsageTable.HomoSapiens,
    Constraints = SequenceConstraints.ForGeneSynthesis,
};
var recoded = CodonOptimizer.OptimizeCoding(result.Dna, humanOptions);
Console.WriteLine($"codons changed: {recoded.CodonsChanged}");
```

## 2. Primer design

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;
using Biotech.Net.SyntheticBiology.Primers;

string template = RandomSequence.Dna(1_200, seed: 8, gcContent: 0.5).ToString();
var pairs = PrimerDesigner.DesignPairs(template, new PrimerDesignOptions
{
    OptTm = 60,
    ProductSizeMin = 300,
    ProductSizeMax = 600,
    MaxPairsReturned = 3,
});
foreach (var p in pairs)
    Console.WriteLine($"{p.Forward.Sequence} / {p.Reverse.Sequence}  {p.ProductSize} bp  Tm {p.Forward.Tm:F1}/{p.Reverse.Tm:F1}  penalty {p.Penalty:F2}");

// Evaluate an existing oligo
var primer = PrimerDesigner.Evaluate("AGCGGATAACAATTTCACACAGG");
Console.WriteLine($"Tm {primer.Tm:F1} °C, GC {primer.GcFraction:P0}, hairpin dG {primer.HairpinDeltaG:F2} kcal/mol");
Console.WriteLine($"nearest-neighbour Tm: {NearestNeighbourTm.Tm("AGCGGATAACAATTTCACACAGG"):F1} °C");
```

## 3. Assembly simulation

```csharp
using Biotech.Net.SyntheticBiology.Cloning;

// Gibson: design overlaps between fragments, then simulate the join
string a = "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTC";
string b = "AGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACA";
var design = GibsonAssembly.DesignOverlaps([a, b]);
Console.WriteLine($"{design.Junctions.Count} junction(s), ideal: {design.IsIdeal}");
var gibson = GibsonAssembly.Simulate(design.TailedFragments);
Console.WriteLine(gibson);

// Golden Gate with a Type IIS enzyme
var bsaI = RestrictionEnzyme.FromName("BsaI");
string part1 = GoldenGate.AddFlanks(a, bsaI, "AATG", "GCTT");
string part2 = GoldenGate.AddFlanks(b, bsaI, "GCTT", "AATG");
var gg = GoldenGate.Simulate([part1, part2], bsaI);
Console.WriteLine(gg);

// Plasmids: unique cutters, digests and a GenBank-style summary
var plasmid = Plasmid.DemonstrationBackbone;
Console.WriteLine($"{plasmid}  unique cutters: {string.Join(", ", plasmid.UniqueCutters().Take(5).Select(e => e.Name))}");
```

## 4. RBS strength and genetic circuits

```csharp
using Biotech.Net.SyntheticBiology.Circuits;

var rbs = RbsStrengthModel.Predict("TAAGGAGGTAAAAA", "ATGAAAGAAGCAGAAGCAAAAGAAGCATAA");
Console.WriteLine($"initiation rate {rbs.TranslationInitiationRate:F0} (a.u.), spacing {rbs.SpacingNucleotides} nt, SD match {rbs.ShineDalgarnoMatch}");

// Parts registry and transcriptional units
var unit = PartsRegistry.BuildUnit("reporter", "J23100", "B0034", "GFP", "B0015");
Console.WriteLine($"{unit}  GC {unit.ToSequence().GcContent():P1}");

// Elowitz repressilator, integrated with RK4
var solution = CircuitModels.SimulateRepressilator(tEnd: 300);
double period = CircuitModels.EstimateOscillationPeriod(solution, stateIndex: 3);
Console.WriteLine($"{solution.StepCount} output steps, period ≈ {period:F1} time units");

// Toggle switch and inducible expression are also provided
var toggle = CircuitModels.SimulateToggleSwitch(alpha1: 50, alpha2: 16, beta: 2.5, gamma: 1,
    initial: [0.1, 0.2], tEnd: 50);
Console.WriteLine($"final state: {string.Join(", ", toggle.Final.Select(v => v.ToString("F3")))}");
```

## 5. CRISPR analysis types

The `Biotech.Net.SyntheticBiology.Crispr` namespace models the published scoring schemes so that guide sets can be analysed and compared in .NET:

| Type | Purpose |
|---|---|
| `CasNuclease` | PAM, guide length and cut geometry for SpCas9, SaCas9, Cas12a, SpG, SpRY, NmeCas9 and engineered variants (`CasNuclease.All`) |
| `GuideFinder`, `Guide` | enumerates PAM-adjacent protospacers in a sequence, with GC, homopolymer and poly-T checks |
| `OnTargetScorer` | Doench et al. 2014 Rule Set 1 on the 30-nt context |
| `CrisprSiteIndex`, `OffTargetSearch` | seeded mismatch search over indexed sequences |
| `CfdScorer`, `MitSpecificityScorer` | CFD (Doench et al. 2016) and MIT specificity aggregation |
| `GuideDesigner`, `GuideDesign` | combines the above into a ranked table |
| `BaseEditing`, `BaseEditor` | editing-window models for published cytosine and adenine base editors |
| `PrimeEditor`, `PegRna` | pegRNA component layout (spacer, scaffold, PBS, RTT) and nicking-guide candidates |

See the [API reference](../api-reference/synthetic-biology.md) for the member lists, and the CRISPR Designer page of the [Gallery](../samples.md) for an interactive view.

## Limitations

- On-target scoring is Rule Set 1 (the 2014 model). Rule Set 2 / Azimuth and DeepCRISPR-style models are not included, and chromatin accessibility is not modelled at all.
- Off-target search is sequence-based over the sequences you index. It is not a whole-genome search service, and `MaxBulges > 0` is slower because bulged alignments are enumerated explicitly.
- Base and prime editor models are window-weight and layout models from the primary literature, not trained efficiency predictors; use them to compare candidates, not to predict absolute editing rates.
- The RBS model is a simplified free-energy model in the spirit of the RBS Calculator, reported in arbitrary units. It is not calibrated against a specific host data set.
- Circuit simulation covers the classical ODE models (toggle switch, repressilator, inducible expression) with RK4/RK45. There is no stochastic (Gillespie) simulation.

---
Previous: [Drug discovery](drug-discovery.md) · Next: [AI](ai.md) · [Tutorial index](index.md)
