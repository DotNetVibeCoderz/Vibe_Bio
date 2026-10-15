# Tutorial: Alignments and phylogenetics

Package: `Gravicode.Biotech.Net.Alignments`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

The module provides substitution matrices, Biopython-style pairwise alignment (Gotoh affine gaps, run on the Rust kernel), progressive multiple alignment, distance-based trees (UPGMA, neighbour joining) and a seeded database search with Karlin–Altschul E-values.

## 1. Pairwise alignment

`PairwiseAligner` follows Biopython's `PairwiseAligner`: scores are rewards and gap scores are zero or negative. `AlignmentMode` lives in `Biotech.Net.Native`, because the dynamic programming runs in the native alignment kernel.

```csharp
using Biotech.Net.Alignments;
using Biotech.Net.Native;

var aligner = new PairwiseAligner
{
    Mode = AlignmentMode.Local,
    SubstitutionMatrix = SubstitutionMatrix.Blosum62,
    OpenGapScore = -11,
    ExtendGapScore = -1,
};
var aln = aligner.Align("HEAGAWGHEE", "PAWHEAE");
Console.WriteLine(aln.Score);            // 17
Console.WriteLine(aln.Format());        // three-line view with match markers
Console.WriteLine(aln.Cigar());
Console.WriteLine($"identity {aln.Identity:P1}, gaps {aln.Gaps}");

// Nucleotide global alignment with a simple match/mismatch scheme
var dnaAligner = new PairwiseAligner { MatchScore = 2, MismatchScore = -3, OpenGapScore = -5, ExtendGapScore = -2 };
int score = dnaAligner.Score("GATTACAGATTACA", "GATCACAGATACA");
```

Built-in matrices: `BLOSUM45/50/62/80/90`, `PAM30/70/250` and `NUC.4.4`. `SubstitutionMatrix.Load(name)` loads one by name, and `SubstitutionMatrix.Parse(text)` reads NCBI-format text.

## 2. Multiple sequence alignment and trees

```csharp
using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Core;

Sequence[] globins =
[
    Sequence.Protein("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF", "HBA_HUMAN"),
    Sequence.Protein("VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESF", "HBB_HUMAN"),
    Sequence.Protein("VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRF", "MYG_HUMAN"),
    Sequence.Protein("VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHF", "HBA_MOUSE"),
];

var msa = new ProgressiveAligner { SubstitutionMatrix = SubstitutionMatrix.Blosum62 }.Align(globins);
Console.WriteLine(msa.ToClustal());
Console.WriteLine($"consensus: {msa.Consensus()}");

var distances = DistanceMatrix.FromAlignment(msa, DistanceModel.PDistance);
var tree = TreeBuilder.NeighborJoining(distances).RootAtMidpoint();
Console.WriteLine(tree.ToNewick(decimals: 3));
```

`PhylogeneticTree.ParseNewick`, `RobinsonFoulds` and `Distance(leafA, leafB)` cover the usual tree operations.

## 3. Database search

`SequenceDatabase` indexes words of each subject, uses them as seeds and scores candidates with the batched local-score kernel. Hits carry raw scores, bit scores and E-values.

```csharp
using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Search;
using Biotech.Net.Core;

var db = new SequenceDatabase(
[
    Sequence.Protein("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHV", "HBA_HUMAN"),
    Sequence.Protein("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGL", "HBB_HUMAN"),
    Sequence.Protein("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKW", "UNRELATED"),
]);
var query = Sequence.Protein("VLSPADKTNVKAAWGKVGAHAGEYGAEAL", "query");
foreach (var hit in db.Search(query, new SearchOptions { SubstitutionMatrix = SubstitutionMatrix.Blosum62, MaxEValue = 1 }))
    Console.WriteLine(hit);
```

## Limitations

- `SequenceDatabase` is a seed-and-score search for databases that fit in memory, not a full BLAST replacement (no two-hit extension or composition-based statistics).
- E-values use published Karlin–Altschul parameters for BLOSUM62 (11/1) and blastn (2/−3, 5/2). For other scoring schemes the bit score and E-value are `NaN` unless you set `SearchOptions.Statistics` yourself, so rank those hits by raw score.
- Progressive alignment does no iterative refinement, so expect quality closer to ClustalW than to MAFFT L-INS-i.

---
Previous: [Sequences](sequences.md) · Next: [Genomics](genomics.md) · [Tutorial index](index.md)
