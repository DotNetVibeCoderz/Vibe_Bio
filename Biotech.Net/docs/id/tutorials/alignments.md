# Tutorial: Penjajaran dan filogenetika

Paket: `Gravicode.Biotech.Net.Alignments`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Modul ini menyediakan matriks substitusi, penjajaran berpasangan gaya Biopython (Gotoh dengan gap afin, dijalankan di kernel Rust), penjajaran ganda progresif, pohon berbasis jarak (UPGMA, neighbour joining), serta pencarian basis data berumpan dengan E-value Karlin–Altschul.

## 1. Penjajaran berpasangan

`PairwiseAligner` mengikuti `PairwiseAligner` milik Biopython: skor berupa imbalan dan skor gap bernilai nol atau negatif. `AlignmentMode` berada di `Biotech.Net.Native` karena pemrograman dinamisnya berjalan di kernel penjajaran native.

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
Console.WriteLine(aln.Format());         // tampilan tiga baris dengan penanda kecocokan
Console.WriteLine(aln.Cigar());
Console.WriteLine($"identitas {aln.Identity:P1}, gap {aln.Gaps}");

// Penjajaran global nukleotida dengan skema match/mismatch sederhana
var dnaAligner = new PairwiseAligner { MatchScore = 2, MismatchScore = -3, OpenGapScore = -5, ExtendGapScore = -2 };
int score = dnaAligner.Score("GATTACAGATTACA", "GATCACAGATACA");
```

Matriks bawaan: `BLOSUM45/50/62/80/90`, `PAM30/70/250`, dan `NUC.4.4`. `SubstitutionMatrix.Load(nama)` memuat berdasarkan nama, dan `SubstitutionMatrix.Parse(teks)` membaca teks berformat NCBI.

## 2. Penjajaran ganda dan pohon

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
Console.WriteLine($"konsensus: {msa.Consensus()}");

var distances = DistanceMatrix.FromAlignment(msa, DistanceModel.PDistance);
var tree = TreeBuilder.NeighborJoining(distances).RootAtMidpoint();
Console.WriteLine(tree.ToNewick(decimals: 3));
```

`PhylogeneticTree.ParseNewick`, `RobinsonFoulds`, dan `Distance(leafA, leafB)` menutup operasi pohon yang biasa dipakai.

## 3. Pencarian basis data

`SequenceDatabase` mengindeks kata pada tiap subjek, memakainya sebagai umpan, lalu menilai kandidat dengan kernel skor lokal berkelompok. Setiap hit membawa skor mentah, bit score, dan E-value.

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

## Batasan

- `SequenceDatabase` adalah pencarian seed-and-score untuk basis data yang muat di memori, bukan pengganti BLAST penuh (tanpa perluasan two-hit maupun statistik berbasis komposisi).
- E-value memakai parameter Karlin–Altschul terpublikasi untuk BLOSUM62 (11/1) dan blastn (2/−3, 5/2). Untuk skema penilaian lain, bit score dan E-value bernilai `NaN` kecuali Anda mengisi `SearchOptions.Statistics` sendiri — urutkan hit tersebut berdasarkan skor mentah.
- Penjajaran progresif tidak melakukan penghalusan iteratif, jadi kualitasnya lebih dekat ke ClustalW daripada MAFFT L-INS-i.

---
Sebelumnya: [Sekuens](sequences.md) · Selanjutnya: [Genomika](genomics.md) · [Indeks tutorial](index.md)
