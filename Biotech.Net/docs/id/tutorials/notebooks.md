# Tutorial: Notebook

Paket: `Gravicode.Biotech.Net.Notebooks`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Formatter HTML kaya untuk .NET Interactive, sehingga `Sequence`, `Dataset`, tabel ekspresi diferensial, molekul, struktur, penjajaran, pohon, atau plot apa pun tampil sebagai kartu yang mudah dibaca, bukan sekadar keluaran `ToString()`. Bekerja di Polyglot Notebooks (VS Code), Visual Studio, Jupyter dengan kernel .NET, dan di mana pun `Microsoft.DotNet.Interactive.Formatting` dimuat.

## 1. Persiapan di notebook

```text
#r "nuget: Gravicode.Biotech.Net.Notebooks, 0.2.0"
```

```csharp
using Biotech.Net.Core;
using Biotech.Net.Notebooks;

BiotechFormatters.Register();
Console.WriteLine(BiotechInfo.Banner);
Console.WriteLine($"{BiotechFormatters.FormattedTypes.Count} tipe mendapat keluaran kaya: " +
                  string.Join(", ", BiotechFormatters.FormattedTypes.Select(t => t.Name)));
```

Setelah `Register()`, ekspresi terakhir sebuah sel akan dirender sebagai HTML:

```csharp
using Biotech.Net.Core;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
// Di sel notebook, menulis `dna` sendirian kini menampilkan kartu residu berwarna yang terbungkus rapi.
Console.WriteLine(dna.Preview());
```

`BiotechFormatters.Reset()` mengembalikan formatter bawaan .NET Interactive.

## 2. Menghasilkan HTML sendiri

`BiotechHtml` dan `AlignmentHtml` membangun markup yang sama di luar notebook — berguna untuk laporan, halaman Blazor, atau surel.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Notebooks;
using Biotech.Net.Visualization;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
string card = BiotechHtml.Sequence(dna);
string footer = BiotechHtml.Footer();          // baris kredit Gravicode

var logo = new SequenceLogo(["ACGTG", "ACGTA", "ACGTG", "TCGTG"]);
string svgPlot = BiotechHtml.Plot(logo, PlotEmbedding.Svg);
string pngPlot = BiotechHtml.Plot(logo, PlotEmbedding.Png);

BiotechHtml.LineWidth = 80;        // residu per baris
BiotechHtml.MaxResidues = 5_000;   // memotong sekuens yang lebih panjang di kartu
File.WriteAllText(Path.Combine(Path.GetTempPath(), "report.html"),
    $"<!doctype html><meta charset=\"utf-8\">{BiotechHtml.Styles}{card}{svgPlot}{footer}");
Console.WriteLine($"{card.Length + svgPlot.Length} karakter HTML");
```

Untuk penjajaran:

```csharp
using Biotech.Net.Alignments;
using Biotech.Net.Native;
using Biotech.Net.Notebooks;

var aln = new PairwiseAligner
{
    Mode = AlignmentMode.Local,
    SubstitutionMatrix = SubstitutionMatrix.Blosum62,
    OpenGapScore = -11,
    ExtendGapScore = -1,
}.Align("HEAGAWGHEE", "PAWHEAE");

string html = AlignmentHtml.Pairwise(aln, width: 60);
Console.WriteLine($"{html.Length} karakter HTML penjajaran");  // residu dipancarkan sebagai span berwarna satu per satu
```

## 3. Notebook contoh

Repositori menyertakan enam notebook `.dib` di `notebooks/`, semuanya dwibahasa (Inggris dan Bahasa Indonesia):

| Notebook | Isi |
|---|---|
| `01-getting-started.dib` | paket, formatter, tipe inti |
| `02-sequences-and-io.dib` | FASTA/FASTQ/GenBank, ORF, motif, digesti restriksi |
| `03-genome-browser-and-plots.dib` | genome browser, heatmap, sequence logo, ekspor statis vs interaktif |
| `04-rnaseq-differential-expression.dib` | eksekusi DESeq, volcano, MA, PCA, heatmap terklaster |
| `05-protein-structure.dib` | ProtParam, DSSP, SASA, superposisi, penampil 3-D |
| `06-cheminformatics.dib` | SMILES, deskriptor, sidik jari, kemiripan, penggambaran |

Bukalah di VS Code dengan ekstensi Polyglot Notebooks, atau jalankan dari Jupyter dengan kernel .NET. Setiap notebook diawali baris `#r "nuget: ..."` beserta alternatif terkomentar yang merujuk hasil build lokal.

## Batasan

- Keluaran kaya memerlukan `Microsoft.DotNet.Interactive.Formatting`; di aplikasi konsol biasa formatter tidak dipakai (panggil `BiotechHtml` langsung seperti contoh di atas).
- Kartu HTML bersifat mandiri kecuali penampil 3-D yang memuat 3Dmol.js dari CDN.
- Sekuens yang sangat panjang dipotong di kartu sesuai `BiotechHtml.MaxResidues`, dan penjajaran sesuai `AlignmentHtml.MaxColumns`.

---
Sebelumnya: [Visualisasi](visualization.md) · Selanjutnya: [Data cloud](data.md) · [Indeks tutorial](index.md)
