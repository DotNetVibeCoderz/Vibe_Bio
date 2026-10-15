# Tutorial: Notebooks

Package: `Gravicode.Biotech.Net.Notebooks`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Rich HTML formatters for .NET Interactive, so that a `Sequence`, a `Dataset`, a differential-expression table, a molecule, a structure, an alignment, a tree or any plot renders as a readable card instead of `ToString()` output. Works in Polyglot Notebooks (VS Code), Visual Studio, Jupyter with the .NET kernel, and anywhere else `Microsoft.DotNet.Interactive.Formatting` is loaded.

## 1. Set-up in a notebook

```text
#r "nuget: Gravicode.Biotech.Net.Notebooks, 0.1.0"
```

```csharp
using Biotech.Net.Core;
using Biotech.Net.Notebooks;

BiotechFormatters.Register();
Console.WriteLine(BiotechInfo.Banner);
Console.WriteLine($"{BiotechFormatters.FormattedTypes.Count} types get rich output: " +
                  string.Join(", ", BiotechFormatters.FormattedTypes.Select(t => t.Name)));
```

After `Register()`, the last expression of a cell renders as HTML:

```csharp
using Biotech.Net.Core;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
// In a notebook cell, `dna` on its own line now shows a coloured, wrapped residue card.
Console.WriteLine(dna.Preview());
```

`BiotechFormatters.Reset()` restores the default .NET Interactive formatters.

## 2. Producing the HTML yourself

`BiotechHtml` and `AlignmentHtml` build the same markup outside a notebook — handy for reports, Blazor pages or e-mail.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Notebooks;
using Biotech.Net.Visualization;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", "demo");
string card = BiotechHtml.Sequence(dna);
string footer = BiotechHtml.Footer();          // the Gravicode credit line

var logo = new SequenceLogo(["ACGTG", "ACGTA", "ACGTG", "TCGTG"]);
string svgPlot = BiotechHtml.Plot(logo, PlotEmbedding.Svg);
string pngPlot = BiotechHtml.Plot(logo, PlotEmbedding.Png);

BiotechHtml.LineWidth = 80;        // residues per line
BiotechHtml.MaxResidues = 5_000;   // truncate longer sequences in the card
File.WriteAllText(Path.Combine(Path.GetTempPath(), "report.html"),
    $"<!doctype html><meta charset=\"utf-8\">{BiotechHtml.Styles}{card}{svgPlot}{footer}");
Console.WriteLine($"{card.Length + svgPlot.Length} characters of HTML");
```

For alignments:

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
Console.WriteLine($"{html.Length} characters of alignment HTML");  // residues are emitted as individually coloured spans
```

## 3. The sample notebooks

The repository ships six `.dib` notebooks under `notebooks/`, each bilingual (English and Bahasa Indonesia):

| Notebook | Contents |
|---|---|
| `01-getting-started.dib` | packages, formatters, the core types |
| `02-sequences-and-io.dib` | FASTA/FASTQ/GenBank, ORFs, motifs, restriction digests |
| `03-genome-browser-and-plots.dib` | genome browser, heatmap, sequence logo, static vs interactive export |
| `04-rnaseq-differential-expression.dib` | DESeq run, volcano, MA, PCA, clustered heatmap |
| `05-protein-structure.dib` | ProtParam, DSSP, SASA, superposition, 3-D viewer |
| `06-cheminformatics.dib` | SMILES, descriptors, fingerprints, similarity, depiction |

Open them in VS Code with the Polyglot Notebooks extension, or run them from Jupyter with the .NET kernel. Each one starts with a `#r "nuget: ..."` line and a commented alternative that references the local build instead.

## Limitations

- Rich output needs `Microsoft.DotNet.Interactive.Formatting`; in a plain console app the formatters are not used (call `BiotechHtml` directly instead, as above).
- The HTML cards are self-contained except for the 3-D viewers, which load 3Dmol.js from a CDN.
- Very long sequences are truncated in the card at `BiotechHtml.MaxResidues`; alignments at `AlignmentHtml.MaxColumns`.

---
Previous: [Visualization](visualization.md) · [Tutorial index](index.md)
