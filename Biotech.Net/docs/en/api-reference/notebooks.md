# API reference: Notebooks

Package `Gravicode.Biotech.Net.Notebooks` · namespace `Biotech.Net.Notebooks`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## BiotechFormatters

| Member | Purpose |
|---|---|
| `Register()` | registers HTML formatters with .NET Interactive and makes HTML the preferred MIME type for the formatted types (idempotent) |
| `Reset()` | restores the default formatters |
| `IsRegistered` | whether `Register()` has run |
| `FormattedTypes` | the 14 types with rich output: `Sequence`, `SequenceRead`, `Protein`, `Variant`, `GenomicInterval`, `Gene`, `Dataset`, `DifferentialExpressionResult`, `Molecule`, `Structure`, `BioPlot`, `PairwiseAlignment`, `MultipleSequenceAlignment`, `PhylogeneticTree` |
| `BannerHtml()` | the product banner as an HTML card |

## BiotechHtml

Produces the same markup outside a notebook.

`Styles` (the shared CSS), `Footer()` (the Gravicode credit line), `Sequence`, `Read`, `Protein`, `Variant`, `Interval`, `Gene`, `Dataset(ds, maxRows, maxColumns)`, `DifferentialExpression(result, top)`, `Molecule`, `Structure(structure, options)`, `Residues(residues, nucleotide, lineWidth, maxResidues)`, `Plot(plot, embedding)`.

Settings: `DefaultPlotEmbedding` (`PlotEmbedding.Auto`, `.Svg`, `.Png`), `LineWidth`, `MaxResidues`.

## AlignmentHtml

`Pairwise(alignment, width)`, `Multiple(msa, width, maxSequences)`, `Tree(tree, width)`, and `MaxColumns` to cap very wide alignments.

## Sample notebooks

`notebooks/*.dib`, bilingual (EN/ID): getting started, sequences and I/O, genome browser and plots, RNA-seq differential expression, protein structure, cheminformatics. See [Samples](../samples.md).

---
[API reference index](index.md) · [Notebooks tutorial](../tutorials/notebooks.md)
