# Referensi API: Visualization

Paket `Gravicode.Biotech.Net.Visualization` · namespace `Biotech.Net.Visualization`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## BioPlot (kelas dasar)

`Width`, `Height`, `Title`, `Theme`, `Render(renderer)`, `ToPng(scale)`, `ToSvg()`, `ToPdf()`, `ToPngDataUri(scale)`, `ToHtml(interactive)`, `ToHtmlDocument(interactive)`, `ToPlotly()`, `SupportsInteractive`, `Save(path)`, `Save(path, format)`; `PlotFormat`. `FigurePlot` menambahkan `XLabel`, `YLabel`, dan `ToFigure()` untuk plot berbasis Matplotlib.Net.

## Tema

`BioTheme.Light`, `BioTheme.Dark`, `BioTheme.Default` (dapat diubah), `With(change)`; `Background`, `Panel`, `Foreground`, `Muted`, `Grid`, `Border`, `Accent`, `Palette`, `Up`, `Down`, `Neutral`, `BaseA/C/G/T/Other`, `FontFamily`, `FontSize`, `DivergingColormap`, `SequentialColormap`, `PlotlyTemplate`, `CategoryColor(index)`, `NucleotideColor(symbol)`, `AminoAcidColor(residue)`, `ResidueColor`, `ElementColor`, `Font(size, bold)`, `TitleFont`, `SmallFont`, `ToPlotStyle()`; `BiotechColormaps.EnsureRegistered()`.

## Plot ekspresi

- `VolcanoPlot(genes, log2FoldChange, pValue, pAdj)` atau `VolcanoPlot(DifferentialExpressionResult)`; `UseAdjustedPValue`.
- `MaPlot(genes, meanExpression, log2FoldChange, pValue, pAdj, meanIsLog)` atau `MaPlot(DifferentialExpressionResult)`.
- Bersama di `DifferentialPlot`: `Alpha`, `FoldChangeThreshold`, `LabelTop`, `Highlight`, `PointSize`, `Classify()`, `Counts()`, `LabelledGenes()`; enum `Regulation`.
- `EmbeddingPlot(x, y, groups, labels)`, `FromPca(pca, groups, sampleNames)`, `FromEmbedding(embedding, groups, method)`, `FromSingleCell(sc, embedding, colorBy)`; `Values`, `ValueLabel`, `ShowLabels`, `PointSize`, `GroupLevels()`.
- `Heatmap(values, rowLabels, columnLabels)`, `Heatmap.FromDataset(dataset, groupAnnotation)`: `Scale` (`HeatmapScale`), `ClusterRows`, `ClusterColumns`, `ShowDendrograms`, `Metric`, `Linkage`, `Colormap`, `Min`, `Max`, `ZClip`, `ColumnGroups`, `MaxRowLabels`, `MaxColumnLabels`, `ColorbarLabel`, `ScaledValues()`, `ComputeLayout()`.
- `HierarchicalClustering.Cluster(observations atau distances, metric, linkage)`, `DistanceMatrix`; `Dendrogram` dengan `Merges`, `Order`, `Cut(k)`, `Segments()`; `DistanceMetric`, `LinkageMethod`.

## Genome browser

`GenomeBrowser(region, width)`: `AddRuler()`, `AddGenes(genes, name, showTranscripts)`, `AddTranscripts`, `AddFeatures`, `AddIntervals`, `AddCoverage(values, start, name, color)`, `AddVariants(variants, name)`, `AddSequence(sequence, offset, name)`, `AddTrack(track)`, `Tracks`, `NameGutter`, `Padding`, `TrackGap`, `BasesPerPixel`, `ToX(position, area)`, `ComputeLayout()`.

Track: `RulerTrack`, `GeneTrack` (`ShowTranscripts`, `RowHeight`, `Models()`, `PackRows`), `FeatureTrack` (`FromIntervals`, `TypeColor`), `CoverageTrack` (`TrackHeight`, `Color`, `Max`, `Bin`), `VariantTrack` (`TrackHeight`, `ShowLabels`, `StemHeights`, `TypeColor`), `SequenceTrack`. Turunkan `GenomeTrack` (`MeasureHeight`, `Draw`) untuk track kustom; `TrackContext` menyediakan renderer, area, tema, dan `X(position)`.

## Plot sekuens dan struktur

- `SequenceLogo(alignedSequences)` atau `FromMatrix(matrix, symbols, sampleSize)`: `Mode` (`LogoMode`), `SmallSampleCorrection`, `FirstPosition`, `MaxBits`, `ComputeColumns()`; `LogoColumn`, `LogoLetter`.
- `DotPlot(first, second, window, minMatches)`: `IncludeReverseComplement`, `MaxDots`, `FindMatches()`; `DotMatch`.
- `ProteinViewer.CaTrace(structure, chainId)`, `ContactMap(structure, chainId, cutoff)`, `To3DHtml(structure, options)`, `Html(modelText, format, options, title)`, `ThreeDmolVersion`, `ThreeDmolUrl`; `CaTracePlot` (`Coloring`, `LabelEvery`, `Project()`), `ContactMapPlot` (`Cutoff`, `ShowDistances`, `Matrix()`); `Viewer3DOptions`, `MolecularStyle`, `ResidueColoring`.
- `MoleculeViewer(molecule, width, height)`: `ShowCarbonLabels`, `ShowImplicitHydrogens`, `ShowAtomIndices`, `ColorByElement`, `ShowFormula`, `HighlightAtoms`, `Coordinates2D()`, `To3DHtml(options)`.

---
[Indeks referensi API](index.md) · [Tutorial visualisasi](../tutorials/visualization.md)
