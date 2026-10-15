# API reference: Visualization

Package `Gravicode.Biotech.Net.Visualization` · namespace `Biotech.Net.Visualization`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## BioPlot (base class)

`Width`, `Height`, `Title`, `Theme`, `Render(renderer)`, `ToPng(scale)`, `ToSvg()`, `ToPdf()`, `ToPngDataUri(scale)`, `ToHtml(interactive)`, `ToHtmlDocument(interactive)`, `ToPlotly()`, `SupportsInteractive`, `Save(path)`, `Save(path, format)`; `PlotFormat`. `FigurePlot` adds `XLabel`, `YLabel` and `ToFigure()` for the Matplotlib.Net-backed plots.

## Theme

`BioTheme.Light`, `BioTheme.Dark`, `BioTheme.Default` (settable), `With(change)`; `Background`, `Panel`, `Foreground`, `Muted`, `Grid`, `Border`, `Accent`, `Palette`, `Up`, `Down`, `Neutral`, `BaseA/C/G/T/Other`, `FontFamily`, `FontSize`, `DivergingColormap`, `SequentialColormap`, `PlotlyTemplate`, `CategoryColor(index)`, `NucleotideColor(symbol)`, `AminoAcidColor(residue)`, `ResidueColor`, `ElementColor`, `Font(size, bold)`, `TitleFont`, `SmallFont`, `ToPlotStyle()`; `BiotechColormaps.EnsureRegistered()`.

## Expression plots

- `VolcanoPlot(genes, log2FoldChange, pValue, pAdj)` or `VolcanoPlot(DifferentialExpressionResult)`; `UseAdjustedPValue`.
- `MaPlot(genes, meanExpression, log2FoldChange, pValue, pAdj, meanIsLog)` or `MaPlot(DifferentialExpressionResult)`.
- Shared on `DifferentialPlot`: `Alpha`, `FoldChangeThreshold`, `LabelTop`, `Highlight`, `PointSize`, `Classify()`, `Counts()`, `LabelledGenes()`; enum `Regulation`.
- `EmbeddingPlot(x, y, groups, labels)`, `FromPca(pca, groups, sampleNames)`, `FromEmbedding(embedding, groups, method)`, `FromSingleCell(sc, embedding, colorBy)`; `Values`, `ValueLabel`, `ShowLabels`, `PointSize`, `GroupLevels()`.
- `Heatmap(values, rowLabels, columnLabels)`, `Heatmap.FromDataset(dataset, groupAnnotation)`: `Scale` (`HeatmapScale`), `ClusterRows`, `ClusterColumns`, `ShowDendrograms`, `Metric`, `Linkage`, `Colormap`, `Min`, `Max`, `ZClip`, `ColumnGroups`, `MaxRowLabels`, `MaxColumnLabels`, `ColorbarLabel`, `ScaledValues()`, `ComputeLayout()`.
- `HierarchicalClustering.Cluster(observations or distances, metric, linkage)`, `DistanceMatrix`; `Dendrogram` with `Merges`, `Order`, `Cut(k)`, `Segments()`; `DistanceMetric`, `LinkageMethod`.

## Genome browser

`GenomeBrowser(region, width)`: `AddRuler()`, `AddGenes(genes, name, showTranscripts)`, `AddTranscripts`, `AddFeatures`, `AddIntervals`, `AddCoverage(values, start, name, color)`, `AddVariants(variants, name)`, `AddSequence(sequence, offset, name)`, `AddTrack(track)`, `Tracks`, `NameGutter`, `Padding`, `TrackGap`, `BasesPerPixel`, `ToX(position, area)`, `ComputeLayout()`.

Tracks: `RulerTrack`, `GeneTrack` (`ShowTranscripts`, `RowHeight`, `Models()`, `PackRows`), `FeatureTrack` (`FromIntervals`, `TypeColor`), `CoverageTrack` (`TrackHeight`, `Color`, `Max`, `Bin`), `VariantTrack` (`TrackHeight`, `ShowLabels`, `StemHeights`, `TypeColor`), `SequenceTrack`. Derive from `GenomeTrack` (`MeasureHeight`, `Draw`) for custom tracks; `TrackContext` gives the renderer, area, theme and `X(position)`.

## Sequence and structure plots

- `SequenceLogo(alignedSequences)` or `FromMatrix(matrix, symbols, sampleSize)`: `Mode` (`LogoMode`), `SmallSampleCorrection`, `FirstPosition`, `MaxBits`, `ComputeColumns()`; `LogoColumn`, `LogoLetter`.
- `DotPlot(first, second, window, minMatches)`: `IncludeReverseComplement`, `MaxDots`, `FindMatches()`; `DotMatch`.
- `ProteinViewer.CaTrace(structure, chainId)`, `ContactMap(structure, chainId, cutoff)`, `To3DHtml(structure, options)`, `Html(modelText, format, options, title)`, `ThreeDmolVersion`, `ThreeDmolUrl`; `CaTracePlot` (`Coloring`, `LabelEvery`, `Project()`), `ContactMapPlot` (`Cutoff`, `ShowDistances`, `Matrix()`); `Viewer3DOptions`, `MolecularStyle`, `ResidueColoring`.
- `MoleculeViewer(molecule, width, height)`: `ShowCarbonLabels`, `ShowImplicitHydrogens`, `ShowAtomIndices`, `ColorByElement`, `ShowFormula`, `HighlightAtoms`, `Coordinates2D()`, `To3DHtml(options)`.

---
[API reference index](index.md) · [Visualization tutorial](../tutorials/visualization.md)
