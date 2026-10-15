# API reference: Transcriptomics

Package `Gravicode.Biotech.Net.Transcriptomics` · namespaces `Biotech.Net.Transcriptomics`, `.Enrichment`, `.SingleCell`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Count matrices are `Dataset` objects (genes × samples); sample metadata lives in `Dataset.ColumnAnnotations`.

## Input and normalisation

- `CountMatrixReader.ReadFeatureCounts`, `ReadHtSeq`, `ReadHtSeqFiles`, `ReadTable`, `GeneLengths`, `WithGeneLengths`, `OpenText`.
- `Normalization.LibrarySizes`, `MedianOfRatiosSizeFactors` (DESeq2), `TmmFactors` (edgeR), `UpperQuartileFactors`, `Factors(method)`, `EffectiveLibrarySizes`, `Cpm`, `LogCpm`, `Tpm`, `Fpkm`, `Rpkm`, `NormalizedCounts`; enum `NormalizationMethod`.
- `Transforms.Vst`, `RlogApproximation`; `ExpressionFilter.FilterByExpr`, `MinimumCount`, `Apply`.

## Differential expression

| Type | Key members |
|---|---|
| `DESeq` | `Fit(counts, design, options)`, `Run(counts, factor, reference, treatment, options, covariates)` |
| `DESeqFit` | `SizeFactors`, `BaseMean`, `GeneDispersions`, `TrendDispersions`, `Dispersions`, `DispersionOutliers`, `Coefficients`, `StandardErrors`, `Converged`, `Results(coefficient, options)`, `VarianceStabilizingTransform()` |
| `DESeqOptions` | `FitType` (`Parametric` or `Mean`), `ShrinkDispersions`, `OutlierSD`, `MaxIterations`, `IndependentFiltering`, `Alpha`, `SizeFactors` |
| `Limma` | `LmFit`, `EBayes`, `FitFDist`, `TopTable`, `Voom`, `VoomTest`, `ModeratedTTest` |
| `DifferentialExpression` | `TTest(logExpression, reference, treatment, welch)` plus an annotation-driven overload |
| `DesignMatrix` | `Intercept`, `TwoGroup`, `FromAnnotations(data, factor, reference, covariates)`, `IndexOf`, `Rank()`, `Leverages()` |
| `DifferentialExpressionResult` | `Genes`, `BaseMean`, `Log2FoldChange`, `LfcSE`, `Stat`, `PValue`, `PAdj`, `Rows()`, `Top(n)`, `Significant(alpha, minAbsLog2FoldChange)`, `CountSignificant`, `Volcano()`, `MaPlot()`, `ToDataset()`, `WriteTsv` |

## Matrices and PCA

`DenseMatrix(rows, columns)` with `FromDataset(dataset, transpose)`, `ToDataset`, `Row`, `Column`, `Transpose`; `Pca.Fit(data, components, solver, seed, powerIterations, center)` returning `PcaResult(Scores, Loadings, ExplainedVariance, ExplainedVarianceRatio, Means)`.

## Enrichment (`.Enrichment`)

- `GeneSet(name, genes, description)`; `Gmt.Read`, `ReadFile`, `Write`.
- `OverRepresentation.Test(hits, universe, sets, minSetSize, maxSetSize)` returning `OverRepresentationResult(Set, SetSize, Overlap, Expected, FoldEnrichment, OddsRatio, PValue, PAdj, Genes)`; `HypergeometricUpperTail`.
- `Gsea.Run(ranking, sets, permutations, weight, minSize, maxSize, seed)` returning `GseaResult(Set, Size, EnrichmentScore, NormalizedEnrichmentScore, PValue, PAdj, FdrQ, LeadingEdge)`; `EnrichmentScore`.

## Single cell (`.SingleCell`)

- `SparseMatrix` (CSR): `FromTriplets`, `FromDense`, `ToDense`, `Transpose`, `SelectRows`, `SelectColumns`, `RowSums`, `ColumnSums`, `ColumnMeanVariance`, `Density`.
- `SingleCellExperiment(x, cellNames, geneIds, geneNames)`: `X`, `Layers`, `CellMetrics`, `CellLabels`, `GeneMetrics`, `GeneFlags`, `Embeddings`, `Neighbors`, `GeneIndex`, `SubsetCells`, `SubsetGenes`, `Copy`.
- `MatrixMarket.Read`, `Write`, `ReadSingleCell`, `Read10x(directory)`.
- `Preprocessing.CalculateQcMetrics`, `FilterCells`, `FilterGenes`, `NormalizeTotal`, `Log1p`, `HighlyVariableGenes(nTopGenes, flavor)`, `Scale`.
- `NeighborGraph.Build(data or sc, nNeighbors, metric or embedding)`: `Indices`, `Distances`, `Connectivities`.
- `CommunityDetection.Cluster(adjacency or sc, resolution, algorithm, seed, key)`, `Modularity`; `CommunityAlgorithm.Leiden`, `.Louvain`.
- `Umap.Fit(graph or sc, options)` with `UmapOptions` and `FindAbParams`; `Tsne.Fit(data or sc, options)` with `TsneOptions`.
- `MarkerGenes.Rank(sc, labels or key, test)` returning `MarkerGene(Gene, GeneIndex, Score, Log2FoldChange, PValue, PAdj, PctInGroup, PctRest)`; enum `MarkerTest`.

Embedding keys follow scanpy: `X_pca`, `X_umap`, `X_tsne`; cluster labels land in `CellLabels["leiden"]` or `["louvain"]`.

---
[API reference index](index.md) · [Transcriptomics tutorial](../tutorials/transcriptomics.md)
