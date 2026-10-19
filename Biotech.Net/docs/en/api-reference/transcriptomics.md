# API reference: Transcriptomics

Package `Gravicode.Biotech.Net.Transcriptomics` · namespaces `Biotech.Net.Transcriptomics`, `.Enrichment`, `.SingleCell`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Count matrices are `Dataset` objects (genes × samples); sample metadata lives in `Dataset.ColumnAnnotations`.

## Input and normalisation

- `CountMatrixReader.ReadFeatureCounts`, `ReadHtSeq`, `ReadHtSeqFiles`, `ReadTable`, `GeneLengths`, `WithGeneLengths`, `OpenText`.
- `Normalization.LibrarySizes`, `MedianOfRatiosSizeFactors` (DESeq2 "ratio"), `PosCountsSizeFactors` (DESeq2 "poscounts"), `TmmFactors` (edgeR), `UpperQuartileFactors`, `Factors(method)`, `EffectiveLibrarySizes`, `Cpm`, `LogCpm`, `Tpm`, `Fpkm`, `Rpkm`, `NormalizedCounts`; enum `NormalizationMethod`.
- `Transforms.Vst`, `RlogApproximation`; `ExpressionFilter.FilterByExpr`, `MinimumCount`, `Apply`.

## Differential expression

| Type | Key members |
|---|---|
| `DESeq` | `Fit(counts, design, options)` (Wald, or LRT when `options.ReducedDesign` is set), `FitLrt(counts, full, reduced, options)`, `Run(counts, factor, reference, treatment, options, covariates)` |
| `DESeqFit` | `SizeFactors`, `BaseMean`, `BaseVar`, `GeneDispersions`, `TrendDispersions`, `MapDispersions`, `Dispersions`, `DispersionOutliers`, `Trend`, `TrendType`, `TrendFunction`, `VarLogDispersionEstimates`, `DispersionPriorVariance`, `Coefficients`, `StandardErrors`, `Converged`, `BetaIterations`, `Deviance`, `Mu`, `HatDiagonals`, `CooksDistances`, `MaxCooks`, `Replaced`, `ReplacementCounts`, `LrtStatistic`, `LrtPValue`, `Test`; `Results(coefficient, options)`, `Results(name, options)`, `Results(factor, numerator, denominator, options)`, `Results(double[] contrast, options)`, `VarianceStabilizingTransform()` |
| `DESeqOptions` | `SizeFactorType` (`Ratio`, `PosCounts`), `FitType` (`Parametric`, `Local`, `Mean`), `ShrinkDispersions`, `OutlierSD`, `MinDispersion`, `MaxIterations`, `Lambda`, `BetaTolerance`, `MinMu`, `UseOptim`, `DispersionPriorVariance`, `ReducedDesign`, `MinReplicatesForReplace` (7), `CooksFiltering`, `CooksCutoff`, `IndependentFiltering`, `Alpha`, `SizeFactors` |
| `LfcShrink` | `Apeglm(fit, coefficient)`, `Normal(fit, coefficient)`, `Shrink(fit, coefficient, LfcShrinkType)` returning `LfcShrinkResult` (`Table`, `PriorVariance`, `PriorScale`, `SValue`, `FalseSignRate`, `Convergence`) |
| `Limma` | `LmFit`, `EBayes`, `FitFDist`, `TopTable`, `Voom`, `VoomTest`, `ModeratedTTest` |
| `DifferentialExpression` | `TTest(logExpression, reference, treatment, welch)` plus an annotation-driven overload |
| `DesignMatrix` | `Intercept`, `TwoGroup`, `FromAnnotations(data, factor, reference, covariates)`, `IndexOf`, `Rank()`, `Leverages()`, `Factors` (`DesignFactor`: name, levels, columns) |
| `DifferentialExpressionResult` | `Genes`, `BaseMean`, `Log2FoldChange`, `LfcSE`, `Stat`, `PValue`, `PAdj`, `Rows()`, `Top(n)`, `Significant(alpha, minAbsLog2FoldChange)`, `CountSignificant`, `Volcano()`, `MaPlot()`, `ToDataset()`, `WriteTsv` |

`DESeq` and `LfcShrink` are ports of DESeq2 1.52.0 and apeglm 1.34.0 tested against R on the pasilla and airway data sets: closed-form stages agree to ~1e-12 relative, iterative stages to 1e-6 for at least 99% of genes (1e-3 for the rest, from line-search stopping points). See the [tutorial](../tutorials/transcriptomics.md) for details and remaining differences.

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

## AnnData / h5ad and atlas scale (`.SingleCell`)

HDF5 access is pure managed code ([PureHDF](https://github.com/Apollo3zehn/PureHDF) 2.2.0, MIT): no native HDF5 library.

- `H5ad.Read(path)` → `AnnData`; `H5ad.ReadSingleCell(path, geneSymbolColumn)`; `H5ad.Write(path, AnnData or SingleCellExperiment, options)`; `H5ad.OpenBacked(path)` → `BackedAnnData`; `H5ad.WriteStreaming(path, x, obs, var, layers, obsm, uns, options, blockRows)` writes X and layers from out-of-core sources block by block.
- `H5adWriteOptions`: `Compression` (`H5adCompression.None` or `.Gzip` = shuffle + deflate), `GzipLevel` (1 or 9), `Shuffle`, `ChunkElements`, `CompressionThreshold`.
- `AnnData(x, obs, var)`: `X`, `Obs`, `Var`, `Obsm`, `Varm`, `Obsp`, `Varp`, `Layers`, `Uns`, `Raw`, `ToSingleCellExperiment(geneSymbolColumn)`, `FromSingleCellExperiment(sc)`.
- `AnnMatrix`: `FromSparse(matrix, format, dataType, values64)`, `FromDense(matrix, dataType)`, `Format` (`AnnMatrixFormat.Dense`, `.Csr`, `.Csc`), `DataType` (numpy dtype name), `Sparse`, `Dense`, `SparseValues64`, `ToSparse`, `ToDense`, indexer.
- `AnnDataFrame(index, indexName)`: `Index`, `IndexName`, `ColumnNames`, indexer, `Set`, `Remove`, `Contains`, `GetDoubles`, `GetStrings`, `SelectRows`. Columns are primitive arrays, `string[]`, `Categorical(codes, categories, ordered)` (`FromLabels`, `ToStrings`) or `MaskedArray(values, mask)` (nullable integer/boolean/string, `NaValue`).
- `AnnArray(data, shape)` (N-d `uns` arrays), `AnnRaw(x, var)` with `Varm`, `AnnUnsupported(Path, EncodingType, Reason)`.
- `BackedAnnData`: `X`, `Obs`, `Var`, `Uns`, `LayerNames`, `Layer(name)`, `ObsmKeys`, `ReadObsm(key)`, `ToMemory()`, `Dispose`.
- `H5adMatrix` (on-disk X or layer, an `ISparseMatrixSource`): `Format`, `DataType`, `ReadRows(start, count)`, `ReadColumns(start, count)`.
- `ISparseMatrixSource` (`Rows`, `Columns`, `NonZeroCount`, `ReadRows`), implemented by `SparseMatrix`, `H5adMatrix` and `MappedSparseMatrix`; extension `RowBlocks(blockRows, maxEntries)`, `ToSparseMatrix()`; constants `SparseMatrixSource.DefaultBlockRows`, `DefaultBlockEntries`.
- `MappedSparseMatrix` (memory-mapped `.bnsparse` directory): `Create(directory, source, blockRows, progress)` (CSC h5ad sources are transposed out of core), `Open(directory, layer)`, `OpenLayer`, `WriteLayer(name, transform, blockRows)`, `LayerNames`, `RowPointer`, `RowIndices`, `RowValues`, `ReadRows`; `MappedSparseMatrixWriter.Create(directory, columns)` with `AppendRows`, `AppendRow`, `Complete`.
- `OutOfCore.CalculateQcMetrics(counts, geneNames, mitoPrefix, percentTop)` → `QcMetrics` (`WriteTo(obs, var)` with scanpy column names, `WriteTo(sc)`); `GeneMeanVariance(x, expm1)`; `NormalizeTotalLog1p(counts, outputDirectory, targetSum, log1p, cellTotals)`; `NormalizeTotalLog1pToLayer(store, layerName, …)`; `HighlyVariableGenes(x, nTopGenes, flavor, bins, span)` → `HvgResult` (`Means`, `Spread`, `Normalized`, `HighlyVariable`, `Rank`, `Selected`, `WriteTo`).
- `OutOfCorePca.Fit(x, genes, components, solver, scale, seed, powerIterations, oversamples)` → `PcaResult`; `OutOfCorePcaSolver.Auto`, `.Covariance` (exact), `.Randomized`.

`Preprocessing.HighlyVariableGenes` now shares the out-of-core implementation, so both match scanpy 1.12. Parity, memory and time figures are in the [tutorial](../tutorials/transcriptomics.md#5-atlas-scale-and-h5ad).

---
[API reference index](index.md) · [Transcriptomics tutorial](../tutorials/transcriptomics.md)
