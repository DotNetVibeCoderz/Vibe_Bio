# Referensi API: Transcriptomics

Paket `Gravicode.Biotech.Net.Transcriptomics` · namespace `Biotech.Net.Transcriptomics`, `.Enrichment`, `.SingleCell`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Matriks hitungan berupa objek `Dataset` (gen × sampel); metadata sampel ada di `Dataset.ColumnAnnotations`.

## Masukan dan normalisasi

- `CountMatrixReader.ReadFeatureCounts`, `ReadHtSeq`, `ReadHtSeqFiles`, `ReadTable`, `GeneLengths`, `WithGeneLengths`, `OpenText`.
- `Normalization.LibrarySizes`, `MedianOfRatiosSizeFactors` (DESeq2 "ratio"), `PosCountsSizeFactors` (DESeq2 "poscounts"), `TmmFactors` (edgeR), `UpperQuartileFactors`, `Factors(method)`, `EffectiveLibrarySizes`, `Cpm`, `LogCpm`, `Tpm`, `Fpkm`, `Rpkm`, `NormalizedCounts`; enum `NormalizationMethod`.
- `Transforms.Vst`, `RlogApproximation`; `ExpressionFilter.FilterByExpr`, `MinimumCount`, `Apply`.

## Ekspresi diferensial

| Tipe | Anggota utama |
|---|---|
| `DESeq` | `Fit(counts, design, options)` (Wald, atau LRT bila `options.ReducedDesign` diisi), `FitLrt(counts, full, reduced, options)`, `Run(counts, factor, reference, treatment, options, covariates)` |
| `DESeqFit` | `SizeFactors`, `BaseMean`, `BaseVar`, `GeneDispersions`, `TrendDispersions`, `MapDispersions`, `Dispersions`, `DispersionOutliers`, `Trend`, `TrendType`, `TrendFunction`, `VarLogDispersionEstimates`, `DispersionPriorVariance`, `Coefficients`, `StandardErrors`, `Converged`, `BetaIterations`, `Deviance`, `Mu`, `HatDiagonals`, `CooksDistances`, `MaxCooks`, `Replaced`, `ReplacementCounts`, `LrtStatistic`, `LrtPValue`, `Test`; `Results(coefficient, options)`, `Results(name, options)`, `Results(factor, numerator, denominator, options)`, `Results(double[] contrast, options)`, `VarianceStabilizingTransform()` |
| `DESeqOptions` | `SizeFactorType` (`Ratio`, `PosCounts`), `FitType` (`Parametric`, `Local`, `Mean`), `ShrinkDispersions`, `OutlierSD`, `MinDispersion`, `MaxIterations`, `Lambda`, `BetaTolerance`, `MinMu`, `UseOptim`, `DispersionPriorVariance`, `ReducedDesign`, `MinReplicatesForReplace` (7), `CooksFiltering`, `CooksCutoff`, `IndependentFiltering`, `Alpha`, `SizeFactors` |
| `LfcShrink` | `Apeglm(fit, coefficient)`, `Normal(fit, coefficient)`, `Shrink(fit, coefficient, LfcShrinkType)` yang mengembalikan `LfcShrinkResult` (`Table`, `PriorVariance`, `PriorScale`, `SValue`, `FalseSignRate`, `Convergence`) |
| `Limma` | `LmFit`, `EBayes`, `FitFDist`, `TopTable`, `Voom`, `VoomTest`, `ModeratedTTest` |
| `DifferentialExpression` | `TTest(logExpression, reference, treatment, welch)` plus kelebihan beban berbasis anotasi |
| `DesignMatrix` | `Intercept`, `TwoGroup`, `FromAnnotations(data, factor, reference, covariates)`, `IndexOf`, `Rank()`, `Leverages()`, `Factors` (`DesignFactor`: name, levels, columns) |
| `DifferentialExpressionResult` | `Genes`, `BaseMean`, `Log2FoldChange`, `LfcSE`, `Stat`, `PValue`, `PAdj`, `Rows()`, `Top(n)`, `Significant(alpha, minAbsLog2FoldChange)`, `CountSignificant`, `Volcano()`, `MaPlot()`, `ToDataset()`, `WriteTsv` |

`DESeq` dan `LfcShrink` adalah port DESeq2 1.52.0 dan apeglm 1.34.0 yang diuji terhadap R pada data pasilla dan airway: tahap berbentuk tertutup sama hingga ~1e-12 relatif, tahap iteratif hingga 1e-6 untuk sedikitnya 99% gen (1e-3 untuk sisanya, akibat titik berhenti pencarian garis). Lihat [tutorial](../tutorials/transcriptomics.md) untuk rincian dan perbedaan yang tersisa.

## Matriks dan PCA

`DenseMatrix(rows, columns)` dengan `FromDataset(dataset, transpose)`, `ToDataset`, `Row`, `Column`, `Transpose`; `Pca.Fit(data, components, solver, seed, powerIterations, center)` mengembalikan `PcaResult(Scores, Loadings, ExplainedVariance, ExplainedVarianceRatio, Means)`.

## Pengayaan (`.Enrichment`)

- `GeneSet(name, genes, description)`; `Gmt.Read`, `ReadFile`, `Write`.
- `OverRepresentation.Test(hits, universe, sets, minSetSize, maxSetSize)` mengembalikan `OverRepresentationResult(Set, SetSize, Overlap, Expected, FoldEnrichment, OddsRatio, PValue, PAdj, Genes)`; `HypergeometricUpperTail`.
- `Gsea.Run(ranking, sets, permutations, weight, minSize, maxSize, seed)` mengembalikan `GseaResult(Set, Size, EnrichmentScore, NormalizedEnrichmentScore, PValue, PAdj, FdrQ, LeadingEdge)`; `EnrichmentScore`.

## Sel tunggal (`.SingleCell`)

- `SparseMatrix` (CSR): `FromTriplets`, `FromDense`, `ToDense`, `Transpose`, `SelectRows`, `SelectColumns`, `RowSums`, `ColumnSums`, `ColumnMeanVariance`, `Density`.
- `SingleCellExperiment(x, cellNames, geneIds, geneNames)`: `X`, `Layers`, `CellMetrics`, `CellLabels`, `GeneMetrics`, `GeneFlags`, `Embeddings`, `Neighbors`, `GeneIndex`, `SubsetCells`, `SubsetGenes`, `Copy`.
- `MatrixMarket.Read`, `Write`, `ReadSingleCell`, `Read10x(directory)`.
- `Preprocessing.CalculateQcMetrics`, `FilterCells`, `FilterGenes`, `NormalizeTotal`, `Log1p`, `HighlyVariableGenes(nTopGenes, flavor)`, `Scale`.
- `NeighborGraph.Build(data atau sc, nNeighbors, metric atau embedding)`: `Indices`, `Distances`, `Connectivities`.
- `CommunityDetection.Cluster(adjacency atau sc, resolution, algorithm, seed, key)`, `Modularity`; `CommunityAlgorithm.Leiden`, `.Louvain`.
- `Umap.Fit(graph atau sc, options)` dengan `UmapOptions` dan `FindAbParams`; `Tsne.Fit(data atau sc, options)` dengan `TsneOptions`.
- `MarkerGenes.Rank(sc, labels atau key, test)` mengembalikan `MarkerGene(Gene, GeneIndex, Score, Log2FoldChange, PValue, PAdj, PctInGroup, PctRest)`; enum `MarkerTest`.

Kunci embedding mengikuti scanpy: `X_pca`, `X_umap`, `X_tsne`; label klaster disimpan di `CellLabels["leiden"]` atau `["louvain"]`.

## AnnData / h5ad dan skala atlas (`.SingleCell`)

Akses HDF5 sepenuhnya kode terkelola ([PureHDF](https://github.com/Apollo3zehn/PureHDF) 2.2.0, MIT): tidak perlu pustaka HDF5 native.

- `H5ad.Read(path)` → `AnnData`; `H5ad.ReadSingleCell(path, geneSymbolColumn)`; `H5ad.Write(path, AnnData atau SingleCellExperiment, options)`; `H5ad.OpenBacked(path)` → `BackedAnnData`; `H5ad.WriteStreaming(path, x, obs, var, layers, obsm, uns, options, blockRows)` menulis X dan layer dari sumber out-of-core blok demi blok.
- `H5adWriteOptions`: `Compression` (`H5adCompression.None` atau `.Gzip` = shuffle + deflate), `GzipLevel` (1 atau 9), `Shuffle`, `ChunkElements`, `CompressionThreshold`.
- `AnnData(x, obs, var)`: `X`, `Obs`, `Var`, `Obsm`, `Varm`, `Obsp`, `Varp`, `Layers`, `Uns`, `Raw`, `ToSingleCellExperiment(geneSymbolColumn)`, `FromSingleCellExperiment(sc)`.
- `AnnMatrix`: `FromSparse(matrix, format, dataType, values64)`, `FromDense(matrix, dataType)`, `Format` (`AnnMatrixFormat.Dense`, `.Csr`, `.Csc`), `DataType` (nama dtype numpy), `Sparse`, `Dense`, `SparseValues64`, `ToSparse`, `ToDense`, indexer.
- `AnnDataFrame(index, indexName)`: `Index`, `IndexName`, `ColumnNames`, indexer, `Set`, `Remove`, `Contains`, `GetDoubles`, `GetStrings`, `SelectRows`. Kolom berupa array primitif, `string[]`, `Categorical(codes, categories, ordered)` (`FromLabels`, `ToStrings`) atau `MaskedArray(values, mask)` (integer/boolean/string nullable, `NaValue`).
- `AnnArray(data, shape)` (array N-dimensi di `uns`), `AnnRaw(x, var)` dengan `Varm`, `AnnUnsupported(Path, EncodingType, Reason)`.
- `BackedAnnData`: `X`, `Obs`, `Var`, `Uns`, `LayerNames`, `Layer(name)`, `ObsmKeys`, `ReadObsm(key)`, `ToMemory()`, `Dispose`.
- `H5adMatrix` (X atau layer di disk, sebuah `ISparseMatrixSource`): `Format`, `DataType`, `ReadRows(start, count)`, `ReadColumns(start, count)`.
- `ISparseMatrixSource` (`Rows`, `Columns`, `NonZeroCount`, `ReadRows`), diimplementasikan oleh `SparseMatrix`, `H5adMatrix` dan `MappedSparseMatrix`; ekstensi `RowBlocks(blockRows, maxEntries)`, `ToSparseMatrix()`; konstanta `SparseMatrixSource.DefaultBlockRows`, `DefaultBlockEntries`.
- `MappedSparseMatrix` (direktori `.bnsparse` yang di-memory-map): `Create(directory, source, blockRows, progress)` (sumber h5ad CSC ditransposisi secara out-of-core), `Open(directory, layer)`, `OpenLayer`, `WriteLayer(name, transform, blockRows)`, `LayerNames`, `RowPointer`, `RowIndices`, `RowValues`, `ReadRows`; `MappedSparseMatrixWriter.Create(directory, columns)` dengan `AppendRows`, `AppendRow`, `Complete`.
- `OutOfCore.CalculateQcMetrics(counts, geneNames, mitoPrefix, percentTop)` → `QcMetrics` (`WriteTo(obs, var)` dengan nama kolom scanpy, `WriteTo(sc)`); `GeneMeanVariance(x, expm1)`; `NormalizeTotalLog1p(counts, outputDirectory, targetSum, log1p, cellTotals)`; `NormalizeTotalLog1pToLayer(store, layerName, …)`; `HighlyVariableGenes(x, nTopGenes, flavor, bins, span)` → `HvgResult` (`Means`, `Spread`, `Normalized`, `HighlyVariable`, `Rank`, `Selected`, `WriteTo`).
- `OutOfCorePca.Fit(x, genes, components, solver, scale, seed, powerIterations, oversamples)` → `PcaResult`; `OutOfCorePcaSolver.Auto`, `.Covariance` (eksak), `.Randomized`.

`Preprocessing.HighlyVariableGenes` kini memakai implementasi out-of-core yang sama, sehingga keduanya cocok dengan scanpy 1.12. Angka kesetaraan, memori, dan waktu ada di [tutorial](../tutorials/transcriptomics.md#5-skala-atlas-dan-h5ad).

---
[Indeks referensi API](index.md) · [Tutorial transkriptomika](../tutorials/transcriptomics.md)
