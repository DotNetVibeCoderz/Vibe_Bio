# Referensi API: Transcriptomics

Paket `Gravicode.Biotech.Net.Transcriptomics` · namespace `Biotech.Net.Transcriptomics`, `.Enrichment`, `.SingleCell`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Matriks hitungan berupa objek `Dataset` (gen × sampel); metadata sampel ada di `Dataset.ColumnAnnotations`.

## Masukan dan normalisasi

- `CountMatrixReader.ReadFeatureCounts`, `ReadHtSeq`, `ReadHtSeqFiles`, `ReadTable`, `GeneLengths`, `WithGeneLengths`, `OpenText`.
- `Normalization.LibrarySizes`, `MedianOfRatiosSizeFactors` (DESeq2), `TmmFactors` (edgeR), `UpperQuartileFactors`, `Factors(method)`, `EffectiveLibrarySizes`, `Cpm`, `LogCpm`, `Tpm`, `Fpkm`, `Rpkm`, `NormalizedCounts`; enum `NormalizationMethod`.
- `Transforms.Vst`, `RlogApproximation`; `ExpressionFilter.FilterByExpr`, `MinimumCount`, `Apply`.

## Ekspresi diferensial

| Tipe | Anggota utama |
|---|---|
| `DESeq` | `Fit(counts, design, options)`, `Run(counts, factor, reference, treatment, options, covariates)` |
| `DESeqFit` | `SizeFactors`, `BaseMean`, `GeneDispersions`, `TrendDispersions`, `Dispersions`, `DispersionOutliers`, `Coefficients`, `StandardErrors`, `Converged`, `Results(coefficient, options)`, `VarianceStabilizingTransform()` |
| `DESeqOptions` | `FitType` (`Parametric` atau `Mean`), `ShrinkDispersions`, `OutlierSD`, `MaxIterations`, `IndependentFiltering`, `Alpha`, `SizeFactors` |
| `Limma` | `LmFit`, `EBayes`, `FitFDist`, `TopTable`, `Voom`, `VoomTest`, `ModeratedTTest` |
| `DifferentialExpression` | `TTest(logExpression, reference, treatment, welch)` plus kelebihan beban berbasis anotasi |
| `DesignMatrix` | `Intercept`, `TwoGroup`, `FromAnnotations(data, factor, reference, covariates)`, `IndexOf`, `Rank()`, `Leverages()` |
| `DifferentialExpressionResult` | `Genes`, `BaseMean`, `Log2FoldChange`, `LfcSE`, `Stat`, `PValue`, `PAdj`, `Rows()`, `Top(n)`, `Significant(alpha, minAbsLog2FoldChange)`, `CountSignificant`, `Volcano()`, `MaPlot()`, `ToDataset()`, `WriteTsv` |

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

---
[Indeks referensi API](index.md) · [Tutorial transkriptomika](../tutorials/transcriptomics.md)
