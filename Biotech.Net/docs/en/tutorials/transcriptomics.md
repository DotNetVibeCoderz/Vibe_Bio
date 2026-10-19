# Tutorial: Transcriptomics

Package: `Gravicode.Biotech.Net.Transcriptomics`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Bulk RNA-seq (normalisation, DESeq2-style and limma-voom-style differential expression, VST, PCA), enrichment (ORA, preranked GSEA) and single-cell analysis (QC, HVG, PCA, neighbour graph, Leiden/Louvain, UMAP, t-SNE, marker genes).

Count matrices are `Dataset` objects (genes in rows, samples in columns). Sample metadata lives in `Dataset.ColumnAnnotations`.

## 1. Load counts

```csharp
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

// From files:
//   CountMatrixReader.ReadFeatureCounts(reader)     featureCounts output (keeps gene lengths)
//   CountMatrixReader.ReadHtSeqFiles(paths)         one HTSeq-count file per sample
//   CountMatrixReader.ReadTable(reader)             any genes x samples table

// A synthetic experiment: 1,000 genes, 4 control + 4 treated samples, 60 true changes
var rng = new Random(42);
int genes = 1000, samples = 8;
var values = new double[genes * samples];
for (int g = 0; g < genes; g++)
{
    double mean = Math.Exp(2 + rng.NextDouble() * 6);
    double fold = g < 30 ? 4 : g < 60 ? 0.25 : 1;
    for (int s = 0; s < samples; s++)
    {
        double mu = mean * (s >= 4 ? fold : 1);
        values[g * samples + s] = Math.Round(mu * Math.Exp(0.25 * (rng.NextDouble() + rng.NextDouble() + rng.NextDouble() - 1.5)));
    }
}
var counts = new Dataset("counts",
    Enumerable.Range(1, genes).Select(i => $"gene{i:D4}").ToArray(),
    ["ctl1", "ctl2", "ctl3", "ctl4", "trt1", "trt2", "trt3", "trt4"], values);
counts.Annotate("condition", "ctl", "ctl", "ctl", "ctl", "trt", "trt", "trt", "trt");

// Normalisation
double[] sizeFactors = Normalization.MedianOfRatiosSizeFactors(counts);  // DESeq2
double[] tmm = Normalization.TmmFactors(counts);                          // edgeR
Dataset logCpm = Normalization.LogCpm(counts, normFactors: tmm);
Console.WriteLine(string.Join(", ", sizeFactors.Select(f => f.ToString("F3"))));
```

## 2. Differential expression

```csharp
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

var rng = new Random(1);
var values = new double[200 * 6];
for (int g = 0; g < 200; g++)
    for (int s = 0; s < 6; s++)
        values[g * 6 + s] = Math.Round(50 * (g < 20 && s >= 3 ? 3 : 1) * (0.8 + 0.4 * rng.NextDouble()));
var counts = new Dataset("counts", Enumerable.Range(0, 200).Select(i => $"g{i}").ToArray(),
    ["a1", "a2", "a3", "b1", "b2", "b3"], values);
counts.Annotate("group", "A", "A", "A", "B", "B", "B");

// DESeq2-style: negative-binomial GLM, dispersion shrinkage, Wald test, BH adjustment
DifferentialExpressionResult deseq = DESeq.Run(counts, "group", "A", "B");
var (up, down) = deseq.CountSignificant(alpha: 0.05);
Console.WriteLine($"DESeq: {up} up, {down} down");
foreach (var row in deseq.Top(3))
    Console.WriteLine($"{row.Gene}  log2FC {row.Log2FoldChange:F2}  padj {row.PAdj:E2}");

// limma-voom-style
var design = DesignMatrix.FromAnnotations(counts, "group", "A");
DifferentialExpressionResult voom = Limma.VoomTest(counts, design);

// Variance-stabilised values and PCA for QC plots
Dataset vst = Transforms.Vst(counts);
PcaResult pca = Pca.Fit(DenseMatrix.FromDataset(vst, transpose: true), components: 2);
Console.WriteLine($"PC1 explains {pca.ExplainedVarianceRatio[0]:P1}");

// Export the result table
var tsv = new StringWriter();
deseq.WriteTsv(tsv);
```


### DESeq2 in depth: designs, contrasts, LRT and shrinkage

`DESeq` is a step-by-step port of DESeq2 1.52: size factors (`Ratio` or `PosCounts`), gene-wise dispersions by DESeq2's
Cox–Reid line search, a parametric / local / mean trend, the dispersion prior variance (including DESeq2's seeded
simulation when there are 3 or fewer residual degrees of freedom), MAP shrinkage with dispersion outliers, the ridge
IRLS GLM, Wald tests or a likelihood-ratio test, Cook's distance filtering and outlier replacement
(`MinReplicatesForReplace`, default 7), and independent filtering on the mean of normalised counts. `LfcShrink` adds
apeglm (Cauchy prior, MAP, Laplace posterior SD, s-values) and DESeq2's original normal prior.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

// 500 genes, 4 control + 4 treated samples sequenced in two batches; genes 0-39 change 3-fold.
var rng = new Random(7);
int genes = 500;
var values = new double[genes * 8];
for (int g = 0; g < genes; g++)
{
    double mean = Math.Exp(3 + 5 * rng.NextDouble());
    for (int s = 0; s < 8; s++)
    {
        double mu = mean * (g < 40 && s >= 4 ? 3 : 1) * (s % 4 >= 2 ? 1.3 : 1);
        double gamma = -Math.Log(rng.NextDouble()) * 0.1 + 0.9;            // mild extra-Poisson noise
        values[g * 8 + s] = Math.Round(mu * gamma);
    }
}
var counts = new Dataset("counts", Enumerable.Range(0, genes).Select(i => $"g{i}").ToArray(),
    ["c1", "c2", "c3", "c4", "t1", "t2", "t3", "t4"], values);
counts.Annotate("condition", "ctl", "ctl", "ctl", "ctl", "trt", "trt", "trt", "trt");
counts.Annotate("batch", "b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2");

// ~ batch + condition
var design = DesignMatrix.FromAnnotations(counts, "condition", "ctl", "batch");
DESeqFit fit = DESeq.Fit(counts, design);
Console.WriteLine(string.Join(", ", design.CoefficientNames));   // Intercept, batch_b2, condition_trt_vs_ctl
Console.WriteLine($"trend {fit.TrendType}: {fit.Trend.AsymptoticDispersion:F4} + {fit.Trend.ExtraPoisson:F3}/mean, prior var {fit.DispersionPriorVariance:F3}");

DifferentialExpressionResult res = fit.Results("condition_trt_vs_ctl");     // Wald test of a coefficient
DifferentialExpressionResult flipped = fit.Results("condition", "ctl", "trt"); // level contrast (sign flipped)
DifferentialExpressionResult custom = fit.Results([0, 1, -1]);               // numeric contrast c'β
Console.WriteLine(res);

// Likelihood-ratio test of the condition effect: full ~ batch + condition vs reduced ~ batch
var reduced = DesignMatrix.FromAnnotations(counts, "batch", "b1");
DifferentialExpressionResult lrt = DESeq.FitLrt(counts, design, reduced).Results();
Console.WriteLine($"LRT: {lrt.CountSignificant(0.1)}");

// Fold-change shrinkage for ranking and plots
LfcShrinkResult ape = LfcShrink.Apeglm(fit, "condition_trt_vs_ctl");
LfcShrinkResult normal = LfcShrink.Normal(fit, "condition_trt_vs_ctl");
Console.WriteLine($"apeglm prior scale {ape.PriorScale:F3}; g0 MLE {res.Log2FoldChange[0]:F3}, apeglm {ape.Table.Log2FoldChange[0]:F3} " +
                  $"(s-value {ape.SValue![0]:E1}), normal {normal.Table.Log2FoldChange[0]:F3}");

// Options mirror DESeq() / results() arguments
var options = new DESeqOptions
{
    SizeFactorType = SizeFactorType.PosCounts,   // sfType = "poscounts"
    FitType = DispersionFitType.Local,           // fitType = "local"
    MinReplicatesForReplace = 7,                 // Cook's outlier replacement (default 7)
    CooksFiltering = true,                       // cooksCutoff (default qf(0.99, p, m - p))
    IndependentFiltering = true, Alpha = 0.05,   // independentFiltering, alpha
};
var res2 = DESeq.Fit(counts, design, options).Results(options: options);
Console.WriteLine(res2);
_ = flipped; _ = custom;
```

**Parity with R.** The test suite compares every stage with R 4.6.1 / DESeq2 1.52.0 / apeglm 1.34.0 on the public
`pasilla` (7 samples) and `airway` (8 samples) data sets. Closed-form steps agree to ~1e-12 relative. Iterative steps
(dispersion line search, IRLS, apeglm's L-BFGS) agree to 1e-6 relative for at least 99% of genes and to 1e-3 for the
rest: DESeq2's line search stops after a log-posterior gain below 1e-6, and on the flat posteriors of high-count genes
the last-bit differences between R and .NET arithmetic move that stopping point. On the full pasilla data (14,599
genes, `~ condition`) the largest relative differences are 7.8e-6 for dispersions, 1.1e-5 for p-values and padj, and
the independent-filtering threshold, the set of filtered genes and the Cook's outliers are identical. Genes whose IRLS
does not converge (5 of 12,359 tested genes for `~ type + condition`) are refitted by R with `optim(L-BFGS-B)` and here
with a projected Newton method; their fold changes differ by up to 2e-4 (log2). The `Local` trend approximates locfit
(median relative difference ~0.5%). On the full airway data (63,677 genes, `~ cell + dex`, 3 residual degrees of
freedom) the simulated dispersion prior variance is identical to R's; apeglm fold changes agree to 1e-5 except for
0.25% of genes - low counts at the maximum dispersion, where the posterior is flat and apeglm's L-BFGS stopping point
moves by up to 1.4e-4 (log2), about as far as apeglm's own two starting points end apart.

`Biotech.Net.Visualization` turns a result into a volcano or MA plot (`new VolcanoPlot(deseq)`, `new MaPlot(deseq)`):

![Volcano plot](../../images/volcano.png)

## 3. Enrichment

```csharp
using Biotech.Net.Transcriptomics.Enrichment;

var sets = new List<GeneSet>
{
    new("CELL_CYCLE", ["CDK1", "CCNB1", "MKI67", "TOP2A", "BUB1"]),
    new("INTERFERON", ["ISG15", "IFIT1", "MX1", "OAS1", "STAT1"]),
};
var universe = sets.SelectMany(s => s.Genes).Concat(Enumerable.Range(0, 200).Select(i => $"BG{i}")).ToList();
var hits = new[] { "CDK1", "CCNB1", "MKI67", "TOP2A", "BG1", "BG2" };

foreach (var r in OverRepresentation.Test(hits, universe, sets))
    Console.WriteLine($"{r.Set}: overlap {r.Overlap}/{r.SetSize}, p {r.PValue:E2}, padj {r.PAdj:E2}");

// Preranked GSEA takes gene -> score (for example log2FC or the Wald statistic)
var ranking = universe.Select((g, i) => (g, score: g.StartsWith("BG") ? (i % 7) - 3.0 : 5.0))
                      .ToDictionary(x => x.g, x => x.score);
var gsea = Gsea.Run(ranking, sets, permutations: 500, minSize: 3);
```

Gene sets can be read from MSigDB `.gmt` files with `Gmt.ReadFile(path)`.

## 4. Single-cell RNA-seq

```csharp
using Biotech.Net.Transcriptomics;
using Biotech.Net.Transcriptomics.SingleCell;

// From disk: MatrixMarket.Read10x("filtered_feature_bc_matrix/")
// Here: 300 cells, 3 populations, 100 genes
var rng = new Random(7);
var rows = new List<int>(); var cols = new List<int>(); var vals = new List<float>();
for (int c = 0; c < 300; c++)
    for (int g = 0; g < 100; g++)
    {
        int population = c / 100;
        double rate = (g / 10 == population) ? 8 : 1;   // genes 0-9, 10-19, 20-29 mark the populations
        int count = (int)Math.Round(rate * rng.NextDouble() * 2);
        if (count > 0) { rows.Add(c); cols.Add(g); vals.Add(count); }
    }
var x = SparseMatrix.FromTriplets(300, 100, rows.ToArray(), cols.ToArray(), vals.ToArray());
var sc = new SingleCellExperiment(x,
    Enumerable.Range(0, 300).Select(i => $"cell{i}").ToArray(),
    Enumerable.Range(0, 100).Select(i => $"GENE{i}").ToArray());

Preprocessing.CalculateQcMetrics(sc);
Preprocessing.FilterGenes(sc, minCells: 3);
Preprocessing.NormalizeTotal(sc);
Preprocessing.Log1p(sc);
int[] hvg = Preprocessing.HighlyVariableGenes(sc, nTopGenes: 50);

sc.Embeddings["X_pca"] = Pca.Fit(Preprocessing.Scale(sc, hvg), components: 10).Scores;
NeighborGraph.Build(sc, nNeighbors: 15);
CommunityDetection.Cluster(sc, resolution: 0.5);          // stores CellLabels["leiden"]
Umap.Fit(sc, new UmapOptions { Seed = 1 });                // stores Embeddings["X_umap"]

var markers = MarkerGenes.Rank(sc, "leiden");
foreach (var (cluster, list) in markers.OrderBy(kv => kv.Key))
    Console.WriteLine($"cluster {cluster}: {string.Join(", ", list.Take(3).Select(m => m.Gene))}");
```

`EmbeddingPlot.FromSingleCell(sc, "X_umap", colorBy: "leiden")` plots the embedding.

## 5. Atlas scale and h5ad

AnnData `.h5ad` files are read and written with pure managed code (PureHDF 2.2.0, no native HDF5). The reader
follows the AnnData on-disk encodings (dense arrays, `csr_matrix`/`csc_matrix`, dataframes with categorical, string,
nullable integer/boolean/string columns, `obsm`/`varm`/`obsp`/`varp`, `layers`, nested `uns`, `raw`, and the pre-0.7
`h5sparse` layout). For data that does not fit in memory, every preprocessing step streams over row blocks of an
`ISparseMatrixSource`: an h5ad file opened in backed mode, or a memory-mapped `.bnsparse` store
(`MappedSparseMatrix`). Results are ordinary arrays and a `PcaResult`, so the existing neighbour graph, Leiden/Louvain,
UMAP and marker-gene code takes over from there.

```csharp
using Biotech.Net.Transcriptomics;
using Biotech.Net.Transcriptomics.SingleCell;

// A small stand-in for a real atlas file: 3 000 cells x 2 000 genes, written as .h5ad.
var rng = new Random(7);
var rows = new List<int>(); var cols = new List<int>(); var vals = new List<float>();
for (int c = 0; c < 3000; c++)
    for (int g = 0; g < 2000; g++)
    {
        double rate = (g / 20 == c % 5) ? 3.0 : 0.15;     // genes 0-19, 20-39, ... mark five populations
        int count = rng.NextDouble() < rate / (1 + rate) ? 1 + (int)(rng.NextDouble() * rate * 2) : 0;
        if (count > 0) { rows.Add(c); cols.Add(g); vals.Add(count); }
    }
var x = SparseMatrix.FromTriplets(3000, 2000, rows.ToArray(), cols.ToArray(), vals.ToArray());
var cellTable = new AnnDataFrame(Enumerable.Range(0, 3000).Select(i => $"cell{i}").ToArray());
cellTable["sample"] = Categorical.FromLabels(Enumerable.Range(0, 3000).Select(i => i < 1500 ? "A" : "B").ToArray());
var geneTable = new AnnDataFrame(Enumerable.Range(0, 2000).Select(j => j < 5 ? $"MT-{j}" : $"GENE{j}").ToArray());
var adata = new AnnData(AnnMatrix.FromSparse(x), cellTable, geneTable);
adata.Uns["description"] = "synthetic example";
H5ad.Write("atlas.h5ad", adata, new H5adWriteOptions { Compression = H5adCompression.Gzip });

// 1. Read it back in memory (works for any anndata-written file) and hand it to the in-memory pipeline.
AnnData loaded = H5ad.Read("atlas.h5ad");
Console.WriteLine(loaded);
SingleCellExperiment small = loaded.ToSingleCellExperiment();

// 2. Backed mode: obs/var in memory, X stays on disk and is read in slices.
using var backed = H5ad.OpenBacked("atlas.h5ad");
SparseMatrix firstCells = backed.X!.ReadRows(0, 100);
Console.WriteLine($"{backed}; first 100 cells hold {firstCells.NonZeroCount} counts");

// 3. Out of core: copy X into a memory-mapped store, then stream every step over it.
using var counts = MappedSparseMatrix.Create("atlas.bnsparse", backed.X);
QcMetrics qc = OutOfCore.CalculateQcMetrics(counts, backed.Var.Index);              // mito prefix "MT-"
HvgResult hvg = OutOfCore.HighlyVariableGenes(counts, nTopGenes: 500, flavor: HvgFlavor.SeuratV3);
using var logNorm = OutOfCore.NormalizeTotalLog1pToLayer(counts, "log1p_norm", cellTotals: qc.TotalCounts);
PcaResult pca = OutOfCorePca.Fit(logNorm, hvg.Selected, components: 30);
Console.WriteLine($"PC1 explains {pca.ExplainedVarianceRatio[0]:P1}; {hvg.Selected.Length} HVGs");

// 4. The PCA scores feed the existing graph, clustering and UMAP code.
var graph = NeighborGraph.Build(pca.Scores, nNeighbors: 15);
int[] clusters = CommunityDetection.Cluster(graph.Connectivities, resolution: 0.5);
var umap = Umap.Fit(graph, new UmapOptions { Seed = 1 });
Console.WriteLine($"{clusters.Distinct().Count()} Leiden clusters");

// 5. Write the result as .h5ad without loading the matrix: log-normalised X, raw counts layer, PCA, UMAP, QC.
qc.WriteTo(backed.Obs, backed.Var);
hvg.WriteTo(backed.Var);
backed.Obs["leiden"] = Categorical.FromLabels(clusters.Select(c => c.ToString()).ToArray());
H5ad.WriteStreaming("atlas_processed.h5ad", logNorm, backed.Obs, backed.Var,
    layers: new Dictionary<string, ISparseMatrixSource> { ["counts"] = counts },
    obsm: new Dictionary<string, AnnMatrix> { ["X_pca"] = AnnMatrix.FromDense(pca.Scores), ["X_umap"] = AnnMatrix.FromDense(umap) },
    options: new H5adWriteOptions { Compression = H5adCompression.Gzip });
```

`OutOfCore.NormalizeTotalLog1p(source, directory)` writes a new store instead of a layer (use it for h5ad-backed
sources); `MappedSparseMatrix.Create` transposes CSC files out of core, because row blocks of a CSC matrix would
otherwise need a full scan each.

### Parity with scanpy

Measured against anndata 0.13.4 / scanpy 1.12.4 / scikit-learn 1.9.1 on a 300-cell × 1 000-gene subset of 10x
pbmc3k, plus synthetic files covering every encoding (`tests/Biotech.Net.Transcriptomics.Tests/Fixtures/AnnData`,
generator `make_fixtures.py`):

| Step | Agreement |
|---|---|
| Reading (dense, CSR gzip, CSC, categoricals, nullable, `uns`, `raw`, legacy) | exact, compared with values dumped by numpy |
| Writing | anndata reads our files back identical to the source (`X`, `obs`/`var` dtypes, `obsm`, `obsp`, `layers`, `uns`, `raw`) |
| `calculate_qc_metrics` (all obs/var columns, `percent_top`) | ≤ 1e-6 relative (scanpy sums in float32) |
| `normalize_total` + `log1p` (1e4 and median target) | ≤ 1e-6 relative, row and column sums and 200 sampled entries |
| `highly_variable_genes(flavor="seurat_v3")` | means 1e-12, variances 1e-10, `variances_norm` 1e-9; same genes and ranks (ties broken differently: numpy's argsort is unstable) |
| `highly_variable_genes(flavor="seurat")` | means/dispersions 1e-6, `dispersions_norm` 1e-5; same genes |
| PCA, `OutOfCorePcaSolver.Covariance` | explained variance 1e-8 relative to scikit-learn; loadings \|cos\| > 1 − 1e-8; scores 1e-6 (up to sign) |
| PCA, `OutOfCorePcaSolver.Randomized` (6 power iterations, 20 oversamples) | variance of PC1–10 within 0.1 %, PC1–20 within 2 %; PC1–5 \|cos\| > 0.999 |

`Preprocessing.HighlyVariableGenes` now uses the same code, so the in-memory pipeline gained the same parity.

### Memory and time at one million cells

A synthetic 1 000 000 cells × 30 000 genes matrix (about 1 100 detected genes per cell, 1.11 × 10⁹ stored entries,
8.9 GB on disk) on a laptop (Intel i7-8650U, 4 cores / 8 threads, 16 GB RAM, local SSD), Release build, generated and
run by `dotnet run -c Release -p:SkipRustBuild=true -- --atlas` in `benchmarks/Biotech.Net.Benchmarks`:

| Step (memory-mapped store) | Time | Peak private memory |
|---|---|---|
| QC metrics (1 pass) | 34 s | 0.64 GB |
| HVG Seurat v3, 2 000 genes (2 passes) | 54 s | 0.82 GB |
| `normalize_total` + `log1p` written as a new layer (4.4 GB) | 71 s | 0.39 GB |
| PCA, covariance solver, 2 000 HVGs, 50 PCs (2 passes) | 102 s | 2.8 GB (of which 0.4 GB are the 1 M × 50 scores) |
| Streamed gzip h5ad export of 200 000 cells (0.45 GB file) | 18 s | 1.8 GB |
| QC from that h5ad in backed mode | 36 s | 2.1 GB |

Private memory stays the same from 100 000 to 1 000 000 cells (0.63 vs 0.64 GB for QC): it is bounded by the
block size (≤ 2 M entries per block, a few blocks in flight) plus O(cells + genes) result arrays, and for PCA the
genes² covariance and the cells × components scores. Peak private memory includes garbage the GC had not yet
collected; the live heap measured after forced collections during the h5ad export was 0.09 GB. The working set
grows up to 12 GB because the OS keeps mapped file pages resident; those pages are clean file cache that it reclaims
under memory pressure, not process allocations.

### Limitations at atlas scale

- The neighbour graph is still exact brute force: 49 s for 50 000 cells, so roughly hours at a million. Build it on a
  subsample or on PCA scores of a subset; an approximate nearest-neighbour index is not implemented yet. t-SNE is
  likewise exact O(n²).
- Marker genes, `Scale` and the dense in-memory `Pca` need a `SingleCellExperiment` in memory.
- PureHDF cannot modify an existing HDF5 file, so new layers go to a `.bnsparse` store or to a new h5ad written by
  `H5ad.WriteStreaming`; there is no in-place append to an h5ad.
- Backed reads through PureHDF are about 3–4× slower than the memory-mapped store (gzip decompression and hyperslab
  bookkeeping); convert once with `MappedSparseMatrix.Create` when several passes are needed.
- Sparse values are float32 in memory (float64 sparse elements keep an exact copy only for round-tripping); dense
  matrices are float64. Awkward arrays, compound/record arrays (older `rank_genes_groups`), the pre-0.7 compound
  `obs`/`var` layout, zarr and float16 writing are not supported; unsupported `uns` entries are read as
  `AnnUnsupported` placeholders and skipped on write.
- Gzip levels are limited to 1 and 9 by PureHDF; LZF and Blosc filters are not supported.

## Limitations

- `DESeq` reproduces DESeq2 1.52 (see "Parity with R" above) but does not support observation weights, normalisation-factor matrices, `betaPrior = TRUE`, the glmGamPoi backend or interaction-aware level contrasts; `LfcShrink` offers apeglm and the normal prior, not ashr.
- The in-memory single-cell pipeline keeps matrices as CSR `SparseMatrix`; for larger data use the out-of-core steps and h5ad I/O of section 5, which has its own limitations list.
- The neighbour graph is found by exact brute-force search (O(n²)), and t-SNE is the exact O(n²) method (Barnes–Hut is not implemented). Both are meant for up to tens of thousands of cells, not millions.

---
Previous: [Genomics](genomics.md) · Next: [Proteomics](proteomics.md) · [Tutorial index](index.md)
