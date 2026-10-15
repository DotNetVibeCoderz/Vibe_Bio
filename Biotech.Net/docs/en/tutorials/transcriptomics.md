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

## Limitations

- `DESeq` follows the published DESeq2 method (median-of-ratios size factors, a parametric or mean dispersion trend, empirical-Bayes dispersion shrinkage with dispersion-outlier handling, a Wald test and independent filtering), but exact numerical parity with R DESeq2 on public datasets is still a roadmap item, and apeglm LFC shrinkage is not implemented.
- The single-cell pipeline keeps matrices in memory as CSR `SparseMatrix`. Atlas-scale data and an AnnData/h5ad reader are on the [roadmap](../roadmap.md).
- The neighbour graph is found by exact brute-force search (O(n²)), and t-SNE is the exact O(n²) method (Barnes–Hut is not implemented). Both are meant for up to a few thousand cells, not millions.

---
Previous: [Genomics](genomics.md) · Next: [Proteomics](proteomics.md) · [Tutorial index](index.md)
