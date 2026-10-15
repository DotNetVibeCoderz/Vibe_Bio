# Tutorial: Transkriptomika

Paket: `Gravicode.Biotech.Net.Transcriptomics`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

RNA-seq bulk (normalisasi, ekspresi diferensial gaya DESeq2 dan limma-voom, VST, PCA), analisis pengayaan (ORA, GSEA preranked), dan analisis sel tunggal (QC, HVG, PCA, graf tetangga, Leiden/Louvain, UMAP, t-SNE, gen penanda).

Matriks hitungan berupa objek `Dataset` (gen di baris, sampel di kolom). Metadata sampel disimpan di `Dataset.ColumnAnnotations`.

## 1. Memuat hitungan

```csharp
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

// Dari berkas:
//   CountMatrixReader.ReadFeatureCounts(reader)     keluaran featureCounts (menyimpan panjang gen)
//   CountMatrixReader.ReadHtSeqFiles(paths)         satu berkas HTSeq-count per sampel
//   CountMatrixReader.ReadTable(reader)             tabel gen x sampel apa pun

// Eksperimen sintetis: 1.000 gen, 4 kontrol + 4 perlakuan, 60 perubahan nyata
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

// Normalisasi
double[] sizeFactors = Normalization.MedianOfRatiosSizeFactors(counts);  // DESeq2
double[] tmm = Normalization.TmmFactors(counts);                          // edgeR
Dataset logCpm = Normalization.LogCpm(counts, normFactors: tmm);
Console.WriteLine(string.Join(", ", sizeFactors.Select(f => f.ToString("F3"))));
```

## 2. Ekspresi diferensial

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

// Gaya DESeq2: GLM binomial negatif, penyusutan dispersi, uji Wald, koreksi BH
DifferentialExpressionResult deseq = DESeq.Run(counts, "group", "A", "B");
var (up, down) = deseq.CountSignificant(alpha: 0.05);
Console.WriteLine($"DESeq: {up} naik, {down} turun");
foreach (var row in deseq.Top(3))
    Console.WriteLine($"{row.Gene}  log2FC {row.Log2FoldChange:F2}  padj {row.PAdj:E2}");

// Gaya limma-voom
var design = DesignMatrix.FromAnnotations(counts, "group", "A");
DifferentialExpressionResult voom = Limma.VoomTest(counts, design);

// Nilai terstabilkan ragam dan PCA untuk plot QC
Dataset vst = Transforms.Vst(counts);
PcaResult pca = Pca.Fit(DenseMatrix.FromDataset(vst, transpose: true), components: 2);
Console.WriteLine($"PC1 menjelaskan {pca.ExplainedVarianceRatio[0]:P1}");

// Ekspor tabel hasil
var tsv = new StringWriter();
deseq.WriteTsv(tsv);
```

`Biotech.Net.Visualization` mengubah hasil menjadi volcano atau MA plot (`new VolcanoPlot(deseq)`, `new MaPlot(deseq)`):

![Volcano plot](../../images/volcano.png)

## 3. Analisis pengayaan

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
    Console.WriteLine($"{r.Set}: irisan {r.Overlap}/{r.SetSize}, p {r.PValue:E2}, padj {r.PAdj:E2}");

// GSEA preranked menerima pemetaan gen -> skor (misalnya log2FC atau statistik Wald)
var ranking = universe.Select((g, i) => (g, score: g.StartsWith("BG") ? (i % 7) - 3.0 : 5.0))
                      .ToDictionary(x => x.g, x => x.score);
var gsea = Gsea.Run(ranking, sets, permutations: 500, minSize: 3);
```

Gene set dapat dibaca dari berkas `.gmt` MSigDB dengan `Gmt.ReadFile(path)`.

## 4. RNA-seq sel tunggal

```csharp
using Biotech.Net.Transcriptomics;
using Biotech.Net.Transcriptomics.SingleCell;

// Dari disk: MatrixMarket.Read10x("filtered_feature_bc_matrix/")
// Di sini: 300 sel, 3 populasi, 100 gen
var rng = new Random(7);
var rows = new List<int>(); var cols = new List<int>(); var vals = new List<float>();
for (int c = 0; c < 300; c++)
    for (int g = 0; g < 100; g++)
    {
        int population = c / 100;
        double rate = (g / 10 == population) ? 8 : 1;   // gen 0-9, 10-19, 20-29 menandai tiap populasi
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
CommunityDetection.Cluster(sc, resolution: 0.5);          // menyimpan CellLabels["leiden"]
Umap.Fit(sc, new UmapOptions { Seed = 1 });                // menyimpan Embeddings["X_umap"]

var markers = MarkerGenes.Rank(sc, "leiden");
foreach (var (cluster, list) in markers.OrderBy(kv => kv.Key))
    Console.WriteLine($"klaster {cluster}: {string.Join(", ", list.Take(3).Select(m => m.Gene))}");
```

`EmbeddingPlot.FromSingleCell(sc, "X_umap", colorBy: "leiden")` memplot embedding-nya.

## Batasan

- `DESeq` mengikuti metode DESeq2 terpublikasi (size factor median-of-ratios, tren dispersi parametrik atau mean, penyusutan dispersi empirical-Bayes dengan penanganan pencilan, uji Wald, dan penyaringan independen), tetapi kesetaraan numerik persis dengan DESeq2 di R pada data publik masih menjadi item roadmap, dan penyusutan LFC apeglm belum ada.
- Pipeline sel tunggal menyimpan matriks di memori sebagai `SparseMatrix` CSR. Data skala atlas dan pembaca AnnData/h5ad masih di [roadmap](../roadmap.md).
- Graf tetangga dicari secara brute force eksak (O(n²)), dan t-SNE memakai metode eksak O(n²) (Barnes–Hut belum ada). Keduanya ditujukan untuk hingga beberapa ribu sel, bukan jutaan.

---
Sebelumnya: [Genomika](genomics.md) · Selanjutnya: [Proteomika](proteomics.md) · [Indeks tutorial](index.md)
