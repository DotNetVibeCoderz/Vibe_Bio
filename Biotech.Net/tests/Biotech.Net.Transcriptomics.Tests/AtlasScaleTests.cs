using Biotech.Net.Transcriptomics.SingleCell;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>
/// Out-of-core QC, normalisation, highly variable genes and PCA against scanpy 1.12.4 / scikit-learn 1.9.1 on a
/// 300-cell × 1 000-gene subset of 10x pbmc3k (Fixtures/AnnData/make_fixtures.py). Every algorithm runs on the
/// h5ad-backed matrix and on the memory-mapped store with small row blocks so the streaming paths are exercised.
/// </summary>
public class AtlasScaleTests
{
    private const int Block = 37; // deliberately not a divisor of 300

    private static readonly Lazy<(AnnData Data, SparseMatrix Counts)> Pbmc = new(() =>
    {
        var a = H5ad.Read(AnnFixtures.File("pbmc3k_subset.h5ad"));
        return (a, a.X!.Sparse!);
    });

    private static void Close(double[] expected, double[] actual, double rel, double abs = 1e-12, string what = "")
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
        {
            if (double.IsNaN(expected[i])) { Assert.True(double.IsNaN(actual[i]), $"{what}[{i}]: expected NaN, got {actual[i]}"); continue; }
            double tol = abs + rel * Math.Abs(expected[i]);
            Assert.True(Math.Abs(expected[i] - actual[i]) <= tol, $"{what}[{i}]: expected {expected[i]:R}, got {actual[i]:R}");
        }
    }

    private static double MaxRel(double[] e, double[] a)
        => e.Zip(a).Where(t => !double.IsNaN(t.First)).Max(t => Math.Abs(t.First - t.Second) / Math.Max(Math.Abs(t.First), 1e-300));

    /// <summary>Runs <paramref name="check"/> on the in-memory, h5ad-backed and memory-mapped versions of the counts.</summary>
    private static void ForEachSource(string tag, Action<ISparseMatrixSource> check)
    {
        check(Pbmc.Value.Counts);
        using var backed = H5ad.OpenBacked(AnnFixtures.File("pbmc3k_subset.h5ad"));
        check(backed.X!);
        using var mapped = MappedSparseMatrix.Create(AnnFixtures.TempPath($"{tag}.bnsparse"), backed.X!, blockRows: 64);
        check(mapped);
    }

    /// <summary>
    /// scanpy ranks Seurat v3 genes with numpy's unstable argsort, so genes with exactly equal normalised variance
    /// (identical count vectors) may be ordered either way. Flags must agree except inside a tie group that straddles
    /// the cut-off, where only the number selected must agree.
    /// </summary>
    private static void AssertSameSelection(bool[] expected, bool[] actual, double[] score)
    {
        var groups = Enumerable.Range(0, score.Length).GroupBy(j => score[j]);
        foreach (var g in groups)
        {
            var idx = g.ToArray();
            int e = idx.Count(j => expected[j]), a = idx.Count(j => actual[j]);
            Assert.Equal(e, a);
            if (e == 0 || e == idx.Length) foreach (int j in idx) Assert.Equal(expected[j], actual[j]);
        }
    }

    [Fact]
    public void QcMetrics_MatchScanpy()
    {
        var (_, obs) = AnnFixtures.Csv("pbmc_qc_obs.csv.gz");
        var (_, var) = AnnFixtures.Csv("pbmc_qc_var.csv.gz");
        var genes = Pbmc.Value.Data.Var.Index;
        ForEachSource("qc", src =>
        {
            var q = OutOfCore.CalculateQcMetrics(src, genes, blockRows: Block);
            var o = new AnnDataFrame(Pbmc.Value.Data.Obs.Index); var v = new AnnDataFrame(genes);
            q.WriteTo(o, v);
            foreach (var col in obs.Keys) Close(AnnFixtures.Numbers(obs[col]), o.GetDoubles(col), 1e-6, 1e-9, col);
            foreach (var col in var.Keys) Close(AnnFixtures.Numbers(var[col]), v.GetDoubles(col), 1e-6, 1e-9, col);
        });
    }

    [Fact]
    public void NormalizeTotalLog1p_MatchesScanpy()
    {
        var (_, sums) = AnnFixtures.Csv("pbmc_lognorm_sums.csv.gz");
        var (_, colSums) = AnnFixtures.Csv("pbmc_lognorm_colsums.csv.gz");
        var (_, entries) = AnnFixtures.Csv("pbmc_lognorm_entries.csv");
        var rowRef = AnnFixtures.Numbers(sums["row_sum"]);
        var medRef = AnnFixtures.Numbers(sums["median_norm_row_sum"]);
        ForEachSource("norm", src =>
        {
            using var lg = src is MappedSparseMatrix m
                ? OutOfCore.NormalizeTotalLog1pToLayer(m, "lognorm", blockRows: Block)
                : OutOfCore.NormalizeTotalLog1p(src, AnnFixtures.TempPath("lognorm.bnsparse"), blockRows: Block);
            var x = lg.ReadRows(0, lg.Rows);
            Close(rowRef, x.RowSums(), 1e-6, 0, "row_sum");
            Close(AnnFixtures.Numbers(colSums["col_sum"]), x.ColumnSums(), 1e-6, 1e-9, "col_sum");
            var r = AnnFixtures.Numbers(entries["col"]); // index column holds the row
            var (rowsIdx, _) = AnnFixtures.Csv("pbmc_lognorm_entries.csv");
            var vals = AnnFixtures.Numbers(entries["value"]);
            for (int k = 0; k < vals.Length; k++)
            {
                int row = int.Parse(rowsIdx[k]), col = (int)r[k];
                Assert.True(Math.Abs(vals[k] - x[row, col]) <= 1e-6 * vals[k], $"entry ({row},{col}): {vals[k]} vs {x[row, col]}");
            }
            using var med = OutOfCore.NormalizeTotalLog1p(src, AnnFixtures.TempPath("median.bnsparse"), targetSum: null, log1p: false, blockRows: Block);
            Close(medRef, med.ReadRows(0, med.Rows).RowSums(), 1e-6, 0, "median_norm_row_sum");
        });
    }

    [Fact]
    public void HighlyVariableGenes_SeuratV3_MatchesScanpy()
    {
        var (_, h) = AnnFixtures.Csv("pbmc_hvg_seurat_v3.csv.gz");
        var flags = h["highly_variable"].Select(s => s == "True").ToArray();
        var norm = AnnFixtures.Numbers(h["variances_norm"]);
        ForEachSource("hvg3", src =>
        {
            var r = OutOfCore.HighlyVariableGenes(src, 250, HvgFlavor.SeuratV3, blockRows: Block);
            Close(AnnFixtures.Numbers(h["means"]), r.Means, 1e-12, 0, "means");
            Close(AnnFixtures.Numbers(h["variances"]), r.Spread, 1e-10, 0, "variances");
            Close(norm, r.Normalized, 1e-9, 0, "variances_norm");
            AssertSameSelection(flags, r.HighlyVariable, norm);
            // Ranks agree except for permutations inside tie groups.
            var rank = AnnFixtures.Numbers(h["highly_variable_rank"]);
            foreach (var g in Enumerable.Range(0, norm.Length).Where(j => flags[j] && r.HighlyVariable[j]).GroupBy(j => norm[j]))
                Assert.Equal(g.Select(j => rank[j]).Order(), g.Select(j => r.Rank![j]).Order());
        });
    }

    [Fact]
    public void HighlyVariableGenes_Seurat_MatchesScanpy()
    {
        var (_, h) = AnnFixtures.Csv("pbmc_hvg_seurat.csv.gz");
        var flags = h["highly_variable"].Select(s => s == "True").ToArray();
        using var lg = OutOfCore.NormalizeTotalLog1p(Pbmc.Value.Counts, AnnFixtures.TempPath("hvg.bnsparse"));
        var r = OutOfCore.HighlyVariableGenes(lg, 250, HvgFlavor.Seurat, blockRows: Block);
        Close(AnnFixtures.Numbers(h["means"]), r.Means, 1e-6, 1e-12, "means");
        Close(AnnFixtures.Numbers(h["dispersions"]), r.Spread, 1e-6, 1e-9, "dispersions");
        Close(AnnFixtures.Numbers(h["dispersions_norm"]), r.Normalized, 1e-5, 1e-7, "dispersions_norm");
        Assert.Equal(flags, r.HighlyVariable);

        // The in-memory API shares the implementation (sums differ only by block order: ~1e-15).
        var sc = new SingleCellExperiment(lg.ReadRows(0, lg.Rows), Pbmc.Value.Data.Obs.Index, Pbmc.Value.Data.Var.Index);
        var sel = Preprocessing.HighlyVariableGenes(sc, 250, HvgFlavor.Seurat);
        Assert.Equal(r.Selected, sel);
        Close(r.Normalized, sc.GeneMetrics["dispersions_norm"], 1e-12, 1e-12, "in-memory dispersions_norm");
    }

    private static (MappedSparseMatrix Log, int[] Hvg) LogHvg(string tag)
    {
        var (_, h) = AnnFixtures.Csv("pbmc_hvg_seurat_v3.csv.gz");
        var hvg = h["highly_variable"].Select((s, j) => (s, j)).Where(t => t.s == "True").Select(t => t.j).ToArray();
        return (OutOfCore.NormalizeTotalLog1p(Pbmc.Value.Counts, AnnFixtures.TempPath($"{tag}.bnsparse")), hvg);
    }

    private static double[,] Table(string name, int rows, int cols)
    {
        var (_, c) = AnnFixtures.Csv(name);
        var t = new double[rows, cols];
        for (int k = 0; k < cols; k++)
        {
            var v = AnnFixtures.Numbers(c[k.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
            for (int i = 0; i < rows; i++) t[i, k] = v[i];
        }
        return t;
    }

    [Fact]
    public void Pca_Covariance_MatchesSklearnUpToSign()
    {
        var (log, hvg) = LogHvg("pca_cov");
        using var _ = log;
        var pca = OutOfCorePca.Fit(log, hvg, components: 20, solver: OutOfCorePcaSolver.Covariance, blockRows: Block);
        var (_, v) = AnnFixtures.Csv("pbmc_pca_variance.csv.gz");
        // The reference ran on scanpy's float32 log values, ours on our own (identical up to 1 ulp), hence 1e-8 rather than 1e-12.
        Close(AnnFixtures.Numbers(v["sklearn_variance"]), pca.ExplainedVariance, 1e-8, 0, "variance");
        Close(AnnFixtures.Numbers(v["sklearn_ratio"]), pca.ExplainedVarianceRatio, 1e-8, 0, "ratio");
        // scanpy computes in float32 (arpack): agreement at its precision.
        Close(AnnFixtures.Numbers(v["scanpy_variance"]), pca.ExplainedVariance, 1e-5, 0, "scanpy variance");
        var loadings = Table("pbmc_pca_loadings.csv.gz", hvg.Length, 10);
        var scores = Table("pbmc_pca_scores.csv.gz", 300, 10);
        for (int c = 0; c < 10; c++)
        {
            double dot = 0;
            for (int j = 0; j < hvg.Length; j++) dot += loadings[j, c] * pca.Loadings[j, c];
            double sign = Math.Sign(dot);
            Assert.True(Math.Abs(dot) > 1 - 1e-8, $"PC{c + 1} |cos| = {Math.Abs(dot)}");
            for (int i = 0; i < 300; i++) Assert.True(Math.Abs(scores[i, c] - sign * pca.Scores[i, c]) < 1e-6, $"score {i},{c}: {scores[i, c]:R} vs {sign * pca.Scores[i, c]:R}");
        }
    }

    [Fact]
    public void Pca_Randomized_AgreesWithinDocumentedTolerance()
    {
        var (log, hvg) = LogHvg("pca_rand");
        using var _ = log;
        var pca = OutOfCorePca.Fit(log, hvg, components: 20, solver: OutOfCorePcaSolver.Randomized, powerIterations: 6, oversamples: 20, blockRows: Block);
        var (_, v) = AnnFixtures.Csv("pbmc_pca_variance.csv.gz");
        var exact = AnnFixtures.Numbers(v["sklearn_variance"]);
        // Leading components: 0.1 %; the full 20: 2 %.
        Close(exact[..10], pca.ExplainedVariance[..10], 1e-3, 0, "variance (PC1-10)");
        Close(exact, pca.ExplainedVariance, 2e-2, 0, "variance (PC1-20)");
        var loadings = Table("pbmc_pca_loadings.csv.gz", hvg.Length, 10);
        for (int c = 0; c < 5; c++)
        {
            double dot = 0;
            for (int j = 0; j < hvg.Length; j++) dot += loadings[j, c] * pca.Loadings[j, c];
            Assert.True(Math.Abs(dot) > 0.999, $"PC{c + 1} |cos| = {Math.Abs(dot)}");
        }
    }

    [Fact]
    public void SubspaceEigen_MatchesTheDenseSolverAbovePca600Genes()
    {
        // A covariance-like PSD matrix with a decaying spectrum, larger than the dense-solver limit.
        int p = 700, r = 900, k = 30;
        var rng = new Random(3);
        var f = new double[r * p];
        for (int i = 0; i < r; i++) for (int j = 0; j < p; j++) f[i * p + j] = (rng.NextDouble() - 0.5) / Math.Sqrt(1 + j * 0.05);
        var a = new double[p * p];
        for (int i = 0; i < r; i++)
            for (int x = 0; x < p; x++)
            {
                double v = f[i * p + x];
                for (int y = 0; y < p; y++) a[x * p + y] += v * f[i * p + y];
            }
        Assert.True(p > OutOfCorePca.DenseEigenLimit);
        var (vs, vecs) = OutOfCorePca.TopEigen(a, p, k, 4);
        var e = Gravicode.Science.GraviNum.Decomposition.SymmetricEigen(new Gravicode.Science.GraviNum.NdArray((double[])a.Clone(), [p, p]));
        var ev = e.Values.ToArray(); var evec = e.Vectors.ToArray();
        var order = Enumerable.Range(0, p).OrderByDescending(i => ev[i]).ToArray();
        for (int c = 0; c < k; c++)
        {
            Assert.True(Math.Abs(vs[c] - ev[order[c]]) <= 1e-9 * ev[order[0]], $"eigenvalue {c}");
            double dot = 0;
            for (int j = 0; j < p; j++) dot += vecs[j * k + c] * evec[j * p + order[c]];
            Assert.True(Math.Abs(dot) > 1 - 1e-8, $"eigenvector {c}: {dot}");
        }
    }

    [Fact]
    public void Pca_MatchesInMemoryPca()
    {
        var (log, hvg) = LogHvg("pca_mem");
        using var _ = log;
        var ooc = OutOfCorePca.Fit(log, hvg, components: 10);
        var dense = log.ReadRows(0, log.Rows).ToDense(hvg);
        var mem = Pca.Fit(dense, 10, PcaSolver.Full);
        Close(mem.ExplainedVariance, ooc.ExplainedVariance, 1e-9, 0, "variance");
        for (int i = 0; i < 300; i++) for (int c = 0; c < 10; c++) Assert.True(Math.Abs(mem.Scores[i, c] - ooc.Scores[i, c]) < 1e-7);
    }

    [Fact]
    public void OutOfCorePipeline_FeedsTheExistingGraphClusteringAndUmap()
    {
        using var backed = H5ad.OpenBacked(AnnFixtures.File("pbmc3k_subset.h5ad"));
        using var store = MappedSparseMatrix.Create(AnnFixtures.TempPath("pipe.bnsparse"), backed.X!);
        var qc = OutOfCore.CalculateQcMetrics(store, backed.Var.Index);
        var hvg = OutOfCore.HighlyVariableGenes(store, 250, HvgFlavor.SeuratV3);
        using var log = OutOfCore.NormalizeTotalLog1pToLayer(store, cellTotals: qc.TotalCounts);
        var pca = OutOfCorePca.Fit(log, hvg.Selected, components: 15);
        var graph = NeighborGraph.Build(pca.Scores, nNeighbors: 10);
        var clusters = CommunityDetection.Cluster(graph.Connectivities, resolution: 0.8);
        Assert.Equal(300, clusters.Length);
        Assert.InRange(clusters.Distinct().Count(), 2, 30);
        var umap = Umap.Fit(graph, new UmapOptions { Epochs = 50, Seed = 1 });
        Assert.Equal(300, umap.Rows);

        // And into a SingleCellExperiment for the in-memory tools (markers etc.) on whatever fits.
        var sc = new SingleCellExperiment(log.ReadRows(0, log.Rows), backed.Obs.Index, backed.Var.Index);
        qc.WriteTo(sc); hvg.WriteTo(sc);
        sc.Embeddings["X_pca"] = pca.Scores;
        sc.CellLabels["leiden"] = clusters.Select(c => c.ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray();
        var markers = MarkerGenes.Rank(sc, "leiden");
        Assert.NotEmpty(markers);
    }
}
