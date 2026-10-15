using Biotech.Net.Transcriptomics.SingleCell;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>
/// The single-cell pipeline end to end: QC, filtering, normalisation, highly variable genes, PCA, the kNN
/// graph, community detection, embeddings and marker genes, on a simulation with planted cluster structure.
/// </summary>
public class SingleCellTests
{
    private const int ClusterCount = 3;
    private const int CellsPerCluster = 60;
    private const int GeneCount = 300;
    private const int MarkersPerCluster = 10;
    private const double MarkerFoldChange = 12.0;

    /// <summary>
    /// Three cell populations of 60 cells over 300 genes. Genes 0–9 are specific to cluster 0, 10–19 to
    /// cluster 1 and 20–29 to cluster 2: each is expressed <see cref="MarkerFoldChange"/>-fold higher in its
    /// own cluster than anywhere else. The remaining 270 genes are plain Poisson noise with the same mean in
    /// every cluster, so any structure a method finds must come from the planted markers.
    /// </summary>
    /// <remarks>
    /// Baseline means are drawn log-uniformly over 0.2–20 counts rather than fixed, so the genes spread across
    /// the whole expression range. That matters for feature selection: the binned dispersion of Seurat (and the
    /// mean–variance loess of Seurat v3) compares each gene with others of similar mean, and a simulation where
    /// every interesting gene shares one expression bin would tell us nothing about either method.
    /// </remarks>
    private static (SingleCellExperiment Sc, int[] Truth) Simulate(int seed = 4242)
    {
        var random = new Random(seed);
        int cells = ClusterCount * CellsPerCluster;
        var truth = new int[cells];

        // Baseline expression level of each gene, log-uniform over 0.2-20 counts.
        var baseline = new double[GeneCount];
        for (int gene = 0; gene < GeneCount; gene++)
            baseline[gene] = Math.Exp(Math.Log(0.2) + random.NextDouble() * (Math.Log(20) - Math.Log(0.2)));

        var rowIdx = new List<int>();
        var colIdx = new List<int>();
        var values = new List<float>();

        for (int cell = 0; cell < cells; cell++)
        {
            int cluster = cell / CellsPerCluster;
            truth[cell] = cluster;
            for (int gene = 0; gene < GeneCount; gene++)
            {
                bool isOwnMarker = gene >= cluster * MarkersPerCluster && gene < (cluster + 1) * MarkersPerCluster;
                double mean = isOwnMarker ? baseline[gene] * MarkerFoldChange : baseline[gene];

                int count = Poisson(random, mean);
                if (count == 0) continue;
                rowIdx.Add(cell);
                colIdx.Add(gene);
                values.Add(count);
            }
        }

        var matrix = SparseMatrix.FromTriplets(cells, GeneCount, rowIdx.ToArray(), colIdx.ToArray(), values.ToArray());
        var cellNames = Enumerable.Range(0, cells).Select(i => $"cell{i}").ToArray();
        var geneIds = Enumerable.Range(0, GeneCount).Select(i => $"G{i}").ToArray();
        // A handful of mitochondrial genes for the QC test.
        var geneNames = geneIds.Select((g, i) => i < 5 ? $"MT-{g}" : g).ToArray();

        return (new SingleCellExperiment(matrix, cellNames, geneIds, geneNames), truth);
    }

    private static int Poisson(Random random, double mean)
    {
        // Knuth's method; the means here are small enough for it.
        double limit = Math.Exp(-mean), product = random.NextDouble();
        int k = 0;
        while (product > limit) { k++; product *= random.NextDouble(); }
        return k;
    }

    /// <summary>Adjusted Rand index between two labellings, the standard clustering-agreement measure.</summary>
    private static double AdjustedRandIndex(IReadOnlyList<int> a, IReadOnlyList<int> b)
    {
        int n = a.Count;
        var table = new Dictionary<(int, int), int>();
        var rowTotals = new Dictionary<int, int>();
        var columnTotals = new Dictionary<int, int>();
        for (int i = 0; i < n; i++)
        {
            table[(a[i], b[i])] = table.GetValueOrDefault((a[i], b[i])) + 1;
            rowTotals[a[i]] = rowTotals.GetValueOrDefault(a[i]) + 1;
            columnTotals[b[i]] = columnTotals.GetValueOrDefault(b[i]) + 1;
        }

        static double Choose2(double x) => x * (x - 1) / 2.0;
        double index = table.Values.Sum(v => Choose2(v));
        double rowSum = rowTotals.Values.Sum(v => Choose2(v));
        double columnSum = columnTotals.Values.Sum(v => Choose2(v));
        double expected = rowSum * columnSum / Choose2(n);
        double max = (rowSum + columnSum) / 2.0;
        return (index - expected) / (max - expected);
    }

    // ---- QC and filtering ----

    [Fact]
    public void QcMetrics_CountGenesTotalsAndMitochondrialFraction()
    {
        var (sc, _) = Simulate();
        Preprocessing.CalculateQcMetrics(sc);

        var nGenes = sc.CellMetrics[Preprocessing.NGenes];
        var totals = sc.CellMetrics[Preprocessing.TotalCounts];
        var pctMito = sc.CellMetrics[Preprocessing.PctMito];
        var mitoTotals = sc.CellMetrics[Preprocessing.TotalMito];

        Assert.Equal(sc.CellCount, nGenes.Length);

        for (int cell = 0; cell < sc.CellCount; cell++)
        {
            // Recount straight from the sparse rows.
            var indices = sc.X.RowIndices(cell);
            var row = sc.X.RowValues(cell);
            Assert.Equal(indices.Length, (int)nGenes[cell]);
            Assert.Equal(SumOf(row), totals[cell], 6);

            double mito = 0;
            for (int k = 0; k < indices.Length; k++) if (indices[k] < 5) mito += row[k];
            Assert.Equal(mito, mitoTotals[cell], 6);
            Assert.Equal(totals[cell] == 0 ? 0 : 100 * mito / totals[cell], pctMito[cell], 6);
            Assert.InRange(pctMito[cell], 0, 100);
        }
    }

    [Fact]
    public void FilterCells_And_FilterGenes_RemoveWhatTheyAreAskedTo()
    {
        var (sc, _) = Simulate();
        Preprocessing.CalculateQcMetrics(sc);

        int before = sc.CellCount;
        var totals = sc.CellMetrics[Preprocessing.TotalCounts].ToArray();
        double threshold = totals.OrderBy(x => x).ElementAt(totals.Length / 2);

        int removed = Preprocessing.FilterCells(sc, minCounts: threshold);
        Assert.Equal(before - sc.CellCount, removed);
        Assert.True(sc.CellCount < before);
        Preprocessing.CalculateQcMetrics(sc);
        Assert.All(sc.CellMetrics[Preprocessing.TotalCounts], t => Assert.True(t >= threshold));

        var (fresh, _) = Simulate();
        int genesBefore = fresh.GeneCount;
        int genesRemoved = Preprocessing.FilterGenes(fresh, minCells: 3);
        Assert.Equal(genesBefore - fresh.GeneCount, genesRemoved);
        // Every surviving gene is seen in at least three cells.
        var counts = fresh.X.ColumnNonZeroCounts();
        Assert.All(counts, c => Assert.True(c >= 3));
    }

    [Fact]
    public void NormalizeTotal_MakesEveryCellSumToTheTargetAndLog1pIsItsOwnInverse()
    {
        var (sc, _) = Simulate();
        Preprocessing.NormalizeTotal(sc, targetSum: 1e4);

        foreach (double total in sc.X.RowSums()) Assert.Equal(1e4, total, 1);

        var before = sc.X.Values.ToArray();
        Preprocessing.Log1p(sc);
        // The values are stored as float, so compare with a tolerance rather than decimal places.
        for (int i = 0; i < before.Length; i++)
            Assert.Equal(Math.Log(1 + before[i]), sc.X.Values[i], tolerance: 1e-5);
    }

    // ---- feature selection ----

    [Fact]
    public void HighlyVariableGenes_PickThePlantedMarkers()
    {
        var (sc, _) = Simulate();
        Preprocessing.NormalizeTotal(sc);
        Preprocessing.Log1p(sc);

        var hvg = Preprocessing.HighlyVariableGenes(sc, nTopGenes: 50);

        Assert.Equal(50, hvg.Length);
        Assert.Equal(hvg.Length, hvg.Distinct().Count());
        Assert.True(sc.GeneFlags.ContainsKey("highly_variable"));
        Assert.Equal(50, sc.GeneFlags["highly_variable"].Count(f => f));

        // The 30 planted marker genes are exactly the genes that vary between clusters, so the top 50 should
        // contain nearly all of them.
        int planted = hvg.Count(g => g < ClusterCount * MarkersPerCluster);
        // Chance alone would put 50/300 * 30 = 5 of them in the list.
        Assert.True(planted >= 18, $"only {planted} of the 30 planted markers were selected");

        // Seurat v3 selects on the mean-variance trend instead and should do at least as well.
        var (fresh, _) = Simulate();
        Preprocessing.NormalizeTotal(fresh);
        Preprocessing.Log1p(fresh);
        var v3 = Preprocessing.HighlyVariableGenes(fresh, nTopGenes: 50, flavor: HvgFlavor.SeuratV3);
        Assert.True(v3.Count(g => g < ClusterCount * MarkersPerCluster) >= 18);
        Assert.True(fresh.GeneMetrics.ContainsKey("variances_norm"));
    }

    [Fact]
    public void Scale_CentresAndClipsTheSelectedGenes()
    {
        var (sc, _) = Simulate();
        Preprocessing.NormalizeTotal(sc);
        Preprocessing.Log1p(sc);
        var hvg = Preprocessing.HighlyVariableGenes(sc, nTopGenes: 50);

        var scaled = Preprocessing.Scale(sc, hvg, maxValue: 10);

        Assert.Equal(sc.CellCount, scaled.Rows);
        Assert.Equal(hvg.Length, scaled.Columns);
        for (int column = 0; column < scaled.Columns; column++)
        {
            var values = scaled.Column(column);
            Assert.Equal(0.0, values.Average(), 6);          // zero-centred
            Assert.All(values, v => Assert.InRange(v, -10.0, 10.0));
        }
    }

    // ---- dimensionality reduction and graph ----

    private static (SingleCellExperiment Sc, int[] Truth, DenseMatrix Pcs) Processed(int components = 15)
    {
        var (sc, truth) = Simulate();
        Preprocessing.NormalizeTotal(sc);
        Preprocessing.Log1p(sc);
        var hvg = Preprocessing.HighlyVariableGenes(sc, nTopGenes: 60);
        var scaled = Preprocessing.Scale(sc, hvg);
        var pca = Pca.Fit(scaled, components);
        sc.Embeddings["X_pca"] = pca.Scores;
        return (sc, truth, pca.Scores);
    }

    [Fact]
    public void Pca_ProducesOrderedComponentsThatSeparateTheClusters()
    {
        var (_, truth, pcs) = Processed();

        Assert.Equal(ClusterCount * CellsPerCluster, pcs.Rows);
        Assert.Equal(15, pcs.Columns);

        // Cluster centroids in the first two PCs must be further apart than the within-cluster spread.
        double between = 0, within = 0;
        for (int a = 0; a < ClusterCount; a++)
        {
            var cellsA = Enumerable.Range(0, truth.Length).Where(i => truth[i] == a).ToArray();
            double cx = cellsA.Average(i => pcs[i, 0]), cy = cellsA.Average(i => pcs[i, 1]);
            within += cellsA.Average(i => Math.Sqrt((pcs[i, 0] - cx) * (pcs[i, 0] - cx) + (pcs[i, 1] - cy) * (pcs[i, 1] - cy)));
            for (int b = a + 1; b < ClusterCount; b++)
            {
                var cellsB = Enumerable.Range(0, truth.Length).Where(i => truth[i] == b).ToArray();
                double dx = cx - cellsB.Average(i => pcs[i, 0]), dy = cy - cellsB.Average(i => pcs[i, 1]);
                between += Math.Sqrt(dx * dx + dy * dy);
            }
        }
        Assert.True(between / 3 > within / ClusterCount,
            "the planted clusters are not separated in the first two principal components");
    }

    [Fact]
    public void NeighborGraph_IsSymmetricAndMostlyWithinCluster()
    {
        var (sc, truth, _) = Processed();
        var graph = NeighborGraph.Build(sc, nNeighbors: 15);

        Assert.Equal(sc.CellCount, graph.CellCount);
        // K counts neighbours excluding the cell itself, so a request for 15 stores 14.
        Assert.Equal(14, graph.K);
        Assert.Equal(sc.CellCount * graph.K, graph.Indices.Length);
        Assert.All(graph.Distances, d => Assert.True(d >= 0));

        // Distances are sorted nearest-first for each cell, and no cell is its own neighbour.
        for (int cell = 0; cell < sc.CellCount; cell++)
        {
            for (int k = 1; k < graph.K; k++)
                Assert.True(graph.Distances[cell * graph.K + k] >= graph.Distances[cell * graph.K + k - 1]);
            for (int k = 0; k < graph.K; k++)
                Assert.NotEqual(cell, graph.Indices[cell * graph.K + k]);
        }

        // The fuzzy-union connectivity graph is symmetric.
        var connectivities = graph.Connectivities;
        Assert.Equal(sc.CellCount, connectivities.Rows);
        for (int row = 0; row < connectivities.Rows; row += 17)
        {
            var indices = connectivities.RowIndices(row);
            var values = connectivities.RowValues(row);
            for (int k = 0; k < indices.Length; k++)
                Assert.Equal(values[k], connectivities[indices[k], row], 6);
        }

        // Nearly every neighbour should come from the same planted cluster.
        int sameCluster = 0;
        for (int cell = 0; cell < sc.CellCount; cell++)
            for (int k = 0; k < graph.K; k++)
                if (truth[graph.Indices[cell * graph.K + k]] == truth[cell]) sameCluster++;
        double purity = (double)sameCluster / (sc.CellCount * graph.K);
        Assert.True(purity > 0.9, $"neighbour purity was only {purity:P0}");
    }

    [Theory]
    [InlineData(CommunityAlgorithm.Leiden)]
    [InlineData(CommunityAlgorithm.Louvain)]
    public void CommunityDetection_RecoversThePlantedClusters(CommunityAlgorithm algorithm)
    {
        var (sc, truth, _) = Processed();
        sc.Neighbors = NeighborGraph.Build(sc, nNeighbors: 15);

        var labels = CommunityDetection.Cluster(sc, resolution: 1.0, algorithm: algorithm);

        Assert.Equal(sc.CellCount, labels.Length);
        double ari = AdjustedRandIndex(truth, labels);
        Assert.True(ari > 0.9, $"{algorithm} recovered the clusters with ARI {ari:F3}");
        Assert.Equal(ClusterCount, labels.Distinct().Count());

        // Labels are a dense 0..k-1 range.
        Assert.Equal(Enumerable.Range(0, labels.Distinct().Count()), labels.Distinct().Order());
    }

    [Fact]
    public void CommunityDetection_IsDeterministicForAGivenSeed()
    {
        var (sc, _, _) = Processed();
        sc.Neighbors = NeighborGraph.Build(sc, nNeighbors: 15);

        var first = CommunityDetection.Cluster(sc.Neighbors.Connectivities);
        var second = CommunityDetection.Cluster(sc.Neighbors.Connectivities);
        Assert.Equal(first, second);
    }

    [Fact]
    public void Modularity_IsHigherForTheTruthThanForAShuffling()
    {
        var (sc, truth, _) = Processed();
        var graph = NeighborGraph.Build(sc, nNeighbors: 15);

        var shuffled = truth.ToArray();
        var random = new Random(1);
        for (int i = shuffled.Length - 1; i > 0; i--)
        {
            int j = random.Next(i + 1);
            (shuffled[i], shuffled[j]) = (shuffled[j], shuffled[i]);
        }

        double truthModularity = CommunityDetection.Modularity(graph.Connectivities, truth);
        double shuffledModularity = CommunityDetection.Modularity(graph.Connectivities, shuffled);

        Assert.True(truthModularity > shuffledModularity);
        Assert.InRange(truthModularity, -1.0, 1.0);
        Assert.True(truthModularity > 0.3, $"modularity of the true partition was only {truthModularity:F3}");
    }

    [Fact]
    public void HigherResolutionGivesAtLeastAsManyClusters()
    {
        var (sc, _, _) = Processed();
        var graph = NeighborGraph.Build(sc, nNeighbors: 15);

        int coarse = CommunityDetection.Cluster(graph.Connectivities, resolution: 0.2).Distinct().Count();
        int fine = CommunityDetection.Cluster(graph.Connectivities, resolution: 3.0).Distinct().Count();

        Assert.True(fine >= coarse, $"resolution 3.0 gave {fine} clusters, resolution 0.2 gave {coarse}");
    }

    // ---- embeddings ----

    [Fact]
    public void Umap_ReturnsTheRightShapeAndKeepsTheClustersApart()
    {
        var (sc, truth, _) = Processed();
        sc.Neighbors = NeighborGraph.Build(sc, nNeighbors: 15);

        var embedding = Umap.Fit(sc, new UmapOptions { Components = 2, Epochs = 100, Seed = 7 });

        Assert.Equal(sc.CellCount, embedding.Rows);
        Assert.Equal(2, embedding.Columns);
        Assert.All(embedding.Data, v => Assert.True(double.IsFinite(v)));

        Assert.True(CentroidSeparation(embedding, truth) > 1.0,
            "UMAP did not separate the planted clusters");
    }

    [Fact]
    public void Umap_IsReproducibleForAGivenSeed()
    {
        var (sc, _, _) = Processed();
        sc.Neighbors = NeighborGraph.Build(sc, nNeighbors: 15);
        var options = new UmapOptions { Components = 2, Epochs = 50, Seed = 11 };

        var a = Umap.Fit(sc, options);
        var b = Umap.Fit(sc, options);

        Assert.Equal(a.Data, b.Data);
    }

    [Fact]
    public void UmapAbParameters_ReproduceTheTargetCurve()
    {
        var (a, b) = Umap.FindAbParams(spread: 1.0, minDist: 0.1);

        Assert.True(a > 0 && b > 0);
        // The fitted curve 1/(1 + a d^(2b)) must match the piecewise target at the knee: at d = minDist the
        // target is 1, and it decays like exp(-(d - minDist)/spread) beyond it.
        foreach (double d in new[] { 0.3, 0.5, 1.0, 2.0 })
        {
            double fitted = 1.0 / (1.0 + a * Math.Pow(d, 2 * b));
            double target = Math.Exp(-(d - 0.1) / 1.0);
            Assert.True(Math.Abs(fitted - target) < 0.1, $"at d={d} the fit gave {fitted:F3}, target {target:F3}");
        }
    }

    [Fact]
    public void Tsne_ReturnsTheRightShapeAndAFiniteDivergence()
    {
        var (sc, truth, pcs) = Processed();

        var (embedding, kl) = Tsne.Fit(pcs, new TsneOptions { Components = 2, Seed = 3 });

        Assert.Equal(pcs.Rows, embedding.Rows);
        Assert.Equal(2, embedding.Columns);
        Assert.All(embedding.Data, v => Assert.True(double.IsFinite(v)));
        Assert.True(kl > 0 && double.IsFinite(kl), $"KL divergence was {kl}");

        Assert.True(CentroidSeparation(embedding, truth) > 1.0,
            "t-SNE did not separate the planted clusters");
    }

    [Fact]
    public void Tsne_IsReproducibleForAGivenSeed()
    {
        var (_, _, pcs) = Processed();
        var options = new TsneOptions { Perplexity = 20, Iterations = 100, Seed = 5 };

        var a = Tsne.Fit(pcs, options);
        var b = Tsne.Fit(pcs, options);

        Assert.Equal(a.Embedding.Data, b.Embedding.Data);
        Assert.Equal(a.KlDivergence, b.KlDivergence, 12);
    }

    /// <summary>Mean distance between cluster centroids divided by the mean within-cluster spread.</summary>
    private static double CentroidSeparation(DenseMatrix embedding, int[] truth)
    {
        var centroids = new double[ClusterCount][];
        var spread = new double[ClusterCount];
        for (int c = 0; c < ClusterCount; c++)
        {
            var cells = Enumerable.Range(0, truth.Length).Where(i => truth[i] == c).ToArray();
            centroids[c] = [cells.Average(i => embedding[i, 0]), cells.Average(i => embedding[i, 1])];
            spread[c] = cells.Average(i => Distance(embedding[i, 0], embedding[i, 1], centroids[c][0], centroids[c][1]));
        }

        double between = 0;
        int pairs = 0;
        for (int a = 0; a < ClusterCount; a++)
            for (int b = a + 1; b < ClusterCount; b++)
            {
                between += Distance(centroids[a][0], centroids[a][1], centroids[b][0], centroids[b][1]);
                pairs++;
            }

        double meanSpread = spread.Average();
        return meanSpread <= 0 ? double.PositiveInfinity : between / pairs / meanSpread;
    }

    private static double Distance(double x1, double y1, double x2, double y2)
        => Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

    // ---- marker genes ----

    [Theory]
    [InlineData(MarkerTest.Wilcoxon)]
    [InlineData(MarkerTest.TTest)]
    public void MarkerGenes_RecoverThePlantedMarkersOfEachCluster(MarkerTest test)
    {
        var (sc, truth) = Simulate();
        Preprocessing.NormalizeTotal(sc);
        Preprocessing.Log1p(sc);

        var labels = truth.Select(t => $"c{t}").ToArray();
        var markers = MarkerGenes.Rank(sc, labels, test);

        Assert.Equal(ClusterCount, markers.Count);
        for (int cluster = 0; cluster < ClusterCount; cluster++)
        {
            var top = markers[$"c{cluster}"].Take(MarkersPerCluster).ToList();
            Assert.Equal(MarkersPerCluster, top.Count);

            // Every top gene must be one of the genes planted for this cluster.
            foreach (var marker in top)
            {
                Assert.InRange(marker.GeneIndex, cluster * MarkersPerCluster, (cluster + 1) * MarkersPerCluster - 1);
                Assert.Equal($"G{marker.GeneIndex}", marker.Gene == $"MT-G{marker.GeneIndex}" ? $"G{marker.GeneIndex}" : marker.Gene);
                Assert.True(marker.Log2FoldChange > 0, $"{marker.Gene} should be up in c{cluster}");
                Assert.True(marker.PValue < 1e-6, $"{marker.Gene} p-value was {marker.PValue}");
                Assert.InRange(marker.PValue, 0, 1);
            }

            // Scores are sorted, best first.
            var scores = markers[$"c{cluster}"].Select(m => m.Score).ToList();
            Assert.Equal(scores.OrderByDescending(s => s), scores);
        }
    }

    [Fact]
    public void MarkerGenes_UseTheStoredClusterLabelsByDefault()
    {
        var (sc, truth) = Simulate();
        Preprocessing.NormalizeTotal(sc);
        Preprocessing.Log1p(sc);
        sc.CellLabels["leiden"] = truth.Select(t => $"c{t}").ToArray();

        var markers = MarkerGenes.Rank(sc);

        Assert.Equal(ClusterCount, markers.Count);
        Assert.All(markers.Values, list => Assert.Equal(sc.GeneCount, list.Count));
    }

    [Fact]
    public void MarkerGenes_RejectBadInput()
    {
        var (sc, truth) = Simulate();
        Assert.Throws<ArgumentException>(() => MarkerGenes.Rank(sc, ["only", "a", "few"]));
        Assert.Throws<ArgumentException>(() => MarkerGenes.Rank(sc, truth.Select(_ => "same").ToArray()));
    }

    // ---- the experiment container ----

    [Fact]
    public void SubsetCellsAndGenes_KeepEverythingConsistent()
    {
        var (sc, _) = Simulate();
        var copy = sc.Copy();

        copy.SubsetCells([0, 1, 2, 3, 4]);
        Assert.Equal(5, copy.CellCount);
        Assert.Equal(GeneCount, copy.GeneCount);
        Assert.Equal(["cell0", "cell1", "cell2", "cell3", "cell4"], copy.CellNames);

        copy.SubsetGenes([10, 11, 12]);
        Assert.Equal(3, copy.GeneCount);
        Assert.Equal(["G10", "G11", "G12"], copy.GeneIds);
        for (int cell = 0; cell < copy.CellCount; cell++)
            for (int gene = 0; gene < 3; gene++)
                Assert.Equal(sc.X[cell, 10 + gene], copy.X[cell, gene]);

        // The original is untouched.
        Assert.Equal(ClusterCount * CellsPerCluster, sc.CellCount);
        Assert.Equal(GeneCount, sc.GeneCount);
    }

    [Fact]
    public void GeneIndex_FindsByIdOrNameAndReportsMisses()
    {
        var (sc, _) = Simulate();
        Assert.Equal(50, sc.GeneIndex("G50"));
        Assert.Equal(0, sc.GeneIndex("MT-G0"));
        Assert.Equal(-1, sc.TryGeneIndex("nothing"));
        Assert.Throws<KeyNotFoundException>(() => sc.GeneIndex("nothing"));
    }

    [Fact]
    public void MatrixMarket_RoundTrips()
    {
        var (sc, _) = Simulate();

        var writer = new StringWriter();
        MatrixMarket.Write(writer, sc.X);
        var read = MatrixMarket.Read(new StringReader(writer.ToString()));

        Assert.Equal(sc.X.Rows, read.Rows);
        Assert.Equal(sc.X.Columns, read.Columns);
        Assert.Equal(sc.X.NonZeroCount, read.NonZeroCount);
        Assert.Equal(sc.X.Values, read.Values);
        Assert.Equal(sc.X.ColumnIndices, read.ColumnIndices);
    }

    // ---- sparse matrix ----

    [Fact]
    public void SparseMatrix_SumsAndCountsMatchTheDenseForm()
    {
        var (sc, _) = Simulate();
        var dense = sc.X.ToDense();

        var rowSums = sc.X.RowSums();
        var columnSums = sc.X.ColumnSums();
        for (int row = 0; row < dense.Rows; row += 13)
            Assert.Equal(dense.Row(row).ToArray().Sum(), rowSums[row], 6);
        for (int column = 0; column < dense.Columns; column += 29)
            Assert.Equal(dense.Column(column).Sum(), columnSums[column], 6);

        var (mean, variance) = sc.X.ColumnMeanVariance();
        for (int column = 0; column < dense.Columns; column += 29)
        {
            var values = dense.Column(column);
            double m = values.Average();
            Assert.Equal(m, mean[column], 6);
            // ColumnMeanVariance documents the sample variance (ddof = 1).
            double sampleVariance = values.Sum(v => (v - m) * (v - m)) / (values.Length - 1);
            Assert.Equal(sampleVariance, variance[column], 4);
        }

        Assert.Equal(sc.X.NonZeroCount, sc.X.RowNonZeroCounts().Sum());
        Assert.Equal(sc.X.NonZeroCount, sc.X.ColumnNonZeroCounts().Sum());
    }

    [Fact]
    public void SparseMatrix_TransposeIsAnInvolution()
    {
        var (sc, _) = Simulate();
        var twice = sc.X.Transpose().Transpose();

        Assert.Equal(sc.X.Rows, twice.Rows);
        Assert.Equal(sc.X.Columns, twice.Columns);
        Assert.Equal(sc.X.Values, twice.Values);
        Assert.Equal(sc.X.ColumnIndices, twice.ColumnIndices);
        Assert.Equal(sc.X.RowPointers, twice.RowPointers);
    }

    [Fact]
    public void SparseMatrix_FromDenseAndBackIsLossless()
    {
        var dense = new DenseMatrix(3, 4, [0, 1, 0, 2, 0, 0, 0, 0, 3, 0, 4, 0]);
        var sparse = SparseMatrix.FromDense(dense);

        Assert.Equal(4, sparse.NonZeroCount);
        Assert.Equal(dense.Data, sparse.ToDense().Data);
        Assert.Equal(2.0f, sparse[0, 3]);
        Assert.Equal(0.0f, sparse[1, 1]);
    }

    private static double SumOf(ReadOnlySpan<float> values)
    {
        double sum = 0;
        foreach (float v in values) sum += v;
        return sum;
    }
}
