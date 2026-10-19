using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using Biotech.Net.Transcriptomics.SingleCell;

namespace Biotech.Net.Transcriptomics.Tests;

/// <summary>Shared access to the AnnData fixtures written by Fixtures/AnnData/make_fixtures.py (anndata 0.13.4).</summary>
internal static class AnnFixtures
{
    public static string Dir => Path.Combine(AppContext.BaseDirectory, "Fixtures", "AnnData");

    public static string File(string name) => Path.Combine(Dir, name);

    private static JsonElement? _expected;

    public static JsonElement Expected => _expected ??= JsonDocument.Parse(System.IO.File.ReadAllText(File("expected_values.json"))).RootElement;

    public static double[] Doubles(string key) => Expected.GetProperty(key).EnumerateArray().Select(e => e.GetDouble()).ToArray();

    /// <summary>Reads a gzipped pandas CSV: header row, then index + columns. Empty cells become NaN.</summary>
    public static (string[] Index, Dictionary<string, string[]> Columns) Csv(string name)
    {
        using var fs = System.IO.File.OpenRead(File(name));
        using Stream s = name.EndsWith(".gz", StringComparison.Ordinal) ? new GZipStream(fs, CompressionMode.Decompress) : fs;
        using var r = new StreamReader(s);
        var header = SplitCsv(r.ReadLine()!);
        var rows = new List<string[]>();
        string? line;
        while ((line = r.ReadLine()) is not null) if (line.Length > 0) rows.Add(SplitCsv(line));
        var cols = new Dictionary<string, string[]>(StringComparer.Ordinal);
        for (int c = 1; c < header.Length; c++) cols[header[c]] = rows.Select(x => x[c]).ToArray();
        return (rows.Select(x => x[0]).ToArray(), cols);
    }

    public static double[] Numbers(string[] column) => column.Select(v => v.Length == 0 ? double.NaN
        : v == "True" ? 1 : v == "False" ? 0 : double.Parse(v, CultureInfo.InvariantCulture)).ToArray();

    private static string[] SplitCsv(string line)
    {
        var parts = new List<string>();
        var cur = new System.Text.StringBuilder();
        bool quoted = false;
        foreach (char ch in line)
        {
            if (ch == '"') quoted = !quoted;
            else if (ch == ',' && !quoted) { parts.Add(cur.ToString()); cur.Clear(); }
            else cur.Append(ch);
        }
        parts.Add(cur.ToString());
        return parts.ToArray();
    }

    public static string TempPath(string name)
    {
        var dir = Path.Combine(Path.GetTempPath(), "biotechnet-anndata-tests", Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        return Path.Combine(dir, name);
    }
}

/// <summary>
/// h5ad reading against values dumped by anndata/numpy (not by an HDF5 reader), writing round trips, backed mode
/// and the memory-mapped store.
/// </summary>
public class AnnDataTests
{
    private static void AssertMatrix(double[] expected, AnnMatrix m)
    {
        Assert.Equal(expected.Length, m.Rows * m.Columns);
        for (int i = 0; i < m.Rows; i++)
            for (int j = 0; j < m.Columns; j++)
                Assert.Equal(expected[i * m.Columns + j], m[i, j]);
    }

    private static void AssertSame(SparseMatrix a, SparseMatrix b)
    {
        Assert.Equal(a.Rows, b.Rows);
        Assert.Equal(a.Columns, b.Columns);
        Assert.Equal(a.RowPointers, b.RowPointers);
        Assert.Equal(a.ColumnIndices.AsSpan(0, a.NonZeroCount).ToArray(), b.ColumnIndices.AsSpan(0, b.NonZeroCount).ToArray());
        Assert.Equal(a.Values.AsSpan(0, a.NonZeroCount).ToArray(), b.Values.AsSpan(0, b.NonZeroCount).ToArray());
    }

    // ---- reading ----

    [Fact]
    public void KitchenSink_MatricesMatchNumpyExactly()
    {
        var a = H5ad.Read(AnnFixtures.File("kitchen_sink.h5ad"));
        Assert.Equal(24, a.ObsCount);
        Assert.Equal(15, a.VarCount);
        Assert.Equal(AnnMatrixFormat.Dense, a.X!.Format);
        Assert.Equal("float32", a.X.DataType);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X"), a.X);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X_pca"), (AnnMatrix)a.Obsm["X_pca"]);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X_umap"), (AnnMatrix)a.Obsm["X_umap"]);
        Assert.Equal("float32", ((AnnMatrix)a.Obsm["X_umap"]).DataType);
        AssertMatrix(AnnFixtures.Doubles("kitchen_PCs"), (AnnMatrix)a.Varm["PCs"]);
        AssertMatrix(AnnFixtures.Doubles("kitchen_distances"), a.Obsp["distances"]);
        Assert.Equal(115, a.Obsp["distances"].Sparse!.NonZeroCount);
        Assert.Equal(68, a.Varp["corr"].Sparse!.NonZeroCount);
        // layers: counts (CSR int32) equals X, dense_layer = 2X; raw.X equals X.
        Assert.Equal("int32", a.Layers["counts"].DataType);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X"), a.Layers["counts"]);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X").Select(v => 2 * v).ToArray(), a.Layers["dense_layer"]);
        Assert.NotNull(a.Raw);
        AssertMatrix(AnnFixtures.Doubles("kitchen_X"), a.Raw!.X);
        Assert.Equal(15, a.Raw.Var.RowCount);
        Assert.True(a.Raw.Varm.ContainsKey("PCs"));
    }

    [Fact]
    public void KitchenSink_DataFramesKeepEveryColumnType()
    {
        var a = H5ad.Read(AnnFixtures.File("kitchen_sink.h5ad"));
        var obs = a.Obs;
        Assert.Equal("barcode", obs.IndexName);
        Assert.Equal("cell0", obs.Index[0]);
        Assert.Equal(["cluster", "stage", "batch_num", "n_counts", "n_genes", "score", "big", "is_doublet", "donor", "nullable_int", "nullable_bool", "with_missing_cat"], obs.ColumnNames);

        var cluster = Assert.IsType<Categorical>(obs["cluster"]);
        Assert.Equal(["B", "NK", "T"], (string[])cluster.Categories);
        Assert.False(cluster.Ordered);
        Assert.Equal(AnnFixtures.Expected.GetProperty("kitchen_cluster").EnumerateArray().Select(e => e.GetString()), cluster.ToStrings());
        var stage = Assert.IsType<Categorical>(obs["stage"]);
        Assert.True(stage.Ordered);
        Assert.Equal(["early", "mid", "late"], (string[])stage.Categories);
        var batch = Assert.IsType<Categorical>(obs["batch_num"]);
        Assert.Equal([1L, 2L], (long[])batch.Categories);
        var missing = Assert.IsType<Categorical>(obs["with_missing_cat"]);
        Assert.Equal(-1, missing.Codes[1]);
        Assert.Null(missing[1]);
        Assert.Equal("a", missing[0]);

        var rowSums = Enumerable.Range(0, 24).Select(i => Enumerable.Range(0, 15).Sum(j => a.X![i, j])).ToArray();
        Assert.Equal(rowSums, Assert.IsType<double[]>(obs["n_counts"]));
        Assert.IsType<int[]>(obs["n_genes"]);
        Assert.Equal(AnnFixtures.Doubles("kitchen_score"), ((float[])obs["score"]).Select(f => (double)f));
        Assert.Equal(AnnFixtures.Expected.GetProperty("kitchen_big").EnumerateArray().Select(e => e.GetInt64()), (long[])obs["big"]);
        Assert.Equal(AnnFixtures.Expected.GetProperty("kitchen_is_doublet").EnumerateArray().Select(e => e.GetBoolean()), (bool[])obs["is_doublet"]);
        Assert.IsType<Categorical>(obs["donor"]);

        var ni = Assert.IsType<MaskedArray>(obs["nullable_int"]);
        var expectedNi = AnnFixtures.Expected.GetProperty("kitchen_nullable_int").EnumerateArray().ToArray();
        for (int i = 0; i < 24; i++)
            if (expectedNi[i].ValueKind == JsonValueKind.Null) Assert.Null(ni[i]);
            else Assert.Equal(expectedNi[i].GetInt64(), (long)ni[i]!);
        var nb = Assert.IsType<MaskedArray>(obs["nullable_bool"]);
        var expectedNb = AnnFixtures.Expected.GetProperty("kitchen_nullable_bool").EnumerateArray().ToArray();
        for (int i = 0; i < 24; i++)
            if (expectedNb[i].ValueKind == JsonValueKind.Null) Assert.Null(nb[i]);
            else Assert.Equal(expectedNb[i].GetBoolean(), (bool)nb[i]!);
        Assert.True(double.IsNaN(obs.GetDoubles("nullable_int")[1]));

        // var: gene_symbol is a pandas-3 nullable string column.
        Assert.Equal("MT-SYM0", a.Var.GetStrings("gene_symbol")[0]);
        Assert.Equal("GENE3", a.Var.Index[3]);
        Assert.IsType<bool[]>(a.Var["highly_variable"]);
    }

    [Fact]
    public void KitchenSink_UnsNestedDictionaries()
    {
        var u = H5ad.Read(AnnFixtures.File("kitchen_sink.h5ad")).Uns;
        Assert.Equal("kitchen sink", u["title"]);
        Assert.Equal(3L, u["version"]);
        Assert.Equal(3.25, u["pi"]);
        Assert.Equal(true, u["flag"]);
        Assert.Equal([1L, 2L, 3L], (long[])u["ints"]!);
        Assert.Equal(["x", "y", "zed"], (string[])u["names"]!);
        Assert.Equal(["#1f77b4", "#ff7f0e", "#2ca02c"], (string[])u["cluster_colors"]!);
        var m = Assert.IsType<AnnArray>(u["matrix"]);
        Assert.Equal([2L, 3L], m.Shape);
        Assert.Equal([0.0, 1, 2, 3, 4, 5], (double[])m.Data);
        var neighbors = Assert.IsType<Dictionary<string, object?>>(u["neighbors"]);
        Assert.Equal("connectivities", neighbors["connectivities_key"]);
        var p = Assert.IsType<Dictionary<string, object?>>(neighbors["params"]);
        Assert.Equal(15L, p["n_neighbors"]);
        Assert.Equal("umap", p["method"]);
        var l2 = (Dictionary<string, object?>)((Dictionary<string, object?>)((Dictionary<string, object?>)u["nested"]!)["level1"]!)["level2"]!;
        Assert.Equal(1.5, l2["value"]);
        Assert.Equal("deep", l2["text"]);
    }

    [Fact]
    public void CompressedCsrAndCscMatchTheDenseValues()
    {
        var expected = AnnFixtures.Doubles("csr_dense");
        var csr = H5ad.Read(AnnFixtures.File("csr_gzip.h5ad"));
        Assert.Equal(AnnMatrixFormat.Csr, csr.X!.Format);
        AssertMatrix(expected, csr.X);
        Assert.Equal(["g1", "g2"], (string[])((Categorical)csr.Obs["group"]).Categories);
        var csc = H5ad.Read(AnnFixtures.File("csc.h5ad"));
        Assert.Equal(AnnMatrixFormat.Csc, csc.X!.Format);
        Assert.Equal("float64", csc.X.DataType);
        AssertMatrix(expected, csc.X);
    }

    [Fact]
    public void LegacyLayout_H5sparseAndCategoriesGroup()
    {
        var a = H5ad.Read(AnnFixtures.File("legacy_v01.h5ad"));
        Assert.Equal(6, a.ObsCount);
        Assert.Equal(5, a.VarCount);
        Assert.Equal(new string?[] { "u", "v", "u", "w", "v", "u" }, a.Obs.GetStrings("label"));
        var full = AnnFixtures.Doubles("csr_dense");
        for (int i = 0; i < 6; i++) for (int j = 0; j < 5; j++) Assert.Equal(full[i * 30 + j], a.X![i, j]);
    }

    [Fact]
    public void Pbmc3kSubset_ConvertsToSingleCellExperiment()
    {
        var sc = H5ad.ReadSingleCell(AnnFixtures.File("pbmc3k_subset.h5ad"));
        Assert.Equal(300, sc.CellCount);
        Assert.Equal(1000, sc.GeneCount);
        Assert.StartsWith("ENSG", H5ad.Read(AnnFixtures.File("pbmc3k_subset.h5ad")).Var.GetStrings("gene_ids")[0]);
        Assert.Contains(sc.GeneNames, g => g.StartsWith("MT-", StringComparison.Ordinal));
        Assert.All(sc.X.Values, v => Assert.Equal(Math.Round(v), v)); // raw UMI counts
    }

    // ---- writing ----

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Write_RoundTripsEveryElement(bool gzip)
    {
        var a = H5ad.Read(AnnFixtures.File("kitchen_sink.h5ad"));
        var path = AnnFixtures.TempPath("roundtrip.h5ad");
        H5ad.Write(path, a, new H5adWriteOptions { Compression = gzip ? H5adCompression.Gzip : H5adCompression.None, CompressionThreshold = 16, ChunkElements = 64 });
        var b = H5ad.Read(path);
        Assert.Equal(a.X!.ToDense().Data, b.X!.ToDense().Data);
        Assert.Equal(a.X.DataType, b.X.DataType);
        Assert.Equal(a.Obs.ColumnNames, b.Obs.ColumnNames);
        Assert.Equal(a.Obs.Index, b.Obs.Index);
        Assert.Equal("barcode", b.Obs.IndexName);
        foreach (var c in a.Obs.ColumnNames) Assert.Equal(a.Obs.GetStrings(c), b.Obs.GetStrings(c));
        Assert.Equal(((Categorical)a.Obs["stage"]).Ordered, ((Categorical)b.Obs["stage"]).Ordered);
        Assert.Equal(((MaskedArray)a.Obs["nullable_int"]).Mask, ((MaskedArray)b.Obs["nullable_int"]).Mask);
        foreach (var c in a.Var.ColumnNames) Assert.Equal(a.Var.GetStrings(c), b.Var.GetStrings(c));
        foreach (var k in a.Obsm.Keys) Assert.Equal(((AnnMatrix)a.Obsm[k]).ToDense().Data, ((AnnMatrix)b.Obsm[k]).ToDense().Data);
        Assert.Equal("float32", ((AnnMatrix)b.Obsm["X_umap"]).DataType);
        AssertSame(a.Obsp["distances"].Sparse!, b.Obsp["distances"].Sparse!);
        AssertSame(a.Layers["counts"].Sparse!, b.Layers["counts"].Sparse!);
        Assert.Equal("int32", b.Layers["counts"].DataType);
        Assert.Equal(a.Raw!.X.ToDense().Data, b.Raw!.X.ToDense().Data);
        Assert.Equal(a.Uns["title"], b.Uns["title"]);
        Assert.Equal(a.Uns["version"], b.Uns["version"]);
        Assert.Equal(true, b.Uns["flag"]);
        Assert.Equal(((AnnArray)a.Uns["matrix"]!).Data, ((AnnArray)b.Uns["matrix"]!).Data);
        var l2 = (Dictionary<string, object?>)((Dictionary<string, object?>)((Dictionary<string, object?>)b.Uns["nested"]!)["level1"]!)["level2"]!;
        Assert.Equal("deep", l2["text"]);
    }

    [Fact]
    public void Write_CscStaysCsc()
    {
        var a = H5ad.Read(AnnFixtures.File("csc.h5ad"));
        var path = AnnFixtures.TempPath("csc.h5ad");
        H5ad.Write(path, a);
        var b = H5ad.Read(path);
        Assert.Equal(AnnMatrixFormat.Csc, b.X!.Format);
        Assert.Equal("float64", b.X.DataType);
        AssertSame(a.X!.Sparse!, b.X.Sparse!);
    }

    [Fact]
    public void SingleCellExperiment_RoundTripsThroughH5ad()
    {
        var sc = H5ad.ReadSingleCell(AnnFixtures.File("pbmc3k_subset.h5ad"));
        Preprocessing.CalculateQcMetrics(sc);
        sc.CellLabels["half"] = Enumerable.Range(0, sc.CellCount).Select(i => i < 150 ? "first" : "second").ToArray();
        sc.Embeddings["X_test"] = new DenseMatrix(sc.CellCount, 2, Enumerable.Range(0, sc.CellCount * 2).Select(i => i * 0.5).ToArray());
        var path = AnnFixtures.TempPath("sc.h5ad");
        H5ad.Write(path, sc, new H5adWriteOptions { Compression = H5adCompression.Gzip });
        var back = H5ad.ReadSingleCell(path);
        AssertSame(sc.X, back.X);
        Assert.Equal(sc.CellNames, back.CellNames);
        Assert.Equal(sc.GeneIds, back.GeneIds);
        Assert.Equal(sc.CellMetrics[Preprocessing.PctMito], back.CellMetrics[Preprocessing.PctMito]);
        Assert.Equal(sc.CellLabels["half"], back.CellLabels["half"]);
        Assert.Equal(sc.GeneFlags["mt"], back.GeneFlags["mt"]);
        Assert.Equal(sc.Embeddings["X_test"].Data, back.Embeddings["X_test"].Data);
    }

    // ---- backed mode ----

    [Fact]
    public void Backed_CsrRowSlicesMatchTheInMemoryMatrix()
    {
        var mem = H5ad.Read(AnnFixtures.File("pbmc3k_subset.h5ad")).X!.Sparse!;
        using var backed = H5ad.OpenBacked(AnnFixtures.File("pbmc3k_subset.h5ad"));
        Assert.Equal(AnnMatrixFormat.Csr, backed.X!.Format);
        Assert.Equal(mem.NonZeroCount, backed.X.NonZeroCount);
        Assert.Equal(300, backed.ObsCount);
        AssertSame(mem.ReadRows(37, 101), backed.X.ReadRows(37, 101));
        AssertSame(mem.ReadRows(299, 1), backed.X.ReadRows(299, 1));
        AssertSame(mem.Transpose().ReadRows(10, 25).Transpose(), backed.X.ReadColumns(10, 25));
        var blocks = backed.X.RowBlocks(64).Select(b => b.Block).ToList();
        Assert.Equal(5, blocks.Count);
        Assert.Equal(mem.NonZeroCount, blocks.Sum(b => b.NonZeroCount));
    }

    [Fact]
    public void Backed_CscAndDenseSlices()
    {
        var expected = AnnFixtures.Doubles("csr_dense");
        using (var csc = H5ad.OpenBacked(AnnFixtures.File("csc.h5ad")))
        {
            var cols = csc.X!.ReadColumns(3, 4);
            var rows = csc.X.ReadRows(5, 7);
            for (int i = 0; i < 40; i++) for (int j = 0; j < 4; j++) Assert.Equal(expected[i * 30 + 3 + j], cols[i, j]);
            for (int i = 0; i < 7; i++) for (int j = 0; j < 30; j++) Assert.Equal(expected[(5 + i) * 30 + j], rows[i, j]);
        }
        using var dense = H5ad.OpenBacked(AnnFixtures.File("kitchen_sink.h5ad"));
        var x = AnnFixtures.Doubles("kitchen_X");
        var r = dense.X!.ReadRows(4, 3);
        for (int i = 0; i < 3; i++) for (int j = 0; j < 15; j++) Assert.Equal(x[(4 + i) * 15 + j], r[i, j]);
        Assert.Contains("counts", dense.LayerNames);
        Assert.Equal(x[16], dense.Layer("counts").ReadRows(1, 1)[0, 1]);
        Assert.Equal("kitchen sink", dense.Uns["title"]);
        Assert.Equal(AnnFixtures.Doubles("kitchen_X_pca"), ((AnnMatrix)dense.ReadObsm("X_pca")).ToDense().Data);
    }

    [Fact]
    public void WriteStreaming_ExportsAnOutOfCoreMatrix()
    {
        var a = H5ad.Read(AnnFixtures.File("pbmc3k_subset.h5ad"));
        var x = a.X!.Sparse!;
        var path = AnnFixtures.TempPath("streamed.h5ad");
        var log = x.Clone();
        for (int k = 0; k < log.NonZeroCount; k++) log.Values[k] = MathF.Log(1 + log.Values[k]);
        H5ad.WriteStreaming(path, x, a.Obs, a.Var, new Dictionary<string, ISparseMatrixSource> { ["log"] = log },
            uns: new Dictionary<string, object?> { ["note"] = "streamed" },
            options: new H5adWriteOptions { Compression = H5adCompression.Gzip, ChunkElements = 4096 }, blockRows: 37);
        var b = H5ad.Read(path);
        AssertSame(x, b.X!.Sparse!);
        AssertSame(log, b.Layers["log"].Sparse!);
        Assert.Equal(a.Obs.Index, b.Obs.Index);
        Assert.Equal("streamed", b.Uns["note"]);
    }

    // ---- memory-mapped store ----

    [Fact]
    public void MappedStore_RoundTripsAndWritesLayers()
    {
        var x = H5ad.Read(AnnFixtures.File("pbmc3k_subset.h5ad")).X!.Sparse!;
        var dir = AnnFixtures.TempPath("store.bnsparse");
        using (var m = MappedSparseMatrix.Create(dir, x, blockRows: 50))
        {
            Assert.Equal(x.NonZeroCount, m.NonZeroCount);
            AssertSame(x, m.ReadRows(0, m.Rows));
            AssertSame(x.ReadRows(123, 45), m.ReadRows(123, 45));
            Assert.Equal(x.RowIndices(7).ToArray(), m.RowIndices(7).ToArray());
            using var doubled = m.WriteLayer("doubled", (_, b) => { for (int k = 0; k < b.NonZeroCount; k++) b.Values[k] *= 2; }, blockRows: 33);
            Assert.Equal(x.Values.Select(v => 2 * v), doubled.ReadRows(0, doubled.Rows).Values);
        }
        using var reopened = MappedSparseMatrix.Open(dir, "doubled");
        Assert.Equal(["X", "doubled"], reopened.LayerNames);
        Assert.Equal(2 * x[5, x.RowIndices(5)[0]], reopened.ReadRows(5, 1)[0, x.RowIndices(5)[0]]);
    }

    [Fact]
    public void MappedStore_TransposesCscOutOfCore()
    {
        var expected = AnnFixtures.Doubles("csr_dense");
        using var csc = H5ad.OpenBacked(AnnFixtures.File("csc.h5ad"));
        using var m = MappedSparseMatrix.Create(AnnFixtures.TempPath("csc.bnsparse"), csc.X!);
        var all = m.ReadRows(0, m.Rows);
        for (int i = 0; i < 40; i++) for (int j = 0; j < 30; j++) Assert.Equal(expected[i * 30 + j], all[i, j]);
        for (int i = 0; i < 40; i++) Assert.True(all.RowIndices(i).ToArray().SequenceEqual(all.RowIndices(i).ToArray().Order()));
    }
}
