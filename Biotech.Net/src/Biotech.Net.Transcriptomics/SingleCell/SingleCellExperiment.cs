using System.Globalization;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>
/// An annotated single-cell matrix in the AnnData layout: cells in rows, genes in columns. Holds the current
/// expression matrix <see cref="X"/>, extra layers (e.g. raw counts), per-cell and per-gene annotations,
/// embeddings (PCA, UMAP, t-SNE) and the neighbour graph.
/// </summary>
public sealed class SingleCellExperiment
{
    /// <summary>Creates an experiment; gene names default to the gene ids.</summary>
    public SingleCellExperiment(SparseMatrix x, IReadOnlyList<string> cellNames, IReadOnlyList<string> geneIds, IReadOnlyList<string>? geneNames = null)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(cellNames);
        ArgumentNullException.ThrowIfNull(geneIds);
        if (cellNames.Count != x.Rows) throw new ArgumentException($"Need {x.Rows} cell names, got {cellNames.Count}.");
        if (geneIds.Count != x.Columns) throw new ArgumentException($"Need {x.Columns} gene ids, got {geneIds.Count}.");
        if (geneNames is not null && geneNames.Count != x.Columns) throw new ArgumentException("Gene names must match gene ids.");
        X = x;
        CellNames = cellNames;
        GeneIds = geneIds;
        GeneNames = geneNames ?? geneIds;
    }

    /// <summary>The current expression matrix (cells × genes). Preprocessing steps replace it.</summary>
    public SparseMatrix X { get; set; }

    /// <summary>Number of cells.</summary>
    public int CellCount => X.Rows;

    /// <summary>Number of genes.</summary>
    public int GeneCount => X.Columns;

    /// <summary>Cell barcodes.</summary>
    public IReadOnlyList<string> CellNames { get; private set; }

    /// <summary>Gene identifiers (e.g. Ensembl ids).</summary>
    public IReadOnlyList<string> GeneIds { get; private set; }

    /// <summary>Gene symbols.</summary>
    public IReadOnlyList<string> GeneNames { get; private set; }

    /// <summary>Additional matrices with the same shape as <see cref="X"/>, e.g. <c>counts</c>.</summary>
    public Dictionary<string, SparseMatrix> Layers { get; } = new(StringComparer.Ordinal);

    /// <summary>Numeric per-cell annotations (QC metrics, scores).</summary>
    public Dictionary<string, double[]> CellMetrics { get; } = new(StringComparer.Ordinal);

    /// <summary>Categorical per-cell annotations (cluster, cell type, sample).</summary>
    public Dictionary<string, string[]> CellLabels { get; } = new(StringComparer.Ordinal);

    /// <summary>Numeric per-gene annotations (means, dispersions).</summary>
    public Dictionary<string, double[]> GeneMetrics { get; } = new(StringComparer.Ordinal);

    /// <summary>Boolean per-gene flags (<c>highly_variable</c>, <c>mt</c>).</summary>
    public Dictionary<string, bool[]> GeneFlags { get; } = new(StringComparer.Ordinal);

    /// <summary>Cell embeddings by name (<c>X_pca</c>, <c>X_umap</c>, <c>X_tsne</c>), each cells × dimensions.</summary>
    public Dictionary<string, DenseMatrix> Embeddings { get; } = new(StringComparer.Ordinal);

    /// <summary>The k-nearest-neighbour graph, once computed.</summary>
    public NeighborGraph? Neighbors { get; set; }

    /// <summary>Index of a gene by symbol or id.</summary>
    public int GeneIndex(string gene)
    {
        for (int j = 0; j < GeneCount; j++) if (GeneNames[j] == gene || GeneIds[j] == gene) return j;
        throw new KeyNotFoundException($"No gene '{gene}'.");
    }

    /// <summary>Index of a gene, or −1.</summary>
    public int TryGeneIndex(string gene)
    {
        for (int j = 0; j < GeneCount; j++) if (GeneNames[j] == gene || GeneIds[j] == gene) return j;
        return -1;
    }

    /// <summary>Keeps only the given cells (in order), subsetting every per-cell structure.</summary>
    public void SubsetCells(ReadOnlySpan<int> cells)
    {
        X = X.SelectRows(cells);
        foreach (var k in Layers.Keys.ToList()) Layers[k] = Layers[k].SelectRows(cells);
        var idx = cells.ToArray();
        CellNames = idx.Select(i => CellNames[i]).ToArray();
        foreach (var k in CellMetrics.Keys.ToList()) CellMetrics[k] = idx.Select(i => CellMetrics[k][i]).ToArray();
        foreach (var k in CellLabels.Keys.ToList()) CellLabels[k] = idx.Select(i => CellLabels[k][i]).ToArray();
        foreach (var k in Embeddings.Keys.ToList())
        {
            var e = Embeddings[k];
            var n = new DenseMatrix(idx.Length, e.Columns);
            for (int r = 0; r < idx.Length; r++) e.Row(idx[r]).CopyTo(n.Row(r));
            Embeddings[k] = n;
        }
        Neighbors = null;
    }

    /// <summary>Keeps only the given genes (ascending indices), subsetting every per-gene structure.</summary>
    public void SubsetGenes(ReadOnlySpan<int> genes)
    {
        X = X.SelectColumns(genes);
        foreach (var k in Layers.Keys.ToList()) Layers[k] = Layers[k].SelectColumns(genes);
        var idx = genes.ToArray();
        GeneIds = idx.Select(i => GeneIds[i]).ToArray();
        GeneNames = idx.Select(i => GeneNames[i]).ToArray();
        foreach (var k in GeneMetrics.Keys.ToList()) GeneMetrics[k] = idx.Select(i => GeneMetrics[k][i]).ToArray();
        foreach (var k in GeneFlags.Keys.ToList()) GeneFlags[k] = idx.Select(i => GeneFlags[k][i]).ToArray();
    }

    /// <summary>A deep copy of the matrix and annotations (embeddings and graph are shared references).</summary>
    public SingleCellExperiment Copy()
    {
        var c = new SingleCellExperiment(X.Clone(), CellNames.ToArray(), GeneIds.ToArray(), GeneNames.ToArray());
        foreach (var kv in Layers) c.Layers[kv.Key] = kv.Value.Clone();
        foreach (var kv in CellMetrics) c.CellMetrics[kv.Key] = (double[])kv.Value.Clone();
        foreach (var kv in CellLabels) c.CellLabels[kv.Key] = (string[])kv.Value.Clone();
        foreach (var kv in GeneMetrics) c.GeneMetrics[kv.Key] = (double[])kv.Value.Clone();
        foreach (var kv in GeneFlags) c.GeneFlags[kv.Key] = (bool[])kv.Value.Clone();
        foreach (var kv in Embeddings) c.Embeddings[kv.Key] = kv.Value;
        c.Neighbors = Neighbors;
        return c;
    }

    /// <inheritdoc />
    public override string ToString()
        => $"SingleCellExperiment: {CellCount} cells × {GeneCount} genes; labels [{string.Join(", ", CellLabels.Keys)}]; embeddings [{string.Join(", ", Embeddings.Keys)}]";
}

/// <summary>
/// Matrix Market (<c>.mtx</c>) coordinate-format I/O, including the 10x Genomics Cell Ranger layout
/// (<c>matrix.mtx</c> genes × cells with <c>barcodes.tsv</c> and <c>features.tsv</c>/<c>genes.tsv</c>, optionally gzipped).
/// </summary>
public static class MatrixMarket
{
    /// <summary>
    /// Reads a coordinate Matrix Market matrix (<c>integer</c>, <c>real</c> or <c>pattern</c>; <c>general</c> or
    /// <c>symmetric</c>) into CSR with the file's own orientation.
    /// </summary>
    public static SparseMatrix Read(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        string header = reader.ReadLine() ?? throw new InvalidDataException("Empty Matrix Market file.");
        var h = header.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        if (h.Length < 4 || !h[0].Equals("%%MatrixMarket", StringComparison.OrdinalIgnoreCase) || !h[1].Equals("matrix", StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("Missing '%%MatrixMarket matrix' header.");
        if (!h[2].Equals("coordinate", StringComparison.OrdinalIgnoreCase))
            throw new InvalidDataException("Only the coordinate (sparse) Matrix Market format is supported.");
        string field = h[3].ToLowerInvariant();
        if (field == "complex") throw new InvalidDataException("Complex matrices are not supported.");
        bool pattern = field == "pattern";
        bool symmetric = h.Length > 4 && h[4].Equals("symmetric", StringComparison.OrdinalIgnoreCase);
        string? line;
        do line = reader.ReadLine(); while (line is not null && (line.Length == 0 || line[0] == '%'));
        if (line is null) throw new InvalidDataException("Matrix Market size line missing.");
        var size = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        int rows = int.Parse(size[0], CultureInfo.InvariantCulture), cols = int.Parse(size[1], CultureInfo.InvariantCulture);
        long nnz = long.Parse(size[2], CultureInfo.InvariantCulture);
        int cap = (int)Math.Min(nnz * (symmetric ? 2 : 1), int.MaxValue);
        var ri = new List<int>(cap); var ci = new List<int>(cap); var vs = new List<float>(cap);
        long read = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '%') continue;
            var span = line.AsSpan();
            int r = ParseIntToken(ref span) - 1, c = ParseIntToken(ref span) - 1;
            float v = pattern ? 1f : float.Parse(NextToken(ref span), NumberStyles.Float, CultureInfo.InvariantCulture);
            ri.Add(r); ci.Add(c); vs.Add(v);
            if (symmetric && r != c) { ri.Add(c); ci.Add(r); vs.Add(v); }
            read++;
        }
        if (read != nnz) throw new InvalidDataException($"Matrix Market header announces {nnz} entries but {read} were read.");
        return SparseMatrix.FromTriplets(rows, cols, ri.ToArray(), ci.ToArray(), vs.ToArray());
    }

    /// <summary>Writes a matrix in coordinate <c>real general</c> format (1-based indices).</summary>
    public static void Write(TextWriter writer, SparseMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(writer);
        ArgumentNullException.ThrowIfNull(matrix);
        writer.WriteLine("%%MatrixMarket matrix coordinate real general");
        writer.WriteLine($"{matrix.Rows} {matrix.Columns} {matrix.NonZeroCount}");
        for (int i = 0; i < matrix.Rows; i++)
        {
            var idx = matrix.RowIndices(i); var v = matrix.RowValues(i);
            for (int k = 0; k < idx.Length; k++)
                writer.WriteLine($"{i + 1} {idx[k] + 1} {v[k].ToString("R", CultureInfo.InvariantCulture)}");
        }
    }

    /// <summary>
    /// Reads a 10x-style triple. With <paramref name="genesInRows"/> (Cell Ranger convention) the matrix is
    /// transposed to cells × genes. Features: first column id, second column symbol (if present).
    /// </summary>
    public static SingleCellExperiment ReadSingleCell(TextReader matrix, TextReader barcodes, TextReader features, bool genesInRows = true)
    {
        var m = Read(matrix);
        if (genesInRows) m = m.Transpose();
        var cells = ReadLines(barcodes).Select(l => l.Split('\t')[0]).ToArray();
        var ids = new List<string>(); var names = new List<string>();
        foreach (var l in ReadLines(features))
        {
            var p = l.Split('\t');
            ids.Add(p[0]); names.Add(p.Length > 1 ? p[1] : p[0]);
        }
        var sc = new SingleCellExperiment(m, cells, ids, names);
        sc.Layers["counts"] = m;
        return sc;
    }

    /// <summary>
    /// Reads a Cell Ranger output directory (<c>matrix.mtx[.gz]</c>, <c>barcodes.tsv[.gz]</c>, and
    /// <c>features.tsv[.gz]</c> or the older <c>genes.tsv</c>).
    /// </summary>
    public static SingleCellExperiment Read10x(string directory)
    {
        string Find(params string[] names)
        {
            foreach (var n in names)
            {
                var p = Path.Combine(directory, n);
                if (File.Exists(p)) return p;
            }
            throw new FileNotFoundException($"None of {string.Join(", ", names)} found in {directory}.");
        }
        using var mtx = CountMatrixReader.OpenText(Find("matrix.mtx.gz", "matrix.mtx"));
        using var bc = CountMatrixReader.OpenText(Find("barcodes.tsv.gz", "barcodes.tsv"));
        using var ft = CountMatrixReader.OpenText(Find("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"));
        return ReadSingleCell(mtx, bc, ft);
    }

    private static IEnumerable<string> ReadLines(TextReader r)
    {
        string? l;
        while ((l = r.ReadLine()) is not null) if (l.Length > 0) yield return l;
    }

    private static ReadOnlySpan<char> NextToken(ref ReadOnlySpan<char> s)
    {
        s = s.TrimStart();
        int e = 0;
        while (e < s.Length && !char.IsWhiteSpace(s[e])) e++;
        var tok = s[..e];
        s = s[e..];
        if (tok.IsEmpty) throw new InvalidDataException("Truncated Matrix Market entry.");
        return tok;
    }

    private static int ParseIntToken(ref ReadOnlySpan<char> s) => int.Parse(NextToken(ref s), NumberStyles.Integer, CultureInfo.InvariantCulture);
}
