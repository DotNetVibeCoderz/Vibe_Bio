using System.Diagnostics;
using Biotech.Net.Transcriptomics;
using Biotech.Net.Transcriptomics.SingleCell;

// Atlas-scale single-cell run (not a BenchmarkDotNet job: one pass over a synthetic matrix far larger than the
// managed heap, timing each out-of-core step and sampling the process memory while it runs).
//
//   dotnet run -c Release -p:SkipRustBuild=true -- --atlas [cells=1000000] [genes=30000] [umis=1500] [dir]
//
// The synthetic matrix is written to a .bnsparse store under [dir] (default: %TEMP%/biotechnet-atlas) and kept
// for re-runs with the same shape. Nothing is written to the repository.
// Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.

internal static class AtlasScaleBenchmark
{
    public static void Run(string[] args)
    {
        int cells = args.Length > 0 ? int.Parse(args[0]) : 1_000_000;
        int genes = args.Length > 1 ? int.Parse(args[1]) : 30_000;
        int umis = args.Length > 2 ? int.Parse(args[2]) : 1500;
        string root = args.Length > 3 ? args[3] : Path.Combine(Path.GetTempPath(), "biotechnet-atlas");
        string store = Path.Combine(root, $"synthetic_{cells}x{genes}_{umis}.bnsparse");
        Console.WriteLine($"Atlas-scale run: {cells:N0} cells x {genes:N0} genes, ~{umis} UMIs/cell, {Environment.ProcessorCount} threads, store {store}");

        using var sampler = new MemorySampler();
        MappedSparseMatrix counts;
        if (File.Exists(Path.Combine(store, "matrix.json")))
        {
            counts = MappedSparseMatrix.Open(store);
            Console.WriteLine($"reusing existing store ({counts.NonZeroCount:N0} stored entries)");
        }
        else counts = Step(sampler, "generate synthetic counts", () => Generate(store, cells, genes, umis));
        Console.WriteLine($"stored entries: {counts.NonZeroCount:N0} ({counts.NonZeroCount * 8 / 1e9:F1} GB of indices + values)");

        var names = Enumerable.Range(0, genes).Select(j => j < 13 ? $"MT-G{j}" : $"G{j}").ToArray();
        var qc = Step(sampler, "QC metrics (1 pass)", () => OutOfCore.CalculateQcMetrics(counts, names));
        var hvg = Step(sampler, "HVG seurat_v3, 2000 genes (2 passes)", () => OutOfCore.HighlyVariableGenes(counts, 2000, HvgFlavor.SeuratV3));
        var log = Step(sampler, "normalize_total + log1p -> new layer (1 pass)",
            () => OutOfCore.NormalizeTotalLog1pToLayer(counts, "log1p_norm", cellTotals: qc.TotalCounts));
        var pca = Step(sampler, "PCA covariance, 2000 HVG, 50 PCs (2 passes)", () => OutOfCorePca.Fit(log, hvg.Selected, components: 50));
        Console.WriteLine($"  PC1-5 variance ratio: {string.Join(", ", pca.ExplainedVarianceRatio.Take(5).Select(v => v.ToString("F4")))}");
        int sub = Math.Min(cells, 50_000);
        var subScores = new DenseMatrix(sub, 30);
        for (int i = 0; i < sub; i++) pca.Scores.Row(i)[..30].CopyTo(subScores.Row(i));
        var graph = Step(sampler, $"kNN graph (exact, brute force) on the first {sub:N0} cells, 30 PCs", () => NeighborGraph.Build(subScores, 15));
        var leiden = Step(sampler, $"Leiden on {sub:N0} cells", () => CommunityDetection.Cluster(graph.Connectivities));
        Console.WriteLine($"  {leiden.Distinct().Count()} clusters");

        int h5Cells = Math.Min(cells, 200_000);
        string h5 = Path.Combine(root, $"synthetic_{h5Cells}x{genes}.h5ad");
        var head = new RowRange(counts, h5Cells);
        var obs = new AnnDataFrame(Enumerable.Range(0, h5Cells).Select(i => $"cell{i}").ToArray());
        var var = new AnnDataFrame(names);
        Step(sampler, $"write h5ad (gzip, streamed) of the first {h5Cells:N0} cells", () =>
        {
            H5ad.WriteStreaming(h5, head, obs, var, options: new H5adWriteOptions { Compression = H5adCompression.Gzip });
            return 0;
        });
        Console.WriteLine($"  h5ad size {new FileInfo(h5).Length / 1e9:F2} GB");
        Step(sampler, $"QC metrics from the backed h5ad ({h5Cells:N0} cells)", () =>
        {
            using var backed = H5ad.OpenBacked(h5);
            return OutOfCore.CalculateQcMetrics(backed.X!, names);
        });
        log.Dispose();
        counts.Dispose();
        Console.WriteLine($"process peak working set {Process.GetCurrentProcess().PeakWorkingSet64 / 1e9:F2} GB, peak private (commit) {Process.GetCurrentProcess().PeakPagedMemorySize64 / 1e9:F2} GB");
    }

    private static T Step<T>(MemorySampler sampler, string name, Func<T> f)
    {
        GC.Collect();
        sampler.Reset();
        var sw = Stopwatch.StartNew();
        var r = f();
        sw.Stop();
        var (ws, priv, heap) = sampler.Peaks();
        Console.WriteLine($"{name,-70} {sw.Elapsed.TotalSeconds,8:F1} s   peak private {priv / 1e9,5:F2} GB   peak managed heap {heap / 1e9,5:F2} GB   peak working set {ws / 1e9,5:F2} GB");
        return r;
    }

    /// <summary>
    /// Ten cell types; each cell draws its library size from a log-normal and its UMIs from its type's gene
    /// distribution (Zipf-like baseline with 200 type-specific genes up-weighted 8x), via an alias table.
    /// </summary>
    private static MappedSparseMatrix Generate(string dir, int cells, int genes, int umis)
    {
        const int types = 10, block = 8192;
        var baseline = Enumerable.Range(0, genes).Select(j => 1.0 / Math.Pow(j + 10, 0.9)).ToArray();
        var perm = Enumerable.Range(0, genes).ToArray();
        new Random(7).Shuffle(perm);
        var tables = new Alias[types];
        for (int t = 0; t < types; t++)
        {
            var w = new double[genes];
            for (int j = 0; j < genes; j++) w[perm[j]] = baseline[j];
            var rnd = new Random(100 + t);
            for (int m = 0; m < 200; m++) w[rnd.Next(genes)] *= 8;
            tables[t] = new Alias(w);
        }
        using var writer = MappedSparseMatrixWriter.Create(dir, genes);
        for (int start = 0; start < cells; start += block)
        {
            int n = Math.Min(block, cells - start);
            var rows = new (int[] Cols, float[] Vals)[n];
            Parallel.For(0, n, () => new int[genes], (i, _, acc) =>
            {
                int cell = start + i;
                var rnd = new Random(cell * 31 + 17);
                int type = cell % types;
                double z = Math.Sqrt(-2 * Math.Log(1 - rnd.NextDouble())) * Math.Cos(2 * Math.PI * rnd.NextDouble());
                int lib = Math.Max(200, (int)(umis * Math.Exp(0.4 * z - 0.08)));
                var touched = new List<int>(lib);
                for (int u = 0; u < lib; u++)
                {
                    int g = tables[type].Sample(rnd);
                    if (acc[g]++ == 0) touched.Add(g);
                }
                touched.Sort();
                var vals = new float[touched.Count];
                for (int k = 0; k < touched.Count; k++) { vals[k] = acc[touched[k]]; acc[touched[k]] = 0; }
                rows[i] = (touched.ToArray(), vals);
                return acc;
            }, _ => { });
            foreach (var (c, v) in rows) writer.AppendRow(c, v);
        }
        return writer.Complete();
    }

    private sealed class Alias
    {
        private readonly double[] _prob;
        private readonly int[] _alias;

        public Alias(double[] w)
        {
            int n = w.Length;
            double s = w.Sum();
            _prob = new double[n]; _alias = new int[n];
            var p = w.Select(x => x * n / s).ToArray();
            var small = new Stack<int>(); var large = new Stack<int>();
            for (int i = 0; i < n; i++) (p[i] < 1 ? small : large).Push(i);
            while (small.Count > 0 && large.Count > 0)
            {
                int l = small.Pop(), g = large.Pop();
                _prob[l] = p[l]; _alias[l] = g;
                p[g] = p[g] + p[l] - 1;
                (p[g] < 1 ? small : large).Push(g);
            }
            while (large.Count > 0) _prob[large.Pop()] = 1;
            while (small.Count > 0) _prob[small.Pop()] = 1;
        }

        public int Sample(Random r)
        {
            int i = r.Next(_prob.Length);
            return r.NextDouble() < _prob[i] ? i : _alias[i];
        }
    }

    /// <summary>The first <c>rows</c> rows of a source.</summary>
    private sealed class RowRange(ISparseMatrixSource inner, int rows) : ISparseMatrixSource
    {
        public int Rows => rows;
        public int Columns => inner.Columns;
        public long NonZeroCount { get; } = inner is MappedSparseMatrix m ? m.RowPointer(rows) : inner.ReadRows(0, rows).NonZeroCount;
        public SparseMatrix ReadRows(int start, int count) => inner.ReadRows(start, count);
    }

    /// <summary>Samples working set, private bytes and the managed heap every 50 ms.</summary>
    private sealed class MemorySampler : IDisposable
    {
        private readonly Thread _thread;
        private volatile bool _stop;
        private long _ws, _priv, _heap;

        public MemorySampler()
        {
            _thread = new Thread(() =>
            {
                var p = Process.GetCurrentProcess();
                while (!_stop)
                {
                    p.Refresh();
                    Max(ref _ws, p.WorkingSet64);
                    Max(ref _priv, p.PrivateMemorySize64);
                    Max(ref _heap, GC.GetTotalMemory(false));
                    Thread.Sleep(50);
                }
            }) { IsBackground = true };
            _thread.Start();
        }

        private static void Max(ref long field, long v)
        {
            long cur;
            while (v > (cur = Interlocked.Read(ref field)) && Interlocked.CompareExchange(ref field, v, cur) != cur) { }
        }

        public void Reset() { Interlocked.Exchange(ref _ws, 0); Interlocked.Exchange(ref _priv, 0); Interlocked.Exchange(ref _heap, 0); }

        public (long WorkingSet, long Private, long Heap) Peaks()
        {
            Thread.Sleep(60);
            return (Interlocked.Read(ref _ws), Interlocked.Read(ref _priv), Interlocked.Read(ref _heap));
        }

        public void Dispose() { _stop = true; _thread.Join(); }
    }
}
