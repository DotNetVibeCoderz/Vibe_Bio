using BenchmarkDotNet.Attributes;
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;

namespace Biotech.Net.Benchmarks;

/// <summary>DESeq-style median-of-ratios normalisation, TMM and log-CPM on a 20,000 gene x 12 sample count matrix.</summary>
public class TranscriptomicsBenchmarks
{
    public const int Genes = 20_000, Samples = 12;
    private Dataset _counts = null!;
    private double[] _sizeFactors = [];

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(3);
        var genes = Enumerable.Range(0, Genes).Select(i => $"gene{i}").ToArray();
        var samples = Enumerable.Range(0, Samples).Select(i => $"S{i}").ToArray();
        var values = new double[Genes * Samples];
        for (int g = 0; g < Genes; g++)
        {
            double mean = Math.Exp(rng.NextDouble() * 8);
            for (int s = 0; s < Samples; s++)
            {
                double depth = 0.6 + 0.08 * s;
                values[g * Samples + s] = Math.Round(mean * depth * (0.7 + 0.6 * rng.NextDouble()));
            }
        }
        _counts = new Dataset("counts", genes, samples, values);
        _sizeFactors = Normalization.MedianOfRatiosSizeFactors(_counts);
    }

    [Benchmark(Description = "Median-of-ratios size factors")]
    public double[] SizeFactors() => Normalization.MedianOfRatiosSizeFactors(_counts);

    [Benchmark(Description = "Normalized counts (DESeq)")]
    public Dataset NormalizedCounts() => Normalization.NormalizedCounts(_counts, _sizeFactors);

    [Benchmark(Description = "TMM factors (edgeR)")]
    public double[] Tmm() => Normalization.TmmFactors(_counts);

    [Benchmark(Description = "log-CPM")]
    public Dataset LogCpm() => Normalization.LogCpm(_counts);
}
