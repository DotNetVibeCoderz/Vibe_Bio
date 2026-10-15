using BenchmarkDotNet.Attributes;
using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Benchmarks;

/// <summary>
/// Rust (<c>NativeBridge</c>) versus managed C# (<c>ManagedKernels</c>) for every hot-path kernel.
/// The managed method of each category is the baseline, so Ratio reads as "Rust time / managed time".
/// </summary>
public class KernelBenchmarks
{
    public const int DnaLength = 1_000_000;
    public const int Genes = 20_000, Samples = 12;
    public const int PValueCount = 100_000;
    public const int ProteinLength = 10_000;
    public const int AlignLength = 2_000;
    public const int BatchTargets = 500, BatchLength = 300;
    public const int KabschPoints = 10_000;

    private byte[] _dna = [], _dnaOut = [], _protOut = [], _codonTable = [];
    private int[] _kmerCounts = [];
    private double[] _matrix = [], _t = [], _p = [], _diff = [], _pValues = [], _qValues = [];
    private int[] _groupA = [], _groupB = [];
    private byte[] _protein = [];
    private byte[] _alignA = [], _alignB = [], _ops = [];
    private int[] _dnaScoring = [], _proteinScoring = [];
    private readonly int[] _summary = new int[6];
    private byte[] _query = [], _targets = [];
    private int[] _offsets = [], _scores = [];
    private double[] _coordsA = [], _coordsB = [];
    private readonly double[] _rotation = new double[9];

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(42);
        const string acgt = "ACGT";
        const string aa = "ACDEFGHIKLMNPQRSTVWY";

        _dna = new byte[DnaLength];
        for (int i = 0; i < _dna.Length; i++) _dna[i] = (byte)acgt[rng.Next(4)];
        _dnaOut = new byte[DnaLength];
        _protOut = new byte[DnaLength / 3];
        _codonTable = GeneticCode.Standard.CodonTable.ToArray();
        _kmerCounts = new int[1 << 16];

        _matrix = new double[Genes * Samples];
        for (int i = 0; i < _matrix.Length; i++)
            _matrix[i] = 5 + rng.NextDouble() * 3 + ((i % Samples) >= Samples / 2 && (i / Samples) % 10 == 0 ? 1.5 : 0);
        _groupA = Enumerable.Range(0, Samples / 2).ToArray();
        _groupB = Enumerable.Range(Samples / 2, Samples / 2).ToArray();
        _t = new double[Genes]; _p = new double[Genes]; _diff = new double[Genes];

        _pValues = new double[PValueCount];
        for (int i = 0; i < _pValues.Length; i++) _pValues[i] = rng.NextDouble() * (i % 20 == 0 ? 1e-4 : 1);
        _qValues = new double[PValueCount];

        _protein = new byte[ProteinLength];
        for (int i = 0; i < _protein.Length; i++) _protein[i] = (byte)aa[rng.Next(aa.Length)];

        _alignA = new byte[AlignLength];
        for (int i = 0; i < AlignLength; i++) _alignA[i] = (byte)acgt[rng.Next(4)];
        _alignB = (byte[])_alignA.Clone();
        for (int i = 0; i < AlignLength / 10; i++) _alignB[rng.Next(AlignLength)] = (byte)acgt[rng.Next(4)];
        _ops = new byte[AlignLength * 2];
        _dnaScoring = BioKernels.MatchMismatchMatrix(2, -3);
        _proteinScoring = BioKernels.MatchMismatchMatrix(5, -2);

        _query = new byte[BatchLength];
        for (int i = 0; i < BatchLength; i++) _query[i] = (byte)aa[rng.Next(aa.Length)];
        _targets = new byte[BatchTargets * BatchLength];
        for (int i = 0; i < _targets.Length; i++) _targets[i] = (byte)aa[rng.Next(aa.Length)];
        _offsets = Enumerable.Range(0, BatchTargets + 1).Select(k => k * BatchLength).ToArray();
        _scores = new int[BatchTargets];

        _coordsA = new double[KabschPoints * 3];
        for (int i = 0; i < _coordsA.Length; i++) _coordsA[i] = rng.NextDouble() * 50;
        _coordsB = new double[_coordsA.Length];
        double c = Math.Cos(0.7), s = Math.Sin(0.7);
        for (int i = 0; i < KabschPoints; i++)
        {
            double x = _coordsA[3 * i], y = _coordsA[3 * i + 1], z = _coordsA[3 * i + 2];
            _coordsB[3 * i] = c * x - s * y + 3 + rng.NextDouble() * 0.1;
            _coordsB[3 * i + 1] = s * x + c * y - 2;
            _coordsB[3 * i + 2] = z + 1;
        }

        if (!NativeBridge.IsAvailable)
            throw new InvalidOperationException($"Rust kernels unavailable: {NativeBridge.UnavailableReason}");
    }

    // ---- Sequence engine ---------------------------------------------------------------------------

    [BenchmarkCategory("ReverseComplement 1 Mb"), Benchmark(Baseline = true)]
    public void ReverseComplement_Managed() => ManagedKernels.ReverseComplement(_dna, _dnaOut);

    [BenchmarkCategory("ReverseComplement 1 Mb"), Benchmark]
    public void ReverseComplement_Rust() => NativeBridge.ReverseComplement(_dna, _dnaOut);

    [BenchmarkCategory("Translate 1 Mb"), Benchmark(Baseline = true)]
    public int Translate_Managed() => ManagedKernels.Translate(_dna, _codonTable, _protOut);

    [BenchmarkCategory("Translate 1 Mb"), Benchmark]
    public int Translate_Rust() => NativeBridge.Translate(_dna, _codonTable, _protOut);

    [BenchmarkCategory("CountKmers k=8, 1 Mb"), Benchmark(Baseline = true)]
    public long CountKmers_Managed() { Array.Clear(_kmerCounts); return ManagedKernels.CountKmers(_dna, 8, _kmerCounts); }

    [BenchmarkCategory("CountKmers k=8, 1 Mb"), Benchmark]
    public long CountKmers_Rust() { Array.Clear(_kmerCounts); return NativeBridge.CountKmers(_dna, 8, _kmerCounts); }

    // ---- Statistics engine -------------------------------------------------------------------------

    [BenchmarkCategory("Welch t-test 20k x 12"), Benchmark(Baseline = true)]
    public void WelchTTest_Managed() => ManagedKernels.WelchTTest(_matrix, Genes, Samples, _groupA, _groupB, _t, _p, _diff);

    [BenchmarkCategory("Welch t-test 20k x 12"), Benchmark]
    public void WelchTTest_Rust() => NativeBridge.WelchTTest(_matrix, Genes, Samples, _groupA, _groupB, _t, _p, _diff);

    [BenchmarkCategory("Benjamini-Hochberg 100k"), Benchmark(Baseline = true)]
    public void BenjaminiHochberg_Managed() => ManagedKernels.BenjaminiHochberg(_pValues, _qValues);

    [BenchmarkCategory("Benjamini-Hochberg 100k"), Benchmark]
    public void BenjaminiHochberg_Rust() => NativeBridge.BenjaminiHochberg(_pValues, _qValues);

    // ---- Protein engine ----------------------------------------------------------------------------

    [BenchmarkCategory("Protein mass 10k aa"), Benchmark(Baseline = true)]
    public double ProteinMass_Managed() => ManagedKernels.ProteinMass(_protein, false);

    [BenchmarkCategory("Protein mass 10k aa"), Benchmark]
    public double ProteinMass_Rust() => NativeBridge.ProteinMass(_protein, false);

    [BenchmarkCategory("Isoelectric point 10k aa"), Benchmark(Baseline = true)]
    public double IsoelectricPoint_Managed() => ManagedKernels.IsoelectricPoint(_protein);

    [BenchmarkCategory("Isoelectric point 10k aa"), Benchmark]
    public double IsoelectricPoint_Rust() => NativeBridge.IsoelectricPoint(_protein);

    // ---- Alignment engine --------------------------------------------------------------------------

    [BenchmarkCategory("Gotoh global 2 kb x 2 kb"), Benchmark(Baseline = true)]
    public int GotohGlobal_Managed()
        => ManagedKernels.Align(_alignA, _alignB, _dnaScoring, 5, 2, (int)AlignmentMode.Global, -1, _ops, _summary);

    [BenchmarkCategory("Gotoh global 2 kb x 2 kb"), Benchmark]
    public int GotohGlobal_Rust()
        => NativeBridge.Align(_alignA, _alignB, _dnaScoring, 5, 2, (int)AlignmentMode.Global, -1, _ops, _summary);

    [BenchmarkCategory("Local score batch 300 aa x 500"), Benchmark(Baseline = true)]
    public void LocalScoreBatch_Managed()
        => ManagedKernels.LocalScoreBatch(_query, _targets, _offsets, _proteinScoring, 11, 1, _scores);

    [BenchmarkCategory("Local score batch 300 aa x 500"), Benchmark]
    public void LocalScoreBatch_Rust()
        => NativeBridge.LocalScoreBatch(_query, _targets, _offsets, _proteinScoring, 11, 1, _scores);

    // ---- Structure engine --------------------------------------------------------------------------

    [BenchmarkCategory("Kabsch RMSD 10k atoms"), Benchmark(Baseline = true)]
    public double Kabsch_Managed() => ManagedKernels.KabschRmsd(_coordsA, _coordsB, KabschPoints, _rotation);

    [BenchmarkCategory("Kabsch RMSD 10k atoms"), Benchmark]
    public double Kabsch_Rust() => NativeBridge.KabschRmsd(_coordsA, _coordsB, KabschPoints, _rotation);
}
