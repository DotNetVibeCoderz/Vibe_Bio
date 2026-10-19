using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Search;
using Biotech.Net.Core;
using Biotech.Net.Gpu.Alignment;
using Biotech.Net.Native;

namespace Biotech.Net.Gpu.Tests;

public class AllVsAllTests
{
    private const string Amino = "ACDEFGHIKLMNPQRSTVWY";

    private static List<Sequence> RandomProteins(int count, int minLen, int maxLen, int seed, string prefix = "p")
    {
        var rng = new Random(seed);
        var list = new List<Sequence>();
        for (int i = 0; i < count; i++)
        {
            int len = rng.Next(minLen, maxLen + 1);
            var chars = new char[len];
            for (int j = 0; j < len; j++) chars[j] = Amino[rng.Next(Amino.Length)];
            list.Add(Sequence.Protein(new string(chars), $"{prefix}{i}"));
        }
        return list;
    }

    /// <summary>Mutated copies of a parent, so scores span weak to strong similarity.</summary>
    private static List<Sequence> Family(Sequence parent, int count, int seed)
    {
        var rng = new Random(seed);
        var list = new List<Sequence>();
        string p = parent.ToString();
        for (int i = 0; i < count; i++)
        {
            var sb = new System.Text.StringBuilder();
            foreach (char c in p)
            {
                double u = rng.NextDouble();
                if (u < 0.1) sb.Append(Amino[rng.Next(Amino.Length)]);
                else if (u < 0.13) continue;                             // deletion
                else if (u < 0.16) sb.Append(c).Append(Amino[rng.Next(Amino.Length)]).Append(Amino[rng.Next(Amino.Length)]);
                else sb.Append(c);
            }
            list.Add(Sequence.Protein(sb.ToString(), $"{parent.Id}_m{i}"));
        }
        return list;
    }

    private static void AssertMatchesCpu(LocalScoreMatrix m, IReadOnlyList<Sequence> q, IReadOnlyList<Sequence> t, SubstitutionMatrix matrix, int open, int extend)
    {
        Assert.DoesNotContain("CPU (", m.Backend);   // the kernel ran; no silent fallback
        for (int i = 0; i < q.Count; i++)
            for (int j = 0; j < t.Count; j++)
                Assert.Equal(BioKernels.LocalScore(q[i].Span, t[j].Span, matrix.KernelMatrix, open, extend), m[i, j]);
    }

    [Fact]
    public void Device_enumeration_always_lists_the_cpu_accelerator()
    {
        var devices = GpuDevice.Enumerate();
        Assert.Contains(devices, d => d.Backend == GpuBackend.Cpu);
        Assert.False(GpuTestDevices.Cpu.IsGpu);
        Assert.Equal(GpuBackend.Cpu, GpuTestDevices.Cpu.Backend);
        Assert.Null(GpuDevice.TryCreate(GpuBackend.Off));
    }

    [Theory]
    [InlineData("cuda", GpuBackend.Cuda)]
    [InlineData("OpenCL", GpuBackend.OpenCL)]
    [InlineData(" cpu ", GpuBackend.Cpu)]
    [InlineData("off", GpuBackend.Off)]
    [InlineData("0", GpuBackend.Off)]
    [InlineData("auto", GpuBackend.Auto)]
    [InlineData(null, GpuBackend.Auto)]
    [InlineData("vulkan", GpuBackend.Auto)]
    public void Backend_names_parse(string? value, GpuBackend expected) => Assert.Equal(expected, GpuDevice.ParseBackend(value));

    [Fact]
    public void Disposed_device_rejects_work()
    {
        var d = GpuDevice.CreateCpu();
        d.Dispose();
        d.Dispose();
        var aligner = new GpuLocalAligner(d);
        Assert.Throws<ObjectDisposedException>(() => aligner.ScoreAll(RandomProteins(2, 5, 9, 1), RandomProteins(2, 5, 9, 2)));
    }

    [Fact]
    public void Protein_all_vs_all_is_bit_identical_to_the_cpu_kernel()
    {
        var parent = RandomProteins(3, 80, 160, 11);
        var seqs = parent.SelectMany((p, i) => Family(p, 6, 100 + i).Prepend(p)).Concat(RandomProteins(10, 1, 220, 12, "r")).ToList();
        var aligner = new GpuLocalAligner(GpuTestDevices.Cpu, SubstitutionMatrix.Blosum62, 11, 1);
        var m = aligner.ScoreAll(seqs);
        Assert.Null(aligner.LastFallbackError?.ToString());
        Assert.StartsWith("Cpu:", m.Backend);
        AssertMatchesCpu(m, seqs, seqs, SubstitutionMatrix.Blosum62, 11, 1);
        // family members score far above unrelated sequences
        Assert.True(m[0, 1] > 3 * m[0, seqs.Count - 1]);
    }

    [Theory]
    [InlineData("BLOSUM45", 15, 2)]
    [InlineData("PAM30", 9, 1)]
    [InlineData("BLOSUM80", 3, 3)]
    [InlineData("BLOSUM62", 0, 0)]
    public void Other_matrices_and_gap_costs_agree(string matrix, int open, int extend)
    {
        var q = RandomProteins(7, 0, 90, 21);
        var t = RandomProteins(13, 0, 120, 22, "t");
        var sm = SubstitutionMatrix.Load(matrix);
        var m = new GpuLocalAligner(GpuTestDevices.Cpu, sm, open, extend).ScoreAll(q, t);
        AssertMatchesCpu(m, q, t, sm, open, extend);
    }

    [Fact]
    public void Nucleotides_mixed_case_and_odd_lengths_agree_with_managed_and_native()
    {
        var rng = new Random(5);
        string Dna(int n) => new([.. Enumerable.Range(0, n).Select(_ => "ACGTacgtN"[rng.Next(9)])]);
        var q = Enumerable.Range(0, 9).Select(i => Sequence.Dna(Dna(rng.Next(0, 70)), $"q{i}")).ToList();
        var t = Enumerable.Range(0, 11).Select(i => Sequence.Dna(Dna(rng.Next(1, 83)), $"t{i}")).ToList();
        var sm = SubstitutionMatrix.MatchMismatch(2, -3);
        var m = new GpuLocalAligner(GpuTestDevices.Cpu, sm, 5, 2).ScoreAll(q, t);
        bool before = BioKernels.ForceManaged;
        try
        {
            foreach (bool managed in new[] { true, false })
            {
                BioKernels.ForceManaged = managed;
                AssertMatchesCpu(m, q, t, sm, 5, 2);
            }
        }
        finally { BioKernels.ForceManaged = before; }
    }

    [Fact]
    public void Small_scratch_budget_splits_the_launch_into_chunks()
    {
        var q = RandomProteins(5, 30, 60, 31);
        var t = RandomProteins(40, 10, 90, 32, "t");
        long saved = LocalScoreRunner.ScratchBudget;
        try
        {
            LocalScoreRunner.ScratchBudget = 7 * 60;   // 7 lanes per launch, 29 launches
            var m = new GpuLocalAligner(GpuTestDevices.Cpu).ScoreAll(q, t);
            AssertMatchesCpu(m, q, t, SubstitutionMatrix.Blosum62, 11, 1);
        }
        finally { LocalScoreRunner.ScratchBudget = saved; }
    }

    [Fact]
    public void Matches_the_pairwise_aligner_score_path()
    {
        var aligner = new PairwiseAligner { Mode = AlignmentMode.Local, MatchScore = 3, MismatchScore = -2, OpenGapScore = -5, ExtendGapScore = -1 };
        var q = RandomProteins(4, 20, 50, 41);
        var t = RandomProteins(6, 20, 50, 42, "t");
        var m = GpuLocalAligner.FromPairwiseAligner(aligner, GpuTestDevices.Cpu).ScoreAll(q, t);
        for (int i = 0; i < q.Count; i++)
            for (int j = 0; j < t.Count; j++) Assert.Equal(aligner.Score(q[i], t[j]), m[i, j]);
        Assert.Throws<ArgumentException>(() => GpuLocalAligner.FromPairwiseAligner(new PairwiseAligner(), null));
    }

    [Fact]
    public void Top_k_is_sorted_and_breaks_ties_by_index()
    {
        var q = RandomProteins(3, 40, 60, 51);
        var t = Family(q[0], 5, 52).Concat(RandomProteins(20, 40, 60, 53, "t")).Append(q[0]).Append(q[0]).ToList();
        var m = new GpuLocalAligner(GpuTestDevices.Cpu).ScoreAll(q, t);
        var top = m.Top(0, 4);
        Assert.Equal(4, top.Count);
        Assert.Equal(t.Count - 2, top[0].Target);     // identical copies first, lower index wins the tie
        Assert.Equal(t.Count - 1, top[1].Target);
        Assert.True(top.Zip(top.Skip(1)).All(p => p.First.Score >= p.Second.Score));
        var expected = Enumerable.Range(0, t.Count).OrderByDescending(j => m[0, j]).ThenBy(j => j).Take(4).ToList();
        Assert.Equal(expected, top.Select(h => h.Target));
        Assert.Empty(m.Top(1, 0));
        Assert.Equal(t.Count, m.Top(2, 1000).Count);
        Assert.Equal(3, m.TopAll(2).Count);
    }

    [Fact]
    public void Null_device_and_kernel_failure_fall_back_to_the_cpu()
    {
        var q = RandomProteins(3, 20, 40, 61);
        var t = RandomProteins(4, 20, 40, 62, "t");
        var cpu = new GpuLocalAligner(null).ScoreAll(q, t);
        Assert.StartsWith("CPU", cpu.Backend);
        var failing = new GpuLocalAligner(GpuTestDevices.Cpu) { InjectFailure = true };
        var m = failing.ScoreAll(q, t);
        Assert.StartsWith("CPU", m.Backend);
        Assert.IsType<InvalidOperationException>(failing.LastFallbackError);
        Assert.Equal(cpu.Values.ToArray(), m.Values.ToArray());
    }

    [Fact]
    public void Empty_inputs_give_empty_matrices()
    {
        var a = new GpuLocalAligner(GpuTestDevices.Cpu);
        Assert.Equal(0, a.ScoreAll([], RandomProteins(3, 4, 8, 1)).QueryCount);
        var m = a.ScoreAll(RandomProteins(2, 4, 8, 1), []);
        Assert.Equal(2, m.QueryCount);
        Assert.Equal(0, m.TargetCount);
        Assert.Throws<ArgumentException>(() => a.Score([1, 2], [0, 3], [1], [0, 1], SubstitutionMatrix.Blosum62.KernelMatrix, 1, 1, new int[1]));
    }

    [Fact]
    public void Database_search_with_the_gpu_backend_matches_the_cpu_search()
    {
        var parents = RandomProteins(4, 100, 150, 71);
        var db = new SequenceDatabase(parents.SelectMany((p, i) => Family(p, 5, 200 + i)).Concat(RandomProteins(40, 60, 200, 72, "d")));
        var query = parents[2];
        foreach (bool seeds in new[] { true, false })
        {
            var cpu = db.Search(query, new SearchOptions { UseSeeds = seeds, IncludeAlignments = false });
            var backend = new GpuLocalAligner(GpuTestDevices.Cpu);
            var gpu = db.Search(query, new SearchOptions { UseSeeds = seeds, IncludeAlignments = false, ScoreBackend = backend });
            Assert.Null(backend.LastFallbackError);
            Assert.StartsWith("Cpu:", backend.LastBackend);
            Assert.Equal(cpu.Select(h => (h.SubjectIndex, h.Score, h.EValue)), gpu.Select(h => (h.SubjectIndex, h.Score, h.EValue)));
            Assert.NotEmpty(gpu);
        }
    }

    [Fact]
    public void Real_gpu_matches_the_cpu_when_present()
    {
        if (GpuTestDevices.Gpu is not { } gpu) return;   // GPU-gated: no CUDA/OpenCL device here
        var seqs = RandomProteins(60, 1, 300, 81);
        var aligner = new GpuLocalAligner(gpu);
        var m = aligner.ScoreAll(seqs);
        Assert.Null(aligner.LastFallbackError);
        Assert.DoesNotContain("CPU", m.Backend);
        AssertMatchesCpu(m, seqs, seqs, SubstitutionMatrix.Blosum62, 11, 1);
    }
}
