namespace Biotech.Net.Workflows.Tests;

public class SchedulingTests
{
    /// <summary>
    /// <paramref name="count"/> independent jobs whose bodies block on <paramref name="gate"/>; records how many run at once.
    /// </summary>
    private sealed class Gated
    {
        private int _running;
        public int Max;
        public readonly TaskCompletionSource Gate = new(TaskCreationOptions.RunContinuationsAsynchronously);
        public readonly SemaphoreSlim Started = new(0);

        public async Task Body(JobContext ctx)
        {
            int now = Interlocked.Increment(ref _running);
            int seen;
            do seen = Max; while (now > seen && Interlocked.CompareExchange(ref Max, now, seen) != seen);
            Started.Release();
            await Gate.Task;
            Interlocked.Decrement(ref _running);
            await File.WriteAllTextAsync(ctx.Output("o"), ctx.JobId);
        }
    }

    private static (Workflow Wf, Gated Gated) Independent(TempDir dir, int count, Action<Rule>? configure = null)
    {
        var wf = new Workflow(dir.Path);
        var gated = new Gated();
        var rule = wf.Rule("work").Output("o", "out/{i}.txt").Run(gated.Body);
        configure?.Invoke(rule);
        wf.Target(Expand.Each("out/{i}.txt", "i", Enumerable.Range(0, count).Select(i => i.ToString())));
        return (wf, gated);
    }

    private static async Task<WorkflowRunReport> RunUntilStarted(Workflow wf, Gated gated, RunOptions options, int expectStarted)
    {
        var run = wf.RunAsync(options);
        for (int i = 0; i < expectStarted; i++)
            Assert.True(await gated.Started.WaitAsync(TimeSpan.FromSeconds(30)), "jobs did not start");
        gated.Gate.SetResult();
        return await run;
    }

    [Fact]
    public async Task MaxParallelism_IsRespected()
    {
        using var dir = new TempDir();
        var (wf, gated) = Independent(dir, 6);
        var report = await RunUntilStarted(wf, gated, new RunOptions { MaxParallelism = 2, MaxThreads = 16 }, 2);
        Assert.True(report.Succeeded);
        Assert.Equal(2, gated.Max);
        Assert.Equal(6, report.Count(JobStatus.Succeeded));
    }

    [Fact]
    public async Task ThreadAndCustomResources_AreRespected()
    {
        using var dir = new TempDir();
        var (wf, gated) = Independent(dir, 4, r => r.Threads(2));
        var report = await RunUntilStarted(wf, gated, new RunOptions { MaxParallelism = 8, MaxThreads = 5 }, 2);
        Assert.True(report.Succeeded);
        Assert.Equal(2, gated.Max);                         // 2 + 2 threads fit in 5, a third does not

        using var dir2 = new TempDir();
        var (wf2, gated2) = Independent(dir2, 3, r => r.Resource("gpu", 1).Memory(100));
        var report2 = await RunUntilStarted(wf2, gated2,
            new RunOptions { MaxParallelism = 8, MaxThreads = 8, Resources = new Dictionary<string, int> { ["gpu"] = 1 } }, 1);
        Assert.True(report2.Succeeded);
        Assert.Equal(1, gated2.Max);

        using var dir3 = new TempDir();
        var (wf3, gated3) = Independent(dir3, 3, r => r.Memory(600));
        var report3 = await RunUntilStarted(wf3, gated3, new RunOptions { MaxParallelism = 8, MaxThreads = 8, MaxMemoryMb = 1000 }, 1);
        Assert.True(report3.Succeeded);
        Assert.Equal(1, gated3.Max);
    }

    [Fact]
    public async Task OversizedRequests_AreCappedSoTheJobStillRuns()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        int threads = 0;
        wf.Rule("big").Output("o", "big.txt").Threads(64).Run(ctx =>
        {
            threads = ctx.Threads;
            File.WriteAllText(ctx.Output("o"), "x");
        });
        var report = await wf.RunAsync(new RunOptions { MaxThreads = 3 });
        Assert.True(report.Succeeded);
        Assert.Equal(3, threads);
    }

    [Fact]
    public async Task Dependencies_RunBeforeDependents_InTopologicalOrder()
    {
        using var dir = new TempDir();
        dir.Write("in.txt", "seed");
        var order = new List<string>();
        var wf = new Workflow(dir.Path);
        void Step(string name, string input, string output) => wf.Rule(name).Input("i", input).Output("o", output).Run(ctx =>
        {
            lock (order) order.Add(ctx.JobId);
            Assert.True(File.Exists(ctx.Input("i")), $"{ctx.JobId} started before its input existed");
            File.WriteAllText(ctx.Output("o"), File.ReadAllText(ctx.Input("i")) + "+" + name);
        });
        Step("c", "b.txt", "c.txt");
        Step("a", "in.txt", "a.txt");
        Step("b", "a.txt", "b.txt");
        var report = await wf.RunAsync(new RunOptions { Targets = ["c.txt"], MaxParallelism = 4 });
        Assert.True(report.Succeeded);
        Assert.Equal(["a", "b", "c"], order);
        Assert.Equal("seed+a+b+c", dir.Read("c.txt"));
        Assert.Equal(["a", "b", "c"], report.Jobs.Select(j => j.Id));
    }

    [Fact]
    public async Task Diamond_JoinWaitsForBothBranches()
    {
        using var dir = new TempDir();
        var left = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var leftStarted = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var wf = new Workflow(dir.Path);
        wf.Rule("root").Output("o", "root.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "r"));
        wf.Rule("left").Input("i", "root.txt").Output("o", "left.txt").Run(async ctx =>
        {
            leftStarted.SetResult();
            await left.Task;
            await File.WriteAllTextAsync(ctx.Output("o"), "l");
        });
        wf.Rule("right").Input("i", "root.txt").Output("o", "right.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "r"));
        bool joinSawBoth = false;
        wf.Rule("join").Input("l", "left.txt").Input("r", "right.txt").Output("o", "join.txt").Run(ctx =>
        {
            joinSawBoth = File.Exists(ctx.Input("l")) && File.Exists(ctx.Input("r"));
            File.WriteAllText(ctx.Output("o"), "j");
        });
        var events = new EventLog();
        var run = wf.RunAsync(new RunOptions { Targets = ["join.txt"] }, events);
        await leftStarted.Task.WaitAsync(TimeSpan.FromSeconds(30));
        Assert.DoesNotContain(events.Events, e => e is JobStartedEvent { JobId: "join" });
        left.SetResult();
        var report = await run;
        Assert.True(report.Succeeded);
        Assert.True(joinSawBoth);
    }
}
