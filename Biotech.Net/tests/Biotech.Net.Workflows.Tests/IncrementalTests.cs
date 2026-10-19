namespace Biotech.Net.Workflows.Tests;

public class IncrementalTests
{
    /// <summary>
    /// Two independent chains: data/{s}.txt → upper → work/{s}.up → count → results/{s}.n (length of the text),
    /// gathered by <c>summary</c>. Each body counts its runs in <paramref name="runs"/>.
    /// </summary>
    private static Workflow Pipeline(TempDir dir, Dictionary<string, int> runs, int suffixLength = 0, string version = "1")
    {
        void Count(string id)
        {
            lock (runs) runs[id] = runs.GetValueOrDefault(id) + 1;
        }

        var wf = new Workflow(dir.Path, "incremental");
        wf.Rule("upper").Input("raw", "data/{s}.txt").Output("up", "work/{s}.up").Version(version)
            .Run(ctx =>
            {
                Count(ctx.JobId);
                File.WriteAllText(ctx.Output("up"), File.ReadAllText(ctx.Input("raw")).ToUpperInvariant());
            });
        wf.Rule("count").Input("up", "work/{s}.up").Output("n", "results/{s}.n").Params(new { Suffix = suffixLength })
            .Run(ctx =>
            {
                Count(ctx.JobId);
                int suffix = ctx.Params<CountParams>().Suffix;
                File.WriteAllText(ctx.Output("n"), (File.ReadAllText(ctx.Input("up")).Length + suffix).ToString());
            });
        wf.Rule("summary").Input("all", Expand.Each("results/{s}.n", "s", ["a", "b"])).Output("o", "results/summary.txt")
            .Run(ctx =>
            {
                Count(ctx.JobId);
                File.WriteAllText(ctx.Output("o"), string.Join(",", ctx.Inputs("all").Select(File.ReadAllText)));
            });
        wf.Target("results/summary.txt");
        return wf;
    }

    private sealed record CountParams(int Suffix);

    private static TempDir Seed()
    {
        var dir = new TempDir();
        dir.Write("data/a.txt", "abc");
        dir.Write("data/b.txt", "hello");
        return dir;
    }

    [Fact]
    public async Task SecondRun_SkipsEverything()
    {
        using var dir = Seed();
        var runs = new Dictionary<string, int>();
        Assert.True((await Pipeline(dir, runs).RunAsync()).Succeeded);
        Assert.Equal(5, runs.Count);
        Assert.Equal("3,5", dir.Read("results/summary.txt"));

        var again = await Pipeline(dir, runs).RunAsync();
        Assert.True(again.Succeeded);
        Assert.Equal(5, again.Count(JobStatus.Skipped));
        Assert.All(runs.Values, n => Assert.Equal(1, n));
        Assert.All(again.Jobs, j => Assert.Equal("up to date", j.Reason));
    }

    [Fact]
    public async Task ChangedInput_RerunsOnlyItsDownstream_AndStopsWhenAnOutputIsUnchanged()
    {
        using var dir = Seed();
        var runs = new Dictionary<string, int>();
        await Pipeline(dir, runs).RunAsync();

        // New content, same length: upper reruns, count reruns (its input changed) but writes the same number,
        // so summary is skipped (early cutoff). Sample b is untouched.
        dir.Write("data/a.txt", "xyz");
        var plan = Pipeline(dir, runs).Plan();
        Assert.Equal("input changed: data/a.txt", plan["upper[s=a]"].Reason);
        Assert.Equal("upstream upper[s=a] will run", plan["count[s=a]"].Reason);
        Assert.False(plan["upper[s=b]"].WillRun);

        var r = await Pipeline(dir, runs).RunAsync();
        Assert.Equal(JobStatus.Succeeded, r["upper[s=a]"].Status);
        Assert.Equal(JobStatus.Succeeded, r["count[s=a]"].Status);
        Assert.Equal(JobStatus.Skipped, r["summary"].Status);
        Assert.Equal(JobStatus.Skipped, r["upper[s=b]"].Status);
        Assert.Equal(JobStatus.Skipped, r["count[s=b]"].Status);
        Assert.Equal("XYZ", dir.Read("work/a.up"));

        // A length change propagates all the way.
        dir.Write("data/a.txt", "abcdefg");
        var r2 = await Pipeline(dir, runs).RunAsync();
        Assert.Equal(JobStatus.Succeeded, r2["summary"].Status);
        Assert.Equal("input changed: results/a.n", r2["summary"].Reason);
        Assert.Equal("7,5", dir.Read("results/summary.txt"));
        Assert.Equal(1, runs["upper[s=b]"]);
    }

    [Fact]
    public async Task TouchingAFile_DoesNotRerun_InContentMode()
    {
        using var dir = Seed();
        var runs = new Dictionary<string, int>();
        await Pipeline(dir, runs).RunAsync();
        File.SetLastWriteTimeUtc(dir.File("data/a.txt"), DateTime.UtcNow.AddMinutes(5));
        var r = await Pipeline(dir, runs).RunAsync();
        Assert.Equal(5, r.Count(JobStatus.Skipped));

        var timestamp = Pipeline(dir, runs);
        timestamp.HashMode = HashMode.Timestamp;
        await timestamp.RunAsync();                                    // records timestamps
        File.SetLastWriteTimeUtc(dir.File("data/a.txt"), DateTime.UtcNow.AddMinutes(10));
        var plan = Pipeline(dir, runs);
        plan.HashMode = HashMode.Timestamp;
        Assert.True(plan.Plan()["upper[s=a]"].WillRun);
    }

    [Fact]
    public async Task ParameterOrVersionChange_InvalidatesThatRule()
    {
        using var dir = Seed();
        var runs = new Dictionary<string, int>();
        await Pipeline(dir, runs).RunAsync();

        var plan = Pipeline(dir, runs, suffixLength: 10).Plan();
        Assert.False(plan["upper[s=a]"].WillRun);
        Assert.StartsWith("rule definition changed", plan["count[s=a]"].Reason);
        var r = await Pipeline(dir, runs, suffixLength: 10).RunAsync();
        Assert.Equal(2, r.Jobs.Count(j => j.Rule == "count" && j.Status == JobStatus.Succeeded));
        Assert.Equal("13,15", dir.Read("results/summary.txt"));

        var v = await Pipeline(dir, runs, suffixLength: 10, version: "2").RunAsync();
        Assert.Equal(2, v.Jobs.Count(j => j.Rule == "upper" && j.Status == JobStatus.Succeeded));
        Assert.Equal(2, v.Jobs.Count(j => j.Rule == "count" && j.Status == JobStatus.Skipped)); // same upper output
    }

    [Fact]
    public async Task MissingOrEditedOutputs_AndForce_Rerun()
    {
        using var dir = Seed();
        var runs = new Dictionary<string, int>();
        await Pipeline(dir, runs).RunAsync();

        File.Delete(dir.File("work/b.up"));
        dir.Write("results/a.n", "999");
        var plan = Pipeline(dir, runs).Plan();
        Assert.Equal("missing output work/b.up", plan["upper[s=b]"].Reason);
        Assert.Equal("output modified outside the workflow: results/a.n", plan["count[s=a]"].Reason);
        Assert.Contains("upper[s=b]", plan.ToText());
        Assert.DoesNotContain("upper[s=a]", plan.ToText());
        Assert.Contains("skip  upper[s=a]", plan.ToText(includeUpToDate: true));
        await Pipeline(dir, runs).RunAsync();
        Assert.Equal("3", dir.Read("results/a.n"));

        var forced = await Pipeline(dir, runs).RunAsync(new RunOptions { ForceRules = ["upper"] });
        Assert.Equal(2, forced.Jobs.Count(j => j.Rule == "upper" && j.Status == JobStatus.Succeeded));
        Assert.Equal(3, forced.Count(JobStatus.Skipped));            // outputs unchanged, so downstream stays
        var all = await Pipeline(dir, runs).RunAsync(new RunOptions { ForceAll = true });
        Assert.Equal(5, all.Count(JobStatus.Succeeded));
    }

    [Fact]
    public async Task Values_ArePassedTyped_PersistedAndSkipped()
    {
        using var dir = Seed();
        int producerRuns = 0;
        Workflow Build()
        {
            var wf = new Workflow(dir.Path, "values");
            var stats = wf.Value<Stats>("stats", "stats/{s}");
            wf.Rule("measure").Input("raw", "data/{s}.txt").Output(stats).Run(ctx =>
            {
                Interlocked.Increment(ref producerRuns);
                ctx.Set(stats, new Stats(ctx.Wildcards["s"], File.ReadAllText(ctx.Input("raw")).Length));
            });
            wf.Rule("report").Input(stats, SampleSheet.FromColumn("s", "a", "b")).Output("o", "report.tsv").Run(ctx =>
                File.WriteAllLines(ctx.Output("o"), ctx.GetAll(stats).Select(x => $"{x.Sample}\t{x.Length}")));
            wf.Rule("one").Input(stats).Output("o", "one/{s}.txt").Run(ctx =>
                File.WriteAllText(ctx.Output("o"), ctx.Get(stats).Length.ToString()));
            return wf;
        }

        var r = await Build().RunAsync(new RunOptions { Targets = ["report.tsv", "one/b.txt"] });
        Assert.True(r.Succeeded, r.ToMarkdown());
        Assert.Equal(["a\t3", "b\t5"], File.ReadAllLines(dir.File("report.tsv")));
        Assert.Equal("5", dir.Read("one/b.txt"));
        Assert.Equal(2, producerRuns);
        Assert.Equal("value:stats/a", r["measure[s=a]"].Outputs[0]);

        var again = await Build().RunAsync(new RunOptions { Targets = ["report.tsv", "one/b.txt"] });
        Assert.Equal(4, again.Count(JobStatus.Skipped));
        Assert.Equal(2, producerRuns);

        var valueTarget = await Build().RunAsync(new RunOptions { Targets = ["value:stats/a"] });
        Assert.Equal(JobStatus.Skipped, Assert.Single(valueTarget.Jobs).Status);
    }

    private sealed record Stats(string Sample, int Length);
}
