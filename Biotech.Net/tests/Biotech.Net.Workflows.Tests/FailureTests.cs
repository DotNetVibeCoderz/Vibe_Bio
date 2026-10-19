namespace Biotech.Net.Workflows.Tests;

public class FailureTests
{
    private static Workflow Chain(TempDir dir, Func<bool> shouldFail, List<string>? ran = null)
    {
        dir.Write("in.txt", "x");
        var wf = new Workflow(dir.Path, "chain");
        void Note(string id)
        {
            if (ran is not null) lock (ran) ran.Add(id);
        }

        wf.Rule("first").Input("i", "in.txt").Output("o", "first.txt").Run(ctx =>
        {
            Note(ctx.JobId);
            File.WriteAllText(ctx.Output("o"), "1");
        });
        wf.Rule("flaky").Input("i", "first.txt").Output("o", "flaky.txt").Output("side", "flaky.side").Run(ctx =>
        {
            Note(ctx.JobId);
            File.WriteAllText(ctx.Output("side"), "complete");
            File.WriteAllText(ctx.Output("o"), "half-written");     // partial output, then the step dies
            if (shouldFail()) throw new InvalidOperationException("tool crashed");
            File.WriteAllText(ctx.Output("o"), "2");
        });
        wf.Rule("last").Input("i", "flaky.txt").Output("o", "last.txt").Run(ctx =>
        {
            Note(ctx.JobId);
            File.WriteAllText(ctx.Output("o"), "3");
        });
        wf.Target("last.txt");
        return wf;
    }

    [Fact]
    public async Task FailureMidRun_LeavesNoOutputs_AndResumeRunsOnlyWhatIsLeft()
    {
        using var dir = new TempDir();
        var ran = new List<string>();
        var r = await Chain(dir, () => true, ran).RunAsync();
        Assert.Equal(RunStatus.Failed, r.Status);
        Assert.Equal(JobStatus.Succeeded, r["first"].Status);
        Assert.Equal(JobStatus.Failed, r["flaky"].Status);
        Assert.Equal("InvalidOperationException: tool crashed", r["flaky"].Error);
        Assert.Equal(JobStatus.NotRun, r["last"].Status);
        Assert.False(dir.Exists("flaky.txt"));                      // atomic: the half-written output never appeared
        Assert.False(dir.Exists("flaky.side"));                     // nor the complete one: outputs commit together
        Assert.Empty(dir.TempFolders());
        Assert.Throws<WorkflowFailedException>(() => r.EnsureSuccess());

        var plan = Chain(dir, () => false).Plan();
        Assert.False(plan["first"].WillRun);
        Assert.StartsWith("incomplete", plan["flaky"].Reason);

        ran.Clear();
        var resumed = await Chain(dir, () => false, ran).RunAsync();
        Assert.True(resumed.Succeeded);
        Assert.Equal(["flaky", "last"], ran);
        Assert.Equal("2", dir.Read("flaky.txt"));
        Assert.Equal("3", dir.Read("last.txt"));
    }

    [Fact]
    public async Task KeepGoing_RunsIndependentJobsAfterAFailure()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("bad").Output("o", "bad.txt").Run(_ => throw new InvalidOperationException("nope"));
        wf.Rule("good").Output("o", "good.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "ok"));
        var stop = await wf.RunAsync(new RunOptions { Targets = ["bad.txt", "good.txt"], MaxParallelism = 1 });
        Assert.Equal(JobStatus.NotRun, stop["good"].Status);
        var go = await wf.RunAsync(new RunOptions { Targets = ["bad.txt", "good.txt"], MaxParallelism = 1, KeepGoing = true });
        Assert.Equal(JobStatus.Succeeded, go["good"].Status);
        Assert.Equal(RunStatus.Failed, go.Status);
    }

    [Fact]
    public async Task UndeclaredOrMissingOutputs_AndMissingSources_FailTheJob()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("lazy").Output("o", "lazy.txt").Run(_ => { });
        wf.Rule("needs").Input("i", "absent.txt").Output("o", "needs.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "x"));
        var r = await wf.RunAsync(new RunOptions { Targets = ["lazy.txt", "needs.txt"], KeepGoing = true });
        Assert.Contains("did not write output 'o'", r["lazy"].Error);
        Assert.Equal("missing input absent.txt (no rule produces it and it does not exist)", r["needs"].Error);
    }

    [Fact]
    public async Task Retries_BackOffExponentially_OnTheInjectedClock()
    {
        using var dir = new TempDir();
        var clock = new ManualTimeProvider();
        int calls = 0;
        var wf = new Workflow(dir.Path);
        wf.Rule("flaky").Output("o", "o.txt").Retries(3, TimeSpan.FromSeconds(10), backoffFactor: 3).Run(ctx =>
        {
            if (Interlocked.Increment(ref calls) < 3) throw new IOException($"transient {ctx.Attempt}");
            File.WriteAllText(ctx.Output("o"), ctx.Attempt.ToString());
        });
        var events = new EventLog();
        var r = await wf.RunAsync(new RunOptions { TimeProvider = clock }, events);
        Assert.True(r.Succeeded);
        Assert.Equal(3, r["flaky"].Attempts);
        Assert.Equal("3", dir.Read("o.txt"));
        Assert.Equal([TimeSpan.FromSeconds(10), TimeSpan.FromSeconds(30)], clock.Delays);
        var retries = events.Events.OfType<JobRetryingEvent>().ToList();
        Assert.Equal([1, 2], retries.Select(e => e.FailedAttempt));
        Assert.Equal("IOException: transient 1", retries[0].Error);
        Assert.Equal(TimeSpan.FromSeconds(40), r.FinishedAt - r.StartedAt);   // the fake clock advanced by the backoff only

        int always = 0;
        var wf2 = new Workflow(dir.Path);
        wf2.Rule("dead").Output("o", "dead.txt").Retries(2, TimeSpan.FromSeconds(1)).Run(_ =>
        {
            Interlocked.Increment(ref always);
            throw new IOException("still broken");
        });
        var r2 = await wf2.RunAsync(new RunOptions { TimeProvider = new ManualTimeProvider() });
        Assert.Equal(3, always);
        Assert.Equal(JobStatus.Failed, r2["dead"].Status);
        Assert.Equal(3, r2["dead"].Attempts);
    }

    [Fact]
    public async Task Timeout_FailsTheAttempt()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("slow").Output("o", "slow.txt").Timeout(TimeSpan.FromMilliseconds(100)).Run(async ctx =>
        {
            await File.WriteAllTextAsync(ctx.Output("o"), "partial");
            await Task.Delay(Timeout.Infinite, ctx.CancellationToken);
        });
        var r = await wf.RunAsync();
        Assert.Equal(JobStatus.Failed, r["slow"].Status);
        Assert.Equal("timed out after 0.1 s", r["slow"].Error);
        Assert.False(dir.Exists("slow.txt"));
        Assert.Empty(dir.TempFolders());
    }

    [Fact]
    public async Task Cancellation_StopsRunningJobs_AndTheNextRunResumes()
    {
        using var dir = new TempDir();
        var started = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        bool block = true;
        Workflow Build()
        {
            var wf = new Workflow(dir.Path);
            wf.Rule("long").Output("o", "long.txt").Run(async ctx =>
            {
                await File.WriteAllTextAsync(ctx.Output("o"), "partial");
                if (block)
                {
                    started.TrySetResult();
                    await Task.Delay(Timeout.Infinite, ctx.CancellationToken);
                }
            });
            wf.Rule("after").Input("i", "long.txt").Output("o", "after.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "a"));
            wf.Target("after.txt");
            return wf;
        }

        using var cts = new CancellationTokenSource();
        var run = Build().RunAsync(cancellationToken: cts.Token);
        await started.Task.WaitAsync(TimeSpan.FromSeconds(30));
        await cts.CancelAsync();
        var r = await run;
        Assert.Equal(RunStatus.Cancelled, r.Status);
        Assert.Equal(JobStatus.Cancelled, r["long"].Status);
        Assert.Equal(JobStatus.NotRun, r["after"].Status);
        Assert.False(dir.Exists("long.txt"));
        Assert.Empty(dir.TempFolders());

        block = false;
        Assert.StartsWith("incomplete", Build().Plan()["long"].Reason);
        Assert.True((await Build().RunAsync()).Succeeded);
        Assert.Equal("a", dir.Read("after.txt"));
    }

    [Fact]
    public async Task ConcurrentRuns_OfOneWorkingDirectory_AreRefused()
    {
        using var dir = new TempDir();
        var gate = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var started = new TaskCompletionSource(TaskCreationOptions.RunContinuationsAsynchronously);
        var wf = new Workflow(dir.Path);
        wf.Rule("hold").Output("o", "hold.txt").Run(async ctx =>
        {
            started.SetResult();
            await gate.Task;
            await File.WriteAllTextAsync(ctx.Output("o"), "x");
        });
        var first = wf.RunAsync();
        await started.Task.WaitAsync(TimeSpan.FromSeconds(30));
        await Assert.ThrowsAsync<InvalidOperationException>(() => wf.RunAsync());
        gate.SetResult();
        Assert.True((await first).Succeeded);
    }
}
