namespace Biotech.Net.Workflows.Tests;

public class ProcessTests
{
    [Fact]
    public async Task ShellStep_WritesOutputs_RendersPlaceholders_AndLogs()
    {
        using var dir = new TempDir();
        dir.Write("in dir/greeting.txt", "hello");
        var wf = new Workflow(dir.Path, "shell");
        // Both shells understand "type/cat <in> > <out>"; the command differs only in the tool name.
        string copy = OperatingSystem.IsWindows() ? "type" : "cat";
        wf.Rule("copy").Input("src", "in dir/{name}.txt").Output("dst", "out/{name}.copy").Shell($"{copy} {{input.src}} > {{output.dst}}");
        var r = await wf.RunAsync(new RunOptions { Targets = ["out/greeting.copy"] });
        Assert.True(r.Succeeded, r.ToMarkdown());
        Assert.Equal("hello", dir.Read("out/greeting.copy").Trim());
        Assert.Equal(0, r["copy[name=greeting]"].ExitCode);
        string log = File.ReadAllText(dir.File(r["copy[name=greeting]"].Log!));
        Assert.Contains("\"in dir", log);                              // a path with a space is quoted
        Assert.Contains("exit code 0", log);

        // The rendered command is part of the definition: editing it reruns the rule.
        var wf2 = new Workflow(dir.Path, "shell");
        wf2.Rule("copy").Input("src", "in dir/{name}.txt").Output("dst", "out/{name}.copy").Shell($"{copy}  {{input.src}} > {{output.dst}}");
        Assert.StartsWith("rule definition changed", wf2.Plan(new RunOptions { Targets = ["out/greeting.copy"] }).Jobs[0].Reason);
    }

    [Fact]
    public async Task ShellStep_NonZeroExit_FailsWithCodeAndStderr()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("boom").Output("o", "never.txt").Shell("echo bad-thing 1>&2 && exit 3");
        var r = await wf.RunAsync();
        var job = r["boom"];
        Assert.Equal(JobStatus.Failed, job.Status);
        Assert.Equal(3, job.ExitCode);
        Assert.Contains("exited with code 3", job.Error);
        Assert.Contains("bad-thing", job.Error);
        Assert.False(dir.Exists("never.txt"));
    }

    [Fact]
    public async Task ProcessStep_AndInProcessToolCalls_CaptureOutput()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        var events = new EventLog();
        wf.Rule("version").Output("o", "dotnet-version.txt").Run(async ctx =>
        {
            int code = await ctx.RunProcessAsync("dotnet", ["--version"]);
            await File.WriteAllTextAsync(ctx.Output("o"), code.ToString());
        });
        wf.Rule("direct").Input("i", "dotnet-version.txt").Output("o", "direct.txt").Process("dotnet", "--version");
        // A process step that does not write its declared output fails the contract check.
        var r = await wf.RunAsync(new RunOptions { Targets = ["direct.txt"] }, events);
        Assert.Equal(JobStatus.Succeeded, r["version"].Status);
        Assert.Equal("0", dir.Read("dotnet-version.txt"));
        Assert.Contains(events.Events.OfType<JobLogEvent>(), e => e.JobId == "version" && e.Message.StartsWith("[out] 1"));
        Assert.Equal(JobStatus.Failed, r["direct"].Status);
        Assert.Equal(0, r["direct"].ExitCode);
        Assert.Contains("did not write output", r["direct"].Error);

        var missing = new Workflow(dir.Path);
        missing.Rule("nope").Output("o", "x.txt").Process("definitely-not-a-real-tool-xyz");
        var m = await missing.RunAsync();
        Assert.Contains("could not start", m["nope"].Error);
    }
}

public class EventAndReportTests
{
    private static Workflow Small(TempDir dir)
    {
        dir.Write("data/a.txt", "a");
        dir.Write("data/b.txt", "b");
        var wf = new Workflow(dir.Path, "small");
        wf.Rule("dup").Input("i", "data/{s}.txt").Output("o", "work/{s}.txt").Describe("doubles the text").Run(ctx =>
        {
            ctx.Log("doubling " + ctx.Wildcards["s"]);
            File.WriteAllText(ctx.Output("o"), File.ReadAllText(ctx.Input("i")) + File.ReadAllText(ctx.Input("i")));
        });
        wf.Rule("all").Input("parts", Expand.Each("work/{s}.txt", "s", ["a", "b"]));
        wf.Target("all");
        return wf;
    }

    [Fact]
    public async Task Events_AreRaisedInOrder_ThroughProgressAndStreaming()
    {
        using var dir = new TempDir();
        var events = new EventLog();
        var report = await Small(dir).RunAsync(new RunOptions { MaxParallelism = 1 }, events);
        Assert.True(report.Succeeded);
        var e = events.Events;
        Assert.IsType<RunStartedEvent>(e[0]);
        Assert.IsType<RunCompletedEvent>(e[^1]);
        Assert.Equal(3, ((RunStartedEvent)e[0]).JobCount);
        Assert.Equal(2, e.OfType<JobSucceededEvent>().Count());
        Assert.Contains(e, x => x is JobSkippedEvent { JobId: "all", Reason: "aggregate rule (no body)" });
        Assert.Contains(e, x => x is JobLogEvent { JobId: "dup[s=a]", Message: "doubling a" });
        Assert.Equal([1, 2, 3], e.OfType<RunProgressEvent>().Select(p => p.Finished));
        int started = e.ToList().FindIndex(x => x is JobStartedEvent { JobId: "dup[s=a]" });
        int succeeded = e.ToList().FindIndex(x => x is JobSucceededEvent { JobId: "dup[s=a]" });
        Assert.True(started >= 0 && started < succeeded);

        var streamed = new List<WorkflowEvent>();
        await foreach (var ev in Small(dir).RunStreamingAsync()) streamed.Add(ev);
        var last = Assert.IsType<RunCompletedEvent>(streamed[^1]);
        Assert.Equal(3, last.Report.Count(JobStatus.Skipped));
    }

    [Fact]
    public async Task Report_IsWrittenAsJsonMarkdownHtml_AndRoundTrips()
    {
        using var dir = new TempDir();
        var report = await Small(dir).RunAsync();
        string runDir = Path.Combine(dir.Path, Workflow.StateDirectoryName, "runs", report.RunId);
        foreach (var f in new[] { "report.json", "report.md", "report.html", "dag.dot" }) Assert.True(File.Exists(Path.Combine(runDir, f)), f);

        var back = WorkflowRunReport.FromJson(File.ReadAllText(Path.Combine(runDir, "report.json")));
        Assert.Equal(report.RunId, back.RunId);
        Assert.Equal(JobStatus.Succeeded, back["dup[s=a]"].Status);
        Assert.Equal("a", back["dup[s=a]"].Wildcards["s"]);
        Assert.Equal(report.RunId, Small(dir).LastReport()!.RunId);
        Assert.Contains("\"Status\": \"Succeeded\"", report.ToJson());

        string md = report.ToMarkdown();
        Assert.Contains("# Workflow run: small", md);
        Assert.Contains("| `dup[s=b]` | Succeeded | 1 |", md);
        Assert.Contains(Biotech.Net.Core.BiotechInfo.CreditsId, md);
        string html = report.ToHtml();
        Assert.Contains("<td class=\"Succeeded\">Succeeded</td>", html);
        Assert.Contains(Biotech.Net.Core.BiotechInfo.CreditsId, html);

        string log = File.ReadAllText(Path.Combine(dir.Path, report["dup[s=a]"].Log!));
        Assert.Contains("doubling a", log);
        Assert.Contains("succeeded", log);
    }

    [Fact]
    public void Plan_ExportsColouredGraphs()
    {
        using var dir = new TempDir();
        var plan = Small(dir).Plan();
        Assert.Equal(2, plan.ToRun.Count);
        Assert.Contains("fillcolor=\"#c9dcf2\"", plan.ToDot());
        Assert.Contains("class n0 willrun", plan.ToMermaid());
        Assert.Contains("    2 to run,     0 up to date  dup", plan.ToText());
    }
}
