namespace Biotech.Net.Workflows.Tests;

public class PatternTests
{
    [Fact]
    public void Pattern_MatchesAndFormatsWildcards()
    {
        var p = new WildcardPattern("work/{sample}.{chunk,\\d+}.fastq");
        Assert.Equal(["sample", "chunk"], p.WildcardNames);
        var w = p.Match("work/A1.07.fastq")!;
        Assert.Equal("A1", w["sample"]);
        Assert.Equal("07", w["chunk"]);
        Assert.Null(p.Match("work/A1.x7.fastq"));       // constraint
        Assert.Null(p.Match("work/a/b.1.fastq"));        // default wildcard excludes '/'
        Assert.Equal("work/B.2.fastq", p.Format(Wildcards.Of(("sample", "B"), ("chunk", "2"))));
        Assert.Throws<KeyNotFoundException>(() => p.Format(Wildcards.Of(("sample", "B"))));
    }

    [Fact]
    public void Pattern_RepeatedWildcardMustAgree_AndBracesEscape()
    {
        var p = new WildcardPattern("{s}/{s}.txt");
        Assert.NotNull(p.Match("a/a.txt"));
        Assert.Null(p.Match("a/b.txt"));
        var lit = new WildcardPattern("x{{1}}/{s}");
        Assert.Equal("x{1}/q", lit.Format(Wildcards.Of(("s", "q"))));
        Assert.NotNull(lit.Match("x{1}/q"));
        Assert.Throws<FormatException>(() => new WildcardPattern("a/{s"));
        Assert.Throws<FormatException>(() => new WildcardPattern("a/{9x}"));
    }

    [Fact]
    public void Expand_ProductEachAndRows()
    {
        Assert.Equal(["c/A.1", "c/A.2", "c/B.1", "c/B.2"], Expand.Product("c/{s}.{n}", ("s", ["A", "B"]), ("n", ["1", "2"])));
        Assert.Equal(["v/A.vcf", "v/B.vcf"], Expand.Each("v/{sample}.vcf", "sample", ["A", "B"]));
        var sheet = SampleSheet.Parse("sample\tfastq\n# comment\nA\tdata/a.fq\nB\tdata/b.fq\n");
        Assert.Equal(["sample", "fastq"], sheet.Columns);
        Assert.Equal(["A", "B"], sheet.Column("sample"));
        Assert.Equal("data/b.fq", sheet.Find("sample", "B")["fastq"]);
        Assert.Equal(["out/A.txt", "out/B.txt"], Expand.Rows("out/{sample}.txt", sheet));
        Assert.Throws<FormatException>(() => SampleSheet.Parse("a,b\n1\n"));
        Assert.Equal(2, SampleSheet.Parse("sample,group\nA,x\nB,y").Count);
    }
}

public class GraphTests
{
    private static Workflow ScatterGather(TempDir dir, params string[] samples)
    {
        var wf = new Workflow(dir.Path, "sg");
        foreach (var s in samples) dir.Write($"data/{s}.txt", s);
        wf.Rule("upper").Input("raw", "data/{sample}.txt").Output("up", "work/{sample}.up")
            .Run(ctx => File.WriteAllText(ctx.Output("up"), File.ReadAllText(ctx.Input("raw")).ToUpperInvariant()));
        wf.Rule("merge").Input("parts", Expand.Each("work/{sample}.up", "sample", samples)).Output("all", "results/all.txt")
            .Run(ctx => File.WriteAllText(ctx.Output("all"), string.Join(",", ctx.Inputs("parts").Select(File.ReadAllText))));
        wf.Target("results/all.txt");
        return wf;
    }

    [Fact]
    public void Build_ScatterGather_TopologicalOrderAndEdges()
    {
        using var dir = new TempDir();
        var g = ScatterGather(dir, "a", "b", "c").Build();
        Assert.Equal(4, g.Jobs.Count);
        Assert.Equal(["upper[sample=a]", "upper[sample=b]", "upper[sample=c]", "merge"], g.Jobs.Select(j => j.Id));
        var merge = g["merge"];
        Assert.Equal(3, merge.Dependencies.Count);
        Assert.All(g.Jobs, j => Assert.All(j.Dependencies, d => Assert.True(d.Index < j.Index)));
        Assert.Equal(3, g.SourceArtifacts.Count);
        Assert.Equal("work/b.up", g["upper[sample=b]"].Outputs["up"].Id);
    }

    [Fact]
    public void Build_DetectsCyclesWithTheirPath()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("a").Input("x", "b.txt").Output("o", "a.txt").Run(_ => { });
        wf.Rule("b").Input("x", "a.txt").Output("o", "b.txt").Run(_ => { });
        var ex = Assert.Throws<WorkflowDefinitionException>(() => wf.Build("a.txt"));
        Assert.Contains("cycle", ex.Message);
        Assert.Contains("a -> b -> a", ex.Message);
    }

    [Fact]
    public void Build_RejectsAmbiguousProducersMissingInputsAndBodilessOutputs()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("one").Output("o", "out/{x}.txt").Run(_ => { });
        wf.Rule("two").Output("o", "out/{y}.txt").Run(_ => { });
        Assert.Contains("more than one rule", Assert.Throws<WorkflowDefinitionException>(() => wf.Build("out/q.txt")).Message);

        var wf2 = new Workflow(dir.Path);
        wf2.Rule("r").Input("i", "missing/{s}.txt").Output("o", "o/{s}.txt").Run(_ => { });
        Assert.Contains("not produced by any rule", Assert.Throws<WorkflowDefinitionException>(() => wf2.Build("nothing.txt")).Message);
        var g = wf2.Build("o/z.txt");                 // the missing source is reported at run time, per job
        Assert.Equal("missing/z.txt", Assert.Single(g.SourceArtifacts).Id);

        var wf3 = new Workflow(dir.Path);
        wf3.Rule("nobody").Output("o", "x.txt");
        Assert.Contains("no body", Assert.Throws<WorkflowDefinitionException>(() => wf3.Build("x.txt")).Message);

        var wf4 = new Workflow(dir.Path);
        wf4.Rule("bad").Input("i", "in/{other}.txt").Output("o", "out/{s}.txt").Run(_ => { });
        Assert.Contains("other", Assert.Throws<WorkflowDefinitionException>(() => wf4.Build("out/a.txt")).Message);
        Assert.Throws<ArgumentException>(() => new Workflow(dir.Path).Rule("m").Output("a", "{s}.a").Output("b", "{t}.b"));
    }

    [Fact]
    public void Export_DotAndMermaid()
    {
        using var dir = new TempDir();
        var wf = ScatterGather(dir, "a", "b");
        var g = wf.Build();
        string dot = g.ToDot();
        Assert.StartsWith("digraph \"sg\" {", dot);
        Assert.Contains("\"upper[sample=a]\" -> \"merge\";", dot);
        Assert.Contains("label=\"upper\\nsample=a\"", dot);
        Assert.EndsWith("}\n", dot);

        string mermaid = g.ToMermaid();
        Assert.StartsWith("flowchart LR\n", mermaid);
        Assert.Contains("n0[\"upper<br/>sample=a\"]", mermaid);
        Assert.Contains("n0 --> n2", mermaid);
        Assert.Contains("n1 --> n2", mermaid);

        string rules = wf.ToDot(rules: true);
        Assert.Contains("\"upper\" -> \"merge\";", rules);
        Assert.Contains("×2", rules);
        Assert.Contains("n0 --> n1", wf.ToMermaid(rules: true));
    }

    [Fact]
    public void Paths_AreNormalisedRelativeToTheWorkingDirectory()
    {
        using var dir = new TempDir();
        var wf = new Workflow(dir.Path);
        wf.Rule("w").Output("o", "./out\\{s}.txt").Run(ctx => File.WriteAllText(ctx.Output("o"), "x"));
        var g = wf.Build(Path.Combine(dir.Path, "out", "k.txt"));
        Assert.Equal("out/k.txt", Assert.Single(g.Jobs).Outputs["o"].Id);
    }
}
