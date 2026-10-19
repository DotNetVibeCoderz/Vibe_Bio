using Biotech.Net.Core;
using Biotech.Net.Genomics.Calling;
using Biotech.Net.Genomics.Cram;
using Biotech.Net.Genomics.IO;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Sequences;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

public class AlignmentInteropTests
{
    public static IEnumerable<(string Args, string Region, List<int> Depth)> DepthExpectations()
    {
        string? key = null;
        var cur = new List<int>();
        foreach (var line in File.ReadAllLines(Hts("reads.depth.tsv")))
        {
            if (line.StartsWith("## ", StringComparison.Ordinal))
            {
                if (key is not null) yield return (key.Split('|')[0], key.Split('|')[1], cur);
                key = line[3..];
                cur = [];
            }
            else if (line.Length > 0) cur.Add(int.Parse(line.Split('\t')[2]));
        }
        if (key is not null) yield return (key.Split('|')[0], key.Split('|')[1], cur);
    }

    [Theory]
    [InlineData("bam")]
    [InlineData("cram")]
    public void DepthMatchesSamtoolsDepth(string format)
    {
        using IAlignmentReader reader = format == "bam" ? BamReader.Open(Hts("reads.bam")) : CramReader.Open(Hts("reads.31.cram"), Hts("ref.fa"));
        int checkedSets = 0;
        foreach (var (args, region, expected) in DepthExpectations())
        {
            var options = args.Length == 0 ? new DepthOptions() : new DepthOptions { CountDeletions = true, MinMapQ = 30, MinBaseQuality = 20 };
            var r = HtsRegion.Parse(region);
            var depth = AlignmentDepth.Compute(reader, r.Name, r.Start, r.End, options);
            int bad = expected.Select((v, i) => v == depth[i] ? -1 : i).FirstOrDefault(i => i >= 0, -1);
            Assert.True(bad < 0, $"{args}|{region}: first difference at +{bad}: samtools {(bad < 0 ? 0 : expected[bad])}, ours {(bad < 0 ? 0 : depth[bad])}");
            checkedSets++;
        }
        Assert.Equal(8, checkedSets);
    }

    [Fact]
    public void MapperOutputGoesThroughSortedIndexedBamIntoTheCaller()
    {
        var reference = RandomSequence.Dna(3_000, seed: 21, id: "ref");
        var mutated = reference.ToString().ToCharArray();
        mutated[1500] = mutated[1500] == 'A' ? 'G' : 'A';
        var sample = new string(mutated);
        var rng = new Random(3);
        var reads = Enumerable.Range(0, 500).Select(i =>
        {
            int s = rng.Next(0, sample.Length - 100);
            return SequenceRead.FromAscii($"r{i}", sample.Substring(s, 100), new string('I', 100));
        }).ToList();
        var header = SamHeader.Create([reference], sortOrder: "coordinate");
        var records = CoordinateComparer.Sort(new ReadMapper([reference]).Map(reads).Select(a => a.ToAlignmentRecord(header)));

        var dir = Directory.CreateTempSubdirectory();
        try
        {
            string bam = Path.Combine(dir.FullName, "mapped.bam");
            using (var w = BamWriter.Create(bam, header, new BamWriterOptions { Index = BinningIndexFormat.Bai }))
                w.Write(records);
            using var reader = BamReader.Open(bam);
            Assert.True(reader.IsIndexed);
            var window = reader.Query("ref:1450-1550").ToList();
            Assert.NotEmpty(window);
            Assert.All(window, r => Assert.True(r.Overlaps(1449, 1550)));

            var pileup = new Pileup([reference]).AddRange(reader.Query("ref"));
            var calls = SnvCaller.Call(pileup);
            var call = Assert.Single(calls, v => v.Variant.Position == 1501);
            Assert.Equal(mutated[1500].ToString(), call.Variant.Alternates[0]);
        }
        finally { dir.Delete(true); }
    }
}
