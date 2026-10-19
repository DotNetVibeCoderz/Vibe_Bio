using Biotech.Net.Genomics.Cram;
using Biotech.Net.Genomics.Mapping;
using static Biotech.Net.Genomics.Tests.HtsFixtures;

namespace Biotech.Net.Genomics.Tests;

public class CramTests
{
    private static string Conformance(string name) => Spec(Path.Combine("cram30", name));

    public static TheoryData<string> ConformanceFiles()
    {
        var data = new TheoryData<string>();
        foreach (var f in Directory.GetFiles(Path.Combine(AppContext.BaseDirectory, "data", "hts-specs", "cram30"), "*.cram").Order())
        {
            var name = Path.GetFileName(f);
            if (name.StartsWith("0903", StringComparison.Ordinal)) continue; // LZMA, see below
            data.Add(name);
        }
        return data;
    }

    private static void AssertSameLines(IReadOnlyList<string> expected, IReadOnlyList<string> actual, string what)
    {
        Assert.True(expected.Count == actual.Count, $"{what}: expected {expected.Count} records, got {actual.Count}");
        for (int i = 0; i < expected.Count; i++)
        {
            if (expected[i] == actual[i]) continue;
            var e = expected[i].Split('\t');
            var a = actual[i].Split('\t');
            for (int k = 0; k < Math.Max(e.Length, a.Length); k++)
                Assert.True(k < e.Length && k < a.Length && e[k] == a[k],
                    $"{what} record {i} ({e[0]}) column {k + 1}: expected '{(k < e.Length ? e[k] : "<none>")}' got '{(k < a.Length ? a[k] : "<none>")}'");
        }
    }

    [Theory]
    [MemberData(nameof(ConformanceFiles))]
    public void DecodesTheHtsSpecsConformanceFiles(string file)
    {
        string samPath = Conformance(Path.ChangeExtension(file, ".sam"));
        if (!File.Exists(samPath)) samPath += ".gz";
        var (_, expected) = ReadSam(samPath);
        using var reader = CramReader.Open(Conformance(file), Conformance("ce_trunc.fa"), new CramReaderOptions { RegenerateMdNm = false });
        var actual = reader.Records().Select(r => r.ToString()).ToList();
        AssertSameLines(expected, actual, file);
    }

    [Fact]
    public void LzmaBlocksAreReportedAsUnsupported()
    {
        using var reader = CramReader.Open(Conformance("0903_comp_lzma.cram"), Conformance("ce_trunc.fa"));
        var ex = Assert.Throws<NotSupportedException>(() => reader.Records().ToList());
        Assert.Contains("LZMA", ex.Message);
    }

    public static TheoryData<string, string, string> Profiles()
    {
        var data = new TheoryData<string, string, string>();
        foreach (var line in File.ReadAllLines(Hts("cram.manifest.tsv")))
        {
            var f = line.Split('\t');
            if (f[1].Contains("use_lzma")) continue;
            data.Add(f[0], f[1], f[2]);
        }
        return data;
    }

    [Theory]
    [MemberData(nameof(Profiles))]
    public void DecodesSamtoolsCramProfilesLikeSamtoolsView(string file, string options, string expectedFile)
    {
        var (header, expected) = ReadSam(Hts(expectedFile));
        using var reader = CramReader.Open(Hts(file), Hts("ref.fa"));
        Assert.Equal(string.Join('\n', header) + "\n", reader.Header.Text.EndsWith('\n') ? reader.Header.Text : reader.Header.Text + "\n");
        var actual = reader.Records().Select(r => r.ToString()).ToList();
        AssertSameLines(expected, actual, $"{file} ({options})");
    }

    [Fact]
    public void LzmaProfileIsReportedAsUnsupported()
    {
        using var reader = CramReader.Open(Hts("reads.30_lzma.cram"), Hts("ref.fa"));
        Assert.Throws<NotSupportedException>(() => reader.Records().ToList());
    }

    [Theory]
    [InlineData("reads.30.cram", "reads.30.cram.regions.tsv")]
    [InlineData("reads.30_multi.cram", "reads.30_multi.cram.regions.tsv")]
    public void RegionQueriesMatchSamtoolsView(string file, string regions)
    {
        using var reader = CramReader.Open(Hts(file), Hts("ref.fa"));
        Assert.True(reader.IsIndexed);
        foreach (var (region, expected) in ReadRegions(Hts(regions)))
        {
            var actual = reader.Query(region).Select(r => Columns(r.ToString())).ToList();
            Assert.True(expected.SequenceEqual(actual), $"{file} {region}: expected {expected.Count} records, got {actual.Count}");
        }
        // the same without the .crai (container-header pruning)
        reader.Index = null;
        foreach (var (region, expected) in ReadRegions(Hts(regions)))
            Assert.Equal(expected, reader.Query(region).Select(r => Columns(r.ToString())));
    }

    [Theory]
    [InlineData("reads.30.cram")]
    [InlineData("reads.30_multi.cram")]
    [InlineData("reads.31_archive.cram")]
    public void BuildsTheSameIndexAsSamtools(string file)
    {
        var expected = CramIndex.Load(Hts(file + ".crai")).Entries;
        using var reader = CramReader.Open(Hts(file), Hts("ref.fa"));
        Assert.Equal(expected, reader.BuildIndex().Entries);
    }

    [Theory]
    [InlineData("1400_index_simple.cram")]
    [InlineData("1402_index_3ref.cram")]
    [InlineData("1403_index_multiref.cram")]
    [InlineData("1404_index_multislice.cram")]
    public void ReadsConformanceIndexes(string file)
    {
        var crai = CramIndex.Load(Conformance(file + ".crai"));
        using var reader = CramReader.Open(Conformance(file), Conformance("ce_trunc.fa"), new CramReaderOptions { RegenerateMdNm = false });
        Assert.Equal(crai.Entries, reader.BuildIndex().Entries);
        var all = reader.Records().ToList();
        // every placed record is found by a query over its own reference
        foreach (var name in all.Where(r => r.ReferenceId >= 0).Select(r => r.ReferenceName).Distinct())
        {
            var expected = all.Where(r => r.ReferenceName == name).Select(r => r.ToString()).ToList();
            Assert.Equal(expected, reader.Query(name).Select(r => r.ToString()));
        }
    }

    [Fact]
    public void OpensCramByContent()
    {
        using var a = AlignmentFile.Open(Hts("reads.31.cram"), Hts("ref.fa"));
        Assert.IsType<CramReader>(a);
        Assert.Equal(3, ((CramReader)a).MajorVersion);
        Assert.Equal(1, ((CramReader)a).MinorVersion);
    }

    [Fact]
    public void MissingReferenceIsAClearError()
    {
        using var reader = CramReader.Open(Hts("reads.30.cram"));
        var ex = Assert.Throws<InvalidOperationException>(() => reader.Records().ToList());
        Assert.Contains("reference", ex.Message);
    }
}
