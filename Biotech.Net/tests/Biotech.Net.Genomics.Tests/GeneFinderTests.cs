using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Genomics.Annotation;
using Biotech.Net.Genomics.GeneFinding;

namespace Biotech.Net.Genomics.Tests;

public class GeneFinderTests
{
    private sealed record TrueGene(int Start, int End, Strand Strand);

    /// <summary>
    /// A synthetic bacterial genome: genes built from a skewed codon-usage table, each preceded by an AGGAGG
    /// Shine–Dalgarno site 7 bases upstream of ATG, separated by AT-rich intergenic DNA; about half on the minus strand.
    /// </summary>
    private static (Sequence Genome, List<TrueGene> Genes) SyntheticGenome(int seed, int geneCount)
    {
        var rnd = new Random(seed);
        var code = GeneticCode.Bacterial;
        var sense = Enumerable.Range(0, 64).Where(i => code.AminoAcidAt(i) != '*').ToArray();
        var weights = sense.Select(_ => Math.Pow(rnd.NextDouble(), 3) + 0.01).ToArray();
        double total = weights.Sum();
        string Codon()
        {
            double x = rnd.NextDouble() * total;
            for (int i = 0; i < sense.Length; i++) { x -= weights[i]; if (x <= 0) return GeneticCode.CodonAt(sense[i]); }
            return GeneticCode.CodonAt(sense[^1]);
        }
        var sb = new StringBuilder();
        var genes = new List<TrueGene>();
        sb.Append(Simulation.RandomDna(rnd, 200, 0.3));
        for (int g = 0; g < geneCount; g++)
        {
            int codons = rnd.Next(120, 450);
            var gene = new StringBuilder("ATG");
            for (int c = 0; c < codons; c++) gene.Append(Codon());
            gene.Append(new[] { "TAA", "TAG", "TGA" }[rnd.Next(3)]);
            string unit = "AGGAGG" + Simulation.RandomDna(rnd, 7, 0.2) + gene;
            bool reverse = rnd.Next(2) == 1;
            int start = sb.Length;
            if (reverse)
            {
                sb.Append(Simulation.ReverseComplement(unit));
                genes.Add(new TrueGene(start, start + gene.Length, Strand.Reverse));
            }
            else
            {
                sb.Append(unit);
                genes.Add(new TrueGene(start + 13, start + unit.Length, Strand.Forward));
            }
            sb.Append(Simulation.RandomDna(rnd, rnd.Next(60, 250), 0.3));
        }
        return (Sequence.Dna(sb.ToString(), "synth"), genes);
    }

    [Fact]
    public void FindsGenesInASyntheticGenome()
    {
        var (genome, truth) = SyntheticGenome(21, 60);
        var finder = new ProkaryoticGeneFinder();
        var predicted = finder.Predict(genome);
        Assert.True(finder.IsTrained);

        static int StopEnd(int start, int end, Strand s) => s == Strand.Forward ? end : start;
        var predictedStops = predicted.Select(p => (StopEnd(p.Start, p.End, p.Strand), p.Strand)).ToHashSet();
        int found = truth.Count(t => predictedStops.Contains((StopEnd(t.Start, t.End, t.Strand), t.Strand)));
        int exact = truth.Count(t => predicted.Any(p => p.Start == t.Start && p.End == t.End && p.Strand == t.Strand));
        var truthStops = truth.Select(t => (StopEnd(t.Start, t.End, t.Strand), t.Strand)).ToHashSet();
        int falsePositives = predicted.Count(p => !truthStops.Contains((StopEnd(p.Start, p.End, p.Strand), p.Strand)));

        Assert.True(found >= 0.95 * truth.Count, $"sensitivity {found}/{truth.Count}");
        var wrong = truth.Where(t => !predicted.Any(p => p.Start == t.Start && p.End == t.End && p.Strand == t.Strand))
            .Select(t => (t, predicted.FirstOrDefault(p => p.Strand == t.Strand && StopEnd(p.Start, p.End, p.Strand) == StopEnd(t.Start, t.End, t.Strand))))
            .Select(x => $"{x.t.Start}-{x.t.End}{x.t.Strand}:{x.Item2}").Take(8);
        Assert.True(exact >= 0.8 * truth.Count, $"exact starts {exact}/{truth.Count} " + string.Join("; ", wrong));
        Assert.True(falsePositives <= 0.1 * predicted.Count, $"false positives {falsePositives}/{predicted.Count}");

        var withRbs = predicted.Where(p => p.RbsMotif is not null).ToList();
        Assert.True(withRbs.Count >= 0.8 * predicted.Count);
        Assert.Contains(withRbs, p => p.RbsMotif == "AGGAGG" && p.RbsSpacer == 7);

        // Proteins start with Met and contain no internal stops.
        var first = predicted[0];
        var protein = first.Translate(genome).ToString();
        Assert.StartsWith("M", protein);
        Assert.DoesNotContain('*', protein);
    }

    [Fact]
    public void WritesFeaturesAndGff3()
    {
        var (genome, _) = SyntheticGenome(22, 25);
        var genes = new ProkaryoticGeneFinder().Predict(genome);
        var feature = genes[0].ToFeature("SYN_0001");
        Assert.Equal("CDS", feature.Type);
        Assert.Equal("SYN_0001", feature["locus_tag"]);
        Assert.Equal(genes[0].Length, feature.Location.Length);

        string gff = GffWriter.ToString(genes);
        Assert.StartsWith("##gff-version 3\n", gff);
        var records = GffReader.Read(new StringReader(gff)).ToList();
        Assert.Equal(2 * genes.Count, records.Count);
        var cds = records.First(r => r.Type == "CDS");
        Assert.Equal(genes[0].Start, cds.Start);
        Assert.Equal(genes[0].End, cds.End);
        Assert.Equal("gene_1", cds.Parent);

        var sw = new StringWriter();
        GffWriter.Write(sw, genome.Id, genes.Select(g => g.ToFeature()));
        Assert.Equal(genes.Count, sw.ToString().Split('\n', StringSplitOptions.RemoveEmptyEntries).Length);
        var models = GffReader.BuildGenes(GffReader.Read(new StringReader(gff)));
        Assert.Equal(genes.Count, models.Count);
    }
}
