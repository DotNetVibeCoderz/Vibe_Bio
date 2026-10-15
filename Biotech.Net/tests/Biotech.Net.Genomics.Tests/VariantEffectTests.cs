using Biotech.Net.Core;
using Biotech.Net.Genomics.Annotation;
using Biotech.Net.Genomics.Vcf;

namespace Biotech.Net.Genomics.Tests;

/// <summary>
/// A synthetic two-exon gene: exon1 [50,100), intron [100,150) (GT…AG), exon2 [150,200), CDS [60,191) = 27 codons.
/// chr1 carries it on the plus strand; chr2 is the reverse complement of chr1, so the same gene sits on the minus
/// strand and every variant must get identical transcript-level (HGVS) descriptions.
/// </summary>
public class VariantEffectTests
{
    private const string Cds = "ATGGCTAAACGGTTAGAATGGCCCAAGGTTCATTGCAGGGAGTGCGATTACCTGAAAGCCGTTGGCACCTTCATCAGGTAA";
    private const int L = 300;
    private readonly string _chr1;
    private readonly string _chr2;
    private readonly VariantEffectPredictor _vep;

    public VariantEffectTests()
    {
        Assert.Equal(81, Cds.Length);
        var rnd = new Random(7);
        var c = new char[L];
        for (int i = 0; i < L; i++) c[i] = "ACGT"[rnd.Next(4)];
        for (int i = 0; i < 40; i++) c[60 + i] = Cds[i];
        for (int i = 40; i < 81; i++) c[150 + i - 40] = Cds[i];
        c[100] = 'G'; c[101] = 'T'; c[148] = 'A'; c[149] = 'G';
        _chr1 = new string(c);
        _chr2 = Sequence.Dna(_chr1).ReverseComplement().ToString();

        var genome = new Genome("test");
        genome.Add(Chromosome.FromSequence(Sequence.Dna(_chr1, "chr1")));
        genome.Add(Chromosome.FromSequence(Sequence.Dna(_chr2, "chr2")));
        genome.Add(MakeGene("chr1", Strand.Forward, x => x));
        genome.Add(MakeGene("chr2", Strand.Reverse, x => x with { Start = L - x.End, End = L - x.Start, Strand = Strand.Reverse }));
        _vep = new VariantEffectPredictor(genome);
    }

    private static Gene MakeGene(string chrom, Strand strand, Func<GenomicInterval, GenomicInterval> map)
    {
        var gene = new Gene($"G_{chrom}", map(new GenomicInterval(chrom, 50, 200, Strand.Forward)), "TEST1");
        gene.Transcripts.Add(new Transcript($"T_{chrom}", gene.Id,
            [map(new GenomicInterval(chrom, 50, 100, Strand.Forward)), map(new GenomicInterval(chrom, 150, 200, Strand.Forward))])
        {
            CodingRegion = map(new GenomicInterval(chrom, 60, 191, Strand.Forward)),
        });
        return gene;
    }

    private static int CdsToGenomic(int cdsIndex) => cdsIndex < 40 ? 60 + cdsIndex : 150 + cdsIndex - 40;

    private static char Comp(char b) => b switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', _ => 'A' };

    /// <summary>Effects of an SNV at 0-based chr1 position <paramref name="g"/> on both strands' genes; asserts they agree.</summary>
    private VariantEffect Snv(int g, char alt)
    {
        var plus = _vep.Predict("chr1", g + 1, _chr1[g].ToString(), alt.ToString()).Single();
        int g2 = L - 1 - g;
        Assert.Equal(Comp(_chr1[g]), _chr2[g2]);
        var minus = _vep.Predict("chr2", g2 + 1, _chr2[g2].ToString(), Comp(alt).ToString()).Single();
        Assert.Equal(plus.Consequences, minus.Consequences);
        Assert.Equal(plus.HgvsC, minus.HgvsC);
        Assert.Equal(plus.HgvsP, minus.HgvsP);
        return plus;
    }

    [Fact]
    public void Missense()
    {
        var e = Snv(CdsToGenomic(9), 'T'); // CGG → TGG, codon 4
        Assert.Equal(Consequence.MissenseVariant, e.MostSevere);
        Assert.Equal(VariantImpact.Moderate, e.Impact);
        Assert.Equal("c.10C>T", e.HgvsC);
        Assert.Equal("p.Arg4Trp", e.HgvsP);
        Assert.Equal("R/W", e.AminoAcids);
        Assert.Equal("Cgg/Tgg", e.Codons);
        Assert.Equal(10, e.CdsPosition);
        Assert.Equal(4, e.ProteinPosition);
        Assert.Equal("TEST1", e.GeneName);
    }

    [Fact]
    public void Synonymous() // GCT → GCC (Ala2)
    {
        var e = Snv(CdsToGenomic(5), 'C');
        Assert.Equal([Consequence.SynonymousVariant], e.Consequences);
        Assert.Equal("c.6T>C", e.HgvsC);
        Assert.Equal("p.Ala2=", e.HgvsP);
    }

    [Fact]
    public void StopGained() // AAA → TAA (Lys3)
    {
        var e = Snv(CdsToGenomic(6), 'T');
        Assert.Equal(Consequence.StopGained, e.MostSevere);
        Assert.Equal(VariantImpact.High, e.Impact);
        Assert.Equal("c.7A>T", e.HgvsC);
        Assert.Equal("p.Lys3Ter", e.HgvsP);
    }

    [Fact]
    public void StartLost() // ATG → ACG
    {
        var e = Snv(CdsToGenomic(1), 'C');
        Assert.Equal(Consequence.StartLost, e.MostSevere);
        Assert.Equal("c.2T>C", e.HgvsC);
        Assert.Equal("p.Met1?", e.HgvsP);
    }

    [Fact]
    public void StopLostAndStopRetained()
    {
        var lost = Snv(CdsToGenomic(78), 'C'); // TAA → CAA (Gln)
        Assert.Equal(Consequence.StopLost, lost.MostSevere);
        Assert.Equal("c.79T>C", lost.HgvsC);
        Assert.StartsWith("p.Ter27GlnextTer", lost.HgvsP);

        var kept = Snv(CdsToGenomic(80), 'G'); // TAA → TAG
        Assert.Equal([Consequence.StopRetainedVariant], kept.Consequences);
        Assert.Equal("p.Ter27=", kept.HgvsP);
    }

    [Fact]
    public void SpliceSitesIntronAndUtrs()
    {
        var donor = Snv(100, 'A');
        Assert.Equal([Consequence.SpliceDonorVariant], donor.Consequences);
        Assert.Equal("c.40+1G>A", donor.HgvsC);

        var acceptor = Snv(149, 'A');
        Assert.Equal([Consequence.SpliceAcceptorVariant], acceptor.Consequences);
        Assert.Equal("c.41-1G>A", acceptor.HgvsC);

        var region = Snv(104, _chr1[104] == 'A' ? 'C' : 'A'); // +5
        Assert.Equal([Consequence.SpliceRegionVariant, Consequence.IntronVariant], region.Consequences);
        Assert.StartsWith("c.40+5", region.HgvsC);

        var exonRegion = Snv(98, _chr1[98] == 'A' ? 'C' : 'A'); // second-last base of exon 1
        Assert.Contains(Consequence.SpliceRegionVariant, exonRegion.Consequences);

        var intron = Snv(124, _chr1[124] == 'A' ? 'C' : 'A');
        Assert.Equal([Consequence.IntronVariant], intron.Consequences);
        Assert.StartsWith("c.40+25", intron.HgvsC);
        var intron2 = Snv(126, _chr1[126] == 'A' ? 'C' : 'A');
        Assert.StartsWith("c.41-24", intron2.HgvsC);

        var utr5 = Snv(55, _chr1[55] == 'A' ? 'C' : 'A');
        Assert.Equal([Consequence.FivePrimeUtrVariant], utr5.Consequences);
        Assert.StartsWith("c.-5", utr5.HgvsC);

        var utr3 = Snv(195, _chr1[195] == 'A' ? 'C' : 'A');
        Assert.Equal([Consequence.ThreePrimeUtrVariant], utr3.Consequences);
        Assert.StartsWith("c.*5", utr3.HgvsC);
    }

    [Fact]
    public void UpstreamDownstreamAndIntergenic()
    {
        var up = _vep.Predict("chr1", 20, _chr1[19].ToString(), "N").Single();
        Assert.Equal(Consequence.UpstreamGeneVariant, up.MostSevere);
        Assert.Equal(30, up.Distance);
        var down = _vep.Predict("chr1", 250, _chr1[249].ToString(), "N").Single();
        Assert.Equal(Consequence.DownstreamGeneVariant, down.MostSevere);
        // On the minus strand the 5' side is to the right.
        var upMinus = _vep.Predict("chr2", 280, _chr2[279].ToString(), "N").Single();
        Assert.Equal(Consequence.UpstreamGeneVariant, upMinus.MostSevere);

        var genome = new Genome("g");
        genome.Add(Chromosome.FromSequence(Sequence.Dna(_chr1, "chr1")));
        var gene = new Gene("G", new GenomicInterval("chr1", 50, 200, Strand.Forward));
        genome.Add(gene);
        var near = new VariantEffectPredictor(genome, new VariantEffectOptions { UpstreamDistance = 10, DownstreamDistance = 10 });
        var inter = near.Predict("chr1", 250, _chr1[249].ToString(), "N").Single();
        Assert.Equal([Consequence.IntergenicVariant], inter.Consequences);
        Assert.Null(inter.GeneId);
        // A gene without transcripts is treated as one non-coding exon.
        Assert.Equal(Consequence.NonCodingTranscriptExonVariant, near.Predict("chr1", 101, _chr1[100].ToString(), "C").Single().MostSevere);
        Assert.StartsWith("n.51", near.Predict("chr1", 101, _chr1[100].ToString(), "C").Single().HgvsC);
    }

    [Fact]
    public void FrameshiftDeletion()
    {
        // Delete one G of codon 7 (TGG); anchored VCF style at the previous base. HGVS shifts 3' → c.21del.
        int g = CdsToGenomic(18);
        string refAllele = _chr1.Substring(g, 2);
        var e = _vep.Predict("chr1", g + 1, refAllele, refAllele[..1]).Single();
        Assert.Equal(Consequence.FrameshiftVariant, e.MostSevere);
        Assert.Equal("c.21del", e.HgvsC);
        Assert.StartsWith("p.Trp7Cysfs", e.HgvsP);

        // Same deletion described on the minus-strand copy.
        int g2 = L - 1 - CdsToGenomic(20);
        string ref2 = _chr2.Substring(g2 - 1, 2);
        var e2 = _vep.Predict("chr2", g2, ref2, ref2[..1]).Single();
        Assert.Equal("c.21del", e2.HgvsC);
        Assert.Equal(e.HgvsP, e2.HgvsP);
    }

    [Fact]
    public void InframeDeletionAndDuplication()
    {
        // Delete codon 5 (TTA, CDS 13–15).
        int g = CdsToGenomic(11);
        string refAllele = _chr1.Substring(g, 4);
        var del = _vep.Predict("chr1", g + 1, refAllele, refAllele[..1]).Single();
        Assert.Equal(Consequence.InframeDeletion, del.MostSevere);
        Assert.Equal("c.13_15del", del.HgvsC);
        Assert.Equal("p.Leu5del", del.HgvsP);

        // Insert another TTA right after codon 5 → duplication.
        int gi = CdsToGenomic(14);
        var dup = _vep.Predict("chr1", gi + 1, _chr1[gi].ToString(), _chr1[gi] + "TTA").Single();
        Assert.Equal(Consequence.InframeInsertion, dup.MostSevere);
        Assert.Equal("c.13_15dup", dup.HgvsC);
        Assert.Equal("p.Leu5dup", dup.HgvsP);

        // Insert GGG (Gly) between codons 6 and 7.
        int gj = CdsToGenomic(17);
        var ins = _vep.Predict("chr1", gj + 1, _chr1[gj].ToString(), _chr1[gj] + "GGG").Single();
        Assert.Equal(Consequence.InframeInsertion, ins.MostSevere);
        Assert.Equal("c.18_19insGGG", ins.HgvsC);
        Assert.Equal("p.Glu6_Trp7insGly", ins.HgvsP);
    }

    [Fact]
    public void AnnotatesVcfRecordsWithCsq()
    {
        int g = CdsToGenomic(9);
        var v = new Variant("chr1", g + 1, _chr1[g].ToString(), ["T", "A"]);
        var rec = new VcfRecord(v);
        var effects = _vep.Annotate(rec);
        Assert.Equal(2, effects.Count);
        Assert.Contains("missense_variant", rec.Info["CSQ"]);
        Assert.Contains("p.Arg4Trp", rec.Info["CSQ"]);
        var header = new VcfHeader();
        VariantEffectPredictor.AddCsqHeader(header);
        Assert.NotNull(header.GetInfo("CSQ"));
        Assert.Equal("SO:0001583", Consequence.MissenseVariant.SoAccession());
        Assert.Equal(Consequence.StopGained, ConsequenceExtensions.ParseSoTerm("stop_gained"));
    }

    [Fact]
    public void BuildsGeneModelsFromGff3()
    {
        string gff = """
            ##gff-version 3
            chr1	test	gene	51	200	.	+	.	ID=G1;Name=TEST1
            chr1	test	mRNA	51	200	.	+	.	ID=T1;Parent=G1
            chr1	test	exon	51	100	.	+	.	Parent=T1
            chr1	test	exon	151	200	.	+	.	Parent=T1
            chr1	test	CDS	61	100	.	+	0	Parent=T1
            chr1	test	CDS	151	191	.	+	2	Parent=T1
            """;
        var genes = GffReader.BuildGenes(GffReader.Read(new StringReader(gff)));
        var gene = Assert.Single(genes);
        Assert.Equal("TEST1", gene.Name);
        var t = Assert.Single(gene.Transcripts);
        Assert.Equal(2, t.Exons.Count);
        Assert.Equal(new GenomicInterval("chr1", 60, 191, Strand.Forward), t.CodingRegion);

        var genome = new Genome("g");
        genome.Add(Chromosome.FromSequence(Sequence.Dna(_chr1, "chr1")));
        genome.Add(gene);
        int g = CdsToGenomic(9);
        Assert.Equal("p.Arg4Trp", new VariantEffectPredictor(genome).Predict("chr1", g + 1, "C", "T").Single().HgvsP);
    }

    [Fact]
    public void BuildsGeneModelsFromGtf()
    {
        string gtf = "chr1\ttest\ttranscript\t51\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
            + "chr1\ttest\texon\t51\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
            + "chr1\ttest\texon\t151\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
            + "chr1\ttest\tCDS\t61\t100\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n"
            + "chr1\ttest\tCDS\t151\t191\t.\t+\t2\tgene_id \"G1\"; transcript_id \"T1\";\n";
        var gene = Assert.Single(GffReader.BuildGenes(GffReader.Read(new StringReader(gtf))));
        Assert.Equal("G1", gene.Id);
        var t = Assert.Single(gene.Transcripts);
        Assert.Equal("T1", t.Id);
        Assert.Equal(2, t.Exons.Count);
        Assert.Equal(new GenomicInterval("chr1", 60, 191, Strand.Forward), t.CodingRegion);
    }
}
