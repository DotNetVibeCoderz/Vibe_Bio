using Biotech.Net.Core;
using Biotech.Net.Genomics.Calling;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Genomics.Vcf;

namespace Biotech.Net.Genomics.Tests;

public class MappingAndCallingTests
{
    [Fact]
    public void MapsSimulatedReadsBackToTheirOrigin()
    {
        var rnd = new Random(31);
        string chrA = Simulation.RandomDna(rnd, 12000), chrB = Simulation.RandomDna(rnd, 8000);
        var reference = new[] { Sequence.Dna(chrA, "chrA"), Sequence.Dna(chrB, "chrB") };
        var mapper = new ReadMapper(reference);
        var simA = Simulation.Reads(rnd, chrA, 300, 100, 0.01, prefix: "a");
        var simB = Simulation.Reads(rnd, chrB, 200, 100, 0.01, prefix: "b");
        var all = simA.Select(s => (s, "chrA")).Concat(simB.Select(s => (s, "chrB"))).ToList();
        var alns = mapper.Map(all.Select(x => x.s.Read)).ToList();
        Assert.Equal(all.Count, alns.Count);
        int correct = 0;
        for (int i = 0; i < all.Count; i++)
        {
            var (sim, chrom) = all[i];
            var a = alns[i];
            Assert.Equal(sim.Read.Id, a.ReadId);
            if (a.IsMapped && a.ReferenceName == chrom && a.Position == sim.Origin && a.IsReverse == sim.Reverse) correct++;
        }
        Assert.True(correct >= 0.99 * all.Count, $"{correct}/{all.Count} placed correctly");
        Assert.True(alns.Where(a => a.IsMapped).Average(a => a.MapQ) > 50);
        var a0 = alns.First(a => a.IsMapped && a.EditDistance > 0);
        Assert.Equal(a0.EditDistance, a0.Md.Count(char.IsLetter));
    }

    [Fact]
    public void AlignsReadsWithIndelsUsingBandedDp()
    {
        var rnd = new Random(32);
        string genome = Simulation.RandomDna(rnd, 5000);
        var mapper = new ReadMapper([Sequence.Dna(genome, "chr")]);

        // 3-base deletion in the read at offset 60.
        string del = genome.Substring(1000, 60) + genome.Substring(1063, 40);
        var d = mapper.Map(Sequence.Dna(del, "del"));
        Assert.True(d.IsMapped);
        Assert.Equal(1000, d.Position);
        Assert.Equal("60M3D40M", Cigar.Format(d.Cigar));
        Assert.Equal(3, d.EditDistance);
        Assert.Equal($"60^{genome.Substring(1060, 3)}40", d.Md);

        // 2-base insertion, read on the reverse strand.
        string ins = genome.Substring(2000, 50) + "TT" + genome.Substring(2050, 50);
        var r = mapper.Map(Sequence.Dna(Simulation.ReverseComplement(ins), "ins"));
        Assert.True(r.IsMapped);
        Assert.True(r.IsReverse);
        Assert.Equal(2000, r.Position);
        Assert.Equal(2, r.Cigar.Where(c => c.Op == 'I').Sum(c => c.Length));
        Assert.Equal(100, Cigar.ReferenceLength(r.Cigar));
        Assert.Equal(102, Cigar.QueryLength(r.Cigar));

        // Unrelated sequence stays unmapped.
        Assert.False(mapper.Map(Sequence.Dna(Simulation.RandomDna(new Random(99), 100), "junk")).IsMapped);
    }

    [Fact]
    public void RepeatsGetLowMappingQuality()
    {
        var rnd = new Random(33);
        string repeat = Simulation.RandomDna(rnd, 400);
        string genome = Simulation.RandomDna(rnd, 3000) + repeat + Simulation.RandomDna(rnd, 3000) + repeat + Simulation.RandomDna(rnd, 3000);
        var mapper = new ReadMapper([Sequence.Dna(genome, "chr")]);
        var inRepeat = mapper.Map(Sequence.Dna(repeat.Substring(100, 100), "rep"));
        Assert.True(inRepeat.IsMapped);
        Assert.Equal(0, inRepeat.MapQ);
        var unique = mapper.Map(Sequence.Dna(genome.Substring(500, 100), "uniq"));
        Assert.Equal(60, unique.MapQ);
    }

    [Fact]
    public void WritesAndReadsSam()
    {
        var rnd = new Random(34);
        string genome = Simulation.RandomDna(rnd, 3000);
        var reference = new[] { Sequence.Dna(genome, "chr1") };
        var mapper = new ReadMapper(reference);
        var reads = Simulation.Reads(rnd, genome, 20, 80, 0.0);
        var sw = new StringWriter();
        using (var w = new SamWriter(sw, reference, leaveOpen: true)) w.Write(mapper.Map(reads.Select(r => r.Read)));
        string sam = sw.ToString();
        Assert.StartsWith("@HD\tVN:1.6", sam);
        Assert.Contains("@SQ\tSN:chr1\tLN:3000", sam);
        var header = SamReader.ReadHeader(new StringReader(sam));
        Assert.Equal(("chr1", 3000L), header.Single());
        var recs = SamReader.Read(new StringReader(sam)).ToList();
        Assert.Equal(20, recs.Count);
        for (int i = 0; i < recs.Count; i++)
        {
            Assert.Equal(reads[i].Origin + 1, recs[i].Position); // SAM is 1-based
            Assert.Equal(reads[i].Reverse, recs[i].IsReverse);
            Assert.Equal("80M", recs[i].Cigar);
            Assert.Equal("0", recs[i].Tag("NM"));
            Assert.Equal(genome.Substring(reads[i].Origin, 80), recs[i].Sequence); // SEQ in reference orientation
            Assert.Equal(new string('?', 80), recs[i].Quality);                    // Q30 + 33
        }
        Assert.Equal(recs[0].ToString(), SamRecord.Parse(recs[0].ToString()).ToString());
    }

    [Fact]
    public void CallsPlantedHomozygousAndHeterozygousSnvs()
    {
        var rnd = new Random(35);
        string genome = Simulation.RandomDna(rnd, 10000);
        var reference = new[] { Sequence.Dna(genome, "chr1") };
        char Alt(char b) => b == 'A' ? 'G' : b == 'C' ? 'T' : b == 'G' ? 'A' : 'C';
        int homPos = 3000, hetPos = 6000;
        string hap1 = Simulation.Mutate(genome, (homPos, Alt(genome[homPos])), (hetPos, Alt(genome[hetPos])));
        string hap2 = Simulation.Mutate(genome, (homPos, Alt(genome[homPos])));
        // 2 × 1,500 reads of 100 bp ≈ 30× diploid coverage, 0.5% substitution errors at Q30.
        var reads = Simulation.Reads(rnd, hap1, 1500, 100, 0.005, prefix: "h1_")
            .Concat(Simulation.Reads(rnd, hap2, 1500, 100, 0.005, prefix: "h2_"))
            .Select(s => s.Read).ToList();

        var mapper = new ReadMapper(reference);
        var alignments = mapper.Map(reads).ToList();
        var pileup = new Pileup(reference).AddRange(alignments);
        var homColumn = pileup.Column("chr1", homPos);
        Assert.True(homColumn[Alt(genome[homPos])] >= 15);

        var calls = SnvCaller.Call(pileup).ToList();
        Assert.Equal(2, calls.Count);
        var hom = calls.Single(c => c.Position == homPos + 1);
        Assert.Equal(genome[homPos].ToString(), hom.Reference);
        Assert.Equal(Alt(genome[homPos]).ToString(), hom.Alternates.Single());
        Assert.True(hom.GetGenotype(0).IsHomozygousAlternate);
        Assert.True(hom.Quality > 100);
        Assert.Equal(0, hom.PhredLikelihoods(0)![2]);
        Assert.Equal(hom.Depth(0), hom.AlleleDepths(0)!.Sum());

        var het = calls.Single(c => c.Position == hetPos + 1);
        Assert.True(het.GetGenotype(0).IsHeterozygous);
        Assert.Equal(0, het.PhredLikelihoods(0)![1]);
        Assert.True(het.GenotypeQuality(0) > 20);
        var ad = het.AlleleDepths(0)!;
        Assert.InRange(ad[1] / (double)(ad[0] + ad[1]), 0.25, 0.75);

        // VCF round trip of the calls.
        var header = SnvCaller.CreateHeader(reference);
        string vcf = VcfWriter.ToString(header, calls);
        var (h2, back) = VcfReader.Parse(vcf);
        Assert.Equal(["sample"], h2.Samples);
        Assert.Equal(10000, h2.Contigs.Single().Length);
        Assert.Equal(calls.Select(c => c.ToString()), back.Select(c => c.ToString()));
    }

    [Fact]
    public void PileupHonoursBaseAndMappingQualityFilters()
    {
        var reference = new[] { Sequence.Dna("ACGTACGTAC", "r") };
        var p = new Pileup(reference, new PileupOptions { MinBaseQuality = 20, MinMapQ = 10 });
        CigarElement[] cig = [new(10, 'M')];
        p.Add("r", 0, cig, "ACGTACGTAC"u8, [30, 30, 30, 30, 30, 30, 30, 30, 30, 30], 60);
        p.Add("r", 0, cig, "ACGTTCGTAC"u8, [30, 30, 30, 30, 5, 30, 30, 30, 30, 30], 60); // low-quality mismatch
        p.Add("r", 0, cig, "ACGTTCGTAC"u8, [30, 30, 30, 30, 30, 30, 30, 30, 30, 30], 5); // low MAPQ read
        Assert.Equal(1, p.Depth("r", 4));
        Assert.Equal(2, p.Depth("r", 3));
        Assert.Empty(p.VariableColumns());
        var sam = SamRecord.Parse("q\t0\tr\t3\t60\t2M1D3M\t*\t0\t0\tGTCGT\t?????");
        p.Add(sam);
        Assert.Equal(3, p.Depth("r", 2));
        Assert.Equal(1, p.Depth("r", 4)); // deleted base not counted
    }
}
