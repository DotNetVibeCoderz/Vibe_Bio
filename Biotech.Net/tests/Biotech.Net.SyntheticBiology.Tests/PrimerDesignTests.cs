using System.Diagnostics;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Tests;

public class PrimerDesignTests
{
    private static string RandomTemplate(int length, int seed)
    {
        var rng = new Random(seed);
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) sb.Append("ACGT"[rng.Next(4)]);
        return sb.ToString();
    }

    private static readonly string Template = RandomTemplate(1500, 7);

    // ---- OligoThermodynamics ---------------------------------------------------------------------------------

    [Fact]
    public void SelfDimer_OfSelfComplementaryOligo_IsTheFullDuplex()
    {
        const string s = "GCGCATATGCGC"; // its own reverse complement
        var d = OligoThermodynamics.SelfDimer(s);
        Assert.True(d.DeltaG37 < -10, $"{d.DeltaG37}");
        Assert.Equal(s.Length, d.Length);
        Assert.Equal(0, d.TopStart);
        Assert.Equal(0, d.BottomStart);
        // Stacks + both initiations + symmetry term: identical to the nearest-neighbour duplex.
        Assert.Equal(NearestNeighbourTm.Duplex(s).FreeEnergy37, d.DeltaG37, 9);
    }

    [Fact]
    public void SelfDimer_OfPolyA_IsZero()
    {
        var d = OligoThermodynamics.SelfDimer(new string('A', 20));
        Assert.Equal(0.0, d.DeltaG37);
        Assert.False(d.Found);
        Assert.Same(OligoDuplex.None, d);
    }

    [Fact]
    public void HeteroDimer_WithReverseComplement_EqualsPerfectDuplexFreeEnergy()
    {
        // The full-length register is one helix scored as stacks + one initiation per end — exactly the terms of
        // NearestNeighbourTm.Duplex, so the two agree to rounding (no terminal-term difference for a
        // non-self-complementary oligo; the symmetry term only applies to self-dimers).
        const string s = "AGCTTGCATGCCTGCAGGTC";
        var rc = DnaUtility.ReverseComplement(s);
        var d = OligoThermodynamics.HeteroDimer(s, rc);
        Assert.Equal(s.Length, d.Length);
        Assert.Equal(NearestNeighbourTm.Duplex(s).FreeEnergy37, d.DeltaG37, 9);
        Assert.True(d.DeltaG37 < -20);
    }

    [Fact]
    public void HeteroDimer_Diagram_HasThreeAlignedLines()
    {
        const string a = "AAAAGGATCCAAAA";
        const string b = "TTTTGGATCCTTTT";
        var d = OligoThermodynamics.HeteroDimer(a, b);
        Assert.True(d.Found);
        var lines = d.Diagram.Split('\n');
        Assert.Equal(3, lines.Length);
        Assert.StartsWith("5' ", lines[0]);
        Assert.StartsWith("3' ", lines[2]);
        Assert.Equal(d.Length, lines[1].Count(c => c == '|'));
        // Every '|' column holds a Watson-Crick pair.
        for (int col = 0; col < lines[1].Length; col++)
        {
            if (lines[1][col] != '|') continue;
            Assert.Equal(3, NearestNeighbourTm.CodeOf(lines[0][col]) + NearestNeighbourTm.CodeOf(lines[2][col]));
        }
    }

    [Fact]
    public void HeteroDimer_ReportsPositionsConsistently()
    {
        const string a = "AAAAGGATCCAAAA";
        const string b = "TTTTGGATCCTTTT";
        var d = OligoThermodynamics.HeteroDimer(a, b);
        for (int k = 0; k < d.Length; k++)
        {
            int i = d.TopStart + k, j = d.Offset - i;
            Assert.Equal(3, NearestNeighbourTm.CodeOf(a[i]) + NearestNeighbourTm.CodeOf(b[j]));
        }
        Assert.Equal(d.Offset - d.TopStart - d.Length + 1, d.BottomStart);
    }

    [Fact]
    public void ThreePrimeDimer_OnlyCountsHelicesAtA3End()
    {
        // GCGC is self-complementary. At the 3' end it makes a 3' self-dimer; at the 5' end it does not.
        const string at3 = "AAAAAAAAAAGCGC";
        const string at5 = "GCGCAAAAAAAAAA";

        var three = OligoThermodynamics.ThreePrimeDimer(at3);
        Assert.True(three.Found);
        Assert.Equal(at3.Length - 4, three.TopStart);
        Assert.Equal(4, three.Length);
        Assert.Equal(NearestNeighbourTm.Duplex("GCGC").FreeEnergy37, three.DeltaG37, 9);

        Assert.True(OligoThermodynamics.SelfDimer(at5).Found);
        Assert.False(OligoThermodynamics.ThreePrimeDimer(at5).Found);
    }

    [Fact]
    public void ThreePrimeDimer_IsNeverMoreStableThanTheUnanchoredDimer()
    {
        const string a = "ACGTTGCAGGCTAGCATCGA";
        const string b = "TCGATGCTAGCCTGCAACGT";
        Assert.True(OligoThermodynamics.ThreePrimeDimer(a, b).DeltaG37 >= OligoThermodynamics.HeteroDimer(a, b).DeltaG37);
    }

    [Fact]
    public void Hairpin_FoundInDesignedStemLoop()
    {
        // Stem CCGAGG / CCTCGG (reverse complements), loop TTTT.
        const string s = "CCGAGGTTTTCCTCGG";
        var h = OligoThermodynamics.Hairpin(s);
        Assert.True(h.Found);
        Assert.Equal(0, h.StemStart);
        Assert.Equal(6, h.StemLength);
        Assert.Equal(4, h.LoopLength);
        Assert.Equal("((((((....))))))", h.Structure);

        // Stacks of CCGAGG (CC=GG, CG, GA, AG=CT, GG) + 4-nt loop 3.5 kcal/mol.
        double stacks = NearestNeighbourTm.Duplex("CCGAGG").FreeEnergy37
                        - NearestNeighbourTm.InitiationDeltaG37(1) - NearestNeighbourTm.InitiationDeltaG37(2);
        Assert.Equal(stacks + 3.5, h.DeltaG37, 9);
        Assert.True(h.DeltaG37 < -3);
    }

    [Fact]
    public void Hairpin_GgggccccStemLoop_IsStable()
    {
        var h = OligoThermodynamics.Hairpin("GGGGCCCCTTTTGGGGCCCC");
        Assert.True(h.DeltaG37 < -4, $"{h.DeltaG37}");
        Assert.True(h.LoopLength >= 3);
    }

    [Theory]
    [InlineData("ACACACACACACACACACAC")] // no base can pair with any other
    [InlineData("TTTTTTTTTTCCCCCCCCCC")]
    public void Hairpin_AbsentInUnstructuredOligo(string s)
    {
        var h = OligoThermodynamics.Hairpin(s);
        Assert.False(h.Found);
        Assert.Equal(0.0, h.DeltaG37);
        Assert.Same(OligoHairpin.None, h);
    }

    [Fact]
    public void Hairpin_RespectsMinimumLoop()
    {
        var h = OligoThermodynamics.Hairpin("CCGAGGTTTTCCTCGG", minLoop: 6);
        Assert.True(!h.Found || h.LoopLength >= 6);
        Assert.Throws<ArgumentOutOfRangeException>(() => OligoThermodynamics.Hairpin("ACGT", minLoop: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => OligoThermodynamics.Hairpin("ACGT", minStem: 1));
    }

    [Fact]
    public void HairpinLoopPenalty_FollowsTableAndJacobsonStockmayer()
    {
        // SantaLucia & Hicks (2004) Annu. Rev. Biophys. Biomol. Struct. 33:415, Table 4 (hairpin column).
        Assert.Equal(3.5, OligoThermodynamics.HairpinLoopDeltaG37(4));
        Assert.Equal(6.3, OligoThermodynamics.HairpinLoopDeltaG37(30));
        // Extrapolation: dG(n) = dG(x) + 2.44 R T ln(n/x).
        double rt = 2.44 * 1.987 * 310.15 / 1000.0;
        Assert.Equal(6.3 + rt * Math.Log(2.0), OligoThermodynamics.HairpinLoopDeltaG37(60), 9);
        Assert.Equal(4.6 + rt * Math.Log(11.0 / 10.0), OligoThermodynamics.HairpinLoopDeltaG37(11), 9);
    }

    // ---- Evaluate ---------------------------------------------------------------------------------------------

    [Fact]
    public void Evaluate_ReportsTheOligoProperties()
    {
        const string s = "AGCTTGCATGCCTGCAGGTC";
        var p = PrimerDesigner.Evaluate(s);
        Assert.Equal(s, p.Sequence);
        Assert.Equal(-1, p.Start);
        Assert.Equal(Strand.Unknown, p.Strand);
        Assert.Equal(20, p.Length);
        Assert.Equal(NearestNeighbourTm.Tm(s, SaltConditions.StandardPcr), p.Tm, 9);
        Assert.Equal(0.6, p.GcFraction, 9);
        Assert.Equal(3, p.GcClamp);                   // last 5 nt: GGTC
        Assert.Equal(2, p.LongestHomopolymer);
        Assert.Equal(NearestNeighbourTm.Duplex("AGGTC").FreeEnergy37, p.EndStability, 9);
        Assert.Equal(OligoThermodynamics.SelfDimer(s).DeltaG37, p.SelfDimerDeltaG, 9);
        Assert.Equal(OligoThermodynamics.Hairpin(s).DeltaG37, p.HairpinDeltaG, 9);
        Assert.True(p.Penalty >= 0);

        var low = PrimerDesigner.Evaluate(s, new SaltConditions(0.010));
        Assert.True(low.Tm < p.Tm);
        Assert.Throws<ArgumentException>(() => PrimerDesigner.Evaluate("ACGNNT"));
    }

    // ---- Design -----------------------------------------------------------------------------------------------

    [Fact]
    public void DesignPairs_RespectsEveryConstraint()
    {
        var o = PrimerDesignOptions.Default;
        var pairs = PrimerDesigner.DesignPairs(Template, o);
        Assert.NotEmpty(pairs);
        Assert.True(pairs.Count <= o.MaxPairsReturned);

        foreach (var pair in pairs)
        {
            foreach (var p in new[] { pair.Forward, pair.Reverse })
            {
                Assert.InRange(p.Length, o.MinLength, o.MaxLength);
                Assert.InRange(p.Tm, o.MinTm, o.MaxTm);
                Assert.InRange(p.GcFraction, o.MinGc, o.MaxGc);
                Assert.InRange(p.GcClamp, o.MinGcClamp, o.MaxGcClamp);
                Assert.True(p.LongestHomopolymer <= o.MaxPolyX);
                Assert.True(p.SelfDimerDeltaG >= o.MaxSelfDimerDeltaG);
                Assert.True(p.ThreePrimeSelfDimerDeltaG >= o.MaxThreePrimeDimerDeltaG);
                Assert.True(p.HairpinDeltaG >= o.MaxHairpinDeltaG);
                Assert.True(p.EndStability >= o.MaxEndStabilityDeltaG);

                // Independently recomputed properties agree.
                Assert.Equal(NearestNeighbourTm.Tm(p.Sequence, o.SaltConditions), p.Tm, 9);
                Assert.Equal(DnaUtility.GcFraction(p.Sequence), p.GcFraction, 9);
                Assert.Equal(DnaUtility.LongestHomopolymer(DnaUtility.Ascii(p.Sequence)), p.LongestHomopolymer);
                Assert.Equal(p.Sequence[^5..].Count(c => c is 'G' or 'C'), p.GcClamp);
                Assert.Equal(OligoThermodynamics.SelfDimer(p.Sequence).DeltaG37, p.SelfDimerDeltaG, 9);
                Assert.Equal(OligoThermodynamics.Hairpin(p.Sequence).DeltaG37, p.HairpinDeltaG, 9);
            }

            Assert.Equal(Strand.Forward, pair.Forward.Strand);
            Assert.Equal(Strand.Reverse, pair.Reverse.Strand);
            Assert.InRange(pair.ProductSize, o.ProductSizeMin, o.ProductSizeMax);
            Assert.True(pair.TmDifference <= o.MaxTmDifference);
            Assert.Equal(Math.Abs(pair.Forward.Tm - pair.Reverse.Tm), pair.TmDifference, 9);
            Assert.True(pair.HeteroDimerDeltaG >= o.MaxHeteroDimerDeltaG);
            Assert.True(pair.ThreePrimeHeteroDimerDeltaG >= o.MaxThreePrimeDimerDeltaG);
            Assert.Equal(OligoThermodynamics.HeteroDimer(pair.Forward.Sequence, pair.Reverse.Sequence).DeltaG37, pair.HeteroDimerDeltaG, 9);
            Assert.True(pair.Reverse.Start >= pair.Forward.End);
            Assert.True(pair.Penalty >= pair.Forward.Penalty + pair.Reverse.Penalty);
        }
    }

    [Fact]
    public void DesignPairs_ProductAndPrimerSequencesMatchTheTemplate()
    {
        foreach (var pair in PrimerDesigner.DesignPairs(Template))
        {
            var f = pair.Forward; var r = pair.Reverse;
            Assert.Equal(Template.Substring(f.Start, f.Length), f.Sequence);
            Assert.Equal(DnaUtility.ReverseComplement(Template.Substring(r.Start, r.Length)), r.Sequence);
            Assert.Equal(r.End - f.Start, pair.ProductSize);
            Assert.Equal(Template[f.Start..r.End], pair.ProductSequence);
            Assert.StartsWith(f.Sequence, pair.ProductSequence);
            Assert.EndsWith(DnaUtility.ReverseComplement(r.Sequence), pair.ProductSequence);
        }
    }

    [Fact]
    public void DesignPairs_AreRankedByAscendingPenalty()
    {
        var pairs = PrimerDesigner.DesignPairs(Template, PrimerDesignOptions.Default with { MaxPairsReturned = 20 });
        Assert.True(pairs.Count > 1);
        for (int i = 1; i < pairs.Count; i++) Assert.True(pairs[i - 1].Penalty <= pairs[i].Penalty);
    }

    [Fact]
    public void DesignPairs_TopPairIsTheBestOverAllValidPairs()
    {
        // Exhaustive check of the branch-and-bound on a short template.
        string template = RandomTemplate(400, 11);
        var o = PrimerDesignOptions.Default with { MaxPairsReturned = 1, MaxPrimersReturned = int.MaxValue };
        var best = Assert.Single(PrimerDesigner.DesignPairs(template, o));
        var fs = PrimerDesigner.DesignForward(template, o);
        var rs = PrimerDesigner.DesignReverse(template, o);
        double bruteForce = double.PositiveInfinity;
        foreach (var f in fs)
        foreach (var r in rs)
        {
            int size = r.End - f.Start;
            if (r.Start < f.End || size < o.ProductSizeMin || size > o.ProductSizeMax) continue;
            double dTm = Math.Abs(f.Tm - r.Tm);
            if (dTm > o.MaxTmDifference) continue;
            double het = OligoThermodynamics.HeteroDimer(f.Sequence, r.Sequence).DeltaG37;
            double het3 = OligoThermodynamics.ThreePrimeDimer(f.Sequence, r.Sequence).DeltaG37;
            if (het < o.MaxHeteroDimerDeltaG || het3 < o.MaxThreePrimeDimerDeltaG) continue;
            double penalty = f.Penalty + r.Penalty + o.TmDifferenceWeight * dTm
                             + o.HeteroDimerWeight * -het + o.ThreePrimeDimerWeight * -het3;
            bruteForce = Math.Min(bruteForce, penalty);
        }
        Assert.Equal(bruteForce, best.Penalty, 9);
    }

    [Fact]
    public void DesignForwardAndReverse_ReturnRankedPrimersOnTheRightStrand()
    {
        var fs = PrimerDesigner.DesignForward(Template);
        var rs = PrimerDesigner.DesignReverse(Template);
        Assert.NotEmpty(fs);
        Assert.NotEmpty(rs);
        for (int i = 1; i < fs.Count; i++) Assert.True(fs[i - 1].Penalty <= fs[i].Penalty);
        for (int i = 1; i < rs.Count; i++) Assert.True(rs[i - 1].Penalty <= rs[i].Penalty);
        foreach (var f in fs) Assert.Equal(Template.Substring(f.Start, f.Length), f.Sequence);
        foreach (var r in rs) Assert.Equal(DnaUtility.ReverseComplement(Template.Substring(r.Start, r.Length)), r.Sequence);
    }

    [Fact]
    public void DesignPairs_TargetAndIncludedRegionAreRespected()
    {
        var target = new TemplateRegion(700, 50);
        var included = new TemplateRegion(400, 600);
        var o = PrimerDesignOptions.Default with { Target = target, IncludedRegion = included, MaxPairsReturned = 10 };
        var pairs = PrimerDesigner.DesignPairs(Template, o);
        Assert.NotEmpty(pairs);
        foreach (var pair in pairs)
        {
            Assert.True(pair.Forward.End <= target.Start);
            Assert.True(pair.Reverse.Start >= target.End);
            Assert.True(pair.Forward.Start >= included.Start);
            Assert.True(pair.Reverse.End <= included.End);
        }
    }

    [Fact]
    public void ImpossibleConstraints_ReturnEmpty()
    {
        Assert.Empty(PrimerDesigner.DesignPairs(Template, PrimerDesignOptions.Default with { MinTm = 90, MaxTm = 95, OptTm = 92 }));
        Assert.Empty(PrimerDesigner.DesignPairs(Template, PrimerDesignOptions.Default with { ProductSizeMin = 5000, ProductSizeMax = 6000 }));
        Assert.Empty(PrimerDesigner.DesignPairs("ACGTACGT"));
        Assert.Empty(PrimerDesigner.DesignPairs(string.Empty));
        Assert.Empty(PrimerDesigner.DesignPairs(new string('N', 1000)));
        Assert.Empty(PrimerDesigner.DesignPairs(new string('A', 1000)));
        Assert.Empty(PrimerDesigner.DesignForward(Template, PrimerDesignOptions.Default with { MinGc = 0.95 }));
        Assert.Empty(PrimerDesigner.DesignPairs(Template, PrimerDesignOptions.Default with { MinLength = 30, MaxLength = 25 }));
    }

    [Fact]
    public void DesignPairs_HandlesLowerCaseAndAmbiguousTemplate()
    {
        string template = Template[..600].ToLowerInvariant() + new string('N', 30) + Template[600..];
        var pairs = PrimerDesigner.DesignPairs(template);
        Assert.NotEmpty(pairs);
        foreach (var pair in pairs)
        {
            Assert.DoesNotContain('N', pair.Forward.Sequence);
            Assert.DoesNotContain('N', pair.Reverse.Sequence);
            Assert.Equal(template.ToUpperInvariant()[pair.Forward.Start..pair.Reverse.End], pair.ProductSequence);
        }
    }

    [Fact]
    public void DesignPairs_TwoKilobaseTemplate_IsFast()
    {
        string template = RandomTemplate(2000, 3);
        PrimerDesigner.DesignPairs(template); // warm-up (JIT, static tables)
        var sw = Stopwatch.StartNew();
        var pairs = PrimerDesigner.DesignPairs(template, PrimerDesignOptions.Default with { ProductSizeMin = 150, ProductSizeMax = 1500 });
        sw.Stop();
        Assert.NotEmpty(pairs);
        Assert.True(sw.Elapsed < TimeSpan.FromSeconds(1), $"took {sw.Elapsed.TotalMilliseconds:F0} ms");
    }
}
