using Biotech.Net.Core;
using Biotech.Net.Proteomics.MassSpec;

namespace Biotech.Net.Proteomics.Tests;

public class SearchTests
{
    private const string Ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    private static readonly Protein[] Database =
    [
        Protein.FromString("P0CG48", Ubiquitin),
        Protein.FromString("SYN1", "MSAMPLERKGLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKK"),
    ];

    /// <summary>A clean b/y (1+) spectrum of <paramref name="sequence"/> with a few unrelated noise peaks.</summary>
    private static Spectrum Synthetic(string id, string sequence, int charge = 2, bool chargeKnown = true)
    {
        var pep = new Peptide(sequence);
        var mz = new List<double>();
        var it = new List<double>();
        foreach (var f in FragmentCalculator.Fragments(pep, [IonType.B, IonType.Y]))
        {
            mz.Add(f.Mz);
            it.Add(f.Type == IonType.Y ? 100 : 60);
        }
        foreach (double noise in new[] { 150.5, 333.3, 612.7, 1001.1 }) { mz.Add(noise); it.Add(10); }
        return new Spectrum(id, mz.ToArray(), it.ToArray())
        {
            PrecursorMz = pep.Mz(charge),
            PrecursorCharge = chargeKnown ? charge : 0,
            PrecursorIntensity = 1e6,
        };
    }

    [Fact]
    public void SyntheticSpectrumIdentifiesTargetPeptide()
    {
        var engine = new PeptideSpectrumMatcher(Database);
        var psm = engine.Match(Synthetic("s1", "TITLEVEPSDTIENVK"));
        Assert.NotNull(psm);
        Assert.Equal("TITLEVEPSDTIENVK", psm.Peptide.Sequence);
        Assert.False(psm.IsDecoy);
        Assert.Equal(["P0CG48"], psm.ProteinIds);
        Assert.Equal(2, psm.Charge);
        // Precursor charge 2 → fragment charge 1 only: all 15 b and 15 y ions are present.
        Assert.Equal(15, psm.MatchedNTerminalIons);
        Assert.Equal(15, psm.MatchedCTerminalIons);
        Assert.Equal(30, psm.TheoreticalIons);
        Assert.InRange(psm.PrecursorErrorPpm, -1e-3, 1e-3);
        // Hyperscore = ln(ΣI) + ln(15!) + ln(15!) with intensities scaled to a base peak of 100.
        double sumI = 15 * 100 + 15 * 60;
        Assert.Equal(Math.Log(sumI) + 2 * LogFactorial(15), psm.Score, 6);
        Assert.True(double.IsNaN(psm.QValue));
    }

    [Fact]
    public void ReversedDecoyWithSameCompositionIsDistinguishedByFragments()
    {
        // Reversed ubiquitin yields VNEITDSPEVELTITK: same composition (hence same precursor mass) as the target
        // TITLEVEPSDTIENVK, so only the fragment ladder can tell them apart.
        var engine = new PeptideSpectrumMatcher(Database);
        Assert.Equal(new Peptide("TITLEVEPSDTIENVK").MonoisotopicMass, new Peptide("VNEITDSPEVELTITK").MonoisotopicMass, 9);

        var decoy = engine.Match(Synthetic("d1", "VNEITDSPEVELTITK"));
        Assert.NotNull(decoy);
        Assert.Equal("VNEITDSPEVELTITK", decoy.Peptide.Sequence);
        Assert.True(decoy.IsDecoy);
        Assert.Equal(["DECOY_P0CG48"], decoy.ProteinIds);

        var target = engine.Match(Synthetic("t1", "TITLEVEPSDTIENVK"));
        Assert.NotNull(target);
        Assert.False(target.IsDecoy);
        Assert.Equal(decoy.Score, target.Score, 9); // mirror-image spectra score identically
    }

    [Fact]
    public void DecoyGenerationCanBeDisabledAndUsesPrefix()
    {
        var withDecoys = new PeptideSpectrumMatcher(Database);
        var noDecoys = new PeptideSpectrumMatcher(Database, new SearchParameters { GenerateDecoys = false });
        Assert.True(withDecoys.IndexedPeptides > noDecoys.IndexedPeptides);

        // Without decoys the decoy-derived spectrum can only hit the target of equal mass, with fewer matched ions.
        var psm = noDecoys.Match(Synthetic("d1", "VNEITDSPEVELTITK"));
        Assert.NotNull(psm);
        Assert.False(psm.IsDecoy);
        Assert.Equal("TITLEVEPSDTIENVK", psm.Peptide.Sequence);
        Assert.True(psm.MatchedNTerminalIons + psm.MatchedCTerminalIons < 30);

        var custom = new PeptideSpectrumMatcher(Database, new SearchParameters { DecoyPrefix = "REV_" });
        var d = custom.Match(Synthetic("d1", "VNEITDSPEVELTITK"));
        Assert.NotNull(d);
        Assert.Equal(["REV_P0CG48"], d.ProteinIds);
    }

    [Fact]
    public void SearchSkipsUnmatchableSpectraAndTriesDefaultCharges()
    {
        var engine = new PeptideSpectrumMatcher(Database);
        var empty = new Spectrum("none", [200.0, 300.0], [1.0, 1.0]) { PrecursorMz = 4321.0, PrecursorCharge = 1 };
        var noPrecursor = new Spectrum("ms1", [200.0], [1.0]) { MsLevel = 1 };
        Assert.Null(engine.Match(empty));
        Assert.Null(engine.Match(noPrecursor));

        var unknownCharge = Synthetic("u", "EGIPPDQQR", charge: 2, chargeKnown: false);
        var results = engine.Search([empty, Synthetic("a", "ESTLHLVLR"), noPrecursor, unknownCharge]);
        Assert.Equal(["a", "u"], results.Select(r => r.Spectrum.Id));
        Assert.Equal("ESTLHLVLR", results[0].Peptide.Sequence);
        Assert.Equal("EGIPPDQQR", results[1].Peptide.Sequence);
        Assert.Equal(2, results[1].Charge);
    }

    [Fact]
    public void FixedCarbamidomethylIsAppliedToCysteinePeptides()
    {
        var protein = Protein.FromString("CYS1", "MKAGDCLLPEKWTTSCMEAPKR");
        var engine = new PeptideSpectrumMatcher([protein], new SearchParameters
        {
            Digestion = new DigestionOptions(MissedCleavages: 0, MinLength: 4, MaxLength: 40),
            VariableModifications = [],
        });
        var mods = new Modification?[] { null, null, null, null, Modification.Carbamidomethyl, null, null, null, null };
        var pep = new Peptide("AGDCLLPEK", mods);
        var fragments = FragmentCalculator.Fragments(pep);
        var s = new Spectrum("c", fragments.Select(f => f.Mz).ToArray(), fragments.Select(_ => 50.0).ToArray())
        {
            PrecursorMz = pep.Mz(2), PrecursorCharge = 2,
        };
        var psm = engine.Match(s);
        Assert.NotNull(psm);
        Assert.Equal("AGDC[+57.021]LLPEK", psm.Peptide.ModifiedSequence);
    }

    private static PeptideSpectrumMatch Psm(double score, bool decoy, string id = "s")
        => new(new Spectrum(id, [], []), new Peptide("PEPTIDE"), [decoy ? "DECOY_P" : "P"], decoy, 2, score, 0, 0, 0, 0);

    [Fact]
    public void QValuesOnConstructedList()
    {
        // Scores (desc): T10 T9 D8 T7 T6 D5 D4 T3
        // FDR = D/T:      0  0 .5 1/3 .25 .5 .75 .6
        // q (running min from the bottom): 0 0 .25 .25 .25 .5 .6 .6
        var psms = new List<PeptideSpectrumMatch>
        {
            Psm(7, false), Psm(3, false), Psm(10, false), Psm(5, true), Psm(9, false), Psm(4, true), Psm(8, true), Psm(6, false),
        };
        TargetDecoy.AssignQValues(psms);
        var byScore = psms.OrderByDescending(p => p.Score).Select(p => p.QValue).ToArray();
        double[] expected = [0, 0, 0.25, 0.25, 0.25, 0.5, 0.6, 0.6];
        for (int i = 0; i < expected.Length; i++) Assert.Equal(expected[i], byScore[i], 12);

        var at1 = TargetDecoy.Filter(psms, 0.01);
        Assert.Equal([10.0, 9.0], at1.Select(p => p.Score));
        var at25 = TargetDecoy.Filter(psms, 0.25);
        Assert.Equal([10.0, 9.0, 7.0, 6.0], at25.Select(p => p.Score));
        Assert.All(at25, p => Assert.False(p.IsDecoy));
    }

    [Fact]
    public void TiedScoresShareQValueAndAllDecoysGiveOne()
    {
        var tied = new List<PeptideSpectrumMatch> { Psm(5, false), Psm(5, true), Psm(5, false), Psm(1, false) };
        TargetDecoy.AssignQValues(tied);
        // The three tied PSMs form one threshold: FDR = 1/2; below it, 1/3 → q = 1/3 for all of them.
        Assert.Equal(1.0 / 3, tied[0].QValue, 12);
        Assert.Equal(tied[0].QValue, tied[1].QValue);
        Assert.Equal(tied[0].QValue, tied[2].QValue);

        var decoys = new List<PeptideSpectrumMatch> { Psm(3, true), Psm(2, true) };
        TargetDecoy.AssignQValues(decoys);
        Assert.All(decoys, p => Assert.Equal(1.0, p.QValue));
        Assert.Empty(TargetDecoy.Filter(decoys, 1.0));
    }

    [Fact]
    public void QValuesAreMonotoneInScoreAndBoundedByFdr()
    {
        var rng = new Random(7);
        var psms = new List<PeptideSpectrumMatch>();
        for (int i = 0; i < 400; i++)
        {
            bool decoy = rng.NextDouble() < 0.3;
            // Targets are a mixture of correct (high-scoring) and random matches; decoys model random matches only.
            double score = decoy || rng.NextDouble() < 0.4 ? rng.NextDouble() * 20 : 15 + rng.NextDouble() * 20;
            psms.Add(Psm(Math.Round(score, 1), decoy));
        }
        TargetDecoy.AssignQValues(psms);
        var sorted = psms.OrderByDescending(p => p.Score).ToList();
        int t = 0, d = 0;
        for (int i = 0; i < sorted.Count; i++)
        {
            Assert.InRange(sorted[i].QValue, 0, 1);
            if (i > 0) Assert.True(sorted[i].QValue >= sorted[i - 1].QValue, "q-values must not decrease as the score falls");
            if (sorted[i].IsDecoy) d++; else t++;
            // At the end of each tie block, q never exceeds the FDR at that threshold.
            if (i == sorted.Count - 1 || sorted[i + 1].Score != sorted[i].Score)
                Assert.True(sorted[i].QValue <= (t == 0 ? 1 : Math.Min(1, (double)d / t)) + 1e-12);
        }
        // Filtering at 5% returns only targets, best first, all with q ≤ 0.05.
        var kept = TargetDecoy.Filter(psms, 0.05);
        Assert.NotEmpty(kept);
        Assert.All(kept, p => { Assert.False(p.IsDecoy); Assert.True(p.QValue <= 0.05); });
        Assert.True(kept.Zip(kept.Skip(1)).All(x => x.First.Score >= x.Second.Score));
    }

    [Fact]
    public void EndToEndSearchFdrAndInferenceFromMgf()
    {
        string[] targets = ["TITLEVEPSDTIENVK", "EGIPPDQQR", "TLSDYNIQK", "ESTLHLVLR", "LIFAGKQLEDGR"];
        var sw = new StringWriter();
        int n = 0;
        foreach (var seq in targets)
        {
            var s = Synthetic($"scan{++n}", seq);
            sw.WriteLine("BEGIN IONS");
            sw.WriteLine($"TITLE={s.Id}");
            sw.WriteLine(FormattableString.Invariant($"PEPMASS={s.PrecursorMz:R} 1000"));
            sw.WriteLine("CHARGE=2+");
            for (int i = 0; i < s.Count; i++) sw.WriteLine(FormattableString.Invariant($"{s.Mz[i]:R} {s.Intensity[i]:R}"));
            sw.WriteLine("END IONS");
        }
        var spectra = MgfReader.Read(new StringReader(sw.ToString())).ToList();
        Assert.Equal(5, spectra.Count);
        Assert.All(spectra, s => Assert.Equal(2, s.PrecursorCharge));

        var engine = new PeptideSpectrumMatcher(Database);
        var psms = engine.Search(spectra);
        Assert.Equal(targets, psms.Select(p => p.Peptide.Sequence));
        var accepted = TargetDecoy.Filter(psms, 0.01);
        Assert.Equal(5, accepted.Count);

        var groups = ProteinInference.Parsimony(accepted);
        var group = Assert.Single(groups);
        Assert.Equal("P0CG48", group.Representative);
        Assert.Equal(5, group.UniquePeptides.Count);

        var quant = LabelFreeQuantification.Quantify(accepted, Database);
        var q = Assert.Single(quant);
        Assert.Equal(5, q.SpectralCount);
        Assert.Equal(1.0, q.Nsaf, 12);
        Assert.Equal(5000, q.Intensity, 9);
        Assert.Equal(7, q.ObservablePeptides);
        Assert.Equal(5000.0 / 7, q.Ibaq, 9);
    }

    [Fact]
    public void ParsimonyMergesIndistinguishableAndDropsSubsumedProteins()
    {
        var map = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["PEPA"] = ["P1"],
            ["PEPB"] = ["P1", "P2"],
            ["PEPC"] = ["P3", "P4"],
            ["PEPD"] = ["P4", "P3"],
        };
        var groups = ProteinInference.Parsimony(map);
        Assert.Equal(2, groups.Count);
        Assert.Equal(["P1"], groups[0].ProteinIds);
        Assert.Equal(["PEPA", "PEPB"], groups[0].Peptides);
        Assert.Equal(["PEPA"], groups[0].UniquePeptides); // PEPB is shared with the subsumed P2
        Assert.Equal(["P3", "P4"], groups[1].ProteinIds);
        Assert.Equal(["PEPC", "PEPD"], groups[1].UniquePeptides);
        Assert.DoesNotContain(groups, g => g.ProteinIds.Contains("P2"));
    }

    [Fact]
    public void NsafNormalisesLengthCorrectedCounts()
    {
        var nsaf = LabelFreeQuantification.Nsaf([2.0, 4.0, 4.0], [100, 100, 200]);
        // SpC/L = 0.02, 0.04, 0.02 → sum 0.08.
        Assert.Equal(0.25, nsaf[0], 12);
        Assert.Equal(0.50, nsaf[1], 12);
        Assert.Equal(0.25, nsaf[2], 12);
        Assert.Throws<ArgumentException>(() => LabelFreeQuantification.Nsaf([1.0], [1, 2]));
    }

    [Fact]
    public void TopPeaksKeepsMostIntenseInMzOrder()
    {
        var s = new Spectrum("x", [500.0, 100.0, 300.0, 200.0], [5.0, 1.0, 9.0, 7.0]) { PrecursorMz = 600, PrecursorCharge = 2 };
        Assert.Equal([100.0, 200.0, 300.0, 500.0], s.Mz); // sorted on construction, intensities follow
        Assert.Equal([1.0, 7.0, 9.0, 5.0], s.Intensity);
        var top = s.TopPeaks(2);
        Assert.Equal([200.0, 300.0], top.Mz);
        Assert.Equal(600, top.PrecursorMz);
        Assert.Equal(2, s.FindPeak(300.004, 0.01));
        Assert.Equal(-1, s.FindPeak(400, 0.01));
        Assert.Equal(22, s.TotalIonCurrent);
    }

    private static double LogFactorial(int n)
    {
        double s = 0;
        for (int i = 2; i <= n; i++) s += Math.Log(i);
        return s;
    }
}
