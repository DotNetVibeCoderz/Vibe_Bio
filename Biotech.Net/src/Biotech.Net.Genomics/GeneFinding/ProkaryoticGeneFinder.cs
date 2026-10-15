using System.Globalization;
using Biotech.Net.Core;
using Biotech.Net.Sequences;

namespace Biotech.Net.Genomics.GeneFinding;

/// <summary>Options for <see cref="ProkaryoticGeneFinder"/>.</summary>
public sealed record GeneFinderOptions
{
    /// <summary>Genetic code (table 11 by default; its starts are ATG, GTG, TTG, ...).</summary>
    public GeneticCode Code { get; init; } = GeneticCode.Bacterial;

    /// <summary>Minimum gene length in bases including the stop codon (default 90).</summary>
    public int MinGeneLength { get; init; } = 90;

    /// <summary>Order of the Markov chains (default 5, GeneMark-like: a base conditioned on the previous five).</summary>
    public int MarkovOrder { get; init; } = 5;

    /// <summary>Pseudocount added to every Markov transition (default 0.5).</summary>
    public double Pseudocount { get; init; } = 0.5;

    /// <summary>Maximum overlap in bases allowed between two predicted genes (default 60).</summary>
    public int MaxOverlap { get; init; } = 60;

    /// <summary>Self-training iterations (default 3).</summary>
    public int TrainingIterations { get; init; } = 3;

    /// <summary>Minimum ORF length used to seed training (default 300).</summary>
    public int TrainingMinLength { get; init; } = 300;

    /// <summary>Minimum total log-odds score (natural log) for a gene to be reported (default 5).</summary>
    public double MinScore { get; init; } = 5;

    /// <summary>Shine–Dalgarno consensus searched upstream of starts (default AGGAGG).</summary>
    public string RbsConsensus { get; init; } = "AGGAGG";

    /// <summary>Weight of each matching RBS base (default 1.0 nat).</summary>
    public double RbsWeight { get; init; } = 1.0;

    /// <summary>Log prior of each start codon (ATG favoured over GTG and TTG).</summary>
    public IReadOnlyDictionary<string, double> StartCodonLogPrior { get; init; } = new Dictionary<string, double>
    {
        ["ATG"] = 0.0, ["GTG"] = -1.5, ["TTG"] = -2.0, ["CTG"] = -3.0, ["ATT"] = -3.0, ["ATC"] = -3.0, ["ATA"] = -3.0,
    };
}

/// <summary>A gene predicted by <see cref="ProkaryoticGeneFinder"/>.</summary>
/// <param name="SequenceId">Contig the gene is on.</param>
/// <param name="Start">0-based start (leftmost base, inclusive).</param>
/// <param name="End">0-based end (exclusive), including the stop codon.</param>
/// <param name="Strand">Coding strand.</param>
/// <param name="StartCodon">The initiation codon (e.g. ATG).</param>
/// <param name="Score">Total score (coding log-odds + RBS + start prior).</param>
/// <param name="CodingScore">Markov log-likelihood ratio, coding vs non-coding.</param>
/// <param name="RbsScore">Ribosome binding site score (0 when none found).</param>
/// <param name="RbsMotif">Matched RBS bases, or null.</param>
/// <param name="RbsSpacer">Bases between the RBS and the start codon (−1 when none).</param>
public sealed record PredictedGene(string SequenceId, int Start, int End, Strand Strand, string StartCodon,
    double Score, double CodingScore, double RbsScore, string? RbsMotif, int RbsSpacer)
{
    /// <summary>Length in bases (including the stop codon).</summary>
    public int Length => End - Start;

    /// <summary>The gene as a 0-based half-open genomic interval.</summary>
    public GenomicInterval Interval => new(SequenceId, Start, End, Strand);

    /// <summary>The coding sequence (reverse-complemented for minus-strand genes).</summary>
    public Sequence Extract(Sequence contig) => ToLocation().Extract(contig, SequenceId);

    /// <summary>The protein (translation of the CDS, first codon read as Met, stop removed).</summary>
    public Sequence Translate(Sequence contig, GeneticCode? code = null)
    {
        var p = Extract(contig).Translate(code ?? GeneticCode.Bacterial, toStop: true).ToString();
        if (p.Length > 0) p = "M" + p[1..];
        return Sequence.Protein(p, $"{SequenceId}_{Start + 1}_{End}");
    }

    private FeatureLocation ToLocation() => FeatureLocation.Span(Start, End, Strand);

    /// <summary>A <c>CDS</c> feature with score, start codon and RBS qualifiers.</summary>
    public SequenceFeature ToFeature(string? locusTag = null)
    {
        var f = new SequenceFeature("CDS", ToLocation())
            .With("score", Score.ToString("F2", CultureInfo.InvariantCulture))
            .With("start_codon", StartCodon)
            .With("inference", "ab initio prediction:Biotech.Net ProkaryoticGeneFinder");
        if (locusTag is not null) f.With("locus_tag", locusTag);
        if (RbsMotif is not null) f.With("rbs_motif", RbsMotif).With("rbs_spacer", RbsSpacer.ToString(CultureInfo.InvariantCulture));
        return f;
    }

    public override string ToString() => $"{SequenceId}:{Start + 1}..{End} ({(Strand == Strand.Reverse ? '-' : '+')}) {StartCodon} score {Score:F1}";
}

/// <summary>
/// Self-training ab initio gene finder for prokaryotic genomes, in the spirit of GeneMark/Glimmer: candidate ORFs
/// (stop to stop, every allowed start) are scored with a three-periodic (codon-position specific) Markov chain of
/// order <see cref="GeneFinderOptions.MarkovOrder"/> against a homogeneous non-coding chain; the start is chosen by
/// coding score plus a Shine–Dalgarno (RBS) score and a start-codon prior; overlapping candidates are resolved
/// greedily by score. Models are first trained on long ORFs and then re-estimated from the predictions.
/// </summary>
public sealed class ProkaryoticGeneFinder
{
    private double[][]? _coding;   // [codon position][ctx * 4 + base] log probabilities
    private double[]? _noncoding;  // [ctx * 4 + base]
    private readonly int _contexts;

    /// <summary>Creates a gene finder.</summary>
    public ProkaryoticGeneFinder(GeneFinderOptions? options = null)
    {
        Options = options ?? new GeneFinderOptions();
        if (Options.MarkovOrder is < 0 or > 8) throw new ArgumentOutOfRangeException(nameof(options), "Markov order must be 0–8.");
        _contexts = 1 << (2 * Options.MarkovOrder);
    }

    /// <summary>The options in use.</summary>
    public GeneFinderOptions Options { get; }

    /// <summary>Whether models have been trained.</summary>
    public bool IsTrained => _coding is not null;

    /// <summary>
    /// Trains the coding model on <paramref name="codingSequences"/> (in frame, start to stop) and the non-coding model on
    /// <paramref name="noncoding"/> (both strands are used).
    /// </summary>
    public void Train(IEnumerable<Sequence> codingSequences, IEnumerable<Sequence> noncoding)
    {
        var c = new double[3][];
        for (int f = 0; f < 3; f++) c[f] = new double[_contexts * 4];
        foreach (var g in codingSequences) CountInto(g.Span, c, periodic: true);
        var n = new double[_contexts * 4];
        var nn = new[] { n };
        foreach (var s in noncoding)
        {
            CountInto(s.Span, nn, periodic: false);
            CountInto(s.ReverseComplement().Span, nn, periodic: false);
        }
        for (int f = 0; f < 3; f++) ToLog(c[f]);
        ToLog(n);
        _coding = c;
        _noncoding = n;
    }

    private void CountInto(ReadOnlySpan<byte> s, double[][] counts, bool periodic)
    {
        int m = Options.MarkovOrder;
        int mask = _contexts - 1;
        int ctx = 0, valid = 0;
        for (int j = 0; j < s.Length; j++)
        {
            int b = Kmer.Code(s[j]);
            if (b < 0) { valid = 0; ctx = 0; continue; }
            if (valid >= m) counts[periodic ? j % 3 : 0][ctx * 4 + b]++;
            ctx = m == 0 ? 0 : ((ctx << 2) | b) & mask;
            valid++;
        }
    }

    private void ToLog(double[] t)
    {
        double pc = Options.Pseudocount;
        for (int ctx = 0; ctx < _contexts; ctx++)
        {
            double sum = 4 * pc;
            for (int b = 0; b < 4; b++) sum += t[ctx * 4 + b];
            for (int b = 0; b < 4; b++) t[ctx * 4 + b] = Math.Log((t[ctx * 4 + b] + pc) / sum);
        }
    }

    /// <summary>Self-trains on <paramref name="genome"/> and predicts its genes.</summary>
    public IReadOnlyList<PredictedGene> Predict(Sequence genome) => Predict([genome]);

    /// <summary>Self-trains on all <paramref name="contigs"/> together and predicts genes on each, sorted by position.</summary>
    public IReadOnlyList<PredictedGene> Predict(IReadOnlyList<Sequence> contigs)
    {
        // Seed: long non-overlapping ORFs as coding, whole contigs as background.
        var seeds = new List<Sequence>();
        foreach (var contig in contigs)
        {
            var orfs = OrfFinder.Find(contig, new OrfOptions
            {
                MinLength = Options.TrainingMinLength, Code = Options.Code, StartMode = OrfStartMode.AlternativeStarts,
            });
            var accepted = new List<Orf>();
            foreach (var o in orfs.OrderByDescending(o => o.Length))
                if (accepted.All(a => Overlap(a.Start, a.End, o.Start, o.End) <= Options.MaxOverlap)) accepted.Add(o);
            seeds.AddRange(accepted.Select(o => o.Extract(contig)));
        }
        Train(seeds, contigs);

        IReadOnlyList<PredictedGene> genes = PredictWithModels(contigs);
        for (int it = 1; it < Options.TrainingIterations && genes.Count > 0; it++)
        {
            var byId = contigs.ToDictionary(c => c.Id);
            var coding = genes.Select(g => g.Extract(byId[g.SequenceId])).ToList();
            var intergenic = new List<Sequence>();
            foreach (var contig in contigs)
            {
                int last = 0;
                foreach (var g in genes.Where(g => g.SequenceId == contig.Id).OrderBy(g => g.Start))
                {
                    if (g.Start - last > Options.MarkovOrder + 1) intergenic.Add(contig.Subsequence(last, g.Start - last));
                    last = Math.Max(last, g.End);
                }
                if (contig.Length - last > Options.MarkovOrder + 1) intergenic.Add(contig.Subsequence(last, contig.Length - last));
            }
            long intergenicBases = intergenic.Sum(s => (long)s.Length);
            Train(coding, intergenicBases >= 2000 ? intergenic : contigs);
            genes = PredictWithModels(contigs);
        }
        return genes;
    }

    /// <summary>Predicts genes with the current models (call <see cref="Train"/> first).</summary>
    public IReadOnlyList<PredictedGene> PredictWithModels(IReadOnlyList<Sequence> contigs)
    {
        if (_coding is null || _noncoding is null) throw new InvalidOperationException("Train the models first.");
        var all = new List<PredictedGene>();
        foreach (var contig in contigs)
        {
            var candidates = new List<PredictedGene>();
            Candidates(contig, contig.Span, Strand.Forward, candidates);
            var rc = contig.ReverseComplement();
            Candidates(contig, rc.Span, Strand.Reverse, candidates);
            var accepted = new List<PredictedGene>();
            foreach (var g in candidates.Where(g => g.Score >= Options.MinScore).OrderByDescending(g => g.Score))
                if (accepted.All(a => Overlap(a.Start, a.End, g.Start, g.End) <= Options.MaxOverlap)) accepted.Add(g);
            all.AddRange(accepted.OrderBy(g => g.Start));
        }
        return all;
    }

    private static int Overlap(int s1, int e1, int s2, int e2) => Math.Max(0, Math.Min(e1, e2) - Math.Max(s1, s2));

    private void Candidates(Sequence contig, ReadOnlySpan<byte> s, Strand strand, List<PredictedGene> output)
    {
        int n = s.Length, m = Options.MarkovOrder, mask = _contexts - 1;
        var code = Options.Code;
        var prefix = new double[n + 1];
        for (int phase = 0; phase < 3; phase++)
        {
            // Per-position log-odds assuming codons start at positions ≡ phase (mod 3); prefix sums.
            var cod = _coding!;
            var non = _noncoding!;
            int ctx = 0, valid = 0;
            prefix[0] = 0;
            for (int i = 0; i < n; i++)
            {
                int b = Kmer.Code(s[i]);
                double llr = 0;
                if (b < 0) { valid = 0; ctx = 0; }
                else
                {
                    if (valid >= m)
                    {
                        int cp = ((i - phase) % 3 + 3) % 3;
                        int k = ctx * 4 + b;
                        llr = cod[cp][k] - non[k];
                    }
                    ctx = m == 0 ? 0 : ((ctx << 2) | b) & mask;
                    valid++;
                }
                prefix[i + 1] = prefix[i] + llr;
            }

            // Walk codons of this frame: starts since the previous stop, closed by the next stop.
            var starts = new List<int>();
            bool sawStop = false;
            for (int p = phase; p + 3 <= n; p += 3)
            {
                int idx = GeneticCode.Index(s.Slice(p, 3));
                if (idx >= 0 && code.AminoAcidAt(idx) == '*')
                {
                    if (sawStop && starts.Count > 0) Choose(contig, s, strand, starts, p, prefix, output);
                    starts.Clear();
                    sawStop = true;
                    continue;
                }
                if (sawStop && idx >= 0 && code.IsStartIndex(idx) && Options.StartCodonLogPrior.ContainsKey(GeneticCode.CodonAt(idx)))
                    starts.Add(p);
            }
        }
    }

    private void Choose(Sequence contig, ReadOnlySpan<byte> s, Strand strand, List<int> starts, int stop, double[] prefix, List<PredictedGene> output)
    {
        int end = stop + 3;
        PredictedGene? best = null;
        foreach (int st in starts)
        {
            if (end - st < Options.MinGeneLength) break; // later starts are shorter
            string codon = GeneticCode.CodonAt(GeneticCode.Index(s.Slice(st, 3)));
            double prior = Options.StartCodonLogPrior[codon];
            double coding = prefix[stop] - prefix[st];
            var (rbs, motif, spacer) = ScoreRbs(s, st);
            double total = coding + rbs + prior;
            if (best is null || total > best.Score)
            {
                int gs = strand == Strand.Forward ? st : s.Length - end;
                int ge = strand == Strand.Forward ? end : s.Length - st;
                best = new PredictedGene(contig.Id, gs, ge, strand, codon, total, coding, rbs, motif, spacer);
            }
        }
        if (best is not null) output.Add(best);
    }

    /// <summary>Best Shine–Dalgarno match (longest common substring with the consensus, ≥ 3 bases) 4–15 bases upstream of the start.</summary>
    private (double Score, string? Motif, int Spacer) ScoreRbs(ReadOnlySpan<byte> s, int start)
    {
        string cons = Options.RbsConsensus;
        int bestLen = 0, bestSpacer = -1;
        string? bestMotif = null;
        int from = Math.Max(0, start - 21);
        for (int len = cons.Length; len >= 3 && bestLen == 0; len--)
        {
            for (int c0 = 0; c0 + len <= cons.Length; c0++)
            {
                var sub = cons.AsSpan(c0, len);
                for (int p = start - 4 - len; p >= from; p--)
                {
                    int spacer = start - (p + len);
                    if (spacer is < 4 or > 15) continue;
                    bool ok = true;
                    for (int k = 0; k < len && ok; k++) ok = (s[p + k] & 0xDF) == sub[k];
                    if (!ok) continue;
                    // Prefer spacers near the optimal 5–10 bases.
                    if (bestLen == 0 || Math.Abs(spacer - 7) < Math.Abs(bestSpacer - 7))
                    {
                        bestLen = len; bestSpacer = spacer; bestMotif = sub.ToString();
                    }
                }
            }
        }
        if (bestLen == 0) return (0, null, -1);
        double score = Options.RbsWeight * (bestLen - 2) - (bestSpacer is >= 5 and <= 10 ? 0 : 0.5);
        return (score, bestMotif, bestSpacer);
    }
}
