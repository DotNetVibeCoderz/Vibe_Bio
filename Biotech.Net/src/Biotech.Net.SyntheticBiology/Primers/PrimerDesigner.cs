using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Primers;

/// <summary>
/// A Primer3-style PCR primer designer: enumerates candidate primers on both strands of a template, filters them
/// by Tm, composition and secondary structure, scores them with a weighted penalty and pairs them into
/// amplicons of the requested size.
/// </summary>
/// <remarks>
/// <para>
/// Tm is the SantaLucia (1998) nearest-neighbour Tm of <see cref="NearestNeighbourTm.Tm"/> under
/// <see cref="PrimerDesignOptions.SaltConditions"/>; dimers and hairpins come from
/// <see cref="OligoThermodynamics"/> (ungapped screening model, see its remarks). Constraints and penalty weights
/// are documented on <see cref="PrimerDesignOptions"/>.
/// </para>
/// <para>
/// <b>Algorithm.</b> Every window of <c>MinLength..MaxLength</c> nt on each strand is screened once, cheapest
/// tests first (ambiguous bases, GC fraction, GC clamp, poly-X, Tm, end stability, then hairpin and dimers).
/// The best <see cref="PrimerDesignOptions.MaxCandidatesPerStrand"/> survivors per strand are then paired by a
/// branch-and-bound search in ascending penalty order, so the comparatively expensive hetero-dimer check is only
/// run for pairs that could still enter the top <see cref="PrimerDesignOptions.MaxPairsReturned"/>. A 2 kb
/// template designs in a few tens of milliseconds.
/// </para>
/// <para>
/// Templates are upper-cased and U is read as T; windows containing any symbol other than A, C, G, T are
/// skipped. Impossible constraints simply yield an empty list.
/// </para>
/// </remarks>
public static class PrimerDesigner
{
    private const int EndWindow = 5;

    /// <summary>
    /// Measures a user-supplied oligo as a primer (no template): Tm, GC, dimers, hairpin, 3′ end properties and
    /// the penalty under the default weights. Nothing is rejected; compare the fields with your own limits.
    /// </summary>
    /// <param name="oligo">The primer, 5′→3′, at least 2 nt of A/C/G/T/U.</param>
    /// <param name="conditions">Salt and primer concentrations for Tm (default <see cref="SaltConditions.StandardPcr"/>).</param>
    /// <returns>A <see cref="Primer"/> with <c>Start = −1</c> and <see cref="Strand.Unknown"/>.</returns>
    /// <exception cref="ArgumentException">Shorter than 2 nt or containing a symbol other than A, C, G, T, U.</exception>
    public static Primer Evaluate(string oligo, SaltConditions? conditions = null)
        => Evaluate(oligo, conditions is null ? PrimerDesignOptions.Default : PrimerDesignOptions.Default with { SaltConditions = conditions });

    /// <summary>
    /// Measures a user-supplied oligo as a primer, using the salt conditions and penalty weights of
    /// <paramref name="options"/>. Nothing is rejected.
    /// </summary>
    /// <param name="oligo">The primer, 5′→3′, at least 2 nt of A/C/G/T/U.</param>
    /// <param name="options">Salt conditions and penalty weights.</param>
    /// <returns>A <see cref="Primer"/> with <c>Start = −1</c> and <see cref="Strand.Unknown"/>.</returns>
    /// <exception cref="ArgumentException">Shorter than 2 nt or containing a symbol other than A, C, G, T, U.</exception>
    public static Primer Evaluate(string oligo, PrimerDesignOptions options)
    {
        ArgumentNullException.ThrowIfNull(oligo);
        ArgumentNullException.ThrowIfNull(options);
        string seq = DnaUtility.NormalizeDna(oligo);
        if (seq.Length < 2) throw new ArgumentException("A primer needs at least 2 nucleotides.", nameof(oligo));
        foreach (char c in seq) NearestNeighbourTm.RequireCode(c);
        var codes = OligoThermodynamics.Encode(seq);
        Measure(seq, codes, 0, seq.Length, options, reject: false, out var m);
        return m.ToPrimer(seq, -1, Strand.Unknown);
    }

    /// <summary>
    /// Designs forward (left, plus-strand) primers on <paramref name="template"/>, best first.
    /// </summary>
    /// <param name="template">The template, plus strand 5′→3′.</param>
    /// <param name="options">Constraints and weights (default <see cref="PrimerDesignOptions.Default"/>).</param>
    /// <returns>Up to <see cref="PrimerDesignOptions.MaxPrimersReturned"/> primers in ascending penalty order.</returns>
    public static IReadOnlyList<Primer> DesignForward(string template, PrimerDesignOptions? options = null)
        => DesignSingle(template, options, Strand.Forward);

    /// <summary>
    /// Designs reverse (right, minus-strand) primers on <paramref name="template"/>, best first. Each primer's
    /// sequence is the reverse complement of the template at <c>[Start, End)</c>.
    /// </summary>
    /// <param name="template">The template, plus strand 5′→3′.</param>
    /// <param name="options">Constraints and weights (default <see cref="PrimerDesignOptions.Default"/>).</param>
    /// <returns>Up to <see cref="PrimerDesignOptions.MaxPrimersReturned"/> primers in ascending penalty order.</returns>
    public static IReadOnlyList<Primer> DesignReverse(string template, PrimerDesignOptions? options = null)
        => DesignSingle(template, options, Strand.Reverse);

    /// <summary>
    /// Designs PCR primer pairs on <paramref name="template"/>, best (lowest pair penalty) first.
    /// </summary>
    /// <remarks>
    /// Every returned pair satisfies all single-primer constraints, the product size range, the Tm difference cap,
    /// the hetero-dimer and 3′ hetero-dimer limits, the target and included regions, and has non-overlapping
    /// primers.
    /// </remarks>
    /// <param name="template">The template, plus strand 5′→3′.</param>
    /// <param name="options">Constraints and weights (default <see cref="PrimerDesignOptions.Default"/>).</param>
    /// <returns>Up to <see cref="PrimerDesignOptions.MaxPairsReturned"/> pairs; empty when nothing satisfies the constraints.</returns>
    public static IReadOnlyList<PrimerPair> DesignPairs(string template, PrimerDesignOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        options ??= PrimerDesignOptions.Default;
        Validate(options);
        if (options.MaxPairsReturned <= 0) return [];

        var ctx = new TemplateContext(template);
        var forward = Candidates(ctx, options, Strand.Forward, options.MaxCandidatesPerStrand);
        var reverse = Candidates(ctx, options, Strand.Reverse, options.MaxCandidatesPerStrand);
        if (forward.Count == 0 || reverse.Count == 0) return [];

        int k = options.MaxPairsReturned;
        var top = new List<(double Penalty, Candidate F, Candidate R, double TmDiff, double Hetero, double Hetero3)>(k + 1);
        double worst = double.PositiveInfinity;
        double bestReverse = reverse[0].Primer.Penalty;
        var target = options.Target;

        foreach (var f in forward)
        {
            if (f.Primer.Penalty + bestReverse >= worst) break;
            foreach (var r in reverse)
            {
                double basePenalty = f.Primer.Penalty + r.Primer.Penalty;
                if (basePenalty >= worst) break;

                var fp = f.Primer; var rp = r.Primer;
                if (rp.Start < fp.End) continue;
                int size = rp.End - fp.Start;
                if (size < options.ProductSizeMin || size > options.ProductSizeMax) continue;
                if (target is { } t && (fp.End > t.Start || rp.Start < t.End)) continue;
                double tmDiff = Math.Abs(fp.Tm - rp.Tm);
                if (tmDiff > options.MaxTmDifference) continue;

                double cheap = basePenalty + options.TmDifferenceWeight * tmDiff;
                if (cheap >= worst) continue;

                var hetero = OligoThermodynamics.BestDuplex(f.Codes, r.Codes, selfDimer: false, threePrimeAnchored: false).DeltaG;
                if (hetero < options.MaxHeteroDimerDeltaG) continue;
                var hetero3 = OligoThermodynamics.BestDuplex(f.Codes, r.Codes, selfDimer: false, threePrimeAnchored: true).DeltaG;
                if (hetero3 < options.MaxThreePrimeDimerDeltaG) continue;

                double penalty = cheap + options.HeteroDimerWeight * Math.Max(0, -hetero)
                                       + options.ThreePrimeDimerWeight * Math.Max(0, -hetero3);
                if (penalty >= worst) continue;

                int at = top.Count;
                while (at > 0 && top[at - 1].Penalty > penalty) at--;
                top.Insert(at, (penalty, f, r, tmDiff, hetero, hetero3));
                if (top.Count > k) top.RemoveAt(top.Count - 1);
                if (top.Count == k) worst = top[^1].Penalty;
            }
        }

        var result = new PrimerPair[top.Count];
        for (int i = 0; i < top.Count; i++)
        {
            var (penalty, f, r, tmDiff, hetero, hetero3) = top[i];
            int size = r.Primer.End - f.Primer.Start;
            result[i] = new PrimerPair(f.Primer, r.Primer, size, ctx.Plus.Substring(f.Primer.Start, size),
                tmDiff, hetero, hetero3, penalty);
        }
        return result;
    }

    // ------------------------------------------------------------------ internals

    private sealed class Candidate(Primer primer, byte[] codes)
    {
        public Primer Primer { get; } = primer;
        public byte[] Codes { get; } = codes;
    }

    private sealed class TemplateContext
    {
        public TemplateContext(string template)
        {
            Plus = DnaUtility.NormalizeDna(template);
            Minus = DnaUtility.ReverseComplement(Plus);
            PlusCodes = OligoThermodynamics.Encode(Plus);
            MinusCodes = OligoThermodynamics.Encode(Minus);
            AmbiguousPrefix = new int[Plus.Length + 1];
            for (int i = 0; i < Plus.Length; i++)
                AmbiguousPrefix[i + 1] = AmbiguousPrefix[i] + (PlusCodes[i] > 3 ? 1 : 0);
        }

        public string Plus { get; }
        public string Minus { get; }
        public byte[] PlusCodes { get; }
        public byte[] MinusCodes { get; }
        /// <summary>Prefix count of non-ACGT symbols on the plus strand.</summary>
        public int[] AmbiguousPrefix { get; }
        public int Length => Plus.Length;
    }

    private struct Metrics
    {
        public double Tm, Gc, SelfDimer, SelfDimer3, Hairpin, EndStability, Penalty;
        public int GcClamp, PolyX;

        public readonly Primer ToPrimer(string sequence, int start, Strand strand)
            => new(sequence, start, sequence.Length, strand, Tm, Gc, SelfDimer, SelfDimer3, Hairpin, EndStability,
                   GcClamp, PolyX, Penalty);
    }

    private static IReadOnlyList<Primer> DesignSingle(string template, PrimerDesignOptions? options, Strand strand)
    {
        ArgumentNullException.ThrowIfNull(template);
        options ??= PrimerDesignOptions.Default;
        Validate(options);
        if (options.MaxPrimersReturned <= 0) return [];
        var candidates = Candidates(new TemplateContext(template), options, strand, options.MaxPrimersReturned);
        var result = new Primer[candidates.Count];
        for (int i = 0; i < result.Length; i++) result[i] = candidates[i].Primer;
        return result;
    }

    private static void Validate(PrimerDesignOptions o)
    {
        ArgumentNullException.ThrowIfNull(o.SaltConditions, nameof(PrimerDesignOptions.SaltConditions));
        if (o.TmWeight < 0 || o.SizeWeight < 0 || o.GcWeight < 0 || o.SelfDimerWeight < 0 || o.ThreePrimeDimerWeight < 0
            || o.HairpinWeight < 0 || o.TmDifferenceWeight < 0 || o.HeteroDimerWeight < 0)
            throw new ArgumentException("Penalty weights must be non-negative.", nameof(o));
    }

    /// <summary>All windows on one strand passing the single-primer constraints, best first, at most <paramref name="keep"/>.</summary>
    private static List<Candidate> Candidates(TemplateContext ctx, PrimerDesignOptions o, Strand strand, int keep)
    {
        var list = new List<Candidate>();
        int n = ctx.Length;
        int minLen = Math.Max(2, o.MinLength), maxLen = Math.Min(o.MaxLength, n);
        if (minLen > maxLen || keep <= 0) return list;

        int regionStart = 0, regionEnd = n;
        if (o.IncludedRegion is { } inc)
        {
            regionStart = Math.Max(0, inc.Start);
            regionEnd = Math.Min(n, inc.End);
        }

        bool plus = strand == Strand.Forward;
        string strandSeq = plus ? ctx.Plus : ctx.Minus;
        byte[] strandCodes = plus ? ctx.PlusCodes : ctx.MinusCodes;

        for (int plusStart = regionStart; plusStart + minLen <= regionEnd; plusStart++)
        {
            for (int len = minLen; len <= maxLen; len++)
            {
                int plusEnd = plusStart + len;
                if (plusEnd > regionEnd) break;
                if (o.Target is { } t && (plus ? plusEnd > t.Start : plusStart < t.End)) continue;
                if (ctx.AmbiguousPrefix[plusEnd] != ctx.AmbiguousPrefix[plusStart]) break; // longer windows contain it too

                int s = plus ? plusStart : n - plusEnd; // window start on the strand being read
                if (!Measure(strandSeq, strandCodes, s, len, o, reject: true, out var m)) continue;
                string seq = strandSeq.Substring(s, len);
                list.Add(new Candidate(m.ToPrimer(seq, plusStart, strand), strandCodes.AsSpan(s, len).ToArray()));
            }
        }

        list.Sort(static (x, y) =>
        {
            int c = x.Primer.Penalty.CompareTo(y.Primer.Penalty);
            if (c != 0) return c;
            c = x.Primer.Start.CompareTo(y.Primer.Start);
            return c != 0 ? c : x.Primer.Length.CompareTo(y.Primer.Length);
        });
        if (list.Count > keep) list.RemoveRange(keep, list.Count - keep);
        return list;
    }

    /// <summary>
    /// Computes the primer metrics of <c>seq[start .. start+len)</c>, cheapest first. With <paramref name="reject"/>
    /// it stops (returning false) at the first violated constraint.
    /// </summary>
    private static bool Measure(string seq, byte[] codes, int start, int len, PrimerDesignOptions o, bool reject, out Metrics m)
    {
        m = default;
        var window = codes.AsSpan(start, len);

        int gc = 0;
        foreach (byte c in window) if (c is 1 or 2) gc++;
        m.Gc = (double)gc / len;
        if (reject && (m.Gc < o.MinGc || m.Gc > o.MaxGc)) return false;

        int endLen = Math.Min(EndWindow, len);
        for (int i = len - endLen; i < len; i++) if (window[i] is 1 or 2) m.GcClamp++;
        if (reject && (m.GcClamp < o.MinGcClamp || m.GcClamp > o.MaxGcClamp)) return false;

        int run = 1, longest = 1;
        for (int i = 1; i < len; i++)
        {
            run = window[i] == window[i - 1] ? run + 1 : 1;
            if (run > longest) longest = run;
        }
        m.PolyX = longest;
        if (reject && m.PolyX > o.MaxPolyX) return false;

        var chars = seq.AsSpan(start, len);
        var duplex = NearestNeighbourTm.Duplex(chars);
        m.Tm = NearestNeighbourTm.TmFromThermodynamics(duplex, o.SaltConditions, len, m.Gc, NearestNeighbourTm.IsSelfComplementary(chars));
        if (reject && (m.Tm < o.MinTm || m.Tm > o.MaxTm)) return false;

        m.EndStability = endLen >= 2 ? NearestNeighbourTm.Duplex(chars[(len - endLen)..]).FreeEnergy37 : 0.0;
        if (reject && m.EndStability < o.MaxEndStabilityDeltaG) return false;

        m.Hairpin = OligoThermodynamics.BestHairpin(window, 3, 2).DeltaG;
        if (reject && m.Hairpin < o.MaxHairpinDeltaG) return false;

        m.SelfDimer3 = OligoThermodynamics.BestDuplex(window, window, selfDimer: true, threePrimeAnchored: true).DeltaG;
        if (reject && m.SelfDimer3 < o.MaxThreePrimeDimerDeltaG) return false;

        m.SelfDimer = OligoThermodynamics.BestDuplex(window, window, selfDimer: true, threePrimeAnchored: false).DeltaG;
        if (reject && m.SelfDimer < o.MaxSelfDimerDeltaG) return false;

        m.Penalty = o.TmWeight * Math.Abs(m.Tm - o.OptTm)
                    + o.SizeWeight * Math.Abs(len - o.OptLength)
                    + o.GcWeight * 100.0 * Math.Abs(m.Gc - o.OptGc)
                    + o.SelfDimerWeight * Math.Max(0, -m.SelfDimer)
                    + o.ThreePrimeDimerWeight * Math.Max(0, -m.SelfDimer3)
                    + o.HairpinWeight * Math.Max(0, -m.Hairpin);
        return true;
    }
}
