namespace Biotech.Net.SyntheticBiology.Codons;

/// <summary>
/// A single place where a designed sequence breaks a manufacturing or cloning rule.
/// </summary>
/// <param name="Kind">The rule that was broken; one of the <c>Kind*</c> constants on <see cref="SequenceConstraints"/>.</param>
/// <param name="Start">Zero-based offset of the offending stretch on the given (sense) strand.</param>
/// <param name="Length">Length of the offending stretch in bases.</param>
/// <param name="Detail">Human-readable explanation, for instance the site that matched or the measured GC.</param>
public sealed record ConstraintViolation(string Kind, int Start, int Length, string Detail);

/// <summary>
/// The rules an optimised coding sequence has to satisfy before it can be ordered: no cloning sites, no long
/// homopolymers, GC inside the window a synthesiser can handle, no strong hairpins and no long direct repeats.
/// </summary>
/// <remarks>
/// <para>
/// These are the constraints commercial gene synthesis quotes are scored against. Homopolymer runs and extreme
/// local GC cause polymerase slippage and assembly failure; perfect inverted repeats fold into hairpins that
/// stall both synthesis and reverse transcription; long direct repeats cause mis-assembly and recombination;
/// and the forbidden sites are simply the restriction sites the downstream cloning scheme needs to stay unique.
/// </para>
/// <para>
/// All checks are made on the sense strand and, where a site is not palindromic, on its reverse complement too,
/// because a restriction enzyme cuts double-stranded DNA and does not care which strand carries the recognition
/// sequence. Create an instance with an object initialiser and keep it: it is immutable and thread safe.
/// </para>
/// </remarks>
public sealed class SequenceConstraints
{
    /// <summary><see cref="ConstraintViolation.Kind"/> for a forbidden (usually restriction) site.</summary>
    public const string KindForbiddenSite = "ForbiddenSite";

    /// <summary><see cref="ConstraintViolation.Kind"/> for a run of one base longer than <see cref="MaxHomopolymer"/>.</summary>
    public const string KindHomopolymer = "Homopolymer";

    /// <summary><see cref="ConstraintViolation.Kind"/> for a window whose GC is outside <see cref="MinGc"/>–<see cref="MaxGc"/>.</summary>
    public const string KindGcWindow = "GcWindow";

    /// <summary><see cref="ConstraintViolation.Kind"/> for an inverted repeat able to fold into a hairpin.</summary>
    public const string KindHairpin = "Hairpin";

    /// <summary><see cref="ConstraintViolation.Kind"/> for a direct repeat longer than <see cref="MaxDirectRepeat"/>.</summary>
    public const string KindDirectRepeat = "DirectRepeat";

    /// <summary>
    /// Recognition sequences of the enzymes most often kept out of a synthetic gene: EcoRI (GAATTC),
    /// BamHI (GGATCC), HindIII (AAGCTT), XhoI (CTCGAG), XbaI (TCTAGA), SpeI (ACTAGT), PstI (CTGCAG),
    /// NdeI (CATATG), NheI (GCTAGC), SalI (GTCGAC), KpnI (GGTACC), SacI (GAGCTC), NotI (GCGGCCGC),
    /// and the Golden Gate type IIS enzymes BsaI (GGTCTC), BsmBI/Esp3I (CGTCTC) and BbsI (GAAGAC).
    /// The type IIS sites are not palindromic, so <see cref="Check(string)"/> looks for their reverse
    /// complements as well.
    /// </summary>
    public static IReadOnlyList<string> CommonRestrictionSites { get; } =
    [
        "GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "TCTAGA", "ACTAGT", "CTGCAG", "CATATG",
        "GCTAGC", "GTCGAC", "GGTACC", "GAGCTC", "GCGGCCGC", "GGTCTC", "CGTCTC", "GAAGAC",
    ];

    /// <summary>
    /// Sequences that must not occur, written as IUPAC patterns (so <c>GGNCC</c> and the like are allowed).
    /// Both strands are searched. Empty by default.
    /// </summary>
    public IReadOnlyList<string> ForbiddenSites { get; init; } = [];

    /// <summary>Longest run of a single base that is still acceptable; runs longer than this are reported.</summary>
    public int MaxHomopolymer { get; init; } = 6;

    /// <summary>Lowest acceptable GC fraction in any window (0–1).</summary>
    public double MinGc { get; init; } = 0.30;

    /// <summary>Highest acceptable GC fraction in any window (0–1).</summary>
    public double MaxGc { get; init; } = 0.70;

    /// <summary>
    /// Width in bases of the sliding window the GC limits are applied to. A sequence shorter than the window is
    /// evaluated as a single window, so the limits then apply to the whole sequence.
    /// </summary>
    public int GcWindow { get; init; } = 50;

    /// <summary>Whether to look for inverted repeats that can fold back into a hairpin.</summary>
    public bool AvoidHairpins { get; init; } = true;

    /// <summary>Shortest perfectly paired stem that counts as a hairpin.</summary>
    public int MinHairpinStem { get; init; } = 10;

    /// <summary>Shortest loop between the two arms of a stem; closer arms cannot pair sterically.</summary>
    public int MinHairpinLoop { get; init; } = 3;

    /// <summary>Whether to look for long exact direct repeats.</summary>
    public bool AvoidRepeats { get; init; } = true;

    /// <summary>Longest exact direct repeat that is still acceptable; anything longer is reported once.</summary>
    public int MaxDirectRepeat { get; init; } = 15;

    /// <summary>
    /// The everyday default: no forbidden sites, homopolymers up to 6, GC between 30% and 70% over 50-base
    /// windows, hairpin stems up to 9 bases and direct repeats up to 15.
    /// </summary>
    public static SequenceConstraints Default { get; } = new();

    /// <summary>
    /// The stricter set a gene synthesis vendor typically applies: all of <see cref="CommonRestrictionSites"/>
    /// forbidden, homopolymers capped at 5, GC held between 35% and 65% over 50-base windows, hairpin stems
    /// capped at 8 and direct repeats at 12.
    /// </summary>
    public static SequenceConstraints ForGeneSynthesis { get; } = new()
    {
        ForbiddenSites = CommonRestrictionSites,
        MaxHomopolymer = 5,
        MinGc = 0.35,
        MaxGc = 0.65,
        GcWindow = 50,
        AvoidHairpins = true,
        MinHairpinStem = 9,
        MinHairpinLoop = 3,
        AvoidRepeats = true,
        MaxDirectRepeat = 12,
    };

    /// <summary>
    /// Finds every place the sequence breaks one of the rules, ordered by position. An empty list means the
    /// sequence is ready to order.
    /// </summary>
    /// <remarks>
    /// Cost: the site, homopolymer, GC and repeat scans are linear in the sequence length (the repeat scan
    /// hashes one k-mer per position). The hairpin scan is the expensive one — it tries every pair of positions
    /// that could close a stem, so it is O(n²) in the worst case, with an early exit on the first unpaired base
    /// of each candidate, which makes the constant small in practice. Switch <see cref="AvoidHairpins"/> off
    /// for sequences of many kilobases if the check dominates the design loop.
    /// </remarks>
    /// <param name="dna">The designed sequence; U is read as T and case is ignored.</param>
    public IReadOnlyList<ConstraintViolation> Check(string dna)
    {
        ArgumentNullException.ThrowIfNull(dna);
        string sequence = DnaUtility.NormalizeDna(dna);
        var bases = DnaUtility.Ascii(sequence);
        var violations = new List<ConstraintViolation>();
        if (bases.Length == 0) return violations;

        CheckForbiddenSites(bases, violations);
        CheckHomopolymers(bases, violations);
        CheckGcWindows(bases, violations);
        if (AvoidHairpins) CheckHairpins(bases, violations);
        if (AvoidRepeats) CheckDirectRepeats(sequence, violations);

        violations.Sort(static (a, b) =>
        {
            int byStart = a.Start.CompareTo(b.Start);
            if (byStart != 0) return byStart;
            int byKind = string.CompareOrdinal(a.Kind, b.Kind);
            return byKind != 0 ? byKind : a.Length.CompareTo(b.Length);
        });
        return violations;
    }

    /// <summary>Whether <see cref="Check(string)"/> finds nothing to complain about.</summary>
    public bool IsSatisfiedBy(string dna) => Check(dna).Count == 0;

    private void CheckForbiddenSites(byte[] bases, List<ConstraintViolation> violations)
    {
        for (int s = 0; s < ForbiddenSites.Count; s++)
        {
            string site = DnaUtility.NormalizeDna(ForbiddenSites[s]);
            if (site.Length == 0 || site.Length > bases.Length) continue;
            var pattern = DnaUtility.Ascii(site);
            string reverse = DnaUtility.ReverseComplement(pattern);
            var reversePattern = DnaUtility.Ascii(reverse);
            bool palindromic = string.Equals(site, reverse, StringComparison.Ordinal);

            for (int i = 0; i + pattern.Length <= bases.Length; i++)
            {
                var window = bases.AsSpan(i, pattern.Length);
                if (DnaUtility.MatchesIupac(window, pattern))
                {
                    violations.Add(new ConstraintViolation(
                        KindForbiddenSite, i, pattern.Length, $"forbidden site {site} on the sense strand"));
                }
                else if (!palindromic && DnaUtility.MatchesIupac(window, reversePattern))
                {
                    violations.Add(new ConstraintViolation(
                        KindForbiddenSite, i, pattern.Length, $"forbidden site {site} on the antisense strand"));
                }
            }
        }
    }

    private void CheckHomopolymers(byte[] bases, List<ConstraintViolation> violations)
    {
        if (MaxHomopolymer <= 0) return;
        int runStart = 0;
        for (int i = 1; i <= bases.Length; i++)
        {
            bool same = i < bases.Length && (bases[i] & 0xDF) == (bases[runStart] & 0xDF);
            if (same) continue;
            int length = i - runStart;
            if (length > MaxHomopolymer)
            {
                violations.Add(new ConstraintViolation(
                    KindHomopolymer, runStart, length,
                    $"{length}×{(char)bases[runStart]} run exceeds the limit of {MaxHomopolymer}"));
            }
            runStart = i;
        }
    }

    private void CheckGcWindows(byte[] bases, List<ConstraintViolation> violations)
    {
        int window = Math.Min(GcWindow <= 0 ? bases.Length : GcWindow, bases.Length);
        if (window <= 0) return;

        // One rolling GC count over the sequence; consecutive offending windows are merged into one violation so
        // that a single GC-rich stretch is reported once rather than once per offset.
        int gc = DnaUtility.CountGc(bases.AsSpan(0, window));
        int openStart = -1, openEnd = -1;
        bool openHigh = false;
        double extreme = 0;

        for (int i = 0; ; i++)
        {
            double fraction = (double)gc / window;
            bool high = fraction > MaxGc, low = fraction < MinGc;
            if (high || low)
            {
                if (openStart >= 0 && openHigh == high && i <= openEnd)
                {
                    openEnd = i + window;
                    extreme = high ? Math.Max(extreme, fraction) : Math.Min(extreme, fraction);
                }
                else
                {
                    FlushGcWindow(violations, openStart, openEnd, openHigh, extreme);
                    openStart = i;
                    openEnd = i + window;
                    openHigh = high;
                    extreme = fraction;
                }
            }

            if (i + window >= bases.Length) break;
            if ((bases[i] & 0xDF) is (byte)'G' or (byte)'C') gc--;
            if ((bases[i + window] & 0xDF) is (byte)'G' or (byte)'C') gc++;
        }
        FlushGcWindow(violations, openStart, openEnd, openHigh, extreme);
    }

    private void FlushGcWindow(List<ConstraintViolation> violations, int start, int end, bool high, double extreme)
    {
        if (start < 0) return;
        violations.Add(new ConstraintViolation(
            KindGcWindow, start, end - start,
            high
                ? $"GC {extreme:P1} over a {Math.Min(GcWindow, end - start)}-base window exceeds the limit of {MaxGc:P1}"
                : $"GC {extreme:P1} over a {Math.Min(GcWindow, end - start)}-base window is below the limit of {MinGc:P1}"));
    }

    private void CheckHairpins(byte[] bases, List<ConstraintViolation> violations)
    {
        int stem = Math.Max(1, MinHairpinStem);
        int loop = Math.Max(0, MinHairpinLoop);
        int n = bases.Length;
        int span = 2 * stem + loop;
        if (n < span) return;

        for (int i = 0; i + span <= n; i++)
        {
            for (int j = i + span - 1; j < n; j++)
            {
                if (!Pairs(bases[i], bases[j])) continue;
                // Report maximal stems only: if the flanking pair also pairs, this stem is part of a longer one
                // that will be (or was) reported from its own outer start.
                if (i > 0 && j + 1 < n && Pairs(bases[i - 1], bases[j + 1])) continue;

                int k = 1;
                while (i + k < j - k && Pairs(bases[i + k], bases[j - k])) k++;
                int loopLength = j - i + 1 - 2 * k;
                if (k >= stem && loopLength >= loop)
                {
                    violations.Add(new ConstraintViolation(
                        KindHairpin, i, j - i + 1,
                        $"{k}-base stem with a {loopLength}-base loop (limit {stem})"));
                }
            }
        }
    }

    private static bool Pairs(byte left, byte right)
    {
        int l = left & 0xDF, r = right & 0xDF;
        return (l == 'A' && r == 'T') || (l == 'T' && r == 'A')
            || (l == 'G' && r == 'C') || (l == 'C' && r == 'G');
    }

    private void CheckDirectRepeats(string sequence, List<ConstraintViolation> violations)
    {
        int k = MaxDirectRepeat + 1;
        if (k <= 0 || sequence.Length < 2 * k) return;

        var first = new Dictionary<string, int>(sequence.Length, StringComparer.Ordinal);
        int lastReportedEnd = -1;
        for (int i = 0; i + k <= sequence.Length; i++)
        {
            string word = sequence.Substring(i, k);
            if (first.TryGetValue(word, out int origin))
            {
                if (i >= lastReportedEnd)
                {
                    violations.Add(new ConstraintViolation(
                        KindDirectRepeat, i, k,
                        $"repeats the {k} bases at {origin} (limit {MaxDirectRepeat})"));
                    lastReportedEnd = i + k;
                }
            }
            else
            {
                first[word] = i;
            }
        }
    }
}
