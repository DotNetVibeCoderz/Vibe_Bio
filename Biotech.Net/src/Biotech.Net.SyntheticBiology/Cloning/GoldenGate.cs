using System.Text;

namespace Biotech.Net.SyntheticBiology.Cloning;

/// <summary>
/// Golden Gate / MoClo assembly: one-pot cloning in which a type IIS enzyme cuts outside its own recognition
/// site, so every part can be released with a freely chosen 4 nt overhang and the enzyme sites are lost from
/// the product.
/// </summary>
/// <remarks>
/// <para>
/// This class covers the three things a designer actually needs: choosing a mutually orthogonal overhang set
/// (<see cref="DesignOverhangs"/>), scoring how well such a set should ligate
/// (<see cref="OverhangFidelity"/>), adding correctly oriented flanking sites to a part
/// (<see cref="AddFlanks"/>), and simulating the digestion-ligation reaction
/// (<see cref="Simulate(IEnumerable{string}, RestrictionEnzyme, string?, bool)"/>).
/// </para>
/// <para>
/// <b>Geometry.</b> A flanked part reads
/// <c>pad · SITE · spacer · upstream overhang · part · downstream overhang · spacer′ · SITE′ · pad′</c>,
/// where <c>SITE′</c> is the reverse complement of the recognition site and the spacer length is the enzyme's
/// top-strand cut offset (1 nt for BsaI <c>GGTCTC(1/5)</c>). Both sites therefore point <i>inward</i>: the
/// enzyme cuts between its site and the part, the sites stay on the discarded outer pieces, and the released
/// fragment is exactly <c>upstream overhang · part · downstream overhang</c>. Overhang bookkeeping follows the
/// convention documented on <see cref="DnaFragment"/>: an overhang is the top-strand text of the window it
/// occupies, shared by the two fragments that meet there.
/// </para>
/// </remarks>
public static class GoldenGate
{
    /// <summary>
    /// Cuts a linear or circular molecule with one type IIS (or type IIP) enzyme and returns the fragments in
    /// order along the molecule.
    /// </summary>
    /// <remarks>
    /// Overhangs follow the <see cref="DnaFragment"/> convention: the single-stranded window of each cut is
    /// written as top-strand text and appears in both of the fragments that meet at that cut. A linear
    /// molecule with <i>n</i> usable cuts yields <i>n</i>+1 fragments, the first and last of which keep a blunt
    /// outer end; a circle yields <i>n</i>. Cuts whose nicks would fall outside a linear molecule (a type IIS
    /// site too close to an end) are ignored.
    /// </remarks>
    /// <param name="dna">Top strand of the molecule; case and U/T are normalised.</param>
    /// <param name="enzyme">The enzyme to digest with.</param>
    /// <param name="circular">Whether the molecule is a circle.</param>
    /// <param name="id">Identifier stem for the fragments; defaults to <c>frag</c>.</param>
    /// <returns>The fragments, in top-strand order (for a circle, starting at the first cut).</returns>
    public static IReadOnlyList<DnaFragment> Digest(string dna, RestrictionEnzyme enzyme, bool circular = false, string? id = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(dna);
        ArgumentNullException.ThrowIfNull(enzyme);
        string top = DnaUtility.NormalizeDna(dna);
        string stem = id ?? "frag";

        if (circular) return new Plasmid(stem, top).Fragments(enzyme);

        int overhang = Math.Abs(enzyme.OverhangLength);
        var edges = new List<(int Lo, int Hi)>();
        foreach (var site in enzyme.Find(top, circular: false))
        {
            int lo = enzyme.OverhangLength < 0 ? site.TopCutPosition - overhang : site.TopCutPosition;
            int hi = lo + overhang;
            if (lo < 0 || hi > top.Length) continue;
            edges.Add((lo, hi));
        }
        edges.Sort(static (a, b) => a.Lo.CompareTo(b.Lo));

        if (edges.Count == 0) return [new DnaFragment(stem, top)];

        var fragments = new List<DnaFragment>(edges.Count + 1);
        int cursor = 0;
        string carried = string.Empty;
        for (int i = 0; i < edges.Count; i++)
        {
            var (lo, hi) = edges[i];
            if (lo < cursor) continue; // overlapping cuts: keep the first one
            string sequence = top[cursor..hi];
            string three = overhang == 0 ? string.Empty : top[lo..hi];
            fragments.Add(new DnaFragment($"{stem}_{fragments.Count + 1}", sequence, carried, three));
            carried = three;
            cursor = lo;
        }
        fragments.Add(new DnaFragment($"{stem}_{fragments.Count + 1}", top[cursor..], carried, string.Empty));
        return fragments;
    }

    /// <summary>Whether a part is free of internal recognition sites for an enzyme ("domesticated").</summary>
    /// <remarks>Both strands are searched, so a site on the bottom strand also fails the check.</remarks>
    /// <param name="part">The part sequence to test.</param>
    /// <param name="enzyme">The assembly enzyme the part must avoid.</param>
    public static bool IsDomesticated(string part, RestrictionEnzyme enzyme)
    {
        ArgumentException.ThrowIfNullOrEmpty(part);
        ArgumentNullException.ThrowIfNull(enzyme);
        return enzyme.Find(part, circular: false).Count == 0;
    }

    /// <summary>
    /// Flanks a part with two inward-pointing type IIS sites, so that the enzyme releases
    /// <paramref name="upstreamOverhang"/> + part + <paramref name="downstreamOverhang"/> and the sites are
    /// lost with the outer pieces.
    /// </summary>
    /// <remarks>
    /// The spacer between a site and its overhang is <see cref="RestrictionEnzyme.CutTop"/> bases long, which
    /// is what puts the top-strand nick exactly at the start of the overhang. The downstream site is written as
    /// the reverse complement of the recognition site, which is what makes it read inward on the bottom strand.
    /// </remarks>
    /// <param name="part">The (domesticated) part sequence.</param>
    /// <param name="enzyme">A type IIS enzyme with a 5′ overhang, e.g. BsaI or BsmBI.</param>
    /// <param name="upstreamOverhang">Overhang to expose at the 5′ end; must match the enzyme's overhang length.</param>
    /// <param name="downstreamOverhang">Overhang to expose at the 3′ end; must match the enzyme's overhang length.</param>
    /// <param name="spacerBase">Base used to pad between site and overhang; defaults to <c>A</c>.</param>
    /// <param name="pad">Outer padding placed beyond each site, so the construct can be amplified; may be empty.</param>
    /// <returns>The flanked part, ready for a one-pot reaction.</returns>
    /// <exception cref="ArgumentException">
    /// The enzyme is not a 5′-overhang type IIS enzyme, an overhang has the wrong length, or the part is not
    /// domesticated for this enzyme.
    /// </exception>
    public static string AddFlanks(
        string part,
        RestrictionEnzyme enzyme,
        string upstreamOverhang,
        string downstreamOverhang,
        char spacerBase = 'A',
        string pad = "TT")
    {
        ArgumentException.ThrowIfNullOrEmpty(part);
        ArgumentNullException.ThrowIfNull(enzyme);
        ArgumentNullException.ThrowIfNull(upstreamOverhang);
        ArgumentNullException.ThrowIfNull(downstreamOverhang);
        ArgumentNullException.ThrowIfNull(pad);
        if (!enzyme.IsTypeIis || enzyme.OverhangLength <= 0)
            throw new ArgumentException($"{enzyme.Name} is not a type IIS enzyme with a 5' overhang.", nameof(enzyme));

        int overhang = enzyme.OverhangLength;
        string up = DnaUtility.NormalizeDna(upstreamOverhang);
        string down = DnaUtility.NormalizeDna(downstreamOverhang);
        if (up.Length != overhang || down.Length != overhang)
            throw new ArgumentException($"{enzyme.Name} leaves a {overhang} nt overhang.", nameof(upstreamOverhang));

        string body = DnaUtility.NormalizeDna(part);
        if (!IsDomesticated(body, enzyme))
            throw new ArgumentException($"The part contains an internal {enzyme.Name} site and must be domesticated first.", nameof(part));

        string site = DnaUtility.NormalizeDna(enzyme.RecognitionSite);
        string spacer = new(spacerBase, enzyme.CutTop);
        string padTop = DnaUtility.NormalizeDna(pad);

        var builder = new StringBuilder();
        builder.Append(padTop).Append(site).Append(spacer).Append(up);
        builder.Append(body);
        builder.Append(down).Append(DnaUtility.ReverseComplement(spacer)).Append(DnaUtility.ReverseComplement(site))
               .Append(DnaUtility.ReverseComplement(padTop));
        return builder.ToString();
    }

    /// <summary>
    /// Chooses a set of mutually orthogonal overhangs.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Accepted overhangs are pairwise orthogonal under four rules: none is palindromic (it would ligate to
    /// itself), none equals another or the reverse complement of another, every pair differs in at least two
    /// positions (Hamming distance ≥ 2) both directly and against the other's reverse complement, and all-one-base
    /// overhangs are skipped. The search walks the full 4<sup>L</sup> space in an order shuffled with
    /// <paramref name="seed"/>, so a given seed always gives the same set.
    /// </para>
    /// </remarks>
    /// <param name="count">How many overhangs are needed (one per junction).</param>
    /// <param name="overhangLength">Overhang length in nucleotides; 4 for BsaI and friends.</param>
    /// <param name="seed">Seed of the deterministic shuffle.</param>
    /// <param name="reserved">Overhangs that are already in use and must also be avoided; may be <see langword="null"/>.</param>
    /// <returns><paramref name="count"/> overhangs, upper case.</returns>
    /// <exception cref="ArgumentOutOfRangeException">The count or the length is not positive.</exception>
    /// <exception cref="InvalidOperationException">No orthogonal set of that size exists for these constraints.</exception>
    public static IReadOnlyList<string> DesignOverhangs(int count, int overhangLength = 4, int seed = 20250101, IEnumerable<string>? reserved = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(count);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(overhangLength);
        if (overhangLength > 8) throw new ArgumentOutOfRangeException(nameof(overhangLength), "Overhangs longer than 8 nt are not supported.");

        var accepted = new List<string>();
        if (reserved is not null)
            foreach (var item in reserved) accepted.Add(DnaUtility.NormalizeDna(item));

        int reservedCount = accepted.Count;
        int total = 1 << (2 * overhangLength);
        var order = new int[total];
        for (int i = 0; i < total; i++) order[i] = i;
        var random = new Random(seed);
        for (int i = total - 1; i > 0; i--)
        {
            int j = random.Next(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }

        var chosen = new List<string>(count);
        foreach (int code in order)
        {
            if (chosen.Count == count) break;
            string candidate = DnaUtility.Unpack((ulong)code, overhangLength);
            if (IsUniformBase(candidate)) continue;
            if (!IsOrthogonalAgainst(candidate, accepted)) continue;
            accepted.Add(candidate);
            chosen.Add(candidate);
        }

        if (chosen.Count < count)
            throw new InvalidOperationException(
                $"Only {chosen.Count} orthogonal {overhangLength} nt overhangs could be found next to {reservedCount} reserved ones; {count} were requested.");
        return chosen;
    }

    private static bool IsUniformBase(string overhang)
    {
        for (int i = 1; i < overhang.Length; i++)
            if (overhang[i] != overhang[0]) return false;
        return true;
    }

    private static bool IsOrthogonalAgainst(string candidate, IReadOnlyList<string> accepted)
    {
        string selfRc = DnaUtility.ReverseComplement(candidate);
        if (selfRc == candidate) return false; // palindromic: ligates to itself
        foreach (var other in accepted)
        {
            if (other.Length != candidate.Length) continue;
            if (Hamming(candidate, other) < 2) return false;
            if (Hamming(candidate, DnaUtility.ReverseComplement(other)) < 2) return false;
        }
        return true;
    }

    private static int Hamming(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    /// <summary>
    /// Checks a set of overhangs for mutual orthogonality and explains the first violation.
    /// </summary>
    /// <param name="overhangs">The overhangs to test; all must have the same length.</param>
    /// <param name="reason">The first rule that is broken, or <see langword="null"/> when the set is orthogonal.</param>
    /// <returns>True when no rule is broken.</returns>
    public static bool AreOrthogonal(IReadOnlyList<string> overhangs, out string? reason)
    {
        ArgumentNullException.ThrowIfNull(overhangs);
        var normalised = new List<string>(overhangs.Count);
        foreach (var overhang in overhangs) normalised.Add(DnaUtility.NormalizeDna(overhang));

        for (int i = 0; i < normalised.Count; i++)
        {
            string a = normalised[i];
            if (a.Length != normalised[0].Length)
            {
                reason = $"Overhang '{a}' has a different length from '{normalised[0]}'.";
                return false;
            }
            if (DnaUtility.ReverseComplement(a) == a)
            {
                reason = $"Overhang '{a}' is palindromic and would ligate to itself.";
                return false;
            }
            for (int j = i + 1; j < normalised.Count; j++)
            {
                string b = normalised[j];
                if (a == b) { reason = $"Overhang '{a}' is used twice."; return false; }
                if (DnaUtility.ReverseComplement(a) == b)
                {
                    reason = $"Overhangs '{a}' and '{b}' are reverse complements of each other.";
                    return false;
                }
                if (Hamming(a, b) < 2)
                {
                    reason = $"Overhangs '{a}' and '{b}' differ in only one position.";
                    return false;
                }
                if (Hamming(a, DnaUtility.ReverseComplement(b)) < 2)
                {
                    reason = $"Overhang '{a}' is within one mismatch of the complement of '{b}'.";
                    return false;
                }
            }
        }
        reason = null;
        return true;
    }

    /// <summary>
    /// A ligation-fidelity <b>proxy</b> for an overhang set: 1.0 for a perfectly orthogonal set, lower when
    /// overhangs can cross-ligate.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>This is not the published NEB / Potapov–Pryor T4 DNA ligase fidelity dataset.</b> Those values come
    /// from a measured mismatch-frequency matrix that is not reproduced here; nothing in this method is
    /// experimental data and the numbers must not be reported as such. It is a transparent, self-contained
    /// model meant for ranking candidate sets against each other.
    /// </para>
    /// <para>
    /// <b>The model.</b> Overhang <i>a</i> ligates correctly to the complement of <i>a</i>. It can mis-ligate
    /// to the complement of another overhang <i>b</i> (or of <i>b</i>'s reverse complement, i.e. to <i>b</i> in
    /// the other orientation), and the comparison reduces to the Hamming distance between the two overhang
    /// strings. Only near-identical ends are treated as competitors: distance 0 scores a full-strength
    /// competing event (weight 1), distance 1 scores a weak one whose weight rises with GC content
    /// (0.15 + 0.25 × GC fraction, because a G·C-rich mismatched duplex still anneals), and distance ≥ 2 scores
    /// nothing. A palindromic overhang additionally competes with itself. Per overhang the fidelity is
    /// 1/(1 + Σ competing weights), and the set score is the arithmetic mean over all overhangs.
    /// </para>
    /// </remarks>
    /// <param name="overhangs">The overhang set, all of the same length.</param>
    /// <returns>A score in (0, 1]; exactly 1.0 when no overhang competes with any other.</returns>
    public static double OverhangFidelity(IReadOnlyList<string> overhangs)
    {
        ArgumentNullException.ThrowIfNull(overhangs);
        if (overhangs.Count == 0) return 1.0;
        var set = new List<string>(overhangs.Count);
        foreach (var overhang in overhangs) set.Add(DnaUtility.NormalizeDna(overhang));

        double sum = 0;
        foreach (var a in set)
        {
            double competing = MispairWeight(a, a, selfPair: true);
            foreach (var b in set)
            {
                if (ReferenceEquals(a, b)) continue;
                if (a.Length != b.Length) continue;
                competing += MispairWeight(a, b, selfPair: false);
                competing += MispairWeight(a, DnaUtility.ReverseComplement(b), selfPair: false);
            }
            sum += 1.0 / (1.0 + competing);
        }
        return sum / set.Count;
    }

    // Weight of the competing ligation between overhang 'a' and the complement of 'b'. Documented in full on
    // OverhangFidelity: it is a Hamming/GC proxy, not measured ligase data.
    private static double MispairWeight(string a, string b, bool selfPair)
    {
        if (a.Length != b.Length) return 0;
        int distance = selfPair ? Hamming(a, DnaUtility.ReverseComplement(a)) : Hamming(a, b);
        if (distance >= 2) return 0;
        if (distance == 0) return 1.0;
        double gc = 0.5 * (DnaUtility.GcFraction(a) + DnaUtility.GcFraction(b));
        return 0.15 + 0.25 * gc;
    }

    /// <summary>
    /// Simulates a one-pot Golden Gate reaction: every input is digested with the type IIS enzyme, the
    /// released fragments are joined by matching overhangs, and the circular product is returned.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A part must carry exactly two sites, one at each end and both pointing inward (see
    /// <see cref="AddFlanks"/>); more than two means the part is not domesticated. The backbone is treated the
    /// same way, except that it is normally a circle whose two sites sit in the stuffer, so the released
    /// backbone fragment is the one without any remaining site.
    /// </para>
    /// <para>
    /// Fragments are chained on exact overhang matches, starting from the backbone (or from the first part when
    /// no backbone is given). Two fragments offering the same 5′ overhang make the reaction ambiguous and are
    /// reported as a failure. The product is circular, contains every released fragment exactly once, and
    /// carries no recognition site for the enzyme.
    /// </para>
    /// </remarks>
    /// <param name="parts">Flanked part sequences, in any order.</param>
    /// <param name="enzyme">The type IIS enzyme used in the reaction.</param>
    /// <param name="backbone">Destination vector; <see langword="null"/> assembles the parts alone.</param>
    /// <param name="backboneIsCircular">Whether the backbone sequence is a circle.</param>
    /// <returns>The circular product with its junction overhangs, or a failed result carrying the reason.</returns>
    public static GoldenGateResult Simulate(
        IEnumerable<string> parts,
        RestrictionEnzyme enzyme,
        string? backbone = null,
        bool backboneIsCircular = true)
    {
        ArgumentNullException.ThrowIfNull(parts);
        ArgumentNullException.ThrowIfNull(enzyme);

        var released = new List<DnaFragment>();
        if (backbone is not null)
        {
            var outcome = Release(backbone, enzyme, backboneIsCircular, "backbone", released);
            if (outcome is not null) return GoldenGateResult.Failed(outcome);
        }

        int index = 0;
        foreach (var part in parts)
        {
            index++;
            if (string.IsNullOrEmpty(part)) return GoldenGateResult.Failed($"Part {index} is empty.");
            var outcome = Release(part, enzyme, circular: false, $"part{index}", released);
            if (outcome is not null) return GoldenGateResult.Failed(outcome);
        }

        if (released.Count == 0)
            return GoldenGateResult.Failed($"No fragment was released: nothing carried a usable {enzyme.Name} site.");

        var byFivePrime = new Dictionary<string, DnaFragment>(StringComparer.Ordinal);
        foreach (var fragment in released)
        {
            if (byFivePrime.ContainsKey(fragment.FivePrimeOverhang))
                return GoldenGateResult.Failed(
                    $"Ambiguous overhang '{fragment.FivePrimeOverhang}': it is offered by more than one fragment.");
            byFivePrime[fragment.FivePrimeOverhang] = fragment;
        }

        var start = released[0];
        var current = start;
        var junctions = new List<string>();
        var used = new HashSet<string>(StringComparer.Ordinal) { start.Id };
        var chain = new List<DnaFragment> { start };

        while (true)
        {
            if (!byFivePrime.TryGetValue(current.ThreePrimeOverhang, out var next))
                return GoldenGateResult.Failed(
                    $"No fragment accepts the overhang '{current.ThreePrimeOverhang}' left by '{current.Id}'.");
            junctions.Add(current.ThreePrimeOverhang);
            if (ReferenceEquals(next, start)) break;
            if (!used.Add(next.Id))
                return GoldenGateResult.Failed($"Fragment '{next.Id}' would be used twice.");
            chain.Add(next);
            current = next;
        }

        if (chain.Count != released.Count)
            return GoldenGateResult.Failed(
                $"Only {chain.Count} of {released.Count} released fragments form a circle; the overhangs do not close.");

        var product = chain[0];
        for (int i = 1; i < chain.Count; i++) product = product.Ligate(chain[i], $"{product.Id}+{chain[i].Id}");
        product = product.Circularise() with { Id = "goldengate_product" };

        if (enzyme.Find(product.Sequence, circular: true).Count != 0)
            return GoldenGateResult.Failed($"The product still contains a {enzyme.Name} site.");

        return new GoldenGateResult(true, product, junctions, null);
    }

    // Digests one input and appends its usable (two-sided, site-free) fragments to 'released'. Returns a
    // failure reason, or null on success.
    private static string? Release(string dna, RestrictionEnzyme enzyme, bool circular, string id, List<DnaFragment> released)
    {
        int sites = enzyme.Find(DnaUtility.NormalizeDna(dna), circular).Count;
        if (sites == 0) return $"'{id}' contains no {enzyme.Name} site.";
        if (sites > 2) return $"'{id}' contains {sites} {enzyme.Name} sites and is not domesticated.";
        if (sites < 2) return $"'{id}' contains only one {enzyme.Name} site, so nothing is released.";

        var usable = new List<DnaFragment>();
        foreach (var fragment in Digest(dna, enzyme, circular, id))
        {
            if (fragment.FivePrimeOverhang.Length == 0 || fragment.ThreePrimeOverhang.Length == 0) continue;
            if (enzyme.Find(fragment.Sequence, circular: false).Count != 0) continue;
            usable.Add(fragment);
        }
        if (usable.Count == 0) return $"'{id}' released no site-free fragment with two sticky ends.";
        if (usable.Count > 1) return $"'{id}' released {usable.Count} fragments; expected exactly one.";
        released.Add(usable[0]);
        return null;
    }
}

/// <summary>Outcome of <see cref="GoldenGate.Simulate(IEnumerable{string}, RestrictionEnzyme, string?, bool)"/>.</summary>
/// <param name="Success">Whether a single circular product was assembled.</param>
/// <param name="Product">The circular product, or <see langword="null"/> when assembly failed.</param>
/// <param name="JunctionOverhangs">The overhangs consumed, in the order they were used around the circle.</param>
/// <param name="FailureReason">
/// Why the reaction failed — an ambiguous overhang, a missing site, an undomesticated part or an open chain —
/// or <see langword="null"/> on success.
/// </param>
public sealed record GoldenGateResult(
    bool Success,
    DnaFragment? Product,
    IReadOnlyList<string> JunctionOverhangs,
    string? FailureReason)
{
    /// <summary>A failed result carrying only a reason.</summary>
    /// <param name="reason">Why the reaction could not be resolved to one circular product.</param>
    public static GoldenGateResult Failed(string reason) => new(false, null, [], reason);

    /// <summary>Readable summary, e.g. <c>Golden Gate: 2_104 bp circular, 4 junctions</c>.</summary>
    public override string ToString() => Success
        ? $"Golden Gate: {Product!.Length} bp circular, {JunctionOverhangs.Count} junctions"
        : $"Golden Gate failed: {FailureReason}";
}
