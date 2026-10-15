using System.Text;

namespace Biotech.Net.SyntheticBiology.Primers;

/// <summary>
/// The most stable ungapped duplex found between two oligonucleotides (or an oligonucleotide and itself).
/// </summary>
/// <remarks>
/// Positions are 0-based indices into the strands as given (5′→3′). Top-strand base <c>i</c> pairs with
/// bottom-strand base <c>j</c> where <c>i + j = </c><see cref="Offset"/>, so one offset is one antiparallel register.
/// When no duplex with ΔG°37 &lt; 0 exists the record is <see cref="None"/>: ΔG 0, length 0, offset −1, empty diagram.
/// </remarks>
/// <param name="DeltaG37">ΔG°37 of the duplex in kcal/mol (negative = stable); 0 when none forms.</param>
/// <param name="Offset">The antiparallel register <c>i + j</c> of the duplex, or −1 for <see cref="None"/>.</param>
/// <param name="TopStart">Index in the top strand of the 5′-most paired base of the duplex.</param>
/// <param name="BottomStart">Index in the bottom strand of the 5′-most paired base of the duplex.</param>
/// <param name="Length">Number of consecutive Watson–Crick base pairs.</param>
/// <param name="Diagram">
/// A 3-line ASCII alignment: the top strand 5′→3′, a line with <c>|</c> under every scored base pair, and the
/// bottom strand written 3′→5′. Lines are separated by <c>'\n'</c>.
/// </param>
public sealed record OligoDuplex(double DeltaG37, int Offset, int TopStart, int BottomStart, int Length, string Diagram)
{
    /// <summary>No stable duplex (ΔG°37 = 0).</summary>
    public static OligoDuplex None { get; } = new(0.0, -1, -1, -1, 0, string.Empty);

    /// <summary>Whether a duplex with negative ΔG°37 was found.</summary>
    public bool Found => Length > 0;
}

/// <summary>
/// The most stable simple stem-loop (hairpin) of an oligonucleotide.
/// </summary>
/// <remarks>
/// The 5′ arm of the stem is <c>[StemStart, StemStart + StemLength)</c>, the loop follows it, and the 3′ arm is the
/// next <see cref="StemLength"/> bases, paired antiparallel with the 5′ arm. When no hairpin with ΔG°37 &lt; 0 exists
/// the record is <see cref="None"/>.
/// </remarks>
/// <param name="DeltaG37">ΔG°37 in kcal/mol (negative = stable); 0 when none forms.</param>
/// <param name="StemStart">0-based index of the outermost 5′ stem base, or −1 for <see cref="None"/>.</param>
/// <param name="StemLength">Base pairs in the stem.</param>
/// <param name="LoopLength">Unpaired nucleotides in the loop.</param>
/// <param name="Structure">Dot-bracket notation of the hairpin over the whole oligo (empty for <see cref="None"/>).</param>
public sealed record OligoHairpin(double DeltaG37, int StemStart, int StemLength, int LoopLength, string Structure)
{
    /// <summary>No stable hairpin (ΔG°37 = 0).</summary>
    public static OligoHairpin None { get; } = new(0.0, -1, 0, 0, string.Empty);

    /// <summary>Whether a hairpin with negative ΔG°37 was found.</summary>
    public bool Found => StemLength > 0;
}

/// <summary>
/// Secondary-structure screening of oligonucleotides — self-dimers, hetero-dimers, 3′-anchored dimers and
/// hairpins — scored with the SantaLucia (1998) nearest-neighbour ΔG°37 parameters, in the spirit of the
/// "any"/"end" alignment scores of Primer3 and the dimer/hairpin panels of IDT OligoAnalyzer.
/// </summary>
/// <remarks>
/// <para>
/// <b>These are screening heuristics, not a structure prediction.</b> Only perfectly matched, contiguous
/// (ungapped) Watson–Crick helices are considered: there are no mismatches, bulges, internal loops, dangling
/// ends, terminal mismatches, coaxial stacking or multi-branch structures, and no minimum-free-energy folding
/// (no Zuker/McCaskill dynamic programming). Real dimers and hairpins can therefore be somewhat more stable
/// than reported here. The numbers are nonetheless on the same kcal/mol scale as Primer3/OligoAnalyzer and
/// are what primer design needs: a fast, monotone measure of how much an oligo wants to pair with itself or
/// its partner.
/// </para>
/// <para>
/// <b>Dimer model.</b> The two strands are slid antiparallel past each other over every register. In each
/// register every maximal run of ≥ 2 consecutive complementary pairs is scored as an independent duplex:
/// ΔG°37 = Σ nearest-neighbour stacks (SantaLucia 1998 Table 1) + one initiation term per helix end
/// (terminal G·C: +0.98 kcal/mol; terminal A·T: +1.03 kcal/mol at 37 °C). For a self-dimer whose helix is
/// symmetric (the two strands contribute the same bases) the symmetry correction ΔS° = −1.4 e.u.
/// (+0.43 kcal/mol) is added, exactly as in <see cref="NearestNeighbourTm.Duplex(string)"/>. Consequently the
/// hetero-dimer of an oligo with its exact reverse complement equals
/// <c>NearestNeighbourTm.Duplex(oligo).FreeEnergy37</c>. Non-ACGT symbols (N, IUPAC codes) never pair.
/// </para>
/// <para>
/// <b>Hairpin model.</b> A stem of L ≥ <c>minStem</c> consecutive pairs closing a loop of n ≥ <c>minLoop</c>
/// nucleotides scores ΔG°37 = Σ (L − 1) stacks + ΔG°37<sub>loop</sub>(n), with no initiation term (the
/// transition is unimolecular). The hairpin-loop penalties are those of SantaLucia &amp; Hicks (2004)
/// "The thermodynamics of DNA structural motifs", Annu. Rev. Biophys. Biomol. Struct. 33:415–440, Table 4
/// (kcal/mol at 37 °C, 1 M NaCl): n = 3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0, 7: 4.2, 8: 4.3, 9: 4.5, 10: 4.6,
/// 12: 5.0, 14: 5.1, 16: 5.3, 18: 5.5, 20: 5.7, 25: 6.1, 30: 6.3. Sizes between tabulated values, and above 30,
/// use the Jacobson–Stockmayer extrapolation of the same paper from the largest tabulated size x ≤ n:
/// ΔG°37(n) = ΔG°37(x) + 2.44·R·T·ln(n/x), with R = 1.987 cal/(mol·K) and T = 310.15 K
/// (2.44·R·T = 1.504 kcal/mol). Loop-sequence bonuses (tetraloops, triloops), terminal mismatches and the
/// terminal A·T penalty are ignored.
/// </para>
/// <para>All energies are at the 1 M Na⁺ reference state of the nearest-neighbour parameters; no salt correction is applied.</para>
/// </remarks>
public static class OligoThermodynamics
{
    private const byte NonPairing = 4;

    /// <summary>ΔG°37 of the self-complementary symmetry correction: −T·ΔS° with ΔS° = −1.4 e.u.</summary>
    private const double SymmetryDeltaG37 = 310.15 * 1.4 / 1000.0;

    private const double LoopExtrapolation = 2.44 * NearestNeighbourTm.GasConstant * 310.15 / 1000.0;

    private static readonly int[] LoopSizes = [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30];
    private static readonly double[] LoopDeltaG = [3.5, 3.5, 3.3, 4.0, 4.2, 4.3, 4.5, 4.6, 5.0, 5.1, 5.3, 5.5, 5.7, 6.1, 6.3];

    /// <summary>
    /// The most stable ungapped duplex an oligonucleotide forms with a second copy of itself.
    /// </summary>
    /// <remarks>See the class remarks for the model and its limits.</remarks>
    /// <param name="oligo">The oligonucleotide, 5′→3′ (case-insensitive; U read as T).</param>
    /// <returns>The best duplex, or <see cref="OligoDuplex.None"/>.</returns>
    public static OligoDuplex SelfDimer(string oligo)
    {
        ArgumentNullException.ThrowIfNull(oligo);
        var codes = Encode(oligo);
        var best = BestDuplex(codes, codes, selfDimer: true, threePrimeAnchored: false);
        return ToDuplex(best, oligo, oligo);
    }

    /// <summary>
    /// The most stable ungapped duplex formed between two different oligonucleotides.
    /// </summary>
    /// <remarks>See the class remarks for the model and its limits.</remarks>
    /// <param name="a">The top strand, 5′→3′.</param>
    /// <param name="b">The bottom strand, 5′→3′ (it is reversed internally to pair antiparallel).</param>
    /// <returns>The best duplex, or <see cref="OligoDuplex.None"/>.</returns>
    public static OligoDuplex HeteroDimer(string a, string b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        var best = BestDuplex(Encode(a), Encode(b), selfDimer: false, threePrimeAnchored: false);
        return ToDuplex(best, a, b);
    }

    /// <summary>
    /// The most stable 3′-anchored dimer between <paramref name="a"/> and <paramref name="b"/>: only helices that
    /// include the 3′-terminal nucleotide of either strand are considered.
    /// </summary>
    /// <remarks>
    /// A paired 3′ end is what a polymerase extends, so this is the dimer that actually produces primer-dimer
    /// artefacts (Primer3's "end" alignment). Pass the same oligo twice for the 3′ self-dimer, in which case the
    /// symmetry correction applies as in <see cref="SelfDimer"/>.
    /// </remarks>
    /// <param name="a">The first oligonucleotide, 5′→3′.</param>
    /// <param name="b">The second oligonucleotide, 5′→3′.</param>
    /// <returns>The best 3′-anchored duplex, or <see cref="OligoDuplex.None"/>.</returns>
    public static OligoDuplex ThreePrimeDimer(string a, string b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        bool self = string.Equals(a, b, StringComparison.OrdinalIgnoreCase);
        var ca = Encode(a);
        var best = BestDuplex(ca, self ? ca : Encode(b), self, threePrimeAnchored: true);
        return ToDuplex(best, a, b);
    }

    /// <summary>The 3′-anchored self-dimer of <paramref name="oligo"/>; see <see cref="ThreePrimeDimer(string, string)"/>.</summary>
    /// <param name="oligo">The oligonucleotide, 5′→3′.</param>
    /// <returns>The best 3′-anchored self-duplex, or <see cref="OligoDuplex.None"/>.</returns>
    public static OligoDuplex ThreePrimeDimer(string oligo) => ThreePrimeDimer(oligo, oligo);

    /// <summary>
    /// The most stable simple stem-loop of <paramref name="oligo"/>: one perfectly paired, ungapped stem closing
    /// a hairpin loop. See the class remarks for the energy model (SantaLucia &amp; Hicks 2004 loop table with
    /// Jacobson–Stockmayer extrapolation).
    /// </summary>
    /// <param name="oligo">The oligonucleotide, 5′→3′.</param>
    /// <param name="minLoop">Smallest loop considered, in nucleotides (at least 3, the steric minimum).</param>
    /// <param name="minStem">Smallest stem considered, in base pairs (at least 2).</param>
    /// <returns>The best hairpin, or <see cref="OligoHairpin.None"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLoop"/> &lt; 3 or <paramref name="minStem"/> &lt; 2.</exception>
    public static OligoHairpin Hairpin(string oligo, int minLoop = 3, int minStem = 2)
    {
        ArgumentNullException.ThrowIfNull(oligo);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLoop, 3);
        ArgumentOutOfRangeException.ThrowIfLessThan(minStem, 2);
        var best = BestHairpin(Encode(oligo), minLoop, minStem);
        if (best.Length == 0) return OligoHairpin.None;

        var structure = new char[oligo.Length];
        Array.Fill(structure, '.');
        for (int k = 0; k < best.Length; k++)
        {
            structure[best.Start + k] = '(';
            structure[best.Outer3 - k] = ')';
        }
        return new OligoHairpin(best.DeltaG, best.Start, best.Length, best.Loop, new string(structure));
    }

    /// <summary>
    /// ΔG°37 (kcal/mol) of a hairpin loop of <paramref name="loopLength"/> nucleotides: the SantaLucia &amp; Hicks
    /// (2004) table with Jacobson–Stockmayer extrapolation, as documented on the class.
    /// </summary>
    /// <param name="loopLength">Loop size in nucleotides (at least 3).</param>
    /// <returns>The loop penalty in kcal/mol (positive).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="loopLength"/> &lt; 3.</exception>
    public static double HairpinLoopDeltaG37(int loopLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(loopLength, 3);
        return LoopPenalty(loopLength);
    }

    // ---------------------------------------------------------------- internals (shared with PrimerDesigner)

    internal readonly record struct DuplexHit(double DeltaG, int Offset, int TopStart, int Length);

    internal readonly record struct HairpinHit(double DeltaG, int Start, int Length, int Loop, int Outer3);

    /// <summary>2-bit codes A=0, C=1, G=2, T/U=3; anything else becomes a non-pairing code.</summary>
    internal static byte[] Encode(string sequence)
    {
        var codes = new byte[sequence.Length];
        for (int i = 0; i < codes.Length; i++)
        {
            int c = NearestNeighbourTm.CodeOf(sequence[i]);
            codes[i] = c < 0 ? NonPairing : (byte)c;
        }
        return codes;
    }

    /// <summary>
    /// Best (most negative) duplex between <paramref name="a"/> and <paramref name="b"/>; Length 0 when none is
    /// negative. Top base i pairs bottom base j on register s = i + j.
    /// </summary>
    internal static DuplexHit BestDuplex(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, bool selfDimer, bool threePrimeAnchored)
    {
        int la = a.Length, lb = b.Length;
        var best = new DuplexHit(0.0, -1, -1, 0);
        if (la < 2 || lb < 2) return best;

        for (int s = 1; s <= la + lb - 3; s++)
        {
            int iFirst = Math.Max(0, s - (lb - 1));
            int iLast = Math.Min(la - 1, s);
            int runStart = -1;
            double stacks = 0;
            for (int i = iFirst; i <= iLast + 1; i++)
            {
                bool paired = i <= iLast && a[i] + b[s - i] == 3;
                if (paired)
                {
                    if (runStart < 0) { runStart = i; stacks = 0; }
                    else stacks += NearestNeighbourTm.StackDeltaG37(a[i - 1], a[i]);
                    continue;
                }

                if (runStart >= 0 && i - runStart >= 2)
                {
                    int i0 = runStart, i1 = i - 1;
                    // 3' end of a is index la-1; 3' end of b is index lb-1, which pairs with a[s-(lb-1)].
                    if (!threePrimeAnchored || i1 == la - 1 || s - i0 == lb - 1)
                    {
                        double dg = stacks
                                    + NearestNeighbourTm.InitiationDeltaG37(a[i0])
                                    + NearestNeighbourTm.InitiationDeltaG37(a[i1]);
                        if (selfDimer && s == i0 + i1) dg += SymmetryDeltaG37;
                        if (dg < best.DeltaG) best = new DuplexHit(dg, s, i0, i1 - i0 + 1);
                    }
                }
                runStart = -1;
            }
        }
        return best;
    }

    /// <summary>Best (most negative) stem-loop; Length 0 when none is negative.</summary>
    internal static HairpinHit BestHairpin(ReadOnlySpan<byte> x, int minLoop, int minStem)
    {
        int n = x.Length;
        var best = new HairpinHit(0.0, -1, 0, 0, -1);
        // Pair (i, j) with i < j lies on register s = i + j; loop = j - i - 1 for the innermost pair.
        for (int s = minLoop + 1; s <= 2 * n - 3 - minLoop; s++)
        {
            int iFirst = Math.Max(0, s - (n - 1));
            int iLast = (s - minLoop - 1) / 2; // largest i with j - i - 1 >= minLoop
            int i = iFirst;
            while (i <= iLast)
            {
                if (x[i] + x[s - i] != 3) { i++; continue; }

                // Maximal run of consecutive pairs starting at outer pair i, growing inwards.
                int i0 = i;
                double stacks = 0;
                int k = i0;
                while (true)
                {
                    int stem = k - i0 + 1;
                    if (stem >= minStem)
                    {
                        int loop = s - 2 * k - 1;
                        double dg = stacks + LoopPenalty(loop);
                        if (dg < best.DeltaG) best = new HairpinHit(dg, i0, stem, loop, s - i0);
                    }
                    int next = k + 1;
                    if (next > iLast || x[next] + x[s - next] != 3) break;
                    stacks += NearestNeighbourTm.StackDeltaG37(x[k], x[next]);
                    k = next;
                }
                i = k + 1;
            }
        }
        return best;
    }

    /// <summary>Hairpin loop ΔG°37 for loop ≥ 3 (see class remarks).</summary>
    internal static double LoopPenalty(int loop)
    {
        int idx = Array.BinarySearch(LoopSizes, loop);
        if (idx >= 0) return LoopDeltaG[idx];
        int lower = ~idx - 1; // largest tabulated size below loop (loop >= 3 guarantees lower >= 0)
        return LoopDeltaG[lower] + LoopExtrapolation * Math.Log((double)loop / LoopSizes[lower]);
    }

    private static OligoDuplex ToDuplex(DuplexHit hit, string a, string b)
    {
        if (hit.Length == 0) return OligoDuplex.None;
        int s = hit.Offset;
        int bottomStart = s - hit.TopStart - hit.Length + 1;
        return new OligoDuplex(hit.DeltaG, s, hit.TopStart, bottomStart, hit.Length, Diagram(a, b, s, hit.TopStart, hit.Length));
    }

    private static string Diagram(string a, string b, int s, int topStart, int length)
    {
        int lb = b.Length;
        // Top base i at column padA + i; bottom (written 3'->5') base j at column padB + (lb - 1 - j).
        int shift = lb - 1 - s;
        int padA = Math.Max(0, shift), padB = Math.Max(0, -shift);

        var top = new StringBuilder("5' ").Append(' ', padA).Append(a.ToUpperInvariant()).Append(" 3'");
        var bottom = new StringBuilder("3' ").Append(' ', padB);
        for (int j = lb - 1; j >= 0; j--) bottom.Append(char.ToUpperInvariant(b[j]));
        bottom.Append(" 5'");
        var middle = new StringBuilder("   ").Append(' ', padA + topStart).Append('|', length);
        return string.Join('\n', top.ToString(), middle.ToString(), bottom.ToString());
    }
}
