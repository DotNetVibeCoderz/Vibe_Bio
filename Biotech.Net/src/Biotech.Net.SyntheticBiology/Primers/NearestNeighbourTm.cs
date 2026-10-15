namespace Biotech.Net.SyntheticBiology.Primers;

/// <summary>
/// Standard-state enthalpy and entropy of duplex formation, as produced by a nearest-neighbour sum.
/// </summary>
/// <remarks>
/// Values follow the usual sign convention for hybridisation (single strands → duplex): both ΔH° and ΔS° are
/// negative for a stable duplex. The entropy is the value at the 1 M Na⁺ reference state of the
/// SantaLucia (1998) parameter set; salt corrections are applied later by
/// <see cref="NearestNeighbourTm.TmFromThermodynamics"/>.
/// </remarks>
/// <param name="EnthalpyKcalPerMol">ΔH° in kcal/mol.</param>
/// <param name="EntropyCalPerMolK">ΔS° in cal/(mol·K).</param>
public readonly record struct DuplexThermodynamics(double EnthalpyKcalPerMol, double EntropyCalPerMolK)
{
    /// <summary>
    /// ΔG°37 = ΔH° − T·ΔS° at T = 310.15 K, in kcal/mol (entropy converted from cal to kcal).
    /// More negative means a more stable duplex.
    /// </summary>
    public double FreeEnergy37 => FreeEnergyAt(37.0);

    /// <summary>
    /// ΔG° at <paramref name="celsius"/> in kcal/mol, assuming ΔH° and ΔS° are temperature independent
    /// (ΔCp = 0), the standard two-state approximation used by SantaLucia (1998).
    /// </summary>
    /// <param name="celsius">Temperature in degrees Celsius.</param>
    public double FreeEnergyAt(double celsius)
        => EnthalpyKcalPerMol - (celsius + 273.15) * EntropyCalPerMolK / 1000.0;
}

/// <summary>
/// Cation and strand concentrations that set the melting temperature of an oligonucleotide duplex.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="MonovalentMolar"/> is the total monovalent cation concentration, conventionally
/// [Na⁺] + [K⁺] + [Tris⁺] where Tris counts as half its formal concentration because roughly half of it is
/// protonated at PCR pH (von Ahsen et al. 2001, Clin. Chem. 47:1956).
/// </para>
/// <para>
/// Deoxynucleoside triphosphates bind Mg²⁺ approximately 1:1, so the Mg²⁺ available to stabilise the duplex is
/// approximately [Mg²⁺] − [dNTP] (<see cref="FreeMagnesiumMolar"/>); this is the simple stoichiometric
/// correction recommended by von Ahsen et al. (2001) and used by Primer3's <c>dntp_conc</c> handling.
/// </para>
/// <para>
/// <see cref="PrimerMolar"/> is the total single-strand concentration C<sub>T</sub>. For a primer in large excess
/// over its template, the effective concentration in the Tm equation is C<sub>T</sub>/4 (two non-self-complementary
/// strands at equal concentration); see <see cref="NearestNeighbourTm.Tm"/>.
/// </para>
/// </remarks>
/// <param name="MonovalentMolar">Monovalent cations (Na⁺ + K⁺ + Tris/2), mol/L.</param>
/// <param name="MagnesiumMolar">Total Mg²⁺, mol/L.</param>
/// <param name="DntpMolar">Total dNTP, mol/L (sum of the four nucleotides).</param>
/// <param name="PrimerMolar">Oligonucleotide strand concentration, mol/L.</param>
public sealed record SaltConditions(
    double MonovalentMolar = 0.050,
    double MagnesiumMolar = 0.0,
    double DntpMolar = 0.0,
    double PrimerMolar = 0.25e-6)
{
    /// <summary>50 mM Na⁺, no Mg²⁺, 0.25 µM oligo — Primer3's historical defaults and a common benchmark condition.</summary>
    public static SaltConditions Default { get; } = new();

    /// <summary>A typical Taq PCR buffer: 50 mM KCl, 1.5 mM MgCl₂, 0.2 mM total dNTP, 0.25 µM of each primer.</summary>
    public static SaltConditions StandardPcr { get; } = new(0.050, 0.0015, 0.0002, 0.25e-6);

    /// <summary>
    /// A typical hydrolysis-probe qPCR master mix: 50 mM monovalent, 3 mM Mg²⁺, 0.8 mM total dNTP,
    /// 0.2 µM of each primer.
    /// </summary>
    public static SaltConditions Qpcr { get; } = new(0.050, 0.003, 0.0008, 0.2e-6);

    /// <summary>Mg²⁺ not chelated by dNTPs: max(0, [Mg²⁺] − [dNTP]), mol/L.</summary>
    public double FreeMagnesiumMolar => Math.Max(0.0, MagnesiumMolar - DntpMolar);
}

/// <summary>
/// Nearest-neighbour melting temperature of short DNA duplexes, with the SantaLucia (1998) unified parameters,
/// the SantaLucia monovalent salt correction and the Owczarzy et al. (2008) Mg²⁺ correction.
/// </summary>
/// <remarks>
/// <para>
/// <b>Model.</b> The duplex ΔH° and ΔS° are the sums of the ten unique Watson–Crick nearest-neighbour stacks of
/// SantaLucia J. Jr. (1998) "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
/// thermodynamics", PNAS 95:1460–1465 (Table 1), plus the two terminal initiation terms of that table
/// (terminal G·C: ΔH° +0.1 kcal/mol, ΔS° −2.8 e.u.; terminal A·T: ΔH° +2.3, ΔS° +4.1) and, for a
/// self-complementary duplex, the symmetry correction ΔS° −1.4 e.u.
/// </para>
/// <para>
/// <b>Tm.</b> For a two-state bimolecular transition, Tm = ΔH°·1000 / (ΔS° + R·ln(C<sub>T</sub>/x)) − 273.15, with
/// R = 1.987 cal/(mol·K), x = 4 for two non-self-complementary strands at equal concentration and x = 1 for a
/// self-complementary strand (SantaLucia 1998, eq. 3).
/// </para>
/// <para>
/// <b>Salt.</b> Without free Mg²⁺ the entropy is corrected as ΔS° + 0.368·(N − 1)·ln[Na⁺] (SantaLucia 1998, eq. 5).
/// With free Mg²⁺ the Owczarzy et al. (2008) Biochemistry 47:5336–5353 scheme is used, see
/// <see cref="TmFromThermodynamics"/>.
/// </para>
/// <para>
/// Results agree with Biopython's <c>Bio.SeqUtils.MeltingTemp.Tm_NN(table=DNA_NN3, saltcorr=5)</c> without Mg²⁺
/// and with <c>saltcorr=7</c> with Mg²⁺ (except for the dNTP chelation model, which here is the simple
/// stoichiometric subtraction). Expect a ±1–2 °C spread against other tools that use different salt models.
/// </para>
/// </remarks>
public static class NearestNeighbourTm
{
    /// <summary>Gas constant R in cal/(mol·K).</summary>
    public const double GasConstant = 1.987;

    // SantaLucia 1998 Table 1, indexed by 4*code(first) + code(second) with A=0, C=1, G=2, T=3.
    private static readonly double[] StackEnthalpy = new double[16];
    private static readonly double[] StackEntropy = new double[16];
    private static readonly double[] StackFreeEnergy37 = new double[16];

    private const double InitGcEnthalpy = 0.1, InitGcEntropy = -2.8;
    private const double InitAtEnthalpy = 2.3, InitAtEntropy = 4.1;
    private const double SymmetryEntropy = -1.4;

    static NearestNeighbourTm()
    {
        // The ten unique Watson-Crick stacks (5'->3' / 3'->5').
        Set("AA", -7.9, -22.2);
        Set("AT", -7.2, -20.4);
        Set("TA", -7.2, -21.3);
        Set("CA", -8.5, -22.7);
        Set("GT", -8.4, -22.4);
        Set("CT", -7.8, -21.0);
        Set("GA", -8.2, -22.2);
        Set("CG", -10.6, -27.2);
        Set("GC", -9.8, -24.4);
        Set("GG", -8.0, -19.9);

        static void Set(string stack, double dh, double ds)
        {
            int a = CodeOf(stack[0]), b = CodeOf(stack[1]);
            Store(4 * a + b, dh, ds);
            // The same physical stack read on the other strand: reverse complement of the dinucleotide.
            Store(4 * (3 - b) + (3 - a), dh, ds);
        }

        static void Store(int index, double dh, double ds)
        {
            StackEnthalpy[index] = dh;
            StackEntropy[index] = ds;
            StackFreeEnergy37[index] = dh - 310.15 * ds / 1000.0;
        }
    }

    /// <summary>
    /// Nearest-neighbour ΔH° and ΔS° (1 M Na⁺ reference) for the perfectly matched duplex formed by
    /// <paramref name="sequence"/> and its complement.
    /// </summary>
    /// <remarks>
    /// Sums the N − 1 stacks, adds one initiation term per end chosen by that end's base pair (G·C or A·T), and
    /// subtracts 1.4 e.u. when the sequence is self-complementary (SantaLucia 1998). Case-insensitive; U is read
    /// as T.
    /// </remarks>
    /// <param name="sequence">The oligonucleotide, 5′→3′, at least 2 nt of A/C/G/T/U.</param>
    /// <exception cref="ArgumentException">Shorter than 2 nt or containing a symbol other than A, C, G, T, U.</exception>
    public static DuplexThermodynamics Duplex(ReadOnlySpan<char> sequence)
    {
        if (sequence.Length < 2)
            throw new ArgumentException("A nearest-neighbour duplex needs at least 2 nucleotides.", nameof(sequence));

        double dh = 0, ds = 0;
        int previous = RequireCode(sequence[0]);
        for (int i = 1; i < sequence.Length; i++)
        {
            int current = RequireCode(sequence[i]);
            int index = 4 * previous + current;
            dh += StackEnthalpy[index];
            ds += StackEntropy[index];
            previous = current;
        }

        AddInitiation(RequireCode(sequence[0]), ref dh, ref ds);
        AddInitiation(previous, ref dh, ref ds);
        if (IsSelfComplementary(sequence)) ds += SymmetryEntropy;
        return new DuplexThermodynamics(dh, ds);
    }

    /// <inheritdoc cref="Duplex(ReadOnlySpan{char})"/>
    public static DuplexThermodynamics Duplex(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return Duplex(sequence.AsSpan());
    }

    /// <summary>
    /// Melting temperature (°C) of <paramref name="sequence"/> annealed to its perfect complement, under
    /// <paramref name="conditions"/> (default <see cref="SaltConditions.Default"/>: 50 mM Na⁺, 0.25 µM).
    /// </summary>
    /// <remarks>
    /// The primer is assumed to be in large excess or at equal concentration with its target, the convention of
    /// Primer3, OligoAnalyzer and Biopython. See the class remarks for the model.
    /// </remarks>
    /// <param name="sequence">The oligonucleotide, 5′→3′.</param>
    /// <param name="conditions">Salt and strand concentrations.</param>
    /// <exception cref="ArgumentException">Invalid sequence, or no cations at all.</exception>
    public static double Tm(string sequence, SaltConditions? conditions = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var d = Duplex(sequence.AsSpan());
        int gc = 0;
        foreach (char c in sequence)
        {
            char u = char.ToUpperInvariant(c);
            if (u == 'G' || u == 'C') gc++;
        }
        return TmFromThermodynamics(d, conditions ?? SaltConditions.Default, sequence.Length,
            (double)gc / sequence.Length, IsSelfComplementary(sequence));
    }

    /// <summary>
    /// Converts nearest-neighbour ΔH°/ΔS° (1 M Na⁺ reference) into a melting temperature (°C) under the given salt.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Strand term: R·ln(C<sub>T</sub>/x) with x = 1 when <paramref name="selfComplementary"/>, otherwise x = 4.
    /// </para>
    /// <para>
    /// Cation correction (Owczarzy R. et al. 2008, "Predicting stability of DNA duplexes in solutions containing
    /// magnesium and monovalent cations", Biochemistry 47:5336–5353; the same decision tree as Biopython's
    /// <c>MeltingTemp.salt_correction(method=7)</c>):
    /// </para>
    /// <list type="number">
    /// <item><description>No free Mg²⁺ (<see cref="SaltConditions.FreeMagnesiumMolar"/> = 0): SantaLucia (1998)
    /// entropy correction ΔS° + 0.368·(N − 1)·ln[Mon⁺].</description></item>
    /// <item><description>Free Mg²⁺ and monovalent ions with ratio R = √[Mg²⁺]/[Mon⁺] &lt; 0.22: monovalent ions
    /// dominate, so the SantaLucia monovalent correction alone is used.</description></item>
    /// <item><description>0.22 ≤ R &lt; 6.0: Owczarzy eq. 16 with the monovalent-dependent coefficients
    /// a′ = a·(0.843 − 0.352·√[Mon⁺]·ln[Mon⁺]),
    /// d′ = d·(1.279 − 4.03×10⁻³·ln[Mon⁺] − 8.03×10⁻³·ln²[Mon⁺]),
    /// g′ = g·(0.486 − 0.258·ln[Mon⁺] + 5.25×10⁻³·ln³[Mon⁺]) (Owczarzy 2008 eqs. 18–20).</description></item>
    /// <item><description>R ≥ 6.0, or no monovalent ions: Mg²⁺ dominates; Owczarzy eq. 16 with the original
    /// coefficients.</description></item>
    /// </list>
    /// <para>
    /// Owczarzy eq. 16, with Tm in kelvin and the Mg²⁺ term relative to the 1 M Na⁺ reference Tm (the reference
    /// state of the nearest-neighbour parameters, i.e. no entropy salt correction):
    /// 1/Tm(Mg²⁺) = 1/Tm(1 M Na⁺) + a + b·ln[Mg²⁺] + f<sub>GC</sub>·(c + d·ln[Mg²⁺])
    /// + (e + f·ln[Mg²⁺] + g·ln²[Mg²⁺]) / (2(N − 1)),
    /// with a = 3.92×10⁻⁵, b = −9.11×10⁻⁶, c = 6.26×10⁻⁵, d = 1.42×10⁻⁵, e = −4.82×10⁻⁴, f = 5.25×10⁻⁴,
    /// g = 8.31×10⁻⁵ (signs folded into the constants, as in Biopython; the paper writes −b and −e with positive
    /// magnitudes).
    /// </para>
    /// </remarks>
    /// <param name="d">Duplex ΔH°/ΔS° at the 1 M Na⁺ reference state.</param>
    /// <param name="conditions">Salt and strand concentrations.</param>
    /// <param name="length">Duplex length N in base pairs.</param>
    /// <param name="gcFraction">GC fraction (0–1), used by the Mg²⁺ correction.</param>
    /// <param name="selfComplementary">Whether the duplex is formed by one self-complementary strand.</param>
    /// <exception cref="ArgumentException">Length below 2, or neither monovalent nor free Mg²⁺ ions present.</exception>
    public static double TmFromThermodynamics(DuplexThermodynamics d, SaltConditions conditions, int length, double gcFraction, bool selfComplementary)
    {
        ArgumentNullException.ThrowIfNull(conditions);
        if (length < 2) throw new ArgumentException("Length must be at least 2.", nameof(length));
        if (conditions.PrimerMolar <= 0) throw new ArgumentException("Primer concentration must be positive.", nameof(conditions));

        double mon = conditions.MonovalentMolar;
        double mg = conditions.FreeMagnesiumMolar;
        if (mon <= 0 && mg <= 0)
            throw new ArgumentException("At least one of monovalent or free Mg2+ concentration must be positive.", nameof(conditions));

        double x = selfComplementary ? 1.0 : 4.0;
        double strandTerm = GasConstant * Math.Log(conditions.PrimerMolar / x);
        double dh = d.EnthalpyKcalPerMol * 1000.0;

        bool useMonovalentOnly = mg <= 0 || (mon > 0 && Math.Sqrt(mg) / mon < 0.22);
        if (useMonovalentOnly)
        {
            double ds = d.EntropyCalPerMolK + 0.368 * (length - 1) * Math.Log(mon);
            return dh / (ds + strandTerm) - 273.15;
        }

        double tm1M = dh / (d.EntropyCalPerMolK + strandTerm); // kelvin, 1 M Na+ reference
        double a = 3.92e-5, b = -9.11e-6, c = 6.26e-5, dd = 1.42e-5, e = -4.82e-4, f = 5.25e-4, g = 8.31e-5;
        if (mon > 0 && Math.Sqrt(mg) / mon < 6.0)
        {
            double lnMon = Math.Log(mon);
            a *= 0.843 - 0.352 * Math.Sqrt(mon) * lnMon;
            dd *= 1.279 - 4.03e-3 * lnMon - 8.03e-3 * lnMon * lnMon;
            g *= 0.486 - 0.258 * lnMon + 5.25e-3 * lnMon * lnMon * lnMon;
        }

        double lnMg = Math.Log(mg);
        double inverse = 1.0 / tm1M + a + b * lnMg + gcFraction * (c + dd * lnMg)
                         + (e + f * lnMg + g * lnMg * lnMg) / (2.0 * (length - 1));
        return 1.0 / inverse - 273.15;
    }

    /// <summary>
    /// The Wallace rule Tm = 2·(A + T) + 4·(G + C) °C (Suggs et al. 1981; Thein &amp; Wallace 1986), calibrated
    /// for 14–20-mers hybridised to membrane-bound DNA in ~0.9 M Na⁺. Only a rough guide, best below ~14 nt.
    /// </summary>
    /// <param name="sequence">The oligonucleotide (A/C/G/T/U, case-insensitive).</param>
    /// <exception cref="ArgumentException">Contains a symbol other than A, C, G, T, U, or is empty.</exception>
    public static double WallaceTm(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0) throw new ArgumentException("Sequence is empty.", nameof(sequence));
        int at = 0, gc = 0;
        foreach (char ch in sequence)
        {
            int code = RequireCode(ch);
            if (code == 1 || code == 2) gc++; else at++;
        }
        return 2.0 * at + 4.0 * gc;
    }

    /// <summary>
    /// The salt-adjusted GC-content ("basic") formula
    /// Tm = 81.5 + 16.6·log₁₀[Na⁺] + 0.41·(%GC) − 600/N (Howley et al. 1979, J. Biol. Chem. 254:4876;
    /// as used by Primer3's legacy mode and Biopython's <c>Tm_GC</c>). At 50 mM Na⁺ it is close to the
    /// popular 64.9 + 41·(G + C − 16.4)/N form for long oligos, but unlike that formula it accounts for salt.
    /// </summary>
    /// <remarks>Ignores sequence order; use <see cref="Tm"/> for primer design.</remarks>
    /// <param name="sequence">The oligonucleotide (A/C/G/T/U, case-insensitive).</param>
    /// <param name="monovalentMolar">Monovalent cation concentration in mol/L.</param>
    /// <exception cref="ArgumentException">Invalid sequence or non-positive salt.</exception>
    public static double GcContentTm(string sequence, double monovalentMolar = 0.050)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0) throw new ArgumentException("Sequence is empty.", nameof(sequence));
        if (monovalentMolar <= 0) throw new ArgumentException("Salt concentration must be positive.", nameof(monovalentMolar));
        int gc = 0;
        foreach (char ch in sequence)
        {
            int code = RequireCode(ch);
            if (code == 1 || code == 2) gc++;
        }
        double percentGc = 100.0 * gc / sequence.Length;
        return 81.5 + 16.6 * Math.Log10(monovalentMolar) + 0.41 * percentGc - 600.0 / sequence.Length;
    }

    /// <summary>
    /// Whether the sequence equals its own reverse complement (a palindrome in the molecular-biology sense,
    /// e.g. GAATTC). Such a strand forms a duplex with itself, which changes the Tm equation (x = 1) and adds
    /// the symmetry entropy penalty. Odd-length sequences can never be self-complementary.
    /// </summary>
    /// <param name="sequence">The sequence (case-insensitive, U read as T). Ambiguous symbols give false.</param>
    public static bool IsSelfComplementary(ReadOnlySpan<char> sequence)
    {
        int n = sequence.Length;
        if (n == 0 || (n & 1) == 1) return false;
        for (int i = 0; i < n / 2; i++)
        {
            int a = CodeOf(sequence[i]), b = CodeOf(sequence[n - 1 - i]);
            if (a < 0 || b < 0 || a + b != 3) return false;
        }
        return true;
    }

    // ---- internal helpers shared with OligoThermodynamics ----

    /// <summary>2-bit code of a character (A=0, C=1, G=2, T/U=3), or -1.</summary>
    internal static int CodeOf(char c) => c switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' or 'U' or 'u' => 3,
        _ => -1,
    };

    internal static int RequireCode(char c)
    {
        int code = CodeOf(c);
        if (code < 0) throw new ArgumentException($"Symbol '{c}' is not A, C, G, T or U.", "sequence");
        return code;
    }

    /// <summary>Parses to 2-bit codes, throwing on any non-ACGTU symbol.</summary>
    internal static int[] Codes(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var codes = new int[sequence.Length];
        for (int i = 0; i < codes.Length; i++) codes[i] = RequireCode(sequence[i]);
        return codes;
    }

    /// <summary>ΔG°37 (kcal/mol) of the Watson–Crick stack 5′-ab-3′ / 3′-a′b′-5′.</summary>
    internal static double StackDeltaG37(int a, int b) => StackFreeEnergy37[4 * a + b];

    /// <summary>ΔG°37 (kcal/mol) of one duplex-end initiation term for the base pair whose top-strand base is <paramref name="code"/>.</summary>
    internal static double InitiationDeltaG37(int code)
        => code is 1 or 2
            ? InitGcEnthalpy - 310.15 * InitGcEntropy / 1000.0
            : InitAtEnthalpy - 310.15 * InitAtEntropy / 1000.0;

    private static void AddInitiation(int code, ref double dh, ref double ds)
    {
        if (code is 1 or 2) { dh += InitGcEnthalpy; ds += InitGcEntropy; }
        else { dh += InitAtEnthalpy; ds += InitAtEntropy; }
    }
}
