using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Cloning;

/// <summary>
/// A type II restriction endonuclease: a recognition site plus the two phosphodiester bonds it hydrolyses.
/// </summary>
/// <remarks>
/// <para>
/// <b>Cut coordinate convention.</b> REBASE writes a type IIP enzyme as <c>G^AATTC</c> and a type IIS enzyme
/// as <c>GGTCTC(1/5)</c>. This record stores both in the same form:
/// </para>
/// <list type="bullet">
/// <item><description>
/// For a type IIP (and any enzyme with <see cref="IsTypeIis"/> false), <see cref="CutTop"/> and
/// <see cref="CutBottom"/> are offsets <b>from the first base of the recognition site</b> to the cut, counted
/// as "number of bases 5′ of the nick". EcoRI <c>G^AATTC</c> therefore has <c>CutTop = 1</c> (the top strand
/// is nicked after the G) and <c>CutBottom = 5</c> (the bottom strand is nicked after GAATT when read on the
/// top-strand coordinate axis), leaving the 4 nt 5′ overhang <c>AATT</c>.
/// </description></item>
/// <item><description>
/// For a type IIS (<see cref="IsTypeIis"/> true) the cut lies outside the site, and REBASE measures it from
/// the <b>last base of the recognition site</b>. BsaI <c>GGTCTC(1/5)</c> has <c>CutTop = 1</c> and
/// <c>CutBottom = 5</c>, meaning the top strand is cut 1 nt past the end of the 6 bp site and the bottom
/// strand 5 nt past it, leaving a 4 nt 5′ overhang. The properties
/// <see cref="TopCutOffsetFromSiteEnd"/> and <see cref="BottomCutOffsetFromSiteEnd"/> expose that reading, and
/// <see cref="TopCutOffsetFromSiteStart"/>/<see cref="BottomCutOffsetFromSiteStart"/> give the equivalent
/// site-start offsets used internally by <see cref="Find(string, bool)"/>. Both pairs are always consistent:
/// for a type IIS, start offset = <see cref="SiteLength"/> + end offset.
/// </description></item>
/// </list>
/// <para>
/// All positions produced by <see cref="Find(string, bool)"/> are 0-based indices into the top strand of the
/// searched molecule, and a "cut position" <c>p</c> means the bond between base <c>p-1</c> and base <c>p</c>.
/// </para>
/// <para>
/// Recognition sites may contain IUPAC ambiguity symbols (for example BstXI <c>CCANNNNNNTGG</c>); matching
/// goes through <see cref="DnaUtility.MatchesIupac(byte, byte)"/>.
/// </para>
/// <para>
/// Catalogue data follow the published REBASE specificities (Roberts RJ et al., "REBASE — a database for DNA
/// restriction and modification", Nucleic Acids Res; rebase.neb.com). Only a small, hand-checked working set
/// is included here.
/// </para>
/// </remarks>
/// <param name="Name">Enzyme name as used by REBASE, e.g. <c>EcoRI</c>.</param>
/// <param name="RecognitionSite">Recognition sequence, 5′→3′ on the top strand, IUPAC symbols allowed.</param>
/// <param name="CutTop">Top-strand cut offset (see the convention above).</param>
/// <param name="CutBottom">Bottom-strand cut offset (see the convention above).</param>
/// <param name="IsTypeIis">
/// True when the enzyme cuts outside its recognition site, so the offsets are measured from the site end.
/// </param>
public sealed record RestrictionEnzyme(
    string Name,
    string RecognitionSite,
    int CutTop,
    int CutBottom,
    bool IsTypeIis = false)
{
    /// <summary>Length of the recognition site in base pairs.</summary>
    public int SiteLength => RecognitionSite.Length;

    /// <summary>Top-strand cut offset measured from the first base of the recognition site.</summary>
    public int TopCutOffsetFromSiteStart => IsTypeIis ? SiteLength + CutTop : CutTop;

    /// <summary>Bottom-strand cut offset measured from the first base of the recognition site.</summary>
    public int BottomCutOffsetFromSiteStart => IsTypeIis ? SiteLength + CutBottom : CutBottom;

    /// <summary>
    /// Top-strand cut offset measured from the last base of the recognition site; 0 means the cut falls
    /// immediately after the site. Negative for a type IIP enzyme, whose cut lies inside the site.
    /// </summary>
    public int TopCutOffsetFromSiteEnd => IsTypeIis ? CutTop : CutTop - SiteLength;

    /// <summary>Bottom-strand cut offset measured from the last base of the recognition site.</summary>
    public int BottomCutOffsetFromSiteEnd => IsTypeIis ? CutBottom : CutBottom - SiteLength;

    /// <summary>
    /// Signed overhang length: positive for a 5′ (recessed top strand) overhang, negative for a 3′ overhang,
    /// zero for a blunt cut. It is simply bottom cut minus top cut.
    /// </summary>
    public int OverhangLength => BottomCutOffsetFromSiteStart - TopCutOffsetFromSiteStart;

    /// <summary>Whether the enzyme leaves blunt ends.</summary>
    public bool IsBlunt => OverhangLength == 0;

    /// <summary>Human-readable end chemistry: <c>5'</c>, <c>3'</c> or <c>blunt</c>.</summary>
    public string OverhangKind => OverhangLength > 0 ? "5'" : OverhangLength < 0 ? "3'" : "blunt";

    /// <summary>
    /// Whether the recognition site equals its own reverse complement, so that a single search of the top
    /// strand finds every site. Type IIS enzymes are not palindromic and must be searched on both strands.
    /// </summary>
    public bool IsPalindromic
    {
        get
        {
            string site = RecognitionSite;
            int n = site.Length;
            for (int i = 0; i < n; i++)
            {
                char a = site[i];
                char b = (char)DnaUtility.Complement((byte)site[n - 1 - i]);
                if (char.ToUpperInvariant(a) != char.ToUpperInvariant(b)) return false;
            }
            return true;
        }
    }

    /// <summary>REBASE-style rendering, e.g. <c>EcoRI G^AATTC</c> or <c>BsaI GGTCTC(1/5)</c>.</summary>
    public override string ToString() => IsTypeIis
        ? $"{Name} {RecognitionSite}({CutTop}/{CutBottom})"
        : $"{Name} {RecognitionSite[..Math.Clamp(CutTop, 0, SiteLength)]}^{RecognitionSite[Math.Clamp(CutTop, 0, SiteLength)..]}";

    /// <summary>
    /// Finds every occurrence of the recognition site in <paramref name="dna"/>.
    /// </summary>
    /// <param name="dna">Top-strand sequence; case-insensitive, U is read as T.</param>
    /// <param name="circular">
    /// When true the molecule is treated as a circle, so sites that run across the origin are reported as
    /// well. Their <see cref="RestrictionSite.Position"/> is the 0-based start on the linearised top strand
    /// and may place part of the site past the end of the string.
    /// </param>
    /// <returns>Sites ordered by top-strand cut position.</returns>
    /// <remarks>
    /// A palindromic site is reported once, on the forward strand. A non-palindromic site (every type IIS,
    /// plus enzymes such as BstXI with an asymmetric ambiguity pattern) is searched on both strands; a reverse
    /// hit has <see cref="Strand.Reverse"/> and its cut positions lie on the other side of the site, which is
    /// exactly why a type IIS enzyme can be oriented to cut away from itself.
    /// </remarks>
    public IReadOnlyList<RestrictionSite> Find(string dna, bool circular = false)
    {
        ArgumentNullException.ThrowIfNull(dna);
        string top = DnaUtility.NormalizeDna(dna);
        int n = top.Length;
        var results = new List<RestrictionSite>();
        if (n == 0 || SiteLength == 0) return results;

        byte[] forwardPattern = DnaUtility.Ascii(RecognitionSite.ToUpperInvariant());
        byte[] reversePattern = DnaUtility.Ascii(DnaUtility.ReverseComplement(RecognitionSite.ToUpperInvariant()));
        bool palindromic = IsPalindromic;

        // Search window: for a circular molecule a site may start up to SiteLength-1 bases before the origin.
        int lastStart = circular ? n - 1 : n - SiteLength;
        var window = new byte[SiteLength];

        for (int start = 0; start <= lastStart; start++)
        {
            if (!circular && start + SiteLength > n) break;
            for (int i = 0; i < SiteLength; i++)
            {
                int idx = start + i;
                window[i] = (byte)top[circular ? idx % n : idx];
            }

            if (DnaUtility.MatchesIupac(window, forwardPattern))
                results.Add(BuildSite(top, start, Strand.Forward, circular));

            if (!palindromic && DnaUtility.MatchesIupac(window, reversePattern))
                results.Add(BuildSite(top, start, Strand.Reverse, circular));
        }

        results.Sort(static (a, b) => a.TopCutPosition.CompareTo(b.TopCutPosition));
        return results;
    }

    /// <summary>Number of occurrences of the site in a sequence.</summary>
    /// <param name="dna">Top-strand sequence.</param>
    /// <param name="circular">Whether to include origin-spanning sites.</param>
    public int CountSites(string dna, bool circular = false) => Find(dna, circular).Count;

    private RestrictionSite BuildSite(string top, int start, Strand strand, bool circular)
    {
        int n = top.Length;
        int topCut, bottomCut;
        if (strand == Strand.Forward)
        {
            topCut = start + TopCutOffsetFromSiteStart;
            bottomCut = start + BottomCutOffsetFromSiteStart;
        }
        else
        {
            // On the reverse strand the enzyme reads the site right-to-left, so the offsets are mirrored about
            // the site: a cut d bases into the site from its 5' end (bottom-strand reading) lands d bases to
            // the left of the site's right-hand edge on the top-strand axis.
            int siteEnd = start + SiteLength;
            topCut = siteEnd - BottomCutOffsetFromSiteStart;
            bottomCut = siteEnd - TopCutOffsetFromSiteStart;
        }

        string overhang = ExtractOverhang(top, topCut, bottomCut, circular);
        if (circular && n > 0)
        {
            topCut = Mod(topCut, n);
            bottomCut = Mod(bottomCut, n);
        }
        return new RestrictionSite(this, start, strand, topCut, bottomCut, overhang);
    }

    private static string ExtractOverhang(string top, int topCut, int bottomCut, bool circular)
    {
        int lo = Math.Min(topCut, bottomCut);
        int hi = Math.Max(topCut, bottomCut);
        if (hi == lo) return string.Empty;
        int n = top.Length;
        var chars = new char[hi - lo];
        for (int i = 0; i < chars.Length; i++)
        {
            int idx = lo + i;
            if (circular) idx = Mod(idx, n);
            if (idx < 0 || idx >= n) return string.Empty;
            chars[i] = top[idx];
        }
        return new string(chars);
    }

    private static int Mod(int value, int modulus)
    {
        int r = value % modulus;
        return r < 0 ? r + modulus : r;
    }

    /// <summary>
    /// The hand-checked built-in catalogue: common type IIP cloning enzymes plus the type IIS workhorses used
    /// for Golden Gate and MoClo assembly.
    /// </summary>
    /// <remarks>
    /// Specificities as published by REBASE. The type IIP entries are, in REBASE notation:
    /// EcoRI <c>G^AATTC</c>, BamHI <c>G^GATCC</c>, HindIII <c>A^AGCTT</c>, XhoI <c>C^TCGAG</c>,
    /// SalI <c>G^TCGAC</c>, NotI <c>GC^GGCCGC</c>, XbaI <c>T^CTAGA</c>, SpeI <c>A^CTAGT</c>,
    /// BglII <c>A^GATCT</c>, NcoI <c>C^CATGG</c>, NdeI <c>CA^TATG</c>, PstI <c>CTGCA^G</c>,
    /// SacI <c>GAGCT^C</c>, KpnI <c>GGTAC^C</c>, SmaI <c>CCC^GGG</c>, EcoRV <c>GAT^ATC</c>,
    /// BstXI <c>CCANNNNN^NTGG</c>. The type IIS entries are BsaI <c>GGTCTC(1/5)</c>,
    /// BsmBI/Esp3I <c>CGTCTC(1/5)</c>, BbsI <c>GAAGAC(2/6)</c>, SapI <c>GCTCTTC(1/4)</c> and
    /// AarI <c>CACCTGC(4/8)</c>.
    /// </remarks>
    public static IReadOnlyList<RestrictionEnzyme> All { get; } =
    [
        // --- Type IIP, 4 nt 5' overhangs ---
        new("EcoRI", "GAATTC", 1, 5),
        new("BamHI", "GGATCC", 1, 5),
        new("HindIII", "AAGCTT", 1, 5),
        new("XhoI", "CTCGAG", 1, 5),
        new("SalI", "GTCGAC", 1, 5),
        new("XbaI", "TCTAGA", 1, 5),
        new("SpeI", "ACTAGT", 1, 5),
        new("BglII", "AGATCT", 1, 5),
        new("NcoI", "CCATGG", 1, 5),
        new("NotI", "GCGGCCGC", 2, 6),
        // NdeI leaves a 2 nt 5' overhang (TA).
        new("NdeI", "CATATG", 2, 4),

        // --- Type IIP, 4 nt 3' overhangs ---
        new("PstI", "CTGCAG", 5, 1),
        new("SacI", "GAGCTC", 5, 1),
        new("KpnI", "GGTACC", 5, 1),

        // --- Type IIP, blunt ---
        new("SmaI", "CCCGGG", 3, 3),
        new("EcoRV", "GATATC", 3, 3),

        // --- Interrupted / ambiguous site: BstXI CCANNNNN^NTGG, 4 nt 3' overhang ---
        new("BstXI", "CCANNNNNNTGG", 8, 4),

        // --- Type IIS: offsets measured from the end of the site ---
        new("BsaI", "GGTCTC", 1, 5, IsTypeIis: true),
        new("BsmBI", "CGTCTC", 1, 5, IsTypeIis: true),
        new("Esp3I", "CGTCTC", 1, 5, IsTypeIis: true),
        new("BbsI", "GAAGAC", 2, 6, IsTypeIis: true),
        new("SapI", "GCTCTTC", 1, 4, IsTypeIis: true),
        new("AarI", "CACCTGC", 4, 8, IsTypeIis: true),
    ];

    private static readonly Dictionary<string, RestrictionEnzyme> ByName = BuildIndex();

    private static Dictionary<string, RestrictionEnzyme> BuildIndex()
    {
        var map = new Dictionary<string, RestrictionEnzyme>(StringComparer.OrdinalIgnoreCase);
        foreach (var enzyme in All) map[enzyme.Name] = enzyme;
        return map;
    }

    /// <summary>Looks a catalogue enzyme up by name (case-insensitive).</summary>
    /// <exception cref="KeyNotFoundException">The name is not in the built-in catalogue.</exception>
    public static RestrictionEnzyme FromName(string name)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        return ByName.TryGetValue(name, out var enzyme)
            ? enzyme
            : throw new KeyNotFoundException($"No built-in restriction enzyme named '{name}'.");
    }

    /// <summary>Tries to look a catalogue enzyme up by name (case-insensitive).</summary>
    public static bool TryFromName(string name, out RestrictionEnzyme enzyme)
    {
        if (string.IsNullOrEmpty(name))
        {
            enzyme = null!;
            return false;
        }
        return ByName.TryGetValue(name, out enzyme!);
    }
}

/// <summary>
/// One occurrence of a recognition site on a molecule, with the two cut positions it implies.
/// </summary>
/// <param name="Enzyme">The enzyme that recognises this site.</param>
/// <param name="Position">0-based start of the recognition site on the top strand.</param>
/// <param name="Strand">
/// Which strand carries the site as written: <see cref="Biotech.Net.Core.Strand.Forward"/> when the top strand
/// reads as <see cref="RestrictionEnzyme.RecognitionSite"/>, otherwise <see cref="Biotech.Net.Core.Strand.Reverse"/>.
/// </param>
/// <param name="TopCutPosition">
/// 0-based top-strand cut: the bond between base <c>TopCutPosition-1</c> and <c>TopCutPosition</c>.
/// </param>
/// <param name="BottomCutPosition">0-based bottom-strand cut, on the same top-strand coordinate axis.</param>
/// <param name="Overhang">
/// The single-stranded sequence released between the two cuts, written 5′→3′ on the top strand; empty for a
/// blunt cutter.
/// </param>
public sealed record RestrictionSite(
    RestrictionEnzyme Enzyme,
    int Position,
    Strand Strand,
    int TopCutPosition,
    int BottomCutPosition,
    string Overhang)
{
    /// <summary>Signed overhang length at this site: positive 5′, negative 3′, zero blunt.</summary>
    public int OverhangLength => BottomCutPosition - TopCutPosition;

    /// <summary>Readable summary, e.g. <c>EcoRI@12 forward top=13 bottom=17 AATT</c>.</summary>
    public override string ToString()
        => $"{Enzyme.Name}@{Position} {Strand} top={TopCutPosition} bottom={BottomCutPosition} {Overhang}";
}
