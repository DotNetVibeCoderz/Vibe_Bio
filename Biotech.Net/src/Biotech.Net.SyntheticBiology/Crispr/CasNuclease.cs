namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>Which side of the protospacer the PAM sits on.</summary>
public enum PamSide
{
    /// <summary>PAM follows the protospacer, as for SpCas9 (<c>protospacer-NGG</c>).</summary>
    ThreePrime,

    /// <summary>PAM precedes the protospacer, as for Cas12a/Cpf1 (<c>TTTV-protospacer</c>).</summary>
    FivePrime,
}

/// <summary>
/// A CRISPR nuclease described as data: PAM pattern, protospacer length and where it cuts. Everything in the
/// guide-design pipeline is driven by one of these, so adding a new engineered variant needs no new code.
/// </summary>
/// <param name="Name">Display name, e.g. <c>SpCas9</c>.</param>
/// <param name="Pam">PAM as an IUPAC pattern, e.g. <c>NGG</c> or <c>NNGRRT</c>.</param>
/// <param name="PamSide">Whether the PAM is 3′ or 5′ of the protospacer.</param>
/// <param name="GuideLength">Protospacer (spacer) length in nucleotides.</param>
/// <param name="BluntCutOffset">
/// Distance of the cut from the PAM-proximal end of the protospacer, counted into the protospacer. SpCas9 cuts
/// 3 bp from the PAM, so this is 3. Cas12a cuts outside the protospacer on the PAM-distal side, so the value
/// is negative there (see <see cref="StaggerLength"/>).
/// </param>
/// <param name="StaggerLength">
/// 0 for a blunt cutter. Cas12a leaves a 5′ overhang; the value is the number of staggered nucleotides between
/// the two strand cuts.
/// </param>
/// <param name="Organism">Source organism or engineering lineage, for reporting.</param>
public sealed record CasNuclease(
    string Name,
    string Pam,
    PamSide PamSide,
    int GuideLength,
    int BluntCutOffset,
    int StaggerLength = 0,
    string Organism = "")
{
    private readonly byte[] _pamBytes = DnaUtility.Ascii(Pam);

    /// <summary>PAM length in nucleotides.</summary>
    public int PamLength => Pam.Length;

    /// <summary>True when the nuclease leaves staggered (sticky) ends rather than a blunt cut.</summary>
    public bool IsStaggered => StaggerLength != 0;

    /// <summary>Total footprint on the genome: protospacer plus PAM.</summary>
    public int SiteLength => GuideLength + PamLength;

    /// <summary>Whether <paramref name="site"/> (exactly <see cref="PamLength"/> bases) matches the PAM pattern.</summary>
    public bool MatchesPam(ReadOnlySpan<byte> site) => DnaUtility.MatchesIupac(site, _pamBytes);

    /// <summary>The PAM pattern as ASCII bytes, without copying.</summary>
    public ReadOnlySpan<byte> PamSpan => _pamBytes;

    /// <summary>
    /// The 0-based cut position on the forward strand for a protospacer that starts at
    /// <paramref name="protospacerStart"/> on strand <paramref name="strand"/>. The returned coordinate is the
    /// index of the first base <em>after</em> the break, i.e. the break lies between <c>result - 1</c> and
    /// <c>result</c>, the same convention as a half-open <see cref="Core.GenomicInterval"/> boundary.
    /// </summary>
    /// <remarks>
    /// <paramref name="protospacerStart"/> is always a forward-strand coordinate: it is the leftmost base of the
    /// protospacer as it appears on the plus strand, even for a guide that targets the minus strand.
    /// </remarks>
    public int CutSite(int protospacerStart, Core.Strand strand)
    {
        bool forward = strand != Core.Strand.Reverse;
        if (PamSide == PamSide.ThreePrime)
        {
            // PAM-proximal end is the right end on the plus strand, the left end on the minus strand.
            return forward
                ? protospacerStart + GuideLength - BluntCutOffset
                : protospacerStart + BluntCutOffset;
        }
        // 5′ PAM (Cas12a): the PAM-proximal end is the left end of the protospacer on the targeted strand.
        return forward
            ? protospacerStart + GuideLength + BluntCutOffset
            : protospacerStart - BluntCutOffset;
    }

    /// <summary>
    /// <em>Streptococcus pyogenes</em> Cas9: 20 nt spacer, 3′ NGG PAM, blunt cut 3 bp from the PAM.
    /// The workhorse nuclease and the one every published scoring model was trained on.
    /// </summary>
    public static readonly CasNuclease SpCas9 = new("SpCas9", "NGG", PamSide.ThreePrime, 20, 3, 0, "Streptococcus pyogenes");

    /// <summary>SpCas9 with the relaxed NAG PAM, which it cuts at a few percent of NGG efficiency.</summary>
    public static readonly CasNuclease SpCas9Nag = new("SpCas9 (NAG)", "NAG", PamSide.ThreePrime, 20, 3, 0, "Streptococcus pyogenes");

    /// <summary>VQR variant of SpCas9 (Kleinstiver 2015), NGA PAM.</summary>
    public static readonly CasNuclease SpCas9Vqr = new("SpCas9-VQR", "NGA", PamSide.ThreePrime, 20, 3, 0, "engineered SpCas9");

    /// <summary>VRER variant of SpCas9 (Kleinstiver 2015), NGCG PAM.</summary>
    public static readonly CasNuclease SpCas9Vrer = new("SpCas9-VRER", "NGCG", PamSide.ThreePrime, 20, 3, 0, "engineered SpCas9");

    /// <summary>xCas9 3.7 (Hu 2018): broadened NG PAM recognition.</summary>
    public static readonly CasNuclease XCas9 = new("xCas9 3.7", "NG", PamSide.ThreePrime, 20, 3, 0, "engineered SpCas9");

    /// <summary>SpG (Walton 2020): NGN PAM.</summary>
    public static readonly CasNuclease SpG = new("SpG", "NGN", PamSide.ThreePrime, 20, 3, 0, "engineered SpCas9");

    /// <summary>
    /// SpRY (Walton 2020): near-PAMless, NRN &gt; NYN. Modelled with the permissive NNN pattern; expect very many
    /// candidate sites and correspondingly poor specificity, so always rank by off-target score.
    /// </summary>
    public static readonly CasNuclease SpRy = new("SpRY", "NNN", PamSide.ThreePrime, 20, 3, 0, "engineered SpCas9");

    /// <summary><em>Staphylococcus aureus</em> Cas9: 21 nt spacer, 3′ NNGRRT PAM, small enough for single-AAV delivery.</summary>
    public static readonly CasNuclease SaCas9 = new("SaCas9", "NNGRRT", PamSide.ThreePrime, 21, 3, 0, "Staphylococcus aureus");

    /// <summary>KKH variant of SaCas9 (Kleinstiver 2015): relaxed NNNRRT PAM.</summary>
    public static readonly CasNuclease SaCas9Kkh = new("SaCas9-KKH", "NNNRRT", PamSide.ThreePrime, 21, 3, 0, "engineered SaCas9");

    /// <summary>
    /// Cas12a (Cpf1), <em>Acidaminococcus</em>: 23 nt spacer with a 5′ TTTV PAM. It cuts outside the protospacer,
    /// 18/23 nt from the PAM, leaving a 5 nt 5′ overhang.
    /// </summary>
    public static readonly CasNuclease Cas12a = new("AsCas12a", "TTTV", PamSide.FivePrime, 23, -5, 5, "Acidaminococcus sp.");

    /// <summary>Cas12a from <em>Lachnospiraceae bacterium</em>, same PAM and cutting geometry as AsCas12a.</summary>
    public static readonly CasNuclease LbCas12a = new("LbCas12a", "TTTV", PamSide.FivePrime, 23, -5, 5, "Lachnospiraceae bacterium");

    /// <summary><em>Neisseria meningitidis</em> Cas9: 24 nt spacer, long NNNNGATT PAM, very high specificity.</summary>
    public static readonly CasNuclease NmeCas9 = new("NmeCas9", "NNNNGATT", PamSide.ThreePrime, 24, 3, 0, "Neisseria meningitidis");

    /// <summary>Every built-in nuclease.</summary>
    public static IReadOnlyList<CasNuclease> All { get; } =
    [
        SpCas9, SpCas9Nag, SpCas9Vqr, SpCas9Vrer, XCas9, SpG, SpRy,
        SaCas9, SaCas9Kkh, Cas12a, LbCas12a, NmeCas9,
    ];

    /// <summary>Looks a built-in nuclease up by name (case-insensitive).</summary>
    public static CasNuclease FromName(string name)
        => All.FirstOrDefault(n => string.Equals(n.Name, name, StringComparison.OrdinalIgnoreCase))
           ?? throw new ArgumentOutOfRangeException(nameof(name), $"'{name}' is not a built-in nuclease.");

    /// <inheritdoc/>
    public override string ToString() => $"{Name} ({(PamSide == PamSide.FivePrime ? Pam + "-N" : "N-" + Pam)}{GuideLength})";
}
