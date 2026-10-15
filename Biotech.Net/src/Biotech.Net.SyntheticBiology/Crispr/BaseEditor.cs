using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// The chemistry of a base editor: which deaminase it carries and therefore which transition it installs.
/// </summary>
public enum BaseEditorKind
{
    /// <summary>
    /// Cytosine base editor (CBE). A cytidine deaminase converts C to U on the single-stranded
    /// protospacer strand exposed by the nicked Cas9; replication reads the U as T, giving C→T
    /// on the protospacer strand and G→A on the complementary strand.
    /// </summary>
    Cytosine,

    /// <summary>
    /// Adenine base editor (ABE). An evolved tRNA adenosine deaminase converts A to inosine on the
    /// protospacer strand, which is read as G, giving A→G on the protospacer strand and T→C on the
    /// complementary strand.
    /// </summary>
    Adenine,
}

/// <summary>
/// How a single predicted base change reads at the protein level.
/// </summary>
public enum EditConsequence
{
    /// <summary>The base change does not fall inside a complete codon of the supplied coding sequence.</summary>
    NonCoding,

    /// <summary>The codon changes but the amino acid does not.</summary>
    Silent,

    /// <summary>The codon changes to a different sense codon.</summary>
    Missense,

    /// <summary>A sense codon becomes a stop codon.</summary>
    Nonsense,

    /// <summary>A stop codon becomes a sense codon (read-through).</summary>
    StopLoss,

    /// <summary>The initiator codon of the supplied coding sequence is destroyed.</summary>
    StartLoss,
}

/// <summary>
/// A base editor: a catalytically impaired Cas (usually an nCas9 nickase, which has exactly the same
/// protospacer/PAM geometry as <see cref="CasNuclease.SpCas9"/>) fused to a deaminase that acts on a short
/// single-stranded window of the protospacer strand.
/// </summary>
/// <remarks>
/// <para>
/// <b>Position convention.</b> Every position in this type is a <em>1-based protospacer position counted from
/// the PAM-distal end</em>, the convention used throughout the base-editing literature. For a 20 nt
/// SpCas9 protospacer written 5′→3′ on the targeted strand, position 1 is the first character of
/// <see cref="Guide.Spacer"/> and position 20 is the base immediately 5′ of the PAM. A "position 4–8 window"
/// therefore means <c>spacer[3..8]</c> in zero-based string indexing.
/// </para>
/// <para>
/// <b>Strand convention.</b> The deaminase acts on the protospacer strand — the strand whose sequence
/// <see cref="Guide.Spacer"/> reproduces, i.e. the strand <em>displaced</em> by the guide RNA. A CBE therefore
/// writes C→T in the spacer sequence as written, and the same event reads as G→A on the complementary strand.
/// An ABE writes A→G in the spacer as written, and T→C on the complementary strand. When you map an edit back
/// onto a plus-strand reference and the guide is a minus-strand guide, you must complement the change; the
/// coding-consequence helpers on this type do that for you.
/// </para>
/// <para>
/// <b>Honesty note.</b> The window boundaries and per-position weights of the built-in editors below are
/// approximate, qualitative values read off the activity profiles in the primary papers. They are
/// editor-dependent, site-dependent and sequence-context dependent; real efficiencies vary several-fold
/// between loci. Nothing here is a calibrated efficiency model, and these numbers must not be used as one.
/// For calibrated predictions, use a trained model such as BE-Hive (Arbab <i>et al.</i> 2020, Cell 182:463).
/// </para>
/// </remarks>
/// <param name="Name">Display name of the editor, e.g. <c>"BE4max"</c>.</param>
/// <param name="Kind">Whether this is a cytosine or an adenine editor.</param>
/// <param name="Nuclease">
/// The Cas geometry the editor is built on. This supplies the PAM and the protospacer length; for the
/// common nCas9(D10A) editors it is simply <see cref="CasNuclease.SpCas9"/>.
/// </param>
/// <param name="EditingWindowStart">First (lowest) 1-based PAM-distal protospacer position of the activity window.</param>
/// <param name="EditingWindowEnd">Last (highest) 1-based PAM-distal protospacer position of the activity window, inclusive.</param>
/// <param name="ActivityWeights">
/// Optional relative activity per 1-based protospacer position. Values are relative only — the largest is
/// conventionally 1.0 — and are used solely to rank candidate edits, never as an efficiency. When null, every
/// position inside the window is weighted 1.0.
/// </param>
/// <param name="Reference">Primary literature citation for the window of this editor.</param>
public sealed record BaseEditor(
    string Name,
    BaseEditorKind Kind,
    CasNuclease Nuclease,
    int EditingWindowStart,
    int EditingWindowEnd,
    IReadOnlyDictionary<int, double>? ActivityWeights = null,
    string Reference = "")
{
    /// <summary>The base this editor converts, as read on the protospacer strand: <c>C</c> for a CBE, <c>A</c> for an ABE.</summary>
    public char TargetBase => Kind == BaseEditorKind.Cytosine ? 'C' : 'A';

    /// <summary>The base produced on the protospacer strand: <c>T</c> for a CBE, <c>G</c> for an ABE.</summary>
    public char ProductBase => Kind == BaseEditorKind.Cytosine ? 'T' : 'G';

    /// <summary>The edited base as read on the complementary strand before editing: <c>G</c> for a CBE, <c>T</c> for an ABE.</summary>
    public char ComplementaryTargetBase => Kind == BaseEditorKind.Cytosine ? 'G' : 'T';

    /// <summary>The edited base as read on the complementary strand after editing: <c>A</c> for a CBE, <c>C</c> for an ABE.</summary>
    public char ComplementaryProductBase => Kind == BaseEditorKind.Cytosine ? 'A' : 'C';

    /// <summary>Number of protospacer positions in the activity window, inclusive of both ends.</summary>
    public int WindowLength => EditingWindowEnd - EditingWindowStart + 1;

    /// <summary>
    /// The centre of the activity window as a (possibly fractional) 1-based protospacer position. Used when
    /// ranking candidate guides, which prefer the requested base to sit near the middle of the window.
    /// </summary>
    public double WindowCentre => (EditingWindowStart + EditingWindowEnd) / 2.0;

    /// <summary>True when <paramref name="protospacerPosition"/> (1-based, PAM-distal) lies inside the activity window.</summary>
    /// <param name="protospacerPosition">A 1-based protospacer position counted from the PAM-distal end.</param>
    public bool IsInWindow(int protospacerPosition)
        => protospacerPosition >= EditingWindowStart && protospacerPosition <= EditingWindowEnd;

    /// <summary>
    /// The relative activity weight at a 1-based protospacer position: 0 outside the window, otherwise the
    /// value from <see cref="ActivityWeights"/> or 1.0 when no weights were supplied.
    /// </summary>
    /// <param name="protospacerPosition">A 1-based protospacer position counted from the PAM-distal end.</param>
    public double WeightAt(int protospacerPosition)
    {
        if (!IsInWindow(protospacerPosition)) return 0.0;
        if (ActivityWeights is not null && ActivityWeights.TryGetValue(protospacerPosition, out double w)) return w;
        return 1.0;
    }

    /// <summary>
    /// BE3 (Komor <i>et al.</i> 2016, Nature 533:420): rAPOBEC1–nCas9(D10A)–UGI. The reported activity window
    /// is approximately protospacer positions 4–8, peaking around 5–7.
    /// </summary>
    public static readonly BaseEditor Be3 = new(
        "BE3", BaseEditorKind.Cytosine, CasNuclease.SpCas9, 4, 8,
        new Dictionary<int, double> { [4] = 0.6, [5] = 0.9, [6] = 1.0, [7] = 0.9, [8] = 0.5 },
        "Komor et al. 2016, Nature 533:420");

    /// <summary>
    /// BE4max (Koblan <i>et al.</i> 2018, Nat Biotechnol 36:843): codon- and NLS-optimised BE4. Same
    /// approximate window as BE3, positions 4–8, with higher absolute editing.
    /// </summary>
    public static readonly BaseEditor Be4Max = new(
        "BE4max", BaseEditorKind.Cytosine, CasNuclease.SpCas9, 4, 8,
        new Dictionary<int, double> { [4] = 0.7, [5] = 0.95, [6] = 1.0, [7] = 0.9, [8] = 0.6 },
        "Koblan et al. 2018, Nat Biotechnol 36:843");

    /// <summary>
    /// YE1-BE4 (Kim <i>et al.</i> 2017, Nat Biotechnol 35:371): an rAPOBEC1 W90Y/R126E variant whose window is
    /// narrowed to roughly positions 5–6, trading activity for precision when bystanders must be avoided.
    /// </summary>
    public static readonly BaseEditor Ye1 = new(
        "YE1-BE4", BaseEditorKind.Cytosine, CasNuclease.SpCas9, 5, 6,
        new Dictionary<int, double> { [5] = 1.0, [6] = 0.8 },
        "Kim et al. 2017, Nat Biotechnol 35:371");

    /// <summary>
    /// ABE7.10 (Gaudelli <i>et al.</i> 2017, Nature 551:464): TadA–TadA*–nCas9(D10A). Activity window
    /// approximately protospacer positions 4–7, peaking at 5–6.
    /// </summary>
    public static readonly BaseEditor Abe710 = new(
        "ABE7.10", BaseEditorKind.Adenine, CasNuclease.SpCas9, 4, 7,
        new Dictionary<int, double> { [4] = 0.7, [5] = 1.0, [6] = 0.9, [7] = 0.5 },
        "Gaudelli et al. 2017, Nature 551:464");

    /// <summary>
    /// ABE8e (Richter <i>et al.</i> 2020, Nat Biotechnol 38:883): a phage-assisted-evolved TadA-8e with much
    /// faster deamination and a correspondingly broader window, approximately positions 4–8.
    /// </summary>
    public static readonly BaseEditor Abe8E = new(
        "ABE8e", BaseEditorKind.Adenine, CasNuclease.SpCas9, 4, 8,
        new Dictionary<int, double> { [4] = 0.85, [5] = 1.0, [6] = 1.0, [7] = 0.85, [8] = 0.6 },
        "Richter et al. 2020, Nat Biotechnol 38:883");

    /// <summary>All built-in editors, in no particular order.</summary>
    public static IReadOnlyList<BaseEditor> All { get; } = [Be3, Be4Max, Ye1, Abe710, Abe8E];

    /// <summary>Looks a built-in editor up by name, case-insensitively.</summary>
    /// <param name="name">The editor name, e.g. <c>"BE4max"</c> or <c>"ABE8e"</c>.</param>
    /// <exception cref="ArgumentOutOfRangeException">No built-in editor carries that name.</exception>
    public static BaseEditor FromName(string name)
        => All.FirstOrDefault(e => string.Equals(e.Name, name, StringComparison.OrdinalIgnoreCase))
           ?? throw new ArgumentOutOfRangeException(nameof(name), name, "Unknown base editor.");

    /// <inheritdoc/>
    public override string ToString()
        => $"{Name} ({TargetBase}→{ProductBase}, positions {EditingWindowStart}-{EditingWindowEnd})";
}

/// <summary>
/// One editable base inside a base editor's activity window.
/// </summary>
/// <param name="Position">1-based protospacer position, counted from the PAM-distal end.</param>
/// <param name="From">The base before editing, as read on the protospacer strand.</param>
/// <param name="To">The base after editing, as read on the protospacer strand.</param>
/// <param name="RelativeWeight">Relative activity weight at this position; see <see cref="BaseEditor.WeightAt"/>.</param>
/// <param name="IsIntended">True for the single base the caller wanted edited; false for a bystander.</param>
public sealed record BaseEditSite(int Position, char From, char To, double RelativeWeight, bool IsIntended)
{
    /// <inheritdoc/>
    public override string ToString() => $"{From}{Position}{To}{(IsIntended ? "" : " (bystander)")}";
}

/// <summary>
/// What a base editor is predicted to do to one protospacer.
/// </summary>
/// <param name="Editor">The editor that was applied.</param>
/// <param name="Protospacer">The unedited protospacer, 5′→3′ on the targeted strand.</param>
/// <param name="EditedProtospacer">
/// The protospacer with <em>every</em> editable base in the window converted. Base editors are processive on
/// the exposed R-loop, so when several target bases sit in the window the fully converted product is a real
/// and usually major outcome; this is the conservative worst case to inspect.
/// </param>
/// <param name="IntendedOnlyProtospacer">
/// The protospacer with only <see cref="Intended"/> converted, i.e. the product you were aiming for. Equal to
/// <paramref name="EditedProtospacer"/> when there are no bystanders.
/// </param>
/// <param name="Sites">Every editable base inside the window, ordered by position.</param>
/// <param name="Intended">The site the caller asked for, or the highest-weighted site when none was requested; null when the window contains no editable base.</param>
/// <param name="Bystanders">Every site other than <paramref name="Intended"/>.</param>
/// <param name="IsUnambiguous">True when the window contains exactly one editable base, so the edit cannot be confused with a bystander.</param>
public sealed record BaseEditOutcome(
    BaseEditor Editor,
    string Protospacer,
    string EditedProtospacer,
    string IntendedOnlyProtospacer,
    IReadOnlyList<BaseEditSite> Sites,
    BaseEditSite? Intended,
    IReadOnlyList<BaseEditSite> Bystanders,
    bool IsUnambiguous)
{
    /// <summary>True when at least one base in the window can be edited.</summary>
    public bool HasEdit => Sites.Count > 0;

    /// <summary>Number of bystander edits predicted alongside the intended one.</summary>
    public int BystanderCount => Bystanders.Count;

    /// <inheritdoc/>
    public override string ToString()
        => Sites.Count == 0
            ? $"{Editor.Name}: no {Editor.TargetBase} in window {Editor.EditingWindowStart}-{Editor.EditingWindowEnd}"
            : $"{Editor.Name}: {string.Join(", ", Sites)}";
}

/// <summary>
/// The protein-level reading of one predicted base change inside a coding sequence.
/// </summary>
/// <param name="Site">The predicted base change on the protospacer strand.</param>
/// <param name="CdsPosition">0-based index of the changed base in the supplied coding sequence (plus strand of the CDS).</param>
/// <param name="CodonIndex">0-based codon number within the coding sequence.</param>
/// <param name="PositionInCodon">0, 1 or 2 — which base of the codon changed.</param>
/// <param name="ReferenceCodon">The codon before editing, read on the coding strand.</param>
/// <param name="EditedCodon">The codon after editing, read on the coding strand.</param>
/// <param name="ReferenceAminoAcid">The amino acid before editing, <c>*</c> for a stop.</param>
/// <param name="EditedAminoAcid">The amino acid after editing, <c>*</c> for a stop.</param>
/// <param name="Consequence">The classification of the change.</param>
public sealed record ProteinConsequence(
    BaseEditSite Site,
    int CdsPosition,
    int CodonIndex,
    int PositionInCodon,
    string ReferenceCodon,
    string EditedCodon,
    char ReferenceAminoAcid,
    char EditedAminoAcid,
    EditConsequence Consequence)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"{ReferenceCodon}>{EditedCodon} p.{ReferenceAminoAcid}{CodonIndex + 1}{EditedAminoAcid} {Consequence}";
}

/// <summary>
/// A base-editing prediction mapped onto a coding sequence.
/// </summary>
/// <param name="Outcome">The protospacer-level prediction.</param>
/// <param name="Consequences">One entry per editable base in the window, each evaluated as a lone edit.</param>
/// <param name="EditedCds">The coding sequence with every window edit applied at once.</param>
public sealed record BaseEditCodingOutcome(
    BaseEditOutcome Outcome,
    IReadOnlyList<ProteinConsequence> Consequences,
    string EditedCds)
{
    /// <summary>The consequence of the intended edit, or null when the window held no editable base.</summary>
    public ProteinConsequence? IntendedConsequence => Consequences.FirstOrDefault(c => c.Site.IsIntended);
}

/// <summary>
/// A guide whose editing window covers a requested base, together with the predicted outcome and a rank score.
/// </summary>
/// <param name="Guide">The candidate guide.</param>
/// <param name="Editor">The editor it was ranked for.</param>
/// <param name="ProtospacerPosition">1-based PAM-distal position of the requested base within this guide's protospacer.</param>
/// <param name="Outcome">The predicted editing outcome for this guide.</param>
/// <param name="Score">The rank score; see <see cref="BaseEditing.FindEditingGuides"/> for the formula.</param>
public sealed record EditingGuideCandidate(
    Guide Guide,
    BaseEditor Editor,
    int ProtospacerPosition,
    BaseEditOutcome Outcome,
    double Score)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"{Guide.Spacer} {Guide.Pam} pos {ProtospacerPosition} "
           + $"bystanders {Outcome.BystanderCount} score {Score:F3}";
}

/// <summary>
/// Outcome prediction and guide selection for cytosine and adenine base editors.
/// </summary>
/// <remarks>
/// See <see cref="BaseEditor"/> for the position and strand conventions, which this whole class follows.
/// </remarks>
public static class BaseEditing
{
    /// <summary>
    /// Predicts what <paramref name="editor"/> does to a protospacer.
    /// </summary>
    /// <param name="editor">The base editor.</param>
    /// <param name="protospacer">
    /// The protospacer, 5′→3′ on the targeted strand, exactly as <see cref="Guide.Spacer"/> reports it.
    /// Its length must match <see cref="BaseEditor.Nuclease"/>'s guide length, since positions are counted
    /// from the PAM-distal end.
    /// </param>
    /// <param name="intendedPosition">
    /// The 1-based PAM-distal protospacer position the caller wants edited, or 0 to let the highest-weighted
    /// editable position in the window be treated as the intended one.
    /// </param>
    /// <returns>The predicted outcome, including bystanders.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="editor"/> or <paramref name="protospacer"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="protospacer"/> is not the nuclease's guide length.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="intendedPosition"/> is outside the protospacer.</exception>
    public static BaseEditOutcome PredictEdits(BaseEditor editor, string protospacer, int intendedPosition = 0)
    {
        ArgumentNullException.ThrowIfNull(editor);
        ArgumentNullException.ThrowIfNull(protospacer);

        string spacer = DnaUtility.NormalizeDna(protospacer);
        int length = editor.Nuclease.GuideLength;
        if (spacer.Length != length)
            throw new ArgumentException(
                $"Protospacer must be {length} nt for {editor.Nuclease.Name}, got {spacer.Length}.", nameof(protospacer));
        if (intendedPosition < 0 || intendedPosition > length)
            throw new ArgumentOutOfRangeException(nameof(intendedPosition), intendedPosition,
                "Intended position must be a 1-based protospacer position, or 0 for automatic.");

        char target = editor.TargetBase;
        char product = editor.ProductBase;

        var positions = new List<int>();
        for (int pos = editor.EditingWindowStart; pos <= editor.EditingWindowEnd; pos++)
        {
            if (pos < 1 || pos > length) continue;
            if (spacer[pos - 1] == target) positions.Add(pos);
        }

        int intended = intendedPosition;
        if (intended != 0 && !positions.Contains(intended)) intended = 0;
        if (intended == 0 && positions.Count > 0)
        {
            intended = positions[0];
            foreach (int pos in positions)
                if (editor.WeightAt(pos) > editor.WeightAt(intended)) intended = pos;
        }

        var sites = new List<BaseEditSite>(positions.Count);
        foreach (int pos in positions)
            sites.Add(new BaseEditSite(pos, target, product, editor.WeightAt(pos), pos == intended));

        var edited = new StringBuilder(spacer);
        foreach (int pos in positions) edited[pos - 1] = product;

        var intendedOnly = new StringBuilder(spacer);
        if (intended != 0) intendedOnly[intended - 1] = product;

        var intendedSite = sites.FirstOrDefault(s => s.IsIntended);
        var bystanders = sites.Where(s => !s.IsIntended).ToArray();

        return new BaseEditOutcome(
            editor, spacer, edited.ToString(), intendedOnly.ToString(),
            sites, intendedSite, bystanders, sites.Count == 1);
    }

    /// <summary>
    /// Predicts what <paramref name="editor"/> does to the protospacer of <paramref name="guide"/>.
    /// </summary>
    /// <param name="editor">The base editor.</param>
    /// <param name="guide">The guide whose spacer is the protospacer strand sequence.</param>
    /// <param name="intendedPosition">
    /// The 1-based PAM-distal protospacer position the caller wants edited, or 0 for automatic selection.
    /// </param>
    /// <returns>The predicted outcome, including bystanders.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="guide"/> is null.</exception>
    public static BaseEditOutcome PredictEdits(BaseEditor editor, Guide guide, int intendedPosition = 0)
    {
        ArgumentNullException.ThrowIfNull(guide);
        return PredictEdits(editor, guide.Spacer, intendedPosition);
    }

    /// <summary>
    /// The 1-based PAM-distal protospacer position of a plus-strand target index within a guide, or 0 when the
    /// index lies outside the guide's protospacer.
    /// </summary>
    /// <param name="guide">The guide.</param>
    /// <param name="targetIndexZeroBased">A 0-based index on the plus strand of the searched sequence.</param>
    /// <returns>The 1-based protospacer position, or 0 if the index is not covered.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="guide"/> is null.</exception>
    public static int ProtospacerPosition(Guide guide, int targetIndexZeroBased)
    {
        ArgumentNullException.ThrowIfNull(guide);
        if (targetIndexZeroBased < guide.Start || targetIndexZeroBased >= guide.End) return 0;
        return guide.Strand == Strand.Reverse
            ? guide.End - targetIndexZeroBased
            : targetIndexZeroBased - guide.Start + 1;
    }

    /// <summary>
    /// The stretch of the plus strand covered by an editor's activity window for a given guide, as a 0-based
    /// half-open interval. Strand follows the guide, since the deaminase acts on the protospacer strand.
    /// </summary>
    /// <param name="editor">The base editor whose window is projected.</param>
    /// <param name="guide">The guide whose protospacer anchors the window.</param>
    /// <returns>The window as a plus-strand interval on <see cref="Guide.SequenceId"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="editor"/> or <paramref name="guide"/> is null.</exception>
    public static GenomicInterval EditingWindowInterval(BaseEditor editor, Guide guide)
    {
        ArgumentNullException.ThrowIfNull(editor);
        ArgumentNullException.ThrowIfNull(guide);

        int lo = Math.Max(1, editor.EditingWindowStart);
        int hi = Math.Min(guide.Nuclease.GuideLength, editor.EditingWindowEnd);

        return guide.Strand == Strand.Reverse
            ? new GenomicInterval(guide.SequenceId, guide.End - hi, guide.End - lo + 1, guide.Strand)
            : new GenomicInterval(guide.SequenceId, guide.Start + lo - 1, guide.Start + hi, guide.Strand);
    }

    /// <summary>
    /// Reads a predicted editing outcome at the protein level.
    /// </summary>
    /// <remarks>
    /// <paramref name="cds"/> is a coding sequence on its own plus strand, in frame from index 0. The
    /// protospacer is located inside it by <paramref name="protospacerStartInCds"/> and
    /// <paramref name="strand"/>. For a forward-strand protospacer, protospacer position <c>p</c> is
    /// <c>cds[protospacerStartInCds + p - 1]</c> and a CBE writes C→T there directly. For a reverse-strand
    /// protospacer, position <c>p</c> is <c>cds[protospacerStartInCds + guideLength - p]</c> and the same CBE
    /// event reads as G→A on the coding strand. Each editable base is evaluated as a lone edit; the combined
    /// product is returned separately as <see cref="BaseEditCodingOutcome.EditedCds"/>.
    /// </remarks>
    /// <param name="editor">The base editor.</param>
    /// <param name="cds">The coding sequence, in frame from index 0.</param>
    /// <param name="protospacerStartInCds">0-based index of the leftmost protospacer base on the plus strand of <paramref name="cds"/>.</param>
    /// <param name="strand">Which strand the protospacer lies on.</param>
    /// <param name="intendedPosition">The 1-based PAM-distal protospacer position the caller wants edited, or 0 for automatic.</param>
    /// <param name="code">The genetic code to translate with; <see cref="GeneticCode.Standard"/> when null.</param>
    /// <returns>The outcome plus one consequence per editable base in the window.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="editor"/> or <paramref name="cds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The protospacer does not fit inside <paramref name="cds"/>.</exception>
    public static BaseEditCodingOutcome PredictProteinConsequence(
        BaseEditor editor,
        string cds,
        int protospacerStartInCds,
        Strand strand = Strand.Forward,
        int intendedPosition = 0,
        GeneticCode? code = null)
    {
        ArgumentNullException.ThrowIfNull(editor);
        ArgumentNullException.ThrowIfNull(cds);
        code ??= GeneticCode.Standard;

        string coding = DnaUtility.NormalizeDna(cds);
        int guideLength = editor.Nuclease.GuideLength;
        if (protospacerStartInCds < 0 || protospacerStartInCds + guideLength > coding.Length)
            throw new ArgumentOutOfRangeException(nameof(protospacerStartInCds), protospacerStartInCds,
                "The protospacer must lie entirely inside the coding sequence.");

        string protospacer = coding.Substring(protospacerStartInCds, guideLength);
        if (strand == Strand.Reverse) protospacer = DnaUtility.ReverseComplement(protospacer);

        var outcome = PredictEdits(editor, protospacer, intendedPosition);

        // The change written on the coding strand: direct for a forward protospacer, complemented otherwise.
        char codingFrom = strand == Strand.Reverse ? editor.ComplementaryTargetBase : editor.TargetBase;
        char codingTo = strand == Strand.Reverse ? editor.ComplementaryProductBase : editor.ProductBase;

        var consequences = new List<ProteinConsequence>(outcome.Sites.Count);
        var combined = new StringBuilder(coding);

        foreach (var site in outcome.Sites)
        {
            int cdsIndex = strand == Strand.Reverse
                ? protospacerStartInCds + guideLength - site.Position
                : protospacerStartInCds + site.Position - 1;

            combined[cdsIndex] = codingTo;

            int codonIndex = cdsIndex / 3;
            int inCodon = cdsIndex % 3;
            int codonStart = codonIndex * 3;

            if (codonStart + 3 > coding.Length)
            {
                consequences.Add(new ProteinConsequence(
                    site, cdsIndex, codonIndex, inCodon, "", "", '?', '?', EditConsequence.NonCoding));
                continue;
            }

            string refCodon = coding.Substring(codonStart, 3);
            if (refCodon[inCodon] != codingFrom)
            {
                // Cannot happen for a protospacer taken from this CDS; defensive only.
                consequences.Add(new ProteinConsequence(
                    site, cdsIndex, codonIndex, inCodon, refCodon, refCodon, '?', '?', EditConsequence.NonCoding));
                continue;
            }

            char[] alt = refCodon.ToCharArray();
            alt[inCodon] = codingTo;
            string altCodon = new(alt);

            char refAa = code.Translate(refCodon);
            char altAa = code.Translate(altCodon);

            EditConsequence consequence;
            if (refAa == altAa) consequence = EditConsequence.Silent;
            else if (altAa == '*') consequence = EditConsequence.Nonsense;
            else if (refAa == '*') consequence = EditConsequence.StopLoss;
            else if (codonIndex == 0 && code.IsStart(refCodon) && !code.IsStart(altCodon)) consequence = EditConsequence.StartLoss;
            else consequence = EditConsequence.Missense;

            consequences.Add(new ProteinConsequence(
                site, cdsIndex, codonIndex, inCodon, refCodon, altCodon, refAa, altAa, consequence));
        }

        return new BaseEditCodingOutcome(outcome, consequences, combined.ToString());
    }

    /// <summary>
    /// Finds and ranks guides whose activity window covers a requested base of <paramref name="target"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A guide qualifies when the requested index falls inside its protospacer, at a 1-based PAM-distal
    /// position inside the editor's window, and the base there (read on the protospacer strand) is the
    /// editor's target base. A plus-strand <c>C</c> therefore qualifies for a forward-strand CBE guide, and a
    /// plus-strand <c>G</c> qualifies for a reverse-strand CBE guide.
    /// </para>
    /// <para>
    /// The rank score is <c>WeightAt(position) − 0.08 × distance from the window centre − 0.15 × bystanders
    /// + 0.5 × on-target score</c>, with the on-target term dropped when it is not available. It is a
    /// deliberately simple preference ordering — centred intended base, few bystanders, decent on-target
    /// score — not a calibrated efficiency.
    /// </para>
    /// </remarks>
    /// <param name="target">The sequence to scan.</param>
    /// <param name="editor">The base editor whose window must cover the requested base.</param>
    /// <param name="targetPositionZeroBased">0-based index on the plus strand of <paramref name="target"/> of the base to edit.</param>
    /// <param name="options">
    /// Guide search options. When null, <see cref="GuideSearchOptions.Default"/> is used with the editor's own
    /// nuclease substituted in; when supplied, the caller's nuclease is respected as given.
    /// </param>
    /// <param name="sequenceId">Name recorded on the returned guides.</param>
    /// <returns>Qualifying candidates, best first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="target"/> or <paramref name="editor"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="targetPositionZeroBased"/> is outside <paramref name="target"/>.</exception>
    public static IReadOnlyList<EditingGuideCandidate> FindEditingGuides(
        string target,
        BaseEditor editor,
        int targetPositionZeroBased,
        GuideSearchOptions? options = null,
        string sequenceId = "target")
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(editor);

        string sequence = DnaUtility.NormalizeDna(target);
        if (targetPositionZeroBased < 0 || targetPositionZeroBased >= sequence.Length)
            throw new ArgumentOutOfRangeException(nameof(targetPositionZeroBased), targetPositionZeroBased,
                "Requested position is outside the target sequence.");

        options ??= new GuideSearchOptions { Nuclease = editor.Nuclease };

        var candidates = new List<EditingGuideCandidate>();
        foreach (var guide in GuideFinder.Find(sequence, options, sequenceId))
        {
            int position = ProtospacerPosition(guide, targetPositionZeroBased);
            if (position == 0 || !editor.IsInWindow(position)) continue;
            if (guide.Spacer[position - 1] != editor.TargetBase) continue;

            var outcome = PredictEdits(editor, guide, position);
            if (outcome.Intended is null) continue;

            double score = editor.WeightAt(position)
                           - 0.08 * Math.Abs(position - editor.WindowCentre)
                           - 0.15 * outcome.BystanderCount;
            if (!double.IsNaN(guide.OnTargetScore)) score += 0.5 * guide.OnTargetScore;

            candidates.Add(new EditingGuideCandidate(guide, editor, position, outcome, score));
        }

        candidates.Sort(static (a, b) => b.Score.CompareTo(a.Score));
        return candidates;
    }
}
