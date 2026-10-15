using Biotech.Net.Core;
using Biotech.Net.SyntheticBiology.Primers;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// The kind of change a pegRNA encodes.
/// </summary>
public enum PrimeEditKind
{
    /// <summary>Replace <c>n</c> bases with <c>n</c> other bases (<c>n</c> is usually 1).</summary>
    Substitution,

    /// <summary>Insert new bases, replacing nothing.</summary>
    Insertion,

    /// <summary>Remove bases, inserting nothing.</summary>
    Deletion,
}

/// <summary>
/// A sequence change expressed against a plus-strand reference, 0-based and half-open like
/// <see cref="GenomicInterval"/>: the <paramref name="ReferenceAllele"/> occupies
/// <c>[Position, Position + ReferenceAllele.Length)</c> and is replaced by <paramref name="AlternateAllele"/>.
/// </summary>
/// <param name="Kind">Whether the change is a substitution, insertion or deletion.</param>
/// <param name="Position">0-based plus-strand index where the replaced stretch starts.</param>
/// <param name="ReferenceAllele">The bases removed; empty for an insertion.</param>
/// <param name="AlternateAllele">The bases written in their place; empty for a deletion.</param>
public sealed record PrimeEdit(PrimeEditKind Kind, int Position, string ReferenceAllele, string AlternateAllele)
{
    /// <summary>0-based half-open end of the replaced stretch on the reference.</summary>
    public int ReferenceEnd => Position + ReferenceAllele.Length;

    /// <summary>0-based half-open end of the new stretch in the edited sequence (the edit starts at the same index).</summary>
    public int AlternateEnd => Position + AlternateAllele.Length;

    /// <summary>Net change in sequence length: positive for an insertion, negative for a deletion.</summary>
    public int LengthDelta => AlternateAllele.Length - ReferenceAllele.Length;

    /// <summary>
    /// Applies the edit to a reference sequence.
    /// </summary>
    /// <param name="reference">The unedited plus-strand sequence.</param>
    /// <returns>The edited sequence.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> is null.</exception>
    /// <exception cref="ArgumentException">The reference allele does not match <paramref name="reference"/> at <see cref="Position"/>.</exception>
    public string Apply(string reference)
    {
        ArgumentNullException.ThrowIfNull(reference);
        if (Position < 0 || ReferenceEnd > reference.Length)
            throw new ArgumentException("Edit lies outside the reference sequence.", nameof(reference));
        if (!reference.AsSpan(Position, ReferenceAllele.Length)
                .Equals(ReferenceAllele, StringComparison.OrdinalIgnoreCase))
            throw new ArgumentException(
                $"Reference allele '{ReferenceAllele}' does not match the sequence at {Position}.", nameof(reference));

        return string.Concat(reference.AsSpan(0, Position), AlternateAllele, reference.AsSpan(ReferenceEnd));
    }

    /// <summary>
    /// A substitution replacing <c>replacement.Length</c> bases of <paramref name="reference"/> at
    /// <paramref name="position"/>.
    /// </summary>
    /// <param name="reference">The unedited plus-strand sequence, used to record the reference allele.</param>
    /// <param name="position">0-based plus-strand index of the first replaced base.</param>
    /// <param name="replacement">The new bases.</param>
    /// <returns>The encoded edit.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> or <paramref name="replacement"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The replaced stretch runs off the end of <paramref name="reference"/>.</exception>
    public static PrimeEdit Substitution(string reference, int position, string replacement)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(replacement);
        string alt = DnaUtility.NormalizeDna(replacement);
        if (position < 0 || position + alt.Length > reference.Length)
            throw new ArgumentOutOfRangeException(nameof(position), position, "Substitution runs off the reference.");
        return new PrimeEdit(PrimeEditKind.Substitution, position,
            DnaUtility.NormalizeDna(reference.Substring(position, alt.Length)), alt);
    }

    /// <summary>
    /// An insertion of <paramref name="inserted"/> immediately before <paramref name="position"/>.
    /// </summary>
    /// <param name="position">0-based plus-strand index the new bases are inserted in front of.</param>
    /// <param name="inserted">The bases to insert.</param>
    /// <returns>The encoded edit.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="inserted"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="position"/> is negative.</exception>
    public static PrimeEdit Insertion(int position, string inserted)
    {
        ArgumentNullException.ThrowIfNull(inserted);
        ArgumentOutOfRangeException.ThrowIfNegative(position);
        return new PrimeEdit(PrimeEditKind.Insertion, position, "", DnaUtility.NormalizeDna(inserted));
    }

    /// <summary>
    /// A deletion of <paramref name="length"/> bases starting at <paramref name="position"/>.
    /// </summary>
    /// <param name="reference">The unedited plus-strand sequence, used to record the deleted bases.</param>
    /// <param name="position">0-based plus-strand index of the first deleted base.</param>
    /// <param name="length">Number of bases to delete.</param>
    /// <returns>The encoded edit.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The deletion runs off the end of <paramref name="reference"/>.</exception>
    public static PrimeEdit Deletion(string reference, int position, int length)
    {
        ArgumentNullException.ThrowIfNull(reference);
        if (position < 0 || length < 0 || position + length > reference.Length)
            throw new ArgumentOutOfRangeException(nameof(length), length, "Deletion runs off the reference.");
        return new PrimeEdit(PrimeEditKind.Deletion, position,
            DnaUtility.NormalizeDna(reference.Substring(position, length)), "");
    }

    /// <summary>
    /// Builds an edit from the loose <c>(position, edit string, kind)</c> form used by
    /// <see cref="PrimeEditor.Design(string, int, string, PrimeEditorOptions?)"/>.
    /// </summary>
    /// <param name="reference">The unedited plus-strand sequence.</param>
    /// <param name="position">0-based plus-strand index where the edit starts.</param>
    /// <param name="edit">
    /// For a substitution, the replacement bases. For an insertion, the inserted bases. For a deletion, either
    /// the bases to delete (which must match the reference) or a decimal count of bases to delete.
    /// </param>
    /// <param name="kind">Which of the three readings of <paramref name="edit"/> applies.</param>
    /// <returns>The encoded edit.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reference"/> or <paramref name="edit"/> is null.</exception>
    /// <exception cref="ArgumentException">A deletion's bases do not match the reference, or <paramref name="kind"/> is not a known value.</exception>
    public static PrimeEdit Parse(string reference, int position, string edit, PrimeEditKind kind)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(edit);

        switch (kind)
        {
            case PrimeEditKind.Substitution:
                return Substitution(reference, position, edit);

            case PrimeEditKind.Insertion:
                return Insertion(position, edit);

            case PrimeEditKind.Deletion:
                if (int.TryParse(edit, out int count)) return Deletion(reference, position, count);
                string deleted = DnaUtility.NormalizeDna(edit);
                var made = Deletion(reference, position, deleted.Length);
                if (made.ReferenceAllele != deleted)
                    throw new ArgumentException(
                        $"Deletion '{deleted}' does not match the reference ('{made.ReferenceAllele}') at {position}.", nameof(edit));
                return made;

            default:
                throw new ArgumentException($"Unknown edit kind {kind}.", nameof(kind));
        }
    }

    /// <inheritdoc/>
    public override string ToString()
        => Kind switch
        {
            PrimeEditKind.Insertion => $"{Position + 1}ins{AlternateAllele}",
            PrimeEditKind.Deletion => $"{Position + 1}del{ReferenceAllele}",
            _ => $"{Position + 1}{ReferenceAllele}>{AlternateAllele}",
        };
}

/// <summary>
/// Parameters for pegRNA and nicking-sgRNA design.
/// </summary>
/// <remarks>
/// Defaults follow the ranges Anzalone <i>et al.</i> 2019 (Nature 576:149) screened and the practical advice
/// in their supplementary protocol: PBS 8–17 nt scanned around a melting temperature near 38 °C, RTT 10–30 nt
/// with at least 5 nt of homology past the edit, and a nick no more than about 30 nt 5′ of the edit.
/// </remarks>
public sealed class PrimeEditorOptions
{
    /// <summary>The nickase geometry to design against. PE2/PE3 use SpCas9(H840A), which has SpCas9 geometry.</summary>
    public CasNuclease Nuclease { get; init; } = CasNuclease.SpCas9;

    /// <summary>How to read the <c>edit</c> string passed to <see cref="PrimeEditor.Design(string, int, string, PrimeEditorOptions?)"/>.</summary>
    public PrimeEditKind EditKind { get; init; } = PrimeEditKind.Substitution;

    /// <summary>Shortest primer binding site to try, in nt.</summary>
    public int PbsMinLength { get; init; } = 8;

    /// <summary>Longest primer binding site to try, in nt.</summary>
    public int PbsMaxLength { get; init; } = 17;

    /// <summary>Shortest reverse-transcriptase template to try, in nt.</summary>
    public int RttMinLength { get; init; } = 10;

    /// <summary>Longest reverse-transcriptase template to try, in nt.</summary>
    public int RttMaxLength { get; init; } = 30;

    /// <summary>
    /// Largest acceptable number of nucleotides between the nick and the first edited base, measured 3′ of the
    /// nick along the nicked strand. Efficiency falls off sharply beyond roughly 30 nt.
    /// </summary>
    public int MaxNickToEditDistance { get; init; } = 30;

    /// <summary>
    /// Melting temperature the PBS is scanned towards, in °C. Anzalone <i>et al.</i> recommend designing the
    /// PBS to about 38 °C, lengthening it for AT-rich and shortening it for GC-rich primer sites.
    /// </summary>
    public double PbsTargetTm { get; init; } = 38.0;

    /// <summary>Minimum number of templated bases the RTT must extend past the edit, to give the 3′ flap homology to anneal with.</summary>
    public int MinHomologyAfterEdit { get; init; } = 5;

    /// <summary>Smallest acceptable PE3 nick-to-nick distance, in nt.</summary>
    public int NickingDistanceMin { get; init; } = 40;

    /// <summary>Largest acceptable PE3 nick-to-nick distance, in nt.</summary>
    public int NickingDistanceMax { get; init; } = 90;

    /// <summary>Maximum number of ranked results returned by each design method.</summary>
    public int MaxResults { get; init; } = 20;

    /// <summary>Search protospacers on the plus strand.</summary>
    public bool SearchForward { get; init; } = true;

    /// <summary>Search protospacers on the minus strand.</summary>
    public bool SearchReverse { get; init; } = true;

    /// <summary>The sgRNA scaffold placed between spacer and 3′ extension; see <see cref="PegRna.DefaultScaffold"/>.</summary>
    public string Scaffold { get; init; } = PegRna.DefaultScaffold;

    /// <summary>Salt and oligo conditions used for the PBS melting temperature; <see cref="SaltConditions.Default"/> when null.</summary>
    public SaltConditions? SaltConditions { get; init; }

    /// <summary>
    /// Guide search filters. When null, a permissive set is used (no GC filter, poly-T rejected, no on-target
    /// scoring), because a prime edit must be made at a fixed position and there is rarely a spare protospacer.
    /// </summary>
    public GuideSearchOptions? GuideOptions { get; init; }

    /// <summary>A fresh options object with every default in place.</summary>
    public static PrimeEditorOptions Default { get; } = new();

    internal GuideSearchOptions EffectiveGuideOptions()
        => GuideOptions ?? new GuideSearchOptions
        {
            Nuclease = Nuclease,
            SearchForward = SearchForward,
            SearchReverse = SearchReverse,
            MinGcContent = 0.0,
            MaxGcContent = 1.0,
            RejectPolyT = true,
            MaxHomopolymer = 0,
            ScoreOnTarget = false,
        };
}

/// <summary>
/// A designed pegRNA: a prime-editing guide RNA that both targets the nickase and templates the new sequence.
/// </summary>
/// <remarks>
/// <para>
/// <b>Geometry (Anzalone <i>et al.</i> 2019, Nature 576:149).</b> SpCas9(H840A) nicks the PAM-containing
/// strand — the strand whose sequence the spacer reproduces — 3 bp from the PAM, i.e. between protospacer
/// positions 17 and 18. The liberated 3′ end anneals to the pegRNA's primer binding site, and the reverse
/// transcriptase then copies the reverse-transcriptase template into a new 3′ flap carrying the edit.
/// </para>
/// <para>
/// <b>Reading the parts.</b> Write the nicked strand 5′→3′ and call the first base after the nick index
/// <c>n</c>. Then <see cref="Pbs"/> is the reverse complement of the <see cref="PbsLength"/> bases ending at
/// the nick (so it anneals to the 3′ flap), and <see cref="Rtt"/> is the reverse complement of the first
/// <see cref="RttLength"/> bases of the <em>edited</em> strand from <c>n</c> onwards. Templating from the
/// edited sequence is the whole point: the RTT is what writes the new allele.
/// </para>
/// <para>
/// The 3′ extension is <c>RTT + PBS</c> read 5′→3′, which is equivalently the reverse complement of the edited
/// strand across <c>[n − PbsLength, n + RttLength)</c>, and the full oligo is
/// <c>spacer + scaffold + extension</c>.
/// </para>
/// </remarks>
/// <param name="Nuclease">The nickase geometry the pegRNA was designed against.</param>
/// <param name="Spacer">The 20 nt protospacer, 5′→3′ on the nicked strand.</param>
/// <param name="Pam">The PAM as read on the nicked strand.</param>
/// <param name="ProtospacerStart">0-based plus-strand index of the leftmost protospacer base.</param>
/// <param name="Strand">Which strand the protospacer, and therefore the nick, lies on.</param>
/// <param name="NickPosition">
/// 0-based plus-strand index of the first base after the nick, the same half-open convention as
/// <see cref="CasNuclease.CutSite"/>.
/// </param>
/// <param name="Pbs">The primer binding site, 5′→3′, at the 3′ end of the pegRNA.</param>
/// <param name="Rtt">The reverse-transcriptase template, 5′→3′, immediately 5′ of the PBS.</param>
/// <param name="Edit">The edit this pegRNA encodes, expressed against the plus strand of the target.</param>
/// <param name="NickToEditDistance">Nucleotides from the nick to the first edited base, along the nicked strand.</param>
/// <param name="HomologyAfterEdit">Templated nucleotides the RTT extends past the edit.</param>
/// <param name="Scaffold">The sgRNA scaffold between spacer and 3′ extension.</param>
/// <param name="Score">Rank score; see <see cref="PrimeEditor.Design(string, PrimeEdit, PrimeEditorOptions?)"/>.</param>
public sealed record PegRna(
    CasNuclease Nuclease,
    string Spacer,
    string Pam,
    int ProtospacerStart,
    Strand Strand,
    int NickPosition,
    string Pbs,
    string Rtt,
    PrimeEdit Edit,
    int NickToEditDistance,
    int HomologyAfterEdit,
    string Scaffold,
    double Score)
{
    /// <summary>
    /// The standard SpCas9 sgRNA scaffold (the F+E / "flip and extension" form used for pegRNAs). It is a
    /// placeholder: substitute your own scaffold, or an engineered pegRNA scaffold with a 3′ motif such as
    /// evopreQ1, when you order.
    /// </summary>
    public const string DefaultScaffold =
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC";

    /// <summary>Length of the primer binding site, in nt.</summary>
    public int PbsLength => Pbs.Length;

    /// <summary>Length of the reverse-transcriptase template, in nt.</summary>
    public int RttLength => Rtt.Length;

    /// <summary>The pegRNA 3′ extension, <c>RTT + PBS</c> read 5′→3′.</summary>
    public string ThreePrimeExtension => Rtt + Pbs;

    /// <summary>The whole pegRNA as DNA, <c>spacer + scaffold + 3′ extension</c>.</summary>
    public string FullSequence => Spacer + Scaffold + ThreePrimeExtension;

    /// <summary>
    /// The newly synthesised 3′ flap, i.e. the edited strand downstream of the nick that the RTT templates.
    /// Equal to <c>ReverseComplement(Rtt)</c> by construction.
    /// </summary>
    public string EditedFlap => DnaUtility.ReverseComplement(Rtt);

    /// <summary>
    /// Predicted melting temperature of the primer binding site, in °C, from the nearest-neighbour DNA model.
    /// The real duplex is RNA:DNA, which is somewhat more stable, so treat this as a relative ranking value.
    /// </summary>
    public double PbsMeltingTemperature => NearestNeighbourTm.Tm(Pbs);

    /// <summary>
    /// The 5′-most base of the 3′ extension — the far end of the RTT. Anzalone <i>et al.</i> advise against a
    /// C here, which appears to destabilise the extension.
    /// </summary>
    public char RttTerminalBase => Rtt.Length > 0 ? Rtt[0] : '\0';

    /// <summary>True when <see cref="RttTerminalBase"/> is not C, the preferred case.</summary>
    public bool HasPreferredRttEnd => RttTerminalBase != 'C';

    /// <summary>The protospacer as a 0-based half-open plus-strand interval.</summary>
    /// <param name="sequenceId">Name to record on the interval.</param>
    /// <returns>The protospacer interval.</returns>
    public GenomicInterval ProtospacerInterval(string sequenceId = "target")
        => new(sequenceId, ProtospacerStart, ProtospacerStart + Nuclease.GuideLength, Strand);

    /// <inheritdoc/>
    public override string ToString()
        => $"pegRNA {Spacer} {Pam} ({(Strand == Strand.Reverse ? '-' : '+')}) "
           + $"PBS {PbsLength} RTT {RttLength} {Edit} score {Score:F2}";
}

/// <summary>
/// A PE3 or PE3b nicking sgRNA: an ordinary sgRNA for the same nickase that cuts the <em>non-edited</em>
/// strand, so that strand is preferentially repaired using the edited one as template.
/// </summary>
/// <param name="Guide">The nicking guide.</param>
/// <param name="NickPosition">0-based plus-strand index of the first base after this guide's nick.</param>
/// <param name="DistanceFromPegNick">Absolute distance in nt between this nick and the pegRNA nick.</param>
/// <param name="IsPe3B">
/// True when the guide's protospacer overlaps the edit, so the spacer matches only the edited allele. Such a
/// guide can only nick after the edit is installed, which largely removes the indel burden of PE3.
/// </param>
/// <param name="Score">Rank score; see <see cref="PrimeEditor.DesignNickingGuides"/>.</param>
public sealed record NickingGuide(
    Guide Guide,
    int NickPosition,
    int DistanceFromPegNick,
    bool IsPe3B,
    double Score)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"{(IsPe3B ? "PE3b" : "PE3")} {Guide.Spacer} {Guide.Pam} "
           + $"({(Guide.Strand == Strand.Reverse ? '-' : '+')}) {DistanceFromPegNick} nt away score {Score:F2}";
}

/// <summary>
/// pegRNA design for prime editing (Anzalone <i>et al.</i> 2019, Nature 576:149).
/// </summary>
public static class PrimeEditor
{
    /// <summary>
    /// Designs and ranks pegRNAs for an edit described loosely by position, edit string and
    /// <see cref="PrimeEditorOptions.EditKind"/>.
    /// </summary>
    /// <param name="target">The unedited plus-strand target sequence.</param>
    /// <param name="editPositionZeroBased">0-based plus-strand index where the edit starts.</param>
    /// <param name="edit">The replacement, inserted or deleted bases; see <see cref="PrimeEdit.Parse"/>.</param>
    /// <param name="options">Design parameters, or null for <see cref="PrimeEditorOptions.Default"/>.</param>
    /// <returns>Ranked pegRNA candidates, best first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="target"/> or <paramref name="edit"/> is null.</exception>
    public static IReadOnlyList<PegRna> Design(
        string target, int editPositionZeroBased, string edit, PrimeEditorOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(edit);
        options ??= PrimeEditorOptions.Default;
        string reference = DnaUtility.NormalizeDna(target);
        return Design(reference, PrimeEdit.Parse(reference, editPositionZeroBased, edit, options.EditKind), options);
    }

    /// <summary>
    /// Designs and ranks pegRNAs for an explicit edit.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Every protospacer of the configured nuclease is considered on both strands. A protospacer qualifies
    /// when its nick lies 5′ of the edit on the nicked strand, no further than
    /// <see cref="PrimeEditorOptions.MaxNickToEditDistance"/> away. For each qualifying protospacer the method
    /// scans PBS lengths over <c>[PbsMinLength, PbsMaxLength]</c> and RTT lengths over
    /// <c>[RttMinLength, RttMaxLength]</c>, keeping only RTTs that reach at least
    /// <see cref="PrimeEditorOptions.MinHomologyAfterEdit"/> nt past the edit, and keeps the single best
    /// PBS/RTT pair per protospacer.
    /// </para>
    /// <para>
    /// <b>Score.</b> Starting from 10, the method subtracts <c>0.25</c> per °C that the PBS melting
    /// temperature deviates from <see cref="PrimeEditorOptions.PbsTargetTm"/>, subtracts <c>0.04</c> per nt of
    /// nick-to-edit distance and <c>0.05</c> per nt of RTT length beyond the minimum that still gives the
    /// required homology, adds <c>1.0</c> when the far end of the RTT is not a C, subtracts <c>0.5</c> when
    /// the PBS GC fraction is below 0.3 or above 0.7, and adds <c>2 ×</c> the spacer's on-target score when
    /// one is available. These weights express the published design preferences as an ordering; they are not
    /// a calibrated efficiency model.
    /// </para>
    /// </remarks>
    /// <param name="target">The unedited plus-strand target sequence.</param>
    /// <param name="edit">The edit to install, expressed against <paramref name="target"/>.</param>
    /// <param name="options">Design parameters, or null for <see cref="PrimeEditorOptions.Default"/>.</param>
    /// <returns>Ranked pegRNA candidates, best first, at most <see cref="PrimeEditorOptions.MaxResults"/> of them.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="target"/> or <paramref name="edit"/> is null.</exception>
    public static IReadOnlyList<PegRna> Design(string target, PrimeEdit edit, PrimeEditorOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(edit);
        options ??= PrimeEditorOptions.Default;

        string reference = DnaUtility.NormalizeDna(target);
        string edited = edit.Apply(reference);
        var salt = options.SaltConditions;
        int guideLength = options.Nuclease.GuideLength;

        // The design is strand-symmetric, so each strand is handled in its own "working" frame: the nicked
        // strand written 5′→3′. For the minus strand that frame is the reverse complement of the target, and
        // an index i in it maps back to plus-strand index (length - 1 - i).
        string plusRef = reference, plusEdited = edited;
        string minusRef = DnaUtility.ReverseComplement(reference);
        string minusEdited = DnaUtility.ReverseComplement(edited);

        // Start of the edit in each working frame. In the minus frame the bases before the edit are the
        // plus-strand bases after it, which the edit leaves untouched, so the reference and edited frames
        // agree on this index too.
        int plusEditStart = edit.Position;
        int minusEditStart = reference.Length - edit.Position - edit.ReferenceAllele.Length;

        var results = new List<PegRna>();

        foreach (var guide in GuideFinder.Find(reference, options.EffectiveGuideOptions()))
        {
            bool reverse = guide.Strand == Strand.Reverse;
            if (reverse ? !options.SearchReverse : !options.SearchForward) continue;

            string workRef = reverse ? minusRef : plusRef;
            string workEdited = reverse ? minusEdited : plusEdited;
            int editStart = reverse ? minusEditStart : plusEditStart;
            int editLength = reverse ? edit.AlternateAllele.Length : edit.AlternateAllele.Length;

            int plusNick = guide.CutSite;
            int nick = reverse ? reference.Length - plusNick : plusNick;

            int nickToEdit = editStart - nick;
            if (nickToEdit < 0 || nickToEdit > options.MaxNickToEditDistance) continue;
            if (nick < options.PbsMinLength) continue;

            // Shortest RTT that covers the edit and the required homology past it.
            int minRtt = Math.Max(options.RttMinLength, nickToEdit + editLength + options.MinHomologyAfterEdit);
            if (minRtt > options.RttMaxLength) continue;

            PegRna? best = null;
            for (int pbsLength = options.PbsMinLength; pbsLength <= options.PbsMaxLength; pbsLength++)
            {
                if (nick - pbsLength < 0) break;
                string pbs = DnaUtility.ReverseComplement(workRef.Substring(nick - pbsLength, pbsLength));
                double pbsTm = NearestNeighbourTm.Tm(pbs, salt);
                double pbsGc = DnaUtility.GcFraction(pbs);

                for (int rttLength = minRtt; rttLength <= options.RttMaxLength; rttLength++)
                {
                    if (nick + rttLength > workEdited.Length) break;
                    string rtt = DnaUtility.ReverseComplement(workEdited.Substring(nick, rttLength));
                    int homology = rttLength - nickToEdit - editLength;

                    double score = 10.0
                                   - 0.25 * Math.Abs(pbsTm - options.PbsTargetTm)
                                   - 0.04 * nickToEdit
                                   - 0.05 * (rttLength - minRtt);
                    if (rtt[0] != 'C') score += 1.0;
                    if (pbsGc < 0.3 || pbsGc > 0.7) score -= 0.5;
                    if (!double.IsNaN(guide.OnTargetScore)) score += 2.0 * guide.OnTargetScore;

                    if (best is null || score > best.Score)
                    {
                        best = new PegRna(
                            options.Nuclease, guide.Spacer, guide.Pam, guide.Start, guide.Strand, plusNick,
                            pbs, rtt, edit, nickToEdit, homology, options.Scaffold, score);
                    }
                }
            }

            if (best is not null) results.Add(best);
        }

        results.Sort(static (a, b) => b.Score.CompareTo(a.Score));
        return results.Count > options.MaxResults ? results.GetRange(0, options.MaxResults) : results;
    }

    /// <summary>
    /// Finds and ranks PE3 / PE3b nicking sgRNAs for a designed pegRNA.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A nicking guide must sit on the strand opposite the pegRNA, so that it nicks the non-edited strand, and
    /// its nick must be between <see cref="PrimeEditorOptions.NickingDistanceMin"/> and
    /// <see cref="PrimeEditorOptions.NickingDistanceMax"/> nt from the pegRNA nick — Anzalone <i>et al.</i>
    /// report an optimum around +40 to +90 nt, with roughly +50 nt best on average.
    /// </para>
    /// <para>
    /// A guide whose protospacer overlaps the edit is flagged <see cref="NickingGuide.IsPe3B"/>: its spacer
    /// matches the edited allele only, so the second nick is made after the edit is installed and the
    /// characteristic PE3 indels are largely avoided. PE3b guides are strongly preferred by the score, which
    /// is <c>5 − 0.03 × |distance − 65| + 3 if PE3b + 2 × on-target score</c>.
    /// </para>
    /// </remarks>
    /// <param name="target">The unedited plus-strand target sequence the pegRNA was designed on.</param>
    /// <param name="pegRna">The pegRNA whose nick the second nick is placed relative to.</param>
    /// <param name="options">Design parameters, or null for <see cref="PrimeEditorOptions.Default"/>.</param>
    /// <returns>Ranked nicking guides, best first, at most <see cref="PrimeEditorOptions.MaxResults"/> of them.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="target"/> or <paramref name="pegRna"/> is null.</exception>
    public static IReadOnlyList<NickingGuide> DesignNickingGuides(
        string target, PegRna pegRna, PrimeEditorOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(pegRna);
        options ??= PrimeEditorOptions.Default;

        string reference = DnaUtility.NormalizeDna(target);
        var edit = pegRna.Edit;
        int editStart = edit.Position;
        int editEnd = Math.Max(edit.ReferenceEnd, edit.Position + 1); // an insertion still occupies a junction

        var results = new List<NickingGuide>();

        foreach (var guide in GuideFinder.Find(reference, options.EffectiveGuideOptions()))
        {
            if (guide.Strand == pegRna.Strand) continue;

            int nick = guide.CutSite;
            int distance = Math.Abs(nick - pegRna.NickPosition);
            if (distance < options.NickingDistanceMin || distance > options.NickingDistanceMax) continue;

            bool pe3b = guide.Start < editEnd && guide.End > editStart;

            double score = 5.0 - 0.03 * Math.Abs(distance - 65);
            if (pe3b) score += 3.0;
            if (!double.IsNaN(guide.OnTargetScore)) score += 2.0 * guide.OnTargetScore;

            results.Add(new NickingGuide(guide, nick, distance, pe3b, score));
        }

        results.Sort(static (a, b) => b.Score.CompareTo(a.Score));
        return results.Count > options.MaxResults ? results.GetRange(0, options.MaxResults) : results;
    }
}
