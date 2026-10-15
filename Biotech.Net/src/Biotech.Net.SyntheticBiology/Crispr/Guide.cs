using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Crispr;

/// <summary>
/// A candidate guide RNA found on a target sequence: the spacer, its PAM, where it sits and where it cuts.
/// </summary>
/// <remarks>
/// <see cref="Start"/>/<see cref="End"/> are 0-based half-open coordinates of the protospacer on the plus
/// strand of the searched sequence, regardless of which strand the guide targets — the same convention as
/// <see cref="GenomicInterval"/>. <see cref="Spacer"/> is written 5′→3′ on the targeted strand, which is what
/// you order as an oligo.
/// </remarks>
public sealed class Guide
{
    internal Guide(CasNuclease nuclease, string spacer, string pam, string sequenceId, int start, Strand strand, string? context30)
    {
        Nuclease = nuclease;
        Spacer = spacer;
        Pam = pam;
        SequenceId = sequenceId;
        Start = start;
        Strand = strand;
        Context = context30;
    }

    /// <summary>The nuclease this guide was designed for.</summary>
    public CasNuclease Nuclease { get; }

    /// <summary>The protospacer, 5′→3′ on the targeted strand. This is the sequence you clone or synthesise.</summary>
    public string Spacer { get; }

    /// <summary>The PAM as read 5′→3′ on the targeted strand.</summary>
    public string Pam { get; }

    /// <summary>Identifier of the sequence the guide was found on.</summary>
    public string SequenceId { get; }

    /// <summary>0-based start of the protospacer on the plus strand.</summary>
    public int Start { get; }

    /// <summary>Exclusive end of the protospacer on the plus strand.</summary>
    public int End => Start + Nuclease.GuideLength;

    /// <summary>Which strand the guide anneals to.</summary>
    public Strand Strand { get; }

    /// <summary>
    /// The 30-mer scoring context (4 nt upstream + 20 nt spacer + 3 nt PAM + 3 nt downstream, on the targeted
    /// strand) required by <see cref="OnTargetScorer"/>, or null when the flanks ran off the end of the sequence.
    /// </summary>
    public string? Context { get; }

    /// <summary>
    /// The 0-based plus-strand coordinate just 3′ of the double-strand break; the break lies between
    /// <c>CutSite - 1</c> and <c>CutSite</c>.
    /// </summary>
    public int CutSite => Nuclease.CutSite(Start, Strand);

    /// <summary>The protospacer as a half-open interval on the searched sequence.</summary>
    public GenomicInterval Interval => new(SequenceId, Start, End, Strand);

    /// <summary>GC fraction of the spacer (0–1).</summary>
    public double GcContent => DnaUtility.GcFraction(Spacer);

    /// <summary>Length of the longest single-base run in the spacer.</summary>
    public int LongestHomopolymer => DnaUtility.LongestHomopolymer(DnaUtility.Ascii(Spacer));

    /// <summary>
    /// True when the spacer contains <c>TTTT</c>, which terminates Pol III transcription of a U6-driven sgRNA
    /// and usually kills activity.
    /// </summary>
    public bool HasPolyTTerminator => Spacer.Contains("TTTT", StringComparison.OrdinalIgnoreCase);

    /// <summary>Doench 2014 Rule Set 1 on-target score (0–1), or <see cref="double.NaN"/> when it could not be computed.</summary>
    public double OnTargetScore { get; internal set; } = double.NaN;

    /// <summary>Spacer plus PAM in the order they occur on the targeted strand.</summary>
    public string FullSite => Nuclease.PamSide == PamSide.ThreePrime ? Spacer + Pam : Pam + Spacer;

    /// <summary>The oligo to clone into a U6 cassette: the spacer with a leading G added when it does not start with one.</summary>
    /// <remarks>
    /// Human/murine U6 transcription starts efficiently only from G. Adding a 5′ G gives a 21 nt spacer that is
    /// tolerated by SpCas9; the alternative is to substitute position 1, which costs activity.
    /// </remarks>
    public string U6CloningOligo => Spacer.StartsWith('G') ? Spacer : "G" + Spacer;

    /// <inheritdoc/>
    public override string ToString()
        => $"{Spacer} {Pam} {SequenceId}:{Start + 1}-{End}({(Strand == Strand.Reverse ? '-' : '+')})"
           + (double.IsNaN(OnTargetScore) ? "" : $" on-target {OnTargetScore:F3}");
}
