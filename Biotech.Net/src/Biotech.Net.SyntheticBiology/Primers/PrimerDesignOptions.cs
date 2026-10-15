namespace Biotech.Net.SyntheticBiology.Primers;

/// <summary>A region of the template, 0-based and half-open: <c>[Start, Start + Length)</c>.</summary>
/// <param name="Start">0-based first position.</param>
/// <param name="Length">Number of bases.</param>
public readonly record struct TemplateRegion(int Start, int Length)
{
    /// <summary>0-based exclusive end.</summary>
    public int End => Start + Length;
}

/// <summary>
/// Constraints and penalty weights for <see cref="PrimerDesigner"/>, modelled on Primer3's settings.
/// </summary>
/// <remarks>
/// <para>
/// <b>Hard constraints.</b> A primer is rejected when its length, Tm, GC fraction, poly-X run, GC clamp, or any of
/// its ΔG values falls outside the limits below. ΔG limits are the <i>most negative</i> ΔG°37 (kcal/mol)
/// accepted: with <see cref="MaxSelfDimerDeltaG"/> = −9, a primer whose self-dimer is −10 is rejected and one at
/// −6 is accepted. A pair is rejected when its product size, Tm difference, or hetero-dimer ΔG values fall
/// outside the limits, when the primers overlap, or when it does not flank <see cref="Target"/>.
/// </para>
/// <para>
/// <b>Primer penalty</b> (lower is better, 0 = ideal):
/// P = <see cref="TmWeight"/>·|Tm − <see cref="OptTm"/>|
/// + <see cref="SizeWeight"/>·|L − <see cref="OptLength"/>|
/// + <see cref="GcWeight"/>·100·|GC − <see cref="OptGc"/>|
/// + <see cref="SelfDimerWeight"/>·max(0, −ΔG<sub>self-dimer</sub>)
/// + <see cref="ThreePrimeDimerWeight"/>·max(0, −ΔG<sub>3′ self-dimer</sub>)
/// + <see cref="HairpinWeight"/>·max(0, −ΔG<sub>hairpin</sub>).
/// With the defaults (1, 1, 0.05 per % GC, 0.1, 0.2, 0.2 per kcal/mol) a 1 °C Tm miss costs as much as one
/// nucleotide of size deviation, 20 points of GC deviation, a 10 kcal/mol self-dimer or a 5 kcal/mol 3′ dimer or
/// hairpin. The Tm and size weights equal Primer3's defaults (PRIMER_WT_TM_LT/GT, PRIMER_WT_SIZE_LT/GT = 1).
/// </para>
/// <para>
/// <b>Pair penalty</b>: P<sub>pair</sub> = P<sub>forward</sub> + P<sub>reverse</sub>
/// + <see cref="TmDifferenceWeight"/>·|ΔTm|
/// + <see cref="HeteroDimerWeight"/>·max(0, −ΔG<sub>hetero-dimer</sub>)
/// + <see cref="ThreePrimeDimerWeight"/>·max(0, −ΔG<sub>3′ hetero-dimer</sub>)
/// (defaults 0.5 per °C, 0.1 and 0.2 per kcal/mol). All weights must be non-negative.
/// </para>
/// </remarks>
public sealed record PrimerDesignOptions
{
    /// <summary>The defaults below.</summary>
    public static PrimerDesignOptions Default { get; } = new();

    /// <summary>Shortest primer, nt (default 18).</summary>
    public int MinLength { get; init; } = 18;

    /// <summary>Optimal primer length, nt (default 20).</summary>
    public int OptLength { get; init; } = 20;

    /// <summary>Longest primer, nt (default 27).</summary>
    public int MaxLength { get; init; } = 27;

    /// <summary>Lowest accepted Tm, °C (default 57).</summary>
    public double MinTm { get; init; } = 57.0;

    /// <summary>Optimal Tm, °C (default 60).</summary>
    public double OptTm { get; init; } = 60.0;

    /// <summary>Highest accepted Tm, °C (default 63).</summary>
    public double MaxTm { get; init; } = 63.0;

    /// <summary>Largest accepted |Tm(forward) − Tm(reverse)|, °C (default 3).</summary>
    public double MaxTmDifference { get; init; } = 3.0;

    /// <summary>Lowest accepted GC fraction (default 0.3).</summary>
    public double MinGc { get; init; } = 0.3;

    /// <summary>Optimal GC fraction for the penalty (default 0.5).</summary>
    public double OptGc { get; init; } = 0.5;

    /// <summary>Highest accepted GC fraction (default 0.8).</summary>
    public double MaxGc { get; init; } = 0.8;

    /// <summary>Longest accepted single-base run (default 5, Primer3's PRIMER_MAX_POLY_X).</summary>
    public int MaxPolyX { get; init; } = 5;

    /// <summary>Fewest G/C among the 3′-terminal 5 nt (default 1).</summary>
    public int MinGcClamp { get; init; } = 1;

    /// <summary>Most G/C among the 3′-terminal 5 nt (default 3; more makes mispriming likelier).</summary>
    public int MaxGcClamp { get; init; } = 3;

    /// <summary>Most negative self-dimer ΔG°37 accepted, kcal/mol (default −9, the usual OligoAnalyzer threshold).</summary>
    public double MaxSelfDimerDeltaG { get; init; } = -9.0;

    /// <summary>Most negative 3′-anchored self- or hetero-dimer ΔG°37 accepted, kcal/mol (default −5).</summary>
    public double MaxThreePrimeDimerDeltaG { get; init; } = -5.0;

    /// <summary>Most negative hairpin ΔG°37 accepted, kcal/mol (default −3).</summary>
    public double MaxHairpinDeltaG { get; init; } = -3.0;

    /// <summary>Most negative forward/reverse hetero-dimer ΔG°37 accepted, kcal/mol (default −9).</summary>
    public double MaxHeteroDimerDeltaG { get; init; } = -9.0;

    /// <summary>Most negative 3′ end stability (ΔG°37 of the last 5 nt) accepted, kcal/mol (default −9).</summary>
    public double MaxEndStabilityDeltaG { get; init; } = -9.0;

    /// <summary>Shortest accepted product, bp, primers included (default 100).</summary>
    public int ProductSizeMin { get; init; } = 100;

    /// <summary>Longest accepted product, bp, primers included (default 300).</summary>
    public int ProductSizeMax { get; init; } = 300;

    /// <summary>Salt and primer concentrations for Tm (default <see cref="SaltConditions.StandardPcr"/>).</summary>
    public SaltConditions SaltConditions { get; init; } = SaltConditions.StandardPcr;

    /// <summary>Most pairs returned by <see cref="PrimerDesigner.DesignPairs"/> (default 5, as Primer3).</summary>
    public int MaxPairsReturned { get; init; } = 5;

    /// <summary>Most primers returned by <see cref="PrimerDesigner.DesignForward"/>/<see cref="PrimerDesigner.DesignReverse"/> (default 20).</summary>
    public int MaxPrimersReturned { get; init; } = 20;

    /// <summary>Best-penalty candidates kept per strand before pairing (default 1000); bounds pairing work.</summary>
    public int MaxCandidatesPerStrand { get; init; } = 1000;

    /// <summary>Region the product must contain: the forward primer must end at or before its start and the reverse primer begin at or after its end.</summary>
    public TemplateRegion? Target { get; init; }

    /// <summary>Region both primers (and hence the product) must lie in; null = the whole template.</summary>
    public TemplateRegion? IncludedRegion { get; init; }

    /// <summary>Penalty per °C of |Tm − OptTm| (default 1).</summary>
    public double TmWeight { get; init; } = 1.0;

    /// <summary>Penalty per nucleotide of |length − OptLength| (default 1).</summary>
    public double SizeWeight { get; init; } = 1.0;

    /// <summary>Penalty per percentage point of |GC − OptGc| (default 0.05).</summary>
    public double GcWeight { get; init; } = 0.05;

    /// <summary>Penalty per kcal/mol of self-dimer stability (default 0.1).</summary>
    public double SelfDimerWeight { get; init; } = 0.1;

    /// <summary>Penalty per kcal/mol of 3′-anchored dimer stability, self or hetero (default 0.2).</summary>
    public double ThreePrimeDimerWeight { get; init; } = 0.2;

    /// <summary>Penalty per kcal/mol of hairpin stability (default 0.2).</summary>
    public double HairpinWeight { get; init; } = 0.2;

    /// <summary>Pair penalty per °C of Tm difference (default 0.5).</summary>
    public double TmDifferenceWeight { get; init; } = 0.5;

    /// <summary>Pair penalty per kcal/mol of hetero-dimer stability (default 0.1).</summary>
    public double HeteroDimerWeight { get; init; } = 0.1;
}
