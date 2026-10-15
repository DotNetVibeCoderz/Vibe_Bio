using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Primers;

/// <summary>
/// A PCR primer with the properties used to choose it: melting temperature, composition, secondary structure
/// and its penalty score.
/// </summary>
/// <remarks>
/// <see cref="Start"/> is the 0-based position on the <b>plus</b> strand of the template of the leftmost base
/// covered by the primer, whichever strand the primer anneals to (the convention of
/// <see cref="GenomicInterval"/> and Primer3's left coordinate for both primers); the primer occupies
/// <c>[Start, Start + Length)</c>. <see cref="Sequence"/> is always written 5′→3′ as ordered, so for a
/// <see cref="Strand.Reverse"/> primer it is the reverse complement of the template at that location.
/// Primers produced by <see cref="PrimerDesigner.Evaluate(string, SaltConditions?)"/> have no template:
/// <see cref="Start"/> is −1 and <see cref="Strand"/> is <see cref="Strand.Unknown"/>.
/// All ΔG values are ΔG°37 in kcal/mol from <see cref="OligoThermodynamics"/> (negative = more stable; 0 = none).
/// </remarks>
/// <param name="Sequence">The primer, 5′→3′, upper case.</param>
/// <param name="Start">0-based leftmost plus-strand position on the template, or −1 without a template.</param>
/// <param name="Length">Length in nucleotides.</param>
/// <param name="Strand">Template strand the primer is identical to (Forward = left primer, Reverse = right primer).</param>
/// <param name="Tm">Nearest-neighbour melting temperature in °C under the design salt conditions.</param>
/// <param name="GcFraction">G + C fraction, 0–1.</param>
/// <param name="SelfDimerDeltaG">ΔG°37 of the most stable self-dimer anywhere along the primer.</param>
/// <param name="ThreePrimeSelfDimerDeltaG">ΔG°37 of the most stable self-dimer that includes a 3′ end.</param>
/// <param name="HairpinDeltaG">ΔG°37 of the most stable hairpin.</param>
/// <param name="EndStability">ΔG°37 of the duplex formed by the 3′-terminal 5 nt (the whole primer if shorter), stacks plus initiation.</param>
/// <param name="GcClamp">Number of G/C among the 3′-terminal 5 nt.</param>
/// <param name="LongestHomopolymer">Longest single-base run (poly-X).</param>
/// <param name="Penalty">Weighted penalty (lower is better); see <see cref="PrimerDesignOptions"/>.</param>
public sealed record Primer(
    string Sequence,
    int Start,
    int Length,
    Strand Strand,
    double Tm,
    double GcFraction,
    double SelfDimerDeltaG,
    double ThreePrimeSelfDimerDeltaG,
    double HairpinDeltaG,
    double EndStability,
    int GcClamp,
    int LongestHomopolymer,
    double Penalty)
{
    /// <summary>0-based exclusive end on the plus strand (<c>Start + Length</c>).</summary>
    public int End => Start + Length;

    /// <inheritdoc/>
    public override string ToString()
        => $"{Sequence} {Strand} {Start}..{End} Tm={Tm:F1} GC={GcFraction:P0} penalty={Penalty:F2}";
}

/// <summary>
/// A forward/reverse primer pair and the PCR product they amplify.
/// </summary>
/// <remarks>
/// The product spans <c>[Forward.Start, Reverse.End)</c> on the plus strand of the template and includes both
/// primer sites, so <see cref="ProductSize"/> = <c>Reverse.End − Forward.Start</c>.
/// </remarks>
/// <param name="Forward">The left (plus-strand) primer.</param>
/// <param name="Reverse">The right (minus-strand) primer.</param>
/// <param name="ProductSize">Amplicon length in base pairs, primers included.</param>
/// <param name="ProductSequence">Plus-strand sequence of the amplicon.</param>
/// <param name="TmDifference">|Forward.Tm − Reverse.Tm| in °C.</param>
/// <param name="HeteroDimerDeltaG">ΔG°37 of the most stable forward/reverse dimer anywhere.</param>
/// <param name="ThreePrimeHeteroDimerDeltaG">ΔG°37 of the most stable forward/reverse dimer including a 3′ end.</param>
/// <param name="Penalty">Pair penalty (lower is better); see <see cref="PrimerDesignOptions"/>.</param>
public sealed record PrimerPair(
    Primer Forward,
    Primer Reverse,
    int ProductSize,
    string ProductSequence,
    double TmDifference,
    double HeteroDimerDeltaG,
    double ThreePrimeHeteroDimerDeltaG,
    double Penalty)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"{Forward.Sequence} / {Reverse.Sequence} product={ProductSize} bp ΔTm={TmDifference:F1} penalty={Penalty:F2}";
}
