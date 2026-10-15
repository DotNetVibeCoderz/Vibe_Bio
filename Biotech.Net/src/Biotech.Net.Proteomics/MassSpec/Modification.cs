using System.Globalization;
using System.Text;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>Where a modification may sit.</summary>
public enum ModificationPosition
{
    /// <summary>On any matching residue.</summary>
    Anywhere,

    /// <summary>On the peptide N-terminus (optionally restricted to given residues).</summary>
    PeptideNTerm,

    /// <summary>On the peptide C-terminus (optionally restricted to given residues).</summary>
    PeptideCTerm,

    /// <summary>On the protein N-terminus only.</summary>
    ProteinNTerm,

    /// <summary>On the protein C-terminus only.</summary>
    ProteinCTerm,
}

/// <summary>A post-translational or chemical modification: a mass shift on given residues or a terminus.</summary>
/// <param name="Name">Short name, e.g. <c>Oxidation</c> (Unimod PSI-MS name where one exists).</param>
/// <param name="MonoisotopicDelta">Monoisotopic mass shift in Da.</param>
/// <param name="AverageDelta">Average mass shift in Da.</param>
/// <param name="Residues">Residues it applies to; empty means any residue (used for terminal modifications).</param>
/// <param name="Position">Positional constraint.</param>
public sealed record Modification(string Name, double MonoisotopicDelta, double AverageDelta, string Residues, ModificationPosition Position = ModificationPosition.Anywhere)
{
    /// <summary>True for N- or C-terminal modifications (which do not consume a residue site).</summary>
    public bool IsTerminal => Position != ModificationPosition.Anywhere;

    /// <summary>True when <paramref name="residue"/> is an allowed target.</summary>
    public bool AppliesTo(char residue) => Residues.Length == 0 || Residues.Contains(char.ToUpperInvariant(residue));

    /// <summary>Carbamidomethyl (iodoacetamide alkylation) on C, Unimod:4, +57.021464.</summary>
    public static Modification Carbamidomethyl { get; } = new("Carbamidomethyl", 57.021464, 57.0513, "C");

    /// <summary>Oxidation on M, Unimod:35, +15.994915.</summary>
    public static Modification Oxidation { get; } = new("Oxidation", 15.994915, 15.9994, "M");

    /// <summary>Phosphorylation on S/T/Y, Unimod:21, +79.966331.</summary>
    public static Modification Phospho { get; } = new("Phospho", 79.966331, 79.9799, "STY");

    /// <summary>Acetylation of the protein N-terminus, Unimod:1, +42.010565.</summary>
    public static Modification AcetylProteinNTerm { get; } = new("Acetyl", 42.010565, 42.0367, "", ModificationPosition.ProteinNTerm);

    /// <summary>Acetylation of lysine, Unimod:1, +42.010565.</summary>
    public static Modification AcetylK { get; } = new("Acetyl", 42.010565, 42.0367, "K");

    /// <summary>Deamidation of N/Q, Unimod:7, +0.984016.</summary>
    public static Modification Deamidation { get; } = new("Deamidated", 0.984016, 0.9848, "NQ");

    /// <summary>Methylation of K/R, Unimod:34, +14.015650.</summary>
    public static Modification Methyl { get; } = new("Methyl", 14.015650, 14.0266, "KR");

    /// <summary>Ubiquitin remnant (GlyGly) on K, Unimod:121, +114.042927.</summary>
    public static Modification GlyGly { get; } = new("GG", 114.042927, 114.1026, "K");

    /// <summary>Pyro-glutamate from N-terminal Q, Unimod:28, −17.026549.</summary>
    public static Modification PyroGluQ { get; } = new("Gln->pyro-Glu", -17.026549, -17.0305, "Q", ModificationPosition.PeptideNTerm);

    /// <summary>TMT 6/10/11-plex label on K, Unimod:737, +229.162932 (add <see cref="Tmt6plexNTerm"/> for the N-terminus).</summary>
    public static Modification Tmt6plex { get; } = new("TMT6plex", 229.162932, 229.2634, "K");

    /// <summary>TMT 6-plex label on the peptide N-terminus.</summary>
    public static Modification Tmt6plexNTerm { get; } = new("TMT6plex", 229.162932, 229.2634, "", ModificationPosition.PeptideNTerm);

    /// <inheritdoc/>
    public override string ToString() => $"{Name} ({(Residues.Length == 0 ? Position.ToString() : Residues)}, {MonoisotopicDelta:+0.000000;-0.000000})";
}

/// <summary>
/// A peptide with optional modifications: one optional modification per residue plus N- and C-terminal ones.
/// Masses are computed once at construction.
/// </summary>
public sealed class Peptide : IEquatable<Peptide>
{
    private readonly Modification?[] _residueMods;
    private readonly double[] _monoResidues;

    /// <summary>Creates a peptide.</summary>
    /// <param name="sequence">Residues (upper-cased internally).</param>
    /// <param name="residueModifications">Per-residue modification or null; length must equal the sequence length when given.</param>
    /// <param name="nTermModification">N-terminal modification.</param>
    /// <param name="cTermModification">C-terminal modification.</param>
    public Peptide(string sequence, IReadOnlyList<Modification?>? residueModifications = null, Modification? nTermModification = null, Modification? cTermModification = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(sequence);
        Sequence = sequence.ToUpperInvariant();
        if (residueModifications is not null && residueModifications.Count != Sequence.Length)
            throw new ArgumentException("One modification slot per residue is required.", nameof(residueModifications));
        _residueMods = new Modification?[Sequence.Length];
        _monoResidues = new double[Sequence.Length];
        NTermModification = nTermModification;
        CTermModification = cTermModification;
        double mono = MassConstants.WaterMonoisotopic, avg = MassConstants.WaterAverage;
        for (int i = 0; i < Sequence.Length; i++)
        {
            double rm = ResidueMasses.Monoisotopic(Sequence[i]);
            if (double.IsNaN(rm)) throw new ArgumentException($"Residue '{Sequence[i]}' has no defined mass.", nameof(sequence));
            var mod = residueModifications?[i];
            _residueMods[i] = mod;
            _monoResidues[i] = rm + (mod?.MonoisotopicDelta ?? 0);
            mono += _monoResidues[i];
            avg += ResidueMasses.Average(Sequence[i]) + (mod?.AverageDelta ?? 0);
        }
        if (nTermModification is not null) { mono += nTermModification.MonoisotopicDelta; avg += nTermModification.AverageDelta; }
        if (cTermModification is not null) { mono += cTermModification.MonoisotopicDelta; avg += cTermModification.AverageDelta; }
        MonoisotopicMass = mono;
        AverageMass = avg;
    }

    /// <summary>Unmodified residues.</summary>
    public string Sequence { get; }

    /// <summary>Number of residues.</summary>
    public int Length => Sequence.Length;

    /// <summary>Modification on residue <paramref name="index"/>, or null.</summary>
    public Modification? ModificationAt(int index) => _residueMods[index];

    /// <summary>N-terminal modification, or null.</summary>
    public Modification? NTermModification { get; }

    /// <summary>C-terminal modification, or null.</summary>
    public Modification? CTermModification { get; }

    /// <summary>Number of modifications carried.</summary>
    public int ModificationCount
    {
        get
        {
            int n = (NTermModification is null ? 0 : 1) + (CTermModification is null ? 0 : 1);
            foreach (var m in _residueMods) if (m is not null) n++;
            return n;
        }
    }

    /// <summary>Neutral monoisotopic mass (residues + modifications + water).</summary>
    public double MonoisotopicMass { get; }

    /// <summary>Neutral average mass.</summary>
    public double AverageMass { get; }

    /// <summary>Monoisotopic residue masses including residue modifications (terminal modifications excluded).</summary>
    public ReadOnlySpan<double> ResidueMassesMonoisotopic => _monoResidues;

    /// <summary>m/z of the [M + zH]^z+ ion.</summary>
    public double Mz(int charge, bool monoisotopic = true) => MassConstants.ToMz(monoisotopic ? MonoisotopicMass : AverageMass, charge);

    /// <summary>Mass-shift notation, e.g. <c>[+42.011]PEPT[+79.966]IDEM[+15.995]K</c>.</summary>
    public string ModifiedSequence
    {
        get
        {
            var ci = CultureInfo.InvariantCulture;
            var sb = new StringBuilder(Sequence.Length + 16);
            if (NTermModification is not null) sb.Append('[').Append(NTermModification.MonoisotopicDelta.ToString("+0.000;-0.000", ci)).Append(']');
            for (int i = 0; i < Sequence.Length; i++)
            {
                sb.Append(Sequence[i]);
                var m = _residueMods[i];
                if (m is not null) sb.Append('[').Append(m.MonoisotopicDelta.ToString("+0.000;-0.000", ci)).Append(']');
            }
            if (CTermModification is not null) sb.Append("-[").Append(CTermModification.MonoisotopicDelta.ToString("+0.000;-0.000", ci)).Append(']');
            return sb.ToString();
        }
    }

    /// <inheritdoc/>
    public bool Equals(Peptide? other) => other is not null && other.ModifiedSequence == ModifiedSequence;

    /// <inheritdoc/>
    public override bool Equals(object? obj) => obj is Peptide p && Equals(p);

    /// <inheritdoc/>
    public override int GetHashCode() => ModifiedSequence.GetHashCode(StringComparison.Ordinal);

    /// <inheritdoc/>
    public override string ToString() => ModifiedSequence;
}

/// <summary>Applies fixed modifications and enumerates variable-modification isoforms of a peptide.</summary>
public static class PeptideModifier
{
    /// <summary>
    /// Every isoform of <paramref name="sequence"/> with all applicable fixed modifications and 0 to
    /// <paramref name="maxVariableModifications"/> variable ones (at most one modification per residue or terminus;
    /// fixed modifications win their site). Enumeration stops after <paramref name="maxIsoforms"/> isoforms.
    /// </summary>
    public static IReadOnlyList<Peptide> Expand(string sequence, IReadOnlyList<Modification>? fixedModifications = null,
        IReadOnlyList<Modification>? variableModifications = null, int maxVariableModifications = 2,
        bool isProteinNTerm = false, bool isProteinCTerm = false, int maxIsoforms = 1024)
    {
        ArgumentException.ThrowIfNullOrEmpty(sequence);
        sequence = sequence.ToUpperInvariant();
        int n = sequence.Length;
        var residue = new Modification?[n];
        Modification? nTerm = null, cTerm = null;

        bool TerminalAllowed(Modification m, bool nSide)
        {
            char r = nSide ? sequence[0] : sequence[^1];
            if (!m.AppliesTo(r)) return false;
            return m.Position switch
            {
                ModificationPosition.PeptideNTerm => nSide,
                ModificationPosition.ProteinNTerm => nSide && isProteinNTerm,
                ModificationPosition.PeptideCTerm => !nSide,
                ModificationPosition.ProteinCTerm => !nSide && isProteinCTerm,
                _ => false,
            };
        }

        if (fixedModifications is not null)
            foreach (var m in fixedModifications)
            {
                if (m.IsTerminal)
                {
                    if (nTerm is null && TerminalAllowed(m, true)) nTerm = m;
                    else if (cTerm is null && TerminalAllowed(m, false)) cTerm = m;
                    continue;
                }
                for (int i = 0; i < n; i++) if (residue[i] is null && m.AppliesTo(sequence[i])) residue[i] = m;
            }

        // Variable sites: (slot, mod) with slot −1 = N-term, n = C-term.
        var sites = new List<(int Slot, Modification Mod)>();
        if (variableModifications is not null)
            foreach (var m in variableModifications)
            {
                if (m.IsTerminal)
                {
                    if (nTerm is null && TerminalAllowed(m, true)) sites.Add((-1, m));
                    if (cTerm is null && TerminalAllowed(m, false)) sites.Add((n, m));
                    continue;
                }
                for (int i = 0; i < n; i++) if (residue[i] is null && m.AppliesTo(sequence[i])) sites.Add((i, m));
            }
        sites.Sort((a, b) => a.Slot.CompareTo(b.Slot));

        var results = new List<Peptide>();
        var chosen = new List<(int Slot, Modification Mod)>();
        void Emit()
        {
            var r = (Modification?[])residue.Clone();
            Modification? nt = nTerm, ct = cTerm;
            foreach (var (slot, mod) in chosen)
            {
                if (slot == -1) nt = mod;
                else if (slot == n) ct = mod;
                else r[slot] = mod;
            }
            results.Add(new Peptide(sequence, r, nt, ct));
        }
        void Recurse(int from)
        {
            if (results.Count >= maxIsoforms) return;
            Emit();
            if (chosen.Count == maxVariableModifications) return;
            for (int s = from; s < sites.Count && results.Count < maxIsoforms; s++)
            {
                bool slotTaken = false;
                foreach (var c in chosen) if (c.Slot == sites[s].Slot) { slotTaken = true; break; }
                if (slotTaken) continue;
                chosen.Add(sites[s]);
                Recurse(s + 1);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Recurse(0);
        return results;
    }
}
