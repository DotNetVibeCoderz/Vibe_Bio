using Biotech.Net.Core;
using Biotech.Net.Native;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>Unit of a mass tolerance.</summary>
public enum ToleranceUnit
{
    /// <summary>Absolute, in daltons (or Th for m/z).</summary>
    Dalton,

    /// <summary>Relative, parts per million.</summary>
    Ppm,
}

/// <summary>A symmetric mass tolerance.</summary>
/// <param name="Value">Tolerance magnitude.</param>
/// <param name="Unit">Dalton or ppm.</param>
public readonly record struct Tolerance(double Value, ToleranceUnit Unit)
{
    /// <summary>A ppm tolerance.</summary>
    public static Tolerance Ppm(double value) => new(value, ToleranceUnit.Ppm);

    /// <summary>A dalton tolerance.</summary>
    public static Tolerance Da(double value) => new(value, ToleranceUnit.Dalton);

    /// <summary>Half-width of the window around <paramref name="mass"/>.</summary>
    public double HalfWidth(double mass) => Unit == ToleranceUnit.Ppm ? Math.Abs(mass) * Value * 1e-6 : Value;

    /// <summary>True when <paramref name="observed"/> lies within the tolerance of <paramref name="theoretical"/>.</summary>
    public bool Contains(double theoretical, double observed) => Math.Abs(observed - theoretical) <= HalfWidth(theoretical);

    /// <inheritdoc/>
    public override string ToString() => Unit == ToleranceUnit.Ppm ? $"±{Value} ppm" : $"±{Value} Da";
}

/// <summary>Database-search settings.</summary>
public sealed record SearchParameters
{
    /// <summary>Protease.</summary>
    public Enzyme Enzyme { get; init; } = Enzyme.Trypsin;

    /// <summary>Digestion settings.</summary>
    public DigestionOptions Digestion { get; init; } = new(MissedCleavages: 1, MinLength: 6, MaxLength: 40);

    /// <summary>Precursor tolerance (default ±10 ppm).</summary>
    public Tolerance PrecursorTolerance { get; init; } = Tolerance.Ppm(10);

    /// <summary>Fragment tolerance (default ±0.02 Da, high-resolution HCD).</summary>
    public Tolerance FragmentTolerance { get; init; } = Tolerance.Da(0.02);

    /// <summary>Fixed modifications (default carbamidomethyl C).</summary>
    public IReadOnlyList<Modification> FixedModifications { get; init; } = [Modification.Carbamidomethyl];

    /// <summary>Variable modifications (default oxidation M).</summary>
    public IReadOnlyList<Modification> VariableModifications { get; init; } = [Modification.Oxidation];

    /// <summary>Maximum variable modifications per peptide.</summary>
    public int MaxVariableModifications { get; init; } = 2;

    /// <summary>Fragment series to score.</summary>
    public IReadOnlyList<IonType> IonTypes { get; init; } = [IonType.B, IonType.Y];

    /// <summary>Highest fragment charge considered (capped at precursor charge − 1, minimum 1).</summary>
    public int MaxFragmentCharge { get; init; } = 2;

    /// <summary>Charges tried when a spectrum has no precursor charge.</summary>
    public IReadOnlyList<int> DefaultCharges { get; init; } = [2, 3];

    /// <summary>Append reversed-sequence decoys for target–decoy FDR.</summary>
    public bool GenerateDecoys { get; init; } = true;

    /// <summary>Prefix marking decoy protein accessions.</summary>
    public string DecoyPrefix { get; init; } = "DECOY_";

    /// <summary>Most-intense peaks kept per spectrum before scoring (0 = all).</summary>
    public int TopPeaks { get; init; } = 150;
}

/// <summary>A peptide–spectrum match.</summary>
public sealed class PeptideSpectrumMatch
{
    internal PeptideSpectrumMatch(Spectrum spectrum, Peptide peptide, IReadOnlyList<string> proteins, bool isDecoy, int charge,
        double score, int matchedN, int matchedC, int theoretical, double precursorErrorPpm)
    {
        Spectrum = spectrum; Peptide = peptide; ProteinIds = proteins; IsDecoy = isDecoy; Charge = charge;
        Score = score; MatchedNTerminalIons = matchedN; MatchedCTerminalIons = matchedC; TheoreticalIons = theoretical;
        PrecursorErrorPpm = precursorErrorPpm;
    }

    /// <summary>The spectrum.</summary>
    public Spectrum Spectrum { get; }

    /// <summary>The best-scoring peptide.</summary>
    public Peptide Peptide { get; }

    /// <summary>Every protein containing the peptide (targets, or decoys when <see cref="IsDecoy"/>).</summary>
    public IReadOnlyList<string> ProteinIds { get; }

    /// <summary>True when the peptide only occurs in decoy proteins.</summary>
    public bool IsDecoy { get; }

    /// <summary>Precursor charge used.</summary>
    public int Charge { get; }

    /// <summary>
    /// X!Tandem-style hyperscore: <c>ln(N_N! · N_C! · Σ I)</c> with intensities normalised to a base peak of 100, where
    /// N_N and N_C are matched N-terminal (a/b/c) and C-terminal (x/y/z) ions.
    /// </summary>
    public double Score { get; }

    /// <summary>Matched a/b/c ions.</summary>
    public int MatchedNTerminalIons { get; }

    /// <summary>Matched x/y/z ions.</summary>
    public int MatchedCTerminalIons { get; }

    /// <summary>Theoretical ions considered.</summary>
    public int TheoreticalIons { get; }

    /// <summary>Precursor mass error in ppm.</summary>
    public double PrecursorErrorPpm { get; }

    /// <summary>Target–decoy q-value (set by <see cref="TargetDecoy.AssignQValues"/>; NaN before).</summary>
    public double QValue { get; internal set; } = double.NaN;

    /// <inheritdoc/>
    public override string ToString() => $"{Spectrum.Id} → {Peptide} ({Charge}+) score {Score:F2}{(IsDecoy ? " [decoy]" : "")} q={QValue:F4}";
}

/// <summary>
/// A small, complete database search engine: in-silico digestion of targets and reversed decoys, fixed/variable
/// modification expansion, a mass-sorted peptide index with binary-searched precursor windows, fragment matching and
/// hyperscore ranking. Spectra are scored in parallel.
/// </summary>
public sealed class PeptideSpectrumMatcher
{
    private readonly SearchParameters _p;
    private readonly Peptide[] _peptides;
    private readonly double[] _masses;
    private readonly string[][] _proteins;
    private readonly bool[] _decoy;

    /// <summary>Builds the peptide index for <paramref name="proteins"/>.</summary>
    public PeptideSpectrumMatcher(IEnumerable<Protein> proteins, SearchParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(proteins);
        _p = parameters ?? new SearchParameters();
        var index = new Dictionary<string, (Peptide Peptide, List<string> Targets, List<string> Decoys)>(StringComparer.Ordinal);
        foreach (var protein in proteins)
        {
            var residues = protein.Sequence.Span;
            AddProtein(index, residues, protein.Id, decoy: false);
            if (_p.GenerateDecoys)
            {
                var rev = residues.ToArray();
                Array.Reverse(rev);
                AddProtein(index, rev, _p.DecoyPrefix + protein.Id, decoy: true);
            }
        }
        var entries = index.Values.ToArray();
        _masses = new double[entries.Length];
        for (int i = 0; i < entries.Length; i++) _masses[i] = entries[i].Peptide.MonoisotopicMass;
        var order = new int[entries.Length];
        for (int i = 0; i < order.Length; i++) order[i] = i;
        Array.Sort((double[])_masses.Clone(), order);
        _peptides = new Peptide[entries.Length];
        _proteins = new string[entries.Length][];
        _decoy = new bool[entries.Length];
        for (int i = 0; i < order.Length; i++)
        {
            var e = entries[order[i]];
            _peptides[i] = e.Peptide;
            _masses[i] = e.Peptide.MonoisotopicMass;
            // A peptide present in any target counts as a target (standard concatenated-search convention).
            _decoy[i] = e.Targets.Count == 0;
            _proteins[i] = (e.Targets.Count > 0 ? e.Targets : e.Decoys).Distinct().ToArray();
        }
    }

    /// <summary>Number of indexed (modified) peptide forms, targets and decoys.</summary>
    public int IndexedPeptides => _peptides.Length;

    private void AddProtein(Dictionary<string, (Peptide, List<string>, List<string>)> index, ReadOnlySpan<byte> residues, string id, bool decoy)
    {
        foreach (var d in Digestion.Digest(residues, _p.Enzyme, _p.Digestion, id))
        {
            bool hasUnknown = false;
            foreach (char c in d.Sequence) if (double.IsNaN(ResidueMasses.Monoisotopic(c))) { hasUnknown = true; break; }
            if (hasUnknown) continue;
            bool nTerm = d.Start <= 1, cTerm = d.End == residues.Length;
            foreach (var pep in PeptideModifier.Expand(d.Sequence, _p.FixedModifications, _p.VariableModifications, _p.MaxVariableModifications, nTerm, cTerm))
            {
                string key = pep.ModifiedSequence;
                if (!index.TryGetValue(key, out var entry)) index[key] = entry = (pep, [], []);
                (decoy ? entry.Item3 : entry.Item2).Add(id);
            }
        }
    }

    /// <summary>Best match per spectrum (spectra without any candidate are omitted). Results keep input order.</summary>
    public List<PeptideSpectrumMatch> Search(IEnumerable<Spectrum> spectra)
    {
        ArgumentNullException.ThrowIfNull(spectra);
        var list = spectra as IReadOnlyList<Spectrum> ?? spectra.ToList();
        var results = new PeptideSpectrumMatch?[list.Count];
        Parallel.For(0, list.Count, i => results[i] = Match(list[i]));
        var output = new List<PeptideSpectrumMatch>(list.Count);
        foreach (var r in results) if (r is not null) output.Add(r);
        return output;
    }

    /// <summary>Scores one spectrum against all candidates within the precursor window; null when none.</summary>
    public PeptideSpectrumMatch? Match(Spectrum spectrum)
    {
        ArgumentNullException.ThrowIfNull(spectrum);
        if (double.IsNaN(spectrum.PrecursorMz) || spectrum.Count == 0) return null;
        var s = _p.TopPeaks > 0 ? spectrum.TopPeaks(_p.TopPeaks) : spectrum;
        double maxI = 0;
        foreach (double x in s.Intensity) maxI = Math.Max(maxI, x);
        if (maxI <= 0) return null;

        IReadOnlyList<int> charges = spectrum.PrecursorCharge > 0 ? [spectrum.PrecursorCharge] : _p.DefaultCharges;
        PeptideSpectrumMatch? best = null;
        var ladder = new double[256];
        foreach (int z in charges)
        {
            double mass = MassConstants.ToNeutralMass(spectrum.PrecursorMz, z);
            double w = _p.PrecursorTolerance.HalfWidth(mass);
            int lo = LowerBound(_masses, mass - w);
            for (int k = lo; k < _masses.Length && _masses[k] <= mass + w; k++)
            {
                var pep = _peptides[k];
                if (ladder.Length < pep.Length) ladder = new double[pep.Length * 2];
                var (score, nN, nC, theo) = Score(pep, s, maxI, Math.Clamp(z - 1, 1, _p.MaxFragmentCharge), ladder);
                if (best is null || score > best.Score)
                    best = new PeptideSpectrumMatch(spectrum, pep, _proteins[k], _decoy[k], z, score, nN, nC, theo,
                        MassConstants.PpmError(mass, pep.MonoisotopicMass));
            }
        }
        return best;
    }

    private (double Score, int MatchedN, int MatchedC, int Theoretical) Score(Peptide pep, Spectrum s, double maxI, int maxCharge, double[] ladder)
    {
        int nN = 0, nC = 0, theo = 0;
        double sumI = 0;
        var mz = s.Mz;
        foreach (var type in _p.IonTypes)
        {
            bool nTerminal = type is IonType.A or IonType.B or IonType.C;
            for (int z = 1; z <= maxCharge; z++)
            {
                FragmentCalculator.Ladder(pep, type, z, ladder);
                for (int i = 0; i < pep.Length - 1; i++)
                {
                    theo++;
                    double t = ladder[i], tol = _p.FragmentTolerance.HalfWidth(t);
                    int j = LowerBound(mz, t - tol);
                    int bestPeak = -1;
                    for (; j < mz.Length && mz[j] <= t + tol; j++)
                        if (bestPeak < 0 || s.Intensity[j] > s.Intensity[bestPeak]) bestPeak = j;
                    if (bestPeak < 0) continue;
                    sumI += 100.0 * s.Intensity[bestPeak] / maxI;
                    if (nTerminal) nN++; else nC++;
                }
            }
        }
        if (nN + nC == 0) return (0, 0, 0, theo);
        double score = Math.Log(sumI) + SpecialFunctions.LogFactorial(nN) + SpecialFunctions.LogFactorial(nC);
        return (score, nN, nC, theo);
    }

    private static int LowerBound(double[] sorted, double value)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid] < value) lo = mid + 1; else hi = mid;
        }
        return lo;
    }
}

/// <summary>Target–decoy false-discovery-rate estimation (Elias &amp; Gygi, Nat. Methods 4:207, 2007).</summary>
public static class TargetDecoy
{
    /// <summary>
    /// Sorts by descending score and assigns q-values: FDR at each threshold is <c>#decoys / #targets</c> and the q-value is
    /// the minimum FDR at that score or any lower threshold. Ties share a q-value.
    /// </summary>
    public static void AssignQValues(IList<PeptideSpectrumMatch> psms)
    {
        ArgumentNullException.ThrowIfNull(psms);
        var sorted = psms.OrderByDescending(p => p.Score).ToArray();
        var fdr = new double[sorted.Length];
        int targets = 0, decoys = 0;
        for (int i = 0; i < sorted.Length;)
        {
            int j = i;
            while (j < sorted.Length && sorted[j].Score == sorted[i].Score)
            {
                if (sorted[j].IsDecoy) decoys++; else targets++;
                j++;
            }
            double f = targets == 0 ? 1 : Math.Min(1, (double)decoys / targets);
            for (int k = i; k < j; k++) fdr[k] = f;
            i = j;
        }
        double running = 1;
        for (int i = sorted.Length - 1; i >= 0; i--)
        {
            running = Math.Min(running, fdr[i]);
            sorted[i].QValue = running;
        }
    }

    /// <summary>Assigns q-values and returns the target PSMs at or below <paramref name="maxQValue"/>, best first.</summary>
    public static List<PeptideSpectrumMatch> Filter(IList<PeptideSpectrumMatch> psms, double maxQValue = 0.01)
    {
        AssignQValues(psms);
        var kept = new List<PeptideSpectrumMatch>();
        foreach (var p in psms) if (!p.IsDecoy && p.QValue <= maxQValue) kept.Add(p);
        kept.Sort((a, b) => b.Score.CompareTo(a.Score));
        return kept;
    }
}
