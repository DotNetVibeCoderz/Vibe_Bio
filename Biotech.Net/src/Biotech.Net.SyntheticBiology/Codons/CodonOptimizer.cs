using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Codons;

/// <summary>How a synonymous codon is chosen at each position of the design.</summary>
public enum CodonOptimizationMethod
{
    /// <summary>
    /// Always take the host's most used codon for the amino acid. This maximises the codon adaptation index
    /// (CAI = 1 by construction) but produces a sequence with no codon variety at all, which raises the risk of
    /// long repeats, extreme GC and tRNA pool depletion. It is the baseline every vendor's tool offers.
    /// </summary>
    MostFrequent,

    /// <summary>
    /// Sample a codon from the synonymous family in proportion to its usage in the host, so the finished gene
    /// reproduces the host's codon distribution rather than collapsing onto one codon per amino acid. Sampling
    /// is driven by a seeded generator, so a given seed always yields the same sequence.
    /// </summary>
    WeightedRandom,

    /// <summary>
    /// Codon harmonisation: for each position keep the <i>relative</i> usage the codon had in the source host,
    /// choosing the target codon whose relative adaptiveness is closest to it. Rare codons in the original gene
    /// stay rare, so the translational pauses that let a domain fold before the next is synthesised are
    /// preserved — see Angov, Hillier, Kincaid and Lyon, "Heterologous protein expression is enhanced by
    /// harmonizing the codon usage frequencies of the target gene with those of the expression host",
    /// PLoS ONE 3:e2189 (2008). Requires both a source table and an input coding sequence.
    /// </summary>
    Harmonised,
}

/// <summary>
/// Everything the optimiser needs besides the sequence itself: the host to optimise for, how to pick codons,
/// and the manufacturing constraints the result has to satisfy.
/// </summary>
public sealed class CodonOptimizationOptions
{
    /// <summary>The expression host to optimise for. Required.</summary>
    public required CodonUsageTable Target { get; init; }

    /// <summary>
    /// The organism the gene comes from, needed only by <see cref="CodonOptimizationMethod.Harmonised"/>, which
    /// compares each codon's standing in its native host with the candidates in the target host.
    /// </summary>
    public CodonUsageTable? Source { get; init; }

    /// <summary>How to choose among synonymous codons. <see cref="CodonOptimizationMethod.MostFrequent"/> by default.</summary>
    public CodonOptimizationMethod Method { get; init; } = CodonOptimizationMethod.MostFrequent;

    /// <summary>
    /// Seed of the generator used by <see cref="CodonOptimizationMethod.WeightedRandom"/>. The whole pipeline is
    /// deterministic: the same inputs with the same seed always produce the same sequence.
    /// </summary>
    public int Seed { get; init; } = 42;

    /// <summary>Constraints the design must satisfy; no repair is attempted when this is null.</summary>
    public SequenceConstraints? Constraints { get; init; }

    /// <summary>
    /// Genetic code to read the input with. When null — the default — the target table's own code is used, which
    /// is what keeps an <i>E. coli</i> design on the bacterial table and a mammalian design on the standard one.
    /// </summary>
    public GeneticCode? Code { get; init; }

    /// <summary>The code actually used: <see cref="Code"/> when given, otherwise the target table's code.</summary>
    public GeneticCode EffectiveCode => Code ?? Target.Code;

    /// <summary>
    /// Whether a terminal stop codon is copied through untouched. With the default (true) the stop codon of an
    /// input coding sequence is kept exactly as it was and is never altered by constraint repair; with false it
    /// is re-coded like any other position. No stop codon is ever added to a protein that does not end in
    /// <c>*</c> — the optimiser encodes the residues it is given and nothing else.
    /// </summary>
    public bool PreserveStop { get; init; } = true;

    /// <summary>
    /// Upper bound on the number of accepted repair substitutions. Repair stops early when the sequence is clean
    /// or when no single synonymous substitution improves it any further.
    /// </summary>
    public int MaxRepairIterations { get; init; } = 200;
}

/// <summary>The outcome of an optimisation run.</summary>
/// <param name="Dna">The designed coding sequence, upper case.</param>
/// <param name="Cai">Codon adaptation index of the design in the target host (see <see cref="CodonUsageTable.Cai(string)"/>).</param>
/// <param name="GcContent">Overall GC fraction (0–1) of the design.</param>
/// <param name="CodonsChanged">
/// How many codons differ from the starting point — the input coding sequence for
/// <see cref="CodonOptimizer.OptimizeCoding"/>, or the unrepaired initial design for
/// <see cref="CodonOptimizer.Optimize"/>, where it therefore counts the repair substitutions.
/// </param>
/// <param name="RemainingViolations">Constraint violations repair could not resolve; empty when the design is clean.</param>
/// <param name="RepairIterations">How many repair substitutions were accepted.</param>
public sealed record CodonOptimizationResult(
    string Dna,
    double Cai,
    double GcContent,
    int CodonsChanged,
    IReadOnlyList<ConstraintViolation> RemainingViolations,
    int RepairIterations);

/// <summary>
/// Designs a coding sequence for a chosen expression host: back-translates a protein, or re-codes an existing
/// gene, then repairs whatever manufacturing constraints the first draft happens to break.
/// </summary>
/// <remarks>
/// <para>
/// The design is a two-stage process. First every position is given a codon by the chosen
/// <see cref="CodonOptimizationMethod"/>. Then, if constraints were supplied, the sequence is scanned and each
/// violation is attacked by swapping a synonymous codon inside the offending window for the next best one, which
/// cannot change the encoded protein. Both stages are deterministic, so the same options reproduce the same gene.
/// </para>
/// <para>
/// The encoded protein is invariant throughout: every substitution the optimiser makes is synonymous under the
/// genetic code in use, which the tests assert by translating the result back.
/// </para>
/// </remarks>
public static class CodonOptimizer
{
    /// <summary>
    /// Back-translates a protein into a coding sequence for the target host and repairs any constraint
    /// violations in the draft.
    /// </summary>
    /// <param name="protein">One-letter amino acid sequence; a trailing <c>*</c> is encoded as a stop codon.</param>
    /// <param name="options">Host, method, seed and constraints.</param>
    /// <exception cref="ArgumentException">
    /// The protein is empty or carries a symbol the genetic code cannot encode (an ambiguity code such as
    /// <c>X</c>, <c>B</c> or <c>Z</c>, or a non-residue character).
    /// </exception>
    /// <exception cref="InvalidOperationException">
    /// <see cref="CodonOptimizationMethod.Harmonised"/> was requested, which needs the source gene and therefore
    /// <see cref="OptimizeCoding"/>.
    /// </exception>
    public static CodonOptimizationResult Optimize(string protein, CodonOptimizationOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(options.Target);
        var residues = ValidateProtein(protein, options.Target);

        if (options.Method == CodonOptimizationMethod.Harmonised)
            throw new InvalidOperationException(
                "Codon harmonisation compares the codons of an existing gene with the target host, so it needs a " +
                "coding sequence: use OptimizeCoding.");

        var codons = new string[residues.Length];
        var random = new Random(options.Seed);
        for (int i = 0; i < residues.Length; i++)
            codons[i] = Choose(options.Target, residues[i], options.Method, random, null, options.Source);

        var locked = LockSingletons(options.Target, residues, options.PreserveStop);
        var draft = string.Concat(codons);
        var (dna, iterations) = Repair(codons, residues, locked, options);
        int changed = CountCodonDifferences(draft, dna);
        return Summarise(dna, changed, iterations, options);
    }

    /// <summary>
    /// Re-codes an existing coding sequence for the target host. The encoded protein, including its stop codon,
    /// is preserved exactly; only synonymous positions move.
    /// </summary>
    /// <param name="cds">The source coding sequence; its length must be a multiple of three.</param>
    /// <param name="options">Host, method, seed and constraints. <see cref="CodonOptimizationMethod.Harmonised"/>
    /// additionally requires <see cref="CodonOptimizationOptions.Source"/>.</param>
    /// <exception cref="ArgumentException">The length is not a multiple of three, or a codon is ambiguous.</exception>
    /// <exception cref="InvalidOperationException">Harmonisation was requested without a source table.</exception>
    public static CodonOptimizationResult OptimizeCoding(string cds, CodonOptimizationOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(options.Target);
        ArgumentException.ThrowIfNullOrEmpty(cds);
        if (cds.Length % 3 != 0)
            throw new ArgumentException(
                $"A coding sequence must be a whole number of codons; got {cds.Length} bases.", nameof(cds));
        if (options.Method == CodonOptimizationMethod.Harmonised && options.Source is null)
            throw new InvalidOperationException(
                "Codon harmonisation needs the codon usage of the organism the gene comes from: set Source.");

        string source = DnaUtility.NormalizeDna(cds);
        var code = options.EffectiveCode;
        int count = source.Length / 3;
        var residues = new char[count];
        var original = new string[count];
        for (int i = 0; i < count; i++)
        {
            string codon = source.Substring(i * 3, 3);
            int k = GeneticCode.Index(codon.AsSpan());
            if (k < 0)
                throw new ArgumentException(
                    $"Codon {i + 1} ('{codon}') contains a symbol that is not A, C, G, T or U.", nameof(cds));
            original[i] = codon;
            residues[i] = code.AminoAcidAt(k);
        }

        var codons = new string[count];
        var random = new Random(options.Seed);
        for (int i = 0; i < count; i++)
        {
            bool terminalStop = residues[i] == '*' && i == count - 1 && options.PreserveStop;
            codons[i] = terminalStop
                ? original[i]
                : Choose(options.Target, residues[i], options.Method, random, original[i], options.Source);
        }

        var locked = LockSingletons(options.Target, residues, options.PreserveStop);
        var (dna, iterations) = Repair(codons, residues, locked, options);
        int changed = CountCodonDifferences(source, dna);
        return Summarise(dna, changed, iterations, options);
    }

    /// <summary>
    /// Back-translates a protein with no constraint repair — the raw draft, useful as a baseline and as the
    /// building block of the two <c>Optimize</c> entry points.
    /// </summary>
    /// <param name="protein">One-letter amino acid sequence; <c>*</c> is encoded as a stop codon.</param>
    /// <param name="table">The host whose usage decides the codons.</param>
    /// <param name="method">
    /// <see cref="CodonOptimizationMethod.MostFrequent"/> or <see cref="CodonOptimizationMethod.WeightedRandom"/>;
    /// harmonisation is not available here because it needs the original gene.
    /// </param>
    /// <param name="seed">Seed for weighted sampling.</param>
    /// <exception cref="ArgumentException">The protein is empty or carries an unencodable symbol.</exception>
    /// <exception cref="InvalidOperationException">Harmonisation was requested.</exception>
    public static string BackTranslate(string protein, CodonUsageTable table, CodonOptimizationMethod method, int seed)
    {
        ArgumentNullException.ThrowIfNull(table);
        var residues = ValidateProtein(protein, table);
        if (method == CodonOptimizationMethod.Harmonised)
            throw new InvalidOperationException("Harmonisation needs a source gene; use OptimizeCoding.");

        var random = new Random(seed);
        var codons = new string[residues.Length];
        for (int i = 0; i < residues.Length; i++)
            codons[i] = Choose(table, residues[i], method, random, null, null);
        return string.Concat(codons);
    }

    /// <summary>
    /// The codon composition of a coding sequence: each codon's share of all codons in it, which is the same
    /// quantity <see cref="CodonUsageTable.Frequency(string)"/> reports for a host. Comparing the two is the
    /// quickest way to see whether a design really matches its host.
    /// </summary>
    /// <param name="cds">A coding sequence whose length is a multiple of three; ambiguous codons are skipped.</param>
    /// <exception cref="ArgumentException">The length is not a multiple of three.</exception>
    public static IReadOnlyDictionary<string, double> CodonFrequencyProfile(string cds)
    {
        ArgumentNullException.ThrowIfNull(cds);
        if (cds.Length % 3 != 0)
            throw new ArgumentException(
                $"A coding sequence must be a whole number of codons; got {cds.Length} bases.", nameof(cds));

        var counts = new double[64];
        double total = 0;
        for (int i = 0; i + 3 <= cds.Length; i += 3)
        {
            int k = GeneticCode.Index(cds.AsSpan(i, 3));
            if (k < 0) continue;
            counts[k]++;
            total++;
        }
        var profile = new Dictionary<string, double>(64, StringComparer.Ordinal);
        for (int k = 0; k < 64; k++) profile[GeneticCode.CodonAt(k)] = total == 0 ? 0 : counts[k] / total;
        return profile;
    }

    private static char[] ValidateProtein(string protein, CodonUsageTable table)
    {
        ArgumentException.ThrowIfNullOrEmpty(protein);
        var residues = new char[protein.Length];
        for (int i = 0; i < protein.Length; i++)
        {
            char aa = char.ToUpperInvariant(protein[i]);
            if (table.Code.CodonsFor(aa).Count == 0)
                throw new ArgumentException(
                    $"Residue {i + 1} ('{protein[i]}') is not encoded by genetic code {table.Code.Id}; " +
                    "ambiguity codes such as X, B, Z and J cannot be back-translated.", nameof(protein));
            residues[i] = aa;
        }
        return residues;
    }

    private static string Choose(
        CodonUsageTable table,
        char aminoAcid,
        CodonOptimizationMethod method,
        Random random,
        string? sourceCodon,
        CodonUsageTable? sourceTable)
    {
        var family = table.SynonymousCodons(aminoAcid);
        if (family.Count == 1) return family[0];

        switch (method)
        {
            case CodonOptimizationMethod.MostFrequent:
                return family[0];

            case CodonOptimizationMethod.WeightedRandom:
            {
                // Sample proportionally to usage within the family. The draw happens for every position, even a
                // single-codon one above, so the stream of random numbers depends only on the protein and seed.
                double total = 0;
                for (int i = 0; i < family.Count; i++) total += table.PerThousand(family[i]);
                if (total <= 0) return family[0];
                double point = random.NextDouble() * total;
                double running = 0;
                for (int i = 0; i < family.Count; i++)
                {
                    running += table.PerThousand(family[i]);
                    if (point < running) return family[i];
                }
                return family[family.Count - 1];
            }

            case CodonOptimizationMethod.Harmonised:
            {
                // Angov harmonisation: keep the codon's standing within its family rather than its identity.
                double wanted = sourceTable is null || sourceCodon is null
                    ? 1.0
                    : sourceTable.RelativeAdaptiveness(sourceCodon);
                string best = family[0];
                double bestDistance = double.MaxValue;
                for (int i = 0; i < family.Count; i++)
                {
                    double distance = Math.Abs(table.RelativeAdaptiveness(family[i]) - wanted);
                    // Ties go to the more frequent codon, which is the earlier entry of the ordered family.
                    if (distance < bestDistance - 1e-12)
                    {
                        bestDistance = distance;
                        best = family[i];
                    }
                }
                return best;
            }

            default:
                return family[0];
        }
    }

    private static bool[] LockSingletons(CodonUsageTable table, ReadOnlySpan<char> residues, bool preserveStop)
    {
        var locked = new bool[residues.Length];
        for (int i = 0; i < residues.Length; i++)
        {
            if (table.SynonymousCodons(residues[i]).Count == 1) locked[i] = true;
            if (preserveStop && residues[i] == '*' && i == residues.Length - 1) locked[i] = true;
        }
        return locked;
    }

    /// <summary>
    /// Repairs constraint violations by synonymous substitution. Each round re-checks the whole sequence, walks
    /// the violations in positional order and tries the alternative codons of every unlocked position inside the
    /// offending window, best first; the first substitution that strictly reduces the violation count is kept and
    /// the round restarts. Nothing outside the offending window moves, and because only synonymous codons are
    /// offered the protein cannot change. The loop ends when the sequence is clean, when no single substitution
    /// helps any more, or when the iteration budget runs out.
    /// </summary>
    private static (string Dna, int Iterations) Repair(
        string[] codons, char[] residues, bool[] locked, CodonOptimizationOptions options)
    {
        string dna = string.Concat(codons);
        var constraints = options.Constraints;
        if (constraints is null) return (dna, 0);

        var violations = constraints.Check(dna);
        int iterations = 0;
        while (violations.Count > 0 && iterations < options.MaxRepairIterations)
        {
            bool improved = false;
            double score = Score(violations);

            for (int v = 0; v < violations.Count && !improved; v++)
            {
                var violation = violations[v];
                int firstCodon = Math.Max(0, violation.Start / 3);
                int lastCodon = Math.Min(codons.Length - 1, (violation.Start + Math.Max(1, violation.Length) - 1) / 3);

                for (int ci = firstCodon; ci <= lastCodon && !improved; ci++)
                {
                    if (locked[ci]) continue;
                    var family = options.Target.SynonymousCodons(residues[ci]);
                    string current = codons[ci];
                    for (int f = 0; f < family.Count; f++)
                    {
                        if (string.Equals(family[f], current, StringComparison.Ordinal)) continue;
                        codons[ci] = family[f];
                        string candidateDna = string.Concat(codons);
                        var candidate = constraints.Check(candidateDna);
                        if (Score(candidate) < score)
                        {
                            dna = candidateDna;
                            violations = candidate;
                            improved = true;
                            break;
                        }
                        codons[ci] = current;
                    }
                }
            }

            if (!improved) break;
            iterations++;
        }
        return (dna, iterations);
    }

    // Violation count first, total offending length second: a swap that shortens a stretch without removing it
    // still counts as progress, which keeps the loop from stalling on a wide GC window.
    private static double Score(IReadOnlyList<ConstraintViolation> violations)
    {
        double length = 0;
        for (int i = 0; i < violations.Count; i++) length += violations[i].Length;
        return violations.Count * 1_000_000.0 + length;
    }

    private static int CountCodonDifferences(string before, string after)
    {
        int n = Math.Min(before.Length, after.Length) / 3;
        int changed = 0;
        for (int i = 0; i < n; i++)
            if (!before.AsSpan(i * 3, 3).SequenceEqual(after.AsSpan(i * 3, 3))) changed++;
        return changed;
    }

    private static CodonOptimizationResult Summarise(
        string dna, int codonsChanged, int iterations, CodonOptimizationOptions options)
    {
        var remaining = options.Constraints?.Check(dna) ?? [];
        return new CodonOptimizationResult(
            dna,
            options.Target.Cai(dna),
            DnaUtility.GcFraction(dna),
            codonsChanged,
            remaining,
            iterations);
    }
}
