using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.SyntheticBiology.Codons;

/// <summary>
/// Codon usage of a host: how often each of the 64 codons occurs per thousand codons, together with the
/// derived relative adaptiveness weights that the codon adaptation index (CAI) is built from.
/// </summary>
/// <remarks>
/// <para>
/// Synonymous codons are not used equally. Highly expressed genes in a given organism are biased towards the
/// codons read by its abundant tRNAs, and matching a transgene to that bias is the usual first step of gene
/// synthesis. The two statistics exposed here are the classical ones:
/// </para>
/// <list type="bullet">
///   <item><description><see cref="Frequency(string)"/> — the codon's share of all codons in the reference
///   set (the published per-thousand value divided by 1000).</description></item>
///   <item><description><see cref="RelativeAdaptiveness(string)"/> — <c>w = f(codon) / max f(synonymous)</c>,
///   the weight defined by Sharp and Li, "The codon adaptation index — a measure of directional synonymous
///   codon usage bias, and its potential applications", Nucleic Acids Research 15:1281–1295 (1987).</description></item>
/// </list>
/// <para>
/// The built-in tables carry the per-thousand values published by the Kazusa Codon Usage Database
/// (Nakamura, Gojobori and Ikemura, "Codon usage tabulated from international DNA sequence databases: status
/// for the year 2000", Nucleic Acids Research 28:292, 2000), rounded to the one decimal place the database
/// prints and rescaled so that each table sums to exactly 1000 per thousand (the published rounding leaves a
/// residue of a few parts per thousand). They describe whole-genome average usage; a project that has its own
/// expression data — for instance the codon usage of the most abundant proteins of its production strain —
/// should build its own table with <see cref="FromCodingSequences"/> or <see cref="FromCounts"/> rather than
/// rely on a genome average.
/// </para>
/// </remarks>
public sealed class CodonUsageTable
{
    private readonly double[] _perThousand;
    private readonly double[] _weights;
    private readonly Dictionary<char, string[]> _synonyms;

    /// <summary>Floor applied to a zero-frequency codon when taking logarithms in <see cref="Cai(string)"/>.</summary>
    /// <remarks>
    /// A codon that never occurs in the reference set would give <c>log 0</c> and collapse the geometric mean to
    /// zero. Sharp and Li's original treatment assigns such codons a small non-zero weight; 0.001 is the value
    /// used by most implementations (it corresponds to a codon some thousand-fold rarer than the family optimum).
    /// </remarks>
    public const double ZeroWeightFloor = 0.001;

    private CodonUsageTable(string name, string organism, int taxonomyId, GeneticCode code, double[] perThousand)
    {
        Name = name;
        Organism = organism;
        TaxonomyId = taxonomyId;
        Code = code;
        _perThousand = perThousand;
        _weights = new double[64];
        _synonyms = new Dictionary<char, string[]>(24);

        // Group the 64 codons by the amino acid they encode, order each family by decreasing usage and divide
        // by the family maximum to obtain the CAI weights.
        var families = new Dictionary<char, List<int>>(24);
        for (int k = 0; k < 64; k++)
        {
            char aa = code.AminoAcidAt(k);
            if (!families.TryGetValue(aa, out var list)) families[aa] = list = new List<int>(6);
            list.Add(k);
        }

        foreach (var (aa, members) in families)
        {
            double max = 0;
            for (int i = 0; i < members.Count; i++) max = Math.Max(max, _perThousand[members[i]]);
            for (int i = 0; i < members.Count; i++)
                _weights[members[i]] = max > 0 ? _perThousand[members[i]] / max : 0;

            // Stable ordering: decreasing frequency, then kernel index, so ties resolve the same way on every run.
            members.Sort((x, y) =>
            {
                int byFrequency = _perThousand[y].CompareTo(_perThousand[x]);
                return byFrequency != 0 ? byFrequency : x.CompareTo(y);
            });
            var codons = new string[members.Count];
            for (int i = 0; i < members.Count; i++) codons[i] = GeneticCode.CodonAt(members[i]);
            _synonyms[aa] = codons;
        }

        double gc = 0;
        for (int k = 0; k < 64; k++)
        {
            string codon = GeneticCode.CodonAt(k);
            int gcBases = 0;
            for (int i = 0; i < 3; i++) if (codon[i] is 'G' or 'C') gcBases++;
            gc += _perThousand[k] / 1000.0 * gcBases / 3.0;
        }
        GcContent = gc;
    }

    /// <summary>Short identifier of the table, for example <c>"E. coli K-12"</c>.</summary>
    public string Name { get; }

    /// <summary>Full organism (and strain) the usage was tabulated from.</summary>
    public string Organism { get; }

    /// <summary>NCBI taxonomy identifier of the organism, or 0 when the table was built from user data.</summary>
    public int TaxonomyId { get; }

    /// <summary>The genetic code the table is interpreted with; it decides which codons are synonymous.</summary>
    public GeneticCode Code { get; }

    /// <summary>
    /// Expected GC fraction (0–1) of a coding sequence written with this usage, that is the usage-weighted mean
    /// GC content of the 64 codons. A useful sanity check: an optimised gene should land near this value.
    /// </summary>
    public double GcContent { get; }

    /// <summary>The raw published usage of a codon, in occurrences per thousand codons.</summary>
    /// <param name="codon">A three-base DNA or RNA codon, any case.</param>
    /// <exception cref="ArgumentException">The codon is not three unambiguous bases.</exception>
    public double PerThousand(string codon) => _perThousand[RequireIndex(codon)];

    /// <summary>The raw published usage at a kernel codon index (16·b1 + 4·b2 + b3), per thousand codons.</summary>
    public double PerThousandAt(int index) => _perThousand[index];

    /// <summary>
    /// The fraction of all codons that this codon accounts for (the per-thousand value divided by 1000), so the
    /// 64 values sum to one.
    /// </summary>
    /// <exception cref="ArgumentException">The codon is not three unambiguous bases.</exception>
    public double Frequency(string codon) => _perThousand[RequireIndex(codon)] / 1000.0;

    /// <summary>
    /// The relative adaptiveness <c>w</c> of a codon: its frequency divided by the frequency of the most common
    /// codon for the same amino acid (Sharp and Li 1987). The family optimum has <c>w = 1</c> by construction and
    /// every other synonymous codon lies in (0, 1].
    /// </summary>
    /// <exception cref="ArgumentException">The codon is not three unambiguous bases.</exception>
    public double RelativeAdaptiveness(string codon) => _weights[RequireIndex(codon)];

    /// <summary>The relative adaptiveness at a kernel codon index.</summary>
    public double RelativeAdaptivenessAt(int index) => _weights[index];

    /// <summary>
    /// The codons encoding <paramref name="aminoAcid"/> (<c>*</c> for stops) in decreasing order of usage, so the
    /// first entry is the host's preferred codon and the last its rarest. Ties are broken by codon index, which
    /// keeps the order deterministic across runs.
    /// </summary>
    /// <exception cref="ArgumentException">The symbol encodes nothing under <see cref="Code"/>.</exception>
    public IReadOnlyList<string> SynonymousCodons(char aminoAcid)
    {
        char aa = char.ToUpperInvariant(aminoAcid);
        if (!_synonyms.TryGetValue(aa, out var codons))
            throw new ArgumentException($"'{aminoAcid}' is not encoded by genetic code {Code.Id}.", nameof(aminoAcid));
        return codons;
    }

    /// <summary>The host's most used codon for an amino acid — what <c>MostFrequent</c> optimisation picks.</summary>
    /// <exception cref="ArgumentException">The symbol encodes nothing under <see cref="Code"/>.</exception>
    public string MostFrequentCodon(char aminoAcid) => SynonymousCodons(aminoAcid)[0];

    /// <summary>The host's rarest codon for an amino acid, useful for building deliberately poor controls.</summary>
    /// <exception cref="ArgumentException">The symbol encodes nothing under <see cref="Code"/>.</exception>
    public string RarestCodon(char aminoAcid)
    {
        var codons = SynonymousCodons(aminoAcid);
        return codons[codons.Count - 1];
    }

    /// <summary>
    /// Codon adaptation index of a coding sequence: the geometric mean of the relative adaptiveness of its codons,
    /// a number in (0, 1] where 1 means every codon is its family's optimum in this host (Sharp and Li 1987).
    /// </summary>
    /// <remarks>
    /// Following the original definition, stop codons and the single-codon amino acids methionine and tryptophan
    /// are excluded: they have no synonymous alternative, so their weight is always 1 and including them would
    /// merely dilute the index towards 1 in proportion to how many of them the protein happens to contain.
    /// Codons containing an ambiguity code are skipped as well. The result is 0 when nothing is left to score.
    /// </remarks>
    /// <param name="cds">A coding sequence whose length is a multiple of three.</param>
    /// <exception cref="ArgumentException">The length is not a multiple of three.</exception>
    public double Cai(string cds)
    {
        ArgumentNullException.ThrowIfNull(cds);
        if (cds.Length % 3 != 0)
            throw new ArgumentException($"A coding sequence must be a whole number of codons; got {cds.Length} bases.", nameof(cds));

        double logSum = 0;
        int counted = 0;
        for (int i = 0; i + 3 <= cds.Length; i += 3)
        {
            int k = GeneticCode.Index(cds.AsSpan(i, 3));
            if (k < 0) continue;
            char aa = Code.AminoAcidAt(k);
            if (aa is '*' or 'M' or 'W') continue;
            double w = _weights[k];
            logSum += Math.Log(w > 0 ? w : ZeroWeightFloor);
            counted++;
        }
        return counted == 0 ? 0 : Math.Exp(logSum / counted);
    }

    /// <summary>Codon adaptation index of a <see cref="Sequence"/> (see <see cref="Cai(string)"/>).</summary>
    public double Cai(Sequence cds)
    {
        ArgumentNullException.ThrowIfNull(cds);
        return Cai(cds.ToString());
    }

    /// <summary>
    /// GC content of the third codon positions (GC3) of a coding sequence, 0–1. GC3 absorbs most of the
    /// synonymous freedom in the code and correlates with expression level in many hosts, so it is the usual
    /// quick read-out of how far a design has been pushed towards a GC-rich or AT-rich host.
    /// </summary>
    /// <param name="cds">A coding sequence whose length is a multiple of three.</param>
    /// <exception cref="ArgumentException">The length is not a multiple of three.</exception>
    public double Gc3Content(string cds)
    {
        ArgumentNullException.ThrowIfNull(cds);
        if (cds.Length % 3 != 0)
            throw new ArgumentException($"A coding sequence must be a whole number of codons; got {cds.Length} bases.", nameof(cds));

        int gc = 0, counted = 0;
        for (int i = 2; i < cds.Length; i += 3)
        {
            char third = char.ToUpperInvariant(cds[i]);
            if (third is 'G' or 'C') gc++;
            counted++;
        }
        return counted == 0 ? 0 : (double)gc / counted;
    }

    /// <summary>
    /// Root-mean-square difference between the relative adaptiveness vectors of two tables (0 when identical,
    /// larger the more the two hosts disagree about which synonymous codon to prefer). It is a convenient way to
    /// ask whether a gene optimised for one host is already acceptable in another, and because it compares
    /// <c>w</c> rather than raw frequencies it is insensitive to differences in amino acid composition.
    /// </summary>
    public double CodonAdaptationDistance(CodonUsageTable other)
    {
        ArgumentNullException.ThrowIfNull(other);
        double sum = 0;
        for (int k = 0; k < 64; k++)
        {
            double d = _weights[k] - other._weights[k];
            sum += d * d;
        }
        return Math.Sqrt(sum / 64.0);
    }

    /// <summary>Human-readable identity of the table.</summary>
    public override string ToString() => $"{Name} ({Organism}, GC {GcContent:P1})";

    /// <summary>
    /// Builds a table from published or measured per-thousand usage. Codons missing from the dictionary are
    /// taken as zero; values need not sum to exactly 1000 (nothing here depends on the absolute scale).
    /// </summary>
    /// <param name="name">Short identifier for the table.</param>
    /// <param name="organism">Organism and strain the data describe.</param>
    /// <param name="perThousand">Codon (three DNA bases, any case) to occurrences per thousand codons.</param>
    /// <param name="code">Genetic code deciding synonymy; the standard code by default.</param>
    /// <param name="taxonomyId">Optional NCBI taxonomy identifier.</param>
    /// <exception cref="ArgumentException">A key is not a codon, a value is negative, or a key repeats.</exception>
    public static CodonUsageTable FromCounts(
        string name,
        string organism,
        IReadOnlyDictionary<string, double> perThousand,
        GeneticCode? code = null,
        int taxonomyId = 0)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        ArgumentNullException.ThrowIfNull(organism);
        ArgumentNullException.ThrowIfNull(perThousand);

        var values = new double[64];
        var seen = new bool[64];
        foreach (var (codon, value) in perThousand)
        {
            if (codon is not { Length: 3 })
                throw new ArgumentException($"'{codon}' is not a three-base codon.", nameof(perThousand));
            int k = GeneticCode.Index(codon.AsSpan());
            if (k < 0)
                throw new ArgumentException($"'{codon}' contains a symbol that is not A, C, G, T or U.", nameof(perThousand));
            if (value < 0 || double.IsNaN(value))
                throw new ArgumentException($"Usage of '{codon}' must be a non-negative number; got {value}.", nameof(perThousand));
            if (seen[k])
                throw new ArgumentException($"Codon '{codon}' appears twice (DNA and RNA spellings collide).", nameof(perThousand));
            seen[k] = true;
            values[k] = value;
        }
        return new CodonUsageTable(name, organism, taxonomyId, code ?? GeneticCode.Standard, values);
    }

    /// <summary>
    /// Counts codons across real coding sequences and converts the counts to per-thousand usage. This is how a
    /// host table should be built when the project has its own data: feed in the coding sequences of the
    /// organism, or better of its most highly expressed genes, which is the reference set Sharp and Li used.
    /// </summary>
    /// <param name="name">Short identifier for the resulting table (also used as the organism label).</param>
    /// <param name="cdsList">Coding sequences; each must be a whole number of codons. Ambiguous codons are skipped.</param>
    /// <param name="code">Genetic code deciding synonymy; the standard code by default.</param>
    /// <exception cref="ArgumentException">A sequence length is not a multiple of three, or no codon was counted.</exception>
    public static CodonUsageTable FromCodingSequences(string name, IEnumerable<string> cdsList, GeneticCode? code = null)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        ArgumentNullException.ThrowIfNull(cdsList);

        var counts = new double[64];
        double total = 0;
        foreach (var cds in cdsList)
        {
            if (cds is null) throw new ArgumentException("A coding sequence was null.", nameof(cdsList));
            if (cds.Length % 3 != 0)
                throw new ArgumentException($"Coding sequence of {cds.Length} bases is not a whole number of codons.", nameof(cdsList));
            for (int i = 0; i + 3 <= cds.Length; i += 3)
            {
                int k = GeneticCode.Index(cds.AsSpan(i, 3));
                if (k < 0) continue;
                counts[k]++;
                total++;
            }
        }
        if (total == 0) throw new ArgumentException("No unambiguous codon was found in the input.", nameof(cdsList));
        for (int k = 0; k < 64; k++) counts[k] = counts[k] / total * 1000.0;
        return new CodonUsageTable(name, name, 0, code ?? GeneticCode.Standard, counts);
    }

    private static int RequireIndex(string codon)
    {
        ArgumentNullException.ThrowIfNull(codon);
        if (codon.Length != 3)
            throw new ArgumentException($"A codon has three bases; got '{codon}'.", nameof(codon));
        int k = GeneticCode.Index(codon.AsSpan());
        if (k < 0)
            throw new ArgumentException($"'{codon}' contains a symbol that is not A, C, G, T or U.", nameof(codon));
        return k;
    }

    /// <summary>
    /// <i>Escherichia coli</i> K-12, the classic expression host and the organism whose bias Sharp and Li
    /// measured. Note the extreme preferences: CTG takes more than half of all leucine, AAA dominates lysine,
    /// GAA dominates glutamate, arginine runs through CGT and CGC while AGA, AGG and CGA together stay near
    /// 7 per thousand — these rare arginine codons, with ATA, CTA and CCC, are the ones that stall heterologous
    /// expression and that optimisation exists to remove.
    /// </summary>
    public static readonly CodonUsageTable EColiK12 = FromCounts("E. coli K-12", "Escherichia coli K-12", new Dictionary<string, double>(64)
    {
        ["TTT"] = 22.4, ["TTC"] = 16.6, ["TTA"] = 13.9, ["TTG"] = 13.7,
        ["CTT"] = 11.0, ["CTC"] = 11.0, ["CTA"] = 3.9, ["CTG"] = 53.1,
        ["ATT"] = 30.3, ["ATC"] = 25.0, ["ATA"] = 4.4, ["ATG"] = 27.9,
        ["GTT"] = 18.3, ["GTC"] = 15.3, ["GTA"] = 10.9, ["GTG"] = 26.4,
        ["TCT"] = 8.5, ["TCC"] = 8.6, ["TCA"] = 7.2, ["TCG"] = 8.9,
        ["CCT"] = 7.0, ["CCC"] = 5.5, ["CCA"] = 8.4, ["CCG"] = 23.2,
        ["ACT"] = 8.9, ["ACC"] = 23.4, ["ACA"] = 7.1, ["ACG"] = 14.4,
        ["GCT"] = 15.3, ["GCC"] = 25.5, ["GCA"] = 20.2, ["GCG"] = 33.6,
        ["TAT"] = 16.2, ["TAC"] = 12.2, ["TAA"] = 2.0, ["TAG"] = 0.2,
        ["CAT"] = 12.9, ["CAC"] = 9.6, ["CAA"] = 15.3, ["CAG"] = 28.8,
        ["AAT"] = 17.7, ["AAC"] = 21.6, ["AAA"] = 33.6, ["AAG"] = 10.3,
        ["GAT"] = 32.2, ["GAC"] = 19.1, ["GAA"] = 39.5, ["GAG"] = 17.8,
        ["TGT"] = 5.2, ["TGC"] = 6.4, ["TGA"] = 0.9, ["TGG"] = 15.2,
        ["CGT"] = 20.9, ["CGC"] = 22.0, ["CGA"] = 3.6, ["CGG"] = 5.4,
        ["AGT"] = 8.8, ["AGC"] = 16.1, ["AGA"] = 2.1, ["AGG"] = 1.2,
        ["GGT"] = 24.7, ["GGC"] = 29.6, ["GGA"] = 8.0, ["GGG"] = 11.1,
    }, GeneticCode.Bacterial, taxonomyId: 83333);

    /// <summary>
    /// <i>Saccharomyces cerevisiae</i>, baker's yeast. Strongly AT-biased at third positions: TTG for leucine,
    /// AGA for arginine (the opposite of <i>E. coli</i>), GAA for glutamate, AAG for lysine and GCT for alanine,
    /// with the CG-containing arginine codons CGC and CGG almost unused.
    /// </summary>
    public static readonly CodonUsageTable SaccharomycesCerevisiae = FromCounts(
        "S. cerevisiae", "Saccharomyces cerevisiae", new Dictionary<string, double>(64)
    {
        ["TTT"] = 26.1, ["TTC"] = 18.4, ["TTA"] = 26.2, ["TTG"] = 27.2,
        ["CTT"] = 12.3, ["CTC"] = 5.4, ["CTA"] = 13.4, ["CTG"] = 10.5,
        ["ATT"] = 30.1, ["ATC"] = 17.2, ["ATA"] = 17.8, ["ATG"] = 20.9,
        ["GTT"] = 22.1, ["GTC"] = 11.8, ["GTA"] = 11.8, ["GTG"] = 10.8,
        ["TCT"] = 23.5, ["TCC"] = 14.2, ["TCA"] = 18.7, ["TCG"] = 8.6,
        ["CCT"] = 13.5, ["CCC"] = 6.8, ["CCA"] = 18.3, ["CCG"] = 5.3,
        ["ACT"] = 20.3, ["ACC"] = 12.7, ["ACA"] = 17.8, ["ACG"] = 8.0,
        ["GCT"] = 21.2, ["GCC"] = 12.6, ["GCA"] = 16.2, ["GCG"] = 6.2,
        ["TAT"] = 18.8, ["TAC"] = 14.8, ["TAA"] = 1.1, ["TAG"] = 0.5,
        ["CAT"] = 13.6, ["CAC"] = 7.8, ["CAA"] = 27.3, ["CAG"] = 12.1,
        ["AAT"] = 35.7, ["AAC"] = 24.8, ["AAA"] = 41.9, ["AAG"] = 30.8,
        ["GAT"] = 37.6, ["GAC"] = 20.2, ["GAA"] = 45.1, ["GAG"] = 19.2,
        ["TGT"] = 8.1, ["TGC"] = 4.8, ["TGA"] = 0.7, ["TGG"] = 10.4,
        ["CGT"] = 6.4, ["CGC"] = 2.6, ["CGA"] = 3.0, ["CGG"] = 1.7,
        ["AGT"] = 14.2, ["AGC"] = 9.8, ["AGA"] = 21.3, ["AGG"] = 9.2,
        ["GGT"] = 23.9, ["GGC"] = 9.8, ["GGA"] = 10.9, ["GGG"] = 6.0,
    }, GeneticCode.Standard, taxonomyId: 4932);

    /// <summary>
    /// <i>Homo sapiens</i>, the reference for mammalian expression and for synthetic genes destined for human
    /// cells. GC-rich at third positions: CTG for leucine, GAG above GAA, AAG above AAA, CAG above CAA, GCC for
    /// alanine, and arginine spread thinly over AGA, AGG, CGG and CGC.
    /// </summary>
    public static readonly CodonUsageTable HomoSapiens = FromCounts("H. sapiens", "Homo sapiens", new Dictionary<string, double>(64)
    {
        ["TTT"] = 17.6, ["TTC"] = 20.3, ["TTA"] = 7.7, ["TTG"] = 12.9,
        ["CTT"] = 13.2, ["CTC"] = 19.6, ["CTA"] = 7.2, ["CTG"] = 39.4,
        ["ATT"] = 16.0, ["ATC"] = 20.8, ["ATA"] = 7.5, ["ATG"] = 22.0,
        ["GTT"] = 11.0, ["GTC"] = 14.5, ["GTA"] = 7.1, ["GTG"] = 28.1,
        ["TCT"] = 15.2, ["TCC"] = 17.7, ["TCA"] = 12.2, ["TCG"] = 4.4,
        ["CCT"] = 17.5, ["CCC"] = 19.8, ["CCA"] = 16.9, ["CCG"] = 6.9,
        ["ACT"] = 13.1, ["ACC"] = 18.9, ["ACA"] = 15.1, ["ACG"] = 6.1,
        ["GCT"] = 18.4, ["GCC"] = 27.7, ["GCA"] = 15.8, ["GCG"] = 7.4,
        ["TAT"] = 12.2, ["TAC"] = 15.3, ["TAA"] = 1.0, ["TAG"] = 0.8,
        ["CAT"] = 10.9, ["CAC"] = 15.1, ["CAA"] = 12.3, ["CAG"] = 34.2,
        ["AAT"] = 17.0, ["AAC"] = 19.1, ["AAA"] = 24.4, ["AAG"] = 31.9,
        ["GAT"] = 21.8, ["GAC"] = 25.1, ["GAA"] = 29.0, ["GAG"] = 39.6,
        ["TGT"] = 10.6, ["TGC"] = 12.6, ["TGA"] = 1.6, ["TGG"] = 13.2,
        ["CGT"] = 4.5, ["CGC"] = 10.4, ["CGA"] = 6.2, ["CGG"] = 11.4,
        ["AGT"] = 12.1, ["AGC"] = 19.5, ["AGA"] = 12.2, ["AGG"] = 12.0,
        ["GGT"] = 10.8, ["GGC"] = 22.2, ["GGA"] = 16.5, ["GGG"] = 16.5,
    }, GeneticCode.Standard, taxonomyId: 9606);

    /// <summary>
    /// <i>Cricetulus griseus</i>, the Chinese hamster behind the CHO cell lines used for most therapeutic
    /// antibodies. Its bias closely tracks the human one — which is why genes optimised for human cells
    /// usually express well in CHO — but it is worth optimising against the host actually used.
    /// </summary>
    public static readonly CodonUsageTable CricetulusGriseus = FromCounts(
        "CHO (C. griseus)", "Cricetulus griseus", new Dictionary<string, double>(64)
    {
        ["TTT"] = 17.1, ["TTC"] = 21.1, ["TTA"] = 6.4, ["TTG"] = 13.1,
        ["CTT"] = 13.1, ["CTC"] = 19.5, ["CTA"] = 7.6, ["CTG"] = 38.9,
        ["ATT"] = 15.9, ["ATC"] = 21.7, ["ATA"] = 6.8, ["ATG"] = 22.4,
        ["GTT"] = 10.6, ["GTC"] = 15.1, ["GTA"] = 7.0, ["GTG"] = 28.6,
        ["TCT"] = 15.6, ["TCC"] = 18.5, ["TCA"] = 11.6, ["TCG"] = 4.4,
        ["CCT"] = 18.2, ["CCC"] = 19.7, ["CCA"] = 17.2, ["CCG"] = 6.2,
        ["ACT"] = 13.2, ["ACC"] = 19.6, ["ACA"] = 15.6, ["ACG"] = 6.1,
        ["GCT"] = 19.2, ["GCC"] = 27.3, ["GCA"] = 15.7, ["GCG"] = 6.4,
        ["TAT"] = 12.2, ["TAC"] = 15.9, ["TAA"] = 0.9, ["TAG"] = 0.6,
        ["CAT"] = 10.5, ["CAC"] = 15.1, ["CAA"] = 12.1, ["CAG"] = 34.3,
        ["AAT"] = 16.7, ["AAC"] = 20.1, ["AAA"] = 23.2, ["AAG"] = 33.8,
        ["GAT"] = 21.7, ["GAC"] = 26.0, ["GAA"] = 27.8, ["GAG"] = 39.9,
        ["TGT"] = 10.6, ["TGC"] = 12.6, ["TGA"] = 1.3, ["TGG"] = 12.3,
        ["CGT"] = 4.6, ["CGC"] = 9.9, ["CGA"] = 6.5, ["CGG"] = 10.6,
        ["AGT"] = 12.1, ["AGC"] = 19.7, ["AGA"] = 12.2, ["AGG"] = 11.6,
        ["GGT"] = 11.3, ["GGC"] = 21.7, ["GGA"] = 17.2, ["GGG"] = 15.6,
    }, GeneticCode.Standard, taxonomyId: 10029);

    /// <summary>
    /// <i>Komagataella phaffii</i> (<i>Pichia pastoris</i>), the methylotrophic yeast used for secreted protein
    /// production. AT-rich like <i>S. cerevisiae</i>: TTG for leucine, AGA for arginine, GAA for glutamate,
    /// GCT for alanine and GGT for glycine, with CGG and CGC almost absent.
    /// </summary>
    public static readonly CodonUsageTable PichiaPastoris = FromCounts(
        "P. pastoris", "Komagataella phaffii (Pichia pastoris)", new Dictionary<string, double>(64)
    {
        ["TTT"] = 22.7, ["TTC"] = 20.3, ["TTA"] = 15.3, ["TTG"] = 32.2,
        ["CTT"] = 15.9, ["CTC"] = 8.6, ["CTA"] = 10.2, ["CTG"] = 15.3,
        ["ATT"] = 29.8, ["ATC"] = 20.9, ["ATA"] = 10.3, ["ATG"] = 21.0,
        ["GTT"] = 25.9, ["GTC"] = 15.3, ["GTA"] = 9.1, ["GTG"] = 13.5,
        ["TCT"] = 24.2, ["TCC"] = 15.8, ["TCA"] = 14.2, ["TCG"] = 8.2,
        ["CCT"] = 15.4, ["CCC"] = 7.4, ["CCA"] = 18.3, ["CCG"] = 5.6,
        ["ACT"] = 21.3, ["ACC"] = 16.2, ["ACA"] = 13.4, ["ACG"] = 7.2,
        ["GCT"] = 25.6, ["GCC"] = 16.2, ["GCA"] = 13.9, ["GCG"] = 5.1,
        ["TAT"] = 16.4, ["TAC"] = 18.7, ["TAA"] = 1.0, ["TAG"] = 0.6,
        ["CAT"] = 12.2, ["CAC"] = 10.0, ["CAA"] = 25.5, ["CAG"] = 13.5,
        ["AAT"] = 24.5, ["AAC"] = 24.8, ["AAA"] = 29.8, ["AAG"] = 31.5,
        ["GAT"] = 36.5, ["GAC"] = 22.4, ["GAA"] = 41.3, ["GAG"] = 25.9,
        ["TGT"] = 7.1, ["TGC"] = 4.0, ["TGA"] = 0.7, ["TGG"] = 10.4,
        ["CGT"] = 8.2, ["CGC"] = 3.1, ["CGA"] = 5.1, ["CGG"] = 2.5,
        ["AGT"] = 14.1, ["AGC"] = 8.7, ["AGA"] = 20.7, ["AGG"] = 7.5,
        ["GGT"] = 24.7, ["GGC"] = 9.7, ["GGA"] = 18.2, ["GGG"] = 6.4,
    }, GeneticCode.Standard, taxonomyId: 4922);

    /// <summary>All built-in host tables, in the order they are declared.</summary>
    public static IReadOnlyList<CodonUsageTable> All { get; } =
        [EColiK12, SaccharomycesCerevisiae, HomoSapiens, CricetulusGriseus, PichiaPastoris];

    /// <summary>
    /// Looks a built-in table up by a loose name: the table name, the organism, or a common short form such as
    /// <c>"ecoli"</c>, <c>"yeast"</c>, <c>"human"</c>, <c>"cho"</c> or <c>"pichia"</c>. Case and punctuation are ignored.
    /// </summary>
    /// <exception cref="ArgumentException">No built-in table matches.</exception>
    public static CodonUsageTable FromName(string name)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        string key = Simplify(name);
        switch (key)
        {
            case "ecoli" or "ecolik12" or "escherichiacoli" or "escherichiacolik12" or "bacteria": return EColiK12;
            case "yeast" or "scerevisiae" or "saccharomyces" or "saccharomycescerevisiae": return SaccharomycesCerevisiae;
            case "human" or "hsapiens" or "homosapiens" or "hek293": return HomoSapiens;
            case "cho" or "cgriseus" or "cricetulusgriseus" or "chinesehamster" or "chocgriseus": return CricetulusGriseus;
            case "pichia" or "ppastoris" or "pichiapastoris" or "komagataellaphaffii": return PichiaPastoris;
        }
        for (int i = 0; i < All.Count; i++)
        {
            var table = All[i];
            if (Simplify(table.Name) == key || Simplify(table.Organism) == key) return table;
        }
        throw new ArgumentException(
            $"No built-in codon usage table matches '{name}'. Known tables: {string.Join(", ", All.Select(t => t.Name))}.",
            nameof(name));
    }

    private static string Simplify(string value)
    {
        var builder = new StringBuilder(value.Length);
        foreach (char c in value)
            if (char.IsLetterOrDigit(c)) builder.Append(char.ToLowerInvariant(c));
        return builder.ToString();
    }
}
