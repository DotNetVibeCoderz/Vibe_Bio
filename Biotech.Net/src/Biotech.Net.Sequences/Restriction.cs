using System.Collections.Frozen;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>The kind of ends a restriction enzyme leaves.</summary>
public enum OverhangType
{
    Blunt,
    FivePrime,
    ThreePrime,
}

/// <summary>
/// A type II restriction enzyme. Cut positions are measured from the first base of the recognition site on the
/// top strand: <see cref="CutTop"/> is where the top strand is cut (EcoRI G^AATTC → 1) and <see cref="CutBottom"/>
/// where the bottom strand is cut, in top-strand coordinates (EcoRI → 5). Type IIS enzymes cut outside the
/// site, e.g. BsaI GGTCTC(1/5) → 7 and 11.
/// </summary>
public sealed record RestrictionEnzyme(string Name, string Site, int CutTop, int CutBottom)
{
    /// <summary>Recognition site length.</summary>
    public int SiteLength => Site.Length;

    /// <summary>Whether the site equals its own reverse complement (so it is only searched on the top strand).</summary>
    public bool IsPalindromic { get; } = string.Equals(Site, Iupac.ReverseComplement(Site), StringComparison.OrdinalIgnoreCase);

    /// <summary>Overhang length: positive for 5' overhangs, negative for 3', zero for blunt.</summary>
    public int Overhang => CutBottom - CutTop;

    /// <summary>The overhang type.</summary>
    public OverhangType OverhangType => Overhang == 0 ? OverhangType.Blunt : Overhang > 0 ? OverhangType.FivePrime : OverhangType.ThreePrime;

    /// <summary>Whether cuts fall outside the recognition site (type IIS).</summary>
    public bool CutsOutsideSite => CutTop < 0 || CutTop > SiteLength || CutBottom < 0 || CutBottom > SiteLength;

    /// <summary>REBASE-style notation, e.g. <c>G^AATTC</c> or <c>GGTCTC(1/5)</c>.</summary>
    public override string ToString() => CutsOutsideSite
        ? $"{Name} {Site}({CutTop - SiteLength}/{CutBottom - SiteLength})"
        : $"{Name} {Site[..CutTop]}^{Site[CutTop..]}";
}

/// <summary>Built-in restriction enzymes (recognition sites and cut positions from REBASE).</summary>
public static class RestrictionEnzymes
{
    /// <summary>Every built-in enzyme.</summary>
    public static IReadOnlyList<RestrictionEnzyme> All { get; } =
    [
        new("AatII", "GACGTC", 5, 1), new("AflII", "CTTAAG", 1, 5), new("AgeI", "ACCGGT", 1, 5), new("AluI", "AGCT", 2, 2),
        new("ApaI", "GGGCCC", 5, 1), new("ApaLI", "GTGCAC", 1, 5), new("AscI", "GGCGCGCC", 2, 6), new("AvrII", "CCTAGG", 1, 5),
        new("BamHI", "GGATCC", 1, 5), new("BbsI", "GAAGAC", 8, 12), new("BglII", "AGATCT", 1, 5), new("BsaI", "GGTCTC", 7, 11),
        new("BsiWI", "CGTACG", 1, 5), new("BsmBI", "CGTCTC", 7, 11), new("BsrGI", "TGTACA", 1, 5), new("BstEII", "GGTNACC", 1, 6),
        new("BstXI", "CCANNNNNNTGG", 8, 4), new("ClaI", "ATCGAT", 2, 4), new("DpnII", "GATC", 0, 4), new("DraI", "TTTAAA", 3, 3),
        new("EagI", "CGGCCG", 1, 5), new("EcoRI", "GAATTC", 1, 5), new("EcoRV", "GATATC", 3, 3), new("FseI", "GGCCGGCC", 6, 2),
        new("HaeIII", "GGCC", 2, 2), new("HhaI", "GCGC", 3, 1), new("HindIII", "AAGCTT", 1, 5), new("HinfI", "GANTC", 1, 4),
        new("HpaI", "GTTAAC", 3, 3), new("HpaII", "CCGG", 1, 3), new("KasI", "GGCGCC", 1, 5), new("KpnI", "GGTACC", 5, 1),
        new("MboI", "GATC", 0, 4), new("MfeI", "CAATTG", 1, 5), new("MluI", "ACGCGT", 1, 5), new("MseI", "TTAA", 1, 3),
        new("MspI", "CCGG", 1, 3), new("NarI", "GGCGCC", 2, 4), new("NcoI", "CCATGG", 1, 5), new("NdeI", "CATATG", 2, 4),
        new("NheI", "GCTAGC", 1, 5), new("NlaIII", "CATG", 4, 0), new("NotI", "GCGGCCGC", 2, 6), new("NsiI", "ATGCAT", 5, 1),
        new("PacI", "TTAATTAA", 5, 3), new("PmeI", "GTTTAAAC", 4, 4), new("PstI", "CTGCAG", 5, 1), new("PvuI", "CGATCG", 4, 2),
        new("PvuII", "CAGCTG", 3, 3), new("RsaI", "GTAC", 2, 2), new("SacI", "GAGCTC", 5, 1), new("SacII", "CCGCGG", 4, 2),
        new("SalI", "GTCGAC", 1, 5), new("SapI", "GCTCTTC", 8, 11), new("Sau3AI", "GATC", 0, 4), new("ScaI", "AGTACT", 3, 3),
        new("SfiI", "GGCCNNNNNGGCC", 8, 5), new("SmaI", "CCCGGG", 3, 3), new("SpeI", "ACTAGT", 1, 5), new("SphI", "GCATGC", 5, 1),
        new("StuI", "AGGCCT", 3, 3), new("SwaI", "ATTTAAAT", 4, 4), new("TaqI", "TCGA", 1, 3), new("XbaI", "TCTAGA", 1, 5),
        new("XhoI", "CTCGAG", 1, 5), new("XmaI", "CCCGGG", 1, 5), new("XmnI", "GAANNNNTTC", 5, 5),
    ];

    private static readonly FrozenDictionary<string, RestrictionEnzyme> ByName = All.ToFrozenDictionary(e => e.Name, StringComparer.OrdinalIgnoreCase);

    /// <summary>Looks an enzyme up by name (case-insensitive).</summary>
    public static RestrictionEnzyme Get(string name) => ByName.TryGetValue(name, out var e) ? e
        : throw new KeyNotFoundException($"Restriction enzyme '{name}' is not built in.");

    /// <summary>Tries to look an enzyme up by name.</summary>
    public static bool TryGet(string name, out RestrictionEnzyme enzyme) => ByName.TryGetValue(name, out enzyme!);

    /// <summary>The enzymes that cut <paramref name="sequence"/> exactly once.</summary>
    public static IReadOnlyList<RestrictionEnzyme> SingleCutters(Sequence sequence, bool circular = false)
        => All.Where(e => RestrictionDigest.FindCuts(sequence, e, circular).Count == 1).ToList();

    /// <summary>The enzymes with no site in <paramref name="sequence"/>.</summary>
    public static IReadOnlyList<RestrictionEnzyme> NonCutters(Sequence sequence)
        => All.Where(e => MotifSearch.FindIupac(sequence, e.Site).Count == 0).ToList();
}

/// <summary>A cut made by an enzyme.</summary>
/// <param name="Position">0-based top-strand position of the first base after the cut (Biopython's <c>search</c> reports this plus one).</param>
/// <param name="SiteStart">0-based start of the recognition site on the top strand.</param>
/// <param name="SiteStrand">Strand the site was found on (reverse only for non-palindromic sites).</param>
public readonly record struct RestrictionCut(RestrictionEnzyme Enzyme, int Position, int SiteStart, Strand SiteStrand);

/// <summary>A fragment of a digest. For circular digests the last fragment may wrap (<see cref="End"/> &lt; <see cref="Start"/>).</summary>
public readonly record struct RestrictionFragment(int Start, int End, int Length);

/// <summary>The outcome of a digest.</summary>
public sealed record DigestResult(IReadOnlyList<RestrictionCut> Cuts, IReadOnlyList<RestrictionFragment> Fragments, int SequenceLength, bool Circular)
{
    /// <summary>Fragment lengths, largest first (what a gel shows).</summary>
    public IReadOnlyList<int> FragmentSizes => Fragments.Select(f => f.Length).OrderDescending().ToList();

    /// <summary>The fragment sequences (wrapping fragments of a circular molecule are joined across the origin).</summary>
    public IReadOnlyList<Sequence> FragmentSequences(Sequence source)
    {
        var list = new List<Sequence>(Fragments.Count);
        int i = 0;
        foreach (var f in Fragments)
        {
            i++;
            string id = $"{source.Id}_frag{i}";
            if (f.End >= f.Start && f.Length == f.End - f.Start) list.Add(source.Subsequence(f.Start, f.Length).WithId(id));
            else list.Add(source.Subsequence(f.Start, source.Length - f.Start).Concat(source.Subsequence(0, f.End), id));
        }
        return list;
    }
}

/// <summary>
/// Restriction digests. Sites are found with the IUPAC matcher; non-palindromic sites are also searched on the
/// bottom strand, where the cut is mirrored (top-strand cut = site start + length − <see cref="RestrictionEnzyme.CutBottom"/>).
/// </summary>
public static class RestrictionDigest
{
    /// <summary>The cuts one enzyme makes, sorted by position. Cuts that would fall outside a linear molecule are dropped.</summary>
    public static IReadOnlyList<RestrictionCut> FindCuts(Sequence sequence, RestrictionEnzyme enzyme, bool circular = false)
    {
        int n = sequence.Length;
        var target = sequence;
        if (circular && n > 0)
        {
            // Sites spanning the origin: search the sequence extended by (site length − 1) bases.
            int extra = Math.Min(enzyme.SiteLength - 1, n);
            target = sequence.Concat(sequence.Subsequence(0, extra));
        }
        var cuts = new List<RestrictionCut>();
        var seen = new HashSet<int>();
        foreach (var m in MotifSearch.FindIupac(target, enzyme.Site, bothStrands: !enzyme.IsPalindromic))
        {
            if (m.Start >= n) continue;
            int cut = m.Strand == Strand.Forward ? m.Start + enzyme.CutTop : m.Start + enzyme.SiteLength - enzyme.CutBottom;
            if (circular) cut = ((cut % n) + n) % n;
            else if (cut <= 0 || cut >= n) continue;
            if (seen.Add(cut)) cuts.Add(new RestrictionCut(enzyme, cut, m.Start, m.Strand));
        }
        cuts.Sort((a, b) => a.Position.CompareTo(b.Position));
        return cuts;
    }

    /// <summary>Digests <paramref name="sequence"/> with one or more enzymes.</summary>
    public static DigestResult Digest(Sequence sequence, IEnumerable<RestrictionEnzyme> enzymes, bool circular = false)
    {
        var cuts = enzymes.SelectMany(e => FindCuts(sequence, e, circular)).OrderBy(c => c.Position).ToList();
        int n = sequence.Length;
        var positions = cuts.Select(c => c.Position).Distinct().ToList();
        var fragments = new List<RestrictionFragment>();
        if (!circular)
        {
            int prev = 0;
            foreach (int p in positions) { fragments.Add(new(prev, p, p - prev)); prev = p; }
            fragments.Add(new(prev, n, n - prev));
        }
        else if (positions.Count == 0) fragments.Add(new(0, n, n));
        else
        {
            for (int i = 0; i + 1 < positions.Count; i++) fragments.Add(new(positions[i], positions[i + 1], positions[i + 1] - positions[i]));
            int last = positions[^1], first = positions[0];
            fragments.Add(new(last, first, n - last + first));
        }
        return new DigestResult(cuts, fragments, n, circular);
    }

    /// <summary>Digests with enzymes given by name, e.g. <c>Digest(plasmid, true, "EcoRI", "BamHI")</c>.</summary>
    public static DigestResult Digest(Sequence sequence, bool circular, params string[] enzymeNames)
        => Digest(sequence, enzymeNames.Select(RestrictionEnzymes.Get), circular);
}
