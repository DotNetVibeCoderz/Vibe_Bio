using Biotech.Net.Core;

namespace Biotech.Net.Sequences;

/// <summary>Which codons may open an ORF.</summary>
public enum OrfStartMode
{
    /// <summary>Only ATG.</summary>
    AtgOnly,

    /// <summary>Every start codon of the genetic code (e.g. ATG, GTG, TTG in table 11).</summary>
    AlternativeStarts,

    /// <summary>Stop-to-stop: an ORF begins right after the previous in-frame stop (or at the sequence start).</summary>
    AnyCodon,
}

/// <summary>Options for <see cref="OrfFinder"/>.</summary>
public sealed record OrfOptions
{
    /// <summary>Minimum ORF length in nucleotides, stop codon included (default 75 nt ≈ 24 aa + stop).</summary>
    public int MinLength { get; init; } = 75;

    /// <summary>Genetic code (standard by default).</summary>
    public GeneticCode Code { get; init; } = GeneticCode.Standard;

    /// <summary>Start-codon rule.</summary>
    public OrfStartMode StartMode { get; init; } = OrfStartMode.AtgOnly;

    /// <summary>Report ORFs that run off the end of the sequence without a stop codon.</summary>
    public bool IncludePartial { get; init; }

    /// <summary>Report every in-frame start upstream of a stop (nested ORFs) rather than only the longest.</summary>
    public bool Nested { get; init; }

    /// <summary>Search the forward strand.</summary>
    public bool Forward { get; init; } = true;

    /// <summary>Search the reverse strand.</summary>
    public bool Reverse { get; init; } = true;
}

/// <summary>
/// An open reading frame. <see cref="Start"/>/<see cref="End"/> are 0-based half-open on the forward strand,
/// whatever the strand, so the ORF's bases are <c>parent[Start..End]</c> (reverse-complemented when on the reverse strand).
/// </summary>
/// <param name="Frame">+1..+3 on the forward strand, −1..−3 on the reverse strand (offset from the 5' end of that strand).</param>
/// <param name="HasStop">Whether the ORF ends with a stop codon (false for partial ORFs).</param>
public sealed record Orf(int Start, int End, Strand Strand, int Frame, bool HasStop)
{
    /// <summary>Length in nucleotides, including the stop codon.</summary>
    public int Length => End - Start;

    /// <summary>Length of the encoded protein (stop excluded).</summary>
    public int ProteinLength => (Length / 3) - (HasStop ? 1 : 0);

    /// <summary>The ORF's nucleotides in coding orientation.</summary>
    public Sequence Extract(Sequence parent)
    {
        var s = parent.Subsequence(Start, Length).WithId($"{parent.Id}_{Start + 1}_{End}_{(Strand == Strand.Reverse ? "-" : "+")}");
        return Strand == Strand.Reverse ? s.ReverseComplement() : s;
    }

    /// <summary>The encoded protein (without the terminal stop).</summary>
    public Sequence Translate(Sequence parent, GeneticCode? code = null)
    {
        var protein = Extract(parent).Translate(code);
        return HasStop && protein.Length > 0 ? protein.Subsequence(0, protein.Length - 1) : protein;
    }

    /// <summary>The ORF as a CDS feature of its parent.</summary>
    public SequenceFeature ToFeature() => new("CDS", FeatureLocation.Span(Start, End, Strand));

    public override string ToString() => $"ORF {Start + 1}..{End} ({(Strand == Strand.Reverse ? "-" : "+")}, frame {Frame:+0;-0}, {ProteinLength} aa)";
}

/// <summary>
/// Six-frame open-reading-frame finder. Each strand is scanned once per frame with the codon table (one lookup per
/// codon), the reverse strand on a single reverse-complemented copy, so the cost is linear in sequence length.
/// </summary>
public static class OrfFinder
{
    /// <summary>Finds ORFs in <paramref name="sequence"/>, sorted by start position.</summary>
    public static IReadOnlyList<Orf> Find(Sequence sequence, OrfOptions? options = null)
    {
        options ??= new OrfOptions();
        if (!sequence.Alphabet.IsNucleotide) throw new ArgumentException("ORF finding needs a nucleotide sequence.", nameof(sequence));
        var result = new List<Orf>();
        int n = sequence.Length;
        if (options.Forward) ScanStrand(sequence.Span, n, Strand.Forward, options, result);
        if (options.Reverse)
        {
            var rc = sequence.ReverseComplement();
            ScanStrand(rc.Span, n, Strand.Reverse, options, result);
        }
        result.Sort((a, b) => a.Start != b.Start ? a.Start.CompareTo(b.Start) : b.End.CompareTo(a.End));
        return result;
    }

    /// <summary>The longest ORF, or null when none satisfies the options.</summary>
    public static Orf? Longest(Sequence sequence, OrfOptions? options = null)
    {
        Orf? best = null;
        foreach (var o in Find(sequence, options)) if (best is null || o.Length > best.Length) best = o;
        return best;
    }

    private static void ScanStrand(ReadOnlySpan<byte> s, int n, Strand strand, OrfOptions o, List<Orf> result)
    {
        var table = o.Code.CodonTable;
        for (int frame = 0; frame < 3; frame++)
        {
            // Open starts in this frame (strand-local coordinates); only the first matters unless Nested.
            int openStart = o.StartMode == OrfStartMode.AnyCodon ? frame : -1;
            List<int>? nested = o.Nested ? [] : null;
            if (nested is not null && openStart >= 0) nested.Add(openStart);
            int i = frame;
            for (; i + 3 <= n; i += 3)
            {
                int k = GeneticCode.Index(s.Slice(i, 3));
                bool stop = k >= 0 && table[k] == (byte)'*';
                if (stop)
                {
                    if (nested is not null) foreach (int st in nested) Emit(st, i + 3, true);
                    else if (openStart >= 0) Emit(openStart, i + 3, true);
                    openStart = o.StartMode == OrfStartMode.AnyCodon ? i + 3 : -1;
                    nested?.Clear();
                    if (nested is not null && openStart >= 0) nested.Add(openStart);
                    continue;
                }
                if (o.StartMode == OrfStartMode.AnyCodon) continue;
                bool isStart = o.StartMode == OrfStartMode.AtgOnly ? k == 0b00_11_10 /* ATG */ : o.Code.IsStartIndex(k);
                if (!isStart) continue;
                if (nested is not null) nested.Add(i);
                else if (openStart < 0) openStart = i;
            }
            if (o.IncludePartial)
            {
                int end = i; // last complete codon boundary
                if (nested is not null) foreach (int st in nested) Emit(st, end, false);
                else if (openStart >= 0) Emit(openStart, end, false);
            }

            void Emit(int localStart, int localEnd, bool hasStop)
            {
                if (localEnd - localStart < o.MinLength || localEnd <= localStart) return;
                int fs = strand == Strand.Forward ? localStart : n - localEnd;
                int fe = strand == Strand.Forward ? localEnd : n - localStart;
                result.Add(new Orf(fs, fe, strand, strand == Strand.Forward ? frame + 1 : -(frame + 1), hasStop));
            }
        }
    }
}

/// <summary>Translation in all six reading frames.</summary>
public static class SixFrame
{
    /// <summary>
    /// The six translations, in the order +1, +2, +3, −1, −2, −3. Ids get a <c>_+1</c> … <c>_-3</c> suffix;
    /// stops appear as <c>*</c>.
    /// </summary>
    public static IReadOnlyList<(int Frame, Sequence Protein)> Translate(Sequence sequence, GeneticCode? code = null)
    {
        var rc = sequence.ReverseComplement();
        var list = new List<(int, Sequence)>(6);
        for (int f = 0; f < 3; f++) list.Add((f + 1, sequence.Translate(code, f).WithId($"{sequence.Id}_+{f + 1}")));
        for (int f = 0; f < 3; f++) list.Add((-(f + 1), rc.Translate(code, f).WithId($"{sequence.Id}_-{f + 1}")));
        return list;
    }
}
