namespace Biotech.Net.Core;

/// <summary>Classification of a sequence variant.</summary>
public enum VariantType
{
    Snv,
    Mnv,
    Insertion,
    Deletion,
    Complex,
    Structural,
}

/// <summary>
/// A sample genotype: allele indices (0 = reference, 1.. = alternates, −1 = missing) and phasing.
/// </summary>
public readonly record struct Genotype(int[] Alleles, bool Phased = false)
{
    public static readonly Genotype Missing = new([-1, -1]);

    /// <summary>Parses VCF GT notation: <c>0/1</c>, <c>1|1</c>, <c>./.</c>, <c>1</c>.</summary>
    public static Genotype Parse(ReadOnlySpan<char> gt)
    {
        bool phased = gt.Contains('|');
        Span<Range> parts = stackalloc Range[8];
        int n = gt.SplitAny(parts, "/|");
        var alleles = new int[n];
        for (int i = 0; i < n; i++)
        {
            var p = gt[parts[i]];
            alleles[i] = p is "." ? -1 : int.Parse(p);
        }
        return new(alleles, phased);
    }

    public bool IsMissing => Alleles.All(a => a < 0);
    public bool IsHomozygousReference => !IsMissing && Alleles.All(a => a == 0);
    public bool IsHeterozygous => Alleles.Where(a => a >= 0).Distinct().Count() > 1;
    public bool IsHomozygousAlternate => !IsMissing && Alleles.All(a => a > 0) && Alleles.Distinct().Count() == 1;

    /// <summary>Number of non-reference alleles (dosage) for bi-allelic sites; missing alleles count 0.</summary>
    public int AlternateCount => Alleles.Count(a => a > 0);

    public override string ToString() => string.Join(Phased ? "|" : "/", Alleles.Select(a => a < 0 ? "." : a.ToString()));
}

/// <summary>A sequence variant at one position (one VCF record).</summary>
public sealed class Variant
{
    public Variant(string chromosome, long position, string reference, IReadOnlyList<string> alternates, string? id = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(chromosome);
        ArgumentException.ThrowIfNullOrEmpty(reference);
        if (position < 1) throw new ArgumentOutOfRangeException(nameof(position), "VCF positions are 1-based.");
        Chromosome = chromosome;
        Position = position;
        Reference = reference;
        Alternates = alternates;
        Id = string.IsNullOrEmpty(id) || id == "." ? null : id;
    }

    /// <summary>Convenience constructor for a bi-allelic variant.</summary>
    public Variant(string chromosome, long position, string reference, string alternate, string? id = null)
        : this(chromosome, position, reference, [alternate], id) { }

    public string Chromosome { get; }

    /// <summary>1-based position of the first reference base.</summary>
    public long Position { get; }

    /// <summary>rsID or other identifier; null when absent.</summary>
    public string? Id { get; }

    public string Reference { get; }
    public IReadOnlyList<string> Alternates { get; }

    /// <summary>Phred-scaled quality, NaN when missing.</summary>
    public double Quality { get; init; } = double.NaN;

    /// <summary>FILTER column; <c>PASS</c> or a list of failed filters.</summary>
    public string Filter { get; init; } = "PASS";

    /// <summary>INFO fields; flags have an empty value.</summary>
    public Dictionary<string, string> Info { get; } = new(StringComparer.Ordinal);

    /// <summary>Per-sample genotypes, in the sample order of the source file.</summary>
    public List<Genotype> Genotypes { get; } = [];

    public bool IsPass => Filter is "PASS" or ".";

    /// <summary>The 0-based half-open span of the reference allele.</summary>
    public GenomicInterval Interval => new(Chromosome, Position - 1, Position - 1 + Reference.Length);

    /// <summary>Type of the first alternate allele.</summary>
    public VariantType Type => Classify(Reference, Alternates.Count > 0 ? Alternates[0] : Reference);

    /// <summary>Whether this is a single-nucleotide purine↔purine or pyrimidine↔pyrimidine change.</summary>
    public bool IsTransition
    {
        get
        {
            if (Type != VariantType.Snv) return false;
            char r = char.ToUpperInvariant(Reference[0]), a = char.ToUpperInvariant(Alternates[0][0]);
            return (r, a) is ('A', 'G') or ('G', 'A') or ('C', 'T') or ('T', 'C');
        }
    }

    /// <summary>Single-nucleotide change that is not a transition.</summary>
    public bool IsTransversion => Type == VariantType.Snv && !IsTransition;

    /// <summary>Alternate allele frequency across called genotypes (first alternate), NaN when no calls.</summary>
    public double AlternateAlleleFrequency()
    {
        int alt = 0, total = 0;
        foreach (var g in Genotypes)
            foreach (var a in g.Alleles)
            {
                if (a < 0) continue;
                total++;
                if (a == 1) alt++;
            }
        return total == 0 ? double.NaN : (double)alt / total;
    }

    public static VariantType Classify(string reference, string alternate)
    {
        if (alternate.StartsWith('<') || alternate.Contains('[') || alternate.Contains(']')) return VariantType.Structural;
        if (reference.Length == alternate.Length)
            return reference.Length == 1 ? VariantType.Snv : VariantType.Mnv;
        if (reference.Length < alternate.Length && alternate.StartsWith(reference, StringComparison.OrdinalIgnoreCase))
            return VariantType.Insertion;
        if (reference.Length > alternate.Length && reference.StartsWith(alternate, StringComparison.OrdinalIgnoreCase))
            return VariantType.Deletion;
        return VariantType.Complex;
    }

    public override string ToString() => $"{Chromosome}:{Position} {Reference}>{string.Join(',', Alternates)}{(Id is null ? "" : $" ({Id})")}";
}
