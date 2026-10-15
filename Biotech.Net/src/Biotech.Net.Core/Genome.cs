namespace Biotech.Net.Core;

/// <summary>A gene on a reference genome, with its transcripts.</summary>
public sealed class Gene
{
    public Gene(string id, GenomicInterval location, string? name = null, string biotype = "protein_coding")
    {
        ArgumentException.ThrowIfNullOrEmpty(id);
        Id = id;
        Location = location;
        Name = name ?? id;
        Biotype = biotype;
    }

    /// <summary>Stable identifier (e.g. ENSG..., locus tag).</summary>
    public string Id { get; }

    /// <summary>Symbol (e.g. <c>TP53</c>).</summary>
    public string Name { get; }

    /// <summary>Genomic span of the gene.</summary>
    public GenomicInterval Location { get; }

    /// <summary>Biotype: protein_coding, lncRNA, tRNA, ...</summary>
    public string Biotype { get; }

    /// <summary>Strand of the gene.</summary>
    public Strand Strand => Location.Strand;

    /// <summary>Transcripts (isoforms).</summary>
    public List<Transcript> Transcripts { get; } = [];

    /// <summary>Extra attributes (from GFF/GTF column 9, etc.).</summary>
    public Dictionary<string, string> Attributes { get; } = new(StringComparer.OrdinalIgnoreCase);

    public override string ToString() => $"{Name} ({Id}) {Location}";
}

/// <summary>A transcript: ordered exons and optional coding region.</summary>
public sealed class Transcript
{
    public Transcript(string id, string geneId, IEnumerable<GenomicInterval> exons)
    {
        Id = id;
        GeneId = geneId;
        Exons = exons.OrderBy(e => e.Start).ToList();
        if (Exons.Count == 0) throw new ArgumentException("A transcript has at least one exon.", nameof(exons));
    }

    public string Id { get; }
    public string GeneId { get; }

    /// <summary>Exons in genomic (ascending) order.</summary>
    public IReadOnlyList<GenomicInterval> Exons { get; }

    /// <summary>Coding region, when the transcript is protein-coding.</summary>
    public GenomicInterval? CodingRegion { get; init; }

    /// <summary>Genomic span from first to last exon.</summary>
    public GenomicInterval Span => Exons[0] with { End = Exons[^1].End };

    /// <summary>Spliced length (sum of exons).</summary>
    public long SplicedLength => Exons.Sum(e => e.Length);

    /// <summary>Introns between consecutive exons.</summary>
    public IEnumerable<GenomicInterval> Introns()
    {
        for (int i = 1; i < Exons.Count; i++)
            yield return Exons[i] with { Start = Exons[i - 1].End, End = Exons[i].Start };
    }

    /// <summary>Spliced transcript sequence from the chromosome sequence (reverse-complemented on the minus strand).</summary>
    public Sequence SplicedSequence(Sequence chromosome)
    {
        var loc = new FeatureLocation(Exons.Select(e => ((int)e.Start, (int)e.End)).ToList(), Exons[0].Strand);
        return loc.Extract(chromosome, Id);
    }
}

/// <summary>A chromosome, contig or scaffold of a genome.</summary>
public sealed class Chromosome
{
    public Chromosome(string name, long length, Sequence? sequence = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        if (sequence is not null && sequence.Length != length)
            throw new ArgumentException($"Chromosome {name}: length {length} but sequence has {sequence.Length} bases.");
        Name = name;
        Length = length;
        Sequence = sequence;
    }

    /// <summary>Creates a chromosome from its sequence.</summary>
    public static Chromosome FromSequence(Sequence sequence) => new(sequence.Id, sequence.Length, sequence);

    public string Name { get; }
    public long Length { get; }

    /// <summary>The bases, when loaded.</summary>
    public Sequence? Sequence { get; }

    /// <summary>Circular molecules (bacterial chromosomes, plasmids, mitochondria).</summary>
    public bool IsCircular { get; init; }

    public override string ToString() => $"{Name} ({Length:N0} bp{(IsCircular ? ", circular" : "")})";
}

/// <summary>A genome assembly: chromosomes plus gene annotation.</summary>
public sealed class Genome
{
    private readonly Dictionary<string, Chromosome> _byName = new(StringComparer.Ordinal);
    private readonly List<Chromosome> _chromosomes = [];
    private readonly List<Gene> _genes = [];

    public Genome(string name, string? species = null, string? assembly = null)
    {
        Name = name;
        Species = species ?? "unknown";
        Assembly = assembly ?? name;
    }

    public string Name { get; }
    public string Species { get; }
    public string Assembly { get; }

    public IReadOnlyList<Chromosome> Chromosomes => _chromosomes;
    public IReadOnlyList<Gene> Genes => _genes;

    /// <summary>Total assembly length.</summary>
    public long Length => _chromosomes.Sum(c => c.Length);

    public Genome Add(Chromosome chromosome)
    {
        if (!_byName.TryAdd(chromosome.Name, chromosome))
            throw new ArgumentException($"Chromosome {chromosome.Name} already exists.");
        _chromosomes.Add(chromosome);
        return this;
    }

    public Genome Add(Gene gene)
    {
        _genes.Add(gene);
        return this;
    }

    public Chromosome this[string name] => _byName.TryGetValue(name, out var c) ? c
        : throw new KeyNotFoundException($"No chromosome named '{name}'.");

    public bool TryGetChromosome(string name, out Chromosome chromosome) => _byName.TryGetValue(name, out chromosome!);

    /// <summary>The bases of <paramref name="region"/> (reverse-complemented for minus-strand regions).</summary>
    public Sequence GetSequence(GenomicInterval region)
    {
        var chrom = this[region.Chromosome];
        if (chrom.Sequence is null) throw new InvalidOperationException($"Chromosome {chrom.Name} has no sequence loaded.");
        var s = chrom.Sequence.Subsequence((int)region.Start, (int)region.Length).WithId(region.ToString());
        return region.Strand == Strand.Reverse ? s.ReverseComplement() : s;
    }

    /// <summary>Genes overlapping <paramref name="region"/>.</summary>
    public IEnumerable<Gene> GenesIn(GenomicInterval region) => _genes.Where(g => g.Location.Overlaps(region));

    /// <summary>Finds a gene by id or name (case-insensitive).</summary>
    public Gene? FindGene(string idOrName) => _genes.FirstOrDefault(g =>
        g.Id.Equals(idOrName, StringComparison.OrdinalIgnoreCase) || g.Name.Equals(idOrName, StringComparison.OrdinalIgnoreCase));

    /// <summary>N50: the contig length at which half the assembly is in contigs at least that long.</summary>
    public long N50()
    {
        var lengths = _chromosomes.Select(c => c.Length).OrderDescending().ToList();
        long half = (lengths.Sum() + 1) / 2, acc = 0;
        foreach (var l in lengths) { acc += l; if (acc >= half) return l; }
        return 0;
    }

    public override string ToString() => $"{Name} ({Species}): {_chromosomes.Count} sequences, {Length:N0} bp, {_genes.Count} genes";
}
