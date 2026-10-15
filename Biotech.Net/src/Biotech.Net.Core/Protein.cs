using Biotech.Net.Native;

namespace Biotech.Net.Core;

/// <summary>A protein: its amino-acid sequence plus identity and annotation.</summary>
public sealed class Protein
{
    public Protein(string id, Sequence sequence, string? name = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(id);
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Type != SequenceType.Protein)
            throw new ArgumentException($"Protein {id} needs an amino-acid sequence, got {sequence.Type}.", nameof(sequence));
        Id = id;
        Sequence = sequence;
        Name = name ?? id;
    }

    /// <summary>Convenience: a protein from a one-letter string.</summary>
    public static Protein FromString(string id, string residues, string? name = null)
        => new(id, Sequence.Protein(residues, id), name);

    /// <summary>Accession (UniProt, RefSeq, ...).</summary>
    public string Id { get; }

    /// <summary>Recommended name.</summary>
    public string Name { get; }

    /// <summary>The amino acids.</summary>
    public Sequence Sequence { get; }

    public string? Organism { get; init; }
    public string? GeneName { get; init; }
    public string? Function { get; init; }

    /// <summary>Identifier of a 3-D structure (PDB id), when known.</summary>
    public string? StructureId { get; init; }

    /// <summary>Additional annotations (EC number, GO terms, ...).</summary>
    public Dictionary<string, string> Annotations { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>Number of residues.</summary>
    public int Length => Sequence.Length;

    /// <summary>Average molecular weight in daltons.</summary>
    public double MolecularWeight => BioKernels.ProteinMass(Sequence.Span);

    /// <summary>Isoelectric point.</summary>
    public double IsoelectricPoint => BioKernels.IsoelectricPoint(Sequence.Span);

    public override string ToString() => $"{Name} ({Id}, {Length} aa, {MolecularWeight / 1000:F1} kDa)";
}
