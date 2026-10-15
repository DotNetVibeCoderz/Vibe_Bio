using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>Sequence file formats understood by <see cref="SequenceFile"/>.</summary>
public enum SequenceFormat
{
    Unknown,
    Fasta,
    Fastq,
    GenBank,
    Embl,
    Pdb,
}

/// <summary>
/// One entry point for reading and writing sequence files whatever their format. The format comes from the
/// extension (ignoring a trailing <c>.gz</c>) and, when that is inconclusive, from the first non-blank line.
/// </summary>
public static class SequenceFile
{
    private static readonly Dictionary<string, SequenceFormat> Extensions = new(StringComparer.OrdinalIgnoreCase)
    {
        [".fa"] = SequenceFormat.Fasta, [".fasta"] = SequenceFormat.Fasta, [".fna"] = SequenceFormat.Fasta,
        [".ffn"] = SequenceFormat.Fasta, [".faa"] = SequenceFormat.Fasta, [".frn"] = SequenceFormat.Fasta,
        [".fas"] = SequenceFormat.Fasta, [".fsa"] = SequenceFormat.Fasta, [".mfa"] = SequenceFormat.Fasta, [".seq"] = SequenceFormat.Fasta,
        [".fq"] = SequenceFormat.Fastq, [".fastq"] = SequenceFormat.Fastq,
        [".gb"] = SequenceFormat.GenBank, [".gbk"] = SequenceFormat.GenBank, [".gbff"] = SequenceFormat.GenBank,
        [".genbank"] = SequenceFormat.GenBank, [".gp"] = SequenceFormat.GenBank, [".gpff"] = SequenceFormat.GenBank,
        [".embl"] = SequenceFormat.Embl, [".emb"] = SequenceFormat.Embl,
        [".pdb"] = SequenceFormat.Pdb, [".ent"] = SequenceFormat.Pdb,
    };

    /// <summary>Detects the format of a file from its extension, then its content.</summary>
    public static SequenceFormat DetectFormat(string path)
    {
        if (Extensions.TryGetValue(BioFile.InnerExtension(path), out var f)) return f;
        using var reader = BioFile.OpenText(path);
        return DetectFormat(reader);
    }

    /// <summary>Detects the format from the first non-blank line of <paramref name="reader"/> (which is consumed).</summary>
    public static SequenceFormat DetectFormat(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Trim().Length == 0) continue;
            if (line[0] == '>') return SequenceFormat.Fasta;
            if (line[0] == '@') return SequenceFormat.Fastq;
            if (line.StartsWith("LOCUS", StringComparison.Ordinal)) return SequenceFormat.GenBank;
            if (line.StartsWith("ID   ", StringComparison.Ordinal)) return SequenceFormat.Embl;
            if (line.StartsWith("HEADER", StringComparison.Ordinal) || line.StartsWith("SEQRES", StringComparison.Ordinal)
                || line.StartsWith("ATOM  ", StringComparison.Ordinal) || line.StartsWith("REMARK", StringComparison.Ordinal)
                || line.StartsWith("COMPND", StringComparison.Ordinal) || line.StartsWith("TITLE ", StringComparison.Ordinal))
                return SequenceFormat.Pdb;
            return SequenceFormat.Unknown;
        }
        return SequenceFormat.Unknown;
    }

    /// <summary>
    /// Reads the sequences of any supported file. FASTQ yields the bases of each read (use <see cref="ReadReads"/>
    /// to keep qualities); PDB yields SEQRES chains, or coordinate-derived chains when there is no SEQRES.
    /// </summary>
    /// <exception cref="InvalidDataException">The format cannot be determined.</exception>
    public static IEnumerable<Sequence> Read(string path, SequenceFormat format = SequenceFormat.Unknown)
    {
        if (format == SequenceFormat.Unknown) format = DetectFormat(path);
        return format switch
        {
            SequenceFormat.Fasta => FastaReader.Read(path),
            SequenceFormat.Fastq => FastqReader.Read(path).Select(r => r.Sequence),
            SequenceFormat.GenBank => GenBankReader.Read(path),
            SequenceFormat.Embl => EmblReader.Read(path),
            SequenceFormat.Pdb => PdbSequenceReader.Read(path),
            _ => throw new InvalidDataException($"Cannot tell the sequence format of '{path}'."),
        };
    }

    /// <summary>Reads the first sequence of a file.</summary>
    public static Sequence ReadFirst(string path) => Read(path).FirstOrDefault()
        ?? throw new InvalidDataException($"'{path}' contains no sequences.");

    /// <summary>Reads sequencing reads with qualities from a FASTQ file.</summary>
    public static IEnumerable<SequenceRead> ReadReads(string path, FastqEncoding encoding = FastqEncoding.Auto) => FastqReader.Read(path, encoding);

    /// <summary>
    /// Writes sequences in <paramref name="format"/> (from the extension when Unknown). PDB is read-only; FASTQ needs
    /// <see cref="FastqWriter"/> because qualities are required.
    /// </summary>
    public static void Write(string path, IEnumerable<Sequence> sequences, SequenceFormat format = SequenceFormat.Unknown)
    {
        if (format == SequenceFormat.Unknown)
            format = Extensions.TryGetValue(BioFile.InnerExtension(path), out var f) ? f : SequenceFormat.Fasta;
        switch (format)
        {
            case SequenceFormat.Fasta: FastaWriter.Write(path, sequences); break;
            case SequenceFormat.GenBank: GenBankWriter.Write(path, sequences); break;
            case SequenceFormat.Embl: EmblWriter.Write(path, sequences); break;
            default: throw new NotSupportedException($"Writing {format} from plain sequences is not supported.");
        }
    }
}
