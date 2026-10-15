namespace Biotech.Net.Core;

/// <summary>The kind of biological polymer a sequence represents.</summary>
public enum SequenceType
{
    Unknown,
    Dna,
    Rna,
    Protein,
}

/// <summary>
/// A set of allowed residue symbols. Validation is a 256-entry lookup, case-insensitive, so checking a
/// sequence costs one table read per byte. Gaps (<c>-</c>, <c>.</c>) are accepted by the IUPAC/gapped alphabets
/// so aligned sequences can share the type.
/// </summary>
public sealed class Alphabet
{
    private readonly bool[] _valid = new bool[256];

    private Alphabet(string name, SequenceType type, string symbols)
    {
        Name = name;
        Type = type;
        Symbols = symbols;
        foreach (char c in symbols)
        {
            _valid[char.ToUpperInvariant(c)] = true;
            _valid[char.ToLowerInvariant(c)] = true;
        }
    }

    /// <summary>Unambiguous DNA: A, C, G, T.</summary>
    public static readonly Alphabet Dna = new("DNA", SequenceType.Dna, "ACGT");

    /// <summary>DNA with IUPAC ambiguity codes and gaps.</summary>
    public static readonly Alphabet DnaIupac = new("DNA (IUPAC)", SequenceType.Dna, "ACGTRYSWKMBDHVN-.");

    /// <summary>Unambiguous RNA: A, C, G, U.</summary>
    public static readonly Alphabet Rna = new("RNA", SequenceType.Rna, "ACGU");

    /// <summary>RNA with IUPAC ambiguity codes and gaps.</summary>
    public static readonly Alphabet RnaIupac = new("RNA (IUPAC)", SequenceType.Rna, "ACGURYSWKMBDHVN-.");

    /// <summary>The 20 standard amino acids.</summary>
    public static readonly Alphabet Protein = new("Protein", SequenceType.Protein, "ACDEFGHIKLMNPQRSTVWY");

    /// <summary>Amino acids plus B, Z, J, U (Sec), O (Pyl), X, stop (<c>*</c>) and gaps.</summary>
    public static readonly Alphabet ProteinExtended = new("Protein (extended)", SequenceType.Protein, "ACDEFGHIKLMNPQRSTVWYBZJUOX*-.");

    /// <summary>The symbols, in upper case.</summary>
    public string Symbols { get; }

    /// <summary>Display name.</summary>
    public string Name { get; }

    /// <summary>The polymer type.</summary>
    public SequenceType Type { get; }

    /// <summary>True for nucleotide alphabets.</summary>
    public bool IsNucleotide => Type is SequenceType.Dna or SequenceType.Rna;

    /// <summary>Whether <paramref name="symbol"/> belongs to the alphabet (case-insensitive).</summary>
    public bool IsValid(byte symbol) => _valid[symbol];

    /// <summary>Whether every symbol of <paramref name="residues"/> belongs to the alphabet.</summary>
    public bool IsValid(ReadOnlySpan<byte> residues) => FirstInvalid(residues) < 0;

    /// <summary>Index of the first symbol not in the alphabet, or −1.</summary>
    public int FirstInvalid(ReadOnlySpan<byte> residues)
    {
        var valid = _valid;
        for (int i = 0; i < residues.Length; i++)
            if (!valid[residues[i]]) return i;
        return -1;
    }

    /// <summary>The IUPAC (ambiguity-tolerant) alphabet for a sequence type.</summary>
    public static Alphabet ForType(SequenceType type) => type switch
    {
        SequenceType.Dna => DnaIupac,
        SequenceType.Rna => RnaIupac,
        SequenceType.Protein => ProteinExtended,
        _ => ProteinExtended,
    };

    /// <summary>
    /// Guesses the type of a residue string: nucleotide when at least 90% of the letters are A/C/G/T/U/N
    /// (RNA if it has U and no T), protein otherwise.
    /// </summary>
    public static Alphabet Detect(ReadOnlySpan<byte> residues)
    {
        long nucleotide = 0, letters = 0, t = 0, u = 0;
        foreach (byte b in residues)
        {
            byte c = (byte)(b & 0xDF);
            if (c is < (byte)'A' or > (byte)'Z') continue;
            letters++;
            switch (c)
            {
                case (byte)'A' or (byte)'C' or (byte)'G' or (byte)'N': nucleotide++; break;
                case (byte)'T': nucleotide++; t++; break;
                case (byte)'U': nucleotide++; u++; break;
            }
        }
        if (letters == 0) return DnaIupac;
        if (nucleotide >= 0.9 * letters) return u > 0 && t == 0 ? RnaIupac : DnaIupac;
        return ProteinExtended;
    }

    public override string ToString() => Name;
}
