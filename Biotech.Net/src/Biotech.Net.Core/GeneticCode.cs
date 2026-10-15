namespace Biotech.Net.Core;

/// <summary>
/// An NCBI genetic code (translation table). Holds the 64-entry codon table in the order the
/// translation kernel expects (A=0, C=1, G=2, T=3; index = 16·b1 + 4·b2 + b3) plus start codons.
/// </summary>
public sealed class GeneticCode
{
    // NCBI order is T, C, A, G for each position; this maps it to the kernel's A, C, G, T order.
    private const string NcbiBases = "TCAG";
    private readonly byte[] _table = new byte[64];
    private readonly bool[] _start = new bool[64];

    private GeneticCode(int id, string name, string aminoAcids, string starts)
    {
        Id = id;
        Name = name;
        for (int i = 0; i < 64; i++)
        {
            int b1 = Code(NcbiBases[i / 16]), b2 = Code(NcbiBases[i / 4 % 4]), b3 = Code(NcbiBases[i % 4]);
            int k = (b1 << 4) | (b2 << 2) | b3;
            _table[k] = (byte)aminoAcids[i];
            _start[k] = starts[i] == 'M';
        }
    }

    private static int Code(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };

    /// <summary>NCBI table number.</summary>
    public int Id { get; }

    /// <summary>NCBI table name.</summary>
    public string Name { get; }

    /// <summary>The 64-byte codon table in kernel order.</summary>
    public ReadOnlySpan<byte> CodonTable => _table;

    /// <summary>Table 1, the standard code.</summary>
    public static readonly GeneticCode Standard = new(1, "Standard",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M---------------M----------------------------");

    /// <summary>Table 2, vertebrate mitochondrial.</summary>
    public static readonly GeneticCode VertebrateMitochondrial = new(2, "Vertebrate Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "----------**--------------------MMMM----------**---M------------");

    /// <summary>Table 4, mold/protozoan/coelenterate mitochondrial and Mycoplasma/Spiroplasma.</summary>
    public static readonly GeneticCode MoldMitochondrial = new(4, "Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM------**-------M------------MMMM---------------M------------");

    /// <summary>Table 11, bacterial, archaeal and plant plastid.</summary>
    public static readonly GeneticCode Bacterial = new(11, "Bacterial, Archaeal and Plant Plastid",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M------------MMMM---------------M------------");

    /// <summary>All built-in tables.</summary>
    public static IReadOnlyList<GeneticCode> All { get; } = [Standard, VertebrateMitochondrial, MoldMitochondrial, Bacterial];

    /// <summary>Looks a table up by NCBI number.</summary>
    public static GeneticCode FromId(int id) => All.FirstOrDefault(c => c.Id == id)
        ?? throw new ArgumentOutOfRangeException(nameof(id), $"Genetic code {id} is not built in.");

    /// <summary>Amino acid for a codon (three nucleotides, DNA or RNA, any case); <c>X</c> if ambiguous.</summary>
    public char Translate(ReadOnlySpan<char> codon)
    {
        if (codon.Length != 3) throw new ArgumentException("A codon has three bases.", nameof(codon));
        int k = Index(codon);
        return k < 0 ? 'X' : (char)_table[k];
    }

    /// <summary>Whether a codon is a start codon in this table.</summary>
    public bool IsStart(ReadOnlySpan<char> codon) { int k = Index(codon); return k >= 0 && _start[k]; }

    /// <summary>Whether a codon is a stop codon in this table.</summary>
    public bool IsStop(ReadOnlySpan<char> codon) { int k = Index(codon); return k >= 0 && _table[k] == (byte)'*'; }

    /// <summary>Kernel-order index of a codon, or −1 if it contains a non-ACGTU symbol.</summary>
    public static int Index(ReadOnlySpan<char> codon)
    {
        int k = 0;
        foreach (char ch in codon)
        {
            int c = char.ToUpperInvariant(ch) switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' or 'U' => 3, _ => -1 };
            if (c < 0) return -1;
            k = (k << 2) | c;
        }
        return k;
    }

    /// <summary>Kernel-order index of a codon given as bytes, or −1.</summary>
    public static int Index(ReadOnlySpan<byte> codon)
    {
        int k = 0;
        foreach (byte b in codon)
        {
            int c = (b & 0xDF) switch { (byte)'A' => 0, (byte)'C' => 1, (byte)'G' => 2, (byte)'T' or (byte)'U' => 3, _ => -1 };
            if (c < 0) return -1;
            k = (k << 2) | c;
        }
        return k;
    }

    /// <summary>Whether the codon at kernel index <paramref name="index"/> is a start codon.</summary>
    public bool IsStartIndex(int index) => index >= 0 && _start[index];

    /// <summary>Amino acid at kernel index <paramref name="index"/>.</summary>
    public char AminoAcidAt(int index) => (char)_table[index];

    /// <summary>All codons (DNA, upper case) that encode <paramref name="aminoAcid"/> (<c>*</c> for stops).</summary>
    public IReadOnlyList<string> CodonsFor(char aminoAcid)
    {
        aminoAcid = char.ToUpperInvariant(aminoAcid);
        var list = new List<string>();
        for (int k = 0; k < 64; k++)
            if (_table[k] == aminoAcid) list.Add(CodonAt(k));
        return list;
    }

    /// <summary>The DNA codon for kernel index <paramref name="index"/>.</summary>
    public static string CodonAt(int index)
    {
        const string bases = "ACGT";
        return string.Create(3, index, static (s, k) => { s[0] = bases[k >> 4]; s[1] = bases[(k >> 2) & 3]; s[2] = bases[k & 3]; });
    }

    public override string ToString() => $"{Id}: {Name}";
}
