using System.Collections.Concurrent;
using System.Globalization;
using System.Text;
using Biotech.Net.Native;

namespace Biotech.Net.Alignments;

/// <summary>
/// A residue substitution matrix (BLOSUM, PAM, NUC.4.4 or custom) with integer scores. Lookups are case-insensitive.
/// </summary>
/// <remarks>
/// Symbols missing from the matrix fall back, in order, to an alias (<c>U</c> → <c>T</c> for nucleotide matrices),
/// the <c>X</c> row/column, the <c>N</c> row/column, and finally the matrix minimum. <see cref="KernelMatrix"/> exposes
/// the 128×128 ASCII-indexed form the alignment kernels consume.
/// </remarks>
public sealed class SubstitutionMatrix
{
    private static readonly string[] BuiltInNames = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250", "NUC.4.4"];
    private static readonly ConcurrentDictionary<string, SubstitutionMatrix> Cache = new(StringComparer.OrdinalIgnoreCase);

    private readonly int[,] _scores;
    private readonly int[] _kernel;

    /// <summary>Creates a matrix from its symbols and a square score table (row/column order = <paramref name="symbols"/>).</summary>
    public SubstitutionMatrix(string name, string symbols, int[,] scores)
    {
        ArgumentNullException.ThrowIfNull(name);
        ArgumentNullException.ThrowIfNull(symbols);
        ArgumentNullException.ThrowIfNull(scores);
        symbols = symbols.ToUpperInvariant();
        if (scores.GetLength(0) != symbols.Length || scores.GetLength(1) != symbols.Length)
            throw new ArgumentException("The score table must be square with one row and column per symbol.", nameof(scores));
        if (symbols.Distinct().Count() != symbols.Length) throw new ArgumentException("Duplicate symbols.", nameof(symbols));
        foreach (char c in symbols) if (c > 127) throw new ArgumentException("Symbols must be ASCII.", nameof(symbols));
        Name = name;
        Symbols = symbols;
        _scores = (int[,])scores.Clone();
        int min = int.MaxValue, max = int.MinValue;
        foreach (int v in _scores) { min = Math.Min(min, v); max = Math.Max(max, v); }
        MinScore = symbols.Length == 0 ? 0 : min;
        MaxScore = symbols.Length == 0 ? 0 : max;
        _kernel = BuildKernel();
    }

    private SubstitutionMatrix(string name, string symbols, int[,] scores, int[] kernel, int min, int max)
    {
        Name = name; Symbols = symbols; _scores = scores; _kernel = kernel; MinScore = min; MaxScore = max;
    }

    /// <summary>Matrix name, e.g. <c>BLOSUM62</c>.</summary>
    public string Name { get; }

    /// <summary>Upper-case symbols in row/column order.</summary>
    public string Symbols { get; }

    /// <summary>Lowest score in the table.</summary>
    public int MinScore { get; }

    /// <summary>Highest score in the table.</summary>
    public int MaxScore { get; }

    /// <summary>True for nucleotide matrices (symbols include A, C, G, T and few others).</summary>
    public bool IsNucleotide => Symbols.Contains('T') && !Symbols.Contains('L') && !Symbols.Contains('E');

    /// <summary>Score for aligning <paramref name="a"/> with <paramref name="b"/> (case-insensitive, with fallbacks).</summary>
    public int this[char a, char b] => Score(a, b);

    /// <summary>Score for aligning <paramref name="a"/> with <paramref name="b"/> (case-insensitive, with fallbacks).</summary>
    public int Score(char a, char b) => a > 127 || b > 127 ? MinScore : _kernel[(a << 7) | b];

    /// <summary>Score for two ASCII residue bytes.</summary>
    public int Score(byte a, byte b) => _kernel[((a & 0x7F) << 7) | (b & 0x7F)];

    /// <summary>The 128×128 ASCII-indexed table (row = first sequence) used by <see cref="BioKernels"/>. Do not modify.</summary>
    public ReadOnlySpan<int> KernelMatrix => _kernel;

    internal int[] KernelArray => _kernel;

    /// <summary>True when <c>M[a,b] = M[b,a]</c> for every pair.</summary>
    public bool IsSymmetric
    {
        get
        {
            for (int i = 0; i < Symbols.Length; i++)
                for (int j = 0; j < i; j++) if (_scores[i, j] != _scores[j, i]) return false;
            return true;
        }
    }

    // ---- Built-ins ------------------------------------------------------------------------------

    /// <summary>BLOSUM45 (distant homologues).</summary>
    public static SubstitutionMatrix Blosum45 => Load("BLOSUM45");
    /// <summary>BLOSUM50 (the matrix of the Durbin et al. textbook examples).</summary>
    public static SubstitutionMatrix Blosum50 => Load("BLOSUM50");
    /// <summary>BLOSUM62, the BLASTP default.</summary>
    public static SubstitutionMatrix Blosum62 => Load("BLOSUM62");
    /// <summary>BLOSUM80 (close homologues).</summary>
    public static SubstitutionMatrix Blosum80 => Load("BLOSUM80");
    /// <summary>BLOSUM90 (very close homologues).</summary>
    public static SubstitutionMatrix Blosum90 => Load("BLOSUM90");
    /// <summary>PAM30 (short, near-identical peptides).</summary>
    public static SubstitutionMatrix Pam30 => Load("PAM30");
    /// <summary>PAM70.</summary>
    public static SubstitutionMatrix Pam70 => Load("PAM70");
    /// <summary>PAM250 (distant homologues).</summary>
    public static SubstitutionMatrix Pam250 => Load("PAM250");
    /// <summary>NUC.4.4 (EDNAFULL), the EMBOSS nucleotide matrix with IUPAC ambiguity codes.</summary>
    public static SubstitutionMatrix Nuc44 => Load("NUC.4.4");

    /// <summary>Names accepted by <see cref="Load"/>.</summary>
    public static IReadOnlyList<string> AvailableNames => BuiltInNames;

    /// <summary>Loads a built-in matrix by name (case-insensitive; <c>NUC44</c> and <c>EDNAFULL</c> alias <c>NUC.4.4</c>).</summary>
    /// <exception cref="KeyNotFoundException">No built-in matrix has that name.</exception>
    public static SubstitutionMatrix Load(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        string key = name.Trim().ToUpperInvariant() switch
        {
            "NUC44" or "NUC4.4" or "EDNAFULL" or "NUC.4.4" => "NUC.4.4",
            var other => other,
        };
        if (!BuiltInNames.Contains(key)) throw new KeyNotFoundException($"Unknown substitution matrix '{name}'. Available: {string.Join(", ", BuiltInNames)}.");
        return Cache.GetOrAdd(key, static k =>
        {
            using var stream = typeof(SubstitutionMatrix).Assembly.GetManifestResourceStream($"Biotech.Net.Alignments.Matrices.{k}.txt")
                ?? throw new InvalidOperationException($"Embedded matrix {k} is missing.");
            using var reader = new StreamReader(stream);
            return Read(reader, k);
        });
    }

    /// <summary>A matrix scoring <paramref name="match"/> for identical symbols and <paramref name="mismatch"/> otherwise, over any residues.</summary>
    public static SubstitutionMatrix MatchMismatch(int match, int mismatch, string symbols = "ACGT")
    {
        symbols = symbols.ToUpperInvariant();
        var t = new int[symbols.Length, symbols.Length];
        for (int i = 0; i < symbols.Length; i++)
            for (int j = 0; j < symbols.Length; j++) t[i, j] = i == j ? match : mismatch;
        return new SubstitutionMatrix($"match={match},mismatch={mismatch}", symbols, t, BioKernels.MatchMismatchMatrix(match, mismatch),
            Math.Min(match, mismatch), Math.Max(match, mismatch));
    }

    // ---- NCBI format ----------------------------------------------------------------------------

    /// <summary>Parses an NCBI/EMBOSS-format matrix: <c>#</c> comments, a header row of symbols, then one row per symbol.</summary>
    public static SubstitutionMatrix Parse(string text, string name = "custom") => Read(new StringReader(text), name);

    /// <summary>Reads an NCBI/EMBOSS-format matrix file.</summary>
    public static SubstitutionMatrix FromFile(string path)
    {
        using var reader = new StreamReader(path);
        return Read(reader, Path.GetFileNameWithoutExtension(path));
    }

    /// <summary>Reads an NCBI/EMBOSS-format matrix.</summary>
    public static SubstitutionMatrix Read(TextReader reader, string name = "custom")
    {
        string? header = null;
        var rows = new List<(char Symbol, int[] Values)>();
        string? line;
        int lineNo = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            string t = line.Trim();
            if (t.Length == 0 || t[0] == '#') continue;
            var parts = t.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (header is null)
            {
                if (parts.Any(p => p.Length != 1)) throw new FormatException($"Line {lineNo}: expected a header row of single-letter symbols.");
                header = string.Concat(parts).ToUpperInvariant();
                continue;
            }
            if (parts[0].Length != 1) throw new FormatException($"Line {lineNo}: a row must start with its symbol.");
            if (parts.Length - 1 != header.Length) throw new FormatException($"Line {lineNo}: expected {header.Length} scores, found {parts.Length - 1}.");
            var values = new int[header.Length];
            for (int k = 0; k < values.Length; k++)
            {
                if (!double.TryParse(parts[k + 1], NumberStyles.Float, CultureInfo.InvariantCulture, out double v) || v != Math.Round(v))
                    throw new FormatException($"Line {lineNo}: '{parts[k + 1]}' is not an integer score.");
                values[k] = (int)v;
            }
            rows.Add((char.ToUpperInvariant(parts[0][0]), values));
        }
        if (header is null) throw new FormatException("No header row found.");
        if (rows.Count != header.Length) throw new FormatException($"Expected {header.Length} rows, found {rows.Count}.");
        var table = new int[header.Length, header.Length];
        for (int i = 0; i < header.Length; i++)
        {
            int r = rows.FindIndex(x => x.Symbol == header[i]);
            if (r < 0) throw new FormatException($"Missing row for symbol '{header[i]}'.");
            for (int j = 0; j < header.Length; j++) table[i, j] = rows[r].Values[j];
        }
        return new SubstitutionMatrix(name, header, table);
    }

    /// <summary>Writes the matrix in NCBI format.</summary>
    public string ToNcbiFormat()
    {
        var sb = new StringBuilder();
        sb.Append("# ").AppendLine(Name);
        sb.Append(' ');
        foreach (char c in Symbols) sb.Append("  ").Append(c);
        sb.AppendLine();
        for (int i = 0; i < Symbols.Length; i++)
        {
            sb.Append(Symbols[i]);
            for (int j = 0; j < Symbols.Length; j++) sb.Append(_scores[i, j].ToString(CultureInfo.InvariantCulture).PadLeft(3));
            sb.AppendLine();
        }
        return sb.ToString();
    }

    /// <inheritdoc />
    public override string ToString() => Name;

    private int[] BuildKernel()
    {
        var index = new int[128];
        for (int c = 0; c < 128; c++) index[c] = Resolve((char)c);
        var k = new int[BioKernels.SubstitutionMatrixLength];
        for (int a = 0; a < 128; a++)
            for (int b = 0; b < 128; b++)
                k[(a << 7) | b] = index[a] < 0 || index[b] < 0 ? MinScore : _scores[index[a], index[b]];
        return k;
    }

    private int Resolve(char c)
    {
        char u = char.ToUpperInvariant(c);
        int i = Symbols.IndexOf(u);
        if (i >= 0) return i;
        if (!char.IsAsciiLetter(u) && u != '*') return -1;
        if (u == 'U' && IsNucleotide && Symbols.Contains('T')) return Symbols.IndexOf('T');
        if (u != '*')
        {
            if ((i = Symbols.IndexOf('X')) >= 0) return i;
            if (IsNucleotide && (i = Symbols.IndexOf('N')) >= 0) return i;
        }
        return -1;
    }
}
