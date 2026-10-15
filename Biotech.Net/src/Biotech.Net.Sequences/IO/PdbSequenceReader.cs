using System.Collections.Frozen;
using System.Globalization;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>Where a chain sequence extracted from a PDB file came from.</summary>
public enum PdbSequenceSource
{
    /// <summary>SEQRES records: the full construct, including residues that are disordered in the model.</summary>
    Seqres,

    /// <summary>ATOM/HETATM coordinates of the first model: only residues actually observed.</summary>
    Atom,
}

/// <summary>
/// Extracts per-chain sequences from PDB-format files, either from SEQRES records or from the residues present
/// in the coordinates (like Biopython's <c>pdb-seqres</c> and <c>pdb-atom</c> parsers). Ids have the form
/// <c>1ABC:A</c>.
/// </summary>
public static class PdbSequenceReader
{
    /// <summary>
    /// Three-letter residue names to one-letter codes: the 20 standard amino acids, selenocysteine and pyrrolysine,
    /// common modified residues mapped to their parent (MSE→M, SEP→S, TPO→T, PTR→Y, ...), and DNA/RNA nucleotides.
    /// </summary>
    public static FrozenDictionary<string, char> ResidueCodes { get; } = new Dictionary<string, char>
    {
        ["ALA"] = 'A', ["ARG"] = 'R', ["ASN"] = 'N', ["ASP"] = 'D', ["CYS"] = 'C', ["GLN"] = 'Q', ["GLU"] = 'E',
        ["GLY"] = 'G', ["HIS"] = 'H', ["ILE"] = 'I', ["LEU"] = 'L', ["LYS"] = 'K', ["MET"] = 'M', ["PHE"] = 'F',
        ["PRO"] = 'P', ["SER"] = 'S', ["THR"] = 'T', ["TRP"] = 'W', ["TYR"] = 'Y', ["VAL"] = 'V',
        ["SEC"] = 'U', ["PYL"] = 'O', ["ASX"] = 'B', ["GLX"] = 'Z', ["XLE"] = 'J', ["UNK"] = 'X',
        // Modified amino acids (parent residue).
        ["MSE"] = 'M', ["SEP"] = 'S', ["TPO"] = 'T', ["PTR"] = 'Y', ["HYP"] = 'P', ["MLY"] = 'K', ["M3L"] = 'K',
        ["MLZ"] = 'K', ["ALY"] = 'K', ["KCX"] = 'K', ["LLP"] = 'K', ["CSO"] = 'C', ["CSD"] = 'C', ["CME"] = 'C',
        ["CSS"] = 'C', ["CAS"] = 'C', ["OCS"] = 'C', ["SCY"] = 'C', ["SMC"] = 'C', ["PCA"] = 'E', ["CGU"] = 'E',
        ["NLE"] = 'L', ["MLE"] = 'L', ["MVA"] = 'V', ["ABA"] = 'A', ["AIB"] = 'A', ["DAL"] = 'A', ["SAR"] = 'G',
        ["FME"] = 'M', ["CXM"] = 'M', ["TYS"] = 'Y', ["TRQ"] = 'W', ["HIC"] = 'H', ["NEP"] = 'H', ["MHS"] = 'H',
        ["AGM"] = 'R', ["DSN"] = 'S', ["DTH"] = 'T', ["DCY"] = 'C', ["DLE"] = 'L', ["DVA"] = 'V', ["DPR"] = 'P',
        ["DGL"] = 'E', ["DAS"] = 'D', ["DLY"] = 'K', ["DAR"] = 'R', ["DHI"] = 'H', ["DPN"] = 'F', ["DTR"] = 'W',
        ["DTY"] = 'Y', ["DIL"] = 'I', ["DGN"] = 'Q', ["DSG"] = 'N', ["MEN"] = 'N', ["IAS"] = 'D', ["BMT"] = 'T',
        ["CRO"] = 'G', ["GYS"] = 'G',
        // Nucleotides (DNA, RNA and common modifications).
        ["DA"] = 'A', ["DC"] = 'C', ["DG"] = 'G', ["DT"] = 'T', ["DU"] = 'U', ["DI"] = 'I',
        ["A"] = 'A', ["C"] = 'C', ["G"] = 'G', ["U"] = 'U', ["T"] = 'T', ["I"] = 'I', ["N"] = 'N',
        ["ADE"] = 'A', ["CYT"] = 'C', ["GUA"] = 'G', ["THY"] = 'T', ["URA"] = 'U',
        ["PSU"] = 'U', ["5MU"] = 'U', ["H2U"] = 'U', ["4SU"] = 'U', ["5MC"] = 'C', ["OMC"] = 'C', ["OMG"] = 'G',
        ["2MG"] = 'G', ["7MG"] = 'G', ["M2G"] = 'G', ["1MA"] = 'A', ["MIA"] = 'A', ["1MG"] = 'G', ["YG"] = 'G',
    }.ToFrozenDictionary(StringComparer.OrdinalIgnoreCase);

    private static readonly FrozenSet<string> DnaResidues = new[] { "DA", "DC", "DG", "DT", "DU", "DI" }.ToFrozenSet(StringComparer.OrdinalIgnoreCase);
    private static readonly FrozenSet<string> RnaResidues = new[] { "A", "C", "G", "U", "I", "N", "PSU", "5MU", "H2U", "4SU", "5MC", "OMC", "OMG", "2MG", "7MG", "M2G", "1MA", "MIA", "1MG", "YG" }.ToFrozenSet(StringComparer.OrdinalIgnoreCase);

    /// <summary>One-letter code for a residue name, or <c>X</c> when unknown.</summary>
    public static char ToOneLetter(string residueName) => ResidueCodes.TryGetValue(residueName.Trim(), out var c) ? c : 'X';

    /// <summary>Sequences from SEQRES records, one per chain in file order.</summary>
    public static IReadOnlyList<Sequence> ReadSeqres(TextReader reader)
    {
        string pdbId = "XXXX";
        var chains = new Dictionary<char, List<string>>();
        var order = new List<char>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("HEADER", StringComparison.Ordinal)) pdbId = ParsePdbId(line);
            else if (line.StartsWith("SEQRES", StringComparison.Ordinal) && line.Length > 19)
            {
                char chain = line.Length > 11 ? line[11] : ' ';
                if (!chains.TryGetValue(chain, out var list)) { chains[chain] = list = []; order.Add(chain); }
                list.AddRange(line[19..Math.Min(line.Length, 70)].Split(' ', StringSplitOptions.RemoveEmptyEntries));
            }
            else if (line.StartsWith("ATOM", StringComparison.Ordinal)) break; // SEQRES precedes coordinates
        }
        return order.Select(c => Build(pdbId, c, chains[c], PdbSequenceSource.Seqres)).ToList();
    }

    /// <summary>
    /// Sequences of the residues present in the coordinates of the first model. HETATM residues are included only
    /// when they are known modified residues (so ligands and water are skipped). With <paramref name="fillGaps"/>,
    /// jumps in residue numbering become <c>X</c> (protein) or <c>N</c> (nucleic acid), as Biopython does.
    /// </summary>
    public static IReadOnlyList<Sequence> ReadAtomSequences(TextReader reader, bool fillGaps = true)
    {
        string pdbId = "XXXX";
        var chains = new Dictionary<char, List<(int Number, char Insertion, string Name)>>();
        var order = new List<char>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("HEADER", StringComparison.Ordinal)) { pdbId = ParsePdbId(line); continue; }
            if (line.StartsWith("ENDMDL", StringComparison.Ordinal)) break;
            bool atom = line.StartsWith("ATOM  ", StringComparison.Ordinal);
            bool het = line.StartsWith("HETATM", StringComparison.Ordinal);
            if (!(atom || het) || line.Length < 27) continue;
            string name = line[17..20].Trim();
            if (het && (!ResidueCodes.ContainsKey(name) || name is "UNK")) continue;
            char chain = line[21];
            int number = int.Parse(line.AsSpan(22, 4).Trim(), CultureInfo.InvariantCulture);
            char ins = line[26];
            if (!chains.TryGetValue(chain, out var list)) { chains[chain] = list = []; order.Add(chain); }
            if (list.Count > 0 && list[^1].Number == number && list[^1].Insertion == ins) continue; // same residue, next atom
            list.Add((number, ins, name));
        }
        var result = new List<Sequence>(order.Count);
        foreach (var c in order)
        {
            var residues = chains[c];
            var names = new List<string>(residues.Count);
            bool nucleic = residues.All(r => DnaResidues.Contains(r.Name) || RnaResidues.Contains(r.Name));
            for (int i = 0; i < residues.Count; i++)
            {
                if (fillGaps && i > 0 && residues[i].Insertion == ' ' && residues[i - 1].Insertion == ' ')
                    for (int gap = residues[i - 1].Number + 1; gap < residues[i].Number; gap++) names.Add(nucleic ? "N" : "UNK");
                names.Add(residues[i].Name);
            }
            result.Add(Build(pdbId, c, names, PdbSequenceSource.Atom));
        }
        return result;
    }

    /// <summary>SEQRES sequences of a PDB file (optionally gzip-compressed).</summary>
    public static IReadOnlyList<Sequence> ReadSeqres(string path)
    {
        using var r = BioFile.OpenText(path);
        return ReadSeqres(r);
    }

    /// <summary>Observed-residue sequences of a PDB file (optionally gzip-compressed).</summary>
    public static IReadOnlyList<Sequence> ReadAtomSequences(string path, bool fillGaps = true)
    {
        using var r = BioFile.OpenText(path);
        return ReadAtomSequences(r, fillGaps);
    }

    /// <summary>SEQRES sequences when the file has them, otherwise coordinate-derived sequences.</summary>
    public static IReadOnlyList<Sequence> Read(string path)
    {
        var seqres = ReadSeqres(path);
        return seqres.Count > 0 ? seqres : ReadAtomSequences(path);
    }

    private static string ParsePdbId(string headerLine)
        => headerLine.Length >= 66 ? headerLine[62..66].Trim() : "XXXX";

    private static Sequence Build(string pdbId, char chain, List<string> names, PdbSequenceSource source)
    {
        bool dna = names.Count > 0 && names.All(n => DnaResidues.Contains(n) || n == "N");
        bool rna = !dna && names.Count > 0 && names.All(n => RnaResidues.Contains(n));
        var bytes = new byte[names.Count];
        for (int i = 0; i < bytes.Length; i++)
        {
            char c = ToOneLetter(names[i]);
            if ((dna || rna) && c == 'X') c = 'N';
            if (dna && c == 'U') c = 'T';
            bytes[i] = (byte)c;
        }
        var alphabet = dna ? Alphabet.DnaIupac : rna ? Alphabet.RnaIupac : Alphabet.ProteinExtended;
        string id = $"{pdbId}:{(chain == ' ' ? "_" : chain.ToString())}";
        var seq = Sequence.FromBuffer(id, bytes, alphabet, $"{pdbId} chain {chain} ({source})", validate: false);
        seq.Annotations["chain"] = chain.ToString();
        seq.Annotations["pdb_id"] = pdbId;
        seq.Annotations["source"] = source.ToString();
        return seq;
    }
}
