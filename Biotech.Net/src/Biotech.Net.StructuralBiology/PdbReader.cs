using System.Globalization;

namespace Biotech.Net.StructuralBiology;

/// <summary>
/// Reads the fixed-column PDB format: ATOM/HETATM (with alternate locations, insertion codes, element and charge
/// columns), MODEL/ENDMDL, TER, HEADER/TITLE, CRYST1, HELIX, SHEET, SSBOND and CONECT. Other records are ignored.
/// </summary>
/// <remarks>
/// Lines shorter than 80 columns are accepted. When the element columns are blank the element is inferred from the
/// atom-name alignment. Residues are keyed by chain, number and insertion code, so waters listed after a TER record
/// with the same chain id join that chain (as in Biopython). For atoms with alternate locations the conformer with
/// the highest occupancy is kept in <see cref="Residue.Atoms"/> and the others are listed in <see cref="Atom.Alternates"/>.
/// </remarks>
public static class PdbReader
{
    /// <summary>Reads a PDB file from disk.</summary>
    public static Structure Read(string path)
    {
        using var reader = new StreamReader(path);
        return Read(reader, Path.GetFileNameWithoutExtension(path));
    }

    /// <summary>Parses PDB-format text.</summary>
    public static Structure Parse(string text, string id = "structure")
    {
        using var reader = new StringReader(text);
        return Read(reader, id);
    }

    /// <summary>Reads PDB records from <paramref name="reader"/>.</summary>
    public static Structure Read(TextReader reader, string id = "structure")
    {
        var structure = new Structure(id);
        Model? model = null;
        var builder = new ModelBuilder();
        int modelCount = 0;
        string? line;
        int lineNo = 0;
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            if (line.Length < 3) continue;
            string record = line.Length >= 6 ? line[..6] : line.PadRight(6);
            switch (record)
            {
                case "ATOM  ":
                case "HETATM":
                    if (model is null)
                    {
                        model = new Model(++modelCount);
                        structure.AddModel(model);
                        builder.Reset(model);
                    }
                    try { var (atom, res) = ParseAtom(line, record == "HETATM"); builder.Add(atom, res); }
                    catch (FormatException ex) { throw new InvalidDataException($"PDB line {lineNo}: {ex.Message}", ex); }
                    break;
                case "MODEL ":
                    {
                        int n = ParseInt(line, 10, 4, ++modelCount);
                        modelCount = Math.Max(modelCount, n);
                        model = new Model(n);
                        structure.AddModel(model);
                        builder.Reset(model);
                        break;
                    }
                case "ENDMDL":
                    model = null;
                    break;
                case "HEADER":
                    structure.Classification = Field(line, 10, 40);
                    break;
                case "TITLE ":
                    {
                        string part = Field(line, 10, 70);
                        structure.Title = structure.Title.Length == 0 ? part : structure.Title + " " + part;
                        break;
                    }
                case "CRYST1":
                    structure.Cell = new UnitCell(ParseDouble(line, 6, 9), ParseDouble(line, 15, 9), ParseDouble(line, 24, 9),
                        ParseDouble(line, 33, 7), ParseDouble(line, 40, 7), ParseDouble(line, 47, 7), Field(line, 55, 11), ParseInt(line, 66, 4, 1));
                    break;
                case "HELIX ":
                    structure.Helices.Add(new HelixRecord(ParseInt(line, 7, 3, 0), Field(line, 11, 3), Field(line, 19, 1), Field(line, 15, 3),
                        ParseInt(line, 21, 4, 0), Field(line, 31, 1), Field(line, 27, 3), ParseInt(line, 33, 4, 0), ParseInt(line, 38, 2, 1), ParseInt(line, 71, 5, 0)));
                    break;
                case "SHEET ":
                    structure.Sheets.Add(new SheetRecord(ParseInt(line, 7, 3, 0), Field(line, 11, 3), ParseInt(line, 14, 2, 0), Field(line, 21, 1),
                        Field(line, 17, 3), ParseInt(line, 22, 4, 0), Field(line, 32, 1), Field(line, 28, 3), ParseInt(line, 33, 4, 0), ParseInt(line, 38, 2, 0)));
                    break;
                case "SSBOND":
                    structure.Disulfides.Add(new SsBondRecord(ParseInt(line, 7, 3, 0), Field(line, 15, 1), ParseInt(line, 17, 4, 0),
                        Field(line, 29, 1), ParseInt(line, 31, 4, 0), ParseDouble(line, 73, 5, double.NaN)));
                    break;
                case "CONECT":
                    {
                        int a = ParseInt(line, 6, 5, -1);
                        if (a < 0) break;
                        for (int col = 11; col + 5 <= Math.Min(line.Length, 31); col += 5)
                        {
                            int b = ParseInt(line, col, 5, -1);
                            if (b >= 0) structure.AddConnection(a, b);
                        }
                        break;
                    }
            }
        }
        return structure;
    }

    internal static (Atom Atom, (string Name, string Chain, int Number, char ICode) Residue) ParseAtom(string line, bool hetero)
    {
        if (line.Length < 54) throw new FormatException("ATOM/HETATM record is shorter than 54 columns.");
        string nameField = Slice(line, 12, 4);
        string element = Field(line, 76, 2);
        if (element.Length == 0 || !char.IsLetter(element[0])) element = ElementData.GuessFromAtomName(nameField, hetero);
        var atom = new Atom(nameField, element, new Vec3(ParseDouble(line, 30, 8), ParseDouble(line, 38, 8), ParseDouble(line, 46, 8)))
        {
            Serial = ParseInt(line, 6, 5, 0),
            AltLoc = line[16],
            Occupancy = ParseDouble(line, 54, 6, 1.0),
            BFactor = ParseDouble(line, 60, 6, 0.0),
            Charge = ParseCharge(Field(line, 78, 2)),
            IsHetero = hetero,
        };
        string resName = Field(line, 17, 3);
        string chain = Slice(line, 21, 1);
        int resSeq = ParseInt(line, 22, 4, 0);
        char icode = line.Length > 26 ? line[26] : ' ';
        return (atom, (resName, chain.Trim(), resSeq, icode));
    }

    private sealed class ModelBuilder
    {
        private Model _model = null!;
        private readonly Dictionary<string, Chain> _chains = [];
        private readonly Dictionary<(string, int, char), Residue> _residues = [];

        public void Reset(Model model)
        {
            _model = model;
            _chains.Clear();
            _residues.Clear();
        }

        public void Add(Atom atom, (string Name, string Chain, int Number, char ICode) res)
        {
            var (resName, chainId, number, icode) = res;
            if (!_chains.TryGetValue(chainId, out var chain))
            {
                chain = new Chain(chainId);
                _chains[chainId] = chain;
                _model.AddChain(chain);
            }
            if (!_residues.TryGetValue((chainId, number, icode), out var residue))
            {
                residue = new Residue(resName, number, icode, atom.IsHetero);
                _residues[(chainId, number, icode)] = residue;
                chain.AddResidue(residue);
            }
            residue.AddAtom(atom);
        }
    }

    internal static string Slice(string line, int start, int length)
    {
        if (start >= line.Length) return string.Empty;
        return line.Substring(start, Math.Min(length, line.Length - start));
    }

    internal static string Field(string line, int start, int length) => Slice(line, start, length).Trim();

    private static double ParseDouble(string line, int start, int length, double? fallback = null)
    {
        string s = Field(line, start, length);
        if (s.Length == 0)
            return fallback ?? throw new FormatException($"Missing number at column {start + 1}.");
        if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out double v)) return v;
        return fallback ?? throw new FormatException($"Bad number '{s}' at column {start + 1}.");
    }

    private static int ParseInt(string line, int start, int length, int fallback)
    {
        string s = Field(line, start, length);
        return int.TryParse(s, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) ? v : fallback;
    }

    private static int ParseCharge(string s)
    {
        if (s.Length == 0) return 0;
        if (s.Length == 2 && char.IsDigit(s[0])) return (s[1] == '-' ? -1 : 1) * (s[0] - '0');
        if (s.Length == 2 && char.IsDigit(s[1])) return (s[0] == '-' ? -1 : 1) * (s[1] - '0');
        return int.TryParse(s, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int v) ? v : 0;
    }
}
