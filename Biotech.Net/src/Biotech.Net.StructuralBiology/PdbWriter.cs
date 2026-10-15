using System.Globalization;
using System.Text;

namespace Biotech.Net.StructuralBiology;

/// <summary>Options for <see cref="PdbWriter"/>.</summary>
public sealed record PdbWriteOptions
{
    /// <summary>Write HEADER/TITLE/CRYST1/HELIX/SHEET/SSBOND records (default true).</summary>
    public bool IncludeHeader { get; init; } = true;

    /// <summary>Write CONECT records (default true).</summary>
    public bool IncludeConnections { get; init; } = true;

    /// <summary>Write alternate-location conformers after the selected atom (default true).</summary>
    public bool IncludeAlternates { get; init; } = true;

    /// <summary>Renumber atom serials from 1 instead of keeping file serials (default false).</summary>
    public bool Renumber { get; init; }
}

/// <summary>Writes structures in fixed-column PDB format (ATOM/HETATM, TER, MODEL/ENDMDL and header records).</summary>
public static class PdbWriter
{
    /// <summary>The structure as PDB text.</summary>
    public static string ToPdb(Structure structure, PdbWriteOptions? options = null)
    {
        var sw = new StringWriter(CultureInfo.InvariantCulture);
        Write(structure, sw, options);
        return sw.ToString();
    }

    /// <summary>Writes the structure to a file.</summary>
    public static void Write(Structure structure, string path, PdbWriteOptions? options = null)
    {
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        Write(structure, w, options);
    }

    /// <summary>Writes the structure to <paramref name="writer"/>.</summary>
    public static void Write(Structure structure, TextWriter writer, PdbWriteOptions? options = null)
    {
        options ??= new PdbWriteOptions();
        var inv = CultureInfo.InvariantCulture;
        if (options.IncludeHeader)
        {
            if (structure.Classification.Length > 0) writer.WriteLine($"HEADER    {structure.Classification,-40}".PadRight(80));
            if (structure.Title.Length > 0) WriteTitle(writer, structure.Title);
            foreach (var h in structure.Helices)
                writer.WriteLine(string.Format(inv, "HELIX  {0,3} {1,3} {2,3} {3,1} {4,4}  {5,3} {6,1} {7,4} {8,2}                               {9,5}    ",
                    h.Serial, h.HelixId, h.StartResidueName, Ch(h.StartChain), h.StartNumber, h.EndResidueName, Ch(h.EndChain), h.EndNumber, h.HelixClass, h.Length));
            foreach (var s in structure.Sheets)
                writer.WriteLine(string.Format(inv, "SHEET  {0,3} {1,3}{2,2} {3,3} {4,1}{5,4}  {6,3} {7,1}{8,4} {9,2}".PadRight(80),
                    s.Strand, s.SheetId, s.StrandCount, s.StartResidueName, Ch(s.StartChain), s.StartNumber, s.EndResidueName, Ch(s.EndChain), s.EndNumber, s.Sense));
            foreach (var b in structure.Disulfides)
                writer.WriteLine(string.Format(inv, "SSBOND {0,3} CYS {1,1} {2,4}    CYS {3,1} {4,4}                          1555   1555 {5,5:F2}  ",
                    b.Serial, Ch(b.ChainA), b.NumberA, Ch(b.ChainB), b.NumberB, double.IsNaN(b.Length) ? 0 : b.Length));
            if (structure.Cell is { } c)
                writer.WriteLine(string.Format(inv, "CRYST1{0,9:F3}{1,9:F3}{2,9:F3}{3,7:F2}{4,7:F2}{5,7:F2} {6,-11}{7,4}          ",
                    c.A, c.B, c.C, c.Alpha, c.Beta, c.Gamma, c.SpaceGroup, c.Z));
        }

        bool multi = structure.Models.Count > 1;
        int serial = 0;
        var serialMap = new Dictionary<int, int>();
        foreach (var model in structure.Models)
        {
            if (multi) writer.WriteLine($"MODEL     {model.Number,4}".PadRight(80));
            foreach (var chain in model.Chains)
            {
                Residue? lastPolymer = null;
                bool terWritten = false;
                foreach (var residue in chain.Residues)
                {
                    if (residue.IsHetero && lastPolymer is not null && !terWritten)
                    {
                        WriteTer(writer, ++serial, lastPolymer);
                        terWritten = true;
                    }
                    foreach (var atom in residue.Atoms)
                    {
                        WriteAtom(writer, atom, residue, chain.Id, Serial(atom));
                        if (options.IncludeAlternates)
                            foreach (var alt in atom.Alternates) WriteAtom(writer, alt, residue, chain.Id, Serial(alt));
                    }
                    if (!residue.IsHetero) lastPolymer = residue;
                }
                if (lastPolymer is not null && !terWritten) WriteTer(writer, ++serial, lastPolymer);
            }
            if (multi) writer.WriteLine("ENDMDL".PadRight(80));
        }

        if (options.IncludeConnections)
        {
            foreach (var (a, partners) in structure.Connections.OrderBy(kv => kv.Key))
            {
                int sa = Map(a);
                for (int i = 0; i < partners.Count; i += 4)
                {
                    var sb = new StringBuilder(string.Format(inv, "CONECT{0,5}", sa));
                    for (int j = i; j < Math.Min(i + 4, partners.Count); j++) sb.Append(string.Format(inv, "{0,5}", Map(partners[j])));
                    writer.WriteLine(sb.ToString().PadRight(80));
                }
            }
        }
        writer.WriteLine("END".PadRight(80));

        int Serial(Atom a)
        {
            int s = options.Renumber || a.Serial <= 0 ? ++serial : a.Serial;
            if (!options.Renumber) serial = Math.Max(serial, s);
            if (a.Serial > 0) serialMap.TryAdd(a.Serial, s);
            return s;
        }

        int Map(int s) => serialMap.TryGetValue(s, out int m) ? m : s;
    }

    /// <summary>Formats one ATOM/HETATM line (80 columns).</summary>
    public static string FormatAtom(Atom atom, Residue residue, string chainId, int serial)
    {
        var inv = CultureInfo.InvariantCulture;
        string name = atom.Name.Length >= 4 || atom.Element.Length == 2 ? atom.Name.PadRight(4) : " " + atom.Name.PadRight(3);
        string charge = atom.Charge == 0 ? "  " : $"{Math.Abs(atom.Charge)}{(atom.Charge > 0 ? '+' : '-')}";
        return string.Format(inv, "{0,-6}{1,5} {2}{3}{4,3} {5}{6,4}{7}   {8,8:F3}{9,8:F3}{10,8:F3}{11,6:F2}{12,6:F2}          {13,2}{14}",
            atom.IsHetero ? "HETATM" : "ATOM", serial % 100000, name, atom.AltLoc, residue.Name, Ch(chainId), residue.Number, residue.InsertionCode,
            atom.Position.X, atom.Position.Y, atom.Position.Z, atom.Occupancy, atom.BFactor, atom.Element, charge);
    }

    private static void WriteAtom(TextWriter w, Atom atom, Residue residue, string chainId, int serial)
        => w.WriteLine(FormatAtom(atom, residue, chainId, serial));

    private static void WriteTer(TextWriter w, int serial, Residue last)
        => w.WriteLine(string.Format(CultureInfo.InvariantCulture, "TER   {0,5}      {1,3} {2}{3,4}{4}", serial % 100000, last.Name,
            Ch(last.Chain?.Id ?? ""), last.Number, last.InsertionCode).PadRight(80));

    private static void WriteTitle(TextWriter w, string title)
    {
        // Word-wrap so that the reader's "join continuation lines with a space" restores the original text.
        int cont = 1, i = 0;
        while (i < title.Length)
        {
            int max = cont == 1 ? 70 : 69;
            int len = Math.Min(max, title.Length - i);
            if (i + len < title.Length)
            {
                int space = title.LastIndexOf(' ', i + len, len);
                if (space > i) len = space - i;
            }
            string part = title.Substring(i, len).Trim();
            w.WriteLine((cont == 1 ? $"TITLE     {part}" : $"TITLE   {cont,2} {part}").PadRight(80));
            i += len;
            while (i < title.Length && title[i] == ' ') i++;
            cont++;
        }
    }

    private static char Ch(string chainId) => chainId.Length == 0 ? ' ' : chainId[0];
}
