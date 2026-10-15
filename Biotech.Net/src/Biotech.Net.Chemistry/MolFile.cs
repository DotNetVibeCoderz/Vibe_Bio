using System.Globalization;
using System.Text;

namespace Biotech.Net.Chemistry;

/// <summary>Options for reading MOL/SDF records.</summary>
public sealed record MolReaderOptions
{
    /// <summary>Fold explicit hydrogen atoms into implicit counts after reading (default true, as RDKit does).</summary>
    public bool RemoveHydrogens { get; init; } = true;

    /// <summary>Re-derive aromaticity from the Kekulé structure (default true).</summary>
    public bool PerceiveAromaticity { get; init; } = true;

    /// <summary>SDF only: skip malformed records instead of throwing (default false).</summary>
    public bool SkipInvalid { get; init; }

    /// <summary>Default options.</summary>
    public static MolReaderOptions Default { get; } = new();
}

/// <summary>
/// Reader and writer for MDL MOL files (V2000 connection tables). Reads the atom block (coordinates, symbol, mass
/// difference, charge code, valence), bond block (types 1–4, wedge/hash stereo) and the <c>M  CHG</c>,
/// <c>M  ISO</c> and <c>M  RAD</c> property lines; implicit hydrogens follow the MDL valence model.
/// </summary>
/// <remarks>V3000 connection tables are not supported. Tetrahedral stereo is not derived from wedges or atom parities;
/// wedge flags are kept on <see cref="Bond.Stereo"/> and written back unchanged.</remarks>
public static class MolFile
{
    /// <summary>Parses a MOL block.</summary>
    /// <exception cref="FormatException">The block is malformed or uses V3000.</exception>
    public static Molecule Parse(string molBlock, MolReaderOptions? options = null)
    {
        using var reader = new StringReader(molBlock);
        var lines = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            lines.Add(line);
            if (line.StartsWith("M  END", StringComparison.Ordinal)) break;
        }
        return Parse(lines, options ?? MolReaderOptions.Default);
    }

    /// <summary>Reads a MOL file from disk.</summary>
    public static Molecule Read(string path, MolReaderOptions? options = null) => Parse(File.ReadAllText(path), options);

    internal static Molecule Parse(IReadOnlyList<string> lines, MolReaderOptions options)
    {
        if (lines.Count < 4) throw new FormatException("MOL block too short (needs header and counts line).");
        string name = lines[0].Trim();
        string counts = lines[3];
        if (counts.Contains("V3000", StringComparison.OrdinalIgnoreCase)) throw new FormatException("V3000 MOL files are not supported.");
        int nAtoms = Int(counts, 0, 3), nBonds = Int(counts, 3, 3);
        if (lines.Count < 4 + nAtoms + nBonds) throw new FormatException($"MOL block truncated: expected {nAtoms} atoms and {nBonds} bonds.");
        var mol = new Molecule(name, nAtoms);
        var coords = new Point3[nAtoms];
        var valenceField = new int[nAtoms];
        for (int i = 0; i < nAtoms; i++)
        {
            string l = lines[4 + i];
            double x = Dbl(l, 0, 10), y = Dbl(l, 10, 10), z = Dbl(l, 20, 10);
            string sym = Field(l, 31, 3).Trim();
            int z0;
            if (sym is "D" or "T") z0 = 1;
            else if (sym is "*" or "R" or "R#" or "A" or "Q" or "L" or "X" or "LP") z0 = 0;
            else if (PeriodicTable.TryGet(sym, out var e) || PeriodicTable.TryGetIgnoreCase(sym, out e)) z0 = e.AtomicNumber;
            else throw new FormatException($"Unknown element '{sym}' on atom {i + 1}.");
            var atom = new Atom(z0) { NoImplicit = z0 == 0 };
            if (sym == "D") atom.Isotope = 2;
            else if (sym == "T") atom.Isotope = 3;
            int massDiff = Int(l, 34, 2);
            if (massDiff != 0 && z0 > 0) atom.Isotope = PeriodicTable.Get(z0).MostCommonIsotope + massDiff;
            int chg = Int(l, 36, 3);
            atom.FormalCharge = chg switch { 1 => 3, 2 => 2, 3 => 1, 5 => -1, 6 => -2, 7 => -3, _ => 0 };
            valenceField[i] = Int(l, 48, 3);
            atom.MapNumber = Int(l, 60, 3);
            mol.AddAtom(atom);
            coords[i] = new Point3(x, y, z);
        }
        mol.Coordinates = coords;
        for (int i = 0; i < nBonds; i++)
        {
            string l = lines[4 + nAtoms + i];
            int a = Int(l, 0, 3) - 1, b = Int(l, 3, 3) - 1, type = Int(l, 6, 3), stereo = Int(l, 9, 3);
            if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) throw new FormatException($"Bond {i + 1} references a missing atom.");
            var order = type switch
            {
                1 => BondOrder.Single,
                2 => BondOrder.Double,
                3 => BondOrder.Triple,
                4 => BondOrder.Aromatic,
                _ => BondOrder.Single, // query bond types 5–8
            };
            int bi = mol.AddBond(a, b, order);
            ref var bond = ref mol.BondAt(bi);
            if (order == BondOrder.Single)
                bond.Stereo = stereo switch { 1 => BondStereo.Up, 6 => BondStereo.Down, 4 => BondStereo.Either, _ => BondStereo.None };
            else if (order == BondOrder.Double && stereo == 3) bond.Stereo = BondStereo.Either;
            if (type == 4)
            {
                mol.AtomAt(a).IsAromatic = true;
                mol.AtomAt(b).IsAromatic = true;
            }
        }
        // Property block.
        bool chargeReset = false, isoReset = false;
        for (int li = 4 + nAtoms + nBonds; li < lines.Count; li++)
        {
            string l = lines[li];
            if (l.StartsWith("M  END", StringComparison.Ordinal)) break;
            if (l.StartsWith("M  CHG", StringComparison.Ordinal) || l.StartsWith("M  ISO", StringComparison.Ordinal) || l.StartsWith("M  RAD", StringComparison.Ordinal))
            {
                bool chg = l[3] == 'C', iso = l[3] == 'I';
                if (chg && !chargeReset) { for (int i = 0; i < nAtoms; i++) mol.AtomAt(i).FormalCharge = 0; chargeReset = true; }
                if (iso && !isoReset) { for (int i = 0; i < nAtoms; i++) if (mol.Atoms[i].AtomicNumber != 1 || mol.Atoms[i].Isotope is not (2 or 3)) mol.AtomAt(i).Isotope = 0; isoReset = true; }
                int n = Int(l, 6, 3);
                for (int k = 0; k < n; k++)
                {
                    int atom = Int(l, 10 + 8 * k, 3) - 1, value = Int(l, 14 + 8 * k, 3);
                    if ((uint)atom >= (uint)nAtoms) throw new FormatException($"Property line references atom {atom + 1}.");
                    if (chg) mol.AtomAt(atom).FormalCharge = value;
                    else if (iso) mol.AtomAt(atom).Isotope = value;
                    else if (value != 0) mol.AtomAt(atom).NoImplicit = true; // radical: no implicit H
                }
            }
        }
        // Hydrogens: valence field fixes the total valence (15 = zero).
        mol.UpdateImplicitHydrogens();
        for (int i = 0; i < nAtoms; i++)
        {
            int v = valenceField[i];
            if (v == 0) continue;
            int total = v == 15 ? 0 : v;
            int explicitValence = 0;
            foreach (int b in mol.NeighborBonds(i)) explicitValence += mol.KekuleOrder(b);
            ref var a = ref mol.AtomAt(i);
            a.NoImplicit = true;
            a.ImplicitHydrogens = 0;
            a.ExplicitHydrogens = Math.Max(0, total - explicitValence);
        }
        if (options.RemoveHydrogens) mol = RemoveHydrogensKeepingName(mol);
        if (options.PerceiveAromaticity) Aromaticity.Perceive(mol);
        return mol;
    }

    private static Molecule RemoveHydrogensKeepingName(Molecule mol)
    {
        bool anyH = false;
        foreach (ref readonly var a in mol.Atoms) if (a.AtomicNumber == 1) { anyH = true; break; }
        if (!anyH) return mol;
        var r = mol.RemoveHydrogens();
        // Atoms whose H were removed keep the right count through ExplicitHydrogens; make implicit H consistent.
        for (int i = 0; i < r.AtomCount; i++)
        {
            ref var a = ref r.AtomAt(i);
            if (!a.NoImplicit && a.ExplicitHydrogens > 0)
            {
                // The valence model would already add these implicitly; keep the explicit count and stop recomputation.
                a.NoImplicit = true;
                a.ImplicitHydrogens = 0;
            }
        }
        return r;
    }

    /// <summary>Writes a V2000 MOL block. Molecules without coordinates get a 2-D depiction (<see cref="Depiction"/>).</summary>
    public static string Write(Molecule mol, string? program = "BiotechNet")
    {
        var sb = new StringBuilder();
        Write(mol, new StringWriter(sb, CultureInfo.InvariantCulture), program);
        return sb.ToString();
    }

    /// <summary>Writes a V2000 MOL block to <paramref name="writer"/>.</summary>
    public static void Write(Molecule mol, TextWriter writer, string? program = "BiotechNet")
    {
        if (mol.AtomCount > 999 || mol.BondCount > 999) throw new NotSupportedException("V2000 supports at most 999 atoms and bonds.");
        var coords = mol.Coordinates ?? Depiction.Compute2DCoordinates(mol.Clone());
        bool is3D = coords.Any(p => Math.Abs(p.Z) > 1e-4);
        var inv = CultureInfo.InvariantCulture;
        writer.Write(mol.Name.Replace('\n', ' ').Replace('\r', ' '));
        writer.Write('\n');
        writer.Write($"  {(program ?? "").PadRight(8)[..8]}{DateTime.UtcNow.ToString("MMddyyHHmm", inv)}{(is3D ? "3D" : "2D")}\n");
        writer.Write('\n');
        writer.Write(string.Create(inv, $"{mol.AtomCount,3}{mol.BondCount,3}  0  0  0  0  0  0  0  0999 V2000\n"));
        var chg = new List<(int, int)>();
        var iso = new List<(int, int)>();
        for (int i = 0; i < mol.AtomCount; i++)
        {
            var a = mol.Atoms[i];
            var p = coords[i];
            string sym = a.AtomicNumber == 0 ? "*" : a.Symbol;
            int valence = 0;
            if (a.AtomicNumber > 0 && a.HydrogenCount != DefaultHydrogens(mol, i))
            {
                int total = mol.ExplicitValence(i) + a.HydrogenCount;
                valence = total == 0 ? 15 : total;
            }
            writer.Write(string.Create(inv, $"{p.X,10:F4}{p.Y,10:F4}{p.Z,10:F4} {sym,-3} 0{0,3}  0  0  0{valence,3}  0  0  0{a.MapNumber,3}  0  0\n"));
            if (a.FormalCharge != 0) chg.Add((i + 1, a.FormalCharge));
            if (a.Isotope != 0) iso.Add((i + 1, a.Isotope));
        }
        for (int b = 0; b < mol.BondCount; b++)
        {
            var bond = mol.Bonds[b];
            int type = mol.KekuleOrder(b);
            int stereo = bond.Stereo switch
            {
                BondStereo.Up when bond.Order == BondOrder.Single => 1,
                BondStereo.Down when bond.Order == BondOrder.Single => 6,
                BondStereo.Either when bond.Order == BondOrder.Single => 4,
                BondStereo.Either when bond.Order == BondOrder.Double => 3,
                _ => 0,
            };
            writer.Write(string.Create(inv, $"{bond.Begin + 1,3}{bond.End + 1,3}{type,3}{stereo,3}\n"));
        }
        WriteProperty(writer, "CHG", chg);
        WriteProperty(writer, "ISO", iso);
        writer.Write("M  END\n");
    }

    private static void WriteProperty(TextWriter w, string tag, List<(int Atom, int Value)> items)
    {
        for (int start = 0; start < items.Count; start += 8)
        {
            var chunk = items.Skip(start).Take(8).ToList();
            w.Write(string.Create(CultureInfo.InvariantCulture, $"M  {tag}{chunk.Count,3}"));
            foreach (var (atom, value) in chunk) w.Write(string.Create(CultureInfo.InvariantCulture, $" {atom,3} {value,3}"));
            w.Write('\n');
        }
    }

    /// <summary>Hydrogens the MDL valence model gives atom <paramref name="atom"/> from its charge and Kekulé bonding.</summary>
    internal static int DefaultHydrogens(Molecule mol, int atom)
    {
        var a = mol.Atoms[atom];
        int v = mol.ExplicitValence(atom);
        foreach (int val in Molecule.AllowedValences(a.AtomicNumber, a.FormalCharge))
            if (val >= v) return val - v;
        return 0;
    }

    internal static string Field(string line, int start, int length)
    {
        if (start >= line.Length) return string.Empty;
        return line.Substring(start, Math.Min(length, line.Length - start));
    }

    internal static int Int(string line, int start, int length)
    {
        var f = Field(line, start, length).Trim();
        return f.Length == 0 ? 0 : int.TryParse(f, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) ? v : throw new FormatException($"Expected an integer in '{line}' at column {start + 1}.");
    }

    private static double Dbl(string line, int start, int length)
    {
        var f = Field(line, start, length).Trim();
        return f.Length == 0 ? 0 : double.TryParse(f, NumberStyles.Float, CultureInfo.InvariantCulture, out double v) ? v : throw new FormatException($"Expected a number in '{line}' at column {start + 1}.");
    }
}

/// <summary>
/// Streaming SD-file reader: records are parsed one at a time (<c>$$$$</c>-terminated MOL blocks followed by
/// <c>&gt;  &lt;TAG&gt;</c> data items, stored in <see cref="Molecule.Properties"/>).
/// </summary>
public static class SdfReader
{
    /// <summary>Lazily reads every record from <paramref name="reader"/>.</summary>
    public static IEnumerable<Molecule> Read(TextReader reader, MolReaderOptions? options = null)
    {
        options ??= MolReaderOptions.Default;
        var block = new List<string>();
        string? line;
        int record = 0;
        while (true)
        {
            block.Clear();
            bool any = false;
            while ((line = reader.ReadLine()) is not null)
            {
                if (line.StartsWith("$$$$", StringComparison.Ordinal)) break;
                block.Add(line);
                any = true;
            }
            if (!any && line is null) yield break;
            record++;
            Molecule? mol = null;
            if (!block.TrueForAll(string.IsNullOrWhiteSpace))
            {
                try
                {
                    mol = ParseRecord(block, options);
                }
                catch (Exception ex) when (ex is FormatException or InvalidOperationException or ArgumentException)
                {
                    if (!options.SkipInvalid) throw new FormatException($"SDF record {record}: {ex.Message}", ex);
                }
            }
            if (mol is not null) yield return mol;
            if (line is null) yield break;
        }
    }

    /// <summary>Lazily reads an SD file from disk (gzip-compressed when the name ends in <c>.gz</c>).</summary>
    public static IEnumerable<Molecule> ReadFile(string path, MolReaderOptions? options = null)
    {
        using var stream = File.OpenRead(path);
        using Stream s = path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
            ? new System.IO.Compression.GZipStream(stream, System.IO.Compression.CompressionMode.Decompress)
            : stream;
        using var reader = new StreamReader(s);
        foreach (var m in Read(reader, options)) yield return m;
    }

    /// <summary>Parses all records of an SDF string.</summary>
    public static IReadOnlyList<Molecule> Parse(string sdf, MolReaderOptions? options = null) => Read(new StringReader(sdf), options).ToList();

    private static Molecule ParseRecord(List<string> lines, MolReaderOptions options)
    {
        int end = lines.FindIndex(l => l.StartsWith("M  END", StringComparison.Ordinal));
        if (end < 0) throw new FormatException("Missing 'M  END'.");
        var mol = MolFile.Parse(lines.GetRange(0, end + 1), options);
        // Data items.
        int i = end + 1;
        while (i < lines.Count)
        {
            string l = lines[i];
            if (l.StartsWith('>'))
            {
                int lt = l.IndexOf('<'), gt = l.IndexOf('>', lt + 1);
                string tag = lt >= 0 && gt > lt ? l[(lt + 1)..gt] : l[1..].Trim();
                var value = new StringBuilder();
                i++;
                while (i < lines.Count && lines[i].Length > 0 && !lines[i].StartsWith('>'))
                {
                    if (value.Length > 0) value.Append('\n');
                    value.Append(lines[i]);
                    i++;
                }
                mol.Properties[tag] = value.ToString();
            }
            else i++;
        }
        return mol;
    }
}

/// <summary>SD-file writer: MOL blocks (<see cref="MolFile.Write(Molecule, TextWriter, string?)"/>) plus data items from <see cref="Molecule.Properties"/>.</summary>
public static class SdfWriter
{
    /// <summary>Writes molecules to <paramref name="writer"/>; returns the number written.</summary>
    /// <param name="writer">Destination.</param>
    /// <param name="molecules">Molecules to write (streamed).</param>
    /// <param name="properties">Property names to write (default: each molecule's own properties).</param>
    public static int Write(TextWriter writer, IEnumerable<Molecule> molecules, IReadOnlyList<string>? properties = null)
    {
        int n = 0;
        foreach (var m in molecules)
        {
            WriteRecord(writer, m, properties);
            n++;
        }
        return n;
    }

    /// <summary>Writes one record.</summary>
    public static void WriteRecord(TextWriter writer, Molecule mol, IReadOnlyList<string>? properties = null)
    {
        MolFile.Write(mol, writer);
        IEnumerable<KeyValuePair<string, string>> items = properties is null
            ? mol.Properties
            : properties.Where(p => mol.Properties.ContainsKey(p)).Select(p => new KeyValuePair<string, string>(p, mol.Properties[p]));
        foreach (var (key, value) in items)
        {
            writer.Write($">  <{key}>\n");
            foreach (var line in value.Replace("\r\n", "\n").Split('\n'))
            {
                writer.Write(line.Length == 0 ? " " : line);
                writer.Write('\n');
            }
            writer.Write('\n');
        }
        writer.Write("$$$$\n");
    }

    /// <summary>Writes molecules to an SD file.</summary>
    public static int WriteFile(string path, IEnumerable<Molecule> molecules, IReadOnlyList<string>? properties = null)
    {
        using var w = new StreamWriter(path, false, new UTF8Encoding(false));
        return Write(w, molecules, properties);
    }

    /// <summary>The SDF text for a set of molecules.</summary>
    public static string ToSdf(IEnumerable<Molecule> molecules)
    {
        var sw = new StringWriter(CultureInfo.InvariantCulture);
        Write(sw, molecules);
        return sw.ToString();
    }
}
