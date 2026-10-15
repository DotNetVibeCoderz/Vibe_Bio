using System.Globalization;
using System.Text;

namespace Biotech.Net.StructuralBiology;

/// <summary>A <c>loop_</c> table of a CIF data block.</summary>
public sealed class CifLoop
{
    internal CifLoop(List<string> tags) => Tags = tags;

    /// <summary>Full tag names (e.g. <c>_atom_site.Cartn_x</c>).</summary>
    public IReadOnlyList<string> Tags { get; }

    /// <summary>Rows of values, one entry per tag. Missing values are <c>?</c> or <c>.</c>.</summary>
    public List<string[]> Rows { get; } = [];

    /// <summary>Column index of <paramref name="tag"/> (case-insensitive), or −1.</summary>
    public int IndexOf(string tag)
    {
        for (int i = 0; i < Tags.Count; i++) if (Tags[i].Equals(tag, StringComparison.OrdinalIgnoreCase)) return i;
        return -1;
    }

    /// <summary>The category of the loop (e.g. <c>_atom_site</c>).</summary>
    public string Category => Tags.Count == 0 ? "" : Tags[0][..Math.Max(0, Tags[0].IndexOf('.'))];
}

/// <summary>One <c>data_</c> block of a CIF file: single-valued items and loops.</summary>
public sealed class CifBlock
{
    internal CifBlock(string name) => Name = name;

    /// <summary>Block name (text after <c>data_</c>).</summary>
    public string Name { get; }

    /// <summary>Single-valued items by full tag (case-insensitive).</summary>
    public Dictionary<string, string> Items { get; } = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>Loops in file order.</summary>
    public List<CifLoop> Loops { get; } = [];

    /// <summary>The loop whose tags belong to <paramref name="category"/> (e.g. <c>_atom_site</c>), or null.</summary>
    public CifLoop? GetLoop(string category) =>
        Loops.FirstOrDefault(l => l.Category.Equals(category, StringComparison.OrdinalIgnoreCase));

    /// <summary>
    /// The rows of a category whether it was written as a loop or as single items (a one-row category), as
    /// tag → column index plus rows.
    /// </summary>
    public CifLoop? GetCategory(string category)
    {
        var loop = GetLoop(category);
        if (loop is not null) return loop;
        var tags = Items.Keys.Where(k => k.StartsWith(category + ".", StringComparison.OrdinalIgnoreCase)).ToList();
        if (tags.Count == 0) return null;
        var single = new CifLoop(tags);
        single.Rows.Add(tags.Select(t => Items[t]).ToArray());
        return single;
    }

    /// <summary>An item value, or null when absent or missing (<c>?</c>/<c>.</c>).</summary>
    public string? Get(string tag) => Items.TryGetValue(tag, out var v) && !MmCifReader.IsMissing(v) ? v : null;
}

/// <summary>
/// Reads PDBx/mmCIF files: a general CIF tokenizer (quoted values, semicolon text fields, loops) plus a builder
/// that turns the <c>_atom_site</c> loop into a <see cref="Structure"/>. Author numbering/naming
/// (<c>auth_*</c>) is preferred over label fields so results match PDB-format files. Cell, title, helices
/// (<c>_struct_conf</c>) and disulfides (<c>_struct_conn</c>) are read when present.
/// </summary>
public static class MmCifReader
{
    /// <summary>Reads an mmCIF file from disk.</summary>
    public static Structure Read(string path) => Parse(File.ReadAllText(path), Path.GetFileNameWithoutExtension(path));

    /// <summary>Parses mmCIF text; <paramref name="id"/> defaults to the data block name.</summary>
    public static Structure Parse(string text, string? id = null)
    {
        var blocks = ParseCif(text);
        if (blocks.Count == 0) throw new InvalidDataException("No data_ block found.");
        return BuildStructure(blocks[0], id);
    }

    /// <summary>True for the CIF missing-value markers <c>?</c> and <c>.</c>.</summary>
    public static bool IsMissing(string v) => v is "?" or ".";

    /// <summary>Tokenises CIF text into data blocks.</summary>
    public static List<CifBlock> ParseCif(string text)
    {
        var blocks = new List<CifBlock>();
        CifBlock? block = null;
        var tokens = Tokenize(text).GetEnumerator();
        Token? pushed = null;
        while (Next(out var t))
        {
            if (!t.Quoted && t.Text.StartsWith("data_", StringComparison.OrdinalIgnoreCase))
            {
                block = new CifBlock(t.Text[5..]);
                blocks.Add(block);
            }
            else if (!t.Quoted && t.Text.Equals("loop_", StringComparison.OrdinalIgnoreCase))
            {
                block ??= Add();
                var tags = new List<string>();
                while (Next(out var tag))
                {
                    if (!tag.Quoted && tag.Text.StartsWith('_')) tags.Add(tag.Text);
                    else { pushed = tag; break; }
                }
                var loop = new CifLoop(tags);
                block.Loops.Add(loop);
                var row = new List<string>(tags.Count);
                while (Next(out var v))
                {
                    if (!v.Quoted && (v.Text.StartsWith('_') || v.Text.Equals("loop_", StringComparison.OrdinalIgnoreCase)
                        || v.Text.StartsWith("data_", StringComparison.OrdinalIgnoreCase)))
                    {
                        pushed = v;
                        break;
                    }
                    row.Add(v.Text);
                    if (row.Count == tags.Count)
                    {
                        loop.Rows.Add(row.ToArray());
                        row.Clear();
                    }
                }
                if (row.Count != 0) throw new InvalidDataException($"Loop {loop.Category} has a partial row.");
            }
            else if (!t.Quoted && t.Text.StartsWith('_'))
            {
                block ??= Add();
                if (!Next(out var value)) throw new InvalidDataException($"Tag {t.Text} has no value.");
                block.Items[t.Text] = value.Text;
            }
        }
        return blocks;

        CifBlock Add()
        {
            var b = new CifBlock("");
            blocks.Add(b);
            return b;
        }

        bool Next(out Token token)
        {
            if (pushed is { } p) { token = p; pushed = null; return true; }
            if (tokens.MoveNext()) { token = tokens.Current; return true; }
            token = default;
            return false;
        }
    }

    private readonly record struct Token(string Text, bool Quoted);

    private static IEnumerable<Token> Tokenize(string text)
    {
        int i = 0, n = text.Length;
        bool lineStart = true;
        while (i < n)
        {
            char c = text[i];
            if (c == '\n') { lineStart = true; i++; continue; }
            if (c is ' ' or '\t' or '\r') { lineStart = false; i++; continue; }
            if (c == ';' && lineStart)
            {
                // Semicolon text field: runs until a line starting with ';'.
                int start = i + 1;
                int end = text.IndexOf("\n;", start, StringComparison.Ordinal);
                if (end < 0) throw new InvalidDataException("Unterminated semicolon text field.");
                yield return new Token(text[start..end].TrimEnd('\r').Trim('\n', '\r'), true);
                i = end + 2;
                lineStart = false;
                continue;
            }
            lineStart = false;
            if (c == '#')
            {
                while (i < n && text[i] != '\n') i++;
                continue;
            }
            if (c is '\'' or '"')
            {
                // A quote closes only when followed by whitespace or end of input (so "O5'" style names survive).
                int j = i + 1;
                while (j < n && !(text[j] == c && (j + 1 == n || char.IsWhiteSpace(text[j + 1])))) j++;
                yield return new Token(text[(i + 1)..Math.Min(j, n)], true);
                i = j + 1;
                continue;
            }
            int k = i;
            while (k < n && !char.IsWhiteSpace(text[k])) k++;
            yield return new Token(text[i..k], false);
            i = k;
        }
    }

    private static Structure BuildStructure(CifBlock block, string? id)
    {
        var structure = new Structure(id ?? (block.Name.Length > 0 ? block.Name : "structure"))
        {
            Title = block.Get("_struct.title") ?? string.Empty,
            Classification = block.Get("_struct_keywords.pdbx_keywords") ?? string.Empty,
        };
        if (block.Get("_cell.length_a") is { } a)
        {
            structure.Cell = new UnitCell(D(a), D(block.Get("_cell.length_b")), D(block.Get("_cell.length_c")),
                D(block.Get("_cell.angle_alpha")), D(block.Get("_cell.angle_beta")), D(block.Get("_cell.angle_gamma")),
                block.Get("_symmetry.space_group_name_H-M") ?? "", int.TryParse(block.Get("_cell.Z_PDB"), out int z) ? z : 1);
        }

        var atoms = block.GetCategory("_atom_site") ?? throw new InvalidDataException("No _atom_site category.");
        int Col(string name) => atoms.IndexOf("_atom_site." + name);
        int cGroup = Col("group_PDB"), cId = Col("id"), cType = Col("type_symbol"), cAlt = Col("label_alt_id"),
            cX = Col("Cartn_x"), cY = Col("Cartn_y"), cZ = Col("Cartn_z"), cOcc = Col("occupancy"), cB = Col("B_iso_or_equiv"),
            cCharge = Col("pdbx_formal_charge"), cIns = Col("pdbx_PDB_ins_code"), cModel = Col("pdbx_PDB_model_num");
        int cAtom = First(Col("auth_atom_id"), Col("label_atom_id")), cComp = First(Col("auth_comp_id"), Col("label_comp_id")),
            cAsym = First(Col("auth_asym_id"), Col("label_asym_id")), cSeq = First(Col("auth_seq_id"), Col("label_seq_id"));
        int cLabelSeq = Col("label_seq_id");
        if (cX < 0 || cY < 0 || cZ < 0 || cAtom < 0) throw new InvalidDataException("_atom_site lacks coordinates or atom names.");

        Model? model = null;
        int currentModel = int.MinValue;
        var chains = new Dictionary<string, Chain>();
        var residues = new Dictionary<(string, int, char, bool), Residue>();
        foreach (var row in atoms.Rows)
        {
            int modelNo = cModel >= 0 && int.TryParse(row[cModel], out int m) ? m : 1;
            if (model is null || modelNo != currentModel)
            {
                model = structure.Models.FirstOrDefault(x => x.Number == modelNo);
                if (model is null)
                {
                    model = new Model(modelNo);
                    structure.AddModel(model);
                }
                currentModel = modelNo;
                chains.Clear();
                foreach (var ch in model.Chains) chains[ch.Id] = ch;
                residues.Clear();
            }
            bool het = cGroup >= 0 && row[cGroup].Equals("HETATM", StringComparison.OrdinalIgnoreCase);
            string name = row[cAtom];
            string element = cType >= 0 && !IsMissing(row[cType]) ? row[cType] : ElementData.GuessFromAtomName(name.PadLeft(2).PadRight(4), het);
            string alt = cAlt >= 0 ? row[cAlt] : ".";
            string charge = cCharge >= 0 ? row[cCharge] : "?";
            var atom = new Atom(name, element, new Vec3(D(row[cX]), D(row[cY]), D(row[cZ])))
            {
                Serial = cId >= 0 && int.TryParse(row[cId], out int s) ? s : 0,
                AltLoc = IsMissing(alt) ? ' ' : alt[0],
                Occupancy = cOcc >= 0 && !IsMissing(row[cOcc]) ? D(row[cOcc]) : 1.0,
                BFactor = cB >= 0 && !IsMissing(row[cB]) ? D(row[cB]) : 0.0,
                Charge = !IsMissing(charge) && int.TryParse(charge, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int q) ? q : 0,
                IsHetero = het,
            };
            string chainId = cAsym >= 0 ? row[cAsym] : "A";
            string seqText = cSeq >= 0 ? row[cSeq] : "?";
            if (IsMissing(seqText) && cLabelSeq >= 0) seqText = row[cLabelSeq];
            int seq = int.TryParse(seqText, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int sq) ? sq : 0;
            char ins = cIns >= 0 && !IsMissing(row[cIns]) ? row[cIns][0] : ' ';
            string comp = cComp >= 0 ? row[cComp] : "UNK";
            if (!chains.TryGetValue(chainId, out var chain))
            {
                chain = new Chain(chainId);
                chains[chainId] = chain;
                model.AddChain(chain);
            }
            // Non-polymer entities without a sequence number are kept apart from each other by their component name.
            var key = (chainId, seq, ins, het && seq == 0);
            if (!residues.TryGetValue(key, out var residue) || (key.Item4 && residue.Name != comp))
            {
                residue = new Residue(comp, seq, ins, het);
                residues[key] = residue;
                chain.AddResidue(residue);
            }
            residue.AddAtom(atom);
        }

        if (block.GetCategory("_struct_conf") is { } conf)
        {
            int ty = conf.IndexOf("_struct_conf.conf_type_id"), bc = conf.IndexOf("_struct_conf.beg_auth_asym_id"),
                bs = conf.IndexOf("_struct_conf.beg_auth_seq_id"), bn = conf.IndexOf("_struct_conf.beg_auth_comp_id"),
                ec = conf.IndexOf("_struct_conf.end_auth_asym_id"), es = conf.IndexOf("_struct_conf.end_auth_seq_id"),
                en = conf.IndexOf("_struct_conf.end_auth_comp_id"), cls = conf.IndexOf("_struct_conf.pdbx_PDB_helix_class"),
                len = conf.IndexOf("_struct_conf.pdbx_PDB_helix_length"), hid = conf.IndexOf("_struct_conf.pdbx_PDB_helix_id");
            int serial = 0;
            if (ty >= 0 && bc >= 0 && bs >= 0 && ec >= 0 && es >= 0)
                foreach (var r in conf.Rows.Where(r => r[ty].StartsWith("HELX", StringComparison.OrdinalIgnoreCase)))
                    structure.Helices.Add(new HelixRecord(++serial, hid >= 0 ? r[hid] : serial.ToString(CultureInfo.InvariantCulture), r[bc],
                        bn >= 0 ? r[bn] : "", I(r[bs]), r[ec], en >= 0 ? r[en] : "", I(r[es]), cls >= 0 ? I(r[cls], 1) : 1, len >= 0 ? I(r[len]) : 0));
        }
        if (block.GetCategory("_struct_conn") is { } conn)
        {
            int ty = conn.IndexOf("_struct_conn.conn_type_id"), c1 = conn.IndexOf("_struct_conn.ptnr1_auth_asym_id"),
                s1 = conn.IndexOf("_struct_conn.ptnr1_auth_seq_id"), c2 = conn.IndexOf("_struct_conn.ptnr2_auth_asym_id"),
                s2 = conn.IndexOf("_struct_conn.ptnr2_auth_seq_id"), d = conn.IndexOf("_struct_conn.pdbx_dist_value");
            int serial = 0;
            if (ty >= 0 && c1 >= 0 && s1 >= 0 && c2 >= 0 && s2 >= 0)
                foreach (var r in conn.Rows.Where(r => r[ty].Equals("disulf", StringComparison.OrdinalIgnoreCase)))
                    structure.Disulfides.Add(new SsBondRecord(++serial, r[c1], I(r[s1]), r[c2], I(r[s2]), d >= 0 && !IsMissing(r[d]) ? D(r[d]) : double.NaN));
        }
        return structure;

        static int First(int a, int b) => a >= 0 ? a : b;
        static double D(string? s) => s is null || IsMissing(s) ? 0 : double.Parse(s, NumberStyles.Float, CultureInfo.InvariantCulture);
        static int I(string s, int fallback = 0) => int.TryParse(s, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int v) ? v : fallback;
    }

    /// <summary>Writes the first model's atoms as a minimal mmCIF (<c>_atom_site</c> loop only).</summary>
    public static string ToMmCif(Structure structure)
    {
        var sb = new StringBuilder();
        var inv = CultureInfo.InvariantCulture;
        sb.Append("data_").AppendLine(structure.Id.Replace(' ', '_'));
        sb.AppendLine("#");
        sb.AppendLine("loop_");
        foreach (var tag in new[] { "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
                     "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
                     "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num" })
            sb.Append("_atom_site.").AppendLine(tag);
        int serial = 0;
        foreach (var model in structure.Models)
            foreach (var chain in model.Chains)
                foreach (var res in chain.Residues)
                    foreach (var atom in res.Atoms.SelectMany(a => a.Alternates.Prepend(a)))
                    {
                        string name = atom.Name.Contains('\'') ? $"\"{atom.Name}\"" : atom.Name;
                        string ch = chain.Id.Length == 0 ? "." : chain.Id;
                        sb.AppendLine(string.Format(inv, "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9:F3} {10:F3} {11:F3} {12:F2} {13:F2} {14} {15} {16} {17}",
                            atom.IsHetero ? "HETATM" : "ATOM", ++serial, atom.Element, name, atom.AltLoc == ' ' ? "." : atom.AltLoc.ToString(),
                            res.Name, ch, res.Number, res.InsertionCode == ' ' ? "?" : res.InsertionCode.ToString(), atom.Position.X, atom.Position.Y,
                            atom.Position.Z, atom.Occupancy, atom.BFactor, atom.Charge == 0 ? "?" : atom.Charge.ToString(inv), res.Number, ch, model.Number));
                    }
        sb.AppendLine("#");
        return sb.ToString();
    }
}
