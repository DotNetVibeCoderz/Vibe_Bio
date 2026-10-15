using System.Globalization;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>Annotation keys used by the GenBank and EMBL readers/writers in <see cref="Sequence.Annotations"/>.</summary>
public static class RecordKeys
{
    public const string Locus = "locus";
    public const string MoleculeType = "molecule_type";
    public const string Topology = "topology";
    public const string Division = "division";
    public const string Date = "date";
    public const string Accession = "accession";
    public const string Version = "version";
    public const string Keywords = "keywords";
    public const string Source = "source";
    public const string Organism = "organism";
    public const string Taxonomy = "taxonomy";
    public const string Comment = "comment";
    public const string Definition = "definition";

    /// <summary>Key of a reference sub-field, e.g. <c>reference1.title</c>.</summary>
    public static string Reference(int number, string field) => $"reference{number}.{field}";
}

/// <summary>
/// Streaming GenBank flat-file reader: header (LOCUS, DEFINITION, ACCESSION, VERSION, KEYWORDS, SOURCE/ORGANISM,
/// REFERENCE, COMMENT), the FEATURES table with qualifiers and full location syntax, and ORIGIN.
/// </summary>
/// <remarks>
/// The sequence id is the VERSION (accession.version) when present, else the accession, else the LOCUS name,
/// matching Biopython. The definition, without its final period, becomes the description. ORIGIN residues are
/// upper-cased: flat files write them in lower case, and Biopython's <c>_FeatureConsumer.sequence</c> upper-cases
/// them too, so <c>SeqIO.read(..., "genbank").seq</c> and this reader agree.
/// </remarks>
public static class GenBankReader
{
    /// <summary>Reads every record (<c>//</c>-terminated) of <paramref name="reader"/>.</summary>
    public static IEnumerable<Sequence> Read(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (!line.StartsWith("LOCUS", StringComparison.Ordinal)) continue;
            yield return ReadRecord(line, reader);
        }
    }

    /// <summary>Reads a GenBank file (optionally gzip-compressed) lazily.</summary>
    public static IEnumerable<Sequence> Read(string path)
    {
        using var reader = BioFile.OpenText(path);
        foreach (var s in Read(reader)) yield return s;
    }

    /// <summary>Parses GenBank text held in memory.</summary>
    public static IReadOnlyList<Sequence> Parse(string text) => Read(new StringReader(text)).ToList();

    private static Sequence ReadRecord(string locusLine, TextReader reader)
    {
        var ann = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        ParseLocus(locusLine, ann);
        var residues = new ByteBuilder(4096);
        var features = new FeatureTableBuilder();
        string? keyword = null;
        var value = new StringBuilder();
        int referenceNo = 0;
        bool inFeatures = false, inOrigin = false;
        string? line;

        void Commit()
        {
            if (keyword is null) return;
            string v = value.ToString().Trim();
            switch (keyword)
            {
                case "DEFINITION": ann[RecordKeys.Definition] = v; break;
                case "ACCESSION": ann[RecordKeys.Accession] = v.Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault() ?? v; break;
                case "VERSION": ann[RecordKeys.Version] = v.Split(' ', StringSplitOptions.RemoveEmptyEntries).FirstOrDefault() ?? v; break;
                case "KEYWORDS": ann[RecordKeys.Keywords] = v.TrimEnd('.'); break;
                case "SOURCE": ann[RecordKeys.Source] = v; break;
                case "ORGANISM":
                    {
                        // First line is the organism, following lines the ';'-separated lineage.
                        var lines = value.ToString().Split('\n');
                        ann[RecordKeys.Organism] = lines[0].Trim();
                        string lineage = string.Join(" ", lines.Skip(1).Select(l => l.Trim())).Trim().TrimEnd('.');
                        if (lineage.Length > 0) ann[RecordKeys.Taxonomy] = lineage;
                        break;
                    }
                case "COMMENT": ann[RecordKeys.Comment] = value.ToString().TrimEnd(); break;
                case "REFERENCE": ann[RecordKeys.Reference(referenceNo, "location")] = v; break;
                case "AUTHORS" or "CONSRTM" or "TITLE" or "JOURNAL" or "PUBMED" or "REMARK" or "MEDLINE":
                    if (referenceNo > 0) ann[RecordKeys.Reference(referenceNo, keyword.ToLowerInvariant())] = v;
                    break;
                default:
                    ann[keyword.ToLowerInvariant()] = v;
                    break;
            }
            keyword = null;
            value.Clear();
        }

        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("//", StringComparison.Ordinal)) break;
            if (inOrigin)
            {
                residues.AppendResidues(line, skipDigits: true, upperCase: true);
                continue;
            }
            if (inFeatures)
            {
                if (line.Length > 0 && line[0] != ' ')
                {
                    inFeatures = false; // a new top-level keyword ends the table
                }
                else
                {
                    features.AddLine(line.Length > 5 ? line.AsSpan(5) : ReadOnlySpan<char>.Empty);
                    continue;
                }
            }
            if (line.Length == 0) continue;
            bool topLevel = line[0] != ' ';
            bool subKeyword = !topLevel && !line.AsSpan(0, Math.Min(12, line.Length)).IsWhiteSpace();
            if (topLevel || subKeyword)
            {
                Commit();
                string kw = (line.Length >= 12 ? line[..12] : line).Trim();
                int sp = kw.IndexOf(' ');
                if (sp > 0) kw = kw[..sp];
                string rest = line.Length > 12 ? line[12..] : string.Empty;
                switch (kw)
                {
                    case "FEATURES": inFeatures = true; continue;
                    case "ORIGIN": inOrigin = true; continue;
                    case "REFERENCE": referenceNo++; break;
                    case "BASE": continue; // BASE COUNT is derived data
                }
                keyword = kw;
                value.Append(rest.Trim());
                continue;
            }
            // Continuation line of the current keyword.
            if (keyword is null) continue;
            value.Append(keyword is "ORGANISM" or "COMMENT" ? '\n' : ' ');
            value.Append(keyword == "COMMENT" ? line.Length > 12 ? line[12..].TrimEnd() : line.Trim() : line.Trim());
        }
        Commit();

        var bytes = residues.ToArray();
        bool protein = ann.TryGetValue("_unit", out var unit) && unit == "aa";
        ann.Remove("_unit");
        Alphabet alphabet = protein ? Alphabet.ProteinExtended
            : ann.TryGetValue(RecordKeys.MoleculeType, out var mol) && mol.Contains("RNA", StringComparison.OrdinalIgnoreCase) && !Alphabet.DnaIupac.IsValid(bytes)
                ? Alphabet.RnaIupac
                : Alphabet.DnaIupac;
        if (!alphabet.IsValid(bytes)) alphabet = Alphabet.Detect(bytes);

        string id = ann.GetValueOrDefault(RecordKeys.Version) ?? ann.GetValueOrDefault(RecordKeys.Accession) ?? ann.GetValueOrDefault(RecordKeys.Locus) ?? "unknown";
        string description = ann.GetValueOrDefault(RecordKeys.Definition) ?? string.Empty;
        if (description.EndsWith('.')) description = description[..^1];
        var seq = Sequence.FromBuffer(id, bytes, alphabet, description, validate: false);
        foreach (var kv in ann) seq.Annotations[kv.Key] = kv.Value;
        foreach (var f in features.Finish()) seq.Features.Add(f);
        return seq;
    }

    private static readonly HashSet<string> Divisions =
    [
        "PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV", "CON", "TSA", "UNK",
    ];

    private static void ParseLocus(string line, Dictionary<string, string> ann)
    {
        var tokens = line[5..].Split(' ', StringSplitOptions.RemoveEmptyEntries);
        if (tokens.Length == 0) return;
        ann[RecordKeys.Locus] = tokens[0];
        int unitAt = Array.FindIndex(tokens, t => t is "bp" or "aa" or "rc");
        if (unitAt > 0) ann["_unit"] = tokens[unitAt];
        for (int i = unitAt + 1; i > 0 && i < tokens.Length; i++)
        {
            string t = tokens[i];
            if (t is "linear" or "circular") ann[RecordKeys.Topology] = t;
            else if (DateTime.TryParseExact(t, "dd-MMM-yyyy", CultureInfo.InvariantCulture, DateTimeStyles.None, out _)) ann[RecordKeys.Date] = t;
            else if (Divisions.Contains(t)) ann[RecordKeys.Division] = t;
            else if (!ann.ContainsKey(RecordKeys.MoleculeType)) ann[RecordKeys.MoleculeType] = t;
        }
    }
}

/// <summary>GenBank flat-file writer (79-column layout, INSDC feature table, ORIGIN in blocks of ten).</summary>
public static class GenBankWriter
{
    /// <summary>Writes one record. Missing header values get sensible defaults (DNA, linear, UNK division, today's date).</summary>
    public static void Write(TextWriter w, Sequence sequence)
    {
        var a = sequence.Annotations;
        bool protein = sequence.Type == SequenceType.Protein;
        string locus = a.TryGetValue(RecordKeys.Locus, out var l) ? l : sequence.Id.Split('.')[0];
        string mol = a.TryGetValue(RecordKeys.MoleculeType, out var m) ? m : protein ? "" : sequence.Type == SequenceType.Rna ? "RNA" : "DNA";
        string topology = a.TryGetValue(RecordKeys.Topology, out var t) ? t : protein ? "" : "linear";
        string division = a.TryGetValue(RecordKeys.Division, out var d) ? d : "UNK";
        string date = a.TryGetValue(RecordKeys.Date, out var dt) ? dt : DateTime.UtcNow.ToString("dd-MMM-yyyy", CultureInfo.InvariantCulture).ToUpperInvariant();

        // Columns per the GenBank release notes: name 13-28, length 30-40, unit 42-43, molecule 45-53, topology 56-63, division 65-67, date 69-79.
        var sb = new StringBuilder("LOCUS       ");
        sb.Append(locus.PadRight(16)).Append(' ');
        sb.Append(sequence.Length.ToString(CultureInfo.InvariantCulture).PadLeft(11)).Append(protein ? " aa " : " bp ");
        sb.Append("   ").Append(mol.PadRight(6)).Append("  ").Append(topology.PadRight(8)).Append(' ').Append(division).Append(' ').Append(date);
        w.Write(sb.ToString().TrimEnd()); w.Write('\n');

        string definition = a.TryGetValue(RecordKeys.Definition, out var def) ? def : sequence.Description;
        if (definition.Length == 0) definition = ".";
        else if (!definition.EndsWith('.')) definition += ".";
        FeatureTableWriter.WrapKeyword(w, "DEFINITION  ", "            ", definition);
        string accession = a.TryGetValue(RecordKeys.Accession, out var acc) ? acc : sequence.Id.Split('.')[0];
        FeatureTableWriter.WrapKeyword(w, "ACCESSION   ", "            ", accession);
        FeatureTableWriter.WrapKeyword(w, "VERSION     ", "            ", a.TryGetValue(RecordKeys.Version, out var ver) ? ver : sequence.Id);
        string keywords = a.TryGetValue(RecordKeys.Keywords, out var kw) && kw.Length > 0 ? kw + "." : ".";
        FeatureTableWriter.WrapKeyword(w, "KEYWORDS    ", "            ", keywords);
        string organism = a.TryGetValue(RecordKeys.Organism, out var org) ? org : ".";
        FeatureTableWriter.WrapKeyword(w, "SOURCE      ", "            ", a.TryGetValue(RecordKeys.Source, out var src) ? src : organism);
        FeatureTableWriter.WrapKeyword(w, "  ORGANISM  ", "            ", organism);
        if (a.TryGetValue(RecordKeys.Taxonomy, out var tax) && tax.Length > 0)
            FeatureTableWriter.WrapKeyword(w, "            ", "            ", tax + ".");
        for (int r = 1; a.ContainsKey(RecordKeys.Reference(r, "location")); r++)
        {
            FeatureTableWriter.WrapKeyword(w, "REFERENCE   ", "            ", a[RecordKeys.Reference(r, "location")]);
            foreach (var (field, label) in new[] { ("authors", "  AUTHORS   "), ("consrtm", "  CONSRTM   "), ("title", "  TITLE     "), ("journal", "  JOURNAL   "), ("pubmed", "   PUBMED   "), ("remark", "  REMARK    ") })
                if (a.TryGetValue(RecordKeys.Reference(r, field), out var v))
                    FeatureTableWriter.WrapKeyword(w, label, "            ", v);
        }
        if (a.TryGetValue(RecordKeys.Comment, out var comment))
        {
            var lines = comment.Split('\n');
            for (int i = 0; i < lines.Length; i++)
                FeatureTableWriter.WrapKeyword(w, i == 0 ? "COMMENT     " : "            ", "            ", lines[i]);
        }
        w.Write("FEATURES             Location/Qualifiers\n");
        FeatureTableWriter.Write(w, "     ", sequence.Features);
        w.Write("ORIGIN\n");
        WriteOrigin(w, sequence.Span);
        w.Write("//\n");
    }

    /// <summary>Writes many records.</summary>
    public static void Write(TextWriter w, IEnumerable<Sequence> sequences)
    {
        foreach (var s in sequences) Write(w, s);
    }

    /// <summary>Writes records to a file (gzip when it ends in <c>.gz</c>).</summary>
    public static void Write(string path, IEnumerable<Sequence> sequences)
    {
        using var w = BioFile.CreateText(path);
        Write(w, sequences);
    }

    /// <summary>GenBank text for one record.</summary>
    public static string ToString(Sequence sequence)
    {
        var sw = new StringWriter { NewLine = "\n" };
        Write(sw, sequence);
        return sw.ToString();
    }

    // "        1 gatcctccat atacaacggt ..." — position right-aligned in nine columns, six blocks of ten, lower case.
    private static void WriteOrigin(TextWriter w, ReadOnlySpan<byte> residues)
    {
        Span<char> line = stackalloc char[80];
        for (int i = 0; i < residues.Length; i += 60)
        {
            int n = 0;
            string pos = (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(9);
            pos.AsSpan().CopyTo(line); n = 9;
            int end = Math.Min(residues.Length, i + 60);
            for (int j = i; j < end; j++)
            {
                if ((j - i) % 10 == 0) line[n++] = ' ';
                line[n++] = (char)(residues[j] | 0x20);
            }
            w.Write(line[..n]);
            w.Write('\n');
        }
    }
}
