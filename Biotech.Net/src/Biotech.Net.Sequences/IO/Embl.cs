using System.Globalization;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Sequences.IO;

/// <summary>
/// Streaming EMBL flat-file reader (ID, AC, SV, DT, DE, KW, OS, OC, RN/RA/RT/RL, CC, FT, SQ). Both the current
/// <c>ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.</c> and the pre-2006 ID line are understood.
/// Annotations use the same <see cref="RecordKeys"/> as GenBank, so records convert between the two formats.
/// SQ residues are upper-cased, as Biopython's EMBL parser does.
/// </summary>
public static class EmblReader
{
    /// <summary>Reads every record of <paramref name="reader"/>.</summary>
    public static IEnumerable<Sequence> Read(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.StartsWith("ID   ", StringComparison.Ordinal)) yield return ReadRecord(line, reader);
    }

    /// <summary>Reads an EMBL file (optionally gzip-compressed) lazily.</summary>
    public static IEnumerable<Sequence> Read(string path)
    {
        using var reader = BioFile.OpenText(path);
        foreach (var s in Read(reader)) yield return s;
    }

    /// <summary>Parses EMBL text held in memory.</summary>
    public static IReadOnlyList<Sequence> Parse(string text) => Read(new StringReader(text)).ToList();

    private static Sequence ReadRecord(string idLine, TextReader reader)
    {
        var ann = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        string? svFromId = ParseId(idLine, ann);
        var residues = new ByteBuilder(4096);
        var features = new FeatureTableBuilder();
        var de = new StringBuilder();
        var kw = new StringBuilder();
        var oc = new StringBuilder();
        var cc = new StringBuilder();
        var os = new StringBuilder();
        var accessions = new List<string>();
        int rn = 0;
        string? lastRefCode = null;
        bool inSequence = false;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("//", StringComparison.Ordinal)) break;
            if (inSequence)
            {
                residues.AppendResidues(line, skipDigits: true, upperCase: true);
                continue;
            }
            if (line.Length < 2) continue;
            string code = line[..2];
            string content = line.Length > 5 ? line[5..] : string.Empty;
            if (code != lastRefCode && code is not ("RA" or "RT" or "RL" or "RX" or "RG" or "RC" or "RP")) lastRefCode = null;
            switch (code)
            {
                case "AC":
                    foreach (var a in content.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)) accessions.Add(a);
                    break;
                case "SV": svFromId ??= content.Trim(); break;
                case "DT": ann[RecordKeys.Date] = content.Split(' ')[0]; break;
                case "DE": Append(de, content); break;
                case "KW": Append(kw, content); break;
                case "OS": Append(os, content); break;
                case "OC": Append(oc, content); break;
                case "CC": if (cc.Length > 0) cc.Append('\n'); cc.Append(content.TrimEnd()); break;
                case "RN": rn++; break;
                case "RA" or "RT" or "RL" or "RX" or "RG" or "RC" or "RP":
                    {
                        string field = code switch { "RA" => "authors", "RT" => "title", "RL" => "journal", "RX" => "xref", "RG" => "consrtm", "RC" => "remark", _ => "location" };
                        string key = RecordKeys.Reference(Math.Max(rn, 1), field);
                        ann[key] = ann.TryGetValue(key, out var prev) ? prev + " " + content.Trim() : content.Trim();
                        lastRefCode = code;
                        break;
                    }
                case "FT": features.AddLine(content.Length > 0 || line.Length > 5 ? line.AsSpan(5) : ReadOnlySpan<char>.Empty); break;
                case "SQ": inSequence = true; break;
            }
        }
        if (accessions.Count > 0) ann[RecordKeys.Accession] = accessions[0];
        if (de.Length > 0) ann[RecordKeys.Definition] = de.ToString();
        if (kw.Length > 0) ann[RecordKeys.Keywords] = kw.ToString().TrimEnd('.').Trim();
        if (os.Length > 0) { ann[RecordKeys.Organism] = os.ToString(); ann[RecordKeys.Source] = os.ToString(); }
        if (oc.Length > 0) ann[RecordKeys.Taxonomy] = oc.ToString().TrimEnd('.');
        if (cc.Length > 0) ann[RecordKeys.Comment] = cc.ToString();
        for (int r = 1; r <= rn; r++)
            foreach (var f in new[] { "title", "authors" })
                if (ann.TryGetValue(RecordKeys.Reference(r, f), out var v)) ann[RecordKeys.Reference(r, f)] = v.Trim('"', ';', ' ');

        string acc = ann.GetValueOrDefault(RecordKeys.Accession) ?? ann.GetValueOrDefault(RecordKeys.Locus) ?? "unknown";
        if (svFromId is not null)
        {
            string version = svFromId.Contains('.') ? svFromId : $"{acc}.{svFromId}";
            ann[RecordKeys.Version] = version;
        }
        if (!ann.ContainsKey(RecordKeys.Locus)) ann[RecordKeys.Locus] = acc;

        var bytes = residues.ToArray();
        bool protein = ann.TryGetValue("_unit", out var unit) && unit == "AA";
        ann.Remove("_unit");
        var alphabet = protein ? Alphabet.ProteinExtended : Alphabet.DnaIupac.IsValid(bytes) ? Alphabet.DnaIupac : Alphabet.Detect(bytes);
        string id = ann.GetValueOrDefault(RecordKeys.Version) ?? acc;
        var seq = Sequence.FromBuffer(id, bytes, alphabet, ann.GetValueOrDefault(RecordKeys.Definition) ?? string.Empty, validate: false);
        foreach (var kv in ann) seq.Annotations[kv.Key] = kv.Value;
        foreach (var f in features.Finish()) seq.Features.Add(f);
        return seq;
    }

    private static void Append(StringBuilder sb, string content)
    {
        if (sb.Length > 0) sb.Append(' ');
        sb.Append(content.Trim());
    }

    // Returns the sequence version when the ID line carries one ("SV 1").
    private static string? ParseId(string line, Dictionary<string, string> ann)
    {
        var body = line[5..].Trim().TrimEnd('.');
        var fields = body.Split(';', StringSplitOptions.TrimEntries);
        string? sv = null;
        if (fields.Length >= 7 && fields[1].StartsWith("SV", StringComparison.Ordinal))
        {
            // Current format: accession; SV version; topology; molecule; data class; division; length BP.
            ann[RecordKeys.Accession] = fields[0];
            ann[RecordKeys.Locus] = fields[0];
            sv = fields[1][2..].Trim();
            ann[RecordKeys.Topology] = fields[2];
            ann[RecordKeys.MoleculeType] = fields[3];
            ann["data_class"] = fields[4];
            ann[RecordKeys.Division] = fields[5];
            ann["_unit"] = fields[6].Split(' ', StringSplitOptions.RemoveEmptyEntries).LastOrDefault() ?? "BP";
        }
        else
        {
            // Pre-2006: "entryname  dataclass; [circular] molecule; division; length BP."
            var first = fields[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            ann[RecordKeys.Locus] = first[0];
            if (first.Length > 1) ann["data_class"] = first[1];
            if (fields.Length > 1)
            {
                var mol = fields[1].Split(' ', StringSplitOptions.RemoveEmptyEntries);
                if (mol.Length > 1 && mol[0] == "circular") { ann[RecordKeys.Topology] = "circular"; ann[RecordKeys.MoleculeType] = mol[1]; }
                else { ann[RecordKeys.Topology] = "linear"; if (mol.Length > 0) ann[RecordKeys.MoleculeType] = mol[0]; }
            }
            if (fields.Length > 2) ann[RecordKeys.Division] = fields[2];
            if (fields.Length > 3) ann["_unit"] = fields[^1].Split(' ', StringSplitOptions.RemoveEmptyEntries).LastOrDefault() ?? "BP";
        }
        return sv;
    }
}

/// <summary>EMBL flat-file writer (current ID-line format, 80-column sequence block with running counts).</summary>
public static class EmblWriter
{
    /// <summary>Writes one record.</summary>
    public static void Write(TextWriter w, Sequence sequence)
    {
        var a = sequence.Annotations;
        bool protein = sequence.Type == SequenceType.Protein;
        string accession = a.TryGetValue(RecordKeys.Accession, out var acc) ? acc : sequence.Id.Split('.')[0];
        string version = a.TryGetValue(RecordKeys.Version, out var ver) ? ver : sequence.Id;
        string sv = version.Contains('.') ? version[(version.LastIndexOf('.') + 1)..] : "1";
        string topology = a.TryGetValue(RecordKeys.Topology, out var t) ? t : "linear";
        string mol = a.TryGetValue(RecordKeys.MoleculeType, out var m) ? m : protein ? "protein" : sequence.Type == SequenceType.Rna ? "RNA" : "genomic DNA";
        string dataClass = a.TryGetValue("data_class", out var dc) ? dc : "STD";
        string division = a.TryGetValue(RecordKeys.Division, out var d) ? d : "UNC";
        w.Write($"ID   {accession}; SV {sv}; {topology}; {mol}; {dataClass}; {division}; {sequence.Length} {(protein ? "AA" : "BP")}.\n");
        w.Write("XX\n");
        w.Write($"AC   {accession};\n");
        w.Write("XX\n");
        if (a.TryGetValue(RecordKeys.Date, out var date)) { w.Write($"DT   {date}\n"); w.Write("XX\n"); }
        string de = a.TryGetValue(RecordKeys.Definition, out var def) ? def : sequence.Description;
        FeatureTableWriter.WrapKeyword(w, "DE   ", "DE   ", de.Length == 0 ? "." : de, 80);
        w.Write("XX\n");
        string keywords = a.TryGetValue(RecordKeys.Keywords, out var kw) && kw.Length > 0 ? kw + "." : ".";
        FeatureTableWriter.WrapKeyword(w, "KW   ", "KW   ", keywords, 80);
        w.Write("XX\n");
        if (a.TryGetValue(RecordKeys.Organism, out var org))
        {
            FeatureTableWriter.WrapKeyword(w, "OS   ", "OS   ", org, 80);
            if (a.TryGetValue(RecordKeys.Taxonomy, out var tax)) FeatureTableWriter.WrapKeyword(w, "OC   ", "OC   ", tax + ".", 80);
            w.Write("XX\n");
        }
        for (int r = 1; a.ContainsKey(RecordKeys.Reference(r, "authors")) || a.ContainsKey(RecordKeys.Reference(r, "title")); r++)
        {
            w.Write($"RN   [{r}]\n");
            if (a.TryGetValue(RecordKeys.Reference(r, "authors"), out var au)) FeatureTableWriter.WrapKeyword(w, "RA   ", "RA   ", au.TrimEnd(';') + ";", 80);
            string title = a.TryGetValue(RecordKeys.Reference(r, "title"), out var ti) ? $"\"{ti}\";" : ";";
            FeatureTableWriter.WrapKeyword(w, "RT   ", "RT   ", title, 80);
            if (a.TryGetValue(RecordKeys.Reference(r, "journal"), out var jl)) FeatureTableWriter.WrapKeyword(w, "RL   ", "RL   ", jl, 80);
            w.Write("XX\n");
        }
        if (a.TryGetValue(RecordKeys.Comment, out var cc))
        {
            foreach (var l in cc.Split('\n')) FeatureTableWriter.WrapKeyword(w, "CC   ", "CC   ", l, 80);
            w.Write("XX\n");
        }
        if (sequence.HasFeatures)
        {
            w.Write("FH   Key             Location/Qualifiers\nFH\n");
            FeatureTableWriter.Write(w, "FT   ", sequence.Features);
            w.Write("XX\n");
        }
        var counts = sequence.NucleotideCounts();
        if (protein) w.Write($"SQ   Sequence {sequence.Length} AA;\n");
        else w.Write($"SQ   Sequence {sequence.Length} BP; {counts.A} A; {counts.C} C; {counts.G} G; {counts.T} T; {counts.N + counts.Other} other;\n");
        WriteSequence(w, sequence.Span);
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

    /// <summary>EMBL text for one record.</summary>
    public static string ToString(Sequence sequence)
    {
        var sw = new StringWriter { NewLine = "\n" };
        Write(sw, sequence);
        return sw.ToString();
    }

    // "     aaacaaacca aatatggatt ...        60": five spaces, six blocks of ten, running count ending in column 80.
    private static void WriteSequence(TextWriter w, ReadOnlySpan<byte> residues)
    {
        const int Width = 80;
        var sb = new StringBuilder(Width + 1);
        for (int i = 0; i < residues.Length; i += 60)
        {
            sb.Clear().Append("    ");
            int end = Math.Min(residues.Length, i + 60);
            for (int j = i; j < end; j++)
            {
                if ((j - i) % 10 == 0) sb.Append(' ');
                sb.Append((char)(residues[j] | 0x20));
            }
            // A full line holds 4 + 6 separators + 60 = 70 characters, so the count is padded into columns 71-80.
            sb.Append(' ');
            sb.Append(end.ToString(CultureInfo.InvariantCulture).PadLeft(Math.Max(0, Width - sb.Length)));
            w.Write(sb.ToString()); w.Write('\n');
        }
    }
}
