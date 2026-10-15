using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Genomics.GeneFinding;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Genomics.Annotation;

/// <summary>Writes features as GFF3 (1-based closed coordinates, column 9 percent-encoded).</summary>
public static class GffWriter
{
    /// <summary>Writes the <c>##gff-version 3</c> pragma.</summary>
    public static void WriteHeader(TextWriter w) => w.Write("##gff-version 3\n");

    /// <summary>
    /// Writes one line per segment of each feature on <paramref name="sequenceId"/>. Qualifiers become attributes;
    /// a <c>score</c> qualifier fills column 6.
    /// </summary>
    public static void Write(TextWriter w, string sequenceId, IEnumerable<SequenceFeature> features, string source = "Biotech.Net")
    {
        foreach (var f in features)
        {
            char strand = f.Location.Strand switch { Strand.Forward => '+', Strand.Reverse => '-', _ => '.' };
            string score = f["score"] ?? ".";
            var attrs = new StringBuilder();
            foreach (var kv in f.Qualifiers)
            {
                if (kv.Key.Equals("score", StringComparison.OrdinalIgnoreCase)) continue;
                if (attrs.Length > 0) attrs.Append(';');
                attrs.Append(Escape(kv.Key)).Append('=').Append(string.Join(',', kv.Value.Select(Escape)));
            }
            int phase = 0;
            var segs = f.Location.Segments;
            foreach (var (s, e) in segs)
            {
                string ph = f.Type == "CDS" ? phase.ToString(CultureInfo.InvariantCulture) : ".";
                w.Write($"{sequenceId}\t{source}\t{f.Type}\t{s + 1}\t{e}\t{score}\t{strand}\t{ph}\t{(attrs.Length == 0 ? "." : attrs.ToString())}\n");
                phase = (3 - ((e - s - phase) % 3 + 3) % 3) % 3;
            }
        }
    }

    /// <summary>Writes predicted genes (gene + CDS lines with ID/Parent) as a complete GFF3 document.</summary>
    public static void Write(TextWriter w, IEnumerable<PredictedGene> genes, string source = "Biotech.Net")
    {
        WriteHeader(w);
        int n = 0;
        foreach (var g in genes)
        {
            n++;
            string id = $"gene_{n}";
            char strand = g.Strand == Strand.Reverse ? '-' : '+';
            string score = g.Score.ToString("F2", CultureInfo.InvariantCulture);
            w.Write($"{g.SequenceId}\t{source}\tgene\t{g.Start + 1}\t{g.End}\t{score}\t{strand}\t.\tID={id}\n");
            string rbs = g.RbsMotif is null ? "" : $";rbs_motif={g.RbsMotif};rbs_spacer={g.RbsSpacer}";
            w.Write($"{g.SequenceId}\t{source}\tCDS\t{g.Start + 1}\t{g.End}\t{score}\t{strand}\t0\tID=cds_{n};Parent={id};start_codon={g.StartCodon}{rbs}\n");
        }
    }

    /// <summary>Predicted genes as GFF3 text.</summary>
    public static string ToString(IEnumerable<PredictedGene> genes)
    {
        var sw = new StringWriter { NewLine = "\n" };
        Write(sw, genes);
        return sw.ToString();
    }

    private static string Escape(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
        {
            if (c is ';' or '=' or '&' or ',' or '%' or '\t' or '\n' or '\r' || c < ' ') sb.Append('%').Append(((int)c).ToString("X2"));
            else sb.Append(c);
        }
        return sb.ToString();
    }
}

/// <summary>One GFF3 line.</summary>
public sealed record GffRecord(string SequenceId, string Source, string Type, long Start, long End, double? Score,
    Strand Strand, int? Phase, IReadOnlyDictionary<string, string> Attributes)
{
    /// <summary>0-based half-open interval of the line.</summary>
    public GenomicInterval Interval => new(SequenceId, Start, End, Strand);

    /// <summary>The <c>ID</c> attribute.</summary>
    public string? Id => Attributes.GetValueOrDefault("ID");

    /// <summary>The <c>Parent</c> attribute (first parent).</summary>
    public string? Parent => Attributes.TryGetValue("Parent", out var p) ? p.Split(',')[0] : null;
}

/// <summary>GFF3 / GTF reader, and conversion of gene → mRNA → exon/CDS hierarchies into <see cref="Gene"/> models.</summary>
public static class GffReader
{
    /// <summary>Streams the records of GFF3 (or GTF, detected per line) text; coordinates become 0-based half-open.</summary>
    public static IEnumerable<GffRecord> Read(TextReader reader)
    {
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("##FASTA", StringComparison.Ordinal)) yield break;
            if (line.Length == 0 || line[0] == '#') continue;
            var c = line.Split('\t');
            if (c.Length < 9) continue;
            var attrs = ParseAttributes(c[8]);
            yield return new GffRecord(c[0], c[1], c[2],
                long.Parse(c[3], CultureInfo.InvariantCulture) - 1, long.Parse(c[4], CultureInfo.InvariantCulture),
                c[5] == "." ? null : double.Parse(c[5], CultureInfo.InvariantCulture),
                c[6] switch { "+" => Strand.Forward, "-" => Strand.Reverse, _ => Strand.Unknown },
                c[7] == "." ? null : int.Parse(c[7], CultureInfo.InvariantCulture), attrs);
        }
    }

    /// <summary>Reads a GFF3/GTF file (gzip-aware).</summary>
    public static IEnumerable<GffRecord> Read(string path)
    {
        using var r = BioFile.OpenText(path);
        foreach (var rec in Read(r)) yield return rec;
    }

    private static Dictionary<string, string> ParseAttributes(string col)
    {
        var d = new Dictionary<string, string>(StringComparer.Ordinal);
        if (col == ".") return d;
        foreach (var raw in col.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            int eq = raw.IndexOf('=');
            if (eq > 0) { d[raw[..eq]] = Uri.UnescapeDataString(raw[(eq + 1)..]); continue; }
            int sp = raw.IndexOf(' '); // GTF: key "value"
            if (sp > 0) d[raw[..sp]] = raw[(sp + 1)..].Trim().Trim('"');
        }
        return d;
    }

    /// <summary>
    /// Builds gene models: <c>gene</c> lines become <see cref="Gene"/>s, <c>mRNA</c>/<c>transcript</c> lines their
    /// <see cref="Transcript"/>s with exons from <c>exon</c> lines (or CDS/UTR lines when exons are absent) and the coding
    /// region spanning the transcript's <c>CDS</c> lines. GTF input (gene_id/transcript_id) is supported too.
    /// </summary>
    public static IReadOnlyList<Gene> BuildGenes(IEnumerable<GffRecord> records)
    {
        var recs = records.ToList();
        var genes = new Dictionary<string, Gene>(StringComparer.Ordinal);
        var order = new List<Gene>();
        var txGene = new Dictionary<string, string>(StringComparer.Ordinal);
        var exons = new Dictionary<string, List<GenomicInterval>>(StringComparer.Ordinal);
        var cds = new Dictionary<string, List<GenomicInterval>>(StringComparer.Ordinal);
        var other = new Dictionary<string, List<GenomicInterval>>(StringComparer.Ordinal);

        foreach (var r in recs.Where(r => r.Type == "gene"))
        {
            string id = r.Id ?? r.Attributes.GetValueOrDefault("gene_id") ?? $"gene_{order.Count + 1}";
            var g = new Gene(id, r.Interval, r.Attributes.GetValueOrDefault("Name") ?? r.Attributes.GetValueOrDefault("gene_name"),
                r.Attributes.GetValueOrDefault("biotype") ?? r.Attributes.GetValueOrDefault("gene_biotype") ?? "protein_coding");
            foreach (var kv in r.Attributes) g.Attributes[kv.Key] = kv.Value;
            if (genes.TryAdd(id, g)) order.Add(g);
        }
        foreach (var r in recs.Where(r => r.Type is "mRNA" or "transcript"))
        {
            string tid = r.Id ?? r.Attributes.GetValueOrDefault("transcript_id") ?? "";
            string gid = r.Attributes.GetValueOrDefault("Parent") ?? r.Attributes.GetValueOrDefault("gene_id") ?? tid;
            if (tid.Length > 0) txGene[tid] = gid;
        }
        foreach (var r in recs)
        {
            string? parent = r.Type is "exon" or "CDS" or "five_prime_UTR" or "three_prime_UTR" or "UTR" or "start_codon" or "stop_codon"
                ? r.Parent ?? r.Attributes.GetValueOrDefault("transcript_id") : null;
            if (parent is null) continue;
            var target = r.Type == "exon" ? exons : r.Type == "CDS" ? cds : r.Type.EndsWith("UTR") ? other : null;
            if (target is null) continue;
            if (!target.TryGetValue(parent, out var list)) target[parent] = list = [];
            list.Add(r.Interval);
            if (!txGene.ContainsKey(parent))
                txGene[parent] = r.Attributes.GetValueOrDefault("gene_id") ?? parent;
        }
        foreach (var (tid, gid) in txGene)
        {
            var ex = exons.GetValueOrDefault(tid);
            var cd = cds.GetValueOrDefault(tid);
            if (ex is null)
            {
                var parts = (cd ?? []).Concat(other.GetValueOrDefault(tid) ?? []).ToList();
                if (parts.Count == 0) continue;
                ex = GenomicInterval.Merge(parts).Select(i => i with { Strand = parts[0].Strand }).ToList();
            }
            if (!genes.TryGetValue(gid, out var gene))
            {
                var span = ex[0] with { Start = ex.Min(e => e.Start), End = ex.Max(e => e.End) };
                gene = new Gene(gid, span);
                genes[gid] = gene;
                order.Add(gene);
            }
            GenomicInterval? coding = cd is { Count: > 0 } ? cd[0] with { Start = cd.Min(c => c.Start), End = cd.Max(c => c.End) } : null;
            gene.Transcripts.Add(new Transcript(tid, gid, ex) { CodingRegion = coding });
        }
        return order;
    }

    /// <summary>Reads a GFF3/GTF file and adds its genes to <paramref name="genome"/>.</summary>
    public static Genome LoadGenes(Genome genome, string path)
    {
        foreach (var g in BuildGenes(Read(path))) genome.Add(g);
        return genome;
    }
}
