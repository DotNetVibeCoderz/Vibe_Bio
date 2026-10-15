using System.Globalization;
using System.Text;
using Biotech.Net.Chemistry;
using Biotech.Net.Core;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Transcriptomics;
using Biotech.Net.Visualization;

namespace Biotech.Net.Notebooks;

/// <summary>How a plot is embedded in notebook output.</summary>
public enum PlotEmbedding
{
    /// <summary>Interactive Plotly HTML when the plot supports it, otherwise inline SVG.</summary>
    Auto,

    /// <summary>Inline SVG.</summary>
    Svg,

    /// <summary>A base64 PNG &lt;img&gt;.</summary>
    Png,
}

/// <summary>
/// Rich HTML renderings of Biotech.Net objects for Polyglot / .NET Interactive / Jupyter notebooks. Every method
/// returns a self-contained fragment (scoped CSS included, light and dark aware) and can be used without the
/// formatter registry, e.g. in a Blazor page.
/// </summary>
public static class BiotechHtml
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Plot embedding used by <see cref="Plot(BioPlot, PlotEmbedding?)"/> when none is given.</summary>
    public static PlotEmbedding DefaultPlotEmbedding { get; set; } = PlotEmbedding.Auto;

    /// <summary>Residues per line in sequence views.</summary>
    public static int LineWidth { get; set; } = 60;

    /// <summary>Residues shown before a sequence is truncated.</summary>
    public static int MaxResidues { get; set; } = 3000;

    /// <summary>Scoped stylesheet shared by all fragments.</summary>
    public const string Styles = """
<style>
.btn-card{font-family:'Segoe UI',system-ui,-apple-system,sans-serif;font-size:13px;line-height:1.45;color:var(--vscode-editor-foreground,#1C2430);
 border:1px solid rgba(127,137,150,.35);border-radius:10px;padding:10px 14px;margin:4px 0;max-width:980px;background:transparent}
.btn-card .btn-h{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap;margin-bottom:6px}
.btn-card .btn-t{font-weight:600;font-size:14px}
.btn-card .btn-k{font-size:11px;text-transform:uppercase;letter-spacing:.06em;opacity:.65}
.btn-card .btn-chips{display:flex;gap:6px;flex-wrap:wrap;margin:4px 0 8px}
.btn-card .btn-chip{font-size:11.5px;padding:1px 8px;border-radius:999px;background:rgba(47,111,219,.12);white-space:nowrap}
.btn-card .btn-seq{font-family:'Cascadia Mono',Consolas,'SFMono-Regular',monospace;font-size:12.5px;white-space:pre;overflow-x:auto;margin:0}
.btn-card .btn-seq .n{opacity:.5;user-select:none}
.btn-card .btn-gc{height:6px;border-radius:3px;background:rgba(127,137,150,.25);overflow:hidden;width:180px;display:inline-block;vertical-align:middle}
.btn-card .btn-gc>span{display:block;height:100%;background:#2F6FDB}
.btn-card table.btn-tbl{border-collapse:collapse;font-size:12px;margin:2px 0}
.btn-card table.btn-tbl th,.btn-card table.btn-tbl td{padding:3px 9px;border-bottom:1px solid rgba(127,137,150,.25);text-align:right}
.btn-card table.btn-tbl th{font-weight:600;opacity:.8}
.btn-card table.btn-tbl td.l,.btn-card table.btn-tbl th.l{text-align:left}
.btn-card .btn-foot{font-size:10.5px;opacity:.55;margin-top:8px}
.btn-card .up{color:#D1414A}.btn-card .down{color:#2F6FDB}.btn-card .muted{opacity:.6}
.btn-A{color:#1E9E5A}.btn-C{color:#2F6FDB}.btn-G{color:#C98C00}.btn-T{color:#D1414A}.btn-X{color:#8A93A0}
.btn-ph{color:#2F6FDB}.btn-pos{color:#D1414A}.btn-neg{color:#B23CC8}.btn-pol{color:#1E9E5A}.btn-gly{color:#E07A1F}.btn-pro{color:#B59A00}.btn-aro{color:#1597A5}
</style>
""";

    /// <summary>The credit line shown in fragment footers.</summary>
    public static string Footer() => $"<div class=\"btn-foot\">{Enc(BiotechInfo.Product)} · {Enc(BiotechInfo.CreditsId)}</div>";

    /// <summary>Escapes the HTML-significant characters (&amp;, &lt;, &gt;, quotes) and leaves other Unicode as is.</summary>
    internal static string Enc(string? s)
    {
        if (string.IsNullOrEmpty(s)) return string.Empty;
        if (s.AsSpan().IndexOfAny("&<>\"'") < 0) return s;
        var sb = new StringBuilder(s.Length + 16);
        foreach (char c in s)
            sb.Append(c switch { '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&#39;", _ => c.ToString() });
        return sb.ToString();
    }

    internal static string Card(string kind, string title, string chips, string body, bool footer = true)
        => $"{Styles}<div class=\"btn-card\"><div class=\"btn-h\"><span class=\"btn-k\">{Enc(kind)}</span><span class=\"btn-t\">{Enc(title)}</span></div>"
           + (chips.Length > 0 ? $"<div class=\"btn-chips\">{chips}</div>" : "") + body + (footer ? Footer() : "") + "</div>";

    internal static string Chip(string text) => $"<span class=\"btn-chip\">{Enc(text)}</span>";

    internal static string ResidueClass(char c, bool nucleotide)
    {
        c = char.ToUpperInvariant(c);
        if (nucleotide) return c switch { 'A' => "btn-A", 'C' => "btn-C", 'G' => "btn-G", 'T' or 'U' => "btn-T", _ => "btn-X" };
        return c switch
        {
            'A' or 'V' or 'L' or 'I' or 'M' or 'F' or 'W' or 'C' => "btn-ph",
            'K' or 'R' => "btn-pos",
            'D' or 'E' => "btn-neg",
            'S' or 'T' or 'N' or 'Q' => "btn-pol",
            'G' => "btn-gly",
            'P' => "btn-pro",
            'H' or 'Y' => "btn-aro",
            _ => "btn-X",
        };
    }

    /// <summary>Coloured, wrapped residues with 1-based line numbers (runs of one colour share a span).</summary>
    public static string Residues(ReadOnlySpan<byte> residues, bool nucleotide, int lineWidth = 0, int maxResidues = 0)
    {
        lineWidth = lineWidth > 0 ? lineWidth : LineWidth;
        maxResidues = maxResidues > 0 ? maxResidues : MaxResidues;
        int shown = Math.Min(residues.Length, maxResidues);
        int numWidth = Math.Max(3, shown.ToString(Inv).Length);
        var sb = new StringBuilder("<pre class=\"btn-seq\">");
        for (int start = 0; start < shown; start += lineWidth)
        {
            sb.Append("<span class=\"n\">").Append((start + 1).ToString(Inv).PadLeft(numWidth)).Append("  </span>");
            int end = Math.Min(shown, start + lineWidth);
            string? cls = null;
            for (int i = start; i < end; i++)
            {
                char c = (char)residues[i];
                string k = ResidueClass(c, nucleotide);
                if (k != cls)
                {
                    if (cls is not null) sb.Append("</span>");
                    sb.Append("<span class=\"").Append(k).Append("\">");
                    cls = k;
                }
                if ((i - start) % 10 == 0 && i > start) sb.Append(' ');
                sb.Append(c is '<' or '>' or '&' ? '?' : c);
            }
            if (cls is not null) sb.Append("</span>");
            sb.Append('\n');
        }
        if (shown < residues.Length)
            sb.Append("<span class=\"n\">… ").Append((residues.Length - shown).ToString("N0", Inv)).Append(" more residues</span>");
        sb.Append("</pre>");
        return sb.ToString();
    }

    /// <summary>A sequence card: id, type, length, GC (nucleotides), coloured wrapped residues.</summary>
    public static string Sequence(Sequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        bool nuc = sequence.Alphabet.IsNucleotide;
        var chips = new StringBuilder();
        chips.Append(Chip(sequence.Type.ToString())).Append(Chip($"{sequence.Length:N0} {(nuc ? "bp" : "aa")}"));
        string extra = "";
        if (nuc && sequence.Length > 0)
        {
            double gc = sequence.GcContent();
            chips.Append(Chip($"GC {gc * 100:F1}%"));
            extra = $"<div class=\"muted\" style=\"font-size:11px;margin-bottom:6px\">GC <span class=\"btn-gc\"><span style=\"width:{gc * 100:F1}%\"></span></span></div>";
        }
        if (sequence.HasFeatures) chips.Append(Chip($"{sequence.Features.Count} features"));
        string desc = string.IsNullOrEmpty(sequence.Description) ? "" : $"<div class=\"muted\" style=\"margin-bottom:4px\">{Enc(sequence.Description)}</div>";
        return Card("sequence", sequence.Id, chips.ToString(), desc + extra + Residues(sequence.Span, nuc));
    }

    /// <summary>A read card: bases coloured, with a per-base quality bar strip and summary statistics.</summary>
    public static string Read(SequenceRead read)
    {
        ArgumentNullException.ThrowIfNull(read);
        var chips = Chip($"{read.Length:N0} bp") + Chip($"mean Q{read.MeanQuality():F1}") + Chip($"expected errors {read.ExpectedErrors():F2}");
        var q = read.Quality;
        var bars = new StringBuilder("<svg role=\"img\" aria-label=\"per-base quality\" width=\"")
            .Append(Math.Min(900, q.Length * 4).ToString(Inv)).Append("\" height=\"34\" viewBox=\"0 0 ")
            .Append((q.Length * 4).ToString(Inv)).Append(" 34\" preserveAspectRatio=\"none\" style=\"display:block;margin:4px 0\">");
        for (int i = 0; i < q.Length && i < 2000; i++)
        {
            int v = q[i];
            string col = v >= 30 ? "#1E9E5A" : v >= 20 ? "#EBA300" : "#D1414A";
            double h = Math.Clamp(v, 0, 42) / 42.0 * 32;
            bars.Append(CultureInfo.InvariantCulture, $"<rect x=\"{i * 4}\" y=\"{34 - h:F1}\" width=\"3\" height=\"{h:F1}\" fill=\"{col}\"/>");
        }
        bars.Append("</svg>");
        return Card("read", read.Id, chips, bars + Residues(read.Sequence.Span, true));
    }

    /// <summary>A protein card: name, organism, mass, pI and coloured residues.</summary>
    public static string Protein(Protein protein)
    {
        ArgumentNullException.ThrowIfNull(protein);
        var chips = new StringBuilder();
        chips.Append(Chip($"{protein.Length:N0} aa")).Append(Chip($"{protein.MolecularWeight / 1000:F2} kDa")).Append(Chip($"pI {protein.IsoelectricPoint:F2}"));
        if (protein.GeneName is not null) chips.Append(Chip("gene " + protein.GeneName));
        if (protein.Organism is not null) chips.Append(Chip(protein.Organism));
        if (protein.StructureId is not null) chips.Append(Chip("PDB " + protein.StructureId));
        string fn = protein.Function is null ? "" : $"<div class=\"muted\" style=\"margin-bottom:4px\">{Enc(protein.Function)}</div>";
        return Card("protein", $"{protein.Name} ({protein.Id})", chips.ToString(), fn + Residues(protein.Sequence.Span, false));
    }

    /// <summary>A variant card (VCF fields, 1-based position, genotypes).</summary>
    public static string Variant(Variant v)
    {
        ArgumentNullException.ThrowIfNull(v);
        var chips = new StringBuilder();
        chips.Append(Chip(v.Type.ToString()));
        if (v.Type == VariantType.Snv) chips.Append(Chip(v.IsTransition ? "transition" : "transversion"));
        chips.Append(Chip(v.IsPass ? "PASS" : "FILTER " + v.Filter));
        if (!double.IsNaN(v.Quality)) chips.Append(Chip($"QUAL {v.Quality:G4}"));
        var rows = new StringBuilder("<table class=\"btn-tbl\"><tr><th class=\"l\">CHROM</th><th>POS</th><th class=\"l\">ID</th><th class=\"l\">REF</th><th class=\"l\">ALT</th></tr>");
        rows.Append(CultureInfo.InvariantCulture, $"<tr><td class=\"l\">{Enc(v.Chromosome)}</td><td>{v.Position:N0}</td><td class=\"l\">{Enc(v.Id ?? ".")}</td><td class=\"l\">{Enc(v.Reference)}</td><td class=\"l\">{Enc(string.Join(",", v.Alternates))}</td></tr></table>");
        if (v.Info.Count > 0)
            rows.Append("<div class=\"muted\" style=\"font-size:11.5px;margin-top:4px\">").Append(Enc(string.Join("; ", v.Info.Take(12).Select(kv => $"{kv.Key}={kv.Value}")))).Append("</div>");
        if (v.Genotypes.Count > 0)
        {
            double af = v.AlternateAlleleFrequency();
            rows.Append("<div style=\"margin-top:4px\">").Append(Enc($"{v.Genotypes.Count} genotypes · AF {af:F3} · ")).Append(Enc(string.Join(" ", v.Genotypes.Take(24)))).Append(v.Genotypes.Count > 24 ? " …" : "").Append("</div>");
        }
        return Card("variant", v.ToString(), chips.ToString(), rows.ToString());
    }

    /// <summary>A genomic interval card with 1-based and 0-based coordinates.</summary>
    public static string Interval(GenomicInterval i)
    {
        var chips = Chip($"{i.Length:N0} bp") + Chip(i.Strand switch { Strand.Forward => "+ strand", Strand.Reverse => "− strand", _ => "unstranded" });
        string body = $"<div class=\"muted\" style=\"font-size:12px\">0-based half-open [{i.Start.ToString("N0", Inv)}, {i.End.ToString("N0", Inv)}) · BED: {Enc(i.Chromosome)}\t{i.Start.ToString(Inv)}\t{i.End.ToString(Inv)}</div>";
        return Card("interval", i.ToString(), chips, body, footer: false);
    }

    /// <summary>A gene card with its transcripts.</summary>
    public static string Gene(Gene gene)
    {
        ArgumentNullException.ThrowIfNull(gene);
        var chips = Chip(gene.Biotype) + Chip(gene.Location.ToString()) + Chip($"{gene.Transcripts.Count} transcripts");
        var t = new StringBuilder();
        if (gene.Transcripts.Count > 0)
        {
            t.Append("<table class=\"btn-tbl\"><tr><th class=\"l\">transcript</th><th>exons</th><th>spliced length</th><th class=\"l\">CDS</th></tr>");
            foreach (var tr in gene.Transcripts.Take(30))
                t.Append(CultureInfo.InvariantCulture, $"<tr><td class=\"l\">{Enc(tr.Id)}</td><td>{tr.Exons.Count}</td><td>{tr.SplicedLength:N0}</td><td class=\"l\">{Enc(tr.CodingRegion?.ToString() ?? "—")}</td></tr>");
            t.Append("</table>");
        }
        return Card("gene", $"{gene.Name} ({gene.Id})", chips, t.ToString());
    }

    /// <summary>A dataset table (first rows / columns) with its shape and annotations.</summary>
    public static string Dataset(Dataset ds, int maxRows = 20, int maxColumns = 12)
    {
        ArgumentNullException.ThrowIfNull(ds);
        var chips = Chip($"{ds.Rows:N0} × {ds.Columns:N0}") + string.Concat(ds.ColumnAnnotations.Keys.Take(6).Select(k => Chip("annotation: " + k)));
        int rows = Math.Min(ds.Rows, maxRows), cols = Math.Min(ds.Columns, maxColumns);
        var sb = new StringBuilder("<div style=\"overflow-x:auto\"><table class=\"btn-tbl\"><tr><th class=\"l\"></th>");
        for (int j = 0; j < cols; j++) sb.Append("<th>").Append(Enc(ds.ColumnNames[j])).Append("</th>");
        if (cols < ds.Columns) sb.Append("<th>…</th>");
        sb.Append("</tr>");
        foreach (var (key, values) in ds.ColumnAnnotations.Take(4))
        {
            sb.Append("<tr class=\"muted\"><td class=\"l\"><i>").Append(Enc(key)).Append("</i></td>");
            for (int j = 0; j < cols; j++) sb.Append("<td>").Append(Enc(j < values.Length ? values[j] : "")).Append("</td>");
            if (cols < ds.Columns) sb.Append("<td></td>");
            sb.Append("</tr>");
        }
        for (int i = 0; i < rows; i++)
        {
            sb.Append("<tr><td class=\"l\">").Append(Enc(ds.RowNames[i])).Append("</td>");
            for (int j = 0; j < cols; j++) sb.Append("<td>").Append(FormatNumber(ds[i, j])).Append("</td>");
            if (cols < ds.Columns) sb.Append("<td>…</td>");
            sb.Append("</tr>");
        }
        if (rows < ds.Rows) sb.Append(CultureInfo.InvariantCulture, $"<tr><td class=\"l muted\" colspan=\"{cols + 2}\">… {ds.Rows - rows:N0} more rows</td></tr>");
        sb.Append("</table></div>");
        return Card("dataset", ds.Name, chips, sb.ToString());
    }

    /// <summary>The top of a differential-expression table with up/down counts.</summary>
    public static string DifferentialExpression(DifferentialExpressionResult result, int top = 15)
    {
        ArgumentNullException.ThrowIfNull(result);
        var (up, down) = result.CountSignificant(0.05);
        var chips = Chip($"{result.Count:N0} genes") + $"<span class=\"btn-chip up\">{up:N0} up</span><span class=\"btn-chip down\">{down:N0} down</span>" + Chip("padj < 0.05");
        var sb = new StringBuilder("<table class=\"btn-tbl\"><tr><th class=\"l\">gene</th><th>baseMean</th><th>log2FC</th><th>lfcSE</th><th>stat</th><th>pvalue</th><th>padj</th></tr>");
        foreach (var r in result.Top(top))
        {
            string cls = r.PAdj < 0.05 ? (r.Log2FoldChange > 0 ? "up" : "down") : "";
            sb.Append(CultureInfo.InvariantCulture,
                $"<tr><td class=\"l\">{Enc(r.Gene)}</td><td>{FormatNumber(r.BaseMean)}</td><td class=\"{cls}\">{r.Log2FoldChange:F3}</td><td>{FormatNumber(r.LfcSE)}</td><td>{FormatNumber(r.Stat)}</td><td>{r.PValue:G3}</td><td>{r.PAdj:G3}</td></tr>");
        }
        sb.Append("</table>");
        return Card("differential expression", $"{result.Method}: {result.Contrast}", chips, sb.ToString());
    }

    /// <summary>A molecule card: SVG depiction, SMILES, formula and key descriptors.</summary>
    public static string Molecule(Molecule mol)
    {
        ArgumentNullException.ThrowIfNull(mol);
        string svg;
        try
        {
            var viewer = new MoleculeViewer(mol, 360, 260) { Title = null, ShowFormula = false };
            svg = StripXml(viewer.ToSvg());
        }
        catch (Exception ex)
        {
            svg = $"<div class=\"muted\">(depiction unavailable: {Enc(ex.Message)})</div>";
        }
        var chips = new StringBuilder();
        string smiles = "";
        try
        {
            smiles = mol.ToSmiles();
            var d = MolecularDescriptors.Compute(mol);
            chips.Append(Chip(d.Formula)).Append(Chip($"MW {d.MolecularWeight:F2}")).Append(Chip($"logP {d.LogP:F2}")).Append(Chip($"TPSA {d.Tpsa:F1}"))
                 .Append(Chip($"HBD {d.HBondDonors} · HBA {d.HBondAcceptors}")).Append(Chip($"rot. bonds {d.RotatableBonds}"));
        }
        catch (Exception)
        {
            chips.Append(Chip($"{mol.AtomCount} atoms"));
        }
        string body = $"<div style=\"display:flex;gap:14px;align-items:flex-start;flex-wrap:wrap\"><div style=\"flex:0 0 auto;border-radius:8px;overflow:hidden\">{svg}</div>"
                    + $"<div style=\"min-width:200px;flex:1\"><div class=\"btn-k\">SMILES</div><pre class=\"btn-seq\" style=\"white-space:pre-wrap;word-break:break-all\">{Enc(smiles)}</pre>"
                    + $"<div class=\"muted\" style=\"font-size:11.5px\">{mol.AtomCount} atoms · {mol.BondCount} bonds · {mol.Rings.Count} rings</div></div></div>";
        return Card("molecule", string.IsNullOrEmpty(mol.Name) ? "molecule" : mol.Name, chips.ToString(), body);
    }

    /// <summary>A structure card with chain summary and an interactive 3-D viewer (3Dmol.js).</summary>
    public static string Structure(Structure structure, Viewer3DOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(structure);
        var chips = new StringBuilder();
        chips.Append(Chip($"{structure.Models.Count} model(s)")).Append(Chip($"{structure.Chains.Count} chain(s)")).Append(Chip($"{structure.AtomCount:N0} atoms"));
        int ligands = structure.Ligands.Count();
        if (ligands > 0) chips.Append(Chip($"{ligands} ligand(s)"));
        if (!string.IsNullOrEmpty(structure.Classification)) chips.Append(Chip(structure.Classification));
        var t = new StringBuilder("<table class=\"btn-tbl\"><tr><th class=\"l\">chain</th><th>residues</th><th class=\"l\">sequence</th></tr>");
        foreach (var c in structure.Chains.Take(12))
        {
            string seq;
            try { seq = c.ToSequence().ToString(); } catch (Exception) { seq = ""; }
            t.Append(CultureInfo.InvariantCulture, $"<tr><td class=\"l\">{Enc(c.Id)}</td><td>{c.Residues.Count}</td><td class=\"l\" style=\"font-family:monospace\">{Enc(seq.Length > 70 ? seq[..70] + "…" : seq)}</td></tr>");
        }
        t.Append("</table>");
        string viewer = ProteinViewer.To3DHtml(structure, options ?? new Viewer3DOptions { Width = 640, Height = 420 });
        string title = string.IsNullOrEmpty(structure.Title) ? structure.Id : $"{structure.Id} · {structure.Title}";
        return Card("structure", title, chips.ToString(), viewer + t);
    }

    /// <summary>Embeds a plot (Plotly HTML, SVG or PNG).</summary>
    public static string Plot(BioPlot plot, PlotEmbedding? embedding = null)
    {
        ArgumentNullException.ThrowIfNull(plot);
        var mode = embedding ?? DefaultPlotEmbedding;
        string content = mode switch
        {
            PlotEmbedding.Png => $"<img alt=\"{Enc(plot.Title ?? plot.GetType().Name)}\" src=\"{plot.ToPngDataUri(2)}\" style=\"width:{plot.Width}px;max-width:100%;height:auto\"/>",
            PlotEmbedding.Svg => plot.ToHtml(interactive: false),
            _ => plot.ToHtml(interactive: true),
        };
        return $"<div class=\"biotech-plot-output\">{content}<div style=\"font:10.5px 'Segoe UI',system-ui,sans-serif;opacity:.55;margin-top:2px\">{Enc(BiotechInfo.CreditsId)}</div></div>";
    }

    private static string StripXml(string svg)
    {
        int i = svg.IndexOf("<svg", StringComparison.Ordinal);
        return i >= 0 ? svg[i..] : svg;
    }

    private static string FormatNumber(double v)
        => double.IsNaN(v) ? "NA" : Math.Abs(v) >= 1e5 || (Math.Abs(v) < 1e-3 && v != 0) ? v.ToString("0.###E+0", Inv) : v.ToString("0.####", Inv);
}
