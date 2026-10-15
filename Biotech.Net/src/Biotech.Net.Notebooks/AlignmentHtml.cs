using System.Globalization;
using System.Text;
using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Phylogenetics;

namespace Biotech.Net.Notebooks;

/// <summary>HTML renderings of alignment results (pairwise, multiple, phylogenetic trees).</summary>
public static class AlignmentHtml
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Maximum alignment columns shown.</summary>
    public static int MaxColumns { get; set; } = 1200;

    /// <summary>A pairwise alignment card: score, identity, gaps, CIGAR and coloured blocks.</summary>
    public static string Pairwise(PairwiseAlignment a, int width = 60)
    {
        ArgumentNullException.ThrowIfNull(a);
        bool nuc = a.Target.Alphabet.IsNucleotide;
        var chips = Chip(a.Mode.ToString()) + Chip($"score {a.Score}") + Chip($"identity {a.Identity * 100:F1}%")
                    + Chip($"{a.Identities}/{a.Length} identical") + Chip($"{a.Gaps} gaps");
        string t = a.AlignedTarget, q = a.AlignedQuery;
        int len = Math.Min(t.Length, MaxColumns);
        var sb = new StringBuilder("<pre class=\"btn-seq\">");
        int tPos = a.TargetStart, qPos = a.QueryStart;
        int nameW = Math.Min(14, Math.Max(a.Target.Id.Length, a.Query.Id.Length));
        for (int start = 0; start < len; start += width)
        {
            int end = Math.Min(len, start + width);
            var mid = new StringBuilder();
            for (int i = start; i < end; i++)
                mid.Append(t[i] == '-' || q[i] == '-' ? ' ' : char.ToUpperInvariant(t[i]) == char.ToUpperInvariant(q[i]) ? '|' : '.');
            int tAdv = t[start..end].Count(c => c != '-'), qAdv = q[start..end].Count(c => c != '-');
            Row(sb, Name(a.Target.Id, nameW), tPos + 1, t, start, end, q, nuc, tPos + tAdv);
            sb.Append("<span class=\"n\">").Append(new string(' ', nameW + 9)).Append("</span>").Append(BiotechHtml.Enc(mid.ToString())).Append('\n');
            Row(sb, Name(a.Query.Id, nameW), qPos + 1, q, start, end, t, nuc, qPos + qAdv);
            sb.Append('\n');
            tPos += tAdv; qPos += qAdv;
        }
        if (len < t.Length) sb.Append("<span class=\"n\">… ").Append((t.Length - len).ToString("N0", Inv)).Append(" more columns</span>");
        sb.Append("</pre>");
        string cigar = a.Cigar();
        string foot = $"<div class=\"muted\" style=\"font-size:11.5px\">target {a.TargetStart + 1}–{a.TargetEnd} · query {a.QueryStart + 1}–{a.QueryEnd} · CIGAR {BiotechHtml.Enc(cigar.Length > 120 ? cigar[..120] + "…" : cigar)}</div>";
        return Card("pairwise alignment", $"{a.Target.Id} × {a.Query.Id}", chips, sb + foot);
    }

    private static void Row(StringBuilder sb, string name, int from, string row, int start, int end, string other, bool nuc, int to)
    {
        sb.Append("<span class=\"n\">").Append(BiotechHtml.Enc(name)).Append(' ').Append(from.ToString(Inv).PadLeft(7)).Append(' ').Append("</span>");
        for (int i = start; i < end; i++)
        {
            char c = row[i];
            bool match = c != '-' && char.ToUpperInvariant(c) == char.ToUpperInvariant(other[i]);
            if (c == '-') sb.Append("<span class=\"muted\">-</span>");
            else if (match) sb.Append("<span class=\"").Append(BiotechHtml.ResidueClass(c, nuc)).Append("\">").Append(Safe(c)).Append("</span>");
            else sb.Append("<span class=\"").Append(BiotechHtml.ResidueClass(c, nuc)).Append("\" style=\"opacity:.55\">").Append(Safe(c)).Append("</span>");
        }
        sb.Append("<span class=\"n\"> ").Append(to.ToString(Inv)).Append("</span>\n");
    }

    /// <summary>A multiple sequence alignment: coloured columns, conservation track and consensus.</summary>
    public static string Multiple(MultipleSequenceAlignment msa, int width = 60, int maxSequences = 60)
    {
        ArgumentNullException.ThrowIfNull(msa);
        bool nuc = msa.Count > 0 && msa.IsNucleotide;
        var chips = Chip($"{msa.Count} sequences") + Chip($"{msa.Length:N0} columns") + Chip(nuc ? "nucleotide" : "protein");
        int rows = Math.Min(msa.Count, maxSequences), len = Math.Min(msa.Length, MaxColumns);
        int nameW = Math.Clamp(msa.Ids.Take(rows).Select(i => i.Length).DefaultIfEmpty(4).Max(), 9, 16);
        string conservation = msa.Length > 0 ? msa.ConservationLine() : "";
        string consensus = msa.Length > 0 ? msa.Consensus() : "";
        var sb = new StringBuilder("<pre class=\"btn-seq\">");
        for (int start = 0; start < len; start += width)
        {
            int end = Math.Min(len, start + width);
            sb.Append("<span class=\"n\">").Append(new string(' ', nameW + 1)).Append((start + 1).ToString(Inv)).Append("</span>\n");
            for (int r = 0; r < rows; r++)
            {
                sb.Append("<span class=\"n\">").Append(BiotechHtml.Enc(Name(msa.Ids[r], nameW))).Append(' ').Append("</span>");
                string row = msa[r];
                for (int i = start; i < end; i++)
                {
                    char c = row[i];
                    if (c == '-') sb.Append("<span class=\"muted\">-</span>");
                    else sb.Append("<span class=\"").Append(BiotechHtml.ResidueClass(c, nuc)).Append("\">").Append(Safe(c)).Append("</span>");
                }
                sb.Append('\n');
            }
            if (conservation.Length >= end)
                sb.Append("<span class=\"n\">").Append(new string(' ', nameW + 1)).Append("</span>").Append(BiotechHtml.Enc(conservation[start..end])).Append('\n');
            if (consensus.Length >= end)
                sb.Append("<span class=\"n\">").Append(Name("consensus", nameW)).Append(' ').Append("</span><b>").Append(BiotechHtml.Enc(consensus[start..end])).Append("</b>\n");
            sb.Append('\n');
        }
        if (rows < msa.Count) sb.Append("<span class=\"n\">… ").Append(msa.Count - rows).Append(" more sequences</span>\n");
        if (len < msa.Length) sb.Append("<span class=\"n\">… ").Append((msa.Length - len).ToString("N0", Inv)).Append(" more columns</span>");
        sb.Append("</pre>");
        return Card("multiple alignment", $"{msa.Count} × {msa.Length}", chips, sb.ToString());
    }

    /// <summary>A phylogenetic tree as an inline SVG phylogram (cladogram when branch lengths are missing) plus Newick.</summary>
    public static string Tree(PhylogeneticTree tree, int width = 640)
    {
        ArgumentNullException.ThrowIfNull(tree);
        var leaves = tree.Leaves;
        int rowH = 22, height = Math.Max(60, leaves.Count * rowH + 20);
        bool lengths = tree.Root.PreOrder().Any(n => n != tree.Root && n.BranchLength is > 0);
        var depth = new Dictionary<TreeNode, double>();
        var y = new Dictionary<TreeNode, double>();
        int leafIndex = 0;
        double maxDepth = 0;
        void Layout(TreeNode n, double d)
        {
            depth[n] = d;
            maxDepth = Math.Max(maxDepth, d);
            if (n.IsLeaf) { y[n] = 10 + (leafIndex++ + 0.5) * rowH; return; }
            foreach (var c in n.Children) Layout(c, d + (lengths ? Math.Max(0, c.BranchLength ?? 0) : 1));
            y[n] = (y[n.Children[0]] + y[n.Children[^1]]) / 2;
        }
        Layout(tree.Root, 0);
        if (!lengths) foreach (var leaf in leaves) depth[leaf] = maxDepth; // cladogram: align tips
        double labelW = Math.Min(220, 12 + 7 * leaves.Select(l => (l.Name ?? "").Length).DefaultIfEmpty(4).Max());
        double scale = (width - labelW - 30) / Math.Max(1e-9, maxDepth);
        double X(TreeNode n) => 12 + depth[n] * scale;
        var svg = new StringBuilder();
        svg.Append(CultureInfo.InvariantCulture, $"<svg role=\"img\" aria-label=\"phylogenetic tree\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\" style=\"max-width:100%;font-family:'Segoe UI',system-ui,sans-serif;font-size:12px\">");
        string stroke = "stroke=\"#5B6573\" stroke-width=\"1.4\" fill=\"none\" stroke-linecap=\"round\"";
        foreach (var n in tree.Root.PreOrder())
        {
            if (n.IsLeaf)
            {
                svg.Append(CultureInfo.InvariantCulture, $"<circle cx=\"{X(n):F1}\" cy=\"{y[n]:F1}\" r=\"2.6\" fill=\"#2F6FDB\"/>");
                svg.Append(CultureInfo.InvariantCulture, $"<text x=\"{X(n) + 7:F1}\" y=\"{y[n] + 4:F1}\" fill=\"currentColor\">{BiotechHtml.Enc(n.Name)}</text>");
                continue;
            }
            svg.Append(CultureInfo.InvariantCulture, $"<path {stroke} d=\"M{X(n):F1},{y[n.Children[0]]:F1}V{y[n.Children[^1]]:F1}\"/>");
            foreach (var c in n.Children)
                svg.Append(CultureInfo.InvariantCulture, $"<path {stroke} d=\"M{X(n):F1},{y[c]:F1}H{X(c):F1}\"/>");
        }
        svg.Append("</svg>");
        string newick = tree.ToNewick(decimals: 4);
        var chips = Chip($"{leaves.Count} leaves") + Chip(lengths ? $"total branch length {tree.TotalBranchLength:G4}" : "cladogram");
        string body = svg + $"<details><summary class=\"muted\" style=\"font-size:11.5px;cursor:pointer\">Newick</summary><pre class=\"btn-seq\" style=\"white-space:pre-wrap;word-break:break-all\">{BiotechHtml.Enc(newick)}</pre></details>";
        return Card("phylogenetic tree", tree.Root.Name ?? "tree", chips, body);
    }

    private static string Card(string kind, string title, string chips, string body) => BiotechHtml.Card(kind, title, chips, body);
    private static string Chip(string text) => BiotechHtml.Chip(text);
    private static string Name(string id, int w) => (id.Length > w ? id[..(w - 1)] + "…" : id).PadRight(w);
    private static string Safe(char c) => c is '<' or '>' or '&' ? "?" : c.ToString();
}
