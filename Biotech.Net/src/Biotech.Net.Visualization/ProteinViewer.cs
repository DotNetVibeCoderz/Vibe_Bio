using System.Globalization;
using System.Net;
using System.Text;
using System.Text.Json;
using Biotech.Net.StructuralBiology;
using Matplotlib.Net;
using Matplotlib.Net.Artists;
using Matplotlib.Net.Rendering;
using Plotly.Net;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>Colouring of residues in protein views.</summary>
public enum ResidueColoring
{
    /// <summary>N-to-C rainbow along each chain.</summary>
    Spectrum,

    /// <summary>One colour per chain.</summary>
    Chain,

    /// <summary>Secondary structure from DSSP (helix, strand, loop).</summary>
    SecondaryStructure,

    /// <summary>Cα B-factor (blue = ordered, red = mobile).</summary>
    BFactor,

    /// <summary>Amino-acid physico-chemical class.</summary>
    ResidueType,
}

/// <summary>Representation in the 3-D viewer.</summary>
public enum MolecularStyle
{
    Cartoon,
    Stick,
    Sphere,
    Line,
}

/// <summary>Options of the 3Dmol.js viewer.</summary>
public sealed record Viewer3DOptions
{
    /// <summary>Viewer width in pixels.</summary>
    public int Width { get; init; } = 640;

    /// <summary>Viewer height in pixels.</summary>
    public int Height { get; init; } = 480;

    /// <summary>Main representation.</summary>
    public MolecularStyle Style { get; init; } = MolecularStyle.Cartoon;

    /// <summary>Colouring of the main representation.</summary>
    public ResidueColoring Coloring { get; init; } = ResidueColoring.Spectrum;

    /// <summary>Show ligands (HETATM, not water) as sticks.</summary>
    public bool ShowLigands { get; init; } = true;

    /// <summary>Show water molecules.</summary>
    public bool ShowWater { get; init; }

    /// <summary>Slowly spin the model.</summary>
    public bool Spin { get; init; }

    /// <summary>Theme for the background (null: <see cref="BioTheme.Default"/>).</summary>
    public BioTheme? Theme { get; init; }
}

/// <summary>
/// Protein structure views: a 2-D Cα trace, a residue contact / distance map, and an interactive 3-D viewer built
/// on 3Dmol.js (loaded from a pinned CDN version).
/// </summary>
public static class ProteinViewer
{
    /// <summary>Pinned 3Dmol.js version.</summary>
    public const string ThreeDmolVersion = "2.4.2";

    /// <summary>CDN URL of 3Dmol.js.</summary>
    public const string ThreeDmolUrl = "https://cdn.jsdelivr.net/npm/3dmol@" + ThreeDmolVersion + "/build/3Dmol-min.js";

    /// <summary>A 2-D Cα trace of the structure (or one chain).</summary>
    public static CaTracePlot CaTrace(Structure structure, string? chainId = null) => new(structure, chainId);

    /// <summary>A residue contact map of one chain (the first protein chain by default).</summary>
    public static ContactMapPlot ContactMap(Structure structure, string? chainId = null, double cutoff = 8.0) => new(structure, chainId, cutoff);

    /// <summary>An HTML fragment with an interactive 3-D viewer of the structure.</summary>
    public static string To3DHtml(Structure structure, Viewer3DOptions? options = null)
        => Viewer3D.Html(PdbWriter.ToPdb(structure, new PdbWriteOptions { IncludeConnections = true }), "pdb", options ?? new(), structure.Id);
}

/// <summary>Shared 3Dmol.js HTML generation (also used for small molecules).</summary>
internal static class Viewer3D
{
    public static string Html(string modelText, string format, Viewer3DOptions o, string title)
    {
        var th = o.Theme ?? BioTheme.Default;
        string id = "biotech3d_" + Guid.NewGuid().ToString("N")[..12];
        string color = o.Coloring switch
        {
            ResidueColoring.Chain => "{colorscheme:'chain'}",
            ResidueColoring.SecondaryStructure => "{colorscheme:'ssJmol'}",
            ResidueColoring.BFactor => "{colorscheme:{prop:'b',gradient:'rwb',min:80,max:0}}",
            ResidueColoring.ResidueType => "{colorscheme:'amino'}",
            _ => "{color:'spectrum'}",
        };
        string style = o.Style switch
        {
            MolecularStyle.Stick => $"{{stick:{color}}}",
            MolecularStyle.Sphere => $"{{sphere:{color}}}",
            MolecularStyle.Line => $"{{line:{color}}}",
            _ => format == "pdb" ? $"{{cartoon:{color}}}" : "{stick:{colorscheme:'Jmol'}}",
        };
        var sb = new StringBuilder();
        sb.Append(CultureInfo.InvariantCulture, $"<div class=\"biotech-3d\" style=\"font-family:'Segoe UI',system-ui,sans-serif\">");
        sb.Append(CultureInfo.InvariantCulture, $"<div id=\"{id}\" style=\"width:{o.Width}px;max-width:100%;height:{o.Height}px;position:relative;border:1px solid {th.Border.ToHex()};border-radius:8px;overflow:hidden\" title=\"{WebUtility.HtmlEncode(title)}\"></div>");
        sb.Append("<script type=\"text/javascript\">(function(){");
        sb.Append("var data=").Append(JsonSerializer.Serialize(modelText)).Append(';');
        sb.Append("function draw($3Dmol){var el=document.getElementById('").Append(id).Append("');if(!el)return;");
        sb.Append("var v=$3Dmol.createViewer(el,{backgroundColor:'").Append(th.Background.ToHex()).Append("'});");
        sb.Append("v.addModel(data,'").Append(format).Append("');");
        sb.Append("v.setStyle({},").Append(style).Append(");");
        if (format == "pdb")
        {
            if (o.ShowLigands) sb.Append("v.setStyle({hetflag:true,not:{resn:['HOH','WAT']}},{stick:{colorscheme:'greenCarbon',radius:0.2}});");
            sb.Append(o.ShowWater ? "v.setStyle({resn:['HOH','WAT']},{sphere:{radius:0.25,color:'#7FB2E5'}});" : "v.setStyle({resn:['HOH','WAT']},{});");
        }
        sb.Append("v.zoomTo();v.render();");
        if (o.Spin) sb.Append("v.spin('y',0.5);");
        sb.Append('}');
        sb.Append("if(window.$3Dmol){draw(window.$3Dmol);}");
        sb.Append("else if(typeof require==='function'&&typeof define==='function'&&define.amd){require.config({paths:{'3Dmol':'")
          .Append(ProteinViewer.ThreeDmolUrl[..^3]).Append("'}});require(['3Dmol'],function(m){draw(m||window.$3Dmol);});}");
        sb.Append("else{var s=document.createElement('script');s.src='").Append(ProteinViewer.ThreeDmolUrl)
          .Append("';s.onload=function(){draw(window.$3Dmol);};document.head.appendChild(s);}");
        sb.Append("})();</script></div>");
        return sb.ToString();
    }
}

/// <summary>
/// A 2-D Cα trace: residues projected onto the structure's two principal axes, connected along each chain, with
/// optional secondary-structure or B-factor colouring and residue numbers.
/// </summary>
public sealed class CaTracePlot : BioPlot
{
    private readonly List<(Residue Residue, int Chain)> _residues;

    public CaTracePlot(Structure structure, string? chainId = null) : base(620, 560)
    {
        Structure = structure ?? throw new ArgumentNullException(nameof(structure));
        var chains = chainId is null ? structure.Chains.ToList() : [structure[chainId]];
        _residues = chains.SelectMany((c, ci) => c.AminoAcids().Select(r => (r, ci))).ToList();
        if (_residues.Count == 0) throw new ArgumentException("The structure has no amino-acid residues with Cα atoms.");
        Title = string.IsNullOrEmpty(structure.Title) ? structure.Id : $"{structure.Id} · {structure.Title}";
    }

    /// <summary>The structure.</summary>
    public Structure Structure { get; }

    /// <summary>Residue colouring.</summary>
    public ResidueColoring Coloring { get; set; } = ResidueColoring.Spectrum;

    /// <summary>Label every n-th residue (0 = none).</summary>
    public int LabelEvery { get; set; } = 10;

    /// <summary>Number of residues drawn.</summary>
    public int ResidueCount => _residues.Count;

    /// <summary>
    /// Cα positions projected on the first two principal axes of the Cα cloud (deterministic: axes are signed so the
    /// first residue lies at negative coordinates).
    /// </summary>
    public (double X, double Y)[] Project()
    {
        var p = _residues.Select(t => t.Residue["CA"].Position).ToArray();
        var c = new Vec3(p.Average(v => v.X), p.Average(v => v.Y), p.Average(v => v.Z));
        double[,] cov = new double[3, 3];
        foreach (var v in p)
        {
            double[] d = [v.X - c.X, v.Y - c.Y, v.Z - c.Z];
            for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cov[i, j] += d[i] * d[j];
        }
        var e1 = PowerIteration(cov, null);
        var e2 = PowerIteration(cov, e1);
        var first = p[0] - c;
        if (Dot(first, e1) > 0) e1 = Neg(e1);
        if (Dot(first, e2) > 0) e2 = Neg(e2);
        return p.Select(v => { var d = v - c; return (Dot(d, e1), Dot(d, e2)); }).ToArray();
    }

    private static double Dot(Vec3 v, double[] e) => v.X * e[0] + v.Y * e[1] + v.Z * e[2];
    private static double[] Neg(double[] e) => [-e[0], -e[1], -e[2]];

    private static double[] PowerIteration(double[,] m, double[]? deflate)
    {
        double[] v = deflate is null ? [1, 0.7, 0.3] : [-deflate[1], deflate[0] + 0.3, 0.5];
        for (int it = 0; it < 200; it++)
        {
            var w = new double[3];
            for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) w[i] += m[i, j] * v[j];
            if (deflate is not null)
            {
                double d = w[0] * deflate[0] + w[1] * deflate[1] + w[2] * deflate[2];
                for (int i = 0; i < 3; i++) w[i] -= d * deflate[i];
            }
            double n = Math.Sqrt(w.Sum(x => x * x));
            if (n < 1e-12) break;
            for (int i = 0; i < 3; i++) v[i] = w[i] / n;
        }
        return v;
    }

    private Color[] ResidueColors()
    {
        var th = Theme;
        var spectrum = Visual.Net.Core.Colors.Colormaps.Get(BiotechColormaps.Spectrum);
        switch (Coloring)
        {
            case ResidueColoring.Chain:
                return _residues.Select(t => th.CategoryColor(t.Chain)).ToArray();
            case ResidueColoring.ResidueType:
                return _residues.Select(t => th.AminoAcidColor(t.Residue.OneLetterCode)).ToArray();
            case ResidueColoring.BFactor:
            {
                var b = _residues.Select(t => t.Residue["CA"].BFactor).ToArray();
                double lo = b.Min(), hi = Math.Max(b.Max(), lo + 1e-9);
                var cmap = Visual.Net.Core.Colors.Colormaps.Get(BiotechColormaps.Diverging);
                return b.Select(v => cmap.At((v - lo) / (hi - lo))).ToArray();
            }
            case ResidueColoring.SecondaryStructure:
            {
                DsspResult? dssp = null;
                try { dssp = Dssp.Assign(Structure); } catch (Exception) { /* fall back to loop colouring */ }
                return _residues.Select(t => (dssp?[t.Residue]?.Structure ?? SecondaryStructure.Loop) switch
                {
                    SecondaryStructure.AlphaHelix or SecondaryStructure.Helix310 or SecondaryStructure.PiHelix => th.Up,
                    SecondaryStructure.Strand or SecondaryStructure.Bridge => th.CategoryColor(1),
                    _ => th.Neutral,
                }).ToArray();
            }
            default:
            {
                var colors = new Color[_residues.Count];
                foreach (var g in _residues.Select((t, i) => (t.Chain, i)).GroupBy(t => t.Chain))
                {
                    var idx = g.Select(t => t.i).ToArray();
                    for (int k = 0; k < idx.Length; k++) colors[idx[k]] = spectrum.At(idx.Length == 1 ? 0 : (double)k / (idx.Length - 1));
                }
                return colors;
            }
        }
    }

    /// <inheritdoc />
    public override void Render(IRenderer r)
    {
        var th = Theme;
        r.Clear(th.Background);
        double top = DrawTitle(r, 10) + 10;
        var pts = Project();
        var colors = ResidueColors();
        double minX = pts.Min(p => p.X), maxX = pts.Max(p => p.X), minY = pts.Min(p => p.Y), maxY = pts.Max(p => p.Y);
        var area = new Rect(24, top, Width - 48, Height - top - 40);
        double scale = Math.Min(area.Width / Math.Max(1e-6, maxX - minX), area.Height / Math.Max(1e-6, maxY - minY));
        double ox = area.X + (area.Width - (maxX - minX) * scale) / 2, oy = area.Y + (area.Height - (maxY - minY) * scale) / 2;
        Pt Map((double X, double Y) p) => new(ox + (p.X - minX) * scale, oy + (p.Y - minY) * scale);

        for (int i = 1; i < pts.Length; i++)
        {
            if (_residues[i].Chain != _residues[i - 1].Chain) continue;
            double dist = Vec3.Distance(_residues[i].Residue["CA"].Position, _residues[i - 1].Residue["CA"].Position);
            var pen = dist > 4.3 ? Paint.Pen(th.Muted, 1, LineStyle.Dashed) : Paint.Pen(Color.Lerp(colors[i - 1], colors[i], 0.5), 3.2);
            r.DrawLine(Map(pts[i - 1]), Map(pts[i]), pen);
        }
        for (int i = 0; i < pts.Length; i++)
        {
            var c = Map(pts[i]);
            r.DrawEllipse(new Rect(c.X - 3.2, c.Y - 3.2, 6.4, 6.4), new Brush(colors[i]), Paint.Pen(th.Background, 0.8));
            if (LabelEvery > 0 && (_residues[i].Residue.Number % LabelEvery == 0 || i == 0 || i == pts.Length - 1))
                r.DrawText($"{_residues[i].Residue.OneLetterCode}{_residues[i].Residue.Number}", new Pt(c.X + 5, c.Y - 5), th.Font(th.FontSize * 0.72), th.Muted, TextAlign.Left, TextBaseline.Bottom);
        }
        string caption = $"Cα trace · {pts.Length} residues · projection on principal axes · colour: {Coloring}";
        r.DrawText(caption, new Pt(Width / 2.0, Height - 14), th.SmallFont, th.Muted, TextAlign.Center, TextBaseline.Bottom);
    }

    /// <inheritdoc />
    public override string ToHtml(bool interactive = true)
        => interactive ? ProteinViewer.To3DHtml(Structure, new Viewer3DOptions { Width = Width, Height = Height, Coloring = Coloring, Theme = Theme }) : base.ToHtml(false);
}

/// <summary>A residue contact map (binary at a cutoff) or Cα distance map of one chain.</summary>
public sealed class ContactMapPlot : FigurePlot
{
    public ContactMapPlot(Structure structure, string? chainId = null, double cutoff = 8.0) : base(600, 560)
    {
        ArgumentNullException.ThrowIfNull(structure);
        Chain = chainId is null
            ? structure.Chains.FirstOrDefault(c => c.AminoAcids().Count > 0) ?? throw new ArgumentException("No protein chain.")
            : structure[chainId];
        Cutoff = cutoff;
        Title = $"{structure.Id} chain {Chain.Id}";
    }

    /// <summary>The chain.</summary>
    public Chain Chain { get; }

    /// <summary>Contact cutoff in ångström.</summary>
    public double Cutoff { get; set; }

    /// <summary>Show Cα–Cα distances instead of binary contacts.</summary>
    public bool ShowDistances { get; set; }

    /// <summary>Contact (1/0) or distance matrix, residues in chain order.</summary>
    public double[,] Matrix()
    {
        var res = Chain.AminoAcids();
        int n = res.Count;
        var m = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double d = Vec3.Distance(res[i]["CA"].Position, res[j]["CA"].Position);
                m[i, j] = ShowDistances ? d : d <= Cutoff ? 1 : 0;
            }
        return m;
    }

    protected override void Build(Figure figure)
    {
        var res = Chain.AminoAcids();
        int n = res.Count;
        var m = Matrix();
        var ax = figure.AddSubplot();
        int first = res.Count > 0 ? res[0].Number : 1;
        ax.ImShow(m, ShowDistances ? Theme.SequentialColormap + "" : BiotechColormaps.Blues,
            extent: (first - 0.5, first + n - 0.5, first + n - 0.5, first - 0.5));
        if (ShowDistances) ax.AddColorbar("Cα–Cα distance (Å)");
        ax.SetAspect(1);
        ax.Grid(false);
        ax.SetXLabel(XLabel ?? "residue");
        ax.SetYLabel(YLabel ?? "residue");
        if (Title is not null) ax.SetTitle(ShowDistances ? Title + " · distance map" : $"{Title} · contacts ≤ {Cutoff:0.#} Å");
    }

    /// <inheritdoc />
    public override bool SupportsInteractive => true;

    /// <inheritdoc />
    public override PlotlyFigure ToPlotly()
    {
        var res = Chain.AminoAcids();
        var m = Matrix();
        int n = res.Count;
        var labels = res.Select(r => $"{r.OneLetterCode}{r.Number}").ToArray();
        var f = NewPlotly();
        f.Add(new Plotly.Net.Traces.Heatmap { Z = m, X = labels, Y = labels, ColorScale = ShowDistances ? "Viridis" : "Blues" });
        f.Layout.YAxis.AutoRangeReversed = true;
        return f;
    }
}
