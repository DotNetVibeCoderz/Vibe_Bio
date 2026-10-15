using Biotech.Net.Chemistry;
using Matplotlib.Net.Rendering;
using Visual.Net.Core.Colors;

namespace Biotech.Net.Visualization;

/// <summary>
/// A 2-D structure depiction of a <see cref="Molecule"/>: Kekulé bonds (ring double bonds drawn inside the ring),
/// heteroatom labels with implicit hydrogens and charges, CPK colours. Coordinates come from the molecule's own
/// 2-D coordinates when present, otherwise from Chemistry's <see cref="Depiction"/> (computed on a copy).
/// </summary>
public sealed class MoleculeViewer : BioPlot
{
    private (double X, double Y)[]? _coords;

    public MoleculeViewer(Molecule molecule, int width = 420, int height = 320) : base(width, height)
    {
        Molecule = molecule ?? throw new ArgumentNullException(nameof(molecule));
        Title = string.IsNullOrEmpty(molecule.Name) ? null : molecule.Name;
    }

    /// <summary>The molecule.</summary>
    public Molecule Molecule { get; }

    /// <summary>Label carbon atoms too.</summary>
    public bool ShowCarbonLabels { get; set; }

    /// <summary>Show implicit hydrogens on labelled atoms (OH, NH₂).</summary>
    public bool ShowImplicitHydrogens { get; set; } = true;

    /// <summary>Print atom indices next to atoms.</summary>
    public bool ShowAtomIndices { get; set; }

    /// <summary>Colour heteroatoms by element (CPK).</summary>
    public bool ColorByElement { get; set; } = true;

    /// <summary>Show the molecular formula under the drawing.</summary>
    public bool ShowFormula { get; set; } = true;

    /// <summary>Atoms to highlight (e.g. a substructure match).</summary>
    public IReadOnlyCollection<int>? HighlightAtoms { get; set; }

    /// <summary>2-D coordinates used for drawing (y down is not applied; molecule units).</summary>
    public (double X, double Y)[] Coordinates2D()
    {
        if (_coords is not null) return _coords;
        Point3[] c;
        if (Molecule.Coordinates is { } own && !Molecule.Has3DCoordinates && own.Any(p => p.X != 0 || p.Y != 0)) c = own;
        else c = Molecule.AtomCount == 0 ? [] : Depiction.Compute2DCoordinates(Molecule.Clone());
        return _coords = c.Select(p => (p.X, p.Y)).ToArray();
    }

    private string Label(int i, out bool show)
    {
        var a = Molecule.Atoms[i];
        string sym = a.Symbol;
        show = a.AtomicNumber != 6 || ShowCarbonLabels || Molecule.Degree(i) == 0 || a.FormalCharge != 0 || a.Isotope != 0;
        if (!show) return string.Empty;
        string s = a.Isotope != 0 ? $"{a.Isotope}{sym}" : sym;
        int h = a.HydrogenCount;
        if (ShowImplicitHydrogens && h > 0 && a.AtomicNumber != 1) s += h == 1 ? "H" : "H" + Sub(h);
        if (a.FormalCharge != 0) s += (Math.Abs(a.FormalCharge) > 1 ? Math.Abs(a.FormalCharge).ToString() : "") + (a.FormalCharge > 0 ? "+" : "−");
        return s;
    }

    private static string Sub(int n) => string.Concat(n.ToString().Select(ch => (char)('₀' + (ch - '0'))));

    /// <inheritdoc />
    public override void Render(IRenderer r)
    {
        var th = Theme;
        r.Clear(th.Background);
        double top = DrawTitle(r, 10);
        var c = Coordinates2D();
        double bottomPad = ShowFormula ? 30 : 14;
        if (c.Length == 0)
        {
            r.DrawText("(empty molecule)", new Pt(Width / 2.0, Height / 2.0), th.Font(), th.Muted, TextAlign.Center, TextBaseline.Middle);
            return;
        }
        double minX = c.Min(p => p.X), maxX = c.Max(p => p.X), minY = c.Min(p => p.Y), maxY = c.Max(p => p.Y);
        var area = new Rect(24, top + 16, Width - 48, Height - top - 16 - bottomPad);
        double meanBond = Molecule.BondCount == 0 ? 1.5 : Enumerable.Range(0, Molecule.BondCount)
            .Average(b => Dist(c[Molecule.Bonds[b].Begin], c[Molecule.Bonds[b].End]));
        double scale = Math.Min(area.Width / Math.Max(1e-6, maxX - minX), area.Height / Math.Max(1e-6, maxY - minY));
        scale = Math.Min(scale, 42 / Math.Max(1e-6, meanBond));
        double ox = area.X + area.Width / 2 - (minX + maxX) / 2 * scale;
        double oy = area.Y + area.Height / 2 + (minY + maxY) / 2 * scale;
        var P = c.Select(p => new Pt(ox + p.X * scale, oy - p.Y * scale)).ToArray();
        double bondPx = meanBond * scale;
        var font = th.Font(Math.Clamp(bondPx * 0.42, 9, 17), bold: false);

        var labels = new string[c.Length];
        var shown = new bool[c.Length];
        for (int i = 0; i < c.Length; i++) labels[i] = Label(i, out shown[i]);

        if (HighlightAtoms is { Count: > 0 } hl)
            foreach (int i in hl.Where(i => i >= 0 && i < P.Length))
                r.DrawEllipse(new Rect(P[i].X - bondPx * 0.32, P[i].Y - bondPx * 0.32, bondPx * 0.64, bondPx * 0.64), new Brush(th.CategoryColor(1).WithAlpha(0.35)), Pen.None);

        var rings = Molecule.Rings;
        double width = Math.Clamp(bondPx * 0.045, 1.2, 2.2);
        for (int b = 0; b < Molecule.BondCount; b++)
        {
            var bond = Molecule.Bonds[b];
            var (p1, p2) = Trim(P[bond.Begin], P[bond.End], shown[bond.Begin] ? font.Size * 0.62 : 0, shown[bond.End] ? font.Size * 0.62 : 0);
            int order = bond.Order == BondOrder.Aromatic ? SafeKekule(b) : (int)bond.Order;
            var col1 = AtomColor(bond.Begin);
            var col2 = AtomColor(bond.End);
            void Segment(Pt a, Pt z)
            {
                var mid = new Pt((a.X + z.X) / 2, (a.Y + z.Y) / 2);
                r.DrawLine(a, mid, Paint.Pen(col1, width));
                r.DrawLine(mid, z, Paint.Pen(col2, width));
            }
            double dx = p2.X - p1.X, dy = p2.Y - p1.Y, len = Math.Sqrt(dx * dx + dy * dy);
            if (len < 1e-6) continue;
            double nx = -dy / len, ny = dx / len, off = bondPx * 0.16;
            switch (order)
            {
                case 2:
                {
                    // Ring double bond: second line inside the ring; acyclic: centred pair.
                    var ring = rings.AtomRings.FirstOrDefault(rg => rg.Contains(bond.Begin) && rg.Contains(bond.End));
                    if (ring is not null)
                    {
                        double cx = ring.Average(i => P[i].X), cy = ring.Average(i => P[i].Y);
                        double side = Math.Sign((cx - p1.X) * nx + (cy - p1.Y) * ny);
                        if (side == 0) side = 1;
                        Segment(p1, p2);
                        double sh = len * 0.12;
                        var a = new Pt(p1.X + dx / len * sh + nx * off * side, p1.Y + dy / len * sh + ny * off * side);
                        var z = new Pt(p2.X - dx / len * sh + nx * off * side, p2.Y - dy / len * sh + ny * off * side);
                        Segment(a, z);
                    }
                    else
                    {
                        Segment(new Pt(p1.X + nx * off / 2, p1.Y + ny * off / 2), new Pt(p2.X + nx * off / 2, p2.Y + ny * off / 2));
                        Segment(new Pt(p1.X - nx * off / 2, p1.Y - ny * off / 2), new Pt(p2.X - nx * off / 2, p2.Y - ny * off / 2));
                    }
                    break;
                }
                case 3:
                    Segment(p1, p2);
                    Segment(new Pt(p1.X + nx * off, p1.Y + ny * off), new Pt(p2.X + nx * off, p2.Y + ny * off));
                    Segment(new Pt(p1.X - nx * off, p1.Y - ny * off), new Pt(p2.X - nx * off, p2.Y - ny * off));
                    break;
                default:
                    Segment(p1, p2);
                    break;
            }
        }

        for (int i = 0; i < P.Length; i++)
        {
            if (shown[i])
            {
                var m = r.MeasureText(labels[i], font);
                r.DrawRect(new Rect(P[i].X - m.Width / 2 - 1, P[i].Y - m.Height / 2, m.Width + 2, m.Height), new Brush(th.Background), Pen.None);
                r.DrawText(labels[i], P[i], font, AtomColor(i), TextAlign.Center, TextBaseline.Middle);
            }
            if (ShowAtomIndices)
                r.DrawText(i.ToString(), new Pt(P[i].X + 6, P[i].Y + 6), th.Font(font.Size * 0.6), th.Muted, TextAlign.Left, TextBaseline.Top);
        }

        if (ShowFormula)
        {
            string formula;
            try { formula = MolecularDescriptors.Formula(Molecule); } catch (Exception) { formula = string.Empty; }
            r.DrawText(formula, new Pt(Width / 2.0, Height - 10), th.SmallFont, th.Muted, TextAlign.Center, TextBaseline.Bottom);
        }
    }

    private int SafeKekule(int bond)
    {
        try { return Molecule.KekuleOrder(bond); } catch (Exception) { return 1; }
    }

    private Color AtomColor(int i) => ColorByElement ? Theme.ElementColor(Molecule.Atoms[i].Symbol) : Theme.Foreground;

    private static double Dist((double X, double Y) a, (double X, double Y) b) => Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));

    private static (Pt, Pt) Trim(Pt a, Pt b, double ta, double tb)
    {
        double dx = b.X - a.X, dy = b.Y - a.Y, len = Math.Sqrt(dx * dx + dy * dy);
        if (len < 1e-9 || ta + tb >= len * 0.9) return (a, b);
        return (new Pt(a.X + dx / len * ta, a.Y + dy / len * ta), new Pt(b.X - dx / len * tb, b.Y - dy / len * tb));
    }

    /// <summary>
    /// An interactive 3Dmol.js viewer (ball-and-stick) of the molecule; uses its 3-D coordinates when present,
    /// otherwise the flat 2-D depiction.
    /// </summary>
    public string To3DHtml(Viewer3DOptions? options = null)
    {
        var o = options ?? new Viewer3DOptions { Width = Width, Height = Height, Style = MolecularStyle.Stick, Theme = Theme };
        var copy = Molecule.Clone();
        if (copy.Coordinates is null) Depiction.Compute2DCoordinates(copy);
        return Viewer3D.Html(MolFile.Write(copy), "sdf", o, Molecule.Name ?? "molecule");
    }
}
