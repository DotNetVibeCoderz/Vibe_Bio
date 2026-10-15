using System.Diagnostics;
using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Media;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Proteomics;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Visualization;

namespace Biotech.Net.Gallery.Pages;

public sealed class ProteinExplorerPage : GalleryPage
{
    public override string Title => "Protein Explorer";
    public override string Modules => "BIOTECH.NET.PROTEOMICS · BIOTECH.NET.STRUCTURALBIOLOGY · VISUALIZATION";
    public override string Lede => "ProtParam-style properties and a Kyte–Doolittle hydropathy profile for any sequence, next to crambin (PDB 1CRN): Cα trace, contact map and DSSP secondary structure.";
    public override char Channel => 'T';

    private Structure _structure = null!;
    private TextBox _input = null!;
    private ComboBox _coloring = null!;
    private readonly ContentControl _analysis = new();
    private readonly ContentControl _props = new();
    private readonly ContentControl _profiles = new();
    private readonly ContentControl _structureView = new();
    private CaTracePlot? _trace;
    private string _resultText = "";

    public override BioPlot? ExportPlot => _trace;
    public override string ResultText => _resultText;

    protected override Control Build()
    {
        _structure = SampleStructures.Crambin();
        string crambin = _structure.Chains[0].ToSequence().ToString();
        _input = Ui.Input(crambin, "seq", 110);
        _coloring = new ComboBox
        {
            ItemsSource = Enum.GetValues<ResidueColoring>(),
            SelectedItem = ResidueColoring.SecondaryStructure,
            MinWidth = 170,
        };
        _coloring.SelectionChanged += (_, _) => RenderStructure();

        var input = Ui.Panel("Sequence", Ui.Column(10, _input,
            Ui.Row(6,
                Ui.Button("Crambin", () => { _input.Text = crambin; Analyze(); }, "chip"),
                Ui.Button("Insulin B", () => { _input.Text = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"; Analyze(); }, "chip"),
                Ui.Button("GPCR TM", () => { _input.Text = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"; Analyze(); }, "chip")),
            Ui.Button("Analyze", Analyze, "primary")));

        var left = Ui.Column(14, input, _props);
        var right = Ui.Column(14, _analysis,
            Ui.Panel("Structure · 1CRN crambin", _structureView, headerExtra: Ui.Row(8, Ui.Eyebrow("Colour"), _coloring)),
            _profiles);
        Analyze();
        RenderStructure();
        return Ui.Columns("340,*", 18, Ui.Scroll(left), Ui.Scroll(right));
    }

    private void Analyze()
    {
        try
        {
            var sw = Stopwatch.StartNew();
            string seq = new string((_input.Text ?? "").Where(char.IsLetter).ToArray()).ToUpperInvariant();
            if (seq.Length < 5) { Status = "Enter at least five residues."; return; }
            var pa = new ProteinAnalysis(seq);
            double mw = pa.MolecularWeight(), pi = pa.IsoelectricPoint(), gravy = pa.Gravy(), ii = pa.InstabilityIndex();
            var ext = pa.ExtinctionCoefficient();
            var ss = pa.SecondaryStructureFraction();
            var half = pa.HalfLife();
            int window = seq.Length >= 60 ? 19 : 9;
            var profile = pa.Profile(AminoAcidScale.KyteDoolittle, window);
            var tm = TransmembranePredictor.Predict(seq);
            var counts = pa.AminoAcidCounts();
            sw.Stop();

            var readouts = Ui.ReadoutRow(
                Ui.Readout("Length", seq.Length.ToString(CultureInfo.InvariantCulture), "aa"),
                Ui.Readout("Mass", (mw / 1000).ToString("0.00", CultureInfo.InvariantCulture), "kDa"),
                Ui.Readout("pI", pi.ToString("0.00", CultureInfo.InvariantCulture), "", Palette.DyeC),
                Ui.Readout("GRAVY", gravy.ToString("+0.000;-0.000", CultureInfo.InvariantCulture), "", gravy > 0 ? Palette.DyeT : Palette.DyeC),
                Ui.Readout("Instability", ii.ToString("0.0", CultureInfo.InvariantCulture), pa.IsStable ? "stable" : "unstable", pa.IsStable ? Palette.DyeA : Palette.DyeT),
                Ui.Readout("ε280", ext.Cystines.ToString("N0", CultureInfo.InvariantCulture), "M⁻¹cm⁻¹"));

            var hydro = new ProfileChart
            {
                Height = 190,
                Values = profile,
                Positive = Palette.DyeT,
                Negative = Palette.DyeC,
                References = [(1.6, "TM threshold 1.6")],
                Segments = tm.Select(t => (t.Start, t.End)).ToList(),
                Format = "0.0",
            };

            const string order = "AVLIMFWCGPSTNQYHKRDE";
            var comp = new BarChart
            {
                Height = 190,
                Bars = order.Select(a => new BarChart.Bar(a.ToString(), counts.TryGetValue(a, out var c) ? c : 0, ToColor(Palette.PlotTheme.AminoAcidColor(a)),
                    (counts.TryGetValue(a, out var c2) ? c2 : 0).ToString(CultureInfo.InvariantCulture))).ToList(),
            };

            var props = Ui.KeyValues(
            [
                ("Aliphatic index", pa.AliphaticIndex().ToString("0.0", CultureInfo.InvariantCulture), null),
                ("Aromaticity", pa.Aromaticity().ToString("0.000", CultureInfo.InvariantCulture), null),
                ("Net charge at pH 7.0", pa.NetCharge(7.0).ToString("+0.00;-0.00", CultureInfo.InvariantCulture), null),
                ("ε280 (reduced Cys)", ext.Reduced.ToString("N0", CultureInfo.InvariantCulture), null),
                ("Abs 0.1% (=1 g/L)", pa.Absorbance01Percent().ToString("0.000", CultureInfo.InvariantCulture), null),
                ("Helix / turn / sheet", $"{ss.Helix:0.00} / {ss.Turn:0.00} / {ss.Sheet:0.00}", null),
                ("Half-life (mammalian)", half.Text, null),
                ("TM segments (KD ≥ 1.6)", tm.Count == 0 ? "none" : string.Join(", ", tm.Select(t => $"{t.Start + 1}–{t.End}")), tm.Count > 0 ? Palette.DyeG : null),
            ]);

            _analysis.Content = Ui.Panel("ProtParam", readouts, $"computed in {GalleryData.Ms(sw)}");
            _props.Content = Ui.Panel("Properties", props);
            _profiles.Content = Ui.Columns("*,*", 14,
                Ui.Panel($"Hydropathy · Kyte–Doolittle", hydro, $"window {window} · shaded = TM"),
                Ui.Panel("Amino-acid composition", comp, "by physico-chemical class"));

            TraceSequence = GalleryData.BackTranslate(seq);
            TraceLabel = $"TRACE · back-translated {Math.Min(seq.Length, 46)} aa";
            Status = $"{seq.Length} aa · {mw:N1} Da · pI {pi:0.00} · GRAVY {gravy:0.000} · {tm.Count} TM segments · {GalleryData.Ms(sw)}";
            _resultText = pa.Report();
            Code = Snippet(seq, window);
        }
        catch (Exception ex)
        {
            Status = "Could not analyse the sequence: " + ex.Message;
        }
    }

    private void RenderStructure()
    {
        try
        {
            var coloring = _coloring.SelectedItem is ResidueColoring c ? c : ResidueColoring.SecondaryStructure;
            var trace = ProteinViewer.CaTrace(_structure);
            trace.Theme = Palette.PlotTheme;
            trace.Coloring = coloring;
            trace.Width = 400; trace.Height = 330;
            trace.LabelEvery = 10;
            var contact = ProteinViewer.ContactMap(_structure, cutoff: 8.0);
            contact.Theme = Palette.PlotTheme;
            contact.Width = 420; contact.Height = 330;
            contact.ShowDistances = true;
            _trace = trace;

            var dssp = Dssp.Assign(_structure);
            string ssString = dssp.ToSecondaryStructureString();
            string residues = string.Concat(dssp.Residues.Select(r => r.Residue.OneLetterCode));
            var strip = new SecondaryStructureStrip { Height = 46 };
            strip.Set(ssString, residues);
            var comp = dssp.Composition();
            int Count(params SecondaryStructure[] s) => s.Sum(x => comp.TryGetValue(x, out var n) ? n : 0);

            var traceImg = Ui.PlotImage(trace);
            traceImg.MaxWidth = 400;
            var contactImg = Ui.PlotImage(contact);
            contactImg.MaxWidth = 400;
            var images = Ui.Columns("*,*", 14,
                Ui.Column(6, Ui.Eyebrow("Cα trace · " + coloring), traceImg),
                Ui.Column(6, Ui.Eyebrow("Cα distance map"), contactImg));

            var ssInfo = Ui.Row(18,
                Ui.Legend(("helix H/G/I", Palette.DyeT), ("strand E/B", Palette.DyeG), ("turn/bend T/S", Palette.Ink3)),
                Ui.Text($"{Count(SecondaryStructure.AlphaHelix, SecondaryStructure.Helix310, SecondaryStructure.PiHelix)} helix · {Count(SecondaryStructure.Strand, SecondaryStructure.Bridge)} strand · {dssp.HydrogenBonds.Count} backbone H-bonds · {_structure.Disulfides.Count} disulfides", "unit"));

            _structureView.Content = Ui.Column(14, images,
                Ui.Column(8, Ui.Eyebrow("DSSP secondary structure"), strip, ssInfo,
                    Ui.Mono(ssString, 11, Palette.Ink2Brush)));
        }
        catch (Exception ex)
        {
            _structureView.Content = Ui.Notice("Structure view failed", ex.Message, Palette.DyeT);
        }
    }

    private static Color ToColor(Visual.Net.Core.Colors.Color c) => Color.Parse(c.ToHex());

    private string Snippet(string seq, int window) => $$"""
        using Biotech.Net.Proteomics;
        using Biotech.Net.StructuralBiology;
        using Biotech.Net.Visualization;

        // ProtParam
        string sequence = {{GalleryData.Literal(seq)}};
        var pa = new ProteinAnalysis(sequence);
        double mw = pa.MolecularWeight();          // Da
        double pI = pa.IsoelectricPoint();         // Bjellqvist pKa set
        double gravy = pa.Gravy();                 // Kyte–Doolittle mean
        double instability = pa.InstabilityIndex();
        var eps = pa.ExtinctionCoefficient();      // (Reduced, Cystines)
        Console.WriteLine(pa.Report());

        // Hydropathy profile and transmembrane segments
        double[] kd = pa.Profile(AminoAcidScale.KyteDoolittle, window: {{window}});
        var tm = TransmembranePredictor.Predict(sequence, window: 19, threshold: 1.6);

        // Structure: crambin is embedded in the library (no download)
        Structure crambin = SampleStructures.Crambin();
        DsspResult dssp = Dssp.Assign(crambin);
        Console.WriteLine(dssp.ToSecondaryStructureString());

        var trace = ProteinViewer.CaTrace(crambin);
        trace.Coloring = ResidueColoring.{{(_coloring?.SelectedItem ?? ResidueColoring.SecondaryStructure)}};
        trace.Save("1crn-trace.svg");
        ProteinViewer.ContactMap(crambin, cutoff: 8.0).Save("1crn-contacts.png");
        File.WriteAllText("1crn-3d.html", ProteinViewer.To3DHtml(crambin));   // interactive 3Dmol.js
        """;
}
