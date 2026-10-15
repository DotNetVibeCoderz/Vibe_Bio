using System.Diagnostics;
using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Media;
using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Visualization;

namespace Biotech.Net.Gallery.Pages;

public sealed class DrugDiscoveryPage : GalleryPage
{
    public override string Title => "Drug Discovery Lab";
    public override string Modules => "BIOTECH.NET.CHEMISTRY · BIOTECH.NET.DRUGDISCOVERY · VISUALIZATION";
    public override string Lede => "Type a SMILES string to depict the molecule, compute RDKit-compatible descriptors, check drug-likeness and ADMET flags, and rank a built-in library by fingerprint similarity.";
    public override char Channel => 'G';

    /// <summary>The built-in screening library: approved drugs across a few target classes.</summary>
    private static readonly (string Name, string Smiles)[] Library =
    [
        ("Imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
        ("Nilotinib", "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1"),
        ("Gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
        ("Erlotinib", "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"),
        ("Dasatinib", "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1"),
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
        ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
        ("Diclofenac", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
        ("Celecoxib", "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"),
        ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
        ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
        ("Sildenafil", "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"),
        ("Metformin", "CN(C)C(=N)N=C(N)N"),
        ("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"),
        ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
        ("Nicotine", "CN1CCCC1c1cccnc1"),
        ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    ];

    private TextBox _smiles = null!;
    private readonly ContentControl _left = new();
    private readonly ContentControl _right = new();
    private MoleculeViewer? _viewer;
    private string _resultText = "";
    private IReadOnlyList<LibraryCompound>? _library;

    public override BioPlot? ExportPlot => _viewer;
    public override string ResultText => _resultText;

    protected override Control Build()
    {
        _smiles = Ui.Input(Library[0].Smiles, "seq");
        _smiles.AcceptsReturn = false;
        var chips = new WrapPanel { Orientation = Orientation.Horizontal };
        foreach (var (name, smi) in Library.Take(10))
        {
            var b = Ui.Button(name, () => { _smiles.Text = smi; Analyze(); }, "chip");
            b.Margin = new Thickness(0, 0, 6, 6);
            chips.Children.Add(b);
        }
        var input = Ui.Panel("SMILES", Ui.Column(10, _smiles, chips, Ui.Row(8, Ui.Button("Analyze", Analyze, "primary"))));
        var body = Ui.Columns("420,*", 18, Ui.Scroll(Ui.Column(14, input, _left)), Ui.Scroll(_right));
        Analyze();
        return body;
    }

    private void Analyze()
    {
        string smi = (_smiles.Text ?? "").Trim();
        string name = Library.FirstOrDefault(x => x.Smiles == smi).Name ?? "";
        var mol = SmilesParser.TryParse(smi, name);
        if (mol is null)
        {
            Status = $"\"{smi}\" is not valid SMILES. Check ring closures and brackets, or pick a preset.";
            _right.Content = Ui.Notice("SMILES could not be parsed", "Check ring-closure digits, brackets and aromatic atoms (lower-case c, n, o, s), or start from one of the presets.", Palette.DyeT);
            return;
        }
        try
        {
            var sw = Stopwatch.StartNew();
            var d = MolecularDescriptors.Compute(mol);
            var admet = Admet.Profile(mol);
            var alerts = StructuralAlerts.Match(mol);
            _library ??= Library.Select(x => LibraryCompound.FromSmiles(x.Name, x.Smiles)).ToList();
            var hits = SimilarityScreen.Run(mol, _library, new SimilarityScreenOptions
            {
                Fingerprints = [FingerprintType.Ecfp4, FingerprintType.Maccs],
                TopK = 8,
            });
            var viewer = new MoleculeViewer(mol, 420, 330) { Theme = Palette.PlotTheme, ShowFormula = true };
            if (alerts.Count > 0)
            {
                var atoms = alerts.SelectMany(a => a.Pattern.Matches(mol).SelectMany(m => m)).ToHashSet();
                viewer.HighlightAtoms = atoms;
            }
            _viewer = viewer;
            sw.Stop();

            var img = Ui.PlotImage(viewer);
            img.MaxWidth = 420;
            _left.Content = Ui.Column(14,
                Ui.Panel("2D depiction", Ui.Column(8, new Border { Child = img, Classes = { "well" }, Padding = new Thickness(0) },
                    Ui.Mono(mol.ToSmiles(), 11, Palette.Ink3Brush)), d.Formula),
                Ui.Panel("Structural alerts · PAINS + Brenk", alerts.Count == 0
                    ? Ui.Text("No PAINS or Brenk alerts matched.", "caption")
                    : Ui.Column(4, alerts.Select(a => (Control)Ui.Text($"{a.Family} · {a.Name}", "td")).ToArray()),
                    alerts.Count == 0 ? "clean" : $"{alerts.Count} matched"));

            var readouts = Ui.ReadoutRow(
                Ui.Readout("MW", d.MolecularWeight.ToString("0.0", CultureInfo.InvariantCulture), "Da"),
                Ui.Readout("cLogP", d.LogP.ToString("0.00", CultureInfo.InvariantCulture), "", Palette.DyeG),
                Ui.Readout("TPSA", d.Tpsa.ToString("0.0", CultureInfo.InvariantCulture), "Å²"),
                Ui.Readout("HBD / HBA", $"{d.LipinskiHbd}/{d.LipinskiHba}", ""),
                Ui.Readout("RotB", d.RotatableBonds.ToString(CultureInfo.InvariantCulture), ""),
                Ui.Readout("QED", admet.Qed.ToString("0.00", CultureInfo.InvariantCulture), "", Palette.DyeA));

            var rules = new WrapPanel();
            foreach (var r in admet.Rules)
            {
                var col = r.Passed ? Palette.DyeA : Palette.DyeT;
                var tag = new Border
                {
                    Classes = { "tag" },
                    BorderBrush = new SolidColorBrush(col),
                    Background = new SolidColorBrush(col, 0.08),
                    Margin = new Thickness(0, 0, 8, 8),
                    Child = new TextBlock
                    {
                        Text = (r.Passed ? "PASS  " : "FAIL  ") + r.Rule + (r.Passed ? "" : $" · {string.Join(", ", r.Violations)}"),
                        FontFamily = Palette.Mono, FontSize = 11, Foreground = new SolidColorBrush(col),
                    },
                };
                rules.Children.Add(tag);
            }

            var admetKv = Ui.KeyValues(
            [
                ("ESOL logS", $"{admet.EsolLogS:0.00}", null),
                ("Solubility", $"{admet.EsolSolubilityMgPerL:0.##} mg/L · {admet.SolubilityClass}", null),
                ("GI absorption (BOILED-Egg)", admet.GiAbsorption ? "high" : "low", admet.GiAbsorption ? Palette.DyeA : Palette.DyeT),
                ("Blood–brain barrier", admet.BloodBrainBarrier ? "permeant" : "not permeant", admet.BloodBrainBarrier ? Palette.DyeC : null),
                ("hERG risk", $"{admet.HergRisk:0.00}", admet.HergRisk >= 0.5 ? Palette.DyeT : null),
                ("P-gp substrate (likely)", admet.LikelyPgpSubstrate ? "yes" : "no", null),
                ("CYP inhibition flags", string.Join(" ", admet.CypFlags.Where(kv => kv.Value).Select(kv => kv.Key)) is { Length: > 0 } cyp ? cyp : "none", null),
                ("Synthetic complexity", $"{admet.Complexity:0.00}", null),
                ("Fraction Csp3 · rings", $"{d.FractionCsp3:0.00} · {d.RingCount} ({d.AromaticRingCount} aromatic)", null),
            ]);

            var simTable = Ui.Table(["#", "Compound", "Score", "ECFP4", "MACCS"], ["Auto", "Auto", "*", "Auto", "Auto"],
                hits.Select(h => new Ui.Cell[]
                {
                    new(h.Rank.ToString(CultureInfo.InvariantCulture), Palette.Ink3),
                    new(h.Id, h.Score > 0.99 ? Palette.DyeG : null),
                    new(h.Score.ToString("0.00", CultureInfo.InvariantCulture), Palette.DyeG, Bar: true, BarValue: h.Score),
                    new(h.PerFingerprint.TryGetValue(FingerprintType.Ecfp4, out var e) ? e.ToString("0.00", CultureInfo.InvariantCulture) : "—", Palette.Ink2),
                    new(h.PerFingerprint.TryGetValue(FingerprintType.Maccs, out var m) ? m.ToString("0.00", CultureInfo.InvariantCulture) : "—", Palette.Ink2),
                }));

            _right.Content = Ui.Column(14,
                Ui.Panel("Descriptors", readouts, $"computed in {GalleryData.Ms(sw)}"),
                Ui.Panel("Drug-likeness rules", rules, $"{admet.Rules.Count(r => r.Passed)}/{admet.Rules.Count} pass"),
                Ui.Columns("*,*", 14,
                    Ui.Panel("ADMET estimates", admetKv),
                    Ui.Panel("Similarity search", simTable, $"Tanimoto · {Library.Length} compounds")));

            TraceSequence = GalleryData.TextToBases(smi);
            TraceLabel = "TRACE · SMILES encoded 2 bases/char";
            Status = $"{d.Formula} · MW {d.MolecularWeight:0.00} · cLogP {d.LogP:0.00} · TPSA {d.Tpsa:0.0} · QED {admet.Qed:0.00} · nearest {hits.FirstOrDefault()?.Id} · {GalleryData.Ms(sw)}";
            _resultText = admet.ToString() + "\n" + string.Join('\n', hits.Select(h => $"{h.Rank}\t{h.Id}\t{h.Score:0.000}"));
            Code = Snippet(smi);
        }
        catch (Exception ex)
        {
            Status = "Analysis failed: " + ex.Message;
            _right.Content = Ui.Notice("Analysis failed", ex.Message, Palette.DyeT);
        }
    }

    private static string Snippet(string smiles) => $$"""
        using Biotech.Net.Chemistry;
        using Biotech.Net.DrugDiscovery;
        using Biotech.Net.Visualization;

        var mol = Molecule.FromSmiles("{{smiles}}");

        // RDKit-compatible descriptors and drug-likeness
        DescriptorSet d = MolecularDescriptors.Compute(mol);
        Console.WriteLine($"{d.Formula} MW={d.MolecularWeight:0.00} cLogP={d.LogP:0.00} TPSA={d.Tpsa:0.0}");
        Console.WriteLine(DrugLikeness.Lipinski(d));

        // ADMET profile (ESOL, BOILED-Egg, hERG, CYP flags, QED) and structural alerts
        AdmetProfile admet = Admet.Profile(mol);
        var alerts = StructuralAlerts.Match(mol);   // PAINS + Brenk

        // Similarity screen against a library (ECFP4 + MACCS, mean fusion)
        var library = new[]
        {
            LibraryCompound.FromSmiles("Imatinib", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
            LibraryCompound.FromSmiles("Gefitinib", "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
            // … {{Library.Length}} compounds in the Gallery
        };
        var hits = SimilarityScreen.Run(mol, library, new SimilarityScreenOptions
        {
            Fingerprints = [FingerprintType.Ecfp4, FingerprintType.Maccs],
            TopK = 8,
        });

        // 2D depiction
        new MoleculeViewer(mol, 420, 330) { Theme = BioTheme.Dark }.Save("molecule.svg");
        """;
}
