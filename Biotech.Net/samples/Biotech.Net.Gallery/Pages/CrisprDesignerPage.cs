using System.Diagnostics;
using System.Globalization;
using System.Text;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Media;
using Biotech.Net.Core;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Sequences;
using Biotech.Net.SyntheticBiology.Crispr;

namespace Biotech.Net.Gallery.Pages;

public sealed class CrisprDesignerPage : GalleryPage
{
    public override string Title => "CRISPR Designer";
    public override string Modules => "BIOTECH.NET.SYNTHETICBIOLOGY · CRISPR";
    public override string Lede => "Paste a target, pick a nuclease and PAM, and get guides ranked by Doench 2014 on-target activity and CFD / MIT specificity against a background genome.";
    public override char Channel => 'A';

    private TextBox _target = null!;
    private ComboBox _nuclease = null!;
    private CheckBox _offTarget = null!;
    private readonly ContentControl _results = new();
    private string _resultText = "";
    private Sequence? _background;

    public override string ResultText => _resultText;

    protected override Control Build()
    {
        _target = Ui.Input(GalleryData.EgfpFasta, "seq", 200);
        _nuclease = new ComboBox
        {
            ItemsSource = CasNuclease.All.Select(n => n.ToString()).ToList(),
            SelectedIndex = 0,
            HorizontalAlignment = HorizontalAlignment.Stretch,
        };
        _nuclease.SelectionChanged += (_, _) => Design();
        _offTarget = new CheckBox { Content = "Profile off-targets (target + 60 kb background)", IsChecked = true, FontSize = 12 };
        _offTarget.IsCheckedChanged += (_, _) => Design();

        var left = Ui.Column(14,
            Ui.Panel("Target", Ui.Column(10, _target,
                Ui.Row(6,
                    Ui.Button("EGFP", () => { _target.Text = GalleryData.EgfpFasta; Design(); }, "chip"),
                    Ui.Button("Random 1 kb", () => { _target.Text = ">random_1kb\n" + GalleryData.Wrap(RandomSequence.Dna(1000, 3).ToString(), 60); Design(); }, "chip")))),
            Ui.Panel("Nuclease · PAM", Ui.Column(10,
                _nuclease,
                _offTarget,
                Ui.Text("Composite = 0.50 × on-target·penalties + 0.35 × CFD + 0.15 × MIT specificity, renormalised over the terms available.", "caption"),
                Ui.Button("Design guides", Design, "primary"))));

        var root = Ui.Columns("340,*", 18, Ui.Scroll(left), Ui.Scroll(_results));
        Design();
        return root;
    }

    private Sequence Background()
    {
        if (_background is not null) return _background;
        // 60 kb synthetic background with a few near-copies of EGFP stretches, so specificity is not trivially 100.
        var sb = new StringBuilder(RandomSequence.Dna(60_000, 11, 0.5).ToString());
        var rng = new Random(5);
        var egfpGuides = GuideFinder.Find(Sequence.Dna(GalleryData.Egfp, "egfp"));
        for (int g = 0; g < egfpGuides.Count; g += 2)
        {
            var piece = egfpGuides[g].FullSite.ToCharArray();
            int mismatches = rng.Next(1, 4);
            for (int m = 0; m < mismatches; m++)
            {
                int p = rng.Next(0, 10);
                piece[p] = piece[p] == 'A' ? 'C' : 'A';
            }
            int at = rng.Next(1000, 59_000);
            sb.Remove(at, piece.Length).Insert(at, piece);
        }
        return _background = Sequence.Dna(sb.ToString(), "background");
    }

    private void Design()
    {
        try
        {
            var sw = Stopwatch.StartNew();
            string text = _target.Text ?? "";
            string raw = text.TrimStart().StartsWith('>')
                ? string.Concat(text.Split('\n').Where(l => !l.StartsWith('>')))
                : text;
            string clean = new string(raw.Where(char.IsLetter).ToArray()).ToUpperInvariant();
            if (clean.Length < 30) { Status = "Paste a target of at least 30 bases."; return; }
            var target = Sequence.Dna(clean, "target");
            var nuclease = CasNuclease.All[Math.Max(0, _nuclease.SelectedIndex)];
            bool offTargets = _offTarget.IsChecked == true;
            CrisprSiteIndex? index = offTargets ? CrisprSiteIndex.Build([target, Background()], nuclease) : null;
            var designs = GuideDesigner.Design(target, new GuideDesignOptions
            {
                Search = new GuideSearchOptions { Nuclease = nuclease },
                OffTargetIndex = index,
                MaxResults = 15,
            });
            sw.Stop();

            if (designs.Count == 0)
            {
                _results.Content = Ui.Notice($"No {nuclease.Name} sites passed the filters",
                    $"The {nuclease.Pam} PAM does not occur on either strand, or every spacer failed the GC / poly-T / homopolymer filters. Try a relaxed PAM such as SpG (NGN) or SpRY (NNN).", Palette.DyeG);
                Status = $"0 guides for {nuclease}";
                return;
            }

            var best = designs[0];
            var readouts = Ui.ReadoutRow(
                Ui.Readout("Guides ranked", designs.Count.ToString(CultureInfo.InvariantCulture), ""),
                Ui.Readout("Best score", best.CompositeScore.ToString("0.000", CultureInfo.InvariantCulture), "", Palette.DyeA),
                Ui.Readout("Best on-target", double.IsNaN(best.OnTargetScore) ? "n/a" : best.OnTargetScore.ToString("0.000", CultureInfo.InvariantCulture), "RS1", Palette.DyeG),
                Ui.Readout("Best CFD spec.", double.IsNaN(best.CfdSpecificity) ? "—" : best.CfdSpecificity.ToString("0.0", CultureInfo.InvariantCulture), "/100", Palette.DyeC),
                Ui.Readout("Target", target.Length.ToString("N0", CultureInfo.InvariantCulture), "bp"));

            var map = new SequenceMap { Height = 82, RowLabels = ["+ strand", "− strand"] };
            map.Set(target.Length, designs.Select((d, i) => new SequenceMap.Mark(d.Guide.Start, d.Guide.End, $"#{i + 1}",
                i == 0 ? Palette.DyeG : d.Guide.Strand == Strand.Reverse ? Palette.DyeC : Palette.DyeA,
                d.Guide.Strand == Strand.Reverse ? 1 : 0, Reverse: d.Guide.Strand == Strand.Reverse)).ToList());

            var table = Ui.Table(
                ["#", "Spacer 5'→3'", "PAM", "Str", "Pos", "Cut", "GC", "On-target", "CFD", "MIT", "Off", "Score"],
                ["Auto", "Auto", "Auto", "Auto", "Auto", "Auto", "Auto", "120", "Auto", "Auto", "Auto", "130"],
                designs.Select((d, i) => new Ui.Cell[]
                {
                    new((i + 1).ToString(CultureInfo.InvariantCulture), i == 0 ? Palette.DyeG : Palette.Ink3),
                    new(d.Guide.Spacer, i == 0 ? Palette.DyeG : null),
                    new(d.Guide.Pam, Palette.DyeC),
                    new(d.Guide.Strand == Strand.Reverse ? "−" : "+"),
                    new((d.Guide.Start + 1).ToString("N0", CultureInfo.InvariantCulture)),
                    new(d.Guide.CutSite.ToString("N0", CultureInfo.InvariantCulture), Palette.Ink2),
                    new($"{d.Guide.GcContent * 100:0}%"),
                    new(double.IsNaN(d.OnTargetScore) ? "n/a" : d.OnTargetScore.ToString("0.00", CultureInfo.InvariantCulture), Palette.DyeG, Bar: !double.IsNaN(d.OnTargetScore), BarValue: d.OnTargetScore),
                    new(double.IsNaN(d.CfdSpecificity) ? "—" : d.CfdSpecificity.ToString("0", CultureInfo.InvariantCulture), Palette.DyeC),
                    new(double.IsNaN(d.MitSpecificity) ? "—" : d.MitSpecificity.ToString("0", CultureInfo.InvariantCulture)),
                    new(d.OffTargets is null ? "—" : d.OffTargetCount.ToString(CultureInfo.InvariantCulture), d.OffTargetCount > 0 ? Palette.DyeT : null),
                    new(d.CompositeScore.ToString("0.000", CultureInfo.InvariantCulture), Palette.DyeA, Bar: true, BarValue: d.CompositeScore),
                }));

            Control offPanel;
            if (best.OffTargets is { } prof && prof.OffTargets.Any())
            {
                offPanel = Ui.Table(["Site", "Mismatches", "Alignment", "CFD"], ["Auto", "Auto", "*", "Auto"],
                    prof.OffTargets.Take(5).Select(o => new Ui.Cell[]
                    {
                        new($"{o.SequenceId}:{o.Start + 1:N0}"),
                        new(o.MismatchCount.ToString(CultureInfo.InvariantCulture), o.MismatchCount <= 2 ? Palette.DyeT : null),
                        new(o.Alignment(best.Guide.Spacer), Palette.Ink2),
                        new(o.CfdScore.ToString("0.000", CultureInfo.InvariantCulture)),
                    }));
            }
            else offPanel = Ui.Text(offTargets ? "No off-target sites within 4 mismatches for the top guide." : "Off-target profiling is off.", "caption");

            _results.Content = Ui.Column(14,
                Ui.Panel($"{nuclease.Name} · {nuclease.Pam} · guide {nuclease.GuideLength} nt", readouts, $"designed in {GalleryData.Ms(sw)}"),
                Ui.Panel("Guide positions", map, "top guide in amber"),
                Ui.Panel("Ranked guides", table, "Doench 2014 · CFD · MIT"),
                Ui.Panel($"Top guide off-targets · {best.Guide.Spacer}", offPanel, best.Guide.U6CloningOligo == best.Guide.Spacer ? "U6-ready" : "U6 oligo: " + best.Guide.U6CloningOligo));

            TraceSequence = best.Guide.FullSite + clean[..Math.Min(clean.Length, 120)];
            TraceLabel = $"TRACE · top guide {best.Guide.Spacer}{best.Guide.Pam}";
            Status = $"{designs.Count} guides · {nuclease} · best {best.Guide.Spacer} {best.Guide.Pam} score {best.CompositeScore:0.000}" +
                     (index is null ? "" : $" · index {index.SiteCount:N0} sites") + $" · {GalleryData.Ms(sw)}";
            _resultText = string.Join('\n', designs.Select((d, i) => $"{i + 1}\t{d}"));
            Code = Snippet(clean, nuclease, offTargets);
        }
        catch (Exception ex)
        {
            Status = "Guide design failed: " + ex.Message;
            _results.Content = Ui.Notice("Guide design failed", ex.Message, Palette.DyeT);
        }
    }

    private static string NucleusField(CasNuclease n) => typeof(CasNuclease).GetFields()
        .FirstOrDefault(f => f.IsStatic && f.FieldType == typeof(CasNuclease) && ReferenceEquals(f.GetValue(null), n))?.Name ?? "SpCas9";

    private static string Snippet(string target, CasNuclease nuclease, bool offTargets) => $$"""
        using Biotech.Net.Core;
        using Biotech.Net.Sequences;
        using Biotech.Net.SyntheticBiology.Crispr;

        var target = Sequence.Dna({{GalleryData.Literal(target)}});
        var nuclease = CasNuclease.{{NucleusField(nuclease)}};        // {{nuclease}}

        {{(offTargets ? """
        // Off-target index over the genome you care about (here: target + background)
        var background = RandomSequence.Dna(60_000, seed: 11);
        var index = CrisprSiteIndex.Build([target, background], nuclease);
        """ : "CrisprSiteIndex? index = null;   // no off-target profiling")}}

        var designs = GuideDesigner.Design(target, new GuideDesignOptions
        {
            Search = new GuideSearchOptions { Nuclease = nuclease },
            OffTargetIndex = index,
            MaxResults = 15,
        });

        foreach (var d in designs)
            Console.WriteLine($"{d.Guide.Spacer} {d.Guide.Pam} {d.Guide.Strand} cut@{d.Guide.CutSite} " +
                              $"on={d.OnTargetScore:0.00} CFD={d.CfdSpecificity:0} MIT={d.MitSpecificity:0} score={d.CompositeScore:0.000}");
        """;
}
