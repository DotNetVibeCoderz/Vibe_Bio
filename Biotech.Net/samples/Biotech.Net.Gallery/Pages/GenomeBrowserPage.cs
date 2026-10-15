using System.Diagnostics;
using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Biotech.Net.Core;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Genomics.Annotation;
using Biotech.Net.Genomics.Vcf;
using Biotech.Net.Sequences;
using Biotech.Net.Sequences.IO;
using Biotech.Net.Visualization;

namespace Biotech.Net.Gallery.Pages;

public sealed class GenomeBrowserPage : GalleryPage
{
    public override string Title => "Genome Browser";
    public override string Modules => "BIOTECH.NET.VISUALIZATION · BIOTECH.NET.GENOMICS";
    public override string Lede => "A bundled 12 kb region with GFF3 gene models and VCF calls. Pan, zoom and toggle tracks; every variant is annotated with its predicted consequence.";
    public override char Channel => 'C';

    private Genome _genome = null!;
    private Sequence _chr = null!;
    private List<Variant> _variants = null!;
    private VariantEffectPredictor _vep = null!;
    private double[] _gc = null!;
    private int _gcStart;
    private long _start, _end;
    private GenomeBrowser? _plot;

    private readonly Image _image = new() { Stretch = Avalonia.Media.Stretch.Uniform, HorizontalAlignment = HorizontalAlignment.Left };
    private readonly TextBlock _regionLabel = Ui.Text("", "td");
    private readonly ContentControl _table = new();
    private CheckBox _ruler = null!, _genes = null!, _gcTrack = null!, _varTrack = null!, _seqTrack = null!;

    public override BioPlot? ExportPlot => _plot;

    public override string ResultText => _plot is null ? "" :
        $"{Region}\n" + string.Join('\n', _variants.Where(v => v.Position - 1 >= _start && v.Position - 1 < _end).Select(v => v.ToString()));

    private GenomicInterval Region => new(_chr.Id, _start, _end);

    protected override Control Build()
    {
        _chr = FastaReader.Parse(GalleryData.Resource("region.fa"))[0];
        _genome = new Genome("Gallery demo", "synthetic", "demo-1").Add(Chromosome.FromSequence(_chr));
        using (var gff = new StringReader(GalleryData.Resource("region.gff3")))
            foreach (var g in GffReader.BuildGenes(GffReader.Read(gff))) _genome.Add(g);
        _variants = VcfReader.Parse(GalleryData.Resource("region.vcf")).Records.Select(r => r.Variant).ToList();
        _vep = new VariantEffectPredictor(_genome);
        var gcWin = SlidingWindow.GcContent(_chr, 50, 1);
        _gc = gcWin.Select(w => w.Value * 100).ToArray();
        _gcStart = gcWin.Count > 0 ? gcWin[0].Start + 25 : 0;
        _start = 0; _end = _chr.Length;

        CheckBox Toggle(string label) { var c = new CheckBox { Content = label, IsChecked = true, FontSize = 12 }; c.IsCheckedChanged += (_, _) => Render(); return c; }
        _ruler = Toggle("Ruler"); _genes = Toggle("Genes"); _gcTrack = Toggle("GC %"); _varTrack = Toggle("Variants"); _seqTrack = Toggle("Sequence");

        var nav = Ui.Row(6,
            Ui.Button("◀ Pan", () => Pan(-0.33)),
            Ui.Button("Pan ▶", () => Pan(0.33)),
            Ui.Button("Zoom in", () => Zoom(0.5)),
            Ui.Button("Zoom out", () => Zoom(2)),
            Ui.Button("Whole region", () => { _start = 0; _end = _chr.Length; Render(); }));
        var geneChips = Ui.Row(6, _genome.Genes.Select(g => (Control)Ui.Button(g.Name, () => Focus(g), "chip")).ToArray());

        var toolbar = new Grid { ColumnDefinitions = new ColumnDefinitions("Auto,*,Auto") };
        var regionBox = Ui.Row(10, Ui.Eyebrow("Region"), _regionLabel);
        regionBox.VerticalAlignment = VerticalAlignment.Center;
        toolbar.Children.Add(regionBox);
        var right = Ui.Row(14, Ui.Row(6, Ui.Eyebrow("Go to"), geneChips), nav);
        right.VerticalAlignment = VerticalAlignment.Center;
        Grid.SetColumn(right, 2);
        toolbar.Children.Add(right);

        var tracks = Ui.Row(18, Ui.Eyebrow("Tracks"), _ruler, _genes, _gcTrack, _varTrack, _seqTrack);
        tracks.VerticalAlignment = VerticalAlignment.Center;

        var browserPanel = Ui.Panel("Tracks", Ui.Column(10, toolbar, tracks, new Border { Child = _image, Classes = { "well" }, Padding = new Thickness(0) }),
            "GFF3 · VCF · FASTA");
        var body = Ui.Column(14, browserPanel, Ui.Panel("Variants in view · predicted effect", _table, "VariantEffectPredictor · Sequence Ontology terms"));
        Render();
        return Ui.Scroll(body);
    }

    private void Focus(Gene g)
    {
        long pad = g.Location.Length / 8;
        _start = Math.Max(0, g.Location.Start - pad);
        _end = Math.Min(_chr.Length, g.Location.End + pad);
        Render();
    }

    private void Pan(double fraction)
    {
        long span = _end - _start, shift = (long)(span * fraction);
        _start = Math.Clamp(_start + shift, 0, _chr.Length - span);
        _end = _start + span;
        Render();
    }

    private void Zoom(double factor)
    {
        long span = _end - _start, mid = (_start + _end) / 2;
        long ns = Math.Clamp((long)(span * factor), 60, _chr.Length);
        _start = Math.Clamp(mid - ns / 2, 0, _chr.Length - ns);
        _end = _start + ns;
        Render();
    }

    private void Render()
    {
        try
        {
            var sw = Stopwatch.StartNew();
            var region = Region;
            var view = new GenomeBrowser(region, width: 1080) { Theme = Palette.PlotTheme, NameGutter = 92 };
            if (_ruler.IsChecked == true) view.AddRuler();
            if (_genes.IsChecked == true) view.AddGenes(_genome.GenesIn(region), "Genes");
            if (_gcTrack.IsChecked == true) view.AddCoverage(_gc, _gcStart, "GC % (50 bp)", Visual.Net.Core.Colors.Color.Parse("#F0B03F"));
            var inView = _variants.Where(v => region.Contains(v.Position - 1)).ToList();
            if (_varTrack.IsChecked == true) view.AddVariants(inView, "Variants");
            if (_seqTrack.IsChecked == true) view.AddSequence(_chr, 0, "Sequence");
            _plot = view;
            _image.Source = Ui.PlotBitmap(view);
            _image.MaxWidth = view.Width;
            sw.Stop();

            _regionLabel.Text = $"{_chr.Id}:{_start + 1:N0}-{_end:N0}  ({_end - _start:N0} bp)";
            var rows = inView.Select(v =>
            {
                var effects = _vep.Predict(v);
                var e = effects.OrderBy(x => x.Impact).FirstOrDefault();
                var impactColor = e?.Impact switch
                {
                    VariantImpact.High => Palette.DyeT,
                    VariantImpact.Moderate => Palette.DyeG,
                    VariantImpact.Low => Palette.DyeA,
                    _ => Palette.Ink3,
                };
                double af = v.Info.TryGetValue("AF", out var s) && double.TryParse(s, CultureInfo.InvariantCulture, out var f) ? f : double.NaN;
                return new Ui.Cell[]
                {
                    new(v.Position.ToString("N0", CultureInfo.InvariantCulture)),
                    new(v.Id ?? "."),
                    new($"{Short(v.Reference)}>{Short(string.Join(',', v.Alternates))}", Palette.DyeC),
                    new(v.Type.ToString()),
                    new(double.IsNaN(af) ? "." : af.ToString("0.00", CultureInfo.InvariantCulture), Palette.DyeG, Bar: true, BarValue: double.IsNaN(af) ? 0 : af),
                    new(e?.GeneName ?? "—"),
                    new(e?.ConsequenceTerms ?? "intergenic_variant"),
                    new(e?.HgvsP ?? e?.HgvsC ?? "", Palette.Ink2),
                    new(e?.Impact.ToString().ToUpperInvariant() ?? "MODIFIER", impactColor),
                };
            });
            _table.Content = Ui.Table(["Pos", "ID", "Ref>Alt", "Type", "AF", "Gene", "Consequence", "HGVS", "Impact"],
                ["Auto", "Auto", "Auto", "Auto", "120", "Auto", "2*", "*", "Auto"], rows);

            TraceSequence = _chr.Subsequence((int)_start, (int)Math.Min(400, _end - _start)).ToString();
            TraceLabel = $"TRACE · {_chr.Id}:{_start + 1:N0}";
            Status = $"{region} · {view.Tracks.Count} tracks · {inView.Count} variants · {_genome.Genes.Count} genes in region file · rendered in {GalleryData.Ms(sw)}";
            Code = Snippet(region);
        }
        catch (Exception ex)
        {
            Status = "Could not render the region: " + ex.Message;
        }
    }

    private static string Short(string s) => s.Length <= 6 ? s : s[..5] + "…";

    private string Snippet(GenomicInterval region) => $$"""
        using Biotech.Net.Core;
        using Biotech.Net.Genomics.Annotation;
        using Biotech.Net.Genomics.Vcf;
        using Biotech.Net.Sequences;
        using Biotech.Net.Sequences.IO;
        using Biotech.Net.Visualization;

        // Reference, gene models (GFF3) and calls (VCF)
        var chr = FastaReader.Read("region.fa").First();
        var genome = new Genome("Gallery demo").Add(Chromosome.FromSequence(chr));
        GffReader.LoadGenes(genome, "region.gff3");
        var variants = VcfReader.ReadRecords("region.vcf").Select(r => r.Variant).ToList();

        // Region in view: 0-based half-open ({{region}})
        var region = new GenomicInterval("{{region.Chromosome}}", {{region.Start}}, {{region.End}});
        var gc = SlidingWindow.GcContent(chr, window: 50, step: 1);

        var view = new GenomeBrowser(region, width: 1080) { Theme = BioTheme.Dark }
        {{(_ruler.IsChecked == true ? "    .AddRuler()\n" : "")}}{{(_genes.IsChecked == true ? "    .AddGenes(genome.GenesIn(region))\n" : "")}}{{(_gcTrack.IsChecked == true ? "    .AddCoverage(gc.Select(w => w.Value * 100).ToArray(), start: 25, name: \"GC % (50 bp)\")\n" : "")}}{{(_varTrack.IsChecked == true ? "    .AddVariants(variants.Where(v => region.Contains(v.Position - 1)))\n" : "")}}{{(_seqTrack.IsChecked == true ? "    .AddSequence(chr)\n" : "")}}    ;
        view.Save("region.svg");   // or .png / .html / .pdf

        // Predicted consequences (Sequence Ontology)
        var vep = new VariantEffectPredictor(genome);
        foreach (var v in variants)
            foreach (var effect in vep.Predict(v))
                Console.WriteLine($"{v} {effect.ConsequenceTerms} {effect.HgvsP} {effect.Impact}");
        """;
}
