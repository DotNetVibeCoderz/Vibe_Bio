using System.Diagnostics;
using System.Globalization;
using System.Text;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Platform.Storage;
using Biotech.Net.Core;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Sequences;
using Biotech.Net.Sequences.IO;

namespace Biotech.Net.Gallery.Pages;

public sealed class DnaExplorerPage : GalleryPage
{
    public override string Title => "DNA Explorer";
    public override string Modules => "BIOTECH.NET.SEQUENCES · BIOTECH.NET.CORE";
    public override string Lede => "Paste or load FASTA to read its composition, translate all six frames, find open reading frames, map restriction sites and estimate primer melting temperature.";
    public override char Channel => 'A';

    private TextBox _input = null!;
    private TextBox _oligo = null!;
    private readonly ContentControl _results = new();
    private string _resultText = "";

    public override string ResultText => _resultText;

    protected override Control Build()
    {
        _input = Ui.Input(GalleryData.EgfpFasta, "seq", 230);
        _oligo = Ui.Input(GalleryData.Egfp[..24], "seq");
        _oligo.AcceptsReturn = false;
        _oligo.TextWrapping = Avalonia.Media.TextWrapping.NoWrap;

        var samples = Ui.Row(6,
            Ui.Button("EGFP", () => { _input.Text = GalleryData.EgfpFasta; Analyze(); }, "chip"),
            Ui.Button("Random 2 kb", () => { _input.Text = ">random_2kb seed 7\n" + GalleryData.Wrap(RandomSequence.Dna(2000, 7, 0.52).ToString(), 60); Analyze(); }, "chip"),
            Ui.Button("Demo region", () => { _input.Text = GalleryData.Resource("region.fa"); Analyze(); }, "chip"));

        var actions = Ui.Row(8, Ui.Button("Analyze", Analyze, "primary"), Ui.Button("Load FASTA…", () => _ = LoadAsync()));

        var left = Ui.Column(14,
            Ui.Panel("Input · FASTA", Ui.Column(10, _input, samples, actions)),
            Ui.Panel("Primer Tm", Ui.Column(8,
                Ui.Labelled("Oligo 5'→3'", _oligo),
                Ui.Text("50 mM Na⁺, 25 nM oligo. Nearest-neighbour (SantaLucia 1998) with Wallace and salt-adjusted GC for comparison.", "caption"),
                Ui.Button("Compute Tm", Analyze))));

        var root = Ui.Columns("340,*", 18, Ui.Scroll(left), Ui.Scroll(_results));
        Analyze();
        return root;
    }

    private async Task LoadAsync()
    {
        var top = TopLevel.GetTopLevel(View);
        if (top is null) return;
        var files = await top.StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Open FASTA",
            AllowMultiple = false,
            FileTypeFilter = [new FilePickerFileType("FASTA") { Patterns = ["*.fa", "*.fasta", "*.fna", "*.ffn", "*.txt"] }],
        });
        if (files.Count == 0) return;
        await using var s = await files[0].OpenReadAsync();
        using var r = new StreamReader(s);
        _input.Text = await r.ReadToEndAsync();
        Analyze();
    }

    private void Analyze()
    {
        try
        {
            var sw = Stopwatch.StartNew();
            string text = _input.Text ?? "";
            var records = text.TrimStart().StartsWith('>') ? FastaReader.Parse(text) : [Sequence.Parse(new string(text.Where(char.IsLetter).ToArray()), "input")];
            if (records.Count == 0 || records[0].Length == 0) { Status = "Paste a DNA sequence or FASTA record to analyse."; return; }
            var seq = records[0];
            if (seq.Type != SequenceType.Dna) seq = Sequence.Dna(seq.ToString().ToUpperInvariant().Replace('U', 'T'), seq.Id, seq.Description);

            var counts = seq.Composition();
            double gc = seq.GcContent();
            double mw = NucleicAcidMass.MolecularWeight(seq, doubleStranded: true);
            var orfs = OrfFinder.Find(seq, new OrfOptions { MinLength = 150 });
            var frames = SixFrame.Translate(seq);
            var cutters = RestrictionEnzymes.All
                .Select(e => (Enzyme: e, Cuts: RestrictionDigest.FindCuts(seq, e)))
                .Where(x => x.Cuts.Count is > 0 and <= 3)
                .OrderBy(x => x.Cuts.Count).ThenBy(x => x.Enzyme.Name).Take(12).ToList();
            int window = Math.Clamp(seq.Length / 20, 30, 500);
            var gcWindow = SlidingWindow.GcContent(seq, window, Math.Max(1, window / 4));

            string oligoText = new string((_oligo.Text ?? "").Where(char.IsLetter).ToArray()).ToUpperInvariant();
            Sequence? oligo = oligoText.Length >= 8 ? Sequence.Dna(oligoText, "primer") : null;
            double tmNn = oligo is null ? double.NaN : MeltingTemperature.NearestNeighbor(oligo);
            double tmW = oligo is null ? double.NaN : MeltingTemperature.Wallace(oligo);
            double tmSalt = oligo is null ? double.NaN : MeltingTemperature.SaltAdjustedGc(oligo);
            double dg = oligo is null ? double.NaN : MeltingTemperature.FreeEnergy37(oligo);
            sw.Stop();

            // ---------------- view ----------------
            var readouts = Ui.ReadoutRow(
                Ui.Readout("Length", seq.Length.ToString("N0", CultureInfo.InvariantCulture), "bp"),
                Ui.Readout("GC content", (gc * 100).ToString("0.0", CultureInfo.InvariantCulture), "%", Palette.DyeG),
                Ui.Readout("dsDNA mass", (mw / 1000).ToString("0.0", CultureInfo.InvariantCulture), "kDa"),
                Ui.Readout("ORFs ≥ 50 aa", orfs.Count.ToString(CultureInfo.InvariantCulture), "", Palette.DyeA),
                Ui.Readout("Primer Tm", double.IsNaN(tmNn) ? "—" : tmNn.ToString("0.0", CultureInfo.InvariantCulture), "°C", Palette.DyeC));

            var comp = new BarChart
            {
                Height = 170,
                Bars = "ACGTN".Select(b => new BarChart.Bar(b.ToString(), counts.TryGetValue(b, out var c) ? c : 0, Palette.Channel(b),
                    $"{(counts.TryGetValue(b, out var c2) ? c2 : 0) * 100.0 / seq.Length:0.0}%")).ToList(),
            };
            var gcChart = new ProfileChart
            {
                Height = 170,
                Values = gcWindow.Select(v => v.Value * 100).ToArray(),
                Baseline = gc * 100,
                Positive = Palette.DyeG,
                Negative = Palette.DyeC,
                XLabel = "window",
                Format = "0",
                References = [(gc * 100, $"mean {gc * 100:0.0}%")],
            };

            var marks = new List<SequenceMap.Mark>();
            foreach (var o in orfs.Take(12))
                marks.Add(new(o.Start, o.End, $"{(o.Strand == Strand.Reverse ? "−" : "+")}{Math.Abs(o.Frame)} · {o.ProteinLength} aa", o.Strand == Strand.Reverse ? Palette.DyeC : Palette.DyeA,
                    o.Strand == Strand.Reverse ? 1 : 0, Reverse: o.Strand == Strand.Reverse));
            int row = 0;
            foreach (var (enz, cuts) in cutters.Take(8))
            {
                foreach (var c in cuts) marks.Add(new(c.Position, c.Position, enz.Name, Palette.DyeT, 2, Tick: true));
                row++;
            }
            var map = new SequenceMap { Height = 166, RowLabels = ["ORF +", "ORF −", "Cut sites"] };
            map.Set(seq.Length, marks);

            var frameText = new StringBuilder();
            foreach (var (frame, protein) in frames)
            {
                string p = protein.ToString();
                frameText.Append(frame > 0 ? $"+{frame} " : $"{frame} ").Append(p.Length > 96 ? p[..96] + "…" : p).Append('\n');
            }
            var framesBlock = Ui.Mono(frameText.ToString().TrimEnd(), 11.5);
            framesBlock.LineHeight = 20;

            var orfTable = Ui.Table(["Strand", "Frame", "Start", "End", "Length", "Protein"], ["Auto", "Auto", "Auto", "Auto", "Auto", "*"],
                orfs.Take(8).Select(o => new Ui.Cell[]
                {
                    new(o.Strand == Strand.Reverse ? "−" : "+", o.Strand == Strand.Reverse ? Palette.DyeC : Palette.DyeA),
                    new(o.Frame.ToString("+0;-0", CultureInfo.InvariantCulture)),
                    new((o.Start + 1).ToString("N0", CultureInfo.InvariantCulture)),
                    new(o.End.ToString("N0", CultureInfo.InvariantCulture)),
                    new($"{o.ProteinLength} aa"),
                    new(Preview(o.Translate(seq).ToString(), 28), Palette.Ink2),
                }));

            var enzTable = Ui.Table(["Enzyme", "Site", "Cuts", "Positions"], ["Auto", "Auto", "Auto", "*"],
                cutters.Select(x => new Ui.Cell[]
                {
                    new(x.Enzyme.Name, Palette.DyeT),
                    new(x.Enzyme.Site),
                    new(x.Cuts.Count.ToString(CultureInfo.InvariantCulture)),
                    new(string.Join(", ", x.Cuts.Select(c => c.Position.ToString("N0", CultureInfo.InvariantCulture))), Palette.Ink2),
                }));

            var tmPanel = oligo is null ? (Control)Ui.Text("Enter at least 8 bases.", "caption") : Ui.KeyValues(
            [
                ("Nearest-neighbour Tm", $"{tmNn:0.0} °C", Palette.DyeC),
                ("Salt-adjusted GC Tm", $"{tmSalt:0.0} °C", null),
                ("Wallace (2AT+4GC)", $"{tmW:0.0} °C", null),
                ("ΔG°37", $"{dg:0.00} kcal/mol", null),
                ("Length / GC", $"{oligo.Length} nt / {oligo.GcContent() * 100:0}%", null),
            ]);

            _results.Content = Ui.Column(14,
                Ui.Panel($"{seq.Id} · {Preview(seq.Description, 60)}", readouts, $"analysed in {GalleryData.Ms(sw)}"),
                Ui.Columns("*,*", 14,
                    Ui.Panel("Base composition", comp, "A C G T N"),
                    Ui.Panel($"GC sliding window", gcChart, $"window {window} bp")),
                Ui.Panel("Linear map", map, $"{orfs.Count} ORFs · {cutters.Count} enzymes cutting ≤3×"),
                Ui.Panel("Six-frame translation", framesBlock, "standard code · * = stop"),
                Ui.Columns("3*,2*", 14,
                    Ui.Panel("Open reading frames", orfTable, "ATG → stop, ≥ 150 nt"),
                    Ui.Column(14, Ui.Panel("Restriction sites", enzTable), Ui.Panel("Melting temperature", tmPanel))));

            TraceSequence = seq.ToString();
            TraceLabel = $"TRACE · {seq.Id} 1–{Math.Min(seq.Length, 140)}";
            Status = $"{seq.Id}: {seq.Length:N0} bp · GC {gc * 100:0.0}% · {orfs.Count} ORFs · {frames.Count} frames translated · {GalleryData.Ms(sw)} · kernels {BiotechInfo.Backend}";
            _resultText = $"{seq.Id}\nlength\t{seq.Length}\nGC\t{gc:P2}\nORFs\t{orfs.Count}\n" + string.Join('\n', orfs.Select(o => o.ToString())) +
                          "\n" + frameText + string.Join('\n', cutters.Select(x => $"{x.Enzyme.Name}\t{string.Join(',', x.Cuts.Select(c => c.Position))}")) +
                          (oligo is null ? "" : $"\nTm(NN)\t{tmNn:0.0}");
            Code = Snippet(seq, oligoText);
        }
        catch (Exception ex)
        {
            Status = "Could not analyse the input: " + ex.Message;
        }
    }

    private static string Preview(string s, int n) => s.Length <= n ? s : s[..n] + "…";

    private static string Snippet(Sequence seq, string oligo) => $$"""
        using Biotech.Net.Core;
        using Biotech.Net.Sequences;
        using Biotech.Net.Sequences.IO;

        // Parse FASTA (or build directly with Sequence.Dna)
        var seq = Sequence.Dna({{GalleryData.Literal(seq.ToString())}}, id: "{{seq.Id}}");

        // Composition and GC
        IReadOnlyDictionary<char, long> counts = seq.Composition();
        double gc = seq.GcContent();                                   // {{seq.GcContent():0.0000}}
        var gcProfile = SlidingWindow.GcContent(seq, window: {{Math.Clamp(seq.Length / 20, 30, 500)}}, step: {{Math.Max(1, Math.Clamp(seq.Length / 20, 30, 500) / 4)}});

        // Reverse complement and six-frame translation
        Sequence rc = seq.ReverseComplement();
        foreach (var (frame, protein) in SixFrame.Translate(seq))
            Console.WriteLine($"{frame,3} {protein.Preview(60)}");

        // Open reading frames (ATG → stop, ≥ 150 nt, both strands)
        var orfs = OrfFinder.Find(seq, new OrfOptions { MinLength = 150 });
        foreach (var orf in orfs) Console.WriteLine($"{orf}  {orf.Translate(seq).Preview(30)}");

        // Restriction map: enzymes that cut 1–3 times
        var cutters = RestrictionEnzymes.All
            .Select(e => (Enzyme: e, Cuts: RestrictionDigest.FindCuts(seq, e)))
            .Where(x => x.Cuts.Count is > 0 and <= 3);
        var digest = RestrictionDigest.Digest(seq, circular: false, "EcoRI", "BamHI");

        // Primer melting temperature
        var primer = Sequence.Dna("{{oligo}}");
        double tm = MeltingTemperature.NearestNeighbor(primer);        // SantaLucia 1998, 50 mM Na+
        double tmWallace = MeltingTemperature.Wallace(primer);
        double dG = MeltingTemperature.FreeEnergy37(primer);
        """;
}
