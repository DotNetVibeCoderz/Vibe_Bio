using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.Threading;
using Biotech.Net.AI.Copilot;
using Biotech.Net.AI.Hub;
using Biotech.Net.AI.Inference;
using Biotech.Net.AI.Models;
using Biotech.Net.AI.Preprocessing;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.StructuralBiology;
using Gravicode.HFNet.GraviTokenizers;

namespace Biotech.Net.Gallery.Pages;

public sealed class AiPlaygroundPage : GalleryPage
{
    public override string Title => "AI Playground";
    public override string Modules => "BIOTECH.NET.AI · HF.NET";
    public override string Lede => "Prepare sequences for biological language models offline (k-mers, ESM tokens), browse the model catalogue, and — when a key is in the environment — score protein variants with ESM-2 or ask the copilot.";
    public override char Channel => 'C';

    /// <summary>The ESM-2 vocabulary in id order (identical for every ESM-2 checkpoint).</summary>
    private static readonly string[] EsmVocabulary =
    [
        "<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E", "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y", "M", "H",
        "W", "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>", "<mask>",
    ];

    private const string EsmRepo = "facebook/esm2_t6_8M_UR50D";

    /// <summary>HF.Net's ESM tokenizer over <see cref="EsmVocabulary"/>, built once (it reads the vocabulary from a file).</summary>
    private static readonly Lazy<HfTokenizer> EsmTokenizer = new(() =>
    {
        string path = Path.Combine(Path.GetTempPath(), "biotechnet-gallery-esm2-vocab.txt");
        File.WriteAllLines(path, EsmVocabulary);
        return HfTokenizer.FromEsmVocabulary(path);
    });

    private TextBox _dna = null!, _protein = null!, _mutations = null!, _question = null!;
    private ComboBox _k = null!;
    private readonly ContentControl _prep = new();
    private readonly ContentControl _esmResult = new();
    private readonly ContentControl _copilotResult = new();
    private string _resultText = "";

    public override string ResultText => _resultText;

    private static string? HfToken() =>
        new[] { "HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", "HUGGINGFACE_TOKEN" }.Select(Environment.GetEnvironmentVariable).FirstOrDefault(v => !string.IsNullOrWhiteSpace(v));

    private static bool EsmCached()
    {
        try { return BioModelHub.IsCached(EsmRepo); } catch { return false; }
    }

    protected override Control Build()
    {
        _dna = Ui.Input(GalleryData.Egfp[..96], "seq", 64);
        _protein = Ui.Input(SampleStructures.Crambin().Chains[0].ToSequence().ToString(), "seq", 48);
        _k = new ComboBox { ItemsSource = new[] { 3, 4, 5, 6, 7, 8 }, SelectedItem = 6, MinWidth = 70 };
        _k.SelectionChanged += (_, _) => Preprocess();
        _mutations = Ui.Input("T2A, P5W, C16S, A27G", "seq");
        _mutations.AcceptsReturn = false;
        _question = Ui.Input("What is the GC content and the longest ORF of ATGGCGTAAGCTAGCATGCCCTGA?", "");
        _question.AcceptsReturn = false;

        var inputs = Ui.Panel("Inputs", Ui.Column(10,
            Ui.Labelled("DNA", _dna),
            Ui.Labelled("Protein", _protein),
            Ui.Row(8, Ui.Eyebrow("k"), _k, Ui.Button("Tokenize", Preprocess, "primary"))));

        var body = Ui.Column(14,
            Ui.Columns("340,*", 14, inputs, _prep),
            Ui.Columns("*,*", 14, EsmPanel(), CopilotPanel()),
            Ui.Panel("Model catalogue", Catalog(), $"{BioModelCatalog.All.Count} checkpoints · BioModelCatalog"));
        Preprocess();
        return Ui.Scroll(body);
    }

    private void Preprocess()
    {
        try
        {
            int k = _k.SelectedItem is int kk ? kk : 6;
            string dna = SequencePreprocessor.CleanDna(_dna.Text ?? "");
            string protein = SequencePreprocessor.CleanProtein(_protein.Text ?? "");
            var kmers = SequencePreprocessor.Kmers(dna, k);
            var tok = EsmTokenizer.Value;
            int[] ids = tok.Encode(protein).ToIdArray();
            string dnabert = BioModelCatalog.DnaBert6.Preprocess(dna);
            string protbert = BioModelCatalog.ProtBert.Preprocess(protein);

            var kmerWrap = new WrapPanel();
            foreach (var (km, i) in kmers.Take(48).Select((x, i) => (x, i)))
                kmerWrap.Children.Add(Token(km, Palette.Channel(km[0]), i.ToString(CultureInfo.InvariantCulture)));
            var esmWrap = new WrapPanel();
            foreach (int id in ids.Take(60))
            {
                string t = tok.IdToToken(id);
                esmWrap.Children.Add(Token(t, t.StartsWith('<') ? Palette.Ink3 : Palette.DyeC, id.ToString(CultureInfo.InvariantCulture)));
            }

            _prep.Content = Ui.Column(14,
                Ui.Panel($"Overlapping {k}-mers · DNABERT input", Ui.Column(8, kmerWrap,
                    Ui.Text($"{kmers.Count} tokens from {dna.Length} nt · {kmers.Distinct().Count()} distinct · DNABERT-6 string starts \"{Preview(dnabert, 48)}\"", "unit")),
                    "SequencePreprocessor.Kmers"),
                Ui.Panel("ESM-2 token ids", Ui.Column(8, esmWrap,
                    Ui.Text($"{ids.Length} ids incl. <cls>/<eos> · vocabulary {tok.VocabularySize} · ProtBERT spaced: \"{Preview(protbert, 40)}\"", "unit")),
                    "HfTokenizer.FromEsmVocabulary"));

            TraceSequence = dna;
            TraceLabel = $"TRACE · DNA input {dna.Length} nt";
            Status = $"{kmers.Count} {k}-mers · {ids.Length} ESM tokens · ESM scoring {(EsmAvailable ? "available" : "offline")} · copilot {(ChatClientFactory.FromEnvironment() is null ? "offline" : "available")}";
            _resultText = string.Join(' ', kmers) + "\n" + string.Join(' ', ids);
            Code = Snippet(dna, protein, k);
        }
        catch (Exception ex)
        {
            Status = "Preprocessing failed: " + ex.Message;
        }
    }

    private static bool EsmAvailable => HfToken() is not null || EsmCached();

    private Control EsmPanel()
    {
        var run = Ui.Button("Score with ESM-2", () => _ = ScoreAsync(), "primary");
        run.IsEnabled = EsmAvailable;
        _esmResult.Content = EsmAvailable
            ? Ui.Text(EsmCached() ? "ESM-2 t6 8M is in the local cache." : "HF token found in the environment; the 8M checkpoint (~30 MB) downloads on first use.", "caption")
            : Ui.Notice("ESM scoring is offline",
                "Set HF_TOKEN (or HUGGING_FACE_HUB_TOKEN) in the environment before starting the Gallery, or cache facebook/esm2_t6_8M_UR50D. The token is read at runtime and never stored.",
                Palette.DyeC);
        return Ui.Panel("ESM-2 variant effect", Ui.Column(10,
            Ui.Labelled("Substitutions (1-based, on the protein above)", _mutations),
            run, _esmResult), "masked-marginal log-likelihood ratio");
    }

    private async Task ScoreAsync()
    {
        string wt = SequencePreprocessor.CleanProtein(_protein.Text ?? "");
        var variants = (_mutations.Text ?? "").Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries).ToList();
        _esmResult.Content = Ui.Text("Loading ESM-2 t6 and scoring…", "caption");
        try
        {
            var scores = await Task.Run(() =>
            {
                if (HfToken() is { } token) BioModelHub.Configure(BioModelHub.Options with { Token = token });
                var scorer = EsmVariantScorer.Load(EsmRepo);
                using var model = scorer.Model;
                return scorer.ScoreMany(wt, variants);
            });
            _esmResult.Content = Ui.Table(["Variant", "LLR", "Reading"], ["Auto", "2*", "*"],
                scores.Select(s => new Ui.Cell[]
                {
                    new(string.Join(':', s.Mutations)),
                    new(s.Score.ToString("+0.00;-0.00", CultureInfo.InvariantCulture), s.Score < 0 ? Palette.DyeT : Palette.DyeA, Bar: true, BarValue: Math.Min(1, Math.Abs(s.Score) / 8)),
                    new(s.Score < -2 ? "likely damaging" : s.Score < 0 ? "mildly disfavoured" : "tolerated", Palette.Ink2),
                }));
            Status = $"ESM-2 scored {scores.Count} variants";
        }
        catch (Exception ex)
        {
            _esmResult.Content = Ui.Notice("ESM scoring failed", ex.GetBaseException().Message, Palette.DyeT);
        }
    }

    private Control CopilotPanel()
    {
        var client = ChatClientFactory.FromEnvironment();
        var ask = Ui.Button("Ask the copilot", () => _ = AskAsync(), "primary");
        ask.IsEnabled = client is not null;
        _copilotResult.Content = client is not null
            ? Ui.Text($"Provider: {client.GetType().Name.Replace("ChatClient", "")} · tools from BioToolbox run locally.", "caption")
            : Ui.Notice("Copilot is offline",
                "Set ANTHROPIC_API_KEY, OPENAI_API_KEY, DEEPSEEK_API_KEY or the AZURE_OPENAI_* variables before starting the Gallery. Keys are read from the environment at runtime and never written anywhere.",
                Palette.DyeG);
        return Ui.Panel("Bio copilot", Ui.Column(10, Ui.Labelled("Question", _question), ask, _copilotResult), "tool-calling LLM");
    }

    private async Task AskAsync()
    {
        var client = ChatClientFactory.FromEnvironment();
        if (client is null) return;
        _copilotResult.Content = Ui.Text("Thinking…", "caption");
        try
        {
            var copilot = new BioCopilot(client);
            var answer = await copilot.AskAsync(_question.Text ?? "");
            await Dispatcher.UIThread.InvokeAsync(() =>
            {
                var col = Ui.Column(8, Ui.Text(answer.Text, "body"));
                if (answer.ToolCalls.Count > 0)
                    col.Children.Add(Ui.Text("Tools: " + string.Join(", ", answer.ToolCalls.Select(t => t.Name)), "unit"));
                _copilotResult.Content = col;
            });
        }
        catch (Exception ex)
        {
            _copilotResult.Content = Ui.Notice("Copilot request failed", ex.GetBaseException().Message, Palette.DyeT);
        }
    }

    private static Control Catalog()
    {
        static Color KindColor(BioSequenceKind k) => k switch
        {
            BioSequenceKind.Dna => Palette.DyeA,
            BioSequenceKind.Protein => Palette.DyeT,
            BioSequenceKind.SingleCellExpression => Palette.DyeG,
            _ => Palette.DyeC,
        };
        return Ui.Table(["Model", "Repository", "Input", "Backend", "Params", "Max tokens", "CPU", "Cached"],
            ["Auto", "*", "Auto", "Auto", "Auto", "Auto", "Auto", "Auto"],
            BioModelCatalog.All.Select(m =>
            {
                bool cached;
                try { cached = BioModelHub.IsCached(m.RepoId); } catch { cached = false; }
                return new Ui.Cell[]
                {
                    new(m.DisplayName),
                    new(m.RepoId, Palette.Ink2),
                    new(m.InputKind.ToString(), KindColor(m.InputKind)),
                    new(m.Backend.ToString(), m.Backend == BioModelBackend.NotSupported ? Palette.Ink3 : null),
                    new(m.ApproximateParameters <= 0 ? "—" : m.ApproximateParameters >= 1_000_000_000 ? $"{m.ApproximateParameters / 1e9:0.#}B" : $"{m.ApproximateParameters / 1e6:0}M"),
                    new(m.MaxTokens.ToString(CultureInfo.InvariantCulture)),
                    new(m.KnownToRunOnCpu ? "yes" : "—", m.KnownToRunOnCpu ? Palette.DyeA : Palette.Ink3),
                    new(cached ? "cached" : "—", cached ? Palette.DyeA : Palette.Ink3),
                };
            }));
    }

    private static Control Token(string text, Color color, string index)
    {
        var sp = new StackPanel { Spacing = 0 };
        sp.Children.Add(new TextBlock { Text = text, FontFamily = Palette.Mono, FontSize = 11.5, Foreground = new SolidColorBrush(color) });
        sp.Children.Add(new TextBlock { Text = index, FontFamily = Palette.Mono, FontSize = 8.5, Foreground = Palette.Ink3Brush });
        return new Border
        {
            Child = sp,
            BorderBrush = new SolidColorBrush(color, 0.35),
            BorderThickness = new Thickness(0, 0, 0, 2),
            Background = new SolidColorBrush(Palette.Raised),
            Padding = new Thickness(5, 2),
            Margin = new Thickness(0, 0, 4, 4),
        };
    }

    private static string Preview(string s, int n) => s.Length <= n ? s : s[..n] + "…";

    private static string Snippet(string dna, string protein, int k) => $$"""
        using Biotech.Net.AI.Copilot;
        using Biotech.Net.AI.Inference;
        using Biotech.Net.AI.Models;
        using Biotech.Net.AI.Preprocessing;

        // Offline preprocessing: no model or network needed
        string dna = SequencePreprocessor.CleanDna({{GalleryData.Literal(dna, 48)}});
        IReadOnlyList<string> kmers = SequencePreprocessor.Kmers(dna, k: {{k}});
        string dnabertInput = BioModelCatalog.DnaBert6.Preprocess(dna);        // "ATGGTG TGGTGA …"

        string protein = SequencePreprocessor.CleanProtein({{GalleryData.Literal(protein, 48)}});
        string protbertInput = BioModelCatalog.ProtBert.Preprocess(protein);   // "T T C C P S …"

        // Model catalogue
        foreach (var m in BioModelCatalog.ForKind(BioSequenceKind.Protein))
            Console.WriteLine($"{m.DisplayName} {m.RepoId} {m.Backend}");

        // ESM-2 variant effect (downloads facebook/esm2_t6_8M_UR50D once; HF_TOKEN read from the environment)
        var scorer = EsmVariantScorer.Load("facebook/esm2_t6_8M_UR50D");   // HF.Net runs ESM-2
        int[] ids = scorer.Model.Tokenizer.Encode(protein).ToIdArray();
        foreach (var s in scorer.ScoreMany(protein, ["T2A", "P5W", "C16S"]))
            Console.WriteLine(s);   // negative LLR = less likely than wild type

        // Copilot with local bioinformatics tools (ANTHROPIC_API_KEY / OPENAI_API_KEY / …)
        if (ChatClientFactory.FromEnvironment() is { } client)
        {
            var copilot = new BioCopilot(client);
            var answer = await copilot.AskAsync("What is the GC content of ATGGCGTAAGC?");
            Console.WriteLine(answer.Text);
        }
        """;
}
