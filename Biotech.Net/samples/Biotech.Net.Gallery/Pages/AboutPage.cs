using System.Runtime.InteropServices;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Layout;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using Avalonia.Platform;
using Biotech.Net.Core;
using Biotech.Net.Gallery.Controls;
using Biotech.Net.Gallery.Shell;
using Biotech.Net.Native;

namespace Biotech.Net.Gallery.Pages;

public sealed class AboutPage : GalleryPage
{
    public override string Title => "About";
    public override string Modules => "BIOTECH.NET · GRAVICODE STUDIOS";
    public override string Lede => BiotechInfo.Tagline + ". Bioconductor + Biopython for .NET, with Rust-accelerated kernels.";
    public override char Channel => 'G';

    private static readonly (string Package, string What, char Channel)[] Packages =
    [
        ("Core", "Sequence, Alphabet, GeneticCode, Genome, Variant, Dataset", 'A'),
        ("Native", "Rust kernels with managed fallback (reverse complement, translation, alignment)", 'G'),
        ("Sequences", "FASTA/FASTQ/GenBank/EMBL/PDB I/O, ORFs, restriction, Tm, k-mers", 'A'),
        ("Alignments", "Pairwise and progressive MSA, substitution matrices, phylogenetics", 'C'),
        ("Genomics", "VCF, GFF3, read mapping, SNV calling, VEP, assembly, population genetics", 'C'),
        ("Transcriptomics", "DESeq-style and limma DE, PCA, single-cell clustering, GSEA", 'C'),
        ("Proteomics", "ProtParam, scales, mass spectrometry, protein classification", 'T'),
        ("StructuralBiology", "PDB/mmCIF, DSSP, SASA, contacts, superposition, pockets", 'T'),
        ("Chemistry", "SMILES, SMARTS, fingerprints, descriptors, depiction", 'G'),
        ("DrugDiscovery", "ADMET, alerts, screening, QSAR, docking, lead optimisation", 'G'),
        ("SyntheticBiology", "CRISPR design, primers, Gibson/Golden Gate, codon optimisation, circuits", 'A'),
        ("AI", "DNABERT/ESM/Geneformer preprocessing and inference, bio copilot", 'C'),
        ("Visualization", "Genome browser, heatmaps, volcano, PCA/UMAP, protein and molecule viewers", 'T'),
        ("Notebooks", ".NET Interactive formatters for every type", 'A'),
    ];

    protected override Control Build()
    {
        var logo = new Image { Width = 120, Height = 120, VerticalAlignment = VerticalAlignment.Top };
        try { logo.Source = new Bitmap(AssetLoader.Open(new Uri("avares://Biotech.Net.Gallery/Assets/icon-512.png"))); } catch { }

        var credits = new TextBlock
        {
            Text = BiotechInfo.CreditsId,
            FontFamily = Palette.Display,
            FontSize = 26,
            FontWeight = FontWeight.SemiBold,
            Foreground = Palette.InkBrush,
            TextWrapping = TextWrapping.Wrap,
        };
        var hero = Ui.Columns("Auto,*", 28, logo, Ui.Column(8,
            Ui.Eyebrow("Credits"),
            credits,
            Ui.Text(BiotechInfo.Credits, "lede"),
            Ui.Mono(BiotechInfo.Banner, 11, Palette.Ink3Brush)));

        var runtime = Ui.KeyValues(
        [
            ("Biotech.Net", BiotechInfo.Version, Palette.DyeG),
            ("Kernel backend", BioKernels.Backend.ToString(), BioKernels.Backend == KernelBackend.Native ? Palette.DyeA : null),
            ("Native kernels", BioKernels.NativeVersion, null),
            (".NET runtime", RuntimeInformation.FrameworkDescription, null),
            ("OS", RuntimeInformation.OSDescription, null),
            ("Architecture", RuntimeInformation.ProcessArchitecture.ToString(), null),
            ("UI", "Avalonia 12 · AvaloniaEdit · CommunityToolkit.Mvvm", null),
        ]);

        var grid = new Grid { ColumnDefinitions = new ColumnDefinitions("*,*") };
        for (int i = 0; i < Packages.Length; i++)
        {
            var (pkg, what, ch) = Packages[i];
            if (i % 2 == 0) grid.RowDefinitions.Add(new RowDefinition(GridLength.Auto));
            var cell = new Border
            {
                BorderBrush = Palette.RuleBrush,
                BorderThickness = new Thickness(0, 0, 0, 1),
                Padding = new Thickness(0, 8, 16, 8),
                Child = Ui.Columns("12,*", 0,
                    new Border { Width = 6, Height = 6, Background = new SolidColorBrush(Palette.Channel(ch)), VerticalAlignment = VerticalAlignment.Top, Margin = new Thickness(0, 6, 0, 0) },
                    Ui.Column(2,
                        new TextBlock { Text = "Gravicode.Biotech.Net." + pkg, FontFamily = Palette.Mono, FontSize = 12, Foreground = Palette.InkBrush },
                        Ui.Text(what, "caption"))),
            };
            Grid.SetRow(cell, i / 2);
            Grid.SetColumn(cell, i % 2);
            grid.Children.Add(cell);
        }

        TraceSequence = GalleryData.TextToBases(BiotechInfo.CreditsId);
        TraceLabel = "TRACE · credits encoded 2 bases/char";
        Status = BiotechInfo.Banner;
        Code = """
            using Biotech.Net.Core;

            Console.WriteLine(BiotechInfo.Banner);
            Console.WriteLine(BiotechInfo.CreditsId);   // Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil
            Console.WriteLine($"Kernels: {BiotechInfo.Backend}");

            // dotnet add package Gravicode.Biotech.Net.Core
            // dotnet run --project samples/Biotech.Net.Gallery -- --screenshots docs/screenshots
            """;

        return Ui.Scroll(Ui.Column(14,
            Ui.Panel("Biotech.Net", hero),
            Ui.Columns("3*,2*", 14,
                Ui.Panel("Packages", grid, $"{Packages.Length} on NuGet · Gravicode.Biotech.Net.*"),
                Ui.Panel("Runtime", runtime))));
    }
}
