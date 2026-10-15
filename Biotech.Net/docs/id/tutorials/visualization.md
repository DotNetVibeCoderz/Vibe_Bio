# Tutorial: Visualisasi

Paket: `Gravicode.Biotech.Net.Visualization` (di atas Visual.Net / Matplotlib.Net / Plotly.Net dan SkiaSharp)

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Genome browser, heatmap dengan dendrogram, volcano dan MA plot, embedding PCA/UMAP, sequence logo, dot plot, jejak Cα protein dan peta kontak, penggambaran molekul 2-D, serta tampilan 3-D interaktif lewat 3Dmol.js.

Semua plot diturunkan dari `BioPlot` dan diekspor dengan cara yang sama:

| Metode | Keluaran |
|---|---|
| `ToPng(scale)` | byte PNG (`scale` untuk DPI tinggi) |
| `ToSvg()` | teks SVG |
| `ToPdf()` | byte PDF |
| `ToHtml(interactive)` | HTML — Plotly bila `SupportsInteractive` bernilai true, selain itu gambar tersemat |
| `Save(path)` | menulis format sesuai ekstensi (`.png`, `.svg`, `.pdf`, `.html`) |
| `ToPngDataUri()` | `data:image/png;base64,...` untuk notebook dan halaman web |

Semuanya headless dan lintas platform: tidak perlu display server.

## 1. Tema

```csharp
using Biotech.Net.Visualization;

BioTheme.Default = BioTheme.Light;        // atau BioTheme.Dark
var custom = BioTheme.Light.With(t => t with { Name = "lab", FontFamily = "Inter", FontSize = 13 });
Console.WriteLine($"{custom.Name}: {custom.Palette.Count} warna kategori, nukleotida A = {custom.NucleotideColor('A')}");
```

Sebuah tema membawa palet kategori, warna naik/turun/netral, warna per nukleotida dan per residu, warna unsur untuk molekul, font, dan colormap yang dipakai heatmap.

## 2. Plot ekspresi

```csharp
using Biotech.Net.Core;
using Biotech.Net.Transcriptomics;
using Biotech.Net.Visualization;

var rng = new Random(42);
int genes = 1200, samples = 8;
var values = new double[genes * samples];
for (int g = 0; g < genes; g++)
{
    double mean = Math.Exp(2 + rng.NextDouble() * 6);
    double fold = g < 40 ? 4 : g < 80 ? 0.25 : 1;
    for (int s = 0; s < samples; s++)
        values[g * samples + s] = Math.Round(mean * (s >= 4 ? fold : 1)
            * Math.Exp(0.25 * (rng.NextDouble() + rng.NextDouble() + rng.NextDouble() - 1.5)));
}
var counts = new Dataset("counts",
    Enumerable.Range(1, genes).Select(i => $"gene{i:D4}").ToArray(),
    ["ctl1", "ctl2", "ctl3", "ctl4", "trt1", "trt2", "trt3", "trt4"], values);
counts.Annotate("condition", "ctl", "ctl", "ctl", "ctl", "trt", "trt", "trt", "trt");
var result = DESeq.Run(counts, "condition", "ctl", "trt");

var volcano = new VolcanoPlot(result) { Title = "Perlakuan vs kontrol", LabelTop = 8, FoldChangeThreshold = 1 };
File.WriteAllBytes(Path.Combine(Path.GetTempPath(), "volcano.png"), volcano.ToPng(scale: 1.5));
Console.WriteLine(volcano.Counts());        // (Naik, Turun, TidakSignifikan)

var ma = new MaPlot(result);

// Heatmap terklaster untuk gen teratas (z-score baris, average linkage)
var vst = Transforms.Vst(counts);
var top = result.Top(40).Select(r => vst.RowIndex(r.Gene)).ToArray();
var heat = Heatmap.FromDataset(vst.SelectRows(top, "40 gen teratas"));
heat.ColumnGroups = counts.ColumnAnnotations["condition"];
Console.WriteLine($"heatmap {heat.Rows}×{heat.Columns}, klaster baris: {heat.ClusterRows}");
```

![Volcano plot](../../images/volcano.png)

![Heatmap](../../images/heatmap.png)

## 3. Genome browser

Track ditambahkan secara berantai dan browser menghitung tingginya sendiri.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Visualization;

var region = GenomicInterval.Parse("chr1:1,001-6,000");
var gene = new Gene("G1", new GenomicInterval("chr1", 1200, 3400, Strand.Forward), "ALPHA");
gene.Transcripts.Add(new Transcript("ALPHA-201", "G1",
    [new("chr1", 1200, 1500, Strand.Forward), new("chr1", 2000, 2300, Strand.Forward), new("chr1", 3000, 3400, Strand.Forward)])
    { CodingRegion = new GenomicInterval("chr1", 1350, 3200, Strand.Forward) });

double[] depth = Enumerable.Range(0, 5000).Select(i => 25 + 12 * Math.Sin(i / 280.0)).ToArray();
var variants = new[] { new Variant("chr1", 1400, "A", "G", "rs1001"), new Variant("chr1", 4000, "AT", "A") };

var browser = new GenomeBrowser(region, width: 1000) { Title = "Lokus ALPHA" }
    .AddRuler()
    .AddGenes([gene])
    .AddCoverage(depth, start: 1000, name: "Kedalaman RNA-seq")
    .AddVariants(variants);
Console.WriteLine($"{browser.Tracks.Count} track, tinggi {browser.Height} px, {browser.BasesPerPixel:F1} bp/px");
```

![Genome browser](../../images/genome-browser.png)

Track yang tersedia: `RulerTrack`, `GeneTrack`, `FeatureTrack`, `CoverageTrack`, `VariantTrack`, dan `SequenceTrack` (huruf basa pada zoom penuh). Turunkan `GenomeTrack` untuk membuat track sendiri.

## 4. Sekuens, struktur, dan molekul

```csharp
using Biotech.Net.Chemistry;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Visualization;

string[] sites =
[
    "TAATGTGAGTTAGCTCACTCAT", "TAATGTGAGTTAGCACACTCAT", "TGATGTGAGTTAGCTCACTCAT",
    "TAATGTGAGTTAGCTCACTCTT", "TAATGCGAGTTAGCTCACTCAT", "TAATGTGAGCTAGCTCACTCAT",
];
var logo = new SequenceLogo(sites) { Title = "Situs mirip CRP" };
Console.WriteLine($"{logo.Length} kolom, maksimum {logo.MaxBits:F0} bit");

var crambin = SampleStructures.Crambin();
var trace = ProteinViewer.CaTrace(crambin);
var contacts = ProteinViewer.ContactMap(crambin);
string html3d = ProteinViewer.To3DHtml(crambin, new Viewer3DOptions { Style = MolecularStyle.Cartoon, Spin = true });

var viewer = new MoleculeViewer(Molecule.FromSmiles("CC(=O)Nc1ccc(O)cc1", "paracetamol"), 520, 380);
Console.WriteLine($"{viewer.Coordinates2D().Length} posisi atom ditata");
```

![Sequence logo](../../images/sequence-logo.png)

![Molekul](../../images/molecule.png)

![Peta kontak protein](../../images/protein-contact-map.png)

## Batasan

- Tampilan 3-D protein dan molekul berupa halaman HTML yang memuat 3Dmol.js dari CDN, jadi perlu peramban dan akses jaringan. Keluaran statis (jejak Cα, peta kontak, gambar 2-D) dirender lokal tanpa jaringan.
- Keluaran Plotly interaktif hanya tersedia untuk plot yang `SupportsInteractive`-nya true (volcano, MA, embedding, heatmap, peta kontak, dot plot). Genome browser, sequence logo, jejak Cα, dan penampil molekul diekspor sebagai gambar atau HTML dengan gambar tersemat.
- Klasterisasi heatmap dilakukan hierarkis di memori; puluhan ribu baris akan lambat dan tidak terbaca. Saring dulu gen yang ingin Anda lihat.
- Di Linux, SkiaSharp memerlukan paket native asset; `Gravicode.Biotech.Net.Visualization` sudah bergantung pada `SkiaSharp.NativeAssets.Linux.NoDependencies`.

---
Sebelumnya: [AI](ai.md) · Selanjutnya: [Notebook](notebooks.md) · [Indeks tutorial](index.md)
