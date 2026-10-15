# Contoh

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Ada dua cara melihat Biotech.Net bekerja: aplikasi desktop **Gallery** (Avalonia) dan **notebook contoh**.

## Biotech.Net Gallery

Aplikasi pamer lintas platform berbasis Avalonia 12 di `samples/Biotech.Net.Gallery`. Setiap halaman menjalankan kode pustaka yang sebenarnya, menampilkan C# yang menghasilkan hasilnya di tab *Live code*, dan bisa mengekspor gambarnya sebagai PNG atau SVG.

```bash
cd samples/Biotech.Net.Gallery
dotnet run -c Release
```

| Halaman | Yang ditunjukkan | Modul yang dipakai |
|---|---|---|
| **DNA Explorer** | statistik sekuens, jendela GC, ORF, digesti restriksi, translasi | Core, Sequences |
| **Genome Browser** | gen, transkrip, cakupan, varian, dan fitur dalam satu tampilan | Core, Genomics, Visualization |
| **Protein Explorer** | sifat ProtParam, struktur sekunder, kontak, jejak Cα, tampilan 3-D | Proteomics, StructuralBiology, Visualization |
| **Drug Discovery Lab** | deskriptor, aturan kemiripan obat, triase ADMET, penapisan kemiripan | Chemistry, DrugDiscovery |
| **AI Playground** | katalog model, praproses, embedding, dan kopilot LLM | AI |
| **CRISPR Designer** | pencacahan guide, skor on-target, dan spesifisitas off-target | SyntheticBiology |

![DNA Explorer](../screenshots/gallery-dna-explorer.png)

![Genome Browser](../screenshots/gallery-genome-browser.png)

![Protein Explorer](../screenshots/gallery-protein-explorer.png)

![Drug Discovery Lab](../screenshots/gallery-drug-discovery-lab.png)

![AI Playground](../screenshots/gallery-ai-playground.png)

![CRISPR Designer](../screenshots/gallery-crispr-designer.png)

Tab *Live code* menampilkan cuplikan di balik hasil yang sedang tampil, dan halaman About memuat kreditnya:

![Live code](../screenshots/gallery-live-code.png)

![About](../screenshots/gallery-about.png)

Data contoh untuk halaman Genome Browser ada di `samples/Biotech.Net.Gallery/Data` (`region.fa`, `region.gff3`, `region.vcf`).

## Notebook contoh

Enam notebook `.dib` dwibahasa di `notebooks/`, untuk Polyglot Notebooks (VS Code), Visual Studio, atau Jupyter dengan kernel .NET:

| Notebook | Isi |
|---|---|
| `01-getting-started.dib` | paket, formatter kaya, tipe inti |
| `02-sequences-and-io.dib` | I/O FASTA/FASTQ/GenBank, ORF, motif, digesti restriksi |
| `03-genome-browser-and-plots.dib` | genome browser, heatmap, sequence logo, ekspor statis dan interaktif |
| `04-rnaseq-differential-expression.dib` | eksekusi DESeq, volcano, MA, PCA, heatmap terklaster |
| `05-protein-structure.dib` | ProtParam, DSSP, SASA, superposisi, penampil 3-D |
| `06-cheminformatics.dib` | SMILES, deskriptor, sidik jari, kemiripan, penggambaran |

Masing-masing diawali `#r "nuget: Gravicode.Biotech.Net.Notebooks, 0.1.0"` berikut alternatif terkomentar yang merujuk build lokal, sehingga bisa dipakai baik dari NuGet maupun dari hasil klon repositori.

## Gambar di dokumentasi ini

Semua berkas di `docs/images/` dihasilkan oleh pustaka ini sendiri — volcano, heatmap, genome browser, sequence logo, gambar molekul, PCA, dot plot, jejak Cα, dan peta kontak. Kodenya ada di [tutorial visualisasi](tutorials/visualization.md).

---
[Indeks dokumentasi](index.md) · English: [samples](../en/samples.md)
