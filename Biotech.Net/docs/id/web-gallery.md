# Web Gallery dan C# Playground

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

`samples/Biotech.Net.WebGallery` adalah aplikasi Blazor WebAssembly (standalone): kembaran [Gallery desktop](samples.md) versi peramban, dengan halaman, data, dan desain "bench instrument" yang sama. Ada juga **C# Playground** yang mengompilasi dan menjalankan cuplikan Biotech.Net langsung di tab peramban. Semua komputasi berjalan di sisi klien dengan pustaka Biotech.Net. Tidak ada server dan tidak ada data yang diunggah.

Versi live: [https://dotnetvibecoderz.github.io/Vibe_Bio/](https://dotnetvibecoderz.github.io/Vibe_Bio/) (diterbitkan oleh `.github/workflows/biotechnet-webgallery.yml` setiap push ke `main` yang menyentuh aplikasi atau pustaka).

![DNA Explorer di peramban](../screenshots/webgallery-dna-explorer.png)

## Halaman

| Halaman | Fungsinya | Modul |
|---|---|---|
| **DNA Explorer** | komposisi, jendela GC, translasi enam frame, ORF, peta restriksi, Tm primer | Core, Sequences |
| **Genome Browser** | region 12 kb bawaan (FASTA + GFF3 + VCF): geser, perbesar, nyalakan/matikan track, prediksi efek varian; gambarnya SVG `GenomeBrowser` | Genomics, Visualization |
| **Protein Explorer** | ProtParam, profil Kyte–Doolittle, segmen TM; crambin (1CRN, sudah ada di pustaka) dalam 3Dmol.js, jejak Cα, peta kontak, DSSP | Proteomics, StructuralBiology, Visualization |
| **Drug Discovery Lab** | SMILES → deskriptor, aturan drug-likeness (Lipinski dan lainnya), penanda ADMET, alert, gambar 2-D (SVG), kemiripan ECFP4 + MACCS terhadap 18 obat yang disetujui | Chemistry, DrugDiscovery, Visualization |
| **CRISPR Designer** | guide pada target yang ditempel untuk nuklease apa pun, skor on-target Rule Set 1, spesifisitas CFD / MIT terhadap latar 60 kb | SyntheticBiology |
| **Transcriptomics** | tabel hitungan bawaan (`wwwroot/data/ifn-demo-counts.tsv`, 400 gen × 6 sampel, data demo sintetis) diproses `DESeq.Run`, dengan plot volcano dan MA | Transcriptomics, Visualization |
| **C# Playground** | sunting dan jalankan cuplikan; sepuluh contoh diambil dari tutorial | semua modul di atas |
| **About** | kredit, paket, runtime, apa yang bisa jalan di peramban | Core |

Setiap halaman kecuali Playground punya tab **Live code** berisi C# yang menghasilkan hasil di layar. **Open in Playground** memuat kode itu ke editor, dan kode itu bisa langsung dijalankan tanpa diubah.

## Menjalankan secara lokal

```bash
dotnet workload install wasm-tools      # sekali saja; lihat "Gambar dan workload wasm-tools" di bawah
cd samples/Biotech.Net.WebGallery
dotnet run                              # server dev, http://localhost:5xxx
```

Untuk menyajikan build rilis seperti di Pages:

```bash
dotnet publish samples/Biotech.Net.WebGallery -c Release -p:SkipRustBuild=true -o publish
cd publish/wwwroot && python -m http.server 8080      # server berkas statis apa pun bisa dipakai
```

`python -m http.server` tidak punya fallback SPA, jadi mulailah dari `http://localhost:8080/` lalu bernavigasi di dalam aplikasi. Tautan langsung seperti `/playground` butuh server yang mengembalikan `index.html` untuk path yang tidak dikenal, dan di Pages tugas itu dikerjakan `404.html`.

### Gambar dan workload wasm-tools

Semua gambar pustaka berformat SVG. Penulis SVG (Matplotlib.Net) mengukur teks dengan SkiaSharp, jadi pustaka native SkiaSharp harus di-link ke runtime WebAssembly (`dotnet.native.wasm`). Untuk itu dibutuhkan workload `wasm-tools`. Proyek ini mereferensikan `SkiaSharp.NativeAssets.WebAssembly` hanya jika workload tersebut terpasang:

- **Dengan wasm-tools:** gambar tampil. Publish memakan beberapa menit karena runtime di-link ulang dengan emscripten.
- **Tanpa wasm-tools:** aplikasi tetap bisa dibangun dan dijalankan. Ini yang terjadi pada `dotnet build Biotech.Net.slnx` biasa, jadi CI tidak butuh workload ini. Halaman yang menggambar gambar pustaka menampilkan pemberitahuan di tempatnya, dan bagian lain tetap berjalan.
- `-p:BiotechWasmSkia=false` mematikannya secara eksplisit.

## Apa yang jalan di peramban dan apa yang tidak

| | Status |
|---|---|
| Sekuens, penjajaran, VCF/GFF3, efek varian, uji gaya DESeq2, ProtParam, DSSP, SASA, SMILES, deskriptor, ADMET, fingerprint, desain CRISPR, optimasi kodon | jalan di peramban |
| Gambar SVG (genome browser, volcano, MA, jejak Cα, peta kontak, gambar molekul) | jalan (butuh publish dengan wasm-tools seperti di atas) |
| Kernel Rust | **tidak dimuat**: kode native tidak bisa berjalan di WebAssembly. Probe di `Biotech.Net.Native` gagal tanpa galat dan setiap kernel memakai versi C# managed (halaman About dan baris status menampilkan `kernels Managed`). Hasilnya identik, tetapi input besar berjalan lebih lambat. |
| Thread | runtime peramban hanya punya satu thread. `Parallel.For` berjalan berurutan, dan komputasi panjang membuat halaman sibuk sampai selesai. |
| Penampil 3-D | memuat 3Dmol.js 2.4.2 dari jsDelivr di dalam iframe sandbox, jadi butuh koneksi jaringan |
| Modul AI (unduhan Hugging Face, inferensi lokal, LLM copilot), kernel GPU, klien data cloud, mesin workflow, path berkas | **tidak tersedia** di peramban; pakai Gallery desktop atau notebook |

## Playground

- **Editor:** text area monospace biasa. Tab membuat indentasi dan Ctrl+Enter (Cmd+Enter di macOS) menjalankan kode. Cuplikan ditulis sebagai top-level statement. `System`, `System.Linq`, `System.Collections.Generic`, `System.IO`, `System.Text`, `System.Threading.Tasks`, dan `Biotech.Net.Core` sudah diimpor, begitu juga `Show(...)`.
- **`Show(x)`** menampilkan `x` di bawah output konsol:
  - `BioPlot` menjadi SVG inline;
  - `Structure` membuka penampil 3Dmol.js;
  - `Molecule` digambar;
  - string SVG ditampilkan apa adanya;
  - selain itu dicetak sebagai teks.
- **Preset** (dari `docs/en/tutorials`): sekuens, penjajaran berpasangan, MSA dan pohon neighbour-joining, ProtParam, analisis struktur crambin (dengan tampilan 3-D), deskriptor aspirin (dengan gambar molekul), ADMET dan PAINS, optimasi kodon, efek varian, serta DESeq dengan plot volcano.

### Cara kerjanya

1. **Memuat compiler.** Membuka Playground memuat Roslyn secara lazy (`Microsoft.CodeAnalysis` dan `Microsoft.CodeAnalysis.CSharp`, 2,7 MB Brotli). Halaman lain tidak mengunduhnya.
2. **Membangun referensi dari memori.** Referensi diambil dari assembly yang sudah dimuat di runtime WebAssembly, dibaca langsung dari memori dengan `AssemblyExtensions.TryGetRawMetadata`. Tidak ada reference assembly yang diunduh. Jumlahnya 52 assembly: framework yang sudah di-trim plus semua modul Biotech.Net.
3. **Mengompilasi dan menjalankan.** Setiap eksekusi dikompilasi di memori dengan `concurrentBuild: false` lalu dimuat dengan `Assembly.Load(byte[])`. Entry point dipanggil selama `Console.Out` dialihkan.
4. **Kecepatan.** Pemanasan sekitar 1 detik. Setelah itu cuplikan biasa terkompilasi dalam 60–700 ms dan berjalan dalam hitungan milidetik.

### Keterbatasan

- **Member framework yang hilang.** Framework di-trim. Aplikasi mempertahankan `System.Runtime`, `System.Linq`, `System.Collections`, dan `System.Console` secara utuh, ditambah daftar tipe umum `System.Private.CoreLib` (`Playground/TrimmerRoots.xml`: string, matematika, span, koleksi, task, interpolated string, dan seterusnya). Cuplikan yang memakai member framework lain yang tidak pernah dipanggil gallery bisa gagal dikompilasi dengan pesan "does not contain a definition". Assembly Biotech.Net tidak di-trim.
- **Eksekusi panjang membekukan halaman.** Cuplikan berjalan di thread UI, jadi loop tak berujung membekukan tab dan hanya bisa dihentikan dengan memuat ulang halaman.
- **Tanpa referensi NuGet.** Cuplikan tidak bisa menambah paket NuGet. Yang tersedia hanya assembly yang dibawa aplikasi.
- **Editor sederhana.** Editor tidak punya IntelliSense. Galat compiler ditampilkan lengkap dengan baris dan kolom.

## Ukuran unduhan

Angka berikut diukur pada publish Release, dengan build wasm-tools:

| | Ukuran |
|---|---|
| `wwwroot` di disk (setiap berkas plus salinan `.br` dan `.gz`) | 63 MB |
| `_framework` tanpa kompresi | 31,4 MB (120 berkas) |
| `_framework` Brotli | 9,5 MB, 2,7 MB di antaranya Roslyn (lazy, hanya untuk Playground) |
| Kunjungan pertama (Brotli, tanpa Roslyn) | sekitar 6,8 MB |

**Penyusun ukurannya:**

- **Trimming.** Publish memakai `TrimMode=partial`: framework di-trim, sedangkan pustaka Biotech.Net dibiarkan utuh karena cuplikan dikompilasi terhadapnya.
- **Runtime.** Bagian tunggal terbesar adalah `dotnet.native.wasm` (5,6 MB mentah), yang sudah termasuk Skia.
- **Dependensi yang tak terpakai.** ILGPU dan Parquet ikut masuk lewat paket `Gravicode.Science`, tidak pernah dipakai di peramban, dan ditandai trimmable.

**Kompresi di GitHub Pages.** Pages tidak menyajikan berkas `.br` yang sudah dikompresi, jadi workflow Pages menghapusnya. Pages mengompresi respons sendiri dengan gzip; total salinan gzip 13 MB termasuk Roslyn.

## Deploy ke GitHub Pages

Workflow disimpan di `build/workflows/biotechnet-webgallery.yml`. Workflow ini tidak melakukan apa pun sampai disalin ke `.github/workflows/` di root monorepo. Langkah-langkahnya:

1. memasang .NET 10 dan workload `wasm-tools`;
2. menjalankan `dotnet publish samples/Biotech.Net.WebGallery -c Release -p:SkipRustBuild=true`;
3. mengganti `<base href="/" />` menjadi `/<repository>/`, di sini `/Vibe_Bio/`;
4. menyalin `index.html` ke `404.html`, supaya tautan langsung sampai ke router Blazor;
5. menambah `.nojekyll`, supaya `_framework/` ikut disajikan;
6. menghapus berkas `.br` dan `.gz`;
7. mengunggah folder dengan `actions/upload-pages-artifact` lalu men-deploy-nya dengan `actions/deploy-pages`.

**Pengaturan repositori:**

- **Sumber Pages.** Settings → Pages → Build and deployment → Source: **GitHub Actions**.
- **Environment.** Environment `github-pages`, yang dibuat GitHub pada deploy pertama, harus mengizinkan branch `main`. Itu aturan bawaannya.
- **Sifat publik.** Situs langsung publik begitu di-deploy.

## Berkas

| Path | Isinya |
|---|---|
| `samples/Biotech.Net.WebGallery/Pages/*.razor` | satu komponen per halaman |
| `Shared/PageFrame.razor` | header, pita kromatogram, tab Result / Live code, baris status |
| `Services/Bench.cs` | token desain, tema plot, penyematan SVG yang aman, kromatogram |
| `Services/Charts.cs` | instrumen SVG inline kecil (batang, profil, peta linear) |
| `Playground/SnippetCompiler.cs` | kompilasi dan eksekusi dengan Roslyn |
| `Playground/Presets.cs` | sepuluh contoh |
| `Playground/Display.cs` | `Show(...)` |
| `wwwroot/css/app.css` | gaya chassis |
| `../Biotech.Net.Gallery/Data/*` | region untuk Genome Browser, dipakai bersama dengan Gallery desktop sebagai embedded resource |

## Tangkapan layar

![Genome Browser](../screenshots/webgallery-genome-browser.png)

![Protein Explorer](../screenshots/webgallery-protein-explorer.png)

![Drug Discovery Lab](../screenshots/webgallery-drug-discovery-lab.png)

![CRISPR Designer](../screenshots/webgallery-crispr-designer.png)

![Transcriptomics](../screenshots/webgallery-transcriptomics.png)

![C# Playground](../screenshots/webgallery-playground.png)

![Playground: Show(volcano)](../screenshots/webgallery-playground-volcano.png)

![Playground: Show(structure)](../screenshots/webgallery-playground-structure.png)

![Live code](../screenshots/webgallery-live-code.png)

![About](../screenshots/webgallery-about.png)

![Lebar ponsel](../screenshots/webgallery-mobile.png)

---
[Indeks dokumentasi](index.md) · [Contoh](samples.md) · English: [Web Gallery](../en/web-gallery.md)
