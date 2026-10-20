# Instalasi

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Kebutuhan

- **.NET 10 SDK** (10.0.100 atau lebih baru; `global.json` memakai `latestFeature`).
- **Toolchain Rust** — opsional, hanya bila ingin membangun ulang kernel native sendiri. Tanpa Rust semuanya tetap berjalan lewat fallback managed.
- Runtime identifier yang didukung untuk kernel native: `win-x64`, `linux-x64`, `osx-arm64`.

## Dari NuGet

Tambahkan modul yang diperlukan saja; tiap paket menarik dependensinya sendiri.

```bash
dotnet add package Gravicode.Biotech.Net.Core            # + Native
dotnet add package Gravicode.Biotech.Net.Sequences
dotnet add package Gravicode.Biotech.Net.Alignments
dotnet add package Gravicode.Biotech.Net.Genomics
dotnet add package Gravicode.Biotech.Net.Transcriptomics
dotnet add package Gravicode.Biotech.Net.Proteomics
dotnet add package Gravicode.Biotech.Net.StructuralBiology
dotnet add package Gravicode.Biotech.Net.Chemistry
dotnet add package Gravicode.Biotech.Net.DrugDiscovery
dotnet add package Gravicode.Biotech.Net.SyntheticBiology
dotnet add package Gravicode.Biotech.Net.Data
dotnet add package Gravicode.Biotech.Net.Gpu
dotnet add package Gravicode.Biotech.Net.AI
dotnet add package Gravicode.Biotech.Net.Visualization
dotnet add package Gravicode.Biotech.Net.Notebooks
dotnet add package Gravicode.Biotech.Net.Workflows
```

Di notebook:

```text
#r "nuget: Gravicode.Biotech.Net.Notebooks, 0.2.0"
```

Biner native ikut dalam `Gravicode.Biotech.Net.Native` di `runtimes/{rid}/native/`, sehingga `dotnet publish` atau `dotnet run` biasa otomatis memakai biner yang tepat.

## Dari sumber

```bash
git clone https://github.com/DotNetVibeCoderz/Vibe_Bio.git
cd Vibe_Bio/Biotech.Net

dotnet build Biotech.Net.slnx                  # ikut membangun crate Rust bila cargo ada di PATH
dotnet build -p:SkipRustBuild=true             # pakai biner native yang sudah ada
dotnet test Biotech.Net.slnx                   # seluruh proyek uji
dotnet test tests/Biotech.Net.Core.Tests       # satu modul saja

cd native/biotech-native && cargo test --release   # menguji kernel Rust tersendiri
```

Solusi memakai format `.slnx`, bukan `.sln`. Versi paket dipusatkan di `Directory.Packages.props` dan metadata NuGet bersama di `Directory.Build.props`.

## Memeriksa kernel mana yang aktif

```csharp
using Biotech.Net.Core;
using Biotech.Net.Native;

Console.WriteLine(BiotechInfo.Banner);
Console.WriteLine($"native tersedia: {BioKernels.IsNativeAvailable} ({BioKernels.NativeVersion})");
Console.WriteLine($"backend dipakai: {BioKernels.Backend}");

BioKernels.ForceManaged = true;                 // mengunci jalur managed (uji, reproduksibilitas)
Console.WriteLine($"backend sekarang: {BioKernels.Backend}");
BioKernels.ForceManaged = false;
```

`BIOTECH_NET_DISABLE_NATIVE=1` di lingkungan melakukan hal yang sama tanpa mengubah kode. Kedua backend menghasilkan nilai yang sama; hanya kecepatannya berbeda ([benchmark](benchmarks.md)).

Bila pustaka native gagal dimuat, tidak ada yang dilempar — implementasi managed langsung mengambil alih, dan string versinya menjadi `managed`.

## Kebutuhan tambahan saat runtime

| Fitur | Yang diperlukan |
|---|---|
| Render plot di Linux | paket native asset SkiaSharp, sudah menjadi dependensi `Gravicode.Biotech.Net.Visualization` |
| Tampilan 3-D protein dan molekul | peramban; halaman memuat 3Dmol.js dari CDN |
| Model AI (`Biotech.Net.AI`) | akses jaringan saat pertama kali; bobot disimpan di `BioHubOptions.CacheDirectory`. Repositori bergerbang perlu token Hugging Face |
| Kopilot LLM | salah satu dari `ANTHROPIC_API_KEY` (plus `ANTHROPIC_WORKSPACE_ID` bila kunci tidak terikat workspace), `AZURE_OPENAI_*`, `OPENAI_API_KEY`, atau `DEEPSEEK_API_KEY` |

Untuk mesin tanpa jaringan: setel `BioHubOptions.Offline = true` setelah menyalin cache model, dan bangun dengan `-p:SkipRustBuild=true` agar tidak ada unduhan crate. Crate Rust-nya sendiri tanpa dependensi dan bisa dibangun offline.

## CI

`.github/workflows/biotechnet-ci.yml` membangun kernel Rust untuk ketiga RID, menjalankan semua uji di Windows dan Linux, mengeksekusi notebook, lalu memaketkan pustaka beserta seluruh biner native. Mendorong tag `biotechnet-v<versi>` akan menerbitkan paket ke nuget.org dari monorepo `Vibe_Bio`.

---
Selanjutnya: [Pendahuluan](getting-started.md) · [Tutorial](tutorials/index.md)
