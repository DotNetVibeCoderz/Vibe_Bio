# Installation

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Requirements

- **.NET 10 SDK** (10.0.100 or newer; `global.json` pins `latestFeature`).
- **Rust toolchain** — optional, only to rebuild the native kernels yourself. Without it everything still works through the managed fallback.
- Supported runtime identifiers for the native kernels: `win-x64`, `linux-x64`, `osx-arm64`.

## From NuGet

Add just the modules you need; each one pulls in what it depends on.

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

In a notebook:

```text
#r "nuget: Gravicode.Biotech.Net.Notebooks, 0.2.0"
```

The native binaries ship inside `Gravicode.Biotech.Net.Native` under `runtimes/{rid}/native/`, so a normal `dotnet publish` or `dotnet run` picks up the right one automatically.

## From source

```bash
git clone https://github.com/DotNetVibeCoderz/Vibe_Bio.git
cd Vibe_Bio/Biotech.Net

dotnet build Biotech.Net.slnx                  # builds the Rust crate too when cargo is on PATH
dotnet build -p:SkipRustBuild=true             # reuse whatever native binary already exists
dotnet test Biotech.Net.slnx                   # every test project
dotnet test tests/Biotech.Net.Core.Tests       # one module

cd native/biotech-native && cargo test --release   # the Rust kernels on their own
```

The solution uses the `.slnx` format, not `.sln`. Package versions are centralised in `Directory.Packages.props` and the shared NuGet metadata in `Directory.Build.props`.

## Checking which kernels are active

```csharp
using Biotech.Net.Core;
using Biotech.Net.Native;

Console.WriteLine(BiotechInfo.Banner);
Console.WriteLine($"native available: {BioKernels.IsNativeAvailable} ({BioKernels.NativeVersion})");
Console.WriteLine($"backend in use:   {BioKernels.Backend}");

BioKernels.ForceManaged = true;                 // pin the managed path (tests, reproducibility)
Console.WriteLine($"backend now:      {BioKernels.Backend}");
BioKernels.ForceManaged = false;
```

`BIOTECH_NET_DISABLE_NATIVE=1` in the environment does the same thing without touching code. Both back-ends produce the same results; only speed differs ([benchmarks](benchmarks.md)).

If the native library cannot load, nothing throws — the managed implementation takes over. `Biotech.Net.Native` reports the reason internally, and the version string becomes `managed`.

## Optional runtime extras

| Feature | What it needs |
|---|---|
| Plot rendering on Linux | the SkiaSharp native asset package, already a dependency of `Gravicode.Biotech.Net.Visualization` |
| 3-D protein and molecule views | a browser; the page loads 3Dmol.js from a CDN |
| AI models (`Biotech.Net.AI`) | network access on first use; weights are cached under `BioHubOptions.CacheDirectory`. Gated repositories need a Hugging Face token |
| The LLM copilot | one of `ANTHROPIC_API_KEY` (plus `ANTHROPIC_WORKSPACE_ID` when the key is not workspace-scoped), `AZURE_OPENAI_*`, `OPENAI_API_KEY` or `DEEPSEEK_API_KEY` |

For air-gapped machines: set `BioHubOptions.Offline = true` after copying a model cache across, and build with `-p:SkipRustBuild=true` so no crate download is attempted. The Rust crate itself has no dependencies and builds offline.

## CI

`.github/workflows/biotechnet-ci.yml` builds the Rust kernels for all three RIDs, runs every test on Windows and Linux, executes the notebooks and packs the libraries with all native binaries inside. Pushing a `biotechnet-v<version>` tag publishes the packages to nuget.org from the `Vibe_Bio` monorepo.

---
Next: [Getting started](getting-started.md) · [Tutorials](tutorials/index.md)
