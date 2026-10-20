# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Biotech.Net is an open-source biotechnology and bioinformatics platform for .NET 10 ("Bioconductor + Biopython for .NET"), with Rust-accelerated kernels.

- `BiotechNet.txt` is the project brief (in Indonesian) and `solution-design.md` is the architecture spec. Both are the specification of record.
- `PLAN.md` holds the roadmap and `Progress.md` tracks development.

## Commands

```bash
dotnet build Biotech.Net.slnx                       # whole solution (.slnx format, not .sln)
dotnet test Biotech.Net.slnx                        # every test project
dotnet test tests/Biotech.Net.Core.Tests            # one module
dotnet test tests/Biotech.Net.Core.Tests --filter "FullyQualifiedName~SequenceTests"
dotnet build -p:SkipRustBuild=true                  # skip cargo; use whatever native binary already exists
cd native/biotech-native && cargo test --release    # Rust kernels
```

`Biotech.Net.Native.csproj` runs `cargo build --release` before building when cargo is on PATH. Without cargo, every kernel falls back to managed C#.

## Layout and dependency direction

```
native/biotech-native        Rust crate (C ABI) -> biotech_native.dll / .so / .dylib
src/Biotech.Net.Native       BioKernels: P/Invoke bridge + ManagedKernels fallback. No dependencies.
src/Biotech.Net.Core         Sequence, Alphabet, GeneticCode, SequenceRead, GenomicInterval, SequenceFeature,
                             Gene/Transcript/Chromosome/Genome, Protein, Variant/Genotype, Dataset. Depends on Native.
src/Biotech.Net.<Module>     Sequences, Alignments, Genomics, Transcriptomics, Proteomics, StructuralBiology,
                             Chemistry, DrugDiscovery, SyntheticBiology, AI, Visualization, Notebooks
tests/Biotech.Net.<Module>.Tests   one xunit.v3 project per module
```

The lowest layer is Native, not Core. Core calls `BioKernels` for reverse complement, translation, composition and protein properties. A new hot-path kernel belongs in Native and needs three parts: a `ManagedKernels` version, a Rust version, and the `BioKernels` dispatch.

Each `.csproj` lists its module-to-module references. Keep the graph acyclic; Visualization and Notebooks sit at the top.

## Conventions

- **Package ids** are `Gravicode.<ProjectName>`, e.g. `Gravicode.Biotech.Net.Core`. Shared NuGet metadata and the icon live in `Directory.Build.props`.
- **Package versions** are managed centrally in `Directory.Packages.props`.
- **Coordinates:** `GenomicInterval` and `FeatureLocation` are 0-based half-open. `Variant.Position` is 1-based, like VCF. `ToString()` prints the 1-based form.
- **Sequence storage:**
  - `Sequence` stores ASCII bytes and keeps the case it was given.
  - Comparisons and kernels ignore case.
  - Parsers should build a `byte[]` and call `Sequence.FromBuffer` (no copy).
- **Codon tables** use kernel order, A=0, C=1, G=2, T=3, with index `16*b1 + 4*b2 + b3`. This is not NCBI's TCAG order; `GeneticCode` converts between them.
- **Bilingual docs:** `docs/en/` and `docs/id/` mirror each other, and both must change together.
- **Credits:** docs and apps must show "Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil" (`BiotechInfo.CreditsId`).
- **UI:** the Gallery and sample apps use the `frontend-design` skill (`.claude/skills/`).

## External credentials (outside the repo, never copy their contents into it)

- NuGet API key: `C:\Users\mifma\Documents\CodeSandbox\PackageCredentials.txt`
- LLM API key for real-model tests: `C:\Users\mifma\Documents\CodeSandbox\testkey.txt`
- Hugging Face token: `C:\Users\mifma\Documents\CodeSandbox\HFToken.txt`

This solution lives in the monorepo `C:\experiment\VibeCoding\Vibe_Bio` (GitHub `DotNetVibeCoderz/Vibe_Bio`) as the `Biotech.Net/` subfolder. Workflows run from the monorepo root `.github/workflows/` (source copies in `build/workflows/`): `biotechnet-ci.yml`, `biotechnet-publish.yml` (push a `biotechnet-v<version>` tag; secret `NUGET_API_KEY`) and `biotechnet-webgallery.yml` (GitHub Pages). nuget.org versions cannot be overwritten, so push and wait for green CI before tagging.
