# Biotech.Net documentation (English)

**Ultimate Biotechnology & Bioinformatics Platform for .NET** — the Bioconductor + Biopython of .NET 10, with Rust-accelerated kernels.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Bahasa Indonesia: [dokumentasi](../id/index.md)

## Contents

| Page | What it covers |
|---|---|
| [Getting started](getting-started.md) | what Biotech.Net is, a five-minute tour, the conventions |
| [Installation](installation.md) | requirements, packages, building from source, native kernels, CI |
| [Tutorials](tutorials/index.md) | one per module, every snippet compiled and run |
| [API reference](api-reference/index.md) | key types and methods per module |
| [Samples](samples.md) | the Avalonia Gallery and the sample notebooks |
| [Benchmarks](benchmarks.md) | measured Rust-vs-managed numbers, parsers, chemistry, normalisation |
| [Roadmap](roadmap.md) | what is shipped, what is next |

## Tutorials at a glance

[Sequences](tutorials/sequences.md) · [Alignments](tutorials/alignments.md) · [Genomics](tutorials/genomics.md) · [Transcriptomics](tutorials/transcriptomics.md) · [Proteomics](tutorials/proteomics.md) · [Structural biology](tutorials/structural-biology.md) · [Chemistry](tutorials/chemistry.md) · [Drug discovery](tutorials/drug-discovery.md) · [Synthetic biology](tutorials/synthetic-biology.md) · [AI](tutorials/ai.md) · [Visualization](tutorials/visualization.md) · [Notebooks](tutorials/notebooks.md)

## House rules for this documentation

- Every C# snippet is extracted into a scratch project, compiled and executed against the real API before publication. If a snippet prints a value, that value came from a real run.
- Every figure is rendered by `Biotech.Net.Visualization` itself.
- Every page ends with an honest **Limitations** section. Heuristics are labelled as heuristics, demo data sets as demos, and unsupported models as unsupported.
- The English (`docs/en/`) and Indonesian (`docs/id/`) trees mirror each other file for file.

## License

MIT. Copyright © Gravicode Studios.
