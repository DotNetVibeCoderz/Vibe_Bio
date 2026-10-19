# API reference

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

An overview of the key types and methods of each module — what exists and where to find it, not a generated dump of every member. For runnable examples see the [tutorials](../tutorials/index.md); for exact signatures use IntelliSense or the XML documentation that ships in every package.

| Module | Namespace root | Page |
|---|---|---|
| Native | `Biotech.Net.Native` | [native.md](native.md) |
| Core | `Biotech.Net.Core` | [core.md](core.md) |
| Workflows | `Biotech.Net.Workflows` | [workflows.md](workflows.md) |
| Sequences | `Biotech.Net.Sequences` | [sequences.md](sequences.md) |
| Alignments | `Biotech.Net.Alignments` | [alignments.md](alignments.md) |
| Genomics | `Biotech.Net.Genomics` | [genomics.md](genomics.md) |
| Transcriptomics | `Biotech.Net.Transcriptomics` | [transcriptomics.md](transcriptomics.md) |
| Proteomics | `Biotech.Net.Proteomics` | [proteomics.md](proteomics.md) |
| Structural biology | `Biotech.Net.StructuralBiology` | [structural-biology.md](structural-biology.md) |
| Chemistry | `Biotech.Net.Chemistry` | [chemistry.md](chemistry.md) |
| Drug discovery | `Biotech.Net.DrugDiscovery` | [drug-discovery.md](drug-discovery.md) |
| Synthetic biology | `Biotech.Net.SyntheticBiology` | [synthetic-biology.md](synthetic-biology.md) |
| Data | `Biotech.Net.Data` | [data.md](data.md) |
| GPU | `Biotech.Net.Gpu` | [gpu.md](gpu.md) |
| AI | `Biotech.Net.AI` | [ai.md](ai.md) |
| Visualization | `Biotech.Net.Visualization` | [visualization.md](visualization.md) |
| Notebooks | `Biotech.Net.Notebooks` | [notebooks.md](notebooks.md) |

## Conventions that apply everywhere

- **Coordinates.** `GenomicInterval` and `FeatureLocation` are 0-based half-open; `Variant.Position` is 1-based like VCF.
- **Spans.** Hot paths take `ReadOnlySpan<byte>`/`Span<T>` and write into caller-provided buffers.
- **Streaming.** File readers return `IEnumerable<T>` and never materialise a whole file.
- **Options objects.** Long parameter lists are `record` option types with `init` properties and sensible defaults (`OrfOptions`, `DESeqOptions`, `SearchParameters`, `MapperOptions`, …). Most take `null` to mean "defaults".
- **Case.** Sequence data keeps the case it was given; comparisons and kernels ignore it.

---
[Documentation index](../index.md) · Bahasa Indonesia: [referensi-api](../../id/api-reference/index.md)
