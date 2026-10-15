# API reference: Native

Package `Gravicode.Biotech.Net.Native` · namespace `Biotech.Net.Native` · no dependencies.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

The lowest layer: a P/Invoke bridge to the Rust crate `biotech-native` plus a managed C# implementation of every entry point. `BioKernels` picks one at call time.

## BioKernels

| Member | Purpose |
|---|---|
| `IsNativeAvailable`, `NativeVersion`, `Backend` | whether the Rust library loaded, its version string, and which back-end calls go to |
| `ForceManaged { get; set; }` | pin the managed path (also `BIOTECH_NET_DISABLE_NATIVE=1`) |
| `ReverseComplement(src, dst)`, `Complement(src, dst)` | IUPAC-aware, case-preserving, into a caller buffer |
| `Translate(nucleotides, codonTable, protein)` | 64-entry codon table in kernel order (A=0, C=1, G=2, T=3; index `16·b1 + 4·b2 + b3`) |
| `CountNucleotides(sequence)` | returns `NucleotideCounts` (A/C/G/T/N/other, `GcFraction`, `GcSkew`) |
| `ByteHistogram(data, counts)` | 256-bin histogram |
| `CountKmers(sequence, k, counts)` | dense counting for k ≤ 12; `counts` has `4^k` entries |
| `HammingDistance(a, b)` | equal-length sequences |
| `WelchTTest(matrix, rows, cols, groupA, groupB, t, p, meanDifference)` | row-wise Welch tests over a dense row-major matrix |
| `BenjaminiHochberg(pValues, qValues)` | BH-FDR |
| `ProteinMass(sequence, monoisotopic)`, `NetCharge(sequence, pH)`, `IsoelectricPoint(sequence)` | protein engine |
| `Hydropathy(sequence, window, output)` | Kyte–Doolittle sliding window |

## Alignment and structure kernels

| Member | Purpose |
|---|---|
| `AlignmentMode` | `Global`, `Local`, `SemiGlobal` |
| `MatchMismatchMatrix(match, mismatch)` | builds the 128×128 kernel scoring matrix |
| `Align(a, b, matrix, gapOpen, gapExtend, mode, band)` | Gotoh affine-gap alignment; returns `KernelAlignment` (score, spans, operation string) |
| `Align(..., Span<byte> operations)` | the allocation-free overload, returns `AlignmentSummary` |
| `LocalScore(query, target, matrix, gapOpen, gapExtend)` | Smith–Waterman score only |
| `LocalScoreBatch(query, targets, offsets, matrix, gapOpen, gapExtend, scores)` | one query against many concatenated targets |
| `KabschRmsd(target, mobile, rotation)` | RMSD of flat `xyz` arrays, optional 3×3 rotation output |

Gap penalties are non-negative costs in the kernel API (`PairwiseAligner` in Alignments exposes the usual negative gap *scores* instead). `SubstitutionMatrixLength` is `128 * 128`.

## Notes

- Both back-ends are tested for parity in `tests/Biotech.Net.Native.Tests/ParityTests.cs`.
- A new hot-path kernel needs three parts: a `ManagedKernels` version, a Rust version, and the `BioKernels` dispatch.
- Buffer sizes are validated before dispatch; a too-small destination throws `ArgumentException`, and an alignment matrix that does not fit throws `OutOfMemoryException` (use a band or a score-only kernel).

---
[API reference index](index.md) · [Benchmarks](../benchmarks.md)
