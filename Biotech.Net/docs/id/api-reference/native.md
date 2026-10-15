# Referensi API: Native

Paket `Gravicode.Biotech.Net.Native` · namespace `Biotech.Net.Native` · tanpa dependensi.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Lapisan terbawah: jembatan P/Invoke ke crate Rust `biotech-native` plus implementasi C# untuk setiap titik masuk. `BioKernels` memilih salah satunya saat pemanggilan.

## BioKernels

| Anggota | Kegunaan |
|---|---|
| `IsNativeAvailable`, `NativeVersion`, `Backend` | apakah pustaka Rust termuat, string versinya, dan backend yang dipakai |
| `ForceManaged { get; set; }` | mengunci jalur managed (juga lewat `BIOTECH_NET_DISABLE_NATIVE=1`) |
| `ReverseComplement(src, dst)`, `Complement(src, dst)` | sadar IUPAC, mempertahankan besar-kecil huruf, menulis ke buffer pemanggil |
| `Translate(nucleotides, codonTable, protein)` | tabel kodon 64 entri dalam urutan kernel (A=0, C=1, G=2, T=3; indeks `16·b1 + 4·b2 + b3`) |
| `CountNucleotides(sequence)` | mengembalikan `NucleotideCounts` (A/C/G/T/N/lainnya, `GcFraction`, `GcSkew`) |
| `ByteHistogram(data, counts)` | histogram 256 bin |
| `CountKmers(sequence, k, counts)` | penghitungan padat untuk k ≤ 12; `counts` berisi `4^k` entri |
| `HammingDistance(a, b)` | untuk sekuens sepanjang sama |
| `WelchTTest(matrix, rows, cols, groupA, groupB, t, p, meanDifference)` | uji Welch per baris pada matriks padat row-major |
| `BenjaminiHochberg(pValues, qValues)` | BH-FDR |
| `ProteinMass(sequence, monoisotopic)`, `NetCharge(sequence, pH)`, `IsoelectricPoint(sequence)` | mesin protein |
| `Hydropathy(sequence, window, output)` | jendela geser Kyte–Doolittle |

## Kernel penjajaran dan struktur

| Anggota | Kegunaan |
|---|---|
| `AlignmentMode` | `Global`, `Local`, `SemiGlobal` |
| `MatchMismatchMatrix(match, mismatch)` | membangun matriks penilaian kernel 128×128 |
| `Align(a, b, matrix, gapOpen, gapExtend, mode, band)` | penjajaran Gotoh gap afin; mengembalikan `KernelAlignment` (skor, rentang, string operasi) |
| `Align(..., Span<byte> operations)` | kelebihan beban tanpa alokasi, mengembalikan `AlignmentSummary` |
| `LocalScore(query, target, matrix, gapOpen, gapExtend)` | skor Smith–Waterman saja |
| `LocalScoreBatch(query, targets, offsets, matrix, gapOpen, gapExtend, scores)` | satu kueri terhadap banyak target yang digabung |
| `KabschRmsd(target, mobile, rotation)` | RMSD dari array `xyz` datar, dengan keluaran rotasi 3×3 opsional |

Penalti gap pada API kernel berupa biaya non-negatif (`PairwiseAligner` di Alignments justru memakai *skor* gap negatif seperti biasa). `SubstitutionMatrixLength` bernilai `128 * 128`.

## Catatan

- Kesetaraan kedua backend diuji di `tests/Biotech.Net.Native.Tests/ParityTests.cs`.
- Kernel jalur kritis baru memerlukan tiga bagian: versi `ManagedKernels`, versi Rust, dan dispatch `BioKernels`.
- Ukuran buffer divalidasi sebelum dispatch; tujuan yang terlalu kecil melempar `ArgumentException`, dan matriks penjajaran yang tidak muat melempar `OutOfMemoryException` (pakailah band atau kernel skor saja).

---
[Indeks referensi API](index.md) · [Benchmark](../benchmarks.md)
