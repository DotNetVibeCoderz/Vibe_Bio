# Referensi API: Notebooks

Paket `Gravicode.Biotech.Net.Notebooks` · namespace `Biotech.Net.Notebooks`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## BiotechFormatters

| Anggota | Kegunaan |
|---|---|
| `Register()` | mendaftarkan formatter HTML ke .NET Interactive dan menjadikan HTML tipe MIME pilihan untuk tipe-tipe yang diformat (idempoten) |
| `Reset()` | mengembalikan formatter bawaan |
| `IsRegistered` | apakah `Register()` sudah dijalankan |
| `FormattedTypes` | 14 tipe dengan keluaran kaya: `Sequence`, `SequenceRead`, `Protein`, `Variant`, `GenomicInterval`, `Gene`, `Dataset`, `DifferentialExpressionResult`, `Molecule`, `Structure`, `BioPlot`, `PairwiseAlignment`, `MultipleSequenceAlignment`, `PhylogeneticTree` |
| `BannerHtml()` | banner produk sebagai kartu HTML |

## BiotechHtml

Menghasilkan markup yang sama di luar notebook.

`Styles` (CSS bersama), `Footer()` (baris kredit Gravicode), `Sequence`, `Read`, `Protein`, `Variant`, `Interval`, `Gene`, `Dataset(ds, maxRows, maxColumns)`, `DifferentialExpression(result, top)`, `Molecule`, `Structure(structure, options)`, `Residues(residues, nucleotide, lineWidth, maxResidues)`, `Plot(plot, embedding)`.

Pengaturan: `DefaultPlotEmbedding` (`PlotEmbedding.Auto`, `.Svg`, `.Png`), `LineWidth`, `MaxResidues`.

## AlignmentHtml

`Pairwise(alignment, width)`, `Multiple(msa, width, maxSequences)`, `Tree(tree, width)`, serta `MaxColumns` untuk membatasi penjajaran yang sangat lebar.

## Notebook contoh

`notebooks/*.dib`, dwibahasa (EN/ID): pendahuluan, sekuens dan I/O, genome browser dan plot, ekspresi diferensial RNA-seq, struktur protein, kemoinformatika. Lihat [Contoh](../samples.md).

---
[Indeks referensi API](index.md) · [Tutorial notebook](../tutorials/notebooks.md)
