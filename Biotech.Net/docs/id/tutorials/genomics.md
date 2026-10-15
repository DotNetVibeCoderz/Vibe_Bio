# Tutorial: Genomika

Paket: `Gravicode.Biotech.Net.Genomics`

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Berkas varian (VCF), prediksi efek varian, pemetaan bacaan dan pemanggilan SNV, perakitan de Bruijn, pencarian gen prokariot, dan genetika populasi.

Koordinat: `GenomicInterval` berbasis 0 dan setengah terbuka, sedangkan `Variant.Position` berbasis 1 seperti VCF.

## 1. Membaca dan menulis VCF

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Vcf;

string vcf = """
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	101	rs1	A	G	50	PASS	DP=30	GT	0/1	1/1	0/0
chr1	250	.	C	T	12	LowQual	DP=8	GT	0/0	0/1	./.
chr1	400	rs3	AT	A	70	PASS	DP=41	GT	0/1	0/1	1/1
""".Replace("\r", "");

var (header, records) = VcfReader.Parse(vcf);
Console.WriteLine($"{header.Samples.Count} sampel, {records.Count} rekaman");
foreach (var r in records.PassOnly())
    Console.WriteLine($"{r.Variant} {r.Variant.Type} AF={r.Variant.AlternateAlleleFrequency():F2} GT(S1)={r.GetGenotype(0)}");

// Menulis kembali (bgzip/tabix belum didukung; path ".gz" ditulis sebagai gzip biasa)
var sw = new StringWriter();
using (var writer = new VcfWriter(sw, header, leaveOpen: true)) writer.Write(records);
```

Untuk berkas di disk, `VcfReader.Open(path)` mengalirkan rekaman dan `reader.Query(GenomicInterval)` menyaring berdasarkan wilayah.

## 2. Prediksi efek varian

`VariantEffectPredictor` memberikan konsekuensi Sequence Ontology dan notasi HGVS memakai gen serta transkrip dari sebuah `Genome`. Gen bisa berasal dari GFF3 (`GffReader.LoadGenes`) atau dibangun dalam kode:

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Annotation;

// Gen satu ekson: ATG AAA CGC ... TAA di chr1 posisi 10..39 (berbasis 0)
var chrSeq = Sequence.Dna("GGGGGGGGGG" + "ATGAAACGCATTAGCACCACCATTACCTAA" + "GGGGGGGGGG", "chr1");
var genome = new Genome("toy").Add(Chromosome.FromSequence(chrSeq));
var gene = new Gene("G1", new GenomicInterval("chr1", 10, 40, Strand.Forward), "TOY1");
gene.Transcripts.Add(new Transcript("TOY1-201", "G1", [new GenomicInterval("chr1", 10, 40, Strand.Forward)])
    { CodingRegion = new GenomicInterval("chr1", 10, 40, Strand.Forward) });
genome.Add(gene);

var vep = new VariantEffectPredictor(genome);
// Posisi 17 (berbasis 1) adalah basa pertama kodon ke-3 (CGC, Arg); C>T menghasilkan TGC (Cys)
foreach (var effect in vep.Predict("chr1", 17, "C", "T"))
    Console.WriteLine($"{effect.ConsequenceTerms} {effect.Impact} {effect.HgvsC} {effect.HgvsP}");
```

## 3. Memetakan bacaan dan memanggil SNV

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Calling;
using Biotech.Net.Genomics.Mapping;
using Biotech.Net.Sequences;

var reference = RandomSequence.Dna(2_000, seed: 5, id: "ref");
// Simulasikan perubahan homozigot pada posisi 1000 (berbasis 0)
var sample = reference.ToString().ToCharArray();
sample[1000] = sample[1000] == 'A' ? 'G' : 'A';
var sampleSeq = new string(sample);

var rng = new Random(9);
var reads = Enumerable.Range(0, 400).Select(i =>
{
    int start = rng.Next(0, sampleSeq.Length - 100);
    return SequenceRead.FromAscii($"r{i}", sampleSeq.Substring(start, 100), new string('I', 100));
}).ToList();

var mapper = new ReadMapper([reference]);
var alignments = mapper.Map(reads).ToList();
Console.WriteLine($"terpetakan {alignments.Count(a => a.IsMapped)}/{alignments.Count}");

var calls = SnvCaller.Call([reference], alignments);
foreach (var call in calls) Console.WriteLine(call);
```

`SamWriter` menulis penjajaran sebagai SAM, dan `SamReader` membaca SAM dari perkakas lain.

## 4. Perakitan dan pencarian gen

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Assembly;
using Biotech.Net.Genomics.GeneFinding;
using Biotech.Net.Sequences;

var target = RandomSequence.Dna(3_000, seed: 11, id: "target");
var rng = new Random(4);
var reads = Enumerable.Range(0, 600)
    .Select(i => target.Subsequence(rng.Next(0, target.Length - 120), 120).WithId($"r{i}"))
    .ToList();

var assembly = DeBruijnAssembler.Assemble(reads, new AssemblyOptions { K = 31, MinKmerCount = 2 });
Console.WriteLine(assembly.Stats);

// Pencari gen prokariot yang melatih dirinya sendiri (model Markov + skor RBS)
var finder = new ProkaryoticGeneFinder();
foreach (var g in finder.Predict(target).Take(5)) Console.WriteLine(g);
```

## 5. Genetika populasi

```csharp
using Biotech.Net.Core;
using Biotech.Net.Genomics.Population;

var v = new Variant("chr1", 100, "A", "G");
foreach (var gt in new[] { "0/0", "0/1", "0/1", "1/1", "0/0", "0/1" }) v.Genotypes.Add(Genotype.Parse(gt));

var counts = GenotypeCounts.From(v);
Console.WriteLine($"AF {counts.AlternateFrequency:F2}, Hobs {counts.ObservedHeterozygosity:F2}, HWE p {counts.HardyWeinbergP():F3}");

// Tajima's D dari site frequency spectrum, Fst Hudson dari hitungan alel
long[] sfs = SiteFrequencySpectrum.Unfolded([1, 1, 2, 5, 1, 3, 1], sampleSize: 10);
Console.WriteLine($"Tajima's D {SiteFrequencySpectrum.TajimasD(sfs):F3}");
double fst = Fst.Hudson([(Ac1: 2, An1: 20, Ac2: 15, An2: 20), (Ac1: 5, An1: 20, Ac2: 6, An2: 20)]);
```

Genome browser di `Biotech.Net.Visualization` menggambar gen, cakupan, dan varian dalam satu tampilan:

![Genome browser](../../images/genome-browser.png)

## Batasan

- BAM/CRAM, BGZF, dan tabix belum diimplementasikan (lihat [roadmap](../roadmap.md)). Konversikan dulu ke SAM atau VCF.
- `ReadMapper` ditujukan untuk bacaan pendek terhadap referensi kecil–menengah yang muat di memori, misalnya genom bakteri dan amplikon. Bukan pengganti BWA atau minimap2 untuk genom manusia utuh.
- `SnvCaller` hanya memanggil SNV diploid; tidak memanggil indel dan tidak melakukan perakitan ulang lokal.
- Perakit de Bruijn cocok untuk bacaan pendek bergalat rendah pada genom kecil. Tidak ada scaffolding dan tidak memanfaatkan informasi paired-end.

---
Sebelumnya: [Penjajaran](alignments.md) · Selanjutnya: [Transkriptomika](transcriptomics.md) · [Indeks tutorial](index.md)
