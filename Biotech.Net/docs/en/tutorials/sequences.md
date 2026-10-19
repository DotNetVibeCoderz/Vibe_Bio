# Tutorial: Sequences and file I/O

Package: `Gravicode.Biotech.Net.Sequences` (depends on `Gravicode.Biotech.Net.Core`).

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

This tutorial covers the core `Sequence` type, the readers and writers for FASTA, FASTQ, GenBank and EMBL, and the everyday analyses: ORFs, motifs, restriction digests, k-mers, melting temperature and codon usage.

## 1. Creating sequences

`Sequence` stores ASCII bytes and keeps the case you give it. Comparisons and kernels ignore case. Reverse complement, translation and base counting run on the Rust kernels when the native library is present, and on the managed fallback when it is not.

```csharp
using Biotech.Net.Core;

var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", id: "demo");
Console.WriteLine(dna.Length);                 // 39
Console.WriteLine(dna.GcContent());            // 0.564...
Console.WriteLine(dna.ReverseComplement());    // CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
Console.WriteLine(dna.Transcribe());           // AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Console.WriteLine(dna.Translate());            // MAIVMGR*KGAR*
Console.WriteLine(dna.Translate(toStop: true)); // MAIVMGR

// Other genetic codes: NCBI tables 1, 2, 4 and 11
var mito = dna.Translate(GeneticCode.VertebrateMitochondrial);

// Alphabet detection and subsequences (0-based start, length)
var unknown = Sequence.Parse("MKTAYIAKQR");
Console.WriteLine(unknown.Type);               // Protein
var cds = dna.Subsequence(0, 21);
Console.WriteLine(cds.IndexOf("GGG"));         // 14
```

## 2. Reading and writing files

Every reader streams records with `IEnumerable<T>`, so a large file is never held in memory all at once. Paths ending in `.gz` are decompressed transparently.

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences.IO;

string fasta = """
>seq1 first record
ATGCGTACGTTAGC
>seq2 second record
GGGCCCAAATTT
""";
foreach (var record in FastaReader.Parse(fasta))
    Console.WriteLine($"{record.Id}: {record.Length} bp, '{record.Description}'");

// Write FASTA (60 columns by default) to any TextWriter
var sw = new StringWriter();
FastaWriter.Write(sw, FastaReader.Parse(fasta));

// FASTQ: the Phred offset is detected automatically
string fastq = "@read1\nACGTACGTAC\n+\nIIIIIHHHGG\n";
foreach (var read in FastqReader.Parse(fastq))
{
    Console.WriteLine($"{read.Id} mean Q{read.MeanQuality():F1}");
    var trimmed = read.TrimByQuality(minQuality: 20);
}

// Format-agnostic access for files on disk:
//   SequenceFile.Read("genome.gbk")          // GenBank, EMBL, FASTA... detected from content
//   SequenceFile.Write("out.fasta", records)
//   using var fai = IndexedFastaReader.Open("hg38.fa");   // samtools-style faidx random access
//   var region = fai.Fetch("chr1:1,000,001-1,000,500");
```

GenBank and EMBL records keep their features (`Sequence.Features`) and header fields (`Sequence.Annotations`, keys in `RecordKeys`):

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences.IO;

var plasmid = Sequence.Dna("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA", "pDEMO");
plasmid.Features.Add(new SequenceFeature("CDS", FeatureLocation.Span(0, 66)).With("gene", "demoA"));
string gbk = GenBankWriter.ToString(plasmid);
var roundTrip = GenBankReader.Parse(gbk)[0];
Console.WriteLine(roundTrip.Features[0]);      // CDS 1..66 demoA
```

## 3. Analyses

```csharp
using Biotech.Net.Core;
using Biotech.Net.Sequences;

var genome = RandomSequence.Dna(5_000, seed: 1, gcContent: 0.52, id: "toy");

// Open reading frames on both strands
var orfs = OrfFinder.Find(genome, new OrfOptions { MinLength = 150 });
Console.WriteLine($"{orfs.Count} ORFs; longest: {OrfFinder.Longest(genome)}");

// IUPAC motif search on both strands, optionally with mismatches
var sites = MotifSearch.FindIupac(genome, "TATAWAWR", bothStrands: true);

// Restriction digest with the built-in enzyme table
var digest = RestrictionDigest.Digest(genome, circular: false, "EcoRI", "BamHI", "HindIII");
Console.WriteLine(string.Join(", ", digest.FragmentSizes));

// k-mer profile (dense counting runs on the k-mer kernel)
var kmers = KmerProfile.Count(genome, k: 6);
foreach (var (kmer, count) in kmers.Top(3)) Console.WriteLine($"{kmer} {count}");

// Sliding-window GC content and skew
var gc = SlidingWindow.GcContent(genome, window: 500, step: 250);

// Primer melting temperature (nearest-neighbour, SantaLucia 1998)
var primer = Sequence.Dna("AGCGGATAACAATTTCACACAGG");
Console.WriteLine($"Tm = {MeltingTemperature.NearestNeighbor(primer):F1} °C");

// Codon usage and codon adaptation index
var table = CodonUsageTable.FromSequences(orfs.Select(o => o.Extract(genome)));
Console.WriteLine($"CAI of first ORF: {table.Cai(orfs[0].Extract(genome)):F3}");
```

## 4. Visualising

`Biotech.Net.Visualization` draws sequence logos and dot plots. See the [visualization tutorial](visualization.md).

![Sequence logo](../../images/sequence-logo.png)

## Limitations

- GenBank/EMBL writing covers the common header lines and the feature table; rarely used lines (for example `CONTIG`) are not written back.
- `IndexedFastaReader` reads plain FASTA files only, not BGZF-compressed ones. BAM, CRAM and tabix live in `Biotech.Net.Genomics` (see the [Genomics tutorial](genomics.md)).
- `PdbSequenceReader` extracts sequences only. Use `Biotech.Net.StructuralBiology` for coordinates.

---
Next: [Alignments](alignments.md) · [Tutorial index](index.md)
