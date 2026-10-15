# Tutorial: AI

Paket: `Gravicode.Biotech.Net.AI` (di atas HF.Net: `Gravicode.HFNet.GraviHub`, `GraviTokenizers`, `GraviTransformers`)

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

Katalog model fondasi biologi, praproses offline untuk masukan DNA, protein, dan sel tunggal, embedding, prediksi masked language model, penilaian varian ESM, pengodean rank-value Geneformer, serta kopilot LLM dengan perkakas bioinformatika.

Semua yang ada di bagian 1 dan 2 berjalan offline dan tidak memerlukan unduhan model.

## 1. Katalog model

```csharp
using Biotech.Net.AI.Models;

foreach (var m in BioModelCatalog.All)
    Console.WriteLine($"{m.Key,-32} {m.InputKind,-22} {m.Backend,-14} {m.ApproximateParameters / 1_000_000.0,6:F0}M  {m.RepoId}");

var esm = BioModelCatalog.Esm2T6;
Console.WriteLine($"{esm.DisplayName}: {esm.Description}");
Console.WriteLine($"model embedding protein terkecil: {BioModelCatalog.Smallest(BioSequenceKind.Protein)?.Key}");
```

`BioModelBackend.NotSupported` menandai model yang sengaja ditolak Biotech.Net, dengan alasannya tertulis di `Description`. Lihat matriks dukungan di bawah.

## 2. Praproses

```csharp
using Biotech.Net.AI.Preprocessing;

Console.WriteLine(SequencePreprocessor.CleanDna("acgt nnryACGT"));         // ACGTNNNNACGT (U->T, ambigu->N; karakter non-IUPAC melempar galat)
Console.WriteLine(SequencePreprocessor.KmerTokenize("ACGTACGTAA", k: 6));  // 6-mer tumpang tindih untuk DNABERT
Console.WriteLine(SequencePreprocessor.SpaceResidues("MKTAYIAKQR"));       // masukan ProtBERT
Console.WriteLine(SequencePreprocessor.EsmInput("mktayiakqr"));            // MKTAYIAKQR

foreach (var chunk in SequencePreprocessor.Chunk(new string('A', 2_500), chunkLength: 1_000, overlap: 100))
    Console.WriteLine($"[{chunk.Start}, {chunk.End})");
```

`BioModelDescriptor.Preprocess(raw)` menerapkan apa pun yang diharapkan model bersangkutan (k-mer, residu berspasi, atau residu polos).

## 3. Embedding dan prediksi masked language model

Bagian ini mengunduh bobot dari Hugging Face saat pertama kali dipakai lalu menyimpannya di cache. Pakailah token untuk repositori bergerbang, atau arahkan cache ke direktori yang sudah Anda isi sendiri.

```csharp
using Biotech.Net.AI.Hub;
using Biotech.Net.AI.Inference;

BioModelHub.Configure(new BioHubOptions
{
    CacheDirectory = Path.Combine(Path.GetTempPath(), "biotech-models"),
    Token = Environment.GetEnvironmentVariable("HF_TOKEN"),
    Offline = false,
});

// ESM-2 t6 (8 juta parameter) adalah model protein terkecil yang tetap berguna
using var embedder = BioEmbedder.Load("Esm2T6");
float[] v1 = embedder.Embed("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ");
float[] v2 = embedder.Embed("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEAQ");
Console.WriteLine($"dimensi {embedder.Dimension}, kosinus {BioEmbedder.CosineSimilarity(v1, v2):F4}");

// Prediksi masked dan prediksi residu tunggal
using var mlm = BioMaskedLanguageModel.Load("Esm2T6");
foreach (var p in mlm.PredictResidue("MKTAYIAKQRQISFVKSHF", position: 5, topK: 3))
    Console.WriteLine($"{p.Token} {p.Score:F3}");

// Efek varian zero-shot (masked marginals)
var scorer = new EsmVariantScorer(mlm.Esm!);
Console.WriteLine(scorer.Score("MKTAYIAKQRQISFVKSHF", "Y5F"));   // tipe liar harus cocok dengan sekuensnya
var scan = scorer.SaturationScan("MKTAYIAKQRQISFVKSHF", from: 1, to: 5);
```

Model sel tunggal memakai pengodean rank-value, bukan tokenizer:

```csharp
using Biotech.Net.AI.Preprocessing;

// Encoder dari kamus Geneformer di Hub (gc104M untuk V2)
// var encoder = GeneformerDictionaries.CreateEncoder(GeneformerDictionaryVersion.V2Gc104M);
// int[] tokens = encoder.Encode(new Dictionary<string, double> { ["ENSG00000141510"] = 12, ["ENSG00000012048"] = 3 });
// using var geneformer = GeneformerModel.Load(BioModelCatalog.GeneformerV2);
// float[] cell = geneformer.EmbedCell(tokens);
Console.WriteLine(GeneformerEncoderOptions.V2.MaxLength);   // 4096
```

## 4. Kopilot LLM

`BioCopilot` menjalankan loop pemanggilan tool: model menjawab pertanyaan dengan memanggil perkakas bioinformatika di `BioToolbox` (komplemen balik, transkripsi, translasi, kandungan GC, sifat protein, pencarian ORF), dan Anda bisa mendaftarkan perkakas sendiri.

```csharp
using Biotech.Net.AI.Copilot;

// Klien dipilih dari lingkungan: ANTHROPIC_API_KEY, AZURE_OPENAI_*, OPENAI_API_KEY, atau DEEPSEEK_API_KEY
var client = ChatClientFactory.FromEnvironment();
if (client is null)
{
    Console.WriteLine("Tidak ada kunci LLM yang dikonfigurasi; contoh kopilot dilewati.");
}
else
{
    var copilot = new BioCopilot(client);
    var answer = await copilot.AskAsync("Translate ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG and tell me the protein length.");
    Console.WriteLine(answer.Text);
    foreach (var call in answer.ToolCalls) Console.WriteLine($"  tool {call.Name} -> {call.Result}");
}

// Konstruksi eksplisit, dengan workspace id Anthropic bila kunci tidak terikat workspace
// var anthropic = new AnthropicChatClient(key, "claude-opus-5-5") { WorkspaceId = "wrkspc_..." };
// var openai    = OpenAiCompatibleChatClient.ForOpenAI(key, "gpt-5-mini");

// Menambahkan tool sendiri
var toolbox = new BioToolbox()
    .Add("melting_temperature", "Nearest-neighbour Tm of a DNA oligo in °C",
        BioToolbox.Schema(("sequence", "string", "DNA oligo", true)),
        args => Biotech.Net.Sequences.MeltingTemperature
            .NearestNeighbor(Biotech.Net.Core.Sequence.Dna(args.GetProperty("sequence").GetString()!))
            .ToString("F2"));
Console.WriteLine($"{toolbox.Tools.Count} tool terdaftar");
```

## Matriks dukungan model

| Model | Backend di sini | Status |
|---|---|---|
| ESM-2 (t6 8M, t12 35M, t30 150M, t33 650M), ESM-1v | runtime ESM managed bawaan | didukung — embedding, masked LM, penilaian varian |
| BioBERT, PubMedBERT | transformer HF.Net | didukung — embedding teks dan masked LM |
| ProtBERT, ProtBERT-BFD | transformer HF.Net | didukung; masukan harus residu berspasi (ditangani `Preprocess`) |
| DNABERT (k = 3…6) | transformer HF.Net | didukung; masukan berupa k-mer tumpang tindih |
| Geneformer V1 10M / V2 104M | Geneformer | didukung lewat `GeneformerModel` plus pengodean rank-value (bukan lewat `BioEmbedder`) |
| Nucleotide Transformer 500M human-ref | runtime ESM managed | dapat dimuat, tetapi belum diverifikasi secara numerik terhadap implementasi Python |
| **DNABERT-2 117M** | — | **tidak didukung.** Kode `bert_layers.py`-nya (BPE + ALiBi) bukan BERT standar, sehingga memuatnya akan menghasilkan angka yang salah; Biotech.Net menolaknya. |
| **Nucleotide Transformer v2 50M** | — | **tidak didukung.** Memerlukan posisi rotari dengan feed-forward bergerbang (SwiGLU) dari kode modelling khusus, yang tidak diimplementasikan runtime ESM bawaan. |

Kedua model yang tidak didukung tetap ada di katalog supaya perkakas bisa menjelaskan *mengapa* keduanya tidak bisa dijalankan: `BioEmbedder.Load` dan `BioMaskedLanguageModel.Load` melempar `NotSupportedException` beserta alasannya.

## Batasan

- Inferensi berjalan di CPU sebagai kode managed. ESM-2 t6/t12 nyaman dipakai; checkpoint 650M lambat dan boros memori. Kernel GPU ada di [roadmap](../roadmap.md).
- Model diunduh dari Hugging Face saat pertama kali dipakai. Repositori bergerbang memerlukan `BioHubOptions.Token`, dan `Offline = true` membuat cache yang kosong melempar `BioModelUnavailableException` alih-alih mengunduh.
- Klien Anthropic memerlukan `anthropic-workspace-id` untuk kunci API yang tidak terikat workspace: setel `AnthropicChatClient.WorkspaceId` atau variabel lingkungan `ANTHROPIC_WORKSPACE_ID`, kalau tidak permintaannya ditolak.
- `EsmVariantScorer` memerlukan checkpoint dengan kepala masked-LM dan melaporkan rasio log masked-marginal. Nilai itu berkorelasi dengan kebugaran eksperimental, tetapi bukan probabilitas patogenisitas yang terkalibrasi.
- Kunci API dibaca dari variabel lingkungan saat pemanggilan dan tidak pernah dicatat atau disimpan oleh pustaka.

---
Sebelumnya: [Biologi sintetik](synthetic-biology.md) · Selanjutnya: [Visualisasi](visualization.md) · [Indeks tutorial](index.md)
