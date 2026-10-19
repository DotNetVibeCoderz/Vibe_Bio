# Tutorial: AI

Package: `Gravicode.Biotech.Net.AI` (built on HF.Net: `Gravicode.HFNet.GraviHub`, `GraviTokenizers`, `GraviTransformers`)

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

A catalogue of biological foundation models, offline preprocessing for DNA, protein and single-cell input, embeddings, masked-language-model predictions, ESM variant scoring, Geneformer rank-value encoding, and an LLM copilot with bioinformatics tools.

Everything in sections 1 and 2 runs offline and needs no model download.

## 1. The model catalogue

```csharp
using Biotech.Net.AI.Models;

foreach (var m in BioModelCatalog.All)
    Console.WriteLine($"{m.Key,-32} {m.InputKind,-22} {m.Backend,-14} {m.ApproximateParameters / 1_000_000.0,6:F0}M  {m.RepoId}");

var esm = BioModelCatalog.Esm2T6;
Console.WriteLine($"{esm.DisplayName}: {esm.Description}");
Console.WriteLine($"smallest protein embedding model: {BioModelCatalog.Smallest(BioSequenceKind.Protein)?.Key}");
```

`BioModelBackend.NotSupported` marks a catalogued model that Biotech.Net deliberately refuses to run, with the reason in `Description`. See the support matrix below.

## 2. Preprocessing

```csharp
using Biotech.Net.AI.Preprocessing;

Console.WriteLine(SequencePreprocessor.CleanDna("acgt nnryACGT"));         // ACGTNNNNACGT (U->T, ambiguity->N; a non-IUPAC character throws)
Console.WriteLine(SequencePreprocessor.KmerTokenize("ACGTACGTAA", k: 6));  // overlapping 6-mers for DNABERT
Console.WriteLine(SequencePreprocessor.SpaceResidues("MKTAYIAKQR"));       // ProtBERT input
Console.WriteLine(SequencePreprocessor.EsmInput("mktayiakqr"));            // MKTAYIAKQR

foreach (var chunk in SequencePreprocessor.Chunk(new string('A', 2_500), chunkLength: 1_000, overlap: 100))
    Console.WriteLine($"[{chunk.Start}, {chunk.End})");
```

`BioModelDescriptor.Preprocess(raw)` applies whatever a given model expects (k-mers, spaced residues, plain residues).

## 3. Embeddings and masked-language-model predictions

These download weights from Hugging Face on first use and cache them. Set a token for gated repositories, or point the cache at a directory you populated yourself.

```csharp
using Biotech.Net.AI.Hub;
using Biotech.Net.AI.Inference;

BioModelHub.Configure(new BioHubOptions
{
    CacheDirectory = Path.Combine(Path.GetTempPath(), "biotech-models"),
    Token = Environment.GetEnvironmentVariable("HF_TOKEN"),
    Offline = false,
});

// ESM-2 t6 (8M parameters) is the smallest useful protein model
using var embedder = BioEmbedder.Load("Esm2T6");
float[] v1 = embedder.Embed("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ");
float[] v2 = embedder.Embed("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEAQ");
Console.WriteLine($"dim {embedder.Dimension}, cosine {BioEmbedder.CosineSimilarity(v1, v2):F4}");

// Masked prediction and single-residue predictions
using var mlm = BioMaskedLanguageModel.Load("Esm2T6");
foreach (var p in mlm.PredictResidue("MKTAYIAKQRQISFVKSHF", position: 5, topK: 3))
    Console.WriteLine($"{p.Token} {p.Score:F3}");

// Zero-shot variant effects (masked marginals)
var scorer = new EsmVariantScorer(mlm.Model);
Console.WriteLine(scorer.Score("MKTAYIAKQRQISFVKSHF", "Y5F"));   // wild type must match the sequence
var scan = scorer.SaturationScan("MKTAYIAKQRQISFVKSHF", from: 1, to: 5);
```

Single-cell models use rank-value encoding instead of a tokenizer:

```csharp
using Biotech.Net.AI.Preprocessing;

// Encoder from the Geneformer dictionaries on the Hub (gc104M for V2)
// var encoder = GeneformerDictionaries.CreateEncoder(GeneformerDictionaryVersion.V2Gc104M);
// int[] tokens = encoder.Encode(new Dictionary<string, double> { ["ENSG00000141510"] = 12, ["ENSG00000012048"] = 3 });
// using var geneformer = GeneformerModel.Load(BioModelCatalog.GeneformerV2);
// float[] cell = geneformer.EmbedCell(tokens);
Console.WriteLine(GeneformerEncoderOptions.V2.MaxLength);   // 4096
```

## 4. The LLM copilot

`BioCopilot` runs a tool-calling loop: the model answers questions by calling the bioinformatics tools in `BioToolbox` (reverse complement, transcription, translation, GC content, protein properties, ORF finding) and you can register your own.

```csharp
using Biotech.Net.AI.Copilot;

// Pick a client from the environment: ANTHROPIC_API_KEY, AZURE_OPENAI_*, OPENAI_API_KEY or DEEPSEEK_API_KEY
var client = ChatClientFactory.FromEnvironment();
if (client is null)
{
    Console.WriteLine("No LLM key configured; skipping the copilot example.");
}
else
{
    var copilot = new BioCopilot(client);
    var answer = await copilot.AskAsync("Translate ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG and tell me the protein length.");
    Console.WriteLine(answer.Text);
    foreach (var call in answer.ToolCalls) Console.WriteLine($"  tool {call.Name} -> {call.Result}");
}

// Explicit construction, with an Anthropic workspace id when the key is not workspace-scoped
// var anthropic = new AnthropicChatClient(key, "claude-opus-5-5") { WorkspaceId = "wrkspc_..." };
// var openai    = OpenAiCompatibleChatClient.ForOpenAI(key, "gpt-5-mini");

// Add a custom tool
var toolbox = new BioToolbox()
    .Add("melting_temperature", "Nearest-neighbour Tm of a DNA oligo in °C",
        BioToolbox.Schema(("sequence", "string", "DNA oligo", true)),
        args => Biotech.Net.Sequences.MeltingTemperature
            .NearestNeighbor(Biotech.Net.Core.Sequence.Dna(args.GetProperty("sequence").GetString()!))
            .ToString("F2"));
Console.WriteLine($"{toolbox.Tools.Count} tools registered");
```

## Model support matrix

| Model | Backend here | Status |
|---|---|---|
| ESM-2 (t6 8M, t12 35M, t30 150M, t33 650M), ESM-1v | HF.Net ESM encoder | supported — embeddings, masked LM, variant scoring |
| BioBERT, PubMedBERT | HF.Net transformer | supported — text embeddings and masked LM |
| ProtBERT, ProtBERT-BFD | HF.Net transformer | supported; input must be space-separated residues (handled by `Preprocess`) |
| DNABERT (k = 3…6) | HF.Net transformer | supported; input is overlapping k-mers |
| Geneformer V1 10M / V2 104M | Geneformer | supported through `GeneformerModel` plus rank-value encoding (not `BioEmbedder`) |
| Nucleotide Transformer 500M human-ref | HF.Net ESM encoder | loads, but not verified numerically against the Python implementation |
| **DNABERT-2 117M** | — | **not supported.** Its custom `bert_layers.py` (BPE + ALiBi) is not a standard BERT, so loading it would give wrong numbers; Biotech.Net refuses it. |
| **Nucleotide Transformer v2 50M** | — | **not supported.** Needs rotary positions with a gated (SwiGLU) feed-forward from custom modelling code, which HF.Net's ESM encoder does not implement. |

Both unsupported models stay in the catalogue so that tooling can explain *why* they cannot run: `BioEmbedder.Load` and `BioMaskedLanguageModel.Load` throw `NotSupportedException` with the reason.

## Limitations

- Inference is CPU-only managed code. ESM-2 t6/t12 are comfortable; the 650M checkpoints are slow and memory-hungry. GPU kernels are on the [roadmap](../roadmap.md).
- Models are downloaded from Hugging Face at first use. Gated repositories need `BioHubOptions.Token`, and `Offline = true` makes a missing cache entry throw `BioModelUnavailableException` instead of downloading.
- The Anthropic client needs `anthropic-workspace-id` for API keys that are not scoped to a workspace: set `AnthropicChatClient.WorkspaceId` or the `ANTHROPIC_WORKSPACE_ID` environment variable, otherwise the request is rejected.
- `EsmVariantScorer` needs a checkpoint with a masked-LM head and reports masked-marginal log-ratios. Those correlate with experimental fitness but are not calibrated probabilities of pathogenicity.
- API keys are read from environment variables at call time and are never logged or persisted by the library.

---
Previous: [Synthetic biology](synthetic-biology.md) · Next: [Visualization](visualization.md) · [Tutorial index](index.md)
