# Referensi API: AI

Paket `Gravicode.Biotech.Net.AI` · namespace `Biotech.Net.AI.Models`, `.Preprocessing`, `.Hub`, `.Inference`, `.Copilot`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Katalog (`.Models`)

- `BioModelCatalog`: `BioBert`, `PubMedBert`, `DnaBert6`, `DnaBert2`, `ProtBert`, `ProtBertBfd`, `Esm2T6`, `Esm2T12`, `Esm2T30`, `Esm2T33`, `Esm1v`, `GeneformerV1`, `GeneformerV2`, `NucleotideTransformer500mHuman`, `NucleotideTransformerV2_50m`; `All`, `Get`, `TryGet`, `ForKind`, `Smallest(kind, task)`.
- `BioModelDescriptor(Key, RepoId, DisplayName, InputKind, Tasks, Description)` dengan `Backend`, `Subfolder`, `MaxTokens`, `KmerSize`, `RequiresSpacedResidues`, `Licence`, `Citation`, `ApproximateParameters`, `KnownToRunOnCpu`, `Preprocess(raw)`, `Supports(task)`.
- Enum `BioModelTask`, `BioSequenceKind`, `BioModelBackend` (`HfTransformer`, `Esm`, `Geneformer`, `NotSupported`).

## Praproses (`.Preprocessing`)

- `SequencePreprocessor.CleanDna`, `CleanProtein`, `IsValidDna`, `IsValidProtein`, `Kmers`, `KmerTokenize`, `FromOverlappingKmers`, `DnaInput`, `SpaceResidues`, `EsmInput`, `Chunk(sequence, chunkLength, overlap)`; `SequenceChunk(Start, Text)`.
- `GeneformerEncoder(tokens, medians, options)`: `Encode(geneIds, counts)`, `Encode(expression)`, `EncodeDetailed`, `EncodeMatrix`, `PadId`, `MaskId`, `ClsId`, `EosId`, `VocabularySize`; `GeneformerEncoderOptions.V1/V2`, `GeneformerEncoding`.
- `GeneformerDictionaries.CreateEncoder(version)`, `FromHub`, `LoadTokenDictionary`, `LoadGeneMedians`, `Paths`; `GeneformerDictionaryVersion`.

## Hub (`.Hub`)

`BioModelHub.Configure(options)`, `Options`, `CacheDirectory`, `CachedDirectory(repoId, revision)`, `IsCached`, `EnsureModel`, `DownloadFile`; `BioHubOptions(CacheDirectory, Token, Offline, Timeout)`; `BioModelUnavailableException(RepoId)`.

## Inferensi (`.Inference`)

| Tipe | Anggota utama |
|---|---|
| `BioEmbedder` | `Load(keyOrRepoId)`, `Embed(input)`, `EmbedBatch(inputs, progress, cancellationToken)`, `TokenEmbeddings`, `Pool`, `Dimension`, `Pooling`, `WindowOverlap`, `MaxDegreeOfParallelism`, `CosineSimilarity(a, b)` |
| `IBioEncoder` | diimplementasikan `HfTransformerEncoder` (keluarga BERT) dan `EsmEncoder`; `Encode(modelInput)` mengembalikan `TokenStates(Hidden, Tokens, HiddenSize, IsSpecial)` |
| `EsmModel` | `Load(repoIdOrDirectory, revision)`, `Open(directory)`, `Forward(ids)`, `Hidden(sequence, out ids)`, `LogitsAt`, `LogSoftmax`, `Config`, `Tokenizer`, `HasLanguageModelHead`, `HiddenSize` |
| `EsmTokenizer` | `FromDirectory`, `Encode`, `Decode`, `Id`, `Token`, `IsSpecial`, `MaskToken`, id token khusus |
| `BioMaskedLanguageModel` | `Load(keyOrRepoId)`, `FillMask(input, topK)`, `PredictResidue(protein, position, topK)`, `MaskToken`, `Esm` |
| `EsmVariantScorer` | `SiteLogProbabilities`, `Score(wildType, mutation)`, `Score(wildType, "A1B:C2D")`, `ScoreMany`, `SaturationScan(wildType, from, to)`; `ProteinMutation`, `VariantEffectScore`, `MutationalScan` |
| `GeneformerModel` | `Load(descriptor)`, `Open(directory)`, `Hidden(tokenIds)`, `EmbedCell(tokenIds, specialIds)`, `Dimension` |

`BioEmbedder.Load` dan `BioMaskedLanguageModel.Load` melempar `NotSupportedException` untuk model ber-`BioModelBackend.NotSupported` (DNABERT-2, Nucleotide Transformer v2), beserta alasannya di pesan galat.

## Kopilot (`.Copilot`)

- `IChatClient` dengan `AnthropicChatClient(apiKey, model)` (`WorkspaceId`, `BaseUri`, `Effort`, model bawaan `claude-opus-5-5`) dan `OpenAiCompatibleChatClient` (`ForOpenAI`, `ForAzureOpenAI`, `ForDeepSeek`).
- `ChatClientFactory.FromEnvironment()`, `DetectProvider(apiKey)`; `ChatProvider`.
- `ChatMessage`, `ChatRequest`, `ChatCompletion`, `ChatToolCall`, `ChatToolDefinition`, `ChatRole`, `ChatClientException(StatusCode)`.
- `BioToolbox(includeBuiltIns)`: tool bawaan `reverse_complement`, `transcribe`, `translate`, `gc_content`, `protein_properties`, `find_orfs`; `Add(name, description, schema, handler)`, `Invoke`, `Definitions`, `Schema(...)`.
- `BioCopilot(client, toolbox, options)`: `AskAsync(question, cancellationToken)` mengembalikan `CopilotAnswer(Text, ToolCalls, Rounds, StopReason)`; `History`, `Reset()`; `BioCopilotOptions(SystemPrompt, MaxRounds, MaxTokens)`.

Kunci API dibaca dari variabel lingkungan saat pemanggilan dan tidak pernah dicatat atau disimpan.

---
[Indeks referensi API](index.md) · [Tutorial AI](../tutorials/ai.md)
