# API reference: AI

Package `Gravicode.Biotech.Net.AI` · namespaces `Biotech.Net.AI.Models`, `.Preprocessing`, `.Hub`, `.Inference`, `.Copilot`.

> Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil

## Catalogue (`.Models`)

- `BioModelCatalog`: `BioBert`, `PubMedBert`, `DnaBert6`, `DnaBert2`, `ProtBert`, `ProtBertBfd`, `Esm2T6`, `Esm2T12`, `Esm2T30`, `Esm2T33`, `Esm1v`, `GeneformerV1`, `GeneformerV2`, `NucleotideTransformer500mHuman`, `NucleotideTransformerV2_50m`; `All`, `Get`, `TryGet`, `ForKind`, `Smallest(kind, task)`.
- `BioModelDescriptor(Key, RepoId, DisplayName, InputKind, Tasks, Description)` with `Backend`, `Subfolder`, `MaxTokens`, `KmerSize`, `RequiresSpacedResidues`, `Licence`, `Citation`, `ApproximateParameters`, `KnownToRunOnCpu`, `Preprocess(raw)`, `Supports(task)`.
- Enums `BioModelTask`, `BioSequenceKind`, `BioModelBackend` (`HfTransformer`, `Esm`, `Geneformer`, `NotSupported`).

## Preprocessing (`.Preprocessing`)

- `SequencePreprocessor.CleanDna`, `CleanProtein`, `IsValidDna`, `IsValidProtein`, `Kmers`, `KmerTokenize`, `FromOverlappingKmers`, `DnaInput`, `SpaceResidues`, `EsmInput`, `Chunk(sequence, chunkLength, overlap)`; `SequenceChunk(Start, Text)`.
- `GeneformerEncoder(tokens, medians, options)`: `Encode(geneIds, counts)`, `Encode(expression)`, `EncodeDetailed`, `EncodeMatrix`, `PadId`, `MaskId`, `ClsId`, `EosId`, `VocabularySize`; `GeneformerEncoderOptions.V1/V2`, `GeneformerEncoding`.
- `GeneformerDictionaries.CreateEncoder(version)`, `FromHub`, `LoadTokenDictionary`, `LoadGeneMedians`, `Paths`; `GeneformerDictionaryVersion`.

## Hub (`.Hub`)

`BioModelHub.Configure(options)`, `Options`, `CacheDirectory`, `CachedDirectory(repoId, revision)`, `IsCached`, `EnsureModel`, `DownloadFile`; `BioHubOptions(CacheDirectory, Token, Offline, Timeout)`; `BioModelUnavailableException(RepoId)`.

## Inference (`.Inference`)

| Type | Key members |
|---|---|
| `BioEmbedder` | `Load(keyOrRepoId)`, `Embed(input)`, `EmbedBatch(inputs, progress, cancellationToken)`, `TokenEmbeddings`, `Pool`, `Dimension`, `Pooling`, `WindowOverlap`, `MaxDegreeOfParallelism`, `CosineSimilarity(a, b)` |
| `IBioEncoder` | implemented by `HfTransformerEncoder` (BERT-family) and `EsmEncoder`; `Encode(modelInput)` returns `TokenStates(Hidden, Tokens, HiddenSize, IsSpecial)` |
| `EsmModel` | `Load(repoIdOrDirectory, revision)`, `Open(directory)`, `Forward(ids)`, `Hidden(sequence, out ids)`, `LogitsAt`, `LogSoftmax`, `Config`, `Tokenizer`, `HasLanguageModelHead`, `HiddenSize` |
| `EsmTokenizer` | `FromDirectory`, `Encode`, `Decode`, `Id`, `Token`, `IsSpecial`, `MaskToken`, special-token ids |
| `BioMaskedLanguageModel` | `Load(keyOrRepoId)`, `FillMask(input, topK)`, `PredictResidue(protein, position, topK)`, `MaskToken`, `Esm` |
| `EsmVariantScorer` | `SiteLogProbabilities`, `Score(wildType, mutation)`, `Score(wildType, "A1B:C2D")`, `ScoreMany`, `SaturationScan(wildType, from, to)`; `ProteinMutation`, `VariantEffectScore`, `MutationalScan` |
| `GeneformerModel` | `Load(descriptor)`, `Open(directory)`, `Hidden(tokenIds)`, `EmbedCell(tokenIds, specialIds)`, `Dimension` |

`BioEmbedder.Load` and `BioMaskedLanguageModel.Load` throw `NotSupportedException` for `BioModelBackend.NotSupported` models (DNABERT-2, Nucleotide Transformer v2), with the reason in the message.

## Copilot (`.Copilot`)

- `IChatClient` with `AnthropicChatClient(apiKey, model)` (`WorkspaceId`, `BaseUri`, `Effort`, default model `claude-opus-5-5`) and `OpenAiCompatibleChatClient` (`ForOpenAI`, `ForAzureOpenAI`, `ForDeepSeek`).
- `ChatClientFactory.FromEnvironment()`, `DetectProvider(apiKey)`; `ChatProvider`.
- `ChatMessage`, `ChatRequest`, `ChatCompletion`, `ChatToolCall`, `ChatToolDefinition`, `ChatRole`, `ChatClientException(StatusCode)`.
- `BioToolbox(includeBuiltIns)`: built-in tools `reverse_complement`, `transcribe`, `translate`, `gc_content`, `protein_properties`, `find_orfs`; `Add(name, description, schema, handler)`, `Invoke`, `Definitions`, `Schema(...)`.
- `BioCopilot(client, toolbox, options)`: `AskAsync(question, cancellationToken)` returning `CopilotAnswer(Text, ToolCalls, Rounds, StopReason)`; `History`, `Reset()`; `BioCopilotOptions(SystemPrompt, MaxRounds, MaxTokens)`.

API keys are read from environment variables at call time and never logged or persisted.

---
[API reference index](index.md) · [AI tutorial](../tutorials/ai.md)
