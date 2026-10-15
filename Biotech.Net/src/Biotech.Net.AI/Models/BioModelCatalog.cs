namespace Biotech.Net.AI.Models;

/// <summary>
/// The model zoo: well-known biological language-model checkpoints on the Hugging Face Hub, with the input
/// preparation and runtime each one needs.
/// </summary>
/// <remarks>
/// Every repository id here was checked against the Hub (October 2026). <see cref="BioModelDescriptor.KnownToRunOnCpu"/>
/// is <see langword="true"/> only for checkpoints whose architecture a Biotech.Net runtime implements;
/// <see cref="BioModelBackend.NotSupported"/> entries are listed so callers get a clear refusal rather than a wrong answer.
/// </remarks>
public static class BioModelCatalog
{
    private static readonly BioModelTask[] TextTasks = [BioModelTask.Embedding, BioModelTask.MaskedLanguageModel];
    private static readonly BioModelTask[] MlmTasks = [BioModelTask.Embedding, BioModelTask.MaskedLanguageModel];

    /// <summary>BioBERT v1.2 (cased), BERT-base further pretrained on PubMed abstracts.</summary>
    public static BioModelDescriptor BioBert { get; } = new(
        "BioBert", "dmis-lab/biobert-base-cased-v1.2", "BioBERT v1.2", BioSequenceKind.BiomedicalText, TextTasks,
        "BERT-base (cased vocabulary) pretrained on PubMed abstracts and PMC full text. The repository ships no tokenizer_config.json, " +
        "so — as in Python transformers' default — the tokenizer lower-cases input.")
    {
        ApproximateParameters = 110_000_000, KnownToRunOnCpu = true,
        Citation = "Lee et al., BioBERT, Bioinformatics 36(4):1234–1240 (2020).",
    };

    /// <summary>PubMedBERT (now BiomedBERT), trained from scratch on PubMed abstracts and PMC full text.</summary>
    public static BioModelDescriptor PubMedBert { get; } = new(
        "PubMedBert", "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext", "PubMedBERT (BiomedBERT) abstract+fulltext",
        BioSequenceKind.BiomedicalText, TextTasks,
        "BERT-base with an in-domain uncased vocabulary, pretrained from scratch on PubMed and PMC. Formerly " +
        "microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext.")
    {
        ApproximateParameters = 110_000_000, KnownToRunOnCpu = true, Licence = "MIT",
        Citation = "Gu et al., Domain-specific language model pretraining for biomedical NLP, ACM HEALTH 3(1) (2021).",
    };

    /// <summary>DNABERT with overlapping 6-mers.</summary>
    public static BioModelDescriptor DnaBert6 { get; } = new(
        "DnaBert6", "zhihan1996/DNA_bert_6", "DNABERT (6-mer)", BioSequenceKind.Dna, MlmTasks,
        "BERT-base pretrained on the human genome with overlapping 6-mer tokens (input: space-separated k-mers). " +
        "Siblings DNA_bert_3/4/5 use k = 3/4/5.")
    {
        KmerSize = 6, ApproximateParameters = 89_000_000, KnownToRunOnCpu = true,
        Citation = "Ji et al., DNABERT, Bioinformatics 37(15):2112–2120 (2021).",
    };

    /// <summary>DNABERT-2 (BPE tokens, ALiBi attention) — listed, not runnable here.</summary>
    public static BioModelDescriptor DnaBert2 { get; } = new(
        "DnaBert2", "zhihan1996/DNABERT-2-117M", "DNABERT-2 117M", BioSequenceKind.Dna, MlmTasks,
        "Multi-species genome model with byte-pair-encoded DNA and ALiBi attention. Its custom modelling code (bert_layers.py) " +
        "is not a standard BERT, so HF.Net's encoder would load it incorrectly; Biotech.Net refuses it until a dedicated runtime exists.")
    {
        Backend = BioModelBackend.NotSupported, ApproximateParameters = 117_000_000,
        Citation = "Zhou et al., DNABERT-2, ICLR 2024.",
    };

    /// <summary>ProtBERT (UniRef100).</summary>
    public static BioModelDescriptor ProtBert { get; } = new(
        "ProtBert", "Rostlab/prot_bert", "ProtBERT", BioSequenceKind.Protein, MlmTasks,
        "30-layer BERT pretrained on UniRef100; residues must be space-separated, with U/Z/O/B mapped to X. About 1.7 GB of weights.")
    {
        RequiresSpacedResidues = true, ApproximateParameters = 420_000_000, KnownToRunOnCpu = true,
        Citation = "Elnaggar et al., ProtTrans, IEEE TPAMI 44(10):7112–7127 (2022).",
    };

    /// <summary>ProtBERT-BFD.</summary>
    public static BioModelDescriptor ProtBertBfd { get; } = new(
        "ProtBertBfd", "Rostlab/prot_bert_bfd", "ProtBERT-BFD", BioSequenceKind.Protein, MlmTasks,
        "ProtBERT pretrained on the Big Fantastic Database (2.1 billion sequences); same input convention as ProtBERT.")
    {
        RequiresSpacedResidues = true, ApproximateParameters = 420_000_000, KnownToRunOnCpu = true,
        Citation = "Elnaggar et al., ProtTrans, IEEE TPAMI 44(10):7112–7127 (2022).",
    };

    /// <summary>ESM-2, 6 layers, 8M parameters — the smallest ESM-2, good for tests and quick embeddings.</summary>
    public static BioModelDescriptor Esm2T6 { get; } = Esm2("Esm2T6", "facebook/esm2_t6_8M_UR50D", "ESM-2 t6 8M", 8_000_000, 320);

    /// <summary>ESM-2, 12 layers, 35M parameters.</summary>
    public static BioModelDescriptor Esm2T12 { get; } = Esm2("Esm2T12", "facebook/esm2_t12_35M_UR50D", "ESM-2 t12 35M", 35_000_000, 480);

    /// <summary>ESM-2, 30 layers, 150M parameters.</summary>
    public static BioModelDescriptor Esm2T30 { get; } = Esm2("Esm2T30", "facebook/esm2_t30_150M_UR50D", "ESM-2 t30 150M", 150_000_000, 640);

    /// <summary>ESM-2, 33 layers, 650M parameters — the size used for most published ESM-2 variant-effect benchmarks.</summary>
    public static BioModelDescriptor Esm2T33 { get; } = Esm2("Esm2T33", "facebook/esm2_t33_650M_UR50D", "ESM-2 t33 650M", 650_000_000, 1280);

    /// <summary>ESM-1v (seed 1), trained specifically for zero-shot variant-effect prediction.</summary>
    public static BioModelDescriptor Esm1v { get; } = new(
        "Esm1v", "facebook/esm1v_t33_650M_UR90S_1", "ESM-1v (seed 1)", BioSequenceKind.Protein,
        [BioModelTask.Embedding, BioModelTask.MaskedLanguageModel],
        "33-layer ESM-1 model trained on UniRef90 for zero-shot mutation-effect prediction (Meier et al. 2021); five seeds _1…_5 " +
        "are published and are usually averaged. Ships pytorch_model.bin only (no safetensors).")
    {
        Backend = BioModelBackend.Esm, MaxTokens = 1024, ApproximateParameters = 650_000_000, KnownToRunOnCpu = true, Licence = "MIT",
        Citation = "Meier et al., Language models enable zero-shot prediction of the effects of mutations on protein function, NeurIPS 2021.",
    };

    /// <summary>Geneformer V1 (10M parameters, Genecorpus-30M, 2,048 genes per cell).</summary>
    public static BioModelDescriptor GeneformerV1 { get; } = new(
        "GeneformerV1", "ctheodoris/Geneformer", "Geneformer V1 10M", BioSequenceKind.SingleCellExpression, MlmTasks,
        "6-layer BERT over rank-value-encoded single-cell transcriptomes (sub-folder Geneformer-V1-10M); use the gc30M dictionaries.")
    {
        Backend = BioModelBackend.Geneformer, Subfolder = "Geneformer-V1-10M", MaxTokens = 2048,
        ApproximateParameters = 10_000_000, KnownToRunOnCpu = true, Licence = "Apache-2.0",
        Citation = "Theodoris et al., Transfer learning enables predictions in network biology, Nature 618:616–624 (2023).",
    };

    /// <summary>Geneformer V2 (104M parameters, Genecorpus-104M, 4,096 tokens with &lt;cls&gt;/&lt;eos&gt;).</summary>
    public static BioModelDescriptor GeneformerV2 { get; } = new(
        "GeneformerV2", "ctheodoris/Geneformer", "Geneformer V2 104M", BioSequenceKind.SingleCellExpression, MlmTasks,
        "12-layer Geneformer V2 (sub-folder Geneformer-V2-104M); use the gc104M dictionaries. The repository root holds the 316M V2 model.")
    {
        Backend = BioModelBackend.Geneformer, Subfolder = "Geneformer-V2-104M", MaxTokens = 4096,
        ApproximateParameters = 104_000_000, KnownToRunOnCpu = true, Licence = "Apache-2.0",
        Citation = "Theodoris et al., Nature 618:616–624 (2023); Chen et al., Geneformer V2 (2024).",
    };

    /// <summary>Nucleotide Transformer v1, 500M parameters, human reference genome.</summary>
    public static BioModelDescriptor NucleotideTransformer500mHuman { get; } = new(
        "NucleotideTransformer500mHumanRef", "InstaDeepAI/nucleotide-transformer-500m-human-ref", "Nucleotide Transformer 500M (human ref)",
        BioSequenceKind.Dna, MlmTasks,
        "ESM-architecture model (absolute positions) over non-overlapping 6-mers of the human reference genome; input is raw DNA, " +
        "which the ESM tokenizer splits into 6-mers. Ships pytorch_model.bin only. Not verified numerically against Python here.")
    {
        Backend = BioModelBackend.Esm, MaxTokens = 1000, ApproximateParameters = 500_000_000, KnownToRunOnCpu = true,
        Licence = "CC-BY-NC-SA-4.0",
        Citation = "Dalla-Torre et al., The Nucleotide Transformer, Nature Methods 22:287–297 (2025).",
    };

    /// <summary>Nucleotide Transformer v2, 50M parameters, multi-species — listed, not runnable here.</summary>
    public static BioModelDescriptor NucleotideTransformerV2_50m { get; } = new(
        "NucleotideTransformerV2_50m", "InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", "Nucleotide Transformer v2 50M (multi-species)",
        BioSequenceKind.Dna, MlmTasks,
        "Uses rotary positions and a gated (SwiGLU) feed-forward without biases from custom modelling code; the built-in ESM runtime does not implement that block.")
    {
        Backend = BioModelBackend.NotSupported, MaxTokens = 2048, ApproximateParameters = 50_000_000, Licence = "CC-BY-NC-SA-4.0",
        Citation = "Dalla-Torre et al., Nature Methods 22:287–297 (2025).",
    };

    /// <summary>Every catalogued model.</summary>
    public static IReadOnlyList<BioModelDescriptor> All { get; } =
    [
        BioBert, PubMedBert, DnaBert6, DnaBert2, ProtBert, ProtBertBfd, Esm2T6, Esm2T12, Esm2T30, Esm2T33, Esm1v,
        GeneformerV1, GeneformerV2, NucleotideTransformer500mHuman, NucleotideTransformerV2_50m,
    ];

    /// <summary>Looks a model up by its catalogue key or Hub repository id (case-insensitive).</summary>
    /// <exception cref="KeyNotFoundException">No such model.</exception>
    public static BioModelDescriptor Get(string keyOrRepoId)
        => TryGet(keyOrRepoId, out var d) ? d
            : throw new KeyNotFoundException($"No catalogued model '{keyOrRepoId}'. Known keys: {string.Join(", ", All.Select(m => m.Key))}.");

    /// <summary>Looks a model up by its catalogue key or Hub repository id (case-insensitive).</summary>
    public static bool TryGet(string keyOrRepoId, out BioModelDescriptor descriptor)
    {
        ArgumentNullException.ThrowIfNull(keyOrRepoId);
        descriptor = All.FirstOrDefault(m => m.Key.Equals(keyOrRepoId, StringComparison.OrdinalIgnoreCase))
                     ?? All.FirstOrDefault(m => m.RepoId.Equals(keyOrRepoId, StringComparison.OrdinalIgnoreCase))!;
        return descriptor is not null;
    }

    /// <summary>Models consuming <paramref name="kind"/>, smallest first.</summary>
    public static IReadOnlyList<BioModelDescriptor> ForKind(BioSequenceKind kind)
        => [.. All.Where(m => m.InputKind == kind).OrderBy(m => m.ApproximateParameters)];

    /// <summary>The smallest runnable model for <paramref name="kind"/> that supports <paramref name="task"/>.</summary>
    public static BioModelDescriptor? Smallest(BioSequenceKind kind, BioModelTask task = BioModelTask.Embedding)
        => ForKind(kind).FirstOrDefault(m => m.Backend != BioModelBackend.NotSupported && m.Supports(task));

    private static BioModelDescriptor Esm2(string key, string repo, string name, long parameters, int hidden) => new(
        key, repo, name, BioSequenceKind.Protein, [BioModelTask.Embedding, BioModelTask.MaskedLanguageModel],
        $"ESM-2 protein language model trained on UniRef50 (hidden size {hidden}). Input is the raw one-letter sequence; " +
        "embeddings are mean-pooled over residues; the LM head scores substitutions (masked marginals).")
    {
        Backend = BioModelBackend.Esm, MaxTokens = 1024, ApproximateParameters = parameters, KnownToRunOnCpu = true, Licence = "MIT",
        Citation = "Lin et al., Evolutionary-scale prediction of atomic-level protein structure with a language model, Science 379:1123–1130 (2023).",
    };
}
