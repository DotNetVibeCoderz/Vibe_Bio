using Biotech.Net.AI.Preprocessing;

namespace Biotech.Net.AI.Models;

/// <summary>
/// A task a biological language model checkpoint can be used for.
/// </summary>
public enum BioModelTask
{
    /// <summary>Producing a fixed-length vector representation of an input.</summary>
    Embedding,

    /// <summary>Predicting masked tokens (residues, nucleotides, k-mers or words).</summary>
    MaskedLanguageModel,

    /// <summary>Assigning one or more labels to a whole input.</summary>
    TextClassification,

    /// <summary>Assigning labels to individual tokens, for example named-entity recognition.</summary>
    TokenClassification,

    /// <summary>Extracting an answer span from a context given a question.</summary>
    QuestionAnswering,
}

/// <summary>
/// The kind of input a biological model consumes.
/// </summary>
public enum BioSequenceKind
{
    /// <summary>Natural-language biomedical text (abstracts, clinical notes).</summary>
    BiomedicalText,

    /// <summary>Nucleotide sequences (DNA; RNA is accepted and converted U to T where relevant).</summary>
    Dna,

    /// <summary>Amino-acid sequences in one-letter code.</summary>
    Protein,

    /// <summary>Single-cell transcriptomes, encoded as ranked gene tokens (Geneformer).</summary>
    SingleCellExpression,
}

/// <summary>
/// The runtime that executes a checkpoint in Biotech.Net.
/// </summary>
public enum BioModelBackend
{
    /// <summary>A BERT-family encoder run by HF.Net's <c>TransformerModel</c>.</summary>
    HfTransformer,

    /// <summary>An ESM-family encoder (ESM-1b, ESM-1v, ESM-2, Nucleotide Transformer v1) run by HF.Net's ESM encoder.</summary>
    Esm,

    /// <summary>A Geneformer BERT checkpoint fed with gene-token ids (<c>GeneformerModel</c>).</summary>
    Geneformer,

    /// <summary>Listed for reference only: the architecture needs custom code that does not run here yet.</summary>
    NotSupported,
}

/// <summary>
/// Describes a pretrained biological model on the Hugging Face Hub: where it lives, what it consumes,
/// how its input must be prepared and what it is known to support.
/// </summary>
/// <param name="Key">Short, stable, case-insensitive catalogue key, for example <c>Esm2T6</c>.</param>
/// <param name="RepoId">Hugging Face repository id, for example <c>facebook/esm2_t6_8M_UR50D</c>.</param>
/// <param name="DisplayName">Human-readable model name.</param>
/// <param name="InputKind">The kind of sequence the model consumes.</param>
/// <param name="Tasks">Tasks the checkpoint is suitable for.</param>
/// <param name="Description">Free-text description, including known limitations.</param>
public sealed record BioModelDescriptor(
    string Key,
    string RepoId,
    string DisplayName,
    BioSequenceKind InputKind,
    IReadOnlyList<BioModelTask> Tasks,
    string Description)
{
    /// <summary>
    /// The k-mer length for k-mer tokenised DNA models (6 for DNABERT-6); <see langword="null"/> when the
    /// model does not use k-mer input.
    /// </summary>
    public int? KmerSize { get; init; }

    /// <summary>Which runtime executes the checkpoint. Defaults to <see cref="BioModelBackend.HfTransformer"/>.</summary>
    public BioModelBackend Backend { get; init; } = BioModelBackend.HfTransformer;

    /// <summary>Sub-folder of the repository holding the checkpoint (Geneformer variants), or <see langword="null"/> for the root.</summary>
    public string? Subfolder { get; init; }

    /// <summary>Maximum number of tokens per input, including special tokens. Defaults to 512.</summary>
    public int MaxTokens { get; init; } = 512;

    /// <summary>
    /// <see langword="true"/> when protein residues must be separated by single spaces (ProtBERT / ProtT5 convention).
    /// </summary>
    public bool RequiresSpacedResidues { get; init; }

    /// <summary>Licence of the checkpoint as stated by its publisher, when known.</summary>
    public string? Licence { get; init; }

    /// <summary>Reference citation for the model.</summary>
    public string? Citation { get; init; }

    /// <summary>Approximate number of parameters, used to pick the smallest model for a kind.</summary>
    public long ApproximateParameters { get; init; }

    /// <summary>
    /// <see langword="true"/> when the checkpoint is known to load and run through HF.Net on a CPU in reasonable time.
    /// </summary>
    public bool KnownToRunOnCpu { get; init; }

    /// <summary>
    /// Converts a raw sequence into the string form this model expects:
    /// biomedical text is returned unchanged, k-mer DNA models receive space-separated overlapping k-mers,
    /// other DNA models receive the cleaned upper-case sequence, ProtBERT-style models receive space-separated
    /// residues and ESM-style models receive the raw upper-case sequence.
    /// </summary>
    /// <param name="rawSequence">The raw input.</param>
    /// <returns>The model-ready input string.</returns>
    /// <exception cref="NotSupportedException">
    /// Thrown for <see cref="BioSequenceKind.SingleCellExpression"/>: Geneformer consumes gene-token ids, which
    /// are produced by <see cref="GeneformerEncoder"/>, not by a string transform.
    /// </exception>
    public string Preprocess(string rawSequence)
    {
        ArgumentNullException.ThrowIfNull(rawSequence);
        return InputKind switch
        {
            BioSequenceKind.BiomedicalText => rawSequence,
            BioSequenceKind.Dna => KmerSize is int k
                ? SequencePreprocessor.KmerTokenize(rawSequence, k)
                : SequencePreprocessor.DnaInput(rawSequence),
            BioSequenceKind.Protein => RequiresSpacedResidues
                ? SequencePreprocessor.SpaceResidues(rawSequence)
                : SequencePreprocessor.EsmInput(rawSequence),
            BioSequenceKind.SingleCellExpression => throw new NotSupportedException(
                $"Model '{Key}' consumes gene-token ids, not text. Use GeneformerEncoder.Encode to build its rank-value encoding."),
            _ => throw new NotSupportedException($"Unknown input kind {InputKind}."),
        };
    }

    /// <summary>Returns whether the model is listed as supporting <paramref name="task"/>.</summary>
    /// <param name="task">The task to test.</param>
    /// <returns><see langword="true"/> when the task is supported.</returns>
    public bool Supports(BioModelTask task)
    {
        for (int i = 0; i < Tasks.Count; i++)
        {
            if (Tasks[i] == task)
            {
                return true;
            }
        }

        return false;
    }
}
