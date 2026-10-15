using Biotech.Net.Alignments;
using Biotech.Net.Alignments.Phylogenetics;
using Biotech.Net.Chemistry;
using Biotech.Net.Core;
using Biotech.Net.StructuralBiology;
using Biotech.Net.Transcriptomics;
using Biotech.Net.Visualization;
using Microsoft.DotNet.Interactive.Formatting;

namespace Biotech.Net.Notebooks;

/// <summary>
/// Registers <c>text/html</c> formatters for Biotech.Net types with .NET Interactive, so that a cell ending in a
/// <see cref="Sequence"/>, <see cref="Molecule"/>, <see cref="Structure"/>, plot, … renders richly in Polyglot
/// Notebooks and Jupyter. Loaded automatically by the package's <c>extension.dib</c>; call
/// <see cref="Register"/> yourself when referencing the project directly.
/// </summary>
/// <example>
/// <code>
/// #r "nuget: Gravicode.Biotech.Net.Notebooks"
/// Biotech.Net.Notebooks.BiotechFormatters.Register();
/// Sequence.Dna("ATGGCGTACGT")   // coloured, wrapped, with length and GC
/// </code>
/// </example>
public static class BiotechFormatters
{
    private static readonly object Gate = new();

    /// <summary>Whether <see cref="Register"/> has run.</summary>
    public static bool IsRegistered { get; private set; }

    /// <summary>Types that get an HTML formatter.</summary>
    public static IReadOnlyList<Type> FormattedTypes { get; } =
    [
        typeof(Sequence), typeof(SequenceRead), typeof(Protein), typeof(Variant), typeof(GenomicInterval), typeof(Gene),
        typeof(Dataset), typeof(DifferentialExpressionResult), typeof(Molecule), typeof(Structure), typeof(BioPlot),
        typeof(PairwiseAlignment), typeof(MultipleSequenceAlignment), typeof(PhylogeneticTree),
    ];

    /// <summary>Registers all formatters (idempotent).</summary>
    public static void Register()
    {
        lock (Gate)
        {
            if (IsRegistered) return;
            const string html = HtmlFormatter.MimeType;
            Formatter.Register<Sequence>((s, w) => w.Write(BiotechHtml.Sequence(s)), html);
            Formatter.Register<SequenceRead>((r, w) => w.Write(BiotechHtml.Read(r)), html);
            Formatter.Register<Protein>((p, w) => w.Write(BiotechHtml.Protein(p)), html);
            Formatter.Register<Variant>((v, w) => w.Write(BiotechHtml.Variant(v)), html);
            Formatter.Register<GenomicInterval>((i, w) => w.Write(BiotechHtml.Interval(i)), html);
            Formatter.Register<Gene>((g, w) => w.Write(BiotechHtml.Gene(g)), html);
            Formatter.Register<Dataset>((d, w) => w.Write(BiotechHtml.Dataset(d)), html);
            Formatter.Register<DifferentialExpressionResult>((d, w) => w.Write(BiotechHtml.DifferentialExpression(d)), html);
            Formatter.Register<Molecule>((m, w) => w.Write(BiotechHtml.Molecule(m)), html);
            Formatter.Register<Structure>((s, w) => w.Write(BiotechHtml.Structure(s)), html);
            Formatter.Register<BioPlot>((p, w) => w.Write(BiotechHtml.Plot(p)), html);
            Formatter.Register<PairwiseAlignment>((a, w) => w.Write(AlignmentHtml.Pairwise(a)), html);
            Formatter.Register<MultipleSequenceAlignment>((a, w) => w.Write(AlignmentHtml.Multiple(a)), html);
            Formatter.Register<PhylogeneticTree>((t, w) => w.Write(AlignmentHtml.Tree(t)), html);
            foreach (var t in FormattedTypes) Formatter.SetPreferredMimeTypesFor(t, html);
            IsRegistered = true;
        }
    }

    /// <summary>Removes all custom formatters (resets .NET Interactive formatting to its defaults).</summary>
    public static void Reset()
    {
        lock (Gate)
        {
            Formatter.ResetToDefault();
            IsRegistered = false;
        }
    }

    /// <summary>The banner shown when the extension loads.</summary>
    public static string BannerHtml()
        => $"{BiotechHtml.Styles}<div class=\"btn-card\"><div class=\"btn-h\"><span class=\"btn-t\">{BiotechHtml.Enc(BiotechInfo.Product)} {BiotechHtml.Enc(BiotechInfo.Version)}</span>"
           + $"<span class=\"btn-k\">notebook formatters loaded</span></div><div class=\"muted\">{BiotechHtml.Enc(BiotechInfo.Tagline)}</div>{BiotechHtml.Footer()}</div>";
}
