using Avalonia.Controls;
using Biotech.Net.Visualization;
using CommunityToolkit.Mvvm.ComponentModel;

namespace Biotech.Net.Gallery.Shell;

/// <summary>
/// One module of the Gallery: navigation metadata, the result view (built lazily in code), the C# snippet that
/// reproduces the current result, the sequence that drives the header chromatogram, and what can be exported.
/// </summary>
public abstract partial class GalleryPage : ObservableObject
{
    private Control? _view;

    /// <summary>Navigation and header title.</summary>
    public abstract string Title { get; }

    /// <summary>The Biotech.Net packages the page exercises (header eyebrow).</summary>
    public abstract string Modules { get; }

    /// <summary>One sentence under the title: what the page lets you do.</summary>
    public abstract string Lede { get; }

    /// <summary>Dye channel (A, C, G or T) that colours the page's navigation mark.</summary>
    public abstract char Channel { get; }

    /// <summary>File-name slug used by the screenshot runner (gallery-&lt;slug&gt;.png).</summary>
    public virtual string Slug => new string(Title.ToLowerInvariant().Select(c => char.IsLetterOrDigit(c) ? c : '-').ToArray()).Trim('-');

    /// <summary>Whether the page is the one on screen (drives the rail highlight).</summary>
    [ObservableProperty] private bool _isActive;

    /// <summary>C# that produces the current result (Live code tab).</summary>
    [ObservableProperty] private string _code = "";

    /// <summary>DNA drawn by the header chromatogram ribbon.</summary>
    [ObservableProperty] private string _traceSequence = "";

    /// <summary>Where the trace comes from, shown beside the ribbon.</summary>
    [ObservableProperty] private string _traceLabel = "";

    /// <summary>Short status line ("computed in 3.1 ms", errors).</summary>
    [ObservableProperty] private string _status = "";

    /// <summary>The page body, built on first use.</summary>
    public Control View => _view ??= Build();

    /// <summary>Whether <see cref="View"/> has been built.</summary>
    public bool IsBuilt => _view is not null;

    /// <summary>Builds the page body.</summary>
    protected abstract Control Build();

    /// <summary>The library plot behind the page's main figure, exportable as PNG and SVG; null when none.</summary>
    public virtual BioPlot? ExportPlot => null;

    /// <summary>The visual saved by "Save PNG" when <see cref="ExportPlot"/> is null.</summary>
    public virtual Control? ExportVisual => _view;

    /// <summary>Plain-text summary of the current result for "Copy result".</summary>
    public virtual string ResultText => Status;

    /// <summary>Lets the screenshot runner wait for asynchronous content.</summary>
    public virtual Task WhenReadyAsync() => Task.CompletedTask;
}
