namespace Biotech.Net.WebGallery.Services;

/// <summary>Hands a page's live-code snippet to the Playground ("Open in Playground").</summary>
public sealed class PlaygroundState
{
    /// <summary>Code waiting to be opened in the Playground, with the page it came from.</summary>
    public (string Title, string Code)? Pending { get; set; }

    /// <summary>The editor's text, kept while the user visits other pages.</summary>
    public string? Draft { get; set; }

    /// <summary>The preset or page title the draft came from.</summary>
    public string? DraftTitle { get; set; }
}
