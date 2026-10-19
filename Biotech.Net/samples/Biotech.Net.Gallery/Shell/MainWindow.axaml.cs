using System.ComponentModel;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input.Platform;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using Avalonia.Platform.Storage;
using AvaloniaEdit.Highlighting;

namespace Biotech.Net.Gallery.Shell;

public partial class MainWindow : Window, IShellServices
{
    private GalleryPage? _watched;

    public MainWindow()
    {
        InitializeComponent();
        ViewModel = new MainViewModel { Services = this };
        DataContext = ViewModel;
        ConfigureEditor();
        ViewModel.PropertyChanged += OnViewModelChanged;
        Watch(ViewModel.SelectedPage);
    }

    public MainViewModel ViewModel { get; }

    /// <summary>Shows a page (used by the screenshot runner).</summary>
    public void Navigate(GalleryPage page) => ViewModel.SelectedPage = page;

    private void OnViewModelChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName == nameof(MainViewModel.SelectedPage)) Watch(ViewModel.SelectedPage);
    }

    private void Watch(GalleryPage page)
    {
        if (_watched is not null) _watched.PropertyChanged -= OnPageChanged;
        _watched = page;
        page.PropertyChanged += OnPageChanged;
        CodeEditor.Text = page.Code;
        Avalonia.Threading.Dispatcher.UIThread.Post(ViewModel.RefreshExports, Avalonia.Threading.DispatcherPriority.Background);
    }

    private void OnPageChanged(object? sender, PropertyChangedEventArgs e)
    {
        if (e.PropertyName == nameof(GalleryPage.Code) && sender is GalleryPage p)
        {
            CodeEditor.Text = p.Code;
            ViewModel.RefreshExports();
        }
    }

    private void ConfigureEditor()
    {
        var def = HighlightingManager.Instance.GetDefinition("C#");
        if (def is null) return;
        foreach (var named in def.NamedHighlightingColors)
        {
            string? hex = named.Name switch
            {
                "Comment" => "#5E6F7F",
                "String" or "Char" => "#3FD69B",
                "NumberLiteral" => "#FF8A80",
                "MethodCall" => "#F0B03F",
                "Punctuation" => "#93A3B2",
                "Preprocessor" => "#5E6F7F",
                _ when named.Name.Contains("Keyword", StringComparison.Ordinal) || named.Name is "Modifiers" or "Visibility"
                    or "TrueFalse" or "ThisOrBaseReference" or "NullOrValueKeywords" or "ParameterModifiers" or "GetSetAddRemove" => "#7FB8FF",
                _ => "#E2E8EE",
            };
            named.Foreground = new SimpleHighlightingBrush(Color.Parse(hex));
            named.FontWeight = null;
            if (named.Name == "Comment") named.FontStyle = FontStyle.Italic;
        }
        CodeEditor.SyntaxHighlighting = def;
        CodeEditor.TextArea.TextView.LinkTextForegroundBrush = new SolidColorBrush(Color.Parse("#4EA3FF"));
        CodeEditor.LineNumbersForeground = new SolidColorBrush(Color.Parse("#3A4A59"));
    }

    // ---------------- IShellServices ----------------

    public async Task CopyTextAsync(string text)
    {
        if (Clipboard is { } clipboard) await clipboard.SetTextAsync(text);
    }

    public async Task<bool> SaveAsync(string suggestedName, string extension, string description, Func<Stream, Task> write)
    {
        var file = await StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            SuggestedFileName = suggestedName,
            DefaultExtension = extension,
            FileTypeChoices = [new FilePickerFileType(description) { Patterns = [$"*.{extension}"] }],
        });
        if (file is null) return false;
        await using var stream = await file.OpenWriteAsync();
        await write(stream);
        return true;
    }

    public byte[]? CaptureVisualPng(Control control)
    {
        var size = control.Bounds.Size;
        if (size.Width < 1 || size.Height < 1) return null;
        const double scale = 2;
        using var rtb = new RenderTargetBitmap(new PixelSize((int)(size.Width * scale), (int)(size.Height * scale)), new Vector(96 * scale, 96 * scale));
        rtb.Render(control);
        using var ms = new MemoryStream();
        rtb.Save(ms, PngBitmapEncoderOptions.Default);
        return ms.ToArray();
    }
}
