using System.Collections.ObjectModel;
using Biotech.Net.Core;
using Biotech.Net.Gallery.Pages;
using Biotech.Net.Native;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;

namespace Biotech.Net.Gallery.Shell;

/// <summary>Services the window provides to the view model (clipboard, save dialog, visual capture).</summary>
public interface IShellServices
{
    Task CopyTextAsync(string text);
    Task<bool> SaveAsync(string suggestedName, string extension, string description, Func<Stream, Task> write);
    byte[]? CaptureVisualPng(Avalonia.Controls.Control control);
}

public sealed partial class MainViewModel : ObservableObject
{
    public MainViewModel()
    {
        Pages =
        [
            new DnaExplorerPage(),
            new GenomeBrowserPage(),
            new ProteinExplorerPage(),
            new DrugDiscoveryPage(),
            new AiPlaygroundPage(),
            new CrisprDesignerPage(),
            new AboutPage(),
        ];
        _selectedPage = Pages[0];
        Pages[0].IsActive = true;
    }

    public IShellServices? Services { get; set; }

    public ObservableCollection<GalleryPage> Pages { get; }

    [ObservableProperty] private GalleryPage _selectedPage;

    /// <summary>True while the Live code tab is shown.</summary>
    [ObservableProperty] private bool _showCode;

    /// <summary>Transient feedback after an export ("Saved crispr.svg").</summary>
    [ObservableProperty] private string _toast = "";

    public string Version => $"v{BiotechInfo.Version}";

    public string Kernels => BioKernels.Backend == KernelBackend.Native ? "Rust kernels" : "Managed kernels";

    public string Credits => BiotechInfo.CreditsId;

    partial void OnSelectedPageChanged(GalleryPage? oldValue, GalleryPage newValue)
    {
        if (oldValue is not null) oldValue.IsActive = false;
        newValue.IsActive = true;
        Toast = "";
        OnPropertyChanged(nameof(CanSaveSvg));
    }

    public bool CanSaveSvg => SelectedPage.ExportPlot is not null;

    [RelayCommand]
    private void Navigate(GalleryPage page) => SelectedPage = page;

    [RelayCommand]
    private void ShowResult() => ShowCode = false;

    [RelayCommand]
    private void ShowLiveCode() => ShowCode = true;

    [RelayCommand]
    private async Task CopyCode()
    {
        if (Services is null) return;
        await Services.CopyTextAsync(SelectedPage.Code);
        Toast = "Snippet copied";
    }

    [RelayCommand]
    private async Task CopyResult()
    {
        if (Services is null) return;
        await Services.CopyTextAsync(SelectedPage.ResultText);
        Toast = "Result copied";
    }

    [RelayCommand]
    private async Task SavePng()
    {
        if (Services is null) return;
        var page = SelectedPage;
        byte[]? png = page.ExportPlot is { } plot ? plot.ToPng(2) : page.ExportVisual is { } v ? Services.CaptureVisualPng(v) : null;
        if (png is null) { Toast = "Nothing to export on this page"; return; }
        if (await Services.SaveAsync($"{page.Slug}.png", "png", "PNG image", s => s.WriteAsync(png).AsTask()))
            Toast = $"Saved {page.Slug}.png";
    }

    [RelayCommand]
    private async Task SaveSvg()
    {
        if (Services is null || SelectedPage.ExportPlot is not { } plot) { Toast = "This page has no vector figure"; return; }
        string svg = plot.ToSvg();
        if (await Services.SaveAsync($"{SelectedPage.Slug}.svg", "svg", "SVG drawing", async s =>
            {
                await using var w = new StreamWriter(s);
                await w.WriteAsync(svg);
            }))
            Toast = $"Saved {SelectedPage.Slug}.svg";
    }

    /// <summary>Lets pages refresh export availability after recomputing.</summary>
    public void RefreshExports() => OnPropertyChanged(nameof(CanSaveSvg));
}
