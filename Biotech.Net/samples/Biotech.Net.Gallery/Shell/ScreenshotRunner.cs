using Avalonia;
using Avalonia.Headless;
using Avalonia.Threading;

namespace Biotech.Net.Gallery.Shell;

/// <summary>
/// Renders every page to PNG through Avalonia.Headless with real Skia drawing (1440×900, independent of the screen):
/// <c>dotnet run --project samples/Biotech.Net.Gallery -- --screenshots docs/screenshots</c>.
/// Writes <c>gallery-&lt;page&gt;.png</c> for each page plus <c>gallery-live-code.png</c>.
/// </summary>
public static class ScreenshotRunner
{
    public static int Run(string directory)
    {
        Directory.CreateDirectory(directory);
        AppBuilder.Configure<App>()
            .UseSkia()
            .UseHeadless(new AvaloniaHeadlessPlatformOptions { UseHeadlessDrawing = false })
            .WithInterFont()
            .SetupWithoutStarting();

        var window = new MainWindow { Width = 1440, Height = 900 };
        window.Show();
        Pump(TimeSpan.FromMilliseconds(300));

        int failures = 0;
        foreach (var page in window.ViewModel.Pages)
        {
            window.ViewModel.ShowCode = false;
            window.Navigate(page);
            Pump(TimeSpan.FromMilliseconds(400));
            var ready = page.WhenReadyAsync();
            var deadline = DateTime.UtcNow + TimeSpan.FromSeconds(120);
            while (!ready.IsCompleted && DateTime.UtcNow < deadline) Pump(TimeSpan.FromMilliseconds(100));
            if (!ready.IsCompleted) { Console.Error.WriteLine($"{page.Title}: timed out waiting for content"); failures++; }
            else if (ready.IsFaulted) { Console.Error.WriteLine($"{page.Title}: {ready.Exception?.GetBaseException().Message}"); failures++; }
            Pump(TimeSpan.FromMilliseconds(500));
            Save(window, Path.Combine(directory, $"gallery-{page.Slug}.png"));
        }

        // The Live code tab, shown on the DNA Explorer.
        window.Navigate(window.ViewModel.Pages[0]);
        window.ViewModel.ShowCode = true;
        Pump(TimeSpan.FromMilliseconds(500));
        Save(window, Path.Combine(directory, "gallery-live-code.png"));

        window.Close();
        return failures == 0 ? 0 : 1;
    }

    private static void Save(MainWindow window, string file)
    {
        var frame = window.CaptureRenderedFrame();
        if (frame is null) { Console.Error.WriteLine($"no frame for {file}"); return; }
        using (var fs = File.Create(file)) frame.Save(fs);
        Console.WriteLine($"saved {file}");
    }

    /// <summary>Runs dispatcher jobs and render ticks for roughly <paramref name="duration"/>.</summary>
    private static void Pump(TimeSpan duration)
    {
        var end = DateTime.UtcNow + duration;
        do
        {
            Dispatcher.UIThread.RunJobs();
            AvaloniaHeadlessPlatform.ForceRenderTimerTick();
            Thread.Sleep(15);
        } while (DateTime.UtcNow < end);
        Dispatcher.UIThread.RunJobs();
    }
}
