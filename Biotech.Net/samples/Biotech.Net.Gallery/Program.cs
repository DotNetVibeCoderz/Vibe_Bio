using Avalonia;
using Biotech.Net.Gallery.Shell;

namespace Biotech.Net.Gallery;

/// <summary>
/// Biotech.Net Gallery. Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.
/// <c>--screenshots &lt;dir&gt;</c> renders every page headlessly (Avalonia.Headless + Skia) to PNG and exits.
/// </summary>
public static class Program
{
    [STAThread]
    public static int Main(string[] args)
    {
        int i = Array.FindIndex(args, a => a is "--screenshots" or "--screenshot");
        if (i >= 0) return ScreenshotRunner.Run(i + 1 < args.Length ? args[i + 1] : "docs/screenshots");
        BuildAvaloniaApp().StartWithClassicDesktopLifetime(args);
        return 0;
    }

    public static AppBuilder BuildAvaloniaApp() => AppBuilder.Configure<App>()
        .UsePlatformDetect()
        .WithInterFont()
        .LogToTrace();
}
