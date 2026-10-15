using System.Reflection;
using Biotech.Net.Native;

namespace Biotech.Net.Core;

/// <summary>Product information and runtime details for Biotech.Net.</summary>
public static class BiotechInfo
{
    /// <summary>Product name.</summary>
    public const string Product = "Biotech.Net";

    /// <summary>Tagline.</summary>
    public const string Tagline = "Ultimate Biotechnology & Bioinformatics Platform for .NET";

    /// <summary>Credits (English).</summary>
    public const string Credits = "Created by Gravicode Studios, led by Kang Fadhil";

    /// <summary>Credits (Bahasa Indonesia).</summary>
    public const string CreditsId = "Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil";

    /// <summary>Informational version of the library.</summary>
    public static string Version { get; } =
        typeof(BiotechInfo).Assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion.Split('+')[0]
        ?? typeof(BiotechInfo).Assembly.GetName().Version?.ToString() ?? "0.0.0";

    /// <summary>Kernel backend in use (Rust or managed).</summary>
    public static KernelBackend Backend => BioKernels.Backend;

    /// <summary>One-line banner for consoles and notebooks.</summary>
    public static string Banner => $"{Product} {Version} — {Tagline} · kernels: {BioKernels.NativeVersion} · {CreditsId}";
}
