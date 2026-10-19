namespace Biotech.Net.Gpu.Tests;

/// <summary>Shared accelerators: ILGPU's CPU accelerator (always) and a real GPU when one is present.</summary>
internal static class GpuTestDevices
{
    private static readonly Lazy<GpuDevice> CpuDevice = new(GpuDevice.CreateCpu);
    private static readonly Lazy<GpuDevice?> RealGpu = new(() =>
    {
        var d = GpuDevice.TryCreate(GpuBackend.Cuda) ?? GpuDevice.TryCreate(GpuBackend.OpenCL);
        return d;
    });

    /// <summary>The CPU accelerator: deterministic, runs the same kernels in CI.</summary>
    public static GpuDevice Cpu => CpuDevice.Value;

    /// <summary>A CUDA or OpenCL device, or null (GPU-gated tests then return early).</summary>
    public static GpuDevice? Gpu => RealGpu.Value;

    /// <summary>The CPU accelerator plus the real GPU when present.</summary>
    public static IEnumerable<GpuDevice> All => Gpu is { } g ? [Cpu, g] : [Cpu];
}
