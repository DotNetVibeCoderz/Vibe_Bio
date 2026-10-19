using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;

namespace Biotech.Net.Gpu;

/// <summary>Which ILGPU accelerator to use.</summary>
public enum GpuBackend
{
    /// <summary>The best GPU present: CUDA, then OpenCL. Never picks the CPU accelerator; no GPU means no device.</summary>
    Auto = 0,

    /// <summary>An NVIDIA GPU through CUDA.</summary>
    Cuda,

    /// <summary>An OpenCL device (AMD, Intel or NVIDIA GPUs).</summary>
    OpenCL,

    /// <summary>ILGPU's CPU accelerator: runs the GPU kernels on CPU threads. Always available and deterministic; used by the tests.</summary>
    Cpu,

    /// <summary>No accelerator: every API runs the existing Rust/managed CPU code.</summary>
    Off,
}

/// <summary>An accelerator ILGPU can see on this machine.</summary>
/// <param name="Backend">CUDA, OpenCL or CPU.</param>
/// <param name="Name">Device name reported by the driver.</param>
/// <param name="MemoryBytes">Device memory in bytes (the CPU accelerator reports <see cref="long.MaxValue"/>).</param>
/// <param name="MaxThreadsPerGroup">Largest thread group the device supports.</param>
public sealed record GpuDeviceInfo(GpuBackend Backend, string Name, long MemoryBytes, int MaxThreadsPerGroup)
{
    /// <inheritdoc/>
    public override string ToString() => Backend == GpuBackend.Cpu
        ? $"{Backend}: {Name}"
        : $"{Backend}: {Name} ({MemoryBytes / (1024.0 * 1024 * 1024):F1} GiB)";
}

/// <summary>
/// An ILGPU accelerator plus the kernels compiled for it. Create one with <see cref="TryCreate"/> (null when the requested
/// backend is absent or switched off) and dispose it when done; every GPU API in this package also accepts a null device
/// and then runs the existing Rust/managed CPU implementation.
/// </summary>
/// <remarks>
/// The environment variable <c>BIOTECH_NET_GPU</c> (<c>cuda</c>, <c>opencl</c>, <c>cpu</c>, <c>off</c> or <c>auto</c>)
/// overrides <see cref="GpuBackend.Auto"/>; an explicit backend argument wins over it. Kernel launches on one device are
/// serialised, so a device can be shared between threads.
/// </remarks>
public sealed class GpuDevice : IDisposable
{
    /// <summary>Name of the environment variable that overrides <see cref="GpuBackend.Auto"/>.</summary>
    public const string EnvironmentVariable = "BIOTECH_NET_GPU";

    private readonly Context _context;
    private readonly Dictionary<string, object> _kernels = [];
    private bool _disposed;

    private GpuDevice(Context context, Accelerator accelerator, GpuBackend backend)
    {
        _context = context;
        Accelerator = accelerator;
        Backend = backend;
        Info = new GpuDeviceInfo(backend, accelerator.Name, accelerator.MemorySize, accelerator.MaxNumThreadsPerGroup);
    }

    /// <summary>The backend of this device.</summary>
    public GpuBackend Backend { get; }

    /// <summary>Name, memory and limits.</summary>
    public GpuDeviceInfo Info { get; }

    /// <summary>Device name.</summary>
    public string Name => Info.Name;

    /// <summary>True for CUDA and OpenCL devices, false for the CPU accelerator.</summary>
    public bool IsGpu => Backend is GpuBackend.Cuda or GpuBackend.OpenCL;

    internal Accelerator Accelerator { get; }

    /// <summary>Serialises launches and buffer traffic on this device.</summary>
    internal object Sync { get; } = new();

    /// <summary>The backend requested by <c>BIOTECH_NET_GPU</c>, or <see cref="GpuBackend.Auto"/> when unset or unrecognised.</summary>
    public static GpuBackend EnvironmentBackend => ParseBackend(Environment.GetEnvironmentVariable(EnvironmentVariable));

    /// <summary>Parses <c>cuda</c>, <c>opencl</c>, <c>cpu</c>, <c>off</c> (also <c>none</c>, <c>0</c>, <c>false</c>) or <c>auto</c>; anything else is Auto.</summary>
    public static GpuBackend ParseBackend(string? value) => value?.Trim().ToLowerInvariant() switch
    {
        "cuda" => GpuBackend.Cuda,
        "opencl" or "cl" => GpuBackend.OpenCL,
        "cpu" => GpuBackend.Cpu,
        "off" or "none" or "0" or "false" => GpuBackend.Off,
        _ => GpuBackend.Auto,
    };

    /// <summary>Every accelerator ILGPU can see (CUDA and OpenCL devices plus the CPU accelerator).</summary>
    public static IReadOnlyList<GpuDeviceInfo> Enumerate()
    {
        try
        {
            using var context = CreateContext();
            return context.Devices.Select(d => new GpuDeviceInfo(ToBackend(d.AcceleratorType), d.Name, d.MemorySize, d.MaxNumThreadsPerGroup)).ToList();
        }
        catch (Exception)
        {
            return [];
        }
    }

    /// <summary>True when a CUDA or OpenCL device is present.</summary>
    public static bool IsGpuAvailable => Enumerate().Any(d => d.Backend is GpuBackend.Cuda or GpuBackend.OpenCL);

    /// <summary>
    /// Opens an accelerator. <see cref="GpuBackend.Auto"/> honours <c>BIOTECH_NET_GPU</c> and otherwise picks CUDA, then the
    /// OpenCL device with most memory. Returns null for <see cref="GpuBackend.Off"/>, when no matching device exists, or
    /// when the driver fails to initialise — callers then use the CPU path.
    /// </summary>
    public static GpuDevice? TryCreate(GpuBackend backend = GpuBackend.Auto)
    {
        if (backend == GpuBackend.Auto) backend = EnvironmentBackend;
        if (backend == GpuBackend.Off) return null;
        Context? context = null;
        try
        {
            context = CreateContext();
            Device? device = backend switch
            {
                GpuBackend.Cuda => context.Devices.Where(d => d.AcceleratorType == AcceleratorType.Cuda).MaxBy(d => d.MemorySize),
                GpuBackend.OpenCL => context.Devices.Where(d => d.AcceleratorType == AcceleratorType.OpenCL).MaxBy(d => d.MemorySize),
                GpuBackend.Cpu => context.Devices.FirstOrDefault(d => d.AcceleratorType == AcceleratorType.CPU),
                _ => context.Devices.Where(d => d.AcceleratorType == AcceleratorType.Cuda).MaxBy(d => d.MemorySize)
                     ?? context.Devices.Where(d => d.AcceleratorType == AcceleratorType.OpenCL).MaxBy(d => d.MemorySize),
            };
            if (device is null) { context.Dispose(); return null; }
            var accelerator = device.CreateAccelerator(context);
            return new GpuDevice(context, accelerator, ToBackend(device.AcceleratorType));
        }
        catch (Exception)
        {
            context?.Dispose();
            return null;
        }
    }

    /// <summary>Opens the CPU accelerator (always succeeds). Useful for tests and for checking the GPU kernels without a GPU.</summary>
    public static GpuDevice CreateCpu() => TryCreate(GpuBackend.Cpu) ?? throw new InvalidOperationException("ILGPU's CPU accelerator is unavailable.");

    private static Context CreateContext() => Context.Create(b => b.Default().EnableAlgorithms().Math(MathMode.Default));

    private static GpuBackend ToBackend(AcceleratorType t) => t switch
    {
        AcceleratorType.Cuda => GpuBackend.Cuda,
        AcceleratorType.OpenCL => GpuBackend.OpenCL,
        _ => GpuBackend.Cpu,
    };

    /// <summary>A kernel compiled once per device and cached.</summary>
    internal T Kernel<T>(string key, Func<Accelerator, T> load) where T : class
    {
        lock (_kernels)
        {
            ObjectDisposedException.ThrowIf(_disposed, this);
            if (_kernels.TryGetValue(key, out var k)) return (T)k;
            var kernel = load(Accelerator);
            _kernels[key] = kernel;
            return kernel;
        }
    }

    /// <inheritdoc/>
    public override string ToString() => Info.ToString();

    /// <summary>Releases the accelerator and its ILGPU context.</summary>
    public void Dispose()
    {
        lock (_kernels)
        {
            if (_disposed) return;
            _disposed = true;
            _kernels.Clear();
        }
        Accelerator.Dispose();
        _context.Dispose();
    }

    internal void ThrowIfDisposed() => ObjectDisposedException.ThrowIf(_disposed, this);
}
