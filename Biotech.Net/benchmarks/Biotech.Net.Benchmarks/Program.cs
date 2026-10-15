using BenchmarkDotNet.Columns;
using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Diagnosers;
using BenchmarkDotNet.Exporters;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Loggers;
using BenchmarkDotNet.Order;
using BenchmarkDotNet.Running;
using BenchmarkDotNet.Toolchains.InProcess.Emit;
using Biotech.Net.Native;

// Biotech.Net benchmarks — dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.
//
//   dotnet run -c Release -p:SkipRustBuild=true -- --filter *           short in-process job (a few minutes)
//   dotnet run -c Release -p:SkipRustBuild=true -- --filter *Kernel* --full   BenchmarkDotNet default job
//
// Kernel benchmarks call the Rust bridge and the managed fallback directly (both are internal and
// visible to this assembly), so one process measures both back-ends side by side.

bool full = args.Contains("--full");
args = args.Where(a => a != "--full").ToArray();

Console.WriteLine($"Rust kernels available: {BioKernels.IsNativeAvailable} ({BioKernels.NativeVersion})");
if (!BioKernels.IsNativeAvailable)
    Console.WriteLine("WARNING: the native library was not found, so the *_Rust benchmarks will fail. Build native/biotech-native first.");

var job = full
    ? Job.Default
    : Job.ShortRun.WithToolchain(InProcessEmitToolchain.Instance).WithId("ShortRun-InProcess");

var config = ManualConfig.CreateEmpty()
    .AddJob(job)
    .AddColumnProvider(DefaultColumnProviders.Instance)
    .AddColumn(CategoriesColumn.Default)
    .AddDiagnoser(MemoryDiagnoser.Default)
    .AddExporter(MarkdownExporter.GitHub)
    .AddLogger(ConsoleLogger.Default)
    .WithOrderer(new DefaultOrderer(SummaryOrderPolicy.Declared))
    .AddLogicalGroupRules(BenchmarkLogicalGroupRule.ByCategory);

BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args, config);
