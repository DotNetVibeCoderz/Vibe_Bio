```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  Categories=All-vs-all SW 64 x 128 (100-300 aa)  

```
| Method   | Backend      | Mean        | Error       | StdDev    | Gen0      | Gen1      | Gen2      | Allocated  |
|--------- |------------- |------------:|------------:|----------:|----------:|----------:|----------:|-----------:|
| AllVsAll | ILGPU-CPU    | 3,362.28 ms | 1,307.33 ms | 71.659 ms |         - |         - |         - |   75.98 KB |
| AllVsAll | ILGPU-OpenCL |   112.99 ms |    52.53 ms |  2.880 ms |         - |         - |         - |    76.2 KB |
| AllVsAll | Managed      | 1,008.20 ms | 1,652.77 ms | 90.594 ms | 1000.0000 | 1000.0000 | 1000.0000 | 6346.73 KB |
| AllVsAll | Rust         |    33.69 ms |    68.18 ms |  3.737 ms |         - |         - |         - |    3.27 KB |
