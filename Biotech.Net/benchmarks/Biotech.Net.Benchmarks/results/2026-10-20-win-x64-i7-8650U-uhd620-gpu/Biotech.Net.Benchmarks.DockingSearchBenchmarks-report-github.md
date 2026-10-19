```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  Categories=Vina docking 3PTB (default budget)  

```
| Method | Backend      | Mean    | Error    | StdDev   | Gen0       | Gen1      | Gen2      | Allocated |
|------- |------------- |--------:|---------:|---------:|-----------:|----------:|----------:|----------:|
| Dock   | CPU-Vina     | 1.550 s | 0.3579 s | 0.0196 s | 49000.0000 | 1000.0000 | 1000.0000 |  197.8 MB |
| Dock   | ILGPU-OpenCL | 1.551 s | 0.1973 s | 0.0108 s |          - |         - |         - |   9.36 MB |
