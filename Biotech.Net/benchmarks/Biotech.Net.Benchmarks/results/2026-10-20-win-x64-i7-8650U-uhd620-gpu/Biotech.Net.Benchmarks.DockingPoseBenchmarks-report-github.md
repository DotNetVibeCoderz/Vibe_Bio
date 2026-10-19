```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  Categories=Vina pose energy+gradient x4096 (1STP)  

```
| Method   | Backend          | Mean      | Error      | StdDev     | Gen0     | Gen1    | Gen2    | Allocated |
|--------- |----------------- |----------:|-----------:|-----------:|---------:|--------:|--------:|----------:|
| Evaluate | ILGPU-CPU        | 54.991 ms | 449.675 ms | 24.6482 ms |        - |       - |       - | 420.53 KB |
| Evaluate | ILGPU-OpenCL     |  5.502 ms |   7.671 ms |  0.4205 ms |  46.8750 | 39.0625 | 39.0625 | 418.18 KB |
| Evaluate | Managed-1T       | 22.601 ms |  31.051 ms |  1.7020 ms | 125.0000 |       - |       - | 609.14 KB |
| Evaluate | Managed-Parallel |  8.347 ms |  15.379 ms |  0.8430 ms | 156.2500 | 31.2500 |       - | 650.54 KB |
