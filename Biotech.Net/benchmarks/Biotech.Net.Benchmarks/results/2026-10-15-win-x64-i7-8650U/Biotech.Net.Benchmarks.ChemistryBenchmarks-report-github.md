```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  

```
| Method                                     | Mean     | Error    | StdDev   | Gen0     | Gen1    | Allocated |
|------------------------------------------- |---------:|---------:|---------:|---------:|--------:|----------:|
| &#39;Parse 20 drug SMILES&#39;                     | 429.4 μs | 404.2 μs | 22.16 μs | 100.0977 |       - |  410.3 KB |
| &#39;Morgan r=2 2048-bit x 20 molecules&#39;       | 175.4 μs | 169.8 μs |  9.31 μs |  32.7148 |       - | 133.79 KB |
| &#39;Tanimoto top-10 over 10,000 fingerprints&#39; | 633.3 μs | 623.9 μs | 34.20 μs |  67.3828 | 22.4609 | 277.56 KB |
| &#39;Bulk Tanimoto, 10,000 fingerprints&#39;       | 238.9 μs | 189.7 μs | 10.40 μs |  39.0625 | 10.0098 |  159.2 KB |
