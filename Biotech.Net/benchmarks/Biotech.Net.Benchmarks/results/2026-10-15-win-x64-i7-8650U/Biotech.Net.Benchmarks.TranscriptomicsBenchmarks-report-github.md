```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  

```
| Method                          | Mean      | Error      | StdDev    | Gen0      | Gen1      | Gen2      | Allocated   |
|-------------------------------- |----------:|-----------:|----------:|----------:|----------:|----------:|------------:|
| &#39;Median-of-ratios size factors&#39; | 25.086 ms | 43.7113 ms | 2.3960 ms |   31.2500 |   31.2500 |   31.2500 |   312.67 KB |
| &#39;Normalized counts (DESeq)&#39;     |  1.318 ms |  1.4424 ms | 0.0791 ms |  289.0625 |  267.5781 |  267.5781 |  2453.98 KB |
| &#39;TMM factors (edgeR)&#39;           | 62.866 ms | 28.8186 ms | 1.5796 ms | 4625.0000 | 4500.0000 | 4000.0000 | 20344.42 KB |
| log-CPM                         |  5.228 ms |  0.8870 ms | 0.0486 ms |  468.7500 |  453.1250 |  453.1250 |  2456.82 KB |
