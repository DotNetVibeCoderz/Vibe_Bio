```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  

```
| Method                            | Mean     | Error     | StdDev    | Gen0       | Allocated |
|---------------------------------- |---------:|----------:|----------:|-----------:|----------:|
| &#39;FASTA 2,000 x 5 kb (10 Mb)&#39;      | 61.06 ms |  75.46 ms |  4.136 ms |  8222.2222 |  32.98 MB |
| &#39;FASTQ 50,000 x 150 bp&#39;           | 64.07 ms | 221.74 ms | 12.155 ms | 18428.5714 |  73.55 MB |
| &#39;VCF 20,000 records x 10 samples&#39; | 25.53 ms |  58.14 ms |  3.187 ms |  6937.5000 |  27.76 MB |
