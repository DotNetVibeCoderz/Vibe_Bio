```

BenchmarkDotNet v0.15.8, Windows 11 (10.0.26200.9711/25H2/2025Update/HudsonValley2)
Intel Core i7-8650U CPU 1.90GHz (Max: 2.11GHz) (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET SDK 10.0.401
  [Host] : .NET 10.0.12 (10.0.12, 10.0.1226.42308), X64 RyuJIT x86-64-v3

Job=ShortRun-InProcess  Toolchain=InProcessEmitToolchain  IterationCount=3  
LaunchCount=1  WarmupCount=3  

```
| Method                    | Categories                     | Mean           | Error          | StdDev        | Ratio | RatioSD | Gen0     | Gen1     | Gen2     | Allocated | Alloc Ratio |
|-------------------------- |------------------------------- |---------------:|---------------:|--------------:|------:|--------:|---------:|---------:|---------:|----------:|------------:|
| BenjaminiHochberg_Managed | Benjamini-Hochberg 100k        |  10,417.261 μs |   6,840.782 μs |   374.9663 μs |  1.00 |    0.04 | 109.3750 | 109.3750 | 109.3750 | 1200710 B |       1.000 |
| BenjaminiHochberg_Rust    | Benjamini-Hochberg 100k        |  11,864.286 μs |  10,339.949 μs |   566.7675 μs |  1.14 |    0.06 |        - |        - |        - |      40 B |       0.000 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| CountKmers_Managed        | CountKmers k=8, 1 Mb           |   4,698.151 μs |  19,074.002 μs | 1,045.5104 μs |  1.03 |    0.27 |        - |        - |        - |      42 B |        1.00 |
| CountKmers_Rust           | CountKmers k=8, 1 Mb           |   2,013.714 μs |   4,075.002 μs |   223.3646 μs |  0.44 |    0.09 |        - |        - |        - |      10 B |        0.24 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| GotohGlobal_Managed       | Gotoh global 2 kb x 2 kb       |  42,742.585 μs |  67,238.147 μs | 3,685.5497 μs |  1.01 |    0.11 | 153.8462 | 153.8462 | 153.8462 | 4049495 B |       1.000 |
| GotohGlobal_Rust          | Gotoh global 2 kb x 2 kb       |  51,905.561 μs |  27,242.051 μs | 1,493.2288 μs |  1.22 |    0.10 |        - |        - |        - |     458 B |       0.000 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| IsoelectricPoint_Managed  | Isoelectric point 10k aa       |      13.075 μs |       8.219 μs |     0.4505 μs |  1.00 |    0.04 |        - |        - |        - |         - |          NA |
| IsoelectricPoint_Rust     | Isoelectric point 10k aa       |      11.361 μs |      12.936 μs |     0.7091 μs |  0.87 |    0.05 |        - |        - |        - |         - |          NA |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| Kabsch_Managed            | Kabsch RMSD 10k atoms          |     246.479 μs |      13.949 μs |     0.7646 μs |  1.00 |    0.00 |        - |        - |        - |       3 B |        1.00 |
| Kabsch_Rust               | Kabsch RMSD 10k atoms          |      96.131 μs |      42.806 μs |     2.3463 μs |  0.39 |    0.01 |        - |        - |        - |       1 B |        0.33 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| LocalScoreBatch_Managed   | Local score batch 300 aa x 500 | 295,240.517 μs | 162,589.896 μs | 8,912.1009 μs |  1.00 |    0.04 |        - |        - |        - |  158736 B |       1.000 |
| LocalScoreBatch_Rust      | Local score batch 300 aa x 500 |   6,270.059 μs |   1,628.539 μs |    89.2657 μs |  0.02 |    0.00 |        - |        - |        - |      42 B |       0.000 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| ProteinMass_Managed       | Protein mass 10k aa            |      10.720 μs |       2.911 μs |     0.1596 μs |  1.00 |    0.02 |        - |        - |        - |         - |          NA |
| ProteinMass_Rust          | Protein mass 10k aa            |       9.814 μs |       1.072 μs |     0.0588 μs |  0.92 |    0.01 |        - |        - |        - |         - |          NA |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| ReverseComplement_Managed | ReverseComplement 1 Mb         |     625.776 μs |     149.754 μs |     8.2085 μs |  1.00 |    0.02 |        - |        - |        - |       5 B |        1.00 |
| ReverseComplement_Rust    | ReverseComplement 1 Mb         |      80.694 μs |      38.811 μs |     2.1274 μs |  0.13 |    0.00 |        - |        - |        - |       1 B |        0.20 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| Translate_Managed         | Translate 1 Mb                 |     985.394 μs |     199.465 μs |    10.9333 μs |  1.00 |    0.01 |        - |        - |        - |      10 B |        1.00 |
| Translate_Rust            | Translate 1 Mb                 |     718.266 μs |     182.519 μs |    10.0045 μs |  0.73 |    0.01 |        - |        - |        - |       5 B |        0.50 |
|                           |                                |                |                |               |       |         |          |          |          |           |             |
| WelchTTest_Managed        | Welch t-test 20k x 12          |   4,610.205 μs |     774.315 μs |    42.4428 μs |  1.00 |    0.01 |        - |        - |        - |      42 B |        1.00 |
| WelchTTest_Rust           | Welch t-test 20k x 12          |   5,738.520 μs |   1,518.245 μs |    83.2202 μs |  1.24 |    0.02 |        - |        - |        - |      42 B |        1.00 |
