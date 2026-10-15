using System.Runtime.InteropServices;

namespace Biotech.Net.Native;

/// <summary>
/// Base composition of a nucleotide sequence. Case-insensitive; <c>U</c> is counted as <see cref="T"/>.
/// Laid out sequentially because the Rust kernel fills the same struct across the C ABI.
/// </summary>
[StructLayout(LayoutKind.Sequential)]
public readonly record struct NucleotideCounts(long A, long C, long G, long T, long N, long Other)
{
    /// <summary>Number of unambiguous bases (A, C, G, T/U).</summary>
    public long Acgt => A + C + G + T;

    /// <summary>Every symbol counted, including <c>N</c> and other characters.</summary>
    public long Total => A + C + G + T + N + Other;

    /// <summary>G+C over A+C+G+T; ambiguous bases are excluded from the denominator. 0 when there are no ACGT bases.</summary>
    public double GcFraction => Acgt == 0 ? 0 : (double)(G + C) / Acgt;

    /// <summary>(G − C) / (G + C), the GC skew used to locate replication origins.</summary>
    public double GcSkew => G + C == 0 ? 0 : (double)(G - C) / (G + C);
}
