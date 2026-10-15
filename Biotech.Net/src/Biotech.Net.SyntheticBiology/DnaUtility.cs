using System.Runtime.CompilerServices;
using System.Text;

namespace Biotech.Net.SyntheticBiology;

/// <summary>
/// Low-level DNA helpers shared by the CRISPR, codon, primer and assembly designers: 2-bit encoding,
/// IUPAC matching, reverse complement and GC counting over ASCII spans.
/// </summary>
/// <remarks>
/// Everything here works on <see cref="ReadOnlySpan{T}"/> of ASCII bytes and is case-insensitive, so callers
/// can slice a <see cref="Core.Sequence"/> buffer without copying. The 2-bit helpers are the basis of the
/// genome-wide guide index: A=0, C=1, G=2, T=3, most significant base first, matching the kernel order used
/// by <c>Biotech.Net.Native</c>.
/// </remarks>
public static class DnaUtility
{
    /// <summary>Sentinel returned by <see cref="Code(byte)"/> for a symbol that is not A, C, G, T or U.</summary>
    public const int Ambiguous = -1;

    private static readonly byte[] ComplementTable = BuildComplementTable();

    // Bit sets of the four bases an IUPAC symbol can stand for: bit0=A, bit1=C, bit2=G, bit3=T.
    private static readonly byte[] IupacMask = BuildIupacMasks();

    /// <summary>The 2-bit code of a base (A=0, C=1, G=2, T/U=3), or <see cref="Ambiguous"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Code(byte symbol) => (symbol & 0xDF) switch
    {
        (byte)'A' => 0,
        (byte)'C' => 1,
        (byte)'G' => 2,
        (byte)'T' or (byte)'U' => 3,
        _ => Ambiguous,
    };

    /// <summary>The base for a 2-bit code (0→A, 1→C, 2→G, 3→T).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static char Base(int code) => "ACGT"[code & 3];

    /// <summary>
    /// Packs up to 31 unambiguous bases into a 64-bit word, most significant base first. Returns false when
    /// the window contains any non-ACGT symbol (the caller should then skip the window).
    /// </summary>
    public static bool TryPack(ReadOnlySpan<byte> window, out ulong packed)
    {
        ulong value = 0;
        for (int i = 0; i < window.Length; i++)
        {
            int c = Code(window[i]);
            if (c < 0) { packed = 0; return false; }
            value = (value << 2) | (uint)c;
        }
        packed = value;
        return true;
    }

    /// <summary>Unpacks <paramref name="length"/> bases from a 2-bit word produced by <see cref="TryPack"/>.</summary>
    public static string Unpack(ulong packed, int length)
    {
        return string.Create(length, (packed, length), static (span, state) =>
        {
            var (value, n) = state;
            for (int i = n - 1; i >= 0; i--)
            {
                span[i] = Base((int)(value & 3));
                value >>= 2;
            }
        });
    }

    /// <summary>Number of differing bases between two equal-length 2-bit words of <paramref name="length"/> bases.</summary>
    /// <remarks>
    /// Branch-free: XOR, fold the two bits of every base together, then pop-count. This is the inner loop of
    /// the off-target scan, so it must stay allocation and branch free.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int MismatchCount(ulong a, ulong b, int length)
    {
        ulong diff = a ^ b;
        ulong folded = (diff | (diff >> 1)) & 0x5555_5555_5555_5555UL;
        if (length < 32) folded &= (1UL << (2 * length)) - 1;
        return System.Numerics.BitOperations.PopCount(folded);
    }

    /// <summary>The complement of an ASCII IUPAC symbol, preserving case.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static byte Complement(byte symbol) => ComplementTable[symbol];

    /// <summary>Reverse complement of an ASCII span into a newly allocated string (upper case preserved as given).</summary>
    public static string ReverseComplement(ReadOnlySpan<byte> sequence)
    {
        int n = sequence.Length;
        var buffer = n <= 512 ? stackalloc byte[n] : new byte[n];
        for (int i = 0; i < n; i++) buffer[i] = ComplementTable[sequence[n - 1 - i]];
        return Encoding.ASCII.GetString(buffer);
    }

    /// <summary>Reverse complement of a string.</summary>
    public static string ReverseComplement(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return string.Create(sequence.Length, sequence, static (span, src) =>
        {
            for (int i = 0; i < span.Length; i++) span[i] = (char)ComplementTable[(byte)src[src.Length - 1 - i]];
        });
    }

    /// <summary>Whether an observed base matches an IUPAC pattern symbol (N matches everything).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool MatchesIupac(byte observed, byte pattern)
    {
        int c = Code(observed);
        if (c < 0) return false;
        return (IupacMask[pattern] & (1 << c)) != 0;
    }

    /// <summary>Whether a window matches an IUPAC pattern of the same length.</summary>
    public static bool MatchesIupac(ReadOnlySpan<byte> window, ReadOnlySpan<byte> pattern)
    {
        if (window.Length != pattern.Length) return false;
        for (int i = 0; i < pattern.Length; i++)
            if (!MatchesIupac(window[i], pattern[i])) return false;
        return true;
    }

    /// <summary>Number of G and C bases (case-insensitive, ambiguity codes ignored).</summary>
    public static int CountGc(ReadOnlySpan<byte> sequence)
    {
        int gc = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            int b = sequence[i] & 0xDF;
            if (b == 'G' || b == 'C') gc++;
        }
        return gc;
    }

    /// <summary>GC fraction (0–1); 0 for an empty span.</summary>
    public static double GcFraction(ReadOnlySpan<byte> sequence)
        => sequence.Length == 0 ? 0 : (double)CountGc(sequence) / sequence.Length;

    /// <summary>GC fraction (0–1) of a string.</summary>
    public static double GcFraction(string sequence) => GcFraction(Ascii(sequence));

    /// <summary>Length of the longest run of one base (homopolymer).</summary>
    public static int LongestHomopolymer(ReadOnlySpan<byte> sequence)
    {
        int best = 0, run = 0;
        int previous = -1;
        for (int i = 0; i < sequence.Length; i++)
        {
            int b = sequence[i] & 0xDF;
            if (b == previous) run++;
            else { previous = b; run = 1; }
            if (run > best) best = run;
        }
        return best;
    }

    /// <summary>ASCII bytes of a string (one allocation; used where a span over a literal is needed).</summary>
    public static byte[] Ascii(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return Encoding.ASCII.GetBytes(sequence);
    }

    /// <summary>Upper-cased copy with U rewritten as T, for inputs that may be RNA or soft-masked.</summary>
    public static string NormalizeDna(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return string.Create(sequence.Length, sequence, static (span, src) =>
        {
            for (int i = 0; i < span.Length; i++)
            {
                char c = char.ToUpperInvariant(src[i]);
                span[i] = c == 'U' ? 'T' : c;
            }
        });
    }

    /// <summary>Whether every symbol is A, C, G or T (case-insensitive).</summary>
    public static bool IsUnambiguous(ReadOnlySpan<byte> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
            if (Code(sequence[i]) < 0) return false;
        return true;
    }

    private static byte[] BuildComplementTable()
    {
        var table = new byte[256];
        for (int i = 0; i < 256; i++) table[i] = (byte)i;
        const string from = "ACGTURYSWKMBDHVN";
        const string to = "TGCAAYRSWMKVHDBN";
        for (int i = 0; i < from.Length; i++)
        {
            table[from[i]] = (byte)to[i];
            table[char.ToLowerInvariant(from[i])] = (byte)char.ToLowerInvariant(to[i]);
        }
        return table;
    }

    private static byte[] BuildIupacMasks()
    {
        var masks = new byte[256];
        void Set(char symbol, int mask)
        {
            masks[char.ToUpperInvariant(symbol)] = (byte)mask;
            masks[char.ToLowerInvariant(symbol)] = (byte)mask;
        }
        const int a = 1, c = 2, g = 4, t = 8;
        Set('A', a); Set('C', c); Set('G', g); Set('T', t); Set('U', t);
        Set('R', a | g); Set('Y', c | t); Set('S', c | g); Set('W', a | t);
        Set('K', g | t); Set('M', a | c);
        Set('B', c | g | t); Set('D', a | g | t); Set('H', a | c | t); Set('V', a | c | g);
        Set('N', a | c | g | t);
        return masks;
    }
}
