namespace Biotech.Net.Genomics.IO;

/// <summary>
/// A BGZF virtual file offset: the file offset of a compressed block in the upper 48 bits and the offset inside the
/// decompressed block in the lower 16 bits. Ordering follows file order, as BAI/CSI/tabix indexes require.
/// </summary>
public readonly record struct VirtualOffset(ulong Value) : IComparable<VirtualOffset>
{
    /// <summary>Creates an offset from a block address and an in-block offset.</summary>
    public VirtualOffset(long blockAddress, int blockOffset) : this(((ulong)blockAddress << 16) | (ushort)blockOffset)
    {
        if (blockAddress < 0 || blockAddress >= 1L << 48) throw new ArgumentOutOfRangeException(nameof(blockAddress));
        if (blockOffset is < 0 or > 0xFFFF) throw new ArgumentOutOfRangeException(nameof(blockOffset));
    }

    /// <summary>File offset of the BGZF block.</summary>
    public long BlockAddress => (long)(Value >> 16);

    /// <summary>Offset inside the decompressed block.</summary>
    public int BlockOffset => (int)(Value & 0xFFFF);

    /// <inheritdoc/>
    public int CompareTo(VirtualOffset other) => Value.CompareTo(other.Value);

    public static bool operator <(VirtualOffset a, VirtualOffset b) => a.Value < b.Value;
    public static bool operator >(VirtualOffset a, VirtualOffset b) => a.Value > b.Value;
    public static bool operator <=(VirtualOffset a, VirtualOffset b) => a.Value <= b.Value;
    public static bool operator >=(VirtualOffset a, VirtualOffset b) => a.Value >= b.Value;

    public override string ToString() => $"{BlockAddress}:{BlockOffset}";
}
