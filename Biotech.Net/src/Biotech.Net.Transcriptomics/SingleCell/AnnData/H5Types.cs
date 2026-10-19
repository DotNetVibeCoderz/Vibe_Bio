using PureHDF;
using PureHDF.Selections;

namespace Biotech.Net.Transcriptomics.SingleCell;

/// <summary>The HDF5 boolean type h5py and numpy use: a 1-byte enum {FALSE = 0, TRUE = 1}.</summary>
internal enum H5Bool : byte
{
    FALSE = 0,
    TRUE = 1,
}

/// <summary>
/// Typed reading over PureHDF (which requires the .NET element type to match the stored type exactly) and the
/// mapping between HDF5 element types and numpy dtype names.
/// </summary>
internal static class H5Types
{
    private static readonly string[] NumericTypes = ["bool", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64", "float16", "float32", "float64"];

    public static string CheckNumericType(string dtype)
    {
        if (!NumericTypes.Contains(dtype)) throw new ArgumentException($"Unsupported element type '{dtype}'. Use one of {string.Join(", ", NumericTypes)}.");
        return dtype;
    }

    /// <summary>numpy dtype name of an HDF5 type, or null for strings and unsupported classes.</summary>
    public static string? DTypeName(IH5DataType t) => t.Class switch
    {
        H5DataTypeClass.FloatingPoint => t.Size switch { 2 => "float16", 4 => "float32", 8 => "float64", _ => null },
        H5DataTypeClass.FixedPoint => (t.FixedPoint.IsSigned, t.Size) switch
        {
            (true, 1) => "int8", (true, 2) => "int16", (true, 4) => "int32", (true, 8) => "int64",
            (false, 1) => "uint8", (false, 2) => "uint16", (false, 4) => "uint32", (false, 8) => "uint64",
            _ => null,
        },
        H5DataTypeClass.Enumerated when IsBoolEnum(t) => "bool",
        H5DataTypeClass.BitField when t.Size == 1 => "uint8",
        _ => null,
    };

    public static bool IsString(IH5DataType t) => t.Class == H5DataTypeClass.String
        || (t.Class == H5DataTypeClass.VariableLength && (t.VariableLength.BaseType.Class is H5DataTypeClass.FixedPoint or H5DataTypeClass.String) && t.VariableLength.BaseType.Size == 1);

    private static bool IsBoolEnum(IH5DataType t) => t.Size == 1 && t.Enumeration.BaseType.Class == H5DataTypeClass.FixedPoint;

    /// <summary>Reads a dataset (or a selection of it) into an array of its natural .NET element type.</summary>
    public static Array ReadArray(IH5Dataset ds, Selection? selection = null)
    {
        var t = ds.Type;
        if (IsString(t)) return ds.Read<string[]>(fileSelection: selection);
        return DTypeName(t) switch
        {
            "float32" => ds.Read<float[]>(fileSelection: selection),
            "float64" => ds.Read<double[]>(fileSelection: selection),
            "float16" => ds.Read<ushort[]>(fileSelection: selection).Select(h => (float)BitConverter.UInt16BitsToHalf(h)).ToArray(),
            "int8" => ds.Read<sbyte[]>(fileSelection: selection),
            "int16" => ds.Read<short[]>(fileSelection: selection),
            "int32" => ds.Read<int[]>(fileSelection: selection),
            "int64" => ds.Read<long[]>(fileSelection: selection),
            "uint8" => ds.Read<byte[]>(fileSelection: selection),
            "uint16" => ds.Read<ushort[]>(fileSelection: selection),
            "uint32" => ds.Read<uint[]>(fileSelection: selection),
            "uint64" => ds.Read<ulong[]>(fileSelection: selection),
            "bool" => ds.Read<byte[]>(fileSelection: selection).Select(b => b != 0).ToArray(),
            _ => throw new NotSupportedException($"HDF5 type class {t.Class} (size {t.Size}) of dataset '{ds.Name}' is not supported."),
        };
    }

    /// <summary>Reads a 1-D integer dataset slice as <c>long</c> values (indptr).</summary>
    public static long[] ReadInt64(IH5Dataset ds, ulong start, ulong count)
    {
        var sel = new HyperslabSelection(start, count);
        return DTypeName(ds.Type) switch
        {
            "int64" => ds.Read<long[]>(fileSelection: sel),
            "int32" => Array.ConvertAll(ds.Read<int[]>(fileSelection: sel), x => (long)x),
            "uint32" => Array.ConvertAll(ds.Read<uint[]>(fileSelection: sel), x => (long)x),
            "uint64" => Array.ConvertAll(ds.Read<ulong[]>(fileSelection: sel), x => checked((long)x)),
            "int16" => Array.ConvertAll(ds.Read<short[]>(fileSelection: sel), x => (long)x),
            "uint16" => Array.ConvertAll(ds.Read<ushort[]>(fileSelection: sel), x => (long)x),
            "int8" => Array.ConvertAll(ds.Read<sbyte[]>(fileSelection: sel), x => (long)x),
            "uint8" => Array.ConvertAll(ds.Read<byte[]>(fileSelection: sel), x => (long)x),
            var d => throw new InvalidDataException($"Expected an integer dataset for '{ds.Name}', found {d ?? ds.Type.Class.ToString()}."),
        };
    }

    /// <summary>Reads a 1-D integer dataset slice as <c>int</c> values (column/row indices).</summary>
    public static int[] ReadInt32(IH5Dataset ds, ulong start, ulong count)
    {
        var sel = new HyperslabSelection(start, count);
        return DTypeName(ds.Type) switch
        {
            "int32" => ds.Read<int[]>(fileSelection: sel),
            "int64" => Array.ConvertAll(ds.Read<long[]>(fileSelection: sel), x => checked((int)x)),
            "uint32" => Array.ConvertAll(ds.Read<uint[]>(fileSelection: sel), x => checked((int)x)),
            "uint64" => Array.ConvertAll(ds.Read<ulong[]>(fileSelection: sel), x => checked((int)x)),
            "int16" => Array.ConvertAll(ds.Read<short[]>(fileSelection: sel), x => (int)x),
            "uint16" => Array.ConvertAll(ds.Read<ushort[]>(fileSelection: sel), x => (int)x),
            "int8" => Array.ConvertAll(ds.Read<sbyte[]>(fileSelection: sel), x => (int)x),
            "uint8" => Array.ConvertAll(ds.Read<byte[]>(fileSelection: sel), x => (int)x),
            var d => throw new InvalidDataException($"Expected an integer dataset for '{ds.Name}', found {d ?? ds.Type.Class.ToString()}."),
        };
    }

    /// <summary>Reads a numeric dataset slice as <c>float</c> values (matrix data). The selection may be 1-D or N-D.</summary>
    public static float[] ReadSingle(IH5Dataset ds, Selection? sel)
    {
        return DTypeName(ds.Type) switch
        {
            "float32" => ds.Read<float[]>(fileSelection: sel),
            "float64" => Array.ConvertAll(ds.Read<double[]>(fileSelection: sel), x => (float)x),
            "int32" => Array.ConvertAll(ds.Read<int[]>(fileSelection: sel), x => (float)x),
            "int64" => Array.ConvertAll(ds.Read<long[]>(fileSelection: sel), x => (float)x),
            "uint32" => Array.ConvertAll(ds.Read<uint[]>(fileSelection: sel), x => (float)x),
            "uint64" => Array.ConvertAll(ds.Read<ulong[]>(fileSelection: sel), x => (float)x),
            "int16" => Array.ConvertAll(ds.Read<short[]>(fileSelection: sel), x => (float)x),
            "uint16" => Array.ConvertAll(ds.Read<ushort[]>(fileSelection: sel), x => (float)x),
            "int8" => Array.ConvertAll(ds.Read<sbyte[]>(fileSelection: sel), x => (float)x),
            "uint8" => Array.ConvertAll(ds.Read<byte[]>(fileSelection: sel), x => (float)x),
            "bool" => Array.ConvertAll(ds.Read<byte[]>(fileSelection: sel), x => x != 0 ? 1f : 0f),
            "float16" => Array.ConvertAll(ds.Read<ushort[]>(fileSelection: sel), h => (float)BitConverter.UInt16BitsToHalf(h)),
            var d => throw new InvalidDataException($"Expected a numeric dataset for '{ds.Name}', found {d ?? ds.Type.Class.ToString()}."),
        };
    }

    /// <summary>A numeric array as doubles.</summary>
    public static double[] ToDoubles(Array a) => AnnDataFrame.ToDoubles(a);

    /// <summary>An integer array as ints (categorical codes).</summary>
    public static int[] ToInt32(Array a) => a switch
    {
        int[] i => i,
        sbyte[] s => Array.ConvertAll(s, x => (int)x),
        short[] s => Array.ConvertAll(s, x => (int)x),
        long[] l => Array.ConvertAll(l, x => checked((int)x)),
        byte[] b => Array.ConvertAll(b, x => (int)x),
        ushort[] u => Array.ConvertAll(u, x => (int)x),
        uint[] u => Array.ConvertAll(u, x => checked((int)x)),
        _ => throw new InvalidDataException($"Expected integer codes, found {a.GetType().Name}."),
    };

    // ---- attributes ----

    public static string? StringAttr(IH5Object o, string name)
    {
        if (!o.AttributeExists(name)) return null;
        var a = o.Attribute(name);
        if (!IsString(a.Type)) return null;
        if (a.Space.Type == H5DataspaceType.Scalar) return a.Read<string>();
        var s = a.Read<string[]>();
        return s.Length > 0 ? s[0] : null;
    }

    public static string[] StringArrayAttr(IH5Object o, string name)
    {
        if (!o.AttributeExists(name)) return [];
        var a = o.Attribute(name);
        if (a.Space.Type == H5DataspaceType.Null || (a.Space.Dimensions.Length > 0 && a.Space.Dimensions.Aggregate(1UL, (x, y) => x * y) == 0)) return [];
        if (!IsString(a.Type)) return [];
        if (a.Space.Type == H5DataspaceType.Scalar) return [a.Read<string>()];
        return a.Read<string[]>();
    }

    public static long[]? Int64ArrayAttr(IH5Object o, string name)
    {
        if (!o.AttributeExists(name)) return null;
        var a = o.Attribute(name);
        return DTypeName(a.Type) switch
        {
            "int64" => a.Read<long[]>(),
            "int32" => Array.ConvertAll(a.Read<int[]>(), x => (long)x),
            "uint64" => Array.ConvertAll(a.Read<ulong[]>(), x => (long)x),
            "uint32" => Array.ConvertAll(a.Read<uint[]>(), x => (long)x),
            _ => null,
        };
    }

    public static bool BoolAttr(IH5Object o, string name, bool fallback = false)
    {
        if (!o.AttributeExists(name)) return fallback;
        var a = o.Attribute(name);
        return DTypeName(a.Type) switch
        {
            "bool" or "uint8" => a.Read<byte>() != 0,
            "int8" => a.Read<sbyte>() != 0,
            "int32" => a.Read<int>() != 0,
            "int64" => a.Read<long>() != 0,
            _ => fallback,
        };
    }

    // ---- writing ----

    /// <summary>Converts a numeric array to the .NET array type PureHDF writes for <paramref name="dtype"/>.</summary>
    public static Array ToStorage(ReadOnlySpan<float> values, string dtype)
    {
        var r = dtype switch
        {
            "float32" => (Array)values.ToArray(),
            "float64" => new double[values.Length],
            "int8" => new sbyte[values.Length],
            "int16" => new short[values.Length],
            "int32" => new int[values.Length],
            "int64" => new long[values.Length],
            "uint8" => new byte[values.Length],
            "uint16" => new ushort[values.Length],
            "uint32" => new uint[values.Length],
            "uint64" => new ulong[values.Length],
            "bool" => new H5Bool[values.Length],
            "float16" => new float[values.Length],
            _ => throw new ArgumentException(dtype),
        };
        switch (r)
        {
            case double[] d: for (int i = 0; i < d.Length; i++) d[i] = values[i]; break;
            case sbyte[] d: for (int i = 0; i < d.Length; i++) d[i] = (sbyte)values[i]; break;
            case short[] d: for (int i = 0; i < d.Length; i++) d[i] = (short)values[i]; break;
            case int[] d: for (int i = 0; i < d.Length; i++) d[i] = (int)values[i]; break;
            case long[] d: for (int i = 0; i < d.Length; i++) d[i] = (long)values[i]; break;
            case byte[] d: for (int i = 0; i < d.Length; i++) d[i] = (byte)values[i]; break;
            case ushort[] d: for (int i = 0; i < d.Length; i++) d[i] = (ushort)values[i]; break;
            case uint[] d: for (int i = 0; i < d.Length; i++) d[i] = (uint)values[i]; break;
            case ulong[] d: for (int i = 0; i < d.Length; i++) d[i] = (ulong)values[i]; break;
            case H5Bool[] d: for (int i = 0; i < d.Length; i++) d[i] = values[i] != 0 ? H5Bool.TRUE : H5Bool.FALSE; break;
            case float[] d when dtype == "float16": values.CopyTo(d); break;
        }
        return r;
    }

    /// <summary>Converts doubles to the .NET array type for <paramref name="dtype"/> (dense matrices).</summary>
    public static Array ToStorage(ReadOnlySpan<double> values, string dtype)
    {
        if (dtype == "float64") return values.ToArray();
        var f = new float[values.Length];
        if (dtype is "float32" or "float16")
        {
            for (int i = 0; i < f.Length; i++) f[i] = (float)values[i];
            return f;
        }
        // Integer types: convert through long to keep values above 2^24 exact.
        Array r = dtype switch
        {
            "int8" => new sbyte[values.Length], "int16" => new short[values.Length], "int32" => new int[values.Length],
            "int64" => new long[values.Length], "uint8" => new byte[values.Length], "uint16" => new ushort[values.Length],
            "uint32" => new uint[values.Length], "uint64" => new ulong[values.Length], "bool" => new H5Bool[values.Length],
            _ => throw new ArgumentException(dtype),
        };
        for (int i = 0; i < values.Length; i++)
        {
            double v = values[i];
            switch (r)
            {
                case sbyte[] d: d[i] = (sbyte)v; break;
                case short[] d: d[i] = (short)v; break;
                case int[] d: d[i] = (int)v; break;
                case long[] d: d[i] = (long)v; break;
                case byte[] d: d[i] = (byte)v; break;
                case ushort[] d: d[i] = (ushort)v; break;
                case uint[] d: d[i] = (uint)v; break;
                case ulong[] d: d[i] = (ulong)v; break;
                case H5Bool[] d: d[i] = v != 0 ? H5Bool.TRUE : H5Bool.FALSE; break;
            }
        }
        return r;
    }

    /// <summary>A .NET array ready for PureHDF (bool → the numpy-compatible enum).</summary>
    public static Array ForWrite(Array a) => a is bool[] b ? Array.ConvertAll(b, x => x ? H5Bool.TRUE : H5Bool.FALSE) : a;

    /// <summary>numpy dtype name of a .NET array.</summary>
    public static string DTypeOf(Array a) => a switch
    {
        double[] => "float64", float[] => "float32", int[] => "int32", long[] => "int64", short[] => "int16", sbyte[] => "int8",
        byte[] => "uint8", ushort[] => "uint16", uint[] => "uint32", ulong[] => "uint64", bool[] => "bool",
        _ => throw new ArgumentException($"No numeric dtype for {a.GetType().Name}."),
    };
}
