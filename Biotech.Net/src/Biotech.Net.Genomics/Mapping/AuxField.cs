using System.Buffers.Binary;
using System.Globalization;
using System.Text;
using Biotech.Net.Genomics.IO;

namespace Biotech.Net.Genomics.Mapping;

/// <summary>
/// One optional field (auxiliary tag) of an alignment, typed as in BAM: <c>A</c> (char), <c>c C s S i I</c>
/// (integers, held as <see cref="long"/>), <c>f</c> (float), <c>Z</c> (string), <c>H</c> (hex string) or <c>B</c>
/// (numeric array: <see cref="sbyte"/>[], <see cref="byte"/>[], <see cref="short"/>[], <see cref="ushort"/>[],
/// <see cref="int"/>[], <see cref="uint"/>[] or <see cref="float"/>[]).
/// </summary>
public readonly struct AuxField : IEquatable<AuxField>
{
    /// <summary>Creates a field; <paramref name="value"/> must match <paramref name="type"/> (see the type remarks).</summary>
    public AuxField(string tag, char type, object value)
    {
        if (tag is null || tag.Length != 2) throw new ArgumentException("A tag has exactly two characters.", nameof(tag));
        Tag = tag;
        Type = type;
        Value = type switch
        {
            'A' => value is char ? value : throw new ArgumentException("Type A needs a char value."),
            'c' or 'C' or 's' or 'S' or 'i' or 'I' => Convert.ToInt64(value, CultureInfo.InvariantCulture),
            'f' => Convert.ToSingle(value, CultureInfo.InvariantCulture),
            'Z' or 'H' => value as string ?? throw new ArgumentException($"Type {type} needs a string value."),
            'B' => ArrayCode(value) != '\0' ? value
                : throw new ArgumentException("Type B needs a numeric array value."),
            _ => throw new ArgumentException($"Unknown tag type '{type}'.", nameof(type)),
        };
    }

    /// <summary>Two-character tag name (<c>NM</c>, <c>MD</c>, ...).</summary>
    public string Tag { get; }

    /// <summary>BAM type code.</summary>
    public char Type { get; }

    /// <summary>The value (see the type remarks).</summary>
    public object Value { get; }

    /// <summary>The element type of a <c>B</c> array (<c>c C s S i I f</c>), otherwise <c>'\0'</c>.</summary>
    public char ArrayType => ArrayCode(Value);

    /// <summary>The B-array element code of <paramref name="v"/> by exact runtime type (the CLR lets byte[] pass as sbyte[] in 'is' tests).</summary>
    internal static char ArrayCode(object? v)
    {
        var t = v?.GetType();
        if (t == typeof(sbyte[])) return 'c';
        if (t == typeof(byte[])) return 'C';
        if (t == typeof(short[])) return 's';
        if (t == typeof(ushort[])) return 'S';
        if (t == typeof(int[])) return 'i';
        if (t == typeof(uint[])) return 'I';
        if (t == typeof(float[])) return 'f';
        return '\0';
    }

    /// <summary>An integer field stored with the smallest BAM type that holds it (htslib's choice).</summary>
    public static AuxField Int(string tag, long value) => new(tag, SmallestIntType(value), value);

    /// <summary>A string (<c>Z</c>) field.</summary>
    public static AuxField String(string tag, string value) => new(tag, 'Z', value);

    /// <summary>A character (<c>A</c>) field.</summary>
    public static AuxField Char(string tag, char value) => new(tag, 'A', value);

    /// <summary>A float (<c>f</c>) field.</summary>
    public static AuxField Float(string tag, float value) => new(tag, 'f', value);

    /// <summary>A hex-string (<c>H</c>) field.</summary>
    public static AuxField Hex(string tag, string value) => new(tag, 'H', value);

    /// <summary>A numeric array (<c>B</c>) field.</summary>
    public static AuxField Array(string tag, Array values) => new(tag, 'B', values);

    internal static char SmallestIntType(long x)
    {
        if (x < 0)
        {
            if (x >= sbyte.MinValue) return 'c';
            if (x >= short.MinValue) return 's';
            if (x >= int.MinValue) return 'i';
        }
        else
        {
            if (x <= byte.MaxValue) return 'C';
            if (x <= ushort.MaxValue) return 'S';
            if (x <= uint.MaxValue) return 'I';
        }
        throw new ArgumentOutOfRangeException(nameof(x), $"Integer tag value {x} does not fit in 32 bits.");
    }

    /// <summary>Whether this is an integer field.</summary>
    public bool IsInteger => Type is 'c' or 'C' or 's' or 'S' or 'i' or 'I';

    /// <summary>The integer value (integer types only).</summary>
    public long AsInt64() => IsInteger ? (long)Value : throw new InvalidOperationException($"Tag {Tag} has type {Type}, not an integer.");

    /// <summary>The integer value as <see cref="int"/> (checked).</summary>
    public int AsInt32() => checked((int)AsInt64());

    /// <summary>The numeric value of an integer or float field.</summary>
    public double AsDouble() => Type == 'f' ? (float)Value : AsInt64();

    /// <summary>The string value of <c>Z</c>/<c>H</c> fields (or the character of an <c>A</c> field).</summary>
    public string AsString() => Value switch
    {
        string s => s,
        char c => c.ToString(),
        _ => throw new InvalidOperationException($"Tag {Tag} has type {Type}, not a string."),
    };

    /// <summary>The array of a <c>B</c> field.</summary>
    public T[] AsArray<T>() => Value as T[] ?? throw new InvalidOperationException($"Tag {Tag} is not a {typeof(T).Name} array.");

    /// <summary>The field in SAM text form (<c>NM:i:2</c>, <c>XB:B:c,1,-2</c>, ...).</summary>
    public override string ToString()
    {
        var sb = new StringBuilder();
        AppendSam(sb);
        return sb.ToString();
    }

    internal void AppendSam(StringBuilder sb)
    {
        sb.Append(Tag).Append(':');
        switch (Type)
        {
            case 'A': sb.Append("A:").Append((char)Value); break;
            case 'c' or 'C' or 's' or 'S' or 'i' or 'I': sb.Append("i:").Append(((long)Value).ToString(CultureInfo.InvariantCulture)); break;
            case 'f': sb.Append("f:").Append(SamText.FormatG((float)Value)); break;
            case 'Z' or 'H': sb.Append(Type).Append(':').Append((string)Value); break;
            case 'B':
                sb.Append("B:").Append(ArrayType);
                switch (ArrayType)
                {
                    case 'c': foreach (var v in (sbyte[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 'C': foreach (var v in (byte[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 's': foreach (var v in (short[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 'S': foreach (var v in (ushort[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 'i': foreach (var v in (int[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 'I': foreach (var v in (uint[])Value) sb.Append(',').Append(v.ToString(CultureInfo.InvariantCulture)); break;
                    case 'f': foreach (var v in (float[])Value) sb.Append(',').Append(SamText.FormatG(v)); break;
                }
                break;
        }
    }

    /// <summary>Parses a SAM text field (<c>TG:TYPE:VALUE</c>); integers get the smallest fitting BAM type.</summary>
    public static AuxField Parse(string field)
    {
        var buf = new ByteBuffer(32);
        AuxCodec.AppendSamField(buf, field);
        return AuxCodec.Decode(buf.Span)[0];
    }

    public bool Equals(AuxField other)
    {
        if (Tag != other.Tag || Type != other.Type) return false;
        if (Value is Array a && other.Value is Array b)
        {
            if (a.GetType() != b.GetType() || a.Length != b.Length) return false;
            for (int i = 0; i < a.Length; i++) if (!Equals(a.GetValue(i), b.GetValue(i))) return false;
            return true;
        }
        return Equals(Value, other.Value);
    }

    public override bool Equals(object? obj) => obj is AuxField f && Equals(f);
    public override int GetHashCode() => HashCode.Combine(Tag, Type);
    public static bool operator ==(AuxField a, AuxField b) => a.Equals(b);
    public static bool operator !=(AuxField a, AuxField b) => !a.Equals(b);
}

/// <summary>Encoding of optional fields in the BAM binary layout, and their SAM text form.</summary>
internal static class AuxCodec
{
    public static int ElementSize(char t) => t switch
    {
        'A' or 'c' or 'C' => 1,
        's' or 'S' => 2,
        'i' or 'I' or 'f' => 4,
        _ => throw new InvalidDataException($"Unknown aux array element type '{t}'."),
    };

    /// <summary>Byte length of the field starting at <paramref name="pos"/> (tag and type included).</summary>
    public static int FieldLength(ReadOnlySpan<byte> aux, int pos)
    {
        if (pos + 3 > aux.Length) throw new InvalidDataException("Truncated aux field.");
        char t = (char)aux[pos + 2];
        switch (t)
        {
            case 'A' or 'c' or 'C': return 4;
            case 's' or 'S': return 5;
            case 'i' or 'I' or 'f': return 7;
            case 'd': return 11;
            case 'Z' or 'H':
            {
                int nul = aux[(pos + 3)..].IndexOf((byte)0);
                if (nul < 0) throw new InvalidDataException("Unterminated Z/H aux field.");
                return 3 + nul + 1;
            }
            case 'B':
            {
                if (pos + 8 > aux.Length) throw new InvalidDataException("Truncated B aux field.");
                int n = BinaryPrimitives.ReadInt32LittleEndian(aux[(pos + 4)..]);
                return 8 + n * ElementSize((char)aux[pos + 3]);
            }
            default: throw new InvalidDataException($"Unknown aux field type '{t}'.");
        }
    }

    /// <summary>Offset of tag <paramref name="tag"/> in <paramref name="aux"/>, or -1.</summary>
    public static int Find(ReadOnlySpan<byte> aux, ReadOnlySpan<char> tag)
    {
        for (int p = 0; p + 3 <= aux.Length;)
        {
            if (aux[p] == tag[0] && aux[p + 1] == tag[1]) return p;
            p += FieldLength(aux, p);
        }
        return -1;
    }

    public static List<AuxField> Decode(ReadOnlySpan<byte> aux)
    {
        var list = new List<AuxField>();
        for (int p = 0; p + 3 <= aux.Length;)
        {
            int len = FieldLength(aux, p);
            list.Add(DecodeOne(aux.Slice(p, len)));
            p += len;
        }
        return list;
    }

    public static AuxField DecodeOne(ReadOnlySpan<byte> f)
    {
        string tag = new([(char)f[0], (char)f[1]]);
        char t = (char)f[2];
        var v = f[3..];
        return t switch
        {
            'A' => new AuxField(tag, 'A', (char)v[0]),
            'c' => new AuxField(tag, 'c', (long)(sbyte)v[0]),
            'C' => new AuxField(tag, 'C', (long)v[0]),
            's' => new AuxField(tag, 's', (long)BinaryPrimitives.ReadInt16LittleEndian(v)),
            'S' => new AuxField(tag, 'S', (long)BinaryPrimitives.ReadUInt16LittleEndian(v)),
            'i' => new AuxField(tag, 'i', (long)BinaryPrimitives.ReadInt32LittleEndian(v)),
            'I' => new AuxField(tag, 'I', (long)BinaryPrimitives.ReadUInt32LittleEndian(v)),
            'f' => new AuxField(tag, 'f', BinaryPrimitives.ReadSingleLittleEndian(v)),
            'd' => new AuxField(tag, 'f', (float)BinaryPrimitives.ReadDoubleLittleEndian(v)),
            'Z' or 'H' => new AuxField(tag, t, Encoding.UTF8.GetString(v[..^1])),
            'B' => new AuxField(tag, 'B', DecodeArray((char)v[0], BinaryPrimitives.ReadInt32LittleEndian(v[1..]), v[5..])),
            _ => throw new InvalidDataException($"Unknown aux field type '{t}'."),
        };
    }

    private static Array DecodeArray(char sub, int n, ReadOnlySpan<byte> d)
    {
        switch (sub)
        {
            case 'c': { var a = new sbyte[n]; for (int i = 0; i < n; i++) a[i] = (sbyte)d[i]; return a; }
            case 'C' or 'A': return d[..n].ToArray();
            case 's': { var a = new short[n]; for (int i = 0; i < n; i++) a[i] = BinaryPrimitives.ReadInt16LittleEndian(d[(2 * i)..]); return a; }
            case 'S': { var a = new ushort[n]; for (int i = 0; i < n; i++) a[i] = BinaryPrimitives.ReadUInt16LittleEndian(d[(2 * i)..]); return a; }
            case 'i': { var a = new int[n]; for (int i = 0; i < n; i++) a[i] = BinaryPrimitives.ReadInt32LittleEndian(d[(4 * i)..]); return a; }
            case 'I': { var a = new uint[n]; for (int i = 0; i < n; i++) a[i] = BinaryPrimitives.ReadUInt32LittleEndian(d[(4 * i)..]); return a; }
            case 'f': { var a = new float[n]; for (int i = 0; i < n; i++) a[i] = BinaryPrimitives.ReadSingleLittleEndian(d[(4 * i)..]); return a; }
            default: throw new InvalidDataException($"Unknown aux array element type '{sub}'.");
        }
    }

    public static void Encode(ByteBuffer buf, AuxField f)
    {
        buf.Append((byte)f.Tag[0]);
        buf.Append((byte)f.Tag[1]);
        buf.Append((byte)f.Type);
        Span<byte> tmp = stackalloc byte[8];
        switch (f.Type)
        {
            case 'A': buf.Append((byte)(char)f.Value); break;
            case 'c': buf.Append((byte)(sbyte)(long)f.Value); break;
            case 'C': buf.Append((byte)(long)f.Value); break;
            case 's': BinaryPrimitives.WriteInt16LittleEndian(tmp, checked((short)(long)f.Value)); buf.Append(tmp[..2]); break;
            case 'S': BinaryPrimitives.WriteUInt16LittleEndian(tmp, checked((ushort)(long)f.Value)); buf.Append(tmp[..2]); break;
            case 'i': BinaryPrimitives.WriteInt32LittleEndian(tmp, checked((int)(long)f.Value)); buf.Append(tmp[..4]); break;
            case 'I': BinaryPrimitives.WriteUInt32LittleEndian(tmp, checked((uint)(long)f.Value)); buf.Append(tmp[..4]); break;
            case 'f': BinaryPrimitives.WriteSingleLittleEndian(tmp, (float)f.Value); buf.Append(tmp[..4]); break;
            case 'Z' or 'H': buf.Append(Encoding.UTF8.GetBytes((string)f.Value)); buf.Append(0); break;
            case 'B':
            {
                var a = (Array)f.Value;
                buf.Append((byte)f.ArrayType);
                BinaryPrimitives.WriteInt32LittleEndian(tmp, a.Length);
                buf.Append(tmp[..4]);
                switch (f.ArrayType)
                {
                    case 'c': foreach (var v in (sbyte[])a) buf.Append((byte)v); break;
                    case 'C': buf.Append((byte[])a); break;
                    case 's': foreach (var v in (short[])a) { BinaryPrimitives.WriteInt16LittleEndian(tmp, v); buf.Append(tmp[..2]); } break;
                    case 'S': foreach (var v in (ushort[])a) { BinaryPrimitives.WriteUInt16LittleEndian(tmp, v); buf.Append(tmp[..2]); } break;
                    case 'i': foreach (var v in (int[])a) { BinaryPrimitives.WriteInt32LittleEndian(tmp, v); buf.Append(tmp[..4]); } break;
                    case 'I': foreach (var v in (uint[])a) { BinaryPrimitives.WriteUInt32LittleEndian(tmp, v); buf.Append(tmp[..4]); } break;
                    case 'f': foreach (var v in (float[])a) { BinaryPrimitives.WriteSingleLittleEndian(tmp, v); buf.Append(tmp[..4]); } break;
                }
                break;
            }
        }
    }

    /// <summary>Appends the BAM encoding of one SAM text field (<c>TG:T:VALUE</c>), choosing integer types like htslib.</summary>
    public static void AppendSamField(ByteBuffer buf, ReadOnlySpan<char> f)
    {
        if (f.Length < 5 || f[2] != ':' || f[4] != ':')
            throw new InvalidDataException($"Malformed SAM optional field '{f}'.");
        char t = f[3];
        var v = f[5..];
        var inv = CultureInfo.InvariantCulture;
        string tag = new(f[..2]);
        switch (t)
        {
            case 'A':
                if (v.Length != 1) throw new InvalidDataException($"Type A field '{f}' must hold one character.");
                Encode(buf, new AuxField(tag, 'A', v[0]));
                break;
            case 'i':
            {
                if (!long.TryParse(v, NumberStyles.AllowLeadingSign, inv, out long x)) throw new InvalidDataException($"Invalid integer in '{f}'.");
                Encode(buf, AuxField.Int(tag, x));
                break;
            }
            case 'f':
                if (!float.TryParse(v, NumberStyles.Float, inv, out float fl)) throw new InvalidDataException($"Invalid float in '{f}'.");
                Encode(buf, new AuxField(tag, 'f', fl));
                break;
            case 'Z':
            case 'H':
                Encode(buf, new AuxField(tag, t, v.ToString()));
                break;
            case 'B':
            {
                if (v.IsEmpty) throw new InvalidDataException($"Type B field '{f}' lacks an element type.");
                char sub = v[0];
                var rest = v.Length > 1 ? v[1..] : default;
                if (!rest.IsEmpty && rest[0] == ',') rest = rest[1..];
                int n = rest.IsEmpty ? 0 : rest.Count(',') + 1;
                Array arr = sub switch
                {
                    'c' => new sbyte[n], 'C' => new byte[n], 's' => new short[n], 'S' => new ushort[n],
                    'i' => new int[n], 'I' => new uint[n], 'f' => new float[n],
                    _ => throw new InvalidDataException($"Unknown B array element type in '{f}'."),
                };
                for (int i = 0; i < n; i++)
                {
                    int c = rest.IndexOf(',');
                    var item = c < 0 ? rest : rest[..c];
                    rest = c < 0 ? default : rest[(c + 1)..];
                    switch (sub)
                    {
                        case 'c': ((sbyte[])arr)[i] = sbyte.Parse(item, NumberStyles.AllowLeadingSign, inv); break;
                        case 'C': ((byte[])arr)[i] = byte.Parse(item, NumberStyles.None, inv); break;
                        case 's': ((short[])arr)[i] = short.Parse(item, NumberStyles.AllowLeadingSign, inv); break;
                        case 'S': ((ushort[])arr)[i] = ushort.Parse(item, NumberStyles.None, inv); break;
                        case 'i': ((int[])arr)[i] = int.Parse(item, NumberStyles.AllowLeadingSign, inv); break;
                        case 'I': ((uint[])arr)[i] = uint.Parse(item, NumberStyles.None, inv); break;
                        case 'f': ((float[])arr)[i] = float.Parse(item, NumberStyles.Float, inv); break;
                    }
                }
                Encode(buf, new AuxField(tag, 'B', arr));
                break;
            }
            default:
                throw new InvalidDataException($"Unknown SAM optional field type '{t}' in '{f}'.");
        }
    }

    /// <summary>Appends the SAM text of the encoded fields (each preceded by a tab).</summary>
    public static void AppendSamText(StringBuilder sb, ReadOnlySpan<byte> aux)
    {
        for (int p = 0; p + 3 <= aux.Length;)
        {
            int len = FieldLength(aux, p);
            var f = aux.Slice(p, len);
            sb.Append('\t');
            char t = (char)f[2];
            if (t is 'Z' or 'H')
            {
                sb.Append((char)f[0]).Append((char)f[1]).Append(':').Append(t).Append(':');
                sb.Append(Encoding.UTF8.GetString(f[3..^1]));
            }
            else DecodeOne(f).AppendSam(sb);
            p += len;
        }
    }
}

/// <summary>Number formatting shared by the SAM writers.</summary>
public static class SamText
{
    /// <summary>Formats like C's <c>printf("%g")</c> (6 significant digits), as htslib does for <c>f</c> tags.</summary>
    public static string FormatG(double d)
    {
        if (double.IsNaN(d)) return "nan";
        if (double.IsPositiveInfinity(d)) return "inf";
        if (double.IsNegativeInfinity(d)) return "-inf";
        if (d == 0) return double.IsNegative(d) ? "-0" : "0";
        var inv = CultureInfo.InvariantCulture;
        string e = d.ToString("E5", inv); // [-]d.dddddE+xxx
        int ei = e.IndexOf('E');
        int exp = int.Parse(e.AsSpan(ei + 1), NumberStyles.AllowLeadingSign, inv);
        if (exp < -4 || exp >= 6)
        {
            string m = e[..ei];
            if (m.Contains('.')) m = m.TrimEnd('0').TrimEnd('.');
            return m + "e" + (exp < 0 ? "-" : "+") + Math.Abs(exp).ToString("00", inv);
        }
        string f = d.ToString("F" + (5 - exp), inv);
        if (f.Contains('.')) f = f.TrimEnd('0').TrimEnd('.');
        return f;
    }
}
