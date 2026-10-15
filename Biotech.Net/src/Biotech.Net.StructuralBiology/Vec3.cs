using System.Globalization;
using System.Runtime.CompilerServices;

namespace Biotech.Net.StructuralBiology;

/// <summary>
/// A 3-D vector or point in ångström, stored as three doubles. Value type: arrays of <see cref="Vec3"/> are
/// contiguous, so coordinate sets can be passed around as spans without boxing or pointer chasing.
/// </summary>
public readonly record struct Vec3(double X, double Y, double Z)
{
    /// <summary>The origin.</summary>
    public static readonly Vec3 Zero = new(0, 0, 0);

    /// <summary>Unit vectors along the axes.</summary>
    public static readonly Vec3 UnitX = new(1, 0, 0), UnitY = new(0, 1, 0), UnitZ = new(0, 0, 1);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vec3 operator +(Vec3 a, Vec3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vec3 operator -(Vec3 a, Vec3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

    public static Vec3 operator -(Vec3 a) => new(-a.X, -a.Y, -a.Z);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vec3 operator *(Vec3 a, double s) => new(a.X * s, a.Y * s, a.Z * s);

    public static Vec3 operator *(double s, Vec3 a) => new(a.X * s, a.Y * s, a.Z * s);

    public static Vec3 operator /(Vec3 a, double s) => new(a.X / s, a.Y / s, a.Z / s);

    /// <summary>Dot product.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Dot(Vec3 a, Vec3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

    /// <summary>Cross product <c>a × b</c>.</summary>
    public static Vec3 Cross(Vec3 a, Vec3 b) => new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);

    /// <summary>Euclidean distance between two points.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Distance(Vec3 a, Vec3 b) => Math.Sqrt(DistanceSquared(a, b));

    /// <summary>Squared distance (no square root; use for cut-off tests).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double DistanceSquared(Vec3 a, Vec3 b)
    {
        double dx = a.X - b.X, dy = a.Y - b.Y, dz = a.Z - b.Z;
        return dx * dx + dy * dy + dz * dz;
    }

    /// <summary>Angle between two vectors in degrees (0–180).</summary>
    public static double AngleBetween(Vec3 a, Vec3 b)
    {
        double d = a.Length * b.Length;
        if (d == 0) return double.NaN;
        return Math.Acos(Math.Clamp(Dot(a, b) / d, -1, 1)) * (180 / Math.PI);
    }

    /// <summary>Component-wise minimum.</summary>
    public static Vec3 Min(Vec3 a, Vec3 b) => new(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y), Math.Min(a.Z, b.Z));

    /// <summary>Component-wise maximum.</summary>
    public static Vec3 Max(Vec3 a, Vec3 b) => new(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y), Math.Max(a.Z, b.Z));

    /// <summary>Vector length.</summary>
    public double Length => Math.Sqrt(X * X + Y * Y + Z * Z);

    /// <summary>Squared length.</summary>
    public double LengthSquared => X * X + Y * Y + Z * Z;

    /// <summary>Unit vector in the same direction (zero stays zero).</summary>
    public Vec3 Normalized()
    {
        double l = Length;
        return l == 0 ? Zero : this / l;
    }

    public override string ToString() => string.Create(CultureInfo.InvariantCulture, $"({X:F3}, {Y:F3}, {Z:F3})");
}

/// <summary>A 3×3 matrix (row-major), used for rotations.</summary>
public readonly record struct Matrix3(
    double M11, double M12, double M13,
    double M21, double M22, double M23,
    double M31, double M32, double M33)
{
    /// <summary>The identity matrix.</summary>
    public static readonly Matrix3 Identity = new(1, 0, 0, 0, 1, 0, 0, 0, 1);

    /// <summary>Matrix × column vector.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vec3 operator *(Matrix3 m, Vec3 v) => new(
        m.M11 * v.X + m.M12 * v.Y + m.M13 * v.Z,
        m.M21 * v.X + m.M22 * v.Y + m.M23 * v.Z,
        m.M31 * v.X + m.M32 * v.Y + m.M33 * v.Z);

    /// <summary>Matrix product.</summary>
    public static Matrix3 operator *(Matrix3 a, Matrix3 b) => new(
        a.M11 * b.M11 + a.M12 * b.M21 + a.M13 * b.M31, a.M11 * b.M12 + a.M12 * b.M22 + a.M13 * b.M32, a.M11 * b.M13 + a.M12 * b.M23 + a.M13 * b.M33,
        a.M21 * b.M11 + a.M22 * b.M21 + a.M23 * b.M31, a.M21 * b.M12 + a.M22 * b.M22 + a.M23 * b.M32, a.M21 * b.M13 + a.M22 * b.M23 + a.M23 * b.M33,
        a.M31 * b.M11 + a.M32 * b.M21 + a.M33 * b.M31, a.M31 * b.M12 + a.M32 * b.M22 + a.M33 * b.M32, a.M31 * b.M13 + a.M32 * b.M23 + a.M33 * b.M33);

    /// <summary>Transpose (the inverse of a rotation).</summary>
    public Matrix3 Transpose() => new(M11, M21, M31, M12, M22, M32, M13, M23, M33);

    /// <summary>Determinant (+1 for a proper rotation, −1 for a reflection).</summary>
    public double Determinant =>
        M11 * (M22 * M33 - M23 * M32) - M12 * (M21 * M33 - M23 * M31) + M13 * (M21 * M32 - M22 * M31);

    /// <summary>Rotation of <paramref name="degrees"/> about <paramref name="axis"/> (right-hand rule).</summary>
    public static Matrix3 FromAxisAngle(Vec3 axis, double degrees)
    {
        var u = axis.Normalized();
        double t = degrees * Math.PI / 180, c = Math.Cos(t), s = Math.Sin(t), k = 1 - c;
        return new(
            c + u.X * u.X * k, u.X * u.Y * k - u.Z * s, u.X * u.Z * k + u.Y * s,
            u.Y * u.X * k + u.Z * s, c + u.Y * u.Y * k, u.Y * u.Z * k - u.X * s,
            u.Z * u.X * k - u.Y * s, u.Z * u.Y * k + u.X * s, c + u.Z * u.Z * k);
    }

    /// <summary>Rotation matrix of the unit quaternion (w, x, y, z); the quaternion is normalised first.</summary>
    public static Matrix3 FromQuaternion(double w, double x, double y, double z)
    {
        double n = Math.Sqrt(w * w + x * x + y * y + z * z);
        w /= n; x /= n; y /= n; z /= n;
        return new(
            w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y),
            2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x),
            2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z);
    }

    /// <summary>Rotation angle in degrees (0–180).</summary>
    public double RotationAngle => Math.Acos(Math.Clamp((M11 + M22 + M33 - 1) / 2, -1, 1)) * (180 / Math.PI);
}
