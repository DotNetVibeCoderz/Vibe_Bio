using System.Runtime.CompilerServices;

namespace Biotech.Net.Chemistry;

/// <summary>
/// Geometry primitives with analytic gradients for force-field terms. Coordinates are flat arrays
/// <c>x[dim * atom + k]</c>; <c>dim</c> is 3 for ordinary force fields and 4 for distance-geometry refinement
/// (where only pair distances use the fourth coordinate).
/// </summary>
internal static class FFMath
{
    public const double Deg2Rad = Math.PI / 180.0;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Distance2(ReadOnlySpan<double> x, int i, int j, int dim)
    {
        int a = dim * i, b = dim * j;
        double s = 0;
        for (int k = 0; k < dim; k++) { double d = x[a + k] - x[b + k]; s += d * d; }
        return s;
    }

    /// <summary>Adds <c>dEdr · ∂r/∂x</c> for the distance r between atoms i and j.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void AddDistanceGradient(ReadOnlySpan<double> x, Span<double> g, int i, int j, int dim, double r, double dEdr)
    {
        if (r < 1e-12) return;
        double f = dEdr / r;
        int a = dim * i, b = dim * j;
        for (int k = 0; k < dim; k++)
        {
            double d = (x[a + k] - x[b + k]) * f;
            g[a + k] += d;
            g[b + k] -= d;
        }
    }

    /// <summary>Cosine of the angle i–j–k (first three coordinates).</summary>
    public static double CosAngle(ReadOnlySpan<double> x, int i, int j, int k, int dim)
    {
        Vec(x, i, j, dim, out double ax, out double ay, out double az);
        Vec(x, k, j, dim, out double bx, out double by, out double bz);
        double la = Math.Sqrt(ax * ax + ay * ay + az * az), lb = Math.Sqrt(bx * bx + by * by + bz * bz);
        if (la < 1e-12 || lb < 1e-12) return 1;
        return Math.Clamp((ax * bx + ay * by + az * bz) / (la * lb), -1, 1);
    }

    /// <summary>Adds <c>dEdcos · ∂cosθ/∂x</c> for the angle i–j–k.</summary>
    public static void AddCosAngleGradient(ReadOnlySpan<double> x, Span<double> g, int i, int j, int k, int dim, double dEdcos)
    {
        Vec(x, i, j, dim, out double ax, out double ay, out double az);
        Vec(x, k, j, dim, out double bx, out double by, out double bz);
        double la = Math.Sqrt(ax * ax + ay * ay + az * az), lb = Math.Sqrt(bx * bx + by * by + bz * bz);
        if (la < 1e-12 || lb < 1e-12) return;
        ax /= la; ay /= la; az /= la; bx /= lb; by /= lb; bz /= lb;
        double c = ax * bx + ay * by + az * bz;
        double gix = (bx - c * ax) / la, giy = (by - c * ay) / la, giz = (bz - c * az) / la;
        double gkx = (ax - c * bx) / lb, gky = (ay - c * by) / lb, gkz = (az - c * bz) / lb;
        int pi = dim * i, pj = dim * j, pk = dim * k;
        g[pi] += dEdcos * gix; g[pi + 1] += dEdcos * giy; g[pi + 2] += dEdcos * giz;
        g[pk] += dEdcos * gkx; g[pk + 1] += dEdcos * gky; g[pk + 2] += dEdcos * gkz;
        g[pj] -= dEdcos * (gix + gkx); g[pj + 1] -= dEdcos * (giy + gky); g[pj + 2] -= dEdcos * (giz + gkz);
    }

    /// <summary>Dihedral angle i–j–k–l in radians, in (−π, π].</summary>
    public static double Dihedral(ReadOnlySpan<double> x, int i, int j, int k, int l, int dim)
    {
        Vec(x, i, j, dim, out double fx, out double fy, out double fz);
        Vec(x, j, k, dim, out double gx, out double gy, out double gz);
        Vec(x, l, k, dim, out double hx, out double hy, out double hz);
        Cross(fx, fy, fz, gx, gy, gz, out double ax, out double ay, out double az);
        Cross(hx, hy, hz, gx, gy, gz, out double bx, out double by, out double bz);
        double gl = Math.Sqrt(gx * gx + gy * gy + gz * gz);
        double cosv = ax * bx + ay * by + az * bz;
        Cross(bx, by, bz, ax, ay, az, out double cx, out double cy, out double cz);
        double sinv = gl < 1e-12 ? 0 : (cx * gx + cy * gy + cz * gz) / gl;
        return Math.Atan2(sinv, cosv);
    }

    /// <summary>Adds <c>dEdphi · ∂φ/∂x</c> for the dihedral i–j–k–l (Blondel &amp; Karplus, J. Comput. Chem. 1996, 17, 1132).</summary>
    public static void AddDihedralGradient(ReadOnlySpan<double> x, Span<double> g, int i, int j, int k, int l, int dim, double dEdphi)
    {
        Vec(x, i, j, dim, out double fx, out double fy, out double fz);
        Vec(x, j, k, dim, out double gx, out double gy, out double gz);
        Vec(x, l, k, dim, out double hx, out double hy, out double hz);
        Cross(fx, fy, fz, gx, gy, gz, out double ax, out double ay, out double az);
        Cross(hx, hy, hz, gx, gy, gz, out double bx, out double by, out double bz);
        double a2 = ax * ax + ay * ay + az * az, b2 = bx * bx + by * by + bz * bz;
        double gl = Math.Sqrt(gx * gx + gy * gy + gz * gz);
        if (a2 < 1e-16 || b2 < 1e-16 || gl < 1e-12) return;
        double fg = fx * gx + fy * gy + fz * gz, hg = hx * gx + hy * gy + hz * gz;
        double sa = -gl / a2, sb = gl / b2;
        double ta = fg / (a2 * gl), tb = hg / (b2 * gl);
        int pi = dim * i, pj = dim * j, pk = dim * k, pl = dim * l;
        // ∂φ/∂ri = −|G|/A² A ; ∂φ/∂rl = |G|/B² B
        double dix = sa * ax, diy = sa * ay, diz = sa * az;
        double dlx = sb * bx, dly = sb * by, dlz = sb * bz;
        // ∂φ/∂rj = |G|/A² A + (F·G)/(A²|G|) A − (H·G)/(B²|G|) B
        double djx = -dix + ta * ax - tb * bx, djy = -diy + ta * ay - tb * by, djz = -diz + ta * az - tb * bz;
        // ∂φ/∂rk = (H·G)/(B²|G|) B − (F·G)/(A²|G|) A − |G|/B² B
        double dkx = tb * bx - ta * ax - dlx, dky = tb * by - ta * ay - dly, dkz = tb * bz - ta * az - dlz;
        g[pi] += dEdphi * dix; g[pi + 1] += dEdphi * diy; g[pi + 2] += dEdphi * diz;
        g[pj] += dEdphi * djx; g[pj + 1] += dEdphi * djy; g[pj + 2] += dEdphi * djz;
        g[pk] += dEdphi * dkx; g[pk + 1] += dEdphi * dky; g[pk + 2] += dEdphi * dkz;
        g[pl] += dEdphi * dlx; g[pl + 1] += dEdphi * dly; g[pl + 2] += dEdphi * dlz;
    }

    /// <summary>
    /// Wilson out-of-plane quantity used by UFF inversions: cos Y, where Y is the angle between the bond j→l and
    /// the normal of the plane (i, j, k). Returns 0 for degenerate geometry.
    /// </summary>
    public static double InversionCosY(ReadOnlySpan<double> x, int i, int j, int k, int l, int dim)
    {
        Vec(x, i, j, dim, out double ax, out double ay, out double az);
        Vec(x, k, j, dim, out double bx, out double by, out double bz);
        Vec(x, l, j, dim, out double cx, out double cy, out double cz);
        double la = ax * ax + ay * ay + az * az, lb = bx * bx + by * by + bz * bz, lc = cx * cx + cy * cy + cz * cz;
        if (la < 1e-16 || lb < 1e-16 || lc < 1e-16) return 0;
        Cross(ax, ay, az, bx, by, bz, out double nx, out double ny, out double nz);
        double s = 1 / (Math.Sqrt(la) * Math.Sqrt(lb));
        nx *= s; ny *= s; nz *= s;
        double ln = nx * nx + ny * ny + nz * nz;
        if (ln < 1e-16) return 0;
        return (nx * cx + ny * cy + nz * cz) / (Math.Sqrt(lc) * Math.Sqrt(ln));
    }

    /// <summary>
    /// Adds the gradient of an inversion energy E(Y) given <c>dE/dW</c> where W = π/2 − Y is the out-of-plane
    /// angle (port of the analytic expressions in RDKit's UFF <c>InversionContrib::getGrad</c>).
    /// </summary>
    public static void AddInversionGradient(ReadOnlySpan<double> x, Span<double> g, int i, int j, int k, int l, int dim, double forceConstant, double c1, double c2)
    {
        Vec(x, i, j, dim, out double rjix, out double rjiy, out double rjiz);
        Vec(x, k, j, dim, out double rjkx, out double rjky, out double rjkz);
        Vec(x, l, j, dim, out double rjlx, out double rjly, out double rjlz);
        double dji = Math.Sqrt(rjix * rjix + rjiy * rjiy + rjiz * rjiz);
        double djk = Math.Sqrt(rjkx * rjkx + rjky * rjky + rjkz * rjkz);
        double djl = Math.Sqrt(rjlx * rjlx + rjly * rjly + rjlz * rjlz);
        if (dji < 1e-10 || djk < 1e-10 || djl < 1e-10) return;
        rjix /= dji; rjiy /= dji; rjiz /= dji;
        rjkx /= djk; rjky /= djk; rjkz /= djk;
        rjlx /= djl; rjly /= djl; rjlz /= djl;
        Cross(-rjix, -rjiy, -rjiz, rjkx, rjky, rjkz, out double nx, out double ny, out double nz);
        double nl = Math.Sqrt(nx * nx + ny * ny + nz * nz);
        if (nl < 1e-12) return;
        nx /= nl; ny /= nl; nz /= nl;
        double cosY = Math.Clamp(nx * rjlx + ny * rjly + nz * rjlz, -1, 1);
        double sinY = Math.Max(Math.Sqrt(1 - cosY * cosY), 1e-8);
        double cosTheta = Math.Clamp(rjix * rjkx + rjiy * rjky + rjiz * rjkz, -1, 1);
        double sinThetaSq = 1 - cosTheta * cosTheta;
        double sinTheta = Math.Max(Math.Sqrt(sinThetaSq), 1e-8);
        sinThetaSq = Math.Max(sinThetaSq, 1e-16);
        double dEdW = -forceConstant * (c1 * cosY + 4.0 * c2 * cosY * sinY);
        Cross(rjlx, rjly, rjlz, rjkx, rjky, rjkz, out double t1x, out double t1y, out double t1z);
        Cross(rjix, rjiy, rjiz, rjlx, rjly, rjlz, out double t2x, out double t2y, out double t2z);
        Cross(rjkx, rjky, rjkz, rjix, rjiy, rjiz, out double t3x, out double t3y, out double t3z);
        double term1 = sinY * sinTheta;
        double term2 = cosY / (sinY * sinThetaSq);
        double g1x = (t1x / term1 - (rjix - rjkx * cosTheta) * term2) / dji;
        double g1y = (t1y / term1 - (rjiy - rjky * cosTheta) * term2) / dji;
        double g1z = (t1z / term1 - (rjiz - rjkz * cosTheta) * term2) / dji;
        double g3x = (t2x / term1 - (rjkx - rjix * cosTheta) * term2) / djk;
        double g3y = (t2y / term1 - (rjky - rjiy * cosTheta) * term2) / djk;
        double g3z = (t2z / term1 - (rjkz - rjiz * cosTheta) * term2) / djk;
        double g4x = (t3x / term1 - rjlx * cosY / sinY) / djl;
        double g4y = (t3y / term1 - rjly * cosY / sinY) / djl;
        double g4z = (t3z / term1 - rjlz * cosY / sinY) / djl;
        int pi = dim * i, pj = dim * j, pk = dim * k, pl = dim * l;
        g[pi] += dEdW * g1x; g[pi + 1] += dEdW * g1y; g[pi + 2] += dEdW * g1z;
        g[pk] += dEdW * g3x; g[pk + 1] += dEdW * g3y; g[pk + 2] += dEdW * g3z;
        g[pl] += dEdW * g4x; g[pl + 1] += dEdW * g4y; g[pl + 2] += dEdW * g4z;
        g[pj] -= dEdW * (g1x + g3x + g4x); g[pj + 1] -= dEdW * (g1y + g3y + g4y); g[pj + 2] -= dEdW * (g1z + g3z + g4z);
    }

    /// <summary>Signed volume (p1−p4)·((p2−p4)×(p3−p4)) over the first three coordinates.</summary>
    public static double ChiralVolume(ReadOnlySpan<double> x, int i1, int i2, int i3, int i4, int dim)
    {
        Vec(x, i1, i4, dim, out double ax, out double ay, out double az);
        Vec(x, i2, i4, dim, out double bx, out double by, out double bz);
        Vec(x, i3, i4, dim, out double cx, out double cy, out double cz);
        Cross(bx, by, bz, cx, cy, cz, out double nx, out double ny, out double nz);
        return ax * nx + ay * ny + az * nz;
    }

    /// <summary>Adds <c>dEdV · ∂V/∂x</c> for <see cref="ChiralVolume"/>.</summary>
    public static void AddChiralVolumeGradient(ReadOnlySpan<double> x, Span<double> g, int i1, int i2, int i3, int i4, int dim, double dEdV)
    {
        Vec(x, i1, i4, dim, out double ax, out double ay, out double az);
        Vec(x, i2, i4, dim, out double bx, out double by, out double bz);
        Vec(x, i3, i4, dim, out double cx, out double cy, out double cz);
        // V = a·(b×c): ∂V/∂a = b×c, ∂V/∂b = c×a, ∂V/∂c = a×b, ∂V/∂p4 = −(sum)
        Cross(bx, by, bz, cx, cy, cz, out double da_x, out double da_y, out double da_z);
        Cross(cx, cy, cz, ax, ay, az, out double db_x, out double db_y, out double db_z);
        Cross(ax, ay, az, bx, by, bz, out double dc_x, out double dc_y, out double dc_z);
        int p1 = dim * i1, p2 = dim * i2, p3 = dim * i3, p4 = dim * i4;
        g[p1] += dEdV * da_x; g[p1 + 1] += dEdV * da_y; g[p1 + 2] += dEdV * da_z;
        g[p2] += dEdV * db_x; g[p2 + 1] += dEdV * db_y; g[p2 + 2] += dEdV * db_z;
        g[p3] += dEdV * dc_x; g[p3 + 1] += dEdV * dc_y; g[p3 + 2] += dEdV * dc_z;
        g[p4] -= dEdV * (da_x + db_x + dc_x); g[p4 + 1] -= dEdV * (da_y + db_y + dc_y); g[p4 + 2] -= dEdV * (da_z + db_z + dc_z);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Vec(ReadOnlySpan<double> x, int a, int b, int dim, out double vx, out double vy, out double vz)
    {
        int pa = dim * a, pb = dim * b;
        vx = x[pa] - x[pb]; vy = x[pa + 1] - x[pb + 1]; vz = x[pa + 2] - x[pb + 2];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Cross(double ax, double ay, double az, double bx, double by, double bz, out double cx, out double cy, out double cz)
    {
        cx = ay * bz - az * by;
        cy = az * bx - ax * bz;
        cz = ax * by - ay * bx;
    }
}
