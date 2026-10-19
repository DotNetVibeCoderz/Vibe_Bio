using System.ComponentModel;
using ILGPU;
using ILGPU.Algorithms;

namespace Biotech.Net.Gpu.Docking;

/// <summary>
/// Device-side description of a ligand and its grid maps (single precision). Coordinates live in a frame centred on the
/// box centre so that float rounding stays near 1e-6 Å. Infrastructure: public only because ILGPU's CPU accelerator emits
/// launchers into a dynamic assembly that cannot see internal types.
/// </summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public struct VinaKernelData
{
    /// <summary>Grid maps, map-major; cell (ix, iy, iz) at <c>(ix·Ny + iy)·Nz + iz</c>.</summary>
    public ArrayView<float> Maps;
    /// <summary>Reference coordinates relative to the root origin, xyz per atom.</summary>
    public ArrayView<float> Local;
    /// <summary>Scored (typed) atoms and the map index of each.</summary>
    public ArrayView<int> ScoredAtom, ScoredMap;
    /// <summary>Torsion axis atoms, deepest torsion first, and their moving atoms (CSR).</summary>
    public ArrayView<int> TorFrom, TorTo, MoveOffsets, MoveAtoms;
    /// <summary>Intramolecular pairs; flags bit 0 = both hydrophobic, bit 1 = H-bond donor/acceptor pair.</summary>
    public ArrayView<int> PairI, PairJ, PairFlags;
    /// <summary>Sum of the X-Score radii of each pair.</summary>
    public ArrayView<float> PairRadius;
    /// <summary>Counts.</summary>
    public int Atoms, Scored, Torsions, Pairs;
    /// <summary>Grid points per axis.</summary>
    public int Nx, Ny, Nz;
    /// <summary>Grid origin (frame coordinates) and spacing.</summary>
    public float Ox, Oy, Oz, Spacing;
    /// <summary>Box edge lengths (random starts are drawn inside the box).</summary>
    public float BoxX, BoxY, BoxZ;
    /// <summary>Vina weights and cut-off.</summary>
    public float W1, W2, WRep, WHyd, WHb, CutOff;
    /// <summary>Heavy-atom radius of gyration (orientation mutation amplitude).</summary>
    public float GyrationRadius;
}

/// <summary>Per-lane scratch, interleaved so that neighbouring lanes touch neighbouring words. Infrastructure.</summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public struct VinaScratch
{
    /// <summary>Float scratch: <see cref="VinaLayout.Size"/> values per lane.</summary>
    public ArrayView<float> Data;
    /// <summary>Energy evaluations per lane.</summary>
    public ArrayView<int> Evaluations;
    /// <summary>Lanes in the launch.</summary>
    public int Lanes;
}

/// <summary>Arguments of <see cref="VinaKernels.EvaluatePoses"/>. Infrastructure.</summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public struct VinaEvalArgs
{
    /// <summary>Poses (7 + T values each) in, energies and conformation gradients (6 + T each) out.</summary>
    public ArrayView<float> Poses, Energies, Gradients;
    /// <summary>Curl cap (1000 = none).</summary>
    public float Cap;
    /// <summary>Write gradients.</summary>
    public int WantGradients;
}

/// <summary>Arguments of <see cref="VinaKernels.MonteCarlo"/>. Infrastructure.</summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public struct VinaMcArgs
{
    /// <summary>Per-lane RNG state, best energy and best pose (7 + T values).</summary>
    public ArrayView<ulong> Rng;
    /// <summary>Best (refined, uncapped) energy per chain.</summary>
    public ArrayView<float> BestEnergy;
    /// <summary>Best pose per chain.</summary>
    public ArrayView<float> BestPose;
    /// <summary>Seed and the step range of this launch.</summary>
    public ulong Seed;
    /// <summary>First step of this launch (0 initialises the chains).</summary>
    public int FirstStep;
    /// <summary>Steps in this launch.</summary>
    public int Steps;
    /// <summary>BFGS iterations per local optimisation.</summary>
    public int LocalSteps;
}

/// <summary>Offsets inside a lane's scratch. Infrastructure.</summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public static class VinaLayout
{
    /// <summary>Conformation length: position, quaternion (w, x, y, z), torsions.</summary>
    public static int Conf(int torsions) => 7 + torsions;
    /// <summary>Gradient length: force, torque, torsions.</summary>
    public static int Grad(int torsions) => 6 + torsions;
    /// <summary>Atom coordinates.</summary>
    public static int X(int atoms) => 0;
    /// <summary>Atom gradients.</summary>
    public static int AtomGrad(int atoms) => 3 * atoms;
    /// <summary>Conformation slot k (0 current, 1 candidate, 2 trial, 3 original, 4 best).</summary>
    public static int Slot(int atoms, int torsions, int k) => 6 * atoms + k * Conf(torsions);
    /// <summary>Vector slot k (0 g, 1 g_new, 2 p, 3 y, 4 −Hy).</summary>
    public static int Vec(int atoms, int torsions, int k) => Slot(atoms, torsions, 5) + k * Grad(torsions);
    /// <summary>Inverse Hessian.</summary>
    public static int Hessian(int atoms, int torsions) => Vec(atoms, torsions, 5);
    /// <summary>Floats per lane.</summary>
    public static int Size(int atoms, int torsions) => Hessian(atoms, torsions) + Grad(torsions) * Grad(torsions);
}

/// <summary>
/// Vina scoring and search on the device, mirroring <c>VinaGrid.Evaluate</c>, <c>VinaDocking.Intra</c>,
/// <c>FlexibleLigand.Build/ConfGradient</c> and the CPU engine's BFGS and Monte Carlo in single precision. Infrastructure:
/// use <see cref="GpuVinaScorer"/> and <see cref="GpuVinaDocker"/>.
/// </summary>
[EditorBrowsable(EditorBrowsableState.Never)]
public static class VinaKernels
{
    private const float OutsideSlope = 1e6f;
    private const float Pi = 3.14159265358979f, TwoPi = 6.28318530717959f;

    /// <summary>One lane per pose: energy (inter + intra, capped) and optionally the conformation gradient.</summary>
    public static void EvaluatePoses(Index1D lane, VinaKernelData d, VinaScratch s, VinaEvalArgs a)
    {
        int nc = VinaLayout.Conf(d.Torsions), ng = VinaLayout.Grad(d.Torsions);
        int conf = VinaLayout.Slot(d.Atoms, d.Torsions, 0), g = VinaLayout.Vec(d.Atoms, d.Torsions, 0);
        for (int k = 0; k < nc; k++) Set(s, lane, conf + k, a.Poses[(long)lane.X * nc + k]);
        float e = Evaluate(d, s, lane, conf, a.Cap, g);
        a.Energies[lane] = e;
        if (a.WantGradients != 0)
            for (int k = 0; k < ng; k++) a.Gradients[(long)lane.X * ng + k] = Get(s, lane, g + k);
    }

    /// <summary>One lane per pose: BFGS local optimisation in place (poses in and out), final energy out.</summary>
    public static void OptimizePoses(Index1D lane, VinaKernelData d, VinaScratch s, VinaEvalArgs a, int maxSteps)
    {
        int nc = VinaLayout.Conf(d.Torsions);
        int conf = VinaLayout.Slot(d.Atoms, d.Torsions, 0);
        s.Evaluations[lane] = 0;
        for (int k = 0; k < nc; k++) Set(s, lane, conf + k, a.Poses[(long)lane.X * nc + k]);
        a.Energies[lane] = Bfgs(d, s, lane, conf, maxSteps, a.Cap);
        for (int k = 0; k < nc; k++) a.Poses[(long)lane.X * nc + k] = Get(s, lane, conf + k);
    }

    /// <summary>
    /// One lane per Monte-Carlo chain, <see cref="VinaMcArgs.Steps"/> steps per launch: mutate one degree of freedom,
    /// BFGS with the hunting cap, Metropolis at T = 1.2, and a full-precision-cap refinement of each new best.
    /// </summary>
    public static void MonteCarlo(Index1D lane, VinaKernelData d, VinaScratch s, VinaMcArgs a)
    {
        const float huntCap = 10f, authentic = 1000f, temperature = 1.2f, amplitude = 2f;
        int nc = VinaLayout.Conf(d.Torsions);
        int cur = VinaLayout.Slot(d.Atoms, d.Torsions, 0), cand = VinaLayout.Slot(d.Atoms, d.Torsions, 1);
        int best = VinaLayout.Slot(d.Atoms, d.Torsions, 4);
        // energies of the chain live after the Hessian
        int eCur = VinaLayout.Size(d.Atoms, d.Torsions), eBest = eCur + 1;
        ulong rng;
        if (a.FirstStep == 0)
        {
            rng = SplitMix(a.Seed + 0x9E3779B97F4A7C15UL * (ulong)(lane + 1));
            rng = Next(rng); Set(s, lane, cur + 0, d.Ox + Uniform(rng) * d.BoxX);
            rng = Next(rng); Set(s, lane, cur + 1, d.Oy + Uniform(rng) * d.BoxY);
            rng = Next(rng); Set(s, lane, cur + 2, d.Oz + Uniform(rng) * d.BoxZ);
            rng = Next(rng); float u1 = Uniform(rng);
            rng = Next(rng); float u2 = Uniform(rng) * TwoPi;
            rng = Next(rng); float u3 = Uniform(rng) * TwoPi;
            float qa = XMath.Sqrt(1 - u1), qb = XMath.Sqrt(u1);
            Set(s, lane, cur + 3, qb * XMath.Cos(u3));
            Set(s, lane, cur + 4, qa * XMath.Sin(u2));
            Set(s, lane, cur + 5, qa * XMath.Cos(u2));
            Set(s, lane, cur + 6, qb * XMath.Sin(u3));
            for (int t = 0; t < d.Torsions; t++) { rng = Next(rng); Set(s, lane, cur + 7 + t, (Uniform(rng) * 2 - 1) * Pi); }
            Set(s, lane, eCur, float.MaxValue);
            Set(s, lane, eBest, float.MaxValue);
            Copy(s, lane, cur, best, nc);
            s.Evaluations[lane] = 0;
        }
        else rng = a.Rng[lane];

        float currentE = Get(s, lane, eCur), bestE = Get(s, lane, eBest);
        int entities = 2 + d.Torsions;
        for (int step = a.FirstStep; step < a.FirstStep + a.Steps; step++)
        {
            Copy(s, lane, cur, cand, nc);
            rng = Next(rng);
            int which = (int)(Uniform(rng) * entities);
            if (which >= entities) which = entities - 1;
            if (which == 0)
            {
                float vx, vy, vz;
                do
                {
                    rng = Next(rng); vx = Uniform(rng) * 2 - 1;
                    rng = Next(rng); vy = Uniform(rng) * 2 - 1;
                    rng = Next(rng); vz = Uniform(rng) * 2 - 1;
                } while (vx * vx + vy * vy + vz * vz >= 1);
                Set(s, lane, cand, Get(s, lane, cand) + vx * amplitude);
                Set(s, lane, cand + 1, Get(s, lane, cand + 1) + vy * amplitude);
                Set(s, lane, cand + 2, Get(s, lane, cand + 2) + vz * amplitude);
            }
            else if (which == 1)
            {
                float vx, vy, vz;
                do
                {
                    rng = Next(rng); vx = Uniform(rng) * 2 - 1;
                    rng = Next(rng); vy = Uniform(rng) * 2 - 1;
                    rng = Next(rng); vz = Uniform(rng) * 2 - 1;
                } while (vx * vx + vy * vy + vz * vz >= 1);
                if (d.GyrationRadius > 1e-6f)
                {
                    float f = amplitude / d.GyrationRadius;
                    RotateOrientation(s, lane, cand, vx * f, vy * f, vz * f);
                }
            }
            else
            {
                rng = Next(rng);
                Set(s, lane, cand + 7 + (which - 2), (Uniform(rng) * 2 - 1) * Pi);
            }
            float e = Bfgs(d, s, lane, cand, a.LocalSteps, huntCap);
            bool accept = step == 0 || e < currentE;
            if (!accept)
            {
                rng = Next(rng);
                accept = Uniform(rng) < XMath.Exp((currentE - e) / temperature);
            }
            if (!accept) continue;
            Copy(s, lane, cand, cur, nc);
            currentE = e;
            if (e < bestE)
            {
                float re = Bfgs(d, s, lane, cur, a.LocalSteps, authentic);
                currentE = re;
                if (re < bestE)
                {
                    bestE = re;
                    Copy(s, lane, cur, best, nc);
                }
            }
        }
        Set(s, lane, eCur, currentE);
        Set(s, lane, eBest, bestE);
        a.Rng[lane] = rng;
        a.BestEnergy[lane] = bestE;
        for (int k = 0; k < nc; k++) a.BestPose[(long)lane.X * nc + k] = Get(s, lane, best + k);
    }

    // ---- energy -------------------------------------------------------------------------------------------------

    private static float Evaluate(VinaKernelData d, VinaScratch s, int lane, int conf, float cap, int gOut)
    {
        s.Evaluations[lane] += 1;
        int x = VinaLayout.X(d.Atoms), ag = VinaLayout.AtomGrad(d.Atoms);
        Build(d, s, lane, conf);
        for (int k = 0; k < 3 * d.Atoms; k++) Set(s, lane, ag + k, 0f);
        float e = 0;
        for (int k = 0; k < d.Scored; k++)
        {
            int i = d.ScoredAtom[k];
            var r = Grid(d, d.ScoredMap[k], Get(s, lane, x + 3 * i), Get(s, lane, x + 3 * i + 1), Get(s, lane, x + 3 * i + 2), cap);
            e += r.W;
            Add(s, lane, ag + 3 * i, r.X);
            Add(s, lane, ag + 3 * i + 1, r.Y);
            Add(s, lane, ag + 3 * i + 2, r.Z);
        }
        for (int p = 0; p < d.Pairs; p++)
        {
            int i = d.PairI[p], j = d.PairJ[p];
            float dx = Get(s, lane, x + 3 * i) - Get(s, lane, x + 3 * j);
            float dy = Get(s, lane, x + 3 * i + 1) - Get(s, lane, x + 3 * j + 1);
            float dz = Get(s, lane, x + 3 * i + 2) - Get(s, lane, x + 3 * j + 2);
            float r = XMath.Sqrt(dx * dx + dy * dy + dz * dz);
            if (r >= d.CutOff) continue;
            float dd = r - d.PairRadius[p];
            int flags = d.PairFlags[p];
            float g1 = XMath.Exp(-(dd / 0.5f) * (dd / 0.5f));
            float g2 = XMath.Exp(-((dd - 3f) / 2f) * ((dd - 3f) / 2f));
            float pe = d.W1 * g1 + d.W2 * g2;
            float dEdr = d.W1 * g1 * (-2 * dd / 0.25f) + d.W2 * g2 * (-2 * (dd - 3f) / 4f);
            if (dd < 0) { pe += d.WRep * dd * dd; dEdr += d.WRep * 2 * dd; }
            if ((flags & 1) != 0)
            {
                if (dd < 0.5f) pe += d.WHyd;
                else if (dd <= 1.5f) { pe += d.WHyd * (1.5f - dd); dEdr -= d.WHyd; }
            }
            if ((flags & 2) != 0)
            {
                if (dd < -0.7f) pe += d.WHb;
                else if (dd <= 0) { pe += d.WHb * (-dd / 0.7f); dEdr -= d.WHb / 0.7f; }
            }
            if (pe > 0 && cap < 1000)
            {
                float t = cap / (cap + pe);
                pe *= t;
                dEdr *= t * t;
            }
            e += pe;
            if (r > 1e-9f)
            {
                float f = dEdr / r;
                Add(s, lane, ag + 3 * i, dx * f); Add(s, lane, ag + 3 * i + 1, dy * f); Add(s, lane, ag + 3 * i + 2, dz * f);
                Add(s, lane, ag + 3 * j, -dx * f); Add(s, lane, ag + 3 * j + 1, -dy * f); Add(s, lane, ag + 3 * j + 2, -dz * f);
            }
        }
        ConfGradient(d, s, lane, conf, gOut);
        return e;
    }

    private static Float4 Grid(VinaKernelData d, int map, float px, float py, float pz, float cap)
    {
        float inv = 1f / d.Spacing;
        float sx = (px - d.Ox) * inv, sy = (py - d.Oy) * inv, sz = (pz - d.Oz) * inv;
        float penalty = 0;
        int ax, ay, az, rx, ry, rz;
        if (sx < 0) { penalty += -sx * d.Spacing; rx = -1; ax = 0; sx = 0; }
        else if (sx >= d.Nx - 1) { penalty += (sx - (d.Nx - 1)) * d.Spacing; rx = 1; ax = d.Nx - 2; sx = 1; }
        else { rx = 0; ax = (int)sx; sx -= ax; }
        if (sy < 0) { penalty += -sy * d.Spacing; ry = -1; ay = 0; sy = 0; }
        else if (sy >= d.Ny - 1) { penalty += (sy - (d.Ny - 1)) * d.Spacing; ry = 1; ay = d.Ny - 2; sy = 1; }
        else { ry = 0; ay = (int)sy; sy -= ay; }
        if (sz < 0) { penalty += -sz * d.Spacing; rz = -1; az = 0; sz = 0; }
        else if (sz >= d.Nz - 1) { penalty += (sz - (d.Nz - 1)) * d.Spacing; rz = 1; az = d.Nz - 2; sz = 1; }
        else { rz = 0; az = (int)sz; sz -= az; }
        penalty *= OutsideSlope;
        long cells = (long)d.Nx * d.Ny * d.Nz, b = map * cells;
        long i000 = b + ((long)ax * d.Ny + ay) * d.Nz + az;
        long dxs = (long)d.Ny * d.Nz, dys = d.Nz;
        float f000 = d.Maps[i000], f100 = d.Maps[i000 + dxs];
        float f010 = d.Maps[i000 + dys], f110 = d.Maps[i000 + dxs + dys];
        float f001 = d.Maps[i000 + 1], f101 = d.Maps[i000 + dxs + 1];
        float f011 = d.Maps[i000 + dys + 1], f111 = d.Maps[i000 + dxs + dys + 1];
        float x = sx, y = sy, z = sz, mx = 1 - x, my = 1 - y, mz = 1 - z;
        float f = f000 * mx * my * mz + f100 * x * my * mz + f010 * mx * y * mz + f110 * x * y * mz
                + f001 * mx * my * z + f101 * x * my * z + f011 * mx * y * z + f111 * x * y * z;
        float gx = (f100 - f000) * my * mz + (f110 - f010) * y * mz + (f101 - f001) * my * z + (f111 - f011) * y * z;
        float gy = (f010 - f000) * mx * mz + (f110 - f100) * x * mz + (f011 - f001) * mx * z + (f111 - f101) * x * z;
        float gz = (f001 - f000) * mx * my + (f101 - f100) * x * my + (f011 - f010) * mx * y + (f111 - f110) * x * y;
        float scale = 1;
        if (f > 0 && cap < 1000)
        {
            float t = cap / (cap + f);
            f *= t;
            scale = t * t;
        }
        gx *= scale * inv; gy *= scale * inv; gz *= scale * inv;
        return new Float4(rx == 0 ? gx : OutsideSlope * rx, ry == 0 ? gy : OutsideSlope * ry, rz == 0 ? gz : OutsideSlope * rz, f + penalty);
    }

    private static void Build(VinaKernelData d, VinaScratch s, int lane, int conf)
    {
        int x = VinaLayout.X(d.Atoms);
        for (int k = 0; k < 3 * d.Atoms; k++) Set(s, lane, x + k, d.Local[k]);
        for (int t = 0; t < d.Torsions; t++)
        {
            float angle = Get(s, lane, conf + 7 + t);
            if (angle == 0) continue;
            int from = d.TorFrom[t], to = d.TorTo[t];
            float ax = Get(s, lane, x + 3 * from), ay = Get(s, lane, x + 3 * from + 1), az = Get(s, lane, x + 3 * from + 2);
            float ux = Get(s, lane, x + 3 * to) - ax, uy = Get(s, lane, x + 3 * to + 1) - ay, uz = Get(s, lane, x + 3 * to + 2) - az;
            float len = XMath.Sqrt(ux * ux + uy * uy + uz * uz);
            ux /= len; uy /= len; uz /= len;
            float c = XMath.Cos(angle), sn = XMath.Sin(angle), kk = 1 - c;
            float r00 = c + ux * ux * kk, r01 = ux * uy * kk - uz * sn, r02 = ux * uz * kk + uy * sn;
            float r10 = uy * ux * kk + uz * sn, r11 = c + uy * uy * kk, r12 = uy * uz * kk - ux * sn;
            float r20 = uz * ux * kk - uy * sn, r21 = uz * uy * kk + ux * sn, r22 = c + uz * uz * kk;
            for (int m = d.MoveOffsets[t]; m < d.MoveOffsets[t + 1]; m++)
            {
                int i = d.MoveAtoms[m];
                float vx = Get(s, lane, x + 3 * i) - ax, vy = Get(s, lane, x + 3 * i + 1) - ay, vz = Get(s, lane, x + 3 * i + 2) - az;
                Set(s, lane, x + 3 * i, ax + r00 * vx + r01 * vy + r02 * vz);
                Set(s, lane, x + 3 * i + 1, ay + r10 * vx + r11 * vy + r12 * vz);
                Set(s, lane, x + 3 * i + 2, az + r20 * vx + r21 * vy + r22 * vz);
            }
        }
        float qw = Get(s, lane, conf + 3), qx = Get(s, lane, conf + 4), qy = Get(s, lane, conf + 5), qz = Get(s, lane, conf + 6);
        float n = XMath.Sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
        qw /= n; qx /= n; qy /= n; qz /= n;
        float m00 = 1 - 2 * (qy * qy + qz * qz), m01 = 2 * (qx * qy - qw * qz), m02 = 2 * (qx * qz + qw * qy);
        float m10 = 2 * (qx * qy + qw * qz), m11 = 1 - 2 * (qx * qx + qz * qz), m12 = 2 * (qy * qz - qw * qx);
        float m20 = 2 * (qx * qz - qw * qy), m21 = 2 * (qy * qz + qw * qx), m22 = 1 - 2 * (qx * qx + qy * qy);
        float px = Get(s, lane, conf), py = Get(s, lane, conf + 1), pz = Get(s, lane, conf + 2);
        for (int i = 0; i < d.Atoms; i++)
        {
            float vx = Get(s, lane, x + 3 * i), vy = Get(s, lane, x + 3 * i + 1), vz = Get(s, lane, x + 3 * i + 2);
            Set(s, lane, x + 3 * i, px + m00 * vx + m01 * vy + m02 * vz);
            Set(s, lane, x + 3 * i + 1, py + m10 * vx + m11 * vy + m12 * vz);
            Set(s, lane, x + 3 * i + 2, pz + m20 * vx + m21 * vy + m22 * vz);
        }
    }

    private static void ConfGradient(VinaKernelData d, VinaScratch s, int lane, int conf, int g)
    {
        int x = VinaLayout.X(d.Atoms), ag = VinaLayout.AtomGrad(d.Atoms);
        float px = Get(s, lane, conf), py = Get(s, lane, conf + 1), pz = Get(s, lane, conf + 2);
        float fx = 0, fy = 0, fz = 0, tx = 0, ty = 0, tz = 0;
        for (int i = 0; i < d.Atoms; i++)
        {
            float gx = Get(s, lane, ag + 3 * i), gy = Get(s, lane, ag + 3 * i + 1), gz = Get(s, lane, ag + 3 * i + 2);
            float rx = Get(s, lane, x + 3 * i) - px, ry = Get(s, lane, x + 3 * i + 1) - py, rz = Get(s, lane, x + 3 * i + 2) - pz;
            fx += gx; fy += gy; fz += gz;
            tx += ry * gz - rz * gy; ty += rz * gx - rx * gz; tz += rx * gy - ry * gx;
        }
        Set(s, lane, g, fx); Set(s, lane, g + 1, fy); Set(s, lane, g + 2, fz);
        Set(s, lane, g + 3, tx); Set(s, lane, g + 4, ty); Set(s, lane, g + 5, tz);
        for (int t = 0; t < d.Torsions; t++)
        {
            int from = d.TorFrom[t], to = d.TorTo[t];
            float ax = Get(s, lane, x + 3 * from), ay = Get(s, lane, x + 3 * from + 1), az = Get(s, lane, x + 3 * from + 2);
            float ux = Get(s, lane, x + 3 * to) - ax, uy = Get(s, lane, x + 3 * to + 1) - ay, uz = Get(s, lane, x + 3 * to + 2) - az;
            float len = XMath.Sqrt(ux * ux + uy * uy + uz * uz);
            ux /= len; uy /= len; uz /= len;
            float qx = 0, qy = 0, qz = 0;
            for (int m = d.MoveOffsets[t]; m < d.MoveOffsets[t + 1]; m++)
            {
                int i = d.MoveAtoms[m];
                float rx = Get(s, lane, x + 3 * i) - ax, ry = Get(s, lane, x + 3 * i + 1) - ay, rz = Get(s, lane, x + 3 * i + 2) - az;
                float gx = Get(s, lane, ag + 3 * i), gy = Get(s, lane, ag + 3 * i + 1), gz = Get(s, lane, ag + 3 * i + 2);
                qx += ry * gz - rz * gy; qy += rz * gx - rx * gz; qz += rx * gy - ry * gx;
            }
            Set(s, lane, g + 6 + t, qx * ux + qy * uy + qz * uz);
        }
    }

    // ---- optimisation ---------------------------------------------------------------------------------------------

    /// <summary>The CPU engine's BFGS (dense inverse Hessian, Armijo backtracking), in place on the conformation slot.</summary>
    private static float Bfgs(VinaKernelData d, VinaScratch s, int lane, int conf, int maxSteps, float cap)
    {
        int n = VinaLayout.Grad(d.Torsions), nc = VinaLayout.Conf(d.Torsions);
        int trial = VinaLayout.Slot(d.Atoms, d.Torsions, 2), orig = VinaLayout.Slot(d.Atoms, d.Torsions, 3);
        int g = VinaLayout.Vec(d.Atoms, d.Torsions, 0), gNew = VinaLayout.Vec(d.Atoms, d.Torsions, 1);
        int p = VinaLayout.Vec(d.Atoms, d.Torsions, 2), y = VinaLayout.Vec(d.Atoms, d.Torsions, 3);
        int mhy = VinaLayout.Vec(d.Atoms, d.Torsions, 4), h = VinaLayout.Hessian(d.Atoms, d.Torsions);
        // float literals matter: ILGPU 1.5's OpenCL backend reinterprets the bits of an int ternary (1 -> 1.4e-45)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) Set(s, lane, h + i * n + j, i == j ? 1f : 0f);
        float f0 = Evaluate(d, s, lane, conf, cap, g);
        float fOrig = f0;
        Copy(s, lane, conf, orig, nc);
        for (int step = 0; step < maxSteps; step++)
        {
            for (int i = 0; i < n; i++)
            {
                float acc = 0;
                for (int j = 0; j < n; j++) acc += Get(s, lane, h + i * n + j) * Get(s, lane, g + j);
                Set(s, lane, p + i, -acc);
            }
            float pg = 0;
            for (int i = 0; i < n; i++) pg += Get(s, lane, p + i) * Get(s, lane, g + i);
            float alpha = 1, f1 = 0;
            for (int t = 0; t < 10; t++)
            {
                Copy(s, lane, conf, trial, nc);
                Increment(d, s, lane, trial, p, alpha);
                f1 = Evaluate(d, s, lane, trial, cap, gNew);
                if (f1 - f0 < 1e-4f * alpha * pg) break;
                alpha *= 0.5f;
            }
            for (int i = 0; i < n; i++) Set(s, lane, y + i, Get(s, lane, gNew + i) - Get(s, lane, g + i));
            f0 = f1;
            Copy(s, lane, trial, conf, nc);
            float gnorm = 0;
            for (int i = 0; i < n; i++) gnorm += Get(s, lane, g + i) * Get(s, lane, g + i);
            if (!(XMath.Sqrt(gnorm) >= 1e-5f)) break;
            Copy(s, lane, gNew, g, n);
            if (step == 0)
            {
                float yy = 0, yp = 0;
                for (int i = 0; i < n; i++) { float yi = Get(s, lane, y + i); yy += yi * yi; yp += yi * Get(s, lane, p + i); }
                if (XMath.Abs(yy) > 1e-12f)
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++) Set(s, lane, h + i * n + j, i == j ? alpha * yp / yy : 0f);
            }
            float ypd = 0;
            for (int i = 0; i < n; i++) ypd += Get(s, lane, y + i) * Get(s, lane, p + i);
            if (alpha * ypd < 1e-12f) continue;
            for (int i = 0; i < n; i++)
            {
                float acc = 0;
                for (int j = 0; j < n; j++) acc += Get(s, lane, h + i * n + j) * Get(s, lane, y + j);
                Set(s, lane, mhy + i, -acc);
            }
            float yhy = 0;
            for (int i = 0; i < n; i++) yhy -= Get(s, lane, y + i) * Get(s, lane, mhy + i);
            float r = 1 / (alpha * ypd);
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    float pi = Get(s, lane, p + i), pj = Get(s, lane, p + j);
                    float v = alpha * r * (Get(s, lane, mhy + i) * pj + Get(s, lane, mhy + j) * pi) + alpha * alpha * (r * r * yhy + r) * pi * pj;
                    Add(s, lane, h + i * n + j, v);
                    if (j != i) Add(s, lane, h + j * n + i, v);
                }
        }
        if (!(f0 <= fOrig)) { Copy(s, lane, orig, conf, nc); f0 = fOrig; }
        return f0;
    }

    /// <summary><c>LigandConf.Increment</c>: position += p·α, orientation ← exp(p_rot·α)·orientation, torsions += p·α (wrapped).</summary>
    private static void Increment(VinaKernelData d, VinaScratch s, int lane, int conf, int p, float alpha)
    {
        for (int k = 0; k < 3; k++) Add(s, lane, conf + k, Get(s, lane, p + k) * alpha);
        RotateOrientation(s, lane, conf, Get(s, lane, p + 3) * alpha, Get(s, lane, p + 4) * alpha, Get(s, lane, p + 5) * alpha);
        for (int t = 0; t < d.Torsions; t++)
            Set(s, lane, conf + 7 + t, WrapAngle(Get(s, lane, conf + 7 + t) + alpha * Get(s, lane, p + 6 + t)));
    }

    private static void RotateOrientation(VinaScratch s, int lane, int conf, float vx, float vy, float vz)
    {
        float angle = XMath.Sqrt(vx * vx + vy * vy + vz * vz);
        float aw = 1, ax = 0, ay = 0, az = 0;
        if (angle >= 1e-12f)
        {
            float sn = XMath.Sin(angle / 2) / angle;
            aw = XMath.Cos(angle / 2); ax = vx * sn; ay = vy * sn; az = vz * sn;
        }
        float bw = Get(s, lane, conf + 3), bx = Get(s, lane, conf + 4), by = Get(s, lane, conf + 5), bz = Get(s, lane, conf + 6);
        float w = aw * bw - ax * bx - ay * by - az * bz;
        float x = aw * bx + ax * bw + ay * bz - az * by;
        float y = aw * by - ax * bz + ay * bw + az * bx;
        float z = aw * bz + ax * by - ay * bx + az * bw;
        float n = XMath.Sqrt(w * w + x * x + y * y + z * z);
        Set(s, lane, conf + 3, w / n); Set(s, lane, conf + 4, x / n); Set(s, lane, conf + 5, y / n); Set(s, lane, conf + 6, z / n);
    }

    private static float WrapAngle(float a)
    {
        if (!(a > -1e30f && a < 1e30f)) return 0;   // NaN or infinite
        return a - TwoPi * XMath.Floor(a / TwoPi + 0.5f);
    }

    // ---- scratch and RNG helpers -----------------------------------------------------------------------------------

    private static float Get(VinaScratch s, int lane, int k) => s.Data[(long)k * s.Lanes + lane];

    private static void Set(VinaScratch s, int lane, int k, float v) => s.Data[(long)k * s.Lanes + lane] = v;

    private static void Add(VinaScratch s, int lane, int k, float v) => s.Data[(long)k * s.Lanes + lane] += v;

    private static void Copy(VinaScratch s, int lane, int from, int to, int count)
    {
        for (int k = 0; k < count; k++) Set(s, lane, to + k, Get(s, lane, from + k));
    }

    private static ulong SplitMix(ulong z)
    {
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        z ^= z >> 31;
        return z == 0 ? 0x2545F4914F6CDD1DUL : z;
    }

    /// <summary>xorshift64 step.</summary>
    private static ulong Next(ulong x)
    {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return x;
    }

    /// <summary>Uniform [0, 1) from the top 24 bits of the scrambled state.</summary>
    private static float Uniform(ulong x) => ((x * 0x2545F4914F6CDD1DUL) >> 40) * (1f / 16777216f);

    private readonly struct Float4(float x, float y, float z, float w)
    {
        public readonly float X = x, Y = y, Z = z, W = w;
    }
}
