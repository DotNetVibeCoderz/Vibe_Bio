namespace Biotech.Net.Chemistry;

/// <summary>
/// 2-D depiction coordinates for drawing molecules. Ring systems are built from regular polygons (fused rings
/// share an edge, spiro rings a vertex, bridges are placed on arcs) and laid out as rigid units; chains grow with
/// 120° zig-zags (trans by default), substituents go into the largest free angle, sp centres are linear. Clashes
/// are then removed by flipping or rotating the smaller side of acyclic single bonds, and disconnected fragments
/// are placed side by side. Each fragment is rotated so its long axis is horizontal.
/// </summary>
public static class Depiction
{
    /// <summary>Default bond length (RDKit uses 1.5 for depictions).</summary>
    public const double DefaultBondLength = 1.5;

    /// <summary>Computes 2-D coordinates (Z = 0), stores them in <see cref="Molecule.Coordinates"/> and returns them.</summary>
    public static Point3[] Compute2DCoordinates(Molecule mol, double bondLength = DefaultBondLength)
    {
        var coords = new Layout(mol, bondLength).Run();
        mol.Coordinates = coords;
        return coords;
    }

    /// <summary>
    /// Number of non-bonded atom pairs closer than <paramref name="factor"/> × the mean bond length
    /// (a quick quality check for depictions).
    /// </summary>
    public static int CountClashes(Molecule mol, double factor = 0.5)
    {
        var c = mol.Coordinates ?? throw new InvalidOperationException("The molecule has no coordinates.");
        double mean = mol.BondCount == 0 ? DefaultBondLength : Enumerable.Range(0, mol.BondCount).Average(b => Point3.Distance(c[mol.Bonds[b].Begin], c[mol.Bonds[b].End]));
        double lim = factor * mean;
        int n = 0;
        for (int i = 0; i < mol.AtomCount; i++)
            for (int j = i + 1; j < mol.AtomCount; j++)
                if (mol.BondBetween(i, j) < 0 && Point3.Distance(c[i], c[j]) < lim) n++;
        return n;
    }

    /// <summary>Simple 2-D vector used by the layout.</summary>
    private readonly record struct V(double X, double Y)
    {
        public static V operator +(V a, V b) => new(a.X + b.X, a.Y + b.Y);
        public static V operator -(V a, V b) => new(a.X - b.X, a.Y - b.Y);
        public static V operator *(V a, double s) => new(a.X * s, a.Y * s);
        public double Length => Math.Sqrt(X * X + Y * Y);
        public double Angle => Math.Atan2(Y, X);
        public V Normalized() { double l = Length; return l < 1e-12 ? new V(1, 0) : new V(X / l, Y / l); }
        public V Rotate(double a) { double c = Math.Cos(a), s = Math.Sin(a); return new V(X * c - Y * s, X * s + Y * c); }
        public static V FromAngle(double a) => new(Math.Cos(a), Math.Sin(a));
        public static double Cross(V a, V b) => a.X * b.Y - a.Y * b.X;
    }

    private sealed class Layout(Molecule mol, double L)
    {
        private readonly int _n = mol.AtomCount;
        private V[] _pos = [];
        private bool[] _placed = [];
        private int[] _atomSys = [];
        private List<int[]> _systems = [];

        public Point3[] Run()
        {
            _pos = new V[_n];
            _placed = new bool[_n];
            var rings = mol.Rings;
            _atomSys = Enumerable.Repeat(-1, _n).ToArray();
            _systems = [];
            foreach (var sys in rings.RingSystems)
            {
                var atoms = sys.SelectMany(r => rings.AtomRings[r]).Distinct().ToArray();
                foreach (int a in atoms) _atomSys[a] = _systems.Count;
                _systems.Add(sys);
            }

            double offsetX = 0;
            foreach (var frag in mol.Fragments())
            {
                LayoutFragment(frag);
                Resolve(frag);
                OrientHorizontally(frag);
                double minX = frag.Min(a => _pos[a].X), maxX = frag.Max(a => _pos[a].X);
                double minY = frag.Min(a => _pos[a].Y), maxY = frag.Max(a => _pos[a].Y);
                double dx = offsetX - minX, dy = -(minY + maxY) / 2;
                foreach (int a in frag) _pos[a] = _pos[a] + new V(dx, dy);
                offsetX += (maxX - minX) + 2 * L;
            }
            return _pos.Select(p => new Point3(Math.Round(p.X, 4), Math.Round(p.Y, 4), 0)).ToArray();
        }

        // ---- fragment layout ------------------------------------------------------------------

        private void LayoutFragment(int[] frag)
        {
            int start;
            var sysInFrag = frag.Where(a => _atomSys[a] >= 0).Select(a => _atomSys[a]).Distinct().ToList();
            var queue = new Queue<int>();
            if (sysInFrag.Count > 0)
            {
                int best = sysInFrag.OrderByDescending(s => SystemAtoms(s).Length).First();
                var local = LayoutRingSystem(best);
                foreach (var (a, p) in local) { _pos[a] = p; _placed[a] = true; }
                foreach (var a in local.Keys.OrderBy(x => x)) queue.Enqueue(a);
            }
            else
            {
                start = Peripheral(frag);
                _pos[start] = new V(0, 0);
                _placed[start] = true;
                queue.Enqueue(start);
            }
            while (queue.Count > 0)
            {
                int u = queue.Dequeue();
                var news = new List<int>();
                foreach (int v in mol.Neighbors(u)) if (!_placed[v]) news.Add(v);
                if (news.Count == 0) continue;
                // heavier branches first so they get the "best" directions
                news = news.OrderByDescending(v => BranchSize(u, v)).ToList();
                var angles = ChooseAngles(u, news.Count);
                for (int k = 0; k < news.Count; k++)
                {
                    int v = news[k];
                    if (_placed[v]) continue; // placed as part of a ring system in this loop
                    var pv = _pos[u] + V.FromAngle(angles[k]) * L;
                    if (_atomSys[v] >= 0)
                    {
                        PlaceSystemAt(_atomSys[v], v, pv, _pos[u]);
                        foreach (int a in SystemAtoms(_atomSys[v])) queue.Enqueue(a);
                    }
                    else
                    {
                        _pos[v] = pv;
                        _placed[v] = true;
                        queue.Enqueue(v);
                    }
                }
            }
        }

        private int Peripheral(int[] frag)
        {
            int Far(int s)
            {
                var dist = new Dictionary<int, int> { [s] = 0 };
                var q = new Queue<int>();
                q.Enqueue(s);
                int last = s;
                while (q.Count > 0)
                {
                    int u = q.Dequeue();
                    last = u;
                    foreach (int v in mol.Neighbors(u)) if (dist.TryAdd(v, dist[u] + 1)) q.Enqueue(v);
                }
                return last;
            }
            return Far(Far(frag[0]));
        }

        private int BranchSize(int from, int to)
        {
            var seen = new HashSet<int> { from, to };
            var stack = new Stack<int>();
            stack.Push(to);
            while (stack.Count > 0)
            {
                int u = stack.Pop();
                foreach (int v in mol.Neighbors(u)) if (seen.Add(v)) stack.Push(v);
                if (seen.Count > 64) break;
            }
            return seen.Count;
        }

        private double[] ChooseAngles(int u, int k)
        {
            var placedNbrs = new List<int>();
            foreach (int v in mol.Neighbors(u)) if (_placed[v]) placedNbrs.Add(v);
            var result = new double[k];
            if (placedNbrs.Count == 0)
            {
                if (k == 1) result[0] = 0;
                else if (k == 2) { result[0] = Math.PI / 6; result[1] = 5 * Math.PI / 6; }
                else for (int i = 0; i < k; i++) result[i] = Math.PI / 2 + 2 * Math.PI * i / k;
                return result;
            }
            if (placedNbrs.Count == 1)
            {
                int p = placedNbrs[0];
                double tp = (_pos[p] - _pos[u]).Angle;
                if (IsLinearCentre(u) && k == 1) { result[0] = tp + Math.PI; return result; }
                if (k == 1)
                {
                    // zig-zag: put the new atom trans to p's other placed neighbour
                    double sign = 1;
                    int q = -1;
                    foreach (int w in mol.Neighbors(p)) if (w != u && _placed[w]) { q = w; break; }
                    if (q >= 0)
                    {
                        var axis = _pos[u] - _pos[p];
                        double sideQ = V.Cross(axis, _pos[q] - _pos[p]);
                        var cand = _pos[u] + V.FromAngle(tp + 2 * Math.PI / 3) * L;
                        double sideC = V.Cross(axis, cand - _pos[p]);
                        sign = sideQ * sideC > 0 ? -1 : 1;
                    }
                    result[0] = tp + sign * 2 * Math.PI / 3;
                    return result;
                }
                for (int i = 0; i < k; i++) result[i] = tp + 2 * Math.PI * (i + 1) / (k + 1);
                return result;
            }
            // several placed neighbours: fill the largest angular gap
            var angs = placedNbrs.Select(v => Norm((_pos[v] - _pos[u]).Angle)).OrderBy(a => a).ToList();
            double bestGap = -1, gapStart = 0;
            for (int i = 0; i < angs.Count; i++)
            {
                double a0 = angs[i], a1 = i + 1 < angs.Count ? angs[i + 1] : angs[0] + 2 * Math.PI;
                if (a1 - a0 > bestGap) { bestGap = a1 - a0; gapStart = a0; }
            }
            for (int i = 0; i < k; i++) result[i] = gapStart + bestGap * (i + 1) / (k + 1);
            return result;
        }

        private static double Norm(double a)
        {
            while (a < 0) a += 2 * Math.PI;
            while (a >= 2 * Math.PI) a -= 2 * Math.PI;
            return a;
        }

        private bool IsLinearCentre(int u)
        {
            int doubles = 0;
            foreach (int b in mol.NeighborBonds(u))
            {
                var o = mol.Bonds[b].Order;
                if (o == BondOrder.Triple) return true;
                if (o == BondOrder.Double) doubles++;
            }
            return doubles >= 2;
        }

        // ---- ring systems -------------------------------------------------------------------------

        private int[] SystemAtoms(int sys) => Enumerable.Range(0, _n).Where(a => _atomSys[a] == sys).ToArray();

        private void PlaceSystemAt(int sys, int anchor, V anchorPos, V fromPos)
        {
            var local = LayoutRingSystem(sys);
            var centroid = new V(local.Values.Average(p => p.X), local.Values.Average(p => p.Y));
            var e = centroid - local[anchor];
            var d = anchorPos - fromPos;
            double rot = (e.Length < 1e-6 ? 0 : d.Angle - e.Angle);
            Dictionary<int, V>? best = null;
            double bestScore = double.MaxValue;
            foreach (bool mirror in new[] { false, true })
            {
                var cand = new Dictionary<int, V>();
                foreach (var (a, p) in local)
                {
                    var rel = p - local[anchor];
                    if (mirror)
                    {
                        // mirror across the anchor→centroid axis
                        double ea = e.Length < 1e-6 ? 0 : e.Angle;
                        rel = rel.Rotate(-ea);
                        rel = new V(rel.X, -rel.Y);
                        rel = rel.Rotate(ea);
                    }
                    cand[a] = anchorPos + rel.Rotate(rot);
                }
                double score = 0;
                for (int i = 0; i < _n; i++)
                {
                    if (!_placed[i]) continue;
                    foreach (var p in cand.Values)
                    {
                        double dd = (p - _pos[i]).Length;
                        if (dd < 1.2 * L) score += (1.2 * L - dd) * (1.2 * L - dd);
                    }
                }
                if (score < bestScore - 1e-9) { bestScore = score; best = cand; }
            }
            foreach (var (a, p) in best!) { _pos[a] = p; _placed[a] = true; }
        }

        private Dictionary<int, V> LayoutRingSystem(int sys)
        {
            var rings = mol.Rings.AtomRings;
            var ringIdx = _systems[sys].ToList();
            var pos = new Dictionary<int, V>();
            var bondShare = ringIdx.ToDictionary(r => r, r => ringIdx.Count(s => s != r && rings[s].Intersect(rings[r]).Count() >= 2));
            int first = ringIdx.OrderByDescending(r => bondShare[r]).ThenByDescending(r => rings[r].Length).First();
            var ring0 = rings[first];
            double R0 = L / (2 * Math.Sin(Math.PI / ring0.Length));
            for (int i = 0; i < ring0.Length; i++)
                pos[ring0[i]] = V.FromAngle(Math.PI / 2 + 2 * Math.PI * i / ring0.Length + (ring0.Length % 2 == 0 ? Math.PI / ring0.Length : 0)) * R0;
            var done = new HashSet<int> { first };
            while (done.Count < ringIdx.Count)
            {
                int next = ringIdx.Where(r => !done.Contains(r))
                    .OrderByDescending(r => rings[r].Count(pos.ContainsKey))
                    .ThenBy(r => rings[r].Length).First();
                var ring = rings[next];
                done.Add(next);
                int placedCount = ring.Count(pos.ContainsKey);
                if (placedCount == ring.Length) continue;
                if (placedCount == 0)
                {
                    // disconnected from what is placed so far (should not happen within a system); put it to the right
                    double maxX = pos.Values.Max(p => p.X);
                    double R = L / (2 * Math.Sin(Math.PI / ring.Length));
                    for (int i = 0; i < ring.Length; i++) pos[ring[i]] = new V(maxX + 2 * L + R, 0) + V.FromAngle(2 * Math.PI * i / ring.Length) * R;
                    continue;
                }
                var centroid = new V(pos.Values.Average(p => p.X), pos.Values.Average(p => p.Y));
                if (placedCount == 1)
                {
                    int s = ring.First(pos.ContainsKey);
                    var nbrsPlaced = mol.Neighbors(s).ToArray().Where(pos.ContainsKey).ToList();
                    var away = nbrsPlaced.Count == 0 ? (pos[s] - centroid) : pos[s] - new V(nbrsPlaced.Average(a => pos[a].X), nbrsPlaced.Average(a => pos[a].Y));
                    var dir = away.Normalized();
                    double R = L / (2 * Math.Sin(Math.PI / ring.Length));
                    var c = pos[s] + dir * R;
                    int si = Array.IndexOf(ring, s);
                    double a0 = (pos[s] - c).Angle;
                    for (int i = 1; i < ring.Length; i++)
                        pos[ring[(si + i) % ring.Length]] = c + V.FromAngle(a0 + 2 * Math.PI * i / ring.Length) * R;
                    continue;
                }
                // place each maximal unplaced run on an arc between its placed end atoms
                int nR = ring.Length;
                int startIdx = Enumerable.Range(0, nR).First(i => pos.ContainsKey(ring[i]) && !pos.ContainsKey(ring[(i + 1) % nR]));
                int iIdx = startIdx;
                int guard = 0;
                while (guard++ < nR)
                {
                    // ring[iIdx] placed, ring[iIdx+1] unplaced
                    var run = new List<int>();
                    int j = (iIdx + 1) % nR;
                    while (!pos.ContainsKey(ring[j])) { run.Add(ring[j]); j = (j + 1) % nR; }
                    var arc = Arc(pos[ring[iIdx]], pos[ring[j]], run.Count, centroid);
                    for (int k = 0; k < run.Count; k++) pos[run[k]] = arc[k];
                    // next unplaced run
                    int found = -1;
                    for (int t = 0; t < nR; t++)
                    {
                        int idx = (j + t) % nR;
                        if (pos.ContainsKey(ring[idx]) && !pos.ContainsKey(ring[(idx + 1) % nR])) { found = idx; break; }
                    }
                    if (found < 0) break;
                    iIdx = found;
                }
            }
            return pos;
        }

        /// <summary>k points on an arc from a to b with k+1 equal segments of length L, bulging away from <paramref name="avoid"/>.</summary>
        private List<V> Arc(V a, V b, int k, V avoid)
        {
            var result = new List<V>(k);
            double c = (b - a).Length;
            int segs = k + 1;
            var mid = (a + b) * 0.5;
            var chordDir = (b - a).Normalized();
            var normal = new V(-chordDir.Y, chordDir.X);
            if (Dot(normal, mid - avoid) < 0) normal = normal * -1;
            if (c < 1e-6)
            {
                // closed loop (spiro-like): regular polygon through a
                double R = L / (2 * Math.Sin(Math.PI / segs));
                var centre = a + normal * R;
                double a0 = (a - centre).Angle;
                for (int i = 1; i <= k; i++) result.Add(centre + V.FromAngle(a0 + 2 * Math.PI * i / segs) * R);
                return result;
            }
            if (c >= segs * L * 0.999)
            {
                for (int i = 1; i <= k; i++) result.Add(a + (b - a) * ((double)i / segs));
                return result;
            }
            // solve sin(segs*phi/2)/sin(phi/2) = c/L for phi in (0, 2π/segs)
            double target = c / L;
            double lo = 1e-9, hi = 2 * Math.PI / segs - 1e-9;
            double F(double phi) => Math.Sin(segs * phi / 2) / Math.Sin(phi / 2) - target;
            // F decreases from segs to 0 on (0, 2π/segs)
            for (int it = 0; it < 100; it++)
            {
                double m = (lo + hi) / 2;
                if (F(m) > 0) lo = m; else hi = m;
            }
            double phiS = (lo + hi) / 2;
            double Rr = L / (2 * Math.Sin(phiS / 2));
            double total = segs * phiS; // angle subtended by the arc
            // centre on the opposite side of the bulge if total < π, else on the same side
            double h = Math.Sqrt(Math.Max(0, Rr * Rr - c * c / 4));
            var centre2 = total <= Math.PI ? mid - normal * h : mid + normal * h;
            double angA = (a - centre2).Angle, angB = (b - centre2).Angle;
            // go from a to b through the bulge side
            double step = phiS;
            var test = centre2 + V.FromAngle(angA + step) * Rr;
            var test2 = centre2 + V.FromAngle(angA - step) * Rr;
            if (Dot(test - mid, normal) < Dot(test2 - mid, normal)) step = -step;
            _ = angB;
            for (int i = 1; i <= k; i++) result.Add(centre2 + V.FromAngle(angA + step * i) * Rr);
            return result;
        }

        private static double Dot(V a, V b) => a.X * b.X + a.Y * b.Y;

        // ---- clean-up -------------------------------------------------------------------------------

        private double ClashScore(int[] frag)
        {
            double thr = 0.8 * L, s = 0;
            for (int x = 0; x < frag.Length; x++)
                for (int y = x + 1; y < frag.Length; y++)
                {
                    int i = frag[x], j = frag[y];
                    double d = (_pos[i] - _pos[j]).Length;
                    if (d >= thr) continue;
                    if (mol.BondBetween(i, j) >= 0) continue;
                    s += (thr - d) * (thr - d);
                }
            return s;
        }

        private void Resolve(int[] frag)
        {
            if (frag.Length < 4) return;
            var rings = mol.Rings;
            var rotors = new List<int>();
            for (int b = 0; b < mol.BondCount; b++)
            {
                var bond = mol.Bonds[b];
                if (rings.IsBondInRing(b) || !frag.Contains(bond.Begin)) continue;
                if (mol.Degree(bond.Begin) < 2 || mol.Degree(bond.End) < 2) continue;
                rotors.Add(b);
            }
            double score = ClashScore(frag);
            for (int round = 0; round < 6 && score > 1e-9; round++)
            {
                bool improved = false;
                foreach (int b in rotors)
                {
                    var bond = mol.Bonds[b];
                    var side = SideAtoms(bond.Begin, bond.End);
                    int pivot = bond.Begin;
                    if (side.Count > frag.Length / 2) { side = SideAtoms(bond.End, bond.Begin); pivot = bond.End; }
                    int other = bond.Other(pivot);
                    var saved = side.ToDictionary(a => a, a => _pos[a]);
                    double bestScore = score;
                    Dictionary<int, V>? best = null;
                    // 1. mirror the side across the bond axis
                    var axisP = _pos[pivot];
                    var axisD = (_pos[other] - _pos[pivot]).Normalized();
                    foreach (int a in side)
                    {
                        var rel = saved[a] - axisP;
                        double along = Dot(rel, axisD);
                        var proj = axisD * along;
                        _pos[a] = axisP + proj * 2 - rel;
                    }
                    double sMirror = ClashScore(frag);
                    if (sMirror < bestScore - 1e-9) { bestScore = sMirror; best = side.ToDictionary(a => a, a => _pos[a]); }
                    // 2. rotate the side about the pivot atom
                    foreach (double deg in new[] { 30.0, -30.0, 60.0, -60.0, 90.0, -90.0, 120.0, -120.0 })
                    {
                        double rad = deg * Math.PI / 180;
                        foreach (int a in side) _pos[a] = axisP + (saved[a] - axisP).Rotate(rad);
                        double sr = ClashScore(frag) + 0.02 * Math.Abs(deg) / 30; // mild preference for small moves
                        if (sr < bestScore - 1e-9) { bestScore = sr; best = side.ToDictionary(a => a, a => _pos[a]); }
                    }
                    if (best is not null)
                    {
                        foreach (var (a, p) in best) _pos[a] = p;
                        score = ClashScore(frag);
                        improved = true;
                    }
                    else foreach (var (a, p) in saved) _pos[a] = p;
                    if (score <= 1e-9) break;
                }
                if (!improved) break;
            }
        }

        /// <summary>Atoms reachable from <paramref name="start"/> without passing through <paramref name="block"/> (start included).</summary>
        private List<int> SideAtoms(int block, int start)
        {
            var seen = new HashSet<int> { block, start };
            var list = new List<int> { start };
            var stack = new Stack<int>();
            stack.Push(start);
            while (stack.Count > 0)
            {
                int u = stack.Pop();
                foreach (int v in mol.Neighbors(u))
                    if (seen.Add(v)) { list.Add(v); stack.Push(v); }
            }
            return list;
        }

        private void OrientHorizontally(int[] frag)
        {
            if (frag.Length < 3) return;
            double mx = frag.Average(a => _pos[a].X), my = frag.Average(a => _pos[a].Y);
            double sxx = 0, syy = 0, sxy = 0;
            foreach (int a in frag)
            {
                double dx = _pos[a].X - mx, dy = _pos[a].Y - my;
                sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
            }
            double angle = 0.5 * Math.Atan2(2 * sxy, sxx - syy);
            var c = new V(mx, my);
            foreach (int a in frag) _pos[a] = c + (_pos[a] - c).Rotate(-angle);
        }
    }
}
