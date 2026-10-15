namespace Biotech.Net.StructuralBiology.Tests;

public class SuperpositionTests
{
    private static Vec3[] CrambinCa() => Superposition.CaTrace(SampleStructures.Crambin()["A"]);

    [Fact]
    public void RotatedTranslatedCopyIsRecoveredExactly()
    {
        var target = CrambinCa();
        var rotation = Matrix3.FromAxisAngle(new Vec3(1, 2, 3), 37);
        var shift = new Vec3(12.5, -4.0, 7.25);
        var mobile = target.Select(p => rotation * p + shift).ToArray();
        Assert.True(Superposition.Rmsd(mobile, target) > 5);

        var fit = Superposition.Kabsch(mobile, target);
        Assert.Equal(0, fit.Rmsd, 6);
        Assert.Equal(1, fit.Rotation.Determinant, 9);
        Assert.Equal(37, fit.Rotation.RotationAngle, 6);
        // The fitted rotation is the inverse (transpose) of the applied one.
        var product = fit.Rotation * rotation;
        Assert.Equal(0, product.RotationAngle, 5);

        var moved = mobile.ToArray();
        fit.Apply(moved);
        Assert.Equal(0, Superposition.Rmsd(moved, target), 6);
        for (int i = 0; i < moved.Length; i++) Assert.True(Vec3.Distance(moved[i], target[i]) < 1e-6);
    }

    [Fact]
    public void ScaledSquareHasAnalyticRmsd()
    {
        // Corners (±1, ±1, 0) vs the same square scaled by 1.1: by symmetry the best rotation is the identity and every
        // point is off by 0.1·√2, so RMSD = 0.1·√2.
        Vec3[] square = [new(1, 1, 0), new(-1, 1, 0), new(-1, -1, 0), new(1, -1, 0)];
        var scaled = square.Select(p => p * 1.1 + new Vec3(3, 3, 3)).ToArray();
        var fit = Superposition.Kabsch(scaled, square);
        Assert.Equal(0.1 * Math.Sqrt(2), fit.Rmsd, 9);
        Assert.Equal(0, fit.Rotation.RotationAngle, 4);
        var moved = scaled.ToArray();
        fit.Apply(moved);
        Assert.Equal(fit.Rmsd, Superposition.Rmsd(moved, square), 9);
    }

    [Fact]
    public void PerturbedSetFittedRmsdNeverExceedsDirectRmsd()
    {
        var target = CrambinCa();
        var rng = new Random(11);
        var noisy = target.Select(p => p + new Vec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5)).ToArray();
        double direct = Superposition.Rmsd(noisy, target);
        var fit = Superposition.Kabsch(noisy, target);
        Assert.True(fit.Rmsd <= direct + 1e-12);
        Assert.InRange(fit.Rmsd, 0.2, 0.5); // uniform ±0.5 Å noise: per-axis σ² = 1/12 → RMSD ≈ 0.5
        // The RMSD reported by the eigenvalue formula equals the RMSD after applying the transform.
        var moved = noisy.ToArray();
        fit.Apply(moved);
        Assert.Equal(fit.Rmsd, Superposition.Rmsd(moved, target), 9);
        // Any other rigid transform does worse.
        var other = new SuperpositionResult(Matrix3.FromAxisAngle(Vec3.UnitZ, 1) * fit.Rotation, fit.Translation, 0);
        var moved2 = noisy.Select(other.Apply).ToArray();
        Assert.True(Superposition.Rmsd(moved2, target) > fit.Rmsd);
    }

    [Fact]
    public void MirrorImageIsNotFittedByAReflection()
    {
        var target = CrambinCa();
        var mirror = target.Select(p => new Vec3(-p.X, p.Y, p.Z)).ToArray();
        var fit = Superposition.Kabsch(mirror, target);
        Assert.Equal(1, fit.Rotation.Determinant, 9);
        Assert.True(fit.Rmsd > 1.0, $"a chiral trace cannot be superposed on its mirror image (RMSD {fit.Rmsd:F2})");
    }

    [Fact]
    public void ZeroWeightsIgnoreOutliers()
    {
        var target = CrambinCa();
        var rotation = Matrix3.FromAxisAngle(new Vec3(0, 1, 1), -80);
        var mobile = target.Select(p => rotation * p).ToArray();
        mobile[0] += new Vec3(30, 0, 0);
        mobile[5] += new Vec3(0, -25, 0);
        var weights = Enumerable.Repeat(1.0, target.Length).ToArray();
        weights[0] = weights[5] = 0;
        var fit = Superposition.Kabsch(mobile, target, weights);
        Assert.Equal(0, fit.Rmsd, 6);
        Assert.True(Vec3.Distance(fit.Apply(mobile[10]), target[10]) < 1e-6);
        Assert.True(Superposition.Kabsch(mobile, target).Rmsd > 1);
        Assert.Throws<ArgumentException>(() => Superposition.Kabsch(mobile, target[..10]));
        Assert.Throws<ArgumentException>(() => Superposition.Kabsch(mobile, target, new double[target.Length]));
    }

    [Fact]
    public void ChainSuperpositionAndTmScoreOfSelf()
    {
        var a = SampleStructures.Crambin();
        var b = SampleStructures.Crambin();
        b.FirstModel.Transform(Matrix3.FromAxisAngle(Vec3.UnitX, 120), new Vec3(5, 5, 5));
        var fit = Superposition.SuperposeCa(b["A"], a["A"]);
        Assert.Equal(0, fit.Rmsd, 6);
        fit.Apply(b.FirstModel);
        Assert.True(a.Atoms.Zip(b.Atoms).All(p => Vec3.Distance(p.First.Position, p.Second.Position) < 1e-6));

        var tm = StructuralAlignment.TmScore(b["A"], a["A"]);
        Assert.Equal(1.0, tm.TmScore, 6);
        Assert.Equal(46, tm.PairsWithin5A);
        Assert.Equal(StructuralAlignment.D0(46), tm.D0);
    }
}

public class SasaTests
{
    [Fact]
    public void IsolatedAtomExposesItsWholeSphere()
    {
        double[] area = ShrakeRupley.Compute([new Vec3(1, 2, 3)], [1.7], probe: 1.4, points: 100);
        Assert.Equal(4 * Math.PI * 3.1 * 3.1, area[0], 9); // 120.76 Å²
        var atom = new Atom("CA", "C", Vec3.Zero);
        Assert.Equal(4 * Math.PI * 3.1 * 3.1, ShrakeRupley.Compute([atom])[0], 9);
        Assert.Equal(4 * Math.PI * 1.7 * 1.7, ShrakeRupley.Compute([atom], probe: 0)[0], 9);
    }

    [Fact]
    public void TwoOverlappingSpheresMatchCapFormula()
    {
        // Exposed area of each of two equal spheres (radius R) at distance d is 2πR(R + d/2).
        const double r = 1.5, probe = 1.4, d = 3.0;
        double big = r + probe;
        double[] area = ShrakeRupley.Compute([Vec3.Zero, new Vec3(d, 0, 0)], [r, r], probe, points: 5000);
        double expected = 2 * Math.PI * big * (big + d / 2);
        Assert.Equal(expected, area[0], expected * 0.005);
        Assert.Equal(area[0], area[1], expected * 0.005);
        // Far apart: no occlusion.
        double[] apart = ShrakeRupley.Compute([Vec3.Zero, new Vec3(10, 0, 0)], [r, r], probe);
        Assert.Equal(4 * Math.PI * big * big, apart[0], 9);
    }

    [Fact]
    public void SpherePointsAreUnitVectorsCentredOnOrigin()
    {
        var pts = ShrakeRupley.SpherePoints(500);
        Assert.Equal(500, pts.Length);
        Assert.All(pts, p => Assert.Equal(1.0, p.Length, 12));
        var mean = pts.Aggregate(Vec3.Zero, (s, p) => s + p) / pts.Length;
        Assert.True(mean.Length < 0.01);
    }

    [Fact]
    public void CrambinTotalAndPerResidueAreas()
    {
        var crambin = SampleStructures.Crambin();
        double total = ShrakeRupley.Total(crambin.FirstModel);
        // Crambin (327 heavy atoms, no waters in the deposited model) has an SASA close to 3000 Å² with Bondi-type
        // radii and a 1.4 Å probe; allow a generous window since the exact value depends on the radius set.
        Assert.InRange(total, 2500, 3500);

        var perResidue = ShrakeRupley.PerResidue(crambin.FirstModel);
        Assert.Equal(46, perResidue.Count);
        Assert.Equal(total, perResidue.Sum(r => r.Absolute), 6);
        Assert.All(perResidue, r =>
        {
            Assert.True(r.Absolute >= 0);
            Assert.InRange(r.Relative, 0, 1.6);
        });
        // The disulfide-bonded cysteines of this small, tightly packed protein are mostly buried.
        Assert.True(perResidue.Count(r => r.Residue.Name == "CYS" && r.IsBuried) >= 4);
        // Deterministic despite the parallel loop.
        Assert.Equal(total, ShrakeRupley.Total(crambin.FirstModel), 9);
    }

    [Fact]
    public void SasaIsTranslationInvariant()
    {
        var a = SampleStructures.Crambin();
        var b = SampleStructures.Crambin();
        b.FirstModel.Transform(Matrix3.Identity, new Vec3(100, -50, 25));
        Assert.Equal(ShrakeRupley.Total(a.FirstModel), ShrakeRupley.Total(b.FirstModel), 0.5);
        var rotated = SampleStructures.Crambin();
        rotated.FirstModel.Transform(Matrix3.FromAxisAngle(new Vec3(1, 1, 0), 63), Vec3.Zero);
        // Rotation changes only the dot discretisation.
        Assert.Equal(ShrakeRupley.Total(a.FirstModel), ShrakeRupley.Total(rotated.FirstModel), 60.0);
    }
}

public class ContactTests
{
    private static Vec3[] RandomPoints(int n, double box, int seed)
    {
        var rng = new Random(seed);
        return Enumerable.Range(0, n).Select(_ => new Vec3(rng.NextDouble() * box, rng.NextDouble() * box, rng.NextDouble() * box)).ToArray();
    }

    [Theory]
    [InlineData(1.0, 3.0)]
    [InlineData(5.0, 3.0)]
    [InlineData(2.0, 7.5)]
    public void NeighborSearchMatchesBruteForce(double cellSize, double radius)
    {
        var pts = RandomPoints(600, 30, 3);
        var ns = new NeighborSearch(pts, cellSize);
        Assert.Equal(600, ns.Count);

        var expectedPairs = new HashSet<(int, int)>();
        for (int i = 0; i < pts.Length; i++)
            for (int j = i + 1; j < pts.Length; j++)
                if (Vec3.Distance(pts[i], pts[j]) <= radius) expectedPairs.Add((i, j));
        var pairs = ns.AllPairs(radius);
        Assert.Equal(expectedPairs.Count, pairs.Count); // no duplicates
        Assert.True(expectedPairs.SetEquals(pairs));

        var rng = new Random(5);
        for (int q = 0; q < 50; q++)
        {
            var c = new Vec3(rng.NextDouble() * 34 - 2, rng.NextDouble() * 34 - 2, rng.NextDouble() * 34 - 2);
            var brute = Enumerable.Range(0, pts.Length).Where(i => Vec3.Distance(pts[i], c) <= radius).ToArray();
            Assert.Equal(brute, ns.Query(c, radius));
            Assert.Equal(brute.Length > 0, ns.Any(c, radius));
            var (idx, dist) = ns.Nearest(c);
            double bruteMin = pts.Min(p => Vec3.Distance(p, c));
            Assert.Equal(bruteMin, dist, 12);
            Assert.Equal(bruteMin, Vec3.Distance(pts[idx], c), 12);
        }
    }

    [Fact]
    public void EmptyNeighborSearch()
    {
        var ns = new NeighborSearch(ReadOnlySpan<Vec3>.Empty);
        Assert.Empty(ns.Query(Vec3.Zero, 10));
        Assert.False(ns.Any(Vec3.Zero, 10));
        Assert.Equal(-1, ns.Nearest(Vec3.Zero).Index);
        Assert.Empty(ns.AllPairs(5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NeighborSearch([Vec3.Zero], 0));
    }

    [Theory]
    [InlineData(ContactMode.CAlpha, 8.0)]
    [InlineData(ContactMode.CBeta, 8.0)]
    [InlineData(ContactMode.AnyHeavyAtom, 4.5)]
    public void CrambinContactMapMatchesBruteForce(ContactMode mode, double cutoff)
    {
        var chain = SampleStructures.Crambin()["A"];
        var map = ContactMap.Compute(chain, mode: mode);
        Assert.Equal(cutoff, map.Cutoff);
        Assert.Equal(46, map.Size);
        var res = chain.AminoAcids();

        IEnumerable<Vec3> Repr(Residue r) => mode switch
        {
            ContactMode.CAlpha => [r["CA"].Position],
            ContactMode.CBeta => [(r.GetAtom("CB") ?? r["CA"]).Position],
            _ => r.Atoms.Where(a => a.Element != "H").Select(a => a.Position),
        };

        int count = 0;
        for (int i = 0; i < res.Count; i++)
        {
            Assert.False(map[i, i]);
            for (int j = i + 1; j < res.Count; j++)
            {
                bool expected = Repr(res[i]).Any(p => Repr(res[j]).Any(q => Vec3.Distance(p, q) <= cutoff));
                Assert.Equal(expected, map[i, j]);
                Assert.Equal(map[i, j], map[j, i]);
                if (expected) count++;
            }
        }
        Assert.Equal(count, map.ContactCount);
        Assert.Equal(count, map.Pairs().Count());
        // Consecutive CA atoms (3.8 Å) are always in contact.
        if (mode == ContactMode.CAlpha) for (int i = 0; i + 1 < res.Count; i++) Assert.True(map[i, i + 1]);
    }

    [Fact]
    public void DistanceMatrixAndContactOrder()
    {
        var chain = SampleStructures.Crambin()["A"];
        var d = ContactMap.DistanceMatrix(chain);
        var ca = Superposition.CaTrace(chain);
        Assert.Equal(46, d.GetLength(0));
        Assert.Equal(Vec3.Distance(ca[3], ca[40]), d[3, 40], 12);
        Assert.Equal(d[3, 40], d[40, 3]);
        Assert.InRange(d[0, 1], 3.7, 3.9);
        var map = ContactMap.Compute(chain);
        double rco = map.RelativeContactOrder(minSeparation: 3);
        Assert.InRange(rco, 0.05, 0.5);
        Assert.All(map.Pairs(minSeparation: 3), p => Assert.True(p.J - p.I >= 3));
    }

    [Fact]
    public void CrambinDisulfidesFromGeometryMatchSsbondRecords()
    {
        var pdb = SampleStructures.Crambin();
        var cif = MmCifReader.Read(TestData.PathOf("1crn.cif"));
        foreach (var s in new[] { pdb, cif })
        {
            var found = Disulfides.Find(s);
            Assert.Equal(3, found.Count);
            var pairs = found.Select(x => (Math.Min(x.CysteineA.Number, x.CysteineB.Number), Math.Max(x.CysteineA.Number, x.CysteineB.Number)))
                .Order().ToArray();
            Assert.Equal([(3, 40), (4, 32), (16, 26)], pairs);
            Assert.All(found, x => Assert.InRange(x.Distance, 1.9, 2.15));
            Assert.All(found, x => Assert.InRange(Math.Abs(x.Chi3), 50, 130));
            Assert.Equal(pdb.Disulfides.Select(r => (r.NumberA, r.NumberB)).Order(), pairs);
        }
    }
}

public class PocketTests
{
    [Fact]
    public void TrypsinS1PocketHoldsBenzamidineAndIsLinedByAsp189()
    {
        var trypsin = TestData.Trypsin();
        var benzamidine = trypsin.Ligands.Single(r => r.Name == "BEN");
        var benCentre = Geometry.Centroid(benzamidine.Atoms.Select(a => a.Position).ToArray());

        var pockets = PocketDetector.Detect(trypsin);
        Assert.NotEmpty(pockets);
        Assert.Equal(Enumerable.Range(1, pockets.Count), pockets.Select(p => p.Rank));
        Assert.True(pockets.Zip(pockets.Skip(1)).All(p => p.First.Volume >= p.Second.Volume));

        // HETATMs are ignored by default, so the ligand-free S1 site is found as a pocket.
        var s1 = pockets.FirstOrDefault(p => p.Contains(benCentre, 2.0));
        Assert.NotNull(s1);
        Assert.True(s1.Rank <= 3, $"S1 pocket ranked {s1.Rank}");
        Assert.Contains(s1.LiningResidues, r => r.Name == "ASP" && r.Number == 189);
        Assert.Contains(s1.LiningResidues, r => r.Name == "SER" && r.Number == 190);
        Assert.Contains(s1.LiningResidues, r => r.Name == "GLY" && r.Number == 219);
        Assert.DoesNotContain(s1.LiningResidues, r => r.IsHetero);
        Assert.InRange(s1.MeanBuriedness, 5, 7);
        // Most benzamidine atoms sit on or next to pocket grid points.
        Assert.True(benzamidine.Atoms.Count(a => s1.Contains(a.Position, 1.5)) >= 6);
    }

    [Fact]
    public void PocketDetectionIsDeterministicAndRespectsOptions()
    {
        var trypsin = TestData.Trypsin();
        var a = PocketDetector.Detect(trypsin);
        var b = PocketDetector.Detect(trypsin);
        Assert.Equal(a.Select(p => p.Volume), b.Select(p => p.Volume));
        Assert.Equal(a.Select(p => p.Center), b.Select(p => p.Center));
        var strict = PocketDetector.Detect(trypsin, new PocketOptions { MinBuriedness = 7, MinPocketPoints = 10 });
        Assert.True(strict.Sum(p => p.Volume) < a.Sum(p => p.Volume));
        Assert.Empty(PocketDetector.Detect(Array.Empty<Atom>()));
    }

    [Fact]
    public void SingleExtendedStrandHasNoLargePocket()
    {
        var strand = SampleStructures.AlanineStrand(6);
        Assert.Empty(PocketDetector.Detect(strand, new PocketOptions { MinPocketPoints = 20 }));
    }
}

public class MmCifStructureTests
{
    [Fact]
    public void CrambinMmCifChainAndResidues()
    {
        var cif = MmCifReader.Read(TestData.PathOf("1crn.cif"));
        var chain = Assert.Single(cif.Chains);
        Assert.Equal("A", chain.Id);
        Assert.Equal(46, chain.Residues.Count);
        Assert.Equal(46, chain.AminoAcids().Count);
        Assert.Equal("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", chain.ToSequence().ToString());
        Assert.Equal(327, cif.AtomCount);
        Assert.Equal(321, cif.Atoms.Count(a => a.Element is "C" or "N" or "O")); // heavy atoms only, no hydrogens
        Assert.Equal(6,cif.Atoms.Count(a => a.Element == "S")); // six cysteines, no methionine
        Assert.Equal(Enumerable.Range(1, 46), chain.Residues.Select(r => r.Number));
        Assert.All(cif.Atoms, a => Assert.False(a.IsHetero));
    }
}
