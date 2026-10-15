using Biotech.Net.Chemistry;

namespace Biotech.Net.Chemistry.Tests;

public class BitVectorTests
{
    [Fact]
    public void Set_get_count_and_onbits()
    {
        var v = new BitVector(130, [0, 64, 129]);
        Assert.True(v[64]);
        Assert.False(v[65]);
        Assert.Equal(3, v.Count);
        Assert.Equal([0, 64, 129], v.OnBits().ToArray());
        v.Clear(64);
        Assert.Equal(2, v.Count);
        Assert.Throws<ArgumentOutOfRangeException>(() => v.Set(130));
    }

    [Fact]
    public void Serialization_round_trips()
    {
        var v = new BitVector(100, [1, 5, 99]);
        Assert.Equal(v, BitVector.FromBitString(v.ToBitString()));
        Assert.Equal(v, BitVector.FromBase64(v.ToBase64()));
    }

    [Fact]
    public void Bitwise_operators_and_fold()
    {
        var a = new BitVector(8, [0, 1, 2]);
        var b = new BitVector(8, [2, 3]);
        Assert.Equal([2], (a & b).OnBits().ToArray());
        Assert.Equal([0, 1, 2, 3], (a | b).OnBits().ToArray());
        Assert.Equal([0, 1, 3], (a ^ b).OnBits().ToArray());
        Assert.Equal([0, 1], new BitVector(8, [0, 4, 5]).Fold(4).OnBits().ToArray());
    }

    [Fact]
    public void Similarity_coefficients()
    {
        // a = 4 bits, b = 3 bits, c = 2 common
        var x = new BitVector(16, [0, 1, 2, 3]);
        var y = new BitVector(16, [2, 3, 4]);
        Assert.Equal(2.0 / 5, Similarity.Tanimoto(x, y), 12);
        Assert.Equal(4.0 / 7, Similarity.Dice(x, y), 12);
        Assert.Equal(2 / Math.Sqrt(12), Similarity.Cosine(x, y), 12);
        Assert.Equal(Similarity.Tanimoto(x, y), Similarity.Tversky(x, y, 1, 1), 12);
        Assert.Equal(Similarity.Dice(x, y), Similarity.Tversky(x, y, 0.5, 0.5), 12);
        Assert.Equal(2.0 / 4, Similarity.Tversky(x, y, 1, 0), 12); // fraction of x found in y
        Assert.Equal(2.0 / (8 + 6 - 6), Similarity.Sokal(x, y), 12);
        Assert.Equal(0, Similarity.Tanimoto(new BitVector(8), new BitVector(8)));
        Assert.Throws<ArgumentException>(() => Similarity.Tanimoto(new BitVector(8), new BitVector(16)));
    }

    [Fact]
    public void Count_vector_similarity()
    {
        var a = new CountVector(); a.Add(1, 2); a.Add(2, 1);
        var b = new CountVector(); b.Add(1, 1); b.Add(3, 1);
        Assert.Equal(1.0 / 4, Similarity.Tanimoto(a, b), 12); // min sum 1, max sum 2+1+1
        Assert.Equal(2.0 * 1 / 5, Similarity.Dice(a, b), 12);
    }
}

public class FingerprintTests
{
    private static readonly string[] Drugs =
    [
        "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(=O)Nc1ccc(O)cc1",
        "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    ];

    [Fact]
    public void Benzene_morgan_has_three_bits_like_rdkit()
        => Assert.Equal(3, MorganFingerprint.Ecfp4(Molecule.FromSmiles("c1ccccc1")).Count);

    [Fact]
    public void Aspirin_ecfp4_bits_match_rdkit()
    {
        int[] rdkit = [389, 456, 650, 695, 807, 909, 1017, 1035, 1047, 1057, 1088, 1199, 1380, 1410, 1447, 1468, 1616, 1729, 1750, 1775, 1873, 1917, 1970, 1991];
        Assert.Equal(rdkit, MorganFingerprint.Ecfp4(Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O")).OnBits().ToArray());
    }

    [Theory]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1O", 0.4482758620689655)]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 0.1951219512195122)]
    [InlineData("c1ccccc1", "Cc1ccccc1", 0.2727272727272727)]
    public void Ecfp4_tanimoto_matches_rdkit(string a, string b, double expected)
        => Assert.Equal(expected, Similarity.Tanimoto(MorganFingerprint.Ecfp4(Molecule.FromSmiles(a)), MorganFingerprint.Ecfp4(Molecule.FromSmiles(b))), 9);

    [Fact]
    public void Morgan_bit_info_explains_every_bit()
    {
        var m = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        var fp = MorganFingerprint.Ecfp4(m);
        var info = MorganFingerprint.BitInfo(m);
        Assert.Equal(fp.OnBits().ToHashSet(), info.Keys.ToHashSet());
        Assert.All(info.Values.SelectMany(v => v), e => Assert.InRange(e.Radius, 0, 2));
    }

    [Fact]
    public void Every_fingerprint_type_is_deterministic_and_self_similar()
    {
        foreach (FingerprintType t in Enum.GetValues<FingerprintType>())
            foreach (var s in Drugs)
            {
                var a = Fingerprint.Compute(Molecule.FromSmiles(s), t);
                var b = Fingerprint.Compute(Molecule.FromSmiles(Molecule.FromSmiles(s).ToSmiles()), t);
                Assert.True(a.Count > 0, $"{t} empty for {s}");
                Assert.Equal(a, b);
                Assert.Equal(1.0, Similarity.Tanimoto(a, b));
            }
    }

    [Fact]
    public void Maccs_has_167_bits_and_key_166_marks_mixtures()
    {
        var fp = MaccsKeys.Compute(Molecule.FromSmiles("CC(=O)[O-].[Na+]"));
        Assert.Equal(167, fp.Length);
        Assert.True(fp[166]);
        Assert.False(MaccsKeys.Compute(Molecule.FromSmiles("CCO"))[166]);
    }

    [Fact]
    public void Path_fingerprint_reflects_substructure()
    {
        var benzene = PathFingerprint.Compute(Molecule.FromSmiles("c1ccccc1"));
        var toluene = PathFingerprint.Compute(Molecule.FromSmiles("Cc1ccccc1"));
        // every benzene path is also a toluene path
        Assert.Equal(benzene.Count, BitVector.IntersectionCount(benzene, toluene));
        var linear = PathFingerprint.Compute(Molecule.FromSmiles("CC(C)C"), new PathFingerprint.Options { BranchedPaths = false });
        var branched = PathFingerprint.Compute(Molecule.FromSmiles("CC(C)C"));
        Assert.True(branched.Count > linear.Count);
    }

    [Fact]
    public void Subgraph_enumeration_counts()
    {
        // isobutane: 3 bonds all meeting at the centre: 3 single bonds, 3 pairs, 1 triple
        Assert.Equal(7, Subgraphs.Enumerate(Molecule.FromSmiles("CC(C)C"), 1, 7).Count());
        Assert.Equal(6, Subgraphs.Enumerate(Molecule.FromSmiles("CC(C)C"), 1, 7, branched: false).Count());
        // cyclopropane: 3 + 3 + 1
        Assert.Equal(7, Subgraphs.Enumerate(Molecule.FromSmiles("C1CC1"), 1, 7).Count());
    }

    [Fact]
    public void Fingerprint_index_threshold_and_topk_search()
    {
        var index = new FingerprintIndex();
        index.AddRange(Drugs.Select((s, i) => Molecule.FromSmiles(s, $"m{i}")));
        var q = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        var top = index.TopK(q, 2);
        Assert.Equal("m0", top[0].Id);
        Assert.Equal(1.0, top[0].Score);
        var hits = index.Search(q, threshold: 0.99);
        Assert.Single(hits);
        var all = index.Search(q);
        Assert.Equal(Drugs.Length, all.Count);
        Assert.True(all.Zip(all.Skip(1)).All(p => p.First.Score >= p.Second.Score));
    }

    [Fact]
    public void Bulk_and_matrix_similarity()
    {
        var fps = Drugs.Select(s => MorganFingerprint.Ecfp4(Molecule.FromSmiles(s))).ToList();
        var bulk = Similarity.BulkTanimoto(fps[0], fps);
        Assert.Equal(1.0, bulk[0]);
        var m = Similarity.Matrix(fps);
        for (int i = 0; i < fps.Count; i++)
            for (int j = 0; j < fps.Count; j++)
            {
                Assert.Equal(m[i, j], m[j, i]);
                Assert.Equal(Similarity.Tanimoto(fps[i], fps[j]), m[i, j], 12);
            }
    }
}
