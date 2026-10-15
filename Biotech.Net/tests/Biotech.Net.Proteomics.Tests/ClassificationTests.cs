using System.Text;
using Biotech.Net.Proteomics.Classification;

namespace Biotech.Net.Proteomics.Tests;

public class ClassificationTests
{
    /// <summary>Three well-separated 2-D Gaussian blobs (seeded).</summary>
    private static (List<double[]> X, List<int> Y) Blobs(int perClass = 30, int seed = 1)
    {
        var rng = new Random(seed);
        double[][] centres = [[0, 0], [6, 0], [0, 6]];
        var x = new List<double[]>();
        var y = new List<int>();
        for (int c = 0; c < centres.Length; c++)
            for (int i = 0; i < perClass; i++)
            {
                x.Add([centres[c][0] + Gaussian(rng), centres[c][1] + Gaussian(rng)]);
                y.Add(c);
            }
        return (x, y);
    }

    private static double Gaussian(Random rng) => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    [Fact]
    public void AminoAcidAndDipeptideComposition()
    {
        var aac = ProteinFeatures.AminoAcidComposition("AACX"u8);
        Assert.Equal(20, aac.Length);
        Assert.Equal(2.0 / 3, aac[AminoAcids.Standard.IndexOf('A')], 12);
        Assert.Equal(1.0 / 3, aac[AminoAcids.Standard.IndexOf('C')], 12);
        Assert.Equal(1.0, aac.Sum(), 12); // X is ignored

        var dpc = ProteinFeatures.DipeptideComposition("ACA"u8);
        Assert.Equal(400, dpc.Length);
        Assert.Equal(0.5, dpc[20 * 0 + 1], 12); // AC
        Assert.Equal(0.5, dpc[20 * 1 + 0], 12); // CA
        Assert.Equal(1.0, dpc.Sum(), 12);
        Assert.Equal("AC", ProteinFeatures.DecodeKmer(1, 2));
        Assert.Equal("CA", ProteinFeatures.DecodeKmer(20, 2));
    }

    [Fact]
    public void KmersSkipWindowsWithNonStandardResidues()
    {
        // "ACXDE": valid 2-mers are AC and DE only.
        var counts = ProteinFeatures.KmerCounts("ACXDE"u8, 2);
        Assert.Equal(2, counts.Count);
        Assert.Equal(["AC", "DE"], counts.Keys.Select(k => ProteinFeatures.DecodeKmer(k, 2)).Order(StringComparer.Ordinal));
        var tri = ProteinFeatures.KmerCounts("WWWW"u8, 3);
        Assert.Equal(2, Assert.Single(tri).Value);
        Assert.Equal(8000, ProteinFeatures.KmerSpectrum("ACDEFGHIK"u8, 3).Length);
        Assert.Throws<ArgumentOutOfRangeException>(() => ProteinFeatures.KmerSpectrum("AC"u8, 4));
        // Lower case counts like upper case.
        Assert.Equal(ProteinFeatures.KmerSpectrum("MKVLA"u8, 2), ProteinFeatures.KmerSpectrum("mkvla"u8, 2));
    }

    [Fact]
    public void PhysicochemicalAndCombinedFeatureLengths()
    {
        var seq = Encoding.ASCII.GetBytes("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE");
        var phys = ProteinFeatures.Physicochemical(seq);
        Assert.Equal(7, phys.Length);
        Assert.All(phys, v => Assert.True(double.IsFinite(v)));
        Assert.Equal(Math.Log10(seq.Length), phys[0], 12);
        Assert.Equal(427, ProteinFeatures.Combined(seq).Length);
    }

    [Fact]
    public void KmerSimilarities()
    {
        Assert.Equal(1.0, SequenceSimilarity.KmerJaccard("MKVLAAGIK", "MKVLAAGIK"), 12);
        Assert.Equal(1.0, SequenceSimilarity.KmerCosine("MKVLAAGIK", "MKVLAAGIK"), 12);
        Assert.Equal(0.0, SequenceSimilarity.KmerJaccard("AAAAA", "WWWWW"), 12);
        Assert.Equal(0.0, SequenceSimilarity.KmerCosine("AAAAA", "WWWWW"), 12);
        // ABCDE-style overlap: {ACD, CDE} vs {CDE, DEF} → 1/3.
        Assert.Equal(1.0 / 3, SequenceSimilarity.KmerJaccard("ACDE", "CDEF"), 12);
        Assert.Equal(0.5, SequenceSimilarity.KmerCosine("ACDE", "CDEF"), 12);

        string[] seqs = ["MKVLAAGIK", "MKVLAAGIR", "WWPPGGSS"];
        var m = SequenceSimilarity.CosineMatrix(seqs);
        for (int i = 0; i < 3; i++)
        {
            Assert.Equal(1.0, m[i * 3 + i], 12);
            for (int j = 0; j < 3; j++)
            {
                Assert.Equal(m[i * 3 + j], m[j * 3 + i], 12);
                if (i != j) Assert.Equal(SequenceSimilarity.KmerCosine(seqs[i], seqs[j]), m[i * 3 + j], 12);
            }
        }
        Assert.True(m[1] > m[2]);
    }

    [Theory]
    [InlineData(NeighborMetric.Euclidean)]
    [InlineData(NeighborMetric.Cosine)]
    public void KNearestNeighboursSeparatesBlobs(NeighborMetric metric)
    {
        var (x, y) = Blobs();
        // Shift away from the origin so cosine distance is meaningful for every blob.
        x = x.Select(v => new[] { v[0] + 10, v[1] + 10 }).ToList();
        var knn = new KNearestNeighborsClassifier(5, metric);
        knn.Fit(x, y);
        IProteinClassifier model = knn;
        Assert.Equal(0, model.Predict([10.0, 10.0]));
        Assert.Equal(1, model.Predict([16.0, 10.0]));
        Assert.Equal(2, model.Predict([10.0, 16.0]));
        var p = knn.PredictProbabilities([16.0, 10.0]);
        Assert.Equal(3, p.Length);
        Assert.Equal(1.0, p.Sum(), 5);
        Assert.True(p[1] > 0.99);
    }

    [Fact]
    public void LogisticRegressionFitsSeparableDataDeterministically()
    {
        var (x, y) = Blobs();
        var a = new LogisticRegressionClassifier();
        var b = new LogisticRegressionClassifier();
        a.Fit(x, y);
        b.Fit(x, y);
        int correct = 0;
        for (int i = 0; i < x.Count; i++)
        {
            var pa = a.PredictProbabilities(x[i]);
            Assert.Equal(pa, b.PredictProbabilities(x[i])); // bit-identical: zero init, no randomness
            Assert.Equal(1.0, pa.Sum(), 12);
            Assert.All(pa, v => Assert.InRange(v, 0, 1));
            if (((IProteinClassifier)a).Predict(x[i]) == y[i]) correct++;
        }
        Assert.True(correct >= x.Count - 2, $"training accuracy {correct}/{x.Count}");
        Assert.True(a.PredictProbabilities([6.0, 0.0])[1] > 0.9);
    }

    [Fact]
    public void ClassifiersValidateInput()
    {
        Assert.Throws<InvalidOperationException>(() => new LogisticRegressionClassifier().PredictProbabilities([1.0]));
        Assert.Throws<InvalidOperationException>(() => new KNearestNeighborsClassifier().PredictProbabilities([1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new KNearestNeighborsClassifier(0));
        var lr = new LogisticRegressionClassifier();
        Assert.Throws<ArgumentException>(() => lr.Fit([[1.0], [2.0]], [0]));
        Assert.Throws<ArgumentException>(() => lr.Fit([[1.0], [2.0, 3.0]], [0, 1]));
        Assert.Throws<ArgumentException>(() => lr.Fit([[1.0]], [-1]));
        lr.Fit([[0.0], [1.0]], [0, 1]);
        Assert.Throws<ArgumentException>(() => lr.PredictProbabilities([1.0, 2.0]));
    }

    [Fact]
    public void StratifiedFoldsPreserveClassProportionsAndAreSeeded()
    {
        var labels = Enumerable.Repeat(0, 20).Concat(Enumerable.Repeat(1, 10)).Concat(Enumerable.Repeat(2, 5)).ToArray();
        var folds = Classifiers.StratifiedFolds(labels, 5);
        Assert.Equal(folds, Classifiers.StratifiedFolds(labels, 5));
        Assert.NotEqual(folds, Classifiers.StratifiedFolds(labels, 5, seed: 99));
        for (int f = 0; f < 5; f++)
        {
            Assert.Equal(4, Enumerable.Range(0, labels.Length).Count(i => folds[i] == f && labels[i] == 0));
            Assert.Equal(2, Enumerable.Range(0, labels.Length).Count(i => folds[i] == f && labels[i] == 1));
            Assert.Equal(1, Enumerable.Range(0, labels.Length).Count(i => folds[i] == f && labels[i] == 2));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => Classifiers.StratifiedFolds(labels, 1));
    }

    [Fact]
    public void MacroF1FromConfusionMatrix()
    {
        // Class 0: TP 5, FP 1, FN 0 → F1 10/11. Class 1: TP 4, FP 0, FN 1 → F1 8/9.
        var r = new CrossValidationResult([0.9], 0.9, [5, 0, 1, 4], 2);
        Assert.Equal((10.0 / 11 + 8.0 / 9) / 2, r.MacroF1, 12);
    }

    [Fact]
    public void DemoDatasetIsReproducible()
    {
        var a = ProteinFamilyDemo.Generate(5, seed: 3);
        var b = ProteinFamilyDemo.Generate(5, seed: 3);
        var c = ProteinFamilyDemo.Generate(5, seed: 4);
        Assert.Equal(20, a.Count);
        Assert.Equal(a.Select(e => e.Protein.Sequence.ToString()), b.Select(e => e.Protein.Sequence.ToString()));
        Assert.NotEqual(a.Select(e => e.Protein.Sequence.ToString()), c.Select(e => e.Protein.Sequence.ToString()));
        Assert.All(a, e => Assert.StartsWith("M", e.Protein.Sequence.ToString()));
        Assert.Equal(ProteinFamilyDemo.Families, a.Select(e => e.Family).Distinct().Order(StringComparer.Ordinal));
    }

    [Fact]
    public void FamilyClassifierLearnsSyntheticFamilies()
    {
        var train = ProteinFamilyDemo.Generate(15, seed: 1);
        var test = ProteinFamilyDemo.Generate(10, seed: 2);
        var clf = new ProteinFamilyClassifier().Fit(train);
        Assert.Equal(ProteinFamilyDemo.Families, clf.Classes);
        int correct = test.Count(e => clf.Predict(e.Protein) == e.Family);
        Assert.True(correct >= 0.9 * test.Count, $"held-out accuracy {correct}/{test.Count}");
        var probs = clf.PredictProbabilities(test[0].Protein);
        Assert.Equal(4, probs.Count);
        Assert.Equal(1.0, probs.Values.Sum(), 9);
        Assert.Equal(test[0].Family, probs.MaxBy(kv => kv.Value).Key);
    }

    [Fact]
    public void CrossValidationIsDeterministicAndAccurate()
    {
        var data = ProteinFamilyDemo.Generate(12, seed: 5);
        var r1 = ProteinFamilyClassifier.CrossValidate(data, () => new KNearestNeighborsClassifier(3), b => ProteinFeatures.AminoAcidComposition(b), folds: 4);
        var r2 = ProteinFamilyClassifier.CrossValidate(data, () => new KNearestNeighborsClassifier(3), b => ProteinFeatures.AminoAcidComposition(b), folds: 4);
        Assert.Equal(r1.FoldAccuracies, r2.FoldAccuracies);
        Assert.Equal(r1.ConfusionMatrix, r2.ConfusionMatrix);
        Assert.Equal(4, r1.Classes);
        Assert.Equal(data.Count, r1.ConfusionMatrix.Sum());
        Assert.True(r1.Accuracy >= 0.85, $"CV accuracy {r1.Accuracy:F3}");
        Assert.InRange(r1.MacroF1, 0.8, 1.0);
    }
}
