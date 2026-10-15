namespace Biotech.Net.Transcriptomics.Tests;

public class PcaTests
{
    /// <summary>
    /// X = mean + s1·u1ᵀ + s2·u2ᵀ with zero-mean orthogonal score vectors s1 = (3, 3, −3, −3) and
    /// s2 = (1, −1, 1, −1) and orthonormal axes u1, u2: variances 36/3 = 12 and 4/3, ratios 0.9 and 0.1.
    /// </summary>
    private static DenseMatrix Known(double[] u1, double[] u2, double[] mean)
    {
        double[] s1 = [3, 3, -3, -3], s2 = [1, -1, 1, -1];
        int p = u1.Length;
        var m = new DenseMatrix(4, p);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < p; j++) m[i, j] = mean[j] + s1[i] * u1[j] + s2[i] * u2[j];
        return m;
    }

    private static void AssertKnown(PcaResult r, double[] u1, double[] u2, double[] mean)
    {
        Assert.Equal(12, r.ExplainedVariance[0], 1e-9);
        Assert.Equal(4.0 / 3, r.ExplainedVariance[1], 1e-9);
        Assert.Equal(0.9, r.ExplainedVarianceRatio[0], 1e-9);
        Assert.Equal(0.1, r.ExplainedVarianceRatio[1], 1e-9);
        for (int j = 0; j < mean.Length; j++) Assert.Equal(mean[j], r.Means[j], 1e-12);
        // Axes up to sign; the sign convention makes the largest |loading| positive.
        double[][] axes = [u1, u2];
        double[][] scores = [[3, 3, -3, -3], [1, -1, 1, -1]];
        for (int c = 0; c < 2; c++)
        {
            double dot = 0;
            for (int j = 0; j < u1.Length; j++) dot += r.Loadings[j, c] * axes[c][j];
            Assert.Equal(1, Math.Abs(dot), 1e-9);
            double sign = Math.Sign(dot);
            for (int i = 0; i < 4; i++) Assert.Equal(sign * scores[c][i], r.Scores[i, c], 1e-9);
            int argmax = Enumerable.Range(0, u1.Length).OrderByDescending(j => Math.Abs(r.Loadings[j, c])).First();
            Assert.True(r.Loadings[argmax, c] > 0);
        }
    }

    [Fact]
    public void Full_RecoversKnownAxes_FewFeatures()
    {
        double[] u1 = [1 / 3.0, 2 / 3.0, 2 / 3.0], u2 = [2 / 3.0, 1 / 3.0, -2 / 3.0], mean = [5, -1, 2];
        var r = Pca.Fit(Known(u1, u2, mean), components: 2, solver: PcaSolver.Full);
        Assert.Equal(2, r.Components);
        AssertKnown(r, u1, u2, mean);
    }

    [Fact]
    public void Full_RecoversKnownAxes_MoreFeaturesThanSamples()
    {
        // Gram-matrix branch (p > n).
        double[] u1 = [0.5, 0.5, 0.5, 0.5, 0], u2 = [0.5, -0.5, 0.5, -0.5, 0], mean = [1, 2, 3, 4, 5];
        var r = Pca.Fit(Known(u1, u2, mean), components: 2, solver: PcaSolver.Full);
        AssertKnown(r, u1, u2, mean);
    }

    [Fact]
    public void Randomized_AgreesWithFull()
    {
        var rng = new Random(7);
        int n = 300, p = 80;
        var data = new DenseMatrix(n, p);
        // Three strong latent factors plus noise.
        var load = new double[3, p];
        for (int c = 0; c < 3; c++) for (int j = 0; j < p; j++) load[c, j] = TestData.Normal(rng);
        for (int i = 0; i < n; i++)
        {
            double f0 = 5 * TestData.Normal(rng), f1 = 3 * TestData.Normal(rng), f2 = 2 * TestData.Normal(rng);
            for (int j = 0; j < p; j++) data[i, j] = f0 * load[0, j] + f1 * load[1, j] + f2 * load[2, j] + 0.1 * TestData.Normal(rng);
        }
        var full = Pca.Fit(data, components: 3, solver: PcaSolver.Full);
        var rand = Pca.Fit(data, components: 3, solver: PcaSolver.Randomized, seed: 3);
        for (int c = 0; c < 3; c++)
        {
            Assert.Equal(full.ExplainedVariance[c], rand.ExplainedVariance[c], full.ExplainedVariance[c] * 1e-6);
            double dot = 0;
            for (int j = 0; j < p; j++) dot += full.Loadings[j, c] * rand.Loadings[j, c];
            Assert.Equal(1, dot, 1e-6); // same sign convention
        }
        Assert.True(full.ExplainedVarianceRatio.Sum() > 0.99);
        // Scores are the centred data projected on the loadings, so their variance is the explained variance.
        for (int c = 0; c < 3; c++)
        {
            var col = full.Scores.Column(c);
            double mean = col.Average();
            Assert.Equal(0, mean, 1e-9);
            Assert.Equal(full.ExplainedVariance[c], col.Sum(v => v * v) / (n - 1), 1e-6);
        }
    }
}
