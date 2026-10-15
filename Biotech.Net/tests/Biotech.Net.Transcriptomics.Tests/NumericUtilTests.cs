using Biotech.Net.Native;
using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics.Tests;

public class NumericUtilTests
{
    [Theory]
    // Reference values from scipy.special.digamma / polygamma.
    [InlineData(0.5, -1.9635100260214235, 4.93480220054468, -16.828796644234316)]
    [InlineData(1.0, -0.5772156649015329, 1.6449340668482266, -2.404113806319188)]
    [InlineData(2.5, 0.7031566406452432, 0.4903577561002349, -0.2362040516417274)]
    [InlineData(10.0, 2.251752589066721, 0.10516633568168575, -0.011049834970802069)]
    public void PolygammaFunctions_MatchScipy(double x, double digamma, double trigamma, double tetragamma)
    {
        Assert.Equal(digamma, NumericUtil.Digamma(x), 1e-12);
        Assert.Equal(trigamma, NumericUtil.Trigamma(x), 1e-12);
        Assert.Equal(tetragamma, NumericUtil.Tetragamma(x), 1e-11);
    }

    [Theory]
    [InlineData(0.1, 10.491681821078421)]
    [InlineData(1.0, 1.4262551202150786)]
    [InlineData(10.0, 0.3350810443779569)]
    public void TrigammaInverse_InvertsTrigamma(double x, double expected)
    {
        double y = NumericUtil.TrigammaInverse(x);
        Assert.Equal(expected, y, 1e-6);
        Assert.Equal(x, NumericUtil.Trigamma(y), 1e-8);
    }

    [Fact]
    public void SpecialFunctions_MatchScipy()
    {
        Assert.Equal(0.5723649429247001, SpecialFunctions.LogGamma(0.5), 1e-13);
        Assert.Equal(0.2967529892956665, SpecialFunctions.RegularizedIncompleteBeta(2.5, 3.5, 0.3), 1e-12);
        Assert.Equal(0.15164090963470994, SpecialFunctions.RegularizedIncompleteBeta(10, 0.5, 0.9), 1e-12);
        Assert.Equal(0.10193947882985835, SpecialFunctions.StudentTTwoSided(2.0, 5), 1e-12);
        Assert.Equal(0.004230535234034089, SpecialFunctions.StudentTTwoSided(3.5, 12.3), 1e-12);
        // Normal upper tail must keep relative precision far into the tail.
        Assert.Equal(2.8665157187919333e-07, NumericUtil.NormalUpper(5), 1e-18);
        Assert.Equal(6.22096057427174e-16, NumericUtil.NormalUpper(8), 1e-26);
        Assert.Equal(0.025, NumericUtil.NormalUpper(1.959963984540054), 1e-12);
    }

    [Fact]
    public void MedianQuantileMad_FollowR()
    {
        Assert.Equal(3, NumericUtil.Median([5, 1, 3]));
        Assert.Equal(2.5, NumericUtil.Median([4, 1, 3, 2]));
        Assert.Equal(2, NumericUtil.Median([double.NaN, 1, 3])); // NaN ignored
        // quantile(c(1,2,3,4,10), 0.3) = 2.2 (type 7)
        Assert.Equal(2.2, NumericUtil.Quantile([10, 4, 3, 2, 1], 0.3), 1e-12);
        Assert.Equal(10, NumericUtil.Quantile([10, 4, 3, 2, 1], 1.0));
        // mad(c(1,2,3,4,100)) = 1.4826
        Assert.Equal(1.4826, NumericUtil.Mad([1, 2, 3, 4, 100]), 1e-12);
    }

    [Fact]
    public void Rank_AveragesTies()
    {
        // rank(c(10, 20, 10, 5, 20, 20)) = 2.5 5 2.5 1 5 5
        Assert.Equal([2.5, 5, 2.5, 1, 5, 5], NumericUtil.Rank([10, 20, 10, 5, 20, 20]));
    }

    [Fact]
    public void Lowess_ReproducesLinesAndResistsOutliers()
    {
        var x = Enumerable.Range(1, 30).Select(i => (double)i).ToArray();
        var y = x.Select(v => 3 - 0.5 * v).ToArray();
        var fit = NumericUtil.Lowess(x, y);
        for (int i = 0; i < x.Length; i++) Assert.Equal(y[i], fit[i], 1e-9);

        var yOut = (double[])y.Clone();
        yOut[15] += 50; // gross outlier
        var robust = NumericUtil.Lowess(x, yOut);
        Assert.Equal(y[15], robust[15], 1e-6);
        var nonRobust = NumericUtil.Lowess(x, yOut, iterations: 0);
        Assert.True(nonRobust[15] - y[15] > 1, "Without robustness iterations the outlier pulls the fit.");
    }

    [Fact]
    public void Lowess_MatchesRCars()
    {
        // R: lowess(cars)$y (f = 2/3, iter = 3, delta = 0.01 * diff(range(speed))).
        double[] speed = [4, 4, 7, 7, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 22, 23, 24, 24, 24, 24, 25];
        double[] dist = [2, 10, 4, 22, 16, 10, 18, 26, 34, 17, 28, 14, 20, 24, 28, 26, 34, 34, 46, 26, 36, 60, 80, 20, 26, 54, 32, 40, 32, 40, 50, 42, 56, 76, 84, 36, 46, 68, 32, 48, 52, 56, 64, 66, 54, 70, 92, 93, 120, 85];
        double[] expectedUnique = [4.965459, 13.124495, 15.858633, 18.579691, 21.280313, 24.129277, 27.119549, 30.027276, 32.962506, 36.757728,
            40.435075, 43.463492, 46.885479, 50.793152, 56.491224, 67.585824, 73.079695, 78.643164, 84.328698];
        double[] uniqueSpeed = [4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25];
        var fit = NumericUtil.Lowess(speed, dist);
        for (int u = 0; u < uniqueSpeed.Length; u++)
        {
            int i = Array.IndexOf(speed, uniqueSpeed[u]);
            Assert.Equal(expectedUnique[u], fit[i], 1e-5);
        }
    }

    [Fact]
    public void Loess_IsExactForQuadratics()
    {
        var r = new Random(3);
        var x = Enumerable.Range(0, 200).Select(_ => r.NextDouble() * 10).ToArray();
        var y = x.Select(v => 1 + 2 * v - 0.3 * v * v).ToArray();
        var fit = NumericUtil.Loess(x, y, 0.3, 2);
        for (int i = 0; i < x.Length; i++) Assert.Equal(y[i], fit[i], 1e-8);
    }

    [Fact]
    public void SmallLinearAlgebra_IsCorrect()
    {
        double[] a = [4, 2, 0, 2, 5, 1, 0, 1, 3];
        var x = NumericUtil.SolveSmall((double[])a.Clone(), [2, 9, 7], 3)!;
        // A·x = b
        Assert.Equal(2, 4 * x[0] + 2 * x[1], 1e-12);
        Assert.Equal(9, 2 * x[0] + 5 * x[1] + x[2], 1e-12);
        Assert.Equal(7, x[1] + 3 * x[2], 1e-12);

        var inv = NumericUtil.InvertSpd(a, 3)!;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                double s = 0;
                for (int k = 0; k < 3; k++) s += a[i * 3 + k] * inv[k * 3 + j];
                Assert.Equal(i == j ? 1 : 0, s, 1e-12);
            }
        // det = 4(15 − 1) − 2(6 − 0) = 44
        Assert.Equal(Math.Log(44), NumericUtil.LogDetSpd(a, 3), 1e-12);
        Assert.Null(NumericUtil.InvertSpd([1, 2, 2, 1], 2));
        Assert.True(double.IsNaN(NumericUtil.LogDetSpd([1, 2, 2, 1], 2)));
        Assert.Null(NumericUtil.SolveSmall([1, 2, 2, 4], [1, 2], 2));
    }

    [Fact]
    public void Interpolate_IsLinearWithConstantEnds()
    {
        double[] xk = [0, 1, 3], yk = [10, 20, 0];
        Assert.Equal(10, NumericUtil.Interpolate(xk, yk, -5));
        Assert.Equal(15, NumericUtil.Interpolate(xk, yk, 0.5), 1e-12);
        Assert.Equal(10, NumericUtil.Interpolate(xk, yk, 2), 1e-12);
        Assert.Equal(0, NumericUtil.Interpolate(xk, yk, 9));
    }
}
