using Biotech.Net.Chemistry;
using Biotech.Net.DrugDiscovery;

namespace Biotech.Net.DrugDiscovery.Tests;

public class RidgeRegressionTests
{
    [Fact]
    public void Recovers_an_exact_linear_relationship_without_regularisation()
    {
        // y = 3 + 2*x0 - 1*x1
        double[][] x = [[0, 0], [1, 0], [0, 1], [1, 1], [2, 1], [1, 2], [3, 0], [0, 3]];
        var y = x.Select(r => 3 + 2 * r[0] - r[1]).ToArray();
        var model = new RidgeRegression { Alpha = 0 }.Fit(x, y);
        Assert.Equal(3, model.Intercept, 6);
        Assert.Equal(2, model.Coefficients[0], 6);
        Assert.Equal(-1, model.Coefficients[1], 6);
        Assert.Equal(1.0, RegressionMetrics.R2(y, model.Predict(x)), 9);
    }

    [Fact]
    public void Regularisation_shrinks_the_coefficients()
    {
        double[][] x = [[0, 0], [1, 0], [0, 1], [1, 1], [2, 1], [1, 2]];
        var y = x.Select(r => 3 + 2 * r[0] - r[1]).ToArray();
        var plain = new RidgeRegression { Alpha = 0 }.Fit(x, y);
        var ridge = new RidgeRegression { Alpha = 10 }.Fit(x, y);
        Assert.True(Math.Abs(ridge.Coefficients[0]) < Math.Abs(plain.Coefficients[0]));
    }

    [Fact]
    public void Handles_collinear_features_without_blowing_up()
    {
        double[][] x = [[1, 2], [2, 4], [3, 6], [4, 8]]; // x1 = 2*x0
        double[] y = [1, 2, 3, 4];
        var model = new RidgeRegression { Alpha = 1 }.Fit(x, y);
        Assert.All(model.Predict(x), p => Assert.True(double.IsFinite(p)));
        Assert.True(RegressionMetrics.R2(y, model.Predict(x)) > 0.9);
    }

    [Fact]
    public void Rejects_bad_input_and_unfitted_use()
    {
        Assert.Throws<ArgumentException>(() => new RidgeRegression().Fit([[1.0]], [1.0, 2.0]));
        Assert.Throws<ArgumentException>(() => new RidgeRegression().Fit([], []));
        Assert.Throws<InvalidOperationException>(() => new RidgeRegression().Predict([1.0]));
    }
}

public class RandomForestTests
{
    private static (double[][] X, double[] Y) NonLinearData(int n = 120)
    {
        var rng = new Random(7);
        var x = new double[n][];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double a = rng.NextDouble() * 4 - 2, b = rng.NextDouble() * 4 - 2;
            x[i] = [a, b, rng.NextDouble()];
            y[i] = a * a + (b > 0 ? 2 : -2);
        }
        return (x, y);
    }

    [Fact]
    public void Fits_a_non_linear_function_better_than_a_linear_model()
    {
        var (x, y) = NonLinearData();
        var forest = new RandomForestRegressor { TreeCount = 60, Seed = 1 }.Fit(x, y);
        var ridge = new RidgeRegression { Alpha = 1 }.Fit(x, y);
        double rfR2 = RegressionMetrics.R2(y, forest.Predict(x));
        double ridgeR2 = RegressionMetrics.R2(y, ridge.Predict(x));
        Assert.True(rfR2 > 0.9, $"forest R² {rfR2:F3}");
        Assert.True(rfR2 > ridgeR2);
    }

    [Fact]
    public void Out_of_bag_r2_is_reported_and_reasonable()
    {
        var (x, y) = NonLinearData();
        var forest = new RandomForestRegressor { TreeCount = 60, Seed = 1 }.Fit(x, y);
        Assert.False(double.IsNaN(forest.OutOfBagR2));
        Assert.InRange(forest.OutOfBagR2, 0.5, 1.0);
    }

    [Fact]
    public void Feature_importances_find_the_informative_columns()
    {
        var (x, y) = NonLinearData();
        var forest = new RandomForestRegressor { TreeCount = 60, Seed = 1 }.Fit(x, y);
        var imp = forest.FeatureImportances(3);
        Assert.Equal(1.0, imp.Sum(), 6);
        Assert.True(imp[2] < imp[0]); // the noise column matters least
        Assert.True(imp[2] < imp[1]);
    }

    [Fact]
    public void Is_deterministic_for_a_given_seed()
    {
        var (x, y) = NonLinearData(60);
        var a = new RandomForestRegressor { TreeCount = 20, Seed = 5 }.Fit(x, y);
        var b = new RandomForestRegressor { TreeCount = 20, Seed = 5 }.Fit(x, y);
        Assert.Equal(a.Predict(x), b.Predict(x));
    }
}

public class MetricsAndCvTests
{
    [Fact]
    public void Metrics_on_a_known_pair()
    {
        double[] actual = [1, 2, 3, 4];
        double[] perfect = [1, 2, 3, 4];
        double[] shifted = [2, 3, 4, 5];
        Assert.Equal(1.0, RegressionMetrics.R2(actual, perfect));
        Assert.Equal(0.0, RegressionMetrics.Rmse(actual, perfect));
        Assert.Equal(1.0, RegressionMetrics.Rmse(actual, shifted), 9);
        Assert.Equal(1.0, RegressionMetrics.Mae(actual, shifted), 9);
        Assert.Equal(1.0, RegressionMetrics.Pearson(actual, shifted), 9);
        Assert.Equal(1.0, RegressionMetrics.Spearman(actual, [10, 20, 30, 40]), 9);
        Assert.Equal(-1.0, RegressionMetrics.Spearman(actual, [40, 30, 20, 10]), 9);
        // predicting the mean gives R² = 0
        Assert.Equal(0.0, RegressionMetrics.R2(actual, [2.5, 2.5, 2.5, 2.5]), 9);
    }

    [Fact]
    public void Cross_validation_recovers_a_learnable_relationship()
    {
        var rng = new Random(3);
        var x = Enumerable.Range(0, 60).Select(_ => new[] { rng.NextDouble() * 10, rng.NextDouble() * 10 }).ToArray();
        var y = x.Select(r => 1 + 2 * r[0] - 0.5 * r[1]).ToArray();
        var cv = CrossValidation.KFold(x, y, (a, b) => new RidgeRegression { Alpha = 0.01 }.Fit(a, b), folds: 5);
        Assert.True(cv.R2 > 0.98, cv.ToString());
        Assert.Equal(5, cv.FoldR2.Count);
        Assert.Equal(60, cv.Predictions.Length);
        Assert.True(cv.Rmse < 0.5);
    }

    [Fact]
    public void Cross_validation_does_not_reward_noise()
    {
        var rng = new Random(11);
        var x = Enumerable.Range(0, 40).Select(_ => new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() }).ToArray();
        var y = Enumerable.Range(0, 40).Select(_ => rng.NextDouble()).ToArray();
        var cv = CrossValidation.KFold(x, y, (a, b) => new RidgeRegression { Alpha = 0.001 }.Fit(a, b));
        Assert.True(cv.R2 < 0.5, cv.ToString()); // out-of-fold R² on noise must not look good
    }

    [Fact]
    public void Leave_one_out_predicts_every_row_once()
    {
        var x = Enumerable.Range(0, 12).Select(i => new[] { (double)i }).ToArray();
        var y = x.Select(r => 2 * r[0] + 1).ToArray();
        var cv = CrossValidation.LeaveOneOut(x, y, (a, b) => new RidgeRegression { Alpha = 0 }.Fit(a, b));
        Assert.Equal(12, cv.Predictions.Length);
        Assert.True(cv.R2 > 0.99);
    }

    [Fact]
    public void Fold_count_is_validated()
    {
        var x = new[] { new[] { 1.0 }, new[] { 2.0 } };
        var y = new[] { 1.0, 2.0 };
        Assert.Throws<ArgumentOutOfRangeException>(() => CrossValidation.KFold(x, y, (a, b) => new RidgeRegression().Fit(a, b), folds: 5));
    }
}

public class DescriptorMatrixTests
{
    private static IReadOnlyList<LibraryCompound> Compounds() =>
    [
        LibraryCompound.FromSmiles("aspirin", Drugs.Aspirin),
        LibraryCompound.FromSmiles("ibuprofen", Drugs.Ibuprofen),
        LibraryCompound.FromSmiles("caffeine", Drugs.Caffeine),
        LibraryCompound.FromSmiles("paracetamol", Drugs.Paracetamol),
        LibraryCompound.FromSmiles("naproxen", Drugs.Naproxen),
    ];

    [Fact]
    public void Standard_matrix_has_the_expected_shape_and_ids()
    {
        var m = DescriptorMatrixBuilder.Standard(Compounds());
        Assert.Equal(5, m.RowCount);
        Assert.Equal(DescriptorSet.Names.Count + Enum.GetValues<PharmacophoreFeature>().Length, m.ColumnCount);
        Assert.Equal(["aspirin", "ibuprofen", "caffeine", "paracetamol", "naproxen"], m.RowIds);
        Assert.Equal(MolecularDescriptors.MolecularWeight(Drugs.Parse(Drugs.Aspirin)), m.Rows[0][0], 6);
        Assert.Contains("Id,MW", m.ToCsv());
    }

    [Fact]
    public void Fingerprint_and_combined_matrices()
    {
        var fp = DescriptorMatrixBuilder.Fingerprints(Compounds(), FingerprintType.Ecfp4, 64);
        Assert.Equal(64, fp.ColumnCount);
        Assert.All(fp.Rows.SelectMany(r => r), v => Assert.True(v is 0 or 1));
        var combined = DescriptorMatrixBuilder.Combined(Compounds(), FingerprintType.Ecfp4, 64);
        Assert.Equal(DescriptorMatrixBuilder.Standard(Compounds()).ColumnCount + 64, combined.ColumnCount);
    }

    [Fact]
    public void Standardize_and_low_variance_removal()
    {
        var m = DescriptorMatrixBuilder.Custom(Compounds(),
            [("MW", MolecularDescriptors.MolecularWeight), ("constant", _ => 1.0)]);
        Assert.Equal(1, m.RemoveLowVariance().ColumnCount);
        var (scaled, mean, std) = m.Standardize();
        Assert.Equal(0, scaled.Column(0).Average(), 9);
        Assert.Equal(MolecularDescriptors.MolecularWeight(Drugs.Parse(Drugs.Aspirin)), mean[0] + scaled.Rows[0][0] * std[0], 6);
    }

    [Fact]
    public void Rows_must_match_the_column_names()
        => Assert.Throws<ArgumentException>(() => new DescriptorMatrix([[1.0, 2.0]], ["only-one"]));
}

public class QsarModelTests
{
    [Fact]
    public void Builds_cross_validates_and_predicts_a_descriptor_driven_response()
    {
        // a response that is a known function of descriptors: the model must recover it
        var smiles = new[]
        {
            "CCO", "CCCO", "CCCCO", "CCCCCO", "CCCCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1",
            Drugs.Aspirin, Drugs.Ibuprofen, Drugs.Paracetamol, Drugs.Naproxen, "CC(=O)O", "CCN", "CCCN", "c1ccncc1",
        };
        var compounds = smiles.Select((s, i) => LibraryCompound.FromSmiles($"c{i}", s)).ToList();
        var response = compounds.Select(c => 0.5 * Crippen.LogP(c.Molecule) + 0.01 * MolecularDescriptors.MolecularWeight(c.Molecule)).ToArray();
        var model = QsarModel.Build(compounds, response, folds: 4);
        Assert.True(model.CrossValidation.R2 > 0.8, model.CrossValidation.ToString());
        Assert.Equal(compounds.Count, model.Matrix.RowCount);
        // a held-out molecule inside the training domain is predicted closely; extrapolation is not claimed
        var heldOut = Drugs.Parse("CCCc1ccccc1");
        double expected = 0.5 * Crippen.LogP(heldOut) + 0.01 * MolecularDescriptors.MolecularWeight(heldOut);
        Assert.InRange(model.Predict(heldOut), expected - 0.5, expected + 0.5);
        Assert.True(model.Predict(Drugs.Parse("CCCCc1ccccc1")) > model.Predict(Drugs.Parse("CCO")));
    }

    [Fact]
    public void Works_with_a_random_forest_too()
    {
        var compounds = new[] { "CCO", "CCCO", "CCCCO", "CCCCCO", "c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1" }
            .Select((s, i) => LibraryCompound.FromSmiles($"c{i}", s)).ToList();
        var response = compounds.Select(c => MolecularDescriptors.MolecularWeight(c.Molecule)).ToArray();
        var model = QsarModel.Build(compounds, response, (x, y) => new RandomForestRegressor { TreeCount = 40, Seed = 2 }.Fit(x, y), folds: 4);
        Assert.True(model.Predict(Drugs.Parse("CCCCCCO")) > model.Predict(Drugs.Parse("CCO")));
    }
}
