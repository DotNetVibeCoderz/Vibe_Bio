using Biotech.Net.Chemistry;

namespace Biotech.Net.DrugDiscovery;

/// <summary>A feature matrix with named columns and optional row ids.</summary>
public sealed class DescriptorMatrix
{
    /// <summary>Creates a matrix from rows of equal length.</summary>
    public DescriptorMatrix(double[][] rows, IReadOnlyList<string> columnNames, IReadOnlyList<string>? rowIds = null)
    {
        ArgumentNullException.ThrowIfNull(rows);
        if (rows.Any(r => r.Length != columnNames.Count))
            throw new ArgumentException($"Every row must have {columnNames.Count} values.", nameof(rows));
        Rows = rows;
        ColumnNames = columnNames;
        RowIds = rowIds ?? Enumerable.Range(0, rows.Length).Select(i => i.ToString()).ToList();
    }

    /// <summary>The feature rows.</summary>
    public double[][] Rows { get; }

    /// <summary>Column (descriptor) names.</summary>
    public IReadOnlyList<string> ColumnNames { get; }

    /// <summary>Row identifiers.</summary>
    public IReadOnlyList<string> RowIds { get; }

    /// <summary>Number of rows.</summary>
    public int RowCount => Rows.Length;

    /// <summary>Number of columns.</summary>
    public int ColumnCount => ColumnNames.Count;

    /// <summary>Column <paramref name="j"/> as an array.</summary>
    public double[] Column(int j) => Rows.Select(r => r[j]).ToArray();

    /// <summary>Drops columns whose variance is below <paramref name="threshold"/>.</summary>
    public DescriptorMatrix RemoveLowVariance(double threshold = 1e-8)
    {
        var keep = Enumerable.Range(0, ColumnCount).Where(j =>
        {
            var c = Column(j);
            double mean = c.Average();
            return c.Sum(v => (v - mean) * (v - mean)) / Math.Max(1, c.Length - 1) > threshold;
        }).ToArray();
        return new DescriptorMatrix(Rows.Select(r => keep.Select(j => r[j]).ToArray()).ToArray(),
            keep.Select(j => ColumnNames[j]).ToArray(), RowIds);
    }

    /// <summary>Z-scores every column, returning the scaled matrix with the means and standard deviations used.</summary>
    public (DescriptorMatrix Scaled, double[] Mean, double[] Std) Standardize()
    {
        int p = ColumnCount;
        var mean = new double[p];
        var std = new double[p];
        for (int j = 0; j < p; j++)
        {
            var c = Column(j);
            mean[j] = c.Average();
            double v = c.Sum(x => (x - mean[j]) * (x - mean[j])) / Math.Max(1, c.Length - 1);
            std[j] = Math.Sqrt(v);
            if (std[j] < 1e-12) std[j] = 1;
        }
        var rows = Rows.Select(r => Enumerable.Range(0, p).Select(j => (r[j] - mean[j]) / std[j]).ToArray()).ToArray();
        return (new DescriptorMatrix(rows, ColumnNames, RowIds), mean, std);
    }

    /// <summary>Comma-separated text with a header row.</summary>
    public string ToCsv()
    {
        var sb = new System.Text.StringBuilder();
        sb.Append("Id,").AppendJoin(',', ColumnNames).Append('\n');
        for (int i = 0; i < RowCount; i++)
            sb.Append(RowIds[i]).Append(',').AppendJoin(',', Rows[i].Select(v => v.ToString("G6", System.Globalization.CultureInfo.InvariantCulture))).Append('\n');
        return sb.ToString();
    }
}

/// <summary>Builds QSAR feature matrices from molecules.</summary>
public static class DescriptorMatrixBuilder
{
    /// <summary>The standard 2-D descriptor block plus pharmacophore-feature counts.</summary>
    public static DescriptorMatrix Standard(IReadOnlyList<LibraryCompound> compounds)
    {
        var names = DescriptorSet.Names.Concat(Enum.GetValues<PharmacophoreFeature>().Select(f => "Ph" + f)).ToList();
        var rows = new double[compounds.Count][];
        Parallel.For(0, compounds.Count, i =>
            rows[i] = [.. MolecularDescriptors.Compute(compounds[i].Molecule).ToVector(), .. Pharmacophore.ToVector(compounds[i].Molecule)]);
        return new DescriptorMatrix(rows, names, compounds.Select(c => c.Id).ToList());
    }

    /// <summary>A folded fingerprint matrix (bits as 0/1 columns).</summary>
    public static DescriptorMatrix Fingerprints(IReadOnlyList<LibraryCompound> compounds, FingerprintType type = FingerprintType.Ecfp4, int nBits = 256)
    {
        var rows = new double[compounds.Count][];
        Parallel.For(0, compounds.Count, i =>
        {
            var fp = Fingerprint.Compute(compounds[i].Molecule, type).Fold(nBits);
            var r = new double[nBits];
            foreach (int b in fp.OnBits()) r[b] = 1;
            rows[i] = r;
        });
        return new DescriptorMatrix(rows, Enumerable.Range(0, nBits).Select(b => $"{type}_{b}").ToList(), compounds.Select(c => c.Id).ToList());
    }

    /// <summary>Descriptors and folded fingerprint bits side by side.</summary>
    public static DescriptorMatrix Combined(IReadOnlyList<LibraryCompound> compounds, FingerprintType type = FingerprintType.Ecfp4, int nBits = 256)
    {
        var a = Standard(compounds);
        var b = Fingerprints(compounds, type, nBits);
        var rows = Enumerable.Range(0, a.RowCount).Select(i => a.Rows[i].Concat(b.Rows[i]).ToArray()).ToArray();
        return new DescriptorMatrix(rows, a.ColumnNames.Concat(b.ColumnNames).ToList(), a.RowIds);
    }

    /// <summary>Custom descriptor functions.</summary>
    public static DescriptorMatrix Custom(IReadOnlyList<LibraryCompound> compounds, IReadOnlyList<(string Name, Func<Molecule, double> Value)> descriptors)
    {
        var rows = compounds.Select(c => descriptors.Select(d => d.Value(c.Molecule)).ToArray()).ToArray();
        return new DescriptorMatrix(rows, descriptors.Select(d => d.Name).ToList(), compounds.Select(c => c.Id).ToList());
    }
}

/// <summary>A fitted regression model.</summary>
public interface IRegressionModel
{
    /// <summary>Predicts the response for one feature row.</summary>
    double Predict(double[] features);

    /// <summary>Predicts the response for every row of <paramref name="x"/>.</summary>
    double[] Predict(double[][] x);
}

/// <summary>
/// Ridge (Tikhonov) linear regression: minimises ‖y − Xβ‖² + α‖β‖² with the intercept left unpenalised. Features
/// are centred and scaled internally, and the normal equations are solved by Cholesky decomposition.
/// </summary>
public sealed class RidgeRegression : IRegressionModel
{
    private double[] _coefficients = [];
    private double[] _mean = [], _std = [];

    /// <summary>Regularisation strength (0 gives ordinary least squares).</summary>
    public double Alpha { get; init; } = 1.0;

    /// <summary>The fitted intercept.</summary>
    public double Intercept { get; private set; }

    /// <summary>Coefficients on the original (unscaled) feature scale.</summary>
    public IReadOnlyList<double> Coefficients => _coefficients;

    /// <summary>Fits the model.</summary>
    /// <exception cref="ArgumentException">Shapes do not match or there are no rows.</exception>
    public RidgeRegression Fit(double[][] x, double[] y)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        if (x.Length != y.Length) throw new ArgumentException("x and y must have the same number of rows.");
        if (x.Length == 0) throw new ArgumentException("At least one row is required.", nameof(x));
        int n = x.Length, p = x[0].Length;
        _mean = new double[p];
        _std = new double[p];
        for (int j = 0; j < p; j++)
        {
            double m = 0;
            for (int i = 0; i < n; i++) m += x[i][j];
            m /= n;
            double v = 0;
            for (int i = 0; i < n; i++) v += (x[i][j] - m) * (x[i][j] - m);
            _mean[j] = m;
            _std[j] = Math.Sqrt(v / Math.Max(1, n - 1));
            if (_std[j] < 1e-12) _std[j] = 1;
        }
        double yMean = y.Average();
        // Gram matrix of the standardised features.
        var a = new double[p, p];
        var b = new double[p];
        var z = new double[n][];
        for (int i = 0; i < n; i++)
        {
            z[i] = new double[p];
            for (int j = 0; j < p; j++) z[i][j] = (x[i][j] - _mean[j]) / _std[j];
        }
        for (int j = 0; j < p; j++)
        {
            for (int k = j; k < p; k++)
            {
                double s = 0;
                for (int i = 0; i < n; i++) s += z[i][j] * z[i][k];
                a[j, k] = a[k, j] = s;
            }
            a[j, j] += Alpha;
            double sb = 0;
            for (int i = 0; i < n; i++) sb += z[i][j] * (y[i] - yMean);
            b[j] = sb;
        }
        var beta = SolveSymmetric(a, b);
        _coefficients = new double[p];
        double intercept = yMean;
        for (int j = 0; j < p; j++)
        {
            _coefficients[j] = beta[j] / _std[j];
            intercept -= _coefficients[j] * _mean[j];
        }
        Intercept = intercept;
        return this;
    }

    /// <inheritdoc/>
    public double Predict(double[] features)
    {
        if (_coefficients.Length == 0) throw new InvalidOperationException("The model has not been fitted.");
        double s = Intercept;
        for (int j = 0; j < _coefficients.Length; j++) s += _coefficients[j] * features[j];
        return s;
    }

    /// <inheritdoc/>
    public double[] Predict(double[][] x) => x.Select(Predict).ToArray();

    /// <summary>Solves A·x = b for a symmetric positive-definite A (Cholesky, with a ridge fallback).</summary>
    internal static double[] SolveSymmetric(double[,] a, double[] b)
    {
        int p = b.Length;
        var l = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                if (i == j)
                {
                    if (s <= 1e-12) s = 1e-12; // keep the factorisation defined for singular inputs
                    l[i, i] = Math.Sqrt(s);
                }
                else l[i, j] = s / l[j, j];
            }
        }
        var yv = new double[p];
        for (int i = 0; i < p; i++)
        {
            double s = b[i];
            for (int k = 0; k < i; k++) s -= l[i, k] * yv[k];
            yv[i] = s / l[i, i];
        }
        var x = new double[p];
        for (int i = p - 1; i >= 0; i--)
        {
            double s = yv[i];
            for (int k = i + 1; k < p; k++) s -= l[k, i] * x[k];
            x[i] = s / l[i, i];
        }
        return x;
    }
}

/// <summary>
/// A random-forest regressor: bootstrap-aggregated CART trees that split on a random subset of features
/// (<see cref="MaxFeatures"/>) by variance reduction, as in Breiman (Mach. Learn. 2001, 45, 5).
/// </summary>
public sealed class RandomForestRegressor : IRegressionModel
{
    private Node[] _trees = [];
    private int[][] _oobRows = [];

    /// <summary>Number of trees (default 100).</summary>
    public int TreeCount { get; init; } = 100;

    /// <summary>Maximum depth (default 12).</summary>
    public int MaxDepth { get; init; } = 12;

    /// <summary>Minimum rows required to split a node (default 4).</summary>
    public int MinSamplesSplit { get; init; } = 4;

    /// <summary>Features examined per split; 0 (default) means ⌈p / 3⌉, the usual regression setting.</summary>
    public int MaxFeatures { get; init; }

    /// <summary>Random seed.</summary>
    public int Seed { get; init; } = 42;

    /// <summary>Out-of-bag R² after fitting (NaN when no row was ever out of bag).</summary>
    public double OutOfBagR2 { get; private set; } = double.NaN;

    /// <summary>Fits the forest.</summary>
    public RandomForestRegressor Fit(double[][] x, double[] y)
    {
        ArgumentNullException.ThrowIfNull(x);
        if (x.Length != y.Length) throw new ArgumentException("x and y must have the same number of rows.");
        if (x.Length == 0) throw new ArgumentException("At least one row is required.", nameof(x));
        int n = x.Length, p = x[0].Length;
        int mtry = MaxFeatures > 0 ? Math.Min(MaxFeatures, p) : Math.Max(1, (int)Math.Ceiling(p / 3.0));
        _trees = new Node[TreeCount];
        _oobRows = new int[TreeCount][];
        Parallel.For(0, TreeCount, t =>
        {
            var rng = new Random(Seed + t * 7919);
            var bag = new int[n];
            var inBag = new bool[n];
            for (int i = 0; i < n; i++) { bag[i] = rng.Next(n); inBag[bag[i]] = true; }
            _oobRows[t] = Enumerable.Range(0, n).Where(i => !inBag[i]).ToArray();
            _trees[t] = Build(x, y, bag, 0, mtry, rng);
        });
        ComputeOob(x, y);
        return this;
    }

    private void ComputeOob(double[][] x, double[] y)
    {
        int n = x.Length;
        var sum = new double[n];
        var count = new int[n];
        for (int t = 0; t < _trees.Length; t++)
            foreach (int i in _oobRows[t]) { sum[i] += Evaluate(_trees[t], x[i]); count[i]++; }
        var rows = Enumerable.Range(0, n).Where(i => count[i] > 0).ToArray();
        if (rows.Length < 2) return;
        var pred = rows.Select(i => sum[i] / count[i]).ToArray();
        OutOfBagR2 = RegressionMetrics.R2(rows.Select(i => y[i]).ToArray(), pred);
    }

    /// <inheritdoc/>
    public double Predict(double[] features)
    {
        if (_trees.Length == 0) throw new InvalidOperationException("The model has not been fitted.");
        double s = 0;
        foreach (var t in _trees) s += Evaluate(t, features);
        return s / _trees.Length;
    }

    /// <inheritdoc/>
    public double[] Predict(double[][] x) => x.Select(Predict).ToArray();

    /// <summary>Mean decrease in variance per feature, normalised to sum to 1.</summary>
    public double[] FeatureImportances(int featureCount)
    {
        var imp = new double[featureCount];
        foreach (var t in _trees) Accumulate(t, imp);
        double total = imp.Sum();
        if (total > 0) for (int i = 0; i < imp.Length; i++) imp[i] /= total;
        return imp;
    }

    private static void Accumulate(Node n, double[] imp)
    {
        if (n.IsLeaf) return;
        imp[n.Feature] += n.ImpurityDecrease;
        Accumulate(n.Left!, imp);
        Accumulate(n.Right!, imp);
    }

    private static double Evaluate(Node n, double[] x)
    {
        while (!n.IsLeaf) n = x[n.Feature] <= n.Threshold ? n.Left! : n.Right!;
        return n.Value;
    }

    private Node Build(double[][] x, double[] y, int[] rows, int depth, int mtry, Random rng)
    {
        double mean = 0;
        foreach (int i in rows) mean += y[i];
        mean /= rows.Length;
        if (depth >= MaxDepth || rows.Length < MinSamplesSplit)
            return new Node { Value = mean };
        double sse = 0;
        foreach (int i in rows) sse += (y[i] - mean) * (y[i] - mean);
        if (sse < 1e-12) return new Node { Value = mean };

        int p = x[0].Length;
        var features = Enumerable.Range(0, p).ToArray();
        // partial Fisher–Yates for mtry features
        for (int k = 0; k < mtry; k++)
        {
            int j = k + rng.Next(p - k);
            (features[k], features[j]) = (features[j], features[k]);
        }
        int bestFeature = -1;
        double bestThreshold = 0, bestScore = sse;
        for (int k = 0; k < mtry; k++)
        {
            int f = features[k];
            var sorted = rows.OrderBy(i => x[i][f]).ToArray();
            double total = 0;
            foreach (int i in sorted) total += y[i];
            double leftSum = 0;
            int leftN = 0, rightN = sorted.Length;
            double totalSq = 0;
            foreach (int i in sorted) totalSq += y[i] * y[i];
            for (int t = 0; t < sorted.Length - 1; t++)
            {
                leftSum += y[sorted[t]];
                leftN++; rightN--;
                if (x[sorted[t]][f] >= x[sorted[t + 1]][f]) continue;
                // SSE after the split = Σy² − (Σy_L)²/n_L − (Σy_R)²/n_R
                double rightSum = total - leftSum;
                double residual = totalSq - leftSum * leftSum / leftN - rightSum * rightSum / rightN;
                if (residual < bestScore - 1e-12)
                {
                    bestScore = residual;
                    bestFeature = f;
                    bestThreshold = (x[sorted[t]][f] + x[sorted[t + 1]][f]) / 2;
                }
            }
        }
        if (bestFeature < 0) return new Node { Value = mean };
        var left = rows.Where(i => x[i][bestFeature] <= bestThreshold).ToArray();
        var right = rows.Where(i => x[i][bestFeature] > bestThreshold).ToArray();
        if (left.Length == 0 || right.Length == 0) return new Node { Value = mean };
        return new Node
        {
            Feature = bestFeature,
            Threshold = bestThreshold,
            ImpurityDecrease = (sse - bestScore) / x.Length,
            Left = Build(x, y, left, depth + 1, mtry, rng),
            Right = Build(x, y, right, depth + 1, mtry, rng),
            Value = mean,
        };
    }

    private sealed class Node
    {
        public int Feature = -1;
        public double Threshold;
        public double Value;
        public double ImpurityDecrease;
        public Node? Left, Right;
        public bool IsLeaf => Left is null;
    }
}

/// <summary>Regression quality measures.</summary>
public static class RegressionMetrics
{
    /// <summary>Coefficient of determination 1 − SS<sub>res</sub> / SS<sub>tot</sub>.</summary>
    public static double R2(double[] actual, double[] predicted)
    {
        if (actual.Length != predicted.Length) throw new ArgumentException("Arrays must have the same length.");
        double mean = actual.Average(), ssRes = 0, ssTot = 0;
        for (int i = 0; i < actual.Length; i++)
        {
            double d = actual[i] - predicted[i];
            ssRes += d * d;
            double t = actual[i] - mean;
            ssTot += t * t;
        }
        return ssTot < 1e-12 ? double.NaN : 1 - ssRes / ssTot;
    }

    /// <summary>Root mean squared error.</summary>
    public static double Rmse(double[] actual, double[] predicted)
    {
        double s = 0;
        for (int i = 0; i < actual.Length; i++) { double d = actual[i] - predicted[i]; s += d * d; }
        return Math.Sqrt(s / actual.Length);
    }

    /// <summary>Mean absolute error.</summary>
    public static double Mae(double[] actual, double[] predicted)
    {
        double s = 0;
        for (int i = 0; i < actual.Length; i++) s += Math.Abs(actual[i] - predicted[i]);
        return s / actual.Length;
    }

    /// <summary>Pearson correlation coefficient.</summary>
    public static double Pearson(double[] a, double[] b)
    {
        double ma = a.Average(), mb = b.Average(), sab = 0, sa = 0, sb = 0;
        for (int i = 0; i < a.Length; i++)
        {
            double da = a[i] - ma, db = b[i] - mb;
            sab += da * db; sa += da * da; sb += db * db;
        }
        return sa < 1e-12 || sb < 1e-12 ? double.NaN : sab / Math.Sqrt(sa * sb);
    }

    /// <summary>Spearman rank correlation.</summary>
    public static double Spearman(double[] a, double[] b)
    {
        static double[] Rank(double[] v)
        {
            var idx = Enumerable.Range(0, v.Length).OrderBy(i => v[i]).ToArray();
            var r = new double[v.Length];
            int k = 0;
            while (k < v.Length)
            {
                int j = k;
                while (j + 1 < v.Length && v[idx[j + 1]] == v[idx[k]]) j++;
                double avg = (k + j) / 2.0 + 1;
                for (int t = k; t <= j; t++) r[idx[t]] = avg;
                k = j + 1;
            }
            return r;
        }
        return Pearson(Rank(a), Rank(b));
    }
}

/// <summary>Cross-validation results: per-fold and pooled statistics.</summary>
/// <param name="R2">Pooled R² over the out-of-fold predictions (q² for leave-one-out).</param>
/// <param name="Rmse">Pooled RMSE.</param>
/// <param name="Mae">Pooled MAE.</param>
/// <param name="FoldR2">R² of each fold.</param>
/// <param name="Predictions">Out-of-fold prediction for every row, in input order.</param>
public sealed record CrossValidationResult(double R2, double Rmse, double Mae, IReadOnlyList<double> FoldR2, double[] Predictions)
{
    /// <inheritdoc/>
    public override string ToString() => $"R² {R2:F3}, RMSE {Rmse:F3}, MAE {Mae:F3} ({FoldR2.Count} folds)";
}

/// <summary>k-fold and leave-one-out cross-validation for any model factory.</summary>
public static class CrossValidation
{
    /// <summary>Runs <paramref name="folds"/>-fold cross-validation (shuffled with <paramref name="seed"/>).</summary>
    public static CrossValidationResult KFold(double[][] x, double[] y, Func<double[][], double[], IRegressionModel> fit, int folds = 5, int seed = 42)
    {
        ArgumentNullException.ThrowIfNull(x);
        if (x.Length != y.Length) throw new ArgumentException("x and y must have the same number of rows.");
        if (folds < 2 || folds > x.Length) throw new ArgumentOutOfRangeException(nameof(folds), $"folds must be between 2 and {x.Length}.");
        int n = x.Length;
        var order = Enumerable.Range(0, n).ToArray();
        var rng = new Random(seed);
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        var predictions = new double[n];
        var foldR2 = new List<double>();
        for (int f = 0; f < folds; f++)
        {
            var testRows = order.Where((_, i) => i % folds == f).ToArray();
            var trainRows = order.Where((_, i) => i % folds != f).ToArray();
            if (trainRows.Length == 0 || testRows.Length == 0) continue;
            var model = fit(trainRows.Select(i => x[i]).ToArray(), trainRows.Select(i => y[i]).ToArray());
            foreach (int i in testRows) predictions[i] = model.Predict(x[i]);
            if (testRows.Length > 1)
                foldR2.Add(RegressionMetrics.R2(testRows.Select(i => y[i]).ToArray(), testRows.Select(i => predictions[i]).ToArray()));
        }
        return new CrossValidationResult(RegressionMetrics.R2(y, predictions), RegressionMetrics.Rmse(y, predictions),
            RegressionMetrics.Mae(y, predictions), foldR2, predictions);
    }

    /// <summary>Leave-one-out cross-validation (the pooled R² is the classic q²).</summary>
    public static CrossValidationResult LeaveOneOut(double[][] x, double[] y, Func<double[][], double[], IRegressionModel> fit)
        => KFold(x, y, fit, x.Length, seed: 0);
}

/// <summary>
/// A convenience QSAR workflow: build descriptors, cross-validate a model and fit it on everything.
/// </summary>
public sealed class QsarModel
{
    private QsarModel(IRegressionModel model, DescriptorMatrix matrix, CrossValidationResult cv)
    {
        Model = model;
        Matrix = matrix;
        CrossValidation = cv;
    }

    /// <summary>The model fitted on all rows.</summary>
    public IRegressionModel Model { get; }

    /// <summary>The training feature matrix.</summary>
    public DescriptorMatrix Matrix { get; }

    /// <summary>Cross-validation statistics.</summary>
    public CrossValidationResult CrossValidation { get; }

    /// <summary>How the features were built, so predictions use the same descriptors.</summary>
    public required Func<Molecule, double[]> FeatureFunction { get; init; }

    /// <summary>Predicts the response for a molecule.</summary>
    public double Predict(Molecule mol) => Model.Predict(FeatureFunction(mol));

    /// <summary>
    /// Builds the standard descriptor matrix, cross-validates <paramref name="fit"/> and returns a model fitted
    /// on every compound.
    /// </summary>
    public static QsarModel Build(IReadOnlyList<LibraryCompound> compounds, double[] response,
        Func<double[][], double[], IRegressionModel>? fit = null, int folds = 5)
    {
        fit ??= (x, y) => new RidgeRegression { Alpha = 1.0 }.Fit(x, y);
        var matrix = DescriptorMatrixBuilder.Standard(compounds);
        var cv = Biotech.Net.DrugDiscovery.CrossValidation.KFold(matrix.Rows, response, fit, Math.Min(folds, compounds.Count));
        var model = fit(matrix.Rows, response);
        return new QsarModel(model, matrix, cv)
        {
            FeatureFunction = m => [.. MolecularDescriptors.Compute(m).ToVector(), .. Pharmacophore.ToVector(m)],
        };
    }
}
