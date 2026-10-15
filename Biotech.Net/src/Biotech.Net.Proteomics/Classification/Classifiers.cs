using Biotech.Net.Core;

namespace Biotech.Net.Proteomics.Classification;

/// <summary>A multi-class classifier over dense feature vectors with integer labels 0..C−1.</summary>
public interface IProteinClassifier
{
    /// <summary>Trains on rows of <paramref name="features"/> with labels <paramref name="labels"/>.</summary>
    void Fit(IReadOnlyList<double[]> features, IReadOnlyList<int> labels);

    /// <summary>Class probabilities (or vote fractions) for one sample.</summary>
    double[] PredictProbabilities(ReadOnlySpan<double> features);

    /// <summary>Most probable class.</summary>
    int Predict(ReadOnlySpan<double> features)
    {
        var p = PredictProbabilities(features);
        int best = 0;
        for (int i = 1; i < p.Length; i++) if (p[i] > p[best]) best = i;
        return best;
    }
}

/// <summary>Distance used by <see cref="KNearestNeighborsClassifier"/>.</summary>
public enum NeighborMetric
{
    /// <summary>Euclidean distance.</summary>
    Euclidean,

    /// <summary>1 − cosine similarity.</summary>
    Cosine,
}

/// <summary>k-nearest-neighbour classifier with majority voting (ties broken by summed inverse distance).</summary>
public sealed class KNearestNeighborsClassifier(int k = 5, NeighborMetric metric = NeighborMetric.Cosine) : IProteinClassifier
{
    private double[][] _x = [];
    private int[] _y = [];
    private int _classes;

    /// <summary>Number of neighbours.</summary>
    public int K { get; } = k >= 1 ? k : throw new ArgumentOutOfRangeException(nameof(k));

    /// <summary>Distance metric.</summary>
    public NeighborMetric Metric { get; } = metric;

    /// <inheritdoc/>
    public void Fit(IReadOnlyList<double[]> features, IReadOnlyList<int> labels)
    {
        Classifiers.Validate(features, labels);
        _x = features.ToArray();
        _y = labels.ToArray();
        _classes = _y.Max() + 1;
    }

    /// <inheritdoc/>
    public double[] PredictProbabilities(ReadOnlySpan<double> features)
    {
        if (_x.Length == 0) throw new InvalidOperationException("Call Fit first.");
        int k = Math.Min(K, _x.Length);
        var dist = new double[_x.Length];
        for (int i = 0; i < _x.Length; i++) dist[i] = Distance(features, _x[i]);
        var idx = new int[_x.Length];
        for (int i = 0; i < idx.Length; i++) idx[i] = i;
        Array.Sort(dist, idx);
        var votes = new double[_classes];
        var weight = new double[_classes];
        for (int i = 0; i < k; i++)
        {
            votes[_y[idx[i]]] += 1;
            weight[_y[idx[i]]] += 1 / (dist[i] + 1e-12);
        }
        // Break ties among equal vote counts by proximity, without changing the reported fractions much.
        double wsum = 0;
        foreach (var w in weight) wsum += w;
        for (int c = 0; c < _classes; c++) votes[c] = (votes[c] + 1e-6 * weight[c] / Math.Max(wsum, 1e-300)) / k;
        return votes;
    }

    private double Distance(ReadOnlySpan<double> a, double[] b)
    {
        if (a.Length != b.Length) throw new ArgumentException("Feature length differs from training data.");
        if (Metric == NeighborMetric.Euclidean)
        {
            double s = 0;
            for (int i = 0; i < a.Length; i++) { double d = a[i] - b[i]; s += d * d; }
            return Math.Sqrt(s);
        }
        double dot = 0, na = 0, nb = 0;
        for (int i = 0; i < a.Length; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
        return na == 0 || nb == 0 ? 1 : 1 - dot / Math.Sqrt(na * nb);
    }
}

/// <summary>
/// Multinomial (softmax) logistic regression with L2 regularisation, trained by full-batch gradient descent with
/// Nesterov momentum on z-scored features. The step is <see cref="LearningRate"/> / L, where L bounds the loss
/// curvature (½·λmax(XᵀX/n) + λ, λmax by power iteration), so training is stable for any feature count.
/// Deterministic: weights start at zero.
/// </summary>
public sealed class LogisticRegressionClassifier(double l2 = 1e-3, double learningRate = 1.0, int iterations = 300) : IProteinClassifier
{
    private double[] _w = [];   // classes × (features + 1), last column = bias
    private double[] _mean = [], _scale = [];
    private int _classes, _d;

    /// <summary>L2 penalty.</summary>
    public double L2 { get; } = l2;

    /// <summary>Step size as a fraction of 1/L (1 is the largest provably stable step).</summary>
    public double LearningRate { get; } = learningRate;

    /// <summary>Gradient steps.</summary>
    public int Iterations { get; } = iterations;

    /// <inheritdoc/>
    public void Fit(IReadOnlyList<double[]> features, IReadOnlyList<int> labels)
    {
        Classifiers.Validate(features, labels);
        int n = features.Count;
        _d = features[0].Length;
        _classes = labels.Max() + 1;
        _mean = new double[_d];
        _scale = new double[_d];
        foreach (var row in features) for (int j = 0; j < _d; j++) _mean[j] += row[j];
        for (int j = 0; j < _d; j++) _mean[j] /= n;
        foreach (var row in features) for (int j = 0; j < _d; j++) { double t = row[j] - _mean[j]; _scale[j] += t * t; }
        for (int j = 0; j < _d; j++) { double sd = Math.Sqrt(_scale[j] / n); _scale[j] = sd > 1e-12 ? 1 / sd : 0; }

        int stride = _d + 1;
        var x = new double[n * stride];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < _d; j++) x[i * stride + j] = (features[i][j] - _mean[j]) * _scale[j];
            x[i * stride + _d] = 1;
        }
        double step = LearningRate / (0.5 * LargestEigenvalue(x, n, stride) + L2);
        _w = new double[_classes * stride];
        var velocity = new double[_w.Length];
        var grad = new double[_w.Length];
        var look = new double[_w.Length];
        var p = new double[_classes];
        const double momentum = 0.9;
        for (int it = 0; it < Iterations; it++)
        {
            for (int i = 0; i < look.Length; i++) look[i] = _w[i] + momentum * velocity[i];
            Array.Clear(grad);
            for (int i = 0; i < n; i++)
            {
                var xi = x.AsSpan(i * stride, stride);
                Softmax(look, xi, p);
                p[labels[i]] -= 1;
                for (int c = 0; c < _classes; c++)
                {
                    double g = p[c];
                    if (g == 0) continue;
                    var gc = grad.AsSpan(c * stride, stride);
                    for (int j = 0; j < stride; j++) gc[j] += g * xi[j];
                }
            }
            for (int i = 0; i < _w.Length; i++)
            {
                double g = grad[i] / n + (i % stride == _d ? 0 : L2 * look[i]);
                velocity[i] = momentum * velocity[i] - step * g;
                _w[i] += velocity[i];
            }
        }
    }

    // λmax of XᵀX/n by power iteration (30 steps, deterministic start).
    private static double LargestEigenvalue(double[] x, int n, int stride)
    {
        var v = new double[stride];
        Array.Fill(v, 1 / Math.Sqrt(stride));
        var u = new double[stride];
        double lambda = 1;
        for (int it = 0; it < 30; it++)
        {
            Array.Clear(u);
            for (int i = 0; i < n; i++)
            {
                var xi = x.AsSpan(i * stride, stride);
                double dot = 0;
                for (int j = 0; j < stride; j++) dot += xi[j] * v[j];
                for (int j = 0; j < stride; j++) u[j] += dot * xi[j];
            }
            double norm = 0;
            for (int j = 0; j < stride; j++) { u[j] /= n; norm += u[j] * u[j]; }
            norm = Math.Sqrt(norm);
            if (norm == 0) return 1;
            lambda = norm;
            for (int j = 0; j < stride; j++) v[j] = u[j] / norm;
        }
        return lambda * 1.05; // small safety margin
    }

    /// <inheritdoc/>
    public double[] PredictProbabilities(ReadOnlySpan<double> features)
    {
        if (_w.Length == 0) throw new InvalidOperationException("Call Fit first.");
        if (features.Length != _d) throw new ArgumentException("Feature length differs from training data.");
        var xi = new double[_d + 1];
        for (int j = 0; j < _d; j++) xi[j] = (features[j] - _mean[j]) * _scale[j];
        xi[_d] = 1;
        var p = new double[_classes];
        Softmax(_w, xi, p);
        return p;
    }

    private void Softmax(double[] w, ReadOnlySpan<double> x, double[] p)
    {
        int stride = x.Length;
        double max = double.NegativeInfinity;
        for (int c = 0; c < _classes; c++)
        {
            var wc = w.AsSpan(c * stride, stride);
            double z = 0;
            for (int j = 0; j < stride; j++) z += wc[j] * x[j];
            p[c] = z;
            if (z > max) max = z;
        }
        double s = 0;
        for (int c = 0; c < _classes; c++) { p[c] = Math.Exp(p[c] - max); s += p[c]; }
        for (int c = 0; c < _classes; c++) p[c] /= s;
    }
}

/// <summary>Cross-validation outcome.</summary>
/// <param name="FoldAccuracies">Accuracy of each fold.</param>
/// <param name="Accuracy">Overall accuracy over all held-out predictions.</param>
/// <param name="ConfusionMatrix">Row-major C×C counts, rows = true class, columns = predicted.</param>
/// <param name="Classes">Number of classes.</param>
public sealed record CrossValidationResult(double[] FoldAccuracies, double Accuracy, int[] ConfusionMatrix, int Classes)
{
    /// <summary>Macro-averaged F1 score.</summary>
    public double MacroF1
    {
        get
        {
            double sum = 0;
            for (int c = 0; c < Classes; c++)
            {
                double tp = ConfusionMatrix[c * Classes + c], fp = 0, fn = 0;
                for (int o = 0; o < Classes; o++)
                {
                    if (o == c) continue;
                    fp += ConfusionMatrix[o * Classes + c];
                    fn += ConfusionMatrix[c * Classes + o];
                }
                sum += tp == 0 ? 0 : 2 * tp / (2 * tp + fp + fn);
            }
            return sum / Classes;
        }
    }
}

/// <summary>Validation helpers and stratified k-fold cross-validation.</summary>
public static class Classifiers
{
    internal static void Validate(IReadOnlyList<double[]> features, IReadOnlyList<int> labels)
    {
        ArgumentNullException.ThrowIfNull(features);
        ArgumentNullException.ThrowIfNull(labels);
        if (features.Count == 0) throw new ArgumentException("No training samples.");
        if (features.Count != labels.Count) throw new ArgumentException("Features and labels differ in count.");
        int d = features[0].Length;
        foreach (var f in features) if (f.Length != d) throw new ArgumentException("All feature vectors need the same length.");
        foreach (int y in labels) if (y < 0) throw new ArgumentException("Labels must be ≥ 0.");
    }

    /// <summary>
    /// Stratified fold assignment: samples of each class are shuffled (seeded) and dealt round-robin, so every fold keeps
    /// the class proportions. Returns the fold (0..folds−1) of each sample.
    /// </summary>
    public static int[] StratifiedFolds(IReadOnlyList<int> labels, int folds, int seed = 42)
    {
        ArgumentNullException.ThrowIfNull(labels);
        ArgumentOutOfRangeException.ThrowIfLessThan(folds, 2);
        var rng = new Random(seed);
        var assignment = new int[labels.Count];
        int offset = 0;
        foreach (var group in Enumerable.Range(0, labels.Count).GroupBy(i => labels[i]).OrderBy(g => g.Key))
        {
            var idx = group.ToArray();
            rng.Shuffle(idx);
            for (int i = 0; i < idx.Length; i++) assignment[idx[i]] = (offset + i) % folds;
            offset += idx.Length;
        }
        return assignment;
    }

    /// <summary>Stratified k-fold cross-validation of a classifier produced by <paramref name="factory"/>.</summary>
    public static CrossValidationResult CrossValidate(Func<IProteinClassifier> factory, IReadOnlyList<double[]> features, IReadOnlyList<int> labels, int folds = 5, int seed = 42)
    {
        ArgumentNullException.ThrowIfNull(factory);
        Validate(features, labels);
        int classes = labels.Max() + 1;
        var fold = StratifiedFolds(labels, folds, seed);
        var accuracies = new double[folds];
        var confusion = new int[classes * classes];
        int correct = 0, total = 0;
        for (int f = 0; f < folds; f++)
        {
            var trainX = new List<double[]>();
            var trainY = new List<int>();
            for (int i = 0; i < features.Count; i++) if (fold[i] != f) { trainX.Add(features[i]); trainY.Add(labels[i]); }
            var model = factory();
            model.Fit(trainX, trainY);
            int fc = 0, ft = 0;
            for (int i = 0; i < features.Count; i++)
            {
                if (fold[i] != f) continue;
                int pred = model.Predict(features[i]);
                if (pred < classes) confusion[labels[i] * classes + pred]++;
                if (pred == labels[i]) fc++;
                ft++;
            }
            accuracies[f] = ft == 0 ? double.NaN : (double)fc / ft;
            correct += fc;
            total += ft;
        }
        return new(accuracies, (double)correct / total, confusion, classes);
    }
}

/// <summary>
/// Protein family classifier: a feature extractor plus any <see cref="IProteinClassifier"/>, working on string labels
/// and Core <see cref="Protein"/> objects.
/// </summary>
public sealed class ProteinFamilyClassifier
{
    private readonly Func<byte[], double[]> _extractor;
    private readonly IProteinClassifier _model;
    private string[] _labels = [];

    /// <summary>Creates a classifier (defaults: amino-acid + dipeptide composition with logistic regression).</summary>
    public ProteinFamilyClassifier(IProteinClassifier? model = null, Func<byte[], double[]>? featureExtractor = null)
    {
        _model = model ?? new LogisticRegressionClassifier();
        _extractor = featureExtractor ?? (b =>
        {
            var aac = ProteinFeatures.AminoAcidComposition(b);
            var dpc = ProteinFeatures.DipeptideComposition(b);
            return [.. aac, .. dpc];
        });
    }

    /// <summary>The class names, in label-index order.</summary>
    public IReadOnlyList<string> Classes => _labels;

    /// <summary>Feature vector of one protein.</summary>
    public double[] Features(Protein protein) => _extractor(protein.Sequence.Span.ToArray());

    /// <summary>Trains on labelled proteins.</summary>
    public ProteinFamilyClassifier Fit(IEnumerable<(Protein Protein, string Family)> examples)
    {
        ArgumentNullException.ThrowIfNull(examples);
        var list = examples.ToList();
        _labels = list.Select(e => e.Family).Distinct().Order(StringComparer.Ordinal).ToArray();
        var map = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < _labels.Length; i++) map[_labels[i]] = i;
        var x = new double[list.Count][];
        var y = new int[list.Count];
        Parallel.For(0, list.Count, i => x[i] = Features(list[i].Protein));
        for (int i = 0; i < list.Count; i++) y[i] = map[list[i].Family];
        _model.Fit(x, y);
        return this;
    }

    /// <summary>Predicted family.</summary>
    public string Predict(Protein protein) => _labels[_model.Predict(Features(protein))];

    /// <summary>Family → probability.</summary>
    public IReadOnlyDictionary<string, double> PredictProbabilities(Protein protein)
    {
        var p = _model.PredictProbabilities(Features(protein));
        var d = new Dictionary<string, double>(StringComparer.Ordinal);
        for (int i = 0; i < _labels.Length && i < p.Length; i++) d[_labels[i]] = p[i];
        return d;
    }

    /// <summary>Stratified k-fold cross-validation of this feature/model combination on labelled proteins.</summary>
    public static CrossValidationResult CrossValidate(IEnumerable<(Protein Protein, string Family)> examples, Func<IProteinClassifier> modelFactory,
        Func<byte[], double[]>? featureExtractor = null, int folds = 5, int seed = 42)
    {
        var list = examples.ToList();
        var probe = new ProteinFamilyClassifier(modelFactory(), featureExtractor);
        var labels = list.Select(e => e.Family).Distinct().Order(StringComparer.Ordinal).ToList();
        var x = new double[list.Count][];
        Parallel.For(0, list.Count, i => x[i] = probe.Features(list[i].Protein));
        var y = list.Select(e => labels.IndexOf(e.Family)).ToArray();
        return Classifiers.CrossValidate(modelFactory, x, y, folds, seed);
    }
}
