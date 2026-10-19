using Biotech.Net.Transcriptomics.Internal;

namespace Biotech.Net.Transcriptomics;

/// <summary>Log fold change shrinkage estimator (DESeq2 <c>lfcShrink(type = ...)</c>).</summary>
public enum LfcShrinkType
{
    /// <summary>apeglm (Zhu, Ibrahim &amp; Love 2018): heavy-tailed Cauchy prior, MAP estimate, Laplace posterior SD.</summary>
    Apeglm,

    /// <summary>DESeq2's original normal prior (Love et al. 2014), the ridge fit of <c>betaPrior = TRUE</c>.</summary>
    Normal,
}

/// <summary>Shrunken log2 fold changes for one coefficient.</summary>
public sealed class LfcShrinkResult
{
    internal LfcShrinkResult(LfcShrinkType type, DifferentialExpressionResult table, double priorVariance, double priorScale,
        double[]? sValue, double[]? falseSignRate, int[]? convergence)
    {
        Type = type; Table = table; PriorVariance = priorVariance; PriorScale = priorScale;
        SValue = sValue; FalseSignRate = falseSignRate; Convergence = convergence;
    }

    /// <summary>The shrinkage estimator.</summary>
    public LfcShrinkType Type { get; }

    /// <summary>
    /// The result table: <c>Log2FoldChange</c> is the shrunken (MAP) estimate, <c>LfcSE</c> its posterior SD (apeglm) or
    /// standard error under the prior (normal); <c>PValue</c>/<c>PAdj</c> are those of the unshrunken test and
    /// <c>Stat</c> is NaN, as in DESeq2.
    /// </summary>
    public DifferentialExpressionResult Table { get; }

    /// <summary>Prior variance of the coefficient: apeglm's estimated A (natural-log scale) or the normal prior's
    /// <c>betaPriorVar</c> (log2 scale).</summary>
    public double PriorVariance { get; }

    /// <summary>apeglm's Cauchy prior scale min(√A, 1) (natural-log scale); NaN for the normal prior.</summary>
    public double PriorScale { get; }

    /// <summary>apeglm s-values (Stephens 2016): running mean of the sorted local false sign rates.</summary>
    public double[]? SValue { get; }

    /// <summary>apeglm local false sign rates P(sign error) = Φ(−|MAP| / SD).</summary>
    public double[]? FalseSignRate { get; }

    /// <summary>apeglm convergence codes (<c>fit$diag[, "conv"]</c>): 0 when the C++ prefit was accepted or the R-side BFGS
    /// refit converged, 1 when that refit hit its iteration limit.</summary>
    public int[]? Convergence { get; }
}

/// <summary>
/// DESeq2 <c>lfcShrink</c>: shrinkage of log2 fold changes towards zero for ranking and visualisation, ported from
/// DESeq2 1.52 / apeglm 1.34 (<c>method = "nbinomCR"</c>, the DESeq2 default).
/// </summary>
/// <remarks>
/// apeglm: the prior variance A of the true LFCs is the root of apeglm's <c>priorVar</c> moment equation (R's
/// <c>uniroot</c>), the Cauchy scale is min(√A, 1); per gene the MAP of the NB likelihood times the prior (Cauchy on the
/// coefficient, N(0, 15²) on the others) is found by the same L-BFGS (LBFGSpp) run from two starts, and the posterior SD
/// is the Laplace approximation from R's <c>optimHess</c>. The normal prior uses DESeq2's weighted upper-quantile prior
/// variance and a ridge-penalised refit.
/// </remarks>
public static class LfcShrink
{
    /// <summary>Shrinks the named coefficient (see <see cref="Shrink(DESeqFit, int, LfcShrinkType)"/>).</summary>
    public static LfcShrinkResult Shrink(DESeqFit fit, string coefficient, LfcShrinkType type = LfcShrinkType.Apeglm)
    {
        ArgumentNullException.ThrowIfNull(fit);
        return Shrink(fit, fit.Design.IndexOf(coefficient), type);
    }

    /// <summary>Shrinks the coefficient at <paramref name="coefficient"/> (default: the last one).</summary>
    public static LfcShrinkResult Shrink(DESeqFit fit, int coefficient = -1, LfcShrinkType type = LfcShrinkType.Apeglm)
        => type == LfcShrinkType.Apeglm ? Apeglm(fit, coefficient) : Normal(fit, coefficient);

    /// <summary>apeglm shrinkage of a named coefficient.</summary>
    public static LfcShrinkResult Apeglm(DESeqFit fit, string coefficient)
    {
        ArgumentNullException.ThrowIfNull(fit);
        return Apeglm(fit, fit.Design.IndexOf(coefficient));
    }

    /// <summary>Normal-prior shrinkage of a named coefficient.</summary>
    public static LfcShrinkResult Normal(DESeqFit fit, string coefficient)
    {
        ArgumentNullException.ThrowIfNull(fit);
        return Normal(fit, fit.Design.IndexOf(coefficient));
    }

    /// <summary>apeglm shrinkage (DESeq2 <c>lfcShrink(dds, coef, type = "apeglm")</c>).</summary>
    public static LfcShrinkResult Apeglm(DESeqFit fit, int coefficient = -1)
    {
        ArgumentNullException.ThrowIfNull(fit);
        var x = fit.Design;
        int p = x.Coefficients, m = x.Samples, n = fit.Counts.Rows;
        if (coefficient < 0) coefficient = p - 1;
        if (coefficient < 1 || coefficient >= p) throw new ArgumentOutOfRangeException(nameof(coefficient), "apeglm shrinks a non-intercept coefficient.");
        var res = fit.Results(coefficient);
        double ln2 = Math.Log(2);

        // Prior variance (apeglm::priorVar) from the MLEs on the natural-log scale.
        var X = new List<double>(); var Dv = new List<double>();
        for (int i = 0; i < n; i++)
            if (!double.IsNaN(res.Log2FoldChange[i])) { X.Add(ln2 * res.Log2FoldChange[i]); double s = ln2 * res.LfcSE[i]; Dv.Add(s * s); }
        double priorVar = PriorVar(X.ToArray(), Dv.ToArray());
        double S = Math.Min(Math.Sqrt(priorVar), 1);
        const double sigma = 15;

        var offset = fit.SizeFactors.Select(v => Math.Log(v)).ToArray();
        var refSamples = Enumerable.Range(0, m).Where(j => Enumerable.Range(0, p).Count(k => x[j, k] == 0) == p - 1).ToArray();
        var map = new double[n]; var sd = new double[n]; var fsr = new double[n]; var conv = new int[n];
        Array.Fill(map, double.NaN); Array.Fill(sd, double.NaN); Array.Fill(fsr, double.NaN); Array.Fill(conv, int.MinValue);
        Parallel.For(0, n, i =>
        {
            var yspan = fit.Counts.Row(i);
            double tot = 0;
            foreach (double v in yspan) tot += v;
            if (tot == 0) return;
            var y = yspan.ToArray();
            double size = 1 / fit.Dispersions[i];
            var prob = new NbProblem(x, y, size, offset, sigma, S, coefficient);
            // C++ prefit (nbinomGLM) from two starts.
            double cnst = 0;
            for (int j = 0; j < m; j++) cnst += (y[j] + size) * Math.Log(size + Math.Exp(offset[j]));
            cnst = cnst > 1 ? cnst : 1;
            var b1 = new double[p]; var b2 = new double[p];
            for (int k = 0; k < p; k++) { b1[k] = k % 2 == 0 ? 0.1 : -0.1; b2[k] = -b1[k]; }
            var solver = new Lbfgs(300, 1e-8, 1e-8);
            int status = solver.Minimize((bb, g) => prob.ScaledValueAndGradient(bb, g, cnst), b1, out _);
            solver.Minimize((bb, g) => prob.ScaledValueAndGradient(bb, g, cnst), b2, out _);
            double delta = 0;
            for (int k = 0; k < p; k++) delta = Math.Max(delta, Math.Abs(b1[k] - b2[k]));
            if (delta > .01) status = -1;
            conv[i] = 0; // apeglm.single overwrites diag[, "conv"] with the R-side optimiser code

            // optimNbinomHess (R side).
            double[] init;
            if (b1.Any(double.IsNaN))
            {
                init = new double[p];
                for (int k = 0; k < p; k++) init[k] = k % 2 == 0 ? 1 : -1;
                double bm = refSamples.Length > 0 ? refSamples.Average(j => y[j]) : y.Average();
                init[0] = bm == 0 ? 0 : Math.Log(bm);
            }
            else init = b1;
            double[] par = init;
            double[,]? hess = null;
            double[]? cov = null;
            if (status == 0)
            {
                hess = RMath.OptimHess(init, prob.Gradient);
                cov = InverseDiag(hess, p);
            }
            if (status != 0 || cov is null || cov.Any(v => v <= 0))
            {
                double rc = -prob.Value(init) - 1;
                par = (double[])init.Clone();
                conv[i] = RMath.Vmmin(par, b => prob.Value(b) + rc, prob.Gradient, out _);
                hess = RMath.OptimHess(par, prob.Gradient);
                cov = InverseDiag(hess, p);
            }
            map[i] = par[coefficient];
            if (cov is null || cov.Any(v => v <= 0)) return;
            sd[i] = Math.Sqrt(cov[coefficient]);
            fsr[i] = RMath.Pnorm(-Math.Abs(map[i]), 0, sd[i]);
        });
        var sval = SValue(fsr);
        const double log2e = 1.4426950408889634;
        var lfc = map.Select(v => log2e * v).ToArray();
        var se = sd.Select(v => log2e * v).ToArray();
        var table = new DifferentialExpressionResult("DESeq2-apeglm", res.Contrast, res.Genes, res.BaseMean, lfc, se,
            Enumerable.Repeat(double.NaN, n).ToArray(), res.PValue, res.PAdj);
        return new LfcShrinkResult(LfcShrinkType.Apeglm, table, priorVar, S, sval, fsr,
            conv.Select(c => c == int.MinValue ? 0 : c).ToArray());
    }

    /// <summary>DESeq2's original normal-prior shrinkage (<c>lfcShrink(dds, coef, type = "normal")</c>).</summary>
    public static LfcShrinkResult Normal(DESeqFit fit, int coefficient = -1)
    {
        ArgumentNullException.ThrowIfNull(fit);
        var x = fit.Design;
        int p = x.Coefficients, n = fit.Counts.Rows;
        if (coefficient < 0) coefficient = p - 1;
        if (coefficient < 0 || coefficient >= p) throw new ArgumentOutOfRangeException(nameof(coefficient));
        var res = fit.Results(coefficient);
        var rows = Enumerable.Range(0, n).Where(i => !fit.AllZero[i]).ToArray();
        // estimateBetaPriorVar(betaPriorMethod = "weighted", upperQuantile = 0.05)
        var priorVar = new double[p];
        for (int k = 0; k < p; k++)
        {
            if (x.CoefficientNames[k] == "Intercept") { priorVar[k] = 1e6; continue; }
            if (rows.Length == 1) { double b = fit.Coefficients[rows[0] * p + k]; priorVar[k] = b * b; continue; }
            var vals = new List<double>(); var wts = new List<double>();
            foreach (int i in rows)
            {
                double b = fit.Coefficients[i * p + k];
                if (!(Math.Abs(b) < 10)) continue;
                vals.Add(Math.Abs(b));
                wts.Add(1 / (1 / fit.BaseMean[i] + fit.TrendDispersions[i]));
            }
            if (vals.Count == 0) { priorVar[k] = 1e6; continue; }
            double sdEst = WeightedQuantile(vals.ToArray(), wts.ToArray(), 0.95) / RMath.Qnorm(1 - 0.05 / 2);
            priorVar[k] = sdEst * sdEst;
        }
        if (priorVar.Any(v => v == 0)) throw new InvalidOperationException("beta prior variances are equal to zero for some variables");
        var lambda = priorVar.Select(v => 1 / v).ToArray();
        var cm = new CountMatrix(fit.Counts.Values.ToArray(), n, fit.Counts.Columns);
        var alpha = rows.Select(i => fit.Dispersions[i]).ToArray();
        var glm = DESeqEngine.FitNbinomGlms(cm, rows, fit.SizeFactors, x, alpha, lambda, fit.MinMu, 1e-8, 100, true);
        var lfc = new double[n]; var se = new double[n];
        Array.Fill(lfc, double.NaN); Array.Fill(se, double.NaN);
        for (int r = 0; r < rows.Length; r++)
        {
            int i = rows[r];
            lfc[i] = glm.Beta[r * p + coefficient];
            se[i] = glm.Se[r * p + coefficient];
            if (fit.Replaced[i] && fit.BaseMean[i] == 0) { lfc[i] = 0; se[i] = 0; }
        }
        var table = new DifferentialExpressionResult("DESeq2-normal", res.Contrast, res.Genes, res.BaseMean, lfc, se,
            Enumerable.Repeat(double.NaN, n).ToArray(), res.PValue, res.PAdj);
        return new LfcShrinkResult(LfcShrinkType.Normal, table, priorVar[coefficient], double.NaN, null, null, null);
    }

    /// <summary>apeglm <c>priorVar</c>: root of Â(A) − A with Â(A) = Σ (X² − D) I(A) / Σ I(A), I(A) = 1/(2(A + D)²).</summary>
    internal static double PriorVar(double[] x, double[] d, double minVar = 1e-6, double maxVar = 400)
    {
        double Objective(double a)
        {
            double num = 0, den = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double w = 1 / (2 * (a + d[i]) * (a + d[i]));
                num += (x[i] * x[i] - d[i]) * w;
                den += w;
            }
            return num / den - a;
        }
        if (Objective(minVar) < 0) return minVar;
        return RMath.Uniroot(Objective, minVar, maxVar);
    }

    /// <summary>apeglm <c>svalue</c>: cumulative mean of the sorted local false sign rates (NaN last).</summary>
    internal static double[] SValue(double[] lfsr)
    {
        int n = lfsr.Length;
        var idx = Enumerable.Range(0, n).OrderBy(i => double.IsNaN(lfsr[i]) ? 1 : 0).ThenBy(i => double.IsNaN(lfsr[i]) ? 0 : lfsr[i]).ThenBy(i => i).ToArray();
        var res = new double[n];
        double cum = 0;
        for (int r = 0; r < n; r++)
        {
            cum += lfsr[idx[r]];
            res[idx[r]] = cum / (r + 1);
        }
        return res;
    }

    /// <summary>DESeq2's port of <c>Hmisc::wtd.quantile(x, weights, probs, normwt = TRUE)</c> (type "quantile").</summary>
    internal static double WeightedQuantile(double[] x, double[] w, double prob)
    {
        int n0 = x.Length;
        double sw = w.Sum();
        var order = Enumerable.Range(0, n0).OrderBy(i => x[i]).ToArray();
        var xs = new List<double>(); var ws = new List<double>();
        foreach (int i in order)
        {
            double wi = w[i] * n0 / sw;
            if (xs.Count > 0 && xs[^1] == x[i]) ws[^1] += wi;
            else { xs.Add(x[i]); ws.Add(wi); }
        }
        double n = ws.Sum();
        double ord = 1 + (n - 1) * prob;
        double low = Math.Max(Math.Floor(ord), 1), high = Math.Min(low + 1, n);
        ord %= 1;
        var cum = new double[ws.Count];
        double c = 0;
        for (int i = 0; i < ws.Count; i++) { c += ws[i]; cum[i] = c; }
        double ApproxConstant(double v)
        {
            if (v <= cum[0]) return xs[0];
            if (v >= cum[^1]) return xs[^1];
            int lo = 0, hi = cum.Length - 1;
            while (hi - lo > 1) { int mid = (lo + hi) / 2; if (cum[mid] <= v) lo = mid; else hi = mid; }
            return v == cum[lo] ? xs[lo] : xs[hi];
        }
        return (1 - ord) * ApproxConstant(low) + ord * ApproxConstant(high);
    }

    private static double[]? InverseDiag(double[,] h, int p)
    {
        var a = new double[p * p];
        for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) a[i * p + j] = h[i, j];
        var inv = DESeqEngine.Inverse(a, p);
        if (inv is null) return null;
        var d = new double[p];
        for (int i = 0; i < p; i++) d[i] = inv[i * p + i];
        return d;
    }

    /// <summary>apeglm's negative log posterior (nbinomFn / nbinomGr and the scaled C++ objective).</summary>
    private sealed class NbProblem(DesignMatrix x, double[] y, double size, double[] offset, double sigma, double s, int shrink)
    {
        private readonly int _m = x.Samples, _p = x.Coefficients;

        private double XBeta(double[] beta, int j)
        {
            double v = 0;
            for (int k = 0; k < _p; k++) v += x[j, k] * beta[k];
            return v;
        }

        private double NegPrior(double[] beta)
        {
            double v = 0;
            for (int k = 0; k < _p; k++)
                v += k == shrink ? RMath.Log1p(beta[k] * beta[k] / (s * s)) : beta[k] * beta[k] / (2 * sigma * sigma);
            return v;
        }

        /// <summary>nbinomFn without the constant.</summary>
        public double Value(double[] beta)
        {
            double ll = 0;
            for (int j = 0; j < _m; j++)
            {
                double xb = XBeta(beta, j);
                ll += y[j] * xb - (y[j] + size) * Math.Log(size + Math.Exp(xb + offset[j]));
            }
            return -ll + NegPrior(beta);
        }

        /// <summary>nbinomGr.</summary>
        public void Gradient(double[] beta, double[] g)
        {
            Array.Clear(g);
            for (int j = 0; j < _m; j++)
            {
                double e = Math.Exp(XBeta(beta, j) + offset[j]);
                double r = y[j] - (y[j] + size) * e / (size + e);
                for (int k = 0; k < _p; k++) g[k] -= x[j, k] * r;
            }
            for (int k = 0; k < _p; k++)
                g[k] += k == shrink ? 2 * beta[k] / (s * s + beta[k] * beta[k]) : beta[k] / (sigma * sigma);
        }

        /// <summary>The objective of apeglm's C++ <c>optimFun::f_grad</c>.</summary>
        public double ScaledValueAndGradient(double[] beta, double[] grad, double cnst)
        {
            double d = 0;
            var c = new double[_m];
            for (int j = 0; j < _m; j++)
            {
                double xb = XBeta(beta, j);
                double e = Math.Exp(xb + offset[j]);
                double a = y[j] + size, b = e + size;
                c[j] = y[j] - a * e / b;
                d += y[j] * xb - a * Math.Log(b);
            }
            double negPrior = 0;
            var dPrior = new double[_p];
            for (int k = 0; k < _p; k++)
            {
                if (k == shrink) { negPrior += RMath.Log1p(beta[k] * beta[k] / (s * s)); dPrior[k] = 2.0 * beta[k] / (s * s + beta[k] * beta[k]); }
                else { negPrior += beta[k] * beta[k] / (2.0 * sigma * sigma); dPrior[k] = beta[k] / (sigma * sigma); }
            }
            for (int k = 0; k < _p; k++)
            {
                double dl = 0;
                for (int j = 0; j < _m; j++) dl += x[j, k] * c[j];
                grad[k] = -dl / cnst + dPrior[k] / cnst;
            }
            return -d / cnst + negPrior / cnst + 10.0;
        }
    }
}
