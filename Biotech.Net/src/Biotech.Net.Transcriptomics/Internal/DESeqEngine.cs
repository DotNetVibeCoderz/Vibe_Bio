namespace Biotech.Net.Transcriptomics.Internal;

/// <summary>Row-major genes × samples counts (a mutable copy, so outlier replacement can edit it).</summary>
internal sealed class CountMatrix(double[] values, int rows, int columns)
{
    public double[] Values { get; } = values;
    public int Rows { get; } = rows;
    public int Columns { get; } = columns;
    public ReadOnlySpan<double> Row(int i) => Values.AsSpan(i * Columns, Columns);
    public double this[int i, int j] { get => Values[i * Columns + j]; set => Values[i * Columns + j] = value; }
    public CountMatrix Clone() => new((double[])Values.Clone(), Rows, Columns);
}

/// <summary>Result of <see cref="DESeqEngine.FitNbinomGlms"/> for a set of rows (DESeq2's <c>fitNbinomGLMs</c>).</summary>
internal sealed class GlmFitResult(int n, int m, int p)
{
    /// <summary>log2 coefficients, rows × p.</summary>
    public double[] Beta { get; } = new double[n * p];
    /// <summary>log2 standard errors, rows × p.</summary>
    public double[] Se { get; } = new double[n * p];
    /// <summary>Fitted means (not clamped at minmu), rows × m.</summary>
    public double[] Mu { get; } = new double[n * m];
    /// <summary>Hat-matrix diagonals from the IRLS fit, rows × m.</summary>
    public double[] Hat { get; } = new double[n * m];
    public bool[] Converged { get; } = new bool[n];
    public int[] Iterations { get; } = new int[n];
    public double[] LogLike { get; } = new double[n];
    public int OptimRows { get; set; }
}

/// <summary>
/// Faithful ports of DESeq2's numerical core (R/core.R, R/fitNbinomGLMs.R and src/DESeq2.cpp of DESeq2 1.52): the
/// Armijo line search on the Cox–Reid adjusted log posterior of log α (<c>fitDisp</c>), its grid fallback
/// (<c>fitDispGrid</c>), the ridge-penalised IRLS negative-binomial GLM (<c>fitBeta</c>), the parametric trend
/// (<c>parametricDispersionFit</c> = <c>glm(family = Gamma("identity"))</c>), robust method-of-moments dispersions
/// and Cook's distances.
/// </summary>
internal static class DESeqEngine
{
    public const double Large = 30.0;

    // ---- small dense linear algebra (row-major p × p) ----------------------------------------------------------

    /// <summary>Determinant by LU with partial pivoting.</summary>
    public static double Det(double[] a, int p)
    {
        if (p == 1) return a[0];
        if (p == 2) return a[0] * a[3] - a[1] * a[2];
        var lu = (double[])a.Clone();
        double det = 1;
        for (int c = 0; c < p; c++)
        {
            int piv = c;
            for (int r = c + 1; r < p; r++) if (Math.Abs(lu[r * p + c]) > Math.Abs(lu[piv * p + c])) piv = r;
            if (lu[piv * p + c] == 0) return 0;
            if (piv != c) { for (int k = 0; k < p; k++) (lu[c * p + k], lu[piv * p + k]) = (lu[piv * p + k], lu[c * p + k]); det = -det; }
            det *= lu[c * p + c];
            for (int r = c + 1; r < p; r++)
            {
                double f = lu[r * p + c] / lu[c * p + c];
                for (int k = c; k < p; k++) lu[r * p + k] -= f * lu[c * p + k];
            }
        }
        return det;
    }

    /// <summary>Inverse by Gauss–Jordan elimination with partial pivoting (null when singular).</summary>
    public static double[]? Inverse(double[] a, int p)
    {
        var m = (double[])a.Clone();
        var inv = new double[p * p];
        for (int i = 0; i < p; i++) inv[i * p + i] = 1;
        for (int c = 0; c < p; c++)
        {
            int piv = c;
            for (int r = c + 1; r < p; r++) if (Math.Abs(m[r * p + c]) > Math.Abs(m[piv * p + c])) piv = r;
            if (m[piv * p + c] == 0 || !double.IsFinite(m[piv * p + c])) return null;
            if (piv != c)
                for (int k = 0; k < p; k++)
                {
                    (m[c * p + k], m[piv * p + k]) = (m[piv * p + k], m[c * p + k]);
                    (inv[c * p + k], inv[piv * p + k]) = (inv[piv * p + k], inv[c * p + k]);
                }
            double d = m[c * p + c];
            for (int k = 0; k < p; k++) { m[c * p + k] /= d; inv[c * p + k] /= d; }
            for (int r = 0; r < p; r++)
            {
                if (r == c) continue;
                double f = m[r * p + c];
                if (f == 0) continue;
                for (int k = 0; k < p; k++) { m[r * p + k] -= f * m[c * p + k]; inv[r * p + k] -= f * inv[c * p + k]; }
            }
        }
        return inv;
    }

    /// <summary>Solves A x = b for small p via the inverse (null when singular).</summary>
    public static double[]? Solve(double[] a, double[] b, int p)
    {
        var inv = Inverse(a, p);
        if (inv is null) return null;
        var x = new double[p];
        for (int i = 0; i < p; i++) { double s = 0; for (int k = 0; k < p; k++) s += inv[i * p + k] * b[k]; x[i] = s; }
        return x;
    }

    /// <summary>Least-squares fitted values H y for the design (DESeq2 <c>linearModelMu</c>).</summary>
    public static double[] LinearModelMu(ReadOnlySpan<double> y, DesignMatrix x, double[] gramInverse)
    {
        int m = x.Samples, p = x.Coefficients;
        var xty = new double[p];
        for (int j = 0; j < m; j++) for (int k = 0; k < p; k++) xty[k] += x[j, k] * y[j];
        var b = new double[p];
        for (int a = 0; a < p; a++) { double s = 0; for (int k = 0; k < p; k++) s += gramInverse[a * p + k] * xty[k]; b[a] = s; }
        var mu = new double[m];
        for (int j = 0; j < m; j++) { double s = 0; for (int k = 0; k < p; k++) s += x[j, k] * b[k]; mu[j] = s; }
        return mu;
    }

    /// <summary>Number of distinct rows of the design (DESeq2 <c>nlevels(modelMatrixGroups(x))</c>).</summary>
    public static int[] RowGroups(DesignMatrix x, out int groups)
    {
        int m = x.Samples, p = x.Coefficients;
        var id = new int[m];
        groups = 0;
        for (int j = 0; j < m; j++)
        {
            id[j] = -1;
            for (int k = 0; k < j; k++)
            {
                bool same = true;
                for (int c = 0; c < p && same; c++) same = x[j, c] == x[k, c];
                if (same) { id[j] = id[k]; break; }
            }
            if (id[j] < 0) id[j] = groups++;
        }
        return id;
    }

    /// <summary>DESeq2 <c>nOrMoreInCell</c>: samples whose design row occurs at least <paramref name="n"/> times.</summary>
    public static bool[] NOrMoreInCell(DesignMatrix x, int n)
    {
        var id = RowGroups(x, out int g);
        var count = new int[g];
        foreach (int i in id) count[i]++;
        var res = new bool[id.Length];
        for (int j = 0; j < id.Length; j++) res[j] = count[id[j]] >= n;
        return res;
    }

    // ---- dispersion likelihood (src/DESeq2.cpp) -----------------------------------------------------------------

    private static void CrMatrices(ReadOnlySpan<double> mu, DesignMatrix x, double alpha, double[] b, double[]? db)
    {
        int m = x.Samples, p = x.Coefficients;
        Array.Clear(b);
        if (db is not null) Array.Clear(db);
        for (int j = 0; j < m; j++)
        {
            double s = 1 / mu[j] + alpha;
            double w = 1 / s;
            double dw = -1 / (s * s);
            for (int a = 0; a < p; a++)
            {
                double xa = x[j, a];
                if (xa == 0) continue;
                for (int c = 0; c < p; c++)
                {
                    double xc = x[j, c];
                    b[a * p + c] += xa * w * xc;
                    if (db is not null) db[a * p + c] += xa * dw * xc;
                }
            }
        }
    }

    /// <summary>DESeq2 <c>log_posterior</c> (Cox–Reid adjusted, optional normal prior on log α).</summary>
    public static double LogPosterior(double logAlpha, ReadOnlySpan<double> y, ReadOnlySpan<double> mu, DesignMatrix x,
        double priorMean, double priorVar, bool usePrior, double[] b)
    {
        double alpha = Math.Exp(logAlpha);
        int p = x.Coefficients;
        CrMatrices(mu, x, alpha, b, null);
        double crTerm = -0.5 * Math.Log(Det(b, p));
        double r = 1 / alpha;
        double lgr = RMath.LGamma(r);
        double ll = 0;
        for (int j = 0; j < y.Length; j++)
            ll += RMath.LGamma(y[j] + r) - lgr - y[j] * Math.Log(mu[j] + r) - r * Math.Log(1.0 + mu[j] * alpha);
        double prior = usePrior ? -0.5 * (logAlpha - priorMean) * (logAlpha - priorMean) / priorVar : 0.0;
        return ll + prior + crTerm;
    }

    /// <summary>DESeq2 <c>dlog_posterior</c>: derivative with respect to log α.</summary>
    public static double DLogPosterior(double logAlpha, ReadOnlySpan<double> y, ReadOnlySpan<double> mu, DesignMatrix x,
        double priorMean, double priorVar, bool usePrior, double[] b, double[] db)
    {
        double alpha = Math.Exp(logAlpha);
        int p = x.Coefficients;
        CrMatrices(mu, x, alpha, b, db);
        var binv = Inverse(b, p);
        double tr = 0;
        if (binv is null) tr = double.NaN;
        else
            for (int a = 0; a < p; a++)
                for (int c = 0; c < p; c++) tr += binv[a * p + c] * db[c * p + a];
        double crTerm = -0.5 * tr;
        double r = 1 / alpha;
        double dgr = RMath.Digamma(r);
        double s = 0;
        for (int j = 0; j < y.Length; j++)
        {
            double ma = mu[j] * alpha;
            s += dgr + Math.Log(1 + ma) - ma / (1.0 + ma) - RMath.Digamma(y[j] + r) + y[j] / (mu[j] + r);
        }
        double ll = s / (alpha * alpha);
        double prior = usePrior ? -1.0 * (logAlpha - priorMean) / priorVar : 0.0;
        return (ll + crTerm) * alpha + prior;
    }

    /// <summary>Outcome of one <see cref="FitDisp"/> call.</summary>
    public readonly record struct DispFitResult(double LogAlpha, int Iterations, double InitialLp, double LastLp);

    /// <summary>DESeq2 <c>fitDisp</c>: Armijo backtracking line search for the maximum of the log posterior.</summary>
    public static DispFitResult FitDisp(ReadOnlySpan<double> y, ReadOnlySpan<double> mu, DesignMatrix x, double logAlpha,
        double priorMean, double priorVar, bool usePrior, double minLogAlpha, double kappa0, double tol, int maxit,
        double[] b, double[] db)
    {
        const double epsilon = 1.0e-4;
        double a = logAlpha;
        double lp = LogPosterior(a, y, mu, x, priorMean, priorVar, usePrior, b);
        double dlp = DLogPosterior(a, y, mu, x, priorMean, priorVar, usePrior, b, db);
        double kappa = kappa0;
        double initialLp = lp;
        int iter = 0, iterAccept = 0;
        for (int t = 0; t < maxit; t++)
        {
            iter++;
            double aPropose = a + kappa * dlp;
            if (aPropose < -30.0) kappa = (-30.0 - a) / dlp;
            if (aPropose > 10.0) kappa = (10.0 - a) / dlp;
            double thetaKappa = -1.0 * LogPosterior(a + kappa * dlp, y, mu, x, priorMean, priorVar, usePrior, b);
            double thetaHatKappa = -1.0 * lp - kappa * epsilon * dlp * dlp;
            if (thetaKappa <= thetaHatKappa)
            {
                iterAccept++;
                a += kappa * dlp;
                double lpnew = LogPosterior(a, y, mu, x, priorMean, priorVar, usePrior, b);
                double change = lpnew - lp;
                if (change < tol) { lp = lpnew; break; }
                if (a < minLogAlpha) break;
                lp = lpnew;
                dlp = DLogPosterior(a, y, mu, x, priorMean, priorVar, usePrior, b, db);
                kappa = Math.Min(kappa * 1.1, kappa0);
                if (iterAccept % 5 == 0) kappa /= 2.0;
            }
            else kappa /= 2.0;
        }
        return new DispFitResult(a, iter, initialLp, lp);
    }

    /// <summary>DESeq2 <c>fitDispGridWrapper</c>/<c>fitDispGrid</c>: coarse 20-point grid on log α, then a fine grid.</summary>
    public static double FitDispGrid(ReadOnlySpan<double> y, ReadOnlySpan<double> mu, DesignMatrix x, double priorMean,
        double priorVar, bool usePrior, double[] b)
    {
        const int n = 20;
        double minLog = Math.Log(1e-8), maxLog = Math.Log(Math.Max(10, x.Samples));
        var grid = RSeq(minLog, maxLog, n);
        double delta = grid[1] - grid[0];
        int best = 0; double bestV = double.NegativeInfinity;
        for (int t = 0; t < n; t++)
        {
            double v = LogPosterior(grid[t], y, mu, x, priorMean, priorVar, usePrior, b);
            if (t == 0 || v > bestV) { bestV = v; best = t; }
        }
        double aHat = grid[best];
        var fine = ArmaLinspace(aHat - delta, aHat + delta, n);
        best = 0; bestV = double.NegativeInfinity;
        for (int t = 0; t < n; t++)
        {
            double v = LogPosterior(fine[t], y, mu, x, priorMean, priorVar, usePrior, b);
            if (t == 0 || v > bestV) { bestV = v; best = t; }
        }
        return Math.Exp(fine[best]);
    }

    /// <summary>R's <c>seq(from, to, length.out = n)</c>.</summary>
    public static double[] RSeq(double from, double to, int n)
    {
        var r = new double[n];
        double by = (to - from) / (n - 1);
        for (int i = 0; i < n; i++) r[i] = from + i * by;
        return r;
    }

    private static double[] ArmaLinspace(double start, double end, int n)
    {
        var r = new double[n];
        double delta = (end - start) / (n - 1);
        for (int i = 0; i < n; i++) r[i] = start + i * delta;
        r[n - 1] = end;
        return r;
    }

    // ---- negative-binomial GLM (fitNbinomGLMs / fitBeta / fitNbinomGLMsOptim) ----------------------------------

    /// <summary>
    /// DESeq2 <c>fitNbinomGLMs</c> for the given rows: QR least-squares start on log(normalised + 0.1), ridge-penalised
    /// IRLS (<c>fitBeta</c>), log-likelihood, and the bounded quasi-Newton fallback for rows that did not converge.
    /// <paramref name="lambdaLog2"/> are the ridge penalties on the log2 scale (DESeq2 default 1e-6).
    /// </summary>
    public static GlmFitResult FitNbinomGlms(CountMatrix counts, int[] rows, double[] sf, DesignMatrix x, double[] alpha,
        double[] lambdaLog2, double minmu, double betaTol, int maxit, bool useOptim)
    {
        int n = rows.Length, m = x.Samples, p = x.Coefficients;
        var res = new GlmFitResult(n, m, p);
        const double log2e = 1.4426950408889634; // log2(exp(1))
        double ln2 = Math.Log(2);
        bool justIntercept = p == 1;
        for (int j = 0; j < m && justIntercept; j++) justIntercept = x[j, 0] == 1;
        bool smallLambda = true;
        foreach (double l in lambdaLog2) if (l > 1e-6) smallLambda = false;
        if (justIntercept && smallLambda)
        {
            for (int r = 0; r < n; r++)
            {
                var y = counts.Row(rows[r]);
                double s = 0;
                for (int j = 0; j < m; j++) s += y[j] / sf[j];
                double beta = Math.Log2(s / m);
                double xtwx = 0, ll = 0;
                for (int j = 0; j < m; j++)
                {
                    double mu = sf[j] * Math.Pow(2, beta);
                    res.Mu[r * m + j] = mu;
                    ll += RMath.DnbinomMuLog(y[j], 1 / alpha[r], mu);
                    double w = 1 / (1 / mu + alpha[r]);
                    res.Hat[r * m + j] = w;
                    xtwx += w;
                }
                for (int j = 0; j < m; j++) res.Hat[r * m + j] /= xtwx;
                res.Beta[r] = beta;
                res.Se[r] = log2e * Math.Sqrt(1 / xtwx);
                res.Converged[r] = true;
                res.Iterations[r] = 1;
                res.LogLike[r] = ll;
            }
            return res;
        }

        var lambdaNat = new double[p];
        for (int k = 0; k < p; k++) lambdaNat[k] = lambdaLog2[k] / (ln2 * ln2);
        var gramInv = Inverse(x.Gram(), p);
        int optimRows = 0;
        var lockObj = new object();
        Parallel.For(0, n, () => new double[4 * p * p + 4 * p + 4 * m], (r, _, buf) =>
        {
            var y = counts.Row(rows[r]);
            var beta0 = new double[p];
            if (gramInv is not null)
            {
                var ly = new double[m];
                for (int j = 0; j < m; j++) ly[j] = Math.Log(y[j] / sf[j] + 0.1);
                var xty = new double[p];
                for (int j = 0; j < m; j++) for (int k = 0; k < p; k++) xty[k] += x[j, k] * ly[j];
                for (int a = 0; a < p; a++) { double s = 0; for (int k = 0; k < p; k++) s += gramInv[a * p + k] * xty[k]; beta0[a] = s; }
            }
            else
            {
                double s = 0;
                for (int j = 0; j < m; j++) s += y[j] / sf[j];
                int ic = -1;
                for (int k = 0; k < p; k++) if (x.CoefficientNames[k] == "Intercept") ic = k;
                if (ic >= 0) beta0[ic] = Math.Log(s / m); else Array.Fill(beta0, 1.0);
            }
            var fb = FitBeta(y, sf, x, alpha[r], beta0, lambdaNat, maxit, betaTol, minmu);
            bool stable = true, varPos = true;
            for (int k = 0; k < p; k++)
            {
                if (double.IsNaN(fb.Beta[k])) stable = false;
                if (!(fb.Var[k] > 0)) varPos = false;
            }
            bool conv = fb.Iter < maxit;
            double ll = 0;
            for (int j = 0; j < m; j++)
            {
                double eta = 0;
                for (int k = 0; k < p; k++) eta += x[j, k] * fb.Beta[k];
                double mu = sf[j] * Math.Exp(eta);
                res.Mu[r * m + j] = mu;
                ll += RMath.DnbinomMuLog(y[j], 1 / alpha[r], mu);
                res.Hat[r * m + j] = fb.Hat[j];
            }
            for (int k = 0; k < p; k++)
            {
                res.Beta[r * p + k] = log2e * fb.Beta[k];
                res.Se[r * p + k] = log2e * Math.Sqrt(Math.Max(fb.Var[k], 0));
            }
            res.Converged[r] = conv;
            res.Iterations[r] = fb.Iter;
            res.LogLike[r] = ll;
            bool needOptim = useOptim ? (!conv || !stable || !varPos) : (!stable || !varPos);
            if (needOptim)
            {
                lock (lockObj) optimRows++;
                bool useBetaMatrix = stable;
                for (int k = 0; k < p; k++) if (!(Math.Abs(res.Beta[r * p + k]) < Large)) useBetaMatrix = false;
                var start = new double[p];
                for (int k = 0; k < p; k++) start[k] = useBetaMatrix ? res.Beta[r * p + k] : beta0[k];
                OptimRow(y, sf, x, alpha[r], lambdaLog2, lambdaNat, minmu, start, res, r);
            }
            return buf;
        }, _ => { });
        res.OptimRows = optimRows;
        return res;
    }

    /// <summary>Per-gene output of <see cref="FitBeta"/> (natural-log scale).</summary>
    public sealed record FitBetaResult(double[] Beta, double[] Var, double[] Hat, int Iter, double[] Sigma);

    /// <summary>
    /// DESeq2 <c>fitBeta</c> (src/DESeq2.cpp, QR branch): IRLS with working weights μ/(1+αμ), ridge diag(λ), μ clamped at
    /// <paramref name="minmu"/>, convergence on |dev − dev_old|/(|dev| + 0.1) &lt; tol after the first iteration.
    /// With <paramref name="maxit"/> = 0 only the covariance at <paramref name="beta0"/> is computed (used for contrasts).
    /// </summary>
    public static FitBetaResult FitBeta(ReadOnlySpan<double> y, ReadOnlySpan<double> sf, DesignMatrix x, double alpha,
        double[] beta0, double[] lambdaNat, int maxit, double tol, double minmu)
    {
        int m = x.Samples, p = x.Coefficients;
        var beta = (double[])beta0.Clone();
        var mu = new double[m];
        UpdateMu(beta, x, sf, mu, minmu);
        double dev, devOld = 0.0;
        int iter = 0;
        var w = new double[m];
        var lhs = new double[p * p];
        var rhs = new double[p];
        for (int t = 0; t < maxit; t++)
        {
            iter++;
            Array.Clear(lhs); Array.Clear(rhs);
            for (int j = 0; j < m; j++)
            {
                w[j] = mu[j] / (1.0 + alpha * mu[j]);
                double z = Math.Log(mu[j] / sf[j]) + (y[j] - mu[j]) / mu[j];
                for (int a = 0; a < p; a++)
                {
                    double xa = x[j, a] * w[j];
                    if (xa == 0) continue;
                    rhs[a] += xa * z;
                    for (int c = 0; c < p; c++) lhs[a * p + c] += xa * x[j, c];
                }
            }
            for (int a = 0; a < p; a++) lhs[a * p + a] += lambdaNat[a];
            var nb = Solve(lhs, rhs, p);
            if (nb is null) { for (int k = 0; k < p; k++) beta[k] = double.NaN; iter = maxit; break; }
            Array.Copy(nb, beta, p);
            bool large = false;
            foreach (double v in beta) if (Math.Abs(v) > Large) large = true;
            if (large) { iter = maxit; break; }
            UpdateMu(beta, x, sf, mu, minmu);
            dev = 0.0;
            for (int j = 0; j < m; j++) dev += -2.0 * RMath.DnbinomMuLog(y[j], 1.0 / alpha, mu[j]);
            double conv = Math.Abs(dev - devOld) / (Math.Abs(dev) + 0.1);
            if (double.IsNaN(conv)) { iter = maxit; break; }
            if (t > 0 && conv < tol) break;
            devOld = dev;
        }
        // Covariance and hat diagonals at the final (clamped) mu.
        Array.Clear(lhs);
        for (int j = 0; j < m; j++)
        {
            w[j] = mu[j] / (1.0 + alpha * mu[j]);
            for (int a = 0; a < p; a++)
                for (int c = 0; c < p; c++) lhs[a * p + c] += x[j, a] * w[j] * x[j, c];
        }
        var xtwx = (double[])lhs.Clone();
        for (int a = 0; a < p; a++) lhs[a * p + a] += lambdaNat[a];
        var inv = Inverse(lhs, p);
        var hat = new double[m];
        var varr = new double[p];
        var sigma = new double[p * p];
        if (inv is null)
        {
            Array.Fill(hat, double.NaN); Array.Fill(varr, double.NaN); Array.Fill(sigma, double.NaN);
        }
        else
        {
            for (int j = 0; j < m; j++)
            {
                double s = 0;
                for (int a = 0; a < p; a++)
                    for (int c = 0; c < p; c++) s += x[j, a] * x[j, c] * inv[c * p + a];
                hat[j] = w[j] * s;
            }
            for (int a = 0; a < p; a++)
                for (int c = 0; c < p; c++)
                {
                    double s = 0;
                    for (int u = 0; u < p; u++)
                        for (int v = 0; v < p; v++) s += inv[a * p + u] * xtwx[u * p + v] * inv[v * p + c];
                    sigma[a * p + c] = s;
                }
            for (int a = 0; a < p; a++) varr[a] = sigma[a * p + a];
        }
        return new FitBetaResult(beta, varr, hat, iter, sigma);
    }

    private static void UpdateMu(double[] beta, DesignMatrix x, ReadOnlySpan<double> sf, double[] mu, double minmu)
    {
        int p = x.Coefficients;
        for (int j = 0; j < x.Samples; j++)
        {
            double eta = 0;
            for (int k = 0; k < p; k++) eta += x[j, k] * beta[k];
            mu[j] = Math.Max(sf[j] * Math.Exp(eta), minmu);
        }
    }

    /// <summary>
    /// DESeq2 <c>fitNbinomGLMsOptim</c> for one row: maximise the log2-scale log posterior (NB likelihood plus the
    /// N(0, 1/λ) ridge prior) within [−30, 30]. DESeq2 uses <c>optim(method = "L-BFGS-B")</c> with numerical gradients;
    /// here a projected Newton method with analytic derivatives finds the same optimum (to optimiser tolerance).
    /// </summary>
    private static void OptimRow(ReadOnlySpan<double> y, ReadOnlySpan<double> sf, DesignMatrix x, double alpha,
        double[] lambdaLog2, double[] lambdaNat, double minmu, double[] start, GlmFitResult res, int r)
    {
        int m = x.Samples, p = x.Coefficients;
        double ln2 = Math.Log(2);
        var yy = y.ToArray(); var ss = sf.ToArray();
        double NegLogPost(double[] pp)
        {
            double ll = 0;
            for (int j = 0; j < m; j++)
            {
                double eta = 0;
                for (int k = 0; k < p; k++) eta += x[j, k] * pp[k];
                ll += RMath.DnbinomMuLog(yy[j], 1 / alpha, ss[j] * Math.Pow(2, eta));
            }
            double lp = 0;
            for (int k = 0; k < p; k++)
            {
                double sd = Math.Sqrt(1 / lambdaLog2[k]);
                lp += -0.5 * Math.Log(2 * Math.PI) - Math.Log(sd) - 0.5 * pp[k] * pp[k] / (sd * sd);
            }
            double v = -(ll + lp);
            return double.IsFinite(v) ? v : 1e300;
        }
        var par = (double[])start.Clone();
        for (int k = 0; k < p; k++) par[k] = Math.Clamp(par[k], -Large, Large);
        double f = NegLogPost(par);
        var g = new double[p]; var h = new double[p * p];
        for (int it = 0; it < 500; it++)
        {
            Array.Clear(g); Array.Clear(h);
            for (int j = 0; j < m; j++)
            {
                double eta = 0;
                for (int k = 0; k < p; k++) eta += x[j, k] * par[k];
                double mu = ss[j] * Math.Pow(2, eta);
                double d1 = ln2 * (yy[j] - mu) / (1 + alpha * mu);
                double d2 = ln2 * ln2 * mu * (1 + alpha * yy[j]) / ((1 + alpha * mu) * (1 + alpha * mu));
                for (int a = 0; a < p; a++)
                {
                    g[a] -= x[j, a] * d1;
                    for (int c = 0; c < p; c++) h[a * p + c] += x[j, a] * d2 * x[j, c];
                }
            }
            for (int k = 0; k < p; k++) { g[k] += lambdaLog2[k] * par[k]; h[k * p + k] += lambdaLog2[k]; }
            // Free variables: not at a bound with the gradient pushing outward.
            var free = new bool[p];
            for (int k = 0; k < p; k++) free[k] = !((par[k] <= -Large && g[k] > 0) || (par[k] >= Large && g[k] < 0));
            var hs = (double[])h.Clone(); var gs = (double[])g.Clone();
            for (int k = 0; k < p; k++)
                if (!free[k]) { for (int c = 0; c < p; c++) { hs[k * p + c] = 0; hs[c * p + k] = 0; } hs[k * p + k] = 1; gs[k] = 0; }
            var step = Solve(hs, gs, p);
            if (step is null) break;
            double t = 1, fn = f;
            var np = new double[p];
            bool accepted = false;
            for (int ls = 0; ls < 60; ls++)
            {
                for (int k = 0; k < p; k++) np[k] = Math.Clamp(par[k] - t * step[k], -Large, Large);
                fn = NegLogPost(np);
                if (fn <= f) { accepted = true; break; }
                t /= 2;
            }
            if (!accepted) break;
            double maxStep = 0;
            for (int k = 0; k < p; k++) maxStep = Math.Max(maxStep, Math.Abs(np[k] - par[k]));
            Array.Copy(np, par, p);
            double fOld = f;
            f = fn;
            if (maxStep < 1e-10 || Math.Abs(fOld - f) <= 1e-15 * Math.Max(1, Math.Abs(f))) break;
        }
        res.Converged[r] = true;
        var muRow = new double[m];
        for (int j = 0; j < m; j++)
        {
            double eta = 0;
            for (int k = 0; k < p; k++) eta += x[j, k] * par[k];
            muRow[j] = ss[j] * Math.Pow(2, eta);
            res.Mu[r * m + j] = muRow[j];
            muRow[j] = Math.Max(muRow[j], minmu);
        }
        var xtwx = new double[p * p];
        for (int j = 0; j < m; j++)
        {
            double w = 1 / (1 / muRow[j] + alpha);
            for (int a = 0; a < p; a++) for (int c = 0; c < p; c++) xtwx[a * p + c] += x[j, a] * w * x[j, c];
        }
        var ridged = (double[])xtwx.Clone();
        for (int a = 0; a < p; a++) ridged[a * p + a] += lambdaNat[a];
        var inv = Inverse(ridged, p);
        const double log2e = 1.4426950408889634;
        for (int a = 0; a < p; a++)
        {
            res.Beta[r * p + a] = par[a];
            double s = double.NaN;
            if (inv is not null)
            {
                s = 0;
                for (int u = 0; u < p; u++) for (int v = 0; v < p; v++) s += inv[a * p + u] * xtwx[u * p + v] * inv[v * p + a];
            }
            res.Se[r * p + a] = log2e * Math.Sqrt(Math.Max(s, 0));
        }
        double llr = 0;
        for (int j = 0; j < m; j++) llr += RMath.DnbinomMuLog(yy[j], 1 / alpha, muRow[j]);
        res.LogLike[r] = llr;
    }

    // ---- dispersion trend -------------------------------------------------------------------------------------

    /// <summary>
    /// DESeq2 <c>parametricDispersionFit</c>: repeated <c>glm(disps ~ I(1/means), family = Gamma("identity"),
    /// start = coefs)</c> on genes whose ratio to the current fit lies in (1e-4, 15). Each glm is a port of R's
    /// <c>glm.fit</c> (IRLS from the start values, step halving, |dev − devold|/(|dev| + 0.1) &lt; 1e-8, maxit 25).
    /// </summary>
    public static bool TryParametricFit(double[] means, double[] disps, out double a0, out double a1)
    {
        a0 = 0.1; a1 = 1;
        int iter = 0;
        while (true)
        {
            var gx = new List<double>(); var gy = new List<double>();
            for (int i = 0; i < means.Length; i++)
            {
                double res = disps[i] / (a0 + a1 / means[i]);
                if (res > 1e-4 && res < 15) { gx.Add(1 / means[i]); gy.Add(disps[i]); }
            }
            if (gx.Count < 2) return false;
            if (!GlmGammaIdentity(gx, gy, a0, a1, out double n0, out double n1, out bool converged)) return false;
            double o0 = a0, o1 = a1;
            a0 = n0; a1 = n1;
            if (!(a0 > 0 && a1 > 0)) return false;
            double change = Math.Pow(Math.Log(a0 / o0), 2) + Math.Pow(Math.Log(a1 / o1), 2);
            if (change < 1e-6 && converged) return true;
            iter++;
            if (iter > 10) return false;
        }
    }

    private static double GammaDeviance(List<double> y, double c0, double c1, List<double> x)
    {
        double dev = 0;
        for (int i = 0; i < y.Count; i++)
        {
            double mu = c0 + c1 * x[i];
            dev += -2 * (Math.Log(y[i] == 0 ? 1 : y[i] / mu) - (y[i] - mu) / mu);
        }
        return dev;
    }

    private static bool GlmGammaIdentity(List<double> x, List<double> y, double s0, double s1, out double c0, out double c1, out bool converged)
    {
        const double epsilon = 1e-8;
        const int maxit = 25;
        c0 = s0; c1 = s1;
        converged = false;
        int n = x.Count;
        bool ValidMu(double b0, double b1)
        {
            for (int i = 0; i < n; i++) { double mu = b0 + b1 * x[i]; if (!(double.IsFinite(mu) && mu > 0)) return false; }
            return true;
        }
        if (!ValidMu(c0, c1)) return false; // "cannot find valid starting values"
        double old0 = c0, old1 = c1;
        double devold = GammaDeviance(y, c0, c1, x);
        for (int iter = 1; iter <= maxit; iter++)
        {
            double s00 = 0, s01 = 0, s11 = 0, t0 = 0, t1 = 0;
            for (int i = 0; i < n; i++)
            {
                double mu = c0 + c1 * x[i];
                double z = mu + (y[i] - mu); // eta + (y - mu)/mu.eta, identity link
                double w = 1 / (mu * mu);   // (weights * mu.eta^2) / variance(mu)
                s00 += w; s01 += w * x[i]; s11 += w * x[i] * x[i];
                t0 += w * z; t1 += w * x[i] * z;
            }
            double det = s00 * s11 - s01 * s01;
            if (det == 0 || !double.IsFinite(det)) return false;
            double n0 = (s11 * t0 - s01 * t1) / det;
            double n1 = (s00 * t1 - s01 * t0) / det;
            double dev = GammaDeviance(y, n0, n1, x);
            int halving = 0;
            while (!double.IsFinite(dev) || !ValidMu(n0, n1))
            {
                if (++halving > maxit) return false;
                n0 = (n0 + old0) / 2; n1 = (n1 + old1) / 2;
                dev = GammaDeviance(y, n0, n1, x);
            }
            c0 = n0; c1 = n1;
            if (Math.Abs(dev - devold) / (Math.Abs(dev) + 0.1) < epsilon) { converged = true; return true; }
            devold = dev;
            old0 = c0; old1 = c1;
        }
        return true;
    }

    /// <summary>
    /// Local mean–dispersion trend in the spirit of DESeq2's <c>localDispersionFit</c> (locfit of log dispersion on
    /// log mean, weights = means): a direct local-quadratic weighted fit with a tricube kernel over the nearest 70% of
    /// genes, evaluated exactly at each requested mean. locfit's adaptive tree interpolation is not reproduced, so values
    /// agree with DESeq2 only approximately.
    /// </summary>
    public static Func<double, double> LocalDispersionFit(double[] means, double[] disps, double minDisp)
    {
        bool allSmall = true;
        foreach (double d in disps) if (d >= minDisp * 10) allSmall = false;
        if (allSmall) return _ => minDisp;
        var xs = new List<double>(); var ys = new List<double>(); var ws = new List<double>();
        for (int i = 0; i < means.Length; i++)
            if (disps[i] >= minDisp * 10) { xs.Add(Math.Log(means[i])); ys.Add(Math.Log(disps[i])); ws.Add(means[i]); }
        var x = xs.ToArray(); var y = ys.ToArray(); var w = ws.ToArray();
        int n = x.Length;
        int k = Math.Max(3, Math.Min(n, (int)Math.Floor(0.7 * n)));
        var dist = new double[n];
        return mean =>
        {
            double x0 = Math.Log(mean);
            for (int i = 0; i < n; i++) dist[i] = Math.Abs(x[i] - x0);
            var sorted = (double[])dist.Clone();
            Array.Sort(sorted);
            double h = Math.Max(sorted[k - 1], 1e-10);
            double s00 = 0, s01 = 0, s02 = 0, s03 = 0, s04 = 0, t0 = 0, t1 = 0, t2 = 0;
            for (int i = 0; i < n; i++)
            {
                double u = dist[i] / h;
                if (u >= 1) continue;
                double kw = Math.Pow(1 - u * u * u, 3) * w[i];
                double d = x[i] - x0, d2 = d * d;
                s00 += kw; s01 += kw * d; s02 += kw * d2; s03 += kw * d2 * d; s04 += kw * d2 * d2;
                t0 += kw * y[i]; t1 += kw * d * y[i]; t2 += kw * d2 * y[i];
            }
            var sol = Solve([s00, s01, s02, s01, s02, s03, s02, s03, s04], [t0, t1, t2], 3);
            return Math.Exp(sol is null ? (s00 > 0 ? t0 / s00 : 0) : sol[0]);
        };
    }

    // ---- prior variance by simulation (estimateDispersionsPriorVar, m - p <= 3) ----------------------------------

    /// <summary>
    /// DESeq2's estimate of the dispersion prior variance when the residual degrees of freedom are 1–3: with
    /// <c>set.seed(2)</c>, for 200 candidate variances on [0, 8] simulate log(χ²_df) + N(0, var) − log(df) (10,000 draws),
    /// compare histograms (breaks −10..10 by 0.5) with the observed log dispersion residuals by KL divergence, smooth the
    /// divergences with <c>loess(span = 0.2)</c>, and take the minimiser on a 1,000-point grid (at least 0.25).
    /// </summary>
    public static double SimulatedPriorVar(IReadOnlyList<double> residuals, int df)
    {
        var rng = new RRandom(2);
        var brks = new double[41];
        for (int k = 0; k < 41; k++) brks[k] = (k - 20) / 2.0;
        var obs = residuals.Where(v => v > -10 && v < 10).ToArray();
        var obsDens = HistDensity(obs, brks);
        var grid = RSeq(0, 8, 200);
        var kl = new double[grid.Length];
        var chi = new double[10000];
        double ldf = Math.Log(df);
        for (int g = 0; g < grid.Length; g++)
        {
            for (int i = 0; i < chi.Length; i++) chi[i] = Math.Log(rng.Rchisq(df));
            double sd = Math.Sqrt(grid[g]);
            var rand = new List<double>(chi.Length);
            for (int i = 0; i < chi.Length; i++)
            {
                double v = chi[i] + rng.Rnorm(0, sd) - ldf;
                if (v > -10 && v < 10) rand.Add(v);
            }
            var randDens = HistDensity(rand.ToArray(), brks);
            double small = double.PositiveInfinity;
            foreach (double z in obsDens) if (z > 0 && z < small) small = z;
            foreach (double z in randDens) if (z > 0 && z < small) small = z;
            double sum = 0;
            for (int k = 0; k < obsDens.Length; k++)
                sum += obsDens[k] * (Math.Log(obsDens[k] + small) - Math.Log(randDens[k] + small));
            kl[g] = sum;
        }
        var fine = RSeq(0, 8, 1000);
        var fitted = RLoess.PredictInterpolate(grid, kl, 0.2, fine);
        int best = 0;
        for (int i = 1; i < fitted.Length; i++) if (fitted[i] < fitted[best]) best = i;
        return Math.Max(fine[best], 0.25);
    }

    /// <summary>R's <c>hist(x, breaks, plot = FALSE)$density</c> (right-closed bins, fuzzed breaks, include.lowest).</summary>
    private static double[] HistDensity(double[] x, double[] breaks)
    {
        int nb = breaks.Length;
        var h = new double[nb - 1];
        for (int k = 0; k < nb - 1; k++) h[k] = breaks[k + 1] - breaks[k];
        double diddle = 1e-7 * NumericUtil.Median(h);
        var fb = new double[nb];
        for (int k = 0; k < nb; k++) fb[k] = breaks[k] + (k == 0 ? -diddle : diddle);
        var counts = new int[nb - 1];
        foreach (double v in x)
        {
            if (!double.IsFinite(v)) continue;
            int lo = 0, hi = nb - 1;
            if (fb[lo] <= v && (v < fb[hi] || v == fb[hi]))
            {
                while (hi - lo >= 2)
                {
                    int mid = (hi + lo) / 2;
                    if (v > fb[mid]) lo = mid; else hi = mid;
                }
                counts[lo]++;
            }
        }
        var dens = new double[nb - 1];
        for (int k = 0; k < nb - 1; k++) dens[k] = counts[k] / (x.Length * h[k]);
        return dens;
    }

    // ---- Cook's distances --------------------------------------------------------------------------------------

    /// <summary>DESeq2 <c>robustMethodOfMomentsDisp</c> for one gene's normalised counts.</summary>
    public static double RobustMethodOfMomentsDisp(ReadOnlySpan<double> norm, int[] cellOf, bool[] cellHasThree, int cells, bool anyThree)
    {
        int m = norm.Length;
        double v;
        if (anyThree)
        {
            v = double.NegativeInfinity;
            for (int c = 0; c < cells; c++)
            {
                if (!cellHasThree[c]) continue;
                var vals = new List<double>();
                for (int j = 0; j < m; j++) if (cellOf[j] == c) vals.Add(norm[j]);
                int n = vals.Count;
                int bin = n <= 3.5 ? 0 : n <= 23.5 ? 1 : 2;
                double trim = bin == 0 ? 1.0 / 3 : bin == 1 ? 1.0 / 4 : 1.0 / 8;
                double scale = bin == 0 ? 2.04 : bin == 1 ? 1.86 : 1.51;
                var arr = vals.ToArray();
                double cm = RMath.TrimmedMean(arr, trim);
                var sq = new double[n];
                for (int t = 0; t < n; t++) sq[t] = (arr[t] - cm) * (arr[t] - cm);
                double ve = scale * RMath.TrimmedMean(sq, trim);
                if (ve > v || double.IsNaN(ve)) v = ve;
            }
        }
        else
        {
            var arr = norm.ToArray();
            double rm = RMath.TrimmedMean(arr, 1.0 / 8);
            var sq = new double[m];
            for (int t = 0; t < m; t++) sq[t] = (arr[t] - rm) * (arr[t] - rm);
            v = 1.51 * RMath.TrimmedMean(sq, 1.0 / 8);
        }
        double mean = RMath.Mean(norm);
        double alpha = (v - mean) / (mean * mean);
        return Math.Max(alpha, 0.04);
    }

    /// <summary>Cell structure used by Cook's distance (cells of identical design rows, in order of appearance).</summary>
    public static (int[] CellOf, bool[] HasThree, int Cells, bool AnyThree) CookCells(DesignMatrix x)
    {
        var id = RowGroups(x, out int g);
        var count = new int[g];
        foreach (int i in id) count[i]++;
        var has = new bool[g];
        bool any = false;
        for (int c = 0; c < g; c++) { has[c] = count[c] >= 3; any |= has[c]; }
        return (id, has, g, any);
    }
}
