//! Statistics engine: Welch t-tests over matrix rows and Benjamini–Hochberg FDR.
//!
//! The special functions follow `SpecialFunctions.cs` operation for operation (same Lanczos coefficients,
//! same continued fraction), so native and managed p-values agree to the last few ulps.

const LANCZOS: [f64; 9] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];

/// ln Γ(x) for x > 0 (Lanczos, g = 7).
pub fn log_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - log_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = LANCZOS[0];
    let t = x + 7.5;
    for (i, c) in LANCZOS.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Regularised incomplete beta I_x(a, b).
pub fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_front = log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (1.0 - x).ln();
    if x < (a + 1.0) / (a + b + 2.0) {
        ln_front.exp() * beta_cf(a, b, x) / a
    } else {
        1.0 - ln_front.exp() * beta_cf(b, a, 1.0 - x) / b
    }
}

fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const EPS: f64 = 1e-15;
    const TINY: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=300 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

/// Two-sided Student t p-value.
pub fn student_t_two_sided(t: f64, df: f64) -> f64 {
    if t.is_nan() || df.is_nan() || df <= 0.0 {
        return f64::NAN;
    }
    if t.is_infinite() {
        return 0.0;
    }
    regularized_incomplete_beta(df / 2.0, 0.5, df / (df + t * t))
}

fn mean_var(row: &[f64], idx: &[i32]) -> (f64, f64) {
    let mut s = 0.0;
    for &i in idx {
        s += row[i as usize];
    }
    let mean = s / idx.len() as f64;
    let mut ss = 0.0;
    for &i in idx {
        let d = row[i as usize] - mean;
        ss += d * d;
    }
    let var = if idx.len() > 1 { ss / (idx.len() - 1) as f64 } else { 0.0 };
    (mean, var)
}

/// Row-wise Welch t-test; group column indices must be valid (checked by the caller).
#[allow(clippy::too_many_arguments)]
pub fn welch_t_test(
    matrix: &[f64],
    rows: usize,
    cols: usize,
    ga: &[i32],
    gb: &[i32],
    t: &mut [f64],
    p: &mut [f64],
    diff: &mut [f64],
) {
    for r in 0..rows {
        let row = &matrix[r * cols..(r + 1) * cols];
        let (ma, va) = mean_var(row, ga);
        let (mb, vb) = mean_var(row, gb);
        let (na, nb) = (ga.len() as f64, gb.len() as f64);
        let (sa, sb) = (va / na, vb / nb);
        let se2 = sa + sb;
        diff[r] = mb - ma;
        if se2 <= 0.0 {
            t[r] = if ma == mb { 0.0 } else { f64::INFINITY * sign(mb - ma) };
            p[r] = if ma == mb { 1.0 } else { 0.0 };
            continue;
        }
        let tv = (mb - ma) / se2.sqrt();
        let df = se2 * se2 / (sa * sa / (na - 1.0) + sb * sb / (nb - 1.0));
        t[r] = tv;
        p[r] = student_t_two_sided(tv, df);
    }
}

/// `Math.Sign` semantics: −1, 0 or +1 (NaN never reaches here because `se2 <= 0` excludes it).
fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Benjamini–Hochberg q-values; NaN p-values stay NaN and do not count towards m.
pub fn benjamini_hochberg(p: &[f64], q: &mut [f64]) {
    let mut order: Vec<usize> = Vec::with_capacity(p.len());
    for (i, &v) in p.iter().enumerate() {
        if v.is_nan() {
            q[i] = f64::NAN;
        } else {
            order.push(i);
        }
    }
    let m = order.len();
    if m == 0 {
        return;
    }
    // Ties do not matter: tied p-values always receive the same q-value.
    order.sort_unstable_by(|&x, &y| p[x].partial_cmp(&p[y]).unwrap_or(std::cmp::Ordering::Equal));
    let mut running = 1.0f64;
    for rank in (1..=m).rev() {
        let i = order[rank - 1];
        running = running.min(p[i] * m as f64 / rank as f64);
        q[i] = running.min(1.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_gamma_known_values() {
        // ln Γ(5) = ln 24; ln Γ(0.5) = ln √π.
        assert!((log_gamma(5.0) - 24f64.ln()).abs() < 1e-12);
        assert!((log_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
    }

    #[test]
    fn student_t_matches_r() {
        // R: 2 * pt(-2.0, df = 10) = 0.07338803477074042
        assert!((student_t_two_sided(2.0, 10.0) - 0.07338803477074042).abs() < 1e-12);
        // R: 2 * pt(-1.5, df = 3.7) = 0.2128...  (check symmetry and range instead of an exact value)
        let p = student_t_two_sided(-1.5, 3.7);
        assert!((p - student_t_two_sided(1.5, 3.7)).abs() < 1e-15 && p > 0.0 && p < 1.0);
        assert_eq!(student_t_two_sided(f64::INFINITY, 4.0), 0.0);
        assert!(student_t_two_sided(1.0, 0.0).is_nan());
    }

    #[test]
    fn welch_matches_r_t_test() {
        // R: t.test(c(4,5,6,7,8), c(1,2,3,4,5)) -> t = -3, df = 8, p-value = 0.01707168 (B − A orientation).
        let m = [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let (ga, gb) = ([5, 6, 7, 8, 9], [0, 1, 2, 3, 4]);
        let (mut t, mut p, mut d) = ([0.0], [0.0], [0.0]);
        welch_t_test(&m, 1, 10, &ga, &gb, &mut t, &mut p, &mut d);
        assert!((t[0] + 3.0).abs() < 1e-12);
        assert!((p[0] - 0.017071681233782634).abs() < 1e-10);
        assert_eq!(d[0], -3.0);
    }

    #[test]
    fn bh_matches_r_p_adjust() {
        // R: p.adjust(c(0.01, 0.04, 0.03, 0.005), "BH") = 0.02 0.04 0.04 0.02
        let p = [0.01, 0.04, 0.03, 0.005];
        let mut q = [0.0; 4];
        benjamini_hochberg(&p, &mut q);
        let expect = [0.02, 0.04, 0.04, 0.02];
        for k in 0..4 {
            assert!((q[k] - expect[k]).abs() < 1e-15);
        }
        // NaN is excluded from m: R p.adjust(c(0.01, NA, 0.04), "BH") = 0.02 NA 0.04
        let p = [0.01, f64::NAN, 0.04];
        let mut q = [0.0; 3];
        benjamini_hochberg(&p, &mut q);
        assert!((q[0] - 0.02).abs() < 1e-15 && q[1].is_nan() && (q[2] - 0.04).abs() < 1e-15);
    }
}
