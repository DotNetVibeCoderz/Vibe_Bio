//! Structure engine: optimal-superposition RMSD of two N×3 coordinate sets.
//!
//! Uses Horn's quaternion method (J. Opt. Soc. Am. A 4, 629, 1987): the best rotation is the eigenvector of the
//! largest eigenvalue of a symmetric 4×4 matrix built from the cross-covariance, found with cyclic Jacobi. Unlike
//! a 3×3 SVD it never returns an improper rotation (reflection). `ManagedKernels.KabschRmsd` mirrors this code
//! operation for operation.

/// RMSD after optimally superposing `b` (mobile) onto `a` (target); both are `n` row-major xyz triples.
/// When `rotation` is given it receives the row-major 3×3 matrix R with `R·(b − b̄) + ā ≈ a`.
pub fn kabsch_rmsd(a: &[f64], b: &[f64], n: usize, rotation: Option<&mut [f64; 9]>) -> f64 {
    if n == 0 {
        if let Some(r) = rotation {
            *r = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        }
        return 0.0;
    }
    let mut ca = [0.0f64; 3];
    let mut cb = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            ca[k] += a[3 * i + k];
            cb[k] += b[3 * i + k];
        }
    }
    let nf = n as f64;
    for k in 0..3 {
        ca[k] /= nf;
        cb[k] /= nf;
    }
    // s[p][q] = Σ x_p y_q with x = mobile (b), y = target (a).
    let mut s = [[0.0f64; 3]; 3];
    let (mut ga, mut gb) = (0.0f64, 0.0f64);
    for i in 0..n {
        let x = [b[3 * i] - cb[0], b[3 * i + 1] - cb[1], b[3 * i + 2] - cb[2]];
        let y = [a[3 * i] - ca[0], a[3 * i + 1] - ca[1], a[3 * i + 2] - ca[2]];
        for p in 0..3 {
            for q in 0..3 {
                s[p][q] += x[p] * y[q];
            }
        }
        ga += y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
        gb += x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    }
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    let mut m = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let mut v = [[0.0f64; 4]; 4];
    for (k, row) in v.iter_mut().enumerate() {
        row[k] = 1.0;
    }
    jacobi4(&mut m, &mut v);
    let mut best = 0;
    for k in 1..4 {
        if m[k][k] > m[best][best] {
            best = k;
        }
    }
    let lambda = m[best][best];
    let msd = (ga + gb - 2.0 * lambda) / nf;
    if let Some(r) = rotation {
        let (q0, q1, q2, q3) = (v[0][best], v[1][best], v[2][best], v[3][best]);
        let norm = (q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3).sqrt();
        let (q0, q1, q2, q3) = (q0 / norm, q1 / norm, q2 / norm, q3 / norm);
        *r = [
            q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
            2.0 * (q1 * q2 - q0 * q3),
            2.0 * (q1 * q3 + q0 * q2),
            2.0 * (q1 * q2 + q0 * q3),
            q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
            2.0 * (q2 * q3 - q0 * q1),
            2.0 * (q1 * q3 - q0 * q2),
            2.0 * (q2 * q3 + q0 * q1),
            q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
        ];
    }
    if msd > 0.0 {
        msd.sqrt()
    } else {
        0.0
    }
}

/// Cyclic Jacobi eigen-decomposition of a symmetric 4×4 matrix: on return the diagonal of `a` holds the
/// eigenvalues and the columns of `v` the eigenvectors.
fn jacobi4(a: &mut [[f64; 4]; 4], v: &mut [[f64; 4]; 4]) {
    let mut total = 0.0;
    for row in a.iter() {
        for x in row {
            total += x * x;
        }
    }
    if total == 0.0 {
        return;
    }
    for _sweep in 0..64 {
        let mut off = 0.0;
        for p in 0..3 {
            for q in p + 1..4 {
                off += a[p][q] * a[p][q];
            }
        }
        if off <= 1e-30 * total {
            break;
        }
        for p in 0..3 {
            for q in p + 1..4 {
                let apq = a[p][q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = if theta.abs() > 1e150 {
                    0.5 / theta
                } else {
                    let sgn = if theta >= 0.0 { 1.0 } else { -1.0 };
                    sgn / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rotate(r: &[f64; 9], p: &[f64]) -> Vec<f64> {
        p.chunks(3)
            .flat_map(|x| (0..3).map(move |i| r[3 * i] * x[0] + r[3 * i + 1] * x[1] + r[3 * i + 2] * x[2]))
            .collect()
    }

    #[test]
    fn rotated_and_translated_copy_has_zero_rmsd() {
        let a = [1.0, 2.0, 3.0, -1.0, 0.5, 2.0, 4.0, -2.0, 0.0, 0.0, 0.0, 1.0, 2.5, 1.5, -1.0];
        let (c, s) = (0.6f64.cos(), 0.6f64.sin());
        let rz = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
        let b: Vec<f64> = rotate(&rz, &a).iter().enumerate().map(|(i, x)| x + [5.0, -3.0, 2.0][i % 3]).collect();
        let mut r = [0.0; 9];
        let rmsd = kabsch_rmsd(&a, &b, 5, Some(&mut r));
        assert!(rmsd < 1e-6, "{rmsd}");
        // The recovered rotation is the inverse of rz (transpose).
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[3 * i + j] - rz[3 * j + i]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn analytic_scaling_case() {
        // a = ±1 on x, b = ±2 on y: best superposition leaves each point 1 Å away, so RMSD = 1.
        let a = [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let b = [0.0, -2.0, 0.0, 0.0, 2.0, 0.0];
        assert!((kabsch_rmsd(&a, &b, 2, None) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn mirror_image_is_not_superposable() {
        // A chiral tetrahedron and its mirror image: rotation alone cannot overlay them (RMSD > 0).
        let a = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0];
        let b = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, -3.0];
        assert!(kabsch_rmsd(&a, &b, 4, None) > 0.1);
        assert_eq!(kabsch_rmsd(&[], &[], 0, None), 0.0);
    }
}
