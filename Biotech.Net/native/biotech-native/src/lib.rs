//! biotech-native: the Rust kernels behind `Biotech.Net.Native.BioKernels`.
//!
//! Every `bn_*` function is a C-ABI twin of a `ManagedKernels` method and must return identical results.
//! Callers (the .NET bridge) validate lengths and arguments; this layer additionally tolerates null pointers
//! for empty buffers and catches any panic so it can never unwind into the CLR.
//!
//! Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.

pub mod align;
pub mod protein;
pub mod seq;
pub mod stats;
pub mod structure;

use std::panic::{catch_unwind, AssertUnwindSafe};

/// Bumped whenever a signature or semantic of the C ABI changes; the .NET bridge refuses other versions.
pub const ABI_VERSION: u32 = 1;

/// Status returned when a panic was caught at the boundary.
pub const ERR_PANIC: i32 = -99;

#[inline]
fn guard<T>(fallback: T, f: impl FnOnce() -> T) -> T {
    catch_unwind(AssertUnwindSafe(f)).unwrap_or(fallback)
}

/// `slice::from_raw_parts` requires a non-null pointer even for length 0; .NET passes null for empty spans.
#[inline]
unsafe fn sl<'a, T>(p: *const T, n: usize) -> &'a [T] {
    if n == 0 || p.is_null() {
        &[]
    } else {
        std::slice::from_raw_parts(p, n)
    }
}

#[inline]
unsafe fn sl_mut<'a, T>(p: *mut T, n: usize) -> &'a mut [T] {
    if n == 0 || p.is_null() {
        &mut []
    } else {
        std::slice::from_raw_parts_mut(p, n)
    }
}

#[no_mangle]
pub extern "C" fn bn_abi_version() -> u32 {
    ABI_VERSION
}

/// Static, NUL-terminated version string.
#[no_mangle]
pub extern "C" fn bn_version() -> *const u8 {
    concat!("biotech-native ", env!("CARGO_PKG_VERSION"), "\0").as_ptr()
}

// ---- Sequence engine --------------------------------------------------------------------------------

/// # Safety
/// `src` and `dst` must point to `len` readable/writable bytes (or be null when `len` is 0).
#[no_mangle]
pub unsafe extern "C" fn bn_reverse_complement(src: *const u8, len: usize, dst: *mut u8) {
    guard((), || seq::reverse_complement(sl(src, len), sl_mut(dst, len)))
}

/// # Safety
/// As [`bn_reverse_complement`].
#[no_mangle]
pub unsafe extern "C" fn bn_complement(src: *const u8, len: usize, dst: *mut u8) {
    guard((), || seq::complement(sl(src, len), sl_mut(dst, len)))
}

/// # Safety
/// `table` must hold 64 bytes and `protein` `len / 3` bytes.
#[no_mangle]
pub unsafe extern "C" fn bn_translate(nuc: *const u8, len: usize, table: *const u8, protein: *mut u8) -> usize {
    guard(0, || {
        let t: &[u8; 64] = match sl(table, 64).try_into() {
            Ok(t) => t,
            Err(_) => return 0,
        };
        seq::translate(sl(nuc, len), t, sl_mut(protein, len / 3))
    })
}

/// # Safety
/// `out` must hold 6 i64 values (A, C, G, T, N, other).
#[no_mangle]
pub unsafe extern "C" fn bn_count_nucleotides(s: *const u8, len: usize, out: *mut i64) {
    guard((), || sl_mut(out, 6).copy_from_slice(&seq::count_nucleotides(sl(s, len))))
}

/// # Safety
/// `counts` must hold 256 i64 values; they are added to.
#[no_mangle]
pub unsafe extern "C" fn bn_byte_histogram(data: *const u8, len: usize, counts: *mut i64) {
    guard((), || {
        if let Ok(c) = <&mut [i64; 256]>::try_from(sl_mut(counts, 256)) {
            seq::byte_histogram(sl(data, len), c);
        }
    })
}

/// # Safety
/// `counts` must hold 4^k i32 values; 1 ≤ k ≤ 12.
#[no_mangle]
pub unsafe extern "C" fn bn_count_kmers(s: *const u8, len: usize, k: u32, counts: *mut i32) -> i64 {
    if !(1..=12).contains(&k) {
        return -1;
    }
    guard(-1, || seq::count_kmers(sl(s, len), k, sl_mut(counts, 1usize << (2 * k))))
}

/// # Safety
/// `a` and `b` must each hold `len` bytes.
#[no_mangle]
pub unsafe extern "C" fn bn_hamming(a: *const u8, b: *const u8, len: usize) -> usize {
    guard(0, || seq::hamming(sl(a, len), sl(b, len)))
}

// ---- Statistics engine ------------------------------------------------------------------------------

/// # Safety
/// `matrix` holds rows×cols values, group indices are < cols, outputs hold `rows` values.
#[no_mangle]
pub unsafe extern "C" fn bn_welch_ttest(
    matrix: *const f64,
    rows: usize,
    cols: usize,
    ga: *const i32,
    na: usize,
    gb: *const i32,
    nb: usize,
    t: *mut f64,
    p: *mut f64,
    diff: *mut f64,
) {
    guard((), || {
        stats::welch_t_test(
            sl(matrix, rows * cols),
            rows,
            cols,
            sl(ga, na),
            sl(gb, nb),
            sl_mut(t, rows),
            sl_mut(p, rows),
            sl_mut(diff, rows),
        )
    })
}

/// # Safety
/// `p` and `q` hold `n` values.
#[no_mangle]
pub unsafe extern "C" fn bn_benjamini_hochberg(p: *const f64, n: usize, q: *mut f64) {
    guard((), || stats::benjamini_hochberg(sl(p, n), sl_mut(q, n)))
}

// ---- Protein engine ---------------------------------------------------------------------------------

/// # Safety
/// `s` holds `len` bytes.
#[no_mangle]
pub unsafe extern "C" fn bn_protein_mass(s: *const u8, len: usize, monoisotopic: i32) -> f64 {
    guard(f64::NAN, || protein::protein_mass(sl(s, len), monoisotopic != 0))
}

/// # Safety
/// `s` holds `len` bytes.
#[no_mangle]
pub unsafe extern "C" fn bn_net_charge(s: *const u8, len: usize, ph: f64) -> f64 {
    guard(f64::NAN, || protein::net_charge(sl(s, len), ph))
}

/// # Safety
/// `s` holds `len` bytes.
#[no_mangle]
pub unsafe extern "C" fn bn_isoelectric_point(s: *const u8, len: usize) -> f64 {
    guard(f64::NAN, || protein::isoelectric_point(sl(s, len)))
}

/// # Safety
/// `out` holds `max(0, len − window + 1)` values.
#[no_mangle]
pub unsafe extern "C" fn bn_hydropathy(s: *const u8, len: usize, window: usize, out: *mut f64) -> usize {
    guard(0, || {
        let n = if window == 0 || len < window { 0 } else { len - window + 1 };
        protein::hydropathy(sl(s, len), window, sl_mut(out, n))
    })
}

// ---- Alignment engine -------------------------------------------------------------------------------

/// Pairwise alignment with traceback. `mode`: 0 global, 1 local, 2 semi-global. `band < 0` = full matrix.
/// `ops` (capacity `ops_cap` ≥ alen + blen) receives `M`/`I`/`D`; `result` receives score, a_start, a_end,
/// b_start, b_end, ops_len. Returns 0 or a negative error code.
///
/// # Safety
/// Pointers must be valid for the given lengths; `matrix` holds 128×128 i32 and `result` 6 i32.
#[no_mangle]
pub unsafe extern "C" fn bn_align(
    a: *const u8,
    alen: usize,
    b: *const u8,
    blen: usize,
    matrix: *const i32,
    gap_open: i32,
    gap_extend: i32,
    mode: i32,
    band: i64,
    ops: *mut u8,
    ops_cap: usize,
    result: *mut i32,
) -> i32 {
    guard(ERR_PANIC, || {
        let Some(mode) = align::Mode::from_i32(mode) else { return align::ERR_BAD_ARGUMENT };
        if result.is_null() || matrix.is_null() {
            return align::ERR_BAD_ARGUMENT;
        }
        match align::align(sl(a, alen), sl(b, blen), sl(matrix, align::MATRIX_LEN), gap_open, gap_extend, mode, band, sl_mut(ops, ops_cap)) {
            Ok(s) => {
                let r = sl_mut(result, 6);
                r[0] = s.score;
                r[1] = s.a_start as i32;
                r[2] = s.a_end as i32;
                r[3] = s.b_start as i32;
                r[4] = s.b_end as i32;
                r[5] = s.ops_len as i32;
                0
            }
            Err(code) => code,
        }
    })
}

/// Score-only local alignment (linear memory, AVX2 striped when available). Returns the score, or a
/// negative error code.
///
/// # Safety
/// Pointers valid for their lengths; `matrix` holds 128×128 i32.
#[no_mangle]
pub unsafe extern "C" fn bn_local_score(q: *const u8, qlen: usize, t: *const u8, tlen: usize, matrix: *const i32, gap_open: i32, gap_extend: i32) -> i32 {
    guard(ERR_PANIC, || {
        if matrix.is_null() {
            return align::ERR_BAD_ARGUMENT;
        }
        align::local_score(sl(q, qlen), sl(t, tlen), sl(matrix, align::MATRIX_LEN), gap_open, gap_extend).unwrap_or_else(|c| c)
    })
}

/// Scores one query against `count` targets stored back to back in `targets` (`offsets` holds count + 1
/// ascending positions). Returns 0 or a negative error code.
///
/// # Safety
/// Pointers valid for their lengths; `scores` holds `count` i32.
#[no_mangle]
pub unsafe extern "C" fn bn_local_score_batch(
    q: *const u8,
    qlen: usize,
    targets: *const u8,
    targets_len: usize,
    offsets: *const i32,
    count: usize,
    matrix: *const i32,
    gap_open: i32,
    gap_extend: i32,
    scores: *mut i32,
) -> i32 {
    guard(ERR_PANIC, || {
        if matrix.is_null() || offsets.is_null() {
            return align::ERR_BAD_ARGUMENT;
        }
        match align::local_score_batch(sl(q, qlen), sl(targets, targets_len), sl(offsets, count + 1), sl(matrix, align::MATRIX_LEN), gap_open, gap_extend, sl_mut(scores, count)) {
            Ok(()) => 0,
            Err(c) => c,
        }
    })
}

// ---- Structure engine -------------------------------------------------------------------------------

/// Optimal-superposition RMSD of two N×3 coordinate sets; `rotation` (nullable) receives a row-major 3×3.
///
/// # Safety
/// `a` and `b` hold 3n doubles; `rotation` is null or holds 9 doubles.
#[no_mangle]
pub unsafe extern "C" fn bn_kabsch_rmsd(a: *const f64, b: *const f64, n: usize, rotation: *mut f64) -> f64 {
    guard(f64::NAN, || {
        let rot = if rotation.is_null() { None } else { <&mut [f64; 9]>::try_from(sl_mut(rotation, 9)).ok() };
        structure::kabsch_rmsd(sl(a, 3 * n), sl(b, 3 * n), n, rot)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_and_abi() {
        assert_eq!(bn_abi_version(), ABI_VERSION);
        let v = unsafe { std::ffi::CStr::from_ptr(bn_version() as *const std::ffi::c_char) };
        assert!(v.to_str().unwrap().starts_with("biotech-native "));
    }

    #[test]
    fn ffi_tolerates_null_for_empty_buffers() {
        unsafe {
            bn_reverse_complement(std::ptr::null(), 0, std::ptr::null_mut());
            assert_eq!(bn_hamming(std::ptr::null(), std::ptr::null(), 0), 0);
            assert_eq!(bn_protein_mass(std::ptr::null(), 0, 0), 0.0);
            assert_eq!(bn_count_kmers(std::ptr::null(), 0, 13, std::ptr::null_mut()), -1);
        }
    }

    #[test]
    fn ffi_align_reports_errors() {
        let m = vec![0i32; align::MATRIX_LEN];
        let mut r = [0i32; 6];
        let mut ops = [0u8; 2];
        let status = unsafe { bn_align(b"AC".as_ptr(), 2, b"AC".as_ptr(), 2, m.as_ptr(), 1, 1, 0, -1, ops.as_mut_ptr(), 2, r.as_mut_ptr()) };
        assert_eq!(status, align::ERR_BUFFER_TOO_SMALL);
        let status = unsafe { bn_align(b"AC".as_ptr(), 2, b"AC".as_ptr(), 2, m.as_ptr(), 1, 1, 7, -1, ops.as_mut_ptr(), 2, r.as_mut_ptr()) };
        assert_eq!(status, align::ERR_BAD_ARGUMENT);
    }
}
