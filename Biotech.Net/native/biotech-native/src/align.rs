//! Alignment engine: affine-gap (Gotoh) global, local and semi-global alignment with full traceback, plus a
//! linear-memory score-only Smith–Waterman (AVX2 striped, Farrar 2007) for database search.
//!
//! Conventions (shared with `ManagedAlignment.cs`, which must produce identical output):
//! * `matrix` is 128×128 i32 indexed by ASCII byte (`(a & 0x7F) << 7 | (b & 0x7F)`), rows = sequence `a`.
//! * A gap of length k costs `open + (k − 1)·extend` (EMBOSS/Biopython convention); both are non-negative.
//! * Operations: `M` (a[i] aligned to b[j]), `I` (b residue against a gap), `D` (a residue against a gap).
//! * Ties: H prefers diagonal, then E (`I`), then F (`D`); E/F prefer opening over extending; local and
//!   semi-global keep the first maximum in scan order.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// Entries of an ASCII-indexed substitution matrix.
pub const MATRIX_LEN: usize = 128 * 128;
const NEG: i32 = i32::MIN / 4;

/// Error codes shared with the C ABI.
pub const ERR_OUT_OF_MEMORY: i32 = -2;
pub const ERR_BUFFER_TOO_SMALL: i32 = -3;
pub const ERR_BAD_ARGUMENT: i32 = -4;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    Global,
    Local,
    SemiGlobal,
}

impl Mode {
    pub fn from_i32(v: i32) -> Option<Mode> {
        match v {
            0 => Some(Mode::Global),
            1 => Some(Mode::Local),
            2 => Some(Mode::SemiGlobal),
            _ => None,
        }
    }
}

/// Score and coordinates (0-based, half-open) of an alignment whose operations were written to the caller's buffer.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Summary {
    pub score: i32,
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
    pub ops_len: usize,
}

#[inline]
fn gap(o: i32, e: i32, k: usize) -> i32 {
    o + (k as i32 - 1) * e
}

/// Diagonal band `dmin ≤ j − i ≤ dmax`; widened by |m − n| so the corners stay reachable.
struct Band {
    banded: bool,
    dmin: i64,
    dmax: i64,
}

impl Band {
    fn new(n: usize, m: usize, band: i64) -> Band {
        let diff = m as i64 - n as i64;
        if band < 0 {
            Band { banded: false, dmin: 0, dmax: 0 }
        } else {
            Band { banded: true, dmin: diff.min(0) - band, dmax: diff.max(0) + band }
        }
    }

    #[inline]
    fn row(&self, i: usize, m: usize) -> (usize, usize) {
        if !self.banded {
            return (1, m);
        }
        let lo = (i as i64 + self.dmin).max(1) as usize;
        let hi = (i as i64 + self.dmax).min(m as i64) as usize;
        (lo, hi)
    }

    #[inline]
    fn contains(&self, i: usize, j: usize) -> bool {
        let d = j as i64 - i as i64;
        !self.banded || (self.dmin <= d && d <= self.dmax)
    }
}

fn try_vec<T: Clone>(len: usize, value: T) -> Result<Vec<T>, i32> {
    let mut v = Vec::new();
    v.try_reserve_exact(len).map_err(|_| ERR_OUT_OF_MEMORY)?;
    v.resize(len, value);
    Ok(v)
}

fn trivial(n: usize, m: usize, o: i32, e: i32, mode: Mode, ops: &mut [u8]) -> Summary {
    if mode == Mode::Local {
        return Summary::default();
    }
    let k = n + m;
    ops[..n].fill(b'D');
    ops[n..k].fill(b'I');
    let score = if mode == Mode::Global && k > 0 { -gap(o, e, k) } else { 0 };
    Summary { score, a_start: 0, a_end: n, b_start: 0, b_end: m, ops_len: k }
}

/// Aligns `a` against `b`. `band < 0` means a full matrix; otherwise cells farther than `band` diagonals from
/// the main band are excluded. `ops` needs `a.len() + b.len()` bytes and receives the operations in order.
#[allow(clippy::too_many_arguments)]
pub fn align(a: &[u8], b: &[u8], matrix: &[i32], o: i32, e: i32, mode: Mode, band: i64, ops: &mut [u8]) -> Result<Summary, i32> {
    let (n, m) = (a.len(), b.len());
    if matrix.len() < MATRIX_LEN || o < 0 || e < 0 {
        return Err(ERR_BAD_ARGUMENT);
    }
    if ops.len() < n + m {
        return Err(ERR_BUFFER_TOO_SMALL);
    }
    if n == 0 || m == 0 {
        return Ok(trivial(n, m, o, e, mode, ops));
    }
    let bd = Band::new(n, m, band);
    let local = mode == Mode::Local;
    let boundary = |k: usize| if mode == Mode::Global { -gap(o, e, k) } else { 0 };

    // Traceback rows are stored back to back; row i covers columns lo[i]..=hi[i].
    let mut row_lo = try_vec(n + 1, 0usize)?;
    let mut row_off = try_vec(n + 2, 0usize)?;
    for i in 1..=n {
        let (lo, hi) = bd.row(i, m);
        row_lo[i] = lo;
        row_off[i + 1] = row_off[i] + (hi + 1 - lo);
    }
    let mut tb = try_vec(row_off[n + 1], 0u8)?;
    let mut hprev = try_vec(m + 1, NEG)?;
    let mut hcur = try_vec(m + 1, NEG)?;
    let mut f = try_vec(m + 1, NEG)?;

    hprev[0] = 0;
    for j in 1..=m {
        hprev[j] = if bd.contains(0, j) { boundary(j) } else { NEG };
    }
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    let (mut col_best, mut col_i) = (hprev[m], 0usize);

    for i in 1..=n {
        let lo = row_lo[i];
        let hi = lo + (row_off[i + 1] - row_off[i]) - 1;
        let srow = &matrix[((a[i - 1] & 0x7F) as usize) << 7..][..128];
        let left = if lo == 1 && bd.contains(i, 0) { boundary(i) } else { NEG };
        hcur[lo - 1] = left;
        let mut h_left = left;
        let mut e_run = NEG;
        let mut diag = hprev[lo - 1];
        let trow = &mut tb[row_off[i]..row_off[i + 1]];
        for j in lo..=hi {
            let s = srow[(b[j - 1] & 0x7F) as usize];
            let (e_open, e_ext) = (h_left - o, e_run - e);
            let (ev, eb) = if e_ext > e_open { (e_ext, 4u8) } else { (e_open, 0u8) };
            let up = hprev[j];
            let (f_open, f_ext) = (up - o, f[j] - e);
            let (fv, fb) = if f_ext > f_open { (f_ext, 8u8) } else { (f_open, 0u8) };
            let d = diag + s;
            let (mut h, mut src) = if d >= ev && d >= fv {
                (d, 0u8)
            } else if ev >= fv {
                (ev, 1u8)
            } else {
                (fv, 2u8)
            };
            if local && h <= 0 {
                h = 0;
                src = 3;
            }
            trow[j - lo] = src | eb | fb;
            diag = up;
            hcur[j] = h;
            f[j] = fv;
            e_run = ev;
            h_left = h;
            if local && h > best {
                best = h;
                bi = i;
                bj = j;
            }
        }
        if hi < m {
            hcur[hi + 1] = NEG;
            f[hi + 1] = NEG;
        }
        if mode == Mode::SemiGlobal && i < n && hi == m && hcur[m] > col_best {
            col_best = hcur[m];
            col_i = i;
        }
        std::mem::swap(&mut hprev, &mut hcur);
    }

    // hprev now holds row n.
    let (score, ei, ej) = match mode {
        Mode::Global => (hprev[m], n, m),
        Mode::Local => (best, bi, bj),
        Mode::SemiGlobal => {
            let (mut s, mut ei, mut ej) = (hprev[m], n, m);
            let first = row_lo[n] - 1;
            for j in first..m {
                if hprev[j] > s {
                    s = hprev[j];
                    ej = j;
                }
            }
            if col_best > s {
                s = col_best;
                ei = col_i;
                ej = m;
            }
            (s, ei, ej)
        }
    };

    // Traceback, writing operations in reverse.
    let mut k = 0usize;
    let (mut i, mut j) = (ei, ej);
    if mode == Mode::SemiGlobal {
        for _ in j..m {
            ops[k] = b'I';
            k += 1;
        }
        for _ in i..n {
            ops[k] = b'D';
            k += 1;
        }
    }
    let mut state = 0u8;
    while i > 0 && j > 0 {
        let lo = row_lo[i];
        let width = row_off[i + 1] - row_off[i];
        if j < lo || j >= lo + width {
            break; // only reachable through unreachable (NEG) cells; never on an optimal path
        }
        let t = tb[row_off[i] + j - lo];
        match state {
            0 => match t & 3 {
                0 => {
                    ops[k] = b'M';
                    k += 1;
                    i -= 1;
                    j -= 1;
                }
                1 => state = 1,
                2 => state = 2,
                _ => break,
            },
            1 => {
                ops[k] = b'I';
                k += 1;
                j -= 1;
                state = if t & 4 != 0 { 1 } else { 0 };
            }
            _ => {
                ops[k] = b'D';
                k += 1;
                i -= 1;
                state = if t & 8 != 0 { 2 } else { 0 };
            }
        }
    }
    let (sa, sb) = (i, j);
    if !local {
        for _ in 0..i {
            ops[k] = b'D';
            k += 1;
        }
        for _ in 0..j {
            ops[k] = b'I';
            k += 1;
        }
    }
    ops[..k].reverse();
    Ok(if local {
        Summary { score, a_start: sa, a_end: ei, b_start: sb, b_end: ej, ops_len: k }
    } else {
        Summary { score, a_start: 0, a_end: n, b_start: 0, b_end: m, ops_len: k }
    })
}

// ---- Score-only local alignment -------------------------------------------------------------------

/// 16 × i16 lanes, 32-byte aligned so AVX2 can use aligned loads.
#[derive(Clone, Copy)]
#[repr(C, align(32))]
struct V16([i16; 16]);

/// A query prepared once and scored against many targets.
pub struct LocalProfile {
    qlen: usize,
    o: i32,
    e: i32,
    /// `scalar[c * qlen + j] = matrix[q[j]][c]`: one contiguous row per target residue.
    scalar: Vec<i32>,
    /// Striped 16-bit profile (Farrar); empty when scores do not fit in i16 or AVX2 is unavailable.
    striped: Vec<V16>,
    seg: usize,
    limit: i32,
}

impl LocalProfile {
    pub fn new(q: &[u8], matrix: &[i32], o: i32, e: i32) -> Result<LocalProfile, i32> {
        if matrix.len() < MATRIX_LEN || o < 0 || e < 0 {
            return Err(ERR_BAD_ARGUMENT);
        }
        let qlen = q.len();
        let mut scalar = try_vec(128 * qlen, 0i32)?;
        let (mut lo, mut hi) = (0i32, 0i32);
        for c in 0..128 {
            let row = &mut scalar[c * qlen..(c + 1) * qlen];
            for (j, &r) in q.iter().enumerate() {
                let v = matrix[(((r & 0x7F) as usize) << 7) | c];
                row[j] = v;
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        let mut p = LocalProfile { qlen, o, e, scalar, striped: Vec::new(), seg: 0, limit: 0 };
        let fits = lo >= i16::MIN as i32 / 2 && hi <= i16::MAX as i32 / 2 && o <= i16::MAX as i32 && e <= i16::MAX as i32;
        if qlen > 0 && fits && simd_available() {
            let seg = qlen.div_ceil(16);
            let mut striped = try_vec(128 * seg, V16([0; 16]))?;
            for c in 0..128 {
                let row = &p.scalar[c * qlen..(c + 1) * qlen];
                for s in 0..seg {
                    let v = &mut striped[c * seg + s].0;
                    for (lane, slot) in v.iter_mut().enumerate() {
                        let pos = lane * seg + s;
                        // Padding sits after the end of the query and only ever repeats existing scores.
                        *slot = if pos < qlen { row[pos] as i16 } else { 0 };
                    }
                }
            }
            p.striped = striped;
            p.seg = seg;
            p.limit = i16::MAX as i32 - hi;
        }
        Ok(p)
    }

    /// Best local alignment score of the query against `t`.
    pub fn score(&self, t: &[u8], work: &mut Workspace) -> i32 {
        if self.qlen == 0 || t.is_empty() {
            return 0;
        }
        #[cfg(target_arch = "x86_64")]
        if !self.striped.is_empty() {
            // SAFETY: the striped profile only exists when AVX2 was detected.
            if let Some(s) = unsafe { self.striped_avx2(t, work) } {
                return s;
            }
        }
        self.scalar_score(t, work)
    }

    fn scalar_score(&self, t: &[u8], work: &mut Workspace) -> i32 {
        let q = self.qlen;
        let (o, e) = (self.o, self.e);
        work.h.clear();
        work.h.resize(q, 0);
        work.f.clear();
        work.f.resize(q, 0);
        let (h, f) = (&mut work.h[..q], &mut work.f[..q]);
        let mut best = 0i32;
        // E and F are floored at zero: max(0, E) satisfies the same recurrence and leaves H unchanged.
        for &c in t {
            let row = &self.scalar[(c & 0x7F) as usize * q..][..q];
            let (mut diag, mut h_left, mut ev) = (0i32, 0i32, 0i32);
            for j in 0..q {
                ev = (h_left - o).max(ev - e).max(0);
                let up = h[j];
                let fv = (up - o).max(f[j] - e).max(0);
                f[j] = fv;
                let hv = (diag + row[j]).max(ev).max(fv).max(0);
                diag = up;
                h[j] = hv;
                h_left = hv;
                best = best.max(hv);
            }
        }
        best
    }

    /// Farrar's striped Smith–Waterman on 16 × i16 lanes. Returns `None` if a score could have saturated.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn striped_avx2(&self, t: &[u8], work: &mut Workspace) -> Option<i32> {
        let seg = self.seg;
        work.v.clear();
        work.v.resize(3 * seg, V16([0; 16]));
        let (store, rest) = work.v.split_at_mut(seg);
        let (load, earr) = rest.split_at_mut(seg);
        let mut ps = store.as_mut_ptr() as *mut __m256i;
        let mut pl = load.as_mut_ptr() as *mut __m256i;
        let pe = earr.as_mut_ptr() as *mut __m256i;
        let prof = self.striped.as_ptr() as *const __m256i;
        let zero = _mm256_setzero_si256();
        let vo = _mm256_set1_epi16(self.o as i16);
        let ve = _mm256_set1_epi16(self.e as i16);
        let mut vmax = zero;

        for &c in t {
            let vp = prof.add((c & 0x7F) as usize * seg);
            let mut vf = zero;
            let mut vh = shift_lanes(_mm256_load_si256(ps.add(seg - 1)));
            std::mem::swap(&mut ps, &mut pl);
            for j in 0..seg {
                vh = _mm256_adds_epi16(vh, _mm256_load_si256(vp.add(j)));
                let e_j = _mm256_load_si256(pe.add(j));
                vh = _mm256_max_epi16(_mm256_max_epi16(vh, e_j), _mm256_max_epi16(vf, zero));
                vmax = _mm256_max_epi16(vmax, vh);
                _mm256_store_si256(ps.add(j), vh);
                let vh_o = _mm256_max_epi16(_mm256_subs_epi16(vh, vo), zero);
                let e_next = _mm256_max_epi16(_mm256_max_epi16(_mm256_subs_epi16(e_j, ve), zero), vh_o);
                _mm256_store_si256(pe.add(j), e_next);
                vf = _mm256_max_epi16(_mm256_max_epi16(_mm256_subs_epi16(vf, ve), zero), vh_o);
                vh = _mm256_load_si256(pl.add(j));
            }
            // Lazy-F: carry vertical gaps across segment boundaries until they stop improving anything.
            vf = shift_lanes(vf);
            let mut j = 0;
            loop {
                let vh = _mm256_load_si256(ps.add(j));
                let vh_o = _mm256_max_epi16(_mm256_subs_epi16(vh, vo), zero);
                if _mm256_movemask_epi8(_mm256_cmpgt_epi16(vf, vh_o)) == 0 {
                    break;
                }
                // Only lanes that F actually raised can open a new gap; the others already did in the first pass,
                // and re-seeding F from them would regenerate the same carry forever.
                let raised = _mm256_cmpgt_epi16(vf, vh);
                let vh = _mm256_max_epi16(vh, vf);
                _mm256_store_si256(ps.add(j), vh);
                vmax = _mm256_max_epi16(vmax, vh);
                let vh_o = _mm256_max_epi16(_mm256_subs_epi16(vh, vo), zero);
                _mm256_store_si256(pe.add(j), _mm256_max_epi16(_mm256_load_si256(pe.add(j)), vh_o));
                vf = _mm256_max_epi16(_mm256_max_epi16(_mm256_subs_epi16(vf, ve), zero), _mm256_and_si256(raised, vh_o));
                j += 1;
                if j == seg {
                    vf = shift_lanes(vf);
                    j = 0;
                }
            }
        }
        let mut lanes = [0i16; 16];
        _mm256_storeu_si256(lanes.as_mut_ptr() as *mut __m256i, vmax);
        let best = lanes.iter().copied().max().unwrap_or(0) as i32;
        if best >= self.limit {
            None
        } else {
            Some(best)
        }
    }
}

/// Shifts the 16 lanes up by one (lane k ← lane k−1, lane 0 ← 0) across the two 128-bit halves.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn shift_lanes(v: __m256i) -> __m256i {
    let carry = _mm256_permute2x128_si256::<0x08>(v, v);
    _mm256_alignr_epi8::<14>(v, carry)
}

fn simd_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Scratch buffers reused across targets.
#[derive(Default)]
pub struct Workspace {
    h: Vec<i32>,
    f: Vec<i32>,
    v: Vec<V16>,
}

/// Local alignment score of `q` (rows of the matrix) against `t`.
pub fn local_score(q: &[u8], t: &[u8], matrix: &[i32], o: i32, e: i32) -> Result<i32, i32> {
    let p = LocalProfile::new(q, matrix, o, e)?;
    Ok(p.score(t, &mut Workspace::default()))
}

/// Scores `q` against every target `targets[offsets[k]..offsets[k+1]]`.
pub fn local_score_batch(q: &[u8], targets: &[u8], offsets: &[i32], matrix: &[i32], o: i32, e: i32, out: &mut [i32]) -> Result<(), i32> {
    if offsets.is_empty() || out.len() < offsets.len() - 1 {
        return Err(ERR_BAD_ARGUMENT);
    }
    for w in offsets.windows(2) {
        if w[0] < 0 || w[1] < w[0] || w[1] as usize > targets.len() {
            return Err(ERR_BAD_ARGUMENT);
        }
    }
    let p = LocalProfile::new(q, matrix, o, e)?;
    let mut work = Workspace::default();
    for (k, w) in offsets.windows(2).enumerate() {
        out[k] = p.score(&targets[w[0] as usize..w[1] as usize], &mut work);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_matrix(mat: i32, mis: i32) -> Vec<i32> {
        let mut m = vec![mis; MATRIX_LEN];
        for c in 0..128 {
            m[(c << 7) | c] = mat;
        }
        m
    }

    fn run(a: &[u8], b: &[u8], m: &[i32], o: i32, e: i32, mode: Mode, band: i64) -> (Summary, String) {
        let mut ops = vec![0u8; a.len() + b.len()];
        let s = align(a, b, m, o, e, mode, band, &mut ops).unwrap();
        (s, String::from_utf8(ops[..s.ops_len].to_vec()).unwrap())
    }

    /// Score recomputed from the operations, to check traceback/score consistency.
    fn rescore(a: &[u8], b: &[u8], m: &[i32], o: i32, e: i32, s: &Summary, ops: &str, mode: Mode) -> i32 {
        let (mut i, mut j, mut total) = (s.a_start, s.b_start, 0);
        let bytes = ops.as_bytes();
        let mut k = 0;
        while k < bytes.len() {
            let op = bytes[k];
            if op == b'M' {
                total += m[((a[i] as usize) << 7) | b[j] as usize];
                i += 1;
                j += 1;
                k += 1;
                continue;
            }
            let start = k;
            while k < bytes.len() && bytes[k] == op {
                k += 1;
            }
            let len = k - start;
            let terminal = start == 0 || k == bytes.len();
            if !(mode == Mode::SemiGlobal && terminal) {
                total -= gap(o, e, len);
            }
            if op == b'D' {
                i += len;
            } else {
                j += len;
            }
        }
        total
    }

    #[test]
    fn global_textbook_example() {
        // Durbin et al. (1998) HEAGAWGHEE vs PAWHEAE with BLOSUM50 and linear gap −8 gives score 1.
        // Here: simple +1/−1 with linear gap 1 (o = e = 1) on GATTACA/GCATGCU (Wikipedia NW example) gives 0.
        let m = simple_matrix(1, -1);
        let (s, ops) = run(b"GATTACA", b"GCATGCU", &m, 1, 1, Mode::Global, -1);
        assert_eq!(s.score, 0);
        assert_eq!(rescore(b"GATTACA", b"GCATGCU", &m, 1, 1, &s, &ops, Mode::Global), 0);
    }

    #[test]
    fn local_finds_embedded_match() {
        let m = simple_matrix(2, -3);
        let (s, ops) = run(b"TTTTACGTACGTTTTT", b"GGACGTACGGG", &m, 5, 2, Mode::Local, -1);
        assert_eq!(s.score, 14); // ACGTACG, 7 matches × 2
        assert_eq!(ops, "MMMMMMM");
        assert_eq!((s.a_start, s.a_end, s.b_start, s.b_end), (4, 11, 2, 9));
    }

    #[test]
    fn semiglobal_ignores_end_gaps() {
        let m = simple_matrix(1, -2);
        let (s, ops) = run(b"ACGTACGT", b"GTAC", &m, 3, 1, Mode::SemiGlobal, -1);
        assert_eq!(s.score, 4);
        assert_eq!(ops, "DDMMMMDD");
    }

    #[test]
    fn affine_gap_prefers_one_long_gap() {
        let m = simple_matrix(5, -4);
        let (s, ops) = run(b"AAAACCCCGGGG", b"AAAAGGGG", &m, 10, 1, Mode::Global, -1);
        assert_eq!(s.score, 8 * 5 - (10 + 3));
        assert_eq!(ops, "MMMMDDDDMMMM");
    }

    #[test]
    fn empty_inputs() {
        let m = simple_matrix(1, -1);
        let (s, ops) = run(b"", b"ACG", &m, 4, 1, Mode::Global, -1);
        assert_eq!((s.score, ops.as_str()), (-6, "III"));
        let (s, _) = run(b"", b"", &m, 4, 1, Mode::Global, -1);
        assert_eq!(s.score, 0);
        let (s, ops) = run(b"ACG", b"", &m, 4, 1, Mode::Local, -1);
        assert_eq!((s.score, ops.as_str()), (0, ""));
    }

    #[test]
    fn banded_equals_full_when_band_is_wide() {
        let m = simple_matrix(2, -3);
        let mut seed = 11u64;
        for _ in 0..50 {
            let a = random_seq(&mut seed, 40);
            let b = mutate(&mut seed, &a);
            for mode in [Mode::Global, Mode::Local, Mode::SemiGlobal] {
                let full = run(&a, &b, &m, 5, 2, mode, -1);
                let banded = run(&a, &b, &m, 5, 2, mode, 100);
                assert_eq!(full, banded);
                assert_eq!(rescore(&a, &b, &m, 5, 2, &full.0, &full.1, mode), full.0.score);
            }
            // A narrow band never beats the full matrix.
            assert!(run(&a, &b, &m, 5, 2, Mode::Global, 2).0.score <= run(&a, &b, &m, 5, 2, Mode::Global, -1).0.score);
        }
    }

    fn random_seq(seed: &mut u64, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((*seed >> 33) % 4) as usize]
            })
            .collect()
    }

    fn mutate(seed: &mut u64, a: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for &c in a {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            match (*seed >> 33) % 10 {
                0 => {}
                1 => {
                    out.push(c);
                    out.push(b'A');
                }
                2 => out.push(b'T'),
                _ => out.push(c),
            }
        }
        out
    }

    #[test]
    fn score_only_matches_full_local() {
        let mut seed = 5u64;
        let m = simple_matrix(3, -2);
        for len in [1usize, 15, 16, 17, 33, 100, 257] {
            let q = random_seq(&mut seed, len);
            let base = random_seq(&mut seed, len + 20);
            let t = mutate(&mut seed, &base);
            for (o, e) in [(5, 2), (1, 1), (0, 0), (11, 1)] {
                let full = run(&q, &t, &m, o, e, Mode::Local, -1).0.score;
                let p = LocalProfile::new(&q, &m, o, e).unwrap();
                let mut w = Workspace::default();
                assert_eq!(p.score(&t, &mut w), full, "len {len} o {o} e {e}");
                assert_eq!(p.scalar_score(&t, &mut w), full);
            }
        }
    }

    #[test]
    fn striped_matches_full_local_on_random_matrices() {
        // Random signed matrices over 20 letters and gap costs including extend > open: exercises lazy-F carries.
        let mut seed = 99u64;
        let mut next = |r: u64| {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) % r
        };
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        for round in 0..60 {
            let mut m = vec![0i32; MATRIX_LEN];
            for &x in letters {
                for &y in letters {
                    m[((x as usize) << 7) | y as usize] = next(15) as i32 - 6;
                }
            }
            let ql = 1 + next(70) as usize;
            let tl = next(90) as usize;
            let q: Vec<u8> = (0..ql).map(|_| letters[next(20) as usize]).collect();
            let t: Vec<u8> = (0..tl).map(|_| letters[next(20) as usize]).collect();
            for (o, e) in [(0, 0), (1, 3), (3, 1), (10, 1), (2, 2), (0, 4)] {
                let full = run(&q, &t, &m, o, e, Mode::Local, -1).0.score;
                let p = LocalProfile::new(&q, &m, o, e).unwrap();
                let mut w = Workspace::default();
                assert_eq!(p.score(&t, &mut w), full, "round {round} o {o} e {e}");
                assert_eq!(p.scalar_score(&t, &mut w), full);
            }
        }
    }

    #[test]
    fn striped_falls_back_on_saturation() {
        // 3000 identical residues × 20 = 60000 > i16::MAX, so the scalar path must take over.
        let m = simple_matrix(20, -5);
        let q = vec![b'A'; 3000];
        assert_eq!(local_score(&q, &q, &m, 5, 1).unwrap(), 60000);
    }

    #[test]
    fn batch_scores() {
        let m = simple_matrix(1, -1);
        let targets = b"ACGTTTACG";
        let offsets = [0, 4, 4, 9];
        let mut out = [0; 3];
        local_score_batch(b"ACG", targets, &offsets, &m, 2, 1, &mut out).unwrap();
        assert_eq!(out, [3, 0, 3]);
        assert_eq!(local_score_batch(b"ACG", targets, &[0, 10], &m, 2, 1, &mut out), Err(ERR_BAD_ARGUMENT));
    }
}
