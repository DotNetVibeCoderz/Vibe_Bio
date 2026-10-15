//! Sequence engine: complement, translation, composition, k-mers and Hamming distance.
//!
//! Every function mirrors `ManagedKernels` in `src/Biotech.Net.Native` byte for byte. The AVX2 paths are
//! selected at run time and always fall back to the scalar loops for the tail and on other CPUs.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// IUPAC complement table, case-preserving; anything unrecognised maps to itself.
pub static COMPLEMENT: [u8; 256] = build_complement();

/// 2-bit code for A/C/G/T(U), case-insensitive; 255 for anything else.
pub static TWO_BIT: [u8; 256] = build_two_bit();

const fn build_complement() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    let pairs: [(u8, u8); 6] = [(b'A', b'T'), (b'C', b'G'), (b'R', b'Y'), (b'K', b'M'), (b'B', b'V'), (b'D', b'H')];
    let mut k = 0;
    while k < pairs.len() {
        let (a, b) = pairs[k];
        t[a as usize] = b;
        t[b as usize] = a;
        t[(a | 0x20) as usize] = b | 0x20;
        t[(b | 0x20) as usize] = a | 0x20;
        k += 1;
    }
    t[b'U' as usize] = b'A';
    t[b'u' as usize] = b'a';
    // S, W and N are self-complementary and already map to themselves.
    t
}

const fn build_two_bit() -> [u8; 256] {
    let mut t = [255u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

#[cfg(target_arch = "x86_64")]
#[inline]
fn has_avx2() -> bool {
    std::is_x86_feature_detected!("avx2")
}

/// Reverse complement of `src` into `dst` (same length).
pub fn reverse_complement(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    #[cfg(target_arch = "x86_64")]
    if src.len() >= 32 && has_avx2() {
        // SAFETY: AVX2 was detected at run time; the function bounds every load and store by the slice lengths.
        unsafe { complement_avx2::<true>(src, dst) };
        return;
    }
    let n = src.len();
    for (i, &b) in src.iter().enumerate() {
        dst[n - 1 - i] = COMPLEMENT[b as usize];
    }
}

/// Complement (not reversed) of `src` into `dst`.
pub fn complement(src: &[u8], dst: &mut [u8]) {
    debug_assert_eq!(src.len(), dst.len());
    #[cfg(target_arch = "x86_64")]
    if src.len() >= 32 && has_avx2() {
        // SAFETY: as above.
        unsafe { complement_avx2::<false>(src, dst) };
        return;
    }
    for (d, &s) in dst.iter_mut().zip(src) {
        *d = COMPLEMENT[s as usize];
    }
}

/// Low five bits of the complement of every byte `0x40 | k`. Upper and lower case letters share their low five
/// bits, so one 32-entry table (two `pshufb` lookups) covers both cases; bit 5 (case) and bit 6 are kept as is.
#[cfg(target_arch = "x86_64")]
const fn low5_table() -> [u8; 32] {
    let mut t = [0u8; 32];
    let mut k = 0;
    while k < 32 {
        t[k] = COMPLEMENT[0x40 | k] & 0x1F;
        k += 1;
    }
    t
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn complement_avx2<const REVERSE: bool>(src: &[u8], dst: &mut [u8]) {
    const T: [u8; 32] = low5_table();
    let n = src.len();
    let lo = _mm256_broadcastsi128_si256(_mm_loadu_si128(T.as_ptr() as *const __m128i));
    let hi = _mm256_broadcastsi128_si256(_mm_loadu_si128(T.as_ptr().add(16) as *const __m128i));
    let m0f = _mm256_set1_epi8(0x0F);
    let m10 = _mm256_set1_epi8(0x10);
    let me0 = _mm256_set1_epi8(0xE0u8 as i8);
    let mc0 = _mm256_set1_epi8(0xC0u8 as i8);
    let m40 = _mm256_set1_epi8(0x40);
    let rev = _mm256_setr_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
    );
    let mut i = 0;
    while i + 32 <= n {
        let v = _mm256_loadu_si256(src.as_ptr().add(i) as *const __m256i);
        let idx = _mm256_and_si256(v, m0f);
        let a = _mm256_shuffle_epi8(lo, idx);
        let b = _mm256_shuffle_epi8(hi, idx);
        let upper_half = _mm256_cmpeq_epi8(_mm256_and_si256(v, m10), m10);
        let low5 = _mm256_blendv_epi8(a, b, upper_half);
        // Only bytes 0x40..=0x7F (letters and a few punctuation marks, which the table maps to themselves) change.
        let letterish = _mm256_cmpeq_epi8(_mm256_and_si256(v, mc0), m40);
        let comp = _mm256_or_si256(_mm256_and_si256(v, me0), low5);
        let out = _mm256_blendv_epi8(v, comp, letterish);
        if REVERSE {
            let r = _mm256_permute4x64_epi64::<0x4E>(_mm256_shuffle_epi8(out, rev));
            _mm256_storeu_si256(dst.as_mut_ptr().add(n - i - 32) as *mut __m256i, r);
        } else {
            _mm256_storeu_si256(dst.as_mut_ptr().add(i) as *mut __m256i, out);
        }
        i += 32;
    }
    while i < n {
        let c = COMPLEMENT[src[i] as usize];
        if REVERSE {
            dst[n - 1 - i] = c;
        } else {
            dst[i] = c;
        }
        i += 1;
    }
}

/// Translates whole codons; codons with a non-ACGTU symbol become `X`. Returns the number of residues written.
pub fn translate(nuc: &[u8], table: &[u8; 64], protein: &mut [u8]) -> usize {
    let codons = nuc.len() / 3;
    for (p, c) in protein[..codons].iter_mut().zip(nuc.chunks_exact(3)) {
        let a = TWO_BIT[c[0] as usize];
        let b = TWO_BIT[c[1] as usize];
        let d = TWO_BIT[c[2] as usize];
        *p = if (a | b | d) > 3 { b'X' } else { table[((a as usize) << 4) | ((b as usize) << 2) | d as usize] };
    }
    codons
}

/// A, C, G, T/U, N and other counts (same layout as the managed `NucleotideCounts`).
pub fn count_nucleotides(seq: &[u8]) -> [i64; 6] {
    let mut acgtn = [0u64; 5];
    let mut i = 0;
    #[cfg(target_arch = "x86_64")]
    if seq.len() >= 32 && has_avx2() {
        // SAFETY: AVX2 detected; loads are bounded by the slice length.
        i = unsafe { count_avx2(seq, &mut acgtn) };
    }
    for &b in &seq[i..] {
        match b | 0x20 {
            b'a' => acgtn[0] += 1,
            b'c' => acgtn[1] += 1,
            b'g' => acgtn[2] += 1,
            b't' | b'u' => acgtn[3] += 1,
            b'n' => acgtn[4] += 1,
            _ => {}
        }
    }
    let [a, c, g, t, n] = acgtn.map(|x| x as i64);
    [a, c, g, t, n, seq.len() as i64 - a - c - g - t - n]
}

/// Processes whole 32-byte blocks and returns how many bytes were consumed.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn count_avx2(seq: &[u8], out: &mut [u64; 5]) -> usize {
    // `b | 0x20` folds case; only 'A'/'a' fold to 'a', and so on, so the comparisons stay exact.
    let fold = _mm256_set1_epi8(0x20);
    let targets = [b'a', b'c', b'g', b't', b'n'].map(|c| _mm256_set1_epi8(c as i8));
    let tu = _mm256_set1_epi8(b'u' as i8);
    let zero = _mm256_setzero_si256();
    let blocks = seq.len() / 32;
    let mut done = 0;
    while done < blocks {
        // Byte counters overflow after 255 increments, so flush at most every 255 blocks.
        let chunk = (blocks - done).min(255);
        let mut acc = [zero; 5];
        for k in 0..chunk {
            let v = _mm256_or_si256(_mm256_loadu_si256(seq.as_ptr().add((done + k) * 32) as *const __m256i), fold);
            for s in 0..5 {
                let mut eq = _mm256_cmpeq_epi8(v, targets[s]);
                if s == 3 {
                    eq = _mm256_or_si256(eq, _mm256_cmpeq_epi8(v, tu));
                }
                acc[s] = _mm256_sub_epi8(acc[s], eq);
            }
        }
        for s in 0..5 {
            let sad = _mm256_sad_epu8(acc[s], zero);
            let mut lanes = [0u64; 4];
            _mm256_storeu_si256(lanes.as_mut_ptr() as *mut __m256i, sad);
            out[s] += lanes.iter().sum::<u64>();
        }
        done += chunk;
    }
    blocks * 32
}

/// Adds the occurrences of every byte value to `counts`.
pub fn byte_histogram(data: &[u8], counts: &mut [i64; 256]) {
    // Four interleaved histograms break the store-to-load dependency on runs of the same byte.
    let mut h = vec![0u64; 1024];
    let mut chunks = data.chunks_exact(4);
    for c in &mut chunks {
        h[c[0] as usize] += 1;
        h[256 + c[1] as usize] += 1;
        h[512 + c[2] as usize] += 1;
        h[768 + c[3] as usize] += 1;
    }
    for &b in chunks.remainder() {
        h[b as usize] += 1;
    }
    for k in 0..256 {
        counts[k] += (h[k] + h[256 + k] + h[512 + k] + h[768 + k]) as i64;
    }
}

/// Adds every k-mer (1 ≤ k ≤ 12) to `counts` (4^k entries); windows with non-ACGT symbols are skipped.
pub fn count_kmers(seq: &[u8], k: u32, counts: &mut [i32]) -> i64 {
    let mask: u32 = (1u32 << (2 * k)) - 1;
    let mut kmer = 0u32;
    let mut valid = 0u32;
    let mut total = 0i64;
    for &b in seq {
        let c = TWO_BIT[b as usize];
        if c > 3 {
            valid = 0;
            kmer = 0;
            continue;
        }
        kmer = ((kmer << 2) | c as u32) & mask;
        valid += 1;
        if valid >= k {
            counts[kmer as usize] += 1;
            total += 1;
        }
    }
    total
}

/// Number of byte positions at which `a` and `b` differ (equal lengths).
pub fn hamming(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    let mut d = 0usize;
    #[cfg(target_arch = "x86_64")]
    if n >= 32 && has_avx2() {
        // SAFETY: AVX2 detected; loads bounded by n.
        unsafe {
            let mut same = 0usize;
            while i + 32 <= n {
                let va = _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i);
                let vb = _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i);
                same += (_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb)) as u32).count_ones() as usize;
                i += 32;
            }
            d = i - same;
        }
    }
    for k in i..n {
        d += (a[k] != b[k]) as usize;
    }
    d
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar_rc(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| COMPLEMENT[b as usize]).collect()
    }

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 33
    }

    #[test]
    fn complement_table_is_iupac() {
        assert_eq!(COMPLEMENT[b'A' as usize], b'T');
        assert_eq!(COMPLEMENT[b'g' as usize], b'c');
        assert_eq!(COMPLEMENT[b'R' as usize], b'Y');
        assert_eq!(COMPLEMENT[b'b' as usize], b'v');
        assert_eq!(COMPLEMENT[b'U' as usize], b'A');
        assert_eq!(COMPLEMENT[b'N' as usize], b'N');
        assert_eq!(COMPLEMENT[b'-' as usize], b'-');
    }

    #[test]
    fn reverse_complement_matches_scalar_for_every_byte_and_length() {
        let all: Vec<u8> = (0..=255u8).collect();
        for len in [0usize, 1, 31, 32, 33, 64, 100, 256] {
            let src: Vec<u8> = all.iter().cycle().skip(len * 7).take(len).copied().collect();
            let mut dst = vec![0u8; len];
            reverse_complement(&src, &mut dst);
            assert_eq!(dst, scalar_rc(&src), "len {len}");
            let mut c = vec![0u8; len];
            complement(&src, &mut c);
            let expect: Vec<u8> = src.iter().map(|&b| COMPLEMENT[b as usize]).collect();
            assert_eq!(c, expect);
        }
        let mut d = [0u8; 8];
        reverse_complement(b"ACGTacgu", &mut d);
        assert_eq!(&d, b"acgtACGT");
    }

    #[test]
    fn translate_standard_code() {
        // Standard code in kernel order (A=0, C=1, G=2, T=3).
        let table: [u8; 64] = *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
        let mut p = [0u8; 4];
        assert_eq!(translate(b"ATGGCCtgaNNNAT", &table, &mut p), 4);
        assert_eq!(&p, b"MA*X");
    }

    #[test]
    fn nucleotide_counts_simd_matches_scalar() {
        let mut seed = 7u64;
        let alphabet = b"ACGTUNacgtunRYx-";
        let s: Vec<u8> = (0..10_000).map(|_| alphabet[(lcg(&mut seed) % 16) as usize]).collect();
        let mut expect = [0i64; 6];
        for &b in &s {
            match b {
                b'A' | b'a' => expect[0] += 1,
                b'C' | b'c' => expect[1] += 1,
                b'G' | b'g' => expect[2] += 1,
                b'T' | b't' | b'U' | b'u' => expect[3] += 1,
                b'N' | b'n' => expect[4] += 1,
                _ => expect[5] += 1,
            }
        }
        assert_eq!(count_nucleotides(&s), expect);
        assert_eq!(count_nucleotides(b""), [0; 6]);
    }

    #[test]
    fn histogram_and_kmers() {
        let mut h = [0i64; 256];
        byte_histogram(b"AACGTT", &mut h);
        byte_histogram(b"A", &mut h);
        assert_eq!((h[b'A' as usize], h[b'T' as usize], h[b'C' as usize]), (3, 2, 1));
        let mut counts = vec![0i32; 16];
        // AC, CG, GT; the N breaks the window; then AA.
        assert_eq!(count_kmers(b"ACGTNAA", 2, &mut counts), 4);
        assert_eq!(counts[0b0001], 1); // AC
        assert_eq!(counts[0b0000], 1); // AA
        let mut c12 = vec![0i32; 1 << 24];
        assert_eq!(count_kmers(b"ACGTACGTACGTA", 12, &mut c12), 2);
    }

    #[test]
    fn hamming_simd_matches_scalar() {
        let mut seed = 3u64;
        let a: Vec<u8> = (0..1001).map(|_| b"ACGT"[(lcg(&mut seed) % 4) as usize]).collect();
        let b: Vec<u8> = (0..1001).map(|_| b"ACGT"[(lcg(&mut seed) % 4) as usize]).collect();
        let expect = a.iter().zip(&b).filter(|(x, y)| x != y).count();
        assert_eq!(hamming(&a, &b), expect);
        assert_eq!(hamming(b"", b""), 0);
        assert_eq!(hamming(b"ACGT", b"ACGA"), 1);
    }
}
