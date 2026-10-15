//! Protein engine: mass, net charge, isoelectric point and Kyte–Doolittle hydropathy.
//! Constants and summation order match `ManagedKernels` exactly.

const fn mass_table(mono: bool) -> [f64; 256] {
    let mut t = [0.0f64; 256];
    // Residue (in-chain, water removed) masses: (average, monoisotopic).
    let r: [(u8, f64, f64); 26] = [
        (b'A', 71.0788, 71.03711),
        (b'R', 156.1875, 156.10111),
        (b'N', 114.1038, 114.04293),
        (b'D', 115.0886, 115.02694),
        (b'C', 103.1388, 103.00919),
        (b'E', 129.1155, 129.04259),
        (b'Q', 128.1307, 128.05858),
        (b'G', 57.0519, 57.02146),
        (b'H', 137.1411, 137.05891),
        (b'I', 113.1594, 113.08406),
        (b'L', 113.1594, 113.08406),
        (b'K', 128.1741, 128.09496),
        (b'M', 131.1926, 131.04049),
        (b'F', 147.1766, 147.06841),
        (b'P', 97.1167, 97.05276),
        (b'S', 87.0782, 87.03203),
        (b'T', 101.1051, 101.04768),
        (b'W', 186.2132, 186.07931),
        (b'Y', 163.1760, 163.06333),
        (b'V', 99.1326, 99.06841),
        (b'U', 150.0388, 150.95364),
        (b'O', 237.3018, 237.14773),
        (b'B', 114.5962, 114.53494),
        (b'Z', 128.6231, 128.55059),
        (b'J', 113.1594, 113.08406),
        (b'X', 110.0, 110.0),
    ];
    let mut k = 0;
    while k < r.len() {
        let m = if mono { r[k].2 } else { r[k].1 };
        t[r[k].0 as usize] = m;
        t[(r[k].0 | 0x20) as usize] = m;
        k += 1;
    }
    t
}

static AVERAGE: [f64; 256] = mass_table(false);
static MONO: [f64; 256] = mass_table(true);
const WATER_AVERAGE: f64 = 18.01528;
const WATER_MONO: f64 = 18.0105646863;

/// Peptide mass including one water; 0 for an empty sequence.
pub fn protein_mass(seq: &[u8], mono: bool) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let (table, mut m) = if mono { (&MONO, WATER_MONO) } else { (&AVERAGE, WATER_AVERAGE) };
    // Sequential sum on purpose: re-association would break bit-exact parity with the managed kernel.
    for &b in seq {
        m += table[b as usize];
    }
    m
}

// EMBOSS pKa set.
const PK_NTERM: f64 = 8.6;
const PK_CTERM: f64 = 3.6;
const PK_K: f64 = 10.8;
const PK_R: f64 = 12.5;
const PK_H: f64 = 6.5;
const PK_D: f64 = 3.9;
const PK_E: f64 = 4.1;
const PK_C: f64 = 8.5;
const PK_Y: f64 = 10.1;

fn folded_counts(seq: &[u8]) -> [i64; 256] {
    let mut h = [0i64; 256];
    for &b in seq {
        h[(b & 0xDF) as usize] += 1;
    }
    h
}

#[inline]
fn pos(pk: f64, ph: f64) -> f64 {
    1.0 / (1.0 + 10f64.powf(ph - pk))
}

#[inline]
fn neg(pk: f64, ph: f64) -> f64 {
    -1.0 / (1.0 + 10f64.powf(pk - ph))
}

fn charge_from_counts(h: &[i64; 256], ph: f64) -> f64 {
    let c = |b: u8| h[b as usize] as f64;
    pos(PK_NTERM, ph) + neg(PK_CTERM, ph)
        + c(b'K') * pos(PK_K, ph)
        + c(b'R') * pos(PK_R, ph)
        + c(b'H') * pos(PK_H, ph)
        + c(b'D') * neg(PK_D, ph)
        + c(b'E') * neg(PK_E, ph)
        + c(b'C') * neg(PK_C, ph)
        + c(b'Y') * neg(PK_Y, ph)
}

/// Net charge at `ph` (EMBOSS pKa set).
pub fn net_charge(seq: &[u8], ph: f64) -> f64 {
    charge_from_counts(&folded_counts(seq), ph)
}

/// Isoelectric point by bisection on [0, 14].
pub fn isoelectric_point(seq: &[u8]) -> f64 {
    let h = folded_counts(seq);
    let (mut lo, mut hi) = (0.0f64, 14.0f64);
    let mut i = 0;
    while i < 100 && hi - lo > 1e-4 {
        let mid = (lo + hi) / 2.0;
        if charge_from_counts(&h, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        i += 1;
    }
    (lo + hi) / 2.0
}

const fn kyte_doolittle() -> [f64; 256] {
    let mut t = [0.0f64; 256];
    let v: [(u8, f64); 20] = [
        (b'A', 1.8),
        (b'R', -4.5),
        (b'N', -3.5),
        (b'D', -3.5),
        (b'C', 2.5),
        (b'Q', -3.5),
        (b'E', -3.5),
        (b'G', -0.4),
        (b'H', -3.2),
        (b'I', 4.5),
        (b'L', 3.8),
        (b'K', -3.9),
        (b'M', 1.9),
        (b'F', 2.8),
        (b'P', -1.6),
        (b'S', -0.8),
        (b'T', -0.7),
        (b'W', -0.9),
        (b'Y', -1.3),
        (b'V', 4.2),
    ];
    let mut k = 0;
    while k < v.len() {
        t[v[k].0 as usize] = v[k].1;
        t[(v[k].0 | 0x20) as usize] = v[k].1;
        k += 1;
    }
    t
}

static KD: [f64; 256] = kyte_doolittle();

/// Sliding-window Kyte–Doolittle average; returns the number of values written.
pub fn hydropathy(seq: &[u8], window: usize, out: &mut [f64]) -> usize {
    if window == 0 || seq.len() < window {
        return 0;
    }
    let n = seq.len() - window + 1;
    let w = window as f64;
    let mut s = 0.0;
    for &b in &seq[..window] {
        s += KD[b as usize];
    }
    out[0] = s / w;
    for i in 1..n {
        s += KD[seq[i + window - 1] as usize] - KD[seq[i - 1] as usize];
        out[i] = s / w;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_of_glycine_dipeptide() {
        // GG: 2 × 57.0519 (Gly residue, average) + 18.01528 (water) = 132.11908 Da.
        assert!((protein_mass(b"GG", false) - 132.11908).abs() < 1e-9);
        assert!((protein_mass(b"gg", true) - (2.0 * 57.02146 + WATER_MONO)).abs() < 1e-12);
        assert_eq!(protein_mass(b"", false), 0.0);
    }

    #[test]
    fn charge_and_pi() {
        // Lysine-rich peptides are basic, acidic ones acidic.
        assert!(isoelectric_point(b"KKKKK") > 10.0);
        assert!(isoelectric_point(b"DDDEE") < 4.0);
        assert!(net_charge(b"KKK", 7.0) > 2.0);
        // At the pI the net charge is ~0.
        let pi = isoelectric_point(b"MKVLAAGIDE");
        assert!(net_charge(b"MKVLAAGIDE", pi).abs() < 1e-3);
    }

    #[test]
    fn hydropathy_window() {
        let mut out = [0.0; 3];
        assert_eq!(hydropathy(b"IVLA", 2, &mut out), 3);
        assert!((out[0] - (4.5 + 4.2) / 2.0).abs() < 1e-12);
        assert_eq!(hydropathy(b"IV", 3, &mut out), 0);
    }
}
