//! Exact rank over F_p (A7: only valid when every coboundary entry is an integer).
//! A one-sided exact cross-check (for integer delta, rank over F_p <= rank over R, so betti_p >= betti_R) and a way to exhibit field dependence (e.g. RP^2: H^1 = 0 over R, 1 over F_2).
use faer::Mat;

/// Largest |entry| accepted as an exact integer: 2^53 (every integer up to it is exactly representable in f64).
pub const MAX_EXACT_INT: f64 = 9007199254740992.0;

/// Integer matrix of `a`, or None if any entry is not an integer within ±2^53. ADR-017 round 2 (P-NUM-04): NaN passed
/// the old `(v - round v).abs() > 1e-12` test (comparison with NaN is false) and became 0; ±inf became i64::MAX; 2^63
/// saturated to i64::MAX (odd), so F_2 ranks were computed from a different matrix.
pub fn to_int(a: &Mat<f64>) -> Option<Vec<Vec<i64>>> {
    let mut out = vec![vec![0i64; a.ncols()]; a.nrows()];
    for i in 0..a.nrows() { for j in 0..a.ncols() {
        let v = a[(i, j)];
        if !v.is_finite() || v.abs() > MAX_EXACT_INT || !((v - v.round()).abs() <= 1e-12) { return None; }
        out[i][j] = v.round() as i64;
    }}
    Some(out)
}

/// (a·b) mod p without overflow: i128 product (ADR-017 round 2, P-NUM-03: the i64 product overflowed for p > ~3.04e9 and
/// silently flipped ranks, e.g. Betti [0,0] instead of [1,1] at p = 4294967291). Inputs are in [0, p).
#[inline] fn mulmod(a: i64, b: i64, p: i64) -> i64 { ((a as i128 * b as i128) % p as i128) as i64 }

/// Deterministic Miller–Rabin for 64-bit (bases = first 12 primes), i128 arithmetic.
pub fn is_prime(p: i64) -> bool {
    if p < 2 { return false; }
    for q in [2i64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if p % q == 0 { return p == q; } }
    let (mut d, mut s) = (p - 1, 0); while d % 2 == 0 { d /= 2; s += 1; }
    let pw = |mut b: i64, mut e: i64| { let mut r = 1i64; b %= p; while e > 0 { if e & 1 == 1 { r = mulmod(r, b, p); } b = mulmod(b, b, p); e >>= 1; } r };
    'w: for a in [2i64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pw(a, d); if x == 1 || x == p - 1 { continue; }
        for _ in 1..s { x = mulmod(x, x, p); if x == p - 1 { continue 'w; } }
        return false;
    }
    true
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FpError { NotPrime(i64), NotIntegral }
impl std::fmt::Display for FpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { FpError::NotPrime(p) => write!(f, "F_p rank needs a prime p >= 2, got {p}"), FpError::NotIntegral => write!(f, "coboundary has a non-integer, non-finite or |x| > 2^53 entry") }
    }
}

/// Rank over F_p. `p` must be prime (callers: `betti_mod_p_checked` verifies it).
pub fn rank_mod_p(m: &[Vec<i64>], ncols: usize, p: i64) -> usize {
    let mut a: Vec<Vec<i64>> = m.iter().map(|r| r.iter().map(|x| x.rem_euclid(p)).collect()).collect();
    let inv = |x: i64| { let (mut r, mut b, mut e) = (1i64, x, p - 2); while e > 0 { if e & 1 == 1 { r = mulmod(r, b, p); } b = mulmod(b, b, p); e >>= 1; } r };
    let mut rank = 0;
    for c in 0..ncols {
        let Some(piv) = (rank..a.len()).find(|&r| a[r][c] != 0) else { continue };
        a.swap(rank, piv);
        let iv = inv(a[rank][c]);
        for x in a[rank].iter_mut() { *x = mulmod(*x, iv, p); }
        for r in 0..a.len() {
            if r != rank && a[r][c] != 0 {
                let f = a[r][c];
                for j in 0..ncols { a[r][j] = (a[r][j] - mulmod(f, a[rank][j], p)).rem_euclid(p); }
            }
        }
        rank += 1;
    }
    rank
}

/// Betti numbers over F_p for a complex whose coboundaries are integral. None if not integral OR p is not prime
/// (use `betti_mod_p_checked` for the reason).
pub fn betti_mod_p(delta: &[Mat<f64>], dims: &[usize], p: i64) -> Option<Vec<usize>> { betti_mod_p_checked(delta, dims, p).ok() }

/// Betti numbers over F_p with a typed error: p must be prime (any prime up to i64::MAX; products are i128).
pub fn betti_mod_p_checked(delta: &[Mat<f64>], dims: &[usize], p: i64) -> Result<Vec<usize>, FpError> {
    if !is_prime(p) { return Err(FpError::NotPrime(p)); }
    let mut ranks = Vec::new();
    for d in delta { let m = to_int(d).ok_or(FpError::NotIntegral)?; ranks.push(rank_mod_p(&m, d.ncols(), p)); }
    Ok((0..dims.len()).map(|k| {
        let rin = if k == 0 { 0 } else { ranks[k - 1] };
        let rout = if k < ranks.len() { ranks[k] } else { 0 };
        dims[k] - rin - rout
    }).collect())
}

/// Typed wrapper: Betti numbers in exact F_p arithmetic (`ExactFp`, never confused with the f64 results).
pub fn betti_fp(delta: &[Mat<f64>], dims: &[usize], p: u32) -> Option<sheaf_cell::num::family::ExactFp<Vec<usize>>> {
    betti_mod_p(delta, dims, p as i64).map(|b| sheaf_cell::num::family::ExactFp::new(b, p))
}

#[cfg(test)]
mod mulmod_tests {
    use super::mulmod;
    use proptest::prelude::*;

    fn reference(a: i64, b: i64, p: i64) -> i64 { ((a as u128 * b as u128) % p as u128) as i64 }

    #[test]
    fn mulmod_edge_cases() {
        for p in [2i64, 3, 4294967291, 9223372036854775783, i64::MAX] {
            for a in [0, 1, p / 2, p - 2, p - 1] { for b in [0, 1, p / 2, p - 2, p - 1] {
                if a < 0 || b < 0 { continue; }
                assert_eq!(mulmod(a, b, p), reference(a, b, p), "a={a} b={b} p={p}");
            }}
        }
        // (p-1)^2 = 1 mod p
        assert_eq!(mulmod(4294967290, 4294967290, 4294967291), 1);
        assert_eq!(mulmod(9223372036854775782, 9223372036854775782, 9223372036854775783), 1);
    }

    fn triple() -> impl Strategy<Value = (i64, i64, i64)> {
        prop_oneof![2i64..=i64::MAX, 2i64..=1 << 32, Just(4294967291i64), Just(9223372036854775783i64), Just(i64::MAX)]
            .prop_flat_map(|p| (0..p, 0..p, Just(p)))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(4096))]
        #[test]
        fn mulmod_matches_u128((a, b, p) in triple()) {
            let r = mulmod(a, b, p);
            prop_assert!((0..p).contains(&r));
            prop_assert_eq!(r, reference(a, b, p));
        }
    }
}

/// FIX-BEND-SHEAF-001: Bend's sheaf-is-prime-001 only covered p < 2^16. These check `is_prime` against an independent
/// trial-division reference near 2^32 (exhaustive windows and a proptest), plus known 64-bit primes/composites.
#[cfg(test)]
mod is_prime_tests {
    use super::is_prime;
    use proptest::prelude::*;

    /// Independent reference: trial division by 2 and odd d with d*d <= n (u64, no shared code with Miller–Rabin).
    fn reference(n: i64) -> bool {
        if n < 2 { return false; }
        let n = n as u64;
        if n % 2 == 0 { return n == 2; }
        let mut d = 3u64;
        while d * d <= n { if n % d == 0 { return false; } d += 2; }
        true
    }

    #[test]
    fn is_prime_exhaustive_windows_near_2_32() {
        let t: i64 = 1 << 32;
        for c in [t, (1i64 << 31), (1i64 << 33)] {
            for n in (c - 2000)..=(c + 2000) { assert_eq!(is_prime(n), reference(n), "n={n}"); }
        }
    }

    #[test]
    fn is_prime_known_values() {
        // primes: largest below 2^32, smallest above 2^32, Bend's pmax, largest below 2^63
        for p in [4294967291i64, 4294967311, 2147483647, 9223372036854775783] { assert!(is_prime(p), "{p}"); }
        // composites: 2^32±1, 2^32-1 = 3·5·17·257·65537, F5 = 641·6700417, 65521·65537, strong pseudoprimes
        // to small bases (3215031751 = 151·751·28351, 3825123056546413051)
        for n in [4294967296i64, 4294967295, 4294967297, 65521 * 65537, 3215031751, 3825123056546413051,
                  4294967291 * 2, 4294967291 * 3] {
            assert!(!is_prime(n), "{n}");
        }
        for n in [i64::MIN, -7, -1, 0, 1] { assert!(!is_prime(n), "{n}"); }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(2048))]
        #[test]
        fn is_prime_matches_trial_division_near_2_32(n in (1i64 << 32) - (1 << 24)..=(1i64 << 32) + (1 << 24)) {
            prop_assert_eq!(is_prime(n), reference(n));
        }
    }
}
