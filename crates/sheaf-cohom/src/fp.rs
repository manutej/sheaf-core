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
