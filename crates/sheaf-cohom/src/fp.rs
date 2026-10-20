//! Exact rank over F_p (A7: only valid when every coboundary entry is an integer).
//! A one-sided exact cross-check (for integer delta, rank over F_p <= rank over R, so betti_p >= betti_R) and a way to exhibit field dependence (e.g. RP^2: H^1 = 0 over R, 1 over F_2).
use faer::Mat;

pub fn to_int(a: &Mat<f64>) -> Option<Vec<Vec<i64>>> {
    let mut out = vec![vec![0i64; a.ncols()]; a.nrows()];
    for i in 0..a.nrows() { for j in 0..a.ncols() {
        let v = a[(i, j)];
        if (v - v.round()).abs() > 1e-12 { return None; }
        out[i][j] = v.round() as i64;
    }}
    Some(out)
}

pub fn rank_mod_p(m: &[Vec<i64>], ncols: usize, p: i64) -> usize {
    let mut a: Vec<Vec<i64>> = m.iter().map(|r| r.iter().map(|x| x.rem_euclid(p)).collect()).collect();
    let inv = |x: i64| { let (mut r, mut b, mut e) = (1i64, x, p - 2); while e > 0 { if e & 1 == 1 { r = r * b % p; } b = b * b % p; e >>= 1; } r };
    let mut rank = 0;
    for c in 0..ncols {
        let Some(piv) = (rank..a.len()).find(|&r| a[r][c] != 0) else { continue };
        a.swap(rank, piv);
        let iv = inv(a[rank][c]);
        for x in a[rank].iter_mut() { *x = *x * iv % p; }
        for r in 0..a.len() {
            if r != rank && a[r][c] != 0 {
                let f = a[r][c];
                for j in 0..ncols { a[r][j] = (a[r][j] - f * a[rank][j]).rem_euclid(p); }
            }
        }
        rank += 1;
    }
    rank
}

/// Betti numbers over F_p for a complex whose coboundaries are integral. None if not integral.
pub fn betti_mod_p(delta: &[Mat<f64>], dims: &[usize], p: i64) -> Option<Vec<usize>> {
    let mut ranks = Vec::new();
    for d in delta { let m = to_int(d)?; ranks.push(rank_mod_p(&m, d.ncols(), p)); }
    Some((0..dims.len()).map(|k| {
        let rin = if k == 0 { 0 } else { ranks[k - 1] };
        let rout = if k < ranks.len() { ranks[k] } else { 0 };
        dims[k] - rin - rout
    }).collect())
}

/// Typed wrapper: Betti numbers in exact F_p arithmetic (`ExactFp`, never confused with the f64 results).
pub fn betti_fp(delta: &[Mat<f64>], dims: &[usize], p: u32) -> Option<sheaf_cell::num::family::ExactFp<Vec<usize>>> {
    betti_mod_p(delta, dims, p as i64).map(|b| sheaf_cell::num::family::ExactFp::new(b, p))
}
