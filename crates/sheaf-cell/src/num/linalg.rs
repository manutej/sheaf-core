//! Dense numerical linear algebra over R (faer SVD). Ranks are *numerical*: A6.
//! A6 rank(A) = #{ sigma_i > tol }, tol = max(m,n) * eps * sigma_max (numpy/LAPACK default),
//!    and every rank is reported together with the singular-value gap around the cut.
use faer::Mat;

#[derive(Debug, Clone, PartialEq)]
pub struct RankReport {
    pub rank: usize,
    pub tol: f64,
    pub sigma_max: f64,
    /// smallest singular value kept (None if rank 0)
    pub sigma_kept_min: Option<f64>,
    /// largest singular value dropped (None if full rank)
    pub sigma_dropped_max: Option<f64>,
}
impl RankReport {
    /// sigma_kept_min / max(sigma_dropped_max, tiny): large => the rank decision is well separated.
    pub fn gap_ratio(&self) -> f64 {
        match (self.sigma_kept_min, self.sigma_dropped_max) {
            (Some(k), Some(d)) => k / d.max(f64::MIN_POSITIVE),
            _ => f64::INFINITY,
        }
    }
}

pub struct FullSvd { pub u: Mat<f64>, pub s: Vec<f64>, pub v: Mat<f64> }

pub fn svd_full(a: &Mat<f64>) -> FullSvd {
    let (m, n) = (a.nrows(), a.ncols());
    if m == 0 || n == 0 {
        return FullSvd { u: Mat::identity(m, m), s: vec![], v: Mat::identity(n, n) };
    }
    let svd = faer::linalg::solvers::Svd::new(a.as_ref()).expect("svd converged");
    let k = m.min(n);
    let s: Vec<f64> = (0..k).map(|i| svd.S().column_vector()[i]).collect();
    FullSvd { u: svd.U().to_owned(), s, v: svd.V().to_owned() }
}

pub fn rank_report(s: &[f64], m: usize, n: usize, tol_override: Option<f64>) -> RankReport {
    let smax = s.iter().cloned().fold(0.0, f64::max);
    let tol = tol_override.unwrap_or((m.max(n) as f64) * f64::EPSILON * smax);
    let rank = s.iter().filter(|&&x| x > tol).count();
    let kept = s.iter().cloned().filter(|&x| x > tol).fold(None, |a: Option<f64>, x| Some(a.map_or(x, |y| y.min(x))));
    let dropped = s.iter().cloned().filter(|&x| x <= tol).fold(None, |a: Option<f64>, x| Some(a.map_or(x, |y| y.max(x))));
    RankReport { rank, tol, sigma_max: smax, sigma_kept_min: kept, sigma_dropped_max: dropped }
}

pub fn cols(a: &Mat<f64>, start: usize, end: usize) -> Mat<f64> {
    Mat::from_fn(a.nrows(), end - start, |i, j| a[(i, start + j)])
}
pub fn transpose(a: &Mat<f64>) -> Mat<f64> { a.transpose().to_owned() }
pub fn vstack(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    assert_eq!(a.ncols(), b.ncols());
    Mat::from_fn(a.nrows() + b.nrows(), a.ncols(), |i, j| if i < a.nrows() { a[(i, j)] } else { b[(i - a.nrows(), j)] })
}
pub fn fro(a: &Mat<f64>) -> f64 { a.norm_l2() }
pub fn matvec(a: &Mat<f64>, x: &[f64]) -> Vec<f64> {
    (0..a.nrows()).map(|i| (0..a.ncols()).map(|j| a[(i, j)] * x[j]).sum()).collect()
}
pub fn matvec_t(a: &Mat<f64>, y: &[f64]) -> Vec<f64> {
    (0..a.ncols()).map(|j| (0..a.nrows()).map(|i| a[(i, j)] * y[i]).sum()).collect()
}
/// orthogonal projection of x onto span of the orthonormal columns of q
pub fn project(q: &Mat<f64>, x: &[f64]) -> Vec<f64> { let c = matvec_t(q, x); matvec(q, &c) }
pub fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() }
pub fn norm(a: &[f64]) -> f64 { dot(a, a).sqrt() }
