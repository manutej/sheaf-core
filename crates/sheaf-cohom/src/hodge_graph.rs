//! Scale path for the v1 target (graph H^1 / harmonic part of an edge cochain), matrix-free.
//! exact = d0 x*, x* = argmin ||d0 x - w||  via CG on L0 x = d0^T w  (L0 PSD, rhs in range => CG converges);
//! harmonic = w - exact (graph: no faces, so no coexact part). Certificate: ||d0^T harmonic|| / ||w||.
use sheaf_cell::num::sparse::Csr;

pub struct GraphHodge { pub exact: Vec<f64>, pub harmonic: Vec<f64>, pub potential: Vec<f64>, pub iters: usize, pub cert_coclosed: f64,
    /// true only if the *certificate* ||d0^T h||/||w|| <= cert_tol (not just the CG residual) and no breakdown
    pub converged: bool, pub breakdown: bool }

pub fn graph_hodge_cg(d0: &Csr, w: &[f64], rtol: f64, max_iter: usize) -> GraphHodge { graph_hodge_cg_cert(d0, w, rtol, max_iter, 1e-8) }

/// Grok r1 #1/#2: curvature guard (p^T L p <= eps ||p||^2 => breakdown, stop) and a certificate-based `converged`.
pub fn graph_hodge_cg_cert(d0: &Csr, w: &[f64], rtol: f64, max_iter: usize, cert_tol: f64) -> GraphHodge {
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let b = d0.mul_t(w);
    let bn = dot(&b, &b).sqrt();
    let mut x = vec![0.0; d0.ncols];
    let mut r = b.clone(); let mut p = r.clone();
    let mut rr = dot(&r, &r); let mut it = 0; let mut breakdown = false;
    while it < max_iter && rr.sqrt() > rtol * bn.max(1e-300) {
        let ap = d0.mul_t(&d0.mul(&p));
        let pap = dot(&p, &ap);
        if !(pap > 1e-14 * dot(&p, &p)) { breakdown = true; break; } // direction (numerically) in ker L0
        let alpha = rr / pap;
        for i in 0..x.len() { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
        let rr2 = dot(&r, &r);
        let beta = rr2 / rr; rr = rr2;
        for i in 0..p.len() { p[i] = r[i] + beta * p[i]; }
        it += 1;
    }
    let exact = d0.mul(&x);
    let harmonic: Vec<f64> = w.iter().zip(&exact).map(|(a, b)| a - b).collect();
    let wn = dot(w, w).sqrt().max(1e-300);
    let c = d0.mul_t(&harmonic);
    let cert_coclosed = dot(&c, &c).sqrt() / wn;
    let converged = cert_coclosed.is_finite() && cert_coclosed <= cert_tol && !breakdown;
    GraphHodge { exact, harmonic, potential: x, iters: it, cert_coclosed, converged, breakdown }
}
