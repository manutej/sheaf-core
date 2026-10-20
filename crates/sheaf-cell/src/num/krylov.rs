//! Conjugate gradients on a symmetric positive SEMI-definite operator, matrix-free, deterministic.
//! On a consistent system started from x0, the iterates stay in x0 + K(A, r0); from x0 = 0 that is the minimum-norm
//! solution (no component in ker A). Breakdown guard: pᵀAp ≤ 1e-14·scale·‖p‖² stops (direction numerically in ker A).
//! (sheaf-core does not use faer's matrix-free CG: in faer 0.24.4 `InitialGuessStatus::Zero` misreports.)
use super::{dot, norm};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveTol { pub rtol: f64, pub max_iter: usize }
impl Default for SolveTol { fn default() -> Self { SolveTol { rtol: 1e-12, max_iter: 20_000 } } }

#[derive(Debug, Clone)]
pub struct CgOutcome { pub x: Vec<f64>, pub iters: usize, pub breakdown: bool, pub residual_norm: f64 }

/// Solve A x = b from `x0`; stop when ‖r‖ ≤ rtol·scale. `op_scale` (e.g. ‖A‖_bound) only enters the breakdown guard.
pub fn cg(apply: &dyn Fn(&[f64]) -> Vec<f64>, b: &[f64], x0: Vec<f64>, scale: f64, op_scale: f64, tol: SolveTol) -> CgOutcome {
    let n = b.len();
    let mut x = x0;
    let ax = apply(&x);
    let mut r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
    let mut p = r.clone(); let mut rr = dot(&r, &r);
    let (mut iters, mut breakdown) = (0usize, false);
    while iters < tol.max_iter && rr.sqrt() > tol.rtol * scale {
        let ap = apply(&p); let pap = dot(&p, &ap);
        if !(pap > 1e-14 * op_scale * dot(&p, &p)) { breakdown = true; break; }
        let al = rr / pap;
        for k in 0..n { x[k] += al * p[k]; r[k] -= al * ap[k]; }
        let rr2 = dot(&r, &r); let be = rr2 / rr; rr = rr2;
        for k in 0..n { p[k] = r[k] + be * p[k]; }
        iters += 1;
    }
    let residual_norm = norm(&r);
    CgOutcome { x, iters, breakdown, residual_norm }
}
