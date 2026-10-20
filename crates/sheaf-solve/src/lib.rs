//! sheaf-solve: harmonic extension for cellular sheaves (graphs or any δ⁰). Depends on sheaf-cell only, never on
//! sheaf-cohom (review §5): whether the extension is UNIQUE is a cohomology question (H⁰(X,B;F) = 0, see
//! `sheaf_cohom::relative::h0_relative`), composed by the caller (sheaf-cli `--boundary`), not decided here.
//!
//! Problem. δ: C⁰ → C¹ is the sheaf coboundary (sparse, m × n, any per-edge restriction maps), L = δᵀδ, and the
//! sheaf energy is E(x) = ‖δx‖² = xᵀLx. Given fixed values x_B on a set B of boundary dofs, the harmonic extension is
//! the interior x_I minimising E, i.e. the solution of  L_II x_I = −L_IB x_B  (∂E/∂x_I = 2(Lx)_I = 0).
//! If L_II is singular (an interior part with sections not pinned by the boundary) the system is still consistent;
//! CG started from `x0_I` returns x0_I's component in ker L_II plus the minimum-norm particular solution.
//! With B = ∅ this is Lx = 0 from x0, so the result is the orthogonal projection of x0 onto H⁰ = ker L (global sections).
//!
//! Family `Exact` (ADR-013): a linear solve of the exact energy, not an ordered/interval approximation.
//! Certificate (certified residual, ADR-008): relative residual ‖L_II x_I + L_IB x_B‖ / max(‖L_IB x_B‖,
//! ‖L‖_bound·max(‖x‖, ‖x_start‖), tiny), the energy, iteration count, the route taken (dense Cholesky or CG), and
//! `converged` = residual ≤ tol.
pub use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::family::Exact;
use sheaf_cell::num::krylov::{cg, SolveTol};
use sheaf_cell::num::norm as nrm;
pub use sheaf_cell::ops::energy;

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Route { Dense, Cg, Trivial }

#[derive(Debug, Clone)]
pub struct Options { pub rtol: f64, pub max_iter: usize, pub dense_max: usize }
impl Default for Options { fn default() -> Self { Options { rtol: 1e-12, max_iter: 20_000, dense_max: 400 } } }

#[derive(Debug, Clone)]
pub struct Certificate {
    pub method: Exact,
    pub route: Route,
    /// ‖L_II x_I + L_IB x_B‖ / scale (see module docs)
    pub residual_rel: f64,
    pub residual_abs: f64,
    pub iters: usize,
    pub converged: bool,
    /// CG met a direction with pᵀL_II p ≈ 0 (singular interior); the solution is still valid if `converged`
    pub breakdown: bool,
    pub energy: f64,
}

#[derive(Debug, Clone)]
pub struct HarmonicExtension { pub x: Vec<f64>, pub cert: Certificate }

#[derive(Debug, Clone, PartialEq)]
pub enum HarmonicError { BoundaryIndex(usize), DuplicateBoundary(usize), NonFinite, X0Len { expected: usize, got: usize } }

/// L x = δᵀ(δx) on the full dof space.
pub fn apply_l(delta: &Csr, x: &[f64]) -> Vec<f64> { delta.mul_t(&delta.mul(x)) }
/// Gershgorin-type bound on ‖L‖₂ via ‖δ‖₁‖δ‖∞ ≥ ‖δ‖₂².
fn l_bound(delta: &Csr) -> f64 {
    let rinf = (0..delta.nrows).map(|r| (delta.ptr[r]..delta.ptr[r + 1]).map(|k| delta.val[k].abs()).sum::<f64>()).fold(0.0, f64::max);
    let mut c1 = vec![0.0; delta.ncols];
    for k in 0..delta.val.len() { c1[delta.idx[k]] += delta.val[k].abs(); }
    rinf * c1.into_iter().fold(0.0, f64::max)
}

/// Dof-level boundary from node blocks: `offsets[i]..offsets[i+1]` are node i's dofs.
pub fn node_boundary(offsets: &[usize], nodes: &[(usize, Vec<f64>)]) -> Vec<(usize, f64)> {
    nodes.iter().flat_map(|(i, v)| { assert_eq!(v.len(), offsets[i + 1] - offsets[*i], "boundary value length != stalk dim of node {i}");
        v.iter().enumerate().map(move |(a, &val)| (offsets[*i] + a, val)) }).collect()
}

/// Harmonic extension. `boundary`: (dof, value) pairs; `x0`: optional full-length start (interior part used; default 0).
pub fn harmonic_extension(delta: &Csr, boundary: &[(usize, f64)], x0: Option<&[f64]>, o: &Options) -> Result<HarmonicExtension, HarmonicError> {
    let n = delta.ncols;
    let mut is_b = vec![false; n];
    let mut x = match x0 { Some(v) if v.len() != n => return Err(HarmonicError::X0Len { expected: n, got: v.len() }), Some(v) => v.to_vec(), None => vec![0.0; n] };
    for &(i, v) in boundary {
        if i >= n { return Err(HarmonicError::BoundaryIndex(i)); }
        if is_b[i] { return Err(HarmonicError::DuplicateBoundary(i)); }
        is_b[i] = true; x[i] = v;
    }
    if x.iter().any(|v| !v.is_finite()) || delta.val.iter().any(|v| !v.is_finite()) { return Err(HarmonicError::NonFinite); }
    let interior: Vec<usize> = (0..n).filter(|&i| !is_b[i]).collect();
    let x_start_norm = nrm(&x);
    let ni = interior.len();
    let restrict = |full: &[f64]| -> Vec<f64> { interior.iter().map(|&i| full[i]).collect() };
    let embed = |xi: &[f64], base: &mut Vec<f64>| { for (k, &i) in interior.iter().enumerate() { base[i] = xi[k]; } };
    // b = −L_IB x_B
    let mut xb = vec![0.0; n]; for &(i, v) in boundary { xb[i] = v; }
    let b: Vec<f64> = restrict(&apply_l(delta, &xb)).into_iter().map(|v| -v).collect();
    let lii = |v: &[f64]| -> Vec<f64> { let mut f = vec![0.0; n]; embed(v, &mut f); restrict(&apply_l(delta, &f)) };
    let lb = l_bound(delta);
    let mut route = Route::Trivial; let (mut iters, mut breakdown) = (0usize, false);
    if ni > 0 {
        let mut done = false;
        if ni <= o.dense_max {
            // dense L_II, Cholesky; falls back to CG if L_II is not numerically positive definite
            let mut a = vec![0.0; ni * ni];
            for j in 0..ni { let mut e = vec![0.0; ni]; e[j] = 1.0; let c = lii(&e); for i in 0..ni { a[i * ni + j] = c[i]; } }
            if let Some(xi) = cholesky_solve(&mut a, ni, &b, 1e-12 * lb.max(1e-300)) { embed(&xi, &mut x); route = Route::Dense; done = true; }
        }
        if !done {
            route = Route::Cg;
            let xi0 = restrict(&x);
            let scale = nrm(&b).max(lb * nrm(&xi0)).max(1e-300);
            let out = cg(&lii, &b, xi0, scale, lb, SolveTol { rtol: o.rtol, max_iter: o.max_iter });
            iters = out.iters; breakdown = out.breakdown; let xi = out.x;
            embed(&xi, &mut x);
        }
    }
    let lx = apply_l(delta, &x);
    let res_abs = nrm(&restrict(&lx));
    // scale includes the start: with B = ∅ the answer can be ~0 (no sections) while the start was not
    let scale = nrm(&b).max(lb * nrm(&x)).max(lb * x_start_norm).max(1e-300);
    let residual_rel = if ni == 0 { 0.0 } else { res_abs / scale };
    let converged = residual_rel.is_finite() && residual_rel <= o.rtol.max(1e-10);
    let e = energy(delta, &x);
    Ok(HarmonicExtension { x, cert: Certificate { method: Exact, route, residual_rel, residual_abs: res_abs, iters, converged, breakdown, energy: e } })
}

/// In-place Cholesky of a symmetric n×n (row-major) and solve; None if a pivot ≤ tol (not numerically SPD).
fn cholesky_solve(a: &mut [f64], n: usize, b: &[f64], tol: f64) -> Option<Vec<f64>> {
    for j in 0..n {
        let mut d = a[j * n + j];
        for k in 0..j { d -= a[j * n + k] * a[j * n + k]; }
        if !(d > tol) { return None; }
        let d = d.sqrt(); a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for k in 0..j { s -= a[i * n + k] * a[j * n + k]; }
            a[i * n + j] = s / d;
        }
    }
    let mut y = b.to_vec();
    for i in 0..n { for k in 0..i { y[i] -= a[i * n + k] * y[k]; } y[i] /= a[i * n + i]; }
    for i in (0..n).rev() { for k in i + 1..n { y[i] -= a[k * n + i] * y[k]; } y[i] /= a[i * n + i]; }
    Some(y)
}
