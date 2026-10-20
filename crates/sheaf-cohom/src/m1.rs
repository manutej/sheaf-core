//! H⁰ for the m1 contract (SPEC v1.18 R-H0-08 / R-PERF-18) on faer: `global_sections` (orthonormal basis of ker δ⁰ by
//! per-component SVD) and `h0_dim_from_l` (dim ker L₀ by per-component symmetric eigenvalues). Port of the Python
//! reference `sheafcore.global_sections` / `h0_dim_from_L` (ADR-016).
//!
//! * Zero cut: the single `H0Policy` type, in its contract-named global variants (`H0Policy::m1_svd(1e-9)` =
//!   σ > 1e-9·max(1, σ_max over all components); `H0Policy::m1_eig(1e-10)` = λ ≤ 1e-10·max(1, max|λ|)).
//! * Components: bipartite row/column NONZERO pattern of δ (stored zeros are not couplings), ordered exactly like the
//!   reference (scipy `connected_components` labels in node order with the m row nodes first: components by their
//!   smallest row, row-less columns afterwards by index), so the basis columns come in the reference's block order.
//! * Caps (R-PERF-18) are checked for EVERY component before anything is densified (pass 0), the basis size before it
//!   is allocated; one dense block is alive at a time. No iterative solver is used (no CG / LSMR): dense SVD / EVD only.
//! * The basis of a multi-dimensional kernel is an orthonormal basis of the same subspace as numpy's, not the same
//!   vectors (singular vectors are unique only up to rotation within a repeated σ); parity is on the projector GGᵀ.
use faer::{Mat, Par};
use faer::dyn_stack::{MemBuffer, MemStack};

/// ADR-017 (adversarial panel, P-CONC-01): every m1 dense factorisation runs SEQUENTIALLY. faer's default global
/// parallelism splits the SVD/EVD work by rayon thread count, so the same input gave a different kernel basis G
/// (max |ΔG| 1.66, same projector) and bitwise-different σ/λ on 1 vs 2 vs 8 threads, i.e. a rank that could flip at the
/// R-H0-08 cut depending on the machine. Explicit `Par::Seq` also never reads or mutates faer's global setting.
/// Measured cost on a 4000-column component: none for the SVD (113 s either way), 25 s → 39 s for the eigenvalues.
pub const M1_PAR: Par = Par::Seq;

fn sv_seq(a: faer::MatRef<'_, f64>) -> Result<Vec<f64>, faer::linalg::svd::SvdError> {
    use faer::linalg::svd::{svd, svd_scratch, ComputeSvdVectors};
    let mut s = faer::diag::Diag::<f64>::zeros(a.nrows().min(a.ncols()));
    let mut buf = MemBuffer::new(svd_scratch::<f64>(a.nrows(), a.ncols(), ComputeSvdVectors::No, ComputeSvdVectors::No, M1_PAR, Default::default()));
    svd(a, s.as_mut(), None, None, M1_PAR, MemStack::new(&mut buf), Default::default())?;
    let mut v: Vec<f64> = s.column_vector().iter().copied().collect(); v.sort_by(|x, y| y.total_cmp(x)); Ok(v)
}
/// Eigenvalues ascending, LOWER triangle only (numpy `eigvalsh` UPLO='L'), sequential.
fn eig_seq(a: faer::MatRef<'_, f64>) -> Result<Vec<f64>, faer::linalg::evd::EvdError> {
    use faer::linalg::evd::{self_adjoint_evd, self_adjoint_evd_scratch, ComputeEigenvectors};
    let mut s = faer::diag::Diag::<f64>::zeros(a.nrows());
    let mut buf = MemBuffer::new(self_adjoint_evd_scratch::<f64>(a.nrows(), ComputeEigenvectors::No, M1_PAR, Default::default()));
    self_adjoint_evd(a, s.as_mut(), None, M1_PAR, MemStack::new(&mut buf), Default::default())?;
    Ok(s.column_vector().iter().copied().collect())
}
use sheaf_cell::m1::{adjacency_components, pattern_components};
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::{H0Policy, ZeroCut};

pub const MAX_DENSE_DIM: usize = 5000;
pub const MAX_DENSE_BYTES: u128 = 2 * (1u128 << 30);

#[derive(Debug, Clone, PartialEq)]
pub enum M1H0Error { Perf(String), Policy(String), Numeric(String) }
impl std::fmt::Display for M1H0Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { M1H0Error::Perf(s) => write!(f, "R-PERF-18: {s}"), M1H0Error::Policy(s) => write!(f, "R-H0-08: {s}"),
                     M1H0Error::Numeric(s) => write!(f, "R-FINITE-01: {s}") }
    }
}

/// One coupled block of δ: the rows holding a nonzero in it and its columns, both ascending.
#[derive(Debug, Clone, PartialEq)]
pub struct M1Comp { pub rows: Vec<usize>, pub cols: Vec<usize> }

/// Upper bound on the dense SVD footprint (padded block + U + Vt + workspace), the reference's `_svd_bytes`.
pub fn svd_bytes(m: usize, n: usize) -> u128 { let (m, n) = (m as u128, n as u128); 8 * (2 * m.max(n) * n + 6 * n * n) }

/// Components of δ in the reference's order.
pub fn components(d: &Csr) -> Vec<M1Comp> {
    if d.ncols == 0 { return vec![]; }
    let mut g: Vec<M1Comp> = pattern_components(d).into_iter().map(|(rows, cols)| M1Comp { rows, cols }).collect();
    g.sort_by_key(|c| c.rows.first().copied().unwrap_or(d.nrows + c.cols[0]));
    g
}

/// Pass 0: every component's caps, before any densification (messages = the reference's).
pub fn check_caps(comps: &[M1Comp]) -> Result<(), M1H0Error> {
    for c in comps {
        if c.cols.len() > MAX_DENSE_DIM { return Err(M1H0Error::Perf(format!("component of dim {} exceeds MAX_DENSE_DIM={MAX_DENSE_DIM}", c.cols.len()))); }
        if svd_bytes(c.rows.len(), c.cols.len()) > MAX_DENSE_BYTES { return Err(M1H0Error::Perf(format!("dense SVD of {}x{} block exceeds MAX_DENSE_BYTES", c.rows.len(), c.cols.len()))); }
    }
    Ok(())
}

/// Dense rows×cols block of δ for one component (stored duplicates summed, like `toarray`). Column-major (faer).
pub fn dense_block(d: &Csr, c: &M1Comp) -> Mat<f64> {
    let mut pos = std::collections::HashMap::with_capacity(c.cols.len());
    for (k, &j) in c.cols.iter().enumerate() { pos.insert(j, k); }
    let mut a = Mat::<f64>::zeros(c.rows.len(), c.cols.len());
    for (ri, &i) in c.rows.iter().enumerate() {
        for p in d.ptr[i]..d.ptr[i + 1] { if let Some(&k) = pos.get(&d.idx[p]) { a[(ri, k)] += d.val[p]; } }
    }
    a
}

/// Singular values of one component's block, descending, as the reference sees them: a wide block (rows < cols) is
/// zero-padded to cols × cols there, so the list always has `cols` entries (the padding contributes exact zeros, which
/// matter for the rank count when the caller's tol is negative). Empty for a row-less component (the reference skips it).
pub fn block_singular_values(d: &Csr, c: &M1Comp) -> Result<Vec<f64>, M1H0Error> {
    if c.rows.is_empty() { return Ok(vec![]); }
    let mut s = sv_seq(dense_block(d, c).as_ref()).map_err(|e| M1H0Error::Numeric(format!("SVD did not converge: {e:?}")))?;
    s.resize(c.cols.len().max(s.len()), 0.0);
    Ok(s)
}

fn need(policy: H0Policy, svd: bool) -> Result<(), M1H0Error> {
    match (policy.cut, svd) {
        (ZeroCut::GlobalSigma { .. }, true) | (ZeroCut::GlobalLambda { .. }, false) => Ok(()),
        _ => Err(M1H0Error::Policy(format!("m1 {} needs the contract policy {}, got {}", if svd { "global_sections" } else { "h0_dim_from_L" },
                                            if svd { "m1_global_svd" } else { "m1_global_eig" }, policy.name()))),
    }
}

/// Pass 1: per-component ranks under `policy` (GlobalSigma), plus the global σ_max. One dense block alive at a time.
pub fn ranks(d: &Csr, comps: &[M1Comp], policy: H0Policy) -> Result<(Vec<usize>, f64), M1H0Error> {
    need(policy, true)?;
    let mut svs = Vec::with_capacity(comps.len());
    for c in comps { svs.push(block_singular_values(d, c)?); }
    let smax = svs.iter().flatten().cloned().fold(0.0, f64::max);
    let thr = policy.global_cut(smax).expect("global policy");
    Ok((svs.iter().map(|s| s.iter().filter(|&&x| x > thr).count()).collect(), smax))
}

/// Orthonormal basis (cols × (cols − r)) of the component's kernel: the last right singular vectors of the block;
/// identity when the component has no rows. Only V is computed (no U).
pub fn block_kernel(d: &Csr, c: &M1Comp, r: usize) -> Result<Mat<f64>, M1H0Error> {
    let n = c.cols.len();
    if c.rows.is_empty() { return Ok(Mat::identity(n, n)); }
    use faer::linalg::svd::{svd, svd_scratch, ComputeSvdVectors};
    let a = dense_block(d, c); let m = a.nrows();
    let cv = if m >= n { ComputeSvdVectors::Thin } else { ComputeSvdVectors::Full };   // V is n×n either way
    let par = M1_PAR;
    let mut buf = MemBuffer::new(svd_scratch::<f64>(m, n, ComputeSvdVectors::No, cv, par, Default::default()));
    let mut s = faer::diag::Diag::<f64>::zeros(m.min(n)); let mut v = Mat::<f64>::zeros(n, n);
    svd(a.as_ref(), s.as_mut(), None, Some(v.as_mut()), par, MemStack::new(&mut buf), Default::default())
        .map_err(|e| M1H0Error::Numeric(format!("SVD did not converge: {e:?}")))?;
    drop(a);
    Ok(v.subcols(r.min(n), n - r.min(n)).to_owned())
}

/// Orthonormal basis G (n × N) of ker δ under the m1 policy (`H0Policy::m1_svd`). Column blocks follow `components`.
pub fn global_sections(d: &Csr, policy: H0Policy) -> Result<Mat<f64>, M1H0Error> {
    need(policy, true)?;
    let n = d.ncols;
    if n == 0 { return Ok(Mat::zeros(0, 0)); }
    let comps = components(d); check_caps(&comps)?;
    let (rk, _) = ranks(d, &comps, policy)?;
    let nn: usize = comps.iter().zip(&rk).map(|(c, r)| c.cols.len() - r).sum();
    if 8 * n as u128 * nn as u128 > MAX_DENSE_BYTES { return Err(M1H0Error::Perf(format!("dense basis {n}x{nn} exceeds MAX_DENSE_BYTES"))); }
    let mut g = Mat::<f64>::zeros(n, nn); let mut k = 0;
    for (c, &r) in comps.iter().zip(&rk) {
        let kb = block_kernel(d, c, r)?;
        for (a, &row) in c.cols.iter().enumerate() { for b in 0..kb.ncols() { g[(row, k + b)] = kb[(a, b)]; } }
        k += kb.ncols();
    }
    Ok(g)
}

/// Components of L (undirected nonzero adjacency), any order (only the count matters).
pub fn l_components(l: &Csr) -> Vec<Vec<usize>> { if l.nrows == 0 { vec![] } else { adjacency_components(l) } }

pub fn check_eig_caps(groups: &[Vec<usize>]) -> Result<(), M1H0Error> {
    for g in groups {
        if g.len() > MAX_DENSE_DIM { return Err(M1H0Error::Perf(format!("component of dim {} exceeds MAX_DENSE_DIM={MAX_DENSE_DIM}", g.len()))); }
        if 24 * (g.len() as u128).pow(2) > MAX_DENSE_BYTES { return Err(M1H0Error::Perf(format!("dense eig of dim {} exceeds MAX_DENSE_BYTES", g.len()))); }
    }
    Ok(())
}

/// Eigenvalues (ascending) of the principal block L[cols, cols], read from the LOWER triangle (numpy `eigvalsh`'s
/// default UPLO='L', so a non-symmetric input is treated exactly like the reference does).
pub fn block_eigenvalues(l: &Csr, cols: &[usize]) -> Result<Vec<f64>, M1H0Error> {
    let mut pos = std::collections::HashMap::with_capacity(cols.len());
    for (k, &j) in cols.iter().enumerate() { pos.insert(j, k); }
    let mut a = Mat::<f64>::zeros(cols.len(), cols.len());
    for (ri, &i) in cols.iter().enumerate() { for p in l.ptr[i]..l.ptr[i + 1] { if let Some(&k) = pos.get(&l.idx[p]) { a[(ri, k)] += l.val[p]; } } }
    eig_seq(a.as_ref()).map_err(|e| M1H0Error::Numeric(format!("eigensolver did not converge: {e:?}")))
}

/// dim ker L₀ under the m1 policy (`H0Policy::m1_eig`): λ zero iff λ ≤ rel·max(1, max|λ| over all components).
pub fn h0_dim_from_l(l: &Csr, policy: H0Policy) -> Result<usize, M1H0Error> {
    need(policy, false)?;
    if l.nrows == 0 { return Ok(0); }
    let groups = l_components(l); check_eig_caps(&groups)?;
    let mut ws = Vec::with_capacity(l.nrows);
    for g in &groups { ws.extend(block_eigenvalues(l, g)?); }
    let lmax = ws.iter().fold(0.0f64, |a, &w| a.max(w.abs()));
    let thr = policy.global_cut(lmax).expect("global policy");
    Ok(ws.iter().filter(|&&w| w <= thr).count())
}
