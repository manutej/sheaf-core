//! Scalable spectral decomposition of sheaf Hodge Laplacians (v1).
//!
//! Smallest-k eigenpairs of L_k = d_{k-1} d_{k-1}^T + d_k^T d_k, matrix-free, on top of faer 0.24's
//! thick-restart Lanczos (`matrix_free::eigen::partial_self_adjoint_eigen`), which only returns the
//! LARGEST-magnitude eigenpairs. We therefore run it on the spectral transform
//!     B = P (c I - L) P,   c = 1.01 * ||L||_bound  (c > lambda_max, so B is PSD and its top = L's bottom),
//! where P projects out (a) a *known* kernel (component-constant sections, exact when every edge of a
//! component has F_u == F_v injective) and (b) eigenvectors locked in earlier rounds.
//! Single-vector Krylov cannot resolve a repeated eigenvalue in one run, so we lock and restart until a
//! fresh run finds nothing below the current (k+1)-th value, then do one Rayleigh-Ritz on the locked span.
//!
//! Every returned pair carries its residual r = ||L v - lambda v|| (v unit). For symmetric L, Weyl/Bauer-Fike
//! gives a true eigenvalue in [lambda - r, lambda + r]. What is NOT certified: that no eigenvalue below the
//! reported ones was missed (Krylov completeness). Wording per ADR-008: "certified residual", never "proof".
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::matrix_free::eigen::{partial_self_adjoint_eigen, partial_self_adjoint_eigen_scratch, PartialEigenParams};
use faer::matrix_free::LinOp;
use faer::{Mat, MatMut, MatRef, Par, Side};

/// The zero cut is the single H0 policy (`sheaf_cell::num::tol::H0Policy`, ADR-012/015). `smallest*` apply it to the operator
/// they are given, so they implement the per-component rule when handed ONE component (the v1 driver does that).
pub use sheaf_cell::num::tol::{H0Policy, EIG_TOL};

pub use sheaf_cell::num::{dot, norm};
pub use sheaf_cell::ops::HodgeOp;

/// Known kernel: per-component constant sections. comp[node] = Some(component id) if the node's
/// component is known-constant; dofs of node i are off[i]..off[i+1].
#[derive(Debug, Clone, Default)]
pub struct KnownKernel { pub comp_of_dof: Vec<Option<(u32, u32)>>, pub n_groups: usize, pub group_size: Vec<u32> }
impl KnownKernel {
    pub fn none(n: usize) -> Self { KnownKernel { comp_of_dof: vec![None; n], n_groups: 0, group_size: vec![] } }
    /// groups = (component, coordinate) pairs; group g's kernel vector = indicator / sqrt(size)
    pub fn from_components(node_comp: &[Option<u32>], off: &[usize]) -> Self {
        let mut map = std::collections::BTreeMap::new();
        let n = *off.last().unwrap_or(&0);
        let mut comp_of_dof = vec![None; n];
        for (i, c) in node_comp.iter().enumerate() {
            if let Some(c) = c { for (j, dof) in (off[i]..off[i + 1]).enumerate() { comp_of_dof[dof] = Some((*c, j as u32)); *map.entry((*c, j as u32)).or_insert(0u32) += 1; } }
        }
        let ids: std::collections::BTreeMap<(u32, u32), u32> = map.keys().enumerate().map(|(g, k)| (*k, g as u32)).collect();
        let group_size: Vec<u32> = map.values().cloned().collect();
        let comp_of_dof = comp_of_dof.into_iter().map(|x| x.map(|k| (ids[&k], 0))).collect();
        KnownKernel { comp_of_dof, n_groups: ids.len(), group_size }
    }
    pub fn dim(&self) -> usize { self.n_groups }
    fn project(&self, x: &mut [f64]) {
        if self.n_groups == 0 { return; }
        let mut s = vec![0.0; self.n_groups];
        for (i, g) in self.comp_of_dof.iter().enumerate() { if let Some((g, _)) = g { s[*g as usize] += x[i]; } }
        for (i, g) in self.comp_of_dof.iter().enumerate() { if let Some((g, _)) = g { x[i] -= s[*g as usize] / self.group_size[*g as usize] as f64; } }
    }
    /// orthonormal kernel basis as dense columns (small cases / tests only)
    pub fn basis_vec(&self, g: usize) -> Vec<f64> {
        let f = 1.0 / (self.group_size[g] as f64).sqrt();
        self.comp_of_dof.iter().map(|x| if matches!(x, Some((h, _)) if *h as usize == g) { f } else { 0.0 }).collect()
    }
}

#[derive(Debug)]
struct Shifted<'a> { op: &'a HodgeOp, c: f64, known: &'a KnownKernel, locked: &'a [Vec<f64>] }
impl<'a> Shifted<'a> {
    fn proj(&self, x: &mut [f64]) {
        for _ in 0..2 { // twice: classical Gram-Schmidt re-orthogonalisation
            self.known.project(x);
            for q in self.locked { let a = dot(q, x); x.iter_mut().zip(q).for_each(|(v, b)| *v -= a * b); }
        }
    }
}
impl<'a> LinOp<f64> for Shifted<'a> {
    fn apply_scratch(&self, _rhs_ncols: usize, _par: Par) -> StackReq { StackReq::EMPTY }
    fn nrows(&self) -> usize { self.op.n }
    fn ncols(&self) -> usize { self.op.n }
    fn apply(&self, out: MatMut<'_, f64>, rhs: MatRef<'_, f64>, _par: Par, _stack: &mut MemStack) {
        let mut out = out; let n = self.op.n;
        let mut x = vec![0.0; n]; let mut y = vec![0.0; n];
        for j in 0..rhs.ncols() {
            for i in 0..n { x[i] = rhs[(i, j)]; }
            self.proj(&mut x);
            self.op.apply(&x, &mut y);
            for i in 0..n { y[i] = self.c * x[i] - y[i]; }
            self.proj(&mut y);
            for i in 0..n { out[(i, j)] = y[i]; }
        }
    }
    fn conj_apply(&self, out: MatMut<'_, f64>, rhs: MatRef<'_, f64>, par: Par, stack: &mut MemStack) { self.apply(out, rhs, par, stack) }
}

/// B = P (L + sigma I)^{-1} P, inner solve = Jacobi-preconditioned CG from x0 = 0 (deterministic).
#[derive(Debug)]
struct ShiftInvert<'a> { op: &'a HodgeOp, sigma: f64, diag_inv: &'a [f64], known: &'a KnownKernel, locked: &'a [Vec<f64>], inner_rtol: f64, inner_max: usize,
    stats: std::sync::Mutex<(usize, usize, f64)> /* (solves, total iters, worst rel residual) */ }
impl<'a> ShiftInvert<'a> {
    fn proj(&self, x: &mut [f64]) {
        for _ in 0..2 { self.known.project(x); for q in self.locked { let a = dot(q, x); x.iter_mut().zip(q).for_each(|(v, b)| *v -= a * b); } }
    }
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.op.n; let bn = norm(b).max(1e-300);
        let mut x = vec![0.0; n]; let mut r = b.to_vec();
        let mut z: Vec<f64> = r.iter().zip(self.diag_inv).map(|(a, d)| a * d).collect();
        let mut p = z.clone(); let mut rz = dot(&r, &z); let mut ap = vec![0.0; n]; let mut it = 0;
        while it < self.inner_max && norm(&r) > self.inner_rtol * bn {
            self.op.apply(&p, &mut ap); for i in 0..n { ap[i] += self.sigma * p[i]; }
            let pap = dot(&p, &ap); if !(pap > 0.0) { break; }
            let al = rz / pap;
            for i in 0..n { x[i] += al * p[i]; r[i] -= al * ap[i]; }
            for i in 0..n { z[i] = r[i] * self.diag_inv[i]; }
            let rz2 = dot(&r, &z); let be = rz2 / rz; rz = rz2;
            for i in 0..n { p[i] = z[i] + be * p[i]; }
            it += 1;
        }
        let mut st = self.stats.lock().unwrap(); st.0 += 1; st.1 += it; st.2 = st.2.max(norm(&r) / bn);
        x
    }
}
impl<'a> LinOp<f64> for ShiftInvert<'a> {
    fn apply_scratch(&self, _rhs_ncols: usize, _par: Par) -> StackReq { StackReq::EMPTY }
    fn nrows(&self) -> usize { self.op.n }
    fn ncols(&self) -> usize { self.op.n }
    fn apply(&self, out: MatMut<'_, f64>, rhs: MatRef<'_, f64>, _par: Par, _stack: &mut MemStack) {
        let mut out = out; let n = self.op.n; let mut x = vec![0.0; n];
        for j in 0..rhs.ncols() {
            for i in 0..n { x[i] = rhs[(i, j)]; }
            self.proj(&mut x);
            let mut y = self.solve(&x);
            self.proj(&mut y);
            for i in 0..n { out[(i, j)] = y[i]; }
        }
    }
    fn conj_apply(&self, out: MatMut<'_, f64>, rhs: MatRef<'_, f64>, par: Par, stack: &mut MemStack) { self.apply(out, rhs, par, stack) }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method { ShiftInvert, Polynomial }

#[derive(Debug, Clone)]
pub struct EigPair { pub value: f64, pub residual: f64 }

#[derive(Debug, Clone)]
pub struct SpectralOpts { pub method: Method, pub sigma_rel: f64, pub inner_rtol: f64, pub inner_max: usize, pub k: usize, pub rtol: f64, pub h0: H0Policy, pub max_rounds: usize, pub max_restarts: usize, pub seed: u64, pub dense_max: usize }
impl Default for SpectralOpts {
    fn default() -> Self { SpectralOpts { method: Method::ShiftInvert, sigma_rel: 1e-6, inner_rtol: 1e-12, inner_max: 20_000, k: 8, rtol: 1e-11, h0: H0Policy::adr012(), max_rounds: 8, max_restarts: 2000, seed: 0x9E3779B97F4A7C15, dense_max: 1500 } }
}

#[derive(Debug, Clone)]
pub struct Spectral {
    pub method: &'static str,
    pub n: usize,
    /// dim of the deflated known kernel (exact component-constant sections)
    pub known_kernel_dim: usize,
    /// the k+1 smallest eigenpairs of L on the complement of the known kernel, ascending
    pub pairs: Vec<EigPair>,
    pub vectors: Vec<Vec<f64>>,
    pub c_shift: f64,
    pub lambda_max_bound: f64,
    pub zero_tol: f64,
    /// eigenvalues <= zero_tol among `pairs`
    pub found_zero: usize,
    pub rounds: usize,
    /// shift-invert only: (inner solves, total inner CG iterations, worst inner relative residual)
    pub inner: Option<(usize, usize, f64)>,
    pub sigma: f64,
    /// true if every one of the k+1 requested pairs converged (faer reported) in the last round
    pub converged: bool,
    /// true if all k+1 computed pairs were <= zero_tol: H0 is then only a LOWER bound
    pub h0_lower_bound_only: bool,
}
impl Spectral {
    pub fn h0_dim(&self) -> usize { self.known_kernel_dim + self.found_zero }
    /// largest eigenvalue counted zero (with residual) and smallest counted nonzero (with residual)
    pub fn cut(&self) -> (Option<&EigPair>, Option<&EigPair>) {
        (self.pairs[..self.found_zero].last(), self.pairs.get(self.found_zero))
    }
    /// The zero cut is consistent with the certified residuals (Weyl): (lambda_zero_max + r) <= tol < (lambda_nonzero_min - r).
    /// A policy check, not a kernel certificate (Krylov completeness is not certified).
    pub fn cut_consistent_with_residuals(&self) -> bool {
        let (z, nz) = self.cut();
        z.map_or(true, |p| p.value + p.residual <= self.zero_tol) && nz.map_or(!self.h0_lower_bound_only, |p| p.value - p.residual > self.zero_tol)
    }
    /// Diagnostic flag, NOT a second count (ADR-015): some pair counted zero by the policy is numerically resolvable
    /// from 0, i.e. lambda > max(64 eps ||L||_bound, 10 * residual). Such a pair is a "near section".
    pub fn has_near_sections(&self) -> bool {
        self.pairs[..self.found_zero].iter().any(|p| p.value > (64.0 * f64::EPSILON * self.lambda_max_bound).max(10.0 * p.residual))
    }
    /// spectral gap = smallest eigenvalue above the zero cut
    pub fn spectral_gap(&self) -> Option<f64> { self.pairs.get(self.found_zero).map(|p| p.value) }
}

struct Rng(u64);
impl Rng { fn next(&mut self) -> f64 { let mut s = self.0; s ^= s << 13; s ^= s >> 7; s ^= s << 17; self.0 = s; ((s >> 11) as f64) / ((1u64 << 53) as f64) - 0.5 } }

fn finish(op: &HodgeOp, known: &KnownKernel, mut pairs: Vec<(f64, Vec<f64>)>, k: usize, c: f64, method: &'static str, rounds: usize, converged: bool, h0: H0Policy) -> Spectral {
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    pairs.truncate(k + 1);
    let lam_bound = op.norm_bound();
    let zero_tol = h0.zero_tol(lam_bound);
    let ps: Vec<EigPair> = pairs.iter().map(|(l, v)| EigPair { value: *l, residual: op.residual(v, *l) }).collect();
    let found_zero = ps.iter().filter(|p| p.value <= zero_tol).count();
    let h0_lower_bound_only = found_zero == ps.len() && ps.len() < op.n - known.dim();
    Spectral { method, n: op.n, known_kernel_dim: known.dim(), vectors: pairs.into_iter().map(|p| p.1).collect(), pairs: ps, c_shift: c, lambda_max_bound: lam_bound,
               zero_tol, found_zero, rounds, inner: None, sigma: 0.0, converged, h0_lower_bound_only }
}

/// Dense reference: full EVD of L, known kernel removed by projection (P L P), k+1 smallest of the rest.
pub fn smallest_dense(op: &HodgeOp, known: &KnownKernel, k: usize, h0: H0Policy) -> Spectral {
    let n = op.n;
    let mut l = op.dense();
    // project known kernel: L restricted to complement. Kernel vectors of L are eigvecs with eigenvalue 0;
    // replacing them by eigenvalue +inf-like (shift) keeps the rest unchanged: L + s * K K^T.
    let s = 10.0 * op.norm_bound().max(1.0) + 10.0;
    for g in 0..known.dim() { let q = known.basis_vec(g); for i in 0..n { if q[i] != 0.0 { for j in 0..n { if q[j] != 0.0 { l[(i, j)] += s * q[i] * q[j]; } } } } }
    let evd = l.self_adjoint_eigen(Side::Lower).expect("evd");
    let m = (k + 1).min(n - known.dim());
    let pairs: Vec<(f64, Vec<f64>)> = (0..m).map(|j| (evd.S().column_vector()[j], (0..n).map(|i| evd.U()[(i, j)]).collect())).collect();
    finish(op, known, pairs, k, 0.0, "dense-evd", 1, true, h0)
}

/// Matrix-free smallest k+1 eigenpairs (see module doc). Falls back to dense when n <= dense_max
/// (faer's restarted Lanczos needs n > max_dim >= 64).
pub fn smallest(op: &HodgeOp, known: &KnownKernel, o: &SpectralOpts) -> Spectral {
    let n = op.n; let want = o.k + 1;
    if n <= o.dense_max || n <= 2 * want.max(64) + known.dim() + 2 { return smallest_dense(op, known, o.k, o.h0); }
    let bound = op.norm_bound();
    let c = 1.01 * bound + 1e-12;
    let sigma = o.sigma_rel * bound.max(1.0);
    let diag_inv: Vec<f64> = if o.method == Method::ShiftInvert {
        let d = op.diagonal();
        d.iter().map(|v| 1.0 / (v + sigma)).collect()
    } else { vec![] };
    let mut inner_tot = (0usize, 0usize, 0.0f64);
    let mut locked: Vec<Vec<f64>> = vec![]; let mut lvals: Vec<f64> = vec![];
    let mut rng = Rng(o.seed | 1); let mut rounds = 0; let mut converged = true;
    for round in 0..o.max_rounds {
        rounds = round + 1;
        let nev = want.min(n - known.dim() - locked.len());
        if nev == 0 { break; }
        let a = Shifted { op, c, known, locked: &locked };
        let si = ShiftInvert { op, sigma, diag_inv: &diag_inv, known, locked: &locked, inner_rtol: o.inner_rtol, inner_max: o.inner_max, stats: std::sync::Mutex::new((0, 0, 0.0)) };
        let (lin, scale): (&dyn LinOp<f64>, f64) = if o.method == Method::ShiftInvert { (&si, 1.0 / sigma) } else { (&a, c) };
        let params = PartialEigenParams { max_dim: (3 * nev).max(64), max_restarts: o.max_restarts, ..Default::default() };
        let mut v0 = faer::Col::<f64>::zeros(n);
        { let mut x: Vec<f64> = (0..n).map(|_| rng.next()).collect(); a.proj(&mut x); for i in 0..n { v0[i] = x[i]; } }
        let mut vecs = Mat::<f64>::zeros(n, nev); let mut vals = vec![0.0; nev];
        let mut mem = MemBuffer::new(partial_self_adjoint_eigen_scratch(lin, nev, Par::Seq, params));
        let info = partial_self_adjoint_eigen(vecs.as_mut(), &mut vals, lin, v0.as_ref(), o.rtol * scale, Par::Seq, MemStack::new(&mut mem), params);
        let nc = info.n_converged_eigen;
        { let st = si.stats.lock().unwrap(); inner_tot.0 += st.0; inner_tot.1 += st.1; inner_tot.2 = inner_tot.2.max(st.2); }
        converged = nc >= nev;
        let thr = { let mut s = lvals.clone(); s.sort_by(|a, b| a.partial_cmp(b).unwrap()); s.get(want - 1).cloned().unwrap_or(f64::INFINITY) };
        let mut new_below = 0; let mut fresh: Vec<(f64, Vec<f64>)> = vec![];
        for j in 0..nc {
            let mut v: Vec<f64> = (0..n).map(|i| vecs[(i, j)]).collect();
            a.proj(&mut v); let nv = norm(&v); if nv < 1e-8 { continue; }
            v.iter_mut().for_each(|x| *x /= nv);
            let lam = if o.method == Method::ShiftInvert { 1.0 / vals[j] - sigma } else { c - vals[j] };
            // Grok v1 #2: a further copy AT the current (k+1)-th value is new too (it is orthogonal to the locked span)
            if lam <= thr * (1.0 + 1e-8) + 1e-14 { new_below += 1; }
            fresh.push((lam, v));
        }
        drop(a); drop(si);
        for (lam, v) in fresh { locked.push(v); lvals.push(lam); }
        if round > 0 && new_below == 0 { break; }
        if nc == 0 { break; }
    }
    // Rayleigh-Ritz on span(locked): orthonormalise (MGS x2), H = V^T L V, EVD.
    let mut basis: Vec<Vec<f64>> = vec![];
    for mut v in locked.into_iter() {
        for _ in 0..2 { known.project(&mut v); for q in &basis { let a = dot(q, &v); v.iter_mut().zip(q).for_each(|(x, b)| *x -= a * b); } }
        let nv = norm(&v); if nv > 1e-8 { v.iter_mut().for_each(|x| *x /= nv); basis.push(v); }
    }
    let p = basis.len();
    let lv: Vec<Vec<f64>> = basis.iter().map(|v| { let mut y = vec![0.0; n]; op.apply(v, &mut y); y }).collect();
    let h = Mat::<f64>::from_fn(p, p, |i, j| 0.5 * (dot(&basis[i], &lv[j]) + dot(&basis[j], &lv[i])));
    let evd = h.self_adjoint_eigen(Side::Lower).expect("ritz evd");
    let pairs: Vec<(f64, Vec<f64>)> = (0..p).map(|j| {
        let mut u = vec![0.0; n];
        for (i, b) in basis.iter().enumerate() { let y = evd.U()[(i, j)]; if y != 0.0 { u.iter_mut().zip(b).for_each(|(o, x)| *o += y * x); } }
        let nu = norm(&u); u.iter_mut().for_each(|x| *x /= nu);
        (evd.S().column_vector()[j], u)
    }).collect();
    let name = if o.method == Method::ShiftInvert { "lanczos((L+sigma I)^-1 via Jacobi-PCG, deflated, lock+restart, rayleigh-ritz)" } else { "lanczos(cI-L, deflated, lock+restart, rayleigh-ritz)" };
    let mut sp = finish(op, known, pairs, o.k, c, name, rounds, converged, o.h0);
    if o.method == Method::ShiftInvert { sp.inner = Some(inner_tot); sp.sigma = sigma; if inner_tot.2 > 1e-10 { sp.converged = false; } } // Grok v1 #7
    sp
}

// ---- per-component H0 driver (moved from sheaf-cli v1; the one H0 path for graph specs) ----
use sheaf_cell::graph::{coboundary, components, offsets, GraphEdge};

/// Per-component H0 outcome under the single H0 policy (ADR-012/015).
#[derive(Debug, Clone)]
pub struct H0Components {
    /// merged spectrum: the zero pairs of every component, then the k+1 smallest nonzero ones overall.
    /// `spectral.zero_tol` is the LARGEST per-component cut (see `cut_min`); never re-apply it globally.
    pub spectral: Spectral,
    /// every component's cut is consistent with its certified residuals (Weyl); a policy check, not a proof
    pub cut_consistent: bool,
    /// some pair counted zero is numerically resolvable from 0 (diagnostic flag, not a second count)
    pub near_sections: bool,
    pub cut_min: f64,
    pub cut_max: f64,
    /// min over components of (smallest counted-nonzero eigenvalue / that component's cut)
    pub margin_ratio: Option<f64>,
}

/// Block-diagonal spectral analysis: L0 is block diagonal over connected components, so each component is
/// solved on its own (dense EVD of ALL eigenvalues if its dofs <= dense_max, else Lanczos with k, widened x4 up to
/// twice if all k+1 values are zero) and the results are merged. The zero cut is `so.h0` applied to EACH component's
/// own ||L0|_c||_bound -- the one H0 rule (ADR-012/015). No global or strict recount is made.
pub fn spectral_components(dims: &[usize], edges: &[GraphEdge], so: &SpectralOpts) -> H0Components {
    let (comp, known_ok, nc) = components(dims.len(), edges);
    let on = offsets(dims);
    let (d0, _, _) = coboundary(dims, edges);
    let bound = HodgeOp::l0(&d0).norm_bound();
    let mut nodes_of = vec![vec![]; nc]; for (i, &c) in comp.iter().enumerate() { nodes_of[c].push(i); }
    let mut edges_of = vec![vec![]; nc]; for (j, e) in edges.iter().enumerate() { edges_of[comp[e.u]].push(j); }
    let mut local = vec![0usize; dims.len()];
    let (mut known_dim, mut zeros, mut cands): (usize, Vec<EigPair>, Vec<(f64, f64, usize, Vec<f64>)>) = (0, vec![], vec![]);
    let (mut cut_consistent, mut near_sections, mut cut_min, mut cut_max, mut margin) = (true, false, f64::INFINITY, 0.0f64, None::<f64>);
    let (mut rounds, mut converged, mut lower, mut inner, mut method) = (0, true, false, None::<(usize, usize, f64)>, "dense-evd");
    for c in 0..nc {
        let ns = &nodes_of[c];
        let dofs: usize = ns.iter().map(|&i| dims[i]).sum();
        if dofs == 0 { continue; }
        if edges_of[c].is_empty() { known_dim += dofs; continue; } // isolated node: L block = 0, whole stalk is a section (exact)
        for (li, &i) in ns.iter().enumerate() { local[i] = li; }
        let sd: Vec<usize> = ns.iter().map(|&i| dims[i]).collect();
        let se: Vec<GraphEdge> = edges_of[c].iter().map(|&j| { let e = &edges[j]; GraphEdge { u: local[e.u], v: local[e.v], fu: e.fu.clone(), fv: e.fv.clone(), label: String::new() } }).collect();
        let (sd0, son, _) = coboundary(&sd, &se);
        let op = HodgeOp::l0(&sd0);
        let kk = if known_ok[c] { KnownKernel::from_components(&vec![Some(0); ns.len()], &son) } else { KnownKernel::none(op.n) };
        // Grok v1 #1 (repro: two components at scales 1e5 and 1e-1): the cut is PER COMPONENT (`smallest` applies
        // so.h0 to this component's own bound).
        let mut o = so.clone();
        if dofs <= so.dense_max { o.k = dofs.saturating_sub(1 + kk.dim()).max(1); }
        let mut sp = smallest(&op, &kk, &o);
        for _ in 0..2 { if !sp.h0_lower_bound_only { break; } o.k *= 4; sp = smallest(&op, &kk, &o); }
        let local_tol = sp.zero_tol;
        cut_min = cut_min.min(local_tol); cut_max = cut_max.max(local_tol);
        cut_consistent &= sp.cut_consistent_with_residuals(); near_sections |= sp.has_near_sections();
        if let Some(p) = sp.pairs.get(sp.found_zero) { let r = p.value / local_tol; margin = Some(margin.map_or(r, |m: f64| m.min(r))); }
        known_dim += sp.known_kernel_dim; rounds = rounds.max(sp.rounds); converged &= sp.converged; lower |= sp.h0_lower_bound_only;
        if sp.method != "dense-evd" { method = sp.method; }
        if let Some((a, b, r)) = sp.inner { let t = inner.unwrap_or((0, 0, 0.0)); inner = Some((t.0 + a, t.1 + b, t.2.max(r))); }
        for (p, v) in sp.pairs.into_iter().zip(sp.vectors) {
            if p.value <= local_tol { zeros.push(p); } else { cands.push((p.value, p.residual, c, v)); }
        }
    }
    if cut_min.is_infinite() { cut_min = so.h0.zero_tol(0.0); cut_max = cut_min; }
    zeros.sort_by(|a, b| a.value.partial_cmp(&b.value).unwrap());
    cands.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.2.cmp(&b.2)));
    cands.truncate(so.k + 1);
    let n = *on.last().unwrap();
    let mut pairs: Vec<EigPair> = zeros.clone(); let mut vectors: Vec<Vec<f64>> = zeros.iter().map(|_| vec![]).collect();
    for (val, res, c, v) in cands {
        let mut g = vec![0.0; n]; let mut q = 0;
        for &i in &nodes_of[c] { for t in 0..dims[i] { g[on[i] + t] = v[q]; q += 1; } }
        pairs.push(EigPair { value: val, residual: res }); vectors.push(g);
    }
    let found_zero = zeros.len();
    // zero pairs keep value + residual but NOT vectors (empty Vec): callers must skip the first found_zero vectors (Grok v1 #8)
    let spectral = Spectral { method, n, known_kernel_dim: known_dim, pairs, vectors, c_shift: 0.0, lambda_max_bound: bound, zero_tol: cut_max, found_zero, rounds, inner, sigma: so.sigma_rel * bound.max(1.0),
        converged, h0_lower_bound_only: lower };
    H0Components { spectral, cut_consistent: cut_consistent && !lower, near_sections, cut_min, cut_max, margin_ratio: margin }
}

