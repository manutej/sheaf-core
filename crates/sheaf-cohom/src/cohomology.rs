//! Cochain complex C^0 -> C^1 -> C^2 (k-generic), cohomology dims + harmonic bases, Hodge split.
//! H^0 and every use of rank(delta_0) follow the ONE H0 policy (`sheaf_cell::num::tol::H0Policy`, ADR-012/015): delta_0 is split
//! per connected component of the 1-skeleton and sigma counts as zero iff sigma^2 <= policy.zero_tol(||L0|_c||_bound).
//! rank(delta_k) for k >= 1 (H^1/H^2 on 2-complexes) uses `linalg::rank_report` (numpy default); see tol.rs.
use sheaf_cell::num::family::Numerical;
use sheaf_cell::num::linalg::*;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::H0Policy;
use sheaf_cell::ops::HodgeOp;
use sheaf_cell::sheaf::{Sheaf, SheafError};
use faer::Mat;

#[derive(Debug, Clone)]
pub struct CochainComplex {
    /// delta[k] : C^k -> C^{k+1}, dense in the spike
    pub delta: Vec<Mat<f64>>,
    pub dims: Vec<usize>,
    /// max_k ||d_{k+1} d_k||_F / (||d_{k+1}|| ||d_k||) -- the verified d^2 = 0 residual
    pub d2_residual: f64,
    /// connected component (of the 1-skeleton) of every C^0 dof / every C^1 dof
    pub c0_comp: Vec<usize>,
    pub c1_comp: Vec<usize>,
    pub n_comp: usize,
    /// the H0 zero-cut policy (default ADR-012 per component)
    pub h0: H0Policy,
}

/// delta_0 split per component under the H0 policy: thin SVD on the kept part plus kernel / cokernel bases.
#[derive(Debug, Clone)]
pub struct D0Split {
    /// rank(delta_0) under the policy. `tol` = sqrt(largest per-component lambda cut) (sigma units); sigma_kept_min /
    /// sigma_dropped_max are the global extremes over components.
    pub rank: RankReport,
    /// orthonormal bases: u = im delta_0 (n1 x r), v = im delta_0^T (n0 x r), with delta_0 v_i = s_i u_i
    pub u: Mat<f64>, pub v: Mat<f64>, pub s: Vec<f64>,
    /// orthonormal basis of ker delta_0 = H^0 (n0 x h0) and of ker delta_0^T (n1 x (n1 - r))
    pub kernel: Mat<f64>, pub cokernel: Mat<f64>,
    /// per-component lambda cut actually applied
    pub cuts: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CohomError {
    Sheaf(SheafError),
    /// A5 violated: restriction maps do not compose consistently, so d^2 != 0
    NotAComplex { k: usize, rel_residual: f64 },
    /// Grok r1 #3: rank(d_k) + rank(d_{k-1}) > dim C^k, i.e. the two independent rank cuts are
    /// inconsistent with im d_{k-1} in ker d_k at this tolerance. Never silently wrap.
    RankConflict { k: usize, n: usize, rank_in: usize, rank_out: usize },
}

pub const D2_TOL: f64 = 1e-10;

impl CochainComplex {
    pub fn from_sheaf(s: &Sheaf) -> Result<Self, CohomError> {
        s.validate().map_err(CohomError::Sheaf)?;
        let top = s.complex.dim();
        let dims: Vec<usize> = (0..=top).map(|k| s.cochain_dim(k)).collect();
        let delta: Vec<Mat<f64>> = (0..top).map(|k| s.coboundary_dense(k)).collect();
        let mut worst = 0.0f64;
        for k in 0..delta.len().saturating_sub(1) {
            let p = &delta[k + 1] * &delta[k];
            let denom = (fro(&delta[k + 1]) * fro(&delta[k])).max(1e-300);
            let r = fro(&p) / denom;
            worst = worst.max(r);
            if r > D2_TOL { return Err(CohomError::NotAComplex { k, rel_residual: r }); }
        }
        // components of the 1-skeleton (union-find over vertices via the vertex->edge incidences)
        let nv = s.complex.n_cells[0];
        let mut parent: Vec<usize> = (0..nv).collect();
        fn find(p: &mut [usize], mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
        let mut edge_vert: Vec<Option<usize>> = vec![None; s.complex.n_cells.get(1).copied().unwrap_or(0)];
        if let Some(inc0) = s.complex.inc.first() {
            for inc in inc0 {
                match edge_vert[inc.cofacet] { None => edge_vert[inc.cofacet] = Some(inc.facet),
                    Some(a) => { let (ra, rb) = (find(&mut parent, a), find(&mut parent, inc.facet)); if ra != rb { parent[ra] = rb; } } }
            }
        }
        let mut label = vec![usize::MAX; nv]; let mut n_comp = 0;
        let root: Vec<usize> = (0..nv).map(|v| find(&mut parent, v)).collect();
        for &r in &root { if label[r] == usize::MAX { label[r] = n_comp; n_comp += 1; } }
        let vcomp: Vec<usize> = (0..nv).map(|v| label[root[v]]).collect();
        let c0_comp: Vec<usize> = (0..nv).flat_map(|v| std::iter::repeat(vcomp[v]).take(s.stalk[0][v])).collect();
        // an edge with no incidences cannot occur in a validated complex; fall back to component 0
        let c1_comp: Vec<usize> = edge_vert.iter().enumerate()
            .flat_map(|(e, a)| std::iter::repeat(a.map_or(0, |a| vcomp[a])).take(s.stalk.get(1).map_or(0, |st| st[e]))).collect();
        Ok(CochainComplex { delta, dims, d2_residual: worst, c0_comp, c1_comp, n_comp, h0: H0Policy::adr012() })
    }
    /// Same complex, different H0 zero-cut policy (e.g. a context-chosen absolute permitted delta).
    pub fn with_h0_policy(mut self, h0: H0Policy) -> Self { self.h0 = h0; self }

    /// delta_0 split per component under `self.h0` (the only H0 rule).
    pub fn d0_split(&self) -> D0Split {
        let n0 = self.dims[0];
        let n1 = if self.delta.is_empty() { 0 } else { self.delta[0].nrows() };
        let mut us: Vec<Vec<f64>> = vec![]; let mut vs: Vec<Vec<f64>> = vec![]; let mut ss: Vec<f64> = vec![];
        let mut ker: Vec<Vec<f64>> = vec![]; let mut cok: Vec<Vec<f64>> = vec![]; let mut cuts = vec![];
        let (mut smax, mut kept_min, mut drop_max, mut cut_max) = (0.0f64, None::<f64>, None::<f64>, 0.0f64);
        for c in 0..self.n_comp.max(if n0 == 0 { 0 } else { 1 }) {
            let cols: Vec<usize> = (0..n0).filter(|&i| self.c0_comp[i] == c).collect();
            let rows: Vec<usize> = (0..n1).filter(|&i| self.c1_comp[i] == c).collect();
            let b = Mat::<f64>::from_fn(rows.len(), cols.len(), |i, j| self.delta[0][(rows[i], cols[j])]);
            let mut t = vec![]; for i in 0..b.nrows() { for j in 0..b.ncols() { if b[(i, j)] != 0.0 { t.push((i, j, b[(i, j)])); } } }
            let bound = HodgeOp::l0(&Csr::from_triplets(b.nrows(), b.ncols(), t)).norm_bound();
            let cut = self.h0.zero_tol(bound); cuts.push(cut); cut_max = cut_max.max(cut);
            let sv = svd_full(&b);
            let embed = |src: &Mat<f64>, j: usize, idx: &[usize], n: usize| { let mut g = vec![0.0; n]; for (q, &i) in idx.iter().enumerate() { g[i] = src[(q, j)]; } g };
            let mut r = 0;
            for (i, &x) in sv.s.iter().enumerate() {
                smax = smax.max(x);
                if x * x > cut { r += 1; kept_min = Some(kept_min.map_or(x, |m: f64| m.min(x))); us.push(embed(&sv.u, i, &rows, n1)); vs.push(embed(&sv.v, i, &cols, n0)); ss.push(x); }
                else { drop_max = Some(drop_max.map_or(x, |m: f64| m.max(x))); }
            }
            // svd_full sorts singular values descending, so the kept ones are the leading r columns
            for j in r..cols.len() { ker.push(embed(&sv.v, j, &cols, n0)); }
            for j in r..rows.len() { cok.push(embed(&sv.u, j, &rows, n1)); }
        }
        let to_mat = |cs: &[Vec<f64>], n: usize| Mat::<f64>::from_fn(n, cs.len(), |i, j| cs[j][i]);
        let rank = RankReport { rank: ss.len(), tol: cut_max.sqrt(), sigma_max: smax, sigma_kept_min: kept_min, sigma_dropped_max: drop_max };
        D0Split { rank, u: to_mat(&us, n1), v: to_mat(&vs, n0), s: ss, kernel: to_mat(&ker, n0), cokernel: to_mat(&cok, n1), cuts }
    }
    fn d_in(&self, k: usize) -> Mat<f64> { // delta_{k-1}: C^{k-1} -> C^k
        if k == 0 { Mat::zeros(self.dims[0], 0) } else { self.delta[k - 1].clone() }
    }
    fn d_out(&self, k: usize) -> Mat<f64> { // delta_k: C^k -> C^{k+1}
        if k < self.delta.len() { self.delta[k].clone() } else { Mat::zeros(0, self.dims[k]) }
    }

    /// H^k = ker d_k / im d_{k-1}. dim via rank-nullity on two SVDs; basis = harmonic reps
    /// = ker of [d_k ; d_{k-1}^T] (third independent SVD). Both are reported so they can disagree loudly.
    pub fn cohomology(&self, k: usize) -> Result<Cohomology, CohomError> {
        let none = |_n: usize| RankReport { rank: 0, tol: 0.0, sigma_max: 0.0, sigma_kept_min: None, sigma_dropped_max: None };
        if k == 0 {
            // H^0 = ker delta_0 under the single H0 policy
            let sp = self.d0_split();
            let n = self.dims[0];
            return Ok(Cohomology { k, dim: sp.kernel.ncols(), cochain_dim: n, rank_in: none(n), rank_out: sp.rank.clone(), harmonic_basis: sp.kernel, laplacian_rank: sp.rank });
        }
        if k == 1 && self.delta.len() == 1 {
            // graph: H^1 = ker delta_0^T, rank(delta_0) under the H0 policy (so H^0 and H^1 always satisfy rank-nullity)
            let sp = self.d0_split();
            let n = self.dims[1];
            return Ok(Cohomology { k, dim: sp.cokernel.ncols(), cochain_dim: n, rank_in: sp.rank.clone(), rank_out: none(0), harmonic_basis: sp.cokernel, laplacian_rank: sp.rank });
        }
        let (a, b) = (self.d_in(k), self.d_out(k));
        let sa = svd_full(&a); let sb = svd_full(&b);
        // rank(delta_{k-1}): for k == 1 it is rank(delta_0), which follows the H0 policy
        let ra = if k == 1 { self.d0_split().rank } else { rank_report(&sa.s, a.nrows(), a.ncols(), None) };
        let rb = rank_report(&sb.s, b.nrows(), b.ncols(), None);
        let n = self.dims[k];
        let dim = n.checked_sub(rb.rank + ra.rank)
            .ok_or(CohomError::RankConflict { k, n, rank_in: ra.rank, rank_out: rb.rank })?;
        let m = vstack(&b, &transpose(&a));
        let sm = svd_full(&m);
        let rm = rank_report(&sm.s, m.nrows(), m.ncols(), None);
        let basis = cols(&sm.v, rm.rank, n);
        Ok(Cohomology { k, dim, cochain_dim: n, rank_in: ra, rank_out: rb, harmonic_basis: basis, laplacian_rank: rm })
    }

    /// Orthogonal Hodge decomposition of a k-cochain: w = exact (im d_{k-1}) + harmonic + coexact (im d_k^T).
    pub fn hodge(&self, k: usize, w: &[f64]) -> Hodge {
        assert_eq!(w.len(), self.dims[k]);
        let (a, b) = (self.d_in(k), self.d_out(k));
        // thin SVD of d_{k-1} (in) and d_k (out); any factor that is delta_0 uses the H0 policy split
        let thin = |m: &Mat<f64>, is_d0: bool| -> (Mat<f64>, Mat<f64>, Vec<f64>, RankReport) {
            if is_d0 { let sp = self.d0_split(); return (sp.u, sp.v, sp.s, sp.rank); }
            let sv = svd_full(m); let rr = rank_report(&sv.s, m.nrows(), m.ncols(), None);
            (cols(&sv.u, 0, rr.rank), cols(&sv.v, 0, rr.rank), sv.s[..rr.rank].to_vec(), rr)
        };
        let (ua, va, sa, ra_rep) = thin(&a, k == 1);
        let (_ub, vb, _sb, _rb) = thin(&b, k == 0 && !self.delta.is_empty());
        let qa = ua;                            // orthonormal basis of im d_{k-1}
        let qb = vb;                            // orthonormal basis of im d_k^T = (ker d_k)^perp
        let exact = project(&qa, w);
        let coexact = project(&qb, w);
        let harmonic: Vec<f64> = (0..w.len()).map(|i| w[i] - exact[i] - coexact[i]).collect();
        // potential x with d_{k-1} x = exact (least-squares "best global assignment")
        let potential: Vec<f64> = {
            let c = matvec_t(&qa, w);
            (0..a.ncols()).map(|j| (0..sa.len()).map(|i| va[(j, i)] * c[i] / sa[i]).sum()).collect()
        };
        let wn = norm(w).max(1e-300);
        let cert = HodgeCertificate {
            // Grok r1 #8: e+h+c=w holds by construction, so check an independent identity instead:
            // the potential must reproduce the exact part through the *assembled* coboundary.
            potential_residual: norm(&{ let dx = matvec(&a, &potential); (0..w.len()).map(|i| dx[i] - exact[i]).collect::<Vec<_>>() }) / wn,
            ortho_exact_harm: dot(&exact, &harmonic).abs() / (wn * wn),
            ortho_exact_coex: dot(&exact, &coexact).abs() / (wn * wn),
            ortho_harm_coex: dot(&harmonic, &coexact).abs() / (wn * wn),
            harmonic_closed: norm(&matvec(&b, &harmonic)) / wn,       // d_k h ~ 0
            harmonic_coclosed: norm(&matvec_t(&a, &harmonic)) / wn,   // d_{k-1}^T h ~ 0
        };
        // Grok r1 #10: the potential is only as trustworthy as sigma_max / sigma_kept_min
        let potential_cond = match ra_rep.sigma_kept_min { Some(m) => ra_rep.sigma_max / m, None => 1.0 };
        Hodge { exact, harmonic, coexact, potential, potential_cond, cert }
    }
}

#[derive(Debug, Clone)]
pub struct Cohomology {
    pub k: usize,
    pub dim: usize,
    pub cochain_dim: usize,
    pub rank_in: RankReport,
    pub rank_out: RankReport,
    /// columns = orthonormal harmonic representatives of H^k (ker L_k)
    pub harmonic_basis: Mat<f64>,
    pub laplacian_rank: RankReport,
}
impl Cohomology {
    /// two independent dimension computations agree
    pub fn consistent(&self) -> bool { self.harmonic_basis.ncols() == self.dim }
    pub fn min_gap_ratio(&self) -> f64 { self.rank_in.gap_ratio().min(self.rank_out.gap_ratio()).min(self.laplacian_rank.gap_ratio()) }
}

#[derive(Debug, Clone)]
pub struct HodgeCertificate {
    pub potential_residual: f64,
    pub ortho_exact_harm: f64,
    pub ortho_exact_coex: f64,
    pub ortho_harm_coex: f64,
    pub harmonic_closed: f64,
    pub harmonic_coclosed: f64,
}
impl HodgeCertificate {
    pub fn max(&self) -> f64 {
        [self.potential_residual, self.ortho_exact_harm, self.ortho_exact_coex, self.ortho_harm_coex, self.harmonic_closed, self.harmonic_coclosed]
            .into_iter().fold(0.0, f64::max)
    }
}

#[derive(Debug, Clone)]
pub struct Hodge {
    pub exact: Vec<f64>,
    pub harmonic: Vec<f64>,
    pub coexact: Vec<f64>,
    pub potential: Vec<f64>,
    /// condition number of d_{k-1} on its range (sigma_max / smallest kept sigma)
    pub potential_cond: f64,
    pub cert: HodgeCertificate,
}
impl Hodge {
    pub fn energies(&self) -> (f64, f64, f64) {
        let e = |v: &Vec<f64>| dot(v, v);
        (e(&self.exact), e(&self.harmonic), e(&self.coexact))
    }
}

/// The nearest global section to a 0-cochain x: argmin_{s ∈ H⁰} ‖x − s‖ = P_{ker δ⁰} x (= the harmonic part of the
/// degree-0 Hodge split). Certified residual: ‖δ⁰ s‖ / ‖x‖ (s really is a section) and the distance ‖x − s‖.
#[derive(Debug, Clone)]
pub struct NearestSection { pub section: Vec<f64>, pub distance: f64, pub section_residual: f64, pub h0_dim: usize }

impl CochainComplex {
    pub fn nearest_section(&self, x: &[f64]) -> Numerical<NearestSection> {
        assert_eq!(x.len(), self.dims[0]);
        let k = self.d0_split().kernel;
        let section = project(&k, x);
        let distance = norm(&(0..x.len()).map(|i| x[i] - section[i]).collect::<Vec<_>>());
        let section_residual = if self.delta.is_empty() { 0.0 } else { norm(&matvec(&self.delta[0], &section)) / norm(x).max(1e-300) };
        Numerical::f64(NearestSection { section, distance, section_residual, h0_dim: k.ncols() })
    }
}
