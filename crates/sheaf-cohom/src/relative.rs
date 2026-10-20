//! Relative cohomology in degree 0: H⁰(X,B;F) = sections vanishing on the boundary subcomplex B
//! = ker(δ⁰ restricted to the interior dofs) = ker L_II (L₀ = δ⁰ᵀδ⁰). It is the UNIQUENESS certificate for harmonic
//! extension (Hansen–Ghrist, Prop. 4.1): the extension with boundary data on B is unique iff H⁰(X,B;F) = 0; otherwise
//! it is determined only up to ker L_II (and `sheaf-solve`'s CG returns the minimum-norm / start-projected one).
//!
//! Same single H0 rule as absolute H⁰ (ADR-012/015): the interior block of δ⁰ is split into its connected blocks
//! (interior columns linked through shared nonzero rows) and σ counts as zero iff σ² ≤ policy.zero_tol(‖L_II|_c‖_bound).
//! With B = ∅ the blocks are the connected components and this equals `CochainComplex::cohomology(0).dim`.
use sheaf_cell::num::family::Numerical;
use sheaf_cell::num::linalg::svd_full;
use sheaf_cell::num::sparse::Csr;
use sheaf_cell::num::tol::H0Policy;
use sheaf_cell::ops::HodgeOp;
use faer::Mat;

#[derive(Debug, Clone)]
pub struct RelativeH0 {
    /// dim H⁰(X,B;F); 0 ⇔ the harmonic extension is unique
    pub dim: usize,
    /// orthonormal basis of H⁰(X,B;F) as full-length C⁰ vectors (zero on B)
    pub basis: Mat<f64>,
    pub n_boundary_dofs: usize,
    pub n_blocks: usize,
    /// per-block λ cut applied, and the global extremes of kept / dropped σ (the gap that backs `dim`)
    pub cuts: Vec<f64>,
    pub sigma_kept_min: Option<f64>,
    pub sigma_dropped_max: Option<f64>,
    pub policy: H0Policy,
}
impl RelativeH0 { pub fn unique(&self) -> bool { self.dim == 0 } }

#[derive(Debug, Clone, PartialEq)]
pub enum RelativeError { BoundaryIndex(usize), BlockTooLarge { cols: usize, dense_max: usize } }

/// H⁰(X,B;F) for δ⁰ (`d0`, m × n CSR) and a boundary dof mask. Each interior block is solved by a dense SVD; a block
/// with more than `dense_max` columns is refused (no silent approximation).
pub fn h0_relative(d0: &Csr, boundary: &[usize], policy: H0Policy, dense_max: usize) -> Result<Numerical<RelativeH0>, RelativeError> {
    let n = d0.ncols;
    let mut is_b = vec![false; n];
    for &b in boundary { if b >= n { return Err(RelativeError::BoundaryIndex(b)); } is_b[b] = true; }
    // union interior columns that share a row
    let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, mut a: usize) -> usize { while p[a] != a { p[a] = p[p[a]]; a = p[a]; } a }
    for r in 0..d0.nrows {
        let mut first: Option<usize> = None;
        for k in d0.ptr[r]..d0.ptr[r + 1] {
            let c = d0.idx[k]; if is_b[c] || d0.val[k] == 0.0 { continue; }
            match first { None => first = Some(c), Some(f) => { let (a, b) = (find(&mut p, f), find(&mut p, c)); if a != b { p[a.max(b)] = a.min(b); } } }
        }
    }
    let interior: Vec<usize> = (0..n).filter(|&i| !is_b[i]).collect();
    let mut block_of = vec![usize::MAX; n]; let mut blocks: Vec<Vec<usize>> = vec![];
    for &i in &interior { let r = find(&mut p, i); if block_of[r] == usize::MAX { block_of[r] = blocks.len(); blocks.push(vec![]); } blocks[block_of[r]].push(i); }
    let mut col_pos = vec![usize::MAX; n];
    let mut rows_of: Vec<Vec<usize>> = vec![vec![]; blocks.len()];
    for (bi, cols) in blocks.iter().enumerate() { for (q, &c) in cols.iter().enumerate() { col_pos[c] = q; } let _ = bi; }
    for r in 0..d0.nrows {
        // ADR-017 round 2 (P-NUM-02): a row belongs to the block of its first NONZERO interior entry. Explicitly stored
        // zeros (m1 coboundary on a cancelling self-loop, any hand-built CSR) are not couplings: picking the first stored
        // entry put the row into the wrong block and returned a "kernel" vector with δg ≠ 0.
        if let Some(c) = (d0.ptr[r]..d0.ptr[r + 1]).filter(|&k| d0.val[k] != 0.0).map(|k| d0.idx[k]).find(|&c| !is_b[c]) { rows_of[block_of[find(&mut p, c)]].push(r); }
    }
    let (mut basis_cols, mut cuts) = (vec![], vec![]);
    let (mut kept_min, mut drop_max) = (None::<f64>, None::<f64>);
    for (bi, cols) in blocks.iter().enumerate() {
        if cols.len() > dense_max { return Err(RelativeError::BlockTooLarge { cols: cols.len(), dense_max }); }
        let rows = &rows_of[bi];
        let mut b = Mat::<f64>::zeros(rows.len(), cols.len());
        let mut t = vec![];
        for (i, &r) in rows.iter().enumerate() { for k in d0.ptr[r]..d0.ptr[r + 1] { let c = d0.idx[k]; if !is_b[c] && d0.val[k] != 0.0 { b[(i, col_pos[c])] += d0.val[k]; } } }   // all nonzero interior entries of a row share its block (union above)
        for i in 0..b.nrows() { for j in 0..b.ncols() { if b[(i, j)] != 0.0 { t.push((i, j, b[(i, j)])); } } }
        let cut = policy.zero_tol(HodgeOp::l0(&Csr::from_triplets(b.nrows(), b.ncols(), t)).norm_bound()); cuts.push(cut);
        let sv = svd_full(&b);
        let mut r = 0;
        for &x in &sv.s { if x * x > cut { r += 1; kept_min = Some(kept_min.map_or(x, |m: f64| m.min(x))); } else { drop_max = Some(drop_max.map_or(x, |m: f64| m.max(x))); } }
        for j in r..cols.len() { let mut g = vec![0.0; n]; for (q, &c) in cols.iter().enumerate() { g[c] = sv.v[(q, j)]; } basis_cols.push(g); }
    }
    let basis = Mat::<f64>::from_fn(n, basis_cols.len(), |i, j| basis_cols[j][i]);
    Ok(Numerical::f64(RelativeH0 { dim: basis_cols.len(), basis, n_boundary_dofs: boundary.len(), n_blocks: blocks.len(), cuts, sigma_kept_min: kept_min, sigma_dropped_max: drop_max, policy }))
}
