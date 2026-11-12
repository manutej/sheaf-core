//! The m1 contract math (SPEC v1.18) in Rust: coboundary δ⁰ from per-edge maps with heterogeneous stalks, L₀ = δᵀδ,
//! per-edge discord and energy. Port of the Python reference `sheafcore.py` (frozen in the m1 gate); ADR-016.
//!
//! Semantics pinned to the reference (parity tests in `crates/sheaf-cli/tests/m1_parity.rs`):
//! * (δx)_e = F_v x_v − F_u x_u, edge e owns rows `oe[e]..oe[e]+de[e]`, node i owns columns `ov[i]..ov[i]+dv[i]`.
//! * Only NONZERO map entries are stored (`np.nonzero`); a self-loop (u == v, legal in the core) sums its two blocks and
//!   keeps a cancelled entry as an explicitly stored 0.0 (scipy's COO→CSR `sum_duplicates`), so the pattern matches.
//! * Columns are sorted within each row (canonical CSR).
//! * `laplacian` keeps every structural product entry; `drop_zeros` then removes products that cancel to exactly 0.0,
//!   which gives scipy's SpGEMM pattern (`csr_matmat` stores only nonzero sums; checked by the stored-entry counts in
//!   the parity tests). Columns sorted.
//! * `edge_discord` sums squared residuals over each edge's rows in row order (= `np.bincount`).
//! Input admission (R-DIM-02 / R-FINITE-01 / R-PERF-18 caps with their exact messages) stays with the caller: in the
//! gate that is the Python shim (`python/m1/sheafimpl_rs.py`), in Rust it is `sheaf_contract::SheafSpec::validate`
//! (sheaf-cli's m1 adapter) plus `check_*` below.
use crate::num::sparse::Csr;

/// One m1 edge: maps are row-major, `fu` is de × dv[u], `fv` is de × dv[v].
#[derive(Debug, Clone, Copy)]
pub struct M1Edge<'a> { pub u: usize, pub v: usize, pub de: usize, pub fu: &'a [f64], pub fv: &'a [f64] }

#[derive(Debug, Clone, PartialEq)]
pub enum M1Error { Dim(String), Finite(String) }
impl std::fmt::Display for M1Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { M1Error::Dim(s) => write!(f, "R-DIM-02: {s}"), M1Error::Finite(s) => write!(f, "R-FINITE-01: {s}") }
    }
}

pub fn offsets(d: &[usize]) -> Vec<usize> { let mut o = Vec::with_capacity(d.len() + 1); o.push(0); let mut s = 0; for &x in d { s += x; o.push(s); } o }

/// Shapes, endpoint range and finiteness (the reference's R-DIM-02 / R-ORIENT-01 / R-FINITE-01 checks, in Rust).
pub fn check_edges(dv: &[usize], edges: &[M1Edge]) -> Result<(), M1Error> {
    for (j, e) in edges.iter().enumerate() {
        if e.u >= dv.len() || e.v >= dv.len() { return Err(M1Error::Dim(format!("edge {j} endpoint out of range"))); }
        if e.fu.len() != e.de * dv[e.u] || e.fv.len() != e.de * dv[e.v] { return Err(M1Error::Dim(format!("edge {j}: map shape != ({}, dv)", e.de))); }
        if !e.fu.iter().chain(e.fv).all(|x| x.is_finite()) { return Err(M1Error::Finite(format!("edge {j}: non-finite map"))); }
    }
    Ok(())
}

/// Upper bound on stored entries (nonzero map entries); the exact count is ≤ this (self-loops merge).
pub fn coboundary_nnz_bound(edges: &[M1Edge]) -> usize { edges.iter().map(|e| e.fu.iter().chain(e.fv).filter(|&&x| x != 0.0).count()).sum() }

/// δ⁰ as canonical CSR. Callers must have run `check_edges` (indices are trusted here).
pub fn coboundary(dv: &[usize], edges: &[M1Edge]) -> Csr {
    let ov = offsets(dv);
    let nrows: usize = edges.iter().map(|e| e.de).sum();
    let mut ptr = Vec::with_capacity(nrows + 1); ptr.push(0);
    let cap = coboundary_nnz_bound(edges);
    let (mut idx, mut val) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
    let mut row: Vec<(usize, f64)> = Vec::new();
    for e in edges {
        let (du, dvv) = (dv[e.u], dv[e.v]);
        for r in 0..e.de {
            row.clear();
            for c in 0..du { let a = e.fu[r * du + c]; if a != 0.0 { row.push((ov[e.u] + c, -a)); } }
            for c in 0..dvv { let a = e.fv[r * dvv + c]; if a != 0.0 { row.push((ov[e.v] + c, a)); } }
            row.sort_by_key(|p| p.0);
            let mut k = 0;
            while k < row.len() {
                let (c, mut s) = row[k]; k += 1;
                while k < row.len() && row[k].0 == c { s += row[k].1; k += 1; }   // u == v only: sum, keep an explicit 0
                idx.push(c); val.push(s);
            }
            ptr.push(idx.len());
        }
    }
    Csr { nrows, ncols: *ov.last().unwrap(), ptr, idx, val }
}

/// CSR of Aᵀ (counting sort; columns of the result sorted because rows of A are visited in order).
pub fn transpose(a: &Csr) -> Csr {
    let mut cnt = vec![0usize; a.ncols + 1];
    for &j in &a.idx { cnt[j + 1] += 1; }
    for j in 0..a.ncols { cnt[j + 1] += cnt[j]; }
    let ptr = cnt.clone(); let mut next = cnt; let nnz = a.idx.len();
    let (mut idx, mut val) = (vec![0usize; nnz], vec![0.0; nnz]);
    for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { let j = a.idx[p]; let q = next[j]; idx[q] = i; val[q] = a.val[p]; next[j] += 1; } }
    Csr { nrows: a.ncols, ncols: a.nrows, ptr, idx, val }
}

/// Exact stored-entry count of δᵀδ (symbolic Gustavson), so a caller can admit / allocate the output first.
pub fn laplacian_nnz(d: &Csr) -> usize {
    let t = transpose(d); let mut mark = vec![usize::MAX; d.ncols]; let mut nnz = 0;
    for i in 0..t.nrows { for p in t.ptr[i]..t.ptr[i + 1] { let k = t.idx[p]; for q in d.ptr[k]..d.ptr[k + 1] { let j = d.idx[q]; if mark[j] != i { mark[j] = i; nnz += 1; } } } }
    nnz
}

/// L₀ = δᵀδ (CSR, sorted columns, structural zeros kept; apply `drop_zeros` for the reference's pattern).
pub fn laplacian(d: &Csr) -> Csr {
    let t = transpose(d); let n = d.ncols;
    let mut mark = vec![usize::MAX; n]; let mut acc = vec![0.0; n]; let mut cols: Vec<usize> = Vec::new();
    let mut ptr = Vec::with_capacity(n + 1); ptr.push(0); let (mut idx, mut val) = (Vec::new(), Vec::new());
    for i in 0..n {
        cols.clear();
        for p in t.ptr[i]..t.ptr[i + 1] {
            let (k, a) = (t.idx[p], t.val[p]);
            for q in d.ptr[k]..d.ptr[k + 1] { let j = d.idx[q]; if mark[j] != i { mark[j] = i; acc[j] = 0.0; cols.push(j); } acc[j] += a * d.val[q]; }
        }
        cols.sort_unstable();
        for &j in &cols { idx.push(j); val.push(acc[j]); }
        ptr.push(idx.len());
    }
    Csr { nrows: n, ncols: n, ptr, idx, val }
}

/// discord(e) = ‖(δx)_e‖², edge e owning `de[e]` consecutive rows. Caller checks Σde = rows and len(x) = cols.
pub fn edge_discord(d: &Csr, x: &[f64], de: &[usize]) -> Vec<f64> {
    let r = d.mul(x); let mut out = Vec::with_capacity(de.len()); let mut row = 0;
    for &k in de { let mut s = 0.0; for i in row..row + k { s += r[i] * r[i]; } out.push(s); row += k; }
    out
}
/// Dense row-major variant (an m1 caller may hand a dense δ).
pub fn edge_discord_dense(a: &[f64], nrows: usize, ncols: usize, x: &[f64], de: &[usize]) -> Vec<f64> {
    let r: Vec<f64> = (0..nrows).map(|i| a[i * ncols..(i + 1) * ncols].iter().zip(x).map(|(p, q)| p * q).sum()).collect();
    let mut out = Vec::with_capacity(de.len()); let mut row = 0;
    for &k in de { let mut s = 0.0; for i in row..row + k { s += r[i] * r[i]; } out.push(s); row += k; }
    out
}

/// E(x) = xᵀ L x.
pub fn energy(l: &Csr, x: &[f64]) -> f64 { let lx = l.mul(x); x.iter().zip(&lx).map(|(a, b)| a * b).sum() }
pub fn energy_dense(a: &[f64], n: usize, x: &[f64]) -> f64 {
    (0..n).map(|i| x[i] * a[i * n..(i + 1) * n].iter().zip(x).map(|(p, q)| p * q).sum::<f64>()).sum()
}

/// Column groups coupled through the NONZERO pattern of A (connected components of the bipartite row/column graph;
/// explicitly stored zeros are not couplings). Each group: (rows with a nonzero in it, columns), both ascending.
/// Groups are ordered by their smallest column.
pub fn pattern_components(a: &Csr) -> Vec<(Vec<usize>, Vec<usize>)> {
    let n = a.ncols; let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    let mut row_rep = vec![usize::MAX; a.nrows];
    for (i, rep) in row_rep.iter_mut().enumerate() {
        let mut first = usize::MAX;
        for q in a.ptr[i]..a.ptr[i + 1] {
            if a.val[q] == 0.0 { continue; }
            let c = a.idx[q];
            if first == usize::MAX { first = c; } else { let (x, y) = (find(&mut p, first), find(&mut p, c)); if x != y { p[x.max(y)] = x.min(y); } }
        }
        *rep = first;
    }
    let mut lab = vec![usize::MAX; n]; let mut groups: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    for c in 0..n { let r = find(&mut p, c); if lab[r] == usize::MAX { lab[r] = groups.len(); groups.push((vec![], vec![])); } groups[lab[r]].1.push(c); }
    for i in 0..a.nrows { if row_rep[i] != usize::MAX { let g = lab[find(&mut p, row_rep[i])]; groups[g].0.push(i); } }
    groups
}

/// Connected components of the undirected graph with an edge i–j for every NONZERO L[i,j] (any square L; the pattern
/// need not be symmetric, exactly like scipy `connected_components(directed=False)`). Ordered by smallest index.
pub fn adjacency_components(l: &Csr) -> Vec<Vec<usize>> {
    let n = l.nrows; let mut p: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    for i in 0..n { for q in l.ptr[i]..l.ptr[i + 1] { if l.val[q] != 0.0 { let (x, y) = (find(&mut p, i), find(&mut p, l.idx[q])); if x != y { p[x.max(y)] = x.min(y); } } } }
    let mut lab = vec![usize::MAX; n]; let mut g: Vec<Vec<usize>> = Vec::new();
    for c in 0..n { let r = find(&mut p, c); if lab[r] == usize::MAX { lab[r] = g.len(); g.push(vec![]); } g[lab[r]].push(c); }
    g
}

/// Drop exactly-zero stored entries (scipy's SpGEMM never stores a product that sums to 0.0).
pub fn drop_zeros(a: Csr) -> Csr {
    let mut ptr = Vec::with_capacity(a.nrows + 1); ptr.push(0); let (mut idx, mut val) = (Vec::new(), Vec::new());
    for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { if a.val[p] != 0.0 { idx.push(a.idx[p]); val.push(a.val[p]); } } ptr.push(idx.len()); }
    Csr { nrows: a.nrows, ncols: a.ncols, ptr, idx, val }
}
