//! Graph sheaves with PER-EDGE restriction maps (the v1 / m1 data shape) on `CellComplex::graph`.
//! `coboundary` is the one δ⁰ assembly for graphs: it goes through `Sheaf::from_graph` + `Sheaf::coboundary_triplets`
//! (no second implementation).
use crate::num::linalg::{rank_report, svd_full};
use crate::num::sparse::Csr;
use crate::sheaf::{Sheaf, SheafError};
use faer::Mat;

/// Edge j: (δx)_j = F_v x_v − F_u x_u, F_u: R^{dim u} → R^{d_e}, d_e = rows of F_u. u > v and parallel edges are legal.
#[derive(Debug, Clone)]
pub struct GraphEdge { pub u: usize, pub v: usize, pub fu: Mat<f64>, pub fv: Mat<f64>, pub label: String }

pub fn offsets(d: &[usize]) -> Vec<usize> { let mut o = vec![0]; for &x in d { o.push(o.last().unwrap() + x); } o }

/// The cellular sheaf of a graph spec (edge order, multi-edges and orientation kept).
pub fn to_sheaf(dims: &[usize], edges: &[GraphEdge]) -> Result<Sheaf, SheafError> {
    let t: Vec<(usize, usize, Mat<f64>, Mat<f64>)> = edges.iter().map(|e| (e.u, e.v, e.fu.clone(), e.fv.clone())).collect();
    Sheaf::from_graph(dims, &t)
}

/// Assembled δ⁰ (CSR), node offsets, edge offsets. Panics only on map shapes that `to_sheaf` rejects (callers validate).
pub fn coboundary(dims: &[usize], edges: &[GraphEdge]) -> (Csr, Vec<usize>, Vec<usize>) {
    let s = to_sheaf(dims, edges).expect("graph sheaf: map shapes must be validated by the caller");
    let (r, c, t) = s.coboundary_triplets(0);
    (Csr::from_triplets(r, c, t), offsets(dims), offsets(&edges.iter().map(|e| e.fu.nrows()).collect::<Vec<_>>()))
}

/// Connected components: (component of each node, known-constant flag per component, count). A component is
/// "known constant" iff every edge has F_u == F_v and injective: then its constants are exact sections.
pub fn components(n: usize, edges: &[GraphEdge]) -> (Vec<usize>, Vec<bool>, usize) {
    fn find(p: &mut [usize], mut a: usize) -> usize { while p[a] != a { p[a] = p[p[a]]; a = p[a]; } a }
    let mut p: Vec<usize> = (0..n).collect();
    for e in edges { let (a, b) = (find(&mut p, e.u), find(&mut p, e.v)); if a != b { p[a.max(b)] = a.min(b); } }
    let mut lab = vec![usize::MAX; n]; let mut nc = 0; let mut comp = vec![0; n];
    for (i, c) in comp.iter_mut().enumerate() { let r = find(&mut p, i); if lab[r] == usize::MAX { lab[r] = nc; nc += 1; } *c = lab[r]; }
    let mut ok = vec![true; nc];
    for e in edges {
        let c = comp[e.u];
        if !ok[c] { continue; }
        let same = e.fu.nrows() == e.fv.nrows() && e.fu.ncols() == e.fv.ncols() && (0..e.fu.nrows()).all(|r| (0..e.fu.ncols()).all(|q| e.fu[(r, q)] == e.fv[(r, q)]));
        let inj = same && { let s = svd_full(&e.fu); rank_report(&s.s, e.fu.nrows(), e.fu.ncols(), None).rank == e.fu.ncols() };
        if !inj { ok[c] = false; }
    }
    (comp, ok, nc)
}
