//! Cellular sheaf: stalk dimension per cell, restriction map per codim-1 incidence.
//!
//! A4 stalks are finite-dimensional real vector spaces R^{d(c)};
//! A5 only codimension-1 restriction maps are stored; functoriality on codim-2
//!    is NOT assumed, it is checked via d_{k+1} d_k = 0 on the assembled coboundary.
use crate::complex::CellComplex;
use faer::Mat;

#[derive(Debug, Clone)]
pub struct Sheaf {
    pub complex: CellComplex,
    /// stalk[k][i] = dim of stalk on k-cell i
    pub stalk: Vec<Vec<usize>>,
    /// maps[k][j] = restriction F(facet) -> F(cofacet) for complex.inc[k][j]; shape d(cofacet) x d(facet)
    pub maps: Vec<Vec<Mat<f64>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SheafError {
    MapShape { dim: usize, incidence: usize, expected: (usize, usize), got: (usize, usize) },
    Complex(crate::complex::ComplexError),
}

impl Sheaf {
    /// Graph sheaf with PER-EDGE maps on `CellComplex::graph` (edge order, multi-edges and u > v kept):
    /// edges[j] = (u, v, F_u, F_v) with F_u: R^{dims[u]} -> R^{d_e}, F_v: R^{dims[v]} -> R^{d_e}, d_e = rows of F_u.
    /// (delta x)_j = F_v x_v - F_u x_u, the m1 / v1 convention.
    pub fn from_graph(dims: &[usize], edges: &[(usize, usize, Mat<f64>, Mat<f64>)]) -> Result<Self, SheafError> {
        let pairs: Vec<(usize, usize)> = edges.iter().map(|e| (e.0, e.1)).collect();
        let complex = CellComplex::graph(dims.len(), &pairs).map_err(SheafError::Complex)?;
        let maps: Vec<Mat<f64>> = edges.iter().flat_map(|e| [e.2.clone(), e.3.clone()]).collect();
        let stalk = vec![dims.to_vec(), edges.iter().map(|e| e.2.nrows()).collect()];
        let s = Sheaf { complex, stalk, maps: vec![maps] };
        s.validate()?;
        Ok(s)
    }
    /// Constant sheaf R^d: identity restrictions everywhere.
    pub fn constant(complex: CellComplex, d: usize) -> Self {
        let stalk = complex.n_cells.iter().map(|&n| vec![d; n]).collect();
        let maps = complex.inc.iter().map(|l| l.iter().map(|_| Mat::<f64>::identity(d, d)).collect()).collect();
        Sheaf { complex, stalk, maps }
    }

    pub fn validate(&self) -> Result<(), SheafError> {
        self.complex.validate().map_err(SheafError::Complex)?;
        for (k, list) in self.complex.inc.iter().enumerate() {
            for (j, inc) in list.iter().enumerate() {
                let e = (self.stalk[k + 1][inc.cofacet], self.stalk[k][inc.facet]);
                let m = &self.maps[k][j];
                if (m.nrows(), m.ncols()) != e {
                    return Err(SheafError::MapShape { dim: k, incidence: j, expected: e, got: (m.nrows(), m.ncols()) });
                }
            }
        }
        Ok(())
    }

    /// offsets[k][i] = start of cell i's block in C^k; last entry = dim C^k
    pub fn offsets(&self, k: usize) -> Vec<usize> {
        let mut o = Vec::with_capacity(self.stalk[k].len() + 1);
        let mut s = 0; o.push(0);
        for &d in &self.stalk[k] { s += d; o.push(s); }
        o
    }
    pub fn cochain_dim(&self, k: usize) -> usize {
        if k < self.stalk.len() { self.stalk[k].iter().sum() } else { 0 }
    }

    /// Sparse triplets of delta_k : C^k -> C^{k+1}; block (tau, sigma) = [sigma:tau] F_{sigma<=tau}.
    pub fn coboundary_triplets(&self, k: usize) -> (usize, usize, Vec<(usize, usize, f64)>) {
        let rows = self.cochain_dim(k + 1);
        let cols = self.cochain_dim(k);
        if k >= self.complex.inc.len() { return (rows, cols, vec![]); }
        let (or, oc) = (self.offsets(k + 1), self.offsets(k));
        let mut t = Vec::new();
        for (j, inc) in self.complex.inc[k].iter().enumerate() {
            let m = &self.maps[k][j];
            for r in 0..m.nrows() { for c in 0..m.ncols() {
                let v = m[(r, c)] * inc.sign as f64;
                if v != 0.0 { t.push((or[inc.cofacet] + r, oc[inc.facet] + c, v)); }
            }}
        }
        (rows, cols, t)
    }

    pub fn coboundary_dense(&self, k: usize) -> Mat<f64> {
        let (r, c, t) = self.coboundary_triplets(k);
        let mut m = Mat::<f64>::zeros(r, c);
        for (i, j, v) in t { m[(i, j)] += v; }
        m
    }
}
