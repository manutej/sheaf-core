//! Operators on cochains: matrix-free Hodge Laplacians L_k = d_{k-1} d_{k-1}^T + d_k^T d_k, energy and discord.
use crate::num::sparse::{norm1_inf, par_mul, transpose, Csr};
use faer::Mat;

/// L = A A^T (down, A = d_{k-1}) + B^T B (up, B = d_k), stored as (M, M^T) CSR pairs.
pub struct HodgeOp { pub n: usize, down: Option<(Csr, Csr)>, up: Option<(Csr, Csr)> }
impl std::fmt::Debug for HodgeOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "HodgeOp(n={})", self.n) }
}
impl HodgeOp {
    /// vertex sheaf Laplacian L0 = d0^T d0
    pub fn l0(d0: &Csr) -> Self { HodgeOp { n: d0.ncols, down: None, up: Some((d0.clone_csr(), transpose(d0))) } }
    /// edge Hodge Laplacian of a graph (no 2-cells): L1 = d0 d0^T
    pub fn l1_graph(d0: &Csr) -> Self { HodgeOp { n: d0.nrows, down: Some((d0.clone_csr(), transpose(d0))), up: None } }
    /// general L_k from d_{k-1} (down) and d_k (up)
    pub fn general(n: usize, down: Option<&Csr>, up: Option<&Csr>) -> Self {
        HodgeOp { n, down: down.map(|a| (a.clone_csr(), transpose(a))), up: up.map(|b| (b.clone_csr(), transpose(b))) }
    }
    pub fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.iter_mut().for_each(|v| *v = 0.0);
        if let Some((a, at)) = &self.down {
            let mut t = vec![0.0; a.ncols]; par_mul(at, x, &mut t);
            let mut z = vec![0.0; self.n]; par_mul(a, &t, &mut z);
            y.iter_mut().zip(&z).for_each(|(o, v)| *o += v);
        }
        if let Some((b, bt)) = &self.up {
            let mut t = vec![0.0; b.nrows]; par_mul(b, x, &mut t);
            let mut z = vec![0.0; self.n]; par_mul(bt, &t, &mut z);
            y.iter_mut().zip(&z).for_each(|(o, v)| *o += v);
        }
    }
    /// ||L||_2 <= sum over terms of ||M||_1 ||M||_inf   (||M M^T||_2 = ||M||_2^2 <= ||M||_1 ||M||_inf)
    pub fn norm_bound(&self) -> f64 {
        let mut s = 0.0;
        for (m, _) in self.down.iter().chain(self.up.iter()) { let (a, b) = norm1_inf(m); s += a * b; }
        s
    }
    /// diag(L) from the CSR pairs: (A A^T)_ii = sum_j A_ij^2 ; (B^T B)_ii = sum_r B_ri^2
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.n];
        if let Some((a, _)) = &self.down { for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { d[i] += a.val[p] * a.val[p]; } } }
        if let Some((b, _)) = &self.up { for i in 0..b.nrows { for p in b.ptr[i]..b.ptr[i + 1] { d[b.idx[p]] += b.val[p] * b.val[p]; } } }
        d
    }
    pub fn dense(&self) -> Mat<f64> {
        let mut m = Mat::<f64>::zeros(self.n, self.n);
        let mut e = vec![0.0; self.n]; let mut y = vec![0.0; self.n];
        for j in 0..self.n { e[j] = 1.0; self.apply(&e, &mut y); for i in 0..self.n { m[(i, j)] = y[i]; } e[j] = 0.0; }
        m
    }
    pub fn residual(&self, v: &[f64], lam: f64) -> f64 {
        let mut y = vec![0.0; self.n]; self.apply(v, &mut y);
        y.iter().zip(v).map(|(a, b)| (a - lam * b).powi(2)).sum::<f64>().sqrt()
    }
}


/// Sheaf energy E(x) = ||delta_0 x||^2.
pub fn energy(d0: &Csr, x: &[f64]) -> f64 { let d = d0.mul(x); d.iter().map(|v| v * v).sum() }
/// Per-edge discord ||(delta_0 x)_e||^2, edge e owning rows `edge_offsets[e]..edge_offsets[e+1]`.
pub fn edge_discord(d0: &Csr, x: &[f64], edge_offsets: &[usize]) -> Vec<f64> {
    let d = d0.mul(x);
    edge_offsets.windows(2).map(|w| d[w[0]..w[1]].iter().map(|v| v * v).sum()).collect()
}
