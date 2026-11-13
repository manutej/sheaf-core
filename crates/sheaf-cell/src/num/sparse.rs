//! Minimal CSR (row-major sparse) used for coboundaries. Duplicate (row, col) triplets are KEPT as separate entries
//! (products sum them), so multi-edges never merge: each edge owns its rows anyway (see `graph::coboundary`).
#[derive(Debug, Clone)]
pub struct Csr { pub nrows: usize, pub ncols: usize, pub ptr: Vec<usize>, pub idx: Vec<usize>, pub val: Vec<f64> }

impl Csr {
    pub fn clone_csr(&self) -> Csr { Csr { nrows: self.nrows, ncols: self.ncols, ptr: self.ptr.clone(), idx: self.idx.clone(), val: self.val.clone() } }
    pub fn from_triplets(nrows: usize, ncols: usize, mut t: Vec<(usize, usize, f64)>) -> Self {
        t.sort_by_key(|x| (x.0, x.1));
        let mut ptr = vec![0; nrows + 1];
        for &(i, _, _) in &t { ptr[i + 1] += 1; }
        for i in 0..nrows { ptr[i + 1] += ptr[i]; }
        Csr { nrows, ncols, ptr, idx: t.iter().map(|x| x.1).collect(), val: t.iter().map(|x| x.2).collect() }
    }
    pub fn mul(&self, x: &[f64]) -> Vec<f64> {
        (0..self.nrows).map(|i| (self.ptr[i]..self.ptr[i + 1]).map(|p| self.val[p] * x[self.idx[p]]).sum()).collect()
    }
    pub fn mul_t(&self, y: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.ncols];
        for (i, w) in self.ptr.windows(2).take(self.nrows).enumerate() { for p in w[0]..w[1] { out[self.idx[p]] += self.val[p] * y[i]; } }
        out
    }
}


/// transpose (CSR of A^T)
pub fn transpose(a: &Csr) -> Csr {
    let mut t = Vec::with_capacity(a.val.len());
    for i in 0..a.nrows { for p in a.ptr[i]..a.ptr[i + 1] { t.push((a.idx[p], i, a.val[p])); } }
    Csr::from_triplets(a.ncols, a.nrows, t)
}
/// deterministic row-parallel y = A x (each output row summed sequentially => bitwise reproducible)
pub fn par_mul(a: &Csr, x: &[f64], y: &mut [f64]) {
    use rayon::prelude::*;
    y.par_iter_mut().enumerate().with_min_len(4096).for_each(|(i, o)| {
        let mut s = 0.0; for p in a.ptr[i]..a.ptr[i + 1] { s += a.val[p] * x[a.idx[p]]; } *o = s;
    });
}
/// (||A||_1, ||A||_inf)
pub fn norm1_inf(a: &Csr) -> (f64, f64) {
    let mut col = vec![0.0; a.ncols]; let mut rmax = 0.0f64;
    for i in 0..a.nrows { let mut r = 0.0; for p in a.ptr[i]..a.ptr[i + 1] { r += a.val[p].abs(); col[a.idx[p]] += a.val[p].abs(); } rmax = rmax.max(r); }
    (col.into_iter().fold(0.0, f64::max), rmax)
}
