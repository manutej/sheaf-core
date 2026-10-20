//! Max |lambda_lanczos - lambda_dense| / max(1, lambda_max) and max residual over fixtures + 200 seeded random graph sheaves.
use faer::{Mat, Side};
use sheaf_cell::{complex::CellComplex, fixtures::*, num::sparse::Csr, sheaf::Sheaf};
use sheaf_cohom::{cohomology::CochainComplex, spectral::*};
fn csr(s: &Sheaf, k: usize) -> Csr { let (r, c, t) = s.coboundary_triplets(k); Csr::from_triplets(r, c, t) }
fn dense_eigs(cc: &CochainComplex, k: usize) -> Vec<f64> {
    let n = cc.dims[k]; let mut l = Mat::<f64>::zeros(n, n);
    if k > 0 { let a = &cc.delta[k - 1]; l = &l + a * a.transpose(); }
    if k < cc.delta.len() { let b = &cc.delta[k]; l = &l + b.transpose() * b; }
    l.self_adjoint_eigenvalues(Side::Lower).unwrap()
}
fn main() {
    let mut worst_fix = 0.0f64; let mut worst_res_fix = 0.0f64; let mut h0_ok = 0; let mut h0_n = 0;
    for f in catalogue() { // dense path (tiny n) vs independent dense EVD
        let cc = CochainComplex::from_sheaf(&f.sheaf).unwrap();
        for k in 0..cc.dims.len().min(2) {
            let top = f.sheaf.complex.dim();
            let down = if k > 0 { Some(csr(&f.sheaf, k - 1)) } else { None }; let up = if k < top { Some(csr(&f.sheaf, k)) } else { None };
            let op = HodgeOp::general(f.sheaf.cochain_dim(k), down.as_ref(), up.as_ref());
            let r = dense_eigs(&cc, k); let sp = smallest(&op, &KnownKernel::none(op.n), &SpectralOpts { k: (op.n - 1).min(6).max(1), ..Default::default() });
            let sc = r.last().unwrap().max(1.0);
            for (a, b) in sp.pairs.iter().zip(&r) { worst_fix = worst_fix.max((a.value - b).abs() / sc); worst_res_fix = worst_res_fix.max(a.residual / sc); }
            h0_n += 1; if sp.h0_dim() == f.expect_r[k] { h0_ok += 1; }
        }
    }
    let mut worst = 0.0f64; let mut worst_res = 0.0f64; let mut cases = 0; let mut h0m = 0; let mut st: u64 = 0x1234567;
    let mut rnd = |m: u64| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; st % m };
    for _ in 0..200 {
        let n = 70 + rnd(50) as usize; let d = 2;
        let mut e: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
        for _ in 0..rnd(12) { let a = rnd(n as u64) as usize; let b = rnd(n as u64) as usize; if a != b { e.push((a, b)); } }
        let mut s = Sheaf::constant(CellComplex::simplicial(n, &e, &[]), d);
        match rnd(3) { 0 => {}, 1 => { s.maps[0][1] = rotation(1.0 + rnd(100) as f64 / 50.0); },
            _ => { let mut q = rnd(u64::MAX) | 1; for m in s.maps[0].iter_mut() { *m = Mat::from_fn(d, d, |_, _| { q ^= q << 13; q ^= q >> 7; q ^= q << 17; ((q % 2001) as f64) / 1000.0 - 1.0 }); } } }
        let cc = CochainComplex::from_sheaf(&s).unwrap(); let r = dense_eigs(&cc, 0);
        let op = HodgeOp::l0(&csr(&s, 0));
        let sp = smallest(&op, &KnownKernel::none(op.n), &SpectralOpts { k: 6, dense_max: 0, ..Default::default() });
        assert!(sp.method.starts_with("lanczos"));
        let sc = r.last().unwrap().max(1.0);
        for (a, b) in sp.pairs.iter().zip(&r) { worst = worst.max((a.value - b).abs() / sc); worst_res = worst_res.max(a.residual / sc); }
        if sp.h0_dim() != cc.cohomology(0).unwrap().dim { h0m += 1; let h = cc.cohomology(0).unwrap(); eprintln!("H0 MISMATCH lanczos {} (tol {:e}) cohomology(0) {} ; smallest eigs {:?} ; svd kept_min {:?} dropped_max {:?} svd tol {:e}", sp.h0_dim(), sp.zero_tol, h.dim, sp.pairs.iter().take(3).map(|p| p.value).collect::<Vec<_>>(), h.rank_out.sigma_kept_min, h.rank_out.sigma_dropped_max, h.rank_out.tol); }
        cases += 1;
    }
    let out = serde_json::json!({"fixtures_dense_path": {"max_rel_eig_err": worst_fix, "max_rel_residual": worst_res_fix, "h0_matches": format!("{h0_ok}/{h0_n}")},
        "random_lanczos_path": {"cases": cases, "k_plus_1": 7, "max_rel_eig_err": worst, "max_rel_residual": worst_res, "h0_mismatches_lanczos_vs_cohomology0_same_policy": h0m}});
    std::fs::write("out/spectral_vs_dense.json", serde_json::to_string_pretty(&out).unwrap()).unwrap();
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
