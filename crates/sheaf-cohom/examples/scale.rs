//! Matrix-free graph Hodge split at scale (v1 target). Prints timings; numbers are box-local, single run.
use sheaf_cell::{complex::CellComplex, num::sparse::*, sheaf::Sheaf};
use sheaf_cohom::hodge_graph::*;
use std::time::Instant;
fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(100_000);
    let extra: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(2 * n);
    let d: usize = std::env::args().nth(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    // cycle + `extra` pseudo-random chords => b1 = |E| - |V| + 1 (connected), constant R^d sheaf
    let mut e: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    let mut st: u64 = 0x2545F4914F6CDD1D;
    for _ in 0..extra { st ^= st << 13; st ^= st >> 7; st ^= st << 17; let a = (st % n as u64) as usize; st ^= st << 13; st ^= st >> 7; st ^= st << 17; let b = (st % n as u64) as usize; if a != b { e.push((a, b)); } }
    let t0 = Instant::now();
    let cx = CellComplex::simplicial(n, &e, &[]);
    let s = Sheaf::constant(cx, d);
    let (r, c, t) = s.coboundary_triplets(0);
    let d0 = Csr::from_triplets(r, c, t);
    let t_build = t0.elapsed();
    let w: Vec<f64> = (0..r).map(|i| ((i as f64) * 0.618).sin()).collect();
    let t1 = Instant::now();
    let g = graph_hodge_cg(&d0, &w, 1e-10, 20_000);
    let t_cg = t1.elapsed();
    let ha: f64 = g.harmonic.iter().map(|x| x * x).sum::<f64>();
    let tot: f64 = w.iter().map(|x| x * x).sum::<f64>();
    let ne = s.complex.n_cells[1];
    println!("V={} E={} d={} dimC0={} dimC1={} b1(if connected)={}", n, ne, d, c, r, d * (ne + 1 - n));
    println!("build {:.3}s  cg {:.3}s  iters {}  harmonic_frac {:.4}  random-null b1/E={:.4}  cert ||d0^T h||/||w|| {:.2e} converged {}",
        t_build.as_secs_f64(), t_cg.as_secs_f64(), g.iters, ha / tot, (ne + 1 - n) as f64 / ne as f64, g.cert_coclosed, g.converged);
}
