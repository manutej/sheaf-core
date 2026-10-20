use sheaf_cell::{complex::CellComplex, fixtures::rotation, num::sparse::Csr, sheaf::Sheaf};
use sheaf_cohom::spectral::*;
fn main() {
    // cycle(300) constant R^2 with a 60deg twist + random chords: no known kernel
    let n = 300; let mut e: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    for i in 0..40 { e.push((i * 7 % n, (i * 13 + 50) % n)); }
    let cx = CellComplex::simplicial(n, &e, &[]);
    for twist in [false, true] {
        let mut s = Sheaf::constant(cx.clone(), 2);
        if twist { s.maps[0][1] = rotation(std::f64::consts::PI / 3.0); }
        let (r, c, t) = s.coboundary_triplets(0); let d0 = Csr::from_triplets(r, c, t);
        let op = HodgeOp::l0(&d0); let kk = KnownKernel::none(op.n);
        let o = SpectralOpts { k: 6, dense_max: 0, ..Default::default() };
        let t0 = std::time::Instant::now();
        let a = smallest(&op, &kk, &o); let ta = t0.elapsed();
        let b = smallest_dense(&op, &kk, 6, H0Policy::adr012());
        println!("twist={twist} rounds={} conv={} t={:?} h0 lz={} dense={}", a.rounds, a.converged, ta, a.h0_dim(), b.h0_dim());
        for (x, y) in a.pairs.iter().zip(&b.pairs) { println!("  {:.15e} {:.15e} diff {:.2e} res {:.2e}", x.value, y.value, (x.value - y.value).abs(), x.residual); }
    }
}
