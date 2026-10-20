//! ADR-017 adversarial panel regression tests (sheaf-cohom).
use sheaf_cell::m1::{coboundary, laplacian, M1Edge};
use sheaf_cell::num::tol::H0Policy;
use sheaf_cohom::m1;

fn lcg(s: &mut u64) -> f64 { *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((*s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0 }

/// P-CONC-01: the m1 H⁰ routes must be bitwise identical whatever the rayon pool size (they used faer's global
/// parallelism, which made the kernel basis G differ by up to 1.66 entrywise between 1, 2 and 8 threads).
#[test]
fn m1_h0_bitwise_independent_of_thread_count() {
    let (n, d, de) = (300usize, 4usize, 2usize);
    let mut s = 5u64; let mut maps = Vec::new(); let mut ends = Vec::new();
    for i in 0..n { let j = (i + 1 + ((lcg(&mut s).abs() * 40.0) as usize)) % n; ends.push((i, j));
        maps.push(((0..de * d).map(|_| lcg(&mut s)).collect::<Vec<_>>(), (0..de * d).map(|_| lcg(&mut s)).collect::<Vec<_>>())); }
    let dv = vec![d; n];
    let edges: Vec<M1Edge> = ends.iter().zip(&maps).map(|(&(u, v), (fu, fv))| M1Edge { u, v, de, fu, fv }).collect();
    let dl = coboundary(&dv, &edges); let l = laplacian(&dl);
    let run = |threads: usize| rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| {
        let g = m1::global_sections(&dl, H0Policy::m1_svd(1e-9)).unwrap();
        let comps = m1::components(&dl);
        let sv: Vec<u64> = comps.iter().flat_map(|c| m1::block_singular_values(&dl, c).unwrap()).map(f64::to_bits).collect();
        let ev: Vec<u64> = m1::l_components(&l).iter().flat_map(|c| m1::block_eigenvalues(&l, c).unwrap()).map(f64::to_bits).collect();
        let gb: Vec<u64> = (0..g.nrows()).flat_map(|i| (0..g.ncols()).map(move |j| (i, j))).map(|(i, j)| g[(i, j)].to_bits()).collect();
        (gb, sv, ev, m1::h0_dim_from_l(&l, H0Policy::m1_eig(1e-10)).unwrap())
    });
    let base = run(1);
    for t in [2, 8] { let r = run(t); assert!(r.0 == base.0, "G differs at {t} threads"); assert!(r.1 == base.1, "sigma differs at {t} threads");
        assert!(r.2 == base.2, "lambda differs at {t} threads"); assert_eq!(r.3, base.3); }
}
