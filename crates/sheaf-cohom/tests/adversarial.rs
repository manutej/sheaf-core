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

// ---- ADR-017 round 2 ----
use sheaf_cell::num::sparse::Csr;
use sheaf_cohom::{fp, relative};

fn max_resid(d: &Csr, b: &faer::Mat<f64>) -> f64 {
    (0..b.ncols()).map(|j| { let x: Vec<f64> = (0..b.nrows()).map(|i| b[(i, j)]).collect(); d.mul(&x).iter().fold(0.0f64, |a, v| a.max(v.abs())) }).fold(0.0, f64::max)
}

/// P-NUM-02: explicit stored zeros must not bind a row to the wrong interior block (it returned g = e1 with δg = 1).
#[test]
fn relative_h0_ignores_explicit_zeros() {
    let d = Csr { nrows: 1, ncols: 2, ptr: vec![0, 2], idx: vec![0, 1], val: vec![0.0, 1.0] };
    let r = relative::h0_relative(&d, &[], H0Policy::adr012(), 1500).unwrap().value;
    assert_eq!(r.dim, 1); assert!(max_resid(&d, &r.basis) < 1e-14, "basis not in ker δ");
    assert!((r.basis[(0, 0)].abs() - 1.0).abs() < 1e-14, "kernel must be e0");
    // three columns, zero first: used to index past a width-1 block and panic
    let d3 = Csr { nrows: 1, ncols: 3, ptr: vec![0, 3], idx: vec![0, 1, 2], val: vec![0.0, 1.0, 1.0] };
    let r3 = relative::h0_relative(&d3, &[], H0Policy::adr012(), 1500).unwrap().value;
    assert_eq!(r3.dim, 2); assert!(max_resid(&d3, &r3.basis) < 1e-14);
    // same with a boundary column and a zero in the interior
    let db = Csr { nrows: 2, ncols: 3, ptr: vec![0, 2, 4], idx: vec![0, 1, 1, 2], val: vec![0.0, 1.0, 1.0, -1.0] };
    let rb = relative::h0_relative(&db, &[1], H0Policy::adr012(), 1500).unwrap().value;
    assert_eq!(rb.dim, 1); assert!(max_resid(&db, &rb.basis) < 1e-14);
}

/// Audit (P-NUM-02): the m1 H⁰ routes give the same answer with and without explicitly stored zeros.
#[test]
fn m1_h0_invariant_under_explicit_zeros() {
    let clean = Csr { nrows: 2, ncols: 3, ptr: vec![0, 2, 3], idx: vec![0, 1, 2], val: vec![1.0, -1.0, 2.0] };
    let dirty = Csr { nrows: 2, ncols: 3, ptr: vec![0, 3, 5], idx: vec![0, 1, 2, 0, 2], val: vec![1.0, -1.0, 0.0, 0.0, 2.0] };
    let g = |d: &Csr| { let g = m1::global_sections(d, H0Policy::m1_svd(1e-9)).unwrap(); (g.ncols(), max_resid(d, &g)) };
    let (a, b) = (g(&clean), g(&dirty)); assert_eq!(a.0, b.0); assert!(b.1 < 1e-14);
    let l = sheaf_cell::m1::laplacian(&dirty);   // keeps structural zeros
    assert_eq!(m1::h0_dim_from_l(&l, H0Policy::m1_eig(1e-10)).unwrap(), a.0);
}

/// P-NUM-03: F_p elimination used i64 products and overflowed for p > ~3.04e9 (Betti [0,0] instead of [1,1]).
#[test]
fn betti_mod_p_large_primes() {
    let d = faer::Mat::from_fn(2, 2, |i, j| if i == j { 1.0 } else { -1.0 });
    assert_eq!(fp::betti_mod_p(&[d.clone()], &[2, 2], 4294967291), Some(vec![1, 1]));
    let pmax = 9223372036854775783i64; // largest prime below 2^63 (the cap: p is an i64)
    assert!(fp::is_prime(pmax));
    assert_eq!(fp::betti_mod_p(&[d.clone()], &[2, 2], pmax), Some(vec![1, 1]));
    let e = faer::Mat::from_fn(2, 2, |i, j| [[3.0, 5.0], [7.0, 11.0]][i][j]);   // det = -2: full rank unless p = 2
    assert_eq!(fp::betti_mod_p(&[e.clone()], &[2, 2], pmax), Some(vec![0, 0]));
    assert_eq!(fp::betti_mod_p(&[e], &[2, 2], 2), Some(vec![1, 1]));
    assert_eq!(fp::betti_mod_p_checked(&[d.clone()], &[2, 2], 4), Err(fp::FpError::NotPrime(4)));
    assert_eq!(fp::betti_mod_p_checked(&[d], &[2, 2], 1), Err(fp::FpError::NotPrime(1)));
}

/// P-NUM-04: to_int accepted NaN (as 0), ±inf and 2^63 (saturated to i64::MAX).
#[test]
fn to_int_refuses_nonfinite_and_inexact() {
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 9223372036854775808.0, 9007199254740994.0] {
        assert!(fp::to_int(&faer::Mat::from_fn(1, 1, |_, _| v)).is_none(), "{v} accepted");
        assert_eq!(fp::betti_mod_p_checked(&[faer::Mat::from_fn(1, 1, |_, _| v)], &[1, 1], 2), Err(fp::FpError::NotIntegral));
    }
    assert_eq!(fp::to_int(&faer::Mat::from_fn(1, 1, |_, _| 9007199254740992.0)), Some(vec![vec![9007199254740992i64]]));
}
