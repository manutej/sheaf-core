//! v1 spectral path vs dense EVD (fixtures + proptest), and v1 core behaviour.
use faer::{Mat, Side};
use proptest::prelude::*;
use sheaf_cell::{complex::CellComplex, fixtures::*, num::sparse::Csr, sheaf::Sheaf};
use sheaf_cohom::{cohomology::CochainComplex, spectral::*};

const TOL: f64 = 1e-8;

fn csr(s: &Sheaf, k: usize) -> Csr { let (r, c, t) = s.coboundary_triplets(k); Csr::from_triplets(r, c, t) }
/// independent reference: dense L_k from the spike's dense coboundaries, faer dense EVD
fn dense_eigs(cc: &CochainComplex, k: usize) -> Vec<f64> {
    let n = cc.dims[k];
    let mut l = Mat::<f64>::zeros(n, n);
    if k > 0 { let a = &cc.delta[k - 1]; l = &l + a * a.transpose(); }
    if k < cc.delta.len() { let b = &cc.delta[k]; l = &l + b.transpose() * b; }
    l.self_adjoint_eigenvalues(Side::Lower).unwrap()
}
fn op_for(s: &Sheaf, k: usize) -> HodgeOp {
    let top = s.complex.dim();
    let down = if k > 0 { Some(csr(s, k - 1)) } else { None };
    let up = if k < top { Some(csr(s, k)) } else { None };
    HodgeOp::general(s.cochain_dim(k), down.as_ref(), up.as_ref())
}
fn check(name: &str, op: &HodgeOp, known: &KnownKernel, reference: &[f64], k: usize, dense_max: usize) -> Spectral {
    let o = SpectralOpts { k, dense_max, ..Default::default() };
    let sp = smallest(op, known, &o);
    let scale = reference.last().cloned().unwrap_or(1.0).abs().max(1.0);
    let mut all: Vec<f64> = vec![0.0; sp.known_kernel_dim]; all.extend(sp.pairs.iter().map(|p| p.value));
    for (i, (a, b)) in all.iter().zip(reference).take(k + 1).enumerate() {
        assert!((a - b).abs() <= TOL * scale, "{name}: eig {i}: got {a:e} dense {b:e} ({})", sp.method);
    }
    for p in &sp.pairs { assert!(p.residual <= 1e-8 * scale, "{name}: residual {:e}", p.residual); }
    sp
}

#[test]
fn catalogue_fixtures_l0_and_l1_match_dense_evd_and_h0() {
    for f in catalogue() {
        let cc = CochainComplex::from_sheaf(&f.sheaf).unwrap();
        for k in 0..cc.dims.len().min(2) {
            let op = op_for(&f.sheaf, k);
            let reference = dense_eigs(&cc, k);
            let kk = cc.dims[k].saturating_sub(1).max(1);
            let sp = check(f.name, &op, &KnownKernel::none(op.n), &reference, kk.min(6), 1500);
            let want = f.expect_r[k];
            if !sp.h0_lower_bound_only { assert_eq!(sp.h0_dim(), want, "{}: dim ker L{k}", f.name); }
            assert!(sp.cut_consistent_with_residuals() || sp.h0_lower_bound_only, "{}: cut not certified", f.name);
        }
    }
}

fn graph_sheaf(n: usize, extra: &[(usize, usize)], d: usize, seed: u64, random_maps: bool, twist: bool) -> Sheaf {
    let mut e: Vec<(usize, usize)> = (0..n).map(|i| (i, (i + 1) % n)).collect();
    e.extend_from_slice(extra);
    let mut s = Sheaf::constant(CellComplex::simplicial(n, &e, &[]), d);
    let mut st = seed | 1;
    if random_maps { for m in s.maps[0].iter_mut() { *m = Mat::from_fn(d, d, |_, _| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 }); } }
    if twist { s.maps[0][1] = rotation(std::f64::consts::PI / 3.0); }
    s
}

#[test]
fn lanczos_matches_dense_on_larger_fixtures() {
    // (name, sheaf, use known kernel)
    let cases = vec![
        ("cycle200_const1_known", graph_sheaf(200, &[], 1, 1, false, false), true),
        ("cycle200_const1_unknown", graph_sheaf(200, &[], 1, 1, false, false), false),
        ("cycle150_const2_twist60", graph_sheaf(150, &[(0, 75), (10, 90)], 2, 1, false, true), false),
        ("graph90_random_maps_d2", graph_sheaf(90, &[(0, 45), (3, 60), (7, 33), (20, 80)], 2, 7, true, false), false),
        ("cycle100_const3_known", graph_sheaf(100, &[(5, 50)], 3, 1, false, false), true),
    ];
    for (name, s, use_known) in cases {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        let op = op_for(&s, 0);
        let known = if use_known {
            let d = s.stalk[0][0]; let n = s.stalk[0].len();
            let off: Vec<usize> = (0..=n).map(|i| i * d).collect();
            KnownKernel::from_components(&vec![Some(0); n], &off)
        } else { KnownKernel::none(op.n) };
        let sp = check(name, &op, &known, &dense_eigs(&cc, 0), 6, 0);
        assert!(sp.method.starts_with("lanczos"), "{name} used {}", sp.method);
        assert_eq!(sp.h0_dim(), cc.cohomology(0).unwrap().dim, "{name} h0");
    }
}

#[test]
fn l1_graph_kernel_multiplicity_found_by_lock_and_restart() {
    // tree-like cycle + 2 chords: b1 = 3 => L1 = d0 d0^T has a 3-fold zero eigenvalue the solver must find unaided
    let s = graph_sheaf(160, &[(0, 80), (20, 120)], 1, 1, false, false);
    let cc = CochainComplex::from_sheaf(&s).unwrap();
    let op = HodgeOp::l1_graph(&csr(&s, 0));
    let sp = check("l1_cycle160_chords2", &op, &KnownKernel::none(op.n), &dense_eigs(&cc, 1), 6, 0);
    assert_eq!(sp.found_zero, 3);
    assert!(sp.cut_consistent_with_residuals());
}

#[test]
fn polynomial_method_agrees_on_expander_like_graph() {
    let mut extra = vec![]; for i in 0..200 { extra.push((i, (i * 37 + 11) % 200)); }
    let extra: Vec<(usize, usize)> = extra.into_iter().filter(|(a, b)| a != b).collect();
    let s = graph_sheaf(200, &extra, 1, 1, false, false);
    let cc = CochainComplex::from_sheaf(&s).unwrap();
    let op = op_for(&s, 0);
    let o = SpectralOpts { k: 4, dense_max: 0, method: Method::Polynomial, ..Default::default() };
    let sp = smallest(&op, &KnownKernel::none(op.n), &o);
    let reference = dense_eigs(&cc, 0);
    for (a, b) in sp.pairs.iter().zip(&reference) { assert!((a.value - b).abs() <= TOL * reference.last().unwrap(), "{} vs {}", a.value, b); }
}

fn arb_graph_sheaf() -> impl Strategy<Value = (Sheaf, bool)> {
    (70usize..110, prop::collection::vec((0usize..110, 0usize..110), 0..12), 1usize..3, any::<u64>(), any::<bool>(), any::<bool>()).prop_map(|(n, ch, d, seed, rnd, tw)| {
        let ch: Vec<(usize, usize)> = ch.into_iter().map(|(a, b)| (a % n, b % n)).filter(|(a, b)| a != b).collect();
        let d = if n * d <= 135 { 2 } else { d }; // keep dofs above the dense fallback (> 2*64+2)
        let s = graph_sheaf(n, &ch, d, seed, rnd, tw && d == 2 && !rnd);
        (s, rnd)
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]
    #[test]
    fn lanczos_vs_dense_random_graph_sheaves((s, _rnd) in arb_graph_sheaf()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        let op = op_for(&s, 0);
        prop_assume!(op.n > 140);
        let reference = dense_eigs(&cc, 0);
        let o = SpectralOpts { k: 5, dense_max: 0, ..Default::default() };
        let sp = smallest(&op, &KnownKernel::none(op.n), &o);
        let scale = reference.last().unwrap().max(1.0);
        for (a, b) in sp.pairs.iter().zip(&reference) { prop_assert!((a.value - b).abs() <= TOL * scale, "{:e} vs {:e}", a.value, b); }
        if !sp.h0_lower_bound_only { prop_assert_eq!(sp.h0_dim(), cc.cohomology(0).unwrap().dim); }
    }

    /// ∀ random graph sheaves with more than 140 dofs (so the Lanczos path runs, not the dense fallback):
    /// two `smallest` calls with equal (default shift-invert) options return bit-identical spectra
    /// (values, residuals, vectors, rounds, inner stats, flags).
    #[test]
    fn smallest_is_bitwise_deterministic((s, _rnd) in arb_graph_sheaf()) {
        let op = op_for(&s, 0);
        prop_assume!(op.n > 140);
        let o = SpectralOpts { k: 5, dense_max: 0, ..Default::default() };
        let known = KnownKernel::none(op.n);
        prop_assert_eq!(fingerprint(&smallest(&op, &known, &o)), fingerprint(&smallest(&op, &known, &o)));
    }
}

/// Every numeric field of a `Spectral` as raw bits, so equality is exact (no tolerance, NaN-safe).
fn fingerprint(sp: &Spectral) -> (Vec<u64>, &'static str, usize, bool, bool) {
    let mut b: Vec<u64> = sp.pairs.iter().flat_map(|p| [p.value.to_bits(), p.residual.to_bits()]).collect();
    b.extend(sp.vectors.iter().flatten().map(|x| x.to_bits()));
    b.extend([sp.c_shift, sp.lambda_max_bound, sp.zero_tol, sp.sigma].map(f64::to_bits));
    if let Some((solves, iters, worst)) = sp.inner { b.extend([solves as u64, iters as u64, worst.to_bits()]); }
    b.extend([sp.n, sp.known_kernel_dim, sp.found_zero].map(|x| x as u64));
    (b, sp.method, sp.rounds, sp.converged, sp.h0_lower_bound_only)
}
