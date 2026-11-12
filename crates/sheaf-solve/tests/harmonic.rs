//! Harmonic extension tests: known answer (path interpolation), dense == CG == faer dense solve, energy-minimality
//! proptest, and the B = ∅ reduction to H⁰ (orthogonal projection onto ker L; the kernel basis comes from an
//! independent faer self-adjoint EVD here, so sheaf-solve's tests do not depend on sheaf-cohom either).
use faer::{linalg::solvers::Solve, Mat, Side};
use proptest::prelude::*;
use sheaf_cell::{complex::CellComplex, fixtures::{cycle_graph, rotation, with_twist}, num::family::Exact, sheaf::Sheaf};
use sheaf_solve::*;

/// orthonormal basis of ker L0 = ker δ⁰ via EVD of δᵀδ (eigenvalues ≤ 1e-10·max(1, λ_max) count as zero)
fn kernel_basis(s: &Sheaf) -> Mat<f64> {
    let dm = s.coboundary_dense(0); let l = dm.transpose() * &dm; let n = l.nrows();
    let e = l.self_adjoint_eigen(Side::Lower).unwrap();
    let lam: Vec<f64> = (0..n).map(|i| e.S()[i]).collect();
    let top = lam.iter().cloned().fold(0.0, f64::max).max(1.0);
    let ks: Vec<usize> = (0..n).filter(|&i| lam[i] <= 1e-10 * top).collect();
    Mat::from_fn(n, ks.len(), |i, j| e.U()[(i, ks[j])])
}

fn delta(s: &Sheaf) -> Csr { let (r, c, t) = s.coboundary_triplets(0); Csr::from_triplets(r, c, t) }
fn path(n: usize) -> CellComplex { let e: Vec<_> = (0..n - 1).map(|i| (i, i + 1)).collect(); CellComplex::simplicial(n, &e, &[]) }
fn lcg(seed: &mut u64) -> f64 { *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((*seed >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0 }
fn random_sheaf(n: usize, extra: &[(usize, usize)], d: usize, seed: u64) -> Sheaf {
    let mut e: Vec<(usize, usize)> = (0..n - 1).map(|i| (i, i + 1)).collect(); e.extend_from_slice(extra);
    let mut s = Sheaf::constant(CellComplex::simplicial(n, &e, &[]), d);
    let mut sd = seed;
    for m in s.maps[0].iter_mut() { *m = Mat::from_fn(d, d, |_, _| lcg(&mut sd)); }
    s
}
const CG_ONLY: Options = Options { rtol: 1e-12, max_iter: 20_000, dense_max: 0 };

#[test]
fn constant_sheaf_on_path_is_linear_interpolation() {
    for d in [1usize, 2] {
        let n = 11; let s = Sheaf::constant(path(n), d); let dl = delta(&s);
        let off: Vec<usize> = (0..=n).map(|i| i * d).collect();
        let (a, b): (Vec<f64>, Vec<f64>) = ((0..d).map(|k| k as f64 + 1.0).collect(), (0..d).map(|k| -3.0 * k as f64 + 7.0).collect());
        let bd = node_boundary(&off, &[(0, a.clone()), (n - 1, b.clone())]);
        for o in [Options::default(), CG_ONLY] {
            let h = harmonic_extension(&dl, &bd, None, &o).unwrap();
            assert!(h.cert.converged, "{:?}", h.cert); assert_eq!(h.cert.method, Exact);
            for i in 0..n { for k in 0..d {
                let t = i as f64 / (n - 1) as f64; let want = (1.0 - t) * a[k] + t * b[k];
                assert!((h.x[i * d + k] - want).abs() < 1e-10, "d={d} node {i} comp {k}: {} vs {want} ({:?})", h.x[i * d + k], h.cert.route);
            } }
            // energy of the linear interpolant = |b-a|²/(n-1)
            let want_e: f64 = a.iter().zip(&b).map(|(x, y)| (y - x) * (y - x)).sum::<f64>() / (n - 1) as f64;
            assert!((h.cert.energy - want_e).abs() < 1e-10);
        }
    }
}

#[test]
fn dense_cg_and_faer_dense_solve_agree() {
    for seed in 1..6u64 {
        let (n, d) = (14, 2);
        let s = random_sheaf(n, &[(0, 7), (3, 11), (5, 13), (2, 9)], d, seed); let dl = delta(&s);
        let off: Vec<usize> = (0..=n).map(|i| i * d).collect();
        let bd = node_boundary(&off, &[(0, vec![1.0, -2.0]), (6, vec![0.5, 0.25]), (13, vec![-1.0, 3.0])]);
        let hd = harmonic_extension(&dl, &bd, None, &Options::default()).unwrap();
        let hc = harmonic_extension(&dl, &bd, None, &CG_ONLY).unwrap();
        assert_eq!(hd.cert.route, Route::Dense); assert_eq!(hc.cert.route, Route::Cg);
        assert!(hd.cert.converged && hc.cert.converged, "{:?} {:?}", hd.cert, hc.cert);
        // independent faer dense solve of L_II x_I = -L_IB x_B (LU with partial pivoting)
        let nd = n * d; let l = { let dm = s.coboundary_dense(0); dm.transpose() * &dm };
        let bset: Vec<usize> = bd.iter().map(|p| p.0).collect(); let int: Vec<usize> = (0..nd).filter(|i| !bset.contains(i)).collect();
        let lii = Mat::from_fn(int.len(), int.len(), |i, j| l[(int[i], int[j])]);
        let rhs = Mat::from_fn(int.len(), 1, |i, _| -bd.iter().map(|&(b, v)| l[(int[i], b)] * v).sum::<f64>());
        let xi = lii.partial_piv_lu().solve(&rhs);
        for (k, &i) in int.iter().enumerate() {
            assert!((hd.x[i] - xi[(k, 0)]).abs() < 1e-8 * (1.0 + xi[(k, 0)].abs()), "seed {seed} dof {i}: dense {} vs faer {}", hd.x[i], xi[(k, 0)]);
            assert!((hc.x[i] - xi[(k, 0)]).abs() < 1e-7 * (1.0 + xi[(k, 0)].abs()), "seed {seed} dof {i}: cg {} vs faer {}", hc.x[i], xi[(k, 0)]);
        }
        for &(b, v) in &bd { assert_eq!(hd.x[b], v); assert_eq!(hc.x[b], v); }
    }
}

#[test]
fn no_boundary_reduces_to_h0_projection() {
    // constant sheaf on a cycle: H0 = constants (dim d); twisted cycle (rotation): H0 = 0; random maps on a tree-with-chords
    let cases: Vec<(&str, Sheaf)> = vec![
        ("constant_cycle_d2", Sheaf::constant(cycle_graph(9), 2)),
        ("rotation_twist_cycle", with_twist(Sheaf::constant(cycle_graph(6), 2).complex, 0, 0, rotation(0.7))),
        ("constant_path_d3", Sheaf::constant(path(7), 3)),
        ("random_maps", random_sheaf(8, &[(0, 4)], 2, 42)),
    ];
    for (name, s) in cases {
        let dl = delta(&s); let n = dl.ncols; let mut sd = 9u64;
        let x0: Vec<f64> = (0..n).map(|_| lcg(&mut sd)).collect();
        let h = harmonic_extension(&dl, &[], Some(&x0), &CG_ONLY).unwrap();
        let q = kernel_basis(&s);   // orthonormal basis of ker L0
        let c: Vec<f64> = (0..q.ncols()).map(|j| (0..n).map(|i| q[(i, j)] * x0[i]).sum()).collect();
        let px: Vec<f64> = (0..n).map(|i| (0..q.ncols()).map(|j| q[(i, j)] * c[j]).sum()).collect();
        for (i, &pxi) in px.iter().enumerate() { assert!((h.x[i] - pxi).abs() < 1e-8, "{name} dof {i}: {} vs P_H0 x0 = {}", h.x[i], pxi); }
        assert!(h.cert.energy < 1e-16 * (1.0 + x0.iter().map(|v| v * v).sum::<f64>()), "{name}: energy {}", h.cert.energy);
        assert!(h.cert.converged, "{name}: {:?}", h.cert);
        // default route on B = ∅: if H0 != 0, L is singular, so Cholesky must refuse and fall back to CG; if H0 = 0
        // (twisted cycle) the dense solve is legitimate and gives 0. Same answer either way.
        let hd = harmonic_extension(&dl, &[], Some(&x0), &Options::default()).unwrap();
        assert_eq!(hd.cert.route, if q.ncols() > 0 { Route::Cg } else { Route::Dense }, "{name}");
        for (i, &pxi) in px.iter().enumerate() { assert!((hd.x[i] - pxi).abs() < 1e-8, "{name} dof {i} (default route)"); }
    }
}

#[test]
fn rejects_bad_boundary() {
    let dl = delta(&Sheaf::constant(path(3), 1));
    assert_eq!(harmonic_extension(&dl, &[(5, 1.0)], None, &Options::default()).unwrap_err(), HarmonicError::BoundaryIndex(5));
    assert_eq!(harmonic_extension(&dl, &[(1, 1.0), (1, 2.0)], None, &Options::default()).unwrap_err(), HarmonicError::DuplicateBoundary(1));
    assert_eq!(harmonic_extension(&dl, &[(1, f64::NAN)], None, &Options::default()).unwrap_err(), HarmonicError::NonFinite);
    let all = harmonic_extension(&dl, &[(0, 1.0), (1, 2.0), (2, 4.0)], None, &Options::default()).unwrap();
    assert_eq!(all.cert.route, Route::Trivial); assert_eq!(all.x, vec![1.0, 2.0, 4.0]); assert!((all.cert.energy - 5.0).abs() < 1e-15);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]
    /// The extension minimises the sheaf energy over all interior perturbations (boundary held fixed).
    #[test]
    fn extension_minimises_energy(seed in 1u64..10_000, n in 4usize..12, d in 1usize..3, nb in 1usize..3, eps in 1e-3f64..1.0, cg in any::<bool>()) {
        let s = random_sheaf(n, &[(0, n - 1)], d, seed); let dl = delta(&s);
        let off: Vec<usize> = (0..=n).map(|i| i * d).collect(); let mut sd = seed ^ 0x9e37;
        let bnodes: Vec<(usize, Vec<f64>)> = (0..nb).map(|k| (k * (n - 1) / nb.max(1), (0..d).map(|_| lcg(&mut sd)).collect())).collect();
        let bd = node_boundary(&off, &bnodes);
        let o = if cg { CG_ONLY } else { Options::default() };
        let h = harmonic_extension(&dl, &bd, None, &o).unwrap();
        prop_assert!(h.cert.converged, "{:?}", h.cert);
        let bset: Vec<usize> = bd.iter().map(|p| p.0).collect();
        for _ in 0..4 {
            let mut y = h.x.clone();
            for (i, yi) in y.iter_mut().enumerate() { if !bset.contains(&i) { *yi += eps * lcg(&mut sd); } }
            for &b in &bset { prop_assert_eq!(y[b].to_bits(), h.x[b].to_bits(), "boundary dof {} moved", b); }
            let ey = energy(&dl, &y);
            prop_assert!(ey >= h.cert.energy - 1e-9 * (1.0 + h.cert.energy), "perturbed {} < extension {}", ey, h.cert.energy);
        }
    }

    /// B = ∅ on a constant sheaf over a path: ker δ⁰ is the node-constant vectors, so P_H0 x0 is the per-coordinate
    /// mean over nodes. ∀ n in 2..10, d in 1..4, x0 in [-1, 1)^(n·d), on both routes: every dof (i, c) of the
    /// extension is within 1e-8 of mean_j x0[j·d + c], and the solve reports convergence through CG.
    #[test]
    fn no_boundary_on_constant_path_is_the_coordinate_mean(
        (n, d, x0) in (2usize..10, 1usize..4).prop_flat_map(|(n, d)| (Just(n), Just(d), proptest::collection::vec(-1.0f64..1.0, n * d))),
        cg in any::<bool>(),
    ) {
        let dl = delta(&Sheaf::constant(path(n), d));
        let o = if cg { CG_ONLY } else { Options::default() };
        let h = harmonic_extension(&dl, &[], Some(&x0), &o).unwrap();
        prop_assert!(h.cert.converged, "{:?}", h.cert);
        prop_assert_eq!(h.cert.route, Route::Cg);
        prop_assert_eq!(h.x.len(), n * d);
        let mean: Vec<f64> = (0..d).map(|c| (0..n).map(|j| x0[j * d + c]).sum::<f64>() / n as f64).collect();
        for (k, &xk) in h.x.iter().enumerate() { prop_assert!((xk - mean[k % d]).abs() < 1e-8, "dof {}: {} vs mean {}", k, xk, mean[k % d]); }
    }
}
