use proptest::prelude::*;
use sheaf_cell::{complex::CellComplex, num::linalg::*, sheaf::Sheaf};
use sheaf_cohom::cohomology::*;
use faer::Mat;
use faer::linalg::solvers::DenseSolveCore;

/// random simplicial 2-complex on n vertices: random triangles + random extra edges
fn arb_complex() -> impl Strategy<Value = CellComplex> {
    (3usize..8).prop_flat_map(|n| {
        let tri = prop::collection::vec((0..n, 0..n, 0..n), 0..6);
        let ed = prop::collection::vec((0..n, 0..n), 0..8);
        (Just(n), tri, ed)
    }).prop_map(|(n, tri, ed)| {
        let faces: Vec<[usize; 3]> = tri.into_iter().filter(|&(a, b, c)| a != b && b != c && a != c).map(|(a, b, c)| [a, b, c]).collect();
        let edges: Vec<(usize, usize)> = ed.into_iter().filter(|&(a, b)| a != b).collect();
        CellComplex::simplicial(n, &edges, &faces)
    })
}

/// random sheaf that is guaranteed functorial: vertex/edge stalks R^d with random edge maps on graph
/// part, and face stalks = 0 or (if the complex has faces) a constant sheaf. Two families:
fn arb_sheaf() -> impl Strategy<Value = Sheaf> {
    (arb_complex(), 1usize..3, any::<u64>(), any::<bool>()).prop_map(|(cx, d, seed, twist)| {
        let has_faces = cx.n_cells.len() > 2;
        let mut s = Sheaf::constant(cx, d);
        if twist && !has_faces {
            // graph only: arbitrary per-incidence maps are always a valid sheaf (no codim-2 condition)
            let mut st = seed | 1;
            for m in s.maps[0].iter_mut() {
                *m = Mat::from_fn(d, d, |_, _| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 });
            }
        }
        s
    })
}

fn arb_cochain(n: usize, seed: u64) -> Vec<f64> {
    let mut st = seed | 1;
    (0..n).map(|_| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 20001) as f64) / 10000.0 - 1.0 }).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    #[test]
    fn d_squared_is_zero(s in arb_sheaf()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        prop_assert!(cc.d2_residual <= D2_TOL);
    }

    // NOTE (Grok r1 #4): rank-nullity and Euler are identities given how dim is computed; the non-trivial
    // assertion here is `consistent()` (third, independent SVD of [d_k; d_{k-1}^T]). Golden values live in known.rs.
    #[test]
    fn rank_nullity_and_euler(s in arb_sheaf()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        let mut chi_c = 0i64; let mut chi_h = 0i64;
        for k in 0..cc.dims.len() {
            let h = cc.cohomology(k).unwrap();
            prop_assert!(h.consistent(), "dim {} vs basis {}", h.dim, h.harmonic_basis.ncols());
            // rank-nullity: dim ker d_k = dim C^k - rank d_k, and dim H^k = dim ker d_k - rank d_{k-1}
            prop_assert_eq!(h.dim + h.rank_in.rank + h.rank_out.rank, cc.dims[k]);
            let sgn = if k % 2 == 0 { 1 } else { -1 };
            chi_c += sgn * cc.dims[k] as i64; chi_h += sgn * h.dim as i64;
        }
        prop_assert_eq!(chi_c, chi_h);
    }

    #[test]
    fn hodge_components_orthogonal_and_sum(s in arb_sheaf(), seed in any::<u64>()) {
        let cc = CochainComplex::from_sheaf(&s).unwrap();
        for k in 0..cc.dims.len() {
            if cc.dims[k] == 0 { continue; }
            let w = arb_cochain(cc.dims[k], seed);
            let h = cc.hodge(k, &w);
            // ADR-015: delta_0 is cut by the H0 policy (sigma^2 <= 1e-10 * ||L0|_c||), so a dropped sigma (up to ~1e-5 *
            // sqrt(bound)) legitimately leaves (closed / coclosed) residuals of that size in the delta_0-facing checks;
            // every other identity must hold to 1e-9.
            let sd = if cc.delta.is_empty() { 0.0 } else { cc.d0_split().rank.sigma_dropped_max.unwrap_or(0.0) };
            let c = &h.cert;
            let d0_tol = 1e-9f64.max(1.01 * sd);
            let (closed_tol, coclosed_tol) = match k { 0 => (d0_tol, 1e-9), 1 => (1e-9, d0_tol), _ => (1e-9, 1e-9) };
            prop_assert!(c.potential_residual.max(c.ortho_exact_harm).max(c.ortho_exact_coex).max(c.ortho_harm_coex) < 1e-9
                && c.harmonic_closed <= closed_tol && c.harmonic_coclosed <= coclosed_tol, "k={} sigma_dropped={} cert={:?}", k, sd, h.cert);
            let (a, b, c) = h.energies();
            prop_assert!((a + b + c - dot(&w, &w)).abs() < 1e-9 * dot(&w, &w).max(1.0)); // Pythagoras
            // harmonic part lies in span of the harmonic basis
            let hb = cc.cohomology(k).unwrap().harmonic_basis;
            let p = project(&hb, &h.harmonic);
            prop_assert!(norm(&p.iter().zip(&h.harmonic).map(|(x, y)| x - y).collect::<Vec<_>>()) < 1e-8);
        }
    }
}

/// Gauge-twisted sheaf on an arbitrary 2-complex: F_{s<=t} = A_t A_s^{-1} with random invertible A_c.
/// Always functorial (d^2 = 0) and isomorphic to the constant sheaf, so its Betti numbers must match.
fn gauge(cx: &CellComplex, d: usize, seed: u64) -> Sheaf {
    let mut st = seed | 1;
    let mut rnd = move || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 };
    let a: Vec<Vec<Mat<f64>>> = cx.n_cells.iter().map(|&n| (0..n).map(|_| {
        Mat::from_fn(d, d, |i, j| rnd() + if i == j { 3.0 } else { 0.0 })
    }).collect()).collect();
    let mut s = Sheaf::constant(cx.clone(), d);
    for k in 0..cx.inc.len() {
        for (j, inc) in cx.inc[k].iter().enumerate() {
            let inv = a[k][inc.facet].partial_piv_lu().inverse();
            s.maps[k][j] = &a[k + 1][inc.cofacet] * &inv;
        }
    }
    s
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]
    #[test]
    fn gauge_twist_on_2_complex_preserves_complex_and_betti(cx in arb_complex(), d in 1usize..3, seed in any::<u64>()) {
        let base = CochainComplex::from_sheaf(&Sheaf::constant(cx.clone(), d)).unwrap();
        let tw = CochainComplex::from_sheaf(&gauge(&cx, d, seed)).unwrap();
        prop_assert!(tw.d2_residual <= D2_TOL);
        for k in 0..base.dims.len() {
            prop_assert_eq!(base.cohomology(k).unwrap().dim, tw.cohomology(k).unwrap().dim);
            let w = arb_cochain(tw.dims[k], seed ^ 0x9e37);
            if tw.dims[k] > 0 { prop_assert!(tw.hodge(k, &w).cert.max() < 1e-8); }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]
    /// Grok r1 #5: generic random maps on a complex WITH faces violate functoriality and must be rejected.
    #[test]
    fn random_maps_on_faces_are_rejected(seed in any::<u64>()) {
        use sheaf_cell::fixtures::{triangle_filled, sphere};
        for cx in [triangle_filled(), sphere()] {
            let mut s = Sheaf::constant(cx, 2);
            let mut st = seed | 1;
            for m in s.maps[0].iter_mut() {
                *m = Mat::from_fn(2, 2, |_, _| { st ^= st << 13; st ^= st >> 7; st ^= st << 17; ((st % 2001) as f64) / 1000.0 - 1.0 });
            }
            let rejected = matches!(CochainComplex::from_sheaf(&s), Err(CohomError::NotAComplex { .. }));
            prop_assert!(rejected);
        }
    }
}
