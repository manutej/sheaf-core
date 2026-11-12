use sheaf_cell::{fixtures::*, num::linalg::*, num::sparse::*, sheaf::Sheaf};
use sheaf_cohom::{cohomology::*, fp, hodge_graph::*};

fn betti(s: &Sheaf) -> Vec<usize> {
    let cc = CochainComplex::from_sheaf(s).expect("complex");
    (0..cc.dims.len()).map(|k| { let h = cc.cohomology(k).unwrap(); assert!(h.consistent(), "two dim computations disagree"); h.dim }).collect()
}

#[test]
fn catalogue_matches_known_betti_numbers_over_r() {
    for n in catalogue() {
        let b = betti(&n.sheaf);
        assert_eq!(b, n.expect_r, "{}", n.name);
        let cc = CochainComplex::from_sheaf(&n.sheaf).unwrap();
        let chi_c: i64 = cc.dims.iter().enumerate().map(|(k, &d)| if k % 2 == 0 { d as i64 } else { -(d as i64) }).sum();
        let chi_h: i64 = b.iter().enumerate().map(|(k, &d)| if k % 2 == 0 { d as i64 } else { -(d as i64) }).sum();
        assert_eq!(chi_c, chi_h, "Euler characteristic {}", n.name);
    }
}

#[test]
fn filling_the_face_kills_h1() {
    assert_eq!(betti(&Sheaf::constant(triangle_graph(), 1)), vec![1, 1]);
    assert_eq!(betti(&Sheaf::constant(triangle_filled(), 1)), vec![1, 0, 0]);
}

#[test]
fn field_dependence_rp2() {
    let cc = CochainComplex::from_sheaf(&Sheaf::constant(rp2(), 1)).unwrap();
    assert_eq!(fp::betti_mod_p(&cc.delta, &cc.dims, 2).unwrap(), vec![1, 1, 1]);
    assert_eq!(fp::betti_mod_p(&cc.delta, &cc.dims, 3).unwrap(), vec![1, 0, 0]);
    let cc = CochainComplex::from_sheaf(&Sheaf::constant(torus7(), 1)).unwrap();
    assert_eq!(fp::betti_mod_p(&cc.delta, &cc.dims, 2).unwrap(), vec![1, 2, 1]);
}

#[test]
fn twisted_map_on_filled_face_is_rejected_as_not_a_complex() {
    let s = with_twist(triangle_filled(), 0, 1, rotation(0.5));
    match CochainComplex::from_sheaf(&s) {
        Err(CohomError::NotAComplex { k: 0, rel_residual }) => assert!(rel_residual > 1e-3),
        other => panic!("expected NotAComplex, got {:?}", other.map(|c| c.d2_residual)),
    }
}

#[test]
fn translation_loop_drift_is_purely_harmonic_and_exact_data_has_none() {
    // EN -> ES -> FR -> EN as a 3-cycle, 1-dim "meaning offset" stalk, expected maps = identity.
    // Each hop drifts by +0.1: every edge looks small, the loop sum is 0.3.
    let cc = CochainComplex::from_sheaf(&Sheaf::constant(cycle_graph(3), 1)).unwrap();
    // edge order from builder: (0,1),(0,2),(1,2); orientation low->high, so loop 0->1->2->0 = e01 + e12 - e02
    let w = vec![0.1, -0.1, 0.1];
    let h = cc.hodge(1, &w);
    let (ex, ha, co) = h.energies();
    assert!(h.cert.max() < 1e-12, "{:?}", h.cert);
    assert!(ex < 1e-24 && co == 0.0);
    assert!((ha - 0.03).abs() < 1e-12); // circulation 0.3 spread over 3 edges: 3 * 0.1^2
    // pure node-explainable data: w = d0 x
    let x = vec![0.0, 0.4, 1.0];
    let w2 = matvec(&cc.delta[0], &x);
    let h2 = cc.hodge(1, &w2);
    assert!(h2.energies().1 < 1e-24);
    // potential recovers x up to the constant (H^0) direction
    let shift = x[0] - h2.potential[0];
    for (i, &xi) in x.iter().enumerate() { assert!((h2.potential[i] + shift - xi).abs() < 1e-10); }
}

proptest::proptest! {
    /// Node-explainable data on a cycle: for every node vector x, `hodge(1, d0 x)` has no harmonic part and its
    /// potential equals x up to the constant (H^0) direction, on every node.
    #[test]
    fn hodge_potential_recovers_node_data_up_to_constant(x in proptest::collection::vec(-8.0f64..8.0, 3..=6)) {
        let cc = CochainComplex::from_sheaf(&Sheaf::constant(cycle_graph(x.len()), 1)).unwrap();
        let h = cc.hodge(1, &matvec(&cc.delta[0], &x));
        proptest::prop_assert!(h.energies().1 < 1e-20, "harmonic energy {}", h.energies().1);
        proptest::prop_assert_eq!(h.potential.len(), x.len());
        let shift = x[0] - h.potential[0];
        for (i, &xi) in x.iter().enumerate() { proptest::prop_assert!((h.potential[i] + shift - xi).abs() < 1e-9, "node {}", i); }
    }
}

#[test]
fn holonomy_kills_h0_and_k4_keeps_loop_space() {
    let s = catalogue().into_iter().find(|n| n.name == "k4_rot60").unwrap().sheaf;
    let cc = CochainComplex::from_sheaf(&s).unwrap();
    let h0 = cc.cohomology(0).unwrap();
    assert_eq!(h0.dim, 0);
    // a global-looking assignment (all vertices equal) now has nonzero discord on the twisted edge
    let x: Vec<f64> = (0..4).flat_map(|_| [1.0, 0.0]).collect();
    let w = matvec(&cc.delta[0], &x);
    assert!(norm(&w) > 0.5);
    // a random edge cochain has a nonzero harmonic (loop-level) part
    let w: Vec<f64> = (0..12).map(|i| ((i * 7919 % 13) as f64) / 13.0 - 0.5).collect();
    let hd = cc.hodge(1, &w);
    assert!(hd.energies().1 > 1e-3 && hd.cert.max() < 1e-10);
}

#[test]
fn sparse_cg_matches_dense_and_scales() {
    let s = Sheaf::constant(cycle_graph(50), 1);
    let cc = CochainComplex::from_sheaf(&s).unwrap();
    let w: Vec<f64> = (0..50).map(|i| (i as f64 * 0.37).sin()).collect();
    let dense = cc.hodge(1, &w);
    let (r, c, t) = s.coboundary_triplets(0);
    let g = graph_hodge_cg(&Csr::from_triplets(r, c, t), &w, 1e-12, 10_000);
    for i in 0..50 { assert!((g.harmonic[i] - dense.harmonic[i]).abs() < 1e-8); }
    assert!(g.cert_coclosed < 1e-9 && g.converged && !g.breakdown);
    // max_iter = 0 must not be reported as converged (Grok r1 #2)
    let (r, c, t) = s.coboundary_triplets(0);
    let g0 = graph_hodge_cg(&Csr::from_triplets(r, c, t), &w, 1e-12, 0);
    assert!(!g0.converged);
}

// ---- ADR-015: one H0 policy (per-component ADR-012 cut) for cohomology(0), hodge(0/1) and the spectral path ----
#[test]
fn h0_policy_is_per_component_in_cohomology0() {
    use faer::Mat;
    use sheaf_cell::num::tol::H0Policy;
    // Grok v1 #1 repro: a-b at scale 1e5, c-d at scale 1e-1 -> H0 = 2 (one constant per component)
    let m = |x: f64| Mat::from_fn(1, 1, |_, _| x);
    let s = Sheaf::from_graph(&[1, 1, 1, 1], &[(0, 1, m(1e5), m(1e5)), (2, 3, m(0.1), m(0.1))]).unwrap();
    let cc = CochainComplex::from_sheaf(&s).unwrap();
    assert_eq!(cc.n_comp, 2);
    let h0 = cc.cohomology(0).unwrap();
    assert_eq!(h0.dim, 2); assert!(h0.consistent());
    let h1 = cc.cohomology(1).unwrap();
    assert_eq!(h1.dim, 0);
    let sp = cc.d0_split();
    assert_eq!(h0.dim, 4 - sp.rank.rank);
    assert_eq!(sp.cuts.len(), 2); assert!(sp.cuts[0] > sp.cuts[1]);
    // absolute permitted delta 0.03 > lambda(c-d) = 0.02: that component's edge now counts as a section direction
    let cc2 = CochainComplex::from_sheaf(&s).unwrap().with_h0_policy(H0Policy::absolute(0.03));
    assert_eq!(cc2.cohomology(0).unwrap().dim, 3);
    assert_eq!(cc2.cohomology(1).unwrap().dim, 1, "H1 follows the same rank(d0) (rank-nullity holds)");
    // hodge(1) on the same policy: w on the 0.1-scale edge is exact under the default policy
    let hd = cc.hodge(1, &[0.0, 1.0]);
    assert!(hd.cert.max() < 1e-10); let (e, h, _) = hd.energies(); assert!((e - 1.0).abs() < 1e-12 && h < 1e-20);
}

#[test]
fn h0_policy_matches_spectral_path_on_fixtures() {
    use sheaf_cell::num::sparse::Csr; use sheaf_cohom::spectral::*;
    for f in catalogue() {
        let cc = CochainComplex::from_sheaf(&f.sheaf).unwrap();
        if cc.n_comp != 1 || cc.delta.is_empty() { continue; } // smallest() applies the policy to the whole operator
        let (r, c, t) = f.sheaf.coboundary_triplets(0);
        let op = HodgeOp::l0(&Csr::from_triplets(r, c, t));
        let sp = smallest_dense(&op, &KnownKernel::none(op.n), op.n, H0Policy::adr012());
        assert_eq!(sp.h0_dim(), cc.cohomology(0).unwrap().dim, "{}", f.name);
    }
}
