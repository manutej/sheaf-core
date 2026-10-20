//! Relative H⁰(X,B;F) (uniqueness certificate for harmonic extension), nearest global section, typed F_p Betti numbers.
use faer::Mat;
use sheaf_cell::{complex::CellComplex, fixtures::*, num::{sparse::Csr, tol::H0Policy}, sheaf::Sheaf};
use sheaf_cohom::{cohomology::CochainComplex, fp, relative::{h0_relative, RelativeError}};

fn d0(s: &Sheaf) -> Csr { let (r, c, t) = s.coboundary_triplets(0); Csr::from_triplets(r, c, t) }
fn path(n: usize) -> CellComplex { let e: Vec<_> = (0..n - 1).map(|i| (i, i + 1)).collect(); CellComplex::simplicial(n, &e, &[]) }
fn rel(s: &Sheaf, b: &[usize]) -> usize { h0_relative(&d0(s), b, H0Policy::adr012(), 1500).unwrap().value.dim }

#[test]
fn empty_boundary_equals_absolute_h0_on_the_catalogue() {
    for nm in catalogue() {
        let cx = CochainComplex::from_sheaf(&nm.sheaf).unwrap();
        assert_eq!(rel(&nm.sheaf, &[]), cx.cohomology(0).unwrap().dim, "{}", nm.name);
    }
}

#[test]
fn pinning_kills_sections_and_detects_unpinned_parts() {
    // constant R^2 on a path: pin one end -> unique (H0(X,B) = 0)
    let s = Sheaf::constant(path(6), 2);
    assert_eq!(rel(&s, &[]), 2); assert_eq!(rel(&s, &[0, 1]), 0); assert_eq!(rel(&s, &[0]), 1);
    // twisted cycle (no global sections) is unique even with B = ∅
    let t = with_twist(Sheaf::constant(cycle_graph(5), 2).complex, 0, 0, rotation(0.9));
    assert_eq!(rel(&t, &[]), 0);
    // two components, boundary only on the first: the second keeps its constants
    let cx = CellComplex::graph(6, &[(0, 1), (1, 2), (3, 4), (4, 5)]).unwrap();
    let s = Sheaf::constant(cx, 1);
    let r = h0_relative(&d0(&s), &[0], H0Policy::adr012(), 1500).unwrap().value;
    assert_eq!((r.dim, r.unique()), (1, false));
    for i in 0..3 { assert!(r.basis[(i, 0)].abs() < 1e-12); }
    assert_eq!(h0_relative(&d0(&s), &[9], H0Policy::adr012(), 1500).unwrap_err(), RelativeError::BoundaryIndex(9));
    assert!(matches!(h0_relative(&d0(&s), &[], H0Policy::adr012(), 2), Err(RelativeError::BlockTooLarge { .. })));
}

#[test]
fn nearest_section_is_the_h0_projection() {
    let s = Sheaf::constant(cycle_graph(7), 2);
    let cx = CochainComplex::from_sheaf(&s).unwrap();
    let x: Vec<f64> = (0..14).map(|i| (i as f64 * 0.37).sin()).collect();
    let ns = cx.nearest_section(&x);
    let h = cx.hodge(0, &x);
    for i in 0..14 { assert!((ns.value.section[i] - h.harmonic[i]).abs() < 1e-12); }
    // constants: the section is the per-coordinate mean
    for k in 0..2 { let m = (0..7).map(|v| x[2 * v + k]).sum::<f64>() / 7.0; for v in 0..7 { assert!((ns.value.section[2 * v + k] - m).abs() < 1e-12); } }
    assert!(ns.value.section_residual < 1e-14 && ns.value.h0_dim == 2);
}

#[test]
fn betti_fp_is_typed_exact_fp() {
    let s = Sheaf::constant(rp2(), 1);
    let cx = CochainComplex::from_sheaf(&s).unwrap();
    let b2 = fp::betti_fp(&cx.delta, &cx.dims, 2).unwrap();
    assert_eq!((b2.value.clone(), b2.p()), (vec![1, 1, 1], 2));
    let _: Mat<f64> = cx.delta[0].clone();
}
