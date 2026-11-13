//! Curvature guard of `graph_hodge_cg_cert` on a 1x1 coboundary d0 = [a] with w = [w0].
//! `rtol = -1.0` switches the residual stop off, so the loop is entered whenever rr is not NaN.
//! The first direction is p = a*w0, with pᵀL0p = p*(a*(a*p)) and pᵀp = p*p in f64 arithmetic.
//! CG takes the step iff `pᵀL0p > 1e-14 * pᵀp`; equality and a NaN on either side are a breakdown.
use proptest::prelude::*;
use sheaf_cell::num::sparse::Csr;
use sheaf_cohom::hodge_graph::graph_hodge_cg_cert;

fn any_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        any::<u64>().prop_map(f64::from_bits),
        -1e3f64..1e3,
        Just(f64::NAN),
        Just(f64::INFINITY),
        Just(0.0),
        Just(-0.0),
    ]
}

/// `(a, w0)`: either both free, or `a` below, near and above the threshold |a| = 1e-7 with a moderate `w0`.
fn case() -> impl Strategy<Value = (f64, f64)> {
    prop_oneof![
        (any_f64(), any_f64()),
        (prop_oneof![Just(0.5e-7), Just(1e-7), Just(2e-7)], -1e3f64..1e3),
    ]
}

fn d0(a: f64) -> Csr { Csr::from_triplets(1, 1, vec![(0, 0, a)]) }

proptest! {
    #[test]
    fn first_step_taken_iff_curvature_exceeds_threshold((a, w0) in case()) {
        let out = graph_hodge_cg_cert(&d0(a), &[w0], -1.0, 1, 1e-8);
        let p = a * w0;
        let (pap, pp) = (p * (a * (a * p)), p * p);
        let entered = !pp.is_nan();
        let advances = entered && pap > 1e-14 * pp;
        prop_assert_eq!(out.breakdown, entered && !advances);
        prop_assert_eq!(out.iters, usize::from(advances));
        let x = if advances { 0.0 + (pp / pap) * p } else { 0.0 };
        prop_assert_eq!(out.potential[0].to_bits(), x.to_bits());
        prop_assert!(!(out.breakdown && out.converged));
    }

    /// a = ±2^j and w0 = ±2^-(e+j) give p = ±2^-e, so pᵀp underflows to 0 while pᵀL0p = 2^(2j-2e) > 0.
    /// The first step is taken with alpha = 0 and beta = 0/0, so the second direction is NaN: a breakdown.
    #[test]
    fn nan_direction_is_a_breakdown(j in 100i32..=200, e in 540i32..=600, sa in any::<bool>(), sw in any::<bool>()) {
        let sign = |neg: bool| if neg { -1.0 } else { 1.0 };
        let (a, w0) = (sign(sa) * 2f64.powi(j), sign(sw) * 2f64.powi(-(e + j)));
        let out = graph_hodge_cg_cert(&d0(a), &[w0], -1.0, 2, 1e-8);
        prop_assert!(out.breakdown);
        prop_assert_eq!(out.iters, 1);
        prop_assert_eq!(out.potential[0].to_bits(), 0.0f64.to_bits());
        prop_assert!(!out.converged);
    }
}
