//! Breakdown guard of `krylov::cg` on a 1x1 system.
//! With A = [a], b = [1], x0 = [0], scale = 1 and max_iter = 1 the loop body runs once with p = [1], so pᵀAp is `a`
//! and the guard threshold is `1e-14 * op_scale`. CG takes the step iff `a > 1e-14 * op_scale`; a NaN on either side
//! is a breakdown.
use proptest::prelude::*;
use sheaf_cell::num::{dot, krylov::{cg, SolveTol}};

fn any_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        any::<u64>().prop_map(f64::from_bits),
        -1e3f64..1e3,
        Just(f64::NAN),
        Just(f64::INFINITY),
        Just(f64::NEG_INFINITY),
        Just(0.0),
        Just(-0.0),
    ]
}

/// `(a, op_scale)`: either both free, or `a` placed below, on and above the threshold of a free `op_scale`.
fn case() -> impl Strategy<Value = (f64, f64)> {
    prop_oneof![
        (any_f64(), any_f64()),
        (any_f64(), prop_oneof![Just(0.5), Just(1.0), Just(2.0)]).prop_map(|(s, k)| (1e-14 * s * k, s)),
    ]
}

proptest! {
    #[test]
    fn cg_steps_iff_curvature_exceeds_threshold((a, op_scale) in case()) {
        // A·0 is 0 for every a, so the start residual is b even when a is infinite or NaN.
        let apply = |v: &[f64]| vec![if v[0] == 0.0 { 0.0 } else { a * v[0] }];
        let out = cg(&apply, &[1.0], vec![0.0], 1.0, op_scale, SolveTol { rtol: 1e-12, max_iter: 1 });
        let advances = a > 1e-14 * op_scale;
        prop_assert_eq!(out.breakdown, !advances);
        if advances {
            prop_assert_eq!(out.iters, 1);
            prop_assert_eq!(out.x[0].to_bits(), (1.0 / dot(&[1.0], &[a])).to_bits());
        } else {
            prop_assert_eq!(out.iters, 0);
            prop_assert_eq!(out.x[0].to_bits(), 0.0f64.to_bits());
            prop_assert_eq!(out.residual_norm, 1.0);
        }
    }
}
