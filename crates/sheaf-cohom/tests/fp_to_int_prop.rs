//! `fp::to_int` admission law: a matrix is converted iff every entry is a finite integer (to 1e-12) within ±2^53,
//! and then each output entry is that integer. The reference predicate is the pre-refactor rejection test, verbatim.
use faer::Mat;
use proptest::prelude::*;
use sheaf_cohom::fp::{to_int, MAX_EXACT_INT};

/// The rejection test `to_int` used before it was rewritten as a positive conjunction.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn old_reject(v: f64) -> bool { !v.is_finite() || v.abs() > MAX_EXACT_INT || !((v - v.round()).abs() <= 1e-12) }

fn entry() -> impl Strategy<Value = f64> {
    prop_oneof![
        any::<u64>().prop_map(f64::from_bits), // every bit pattern: all NaNs, ±inf, subnormals, huge
        (-(1i64 << 53)..=(1i64 << 53)).prop_map(|k| k as f64), // exact integers up to the bound
        (-2000i32..=2000, -3e-12f64..3e-12).prop_map(|(k, d)| k as f64 + d), // straddles the 1e-12 tolerance
        prop::sample::select(vec![
            f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.0, 0.5, 1e-12, 1.1e-12,
            MAX_EXACT_INT, -MAX_EXACT_INT, 9007199254740994.0, -9007199254740994.0,
        ]),
    ]
}

proptest! {
    #[test]
    fn to_int_accepts_exactly_the_exact_integer_matrices(
        (nrows, ncols, vals) in (0usize..4, 0usize..4, prop::collection::vec(entry(), 9)),
    ) {
        let a = Mat::from_fn(nrows, ncols, |i, j| vals[i * 3 + j]);
        let rejected = (0..nrows).any(|i| (0..ncols).any(|j| old_reject(a[(i, j)])));
        match to_int(&a) {
            None => prop_assert!(rejected),
            Some(m) => {
                prop_assert!(!rejected);
                prop_assert_eq!(m.len(), nrows);
                for (i, row) in m.iter().enumerate() {
                    prop_assert_eq!(row.len(), ncols);
                    for (j, &x) in row.iter().enumerate() { prop_assert_eq!(x as f64, a[(i, j)].round()); }
                }
            }
        }
    }
}
